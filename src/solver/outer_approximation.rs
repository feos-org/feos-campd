use super::Gradients;
use crate::process::ProcessModel;
use crate::{
    MolecularRepresentation, OptimizationMode, OptimizationProblem, OptimizationResult,
    PropertyModel,
};
use feos::core::EosResult;
use knitro_rs::{Knitro, KnitroError};

impl<
        M: MolecularRepresentation<N>,
        R: PropertyModel<M::ChemicalRecord>,
        P: ProcessModel,
        const N: usize,
    > OptimizationProblem<M, R, P, N>
{
    pub fn solve_outer_approximation(
        &mut self,
        mut y_fixed: Vec<f64>,
        update_lower_bound: bool,
        options_nlp: Option<&str>,
        options_mip: Option<&str>,
    ) -> EosResult<OptimizationResult> {
        let mut f_l = f64::MIN;
        let mut f_u = f64::MAX;

        // Evaluate NLP for initial guess
        let mut res = self.solve_fixed(None, &y_fixed, options_nlp).unwrap();
        let y: Vec<_> = y_fixed.iter().map(|y| *y as usize).collect();
        let smiles = &self.molecules.smiles(&y)[0];
        println!("{:8.5} {smiles}", res.0);
        for _ in 0..500 {
            let (target, [ref x, ref y, ref p]) = res;
            f_u = f_u.min(target);

            let y_usize: Vec<_> = y.iter().map(|&y| y as usize).collect();
            self.solutions.retain(|e| e.y != y_usize);

            // Evaluate gradients
            let mut gradients = Default::default();
            self.solve(
                Some(x),
                Some(y),
                Some(p),
                options_nlp,
                OptimizationMode::Gradients,
                Some(&mut gradients),
            )
            .unwrap();
            let (vars, (f, g, h), lambda) = gradients;
            let y: Vec<_> = y.iter().map(|y| *y as usize).collect();
            let result = OptimizationResult::new(
                target,
                self.molecules.smiles(&y).to_vec(),
                x.clone(),
                y.clone(),
                p.clone(),
            );
            self.solutions.insert(result.clone());

            // Solve master MILP
            let Ok((f, y)) = self.solve_master(vars, (f, g, h), lambda, f_l, f_u, options_mip)
            else {
                return Ok(result);
            };
            y_fixed = y;
            if update_lower_bound {
                f_l = f;
            }

            // Solve NLP subproblem
            let r = self.solve_fixed(None, &y_fixed, options_nlp);
            let y: Vec<_> = y_fixed.iter().map(|y| *y as usize).collect();
            let smiles = self.molecules.smiles(&y);
            if let Ok(r) = r {
                res = r;
                println!("{:8.5} {}", res.0, smiles[0]);
            } else {
                self.solutions.insert(OptimizationResult::new(
                    0.0,
                    smiles.to_vec(),
                    vec![],
                    y,
                    vec![],
                ));
                println!("{} not converged!", smiles[0])
            }
        }
        Err(feos::core::EosError::NotConverged(
            "Outer approximation".into(),
        ))
    }

    pub fn outer_approximation_ranking(
        &mut self,
        y_fixed: &[f64],
        update_lower_bound: bool,
        n_solutions: usize,
        options_nlp: Option<&str>,
        options_mip: Option<&str>,
    ) {
        let y_usize: Vec<_> = y_fixed.iter().map(|&y| y as usize).collect();
        for k in 0..n_solutions {
            self.solutions.retain(|e| e.y != y_usize);
            match self.solve_outer_approximation(
                y_fixed.to_vec(),
                update_lower_bound,
                options_nlp,
                options_mip,
            ) {
                Ok(result) => {
                    let mut solutions: Vec<_> = self.solutions.iter().collect();
                    solutions.sort_by(|s1, s2| s1.target.total_cmp(&s2.target));
                    println!("\nRun {}", k + 1);
                    for (k, solution) in solutions.into_iter().enumerate() {
                        let k = if solution == &result || solution.target < result.target {
                            format!("{:3}", k + 1)
                        } else {
                            "   ".into()
                        };
                        let smiles = match &solution.smiles[..] {
                            [smiles] => smiles.clone(),
                            _ => format!("[{}]", solution.smiles.join(", ")),
                        };
                        println!("{k} {:10.7} {:.5?} {smiles}", solution.target, solution.x);
                    }
                    println!();
                }
                Err(e) => println!("\nRun {}\n{e}", k + 1),
            }
        }
    }

    fn solve_master(
        &self,
        x_k: Vec<f64>,
        (f, g, h): Gradients,
        lambda: Vec<f64>,
        z_l: f64,
        z_u: f64,
        options: Option<&str>,
    ) -> Result<(f64, Vec<f64>), KnitroError> {
        // set up MILP
        let kc = Knitro::new()?;
        let mode = OptimizationMode::MolecularDesign;

        // Declare auxiliary variable
        let mu = kc.add_vars(1)?[0];
        kc.set_var_lobnd(mu, z_l)?;
        kc.set_var_upbnd(mu, z_u)?;

        // Declare process variables
        let x = self.process.variables().setup_knitro(&kc, None, mode)?;

        let [y, p] = self
            .molecules
            .setup_knitro(&kc, None, None, &self.solutions, mode)?;

        let vars = [x, p].concat();

        // Linearized target function
        let con = kc.add_con()?;
        let (coefs, bnd) = f.linearize(&x_k);
        kc.add_con_linear_struct_one(con, &vars, &coefs)?;
        kc.add_con_linear_term(con, mu, -1.0)?;
        kc.set_con_upbnd(con, bnd)?;

        // Linearized equality constraints
        for (h, &l) in h.into_iter().zip(lambda.iter()) {
            if l == 0.0 {
                continue;
            }
            let con = kc.add_con()?;
            let (coefs, bnd) = h.linearize(&x_k);
            kc.add_con_linear_struct_one(con, &vars, &coefs)?;
            if l > 0.0 {
                kc.set_con_upbnd(con, bnd)?;
            } else {
                kc.set_con_lobnd(con, bnd)?;
            }
        }

        // Linearized inequality constraints
        for (g, [lobnd, upbnd]) in g.into_iter().zip(self.process.constraints()) {
            let con = kc.add_con()?;
            let (coefs, bnd) = g.linearize(&x_k);
            kc.add_con_linear_struct_one(con, &vars, &coefs)?;
            if let Some(lobnd) = lobnd {
                kc.set_con_lobnd(con, lobnd + bnd)?;
            }
            if let Some(upbnd) = upbnd {
                kc.set_con_upbnd(con, upbnd + bnd)?;
            }
        }

        // Target
        kc.add_obj_linear_struct(&[mu], &[1.0])?;
        if let Some(options) = options {
            kc.load_param_file(options)?;
        }

        // Solve MILP
        kc.solve()?;
        let f = kc.get_obj_value()?;
        let y = kc.get_var_primal_values(&y)?;

        Ok((f, y.to_vec()))
    }
}
