use super::process_optimization::Gradients;
use crate::process::ProcessModel;
use crate::{MolecularRepresentation, OptimizationProblem, OptimizationResult, PropertyModel};
use feos::core::{EosResult, IdealGas, Residual};
use knitro_rs::{Knitro, KnitroError};
use std::f64;

#[derive(Clone, Copy)]
pub enum OuterApproximationAlgorithm {
    DuranGrossmann(bool),
    FletcherLeyffer,
}

impl<
        E: Residual + IdealGas,
        M: MolecularRepresentation,
        R: PropertyModel<EquationOfState = E>,
        P: ProcessModel<E>,
    > OptimizationProblem<E, M, R, P>
{
    pub fn solve_outer_approximation(
        &mut self,
        y_fixed: Vec<f64>,
        algorithm: OuterApproximationAlgorithm,
        options_nlp: Option<&str>,
        options_mip: Option<&str>,
    ) -> EosResult<OptimizationResult> {
        match algorithm {
            OuterApproximationAlgorithm::DuranGrossmann(update_lower_bound) => self
                .solve_outer_approximation_duran(
                    y_fixed.to_vec(),
                    update_lower_bound,
                    options_nlp,
                    options_mip,
                ),
            OuterApproximationAlgorithm::FletcherLeyffer => {
                self.solve_outer_approximation_fletcher(y_fixed.to_vec(), options_nlp, options_mip)
            }
        }
    }

    fn solve_outer_approximation_duran(
        &mut self,
        mut y_fixed: Vec<f64>,
        update_lower_bound: bool,
        options_nlp: Option<&str>,
        options_mip: Option<&str>,
    ) -> EosResult<OptimizationResult> {
        let mut f_l = f64::MIN;
        let mut f_u = f64::MAX;

        // Evaluate NLP for initial guess
        let mut res = self.solve_fixed(&y_fixed, options_nlp).unwrap();
        println!("{:8.5} {:.5?} {:?}", res.target, res.x, res.smiles);
        for _ in 0..500 {
            f_u = f_u.min(res.target);

            // Evaluate gradients
            let Ok((vars, (f, g, h), lambda)) = self.solve_gradients(&y_fixed, &res.x, options_nlp)
            else {
                return Ok(res);
            };

            // Solve master MILP
            let Ok((f, y)) = self.solve_master(vars, (f, g, h), lambda, f_l, f_u, options_mip)
            else {
                return Ok(res);
            };
            y_fixed = y;
            if update_lower_bound {
                f_l = f;
            }

            // Solve NLP subproblem
            let r = self.solve_fixed(&y_fixed, options_nlp);
            if let Ok(r) = r {
                res = r;
                println!("{:8.5} {:.5?} {:?}", res.target, res.x, res.smiles);
            } else {
                let y_usize: Vec<_> = y_fixed.iter().map(|y| y.round() as usize).collect();
                let smiles = self.molecules.smiles(&y_usize);
                println!("{:?} not converged!", smiles);
                self.solutions.insert(OptimizationResult::new(
                    0.0,
                    smiles,
                    vec![],
                    y_usize,
                    vec![],
                ));
                return Ok(res);
            }
        }
        Err(feos::core::EosError::NotConverged(
            "Outer approximation".into(),
        ))
    }

    fn solve_outer_approximation_fletcher(
        &mut self,
        mut y_fixed: Vec<f64>,
        options_nlp: Option<&str>,
        options_mip: Option<&str>,
    ) -> EosResult<OptimizationResult> {
        let mut ubd = f64::MAX;
        let mut solutions = vec![];

        for _ in 0..500 {
            // Solve NLP subproblem
            let res = self.solve_fixed(&y_fixed, options_nlp);
            let (res, feasible) = if let Ok(res) = res {
                println!("{:8.5} {:.5?} {:?}", res.target, res.x, res.smiles);
                ubd = ubd.min(res.target);
                (res, true)
            } else {
                return Err(feos::core::EosError::NotConverged("NLP subproblem".into()));
            };

            // Evaluate gradients
            let Ok((vars, gradients, _)) = self.solve_gradients(&y_fixed, &res.x, options_nlp)
            else {
                return Ok(res);
            };
            solutions.push((vars, gradients, feasible));

            // Solve master MILP
            let Ok((_, y)) = self.solve_master_fletcher(&solutions, ubd, options_mip) else {
                return Ok(res);
            };
            y_fixed = y;
        }
        Err(feos::core::EosError::NotConverged(
            "Outer approximation".into(),
        ))
    }

    pub fn outer_approximation_ranking(
        &mut self,
        y_fixed: &[f64],
        algorithm: OuterApproximationAlgorithm,
        runs: usize,
        options_nlp: Option<&str>,
        options_mip: Option<&str>,
    ) {
        for k in 0..runs {
            let result = self.solve_outer_approximation(
                y_fixed.to_vec(),
                algorithm,
                options_nlp,
                options_mip,
            );
            match result {
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
                        println!(
                            "{k} {:10.7} {:.5?} {:?}",
                            solution.target, solution.x, solution.smiles
                        );
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

        // Declare auxiliary variable
        let mu = kc.add_vars(1)?[0];
        kc.set_var_lobnd(mu, z_l)?;
        kc.set_var_upbnd(mu, z_u)?;

        // Declare process variables
        let x = self.process.variables().setup_knitro(&kc, None)?;

        let (y, n) = self
            .molecules
            .setup_knitro(&kc, None, &self.solutions, false)?;

        let p = self.property_model.setup_knitro(&kc, &n)?;

        let vars = [x, p.clone()].concat();

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
        for g in g.into_iter() {
            let con = kc.add_con()?;
            let (coefs, bnd) = g.linearize(&x_k);
            kc.add_con_linear_struct_one(con, &vars, &coefs)?;
            kc.set_con_lobnd(con, bnd)?;
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

    fn solve_master_fletcher(
        &self,
        solutions: &[(Vec<f64>, Gradients, bool)],
        ubd: f64,
        options: Option<&str>,
    ) -> Result<(f64, Vec<f64>), KnitroError> {
        // set up MILP
        let kc = Knitro::new()?;

        // Declare auxiliary variable
        let eta = kc.add_vars(1)?[0];
        kc.set_var_upbnd(eta, ubd)?;

        // Declare process variables
        let x = self.process.variables().setup_knitro(&kc, None)?;

        let (y, n) = self
            .molecules
            .setup_knitro(&kc, None, &self.solutions, false)?;

        let p = self.property_model.setup_knitro(&kc, &n)?;

        let vars = [x, p.clone()].concat();

        for (x_i, (f, g, h), feasible) in solutions {
            // Linearized target function
            if *feasible {
                let con = kc.add_con()?;
                let (coefs, bnd) = f.clone().linearize(x_i);
                kc.add_con_linear_struct_one(con, &vars, &coefs)?;
                kc.add_con_linear_term(con, eta, -1.0)?;
                kc.set_con_upbnd(con, bnd)?;
            }

            // Linearized equality constraints
            for h in h {
                let con = kc.add_con()?;
                let (coefs, bnd) = h.clone().linearize(x_i);
                kc.add_con_linear_struct_one(con, &vars, &coefs)?;
                kc.set_con_upbnd(con, bnd)?;
                kc.set_con_lobnd(con, bnd)?;
            }

            // Linearized inequality constraints
            for g in g {
                let con = kc.add_con()?;
                let (coefs, bnd) = g.clone().linearize(x_i);
                kc.add_con_linear_struct_one(con, &vars, &coefs)?;
                kc.set_con_lobnd(con, bnd)?;
            }
        }

        // Target
        kc.add_obj_linear_struct(&[eta], &[1.0])?;
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
