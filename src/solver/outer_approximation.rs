use super::Gradients;
use crate::process::ProcessModel;
use crate::{
    MolecularRepresentation, OptimizationMode, OptimizationProblem, OptimizationResult,
    PropertyModel,
};
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
    ) -> Result<(), KnitroError> {
        let mut f_l = f64::MIN;
        let mut f_u = f64::MAX;

        let mut res = self.solve_fixed(None, &y_fixed, options_nlp)?;
        let y: Vec<_> = y_fixed.iter().map(|y| *y as usize).collect();
        let smiles = &self.molecules.smiles(&y)[0];
        println!("{:8.5} {:?} {smiles}", res.0, y);
        for _ in 0..500 {
            let (target, [ref x, ref y, ref p]) = res;
            f_u = f_u.min(target);
            let mut gradients = Default::default();
            self.solve(
                Some(x),
                Some(y),
                Some(p),
                options_nlp,
                OptimizationMode::Gradients,
                Some(&mut gradients),
            )?;
            let (x, (f, g, h), lambda) = gradients;
            let y: Vec<_> = y.iter().map(|y| *y as usize).collect();
            self.solutions.insert(OptimizationResult::new(
                target,
                vec![],
                x.clone(),
                y.clone(),
                p.clone(),
            ));

            let (f, y) = self.solve_master(x, (f, g, h), lambda, f_l, f_u, options_mip)?;
            y_fixed = y;
            if update_lower_bound {
                f_l = f;
            }

            let r = self.solve_fixed(None, &y_fixed, options_nlp);
            if let Ok(r) = r {
                res = r;
                let (target, [x, y, p]) = res.clone();
                let y: Vec<_> = y.iter().map(|y| *y as usize).collect();
                let smiles = &self.molecules.smiles(&y)[0];
                println!("{:8.5} {:?} {smiles}", res.0, y);

                // println!("{} {:?}", res.1.val, res.0);
            } else {
                println!("{:8.5}", 0.0);
            }
        }
        Ok(())
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

        // // declare composition variables
        // let c = kc.add_vars(self.comps)?;
        // for &i in c.iter() {
        //     kc.set_var_lobnd(i, 0.01)?;
        // }
        // let con = kc.add_con()?;
        // kc.add_con_linear_struct_one(con, &c, &vec![1.0; x.len()])?;
        // kc.set_con_eqbnd(con, 1.0)?;

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
