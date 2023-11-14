use crate::process::ProcessModel;
use crate::*;
use feos_core::EosResult;
use knitro_rs::*;

impl<
        M: MolecularRepresentation<N>,
        R: PropertyModel<M::ChemicalRecord>,
        P: ProcessModel,
        const N: usize,
    > EvalCallback for OptimizationProblem<M, R, P, N>
{
    fn callback(&self, x: &[f64], c: &mut [f64]) -> f64 {
        match self.evaluate(x) {
            Ok((target, constraints)) => {
                c.copy_from_slice(&constraints);
                target
            }
            Err(_) => {
                // println!("{e}");
                c.iter_mut().for_each(|c| *c = -1.0);
                0.0
            }
        }
    }
}

impl<
        M: MolecularRepresentation<N>,
        R: PropertyModel<M::ChemicalRecord>,
        P: ProcessModel,
        const N: usize,
    > OptimizationProblem<M, R, P, N>
{
    fn evaluate(&self, x: &[f64]) -> EosResult<(f64, Vec<f64>)> {
        let n_y = self.molecules.variables().len();
        let (y, x) = x.split_at(n_y);
        let cr = self.molecules.build(y);
        let eos = self.property_model.build_eos(cr.into())?;
        let (_, target, constraints) = self.process.solve(&eos, x)?;
        Ok((target, constraints))
    }

    fn setup_knitro_molecule(
        &self,
        kc: &Knitro,
        y0: Option<&[f64]>,
        vartype: u32,
    ) -> Result<usize, KnitroError> {
        // define variables
        let variables = self.molecules.variables();
        let n_y = variables.len();
        let index_vars = kc.add_vars(variables.len(), Some(vartype))?;
        for (&i, variable) in index_vars.iter().zip(variables) {
            kc.set_var_primal_initial_value(i, variable.init)?;
            kc.set_var_lobnd(i, variable.lobnd)?;
            kc.set_var_upbnd(i, variable.upbnd)?;
        }
        if let Some(y0) = y0 {
            kc.set_var_primal_initial_values(&index_vars, y0)?;
        }

        // define constraints
        let constraints = self.molecules.constraints(&index_vars);
        for constraint in constraints {
            let c = kc.add_con()?;
            kc.add_con_linear_struct_one(c, &constraint.vars, &constraint.coefs)?;
            if let Some(lobnd) = constraint.lobnd {
                kc.set_con_lobnd(c, lobnd)?;
            }
            if let Some(upbnd) = constraint.upbnd {
                kc.set_con_upbnd(c, upbnd)?;
            }
        }

        // integer cuts
        for solution in self.solutions.iter() {
            let y0 = &solution.y;
            let vars = Vec::from_iter(0..y0.len() as i32);
            let c = kc.add_con()?;
            let qcoefs: Vec<_> = y0.iter().map(|_| 1.0).collect();
            let lcoefs: Vec<_> = y0.iter().map(|&n0| -2.0 * n0 as f64).collect();
            let lbond = 1.0 - y0.iter().map(|n0| n0.pow(2) as f64).sum::<f64>();
            kc.add_con_quadratic_struct_one(c, &vars, &vars, &qcoefs)?;
            kc.add_con_linear_struct_one(c, &vars, &lcoefs)?;
            kc.set_con_lobnd(c, lbond)?;
        }
        Ok(n_y)
    }

    fn setup_knitro_process(
        &self,
        kc: &Knitro,
        x0: &[f64],
        vartype: u32,
    ) -> Result<Vec<i32>, KnitroError> {
        // declare continuous variables
        let variables = self.process.variables();
        let index_vars = kc.add_vars(variables.len(), None)?;
        for (&i, [l, u]) in index_vars.iter().zip(variables.into_iter()) {
            if let Some(l) = l {
                kc.set_var_lobnd(i, l)?;
            }
            if let Some(u) = u {
                kc.set_var_upbnd(i, u)?;
            }
        }
        kc.set_var_primal_initial_values(&index_vars, x0)?;

        // declare binary variables
        if self.process.binary_variables() > 0 {
            kc.add_vars(self.process.binary_variables(), Some(vartype))?;
        }

        // add constraints
        let constraints = self.process.constraints();
        let index_cons = kc.add_cons(constraints.len())?;
        for (&i, [l, u, e]) in index_cons.iter().zip(constraints.into_iter()) {
            if let Some(l) = l {
                kc.set_con_lobnd(i, l)?;
            }
            if let Some(u) = u {
                kc.set_con_upbnd(i, u)?;
            }
            if let Some(e) = e {
                kc.set_con_eqbnd(i, e)?;
            }
        }

        Ok(index_cons)
    }

    pub fn solve_target(
        &mut self,
        x0: &[f64],
        options: Option<&str>,
    ) -> Result<(f64, Vec<f64>, Vec<f64>), KnitroError> {
        let kc = Knitro::new()?;

        // Set up CAMD formulation
        let n_y = self.setup_knitro_molecule(&kc, None, KN_VARTYPE_CONTINUOUS)?;

        // Set up process variables and constraints
        let index_cons = self.setup_knitro_process(&kc, x0, KN_VARTYPE_CONTINUOUS)?;

        // Set up function callback
        let cb = kc.add_eval_callback(true, &index_cons, self)?;
        if let Some(options) = options {
            kc.load_param_file(options)?;
        }
        kc.set_cb_user_params(cb, self)?;

        // Solve MINLP
        kc.solve()?;
        let t = kc.get_obj_value()?;
        let x = kc.get_var_primal_values_all()?;
        let (y, x) = x.split_at(n_y);

        Ok((t, x.to_vec(), y.to_vec()))
    }

    pub fn solve_knitro_once(
        &mut self,
        x0: &[f64],
        y0: Option<&[f64]>,
        options: Option<&str>,
    ) -> Result<OptimizationResult, KnitroError> {
        let kc = Knitro::new()?;

        // Set up CAMD formulation
        let n_y = self.setup_knitro_molecule(&kc, y0, KN_VARTYPE_BINARY)?;

        // Set up process variables and constraints
        let index_cons = self.setup_knitro_process(&kc, x0, KN_VARTYPE_BINARY)?;

        // Set up function callback
        let cb = kc.add_eval_callback(true, &index_cons, self)?;
        if let Some(options) = options {
            kc.load_param_file(options)?;
        }
        kc.set_cb_user_params(cb, self)?;

        // Solve MINLP
        kc.solve()?;
        let t = kc.get_obj_value()?;
        let x = kc.get_var_primal_values_all()?;
        let (y, x) = x.split_at(n_y);
        let y: Vec<_> = y.iter().map(|y| *y as usize).collect();
        let smiles = self.molecules.smiles(&y);
        // println!("{t} {y:?} {x:?}");
        Ok(OptimizationResult::new(t, smiles.into(), y, x.to_vec()))
    }

    pub fn solve_knitro(
        &mut self,
        x0: &[f64],
        y0: Option<&[f64]>,
        n_solutions: usize,
        options: Option<&str>,
    ) {
        for k in 0..n_solutions {
            if let Ok(result) = self.solve_knitro_once(x0, y0, options) {
                self.solutions.insert(result.clone());
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
                    println!(
                        "{k} {:.7} {:?} {smiles} {:?}",
                        solution.target, solution.y, solution.x
                    );
                }
            }
        }
    }
}
