use crate::process::ProcessModel;
use crate::*;
use feos::core::{EosResult, IdealGas, Residual};
use knitro_rs::*;

mod outer_approximation;
mod process_optimization;
pub use outer_approximation::OuterApproximationAlgorithm;

struct OptimizationProblemCallback<'a, R, P> {
    property_model: &'a R,
    process: &'a P,
    index_process_vars: &'a [i32],
    index_parameter_vars: &'a [i32],
}

impl<'a, R, P> OptimizationProblemCallback<'a, R, P> {
    fn new<E, M>(
        problem: &'a OptimizationProblem<E, M, R, P>,
        index_process_vars: &'a [i32],
        index_parameter_vars: &'a [i32],
    ) -> Self {
        Self {
            property_model: &problem.property_model,
            process: &problem.process,
            index_process_vars,
            index_parameter_vars,
        }
    }
}

impl<'a, E: Residual + IdealGas, R: PropertyModel<EquationOfState = E>, P: ProcessModel<E>>
    EvalCallback for OptimizationProblemCallback<'a, R, P>
{
    fn callback(&self, x: &[f64], obj: &mut f64, c: &mut [f64]) -> i32 {
        match self.evaluate(x) {
            Ok((target, constraints)) => {
                if constraints.len() != c.len() {
                    println!("Wrong number of constraints returned from process model!");
                    return KN_RC_EVAL_ERR;
                }
                c.copy_from_slice(&constraints);
                *obj = target;
                0
            }
            Err(_) => {
                // println!("{e}");
                // panic!();
                KN_RC_EVAL_ERR
            }
        }
    }
}

impl<'a, E: Residual + IdealGas, R: PropertyModel<EquationOfState = E>, P: ProcessModel<E>>
    OptimizationProblemCallback<'a, R, P>
{
    fn evaluate(&self, vars: &[f64]) -> EosResult<(f64, Vec<f64>)> {
        let x: Vec<_> = self
            .index_process_vars
            .iter()
            .map(|&i| vars[i as usize])
            .collect();
        let p: Vec<_> = self
            .index_parameter_vars
            .iter()
            .map(|&i| vars[i as usize])
            .collect();
        let eos = self.property_model.build_eos(&p);
        let (target, eq_constraints, ineq_constraints) = self.process.solve(&eos, &x)?;
        Ok((target, [eq_constraints, ineq_constraints].concat()))
    }
}

impl<
        E: Residual + IdealGas,
        M: MolecularRepresentation,
        R: PropertyModel<EquationOfState = E>,
        P: ProcessModel<E>,
    > OptimizationProblem<E, M, R, P>
{
    fn solve(
        &mut self,
        x0: Option<&[f64]>,
        y0: Option<&[f64]>,
        options: Option<&str>,
        target: bool,
    ) -> Result<(f64, [Vec<f64>; 3]), KnitroError> {
        let kc = Knitro::new()?;

        // declare process variables
        let index_process_vars = self.process.variables().setup_knitro(&kc, x0)?;

        // add equality constraints
        let index_eq_cons = kc.add_cons(self.process.equality_constraints())?;
        for &i in &index_eq_cons {
            kc.set_con_eqbnd(i, 0.0)?;
        }

        // add inequality constraints
        let index_ineq_cons = kc.add_cons(self.process.inequality_constraints())?;
        for &i in &index_ineq_cons {
            kc.set_con_lobnd(i, 0.0)?;
        }
        let index_cons = [index_eq_cons.clone(), index_ineq_cons.clone()].concat();

        // Set up CAMD formulation
        let (index_structure_vars, index_feature_vars) =
            self.molecules
                .setup_knitro(&kc, y0, &self.solutions, target)?;

        // Set up property model
        let index_parameter_vars = self.property_model.setup_knitro(&kc, &index_feature_vars)?;

        // Set up function callback
        let mut callback =
            OptimizationProblemCallback::new(self, &index_process_vars, &index_parameter_vars);
        let mut cb = kc.add_eval_callback(true, &index_cons, &mut callback)?;
        if let Some(options) = options {
            kc.load_param_file(options)?;
        }

        // only evaluate gradients with respect to explicit variables
        let vars = [index_process_vars.clone(), index_parameter_vars.clone()].concat();
        let jac_vars = vec![&vars as &[i32]; index_cons.len()].concat();
        let jac_cons: Vec<_> = index_cons
            .iter()
            .flat_map(|&i| vec![i; vars.len()])
            .collect();
        kc.set_cb_grad::<OptimizationProblemCallback<R, P>>(
            &mut cb, &vars, &jac_cons, &jac_vars, None,
        )?;

        // Solve (MI)NLP
        kc.solve()?;

        let t = kc.get_obj_value()?;
        let y = kc.get_var_primal_values(&index_structure_vars)?;
        let p = kc.get_var_primal_values(&index_parameter_vars)?;
        let x = kc.get_var_primal_values(&index_process_vars)?;

        Ok((t, [x.to_vec(), y.to_vec(), p.to_vec()]))
    }

    pub fn solve_target(
        &mut self,
        x0: &[f64],
        options: Option<&str>,
    ) -> Result<(f64, [Vec<f64>; 3]), KnitroError> {
        self.solve(Some(x0), None, options, true)
    }

    pub fn solve_knitro_once(
        &mut self,
        x0: &[f64],
        y0: Option<&[f64]>,
        options: Option<&str>,
    ) -> Result<OptimizationResult, KnitroError> {
        let (t, [x, y, p]) = self.solve(Some(x0), y0, options, false)?;

        let y: Vec<_> = y.iter().map(|y| y.round() as usize).collect();
        let smiles = self.molecules.smiles(&y);
        // println!("{t} {y:?} {x:?}");
        Ok(OptimizationResult::new(t, smiles, x, y, p))
    }

    pub fn solve_knitro(
        &mut self,
        x0: &[f64],
        y0: Option<&[f64]>,
        n_solutions: usize,
        options: Option<&str>,
    ) {
        for k in 0..n_solutions {
            match self.solve_knitro_once(x0, y0, options) {
                Ok(result) => {
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
                        println!(
                            "{k} {:.7} {:?} {:?} {:?}",
                            solution.target, solution.y, solution.smiles, solution.x
                        );
                    }
                }
                Err(e) => println!("\nRun {}\n{e}", k + 1),
            }
        }
    }
}
