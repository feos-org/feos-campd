use crate::process::ProcessModel;
use crate::*;
use feos::core::{EosResult, IdealGas, Residual};
use knitro_rs::*;
use std::array;

mod outer_approximation;
mod process_optimization;
pub use outer_approximation::OuterApproximationAlgorithm;

struct OptimizationProblemCallback<'a, R, P, const N: usize> {
    property_model: &'a R,
    process: &'a P,
    process_vars: usize,
    parameter_vars: usize,
}

impl<'a, R, P, const N: usize> OptimizationProblemCallback<'a, R, P, N> {
    fn new<E, M>(
        problem: &'a OptimizationProblem<E, M, R, P, N>,
        process_vars: usize,
        parameter_vars: usize,
    ) -> Self {
        Self {
            property_model: &problem.property_model,
            process: &problem.process,
            process_vars,
            parameter_vars,
        }
    }
}

impl<
        'a,
        E: Residual + IdealGas,
        R: PropertyModel<N, EquationOfState = E>,
        P: ProcessModel<E>,
        const N: usize,
    > EvalCallback for OptimizationProblemCallback<'a, R, P, N>
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

impl<
        'a,
        E: Residual + IdealGas,
        R: PropertyModel<N, EquationOfState = E>,
        P: ProcessModel<E>,
        const N: usize,
    > OptimizationProblemCallback<'a, R, P, N>
{
    fn evaluate(&self, vars: &[f64]) -> EosResult<(f64, Vec<f64>)> {
        let x = &vars[..self.process_vars];
        let p = &vars[vars.len() - self.parameter_vars..];
        let eos = self.property_model.build_eos(p);
        let (target, eq_constraints, ineq_constraints) = self.process.solve(&eos, x)?;
        Ok((target, [eq_constraints, ineq_constraints].concat()))
    }
}

#[allow(clippy::type_complexity)]
impl<
        E: Residual + IdealGas,
        M: MolecularRepresentation,
        R: PropertyModel<N, EquationOfState = E>,
        P: ProcessModel<E>,
        const N: usize,
    > OptimizationProblem<E, M, R, P, N>
{
    fn solve(
        &mut self,
        x0: Option<&[f64]>,
        y0: Option<&[Vec<f64>; N]>,
        options: Option<&str>,
        target: bool,
    ) -> Result<(f64, Vec<f64>, [Vec<f64>; N]), KnitroError> {
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
        let mut index_structure_vars = array::from_fn(|_| Vec::new());
        let index_feature_vars = array::from_fn(|i| {
            let (y, f) = self.molecules[i]
                .setup_knitro(&kc, y0.map(|y0| &y0[i][..]), target)
                .unwrap();
            index_structure_vars[i].extend(y);
            f
        });

        // integer cuts
        let all_structure_vars = index_structure_vars.concat();
        for solution in &self.solutions {
            let y0 = &solution.y;
            let coefs: Vec<_> = y0.iter().map(|&y0| 1.0 - 2.0 * y0 as f64).collect();
            let lobnd = 1.0 - y0.iter().map(|&y0| y0 as f64).sum::<f64>();
            Constraint::new()
                .linear_struct(all_structure_vars.clone(), coefs)
                .lobnd(lobnd)
                .setup_knitro(&kc)?;
        }

        // Set up property model
        let index_parameter_vars = self.property_model.setup_knitro(&kc, &index_feature_vars)?;

        // Set up function callback
        let mut callback = OptimizationProblemCallback::new(
            self,
            index_process_vars.len(),
            index_parameter_vars.len(),
        );
        let mut cb = kc.add_eval_callback(true, &index_cons, &mut callback)?;
        if let Some(options) = options {
            kc.load_param_file(options)?;
        }

        // only evaluate gradients with respect to explicit variables
        let vars = [index_process_vars.clone(), index_parameter_vars].concat();
        let jac_vars = vec![&vars as &[i32]; index_cons.len()].concat();
        let jac_cons: Vec<_> = index_cons
            .iter()
            .flat_map(|&i| vec![i; vars.len()])
            .collect();
        kc.set_cb_grad::<OptimizationProblemCallback<R, P, N>>(
            &mut cb, &vars, &jac_cons, &jac_vars, None,
        )?;

        // Solve (MI)NLP
        kc.solve()?;

        let t = kc.get_obj_value()?;
        let y = index_structure_vars.map(|y| kc.get_var_primal_values(&y).unwrap());
        let x = kc.get_var_primal_values(&index_process_vars)?;

        Ok((t, x, y))
    }

    pub fn solve_target(
        &mut self,
        x0: &[f64],
        options: Option<&str>,
    ) -> Result<(f64, Vec<f64>, [Vec<f64>; N]), KnitroError> {
        self.solve(Some(x0), None, options, true)
    }

    pub fn solve_knitro_once(
        &mut self,
        x0: &[f64],
        y0: Option<&[Vec<f64>; N]>,
        options: Option<&str>,
    ) -> Result<OptimizationResult, KnitroError> {
        let (t, x, y) = self.solve(Some(x0), y0, options, false)?;

        let (y, smiles) = self.smiles(&y);
        // println!("{t} {y:?} {x:?}");
        Ok(OptimizationResult::new(t, smiles, x, y))
    }

    pub fn solve_knitro(
        &mut self,
        x0: &[f64],
        y0: Option<&[Vec<f64>; N]>,
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
