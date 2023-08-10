#![allow(clippy::missing_safety_doc)]
use crate::process::ProcessModel;
use crate::*;
use feos_core::EosResult;
use knitro_rs::*;
use std::sync::Arc;

pub trait MolecularRepresentationKnitro: MolecularRepresentation {
    fn setup_knitro(&self, kc: &Knitro) -> Result<(), KnitroError>;
}

impl MolecularRepresentationKnitro for FixedMolecule {
    fn setup_knitro(&self, _: &Knitro) -> Result<(), KnitroError> {
        Ok(())
    }
}

pub trait ProcessModelKnitro: ProcessModel {
    fn setup_knitro(&self, kc: &Knitro, x0: &[f64]) -> Result<Vec<i32>, KnitroError> {
        // declare continuous variables
        let variables = self.variables();
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
        if self.binary_variables() > 0 {
            kc.add_vars(self.binary_variables(), Some(KN_VARTYPE_BINARY))?;
        }

        // add constraints
        let constraints = self.constraints();
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
}

impl<T: ProcessModel> ProcessModelKnitro for T {}

struct OptimizationProblemCallback<'a, M, R, P> {
    molecules: &'a [M],
    property_model: &'a R,
    process: &'a P,
    solutions: &'a [OptimizationResult],
}

impl<'a, M, R, P> From<&'a OptimizationProblem<M, R, P>>
    for OptimizationProblemCallback<'a, M, R, P>
{
    fn from(problem: &'a OptimizationProblem<M, R, P>) -> Self {
        Self {
            molecules: &problem.molecules,
            property_model: &problem.property_model,
            process: &problem.process,
            solutions: &problem.solutions,
        }
    }
}

impl<M: MolecularRepresentationKnitro, R: PropertyModel<M::ChemicalRecord>, P: ProcessModel>
    Callback for OptimizationProblemCallback<'_, M, R, P>
{
    fn callback(&self, x: &[f64], c: &mut [f64]) -> f64 {
        match self.evaluate(x) {
            Ok((target, constraints)) => {
                c.copy_from_slice(&constraints);
                target
            }
            _ => {
                c.iter_mut().for_each(|c| *c = -1.0);
                0.0
            }
        }
    }
}

impl<M: MolecularRepresentationKnitro, R: PropertyModel<M::ChemicalRecord>, P: ProcessModel>
    OptimizationProblemCallback<'_, M, R, P>
{
    fn evaluate(&self, x: &[f64]) -> EosResult<(f64, Vec<f64>)> {
        let n_y = self.molecules.iter().map(|m| m.variables()).sum();
        let (mut y, x) = x.split_at(n_y);
        let cr = self
            .molecules
            .iter()
            .map(|molecule| {
                let (y1, y2) = y.split_at(molecule.variables());
                y = y2;
                molecule.build(y1.to_vec())
            })
            .collect();
        let eos = Arc::new(self.property_model.build_eos(cr)?);
        let (_, target, constraints) = self.process.solve(&eos, x)?;
        Ok((target, constraints))
    }

    pub fn solve_knitro_once(
        &mut self,
        x0: &[f64],
        options: Option<&str>,
    ) -> Result<OptimizationResult, KnitroError> {
        let n_y = self.molecules.iter().map(|m| m.variables()).sum();

        let kc = Knitro::new()?;
        self.molecules
            .iter()
            .try_for_each(|molecule| molecule.setup_knitro(&kc))?;

        // integer cuts
        for solution in self.solutions {
            let y0 = &solution.y[..n_y];
            let vars = Vec::from_iter(0..n_y as i32);
            let c = kc.add_con()?;
            let qcoefs: Vec<_> = y0.iter().map(|_| 1.0).collect();
            let lcoefs: Vec<_> = y0.iter().map(|&n0| -2.0 * n0 as f64).collect();
            let lbond = 1.0 - y0.iter().map(|n0| n0.pow(2) as f64).sum::<f64>();
            kc.add_con_quadratic_struct_one(c, &vars, &vars, &qcoefs)?;
            kc.add_con_linear_struct_one(c, &vars, &lcoefs)?;
            kc.set_con_lobnd(c, lbond)?;
        }

        let index_cons = self.process.setup_knitro(&kc, x0)?;
        let cb = kc.add_eval_callback(true, &index_cons, self)?;
        if let Some(options) = options {
            kc.load_param_file(options)?;
        }
        kc.set_cb_user_params(cb, self)?;
        kc.solve()?;
        let t = kc.get_obj_value()?;
        let x = kc.get_var_primal_values_all()?;
        let (y, x) = x.split_at(n_y);
        let y: Vec<_> = y.iter().map(|y| *y as usize).collect();
        println!("{t} {y:?} {x:?}");
        Ok(OptimizationResult::new(t, y, x.to_vec()))
    }
}

impl<M: MolecularRepresentationKnitro, R: PropertyModel<M::ChemicalRecord>, P: ProcessModel>
    OptimizationProblem<M, R, P>
{
    pub fn solve_knitro(
        &mut self,
        x0: &[f64],
        n_solutions: usize,
        options: Option<&str>,
    ) -> Result<(), KnitroError> {
        for _ in 0..n_solutions {
            let mut problem = OptimizationProblemCallback::from(&*self);
            let result = problem.solve_knitro_once(x0, options)?;
            self.add_solution(result);
        }
        Ok(())
    }
}
