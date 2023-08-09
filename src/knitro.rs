#![allow(clippy::missing_safety_doc)]
use crate::process::ProcessModel;
use crate::*;
use feos_core::EosResult;
use knitro_rs::*;
use rayon::iter::*;
use std::sync::Arc;

pub trait MolecularRepresentationKnitro: MolecularRepresentation {
    fn setup_knitro(
        &self,
        kc: &Knitro,
        solutions: &[OptimizationResult],
    ) -> Result<(), KnitroError>;
}

impl MolecularRepresentationKnitro for SuperMolecule {
    fn setup_knitro(
        &self,
        kc: &Knitro,
        solutions: &[OptimizationResult],
    ) -> Result<(), KnitroError> {
        let n_y = self.variables();

        // binary variables
        let index_vars = kc.add_vars(n_y, Some(KN_VARTYPE_BINARY))?;

        // maximum size
        let index_con = kc.add_con()?;
        let (coefs, size) = self.size_constraint();
        kc.add_con_linear_struct_one(index_con, &index_vars, &coefs)?;
        kc.set_con_upbnd(index_con, size as f64)?;

        // integer cuts
        for solution in solutions {
            let index_con = kc.add_con()?;
            let coefs: Vec<_> = solution.y.iter().map(|&y| 2.0 * y as f64 - 1.0).collect();
            kc.add_con_linear_struct_one(index_con, &index_vars, &coefs)?;
            kc.set_con_upbnd(index_con, solution.y.iter().sum::<usize>() as f64 - 1.0)?;
        }

        // functional group
        let index_vars = self.functional_group_constraint();
        if !index_vars.is_empty() {
            let index_con = kc.add_con()?;
            kc.add_con_linear_struct_one(index_con, &index_vars, &vec![1.0; index_vars.len()])?;
            kc.set_con_eqbnd(index_con, 1.0)?;
        }

        // minimum size
        let index_con = kc.add_con()?;
        kc.add_con_linear_term(index_con, index_vars.len() as i32, 1.0)?;
        kc.set_con_eqbnd(index_con, 1.0)?;

        // bond constraints
        for bond in self.bond_constraints() {
            let index_con = kc.add_con()?;
            kc.add_con_linear_struct_one(index_con, &bond, &[1.0, -1.0])?;
            kc.set_con_lobnd(index_con, 0.0)?;
        }

        // symmetry constraints
        for (index_vars, coefs) in self.symmetry_constraints() {
            let index_con = kc.add_con()?;
            kc.add_con_linear_struct_one(index_con, &index_vars, &coefs)?;
            kc.set_con_lobnd(index_con, 0.0)?;
        }

        Ok(())
    }
}

impl MolecularRepresentationKnitro for FixedMolecule {
    fn setup_knitro(&self, _: &Knitro, _: &[OptimizationResult]) -> Result<(), KnitroError> {
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
    molecule: &'a M,
    property_model: &'a R,
    process: &'a P,
    solutions: &'a [OptimizationResult],
    chemical: Option<&'a str>,
}

impl<'a, M, R, P> From<&'a OptimizationProblem<M, R, P>>
    for OptimizationProblemCallback<'a, M, R, P>
{
    fn from(problem: &'a OptimizationProblem<M, R, P>) -> Self {
        Self {
            molecule: &problem.molecule,
            property_model: &problem.property_model,
            process: &problem.process,
            solutions: &problem.solutions,
            chemical: None,
        }
    }
}

impl<R, P> MetaOptimizationProblem<R, P> {
    fn callback<'a>(
        &'a self,
        molecule: &'a SuperMolecule,
        solutions: &'a [OptimizationResult],
        chemical: &'a str,
    ) -> OptimizationProblemCallback<'a, SuperMolecule, R, P> {
        OptimizationProblemCallback {
            molecule,
            property_model: &self.property_model,
            process: &self.process,
            solutions,
            chemical: Some(chemical),
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
        let n_y = self.molecule.variables();
        let (y, x) = x.split_at(n_y);
        let cr = self.molecule.build(y.to_vec());
        let eos = Arc::new(self.property_model.build_eos(cr)?);
        let (_, target, constraints) = self.process.solve(&eos, x)?;
        Ok((target, constraints))
    }

    pub fn solve_knitro_once(
        &mut self,
        x0: &[f64],
        options: Option<&str>,
    ) -> Result<OptimizationResult, KnitroError> {
        let n_y = self.molecule.variables();

        let kc = Knitro::new()?;
        self.molecule.setup_knitro(&kc, self.solutions)?;
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
        println!("{} {t} {y:?} {x:?}", self.chemical.unwrap_or(""));
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

impl<R: PropertyModel<SegmentAndBondCount> + Sync, P: ProcessModel + Sync>
    MetaOptimizationProblem<R, P>
{
    pub fn solve_knitro(&mut self, x0: &[f64], n_solutions: usize, options: Option<&str>) {
        let candidates = self
            .molecules
            .par_iter()
            .map(|(chemical, molecule)| {
                let mut problem = self.callback(molecule, &[], chemical);
                (
                    chemical.clone(),
                    vec![problem.solve_knitro_once(x0, options).unwrap()],
                )
            })
            .collect();
        self.candidates = candidates;

        for _ in 0..n_solutions {
            let chemical = self.best_candidate();
            let mut problem = self.callback(
                &self.molecules[&chemical],
                &self.candidates[&chemical],
                &chemical,
            );
            let result = problem.solve_knitro_once(x0, options).unwrap();
            self.update_candidates(&chemical, result);
        }
    }

    pub fn find_best_knitro(&mut self, x0: &[f64], depth: usize, options: Option<&str>) {
        let candidates = self
            .molecules
            .par_iter()
            .map(|(chemical, molecule)| {
                let mut results = vec![];
                for _ in 0..depth {
                    let mut problem = self.callback(molecule, &results, chemical);
                    let result = problem.solve_knitro_once(x0, options).unwrap();
                    results.push(result);
                }
                (chemical.clone(), results)
            })
            .collect();
        self.candidates = candidates;

        let mut best = 0.0;
        let mut best_result = None;
        let mut best_chemical = String::new();
        for (chemical, results) in &self.candidates {
            for result in results {
                if result.target < best {
                    best = result.target;
                    best_result = Some(result.clone());
                    best_chemical = chemical.clone();
                }
            }
        }
        self.solutions = vec![(best_chemical, best_result.unwrap())];
    }
}
