use crate::process::ProcessModel;
use crate::*;
use feos_core::EosResult;
use knitro_rs::*;

pub trait MolecularRepresentationKnitro<const N: usize>: MolecularRepresentation<N> {
    fn setup_knitro(&self, kc: &Knitro) -> Result<(), KnitroError>;
}

impl<C: Clone> MolecularRepresentationKnitro<1> for FixedMolecule<C> {
    fn setup_knitro(&self, _: &Knitro) -> Result<(), KnitroError> {
        Ok(())
    }
}

pub trait MixtureConstraints {
    fn setup_knitro_mixture(&self, kc: &Knitro) -> Result<(), KnitroError>;
}

impl<C, M> MixtureConstraints for (FixedMolecule<C>, M) {
    fn setup_knitro_mixture(&self, _: &Knitro) -> Result<(), KnitroError> {
        Ok(())
    }
}

impl<
        C,
        M1: MolecularRepresentationKnitro<1, ChemicalRecord = C>,
        M2: MolecularRepresentationKnitro<1, ChemicalRecord = C>,
    > MolecularRepresentationKnitro<2> for (M1, M2)
where
    (M1, M2): MixtureConstraints,
{
    fn setup_knitro(&self, kc: &Knitro) -> Result<(), KnitroError> {
        self.0.setup_knitro(kc)?;
        self.1.setup_knitro(kc)?;
        self.setup_knitro_mixture(kc)
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

impl<
        M: MolecularRepresentationKnitro<N>,
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
        M: MolecularRepresentationKnitro<N>,
        R: PropertyModel<M::ChemicalRecord>,
        P: ProcessModel,
        const N: usize,
    > OptimizationProblem<M, R, P, N>
{
    fn evaluate(&self, x: &[f64]) -> EosResult<(f64, Vec<f64>)> {
        let n_y = self.molecules.variables();
        let (y, x) = x.split_at(n_y);
        let cr = self.molecules.build(y);
        let eos = self.property_model.build_eos(cr.into())?;
        let (_, target, constraints) = self.process.solve(&eos, x)?;
        Ok((target, constraints))
    }

    pub fn solve_knitro_once(
        &mut self,
        x0: &[f64],
        options: Option<&str>,
    ) -> Result<OptimizationResult, KnitroError> {
        let n_y = self.molecules.variables();

        let kc = Knitro::new()?;
        self.molecules.setup_knitro(&kc)?;

        // integer cuts
        for solution in self.solutions.iter() {
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
        let smiles = self.molecules.smiles(&y);
        // println!("{t} {y:?} {x:?}");
        Ok(OptimizationResult::new(t, smiles.into(), y, x.to_vec()))
    }

    pub fn solve_knitro(&mut self, x0: &[f64], n_solutions: usize, options: Option<&str>) {
        for k in 0..n_solutions {
            if let Ok(result) = self.solve_knitro_once(x0, options) {
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
