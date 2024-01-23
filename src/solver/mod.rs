use crate::process::ProcessModel;
use crate::*;
use feos::core::EosResult;
use knitro_rs::*;

mod outer_approximation;

impl<
        M: MolecularRepresentation<N>,
        R: PropertyModel<M::ChemicalRecord>,
        P: ProcessModel,
        const N: usize,
    > EvalCallback for OptimizationProblem<M, R, P, N>
{
    fn callback(&self, x: &[f64], obj: &mut f64, c: &mut [f64]) -> i32 {
        match self.evaluate(x) {
            Ok((target, constraints)) => {
                c.copy_from_slice(&constraints);
                *obj = target;
                0
            }
            Err(_) => {
                // println!("{e}");
                KN_RC_EVAL_ERR
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
    fn evaluate(&self, vars: &[f64]) -> EosResult<(f64, Vec<f64>)> {
        let n_p = self.molecules.parameter_variables().len();
        let n_x = self.process.variables().len();
        let (x, vars) = vars.split_at(n_x);
        let (p, y) = vars.split_at(n_p);
        let cr = self.molecules.build(y, p);
        let eos = self.property_model.build_eos(&cr)?;
        let (target, eq_constraints, ineq_constraints) = self.process.solve(&eos, x)?;
        Ok((target, [eq_constraints, ineq_constraints].concat()))
    }

    fn gradients(
        &mut self,
        kc: &Knitro,
        vars: &[i32],
        index_eq_cons: &[i32],
        index_ineq_cons: &[i32],
    ) -> Result<(Vec<f64>, Gradients, Vec<f64>), KnitroError> {
        let x = kc.get_var_primal_values(vars)?;
        let f = kc.get_obj_value()?;
        let h_vec = kc.get_con_values(index_eq_cons)?;
        let g_vec = kc.get_con_values(index_ineq_cons)?;
        let (_, nabla_f) = kc.get_objgrad_values()?;
        let f = Gradient::new(f, nabla_f);
        let (_, _, jac) = kc.get_jacobian_values()?;
        let lambda = kc.get_con_dual_values(index_eq_cons)?;

        let mut h = Vec::with_capacity(index_eq_cons.len());
        let mut g = Vec::with_capacity(index_ineq_cons.len());
        let mut k = 0;
        for val in h_vec {
            let mut nabla = Vec::with_capacity(vars.len());
            for _ in vars {
                nabla.push(jac[k]);
                k += 1;
            }
            h.push(Gradient::new(val, nabla));
        }
        for val in g_vec {
            let mut nabla = Vec::with_capacity(vars.len());
            for _ in vars {
                nabla.push(jac[k]);
                k += 1;
            }
            g.push(Gradient::new(val, nabla));
        }
        Ok((x, (f, g, h), lambda))
    }

    fn solve(
        &mut self,
        x0: Option<&[f64]>,
        y0: Option<&[f64]>,
        p0: Option<&[f64]>,
        options: Option<&str>,
        mode: OptimizationMode,
        gradients: Option<&mut (Vec<f64>, Gradients, Vec<f64>)>,
    ) -> Result<(f64, [Vec<f64>; 3]), KnitroError> {
        let kc = Knitro::new()?;

        // Set up process variables and constraints
        let [index_process_vars, index_eq_cons, index_ineq_cons] =
            self.process.setup_knitro(&kc, x0, mode)?;
        let index_cons = [index_eq_cons.clone(), index_ineq_cons.clone()].concat();

        // Set up CAMD formulation
        let [index_structure_vars, index_parameter_vars] =
            self.molecules
                .setup_knitro(&kc, y0, p0, &self.solutions, mode)?;

        // Set up function callback
        let mut cb = kc.add_eval_callback(true, &index_cons, self)?;
        if let Some(options) = options {
            kc.load_param_file(options)?;
        }

        // only evaluate gradients with respect to explicit variables
        let mut vars = vec![index_process_vars.clone()];
        if let OptimizationMode::MolecularDesign
        | OptimizationMode::Target
        | OptimizationMode::Gradients = mode
        {
            if index_parameter_vars.is_empty() {
                vars.push(index_structure_vars.clone())
            } else {
                vars.push(index_parameter_vars.clone())
            }
        }
        let vars = vars.concat();
        let jac_vars = vec![&vars as &[i32]; index_cons.len()].concat();
        let jac_cons: Vec<_> = index_cons
            .iter()
            .flat_map(|&i| vec![i; vars.len()])
            .collect();
        kc.set_cb_grad::<Self>(&mut cb, &vars, &jac_cons, &jac_vars, None)?;

        // Solve (MI)NLP
        kc.solve()?;
        let t = kc.get_obj_value()?;
        let y = kc.get_var_primal_values(&index_structure_vars)?;
        let p = kc.get_var_primal_values(&index_parameter_vars)?;
        let x = kc.get_var_primal_values(&index_process_vars)?;

        if let Some(gradients) = gradients {
            *gradients = self.gradients(&kc, &vars, &index_eq_cons, &index_ineq_cons)?;
        }

        Ok((t, [x.to_vec(), y.to_vec(), p.to_vec()]))
    }

    pub fn solve_target(
        &mut self,
        x0: &[f64],
        options: Option<&str>,
    ) -> Result<(f64, [Vec<f64>; 3]), KnitroError> {
        self.solve(
            Some(x0),
            None,
            None,
            options,
            OptimizationMode::Target,
            None,
        )
    }

    pub fn solve_fixed(
        &mut self,
        x0: Option<&[f64]>,
        y: &[f64],
        options: Option<&str>,
    ) -> Result<(f64, [Vec<f64>; 3]), KnitroError> {
        self.solve(
            x0,
            Some(y),
            None,
            options,
            OptimizationMode::FixedMolecule,
            None,
        )
    }

    pub fn solve_knitro_once(
        &mut self,
        x0: &[f64],
        y0: Option<&[f64]>,
        p0: Option<&[f64]>,
        options: Option<&str>,
    ) -> Result<OptimizationResult, KnitroError> {
        let (t, [x, y, p]) = self.solve(
            Some(x0),
            y0,
            p0,
            options,
            OptimizationMode::MolecularDesign,
            None,
        )?;

        let y: Vec<_> = y.iter().map(|y| *y as usize).collect();
        let smiles = self.molecules.smiles(&y);
        // println!("{t} {y:?} {x:?}");
        Ok(OptimizationResult::new(t, smiles.into(), x, y, p))
    }

    pub fn solve_knitro(
        &mut self,
        x0: &[f64],
        y0: Option<&[f64]>,
        p0: Option<&[f64]>,
        n_solutions: usize,
        options: Option<&str>,
    ) {
        for k in 0..n_solutions {
            match self.solve_knitro_once(x0, y0, p0, options) {
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
                Err(e) => println!("\nRun {}\n{e}", k + 1),
            }
        }
    }
}

#[derive(Clone, Default)]
struct Gradient {
    val: f64,
    nabla: Vec<f64>,
}

impl Gradient {
    fn new(val: f64, nabla: Vec<f64>) -> Self {
        Self { val, nabla }
    }

    fn linearize(self, x: &[f64]) -> (Vec<f64>, f64) {
        let rhs = self
            .nabla
            .iter()
            .zip(x.iter())
            .map(|(&n, &x)| n * x)
            .sum::<f64>();
        (self.nabla, rhs - self.val)
    }
}

type Gradients = (Gradient, Vec<Gradient>, Vec<Gradient>);
