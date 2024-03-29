use crate::process::ProcessModel;
use crate::*;
use knitro_rs::*;
use std::array;

struct ProcessCallback<'a, R, P, const N: usize> {
    parameters: Option<&'a [f64]>,
    property_model: &'a R,
    process: &'a P,
}

impl<'a, R, P, const N: usize> ProcessCallback<'a, R, P, N> {
    fn new<M>(problem: &'a OptimizationProblem<M, R, P, N>, parameters: Option<&'a [f64]>) -> Self {
        Self {
            parameters,
            property_model: &problem.property_model,
            process: &problem.process,
        }
    }
}

impl<'a, R: PropertyModel<N>, P: ProcessModel<R::EquationOfState>, const N: usize> EvalCallback
    for ProcessCallback<'a, R, P, N>
{
    fn callback(&self, x: &[f64], obj: &mut f64, c: &mut [f64]) -> i32 {
        let (x, p) = x.split_at(self.process.variables().len());
        let p = self.parameters.unwrap_or(p);
        let eos = self.property_model.build_eos(p);
        match self.process.solve(&eos, x) {
            Ok((target, eq_constraints, ineq_constraints)) => {
                let constraints = [eq_constraints, ineq_constraints].concat();
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
        M: MolecularRepresentation,
        R: PropertyModel<N>,
        P: ProcessModel<R::EquationOfState>,
        const N: usize,
    > OptimizationProblem<M, R, P, N>
{
    pub fn solve_fixed(
        &mut self,
        y: &[Vec<f64>; N],
        options: Option<&str>,
    ) -> Result<OptimizationResult, KnitroError> {
        let (y_usize, smiles) = self.smiles(y);
        if let Some(res) = self.solutions.iter().find(|r| r.y == y_usize) {
            Ok(res.clone())
        } else {
            let f = array::from_fn(|i| self.molecules[i].evaluate_feature_variables(&y[i]));
            let p = self.property_model.evaluate_parameter_variables(&f);

            // Set up process variables and constraints
            let kc = Knitro::new()?;

            // declare process variables
            self.process.variables().setup_knitro(&kc, None)?;

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

            // Set up function callback
            let mut callback = ProcessCallback::new(self, Some(&p));
            kc.add_eval_callback(true, &index_cons, &mut callback)?;
            if let Some(options) = options {
                kc.load_param_file(options)?;
            }

            // Solve NLP
            kc.solve()?;

            let target = kc.get_obj_value()?;
            let x = kc.get_var_primal_values_all()?;

            let res = OptimizationResult::new(target, smiles, x, y_usize);
            self.solutions.insert(res.clone());
            Ok(res)
        }
    }

    pub(super) fn solve_gradients(
        &mut self,
        y: &[Vec<f64>; N],
        x: &[f64],
        options: Option<&str>,
    ) -> Result<(Vec<f64>, Gradients, Vec<f64>), KnitroError> {
        let f = array::from_fn(|i| self.molecules[i].evaluate_feature_variables(&y[i]));
        let p = self.property_model.evaluate_parameter_variables(&f);

        let kc = Knitro::new()?;

        // declare process variables
        for &x in x {
            Variable::fixed(x).setup_knitro(&kc, "fixed_process_variable")?;
        }

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

        // Set up property model
        for p in p {
            Variable::fixed(p).setup_knitro(&kc, "fixed_param")?;
        }

        // Set up function callback
        let mut callback = ProcessCallback::new(self, None);
        kc.add_eval_callback(true, &index_cons, &mut callback)?;
        if let Some(options) = options {
            kc.load_param_file(options)?;
        }

        // Solve NLP
        kc.solve()?;

        // Determine gradients
        let x = kc.get_var_primal_values_all()?;
        let n_vars = x.len();
        let f = kc.get_obj_value()?;
        let h_vec = kc.get_con_values(&index_eq_cons)?;
        let g_vec = kc.get_con_values(&index_ineq_cons)?;
        let (_, nabla_f) = kc.get_objgrad_values()?;
        let f = Gradient::new(f, nabla_f);
        let (_, _, jac) = kc.get_jacobian_values()?;
        let lambda = kc.get_con_dual_values(&index_eq_cons)?;

        let mut h = Vec::with_capacity(index_eq_cons.len());
        let mut g = Vec::with_capacity(index_ineq_cons.len());
        let mut k = 0;
        for val in h_vec {
            let mut nabla = Vec::with_capacity(n_vars);
            for _ in 0..n_vars {
                nabla.push(jac[k]);
                k += 1;
            }
            h.push(Gradient::new(val, nabla));
        }
        for val in g_vec {
            let mut nabla = Vec::with_capacity(n_vars);
            for _ in 0..n_vars {
                nabla.push(jac[k]);
                k += 1;
            }
            g.push(Gradient::new(val, nabla));
        }
        Ok((x, (f, g, h), lambda))
    }
}

#[derive(Clone, Default)]
pub(super) struct Gradient {
    val: f64,
    nabla: Vec<f64>,
}

impl Gradient {
    fn new(val: f64, nabla: Vec<f64>) -> Self {
        Self { val, nabla }
    }

    pub(super) fn linearize(self, x: &[f64]) -> (Vec<f64>, f64) {
        let rhs = self
            .nabla
            .iter()
            .zip(x.iter())
            .map(|(&n, &x)| n * x)
            .sum::<f64>();
        (self.nabla, rhs - self.val)
    }
}

pub(super) type Gradients = (Gradient, Vec<Gradient>, Vec<Gradient>);
