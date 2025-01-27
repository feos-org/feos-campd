use super::dual_vec_multiple::DualSVecMult;
use super::{Gradient, MixedIntegerNonLinearProgram, OptimizationResult, OuterApproximation};
use ipopt::{Ipopt, IpoptOption, SolveStatus};
use ipopt_ad::{ADProblem, BasicADProblem, CachedADProblem};
use nalgebra::{Const, DVector, SMatrix, SVector, U1};
use num_dual::{Derivative, DualNum};
use std::fmt::Debug;

struct Nlp<'a, P, const N_Y1: usize, const N_Y2: usize> {
    y: SMatrix<f64, N_Y1, N_Y2>,
    minlp: &'a P,
}

impl<'a, P, const N_Y1: usize, const N_Y2: usize> Nlp<'a, P, N_Y1, N_Y2> {
    pub fn new(y: SMatrix<f64, N_Y1, N_Y2>, minlp: &'a P) -> Self {
        Self { y, minlp }
    }
}

impl<
        P: MixedIntegerNonLinearProgram<N_X, N_Y1, N_Y2>,
        const N_X: usize,
        const N_Y1: usize,
        const N_Y2: usize,
    > BasicADProblem<N_X> for Nlp<'_, P, N_Y1, N_Y2>
{
    fn bounds(&self) -> ([f64; N_X], [f64; N_X]) {
        let vars = self.minlp.x_variables().data.0[0];
        (vars.map(|(l, _, _)| l), vars.map(|(_, u, _)| u))
    }

    fn initial_point(&self) -> [f64; N_X] {
        self.minlp.x_variables().data.0[0].map(|(_, _, i)| i)
    }

    fn constraint_bounds(&self) -> (Vec<f64>, Vec<f64>) {
        self.minlp
            .constraints()
            .into_iter()
            .map(|c| (c.lower_bound(), c.upper_bound()))
            .unzip()
    }
}

impl<
        P: MixedIntegerNonLinearProgram<N_X, N_Y1, N_Y2>,
        const N_X: usize,
        const N_Y1: usize,
        const N_Y2: usize,
    > CachedADProblem<N_X> for Nlp<'_, P, N_Y1, N_Y2>
{
    type Error = P::Error;

    fn evaluate<D: DualNum<f64> + Copy>(&self, x: [D; N_X]) -> Result<(D, Vec<D>), P::Error> {
        self.minlp.evaluate(SVector::from(x), self.y.map(D::from))
    }
}

impl<
        M: MixedIntegerNonLinearProgram<N_X, N_Y1, N_Y2>,
        const N_X: usize,
        const N_Y1: usize,
        const N_Y2: usize,
    > OuterApproximation<'_, M, N_X, N_Y1, N_Y2>
where
    M::Error: Debug,
{
    pub fn solve_nlp(
        &mut self,
        y: SMatrix<f64, N_Y1, N_Y2>,
        s: Vec<f64>,
    ) -> Option<&OptimizationResult<N_X, N_Y1, N_Y2>> {
        self.solve_nlp_with_options::<f64>(y, s, &[])
    }

    pub fn solve_nlp_with_options<'b, O: Into<IpoptOption<'b>> + Copy>(
        &mut self,
        y: SMatrix<f64, N_Y1, N_Y2>,
        s: Vec<f64>,
        options: &[(&str, O)],
    ) -> Option<&OptimizationResult<N_X, N_Y1, N_Y2>> {
        let key = self.minlp.y_to_string(&y);
        if self.known_solutions.contains_key(&key) {
            return self.known_solutions.get(&key);
        }

        let optim = Nlp::new(y, self.minlp);
        let Ok(problem) = ADProblem::new_cached(optim) else {
            return None;
        };
        let mut ipopt = Ipopt::new(problem).unwrap();
        for &(s, o) in options {
            ipopt.set_option(s, o);
        }
        let res = ipopt.solve();
        if let SolveStatus::SolveSucceeded = res.status {
            let x = SVector::from_column_slice(res.solver_data.solution.primal_variables);
            let lambda = res.solver_data.solution.constraint_multipliers;
            let (objective, constraints) = Self::gradients(self.minlp, y, x);
            let result = OptimizationResult::new(
                key.clone(),
                y,
                s,
                objective,
                constraints,
                x,
                DVector::from_column_slice(lambda),
            );
            self.known_solutions.insert(key.clone(), result);
            self.known_solutions.get(&key)
        } else {
            // println!("{}", res.status);
            None
        }
    }

    fn gradients(
        minlp: &M,
        y: SMatrix<f64, N_Y1, N_Y2>,
        x: SVector<f64, N_X>,
    ) -> (Gradient<N_X, N_Y1, N_Y2>, Vec<Gradient<N_X, N_Y1, N_Y2>>) {
        let mut y_dual = y.map(DualSVecMult::from_re);
        let mut x_dual = x.map(DualSVecMult::from_re);
        for (i, y) in y_dual.iter_mut().enumerate() {
            y.eps1 = Derivative::derivative_generic(Const::<N_Y1>, Const::<N_Y2>, i);
        }
        for (i, x) in x_dual.iter_mut().enumerate() {
            x.eps2 = Derivative::derivative_generic(Const::<N_X>, U1, i);
        }
        let (f, con) = minlp
            .evaluate(x_dual, y_dual)
            .expect("Unexpected error ocurred in the calculation of gradients!");
        let f = (
            f.re,
            f.eps1.unwrap_generic(Const::<N_Y1>, Const::<N_Y2>),
            f.eps2.unwrap_generic(Const::<N_X>, U1),
        );
        let con = con
            .into_iter()
            .map(|con| {
                (
                    con.re,
                    con.eps1.unwrap_generic(Const::<N_Y1>, Const::<N_Y2>),
                    con.eps2.unwrap_generic(Const::<N_X>, U1),
                )
            })
            .collect();
        (f, con)
    }
}
