use core::f64;
use good_lp::{
    constraint, variable, Constraint as LinearConstraint, Expression, ProblemVariables, Solution,
    Solver, SolverModel, Variable,
};
use ipopt::IpoptOption;
use nalgebra::{DVector, SMatrix, SVector};
use num_dual::DualNum;
use std::collections::{HashMap, HashSet};
use std::fmt::Debug;

mod dual_vec_multiple;
mod nlp;

#[derive(Clone)]
pub struct OptimizationResult<const N_X: usize, const N_Y1: usize, const N_Y2: usize> {
    pub key: String,
    pub objective: Gradient<N_X, N_Y1, N_Y2>,
    pub constraints: Vec<Gradient<N_X, N_Y1, N_Y2>>,
    pub x: SVector<f64, N_X>,
    pub y: SMatrix<f64, N_Y1, N_Y2>,
    pub s: Vec<f64>,
    pub lambda: DVector<f64>,
}

impl<const N_X: usize, const N_Y1: usize, const N_Y2: usize> OptimizationResult<N_X, N_Y1, N_Y2> {
    pub fn new(
        key: String,
        y: SMatrix<f64, N_Y1, N_Y2>,
        s: Vec<f64>,
        objective: Gradient<N_X, N_Y1, N_Y2>,
        constraints: Vec<Gradient<N_X, N_Y1, N_Y2>>,
        x: SVector<f64, N_X>,
        lambda: DVector<f64>,
    ) -> Self {
        Self {
            key,
            objective,
            constraints,
            x,
            y,
            s,
            lambda,
        }
    }
}

type Gradient<const N_X: usize, const N_Y1: usize, const N_Y2: usize> =
    (f64, SMatrix<f64, N_Y1, N_Y2>, SVector<f64, N_X>);

#[derive(Clone)]
pub struct OptimizationOptions<'a> {
    pub min_iter: usize,
    pub max_iter: usize,
    pub zero_tol: f64,
    pub nlp_options: Vec<(&'a str, IpoptOption<'a>)>,
}

impl Default for OptimizationOptions<'_> {
    fn default() -> Self {
        Self {
            min_iter: 5,
            max_iter: 50,
            zero_tol: 1e-6,
            nlp_options: vec![("print_level", IpoptOption::Int(0))],
        }
    }
}

#[derive(Clone, Copy)]
pub enum GeneralConstraint {
    Equality(f64),
    Inequality(Option<f64>, Option<f64>),
}

impl GeneralConstraint {
    pub fn lower_bound(&self) -> f64 {
        match self {
            Self::Equality(e) => *e,
            Self::Inequality(l, _) => l.unwrap_or(f64::NEG_INFINITY),
        }
    }

    pub fn upper_bound(&self) -> f64 {
        match self {
            Self::Equality(e) => *e,
            Self::Inequality(_, u) => u.unwrap_or(f64::INFINITY),
        }
    }
}

pub trait MixedIntegerNonLinearProgram<const N_X: usize, const N_Y1: usize, const N_Y2: usize> {
    type Error;

    fn x_variables(&self) -> SVector<(f64, f64, f64), N_X>;

    fn y_variables(&self) -> SMatrix<(i32, i32), N_Y1, N_Y2>;

    fn linear_constraints(&self, y: SMatrix<Variable, N_Y1, N_Y2>) -> Vec<LinearConstraint>;

    fn constraints(&self) -> Vec<GeneralConstraint>;

    fn evaluate<D: DualNum<f64> + Copy>(
        &self,
        x: SVector<D, N_X>,
        y: SMatrix<D, N_Y1, N_Y2>,
    ) -> Result<(D, Vec<D>), Self::Error>;

    fn y_to_string(&self, y: &SMatrix<f64, N_Y1, N_Y2>) -> String;

    fn exclude_solutions(&self, s: &[f64]) -> Vec<Vec<f64>> {
        vec![s.to_vec()]
    }
}

pub struct OuterApproximation<
    'a,
    M: MixedIntegerNonLinearProgram<N_X, N_Y1, N_Y2>,
    const N_X: usize,
    const N_Y1: usize,
    const N_Y2: usize,
> {
    minlp: &'a M,
    known_solutions: HashMap<String, OptimizationResult<N_X, N_Y1, N_Y2>>,
    excluded_solutions: Vec<Vec<f64>>,
}

impl<
        'a,
        M: MixedIntegerNonLinearProgram<N_X, N_Y1, N_Y2>,
        const N_X: usize,
        const N_Y1: usize,
        const N_Y2: usize,
    > OuterApproximation<'a, M, N_X, N_Y1, N_Y2>
where
    M::Error: Debug,
{
    pub fn new(minlp: &'a M) -> Self {
        Self {
            minlp,
            known_solutions: HashMap::new(),
            excluded_solutions: vec![],
        }
    }

    fn add_oa_cuts(
        &self,
        constraints: &mut Vec<LinearConstraint>,
        x: SVector<Variable, N_X>,
        y: SMatrix<Variable, N_Y1, N_Y2>,
        mu: Variable,
        result: &OptimizationResult<N_X, N_Y1, N_Y2>,
        zero_tol: f64,
    ) {
        let (f, grad_x, grad_y) = &result.objective;
        let con = &result.constraints;

        let y_expr = y.iter().zip(result.y.data.0[0]).zip(grad_y.data.0[0]);
        let x_expr = x.iter().zip(result.x.data.0[0]).zip(grad_x.data.0[0]);
        let expr: Expression = y_expr.chain(x_expr).map(|((&x, x0), j)| (x - x0) * j).sum();
        constraints.push(constraint!(expr + *f <= mu));

        for ((constraint, (c, jac_y, jac_x)), &l) in self
            .minlp
            .constraints()
            .into_iter()
            .zip(con)
            .zip(result.lambda.iter())
        {
            match constraint {
                GeneralConstraint::Equality(eq) => {
                    let sign = 1f64.copysign(-l);
                    let y_expr = y.iter().zip(result.y.data.0[0]).zip(jac_y.data.0[0]);
                    let x_expr = x.iter().zip(result.x.data.0[0]).zip(jac_x.data.0[0]);
                    let expr: Expression =
                        y_expr.chain(x_expr).map(|((&x, x0), j)| (x - x0) * j).sum();
                    constraints.push(constraint!(expr * sign + *c <= eq));
                }
                GeneralConstraint::Inequality(lo, up) => {
                    if let Some(up) = up {
                        // Check if constraint is active
                        if up - c < zero_tol {
                            let y_expr = y.iter().zip(result.y.data.0[0]).zip(jac_y.data.0[0]);
                            let x_expr = x.iter().zip(result.x.data.0[0]).zip(jac_x.data.0[0]);
                            let expr: Expression =
                                y_expr.chain(x_expr).map(|((&x, x0), j)| (x - x0) * j).sum();
                            constraints.push(constraint!(expr + *c <= up));
                        }
                    }
                    if let Some(lo) = lo {
                        // Check if constraint is active
                        if c - lo < zero_tol {
                            let y_expr = y.iter().zip(result.y.data.0[0]).zip(jac_y.data.0[0]);
                            let x_expr = x.iter().zip(result.x.data.0[0]).zip(jac_x.data.0[0]);
                            let expr: Expression =
                                y_expr.chain(x_expr).map(|((&x, x0), j)| (x - x0) * j).sum();
                            constraints.push(constraint!(expr + *c >= lo));
                        }
                    }
                }
            }
        }
    }

    fn add_integer_cut(constraints: &mut Vec<LinearConstraint>, s: &[Variable], s0: &[f64]) {
        let expr = s.iter().zip(s0).map(|(&s, &s0)| s - 2.0 * s0 * s + s0);
        constraints.push(constraint!(expr.sum::<Expression>() >= 1.0));
    }

    #[expect(clippy::type_complexity)]
    pub fn solve_milp<S: Solver>(
        &mut self,
        solver: S,
        oa_cuts: &[String],
        zero_tol: f64,
    ) -> Result<(SMatrix<f64, N_Y1, N_Y2>, Vec<f64>), <S::Model as SolverModel>::Error> {
        let mut model = ProblemVariables::new();
        let mut constraints = Vec::new();

        // binary variables for integer cuts
        let mut s = Vec::new();

        // discrete variables
        let y = self.minlp.y_variables().map(|(l, u)| {
            let y = model.add(variable().integer().bounds(l..u));

            // add binary variables for integer cuts
            let vars = model.add_vector(variable().binary(), (u - l) as usize);
            for vars in vars.windows(2) {
                constraints.push(constraint!(vars[0] >= vars[1]));
            }
            constraints.push(constraint!(l + vars.iter().sum::<Expression>() == y));
            s.extend_from_slice(&vars);

            y
        });

        // linear constraints
        constraints.append(&mut self.minlp.linear_constraints(y));

        // process variables
        let x = self
            .minlp
            .x_variables()
            .map(|(l, u, _)| model.add(variable().bounds(l..u)));

        // epigraph variable
        let mu = model.add_variable();

        // integer cuts
        for solution in &self.excluded_solutions {
            Self::add_integer_cut(&mut constraints, &s, solution);
        }

        // OA cuts
        for key in oa_cuts {
            let solution = &self.known_solutions[key];
            self.add_oa_cuts(&mut constraints, x, y, mu, solution, zero_tol);
        }

        // setup solver
        let mut model = model.minimise(mu).using(solver);

        // add constraints
        constraints.into_iter().for_each(|c| {
            model.add_constraint(c);
        });

        // solve MILP
        model.solve().map(|solution| {
            let y = y.map(|y| solution.value(y).round());
            let s: Vec<_> = s.iter().map(|s| solution.value(*s).round()).collect();

            (y, s)
        })
    }

    fn calculate_s(&self, y: SMatrix<f64, N_Y1, N_Y2>) -> Vec<f64> {
        let mut s = Vec::new();
        self.minlp
            .y_variables()
            .iter()
            .zip(y.iter())
            .for_each(|(&(l, u), &y)| {
                s.extend_from_slice(&vec![1.0; y as usize - l as usize]);
                s.extend_from_slice(&vec![0.0; u as usize - y as usize]);
            });
        s
    }

    pub fn solve<S>(
        &mut self,
        y_init: SMatrix<f64, N_Y1, N_Y2>,
        solver: &S,
        options: &OptimizationOptions,
    ) -> Vec<String>
    where
        for<'b> &'b S: Solver,
    {
        // Solve the process for the initial structure. Has to converge!
        let s_init = self.calculate_s(y_init);
        self.excluded_solutions
            .extend(self.minlp.exclude_solutions(&s_init));
        let result = self
            .solve_nlp_with_options(y_init, s_init, &options.nlp_options)
            .expect("The optimization did not converge for the initial structure!");
        println!(
            "{:8.5} {:.5?} {}",
            result.objective.0, result.x.data.0[0], result.key
        );

        // Initialize the list of found structures in this run
        let mut new_solutions = vec![result.key.clone()];

        // objective value of the previous structure
        let mut last = result.objective.0;

        for k in 0..options.max_iter {
            // Solve for a new structure
            let (y, s) = match self.solve_milp(solver, &new_solutions, options.zero_tol) {
                Ok(result) => result,
                Err(e) => {
                    // No new structure found -> exit run
                    println!("{e}");
                    return new_solutions;
                }
            };
            // Exclude the found structure and all symmetric structures from future runs
            self.excluded_solutions
                .extend(self.minlp.exclude_solutions(&s));

            // Solve the process for the current structure
            if let Some(result) = self.solve_nlp_with_options(y, s, &options.nlp_options) {
                println!(
                    "{:8.5} {:.5?} {}",
                    result.objective.0, result.x.data.0[0], result.key
                );
                let obj = result.objective.0;
                new_solutions.push(result.key.clone());

                // Exit after at least min_iter iterations and on non-improving objective
                if obj > last && k >= options.min_iter {
                    return new_solutions;
                }
                last = obj;
            } else {
                println!("{} not converged!", self.minlp.y_to_string(&y));
            }
        }

        new_solutions
    }

    pub fn solve_ranking<S>(
        mut self,
        y_init: SMatrix<f64, N_Y1, N_Y2>,
        solver: S,
        runs: usize,
        options: &OptimizationOptions,
    ) -> Vec<OptimizationResult<N_X, N_Y1, N_Y2>>
    where
        for<'b> &'b S: Solver,
    {
        // let mut old_solutions = Vec::new();
        let key_init = self.minlp.y_to_string(&y_init);
        for k in 0..runs {
            println!("\nStarting run {}", k + 1);

            // Calculate a new set of solutions
            let new_solutions = self.solve(y_init, &solver, options);
            let new_solutions: HashSet<_> = new_solutions.into_iter().collect();

            // Calculate a sorted list of all solutions
            let mut all_solutions: Vec<_> = self.known_solutions.values().collect();
            all_solutions.sort_by(|&s1, &s2| s1.objective.0.total_cmp(&s2.objective.0));

            // Print the results
            println!("\nRanking after run {}", k + 1);
            for (k, solution) in all_solutions.into_iter().enumerate() {
                let s = self.minlp.y_to_string(&solution.y);
                let known = if solution.key == key_init {
                    "+"
                } else if new_solutions.contains(&solution.key) {
                    "*"
                } else {
                    " "
                };
                println!(
                    "{:3}{known} {:10.7} {:.5?} {s}",
                    k + 1,
                    solution.objective.0,
                    solution.x.data.0[0]
                );
            }
        }

        self.known_solutions.into_values().collect()
    }
}

#[cfg(test)]
mod test {
    use std::convert::Infallible;

    use good_lp::highs;

    use super::*;

    struct TestMINLP;

    impl MixedIntegerNonLinearProgram<2, 3, 1> for TestMINLP {
        type Error = Infallible;

        fn x_variables(&self) -> SVector<(f64, f64, f64), 2> {
            SVector::from([(0.0, 10.0, 5.0), (0.0, 10.0, 5.0)])
        }

        fn y_variables(&self) -> SVector<(i32, i32), 3> {
            SVector::from([(0, 1); 3])
        }

        fn linear_constraints(&self, y: SVector<Variable, 3>) -> Vec<LinearConstraint> {
            vec![constraint!(-y[0] - y[1] + y[2] <= 0.0)]
        }

        fn constraints(&self) -> Vec<GeneralConstraint> {
            vec![
                GeneralConstraint::Equality(1.25),
                GeneralConstraint::Equality(3.0),
                GeneralConstraint::Inequality(None, Some(1.6)),
                GeneralConstraint::Inequality(None, Some(3.0)),
            ]
        }

        fn evaluate<D: DualNum<f64> + Copy>(
            &self,
            x: SVector<D, 2>,
            y: SVector<D, 3>,
        ) -> Result<(D, Vec<D>), Self::Error> {
            let [x1, x2] = x.data.0[0];
            let [y1, y2, y3] = y.data.0[0];
            let c1 = x1 * x1 + y1;
            let c2 = x2.powf(1.5) + y2 * 1.5;
            let c4 = x1 + y1;
            let c5 = x2 * 1.333 + y2;
            let obj = x1 * 2.0 + x2 * 3.0 + y1 * 1.5 + y2 * 2.0 - y3 * 0.5;
            Ok((obj, vec![c1, c2, c4, c5]))
        }

        fn y_to_string(&self, y: &SMatrix<f64, 3, 1>) -> String {
            format!("{:.1?}", y.data.0[0])
        }
    }

    #[test]
    fn test_minlp_highs() {
        let minlp = OuterApproximation::new(&TestMINLP);
        minlp.solve_ranking(SVector::from([0.0; 3]), &highs, 1, &Default::default());
    }
}
