use super::{LinearConstraint, MolecularRepresentation, Variable};
use std::ops::{Add, Mul};

impl<
        M: MolecularRepresentation<1, ChemicalRecord = C>,
        C: Add<Output = C> + Mul<f64, Output = C>,
        const N: usize,
    > MolecularRepresentation<1> for [M; N]
{
    type ChemicalRecord = C;

    fn variables(&self) -> Vec<Variable> {
        let mut variables = self
            .iter()
            .map(|s| s.variables())
            .max_by_key(|v| v.len())
            .unwrap();
        variables.extend([Variable::binary(Some(1.0 / N as f64)); N]);
        variables
    }

    fn constraints(&self, index_vars: &[i32]) -> Vec<LinearConstraint> {
        // structure variables
        let (y, c) = index_vars.split_at(index_vars.len() - self.len());

        // disjunction
        let mut constraints = Vec::new();
        let coefs = vec![1.0; c.len()];
        constraints.push(LinearConstraint::new(c.to_vec(), coefs).eqbnd(1.0));

        for (i, m) in self.iter().enumerate() {
            let n_y = m.variables().len();
            let (used_vars, unused_vars) = y.split_at(n_y);
            let mut constr = m.constraints(used_vars);

            // unused variable constraint
            if !unused_vars.is_empty() {
                let coefs = vec![1.0; unused_vars.len()];
                constr.push(LinearConstraint::new(unused_vars.to_vec(), coefs).upbnd(0.0));
            }

            // modified constraints with big M formalism
            for constr in constr {
                let mut vars = constr.vars;
                let mut coefs = constr.coefs;
                vars.push(c[i]);
                if let Some(lobnd) = constr.lobnd {
                    let big_m = coefs.iter().map(|c| c.min(0.0)).sum::<f64>() - lobnd;
                    coefs.push(big_m);
                    constraints.push(LinearConstraint::new(vars, coefs).lobnd(lobnd + big_m));
                } else if let Some(upbnd) = constr.upbnd {
                    let big_m = coefs.iter().map(|c| c.max(0.0)).sum::<f64>() - upbnd;
                    coefs.push(big_m);
                    constraints.push(LinearConstraint::new(vars, coefs).upbnd(upbnd + big_m));
                }
            }
        }

        constraints
    }

    fn build(&self, y: &[f64]) -> [C; 1] {
        let (y, c) = y.split_at(y.len() - N);
        [self
            .iter()
            .map(|m| {
                let [cr] = m.build(y);
                cr
            })
            .zip(c.iter())
            .map(|(cr, &c)| cr * c)
            .reduce(|a, b| a + b)
            .unwrap()]
    }

    fn smiles(&self, y: &[usize]) -> [String; 1] {
        let (y, c) = y.split_at(y.len() - N);
        for (&c, m) in c.iter().zip(self.iter()) {
            if c == 1 {
                return m.smiles(y);
            }
        }
        unreachable!();
    }
}
