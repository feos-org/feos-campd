use super::{
    Constraint, MolecularRepresentation, ParameterVariables, StructureVariables, Variable,
};
use std::ops::{Add, Mul};

impl<
        M: MolecularRepresentation<ChemicalRecord = [C; 1]>,
        C: Add<Output = C> + Mul<f64, Output = C>,
        const N: usize,
    > MolecularRepresentation for [M; N]
{
    type ChemicalRecord = [C; 1];

    fn structure_variables(&self) -> StructureVariables {
        let mut variables = self
            .iter()
            .map(|s| s.structure_variables())
            .max_by_key(|v| v.len())
            .unwrap();
        variables.extend([Variable::binary(); N]);
        variables
    }

    fn parameter_variables(&self) -> ParameterVariables {
        self.iter()
            .map(|s| s.parameter_variables())
            .max_by_key(|v| v.len())
            .unwrap()
    }

    fn determine_parameters(&self, y: &[f64]) -> Vec<f64> {
        let (y, c) = y.split_at(y.len() - N);
        self.iter()
            .map(|m| m.determine_parameters(y))
            .zip(c.iter())
            .map(|(p, &c)| p.iter().map(|p| p * c).collect::<Vec<_>>())
            .reduce(|a, b| a.iter().zip(b).map(|(a, b)| a + b).collect())
            .unwrap()
    }

    fn constraints(
        &self,
        index_structure_vars: &[i32],
        index_parameter_vars: Option<&[i32]>,
    ) -> Vec<Constraint> {
        // structure variables
        let (y, c) = index_structure_vars.split_at(index_structure_vars.len() - self.len());

        // disjunction
        let mut constraints = Vec::new();
        let coefs = vec![1.0; c.len()];
        constraints.push(
            Constraint::new()
                .linear_struct(c.to_vec(), coefs)
                .eqbnd(1.0),
        );

        for (i, m) in self.iter().enumerate() {
            let n_y = m.structure_variables().len();
            let (used_vars, unused_vars) = y.split_at(n_y);
            let mut constr = m.constraints(used_vars, index_parameter_vars);

            // unused variable constraint
            if !unused_vars.is_empty() {
                let coefs = vec![1.0; unused_vars.len()];
                constr.push(
                    Constraint::new()
                        .linear_struct(unused_vars.to_vec(), coefs)
                        .upbnd(0.0),
                );
            }

            // modified constraints with big M formalism
            for constr in constr {
                let mut vars = constr.lvars;
                let mut coefs = constr.lcoefs;
                vars.push(c[i]);
                if let Some(lobnd) = constr.lobnd {
                    let big_m = coefs.iter().map(|c| c.min(0.0)).sum::<f64>() - lobnd;
                    coefs.push(big_m);
                    constraints.push(
                        Constraint::new()
                            .linear_struct(vars, coefs)
                            .lobnd(lobnd + big_m),
                    );
                } else if let Some(upbnd) = constr.upbnd {
                    let big_m = coefs.iter().map(|c| c.max(0.0)).sum::<f64>() - upbnd;
                    coefs.push(big_m);
                    constraints.push(
                        Constraint::new()
                            .linear_struct(vars, coefs)
                            .upbnd(upbnd + big_m),
                    );
                }
            }
        }

        constraints
    }

    fn build(&self, y: &[f64], p: &[f64]) -> [C; 1] {
        let (y, c) = y.split_at(y.len() - N);

        [self
            .iter()
            .map(|m| {
                let [cr] = m.build(y, p);
                cr
            })
            .zip(c.iter())
            .map(|(cr, &c)| cr * c)
            .reduce(|a, b| a + b)
            .unwrap()]
    }

    fn smiles(&self, y: &[usize]) -> Vec<String> {
        let (y, c) = y.split_at(y.len() - N);
        for (&c, m) in c.iter().zip(self.iter()) {
            if c == 1 {
                return m.smiles(y);
            }
        }
        unreachable!();
    }
}
