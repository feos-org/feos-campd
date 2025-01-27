use super::{ChemicalRecord, MolecularRepresentation, SuperMolecule};
use crate::GeneralConstraint;
use good_lp::{constraint, Constraint, Expression, Variable};
use num_dual::DualNum;
use std::collections::HashMap;

/// A combination of different molecule superstructures.
#[derive(Clone)]
pub struct Disjunction<const N: usize>(pub [SuperMolecule; N]);

impl<const N: usize> Disjunction<N> {
    pub const fn variables(&self) -> usize {
        let mut variables = 0;
        let mut k = 0;
        while k < N {
            let v = self.0[k].variables();
            if v > variables {
                variables = v;
            }
            k += 1;
        }
        variables + N
    }
}

impl<const N: usize> Disjunction<N> {
    pub fn smiles(&self, s: &[f64]) -> String {
        let (y, c) = s.split_at(s.len() - N);

        self.0
            .iter()
            .zip(c)
            .filter(|(_, &c)| c == 1.0)
            .map(|(m, _)| m.smiles(y))
            .last()
            .unwrap()
    }
}

impl<const N_Y: usize, const N: usize> MolecularRepresentation<N_Y> for Disjunction<N>
where
    SuperMolecule: MolecularRepresentation<N_Y>,
{
    fn structure_variables(&self) -> [(i32, i32); N_Y] {
        [(0, 1); N_Y]
    }

    fn constraints(&self, y: [Variable; N_Y]) -> Vec<Constraint> {
        let (y, c) = y.split_at(y.len() - N);
        let mut constraints = Vec::new();

        // disjunction
        constraints.push(constraint!(c.iter().sum::<Expression>() == 1.0));

        for (m, &c) in self.0.iter().zip(c) {
            let n_y = m.variables();
            let (used_vars, unused_vars) = y.split_at(n_y);
            let mut constr = m.constraints(used_vars);

            // unused variable constraint
            if !unused_vars.is_empty() {
                constr.push((
                    unused_vars.to_vec(),
                    vec![1.0; unused_vars.len()],
                    GeneralConstraint::Equality(0.0),
                ));
            }

            for (vars, coefs, con) in constr {
                let (lhs, rhs) = match con {
                    GeneralConstraint::Equality(eq) => (Some(eq), Some(eq)),
                    GeneralConstraint::Inequality(lhs, rhs) => (lhs, rhs),
                };
                if let Some(lhs) = lhs {
                    let big_m = coefs.iter().map(|c| c.min(0.0)).sum::<f64>() - lhs;
                    let expr: Expression = vars.iter().zip(&coefs).map(|(&v, &c)| v * c).sum();
                    constraints.push(constraint!(expr + big_m * c >= lhs + big_m));
                }
                if let Some(rhs) = rhs {
                    let big_m = coefs.iter().map(|c| c.max(0.0)).sum::<f64>() - rhs;
                    let expr: Expression = vars.iter().zip(&coefs).map(|(&v, &c)| v * c).sum();
                    constraints.push(constraint!(expr + big_m * c <= rhs + big_m));
                }
            }
        }

        constraints
    }

    fn smiles(&self, y: [f64; N_Y]) -> String {
        self.smiles(&y)
    }

    fn build_molecule<D: DualNum<f64> + Copy>(&self, s: [D; N_Y]) -> ChemicalRecord<D> {
        let (y, c) = s.split_at(s.len() - N);

        let mut groups = HashMap::new();
        let mut bonds = HashMap::new();
        for (m, c) in self.0.iter().zip(c) {
            let cr = m.build(y.to_vec());
            cr.groups
                .into_iter()
                .for_each(|(g, v)| *groups.entry(g).or_insert(D::zero()) += v * c);
            cr.bonds
                .into_iter()
                .for_each(|(b, v)| *bonds.entry(b).or_insert(D::zero()) += v * c);
        }
        ChemicalRecord { groups, bonds }
    }
}
