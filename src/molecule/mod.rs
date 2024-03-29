use crate::variables::{Constraint, ExplicitVariable, StructureVariables};
use indexmap::IndexMap;
use itertools::Itertools;
#[cfg(feature = "knitro_rs")]
use knitro_rs::{Knitro, KnitroError};
use ndarray::{Array, Array1};

mod comt_camd;
mod disjunction;
mod polynomial;
mod supermolecule;
pub use comt_camd::CoMTCAMD;
pub use disjunction::Disjunction;
pub use supermolecule::SuperMolecule;

/// A generic molecular representation to be used in an [OptimizationProblem](super::OptimizationProblem).
pub trait MolecularRepresentation {
    fn structure_variables(&self) -> StructureVariables;

    fn feature_variables(&self, index_structure_vars: &[i32])
        -> IndexMap<String, ExplicitVariable>;

    fn evaluate_feature_variables(&self, structure_vars: &[f64]) -> IndexMap<String, f64> {
        let index_structure_vars: Vec<_> = (0..structure_vars.len() as i32).collect();
        self.feature_variables(&index_structure_vars)
            .into_iter()
            .map(|(k, v)| (k, v.evaluate(structure_vars)))
            .collect()
    }

    fn constraints(&self, index_structure_vars: &[i32]) -> Vec<Constraint>;

    fn smiles(&self, y: &[usize]) -> String;

    fn generate_solutions(&self) -> Vec<Array1<i32>> {
        let mut res = Vec::new();
        let n_y = self.structure_variables().len();
        let index_vars: Vec<_> = (0..n_y as i32).collect();
        let constraints = self.constraints(&index_vars);
        let constraints: Vec<_> = constraints
            .into_iter()
            .flat_map(|constraint| {
                let mut con = vec![];
                if let Some(upbnd) = constraint.upbnd {
                    con.push((constraint.lvars.clone(), constraint.lcoefs.clone(), upbnd));
                }
                if let Some(lobnd) = constraint.lobnd {
                    let coefs = constraint.lcoefs.iter().map(|c| -c).collect();
                    con.push((constraint.lvars, coefs, -lobnd));
                }
                con
            })
            .collect();
        let mut a = Array::zeros((constraints.len(), n_y));
        let mut b = Array::zeros(constraints.len());
        for (i, (vars, coefs, upbnd)) in constraints.into_iter().enumerate() {
            for (j, c) in vars.into_iter().zip(coefs) {
                a[[i, j as usize]] = c as i32;
            }
            b[i] = upbnd as i32;
        }

        for y in vec![0..=1; self.structure_variables().len()]
            .into_iter()
            .multi_cartesian_product()
        {
            let y = Array::from_vec(y);
            if b.iter().zip(a.dot(&y)).all(|(&b, x)| x <= b) {
                res.push(y)
            }
        }
        res
    }

    #[cfg(feature = "knitro_rs")]
    fn setup_knitro(
        &self,
        kc: &Knitro,
        y0: Option<&[f64]>,
        target: bool,
    ) -> Result<(Vec<i32>, IndexMap<String, i32>), KnitroError> {
        // define structure variables
        let index_structure_vars = self.structure_variables().setup_knitro(kc, y0, target)?;

        // define feature variables
        let index_feature_vars = self
            .feature_variables(&index_structure_vars)
            .into_iter()
            .map(|(k, f)| Ok((k, f.setup_knitro(kc)?)))
            .collect::<Result<IndexMap<_, _>, _>>()?;

        // define constraints
        self.constraints(&index_structure_vars)
            .iter()
            .try_for_each(|c| c.setup_knitro(kc))?;

        Ok((index_structure_vars, index_feature_vars))
    }
}
