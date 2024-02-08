use crate::variables::{LinearConstraint, ParameterVariables, StructureVariables, Variable};
#[cfg(feature = "knitro_rs")]
use crate::{OptimizationMode, OptimizationResult};
use itertools::Itertools;
#[cfg(feature = "knitro_rs")]
use knitro_rs::{Knitro, KnitroError};
use ndarray::{Array, Array1};
#[cfg(feature = "knitro_rs")]
use std::collections::HashSet;

mod comt_camd;
mod disjunction;
mod mixture;
mod polynomial;
mod supermolecule;
pub use comt_camd::{CoMTCAMD, CoMTCAMDPropertyModel};
pub use supermolecule::{SegmentAndBondCount, SuperMolecule};

/// A generic molecular representation to be used in an [OptimizationProblem](super::OptimizationProblem).
pub trait MolecularRepresentation<const N: usize> {
    type ChemicalRecord;
    fn build(&self, y: &[f64], p: &[f64]) -> [Self::ChemicalRecord; N];

    fn structure_variables(&self) -> StructureVariables;

    fn parameter_variables(&self) -> ParameterVariables;

    fn determine_parameters(&self, y: &[f64]) -> Vec<f64>;

    fn constraints(
        &self,
        index_structure_vars: &[i32],
        index_parameter_vars: Option<&[i32]>,
    ) -> Vec<LinearConstraint>;

    fn smiles(&self, y: &[usize]) -> [String; N];

    fn generate_solutions(&self) -> Vec<Array1<i32>> {
        let mut res = Vec::new();
        let n_y = self.structure_variables().len();
        let index_vars: Vec<_> = (0..n_y as i32).collect();
        let constraints = self.constraints(&index_vars, None);
        let constraints: Vec<_> = constraints
            .into_iter()
            .flat_map(|constraint| {
                let mut con = vec![];
                if let Some(upbnd) = constraint.upbnd {
                    con.push((constraint.vars.clone(), constraint.coefs.clone(), upbnd));
                }
                if let Some(lobnd) = constraint.lobnd {
                    let coefs = constraint.coefs.iter().map(|c| -c).collect();
                    con.push((constraint.vars, coefs, -lobnd));
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
        p0: Option<&[f64]>,
        solutions: &HashSet<OptimizationResult>,
        mode: OptimizationMode,
    ) -> Result<[Vec<i32>; 2], KnitroError> {
        // define parameter variables
        let index_parameter_vars = if let OptimizationMode::FixedMolecule = mode {
            let p0 = self.determine_parameters(y0.unwrap());
            self.parameter_variables().setup_knitro(kc, Some(&p0), mode)
        } else {
            self.parameter_variables().setup_knitro(kc, p0, mode)
        }?;

        // define structure variables
        let index_structure_vars = self.structure_variables().setup_knitro(kc, y0, mode)?;

        // define constraints
        let constraints = self.constraints(&index_structure_vars, Some(&index_parameter_vars));
        for constraint in constraints {
            let c = kc.add_con()?;
            kc.add_con_linear_struct_one(c, &constraint.vars, &constraint.coefs)?;
            if let Some(lobnd) = constraint.lobnd {
                kc.set_con_lobnd(c, lobnd)?;
            }
            if let Some(upbnd) = constraint.upbnd {
                kc.set_con_upbnd(c, upbnd)?;
            }
        }

        // integer cuts
        for solution in solutions.iter() {
            let y0 = &solution.y;
            let c = kc.add_con()?;
            // let qcoefs: Vec<_> = y0.iter().map(|_| 1.0).collect();
            let lcoefs: Vec<_> = y0.iter().map(|&y0| 1.0 - 2.0 * y0 as f64).collect();
            let lbond = 1.0 - y0.iter().map(|&y0| y0 as f64).sum::<f64>();
            // kc.add_con_quadratic_struct_one(
            //     c,
            //     &index_structure_vars,
            //     &index_structure_vars,
            //     &qcoefs,
            // )?;
            kc.add_con_linear_struct_one(c, &index_structure_vars, &lcoefs)?;
            kc.set_con_lobnd(c, lbond)?;
        }
        Ok([index_structure_vars, index_parameter_vars])
    }
}
