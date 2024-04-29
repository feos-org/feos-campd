#![warn(clippy::all)]
#![allow(clippy::too_many_arguments)]
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::hash::{Hash, Hasher};
use std::sync::Arc;
use std::{array, fmt};

mod molecule;
pub mod process;
mod property;
mod variables;
pub use molecule::{CoMTCAMD, MolecularRepresentation, SuperMolecule};
#[doc(inline)]
pub use process::ProcessModel;
pub use property::{GcPcSaftPropertyModel, PcSaftPropertyModel, PropertyModel};
pub use variables::{
    Constraint, ExplicitVariable, ProcessVariables, StructureVariables, Variable, Variables,
};

#[cfg(feature = "knitro_rs")]
mod solver;
#[cfg(feature = "knitro_rs")]
pub use solver::OuterApproximationAlgorithm;

#[cfg(feature = "python")]
mod python;

/// A full optimization problem consisting of a [MolecularRepresentation], a [PropertyModel], and a [ProcessModel](process::ProcessModel).
// #[derive(Serialize, Deserialize)]
pub struct OptimizationProblem<M, R, P, const N: usize> {
    pub molecules: [M; N],
    pub property_model: R,
    pub process: P,
    pub solutions: HashSet<OptimizationResult>,
}

impl<M, R, P, const N: usize> OptimizationProblem<M, R, P, N> {
    pub fn new(molecules: [M; N], property_model: R, process: P) -> Self {
        Self {
            molecules,
            property_model,
            process,
            solutions: HashSet::new(),
        }
    }

    pub fn smiles(&self, y: &[Vec<f64>; N]) -> (Vec<usize>, Vec<String>)
    where
        M: MolecularRepresentation,
    {
        let (y_usize, smiles): (Vec<_>, Vec<_>) = y
            .iter()
            .zip(&self.molecules)
            .map(|(y, m)| {
                let y: Vec<_> = y.iter().map(|y| y.round() as usize).collect();
                let smiles = m.smiles(&y);
                (y, smiles)
            })
            .unzip();
        let y_usize = y_usize.concat();
        (y_usize, smiles)
    }

    pub fn build_eos(&self, result: &OptimizationResult) -> Arc<R::EquationOfState>
    where
        M: MolecularRepresentation,
        R: PropertyModel<N>,
    {
        let mut y_iter = result.y.iter().map(|&y| y as f64);
        let f = array::from_fn(|i| {
            let y = vec![y_iter.next().unwrap(); self.molecules[i].structure_variables().len()];
            self.molecules[i].evaluate_feature_variables(&y)
        });
        let p = self.property_model.evaluate_parameter_variables(&f);
        self.property_model.build_eos(&p)
    }
}

/// The output of a single optimization consisting of the target, structure variables, and process variables.
#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct OptimizationResult {
    pub target: f64,
    pub smiles: Vec<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    #[serde(default)]
    pub x: Vec<f64>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    #[serde(default)]
    pub y: Vec<usize>,
}

impl OptimizationResult {
    pub fn new(target: f64, smiles: Vec<String>, x: Vec<f64>, y: Vec<usize>) -> Self {
        Self {
            target,
            smiles,
            x,
            y,
        }
    }
}

impl fmt::Display for OptimizationResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl PartialEq for OptimizationResult {
    fn eq(&self, other: &Self) -> bool {
        self.y == other.y
    }
}

impl Eq for OptimizationResult {}

impl Hash for OptimizationResult {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.y.hash(state);
    }
}
