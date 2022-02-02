#![warn(clippy::all)]
use feos_core::parameter::ParameterError;
use serde::{Deserialize, Serialize};
use std::fmt;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;

mod molecule;
mod process;
mod property;
use molecule::MolecularRepresentation;
pub use molecule::{FixedMolecule, FunctionalGroup, HomoGc, Molecule, SuperMolecule};
pub use process::{
    Equipment, OrganicRankineCycle, Process, ProcessModel, ProcessPlot, ProcessPoint, ProcessState,
    ProcessStep, StatePoint,
};
pub use property::{EquationsOfState, PropertyModel};

#[cfg(feature = "python")]
mod python;

#[derive(Serialize, Deserialize)]
pub struct OptimizationProblem {
    pub molecule: MolecularRepresentation,
    pub property_model: PropertyModel,
    pub process: ProcessModel,
    pub solutions: Vec<OptimizationResult>,
}

impl OptimizationProblem {
    pub fn new(
        molecule: MolecularRepresentation,
        property_model: PropertyModel,
        process: ProcessModel,
    ) -> Self {
        Self {
            molecule,
            property_model,
            process,
            solutions: vec![],
        }
    }

    pub fn add_solution(&mut self, solutions: OptimizationResult) {
        self.solutions.push(solutions);
    }

    pub fn from_json<P: AsRef<Path>>(file: P) -> Result<Self, ParameterError> {
        Ok(serde_json::from_reader(BufReader::new(File::open(file)?))?)
    }

    pub fn to_json<P: AsRef<Path>>(&self, file: P) -> Result<(), ParameterError> {
        Ok(serde_json::to_writer_pretty(
            BufWriter::new(File::create(file)?),
            self,
        )?)
    }
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct OptimizationResult {
    pub target: f64,
    pub y: Vec<usize>,
    pub x: Vec<f64>,
}

impl OptimizationResult {
    pub fn new(target: f64, y: Vec<usize>, x: Vec<f64>) -> Self {
        Self { target, y, x }
    }
}

impl fmt::Display for OptimizationResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}
