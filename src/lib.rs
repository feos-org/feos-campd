#![warn(clippy::all)]
#![allow(clippy::too_many_arguments)]
use feos_core::parameter::ParameterError;
use serde::{Deserialize, Serialize};
use std::fmt;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;

mod molecule;
pub mod process;
mod property;
pub use molecule::{
    CoMTCAMD, FixedMolecule, GroupCount, MolecularRepresentation, SegmentAndBondCount,
    SuperMolecule, SuperMoleculeDisjunct,
};
pub use property::{GcPcSaftPropertyModel, PcSaftPropertyModel, PropertyModel};

#[cfg(feature = "knitro_rs")]
pub mod knitro;

/// A full optimization problem consisting of a [MolecularRepresentation], a [PropertyModel], and a [ProcessModel](process::ProcessModel).
#[derive(Serialize, Deserialize)]
pub struct OptimizationProblem<M, R, P> {
    pub molecules: Vec<M>,
    pub property_model: R,
    pub process: P,
    pub solutions: Vec<OptimizationResult>,
}

impl<M, R, P> OptimizationProblem<M, R, P> {
    pub fn new(molecules: Vec<M>, property_model: R, process: P) -> Self {
        Self {
            molecules,
            property_model,
            process,
            solutions: vec![],
        }
    }

    pub fn new_pure(molecule: M, property_model: R, process: P) -> Self {
        Self::new(vec![molecule], property_model, process)
    }

    pub fn add_solution(&mut self, solutions: OptimizationResult) {
        self.solutions.push(solutions);
    }

    pub fn from_json<FP: AsRef<Path>>(file: FP) -> Result<Self, ParameterError>
    where
        for<'a> M: Deserialize<'a> + Default,
        for<'a> R: Deserialize<'a> + Default,
        for<'a> P: Deserialize<'a> + Default,
    {
        Ok(serde_json::from_reader(BufReader::new(File::open(file)?))?)
    }

    pub fn to_json<FP: AsRef<Path>>(&self, file: FP) -> Result<(), ParameterError>
    where
        M: Serialize,
        R: Serialize,
        P: Serialize,
    {
        Ok(serde_json::to_writer_pretty(
            BufWriter::new(File::create(file)?),
            self,
        )?)
    }
}

/// The output of a single optimization consisting of the target, structure variables, and process variables.
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
