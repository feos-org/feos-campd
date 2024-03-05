#![warn(clippy::all)]
#![allow(clippy::too_many_arguments)]
use feos::core::parameter::ParameterError;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::fmt;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{BufReader, BufWriter};
use std::marker::PhantomData;
use std::path::Path;

mod molecule;
pub mod process;
mod property;
mod variables;
pub use molecule::{
    CoMTCAMD, CoMTCAMDBinary, MolecularRepresentation, SegmentAndBondCount, SuperMolecule,
};
pub use property::{GcPcSaftPropertyModel, PcSaftPropertyModel, PropertyModel};
pub use variables::{Constraint, ProcessVariables, StructureVariables, Variable, Variables};

#[cfg(feature = "knitro_rs")]
mod solver;
#[cfg(feature = "knitro_rs")]
pub use solver::OuterApproximationAlgorithm;

#[cfg(feature = "python")]
mod python;

#[derive(Clone, Copy)]
pub enum OptimizationMode {
    FixedMolecule,
    Gradients,
    Target,
    MolecularDesign,
    Feasibility,
}

/// A full optimization problem consisting of a [MolecularRepresentation], a [PropertyModel], and a [ProcessModel](process::ProcessModel).
#[derive(Serialize, Deserialize)]
pub struct OptimizationProblem<E, M, R, P> {
    eos: PhantomData<E>,
    pub molecules: M,
    pub property_model: R,
    pub process: P,
    pub solutions: HashSet<OptimizationResult>,
}

impl<E, M, R, P> OptimizationProblem<E, M, R, P> {
    pub fn new(molecules: M, property_model: R, process: P) -> Self {
        Self {
            eos: PhantomData,
            molecules,
            property_model,
            process,
            solutions: HashSet::new(),
        }
    }
}

impl<E, M, R, P> OptimizationProblem<E, M, R, P> {
    pub fn from_json<FP: AsRef<Path>>(file: FP) -> Result<Self, ParameterError>
    where
        for<'a> M: Deserialize<'a>,
        for<'a> R: Deserialize<'a>,
        for<'a> P: Deserialize<'a>,
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
    pub smiles: Vec<String>,
    pub x: Vec<f64>,
    pub y: Vec<usize>,
    pub p: Vec<f64>,
}

impl OptimizationResult {
    pub fn new(target: f64, smiles: Vec<String>, x: Vec<f64>, y: Vec<usize>, p: Vec<f64>) -> Self {
        Self {
            target,
            smiles,
            x,
            y,
            p,
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
