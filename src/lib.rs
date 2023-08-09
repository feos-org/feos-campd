#![warn(clippy::all)]
#![allow(clippy::too_many_arguments)]
use feos_core::parameter::ParameterError;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;

mod molecule;
pub mod process;
mod property;
pub use molecule::{
    CoMTCAMD, FixedMolecule, GroupCount, MolecularRepresentation, SegmentAndBondCount,
    SuperMolecule,
};
pub use property::{
    GcPcSaftPropertyModel, PcSaftFixedPropertyModel, PcSaftPropertyModel, PropertyModel,
};

#[cfg(feature = "knitro_rs")]
pub mod knitro;

/// A full optimization problem consisting of a [MolecularRepresentation], a [PropertyModel], and a [ProcessModel](process::ProcessModel).
#[derive(Serialize, Deserialize)]
pub struct OptimizationProblem<M, R, P> {
    pub molecule: M,
    pub property_model: R,
    pub process: P,
    pub solutions: Vec<OptimizationResult>,
}

impl<M, R, P> OptimizationProblem<M, R, P> {
    pub fn new(molecule: M, property_model: R, process: P) -> Self {
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

/// An optimization problem that can handle multiple molecular representations simultaneously.
#[derive(Serialize, Deserialize)]
pub struct MetaOptimizationProblem<R, P> {
    pub molecules: HashMap<String, SuperMolecule>,
    pub property_model: R,
    pub process: P,
    pub solutions: Vec<(String, OptimizationResult)>,
    pub candidates: HashMap<String, Vec<OptimizationResult>>,
}

impl<R, P> MetaOptimizationProblem<R, P> {
    pub fn new(molecule_size: usize, property_model: R, process: P) -> Self {
        Self {
            molecules: SuperMolecule::all(molecule_size),
            property_model,
            process,
            solutions: Vec::new(),
            candidates: HashMap::new(),
        }
    }

    pub fn best_candidate(&self) -> String {
        let (candidate, _) = self
            .candidates
            .iter()
            .min_by(|(_, c1), (_, c2)| {
                c1.last()
                    .unwrap()
                    .target
                    .partial_cmp(&c2.last().unwrap().target)
                    .unwrap()
            })
            .unwrap();
        candidate.clone()
    }

    pub fn update_candidates(&mut self, chemical: &String, solution: OptimizationResult) {
        let best = self.candidates[chemical].last().unwrap().clone();
        self.candidates.get_mut(chemical).unwrap().push(solution);
        self.solutions.push((chemical.clone(), best))
    }

    pub fn from_json<FP: AsRef<Path>>(file: FP) -> Result<Self, ParameterError>
    where
        for<'a> R: Deserialize<'a>,
        for<'a> P: Deserialize<'a>,
    {
        Ok(serde_json::from_reader(BufReader::new(File::open(file)?))?)
    }

    pub fn to_json<FP: AsRef<Path>>(&self, file: FP) -> Result<(), ParameterError>
    where
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
