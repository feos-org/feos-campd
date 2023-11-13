use serde::{Deserialize, Serialize};

mod comt_camd;
mod disjunction;
mod polynomial;
mod supermolecule;
pub use comt_camd::{CoMTCAMD, GroupCount};
pub use supermolecule::{SegmentAndBondCount, SuperMolecule};

#[derive(Clone, Copy)]
pub struct Variable {
    pub init: f64,
    pub lobnd: f64,
    pub upbnd: f64,
}

impl Variable {
    pub fn binary(init: Option<f64>) -> Self {
        Self {
            init: init.unwrap_or(0.5),
            lobnd: 0.0,
            upbnd: 1.0,
        }
    }

    pub fn integer(init: Option<f64>, lobnd: usize, upbnd: usize) -> Self {
        let lobnd = lobnd as f64;
        let upbnd = upbnd as f64;
        Self {
            init: init.unwrap_or(0.5 * (lobnd + upbnd)),
            lobnd,
            upbnd,
        }
    }
}

pub struct LinearConstraint {
    pub vars: Vec<i32>,
    pub coefs: Vec<f64>,
    pub lobnd: Option<f64>,
    pub upbnd: Option<f64>,
}

impl LinearConstraint {
    pub fn new(vars: Vec<i32>, coefs: Vec<f64>) -> Self {
        Self {
            vars,
            coefs,
            lobnd: None,
            upbnd: None,
        }
    }

    pub fn lobnd(mut self, lobnd: f64) -> Self {
        self.lobnd = Some(lobnd);
        self
    }

    pub fn upbnd(mut self, upbnd: f64) -> Self {
        self.upbnd = Some(upbnd);
        self
    }

    pub fn eqbnd(self, eqbnd: f64) -> Self {
        self.lobnd(eqbnd).upbnd(eqbnd)
    }
}

/// A generic molecular representation to be used in an [OptimizationProblem](super::OptimizationProblem).
pub trait MolecularRepresentation<const N: usize> {
    type ChemicalRecord;
    fn build(&self, y: &[f64]) -> [Self::ChemicalRecord; N];

    fn variables(&self) -> Vec<Variable>;

    fn constraints(&self, index_vars: &[i32]) -> Vec<LinearConstraint>;

    fn smiles(&self, y: &[usize]) -> [String; N];
}

/// Molecular representation for a fixed molecule.
#[derive(Clone, Serialize, Deserialize)]
pub struct FixedMolecule<C> {
    smiles: String,
    chemical_record: C,
}

impl<C> FixedMolecule<C> {
    pub fn new(smiles: String, chemical_record: C) -> Self {
        Self {
            smiles,
            chemical_record,
        }
    }
}

impl<C: Clone> MolecularRepresentation<1> for FixedMolecule<C> {
    type ChemicalRecord = C;

    fn variables(&self) -> Vec<Variable> {
        vec![]
    }

    fn constraints(&self, _: &[i32]) -> Vec<LinearConstraint> {
        vec![]
    }

    fn build(&self, _: &[f64]) -> [C; 1] {
        [self.chemical_record.clone()]
    }

    fn smiles(&self, _: &[usize]) -> [String; 1] {
        [self.smiles.clone()]
    }
}

impl<
        C,
        M1: MolecularRepresentation<1, ChemicalRecord = C>,
        M2: MolecularRepresentation<1, ChemicalRecord = C>,
    > MolecularRepresentation<2> for (M1, M2)
{
    type ChemicalRecord = C;

    fn build(&self, y: &[f64]) -> [Self::ChemicalRecord; 2] {
        let (y1, y2) = y.split_at(self.0.variables().len());
        let [cr1] = self.0.build(y1);
        let [cr2] = self.1.build(y2);
        [cr1, cr2]
    }

    fn variables(&self) -> Vec<Variable> {
        let mut variables = self.0.variables();
        variables.append(&mut self.1.variables());
        variables
    }

    fn constraints(&self, index_vars: &[i32]) -> Vec<LinearConstraint> {
        let (index_vars1, index_vars2) = index_vars.split_at(self.0.variables().len());
        let mut constraints = self.0.constraints(index_vars1);
        constraints.append(&mut self.1.constraints(index_vars2));
        constraints
    }

    fn smiles(&self, y: &[usize]) -> [String; 2] {
        let (y1, y2) = y.split_at(self.0.variables().len());
        let [smiles1] = self.0.smiles(y1);
        let [smiles2] = self.1.smiles(y2);
        [smiles1, smiles2]
    }
}
