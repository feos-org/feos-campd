use crate::CoMTCAMD;
use crate::{molecule::Disjunction, MolecularRepresentation, SuperMolecule};
use pyo3::prelude::*;
use std::collections::HashMap;

#[pyclass(name = "CoMTCAMD")]
#[derive(Clone)]
pub struct PyCoMTCAMD(pub CoMTCAMD);

#[pymethods]
impl PyCoMTCAMD {
    #[staticmethod]
    fn from_json(file: String) -> PyResult<Self> {
        Ok(Self(CoMTCAMD::from_json(file)?))
    }

    #[staticmethod]
    fn from_molecules(molecules: Vec<String>) -> Self {
        Self(CoMTCAMD::from_molecules(molecules))
    }

    #[staticmethod]
    fn from_json_molecules(file: String) -> PyResult<Self> {
        Ok(Self(CoMTCAMD::from_json_molecules(file)?))
    }

    fn get_initial_values(&self, structure: &str, groups: HashMap<&str, usize>) -> Vec<f64> {
        let groups = groups.into_iter().collect();
        self.0.get_initial_values(structure, &groups)
    }

    fn get_initial_values_molecules(&self, molecule: &str) -> Vec<f64> {
        self.0.get_initial_values_molecules(molecule)
    }
}

#[derive(Clone)]
#[allow(clippy::large_enum_variant)]
pub enum SuperMolecules {
    Single(SuperMolecule),
    All(Disjunction<SuperMolecule, 7>),
    NonAssociating(Disjunction<SuperMolecule, 5>),
}

impl MolecularRepresentation for SuperMolecules {
    fn structure_variables(&self) -> crate::StructureVariables {
        match self {
            SuperMolecules::Single(m) => m.structure_variables(),
            SuperMolecules::All(m) => m.structure_variables(),
            SuperMolecules::NonAssociating(m) => m.structure_variables(),
        }
    }

    fn feature_variables(
        &self,
        index_structure_vars: &[i32],
    ) -> indexmap::IndexMap<String, crate::variables::ExplicitVariable> {
        match self {
            SuperMolecules::Single(m) => m.feature_variables(index_structure_vars),
            SuperMolecules::All(m) => m.feature_variables(index_structure_vars),
            SuperMolecules::NonAssociating(m) => m.feature_variables(index_structure_vars),
        }
    }

    fn constraints(&self, index_structure_vars: &[i32]) -> Vec<crate::Constraint> {
        match self {
            SuperMolecules::Single(m) => m.constraints(index_structure_vars),
            SuperMolecules::All(m) => m.constraints(index_structure_vars),
            SuperMolecules::NonAssociating(m) => m.constraints(index_structure_vars),
        }
    }

    fn smiles(&self, y: &[usize]) -> String {
        match self {
            SuperMolecules::Single(m) => m.smiles(y),
            SuperMolecules::All(m) => m.smiles(y),
            SuperMolecules::NonAssociating(m) => m.smiles(y),
        }
    }
}

#[pyclass(name = "SuperMolecule")]
#[derive(Clone)]
pub struct PySuperMolecule(pub SuperMolecules);

#[pymethods]
impl PySuperMolecule {
    #[staticmethod]
    fn alkane(size: usize) -> Self {
        Self(SuperMolecules::Single(SuperMolecule::alkane(size)))
    }

    #[staticmethod]
    fn alkene(size: usize) -> Self {
        Self(SuperMolecules::Single(SuperMolecule::alkene(size)))
    }

    #[staticmethod]
    fn alkyne(size: usize) -> Self {
        Self(SuperMolecules::Single(SuperMolecule::alkyne(size)))
    }

    #[staticmethod]
    fn alcohol(size: usize) -> Self {
        Self(SuperMolecules::Single(SuperMolecule::alcohol(size)))
    }

    #[staticmethod]
    fn methylether(size: usize) -> Self {
        Self(SuperMolecules::Single(SuperMolecule::methylether(size)))
    }

    #[staticmethod]
    fn ketone(size: usize) -> Self {
        Self(SuperMolecules::Single(SuperMolecule::ketone(size)))
    }

    #[staticmethod]
    fn amine(size: usize) -> Self {
        Self(SuperMolecules::Single(SuperMolecule::amine(size)))
    }

    #[staticmethod]
    fn all(size: usize) -> Self {
        Self(SuperMolecules::All(SuperMolecule::all(size)))
    }

    #[staticmethod]
    fn non_associating(size: usize) -> Self {
        Self(SuperMolecules::NonAssociating(
            SuperMolecule::non_associating(size),
        ))
    }

    fn get_initial_values(&self, structure: &str, y0: Vec<f64>) -> Vec<f64> {
        match &self.0 {
            SuperMolecules::Single(_) => y0,
            SuperMolecules::All(molecules) => {
                SuperMolecule::get_initial_values(molecules, structure, &y0)
            }
            SuperMolecules::NonAssociating(molecules) => {
                SuperMolecule::get_initial_values(molecules, structure, &y0)
            }
        }
    }
}
