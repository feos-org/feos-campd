use crate::molecule::{FixedMolecule, MolecularRepresentation, SuperMolecule};
use feos_core::python::parameter::PyChemicalRecord;
use pyo3::prelude::*;
use std::collections::HashMap;

#[pyclass(name = "SuperMolecule")]
#[derive(Clone)]
pub struct PySuperMolecule(pub SuperMolecule);

#[pymethods]
impl PySuperMolecule {
    #[new]
    fn new(size: usize) -> Self {
        Self(SuperMolecule::new(size, None))
    }

    #[staticmethod]
    fn alkane(size: usize) -> Self {
        Self(SuperMolecule::alkane(size))
    }

    #[staticmethod]
    fn alkene(size: usize) -> Self {
        Self(SuperMolecule::alkene(size))
    }

    #[staticmethod]
    fn alkyne(size: usize) -> Self {
        Self(SuperMolecule::alkyne(size))
    }

    #[staticmethod]
    fn alcohol(size: usize) -> Self {
        Self(SuperMolecule::alcohol(size))
    }

    #[staticmethod]
    fn methylether(size: usize) -> Self {
        Self(SuperMolecule::methylether(size))
    }

    #[staticmethod]
    fn ketone(size: usize) -> Self {
        Self(SuperMolecule::ketone(size))
    }

    #[staticmethod]
    fn amine(size: usize) -> Self {
        Self(SuperMolecule::amine(size))
    }

    #[staticmethod]
    fn all(size: usize) -> HashMap<String, Self> {
        SuperMolecule::all(size)
            .into_iter()
            .map(|(s, m)| (s, Self(m)))
            .collect()
    }

    #[getter]
    fn get_variables(&self) -> usize {
        self.0.variables()
    }

    fn build(&self, y: Vec<f64>) -> PyChemicalRecord {
        PyChemicalRecord(self.0.build(y))
    }

    fn smiles(&self, y: Vec<usize>) -> String {
        self.0.smiles(y)
    }
}

#[pyclass(name = "FixedMolecule")]
#[derive(Clone)]
pub struct PyFixedMolecule(pub FixedMolecule);

#[pymethods]
impl PyFixedMolecule {
    #[new]
    fn new(components: Vec<String>) -> Self {
        Self(FixedMolecule::new(components))
    }

    #[getter]
    fn get_variables(&self) -> usize {
        self.0.variables()
    }

    fn build(&self, y: Vec<f64>) -> PyChemicalRecord {
        PyChemicalRecord(self.0.build(y))
    }
}
