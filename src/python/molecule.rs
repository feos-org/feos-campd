use crate::molecule::{FixedMolecule, SuperMolecule};
use feos_core::python::parameter::PyChemicalRecord;
use pyo3::prelude::*;
use std::vec::Vec;

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
    fn all(size: usize) -> Vec<(String, Self)> {
        SuperMolecule::all(size)
            .into_iter()
            .map(|(s, m)| (s, Self(m)))
            .collect()
    }

    #[getter]
    fn get_variables(&self) -> usize {
        self.0.variables()
    }

    #[getter]
    fn get_size_constraint(&self) -> (Vec<usize>, usize) {
        self.0.size_constraint()
    }

    #[getter]
    fn get_functional_group_constraint(&self) -> Vec<usize> {
        self.0.functional_group_constraint()
    }

    #[getter]
    fn get_bond_constraints(&self) -> Vec<(usize, usize)> {
        self.0.bond_constraints()
    }

    #[getter]
    fn get_symmetry_constraints(&self) -> Vec<(Vec<usize>, Vec<isize>)> {
        self.0.symmetry_constraints()
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
