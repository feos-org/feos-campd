use crate::molecule::{FixedMolecule, FunctionalGroup, HomoGc, Molecule, SuperMolecule};
use feos_core::python::parameter::PyChemicalRecord;
use pyo3::prelude::*;
use std::vec::Vec;

#[pyclass(name = "HomoGc")]
#[derive(Clone)]
pub struct PyHomoGc(pub HomoGc);

#[pymethods]
impl PyHomoGc {
    #[new]
    fn new(size: usize, functional_groups: Vec<(String, String)>) -> Self {
        let functional_groups = functional_groups
            .into_iter()
            .map(FunctionalGroup::from)
            .collect();
        Self(HomoGc::new(size, functional_groups))
    }

    #[getter]
    fn get_size(&self) -> usize {
        self.0.size
    }

    #[getter]
    fn get_variables(&self) -> usize {
        self.0.variables()
    }

    fn build(&self, y: Vec<f64>) -> PyChemicalRecord {
        PyChemicalRecord(self.0.build(y))
    }
}

#[pyclass(name = "SuperMolecule")]
#[derive(Clone)]
pub struct PySuperMolecule(pub SuperMolecule);

#[pymethods]
impl PySuperMolecule {
    #[new]
    fn new(size: usize, functional_groups: Vec<(String, String)>) -> Self {
        let functional_groups = functional_groups
            .into_iter()
            .map(FunctionalGroup::from)
            .collect();
        Self(SuperMolecule::new(size, functional_groups))
    }

    #[getter]
    fn get_size(&self) -> usize {
        self.0.size
    }

    #[getter]
    fn get_variables(&self) -> usize {
        self.0.variables()
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

    fn _repr_mimebundle_(
        &self,
        py: Python,
        include: Option<Vec<String>>,
        exclude: Option<Vec<String>>,
    ) -> PyResult<PyObject> {
        let fun: Py<PyAny> = PyModule::from_code(
            py,
            "def f(s, include, exclude):
                import graphviz
                return graphviz.Source(s.replace('\\\\\"', ''))._repr_mimebundle_(include, exclude)",
            "",
            "",
        )?
        .getattr("f")?
        .into();
        fun.call1(py, (self.0.indices(), include, exclude))
    }

    fn molecule(&self, y: Vec<usize>) -> PyMolecule {
        PyMolecule(self.0.molecule(y))
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

#[pyclass(name = "Molecule", unsendable)]
pub struct PyMolecule(Molecule);

#[pymethods]
impl PyMolecule {
    #[getter]
    fn get_smiles(&self) -> String {
        self.0.smiles()
    }

    fn _repr_mimebundle_(
        &self,
        py: Python,
        include: Option<Vec<String>>,
        exclude: Option<Vec<String>>,
    ) -> PyResult<PyObject> {
        let fun: Py<PyAny> = PyModule::from_code(
            py,
            "def f(s, include, exclude):
                import graphviz
                return graphviz.Source(s.replace('\\\\\"', ''))._repr_mimebundle_(include, exclude)",
            "",
            "",
        )?
        .getattr("f")?
        .into();
        fun.call1(py, (self.0.dot(), include, exclude))
    }
}
