use super::{
    MolecularRepresentation, OptimizationProblem, OptimizationResult, OrganicRankineCycle, Process,
    ProcessModel,
};
use feos_core::impl_json_handling;
use feos_core::parameter::ParameterError;
use feos_core::python::parameter::{PyChemicalRecord, PyIdentifier};
use feos_core::python::{PyContributions, PyVerbosity};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::wrap_pymodule;
use quantity::python::PyInit_quantity;

mod molecule;
mod process;
mod property;
use molecule::{PyFixedMolecule, PyHomoGc, PySuperMolecule};
use process::PyOrganicRankineCycle;
use property::{PyPhaseDiagramPure, PyPhaseEquilibrium, PyPropertyModel, PyState};

#[pyclass(name = "OptimizationProblem")]
pub struct PyOptimizationProblem(OptimizationProblem);

#[pymethods]
impl PyOptimizationProblem {
    #[new]
    fn new(molecule: &PyAny, property_model: PyPropertyModel, process: &PyAny) -> PyResult<Self> {
        let molecule = if let Ok(molecule) = molecule.extract::<PySuperMolecule>() {
            MolecularRepresentation::SuperMolecule(molecule.0)
        } else if let Ok(molecule) = molecule.extract::<PyHomoGc>() {
            MolecularRepresentation::HomoGc(molecule.0)
        } else if let Ok(molecule) = molecule.extract::<PyFixedMolecule>() {
            MolecularRepresentation::FixedMolecule(molecule.0)
        } else {
            return Err(PyValueError::new_err(
                "Parameter `molecule` has an invalid type.",
            ));
        };
        let process = if let Ok(process) = process.extract::<PyOrganicRankineCycle>() {
            ProcessModel::OrganicRankineCycle(process.0)
        } else {
            return Err(PyValueError::new_err(
                "Parameter `process` has an invalid type.",
            ));
        };
        Ok(Self(OptimizationProblem::new(
            molecule,
            property_model.0,
            process,
        )))
    }

    fn add_solution(&mut self, solution: PyOptimizationResult) {
        self.0.add_solution(solution.0)
    }

    #[getter]
    fn get_molecule(&self, py: Python) -> PyObject {
        match &self.0.molecule {
            MolecularRepresentation::SuperMolecule(molecule) => {
                PySuperMolecule(molecule.clone()).into_py(py)
            }
            MolecularRepresentation::HomoGc(molecule) => PyHomoGc(molecule.clone()).into_py(py),
            MolecularRepresentation::FixedMolecule(molecule) => {
                PyFixedMolecule(molecule.clone()).into_py(py)
            }
        }
    }

    #[getter]
    fn get_property_model(&self) -> PyPropertyModel {
        PyPropertyModel(self.0.property_model.clone())
    }

    #[getter]
    fn get_process(&self, py: Python) -> PyObject {
        match &self.0.process {
            ProcessModel::OrganicRankineCycle(process) => {
                PyOrganicRankineCycle(process.clone()).into_py(py)
            }
        }
    }

    #[getter]
    fn get_solutions(&self) -> Vec<PyOptimizationResult> {
        self.0
            .solutions
            .iter()
            .map(|s| PyOptimizationResult(s.clone()))
            .collect()
    }

    #[staticmethod]
    fn from_json(file: &str) -> Result<Self, ParameterError> {
        Ok(Self(OptimizationProblem::from_json(file)?))
    }

    fn to_json(&self, file: &str) -> Result<(), ParameterError> {
        self.0.to_json(file)
    }
}

impl_json_handling!(PyOptimizationProblem);

#[pyclass(name = "OptimizationResult")]
#[derive(Clone)]
pub struct PyOptimizationResult(OptimizationResult);

#[pymethods]
impl PyOptimizationResult {
    #[new]
    fn new(target: f64, y: Vec<usize>, x: Vec<f64>) -> Self {
        Self(OptimizationResult::new(target, y, x))
    }

    #[getter]
    fn get_target(&self) -> f64 {
        self.0.target
    }

    #[getter]
    fn get_y(&self) -> Vec<usize> {
        self.0.y.clone()
    }

    #[getter]
    fn get_x(&self) -> Vec<f64> {
        self.0.x.clone()
    }
}

#[pyproto]
impl pyo3::class::basic::PyObjectProtocol for PyOptimizationResult {
    fn __repr__(&self) -> PyResult<String> {
        Ok(self.0.to_string())
    }
}

#[pymodule]
pub fn feos_campd(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyOptimizationProblem>()?;
    m.add_class::<PyOptimizationResult>()?;
    m.add_class::<PyOrganicRankineCycle>()?;
    m.add_class::<PySuperMolecule>()?;
    m.add_class::<PyHomoGc>()?;
    m.add_class::<PyFixedMolecule>()?;
    m.add_class::<PyPropertyModel>()?;

    m.add_class::<PyIdentifier>()?;
    m.add_class::<PyVerbosity>()?;
    m.add_class::<PyContributions>()?;
    m.add_class::<PyChemicalRecord>()?;

    m.add_class::<PyState>()?;
    m.add_class::<PyPhaseDiagramPure>()?;
    m.add_class::<PyPhaseEquilibrium>()?;

    m.add_wrapped(wrap_pymodule!(quantity))?;

    py.run(
        "\
import sys
sys.modules['feos_campd.si'] = quantity
    ",
        None,
        Some(m.dict()),
    )?;
    Ok(())
}
