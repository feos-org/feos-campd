use super::{MolecularRepresentation, OptimizationProblem, OptimizationResult, ProcessModel};
use feos_core::parameter::ParameterError;
use feos_core::python::parameter::{PyChemicalRecord, PyIdentifier};
use feos_core::{impl_json_handling, Contributions, Verbosity};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::wrap_pymodule;
use quantity::python::__PYO3_PYMODULE_DEF_QUANTITY;

mod molecule;
mod orc;
mod process;
mod property;
use molecule::{PyFixedMolecule, PySuperMolecule};
use orc::{PyOrganicRankineCycle, PyOrganicRankineCycleSuperStructure};
use process::PyUtility;
use property::{PyPhaseDiagram, PyPhaseEquilibrium, PyPropertyModel, PyState};

#[pyclass(name = "OptimizationProblem")]
pub struct PyOptimizationProblem(OptimizationProblem);

#[pymethods]
impl PyOptimizationProblem {
    #[new]
    fn new(
        molecule: &PyAny,
        property_model: Option<PyPropertyModel>,
        process: Option<&PyAny>,
    ) -> PyResult<Self> {
        let molecule = if let Ok(molecule) = molecule.extract::<PySuperMolecule>() {
            MolecularRepresentation::SuperMolecule(molecule.0)
        } else if let Ok(molecule) = molecule.extract::<PyFixedMolecule>() {
            MolecularRepresentation::FixedMolecule(molecule.0)
        } else {
            return Err(PyValueError::new_err(
                "Parameter `molecule` has an invalid type.",
            ));
        };
        let process = process
            .map(|process| {
                if let Ok(process) = process.extract::<PyOrganicRankineCycle>() {
                    Ok(ProcessModel::OrganicRankineCycle(process.0))
                } else if let Ok(process) = process.extract::<PyOrganicRankineCycleSuperStructure>()
                {
                    Ok(ProcessModel::OrganicRankineCycleSuperStructure(process.0))
                } else {
                    Err(PyValueError::new_err(
                        "Parameter `process` has an invalid type.",
                    ))
                }
            })
            .transpose()?;
        Ok(Self(OptimizationProblem::new(
            molecule,
            property_model.map(|p| p.0),
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
            MolecularRepresentation::FixedMolecule(molecule) => {
                PyFixedMolecule(molecule.clone()).into_py(py)
            }
        }
    }

    #[getter]
    fn get_property_model(&self) -> Option<PyPropertyModel> {
        self.0.property_model.clone().map(PyPropertyModel)
    }

    #[getter]
    fn get_process(&self, py: Python) -> Option<PyObject> {
        self.0.process.clone().map(|process| match process {
            ProcessModel::OrganicRankineCycle(process) => {
                PyOrganicRankineCycle(process).into_py(py)
            }
            ProcessModel::OrganicRankineCycleSuperStructure(process) => {
                PyOrganicRankineCycleSuperStructure(process).into_py(py)
            }
        })
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

    fn __repr__(&self) -> PyResult<String> {
        Ok(self.0.to_string())
    }
}

#[pymodule]
pub fn feos_campd(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyOptimizationProblem>()?;
    m.add_class::<PyOptimizationResult>()?;
    m.add_class::<PyOrganicRankineCycle>()?;
    m.add_class::<PyOrganicRankineCycleSuperStructure>()?;
    m.add_class::<PyUtility>()?;
    m.add_class::<PySuperMolecule>()?;
    m.add_class::<PyFixedMolecule>()?;
    m.add_class::<PyPropertyModel>()?;

    m.add_class::<PyIdentifier>()?;
    m.add_class::<Verbosity>()?;
    m.add_class::<Contributions>()?;
    m.add_class::<PyChemicalRecord>()?;

    m.add_class::<PyState>()?;
    m.add_class::<PyPhaseDiagram>()?;
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
