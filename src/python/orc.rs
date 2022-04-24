use super::process::{PyProcess, PyUtility};
use super::property::PyEquationsOfState;
use crate::process::{OrganicRankineCycle, ProcessModel};
use feos_core::parameter::ParameterError;
use pyo3::prelude::*;
use quantity::python::*;

#[pyclass(name = "OrganicRankineCycle")]
#[derive(Clone)]
pub struct PyOrganicRankineCycle(pub OrganicRankineCycle);

#[allow(clippy::too_many_arguments)]
#[pymethods]
impl PyOrganicRankineCycle {
    #[new]
    fn new(
        heat_source: PyUtility,
        isentropic_turbine_efficiency: f64,
        isentropic_pump_efficiency: f64,
        min_abs_pressure: PySINumber,
        min_rel_pressure: f64,
        max_abs_pressure: PySINumber,
        max_red_pressure: f64,
        cooling: PyUtility,
    ) -> Self {
        Self(OrganicRankineCycle::new(
            heat_source.0,
            isentropic_turbine_efficiency,
            isentropic_pump_efficiency,
            min_abs_pressure.into(),
            min_rel_pressure,
            max_abs_pressure.into(),
            max_red_pressure,
            cooling.0,
        ))
    }

    #[getter]
    fn get_variables(&self) -> Vec<[Option<f64>; 2]> {
        self.0.variables()
    }

    #[getter]
    fn get_binary_variables(&self) -> usize {
        self.0.binary_variables()
    }

    #[getter]
    fn get_constraints(&self) -> Vec<[Option<f64>; 3]> {
        self.0.constraints()
    }

    fn solve(
        &self,
        eos: PyEquationsOfState,
        x: Vec<f64>,
        py: Python,
    ) -> PyResult<(PyObject, f64, Vec<f64>)> {
        let (process, target, constraints) = self.0.solve(&eos.0, &x)?;
        Ok((PyProcess(process).into_py(py), target, constraints))
    }

    #[staticmethod]
    #[pyo3(text_signature = "(file)")]
    fn from_json(file: &str) -> Result<Self, ParameterError> {
        Ok(Self(OrganicRankineCycle::from_json(file)?))
    }

    #[staticmethod]
    #[pyo3(text_signature = "(json)")]
    fn from_json_str(json: &str) -> Result<Self, ParameterError> {
        Ok(Self(serde_json::from_str(json)?))
    }

    fn to_json_str(&self) -> Result<String, ParameterError> {
        Ok(serde_json::to_string(&self.0)?)
    }
}

// #[pyclass(name = "OrganicRankineCycleSuperStructure")]
// #[derive(Clone)]
// pub struct PyOrganicRankineCycleSuperStructure(pub OrganicRankineCycleSuperStructure);

// #[allow(clippy::too_many_arguments)]
// #[pymethods]
// impl PyOrganicRankineCycleSuperStructure {
//     #[new]
//     fn new(
//         heat_source: PyUtility,
//         min_heat_source_outlet_temperature: PySINumber,
//         isentropic_turbine_efficiency: f64,
//         isentropic_pump_efficiency: f64,
//         min_abs_pressure: PySINumber,
//         min_rel_pressure: f64,
//         max_abs_pressure: PySINumber,
//         max_red_pressure: f64,
//         cooling: PyUtility,
//         min_temperature_difference_recuperator: PySINumber,
//         min_heat_transfer_rate_recuperator: PySINumber,
//     ) -> Self {
//         Self(OrganicRankineCycleSuperStructure::new(
//             heat_source.0,
//             min_heat_source_outlet_temperature.into(),
//             isentropic_turbine_efficiency,
//             isentropic_pump_efficiency,
//             min_abs_pressure.into(),
//             min_rel_pressure,
//             max_abs_pressure.into(),
//             max_red_pressure,
//             cooling.0,
//             min_temperature_difference_recuperator.into(),
//             min_heat_transfer_rate_recuperator.into(),
//         ))
//     }

//     #[getter]
//     fn get_variables(&self) -> Vec<[Option<f64>; 2]> {
//         self.0.variables()
//     }

//     #[getter]
//     fn get_binary_variables(&self) -> usize {
//         self.0.binary_variables()
//     }

//     #[getter]
//     fn get_constraints(&self) -> Vec<[Option<f64>; 3]> {
//         self.0.constraints()
//     }

//     fn solve(
//         &self,
//         eos: PyEquationsOfState,
//         x: Vec<f64>,
//         py: Python,
//     ) -> PyResult<(PyObject, f64, Vec<f64>)> {
//         let (process, target, constraints) = self.0.solve(&eos.0, &x)?;
//         Ok((PyProcess(process).into_py(py), target, constraints))
//     }

//     #[staticmethod]
//     #[pyo3(text_signature = "(file)")]
//     fn from_json(file: &str) -> Result<Self, ParameterError> {
//         Ok(Self(OrganicRankineCycleSuperStructure::from_json(file)?))
//     }

//     #[staticmethod]
//     #[pyo3(text_signature = "(json)")]
//     fn from_json_str(json: &str) -> Result<Self, ParameterError> {
//         Ok(Self(serde_json::from_str(json)?))
//     }

//     fn to_json_str(&self) -> Result<String, ParameterError> {
//         Ok(serde_json::to_string(&self.0)?)
//     }
// }
