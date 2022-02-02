use super::property::PyEquationsOfState;
use super::{OrganicRankineCycle, Process};
use crate::property::EquationsOfState;
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
        heat_source_temperature: PySINumber,
        heat_source_mass_flow_rate: PySINumber,
        heat_source_heat_capacity: PySINumber,
        isentropic_turbine_efficiency: f64,
        isentropic_pump_efficiency: f64,
        min_approach_temperature: PySINumber,
        min_abs_pressure: PySINumber,
        min_rel_pressure: f64,
        max_abs_pressure: PySINumber,
        max_red_pressure: f64,
        min_cooling_temperature: PySINumber,
    ) -> Self {
        Self(OrganicRankineCycle::new(
            heat_source_temperature.into(),
            heat_source_mass_flow_rate.into(),
            heat_source_heat_capacity.into(),
            isentropic_turbine_efficiency,
            isentropic_pump_efficiency,
            min_approach_temperature.into(),
            min_abs_pressure.into(),
            min_rel_pressure,
            max_abs_pressure.into(),
            max_red_pressure,
            min_cooling_temperature.into(),
        ))
    }

    #[getter]
    fn get_variables(&self) -> Vec<[Option<f64>; 2]> {
        self.0.variables()
    }

    #[getter]
    fn get_constraints(&self) -> Vec<[Option<f64>; 3]> {
        self.0.constraints()
    }

    fn solve(
        &self,
        eos: PyEquationsOfState,
        x: [f64; 4],
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

#[pyclass(name = "Process", unsendable)]
struct PyProcess(Process<EquationsOfState>);

#[pymethods]
impl PyProcess {
    #[getter]
    fn net_power(&self) -> PySINumber {
        self.0.net_power().into()
    }

    fn plot(&self) -> PyResult<Vec<(PySIArray1, PySIArray1, Option<PySIArray1>)>> {
        let mut res = self.0.plot()?;
        Ok(res
            .drain(..)
            .map(|plot| {
                (
                    plot.temperature.into(),
                    plot.entropy.into(),
                    plot.utility_temperature.map(PySIArray1::from),
                )
            })
            .collect())
    }
}
