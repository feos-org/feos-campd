use crate::process::{Process, Utility, UtilitySpecification};
use crate::property::EquationsOfState;
use pyo3::prelude::*;
use quantity::python::*;

#[pyclass(name = "UtilitySpecification")]
#[derive(Clone)]
pub struct PyUtilitySpecification(pub UtilitySpecification);

#[pymethods]
#[allow(non_snake_case)]
impl PyUtilitySpecification {
    #[classattr]
    fn ConstantTemperature() -> Self {
        Self(UtilitySpecification::ConstantTemperature)
    }

    #[staticmethod]
    fn HeatCapacityRate(rate: PySINumber) -> Self {
        Self(UtilitySpecification::HeatCapacityRate(rate.into()))
    }

    #[staticmethod]
    fn OutletTemperature(temperature: PySINumber) -> Self {
        Self(UtilitySpecification::OutletTemperature(temperature.into()))
    }
}

#[pyclass(name = "Utility")]
#[derive(Clone)]
pub struct PyUtility(pub Utility);

#[pymethods]
#[allow(non_snake_case)]
impl PyUtility {
    #[new]
    fn new(
        temperature: PySINumber,
        specification: PyUtilitySpecification,
        min_approach_temperature: PySINumber,
    ) -> Self {
        Self(Utility::new(
            temperature.into(),
            specification.0,
            min_approach_temperature.into(),
        ))
    }
}

#[pyclass(name = "Process", unsendable)]
pub struct PyProcess(pub Process<EquationsOfState>);

#[pymethods]
impl PyProcess {
    #[getter]
    fn net_power(&self) -> Option<PySINumber> {
        self.0.net_power().map(|p| p.into())
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
