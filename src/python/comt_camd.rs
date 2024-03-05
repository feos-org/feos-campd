use crate::{CoMTCAMD, CoMTCAMDBinary};
use feos_core::python::parameter::PyBinarySegmentRecord;
use pyo3::prelude::*;
use std::collections::HashMap;

#[pyclass(name = "CoMTCAMD")]
#[derive(Clone)]
pub struct PyCoMTCAMD(pub CoMTCAMD);

#[pymethods]
impl PyCoMTCAMD {
    #[staticmethod]
    fn from_json(camd_records: String) -> PyResult<Self> {
        Ok(Self(CoMTCAMD::from_json(camd_records)?))
    }

    fn get_initial_values(&self, structure: &str, groups: HashMap<&str, usize>) -> Vec<f64> {
        self.0.get_initial_values(structure, &groups)
    }
}

#[pyclass(name = "CoMTCAMDBinary")]
#[derive(Clone)]
pub struct PyCoMTCAMDBinary(pub CoMTCAMDBinary);

#[pymethods]
impl PyCoMTCAMDBinary {
    #[new]
    fn new(
        molecule1: PyCoMTCAMD,
        molecule2: PyCoMTCAMD,
        symmetry_constraints: bool,
        k_ij: Option<Vec<PyBinarySegmentRecord>>,
    ) -> Self {
        Self(CoMTCAMDBinary::new(
            [molecule1.0, molecule2.0],
            symmetry_constraints,
            k_ij.map(|k_ij| k_ij.into_iter().map(|k_ij| k_ij.0).collect()),
        ))
    }
}
