use crate::property::{EquationsOfState, PropertyModel};
use feos_core::parameter::ParameterError;
use feos_core::python::parameter::PyChemicalRecord;
use feos_core::python::{PyContributions, PyVerbosity};
use feos_core::*;
use numpy::convert::ToPyArray;
use numpy::{PyArray1, PyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use quantity::python::{PySIArray1, PySIArray2, PySINumber};
use quantity::si::*;
use std::collections::HashMap;
use std::rc::Rc;

#[pyclass(name = "EquationsOfState", unsendable)]
#[derive(Clone)]
pub struct PyEquationsOfState(pub Rc<EquationsOfState>);

impl_equation_of_state!(PyEquationsOfState);
impl_state!(EquationsOfState, PyEquationsOfState);
impl_state_molarweight!(PcSaft, PyPcSaft);
impl_vle_state!(EquationsOfState, PyEquationsOfState);

#[pyclass(name = "PropertyModel")]
#[derive(Clone)]
pub struct PyPropertyModel(pub PropertyModel);

#[pymethods]
impl PyPropertyModel {
    #[staticmethod]
    fn new_pcsaft(file: &str) -> Result<Self, ParameterError> {
        Ok(Self(PropertyModel::new_pcsaft(file)?))
    }

    #[staticmethod]
    fn new_gc_pcsaft(file: &str) -> Result<Self, ParameterError> {
        Ok(Self(PropertyModel::new_gc_pcsaft(file)?))
    }

    fn build_eos(
        &self,
        chemical_record: PyChemicalRecord,
    ) -> Result<PyEquationsOfState, ParameterError> {
        Ok(PyEquationsOfState(Rc::new(
            self.0.build_eos(chemical_record.0)?,
        )))
    }
}
