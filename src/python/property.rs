use super::{eos::pcsaft::PyPcSaft, PySuperMolecule};
use crate::{GcPcSaftPropertyModel, PcSaftPropertyModel, PropertyModel};
use feos_core::parameter::IdentifierOption;
use indexmap::IndexMap;
use pyo3::prelude::*;

#[pyclass(name = "PcSaftPropertyModel")]
#[derive(Clone)]
pub struct PyPcSaftPropertyModel(pub PcSaftPropertyModel);

#[pymethods]
impl PyPcSaftPropertyModel {
    #[staticmethod]
    #[pyo3(signature = (file_pcsaft, file_ideal_gas, file_binary=None, symmetry_constraints=true))]
    fn from_json_groups(
        file_pcsaft: String,
        file_ideal_gas: String,
        file_binary: Option<String>,
        symmetry_constraints: bool,
    ) -> PyResult<Self> {
        Ok(Self(PcSaftPropertyModel::from_json_groups(
            file_pcsaft,
            file_ideal_gas,
            file_binary,
            symmetry_constraints,
        )?))
    }

    #[staticmethod]
    #[pyo3(signature = (file_pcsaft, file_ideal_gas, file_binary=None, identifier_option=IdentifierOption::Name, symmetry_constraints=true))]
    fn from_json_molecules(
        file_pcsaft: String,
        file_ideal_gas: String,
        file_binary: Option<String>,
        identifier_option: IdentifierOption,
        symmetry_constraints: bool,
    ) -> PyResult<Self> {
        Ok(Self(PcSaftPropertyModel::from_json_molecules(
            file_pcsaft,
            file_ideal_gas,
            file_binary,
            identifier_option,
            symmetry_constraints,
        )?))
    }

    fn build_eos_pure(&self, feature_variables: [IndexMap<String, f64>; 1]) -> PyPcSaft {
        let p = self.0.evaluate_parameter_variables(&feature_variables);
        PyPcSaft(<PcSaftPropertyModel as PropertyModel<1>>::build_eos(
            &self.0, &p,
        ))
    }

    fn build_eos_binary(&self, feature_variables: [IndexMap<String, f64>; 2]) -> PyPcSaft {
        let p = self.0.evaluate_parameter_variables(&feature_variables);
        PyPcSaft(<PcSaftPropertyModel as PropertyModel<2>>::build_eos(
            &self.0, &p,
        ))
    }
}

#[pyclass(name = "GcPcSaftPropertyModel")]
#[derive(Clone)]
pub struct PyGcPcSaftPropertyModel(pub GcPcSaftPropertyModel);

#[pymethods]
impl PyGcPcSaftPropertyModel {
    #[staticmethod]
    pub fn from_json(
        py: Python,
        molecules: Vec<PySuperMolecule>,
        file_residual: String,
        file_ideal_gas: String,
        file_binary: Option<String>,
    ) -> PyResult<PyObject> {
        let molecules: Vec<_> = molecules.into_iter().map(|m| m.0).collect();
        Ok(PyGcPcSaftPropertyModel(GcPcSaftPropertyModel::from_json(
            &molecules,
            file_residual,
            file_ideal_gas,
            file_binary,
        )?)
        .into_py(py))
    }
}
