use super::{PyGcPcSaft, PyPcSaft};
use crate::process::ProcessModel;
use crate::{ProcessVariables, Variable};
use feos::core::{EosResult, EquationOfState};
use feos::gc_pcsaft::GcPcSaft;
use feos::ideal_gas::Joback;
use feos::pcsaft::PcSaft;
use feos_core::EosError;
use pyo3::prelude::*;
use std::sync::Arc;

#[pyclass(name = "ProcessModel")]
#[derive(Clone)]
pub struct PyProcessModel(Py<PyAny>);

#[pymethods]
impl PyProcessModel {
    #[new]
    pub fn new(obj: Py<PyAny>) -> PyResult<Self> {
        Python::with_gil(|py| {
            let attr = obj.as_ref(py).hasattr("variables")?;
            if !attr {
                panic!("Python class must have a method 'variables' with signature:\n\tdef variables(self) -> List[[int, int, int]]")
            }
            let attr = obj.as_ref(py).hasattr("equality_constraints")?;
            if !attr {
                panic!("Python class must have a method 'equality_constraints' with signature:\n\tdef equality_constraints(self) -> int")
            }
            let attr = obj.as_ref(py).hasattr("inequality_constraints")?;
            if !attr {
                panic!("Python class must have a method 'inequality_constraints' with signature:\n\tdef inequality_constraints(self) -> int")
            }
            let attr = obj.as_ref(py).hasattr("solve")?;
            if !attr {
                panic!("Python class must have a method 'solve' with signature:\n\tdef solve(self, eos: EquationOfState, x: List[float]) -> (float, List[float], List[float])")
            }
            Ok(Self(obj))
        })
    }
}

macro_rules! impl_process_model {
    ($py_eos:ident, $eos:ident) => {
        impl ProcessModel<EquationOfState<Joback, $eos>> for PyProcessModel {
            fn variables(&self) -> ProcessVariables {
                Python::with_gil(|py| {
                    let py_result = self.0.as_ref(py).call_method0("variables").unwrap();
                    if py_result.get_type().name().unwrap() != "list" {
                        panic!(
                            "Expected a list for the variables() method signature, got {}",
                            py_result.get_type().name().unwrap()
                        );
                    }
                    let variables = py_result.extract::<Vec<[f64; 3]>>().unwrap();
                    variables
                        .into_iter()
                        .map(|[l, u, i]| Variable::continuous(l, u, i))
                        .collect()
                })
            }

            fn equality_constraints(&self) -> usize {
                Python::with_gil(|py| {
                    let py_result = self
                        .0
                        .as_ref(py)
                        .call_method0("equality_constraints")
                        .unwrap();
                    if py_result.get_type().name().unwrap() != "int" {
                        panic!(
                            "Expected an integer for the equality_constraints() method signature, got {}",
                            py_result.get_type().name().unwrap()
                        );
                    }
                    py_result.extract().unwrap()
                })
            }

            fn inequality_constraints(&self) -> usize {
                Python::with_gil(|py| {
                    let py_result = self
                        .0
                        .as_ref(py)
                        .call_method0("inequality_constraints")
                        .unwrap();
                    if py_result.get_type().name().unwrap() != "int" {
                        panic!(
                            "Expected an integer for the inequality_constraints() method signature, got {}",
                            py_result.get_type().name().unwrap()
                        );
                    }
                    py_result.extract().unwrap()
                })
            }

            fn solve(
                &self,
                eos: &Arc<EquationOfState<Joback, $eos>>,
                x: &[f64],
            ) -> EosResult<(f64, Vec<f64>, Vec<f64>)> {
                Python::with_gil(|py| {
                    let py_eos = $py_eos(eos.clone());
                    self.0
                        .as_ref(py)
                        .call_method1("solve", (py_eos, x.to_vec()))
                        .and_then(|py_result| py_result.extract::<(f64, Vec<f64>, Vec<f64>)>())
                        .map_err(|_| EosError::NotConverged("Process model".into()))
                })
            }
        }
    };
}

impl_process_model!(PyPcSaft, PcSaft);
impl_process_model!(PyGcPcSaft, GcPcSaft);
