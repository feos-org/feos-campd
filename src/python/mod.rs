use crate::{
    CoMTCAMD, CoMTCAMDBinary, OptimizationProblem, OptimizationResult, OuterApproximationAlgorithm,
};
use crate::{MolecularRepresentation, PropertyModel};
use feos::core::EquationOfState;
use feos::ideal_gas::IdealGasModel;
use feos::pcsaft::PcSaft;
use feos_core::EosError;
use pyo3::prelude::*;
use pyo3::wrap_pymodule;
use quantity::python::quantity as quantity_module;
use std::collections::HashMap;

mod comt_camd;
mod dippr;
mod eos;
mod joback;
mod pcsaft;
mod process;

use comt_camd::{PyCoMTCAMD, PyCoMTCAMDBinary};
use dippr::dippr as dippr_module;
use eos::{eos as eos_module, PyEquationOfState};
use joback::joback as joback_module;
use pcsaft::pcsaft as pcsaft_module;
use process::PyProcessModel;

#[pyclass(name = "OuterApproximationAlgorithm")]
#[derive(Clone, Copy)]
pub struct PyOuterApproximationAlgorithm(OuterApproximationAlgorithm);

#[pymethods]
#[allow(non_snake_case)]
impl PyOuterApproximationAlgorithm {
    #[staticmethod]
    fn DuranGrossmann(update_lower_bound: bool) -> Self {
        Self(OuterApproximationAlgorithm::DuranGrossmann(
            update_lower_bound,
        ))
    }

    #[classattr]
    fn FletcherLeyffer() -> Self {
        Self(OuterApproximationAlgorithm::FletcherLeyffer)
    }
}

#[allow(clippy::large_enum_variant)]
enum OptimizationProblems {
    Pure(
        OptimizationProblem<
            EquationOfState<IdealGasModel, PcSaft>,
            CoMTCAMD,
            CoMTCAMD,
            PyProcessModel,
        >,
    ),
    Binary(
        OptimizationProblem<
            EquationOfState<IdealGasModel, PcSaft>,
            CoMTCAMDBinary,
            CoMTCAMDBinary,
            PyProcessModel,
        >,
    ),
}

#[pyclass(name = "OptimizationProblem")]
pub struct PyOptimizationProblem(OptimizationProblems);

#[pymethods]
impl PyOptimizationProblem {
    #[staticmethod]
    fn pure(comt_camd: PyCoMTCAMD, process: PyProcessModel) -> Self {
        Self(OptimizationProblems::Pure(OptimizationProblem::new(
            comt_camd.0.clone(),
            comt_camd.0,
            process,
        )))
    }

    #[staticmethod]
    fn binary(comt_camd: PyCoMTCAMDBinary, process: PyProcessModel) -> Self {
        Self(OptimizationProblems::Binary(OptimizationProblem::new(
            comt_camd.0.clone(),
            comt_camd.0,
            process,
        )))
    }

    fn solve_fixed(
        &mut self,
        y_fixed: Vec<f64>,
        options: Option<&str>,
    ) -> PyResult<PyOptimizationResult> {
        Ok(PyOptimizationResult(match &mut self.0 {
            OptimizationProblems::Pure(p) => p
                .solve_fixed(&y_fixed, options)
                .map_err(|_| EosError::NotConverged("Process optimization".into()))?,
            OptimizationProblems::Binary(b) => b
                .solve_fixed(&y_fixed, options)
                .map_err(|_| EosError::NotConverged("Process optimization".into()))?,
        }))
    }

    fn solve_outer_approximation(
        &mut self,
        y_fixed: Vec<f64>,
        algorithm: PyOuterApproximationAlgorithm,
        options_nlp: Option<&str>,
        options_mip: Option<&str>,
    ) -> PyResult<PyOptimizationResult> {
        Ok(PyOptimizationResult(match &mut self.0 {
            OptimizationProblems::Pure(p) => {
                p.solve_outer_approximation(y_fixed, algorithm.0, options_nlp, options_mip)
            }
            OptimizationProblems::Binary(b) => {
                b.solve_outer_approximation(y_fixed, algorithm.0, options_nlp, options_mip)
            }
        }?))
    }

    fn outer_approximation_ranking(
        &mut self,
        y_fixed: Vec<f64>,
        algorithm: PyOuterApproximationAlgorithm,
        runs: usize,
        options_nlp: Option<&str>,
        options_mip: Option<&str>,
    ) {
        match &mut self.0 {
            OptimizationProblems::Pure(p) => {
                p.outer_approximation_ranking(&y_fixed, algorithm.0, runs, options_nlp, options_mip)
            }
            OptimizationProblems::Binary(b) => {
                b.outer_approximation_ranking(&y_fixed, algorithm.0, runs, options_nlp, options_mip)
            }
        }
    }

    #[getter]
    fn get_solutions(&self) -> Vec<PyOptimizationResult> {
        let mut solutions: Vec<_> = match &self.0 {
            OptimizationProblems::Pure(p) => p.solutions.iter().collect(),
            OptimizationProblems::Binary(b) => b.solutions.iter().collect(),
        };
        solutions.sort_by(|s1, s2| s1.target.total_cmp(&s2.target));
        solutions
            .into_iter()
            .map(|result| PyOptimizationResult(result.clone()))
            .collect()
    }

    fn build_eos(&self, result: PyOptimizationResult) -> PyResult<PyEquationOfState> {
        let y: Vec<_> = result.0.y.into_iter().map(|y| y as f64).collect();
        Ok(PyEquationOfState(match &self.0 {
            OptimizationProblems::Pure(p) => {
                let cr = p.molecules.build(&y, &result.0.p);
                p.property_model.build_eos(cr)?
            }
            OptimizationProblems::Binary(b) => {
                let cr = b.molecules.build(&y, &result.0.p);
                b.property_model.build_eos(cr)?
            }
        }))
    }
}

#[pyclass(name = "PyOptimizationResult")]
#[derive(Clone)]
pub struct PyOptimizationResult(OptimizationResult);

#[pymethods]
impl PyOptimizationResult {
    #[getter]
    fn get_target(&self) -> f64 {
        self.0.target
    }

    #[getter]
    fn get_x(&self) -> Vec<f64> {
        self.0.x.clone()
    }

    #[getter]
    fn get_y(&self) -> Vec<usize> {
        self.0.y.clone()
    }

    #[getter]
    fn get_p(&self) -> Vec<f64> {
        self.0.p.clone()
    }

    #[getter]
    fn get_smiles(&self) -> Vec<String> {
        self.0.smiles.clone()
    }

    #[getter]
    fn get_groups(&self) -> Vec<HashMap<String, usize>> {
        self.0
            .smiles
            .iter()
            .map(|smiles| {
                let mut count = HashMap::new();
                smiles
                    .split(',')
                    .for_each(|g| *count.entry(g.to_string()).or_insert(0) += 1);
                count
            })
            .collect()
    }
}

#[pymodule]
pub fn feos_campd(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_wrapped(wrap_pymodule!(quantity_module))?;
    m.add_wrapped(wrap_pymodule!(eos_module))?;
    m.add_wrapped(wrap_pymodule!(joback_module))?;
    m.add_wrapped(wrap_pymodule!(dippr_module))?;
    m.add_wrapped(wrap_pymodule!(pcsaft_module))?;

    m.add_class::<PyProcessModel>()?;
    m.add_class::<PyCoMTCAMD>()?;
    m.add_class::<PyCoMTCAMDBinary>()?;
    m.add_class::<PyOptimizationProblem>()?;
    m.add_class::<PyOuterApproximationAlgorithm>()?;

    set_path(py, m, "feos_campd.si", "quantity")?;
    set_path(py, m, "feos_campd.eos", "eos")?;
    set_path(py, m, "feos_campd.joback", "joback")?;
    set_path(py, m, "feos_campd.dippr", "dippr")?;
    set_path(py, m, "feos_campd.pcsaft", "pcsaft")?;
    Ok(())
}

fn set_path(py: Python<'_>, m: &PyModule, path: &str, module: &str) -> PyResult<()> {
    py.run(
        &format!(
            "\
import sys
sys.modules['{path}'] = {module}
    "
        ),
        None,
        Some(m.dict()),
    )
}
