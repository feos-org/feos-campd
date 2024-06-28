use crate::GcPcSaftPropertyModel;
use crate::{
    CoMTCAMD, OptimizationProblem, OptimizationResult, OuterApproximationAlgorithm,
    PcSaftPropertyModel,
};
use feos_core::EosError;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::wrap_pymodule;
use quantity::python::quantity as quantity_module;
use std::collections::HashMap;

mod eos;
mod molecule;
mod process;
mod property;

use eos::gc_pcsaft::{gc_pcsaft as gc_pcsaft_module, PyGcPcSaft};
use eos::pcsaft::{pcsaft as pcsaft_module, PyPcSaft};
use molecule::{PyCoMTCAMD, PySuperMolecule, SuperMolecules};
use process::PyProcessModel;
use property::{PyGcPcSaftPropertyModel, PyPcSaftPropertyModel};

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

#[pyclass(name = "OptimizationProblem")]
pub struct PyOptimizationProblem;

#[pymethods]
impl PyOptimizationProblem {
    #[staticmethod]
    fn pure(
        py: Python,
        molecule: &Bound<'_, PyAny>,
        property: &Bound<'_, PyAny>,
        process: PyProcessModel,
    ) -> PyResult<PyObject> {
        if let (Ok(comt_camd), Ok(pcsaft)) = (
            molecule.extract::<PyCoMTCAMD>(),
            property.extract::<PyPcSaftPropertyModel>(),
        ) {
            Ok(
                PyOptimizationProblem1(OptimizationProblem::new([comt_camd.0], pcsaft.0, process))
                    .into_py(py),
            )
        } else if let (Ok(supermolecule), Ok(pcsaft)) = (
            molecule.extract::<PySuperMolecule>(),
            property.extract::<PyPcSaftPropertyModel>(),
        ) {
            Ok(PyOptimizationProblem2(OptimizationProblem::new(
                [supermolecule.0],
                pcsaft.0,
                process,
            ))
            .into_py(py))
        } else if let (Ok(supermolecule), Ok(gc_pcsaft)) = (
            molecule.extract::<PySuperMolecule>(),
            property.extract::<PyGcPcSaftPropertyModel>(),
        ) {
            Ok(PyOptimizationProblem3(OptimizationProblem::new(
                [supermolecule.0],
                gc_pcsaft.0,
                process,
            ))
            .into_py(py))
        } else {
            Err(PyValueError::new_err("Ivalid arguments!"))
        }
    }

    #[staticmethod]
    fn binary(
        py: Python,
        molecule: &Bound<'_, PyAny>,
        property: &Bound<'_, PyAny>,
        process: PyProcessModel,
    ) -> PyResult<PyObject> {
        if let (Ok([m1, m2]), Ok(pcsaft)) = (
            molecule.extract::<[PyCoMTCAMD; 2]>(),
            property.extract::<PyPcSaftPropertyModel>(),
        ) {
            Ok(
                PyOptimizationProblem4(OptimizationProblem::new([m1.0, m2.0], pcsaft.0, process))
                    .into_py(py),
            )
        } else if let (Ok([m1, m2]), Ok(pcsaft)) = (
            molecule.extract::<[PySuperMolecule; 2]>(),
            property.extract::<PyPcSaftPropertyModel>(),
        ) {
            Ok(
                PyOptimizationProblem5(OptimizationProblem::new([m1.0, m2.0], pcsaft.0, process))
                    .into_py(py),
            )
        } else if let (Ok([m1, m2]), Ok(gc_pcsaft)) = (
            molecule.extract::<[PySuperMolecule; 2]>(),
            property.extract::<PyGcPcSaftPropertyModel>(),
        ) {
            Ok(
                PyOptimizationProblem6(OptimizationProblem::new(
                    [m1.0, m2.0],
                    gc_pcsaft.0,
                    process,
                ))
                .into_py(py),
            )
        } else {
            Err(PyValueError::new_err("Ivalid arguments!"))
        }
    }
}

macro_rules! impl_optimization_problem {
    ($py_optimization_problem:ident, $molecule:ident, $property_model:ident, $py_eos:ident, 1) => {
        #[pyclass(name = "OptimizationProblem")]
        pub struct $py_optimization_problem(
            OptimizationProblem<$molecule, $property_model, PyProcessModel, 1>,
        );

        #[pymethods]
        impl $py_optimization_problem {
            fn solve_fixed(
                &mut self,
                y_fixed: Vec<f64>,
                options: Option<&str>,
            ) -> PyResult<PyOptimizationResult> {
                Ok(PyOptimizationResult(
                    self.0
                        .solve_fixed(&[y_fixed], options)
                        .map_err(|_| EosError::NotConverged("Process optimization".into()))?,
                ))
            }

            fn solve_outer_approximation(
                &mut self,
                y_fixed: Vec<f64>,
                algorithm: PyOuterApproximationAlgorithm,
                options_nlp: Option<&str>,
                options_mip: Option<&str>,
            ) -> PyResult<PyOptimizationResult> {
                Ok(PyOptimizationResult(self.0.solve_outer_approximation(
                    [y_fixed],
                    algorithm.0,
                    options_nlp,
                    options_mip,
                )?))
            }

            fn outer_approximation_ranking(
                &mut self,
                y_fixed: Vec<f64>,
                algorithm: PyOuterApproximationAlgorithm,
                runs: usize,
                options_nlp: Option<&str>,
                options_mip: Option<&str>,
            ) {
                self.0.outer_approximation_ranking(
                    [y_fixed],
                    algorithm.0,
                    runs,
                    options_nlp,
                    options_mip,
                )
            }

            #[getter]
            fn get_solutions(&self) -> Vec<PyOptimizationResult> {
                let mut solutions: Vec<_> = self.0.solutions.iter().collect();
                solutions.sort_by(|s1, s2| s1.target.total_cmp(&s2.target));
                solutions
                    .into_iter()
                    .map(|result| PyOptimizationResult(result.clone()))
                    .collect()
            }

            fn build_eos(&self, result: PyOptimizationResult, py: Python) -> PyObject {
                $py_eos(self.0.build_eos(&result.0)).into_py(py)
            }
        }
    };
    ($py_optimization_problem:ident, $molecule:ident, $property_model:ident, $py_eos:ident, 2) => {
        #[pyclass(name = "OptimizationProblem")]
        pub struct $py_optimization_problem(
            OptimizationProblem<$molecule, $property_model, PyProcessModel, 2>,
        );

        #[pymethods]
        impl $py_optimization_problem {
            fn solve_fixed(
                &mut self,
                y_fixed: [Vec<f64>; 2],
                options: Option<&str>,
            ) -> PyResult<PyOptimizationResult> {
                Ok(PyOptimizationResult(
                    self.0
                        .solve_fixed(&y_fixed, options)
                        .map_err(|_| EosError::NotConverged("Process optimization".into()))?,
                ))
            }

            fn solve_outer_approximation(
                &mut self,
                y_fixed: [Vec<f64>; 2],
                algorithm: PyOuterApproximationAlgorithm,
                options_nlp: Option<&str>,
                options_mip: Option<&str>,
            ) -> PyResult<PyOptimizationResult> {
                Ok(PyOptimizationResult(self.0.solve_outer_approximation(
                    y_fixed,
                    algorithm.0,
                    options_nlp,
                    options_mip,
                )?))
            }

            fn outer_approximation_ranking(
                &mut self,
                y_fixed: [Vec<f64>; 2],
                algorithm: PyOuterApproximationAlgorithm,
                runs: usize,
                options_nlp: Option<&str>,
                options_mip: Option<&str>,
            ) {
                self.0.outer_approximation_ranking(
                    y_fixed,
                    algorithm.0,
                    runs,
                    options_nlp,
                    options_mip,
                )
            }

            #[getter]
            fn get_solutions(&self) -> Vec<PyOptimizationResult> {
                let mut solutions: Vec<_> = self.0.solutions.iter().collect();
                solutions.sort_by(|s1, s2| s1.target.total_cmp(&s2.target));
                solutions
                    .into_iter()
                    .map(|result| PyOptimizationResult(result.clone()))
                    .collect()
            }

            fn build_eos(&self, result: PyOptimizationResult, py: Python) -> PyObject {
                $py_eos(self.0.build_eos(&result.0)).into_py(py)
            }
        }
    };
}

impl_optimization_problem!(
    PyOptimizationProblem1,
    CoMTCAMD,
    PcSaftPropertyModel,
    PyPcSaft,
    1
);
impl_optimization_problem!(
    PyOptimizationProblem2,
    SuperMolecules,
    PcSaftPropertyModel,
    PyPcSaft,
    1
);
impl_optimization_problem!(
    PyOptimizationProblem3,
    SuperMolecules,
    GcPcSaftPropertyModel,
    PyGcPcSaft,
    1
);
impl_optimization_problem!(
    PyOptimizationProblem4,
    CoMTCAMD,
    PcSaftPropertyModel,
    PyPcSaft,
    2
);
impl_optimization_problem!(
    PyOptimizationProblem5,
    SuperMolecules,
    PcSaftPropertyModel,
    PyPcSaft,
    2
);
impl_optimization_problem!(
    PyOptimizationProblem6,
    SuperMolecules,
    GcPcSaftPropertyModel,
    PyGcPcSaft,
    2
);

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
pub fn feos_campd(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_wrapped(wrap_pymodule!(quantity_module))?;
    m.add_wrapped(wrap_pymodule!(pcsaft_module))?;
    m.add_wrapped(wrap_pymodule!(gc_pcsaft_module))?;

    m.add_class::<PyProcessModel>()?;
    m.add_class::<PyCoMTCAMD>()?;
    m.add_class::<PySuperMolecule>()?;
    m.add_class::<PyPcSaftPropertyModel>()?;
    m.add_class::<PyGcPcSaftPropertyModel>()?;
    m.add_class::<PyOptimizationProblem>()?;
    m.add_class::<PyOuterApproximationAlgorithm>()?;

    set_path(m, "feos_campd.si", "quantity")?;
    set_path(m, "feos_campd.pcsaft", "pcsaft")?;
    set_path(m, "feos_campd.gc_pcsaft", "gc_pcsaft")?;
    Ok(())
}

fn set_path(m: &Bound<'_, PyModule>, path: &str, module: &str) -> PyResult<()> {
    m.py().run_bound(
        &format!(
            "\
import sys
sys.modules['{path}'] = {module}
    "
        ),
        None,
        Some(&m.dict()),
    )
}
