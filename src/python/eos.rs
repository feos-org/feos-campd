use super::dippr::PyDippr;
use super::joback::PyJoback;
use super::pcsaft::PyPcSaftParameters;
use feos::core::si::*;
use feos::core::*;
use feos::ideal_gas::IdealGasModel;
use feos::pcsaft::{DQVariants, PcSaft, PcSaftOptions};
use numpy::convert::ToPyArray;
use numpy::{PyArray1, PyArray2};
use pyo3::exceptions::{PyIndexError, PyValueError};
use pyo3::prelude::*;
#[cfg(feature = "estimator")]
use pyo3::wrap_pymodule;
use quantity::python::{PySIArray1, PySIArray2, PySINumber};
use std::collections::HashMap;
use std::convert::{TryFrom, TryInto};
use std::sync::Arc;
use typenum::P3;

/// Collection of equations of state.
#[pyclass(name = "EquationOfState")]
#[derive(Clone)]
pub struct PyEquationOfState(pub Arc<EquationOfState<IdealGasModel, PcSaft>>);

#[pymethods]
impl PyEquationOfState {
    /// PC-SAFT equation of state.
    ///
    /// Parameters
    /// ----------
    /// parameters : PcSaftParameters
    ///     The parameters of the PC-SAFT equation of state to use.
    /// max_eta : float, optional
    ///     Maximum packing fraction. Defaults to 0.5.
    /// max_iter_cross_assoc : unsigned integer, optional
    ///     Maximum number of iterations for cross association. Defaults to 50.
    /// tol_cross_assoc : float
    ///     Tolerance for convergence of cross association. Defaults to 1e-10.
    ///
    /// Returns
    /// -------
    /// EquationOfState
    ///     The PC-SAFT equation of state that can be used to compute thermodynamic
    ///     states.
    #[staticmethod]
    #[pyo3(
        signature = (parameters, max_eta=0.5, max_iter_cross_assoc=50, tol_cross_assoc=1e-10),
        text_signature = "(parameters, max_eta=0.5, max_iter_cross_assoc=50, tol_cross_assoc=1e-10)"
    )]
    pub fn pcsaft(
        parameters: PyPcSaftParameters,
        max_eta: f64,
        max_iter_cross_assoc: usize,
        tol_cross_assoc: f64,
    ) -> Self {
        let options = PcSaftOptions {
            max_eta,
            max_iter_cross_assoc,
            tol_cross_assoc,
            dq_variant: DQVariants::DQ35,
        };
        let residual = Arc::new(PcSaft::with_options(parameters.0, options));
        let ideal_gas = Arc::new(IdealGasModel::NoModel(residual.components()));
        Self(Arc::new(EquationOfState::new(ideal_gas, residual)))
    }

    /// Ideal gas model of Joback and Reid.
    ///
    /// Parameters
    /// ----------
    /// joback : Joback
    ///     The parametrized Joback model.
    ///
    /// Returns
    /// -------
    /// EquationOfState
    fn joback(&self, joback: PyJoback) -> Self {
        self.add_ideal_gas(IdealGasModel::Joback(joback.0))
    }

    /// Ideal gas model based on DIPPR equations for the ideal
    /// gas heat capacity.
    ///
    /// Parameters
    /// ----------
    /// dippr : Dippr
    ///     The parametrized Dippr model.
    ///
    /// Returns
    /// -------
    /// EquationOfState
    fn dippr(&self, dippr: PyDippr) -> Self {
        self.add_ideal_gas(IdealGasModel::Dippr(dippr.0))
    }
}

impl PyEquationOfState {
    fn add_ideal_gas(&self, ideal_gas: IdealGasModel) -> Self {
        Self(Arc::new(EquationOfState::new(
            Arc::new(ideal_gas),
            self.0.residual.clone(),
        )))
    }
}

impl_equation_of_state!(PyEquationOfState);
impl_virial_coefficients!(PyEquationOfState);
impl_state!(EquationOfState<IdealGasModel, PcSaft>, PyEquationOfState);
impl_state_entropy_scaling!(EquationOfState<IdealGasModel, PcSaft>, PyEquationOfState);
impl_phase_equilibrium!(EquationOfState<IdealGasModel, PcSaft>, PyEquationOfState);

#[pymodule]
pub fn eos(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_class::<Contributions>()?;
    m.add_class::<Verbosity>()?;

    m.add_class::<PyEquationOfState>()?;
    m.add_class::<PyState>()?;
    m.add_class::<PyStateVec>()?;
    m.add_class::<PyPhaseDiagram>()?;
    m.add_class::<PyPhaseEquilibrium>()?;

    Ok(())
}
