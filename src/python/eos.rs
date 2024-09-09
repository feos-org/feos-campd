use feos::core::si::*;
use feos::core::*;
use feos::gc_pcsaft::GcPcSaft;
use feos::ideal_gas::Joback;
use feos::pcsaft::PcSaft;
use numpy::convert::ToPyArray;
use numpy::{PyArray1, PyArray2};
use numpy::prelude::*;
use pyo3::exceptions::{PyIndexError, PyValueError};
use pyo3::prelude::*;
use quantity::python::{PySIArray1, PySIArray2, PySINumber};
use std::collections::HashMap;
use std::convert::{TryFrom, TryInto};
use std::sync::Arc;
use typenum::P3;

pub mod pcsaft {
    use super::*;

    /// Collection of equations of state.
    #[pyclass(name = "PcSaft")]
    #[derive(Clone)]
    pub struct PyPcSaft(pub Arc<EquationOfState<Joback, PcSaft>>);

    impl_equation_of_state!(PyPcSaft);
    impl_virial_coefficients!(PyPcSaft);
    impl_state!(EquationOfState<Joback, PcSaft>, PyPcSaft);
    impl_state_entropy_scaling!(EquationOfState<Joback, PcSaft>, PyPcSaft);
    impl_phase_equilibrium!(EquationOfState<Joback, PcSaft>, PyPcSaft);

    #[pymodule]
    pub fn pcsaft(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_class::<Contributions>()?;
        m.add_class::<Verbosity>()?;

        m.add_class::<PyState>()?;
        m.add_class::<PyStateVec>()?;
        m.add_class::<PyPhaseDiagram>()?;
        m.add_class::<PyPhaseEquilibrium>()?;

        Ok(())
    }
}

pub mod gc_pcsaft {
    use super::*;

    /// Collection of equations of state.
    #[pyclass(name = "GcPcSaft")]
    #[derive(Clone)]
    pub struct PyGcPcSaft(pub Arc<EquationOfState<Joback, GcPcSaft>>);

    impl_equation_of_state!(PyGcPcSaft);
    impl_virial_coefficients!(PyGcPcSaft);
    impl_state!(EquationOfState<Joback, GcPcSaft>, PyGcPcSaft);
    impl_phase_equilibrium!(EquationOfState<Joback, GcPcSaft>, PyGcPcSaft);

    #[pymodule]
    pub fn gc_pcsaft(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_class::<Contributions>()?;
        m.add_class::<Verbosity>()?;

        m.add_class::<PyState>()?;
        m.add_class::<PyStateVec>()?;
        m.add_class::<PyPhaseDiagram>()?;
        m.add_class::<PyPhaseEquilibrium>()?;

        Ok(())
    }
}
