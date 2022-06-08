use feos_core::parameter::ParameterError;
use feos_core::{EquationOfState, MolarWeight};
use quantity::si::SIUnit;

pub trait PropertyModel<C> {
    type Eos: EquationOfState + MolarWeight<SIUnit>;
    fn build_eos(&self, chemical_record: C) -> Result<Self::Eos, ParameterError>;
}
