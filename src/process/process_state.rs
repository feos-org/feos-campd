use feos_core::{
    Contributions, EosResult, EquationOfState, MolarWeight, PhaseEquilibrium, State, StateBuilder,
};
use quantity::si::{SINumber, SIUnit};
use std::rc::Rc;

pub enum ProcessState<E> {
    SinglePhase(Box<State<SIUnit, E>>),
    TwoPhase(Rc<PhaseEquilibrium<SIUnit, E, 2>>, f64),
}

impl<E> Clone for ProcessState<E> {
    fn clone(&self) -> Self {
        match self {
            Self::SinglePhase(s) => Self::SinglePhase(s.clone()),
            Self::TwoPhase(vle, x) => Self::TwoPhase(vle.clone(), *x),
        }
    }
}

impl<E: EquationOfState + MolarWeight<SIUnit>> ProcessState<E> {
    pub fn new_pure_enthalpy_pressure(
        eos: &Rc<E>,
        enthalpy: SINumber,
        pressure: SINumber,
        vle: &Rc<PhaseEquilibrium<SIUnit, E, 2>>,
        initial_temperature: SINumber,
    ) -> EosResult<Self> {
        let hv = vle.vapor().molar_enthalpy(Contributions::Total);
        let hl = vle.liquid().molar_enthalpy(Contributions::Total);

        Ok(if enthalpy > hl && enthalpy < hv {
            Self::TwoPhase(vle.clone(), (enthalpy - hl).to_reduced(hv - hl)?)
        } else {
            let state = StateBuilder::new(eos)
                .pressure(pressure)
                .enthalpy(enthalpy)
                .initial_temperature(initial_temperature)
                .build()?;
            Self::SinglePhase(Box::new(state))
        })
    }

    pub fn new_pure_entropy_pressure(
        eos: &Rc<E>,
        entropy: SINumber,
        pressure: SINumber,
        vle: &Rc<PhaseEquilibrium<SIUnit, E, 2>>,
        initial_temperature: SINumber,
    ) -> EosResult<Self> {
        let sv = vle.vapor().molar_entropy(Contributions::Total);
        let sl = vle.liquid().molar_entropy(Contributions::Total);
        Ok(if entropy > sl && entropy < sv {
            Self::TwoPhase(vle.clone(), (entropy - sl).to_reduced(sv - sl)?)
        } else {
            let state = StateBuilder::new(eos)
                .pressure(pressure)
                .entropy(entropy)
                .initial_temperature(initial_temperature)
                .build()?;
            Self::SinglePhase(Box::new(state))
        })
    }

    pub fn is_two_phase(&self) -> bool {
        match self {
            Self::SinglePhase(_) => false,
            Self::TwoPhase(_, _) => true,
        }
    }

    pub fn eos(&self) -> &Rc<E> {
        match self {
            Self::SinglePhase(state) => &state.eos,
            Self::TwoPhase(vle, _) => &vle.vapor().eos,
        }
    }

    pub fn temperature(&self) -> SINumber {
        match self {
            Self::SinglePhase(state) => state.temperature,
            Self::TwoPhase(vle, _) => vle.vapor().temperature,
        }
    }

    pub fn density(&self) -> SINumber {
        match self {
            Self::SinglePhase(state) => state.density,
            Self::TwoPhase(vle, _) => vle.vapor().density,
        }
    }

    pub fn pressure(&self) -> SINumber {
        match self {
            Self::SinglePhase(state) => state.pressure(Contributions::Total),
            Self::TwoPhase(vle, _) => vle.vapor().pressure(Contributions::Total),
        }
    }

    pub fn specific_enthalpy(&self) -> SINumber {
        match self {
            Self::SinglePhase(state) => state.specific_enthalpy(Contributions::Total),
            Self::TwoPhase(vle, vapor_fraction) => {
                vle.vapor().specific_enthalpy(Contributions::Total) * *vapor_fraction
                    + vle.liquid().specific_enthalpy(Contributions::Total) * (1.0 - vapor_fraction)
            }
        }
    }

    pub fn specific_entropy(&self) -> SINumber {
        match self {
            Self::SinglePhase(state) => state.specific_entropy(Contributions::Total),
            Self::TwoPhase(vle, vapor_fraction) => {
                vle.vapor().specific_entropy(Contributions::Total) * *vapor_fraction
                    + vle.liquid().specific_entropy(Contributions::Total) * (1.0 - vapor_fraction)
            }
        }
    }

    pub fn molar_enthalpy(&self) -> SINumber {
        match self {
            Self::SinglePhase(state) => state.molar_enthalpy(Contributions::Total),
            Self::TwoPhase(vle, vapor_fraction) => {
                vle.vapor().molar_enthalpy(Contributions::Total) * *vapor_fraction
                    + vle.liquid().molar_enthalpy(Contributions::Total) * (1.0 - vapor_fraction)
            }
        }
    }

    pub fn molar_entropy(&self) -> SINumber {
        match self {
            Self::SinglePhase(state) => state.molar_entropy(Contributions::Total),
            Self::TwoPhase(vle, vapor_fraction) => {
                vle.vapor().molar_entropy(Contributions::Total) * *vapor_fraction
                    + vle.liquid().molar_entropy(Contributions::Total) * (1.0 - vapor_fraction)
            }
        }
    }
}
