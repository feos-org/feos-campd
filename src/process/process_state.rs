use feos_core::{si::*, IdealGas};
use feos_core::{
    Contributions, DensityInitialization, EosError, EosResult, PhaseEquilibrium, Residual, State,
    StateBuilder,
};
use ndarray::arr1;
use std::sync::Arc;

/// Sub- or Supercritical isotherm.
pub enum Isobar<E> {
    SuperCritical(Pressure<f64>),
    SubCritical(Arc<PhaseEquilibrium<E, 2>>),
}

impl<E: Residual> Isobar<E> {
    pub fn new(eos: &Arc<E>, pressure: Pressure<f64>) -> Self {
        if let Ok(vle) = PhaseEquilibrium::pure(eos, pressure, None, Default::default()) {
            Self::SubCritical(Arc::new(vle))
        } else {
            Self::SuperCritical(pressure)
        }
    }

    pub fn pressure(&self) -> Pressure<f64> {
        match self {
            Self::SuperCritical(pressure) => *pressure,
            Self::SubCritical(vle) => vle.vapor().pressure(Contributions::Total),
        }
    }

    pub fn phase_equilibrium(&self) -> EosResult<&Arc<PhaseEquilibrium<E, 2>>> {
        match self {
            Self::SuperCritical(_) => Err(EosError::SuperCritical),
            Self::SubCritical(vle) => Ok(vle),
        }
    }

    pub fn saturated_liquid(&self) -> EosResult<&State<E>> {
        match self {
            Self::SuperCritical(_) => Err(EosError::SuperCritical),
            Self::SubCritical(vle) => Ok(vle.liquid()),
        }
    }

    pub fn saturated_vapor(&self) -> EosResult<&State<E>> {
        match self {
            Self::SuperCritical(_) => Err(EosError::SuperCritical),
            Self::SubCritical(vle) => Ok(vle.vapor()),
        }
    }
}

/// A [State] that can consist of either a single phase or two phases in equilibrium.
pub enum ProcessState<E> {
    SinglePhase(Box<State<E>>, Option<MassFlowRate<f64>>),
    TwoPhase(Arc<PhaseEquilibrium<E, 2>>, f64, Option<MassFlowRate<f64>>),
}

impl<E> Clone for ProcessState<E> {
    fn clone(&self) -> Self {
        match self {
            Self::SinglePhase(s, m) => Self::SinglePhase(s.clone(), *m),
            Self::TwoPhase(vle, x, m) => Self::TwoPhase(vle.clone(), *x, *m),
        }
    }
}

impl<E: Residual + IdealGas> ProcessState<E> {
    pub fn new_pure_temperature_pressure(
        eos: &Arc<E>,
        temperature: Temperature<f64>,
        pressure: Pressure<f64>,
        mass_flow_rate: Option<MassFlowRate<f64>>,
        density_initialization: DensityInitialization,
    ) -> EosResult<Self> {
        Ok(Self::SinglePhase(
            Box::new(State::new_npt(
                eos,
                temperature,
                pressure,
                &(arr1(&[1.0]) * MOL),
                density_initialization,
            )?),
            mass_flow_rate,
        ))
    }

    fn new_pure_enthalpy_pressure(
        eos: &Arc<E>,
        molar_enthalpy: MolarEnergy<f64>,
        isobar: &Isobar<E>,
        initial_temperature: Temperature<f64>,
        mass_flow_rate: Option<MassFlowRate<f64>>,
    ) -> EosResult<Self> {
        if let Isobar::SubCritical(vle) = isobar {
            let hv = vle.vapor().molar_enthalpy(Contributions::Total);
            let hl = vle.liquid().molar_enthalpy(Contributions::Total);

            if molar_enthalpy > hl && molar_enthalpy < hv {
                return Ok(Self::TwoPhase(
                    vle.clone(),
                    (molar_enthalpy - hl).convert_into(hv - hl),
                    mass_flow_rate,
                ));
            }
        }
        let state = StateBuilder::new(eos)
            .pressure(isobar.pressure())
            .molar_enthalpy(molar_enthalpy)
            .initial_temperature(initial_temperature)
            .build()?;
        Ok(Self::SinglePhase(Box::new(state), mass_flow_rate))
    }

    fn new_pure_entropy_pressure(
        eos: &Arc<E>,
        molar_entropy: MolarEntropy<f64>,
        isobar: &Isobar<E>,
        initial_temperature: Temperature<f64>,
        mass_flow_rate: Option<MassFlowRate<f64>>,
    ) -> EosResult<Self> {
        if let Isobar::SubCritical(vle) = isobar {
            let sv = vle.vapor().molar_entropy(Contributions::Total);
            let sl = vle.liquid().molar_entropy(Contributions::Total);
            if molar_entropy > sl && molar_entropy < sv {
                return Ok(Self::TwoPhase(
                    vle.clone(),
                    (molar_entropy - sl).convert_into(sv - sl),
                    mass_flow_rate,
                ));
            }
        }
        let state = StateBuilder::new(eos)
            .pressure(isobar.pressure())
            .molar_entropy(molar_entropy)
            .initial_temperature(initial_temperature)
            .build()?;
        Ok(Self::SinglePhase(Box::new(state), mass_flow_rate))
    }

    pub fn isobaric_temperature_change(
        &self,
        temperature_change: Temperature<f64>,
        density_initialization: DensityInitialization,
    ) -> EosResult<Self> {
        Self::new_pure_temperature_pressure(
            self.eos(),
            self.temperature() + temperature_change,
            self.pressure(),
            self.mass_flow_rate(),
            density_initialization,
        )
    }

    pub fn isobaric_heat_transfer(
        &self,
        heat_flow_rate: Power<f64>,
        isobar: &Isobar<E>,
    ) -> EosResult<Self> {
        ProcessState::new_pure_enthalpy_pressure(
            self.eos(),
            self.molar_enthalpy() + heat_flow_rate / self.mole_flow_rate().unwrap(),
            isobar,
            self.temperature(),
            self.mass_flow_rate(),
        )
    }

    pub fn isentropic_expansion(
        &self,
        isobar: &Isobar<E>,
        initial_temperature: Temperature<f64>,
    ) -> EosResult<Self> {
        ProcessState::new_pure_entropy_pressure(
            self.eos(),
            self.molar_entropy(),
            isobar,
            initial_temperature,
            self.mass_flow_rate(),
        )
    }

    pub fn polytropic_expansion(
        &self,
        isobar: &Isobar<E>,
        molar_work: MolarEnergy<f64>,
        initial_temperature: Temperature<f64>,
    ) -> EosResult<Self> {
        ProcessState::new_pure_enthalpy_pressure(
            self.eos(),
            self.molar_enthalpy() + molar_work,
            isobar,
            initial_temperature,
            self.mass_flow_rate(),
        )
    }

    pub fn split(&self, ratio: f64) -> [Self; 2] {
        match self {
            Self::SinglePhase(state, m) => [
                Self::SinglePhase(state.clone(), Some(ratio * m.unwrap())),
                Self::SinglePhase(state.clone(), Some((1.0 - ratio) * m.unwrap())),
            ],
            Self::TwoPhase(vle, vapor_fraction, m) => [
                Self::TwoPhase(vle.clone(), *vapor_fraction, Some(ratio * m.unwrap())),
                Self::TwoPhase(
                    vle.clone(),
                    *vapor_fraction,
                    Some((1.0 - ratio) * m.unwrap()),
                ),
            ],
        }
    }

    pub fn mix(&self, other: &Self, isobar: &Isobar<E>) -> EosResult<Self> {
        let enthalpy = self.enthalpy_flow_rate().unwrap() + other.enthalpy_flow_rate().unwrap();
        let mole_flow_rate = self.mole_flow_rate().unwrap() + other.mole_flow_rate().unwrap();
        let mass_flow_rate = self.mass_flow_rate().unwrap() + other.mass_flow_rate().unwrap();
        Self::new_pure_enthalpy_pressure(
            self.eos(),
            enthalpy / mole_flow_rate,
            isobar,
            self.temperature(),
            Some(mass_flow_rate),
        )
    }

    pub fn is_two_phase(&self) -> bool {
        match self {
            Self::SinglePhase(_, _) => false,
            Self::TwoPhase(_, _, _) => true,
        }
    }

    pub fn eos(&self) -> &Arc<E> {
        match self {
            Self::SinglePhase(state, _) => &state.eos,
            Self::TwoPhase(vle, _, _) => &vle.vapor().eos,
        }
    }

    pub fn temperature(&self) -> Temperature<f64> {
        match self {
            Self::SinglePhase(state, _) => state.temperature,
            Self::TwoPhase(vle, _, _) => vle.vapor().temperature,
        }
    }

    pub fn density(&self) -> Density<f64> {
        match self {
            Self::SinglePhase(state, _) => state.density,
            Self::TwoPhase(vle, _, _) => vle.vapor().density,
        }
    }

    pub fn pressure(&self) -> Pressure<f64> {
        match self {
            Self::SinglePhase(state, _) => state.pressure(Contributions::Total),
            Self::TwoPhase(vle, _, _) => vle.vapor().pressure(Contributions::Total),
        }
    }

    pub fn specific_enthalpy(&self) -> SpecificEnergy<f64> {
        match self {
            Self::SinglePhase(state, _) => state.specific_enthalpy(Contributions::Total),
            Self::TwoPhase(vle, vapor_fraction, _) => {
                vle.vapor().specific_enthalpy(Contributions::Total) * *vapor_fraction
                    + vle.liquid().specific_enthalpy(Contributions::Total) * (1.0 - vapor_fraction)
            }
        }
    }

    pub fn specific_entropy(&self) -> SpecificEntropy<f64> {
        match self {
            Self::SinglePhase(state, _) => state.specific_entropy(Contributions::Total),
            Self::TwoPhase(vle, vapor_fraction, _) => {
                vle.vapor().specific_entropy(Contributions::Total) * *vapor_fraction
                    + vle.liquid().specific_entropy(Contributions::Total) * (1.0 - vapor_fraction)
            }
        }
    }

    pub fn molar_enthalpy(&self) -> MolarEnergy<f64> {
        match self {
            Self::SinglePhase(state, _) => state.molar_enthalpy(Contributions::Total),
            Self::TwoPhase(vle, vapor_fraction, _) => {
                vle.vapor().molar_enthalpy(Contributions::Total) * *vapor_fraction
                    + vle.liquid().molar_enthalpy(Contributions::Total) * (1.0 - vapor_fraction)
            }
        }
    }

    pub fn molar_entropy(&self) -> MolarEntropy<f64> {
        match self {
            Self::SinglePhase(state, _) => state.molar_entropy(Contributions::Total),
            Self::TwoPhase(vle, vapor_fraction, _) => {
                vle.vapor().molar_entropy(Contributions::Total) * *vapor_fraction
                    + vle.liquid().molar_entropy(Contributions::Total) * (1.0 - vapor_fraction)
            }
        }
    }

    pub fn mass_flow_rate(&self) -> Option<MassFlowRate<f64>> {
        match self {
            Self::SinglePhase(_, m) => *m,
            Self::TwoPhase(_, _, m) => *m,
        }
    }

    pub fn mole_flow_rate(&self) -> Option<MoleFlowRate<f64>> {
        match self {
            Self::SinglePhase(state, m) => m.map(|m| m / state.total_molar_weight()),
            Self::TwoPhase(vle, vapor_fraction, m) => m.map(|m| {
                let recip_molar_weight = *vapor_fraction / vle.vapor().total_molar_weight()
                    + (1.0 - vapor_fraction) / vle.liquid().total_molar_weight();
                m * recip_molar_weight
            }),
        }
    }

    pub fn enthalpy_flow_rate(&self) -> Option<Power<f64>> {
        self.mass_flow_rate().map(|m| m * self.specific_enthalpy())
    }
}
