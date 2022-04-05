use super::{Isobar, Process, ProcessState, ProcessStep, StatePoint};
use feos_core::{Contributions, DensityInitialization, EosResult, EquationOfState, MolarWeight};
use quantity::si::*;

pub struct Equipment {
    pub states: Vec<StatePoint>,
}

impl Equipment {
    pub fn out(&self) -> StatePoint {
        *self.states.last().unwrap()
    }
}

impl<E: EquationOfState + MolarWeight<SIUnit>> Process<E> {
    pub fn total_condenser(
        &mut self,
        feed: StatePoint,
        isobar: &Isobar<E>,
        subcooling: Option<SINumber>,
    ) -> EosResult<Equipment> {
        let mass_flow_rate = self[feed].mass_flow_rate();
        let vle = isobar.phase_equilibrium()?;

        let vapor_state =
            if self[feed].molar_enthalpy() < vle.vapor().molar_enthalpy(Contributions::Total) {
                // Add a dummy step to make sure the number of pinch constraints is always the same
                self[feed].clone()
            } else {
                ProcessState::TwoPhase(vle.clone(), 1.0, mass_flow_rate)
            };
        let liquid_state =
            if self[feed].molar_enthalpy() < vle.liquid().molar_enthalpy(Contributions::Total) {
                // Add a dummy step to make sure the number of pinch constraints is always the same
                self[feed].clone()
            } else {
                ProcessState::TwoPhase(vle.clone(), 0.0, mass_flow_rate)
            };
        let subcooled_state = subcooling
            .map(|dt| liquid_state.isobaric_temperature_change(-dt, DensityInitialization::Liquid))
            .transpose()?;

        let vapor = self.add_step(feed, vapor_state, ProcessStep::Isobaric);
        let liquid = self.add_step(vapor, liquid_state, ProcessStep::PhaseChange);
        let subcooled = subcooled_state.map(|s| self.add_step(liquid, s, ProcessStep::Isobaric));

        let states = [feed, vapor, liquid]
            .into_iter()
            .chain(subcooled.into_iter())
            .collect();

        Ok(Equipment { states })
    }

    pub fn evaporator(
        &mut self,
        feed: StatePoint,
        isobar: &Isobar<E>,
        superheating: Option<SINumber>,
    ) -> EosResult<Equipment> {
        let mass_flow_rate = self[feed].mass_flow_rate();
        let vle = isobar.phase_equilibrium()?;

        let liquid_state =
            if self[feed].molar_enthalpy() > vle.liquid().molar_enthalpy(Contributions::Total) {
                // Add a dummy step to make sure the number of pinch constraints is always the same
                self[feed].clone()
            } else {
                ProcessState::TwoPhase(vle.clone(), 0.0, mass_flow_rate)
            };
        let vapor_state =
            if self[feed].molar_enthalpy() > vle.vapor().molar_enthalpy(Contributions::Total) {
                // Add a dummy step to make sure the number of pinch constraints is always the same
                self[feed].clone()
            } else {
                ProcessState::TwoPhase(vle.clone(), 0.0, mass_flow_rate)
            };
        let superheated_state = superheating
            .map(|dt| vapor_state.isobaric_temperature_change(dt, DensityInitialization::Vapor))
            .transpose()?;

        let liquid = self.add_step(feed, liquid_state, ProcessStep::Isobaric);
        let vapor = self.add_step(liquid, vapor_state, ProcessStep::PhaseChange);
        let superheated = superheated_state.map(|s| self.add_step(vapor, s, ProcessStep::Isobaric));

        let states = [feed, liquid, vapor]
            .into_iter()
            .chain(superheated.into_iter())
            .collect();

        Ok(Equipment { states })
    }

    pub fn pump(
        &mut self,
        feed: StatePoint,
        isobar: &Isobar<E>,
        efficiency: f64,
    ) -> EosResult<Equipment> {
        self.pressure_changer(feed, isobar, Some(1.0 / efficiency))
    }

    pub fn turbine(
        &mut self,
        feed: StatePoint,
        isobar: &Isobar<E>,
        efficiency: f64,
    ) -> EosResult<Equipment> {
        self.pressure_changer(feed, isobar, Some(efficiency))
    }

    pub fn compressor(
        &mut self,
        feed: StatePoint,
        isobar: &Isobar<E>,
        efficiency: f64,
    ) -> EosResult<Equipment> {
        self.pressure_changer(feed, isobar, Some(1.0 / efficiency))
    }

    pub fn throttle(&mut self, feed: StatePoint, isobar: &Isobar<E>) -> EosResult<Equipment> {
        self.pressure_changer(feed, isobar, None)
    }

    fn pressure_changer(
        &mut self,
        feed: StatePoint,
        isobar: &Isobar<E>,
        efficiency: Option<f64>,
    ) -> EosResult<Equipment> {
        let mut initial_temperature = self[feed].temperature();
        let molar_work = if let Some(efficiency) = efficiency {
            let state2s = self[feed].isentropic_expansion(isobar, initial_temperature)?;
            initial_temperature = state2s.temperature();
            let h1 = self[feed].molar_enthalpy();
            let h2s = state2s.molar_enthalpy();
            efficiency * (h2s - h1)
        } else {
            0.0 * WATT
        };
        let out = self[feed].polytropic_expansion(isobar, molar_work, initial_temperature)?;
        let out = self.add_step(feed, out, ProcessStep::Polytropic);

        Ok(Equipment {
            states: vec![feed, out],
        })
    }

    pub fn recuperator(
        &mut self,
        liquid: StatePoint,
        vapor: StatePoint,
        isobar_liquid: &Isobar<E>,
        isobar_vapor: &Isobar<E>,
        transferred_heat: SINumber,
    ) -> EosResult<[Equipment; 2]> {
        let liquid_out = self[liquid].isobaric_heat_transfer(transferred_heat, isobar_liquid)?;
        let vapor_out = self[vapor].isobaric_heat_transfer(-transferred_heat, isobar_vapor)?;

        let liquid_out = self.add_step(liquid, liquid_out, ProcessStep::Isobaric);
        let vapor_out = self.add_step(vapor, vapor_out, ProcessStep::Isobaric);
        Ok([
            Equipment {
                states: vec![liquid, liquid_out],
            },
            Equipment {
                states: vec![vapor, vapor_out],
            },
        ])
    }

    pub fn split(&mut self, feed: StatePoint, ratio: f64) -> [StatePoint; 2] {
        let [out1, out2] = self[feed].split(ratio);
        [self.add_state(out1), self.add_state(out2)]
    }

    pub fn isobaric_mixer(
        &mut self,
        feed1: StatePoint,
        feed2: StatePoint,
        isobar: &Isobar<E>,
    ) -> EosResult<StatePoint> {
        let mix = self[feed1].mix(&self[feed2], isobar)?;
        let out = self.add_step(feed1, mix, ProcessStep::Isobaric);
        self.add_connection(feed2, out, ProcessStep::Isobaric);
        Ok(out)
    }
}
