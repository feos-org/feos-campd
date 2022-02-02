use super::{Process, ProcessState, ProcessStep, StatePoint};
use feos_core::{
    Contributions, EosResult, EquationOfState, MolarWeight, PhaseEquilibrium, StateBuilder,
};
use quantity::si::*;
use std::rc::Rc;

pub struct Equipment {
    pub states: Vec<StatePoint>,
}

impl Equipment {
    pub fn out(&self) -> StatePoint {
        *self.states.last().unwrap()
    }

    pub fn total_condenser<E: EquationOfState + MolarWeight<SIUnit>>(
        process: &mut Process<E>,
        mut feed: StatePoint,
        mass_flow_rate: SINumber,
        vle: &Rc<PhaseEquilibrium<SIUnit, E, 2>>,
    ) -> EosResult<Self> {
        let mut states = Vec::new();
        if !process[feed].is_two_phase() {
            let (vapor, _) = process.add_step(
                feed,
                ProcessState::TwoPhase(vle.clone(), 1.0),
                ProcessStep::isobaric(mass_flow_rate),
            );
            feed = vapor;
            states.push(vapor)
        }

        let (liquid, _) = process.add_step(
            feed,
            ProcessState::TwoPhase(vle.clone(), 0.0),
            ProcessStep::phase_change(mass_flow_rate),
        );
        states.push(liquid);
        Ok(Self { states })
    }

    pub fn evaporator<E: EquationOfState + MolarWeight<SIUnit>>(
        process: &mut Process<E>,
        feed: StatePoint,
        mass_flow_rate: SINumber,
        vle: &Rc<PhaseEquilibrium<SIUnit, E, 2>>,
        superheating: Option<SINumber>,
    ) -> EosResult<Self> {
        let (liquid, _) = process.add_step(
            feed,
            ProcessState::TwoPhase(vle.clone(), 0.0),
            ProcessStep::isobaric(mass_flow_rate),
        );
        let (vapor, _) = process.add_step(
            liquid,
            ProcessState::TwoPhase(vle.clone(), 1.0),
            ProcessStep::phase_change(mass_flow_rate),
        );

        let mut states = vec![feed, liquid, vapor];
        if let Some(dt) = superheating {
            let pressure = process[feed].pressure();
            let state_out = StateBuilder::new(process[feed].eos())
                .pressure(pressure)
                .temperature(vle.vapor().temperature + dt)
                .vapor()
                .build()?;
            let (superheated, _) = process.add_step(
                vapor,
                ProcessState::SinglePhase(Box::new(state_out)),
                ProcessStep::isobaric(mass_flow_rate),
            );
            states.push(superheated);
        }

        Ok(Self { states })
    }

    pub fn pump<E: EquationOfState + MolarWeight<SIUnit>>(
        process: &mut Process<E>,
        feed: StatePoint,
        vle_out: &Rc<PhaseEquilibrium<SIUnit, E, 2>>,
        mass_flow_rate: SINumber,
        efficiency: f64,
    ) -> EosResult<Self> {
        Self::pressure_changer(
            process,
            feed,
            vle_out,
            mass_flow_rate,
            1.0 / efficiency,
            process[feed].temperature(),
        )
    }

    pub fn turbine<E: EquationOfState + MolarWeight<SIUnit>>(
        process: &mut Process<E>,
        feed: StatePoint,
        vle_out: &Rc<PhaseEquilibrium<SIUnit, E, 2>>,
        mass_flow_rate: SINumber,
        efficiency: f64,
    ) -> EosResult<Self> {
        Self::pressure_changer(
            process,
            feed,
            vle_out,
            mass_flow_rate,
            efficiency,
            vle_out.vapor().temperature,
        )
    }

    fn pressure_changer<E: EquationOfState + MolarWeight<SIUnit>>(
        process: &mut Process<E>,
        feed: StatePoint,
        vle_out: &Rc<PhaseEquilibrium<SIUnit, E, 2>>,
        mass_flow_rate: SINumber,
        efficiency: f64,
        initial_temperature: SINumber,
    ) -> EosResult<Self> {
        let pressure_out = vle_out.vapor().pressure(Contributions::Total);
        let h1 = process[feed].molar_enthalpy();
        let s1 = process[feed].molar_entropy();
        let state2s = ProcessState::new_pure_entropy_pressure(
            process[feed].eos(),
            s1,
            pressure_out,
            vle_out,
            initial_temperature,
        )?;
        let h2s = state2s.molar_enthalpy();
        let h2 = h1 - efficiency * (h1 - h2s);
        let out = ProcessState::new_pure_enthalpy_pressure(
            process[feed].eos(),
            h2,
            pressure_out,
            vle_out,
            state2s.temperature(),
        )?;
        let (out, _) = process.add_step(feed, out, ProcessStep::polytropic(mass_flow_rate));

        Ok(Self {
            states: vec![feed, out],
        })
    }
}
