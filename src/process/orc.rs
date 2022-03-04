use super::{Equipment, Process, ProcessState};
use feos_core::parameter::ParameterError;
use feos_core::{Contributions, EosResult, EquationOfState, MolarWeight, PhaseEquilibrium, State};
use quantity::si::*;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::rc::Rc;

#[derive(Serialize, Deserialize)]
pub struct OrganicRankineCycleJSON {
    #[serde(rename = "heat_source_temperature [°C]")]
    heat_source_temperature: f64,
    #[serde(rename = "heat_source_mass_flow_rate [kg/s]")]
    heat_source_mass_flow_rate: f64,
    #[serde(rename = "heat_source_heat_capacity [kJ/kg/K]")]
    heat_source_heat_capacity: f64,
    isentropic_turbine_efficiency: f64,
    isentropic_pump_efficiency: f64,
    #[serde(rename = "min_approach_temperature [K]")]
    min_approach_temperature: f64,
    #[serde(rename = "min_abs_pressure [bar]")]
    min_abs_pressure: f64,
    min_red_pressure: f64,
    #[serde(rename = "max_abs_pressure [bar]")]
    max_abs_pressure: f64,
    max_red_pressure: f64,
    #[serde(rename = "min_cooling_temperature [°C]")]
    min_cooling_temperature: f64,
}

#[derive(Clone, Serialize, Deserialize, Debug)]
#[serde(from = "OrganicRankineCycleJSON")]
#[serde(into = "OrganicRankineCycleJSON")]
pub struct OrganicRankineCycle {
    heat_source_temperature: SINumber,
    heat_source_mass_flow_rate: SINumber,
    heat_source_heat_capacity: SINumber,
    isentropic_turbine_efficiency: f64,
    isentropic_pump_efficiency: f64,
    min_approach_temperature: SINumber,
    min_abs_pressure: SINumber,
    min_red_pressure: f64,
    max_abs_pressure: SINumber,
    max_red_pressure: f64,
    min_cooling_temperature: SINumber,
}

impl From<OrganicRankineCycleJSON> for OrganicRankineCycle {
    fn from(orc: OrganicRankineCycleJSON) -> Self {
        Self {
            heat_source_temperature: orc.heat_source_temperature * CELSIUS,
            heat_source_mass_flow_rate: orc.heat_source_mass_flow_rate * KILOGRAM / SECOND,
            heat_source_heat_capacity: orc.heat_source_heat_capacity * KILO * JOULE
                / KILOGRAM
                / KELVIN,
            isentropic_turbine_efficiency: orc.isentropic_turbine_efficiency,
            isentropic_pump_efficiency: orc.isentropic_pump_efficiency,
            min_approach_temperature: orc.min_approach_temperature * KELVIN,
            min_abs_pressure: orc.min_abs_pressure * BAR,
            min_red_pressure: orc.min_red_pressure,
            max_abs_pressure: orc.max_abs_pressure * BAR,
            max_red_pressure: orc.max_red_pressure,
            min_cooling_temperature: orc.min_cooling_temperature * CELSIUS,
        }
    }
}

impl From<OrganicRankineCycle> for OrganicRankineCycleJSON {
    fn from(orc: OrganicRankineCycle) -> Self {
        Self {
            heat_source_temperature: orc.heat_source_temperature / CELSIUS,
            heat_source_mass_flow_rate: orc
                .heat_source_mass_flow_rate
                .to_reduced(KILOGRAM / SECOND)
                .unwrap(),
            heat_source_heat_capacity: orc
                .heat_source_heat_capacity
                .to_reduced(KILO * JOULE / KILOGRAM / KELVIN)
                .unwrap(),
            isentropic_turbine_efficiency: orc.isentropic_turbine_efficiency,
            isentropic_pump_efficiency: orc.isentropic_pump_efficiency,
            min_approach_temperature: orc.min_approach_temperature.to_reduced(KELVIN).unwrap(),
            min_abs_pressure: orc.min_abs_pressure.to_reduced(BAR).unwrap(),
            min_red_pressure: orc.min_red_pressure,
            max_abs_pressure: orc.max_abs_pressure.to_reduced(BAR).unwrap(),
            max_red_pressure: orc.max_red_pressure,
            min_cooling_temperature: orc.min_cooling_temperature / CELSIUS,
        }
    }
}

#[allow(clippy::too_many_arguments)]
impl OrganicRankineCycle {
    pub fn new(
        heat_source_temperature: SINumber,
        heat_source_mass_flow_rate: SINumber,
        heat_source_heat_capacity: SINumber,
        isentropic_turbine_efficiency: f64,
        isentropic_pump_efficiency: f64,
        min_approach_temperature: SINumber,
        min_abs_pressure: SINumber,
        min_red_pressure: f64,
        max_abs_pressure: SINumber,
        max_red_pressure: f64,
        min_cooling_temperature: SINumber,
    ) -> Self {
        Self {
            heat_source_temperature,
            heat_source_mass_flow_rate,
            heat_source_heat_capacity,
            isentropic_turbine_efficiency,
            isentropic_pump_efficiency,
            min_approach_temperature,
            min_abs_pressure,
            min_red_pressure,
            max_abs_pressure,
            max_red_pressure,
            min_cooling_temperature,
        }
    }

    pub fn from_json<P: AsRef<Path>>(file: P) -> Result<Self, ParameterError> {
        let reader = BufReader::new(File::open(file)?);
        Ok(serde_json::from_reader(reader)?)
    }
}

impl OrganicRankineCycle {
    pub fn variables(&self) -> Vec<[Option<f64>; 2]> {
        vec![
            [Some(0.0), None],
            [
                Some(self.min_red_pressure.ln()),
                Some(self.max_red_pressure.ln()),
            ],
            [
                Some(self.min_red_pressure.ln()),
                Some(self.max_red_pressure.ln()),
            ],
            [Some(0.0), None],
        ]
    }

    pub fn constraints(&self) -> Vec<[Option<f64>; 3]> {
        vec![
            // Pinch constraints
            [Some(1.0), None, None],
            [Some(1.0), None, None],
            [Some(1.0), None, None],
            [Some(1.0), None, None],
            // Min. cooling temperature
            [Some(1.0), None, None],
            // Absolute pressure (bar)
            [
                Some(self.min_abs_pressure.to_reduced(BAR).unwrap()),
                Some(self.max_abs_pressure.to_reduced(BAR).unwrap()),
                None,
            ],
            [
                Some(self.min_abs_pressure.to_reduced(BAR).unwrap()),
                Some(self.max_abs_pressure.to_reduced(BAR).unwrap()),
                None,
            ],
            [Some(1.0), None, None],
        ]
    }

    pub fn solve<E: EquationOfState + MolarWeight<SIUnit>>(
        &self,
        eos: &Rc<E>,
        x: &[f64],
    ) -> EosResult<(Process<E>, f64, Vec<f64>)> {
        // unpack variables
        let mwf = x[0] * self.heat_source_mass_flow_rate;
        let p_cond_red = x[1].exp();
        let p_evap_red = x[2].exp();
        let dt_sh = x[3] * 50.0 * KELVIN;

        // Calculate pressures
        let cp = State::critical_point(eos, None, None, Default::default())?;
        let p_crit = cp.pressure(Contributions::Total);
        let p_cond = p_cond_red * p_crit;
        let p_evap = p_evap_red * p_crit;

        // Calculate phase equilibria
        let vle_cond = PhaseEquilibrium::pure_p(eos, p_cond, None, Default::default())?;
        let vle_cond = Rc::new(vle_cond);
        let vle_evap = PhaseEquilibrium::pure_p(eos, p_evap, None, Default::default())?;
        let vle_evap = Rc::new(vle_evap);

        // Initialize process
        let mut process = Process::new();
        let feed = process.initialize(ProcessState::SinglePhase(Box::new(
            vle_cond.liquid().clone(),
        )));

        // Calculate pump
        let pump = Equipment::pump(
            &mut process,
            feed,
            &vle_evap,
            Some(mwf),
            self.isentropic_pump_efficiency,
        )?;

        // Calculate evaporator
        let evaporator =
            Equipment::evaporator(&mut process, pump.out(), Some(mwf), &vle_evap, Some(dt_sh))?;
        process.add_utility(
            &evaporator,
            self.heat_source_temperature,
            self.heat_source_mass_flow_rate,
            self.heat_source_heat_capacity,
            self.min_approach_temperature,
        );

        // Calculate turbine
        let turbine = Equipment::turbine(
            &mut process,
            evaporator.out(),
            &vle_cond,
            Some(mwf),
            self.isentropic_turbine_efficiency,
        )?;

        // Calculate condenser
        Equipment::total_condenser(&mut process, turbine.out(), Some(mwf), &vle_cond, None)?;

        // Target
        let target = process.net_power().to_reduced(MEGA * WATT)?;

        // Pitch constraints
        let mut constraints = process.pinch_constraints();

        // Min. cooling temperature
        constraints.push(
            vle_cond
                .vapor()
                .temperature
                .to_reduced(self.min_cooling_temperature)?,
        );

        // Absolute pressure constraints
        constraints.push(p_cond.to_reduced(BAR)?);
        constraints.push(p_evap.to_reduced(BAR)?);
        constraints.push(p_evap.to_reduced(p_cond)?);

        Ok((process, target, constraints))
    }
}
