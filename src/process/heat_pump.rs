use super::{Equipment, Process, ProcessState};
use feos_core::parameter::ParameterError;
use feos_core::{
    Contributions, EosResult, EquationOfState, MolarWeight, PhaseEquilibrium, State, StateBuilder,
};
use quantity::si::*;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::rc::Rc;

#[derive(Serialize, Deserialize)]
pub struct HeatPumpJSON {
    #[serde(rename = "heat_sink_temperature_in [°C]")]
    heat_sink_temperature_in: f64,
    #[serde(rename = "heat_sink_temperature_out [°C]")]
    heat_sink_temperature_out: f64,
    #[serde(rename = "heating_power [kW]")]
    heating_power: f64,
    #[serde(rename = "heat_source_temperature [°C]")]
    heat_source_temperature: f64,
    isentropic_compressor_efficiency: f64,
    #[serde(rename = "min_approach_temperature [K]")]
    min_approach_temperature: f64,
    #[serde(rename = "min_abs_pressure [bar]")]
    min_abs_pressure: f64,
    min_red_pressure: f64,
    #[serde(rename = "max_abs_pressure [bar]")]
    max_abs_pressure: f64,
    max_red_pressure: f64,
    #[serde(rename = "min_superheating [K]")]
    min_superheating: f64,
    #[serde(rename = "max_subcooling [K]")]
    max_subcooling: f64,
}

#[derive(Clone, Serialize, Deserialize, Debug)]
#[serde(from = "HeatPumpJSON")]
#[serde(into = "HeatPumpJSON")]
pub struct HeatPump {
    heat_sink_temperature_in: SINumber,
    heat_sink_temperature_out: SINumber,
    heating_power: SINumber,
    heat_source_temperature: SINumber,
    isentropic_compressor_efficiency: f64,
    min_approach_temperature: SINumber,
    min_abs_pressure: SINumber,
    min_red_pressure: f64,
    max_abs_pressure: SINumber,
    max_red_pressure: f64,
    min_superheating: SINumber,
    max_subcooling: SINumber,
}

impl From<HeatPumpJSON> for HeatPump {
    fn from(heat_pump: HeatPumpJSON) -> Self {
        Self {
            heat_sink_temperature_in: heat_pump.heat_sink_temperature_in * CELSIUS,
            heat_sink_temperature_out: heat_pump.heat_sink_temperature_out * CELSIUS,
            heat_source_temperature: heat_pump.heat_source_temperature * CELSIUS,
            heating_power: heat_pump.heating_power * KILO * WATT,
            isentropic_compressor_efficiency: heat_pump.isentropic_compressor_efficiency,
            min_approach_temperature: heat_pump.min_approach_temperature * KELVIN,
            min_abs_pressure: heat_pump.min_abs_pressure * BAR,
            min_red_pressure: heat_pump.min_red_pressure,
            max_abs_pressure: heat_pump.max_abs_pressure * BAR,
            max_red_pressure: heat_pump.max_red_pressure,
            min_superheating: heat_pump.min_superheating * KELVIN,
            max_subcooling: heat_pump.max_subcooling * KELVIN,
        }
    }
}

impl From<HeatPump> for HeatPumpJSON {
    fn from(heat_pump: HeatPump) -> Self {
        Self {
            heat_sink_temperature_in: heat_pump.heat_sink_temperature_in / CELSIUS,
            heat_sink_temperature_out: heat_pump.heat_sink_temperature_out / CELSIUS,
            heat_source_temperature: heat_pump.heat_source_temperature / CELSIUS,
            heating_power: heat_pump.heating_power.to_reduced(KILO * WATT).unwrap(),
            isentropic_compressor_efficiency: heat_pump.isentropic_compressor_efficiency,
            min_approach_temperature: heat_pump
                .min_approach_temperature
                .to_reduced(KELVIN)
                .unwrap(),
            min_abs_pressure: heat_pump.min_abs_pressure.to_reduced(BAR).unwrap(),
            min_red_pressure: heat_pump.min_red_pressure,
            max_abs_pressure: heat_pump.max_abs_pressure.to_reduced(BAR).unwrap(),
            max_red_pressure: heat_pump.max_red_pressure,
            min_superheating: heat_pump.min_superheating.to_reduced(KELVIN).unwrap(),
            max_subcooling: heat_pump.max_subcooling.to_reduced(KELVIN).unwrap(),
        }
    }
}

#[allow(clippy::too_many_arguments)]
impl HeatPump {
    pub fn new(
        heat_sink_temperature_in: SINumber,
        heat_sink_temperature_out: SINumber,
        heating_power: SINumber,
        heat_source_temperature: SINumber,
        isentropic_compressor_efficiency: f64,
        min_approach_temperature: SINumber,
        min_abs_pressure: SINumber,
        min_red_pressure: f64,
        max_abs_pressure: SINumber,
        max_red_pressure: f64,
        min_superheating: SINumber,
        max_subcooling: SINumber,
    ) -> Self {
        Self {
            heat_sink_temperature_in,
            heat_sink_temperature_out,
            heating_power,
            heat_source_temperature,
            isentropic_compressor_efficiency,
            min_approach_temperature,
            min_abs_pressure,
            min_red_pressure,
            max_abs_pressure,
            max_red_pressure,
            min_superheating,
            max_subcooling,
        }
    }

    pub fn from_json<P: AsRef<Path>>(file: P) -> Result<Self, ParameterError> {
        let reader = BufReader::new(File::open(file)?);
        Ok(serde_json::from_reader(reader)?)
    }
}

impl HeatPump {
    pub fn variables(&self) -> Vec<[Option<f64>; 2]> {
        vec![
            // Evaporator temperature: fully determined by source temperature, superheating and minimum approach temperature
            [
                Some(
                    (self.heat_source_temperature
                        - self.min_superheating
                        - self.min_approach_temperature)
                        .to_reduced(KELVIN)
                        .unwrap()
                        / 10.,
                ),
                Some(
                    (self.heat_source_temperature
                        - self.min_superheating
                        - self.min_approach_temperature)
                        .to_reduced(KELVIN)
                        .unwrap()
                        / 10.,
                ),
            ],
            // Condenser temperature
            [
                Some(self.heat_sink_temperature_in.to_reduced(KELVIN).unwrap()),
                None,
            ],
            // Superheating
            [
                Some(self.min_superheating.to_reduced(KELVIN).unwrap()),
                Some(self.min_superheating.to_reduced(KELVIN).unwrap()),
            ],
            [
                Some(0.0),
                Some(self.max_subcooling.to_reduced(KELVIN).unwrap()),
            ],
            [Some(2.5), Some(9.0)],
        ]
    }

    pub fn constraints(&self) -> Vec<[Option<f64>; 3]> {
        vec![
            // "Process constraint": Evaporator temperature lower than condenser temperature: Tevap - Tcond <= 0
            [None, Some(0.0), None],
            // Pinch constraints
            [Some(1.0), None, None],
            [Some(1.0), None, None],
            [Some(1.0), None, None],
            [Some(1.0), None, None],
            // 2-phase compression prohibited: (h_comp - h'') / h'' >= 0
            [Some(0.0), None, None],
            // Maximum reduced pressure
            [None, Some(self.max_red_pressure), None],
            // Minimum pressure
            [Some(0.0), None, None],
            // Fixed heating power
            [Some(-0.02), Some(0.02), None],
        ]
    }

    pub fn solve<E: EquationOfState + MolarWeight<SIUnit>>(
        &self,
        eos: &Rc<E>,
        x: &[f64],
    ) -> EosResult<(Process<E>, f64, Vec<f64>)> {
        // unpack variables
        let t_evap = x[0] * 10.0 * KELVIN;
        let t_cond = x[1] * 10.0 * KELVIN;
        let dt_sh = x[2] * KELVIN;
        let dt_sc = x[3] * KELVIN;
        let d = x[4] * 10.0 * MILLI * METER;

        // Calculate phase equilibria
        let vle_cond = PhaseEquilibrium::pure_t(eos, t_cond, None, Default::default())?;
        let vle_cond = Rc::new(vle_cond);
        let vle_evap = PhaseEquilibrium::pure_t(eos, t_evap, None, Default::default())?;
        let vle_evap = Rc::new(vle_evap);

        // Calculate pressures
        let cp = State::critical_point(eos, None, None, Default::default())?;
        let p_crit = cp.pressure(Contributions::Total);
        let p_cond = vle_cond.liquid().pressure(Contributions::Total);
        let p_evap = vle_evap.liquid().pressure(Contributions::Total);

        // Initialize process
        let mut process = Process::new();
        let feed = process.initialize(ProcessState::SinglePhase(Box::new(
            StateBuilder::new(eos)
                .temperature(t_cond - dt_sc)
                .pressure(p_cond)
                .liquid()
                .build()?,
        )));

        // Calculate throttle
        let throttle = Equipment::throttle(&mut process, feed, &vle_evap, None)?;

        // Calculate evaporator
        let evaporator =
            Equipment::evaporator(&mut process, throttle.out(), None, &vle_evap, Some(dt_sh))?;
        //process.add_utility(
        //    &evaporator,
        //    self.heat_source_temperature,
        //    self.heat_source_mass_flow_rate,
        //    self.heat_source_heat_capacity,
        //    self.min_approach_temperature,
        //);

        // Calculate compressor
        let compressor = Equipment::compressor(
            &mut process,
            evaporator.out(),
            &vle_cond,
            None,
            self.isentropic_compressor_efficiency,
        )?;

        // Calculate condenser
        let condenser = Equipment::total_condenser(
            &mut process,
            compressor.out(),
            None,
            &vle_cond,
            Some(dt_sc),
        )?;
        process.add_utility_heat_sink(
            &condenser,
            self.heat_sink_temperature_in,
            self.heat_sink_temperature_out,
            self.heating_power,
            self.min_approach_temperature,
        );

        // Target
        // let target =
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
