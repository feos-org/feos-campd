use super::process_state::Isobar;
use super::{Process, ProcessState, Utility};
use feos_core::parameter::ParameterError;
use feos_core::{
    Contributions, DensityInitialization, EosResult, EquationOfState, MolarWeight, State,
};
use quantity::si::*;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::rc::Rc;

const BIG_M: f64 = 500.0;

#[derive(Serialize, Deserialize)]
pub struct OrganicRankineCycleSuperStructureJSON {
    heat_source: Utility,
    #[serde(rename = "min_heat_source_outlet_temperature [°C]")]
    min_heat_source_outlet_temperature: f64,
    isentropic_turbine_efficiency: f64,
    isentropic_pump_efficiency: f64,
    #[serde(rename = "min_abs_pressure [bar]")]
    min_abs_pressure: f64,
    min_red_pressure: f64,
    #[serde(rename = "max_abs_pressure [bar]")]
    max_abs_pressure: f64,
    max_red_pressure: f64,
    cooling: Utility,
    #[serde(rename = "min_temperature_difference_recuperator [K]")]
    min_temperature_difference_recuperator: f64,
    #[serde(rename = "min_heat_transfer_rate_recuperator [kW]")]
    min_heat_transfer_rate_recuperator: f64,
}

#[derive(Clone, Serialize, Deserialize, Debug)]
#[serde(from = "OrganicRankineCycleSuperStructureJSON")]
#[serde(into = "OrganicRankineCycleSuperStructureJSON")]
pub struct OrganicRankineCycleSuperStructure {
    heat_source: Utility,
    min_heat_source_outlet_temperature: SINumber,
    isentropic_turbine_efficiency: f64,
    isentropic_pump_efficiency: f64,
    min_abs_pressure: SINumber,
    min_red_pressure: f64,
    max_abs_pressure: SINumber,
    max_red_pressure: f64,
    cooling: Utility,
    min_temperature_difference_recuperator: SINumber,
    min_heat_transfer_rate_recuperator: SINumber,
}

impl From<OrganicRankineCycleSuperStructureJSON> for OrganicRankineCycleSuperStructure {
    fn from(orc: OrganicRankineCycleSuperStructureJSON) -> Self {
        Self {
            heat_source: orc.heat_source,
            min_heat_source_outlet_temperature: orc.min_heat_source_outlet_temperature * CELSIUS,
            isentropic_turbine_efficiency: orc.isentropic_turbine_efficiency,
            isentropic_pump_efficiency: orc.isentropic_pump_efficiency,
            min_abs_pressure: orc.min_abs_pressure * BAR,
            min_red_pressure: orc.min_red_pressure,
            max_abs_pressure: orc.max_abs_pressure * BAR,
            max_red_pressure: orc.max_red_pressure,
            cooling: orc.cooling,
            min_temperature_difference_recuperator: orc.min_temperature_difference_recuperator
                * KELVIN,
            min_heat_transfer_rate_recuperator: orc.min_heat_transfer_rate_recuperator
                * KILO
                * WATT,
        }
    }
}

impl From<OrganicRankineCycleSuperStructure> for OrganicRankineCycleSuperStructureJSON {
    fn from(orc: OrganicRankineCycleSuperStructure) -> Self {
        Self {
            heat_source: orc.heat_source,
            min_heat_source_outlet_temperature: orc.min_heat_source_outlet_temperature / CELSIUS,
            isentropic_turbine_efficiency: orc.isentropic_turbine_efficiency,
            isentropic_pump_efficiency: orc.isentropic_pump_efficiency,
            min_abs_pressure: orc.min_abs_pressure.to_reduced(BAR).unwrap(),
            min_red_pressure: orc.min_red_pressure,
            max_abs_pressure: orc.max_abs_pressure.to_reduced(BAR).unwrap(),
            max_red_pressure: orc.max_red_pressure,
            cooling: orc.cooling,
            min_temperature_difference_recuperator: orc
                .min_temperature_difference_recuperator
                .to_reduced(KELVIN)
                .unwrap(),
            min_heat_transfer_rate_recuperator: orc
                .min_heat_transfer_rate_recuperator
                .to_reduced(KILO * WATT)
                .unwrap(),
        }
    }
}

#[allow(clippy::too_many_arguments)]
impl OrganicRankineCycleSuperStructure {
    pub fn new(
        heat_source: Utility,
        min_heat_source_outlet_temperature: SINumber,
        isentropic_turbine_efficiency: f64,
        isentropic_pump_efficiency: f64,
        min_abs_pressure: SINumber,
        min_red_pressure: f64,
        max_abs_pressure: SINumber,
        max_red_pressure: f64,
        cooling: Utility,
        min_temperature_difference_recuperator: SINumber,
        min_heat_transfer_rate_recuperator: SINumber,
    ) -> Self {
        Self {
            heat_source,
            min_heat_source_outlet_temperature,
            isentropic_turbine_efficiency,
            isentropic_pump_efficiency,
            min_abs_pressure,
            min_red_pressure,
            max_abs_pressure,
            max_red_pressure,
            cooling,
            min_temperature_difference_recuperator,
            min_heat_transfer_rate_recuperator,
        }
    }

    pub fn from_json<P: AsRef<Path>>(file: P) -> Result<Self, ParameterError> {
        let reader = BufReader::new(File::open(file)?);
        Ok(serde_json::from_reader(reader)?)
    }
}

impl OrganicRankineCycleSuperStructure {
    pub fn variables(&self) -> Vec<[Option<f64>; 2]> {
        vec![
            // Mass flow rate
            [Some(0.0), None],
            // pressure condenser
            [
                Some(self.min_red_pressure.ln()),
                Some(self.max_red_pressure.ln()),
            ],
            // pressure evaporator
            [
                Some(self.min_red_pressure.ln()),
                Some(self.max_red_pressure.ln()),
            ],
            // superheating
            [Some(0.0), Some(1.0)],
            // heat flow recuperator
            [Some(0.0), None],
        ]
    }

    pub fn binary_variables(&self) -> usize {
        1
    }

    pub fn constraints(&self) -> Vec<[Option<f64>; 3]> {
        vec![
            // Pinch constraints
            [Some(0.0), None, None],
            [Some(0.0), None, None],
            [Some(0.0), None, None],
            [Some(0.0), None, None],
            [Some(0.0), None, None],
            [Some(0.0), None, None],
            [Some(0.0), None, None],
            // Absolute pressure (bar)
            [Some(0.0), Some(1.0), None],
            [Some(0.0), Some(1.0), None],
            // Recuperator
            [Some(0.0), None, None],
            [Some(0.0), None, None],
            [Some(0.0), None, None],
            [Some(0.0), None, None],
            // Heat source outlet
            [Some(0.0), None, None],
        ]
    }

    pub fn solve<E: EquationOfState + MolarWeight<SIUnit>>(
        &self,
        eos: &Rc<E>,
        x: &[f64],
    ) -> EosResult<(Process<E>, f64, Vec<f64>)> {
        // unpack variables
        let mwf = x[0] * KILOGRAM / SECOND;
        let p_cond_red = x[1].exp();
        let p_evap_red = x[2].exp();
        let dt_sh = x[3] * (self.heat_source.temperature - self.min_heat_source_outlet_temperature);
        let q_rec = x[4] * MEGA * WATT;
        let recuperator = x[5];

        // Calculate pressures
        let cp = State::critical_point(eos, None, None, Default::default())?;
        let p_crit = cp.pressure(Contributions::Total);
        let p_cond = p_cond_red * p_crit;
        let p_evap = p_evap_red * p_crit;

        // Calculate isobars
        let isobar_cond = Isobar::new(eos, p_cond);
        let isobar_evap = Isobar::new(eos, p_evap);

        // Initialize process
        let mut process = Process::new();

        // Calculate pump
        let feed = process.add_state(ProcessState::SinglePhase(
            Box::new(isobar_cond.saturated_liquid()?.clone()),
            Some(mwf),
        ));
        let pump = process.pump(feed, &isobar_evap, self.isentropic_pump_efficiency)?;

        // Calculate turbine
        let feed = process.add_state(ProcessState::new_pure_temperature_pressure(
            eos,
            isobar_evap.saturated_vapor()?.temperature + dt_sh,
            p_evap,
            Some(mwf),
            DensityInitialization::Vapor,
        )?);
        let turbine = process.turbine(feed, &isobar_cond, self.isentropic_turbine_efficiency)?;

        // Calculate Recuperator
        let [recuperator_liquid, recuperator_vapor] =
            process.recuperator(pump.out(), turbine.out(), &isobar_evap, &isobar_cond, q_rec)?;

        // Calculate evaporator
        let evaporator = process.evaporator(recuperator_liquid.out(), &isobar_evap, Some(dt_sh))?;
        process.add_utility(&evaporator, self.heat_source);

        // Calculate condenser
        let condenser = process.total_condenser(recuperator_vapor.out(), &isobar_cond, None)?;
        process.add_utility(&condenser, self.cooling);

        // Target
        let target = process.net_power().unwrap().to_reduced(KILO * WATT)?;

        // Pinch constraints
        let mut constraints = process.pinch_constraints();

        // Absolute pressure constraints
        constraints
            .push((p_cond - self.min_abs_pressure).to_reduced(p_evap - self.min_abs_pressure)?);
        constraints.push((p_evap - p_cond).to_reduced(self.max_abs_pressure - p_cond)?);

        // Recuperator constraints
        constraints.push(
            (process[turbine.out()].temperature()
                - process[recuperator_liquid.out()].temperature())
            .to_reduced(self.min_temperature_difference_recuperator)?
                + (1.0 - recuperator) * BIG_M
                - recuperator,
        );
        constraints.push(
            (process[recuperator_vapor.out()].temperature() - process[pump.out()].temperature())
                .to_reduced(self.min_temperature_difference_recuperator)?
                + (1.0 - recuperator) * BIG_M
                - recuperator,
        );
        constraints.push(q_rec.to_reduced(self.min_heat_transfer_rate_recuperator)? - recuperator);
        constraints
            .push(recuperator * BIG_M - q_rec.to_reduced(self.min_heat_transfer_rate_recuperator)?);

        // Heat source outlet constraint
        constraints.push(
            process.graph[recuperator_liquid.out()]
                .utility_temperature
                .unwrap()[0]
                .to_reduced(self.min_heat_source_outlet_temperature)?
                - 1.0,
        );

        Ok((process, target, constraints))
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use feos_core::parameter::{IdentifierOption, Parameter};
    use feos_pcsaft::{PcSaft, PcSaftParameters};

    #[test]
    fn test_json() {
        let orc: OrganicRankineCycleSuperStructure =
            serde_json::from_reader(BufReader::new(File::open("orc.json").unwrap())).unwrap();
        let params = PcSaftParameters::from_json(
            vec!["propane"],
            "gross2001.json",
            None,
            IdentifierOption::Name,
        )
        .unwrap();
        let eos = Rc::new(PcSaft::new(Rc::new(params)));
        orc.solve(&eos, &[1.15, -4.5, -1.2, 0.8, 2.0, 0.5]).unwrap();
        // let u1 = Utility::ConstantTemperature;
        // let u2 = Utility::HeatCapacityRate(6.0 * WATT / KELVIN);
        // let u3 = Utility::OutletTemperature(25.0 * CELSIUS);
        // println!("{}", serde_json::to_string(&u1).unwrap());
        // println!("{}", serde_json::to_string(&u2).unwrap());
        // println!("{}", serde_json::to_string(&u3).unwrap());
        // println!("{} {}", 0.0f64.signum(), (-0.0f64).signum());
        // assert!(false);
    }
}
