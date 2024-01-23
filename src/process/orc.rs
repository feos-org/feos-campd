use crate::variables::ContinuousVariables;
use crate::Variable;

use super::{Isobar, Process, ProcessModel, ProcessState, Utility};
use feos::core::parameter::ParameterError;
use feos::core::si::*;
use feos::core::{Contributions, EosResult, IdealGas, Residual, State};
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::sync::Arc;

#[derive(Serialize, Deserialize)]
struct OrganicRankineCycleJSON {
    heat_source: Utility,
    isentropic_turbine_efficiency: f64,
    isentropic_pump_efficiency: f64,
    #[serde(rename = "min_abs_pressure [bar]")]
    min_abs_pressure: f64,
    min_red_pressure: f64,
    #[serde(rename = "max_abs_pressure [bar]")]
    max_abs_pressure: f64,
    max_red_pressure: f64,
    cooling: Utility,
}

/// Implementation of a simplified process model for an organic Rankine cycle.
#[derive(Clone, Serialize, Deserialize, Debug)]
#[serde(from = "OrganicRankineCycleJSON")]
#[serde(into = "OrganicRankineCycleJSON")]
pub struct OrganicRankineCycle {
    pub heat_source: Utility,
    pub isentropic_turbine_efficiency: f64,
    pub isentropic_pump_efficiency: f64,
    pub min_abs_pressure: Pressure,
    pub min_red_pressure: f64,
    pub max_abs_pressure: Pressure,
    pub max_red_pressure: f64,
    pub cooling: Utility,
}

impl From<OrganicRankineCycleJSON> for OrganicRankineCycle {
    fn from(orc: OrganicRankineCycleJSON) -> Self {
        Self {
            heat_source: orc.heat_source,
            isentropic_turbine_efficiency: orc.isentropic_turbine_efficiency,
            isentropic_pump_efficiency: orc.isentropic_pump_efficiency,
            min_abs_pressure: orc.min_abs_pressure * BAR,
            min_red_pressure: orc.min_red_pressure,
            max_abs_pressure: orc.max_abs_pressure * BAR,
            max_red_pressure: orc.max_red_pressure,
            cooling: orc.cooling,
        }
    }
}

impl From<OrganicRankineCycle> for OrganicRankineCycleJSON {
    fn from(orc: OrganicRankineCycle) -> Self {
        Self {
            heat_source: orc.heat_source,
            isentropic_turbine_efficiency: orc.isentropic_turbine_efficiency,
            isentropic_pump_efficiency: orc.isentropic_pump_efficiency,
            min_abs_pressure: orc.min_abs_pressure.convert_into(BAR),
            min_red_pressure: orc.min_red_pressure,
            max_abs_pressure: orc.max_abs_pressure.convert_into(BAR),
            max_red_pressure: orc.max_red_pressure,
            cooling: orc.cooling,
        }
    }
}

#[allow(clippy::too_many_arguments)]
impl OrganicRankineCycle {
    pub fn new(
        heat_source: Utility,
        isentropic_turbine_efficiency: f64,
        isentropic_pump_efficiency: f64,
        min_abs_pressure: Pressure,
        min_red_pressure: f64,
        max_abs_pressure: Pressure,
        max_red_pressure: f64,
        cooling: Utility,
    ) -> Self {
        Self {
            heat_source,
            isentropic_turbine_efficiency,
            isentropic_pump_efficiency,
            min_abs_pressure,
            min_red_pressure,
            max_abs_pressure,
            max_red_pressure,
            cooling,
        }
    }

    pub fn from_json<P: AsRef<Path>>(file: P) -> Result<Self, ParameterError> {
        let reader = BufReader::new(File::open(file)?);
        Ok(serde_json::from_reader(reader)?)
    }
}

impl ProcessModel for OrganicRankineCycle {
    fn variables(&self) -> ContinuousVariables {
        vec![
            Variable::continuous(Some(0.0), Some(2.0)),
            Variable::continuous(
                Some(self.min_red_pressure.ln()),
                Some(self.max_red_pressure.ln()),
            ),
            Variable::continuous(
                Some(self.min_red_pressure.ln()),
                Some(self.max_red_pressure.ln()),
            ),
            Variable::continuous(Some(0.0), Some(2.0)),
        ]
        .into()
    }

    fn equality_constraints(&self) -> usize {
        0
    }

    fn constraints(&self) -> Vec<[Option<f64>; 2]> {
        vec![
            // Pinch constraints evaporator
            [Some(0.0), None],
            [Some(0.0), None],
            [Some(0.0), None],
            [Some(0.0), None],
            // Pinch constraints condenser
            [Some(0.0), None],
            [Some(0.0), None],
            [Some(0.0), None],
            // Absolute pressure (bar)
            [
                Some(self.min_abs_pressure.convert_into(BAR)),
                Some(self.max_abs_pressure.convert_into(BAR)),
            ],
            [
                Some(self.min_abs_pressure.convert_into(BAR)),
                Some(self.max_abs_pressure.convert_into(BAR)),
            ],
            [Some(1.0), None],
        ]
    }

    fn solve<E: Residual + IdealGas>(
        &self,
        eos: &Arc<E>,
        x: &[f64],
    ) -> EosResult<(f64, Vec<f64>, Vec<f64>)> {
        // unpack variables
        let mwf = x[0] * 50.0 * KILOGRAM / SECOND;
        let p_cond_red = x[1].exp();
        let p_evap_red = x[2].exp();
        let dt_sh = x[3] * 50.0 * KELVIN;

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

        // Calculate evaporator
        let evaporator = process.evaporator(pump.out(), &isobar_evap, Some(dt_sh))?;
        process.add_utility(&evaporator, self.heat_source);

        // Calculate turbine
        let turbine = process.turbine(
            evaporator.out(),
            &isobar_cond,
            self.isentropic_turbine_efficiency,
        )?;

        // Calculate condenser
        let condenser = process.total_condenser(turbine.out(), &isobar_cond, None)?;
        process.add_utility(&condenser, self.cooling);

        // Target
        let target = process.net_power().unwrap().convert_into(MEGA * WATT);

        // Pinch constraints
        let mut constraints = process.pinch_constraints();

        // Absolute pressure constraints
        constraints.push(p_cond.convert_into(BAR));
        constraints.push(p_evap.convert_into(BAR));
        constraints.push(p_evap.convert_into(p_cond));

        Ok((target, vec![], constraints))
    }
}
