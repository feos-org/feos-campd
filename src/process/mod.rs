use feos_core::{EosResult, EquationOfState, MolarWeight, StateBuilder};
use petgraph::prelude::*;
use quantity::si::{SIArray1, SINumber, SIUnit, CELSIUS, JOULE, KELVIN, KILO, KILOGRAM, WATT};
use serde::{Deserialize, Serialize};
use std::ops::Index;
use std::rc::Rc;

mod equipment;
// mod heat_pump;
mod orc;
mod orc_superstructure;
mod process_state;
pub use equipment::Equipment;
// pub use heat_pump::HeatPump;
pub use orc::OrganicRankineCycle;
pub use orc_superstructure::OrganicRankineCycleSuperStructure;
pub use process_state::{Isobar, ProcessState};

#[derive(Clone, Serialize, Deserialize)]
pub enum ProcessModel {
    OrganicRankineCycle(OrganicRankineCycle),
    OrganicRankineCycleSuperStructure(OrganicRankineCycleSuperStructure),
}

impl ProcessModel {
    pub fn variables(&self) -> Vec<[Option<f64>; 2]> {
        match self {
            Self::OrganicRankineCycle(process) => process.variables(),
            Self::OrganicRankineCycleSuperStructure(process) => process.variables(),
        }
    }

    pub fn binary_variables(&self) -> usize {
        match self {
            Self::OrganicRankineCycle(process) => process.binary_variables(),
            Self::OrganicRankineCycleSuperStructure(process) => process.binary_variables(),
        }
    }

    pub fn constraints(&self) -> Vec<[Option<f64>; 3]> {
        match self {
            Self::OrganicRankineCycle(process) => process.constraints(),
            Self::OrganicRankineCycleSuperStructure(process) => process.constraints(),
        }
    }

    pub fn solve<E: EquationOfState + MolarWeight<SIUnit>>(
        &self,
        eos: &Rc<E>,
        x: &[f64],
    ) -> EosResult<(Process<E>, f64, Vec<f64>)> {
        match self {
            Self::OrganicRankineCycle(process) => process.solve(eos, x),
            Self::OrganicRankineCycleSuperStructure(process) => process.solve(eos, x),
        }
    }
}

pub type StatePoint = NodeIndex;

pub struct Process<E> {
    graph: Graph<ProcessPoint<E>, ProcessStep>,
}

impl<E> Default for Process<E> {
    fn default() -> Self {
        Self {
            graph: Graph::new(),
        }
    }
}

impl<E> Process<E> {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_state(&mut self, state: ProcessState<E>) -> StatePoint {
        self.graph.add_node(state.into())
    }

    pub fn add_step(
        &mut self,
        feed: StatePoint,
        state: ProcessState<E>,
        step: ProcessStep,
    ) -> StatePoint {
        let out = self.graph.add_node(state.into());
        self.graph.add_edge(feed, out, step);
        out
    }

    pub fn add_connection(&mut self, feed: StatePoint, outlet: StatePoint, step: ProcessStep) {
        self.graph.add_edge(feed, outlet, step);
    }
}

impl<E> Index<StatePoint> for Process<E> {
    type Output = ProcessState<E>;
    fn index(&self, index: StatePoint) -> &ProcessState<E> {
        &self.graph[index].state
    }
}

impl<E: EquationOfState + MolarWeight<SIUnit>> Process<E> {
    pub fn add_utility(&mut self, equipment: &Equipment, utility: Utility) {
        let states = &equipment.states;
        let inlet = *states.first().unwrap();
        let outlet = *states.last().unwrap();
        let mut h2 = self[outlet].specific_enthalpy();
        let x = match utility.specification {
            UtilitySpecification::ConstantTemperature => 0.0 * KELVIN / (JOULE / KILOGRAM),
            UtilitySpecification::HeatCapacityRate(c) => self[inlet].mass_flow_rate().unwrap() / c,
            UtilitySpecification::OutletTemperature(t_util_out) => {
                (t_util_out - utility.temperature)
                    / (self[inlet].specific_enthalpy() - self[outlet].specific_enthalpy())
            }
        };
        let mut t2 = utility.temperature;
        self.graph[outlet].utility_temperature = Some([t2, utility.min_approach_temperature]);
        for &s1 in states.iter().rev().skip(1) {
            let h1 = self[s1].specific_enthalpy();
            let t1 = t2 - (h2 - h1) * x;
            self.graph[s1].utility_temperature = Some([t1, utility.min_approach_temperature]);
            h2 = h1;
            t2 = t1;
        }
    }

    // pub fn add_utility_heat_sink(
    //     &mut self,
    //     equipment: &Equipment,
    //     temperature_in: SINumber,
    //     temperature_out: SINumber,
    //     power: SINumber,
    //     min_approach_temperature: SINumber,
    // ) {
    //     let states = &equipment.states;
    //     let out = *states.last().unwrap();
    //     let inlet = *states.first().unwrap();
    //     let mut h2 = self[out].specific_enthalpy();
    //     let mut t2 = temperature_in;
    //     // Calculate necessary mass flow rate
    //     let mwf = power / (self[inlet].specific_enthalpy() - self[out].specific_enthalpy());

    //     let mcp_utility = power / (temperature_out - temperature_in);
    //     self.graph[out].utility_temperature = Some([t2, min_approach_temperature]);
    //     for (&s1, &s2) in states.iter().rev().skip(1).zip(states.iter().rev()) {
    //         let index = self.graph.find_edge(s1, s2).unwrap();
    //         let h1 = self[s1].specific_enthalpy();
    //         let mass_flow_rate_wf = self.graph[index].mass_flow_rate.unwrap_or(mwf);
    //         let t1 = t2 - (h2 - h1) / (mcp_utility) * mass_flow_rate_wf;
    //         self.graph[s1].utility_temperature = Some([t1, min_approach_temperature]);
    //         // self.graph[index].utility = Some(mcp_utility);
    //         h2 = h1;
    //         t2 = t1;
    //     }
    // }

    pub fn net_power(&self) -> Option<SINumber> {
        self.graph
            .edge_references()
            .map(|e| e.weight().power([&self[e.source()], &self[e.target()]]))
            .collect::<Option<Vec<_>>>()
            .map(|x| x.into_iter().reduce(|a, b| a + b).unwrap())
    }

    pub fn pinch_constraints(&self) -> Vec<f64> {
        self.graph
            .node_weights()
            .filter_map(|p| {
                p.utility_temperature.map(|[t, dt]| {
                    let s = dt.to_reduced(KELVIN).unwrap().signum() * KELVIN;
                    (t - p.state.temperature() - dt).to_reduced(dt + s).unwrap()
                })
            })
            .collect()
    }

    pub fn plot(&self) -> EosResult<Vec<ProcessPlot>> {
        self.graph
            .edge_references()
            .map(|e| {
                e.weight()
                    .plot_ts([&self.graph[e.source()], &self.graph[e.target()]])
            })
            .collect()
    }
}

pub struct ProcessPoint<E> {
    state: ProcessState<E>,
    utility_temperature: Option<[SINumber; 2]>,
}

impl<E> From<ProcessState<E>> for ProcessPoint<E> {
    fn from(state: ProcessState<E>) -> Self {
        Self {
            state,
            utility_temperature: None,
        }
    }
}

#[derive(Serialize, Deserialize)]
enum UtilitySpecificationJSON {
    #[serde(rename = "constant")]
    ConstantTemperature,
    #[serde(rename = "heat_capacity_rate [kW/K]")]
    HeatCapacityRate(f64),
    #[serde(rename = "outlet_temperature [°C]")]
    OutletTemperature(f64),
}

#[derive(Clone, Copy, Serialize, Deserialize, Debug)]
#[serde(into = "UtilitySpecificationJSON")]
#[serde(from = "UtilitySpecificationJSON")]
pub enum UtilitySpecification {
    ConstantTemperature,
    HeatCapacityRate(SINumber),
    OutletTemperature(SINumber),
}

impl From<UtilitySpecification> for UtilitySpecificationJSON {
    fn from(utility: UtilitySpecification) -> Self {
        match utility {
            UtilitySpecification::ConstantTemperature => {
                UtilitySpecificationJSON::ConstantTemperature
            }
            UtilitySpecification::HeatCapacityRate(c) => {
                UtilitySpecificationJSON::HeatCapacityRate(
                    c.to_reduced(KILO * WATT / KELVIN).unwrap(),
                )
            }
            UtilitySpecification::OutletTemperature(t) => {
                UtilitySpecificationJSON::OutletTemperature(t / CELSIUS)
            }
        }
    }
}

impl From<UtilitySpecificationJSON> for UtilitySpecification {
    fn from(utility: UtilitySpecificationJSON) -> Self {
        match utility {
            UtilitySpecificationJSON::ConstantTemperature => {
                UtilitySpecification::ConstantTemperature
            }
            UtilitySpecificationJSON::HeatCapacityRate(c) => {
                UtilitySpecification::HeatCapacityRate(c * KILO * WATT / KELVIN)
            }
            UtilitySpecificationJSON::OutletTemperature(t) => {
                UtilitySpecification::OutletTemperature(t * CELSIUS)
            }
        }
    }
}

#[derive(Clone, Copy, Serialize, Deserialize, Debug)]
pub struct UtilityJSON {
    #[serde(rename = "temperature [°C]")]
    temperature: f64,
    specification: UtilitySpecification,
    #[serde(rename = "min_approach_temperature [K]")]
    min_approach_temperature: f64,
}

#[derive(Clone, Copy, Serialize, Deserialize, Debug)]
#[serde(into = "UtilityJSON")]
#[serde(from = "UtilityJSON")]
pub struct Utility {
    temperature: SINumber,
    specification: UtilitySpecification,
    min_approach_temperature: SINumber,
}

impl From<Utility> for UtilityJSON {
    fn from(utility: Utility) -> Self {
        Self {
            temperature: utility.temperature / CELSIUS,
            specification: utility.specification,
            min_approach_temperature: utility.min_approach_temperature.to_reduced(KELVIN).unwrap(),
        }
    }
}

impl From<UtilityJSON> for Utility {
    fn from(utility: UtilityJSON) -> Self {
        Self {
            temperature: utility.temperature * CELSIUS,
            specification: utility.specification,
            min_approach_temperature: utility.min_approach_temperature * KELVIN,
        }
    }
}

impl Utility {
    pub fn new(
        temperature: SINumber,
        specification: UtilitySpecification,
        min_approach_temperature: SINumber,
    ) -> Self {
        Self {
            temperature,
            specification,
            min_approach_temperature,
        }
    }
}

pub enum ProcessStep {
    PhaseChange,
    Isobaric,
    Polytropic,
}

impl ProcessStep {
    pub fn power<E: EquationOfState + MolarWeight<SIUnit>>(
        &self,
        states: [&ProcessState<E>; 2],
    ) -> Option<SINumber> {
        match self {
            Self::Polytropic => states[0]
                .mass_flow_rate()
                .map(|m| m * (states[1].specific_enthalpy() - states[0].specific_enthalpy())),
            _ => Some(0.0 * WATT),
        }
    }

    pub fn plot_ts<E: EquationOfState + MolarWeight<SIUnit>>(
        &self,
        states: [&ProcessPoint<E>; 2],
    ) -> EosResult<ProcessPlot> {
        Ok(match self {
            Self::PhaseChange => ProcessPlot {
                temperature: SIArray1::from_shape_fn(2, |i| states[i].state.temperature()),
                entropy: SIArray1::from_shape_fn(2, |i| states[i].state.specific_entropy()),
                utility_temperature: if let [Some([t1, _]), Some([t2, _])] =
                    [states[0].utility_temperature, states[1].utility_temperature]
                {
                    Some(SIArray1::from_vec(vec![t1, t2]))
                } else {
                    None
                },
            },
            Self::Isobaric => {
                let temperature = SIArray1::linspace(
                    states[0].state.temperature(),
                    states[1].state.temperature(),
                    50,
                )?;
                let pressure = states[0].state.pressure();
                let mut state_vec = vec![states[0].state.clone()];
                state_vec.extend((1..49).map(|i| {
                    ProcessState::SinglePhase(
                        Box::new(
                            StateBuilder::new(states[0].state.eos())
                                .temperature(temperature.get(i))
                                .pressure(pressure)
                                .initial_density(states[0].state.density())
                                .build()
                                .unwrap(),
                        ),
                        None,
                    )
                }));
                state_vec.push(states[1].state.clone());
                let entropy =
                    SIArray1::from_shape_fn(state_vec.len(), |i| state_vec[i].specific_entropy());
                let utility_temperature = if let [Some([t1, _]), Some([t2, _])] =
                    [states[0].utility_temperature, states[1].utility_temperature]
                {
                    let enthalpy = SIArray1::from_shape_fn(state_vec.len(), |i| {
                        state_vec[i].specific_enthalpy()
                    });
                    let h_in = enthalpy.get(0);
                    let h_out = enthalpy.get(49);
                    Some(t2 + (enthalpy - h_out) / (h_in - h_out) * (t1 - t2))
                } else {
                    None
                };
                ProcessPlot {
                    temperature,
                    entropy,
                    utility_temperature,
                }
            }
            Self::Polytropic => {
                let t0 = states[0].state.temperature();
                let t1 = states[1].state.temperature();
                let s0 = states[0].state.specific_entropy();
                let s1 = states[1].state.specific_entropy();

                let k = (s1 - s0) / t1.to_reduced(t0)?.ln();
                let temperature = SIArray1::linspace(t0, t1, 50)?;
                let entropy = s0 + k * temperature.to_reduced(t0)?.mapv(f64::ln);
                ProcessPlot {
                    temperature,
                    entropy,
                    utility_temperature: None,
                }
            }
        })
    }
}

pub struct ProcessPlot {
    pub entropy: SIArray1,
    pub temperature: SIArray1,
    pub utility_temperature: Option<SIArray1>,
}
