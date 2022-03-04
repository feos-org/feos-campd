use feos_core::{EosResult, EquationOfState, MolarWeight, StateBuilder};
use petgraph::prelude::*;
use quantity::si::{SIArray1, SINumber, SIUnit, WATT};
use serde::{Deserialize, Serialize};
use std::ops::Index;
use std::rc::Rc;

mod equipment;
mod heat_pump;
mod orc;
mod process_state;
pub use equipment::Equipment;
pub use heat_pump::HeatPump;
pub use orc::OrganicRankineCycle;
pub use process_state::ProcessState;

#[derive(Clone, Serialize, Deserialize)]
pub enum ProcessModel {
    OrganicRankineCycle(OrganicRankineCycle),
}

impl ProcessModel {
    pub fn variables(&self) -> Vec<[Option<f64>; 2]> {
        match self {
            Self::OrganicRankineCycle(process) => process.variables(),
        }
    }

    pub fn constraints(&self) -> Vec<[Option<f64>; 3]> {
        match self {
            Self::OrganicRankineCycle(process) => process.constraints(),
        }
    }

    pub fn solve<E: EquationOfState + MolarWeight<SIUnit>>(
        &self,
        eos: &Rc<E>,
        x: &[f64],
    ) -> EosResult<(Process<E>, f64, Vec<f64>)> {
        match self {
            Self::OrganicRankineCycle(process) => process.solve(eos, x),
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

    pub fn initialize(&mut self, state: ProcessState<E>) -> StatePoint {
        self.graph.add_node(state.into())
    }
}

impl<E> Index<StatePoint> for Process<E> {
    type Output = ProcessState<E>;
    fn index(&self, index: StatePoint) -> &ProcessState<E> {
        &self.graph[index].state
    }
}

impl<E: EquationOfState + MolarWeight<SIUnit>> Process<E> {
    pub fn add_step(
        &mut self,
        feed: StatePoint,
        state: ProcessState<E>,
        step: ProcessStep,
    ) -> (StatePoint, EdgeIndex) {
        let out = self.graph.add_node(state.into());
        (out, self.graph.add_edge(feed, out, step))
    }

    pub fn add_utility(
        &mut self,
        equipment: &Equipment,
        temperature: SINumber,
        mass_flow_rate: SINumber,
        heat_capacity: SINumber,
        min_approach_temperature: SINumber,
    ) {
        let states = &equipment.states;
        let out = *states.last().unwrap();
        let mut h2 = self[out].specific_enthalpy();
        let mut t2 = temperature;
        self.graph[out].utility_temperature = Some([t2, min_approach_temperature]);
        for (&s1, &s2) in states.iter().rev().skip(1).zip(states.iter().rev()) {
            let index = self.graph.find_edge(s1, s2).unwrap();
            let h1 = self[s1].specific_enthalpy();
            let t1 = t2
                - (h2 - h1) / (mass_flow_rate * heat_capacity)
                    * self.graph[index].mass_flow_rate.unwrap();
            self.graph[s1].utility_temperature = Some([t1, min_approach_temperature]);
            self.graph[index].utility = Some(heat_capacity * mass_flow_rate);
            h2 = h1;
            t2 = t1;
        }
    }

    pub fn add_utility_heat_sink(
        &mut self,
        equipment: &Equipment,
        temperature_in: SINumber,
        temperature_out: SINumber,
        power: SINumber,
        min_approach_temperature: SINumber,
    ) {
        let states = &equipment.states;
        let out = *states.last().unwrap();
        let inlet = *states.first().unwrap();
        let mut h2 = self[out].specific_enthalpy();
        let mut t2 = temperature_in;
        // Calculate necessary mass flow rate
        let mwf = power / (self[inlet].specific_enthalpy() - self[out].specific_enthalpy());

        let mcp_utility = power / (temperature_out - temperature_in);
        self.graph[out].utility_temperature = Some([t2, min_approach_temperature]);
        for (&s1, &s2) in states.iter().rev().skip(1).zip(states.iter().rev()) {
            let index = self.graph.find_edge(s1, s2).unwrap();
            let h1 = self[s1].specific_enthalpy();
            let mass_flow_rate_wf = self.graph[index].mass_flow_rate.unwrap_or(mwf);
            let t1 = t2 - (h2 - h1) / (mcp_utility) * mass_flow_rate_wf;
            self.graph[s1].utility_temperature = Some([t1, min_approach_temperature]);
            self.graph[index].utility = Some(mcp_utility);
            h2 = h1;
            t2 = t1;
        }
    }

    pub fn net_power(&self) -> SINumber {
        self.graph
            .edge_references()
            .map(|e| e.weight().power([&self[e.source()], &self[e.target()]]))
            .reduce(|a, b| a + b)
            .unwrap()
    }

    pub fn pinch_constraints(&self) -> Vec<f64> {
        self.graph
            .node_weights()
            .filter_map(|p| {
                p.utility_temperature
                    .map(|t| (t[0] - p.state.temperature()).to_reduced(t[1]).unwrap())
            })
            .collect()
    }

    pub fn plot(&self) -> EosResult<Vec<ProcessPlot>> {
        Ok(self
            .graph
            .edge_references()
            .map(|e| {
                e.weight()
                    .plot_ts([&self.graph[e.source()], &self.graph[e.target()]])
                    .unwrap()
            })
            .collect())
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

enum StateChange {
    PhaseChange,
    Isobaric,
    Polytropic,
}

pub struct ProcessStep {
    state_change: StateChange,
    mass_flow_rate: Option<SINumber>,
    utility: Option<SINumber>,
}

impl ProcessStep {
    pub fn phase_change(mass_flow_rate: Option<SINumber>) -> Self {
        Self {
            state_change: StateChange::PhaseChange,
            mass_flow_rate,
            utility: None,
        }
    }

    pub fn isobaric(mass_flow_rate: Option<SINumber>) -> Self {
        Self {
            state_change: StateChange::Isobaric,
            mass_flow_rate,
            utility: None,
        }
    }

    pub fn polytropic(mass_flow_rate: Option<SINumber>) -> Self {
        Self {
            state_change: StateChange::Polytropic,
            mass_flow_rate,
            utility: None,
        }
    }

    pub fn power<E: EquationOfState + MolarWeight<SIUnit>>(
        &self,
        states: [&ProcessState<E>; 2],
    ) -> SINumber {
        match self.state_change {
            StateChange::Polytropic => {
                (states[1].specific_enthalpy() - states[0].specific_enthalpy())
                    * self.mass_flow_rate.unwrap()
            }
            _ => 0.0 * WATT,
        }
    }

    pub fn plot_ts<E: EquationOfState + MolarWeight<SIUnit>>(
        &self,
        states: [&ProcessPoint<E>; 2],
    ) -> EosResult<ProcessPlot> {
        Ok(match self.state_change {
            StateChange::PhaseChange => ProcessPlot {
                temperature: SIArray1::from_shape_fn(2, |i| states[i].state.temperature()),
                entropy: SIArray1::from_shape_fn(2, |i| states[i].state.specific_entropy()),
                utility_temperature: self.utility.map(|_| {
                    SIArray1::from_shape_fn(2, |i| states[i].utility_temperature.unwrap()[0])
                }),
            },
            StateChange::Isobaric => {
                let temperature = SIArray1::linspace(
                    states[0].state.temperature(),
                    states[1].state.temperature(),
                    50,
                )?;
                let pressure = states[0].state.pressure();
                let mut state_vec = vec![states[0].state.clone()];
                state_vec.extend((1..49).map(|i| {
                    ProcessState::SinglePhase(Box::new(
                        StateBuilder::new(states[0].state.eos())
                            .temperature(temperature.get(i))
                            .pressure(pressure)
                            .build()
                            .unwrap(),
                    ))
                }));
                state_vec.push(states[1].state.clone());
                let entropy =
                    SIArray1::from_shape_fn(state_vec.len(), |i| state_vec[i].specific_entropy());
                let utility_temperature = self.utility.map(|u| {
                    let enthalpy = SIArray1::from_shape_fn(state_vec.len(), |i| {
                        state_vec[i].specific_enthalpy()
                    });
                    let h_out = enthalpy.get(49);
                    states[1].utility_temperature.unwrap()[0]
                        + (enthalpy - h_out) / u * self.mass_flow_rate.unwrap()
                });
                ProcessPlot {
                    temperature,
                    entropy,
                    utility_temperature,
                }
            }
            StateChange::Polytropic => {
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
