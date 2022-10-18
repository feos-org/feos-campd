//! Simple unit operations and process models.
use feos_core::{EosResult, EquationOfState, MolarWeight, StateBuilder};
use ndarray::arr1;
use petgraph::prelude::*;
use quantity::si::{SIArray1, SINumber, SIUnit, CELSIUS, JOULE, KELVIN, KILO, KILOGRAM, WATT};
use serde::{Deserialize, Serialize};
use std::ops::Index;
use std::rc::Rc;

mod equipment;
mod orc;
mod process_state;
pub use equipment::Equipment;
pub use orc::OrganicRankineCycle;
pub use process_state::{Isobar, ProcessState};

/// Generic process model to be used in an [OptimizationProblem].
pub trait ProcessModel {
    /// For each continuous process variable return the lower and upper bounds.
    fn variables(&self) -> Vec<[Option<f64>; 2]>;

    /// Return the number of binary process variables.
    fn binary_variables(&self) -> usize {
        0
    }

    /// For each constraint return the lower, upper and equality bound.
    fn constraints(&self) -> Vec<[Option<f64>; 3]>;

    /// Solve the process model and return the [Process], the target and the constraint values.
    fn solve<E: EquationOfState + MolarWeight<SIUnit>>(
        &self,
        eos: &Rc<E>,
        x: &[f64],
    ) -> EosResult<(Process<E>, f64, Vec<f64>)>;
}

/// The type used to index a [Process].
pub type StatePoint = NodeIndex;

/// Representation of a simple process as a graph structure.
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
    /// Create a new empty [Process].
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a [ProcessState] to the [Process].
    pub fn add_state(&mut self, state: ProcessState<E>) -> StatePoint {
        self.graph.add_node(state.into())
    }

    /// Connect a given [StatePoint] vie a [ProcessStep] with a new [ProcessState].
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

    /// Add a [ProcessStep] between two existing [StatePoint]s.
    pub fn add_connection(&mut self, feed: StatePoint, outlet: StatePoint, step: ProcessStep) {
        self.graph.add_edge(feed, outlet, step);
    }

    /// Return the utility temperature at a given [StatePoint].
    pub fn utility_temperature(&self, index: StatePoint) -> Option<SINumber> {
        self.graph[index].utility_temperature.map(|[t, _]| t)
    }
}

impl<E> Index<StatePoint> for Process<E> {
    type Output = ProcessState<E>;
    fn index(&self, index: StatePoint) -> &ProcessState<E> {
        &self.graph[index].state
    }
}

impl<E: EquationOfState + MolarWeight<SIUnit>> Process<E> {
    /// Add a [Utility] to a given [Equipment].
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

    /// Return the net power of the entire process.
    pub fn net_power(&self) -> Option<SINumber> {
        self.graph
            .edge_references()
            .map(|e| e.weight().power([&self[e.source()], &self[e.target()]]))
            .collect::<Option<Vec<_>>>()
            .map(|x| x.into_iter().reduce(|a, b| a + b).unwrap())
    }

    /// Return the values of all pinch constraints in the process.
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

    /// Return the necessary data to plot all process steps.
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

/// An edge in the [Process] graph. Contains a [ProcessState] the corresponding utility temperature.
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

/// Different ways of specifying a [Utility] stream.
#[derive(Clone, Copy, Serialize, Deserialize, Debug)]
#[serde(into = "UtilitySpecificationJSON")]
#[serde(from = "UtilitySpecificationJSON")]
pub enum UtilitySpecification {
    /// Constant temperature <-> infinite heat capacity rate.
    ConstantTemperature,
    /// Fixed heat capacity rate.
    HeatCapacityRate(SINumber),
    /// Heat capacity rate inferred from outlet temperature.
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
struct UtilityJSON {
    #[serde(rename = "temperature [°C]")]
    temperature: f64,
    specification: UtilitySpecification,
    #[serde(rename = "min_approach_temperature [K]")]
    min_approach_temperature: f64,
}

/// Data specifying a single utility stream.
#[derive(Clone, Copy, Serialize, Deserialize, Debug)]
#[serde(into = "UtilityJSON")]
#[serde(from = "UtilityJSON")]
pub struct Utility {
    pub temperature: SINumber,
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

/// Possible state changes inbetween [ProcessState]s.
pub enum ProcessStep {
    PhaseChange,
    Isobaric,
    Polytropic,
}

impl ProcessStep {
    /// Return the power in-/output in the step.
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

    /// Return data required for drawing the step in a Ts plot.
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

/// Collection of data required to create a Ts plot of a process.
#[derive(Clone, Serialize, Deserialize)]
#[serde(from = "ProcessPlotJSON", into = "ProcessPlotJSON")]
pub struct ProcessPlot {
    pub entropy: SIArray1,
    pub temperature: SIArray1,
    pub utility_temperature: Option<SIArray1>,
}

#[derive(Serialize, Deserialize)]
struct ProcessPlotJSON {
    #[serde(rename = "entropy [kJ/kg/K]")]
    entropy: Vec<f64>,
    #[serde(rename = "temperature [°C]")]
    temperature: Vec<f64>,
    #[serde(rename = "utility temperature [°C]")]
    utility_temperature: Option<Vec<f64>>,
}

impl From<ProcessPlot> for ProcessPlotJSON {
    fn from(plot: ProcessPlot) -> Self {
        Self {
            entropy: plot
                .entropy
                .to_reduced(KILO * JOULE / KILOGRAM / KELVIN)
                .unwrap()
                .to_vec(),
            temperature: (plot.temperature / CELSIUS).to_vec(),
            utility_temperature: plot
                .utility_temperature
                .map(|utility_temperature| (utility_temperature / CELSIUS).to_vec()),
        }
    }
}

impl From<ProcessPlotJSON> for ProcessPlot {
    fn from(plot: ProcessPlotJSON) -> Self {
        Self {
            entropy: arr1(&plot.entropy) * KILO * JOULE / KILOGRAM / KELVIN,
            temperature: arr1(&plot.temperature) * CELSIUS,
            utility_temperature: plot
                .utility_temperature
                .map(|utility_temperature| arr1(&utility_temperature) * CELSIUS),
        }
    }
}
