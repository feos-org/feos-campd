//! Simple unit operations and process models.
use crate::variables::ProcessVariables;
#[cfg(feature = "knitro_rs")]
use crate::OptimizationMode;
use feos::core::si::*;
use feos::core::{EosResult, IdealGas, Residual, StateBuilder};
#[cfg(feature = "knitro_rs")]
use knitro_rs::{Knitro, KnitroError};
use ndarray::{arr1, Array1};
use petgraph::prelude::*;
use serde::{Deserialize, Serialize};
use std::ops::Index;
use std::sync::Arc;

mod equipment;
mod orc;
mod process_state;
pub use equipment::Equipment;
pub use orc::OrganicRankineCycle;
pub use process_state::{Isobar, ProcessState};

/// Generic process model to be used in an [OptimizationProblem](super::OptimizationProblem).
pub trait ProcessModel<E> {
    /// For each continuous process variable return the lower and upper bounds.
    fn variables(&self) -> ProcessVariables;

    /// Return the number of equality constraints (h(x)=0)
    fn equality_constraints(&self) -> usize;

    /// Return the number of inequality constraints (g(x)>=0)
    fn inequality_constraints(&self) -> usize;

    /// Solve the process model and return the target, and the values of equality and inequality constraints.
    fn solve(&self, eos: &Arc<E>, x: &[f64]) -> EosResult<(f64, Vec<f64>, Vec<f64>)>;

    fn _solve(&self, eos: &Arc<E>, x: &[f64]) -> EosResult<(f64, Vec<f64>, Vec<f64>)>
    where
        E: Residual + IdealGas,
    {
        let (x, u) = x.split_at(self.variables().len());
        match *u {
            [u] => self.solve_infeasible(eos, x, u),
            [] => self.solve(eos, x),
            _ => unreachable!(),
        }
    }

    fn solve_infeasible(
        &self,
        eos: &Arc<E>,
        x: &[f64],
        u: f64,
    ) -> EosResult<(f64, Vec<f64>, Vec<f64>)> {
        let (_, h, mut g) = self.solve(eos, x)?;
        g.iter_mut().for_each(|g| *g += u);
        Ok((u, h, g))
    }

    #[cfg(feature = "knitro_rs")]
    fn setup_knitro(
        &self,
        kc: &Knitro,
        x0: Option<&[f64]>,
        mode: OptimizationMode,
    ) -> Result<[Vec<i32>; 3], KnitroError> {
        // declare continuous variables
        let mut index_vars = self.variables().setup_knitro(kc, x0, mode)?;

        // Declare feasibility bound
        if let OptimizationMode::Feasibility = mode {
            index_vars.push(kc.add_var()?);
        }

        // add equality constraints
        let index_eq_cons = kc.add_cons(self.equality_constraints())?;
        for &i in &index_eq_cons {
            kc.set_con_eqbnd(i, 0.0)?;
        }

        // add inequality constraints
        let index_ineq_cons = kc.add_cons(self.inequality_constraints())?;
        for &i in &index_ineq_cons {
            kc.set_con_lobnd(i, 0.0)?;
        }

        Ok([index_vars, index_eq_cons, index_ineq_cons])
    }
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
    pub fn utility_temperature(&self, index: StatePoint) -> Option<Temperature<f64>> {
        self.graph[index].utility_temperature.map(|[t, _]| t)
    }
}

impl<E> Index<StatePoint> for Process<E> {
    type Output = ProcessState<E>;
    fn index(&self, index: StatePoint) -> &ProcessState<E> {
        &self.graph[index].state
    }
}

impl<E: Residual + IdealGas> Process<E> {
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
    pub fn net_power(&self) -> Option<Power<f64>> {
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
                    let s = dt.signum() * KELVIN;
                    (t - p.state.temperature() - dt).convert_into(dt + s)
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
    utility_temperature: Option<[Temperature<f64>; 2]>,
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
    HeatCapacityRate(HeatCapacityRate<f64>),
    /// Heat capacity rate inferred from outlet temperature.
    OutletTemperature(Temperature<f64>),
}

impl From<UtilitySpecification> for UtilitySpecificationJSON {
    fn from(utility: UtilitySpecification) -> Self {
        match utility {
            UtilitySpecification::ConstantTemperature => {
                UtilitySpecificationJSON::ConstantTemperature
            }
            UtilitySpecification::HeatCapacityRate(c) => {
                UtilitySpecificationJSON::HeatCapacityRate(c.convert_into(KILO * WATT / KELVIN))
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
    pub temperature: Temperature<f64>,
    pub specification: UtilitySpecification,
    pub min_approach_temperature: Temperature<f64>,
}

impl From<Utility> for UtilityJSON {
    fn from(utility: Utility) -> Self {
        Self {
            temperature: utility.temperature / CELSIUS,
            specification: utility.specification,
            min_approach_temperature: utility.min_approach_temperature.convert_into(KELVIN),
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
        temperature: Temperature<f64>,
        specification: UtilitySpecification,
        min_approach_temperature: Temperature<f64>,
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
    pub fn power<E: Residual + IdealGas>(
        &self,
        states: [&ProcessState<E>; 2],
    ) -> Option<Power<f64>> {
        match self {
            Self::Polytropic => states[0]
                .mass_flow_rate()
                .map(|m| m * (states[1].specific_enthalpy() - states[0].specific_enthalpy())),
            _ => Some(0.0 * WATT),
        }
    }

    /// Return data required for drawing the step in a Ts plot.
    pub fn plot_ts<E: Residual + IdealGas>(
        &self,
        states: [&ProcessPoint<E>; 2],
    ) -> EosResult<ProcessPlot> {
        Ok(match self {
            Self::PhaseChange => ProcessPlot {
                temperature: Temperature::from_shape_fn(2, |i| states[i].state.temperature()),
                entropy: SpecificEntropy::from_shape_fn(2, |i| states[i].state.specific_entropy()),
                utility_temperature: if let [Some([t1, _]), Some([t2, _])] =
                    [states[0].utility_temperature, states[1].utility_temperature]
                {
                    Some(Temperature::from_vec(vec![t1, t2]))
                } else {
                    None
                },
            },
            Self::Isobaric => {
                let temperature = Temperature::linspace(
                    states[0].state.temperature(),
                    states[1].state.temperature(),
                    50,
                );
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
                let entropy = SpecificEntropy::from_shape_fn(state_vec.len(), |i| {
                    state_vec[i].specific_entropy()
                });
                let utility_temperature = if let [Some([t1, _]), Some([t2, _])] =
                    [states[0].utility_temperature, states[1].utility_temperature]
                {
                    let enthalpy = SpecificEnergy::from_shape_fn(state_vec.len(), |i| {
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

                let k = (s1 - s0) / t1.convert_into(t0).ln();
                let temperature = Temperature::linspace(t0, t1, 50);
                let entropy = s0 + (&temperature / t0).into_value().mapv(f64::ln) * k;
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
    pub entropy: SpecificEntropy<Array1<f64>>,
    pub temperature: Temperature<Array1<f64>>,
    pub utility_temperature: Option<Temperature<Array1<f64>>>,
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
            entropy: (plot.entropy / (KILO * JOULE / KILOGRAM / KELVIN))
                .into_value()
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
