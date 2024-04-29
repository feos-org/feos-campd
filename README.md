# FeOs-CAMPD

[![crate](https://img.shields.io/crates/v/feos-campd.svg)](https://crates.io/crates/feos-campd)
[![documentation](https://docs.rs/feos-campd/badge.svg)](https://docs.rs/feos-campd)

Computer-aided molecular and process design based on the [`FeOs`](https://github.com/feos-org/feos) framework.

The package provides some basic infrastructure to perform a computer-aided molecular and process design. It consists of
- NLP/MIQCP solver bindings (Artelys Knitro)
- Implementation of custom outer-approximation algorithms to solve the resulting MINLP
- Molecular representations (group counts (CoMT-CAMD) and molecule superstructures)
- Property models (PC-SAFT and (heterosegmented) gc-PC-SAFT)

and the surrounding framework that is used to run the optimization problems.

Currently, solving CAMPD problems is only possible with an installation of Artelys Knitro 12.

## Usage
### Molecular representation
```python
# define the molecular representation using a molecule superstructure
molecule = SuperMolecule.alkane(size)
# or a combination of superstructures
molecule = SuperMolecule.non_associating(size)

# or provide a list of molecules to choose from
molecule = CoMTCAMD.from_molecules(list_of_identifiers)
molecule = CoMTCAMD.from_json_molecules(list_of_identifiers)

# or provide an input file with group and structure definitions for CoMT-CAMD
molecule = CoMTCAMD.from_json(list_of_identifiers)
```
### Property model
```python
# The available property models are compatible with the parameter files in FeOs
property_model = GcPcSaftPropertyModel.from_json([molecule], "sauer2014_hetero.json", "joback1987.json")
property_model = PcSaftPropertyModel.from_json_molecules("gross2001.json", "poling2000.json")
property_model = PcSaftPropertyModel.from_json_groups("sauer2014_homo.json", "joback1987.json")
```
### Process model
```python
#To implement a process model in Python, define a class with the following methods:
class ORC:
    # For each process variable: the lower bound, the upper bound, and the initial value
    def variables(self):
        return [[lb_0, ub_0, init_0], [lb_1, ub_1, init_1], ...]
    
    # The number of equality constraints (h(x) = 0)
    def equality_constraints(self):
        return ...
    
    # The number of inequality constraints (g(x) >= 0)
    def inequality_constraints(self):
        return ...
        
    # For given equation of state and process variables x, return the target, and the values of
    # equality and inequality constraints
    def solve(self, eos, x):
        # eos - the equation of state object as used in FeOs, e.g.,
        state = State(eos, temperature=300*KELVIN, pressure=BAR)

        # x - the list of process degrees of freedom, e.g.,
        T_in, p_des, ... = x

        # The function is called with regular Python data types (floats), so the implementation
        # of the process model can be as flexible as desired and involve external function calls
        ...

        return target, equality_constraints, inequality_constraints

```
### Optimization problem
```python
# combine molecular representation, property model, and process model in an optimization problem
# for a pure component
problem = OptimizationProblem.pure(molecule, property_model, process)

# or a binary mixture
problem = OptimizationProblem.pure([molecule1, molecule2], property_model, process)

# use either of these algorithms
# (The boolean indicates whether the algorithm should update the lower bound of the outer approximation
# True can lead to local optima in non-convex problems)
algorithm = OuterApproximationAlgorithm.DuranGrossmann(True)
algorithm = OuterApproximationAlgorithm.DuranGrossmann(False)
algorithm = OuterApproximationAlgorithm.FletcherLeyffer

# to determine a ranking of the optimal molecules
problem.outer_approximation_ranking(y0, algorithm, num_molecules, "options_NLP.opt", "options_MILP.opt")
```

Molecular representations and property models can be combined according to:
||PC-SAFT|gc-PC-SAFT|
|-|-|-|
| CoMTCAMD | yes | no |
| SuperMolecule | yes | yes |