# FeOs-CAMPD

[![crate](https://img.shields.io/crates/v/feos-campd.svg)](https://crates.io/crates/feos-campd)
[![documentation](https://docs.rs/feos-campd/badge.svg)](https://docs.rs/feos-campd)

Computer-aided molecular and process design based on the [`FeOs`](https://github.com/feos-org/feos) framework.

The package provides infrastructure to perform a computer-aided molecular and process design. It consists of
- Implementation of a custom outer-approximation algorithms to solve the resulting MINLP
- Molecular representations (group counts (CoMT-CAMD) and molecule superstructures)
- Property models (PC-SAFT and (heterosegmented) gc-PC-SAFT)

The underlying NLP subproblems are solved using [IPOPT](https://coin-or.github.io/Ipopt/) via the [`ipopt-ad`](https://github.com/prehner/ipopt-ad) crate. For the solution of th MILP master problems, the solver is linked to the [`good_lp`](https://github.com/rust-or/good_lp) crate that offers a common interface to various open-source MILP solvers.

## Installation
### Rust
Just add the dependency to your `Cargo.toml`
```toml
feos-campd = "0.3"
```

## Usage
The following sections demonstrates the usage of the framework, for the API details, check out the [documentation](https://docs.rs/feos-campd).

### Molecular representation
```rust
// define the molecular representation using a molecule superstructure
let molecule = SuperMolecule::alkane(min_size, size);
// or a combination of superstructures
let molecule = SuperMolecule::non_associating(min_size, size);
// or use group and structure definitions for CoMT-CAMD
// (currently not further customizable)
let molecule = CoMTCAMD;
```
### Property model
```rust
// The property model for heterosegmented gc-PC-SAFT is not customizable. It does not include association.
let property_model = GcPcSaftPropertyModel;
// For homosegmented PC-SAFT, the group parameters can be passes as arguments
// and it is available as the full model
let property_model = PcSaftPropertyModel::full("rehner2023_homo.json", Some("rehner2023_homo_binary.json"));
// or non-associating (ignoring all associating groups)
let property_model = PcSaftPropertyModel::non_associating("rehner2023_homo.json", Some("rehner2023_homo_binary.json"));
```
### Process model
```rust
//To define a process model in Rust, implement the `ProcessModel` trait for your struct.
impl<E: TotalHelmholtzEnergy<N: ...>> ProcessModel<E, N_X: ..., N: ...> for YourModel {
    fn variables(&self) -> [ContinuousVariable; N_X] {
        [
            ContinuousVariable::new(..., ..., ...),
            ...
        ]
    }

    fn constraints(&self) -> Vec<GeneralConstraint> {
        vec![
            GeneralConstraint::Inequality(..., ...),
            GeneralConstraint::Equality(...),
            ...
        ]
    }

    fn evaluate<D: DualNum<f64> + Copy>(
        &self,
        eos: &HelmholtzEnergyWrapper<E, D, N>,
        chemical_records: [&ChemicalRecord<D>; N],
        x: [D; 3],
    ) -> EosResult<(D, Vec<D>)> {
        ...
    }
}
```
### Optimization problem
```rust
// combine molecular representation, property model, and process model in an optimization problem
let campd = IntegratedDesign::new(molecule, property_model, process);

// and pass it to the outer approximation solver
let solver = OuterApproximation::new(&campd);

// to determine a ranking of the optimal molecules for a pur component
solver.solve_ranking(y0, highs, runs, options);
// where highs can be replaced with any of the other MILP solvers provided by the good_lp crate.
```

Molecular representations and property models can be combined according to:
||PC-SAFT|gc-PC-SAFT|
|-|-|-|
| `CoMTCAMD` | yes | no |
| `SuperMolecule` | yes | yes |

### Cite us
If you find FeOs-torch useful for your own research, consider citing our [publication](https://pubs.rsc.org/en/content/articlelanding/2023/me/d2me00230b) from which this library resulted.
```
@article{rehner2023molecule_superstructures,
  author = {Rehner, Philipp and Schilling, Johannes and Bardow, André},
  title = {Molecule superstructures for computer-aided molecular and process design}
  journal = {Molecular Systems Design & Engineering},
  volume = {8},
  issue = {4},
  number = {12},
  pages = {488-499},
  year = {2023}
}
```
