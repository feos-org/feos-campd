# FeOs-CAMPD

[![crate](https://img.shields.io/crates/v/feos-campd.svg)](https://crates.io/crates/feos-campd)
[![documentation](https://docs.rs/feos-campd/badge.svg)](https://docs.rs/feos-campd)

Computer-aided molecular and process design using the [`FeOs`](https://github.com/feos-org/feos) framework.

The crate provides some basic infrastructure to perform a computer-aided molecular and process design. It consists of
- MINLP solver bindings (Artelys Knitro)
- Molecular representations (fixed molecules and molecule superstructures)
- Property models (PC-SAFT and (heterosegmented) gc-PC-SAFT)
- Infrastructure for simple process models including that for an organic Rankine cycle

and the surrounding framework that is used to run the optimization problems.

Actually solving CAMPD problems is only possible with either the `knitro_12` or `knitro_13` feature which require the corresponding version of Artelys Knitro to be installed.

## Installation

Add this to your `Cargo.toml`

```toml
[dependencies]
feos-campd = "0.1"
```