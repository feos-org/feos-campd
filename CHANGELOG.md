# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.3.1] - 2025-01-27
### Added
- Re-release of 0.3.0 with added documentation for docs.rs.

## [0.3.0] - 2025-01-27
### Changed
- Completely revamped the library, by switching from Knitro to open-source solvers and using automatic implicit differentiation for the calculation of Jacobians and Hessians.

## [0.2.1] - 2024-06-28
### Fixed
- Fix bug in OptimizationProblem::build_eos in https://github.com/feos-org/feos-campd/pull/3
- Fix bug for mixtures using GcPcSaft in https://github.com/feos-org/feos-campd/pull/4

## [0.2.0] - 2024-04-29
### Changed
- Complete rewrite including the implementation of outer approximation algorithms and process models in Python.

## [0.1.0] - 2022-10-19
### Added
- Initial release
