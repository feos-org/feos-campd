use feos::pcsaft::{PcSaft, PcSaftParameters};
use feos_campd::process::{OrganicRankineCycle, ProcessModel};
use feos_core::joback::{Joback, JobackParameters, JobackRecord};
use feos_core::parameter::{IdentifierOption, Parameter};
use feos_core::{EosResult, EquationOfState};
use std::fs::File;
use std::io::BufReader;
use std::sync::Arc;

#[test]
pub fn test_orc() -> EosResult<()> {
    let orc: OrganicRankineCycle =
        serde_json::from_reader(BufReader::new(File::open("tests/orc.json").unwrap())).unwrap();
    let params = PcSaftParameters::from_json(
        vec!["propane"],
        "tests/gross2001.json",
        None,
        IdentifierOption::Name,
    )?;
    let pcsaft = Arc::new(PcSaft::new(Arc::new(params)));
    let joback = Arc::new(Joback::new(Arc::new(JobackParameters::from_model_records(
        vec![JobackRecord::new(0.0, 0.0, 0.0, 0.0, 0.0)],
    )?)));
    let eos = Arc::new(EquationOfState::new(joback, pcsaft));
    orc.solve(&eos, &[0.786, -4.5, -1.2, 0.8])?;
    Ok(())
}

#[cfg(feature = "knitro_rs")]
use approx::assert_relative_eq;
#[cfg(feature = "knitro_rs")]
use feos_campd::*;
#[cfg(feature = "knitro_rs")]
use knitro_rs::KnitroError;

#[test]
#[cfg(feature = "knitro_rs")]
fn test_validation() -> Result<(), KnitroError> {
    let molecule = SuperMolecule::alkane(5);
    let orc = OrganicRankineCycle::from_json("tests/orc.json").unwrap();
    let pcsaft =
        PcSaftPropertyModel::new("tests/sauer2014_homo.json", "tests/joback1987.json").unwrap();
    let mut problem = OptimizationProblem::new(molecule, pcsaft, orc);
    #[cfg(feature = "knitro_13")]
    let options = Some("tests/options_13.opt");
    #[cfg(feature = "knitro_12")]
    let options = Some("tests/options_12.opt");
    problem.solve_knitro(&[0.786, -4.5, -1.2, 0.8], 1, options)?;
    assert_relative_eq!(
        -1.5905756277031002,
        problem.solutions[0].target,
        max_relative = 1e-5
    );
    Ok(())
}

#[test]
#[cfg(feature = "knitro_rs")]
fn test_comt_camd() -> Result<(), KnitroError> {
    let orc = OrganicRankineCycle::from_json("tests/orc.json").unwrap();
    let pcsaft =
        PcSaftPropertyModel::new("tests/sauer2014_homo.json", "tests/joback1987.json").unwrap();
    let mut problem = OptimizationProblem::new(CoMTCAMD, pcsaft, orc);
    #[cfg(feature = "knitro_13")]
    let options = Some("tests/options_13.opt");
    #[cfg(feature = "knitro_12")]
    let options = Some("tests/options_12.opt");
    problem.solve_knitro(&[0.786, -4.5, -1.2, 0.8], 1, options)?;
    assert_relative_eq!(
        -1.5905756277031002,
        problem.solutions[0].target,
        max_relative = 1e-5
    );
    Ok(())
}

#[test]
#[cfg(feature = "knitro_rs")]
fn test_supermolecule_disjunct() -> Result<(), KnitroError> {
    let molecule = SuperMoleculeDisjunct::new(5);
    let orc = OrganicRankineCycle::from_json("tests/orc.json").unwrap();
    let pcsaft =
        PcSaftPropertyModel::new("tests/sauer2014_homo.json", "tests/joback1987.json").unwrap();
    let mut problem = OptimizationProblem::new(molecule.clone(), pcsaft, orc);
    #[cfg(feature = "knitro_13")]
    let options = Some("tests/options_13.opt");
    #[cfg(feature = "knitro_12")]
    let options = Some("tests/options_12.opt");
    problem.solve_knitro(&[0.786, -4.5, -1.2, 0.8], 1, options)?;
    assert_relative_eq!(
        -1.5905756277031002,
        problem.solutions[0].target,
        max_relative = 1e-5
    );
    Ok(())
}
