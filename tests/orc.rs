use feos_campd::process::{OrganicRankineCycle, ProcessModel};
use feos_core::parameter::{IdentifierOption, Parameter};
use feos_core::EosResult;
use feos_pcsaft::{PcSaft, PcSaftParameters};
use std::fs::File;
use std::io::BufReader;
use std::rc::Rc;

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
    let eos = Rc::new(PcSaft::new(Rc::new(params)));
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
    let pcsaft = PropertyModel::new_pcsaft("tests/sauer2014_homo_joback.json").unwrap();
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
