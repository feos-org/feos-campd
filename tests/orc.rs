use feos::core::parameter::{IdentifierOption, Parameter};
use feos::core::{EosResult, EquationOfState, IdealGas, Residual};
use feos::ideal_gas::{Joback, JobackRecord};
use feos::pcsaft::{PcSaft, PcSaftParameters};
use feos_campd::process::{OrganicRankineCycle, ProcessModel};
use feos_campd::ContinuousVariables;
#[cfg(feature = "knitro_rs")]
use std::collections::HashMap;
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
    let joback = Arc::new(Joback::from_model_records(vec![JobackRecord::new(
        0.0, 0.0, 0.0, 0.0, 0.0,
    )])?);
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
        PcSaftPropertyModel::new("tests/sauer2014_homo.json", "tests/joback1987.json", None)
            .unwrap();
    let mut problem = OptimizationProblem::new(molecule, pcsaft, orc);
    #[cfg(feature = "knitro_13")]
    let options = Some("tests/options_13.opt");
    #[cfg(feature = "knitro_12")]
    let options = Some("tests/options_12.opt");
    let (_, [x, y, p]) =
        problem.solve_target(&[0.786, -4.5, -1.2, 0.8], Some("tests/options_target.opt"))?;
    println!("{x:?}");
    println!("{y:?}");
    println!("{p:?}");
    let solution = problem.solve_knitro_once(&x, Some(&y), Some(&p), options)?;
    assert_relative_eq!(-1.5905756277031002, solution.target, max_relative = 1e-5);
    println!("{solution}");
    Ok(())
}

#[test]
#[cfg(feature = "knitro_rs")]
fn test_target() -> Result<(), KnitroError> {
    let molecule = SuperMolecule::alkane(5);
    let orc = OrganicRankineCycle::from_json("tests/orc.json").unwrap();
    let pcsaft =
        PcSaftPropertyModel::new("tests/sauer2014_homo.json", "tests/joback1987.json", None)
            .unwrap();
    let mut problem = OptimizationProblem::new(molecule, pcsaft, orc);
    let (target, [x, y, _]) =
        problem.solve_target(&[0.786, -4.5, -1.2, 0.8], Some("tests/options_target.opt"))?;
    println!("{target} {y:?} {x:?}");
    assert_relative_eq!(target, -1.805999183, max_relative = 1e-5);
    Ok(())
}

#[test]
#[cfg(feature = "knitro_rs")]
fn test_fixed() -> Result<(), KnitroError> {
    let molecule = SuperMolecule::alkane(5);
    let orc = OrganicRankineCycle::from_json("tests/orc.json").unwrap();
    let pcsaft =
        PcSaftPropertyModel::new("tests/sauer2014_homo.json", "tests/joback1987.json", None)
            .unwrap();
    let mut problem = OptimizationProblem::new(molecule, pcsaft, orc.clone());
    let (target, [x, y, _]) = problem.solve_fixed(
        Some(&[0.786, -4.5, -1.2, 0.8]),
        &[1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        Some("tests/options_target.opt"),
    )?;
    let mut problem = OptimizationProblem::new(
        CoMTCAMD::from_json("tests/comt_camd_non_assoc.json").unwrap(),
        CoMTCAMDPropertyModel,
        orc,
    );
    println!("{target} {y:?} {x:?}");
    assert_relative_eq!(target, -1.5905756277031002, max_relative = 1e-5);

    let y0 = problem
        .molecules
        .get_initial_values("alkanes", &HashMap::from([("CH3", 2), ("CH2", 1)]));
    let (f, [x, y, p]) = problem.solve_fixed(
        Some(&[0.786, -4.5, -1.2, 0.0]),
        &y0,
        Some("tests/options_target.opt"),
    )?;
    println!("{f}");
    println!("{x:?}");
    println!("{y:?}");
    println!("{p:?}");
    assert_relative_eq!(f, -1.5905756277031002, max_relative = 1e-5);
    Ok(())
}

#[test]
#[cfg(feature = "knitro_rs")]
fn test_outer_approximation() -> EosResult<()> {
    let orc = OrganicRankineCycle::from_json("tests/orc.json")?;
    let pcsaft = CoMTCAMDPropertyModel;
    let mut problem = OptimizationProblem::new(
        CoMTCAMD::from_json("tests/comt_camd_non_assoc.json")?,
        pcsaft,
        orc,
    );
    let y0 = problem.molecules.get_initial_values(
        "alkanes",
        &HashMap::from([("CH3", 3), ("CH2", 1), (">CH", 1)]),
    );
    let solution = problem.solve_outer_approximation(
        y0,
        false,
        Some("tests/options_target.opt"),
        Some("tests/options_MILP.opt"),
    )?;
    assert_relative_eq!(-1.5905756277031002, solution.target, max_relative = 1e-5);
    Ok(())
}

#[test]
#[cfg(feature = "knitro_rs")]
fn test_outer_approximation_molecules() -> EosResult<()> {
    let orc = OrganicRankineCycle::from_json("tests/orc.json")?;
    let pcsaft = CoMTCAMDPropertyModel;
    let mut problem = OptimizationProblem::new(
        CoMTCAMD::from_json("tests/gross_2001_comps.json")?,
        pcsaft,
        orc,
    );
    let y0 = problem
        .molecules
        .get_initial_values("molecule", &HashMap::from([("cyclopentane", 1)]));
    let solution = problem.solve_outer_approximation(
        y0,
        false,
        Some("tests/options_target.opt"),
        Some("tests/options_MILP.opt"),
    )?;
    assert_relative_eq!(-1.6261766497437158, solution.target, max_relative = 1e-5);
    Ok(())
}

#[test]
#[cfg(feature = "knitro_rs")]
fn test_outer_approximation_ranking() -> EosResult<()> {
    let orc = OrganicRankineCycle::from_json("tests/orc.json")?;
    let pcsaft = CoMTCAMDPropertyModel;
    let mut problem = OptimizationProblem::new(
        CoMTCAMD::from_json("tests/comt_camd_non_assoc.json")?,
        pcsaft,
        orc,
    );
    let y0 = problem.molecules.get_initial_values(
        "alkanes",
        &HashMap::from([("CH3", 3), ("CH2", 1), (">CH", 1)]),
    );
    problem.outer_approximation_ranking(
        &y0,
        false,
        10,
        Some("tests/options_target.opt"),
        Some("tests/options_MILP.opt"),
    );
    Ok(())
}

#[test]
#[cfg(feature = "knitro_rs")]
fn test_outer_approximation_ranking_molecules() -> EosResult<()> {
    let orc = OrganicRankineCycle::from_json("tests/orc.json")?;
    let pcsaft = CoMTCAMDPropertyModel;
    let mut problem = OptimizationProblem::new(
        CoMTCAMD::from_json("tests/gross_2001_comps.json")?,
        pcsaft,
        orc,
    );
    let y0 = problem
        .molecules
        .get_initial_values("molecule", &HashMap::from([("cyclopentane", 1)]));
    problem.outer_approximation_ranking(
        &y0,
        false,
        10,
        Some("tests/options_target.opt"),
        Some("tests/options_MILP.opt"),
    );
    Ok(())
}

#[test]
#[cfg(feature = "knitro_rs")]
fn test_supermolecule_disjunct() -> Result<(), KnitroError> {
    let molecule = SuperMolecule::all(5);
    let orc = OrganicRankineCycle::from_json("tests/orc.json").unwrap();
    let pcsaft =
        PcSaftPropertyModel::new("tests/sauer2014_homo.json", "tests/joback1987.json", None)
            .unwrap();
    let mut problem = OptimizationProblem::new(molecule, pcsaft, orc);
    #[cfg(feature = "knitro_13")]
    let options = Some("tests/options_13.opt");
    #[cfg(feature = "knitro_12")]
    let options = Some("tests/options_12.opt");
    let (_, [x, y, p]) =
        problem.solve_target(&[0.786, -4.5, -1.2, 0.8], Some("tests/options_target.opt"))?;
    let solution = problem.solve_knitro_once(&x, Some(&y), Some(&p), options)?;
    assert_relative_eq!(-1.5905756277031002, solution.target, max_relative = 1e-5);
    Ok(())
}

#[test]
#[cfg(feature = "knitro_rs")]
fn test_supermolecule_mix() {
    let molecules = (SuperMolecule::alcohol(5), SuperMolecule::ketone(5));
    let pcsaft =
        PcSaftPropertyModel::new("tests/sauer2014_homo.json", "tests/joback1987.json", None)
            .unwrap();
    let mut problem = OptimizationProblem::new(molecules, pcsaft, NoModel);
    #[cfg(feature = "knitro_13")]
    let options = Some("tests/options_13.opt");
    #[cfg(feature = "knitro_12")]
    let options = Some("tests/options_12.opt");
    problem.solve_knitro(&[], None, None, 10, options)
}

struct NoModel;
impl ProcessModel for NoModel {
    fn variables(&self) -> ContinuousVariables {
        vec![].into()
    }

    fn constraints(&self) -> Vec<[Option<f64>; 2]> {
        vec![]
    }

    fn solve<E: Residual + IdealGas>(
        &self,
        _: &Arc<E>,
        _: &[f64],
    ) -> EosResult<(f64, Vec<f64>, Vec<f64>)> {
        Ok((0.0, vec![], vec![]))
    }
}
