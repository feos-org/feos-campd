use anyhow::Result;
use feos::core::parameter::{IdentifierOption, Parameter};
use feos::core::{EosResult, EquationOfState};
use feos::ideal_gas::{Joback, JobackRecord};
use feos::pcsaft::{PcSaft, PcSaftParameters};
use feos_campd::process::{OrganicRankineCycle, ProcessModel};
use feos_campd::ProcessVariables;
#[cfg(feature = "knitro_rs")]
use indexmap::IndexMap;
use std::fs::File;
use std::io::BufReader;
use std::sync::Arc;

#[test]
pub fn test_orc() -> Result<()> {
    let orc: OrganicRankineCycle =
        serde_json::from_reader(BufReader::new(File::open("tests/orc.json")?))?;
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

#[test]
#[cfg(feature = "knitro_rs")]
fn test_validation() -> Result<()> {
    let molecule = SuperMolecule::alkane(5);
    let orc = OrganicRankineCycle::from_json("tests/orc.json")?;
    let pcsaft = PcSaftPropertyModel::from_json_groups(
        "tests/sauer2014_homo.json",
        "tests/joback1987.json",
        None,
        false,
    )?;
    let mut problem = OptimizationProblem::new([molecule], pcsaft, orc);
    #[cfg(feature = "knitro_13")]
    let options = Some("tests/options_13.opt");
    #[cfg(feature = "knitro_12")]
    let options = Some("tests/options_12.opt");
    let (_, x, y) =
        problem.solve_target(&[0.786, -4.5, -1.2, 0.8], Some("tests/options_target.opt"))?;
    println!("{x:?}");
    println!("{y:?}");
    let solution = problem.solve_knitro_once(&x, Some(&y), options)?;
    assert_relative_eq!(-1.5905756277031002, solution.target, max_relative = 1e-5);
    println!("{solution}");
    Ok(())
}

#[test]
#[cfg(feature = "knitro_rs")]
fn test_validation_hetero() -> Result<()> {
    let molecule = [SuperMolecule::alkane(5)];
    let orc = OrganicRankineCycle::from_json("tests/orc.json")?;
    let pcsaft = GcPcSaftPropertyModel::from_json(
        &molecule,
        "tests/sauer2014_hetero.json",
        "tests/joback1987.json",
        None,
    )?;
    let mut problem = OptimizationProblem::new(molecule.clone(), pcsaft.clone(), orc.clone());
    let solution = problem.solve_outer_approximation(
        [vec![1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0]],
        OuterApproximationAlgorithm::DuranGrossmann(false),
        Some("tests/options_target.opt"),
        Some("tests/options_MILP.opt"),
    )?;
    assert_relative_eq!(-1.5640076533528697, solution.target, max_relative = 1e-5);

    let mut problem = OptimizationProblem::new(molecule, pcsaft, orc);
    #[cfg(feature = "knitro_13")]
    let options = Some("tests/options_13.opt");
    #[cfg(feature = "knitro_12")]
    let options = Some("tests/options_12.opt");
    let (_, x, y) =
        problem.solve_target(&[0.786, -4.5, -1.2, 0.8], Some("tests/options_target.opt"))?;
    println!("{x:?}");
    println!("{y:?}");
    let solution = problem.solve_knitro_once(&x, Some(&y), options)?;
    assert_relative_eq!(-1.5640076533528697, solution.target, max_relative = 1e-5);
    println!("{solution}");
    Ok(())
}

#[test]
#[cfg(feature = "knitro_rs")]
fn test_target() -> Result<()> {
    let molecule = SuperMolecule::alkane(5);
    let orc = OrganicRankineCycle::from_json("tests/orc.json")?;
    let pcsaft = PcSaftPropertyModel::from_json_groups(
        "tests/sauer2014_homo.json",
        "tests/joback1987.json",
        None,
        false,
    )?;
    let mut problem = OptimizationProblem::new([molecule], pcsaft, orc);
    let (target, x, y) =
        problem.solve_target(&[0.786, -4.5, -1.2, 0.8], Some("tests/options_target.opt"))?;
    println!("{target} {y:?} {x:?}");
    assert_relative_eq!(target, -1.8145704740258306, max_relative = 1e-3);
    Ok(())
}

#[test]
#[cfg(feature = "knitro_rs")]
fn test_fixed() -> Result<()> {
    use indexmap::IndexMap;

    let molecule = SuperMolecule::alkane(5);
    let orc = OrganicRankineCycle::from_json("tests/orc.json")?;
    let pcsaft = PcSaftPropertyModel::from_json_groups(
        "tests/sauer2014_homo.json",
        "tests/joback1987.json",
        None,
        false,
    )?;
    let mut problem = OptimizationProblem::new([molecule], pcsaft.clone(), orc.clone());
    let res = problem.solve_fixed(
        &[vec![1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]],
        Some("tests/options_target.opt"),
    )?;
    let camd = CoMTCAMD::from_json("tests/comt_camd_non_assoc.json")?;
    let mut problem = OptimizationProblem::new([camd], pcsaft, orc);
    println!("{} {:?} {:?}", res.target, res.y, res.x);
    assert_relative_eq!(res.target, -1.5905756277031002, max_relative = 1e-5);

    let y0 = problem.molecules[0]
        .get_initial_values("alkanes", &IndexMap::from([("CH3", 2), ("CH2", 1)]));
    let res = problem.solve_fixed(&[y0], Some("tests/options_target.opt"))?;
    println!("{}", res.target);
    println!("{:?}", res.x);
    println!("{:?}", res.y);
    assert_relative_eq!(res.target, -1.5905756277031002, max_relative = 1e-5);
    Ok(())
}

#[test]
#[cfg(feature = "knitro_rs")]
fn test_outer_approximation_duran() -> Result<()> {
    let orc = OrganicRankineCycle::from_json("tests/orc.json")?;
    let camd = CoMTCAMD::from_json("tests/comt_camd_non_assoc.json")?;
    let pcsaft = PcSaftPropertyModel::from_json_groups(
        "tests/sauer2014_homo.json",
        "tests/joback1987.json",
        None,
        false,
    )?;
    let mut problem = OptimizationProblem::new([camd], pcsaft, orc);
    let y0 = problem.molecules[0].get_initial_values(
        "alkanes",
        &IndexMap::from([("CH3", 3), ("CH2", 1), (">CH", 1)]),
    );
    let solution = problem.solve_outer_approximation(
        [y0],
        OuterApproximationAlgorithm::DuranGrossmann(false),
        Some("tests/options_target.opt"),
        Some("tests/options_MILP.opt"),
    )?;
    assert_relative_eq!(-1.5905756277031002, solution.target, max_relative = 1e-5);
    Ok(())
}

#[test]
#[cfg(feature = "knitro_rs")]
fn test_outer_approximation_fletcher() -> Result<()> {
    let orc = OrganicRankineCycle::from_json("tests/orc.json")?;
    let camd = CoMTCAMD::from_json("tests/comt_camd_non_assoc.json")?;
    let pcsaft = PcSaftPropertyModel::from_json_groups(
        "tests/sauer2014_homo.json",
        "tests/joback1987.json",
        None,
        false,
    )?;
    let mut problem = OptimizationProblem::new([camd], pcsaft, orc);
    let y0 = problem.molecules[0].get_initial_values(
        "alkenes",
        &IndexMap::from([("CH3", 1), ("=CH2", 1), ("=CH", 1)]),
    );
    let solution = problem.solve_outer_approximation(
        [y0],
        OuterApproximationAlgorithm::FletcherLeyffer,
        Some("tests/options_target.opt"),
        Some("tests/options_MILP.opt"),
    )?;
    assert_relative_eq!(-1.5905756277031002, solution.target, max_relative = 1e-5);
    Ok(())
}

#[test]
#[cfg(feature = "knitro_rs")]
fn test_outer_approximation_molecules() -> Result<()> {
    let orc = OrganicRankineCycle::from_json("tests/orc.json")?;
    let camd = CoMTCAMD::from_json_molecules("tests/gross_2001_comps.json")?;
    let pcsaft = PcSaftPropertyModel::from_json_molecules(
        "tests/gross2001.json",
        "tests/poling2000.json",
        None,
        IdentifierOption::Name,
        false,
    )?;
    let mut problem = OptimizationProblem::new([camd], pcsaft, orc);
    let y0 = problem.molecules[0]
        .get_initial_values("molecule", &IndexMap::from([("dimethyl ether", 1)]));
    let solution = problem.solve_outer_approximation(
        [y0],
        OuterApproximationAlgorithm::DuranGrossmann(true),
        Some("tests/options_target.opt"),
        Some("tests/options_MILP.opt"),
    )?;
    assert_relative_eq!(-1.6261766509457936, solution.target, max_relative = 1e-5);
    Ok(())
}

#[test]
#[cfg(feature = "knitro_rs")]
fn test_outer_approximation_ranking() -> Result<()> {
    let orc = OrganicRankineCycle::from_json("tests/orc.json")?;
    let camd = CoMTCAMD::from_json("tests/comt_camd_non_assoc.json")?;
    let pcsaft = PcSaftPropertyModel::from_json_groups(
        "tests/sauer2014_homo.json",
        "tests/joback1987.json",
        None,
        false,
    )?;
    let mut problem = OptimizationProblem::new([camd], pcsaft, orc);
    let y0 = problem.molecules[0].get_initial_values(
        "alkanes",
        &IndexMap::from([("CH3", 3), ("CH2", 1), (">CH", 1)]),
    );
    problem.outer_approximation_ranking(
        [y0],
        OuterApproximationAlgorithm::FletcherLeyffer,
        15,
        Some("tests/options_target.opt"),
        Some("tests/options_MILP.opt"),
    );
    Ok(())
}

#[test]
#[cfg(feature = "knitro_rs")]
fn test_outer_approximation_ranking_molecules() -> Result<()> {
    let orc = OrganicRankineCycle::from_json("tests/orc.json")?;
    let camd = CoMTCAMD::from_json_molecules("tests/gross_2001_comps.json")?;
    let pcsaft = PcSaftPropertyModel::from_json_molecules(
        "tests/gross2001.json",
        "tests/poling2000.json",
        None,
        IdentifierOption::Name,
        false,
    )?;
    let mut problem = OptimizationProblem::new([camd], pcsaft, orc);
    let y0 =
        problem.molecules[0].get_initial_values("molecule", &IndexMap::from([("cyclopentane", 1)]));
    // .get_initial_values("molecule", &IndexMap::from([("ethylene", 1)]));
    problem.outer_approximation_ranking(
        [y0],
        OuterApproximationAlgorithm::DuranGrossmann(true),
        10,
        Some("tests/options_target.opt"),
        Some("tests/options_MILP.opt"),
    );
    Ok(())
}

// #[test]
// #[cfg(feature = "knitro_rs")]
// fn test_supermolecule_disjunct() -> Result<()> {
//     let molecule = SuperMolecule::non_associating(5);
//     let orc = OrganicRankineCycle::from_json("tests/orc.json")?;
//     let pcsaft = PcSaftPropertyModel::from_json("tests/comt_camd_non_assoc.json")?;
//     let mut problem = OptimizationProblem::new([molecule], pcsaft, orc);
//     #[cfg(feature = "knitro_13")]
//     let options = Some("tests/options_13.opt");
//     #[cfg(feature = "knitro_12")]
//     let options = Some("tests/options_12.opt");
//     let (_, x, y) =
//         problem.solve_target(&[0.786, -4.5, -1.2, 0.8], Some("tests/options_ms.opt"))?;
//     let solution = problem.solve_knitro_once(&x, Some(&y), options)?;
//     assert_relative_eq!(-1.5905756277031002, solution.target, max_relative = 1e-5);
//     Ok(())
// }

#[test]
#[cfg(feature = "knitro_rs")]
fn test_supermolecule_disjunct_oa() -> Result<()> {
    let molecule = SuperMolecule::non_associating(5);
    let orc = OrganicRankineCycle::from_json("tests/orc.json")?;
    let pcsaft = PcSaftPropertyModel::from_json_groups(
        "tests/sauer2014_homo.json",
        "tests/joback1987.json",
        None,
        false,
    )?;
    let mut problem = OptimizationProblem::new([molecule], pcsaft, orc);
    // println!("{}", molecule.structure_variables().len());
    let y = [1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0];
    let v0 = y;
    let v1 = [0.0; 7];
    let v2 = [0.0; 7];
    let v3 = [0.0; 7];
    let v4 = [0.0; 7];
    let c = [1.0, 0.0, 0.0, 0.0, 0.0];
    let y0 = [&y[..], &v0, &v1, &v2, &v3, &v4, &c].concat();
    problem.solve_fixed(&[y0.clone()], Some("tests/options_target.opt"))?;
    problem.outer_approximation_ranking(
        [y0],
        OuterApproximationAlgorithm::FletcherLeyffer,
        3,
        Some("tests/options_target.opt"),
        Some("tests/options_MILP.opt"),
    );
    let mut solutions: Vec<_> = problem.solutions.iter().collect();
    solutions.sort_by(|s1, s2| s1.target.total_cmp(&s2.target));
    println!("{:?}", solutions[0]);
    assert_relative_eq!(
        -1.5905756277031002,
        solutions[0].target,
        max_relative = 1e-5
    );
    Ok(())
}

struct NoModel;
impl<E> ProcessModel<E> for NoModel {
    fn variables(&self) -> ProcessVariables {
        vec![].into()
    }

    fn equality_constraints(&self) -> usize {
        0
    }

    fn inequality_constraints(&self) -> usize {
        0
    }

    fn solve(&self, _: &Arc<E>, _: &[f64]) -> EosResult<(f64, Vec<f64>, Vec<f64>)> {
        Ok((0.0, vec![], vec![]))
    }
}
