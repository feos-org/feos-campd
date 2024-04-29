#![cfg(feature = "knitro_rs")]
use anyhow::Result;
use approx::assert_relative_eq;
use feos::core::parameter::{IdentifierOption, Parameter, PureRecord};
use feos::core::si::{KELVIN, MOL};
use feos::core::{EosResult, Residual, State};
use feos::pcsaft::{PcSaft, PcSaftBinaryRecord, PcSaftParameters, PcSaftRecord};
use feos_campd::process::ProcessModel;
use feos_campd::{
    CoMTCAMD, OptimizationProblem, OuterApproximationAlgorithm, PcSaftPropertyModel,
    ProcessVariables, Variable,
};
use indexmap::IndexMap;
use ndarray::arr1;
use std::sync::Arc;

#[test]
pub fn test_binary() -> Result<()> {
    let camd = CoMTCAMD::from_molecules(vec![
        "pentane".into(),
        "benzene".into(),
        "diethyl ether".into(),
    ]);
    let pcsaft = PcSaftPropertyModel::from_json_molecules(
        "tests/gross2001.json",
        "tests/poling2000.json",
        Some("tests/mixture_test_comps_binary.json"),
        IdentifierOption::Name,
        true,
    )?;
    let y0 = camd.get_initial_values("molecule", &IndexMap::from([("pentane", 1)]));
    let y1 = camd.get_initial_values("molecule", &IndexMap::from([("diethyl ether", 1)]));
    let y = [y0, y1];
    let camd_binary = [camd.clone(), camd];
    let mut problem = OptimizationProblem::new(camd_binary, pcsaft, CriticalPointModel);
    let result = problem.solve_outer_approximation(
        y,
        OuterApproximationAlgorithm::DuranGrossmann(true),
        Some("tests/options_target.opt"),
        Some("tests/options_MILP.opt"),
    )?;
    println!("{} {}", result.target, result.x[0]);

    let pentane = PureRecord::new(
        Default::default(),
        0.0,
        PcSaftRecord::new(
            2.6896, 3.7729, 231.2, None, None, None, None, None, None, None, None, None, None,
        ),
    );
    let benzene = PureRecord::new(
        Default::default(),
        0.0,
        PcSaftRecord::new(
            2.4653, 3.6478, 287.35, None, None, None, None, None, None, None, None, None, None,
        ),
    );
    let params = PcSaftParameters::new_binary(
        vec![benzene, pentane],
        Some(PcSaftBinaryRecord::new(
            Some(-0.18371738189391876),
            None,
            None,
        )),
    )?;
    let eos = Arc::new(PcSaft::new(Arc::new(params)));
    let (t, _, _) = CriticalPointModel.solve(&eos, &result.x)?;
    println!("{t}");

    assert_relative_eq!(t, result.target, max_relative = 1e-5);

    Ok(())
}

struct CriticalPointModel;

impl<E: Residual> ProcessModel<E> for CriticalPointModel {
    fn variables(&self) -> ProcessVariables {
        vec![("x", Variable::continuous(0.0, 1.0, 0.1))].into()
    }

    fn equality_constraints(&self) -> usize {
        0
    }

    fn inequality_constraints(&self) -> usize {
        0
    }

    fn solve(&self, eos: &Arc<E>, x: &[f64]) -> EosResult<(f64, Vec<f64>, Vec<f64>)> {
        println!("{x:?}");
        let molefracs = arr1(&[x[0], 1.0 - x[0]]);
        let cp = State::critical_point(eos, Some(&(molefracs * MOL)), None, Default::default())?;
        Ok((-cp.temperature.convert_into(KELVIN), vec![], vec![]))
    }
}
