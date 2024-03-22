#![cfg(feature = "knitro_rs")]
use feos::core::parameter::{BinaryRecord, Parameter, PureRecord};
use feos::core::si::{KELVIN, MOL};
use feos::core::{EosResult, Residual, State};
use feos::pcsaft::{PcSaft, PcSaftBinaryRecord, PcSaftParameters, PcSaftRecord};
use feos_campd::process::ProcessModel;
use feos_campd::{
    CoMTCAMD, CoMTCAMDBinary, CoMTCAMDBinaryPropertyModel, CoMTCAMDPropertyModel,
    OptimizationProblem, OuterApproximationAlgorithm, ProcessVariables, Variable,
};
use knitro_rs::KnitroError;
use ndarray::arr1;
use std::collections::HashMap;
use std::sync::Arc;

#[test]
pub fn test_binary() -> Result<(), KnitroError> {
    let camd = CoMTCAMD::from_json("tests/mixture_test_comps.json").unwrap();
    let pcsaft = CoMTCAMDPropertyModel::from_json("tests/mixture_test_comps.json").unwrap();
    let y0 = camd.get_initial_values("molecule", &HashMap::from([("pentane", 1)]));
    let y1 = camd.get_initial_values("molecule", &HashMap::from([("diethyl ether", 1)]));
    let y = [y0, y1].concat();
    let k_ij: Vec<BinaryRecord<String, f64>> =
        BinaryRecord::from_json("tests/mixture_test_comps_binary.json").unwrap();
    let camd_binary = CoMTCAMDBinary::new([camd.clone(), camd], true);
    let pcsaft_binary =
        CoMTCAMDBinaryPropertyModel::new([pcsaft.clone(), pcsaft], Some(k_ij), true, false);
    let mut problem = OptimizationProblem::new(camd_binary, pcsaft_binary, CriticalPointModel);
    let result = problem
        .solve_outer_approximation(
            y,
            OuterApproximationAlgorithm::DuranGrossmann(true),
            Some("tests/options_target.opt"),
            Some("tests/options_MILP.opt"),
        )
        .unwrap();
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
        vec![pentane, benzene],
        Some(PcSaftBinaryRecord::new(
            Some(-0.18371738189391876),
            None,
            None,
        )),
    )
    .unwrap();
    let eos = Arc::new(PcSaft::new(Arc::new(params)));
    let (t, _, _) = CriticalPointModel.solve(&eos, &result.x).unwrap();
    println!("{t}");

    assert_eq!(t, result.target);

    Ok(())
}

struct CriticalPointModel;

impl<E: Residual> ProcessModel<E> for CriticalPointModel {
    fn variables(&self) -> ProcessVariables {
        vec![Variable::continuous("x".into(), 0.0, 1.0, 0.5)].into()
    }

    fn equality_constraints(&self) -> usize {
        0
    }

    fn inequality_constraints(&self) -> usize {
        0
    }

    fn solve(&self, eos: &Arc<E>, x: &[f64]) -> EosResult<(f64, Vec<f64>, Vec<f64>)> {
        let molefracs = arr1(&[x[0], 1.0 - x[0]]);
        let cp = State::critical_point(eos, Some(&(molefracs * MOL)), None, Default::default())?;
        Ok((-cp.temperature.convert_into(KELVIN), vec![], vec![]))
    }
}
