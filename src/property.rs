use crate::variables::ExplicitVariable;
use crate::Constraint;

// use super::molecule::SegmentAndBondCount;
use feos::core::parameter::{
    BinaryRecord, Identifier, Parameter, ParameterError, ParameterHetero, SegmentCount,
    SegmentRecord,
};
use feos::core::{EquationOfState, IdealGas, Residual};
use feos::gc_pcsaft::{GcPcSaft, GcPcSaftChemicalRecord, GcPcSaftEosParameters, GcPcSaftRecord};
use feos::ideal_gas::{Joback, JobackRecord};
use feos::pcsaft::{PcSaft, PcSaftParameters, PcSaftRecord};
#[cfg(feature = "knitro_rs")]
use knitro_rs::{Knitro, KnitroError};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::sync::Arc;

/// A generic property model to be used in an [OptimizationProblem](super::OptimizationProblem).
pub trait PropertyModel {
    type EquationOfState: Residual + IdealGas;

    fn parameter_variables(
        &self,
        index_feature_vars: &HashMap<String, i32>,
    ) -> Vec<ExplicitVariable>;

    fn evaluate_parameter_variables(&self, feature_vars: &HashMap<String, f64>) -> Vec<f64> {
        let (index_feature_vars, feature_vars): (HashMap<_, _>, Vec<_>) = feature_vars
            .iter()
            .enumerate()
            .map(|(i, (k, v))| ((k.clone(), i as i32), v))
            .unzip();
        self.parameter_variables(&index_feature_vars)
            .into_iter()
            .map(|v| v.evaluate(&feature_vars))
            .collect()
    }

    fn symmetry_constraint(&self, index_parameter_vars: &[i32]) -> Option<Constraint>;

    fn build_eos(&self, parameters: &[f64]) -> Arc<Self::EquationOfState>;

    #[cfg(feature = "knitro_rs")]
    fn setup_knitro(
        &self,
        kc: &Knitro,
        index_feature_vars: &HashMap<String, i32>,
    ) -> Result<Vec<i32>, KnitroError> {
        // define parameter variables
        let index_parameter_vars = self
            .parameter_variables(index_feature_vars)
            .iter()
            .map(|f| f.setup_knitro(kc))
            .collect::<Result<Vec<_>, _>>()?;

        // Add symmetry constraint
        self.symmetry_constraint(&index_parameter_vars)
            .map(|c| c.setup_knitro(kc));

        Ok(index_parameter_vars)
    }
}

/// A generic group contribution method
#[derive(Clone, Serialize, Deserialize)]
pub struct GcPropertyModel<I, R> {
    residual: Vec<SegmentRecord<R>>,
    ideal_gas: Vec<SegmentRecord<I>>,
    binary: Option<Vec<BinaryRecord<String, f64>>>,
}

// impl<I, R> GcPropertyModel<I, R>
// where
//     for<'de> R: Deserialize<'de>,
//     for<'de> I: Deserialize<'de>,
// {
//     pub fn new<P: AsRef<Path>>(
//         file_residual: P,
//         file_ideal_gas: P,
//         file_binary: Option<P>,
//     ) -> Result<Self, ParameterError> {
//         Ok(Self {
//             residual: SegmentRecord::from_json(file_residual)?,
//             ideal_gas: SegmentRecord::from_json(file_ideal_gas)?,
//             binary: file_binary
//                 .map(|f| BinaryRecord::<String, f64>::from_json(f))
//                 .transpose()?
//                 .map(|binary| {
//                     binary
//                         .into_iter()
//                         .map(|br| BinaryRecord::new(br.id1, br.id2, br.model_record))
//                         .collect()
//                 }),
//         })
//     }
// }

// /// The (homosegmented) group contribution PC-SAFT model.
// pub type PcSaftPropertyModel = GcPropertyModel<JobackRecord, PcSaftRecord>;

// impl<T: SegmentCount<Count = f64> + Clone, const N: usize> PropertyModel<[T; N]>
//     for PcSaftPropertyModel
// {
//     type EquationOfState = EquationOfState<Joback, PcSaft>;

//     fn build_eos(
//         &self,
//         chemical_records: [T; N],
//     ) -> Result<Arc<EquationOfState<Joback, PcSaft>>, ParameterError> {
//         let pcsaft = PcSaft::new(Arc::new(PcSaftParameters::from_segments(
//             chemical_records.to_vec(),
//             self.residual.clone(),
//             self.binary.clone(),
//         )?));
//         let joback =
//             Joback::from_segments(chemical_records.to_vec(), self.ideal_gas.clone(), None)?;
//         Ok(Arc::new(EquationOfState::new(
//             Arc::new(joback),
//             Arc::new(pcsaft),
//         )))
//     }
// }

// // /// The PC-SAFT equation of state for a fixed molecule.
// // pub struct PcSaftFixedPropertyModel(Arc<EquationOfState<Joback, PcSaft>>);

// // impl PcSaftFixedPropertyModel {
// //     pub fn new(eos: &Arc<EquationOfState<Joback, PcSaft>>) -> Self {
// //         Self(eos.clone())
// //     }
// // }

// // impl PropertyModel<()> for PcSaftFixedPropertyModel {
// //     type EquationOfState = EquationOfState<Joback, PcSaft>;

// //     fn build_eos(
// //         &self,
// //         _: Vec<()>,
// //     ) -> Result<Arc<EquationOfState<Joback, PcSaft>>, ParameterError> {
// //         Ok(self.0.clone())
// //     }
// // }

/// The heterosegmented gc-PC-SAFT equation of state.
pub type GcPcSaftPropertyModel = GcPropertyModel<JobackRecord, GcPcSaftRecord>;

impl PropertyModel for GcPcSaftPropertyModel {
    type EquationOfState = EquationOfState<Joback, GcPcSaft>;

    fn parameter_variables(
        &self,
        index_feature_vars: &HashMap<String, i32>,
    ) -> Vec<ExplicitVariable> {
        todo!()
    }

    fn symmetry_constraint(&self, _: &[i32]) -> Option<Constraint> {
        todo!()
    }

    fn build_eos(&self, parameters: &[f64]) -> Arc<EquationOfState<Joback, GcPcSaft>> {
        todo!()
        // let gc_pcsaft = GcPcSaft::new(Arc::new(GcPcSaftEosParameters::from_segments(
        //     chemical_records.to_vec(),
        //     self.residual.clone(),
        //     self.binary.clone(),
        // )?));
        // let joback =
        //     Joback::from_segments(chemical_records.to_vec(), self.ideal_gas.clone(), None)?;
        // Ok(Arc::new(EquationOfState::new(
        //     Arc::new(joback),
        //     Arc::new(gc_pcsaft),
        // )))
    }
}

// impl From<SegmentAndBondCount> for GcPcSaftChemicalRecord {
//     fn from(segment_and_bond_count: SegmentAndBondCount) -> Self {
//         Self::new(
//             Identifier::default(),
//             segment_and_bond_count.segments,
//             segment_and_bond_count.bonds,
//             1.0,
//         )
//     }
// }
