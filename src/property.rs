use super::molecule::SegmentAndBondCount;
use feos::gc_pcsaft::{GcPcSaft, GcPcSaftChemicalRecord, GcPcSaftEosParameters, GcPcSaftRecord};
use feos::pcsaft::{PcSaft, PcSaftBinaryRecord, PcSaftParameters, PcSaftRecord};
use feos_core::joback::{Joback, JobackParameters, JobackRecord};
use feos_core::parameter::{
    BinaryRecord, Identifier, Parameter, ParameterError, ParameterHetero, SegmentCount,
    SegmentRecord,
};
use feos_core::{EquationOfState, IdealGas, Residual};
use serde::{Deserialize, Serialize};
use std::path::Path;
use std::sync::Arc;

/// A generic property model to be used in an [OptimizationProblem](super::OptimizationProblem).
pub trait PropertyModel<C> {
    type EquationOfState: Residual + IdealGas;

    fn build_eos(&self, chemical_records: Vec<C>) -> Result<Self::EquationOfState, ParameterError>;
}

/// The (homosegmented) group contribution PC-SAFT model.
#[derive(Clone, Serialize, Deserialize)]
pub struct PcSaftPropertyModel {
    pcsaft: Vec<SegmentRecord<PcSaftRecord>>,
    joback: Vec<SegmentRecord<JobackRecord>>,
    binary: Option<Vec<BinaryRecord<String, PcSaftBinaryRecord>>>,
}

impl PcSaftPropertyModel {
    pub fn new<P: AsRef<Path>>(file_pcsaft: P, file_joback: P) -> Result<Self, ParameterError> {
        Ok(Self {
            pcsaft: SegmentRecord::from_json(file_pcsaft)?,
            joback: SegmentRecord::from_json(file_joback)?,
            binary: None,
        })
    }
}

impl<T: SegmentCount<Count = f64> + Clone> PropertyModel<T> for PcSaftPropertyModel {
    type EquationOfState = EquationOfState<Joback, PcSaft>;

    fn build_eos(
        &self,
        chemical_records: Vec<T>,
    ) -> Result<EquationOfState<Joback, PcSaft>, ParameterError> {
        let pcsaft = PcSaft::new(Arc::new(PcSaftParameters::from_segments(
            chemical_records.clone(),
            self.pcsaft.clone(),
            self.binary.clone(),
        )?));
        let joback = Joback::new(Arc::new(JobackParameters::from_segments(
            chemical_records,
            self.joback.clone(),
            None,
        )?));
        Ok(EquationOfState::new(Arc::new(joback), Arc::new(pcsaft)))
    }
}

// /// The PC-SAFT equation of state for a fixed molecule.
// #[derive(Clone, Serialize, Deserialize)]
// pub struct PcSaftFixedPropertyModel(EquationOfState<Joback, PcSaft>);

// impl PcSaftFixedPropertyModel {
//     pub fn new(pcsaft: Vec<PureRecord<PcSaftRecord>>, joback: Vec<JobackRecord>) -> Self {
//         Self { pcsaft, joback }
//     }
// }

// impl PropertyModel<()> for PcSaftFixedPropertyModel {
//     type Residual = PcSaft;
//     type IdealGas = Joback;

//     fn build_eos(&self, _: Vec<()>) -> Result<EquationOfState<Joback, PcSaft>, ParameterError> {
//         let params = PcSaftParameters::from_records(self.pcsaft.clone(), self.binary.clone())?;
//         let pcsaft = PcSaft::new(Arc::new(params));
//         let joback = Joback::new(Arc::new(JobackParameters::from_model_records(
//             self.joback.clone(),
//         )?));
//         Ok(EquationOfState::new(Arc::new(joback), Arc::new(pcsaft)))
//     }
// }

/// The heterosegmented gc-PC-SAFT equation of state.
#[derive(Clone, Serialize, Deserialize)]
pub struct GcPcSaftPropertyModel(
    Vec<SegmentRecord<GcPcSaftRecord>>,
    Vec<SegmentRecord<JobackRecord>>,
);

impl GcPcSaftPropertyModel {
    pub fn new<P: AsRef<Path>>(file_gc_pcsaft: P, file_joback: P) -> Result<Self, ParameterError> {
        Ok(Self(
            SegmentRecord::from_json(file_gc_pcsaft)?,
            SegmentRecord::from_json(file_joback)?,
        ))
    }
}

impl PropertyModel<SegmentAndBondCount> for GcPcSaftPropertyModel {
    type EquationOfState = EquationOfState<Joback, GcPcSaft>;

    fn build_eos(
        &self,
        chemical_records: Vec<SegmentAndBondCount>,
    ) -> Result<EquationOfState<Joback, GcPcSaft>, ParameterError> {
        let gc_pcsaft = GcPcSaft::new(Arc::new(GcPcSaftEosParameters::from_segments(
            chemical_records.clone(),
            self.0.clone(),
            None,
        )?));
        let joback = Joback::new(Arc::new(JobackParameters::from_segments(
            chemical_records,
            self.1.clone(),
            None,
        )?));
        Ok(EquationOfState::new(Arc::new(joback), Arc::new(gc_pcsaft)))
    }
}

impl From<SegmentAndBondCount> for GcPcSaftChemicalRecord {
    fn from(segment_and_bond_count: SegmentAndBondCount) -> Self {
        Self::new(
            Identifier::default(),
            segment_and_bond_count.segments,
            segment_and_bond_count.bonds,
            1.0,
        )
    }
}
