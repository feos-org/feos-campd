use super::molecule::SegmentAndBondCount;
use feos::gc_pcsaft::{GcPcSaft, GcPcSaftChemicalRecord, GcPcSaftEosParameters, GcPcSaftRecord};
use feos::pcsaft::{PcSaft, PcSaftParameters, PcSaftRecord};
use feos_core::joback::{Joback, JobackParameters, JobackRecord};
use feos_core::parameter::{
    Identifier, Parameter, ParameterError, ParameterHetero, PureRecord, SegmentRecord,
};
use feos_core::{EquationOfState, IdealGas, Residual};
use serde::{Deserialize, Serialize};
use std::path::Path;
use std::sync::Arc;

/// A generic property model to be used in an [OptimizationProblem](super::OptimizationProblem).
pub trait PropertyModel<C> {
    type Residual: Residual;
    type IdealGas: IdealGas;
    fn build_eos(
        &self,
        chemical_record: C,
    ) -> Result<EquationOfState<Self::IdealGas, Self::Residual>, ParameterError>;
}

/// The (homosegmented) group contribution PC-SAFT model.
#[derive(Clone, Serialize, Deserialize)]
pub struct PcSaftPropertyModel(
    Vec<SegmentRecord<PcSaftRecord>>,
    Vec<SegmentRecord<JobackRecord>>,
);

impl PcSaftPropertyModel {
    pub fn new<P: AsRef<Path>>(file_pcsaft: P, file_joback: P) -> Result<Self, ParameterError> {
        Ok(Self(
            SegmentRecord::from_json(file_pcsaft)?,
            SegmentRecord::from_json(file_joback)?,
        ))
    }
}

impl PropertyModel<SegmentAndBondCount> for PcSaftPropertyModel {
    type Residual = PcSaft;
    type IdealGas = Joback;

    fn build_eos(
        &self,
        chemical_record: SegmentAndBondCount,
    ) -> Result<EquationOfState<Joback, PcSaft>, ParameterError> {
        let pcsaft = PcSaft::new(Arc::new(PcSaftParameters::from_segments(
            vec![chemical_record.clone()],
            self.0.clone(),
            None,
        )?));
        let joback = Joback::new(Arc::new(JobackParameters::from_segments(
            vec![chemical_record],
            self.1.clone(),
            None,
        )?));
        Ok(EquationOfState::new(Arc::new(joback), Arc::new(pcsaft)))
    }
}

/// The PC-SAFT equation of state for a fixed molecule.
#[derive(Clone, Serialize, Deserialize)]
pub struct PcSaftFixedPropertyModel(PureRecord<PcSaftRecord>, JobackRecord);

impl PcSaftFixedPropertyModel {
    pub fn new(pure_record: PureRecord<PcSaftRecord>, joback_record: JobackRecord) -> Self {
        Self(pure_record, joback_record)
    }
}

impl PropertyModel<()> for PcSaftFixedPropertyModel {
    type Residual = PcSaft;
    type IdealGas = Joback;

    fn build_eos(&self, _: ()) -> Result<EquationOfState<Joback, PcSaft>, ParameterError> {
        let pcsaft = PcSaft::new(Arc::new(PcSaftParameters::new_pure(self.0.clone())?));
        let joback = Joback::new(Arc::new(JobackParameters::from_model_records(vec![self
            .1
            .clone()])?));
        Ok(EquationOfState::new(Arc::new(joback), Arc::new(pcsaft)))
    }
}

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
    type Residual = GcPcSaft;
    type IdealGas = Joback;

    fn build_eos(
        &self,
        chemical_record: SegmentAndBondCount,
    ) -> Result<EquationOfState<Joback, GcPcSaft>, ParameterError> {
        let gc_pcsaft = GcPcSaft::new(Arc::new(GcPcSaftEosParameters::from_segments(
            vec![chemical_record.clone()],
            self.0.clone(),
            None,
        )?));
        let joback = Joback::new(Arc::new(JobackParameters::from_segments(
            vec![chemical_record],
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
