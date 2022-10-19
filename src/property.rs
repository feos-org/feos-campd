use super::molecule::SegmentAndBondCount;
use feos::gc_pcsaft::{GcPcSaft, GcPcSaftChemicalRecord, GcPcSaftEosParameters, GcPcSaftRecord};
use feos::pcsaft::{PcSaft, PcSaftParameters, PcSaftRecord};
use feos_core::joback::JobackRecord;
use feos_core::parameter::{
    Identifier, Parameter, ParameterError, ParameterHetero, PureRecord, SegmentRecord,
};
use feos_core::{EquationOfState, MolarWeight};
use quantity::si::SIUnit;
use serde::{Deserialize, Serialize};
use std::path::Path;
use std::rc::Rc;

/// A generic property model to be used in an [OptimizationProblem](super::OptimizationProblem).
pub trait PropertyModel<C> {
    type Eos: EquationOfState + MolarWeight<SIUnit>;
    fn build_eos(&self, chemical_record: C) -> Result<Self::Eos, ParameterError>;
}

/// The (homosegmented) group contribution PC-SAFT model.
#[derive(Clone, Serialize, Deserialize)]
pub struct PcSaftPropertyModel(Vec<SegmentRecord<PcSaftRecord, JobackRecord>>);

impl PcSaftPropertyModel {
    pub fn new<P: AsRef<Path>>(file: P) -> Result<Self, ParameterError> {
        Ok(Self(SegmentRecord::from_json(file)?))
    }
}

impl PropertyModel<SegmentAndBondCount> for PcSaftPropertyModel {
    type Eos = PcSaft;

    fn build_eos(&self, chemical_record: SegmentAndBondCount) -> Result<PcSaft, ParameterError> {
        Ok(PcSaft::new(Rc::new(PcSaftParameters::from_segments(
            vec![chemical_record],
            self.0.clone(),
            None,
        )?)))
    }
}

/// The PC-SAFT equation of state for a fixed molecule.
#[derive(Clone, Serialize, Deserialize)]
pub struct PcSaftFixedPropertyModel(PureRecord<PcSaftRecord, JobackRecord>);

impl PcSaftFixedPropertyModel {
    pub fn new(pure_record: PureRecord<PcSaftRecord, JobackRecord>) -> Self {
        Self(pure_record)
    }
}

impl PropertyModel<()> for PcSaftFixedPropertyModel {
    type Eos = PcSaft;

    fn build_eos(&self, _: ()) -> Result<PcSaft, ParameterError> {
        Ok(PcSaft::new(Rc::new(PcSaftParameters::new_pure(
            self.0.clone(),
        ))))
    }
}

/// The heterosegmented gc-PC-SAFT equation of state.
#[derive(Clone, Serialize, Deserialize)]
pub struct GcPcSaftPropertyModel(Vec<SegmentRecord<GcPcSaftRecord, JobackRecord>>);

impl GcPcSaftPropertyModel {
    pub fn new<P: AsRef<Path>>(file: P) -> Result<Self, ParameterError> {
        Ok(Self(SegmentRecord::from_json(file)?))
    }
}

impl PropertyModel<SegmentAndBondCount> for GcPcSaftPropertyModel {
    type Eos = GcPcSaft;

    fn build_eos(&self, chemical_record: SegmentAndBondCount) -> Result<GcPcSaft, ParameterError> {
        Ok(GcPcSaft::new(Rc::new(
            GcPcSaftEosParameters::from_segments(vec![chemical_record], self.0.clone(), None)?,
        )))
    }
}

impl From<SegmentAndBondCount> for GcPcSaftChemicalRecord {
    fn from(segment_and_bond_count: SegmentAndBondCount) -> Self {
        Self::new(
            Identifier::default(),
            segment_and_bond_count.segments,
            segment_and_bond_count.bonds,
        )
    }
}
