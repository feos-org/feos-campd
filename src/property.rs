use super::molecule::SegmentAndBondCount;
use feos::gc_pcsaft::{GcPcSaft, GcPcSaftChemicalRecord, GcPcSaftEosParameters, GcPcSaftRecord};
use feos::pcsaft::{PcSaft, PcSaftParameters, PcSaftRecord};
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

    fn build_eos(
        &self,
        chemical_records: Vec<C>,
    ) -> Result<Arc<Self::EquationOfState>, ParameterError>;
}

/// A generic group contribution method
#[derive(Clone, Serialize, Deserialize)]
pub struct GcPropertyModel<I, R> {
    residual: Vec<SegmentRecord<R>>,
    ideal_gas: Vec<SegmentRecord<I>>,
    binary: Option<Vec<BinaryRecord<String, f64>>>,
}

impl<I, R> GcPropertyModel<I, R>
where
    for<'de> R: Deserialize<'de>,
    for<'de> I: Deserialize<'de>,
{
    pub fn new<P: AsRef<Path>>(
        file_residual: P,
        file_ideal_gas: P,
        file_binary: Option<P>,
    ) -> Result<Self, ParameterError> {
        Ok(Self {
            residual: SegmentRecord::from_json(file_residual)?,
            ideal_gas: SegmentRecord::from_json(file_ideal_gas)?,
            binary: file_binary
                .map(|f| BinaryRecord::<String, f64>::from_json(f))
                .transpose()?
                .map(|binary| {
                    binary
                        .into_iter()
                        .map(|br| BinaryRecord::new(br.id1, br.id2, br.model_record))
                        .collect()
                }),
        })
    }
}

/// The (homosegmented) group contribution PC-SAFT model.
pub type PcSaftPropertyModel = GcPropertyModel<JobackRecord, PcSaftRecord>;

impl<T: SegmentCount<Count = f64> + Clone> PropertyModel<T> for PcSaftPropertyModel {
    type EquationOfState = EquationOfState<Joback, PcSaft>;

    fn build_eos(
        &self,
        chemical_records: Vec<T>,
    ) -> Result<Arc<EquationOfState<Joback, PcSaft>>, ParameterError> {
        let pcsaft = PcSaft::new(Arc::new(PcSaftParameters::from_segments(
            chemical_records.clone(),
            self.residual.clone(),
            self.binary.clone(),
        )?));
        let joback = Joback::new(Arc::new(JobackParameters::from_segments(
            chemical_records,
            self.ideal_gas.clone(),
            None,
        )?));
        Ok(Arc::new(EquationOfState::new(
            Arc::new(joback),
            Arc::new(pcsaft),
        )))
    }
}

/// The PC-SAFT equation of state for a fixed molecule.
pub struct PcSaftFixedPropertyModel(Arc<EquationOfState<Joback, PcSaft>>);

impl PcSaftFixedPropertyModel {
    pub fn new(eos: &Arc<EquationOfState<Joback, PcSaft>>) -> Self {
        Self(eos.clone())
    }
}

impl PropertyModel<()> for PcSaftFixedPropertyModel {
    type EquationOfState = EquationOfState<Joback, PcSaft>;

    fn build_eos(
        &self,
        _: Vec<()>,
    ) -> Result<Arc<EquationOfState<Joback, PcSaft>>, ParameterError> {
        Ok(self.0.clone())
    }
}

/// The heterosegmented gc-PC-SAFT equation of state.
pub type GcPcSaftPropertyModel = GcPropertyModel<JobackRecord, GcPcSaftRecord>;

impl PropertyModel<SegmentAndBondCount> for GcPcSaftPropertyModel {
    type EquationOfState = EquationOfState<Joback, GcPcSaft>;

    fn build_eos(
        &self,
        chemical_records: Vec<SegmentAndBondCount>,
    ) -> Result<Arc<EquationOfState<Joback, GcPcSaft>>, ParameterError> {
        let gc_pcsaft = GcPcSaft::new(Arc::new(GcPcSaftEosParameters::from_segments(
            chemical_records.clone(),
            self.residual.clone(),
            self.binary.clone(),
        )?));
        let joback = Joback::new(Arc::new(JobackParameters::from_segments(
            chemical_records,
            self.ideal_gas.clone(),
            None,
        )?));
        Ok(Arc::new(EquationOfState::new(
            Arc::new(joback),
            Arc::new(gc_pcsaft),
        )))
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
