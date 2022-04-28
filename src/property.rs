use feos_core::joback::JobackRecord;
use feos_core::parameter::{ChemicalRecord, Parameter, ParameterError, PureRecord, SegmentRecord};
use feos_core::{EquationOfState, HelmholtzEnergy, IdealGasContribution, MolarWeight};
use feos_gc_pcsaft::{GcPcSaft, GcPcSaftEosParameters, GcPcSaftRecord};
use feos_pcsaft::{PcSaft, PcSaftParameters, PcSaftRecord};
use ndarray::Array1;
use quantity::si::{SIArray1, SIUnit};
use serde::{Deserialize, Serialize};
use std::path::Path;
use std::rc::Rc;

pub enum EquationsOfState {
    PcSaft(PcSaft),
    GcPcSaft(GcPcSaft),
}

impl EquationOfState for EquationsOfState {
    fn components(&self) -> usize {
        match self {
            Self::PcSaft(eos) => eos.components(),
            Self::GcPcSaft(eos) => eos.components(),
        }
    }

    fn subset(&self, component_list: &[usize]) -> Self {
        match self {
            Self::PcSaft(eos) => Self::PcSaft(eos.subset(component_list)),
            Self::GcPcSaft(eos) => Self::GcPcSaft(eos.subset(component_list)),
        }
    }

    fn compute_max_density(&self, moles: &Array1<f64>) -> f64 {
        match self {
            Self::PcSaft(eos) => eos.compute_max_density(moles),
            Self::GcPcSaft(eos) => eos.compute_max_density(moles),
        }
    }

    fn residual(&self) -> &[Box<dyn HelmholtzEnergy>] {
        match self {
            Self::PcSaft(eos) => eos.residual(),
            Self::GcPcSaft(eos) => eos.residual(),
        }
    }

    fn ideal_gas(&self) -> &dyn IdealGasContribution {
        match self {
            Self::PcSaft(eos) => eos.ideal_gas(),
            Self::GcPcSaft(eos) => eos.ideal_gas(),
        }
    }
}

impl MolarWeight<SIUnit> for EquationsOfState {
    fn molar_weight(&self) -> SIArray1 {
        match self {
            Self::PcSaft(eos) => eos.molar_weight(),
            Self::GcPcSaft(eos) => eos.molar_weight(),
        }
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub enum PropertyModel {
    PcSaft(Vec<SegmentRecord<PcSaftRecord, JobackRecord>>),
    GcPcSaft(Vec<SegmentRecord<GcPcSaftRecord, JobackRecord>>),
    PcSaftFixed(Box<PureRecord<PcSaftRecord, JobackRecord>>),
}

impl PropertyModel {
    pub fn new_pcsaft<P: AsRef<Path>>(file: P) -> Result<Self, ParameterError> {
        Ok(Self::PcSaft(SegmentRecord::from_json(file)?))
    }

    pub fn new_gc_pcsaft<P: AsRef<Path>>(file: P) -> Result<Self, ParameterError> {
        Ok(Self::GcPcSaft(SegmentRecord::from_json(file)?))
    }

    pub fn build_eos(
        &self,
        chemical_record: ChemicalRecord,
    ) -> Result<EquationsOfState, ParameterError> {
        Ok(match self {
            Self::PcSaft(segments) => EquationsOfState::PcSaft(PcSaft::new(Rc::new(
                PcSaftParameters::from_segments(vec![chemical_record], segments.clone(), None)?,
            ))),
            Self::GcPcSaft(segments) => EquationsOfState::GcPcSaft(GcPcSaft::new(Rc::new(
                GcPcSaftEosParameters::from_segments(
                    vec![chemical_record],
                    segments.clone(),
                    None,
                )?,
            ))),
            Self::PcSaftFixed(record) => EquationsOfState::PcSaft(PcSaft::new(Rc::new(
                PcSaftParameters::new_pure(*record.clone()),
            ))),
        })
    }
}
