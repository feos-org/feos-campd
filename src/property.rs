use crate::variables::ExplicitVariable;
use crate::{Constraint, MolecularRepresentation};

use super::molecule::SegmentAndBondCount;
use feos::core::parameter::{
    BinaryRecord, Identifier, Parameter, ParameterError, ParameterHetero, SegmentRecord,
};
use feos::core::{EquationOfState, IdealGas, Residual};
use feos::gc_pcsaft::{GcPcSaft, GcPcSaftChemicalRecord, GcPcSaftEosParameters, GcPcSaftRecord};
use feos::ideal_gas::{Joback, JobackRecord};
use indexmap::IndexMap;
#[cfg(feature = "knitro_rs")]
use knitro_rs::{Knitro, KnitroError};
use serde::{Deserialize, Serialize};
use std::path::Path;
use std::sync::Arc;

/// A generic property model to be used in an [OptimizationProblem](super::OptimizationProblem).
pub trait PropertyModel {
    type EquationOfState: Residual + IdealGas;

    fn parameter_variables(
        &self,
        index_feature_vars: &IndexMap<String, i32>,
    ) -> Vec<ExplicitVariable>;

    fn evaluate_parameter_variables(&self, feature_vars: &IndexMap<String, f64>) -> Vec<f64> {
        let (index_feature_vars, feature_vars): (IndexMap<_, _>, Vec<_>) = feature_vars
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
        index_feature_vars: &IndexMap<String, i32>,
    ) -> Result<Vec<i32>, KnitroError> {
        // define parameter variables
        let index_parameter_vars = self
            .parameter_variables(index_feature_vars)
            .into_iter()
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
pub struct GcPcSaftPropertyModel {
    residual: Vec<SegmentRecord<GcPcSaftRecord>>,
    ideal_gas: Vec<SegmentRecord<JobackRecord>>,
    binary: Option<Vec<BinaryRecord<String, f64>>>,
    groups: Vec<(String, usize)>,
    bonds: Vec<([String; 2], usize)>,
}

impl GcPcSaftPropertyModel {
    pub fn new<P: AsRef<Path>, M: MolecularRepresentation>(
        molecule: M,
        file_residual: P,
        file_ideal_gas: P,
        file_binary: Option<P>,
    ) -> Result<Self, ParameterError> {
        let index_structure_vars: Vec<_> =
            (0..molecule.structure_variables().len() as i32).collect();
        let mut groups = Vec::new();
        let mut bonds = Vec::new();
        molecule
            .feature_variables(&index_structure_vars)
            .into_keys()
            .enumerate()
            .for_each(|(i, s)| {
                if let Some(index) = s.find('-') {
                    let (g1, g2) = s.split_at(index);
                    bonds.push(([g1.to_string(), g2[1..].to_string()], i));
                } else {
                    groups.push((s, i))
                }
            });
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
            groups,
            bonds,
        })
    }
}

impl PropertyModel for GcPcSaftPropertyModel {
    type EquationOfState = EquationOfState<Joback, GcPcSaft>;

    fn parameter_variables(
        &self,
        index_feature_vars: &IndexMap<String, i32>,
    ) -> Vec<ExplicitVariable> {
        index_feature_vars
            .iter()
            .map(|(n, &i)| ExplicitVariable::new(n.clone()).linear_struct(vec![i], vec![1.0]))
            .collect()
    }

    fn symmetry_constraint(&self, _: &[i32]) -> Option<Constraint> {
        None
    }

    fn build_eos(&self, parameters: &[f64]) -> Arc<EquationOfState<Joback, GcPcSaft>> {
        let segments = self
            .groups
            .iter()
            .cloned()
            .map(|(s, i)| (s, parameters[i]))
            .collect();
        let bonds = self
            .bonds
            .iter()
            .cloned()
            .map(|(s, i)| (s, parameters[i]))
            .collect();
        let chemical_records = vec![SegmentAndBondCount::new(segments, bonds)];

        let gc_pcsaft = GcPcSaft::new(Arc::new(
            GcPcSaftEosParameters::from_segments(
                chemical_records.clone(),
                self.residual.clone(),
                self.binary.clone(),
            )
            .unwrap(),
        ));
        let joback = Joback::from_segments(chemical_records, self.ideal_gas.clone(), None).unwrap();
        Arc::new(EquationOfState::new(Arc::new(joback), Arc::new(gc_pcsaft)))
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
