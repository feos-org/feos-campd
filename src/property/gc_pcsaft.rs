use super::PropertyModel;
use crate::variables::ExplicitVariable;
use crate::{Constraint, MolecularRepresentation};
use feos::core::parameter::{
    BinaryRecord, Identifier, Parameter, ParameterError, ParameterHetero, SegmentRecord,
};
use feos::core::EquationOfState;
use feos::gc_pcsaft::{GcPcSaft, GcPcSaftChemicalRecord, GcPcSaftEosParameters, GcPcSaftRecord};
use feos::ideal_gas::{Joback, JobackRecord};
use indexmap::IndexMap;
use std::path::Path;
use std::sync::Arc;

/// The heterosegmented gc-PC-SAFT equation of state.
#[derive(Clone)]
pub struct GcPcSaftPropertyModel {
    residual: Vec<SegmentRecord<GcPcSaftRecord>>,
    ideal_gas: Vec<SegmentRecord<JobackRecord>>,
    binary: Option<Vec<BinaryRecord<String, f64>>>,
    groups: Vec<Vec<(String, usize)>>,
    bonds: Vec<Vec<([String; 2], usize)>>,
}

impl GcPcSaftPropertyModel {
    pub fn from_json<P: AsRef<Path>, M: MolecularRepresentation>(
        molecules: &[M],
        file_residual: P,
        file_ideal_gas: P,
        file_binary: Option<P>,
    ) -> Result<Self, ParameterError> {
        let mut i = 0;
        let (groups, bonds) = molecules
            .iter()
            .map(|molecule| {
                let mut groups = Vec::new();
                let mut bonds = Vec::new();
                let index_structure_vars: Vec<_> =
                    (0..molecule.structure_variables().len() as i32).collect();
                molecule
                    .feature_variables(&index_structure_vars)
                    .into_keys()
                    .for_each(|s| {
                        if let Some(index) = s.find('-') {
                            let (g1, g2) = s.split_at(index);
                            bonds.push(([g1.to_string(), g2[1..].to_string()], i));
                        } else {
                            groups.push((s, i))
                        }
                        i += 1;
                    });
                (groups, bonds)
            })
            .unzip();
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

impl<const N: usize> PropertyModel<N> for GcPcSaftPropertyModel {
    type EquationOfState = EquationOfState<Joback, GcPcSaft>;

    fn parameter_variables(
        &self,
        index_feature_vars: &[IndexMap<String, i32>; N],
    ) -> Vec<ExplicitVariable> {
        index_feature_vars
            .iter()
            .flatten()
            .map(|(n, &i)| ExplicitVariable::new(n.clone()).linear_struct(vec![i], vec![1.0]))
            .collect()
    }

    fn symmetry_constraints(&self, _: &[IndexMap<String, i32>; N], _: &[i32]) -> Vec<Constraint> {
        Vec::new()
    }

    fn build_eos(&self, parameters: &[f64]) -> Arc<EquationOfState<Joback, GcPcSaft>> {
        let chemical_records: Vec<_> = self
            .groups
            .iter()
            .zip(&self.bonds)
            .map(|(groups, bonds)| {
                let segments = groups
                    .iter()
                    .cloned()
                    .map(|(s, i)| (s, parameters[i]))
                    .filter(|(_, n)| *n > 0.0)
                    .collect();
                let bonds = bonds
                    .iter()
                    .cloned()
                    .map(|(s, i)| (s, parameters[i]))
                    .filter(|(_, n)| *n > 0.0)
                    .collect();
                GcPcSaftChemicalRecord::new(Identifier::default(), segments, bonds, 1.0)
            })
            .collect();

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
