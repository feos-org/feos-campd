use super::{MolecularRepresentation, SuperMolecule};
use crate::SegmentAndBondCount;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Clone, Serialize, Deserialize)]
pub struct SuperMoleculeDisjunct {
    supermolecules: Vec<SuperMolecule>,
}

impl SuperMoleculeDisjunct {
    pub fn new(size: usize) -> Self {
        let supermolecules = vec![
            SuperMolecule::alkane(size),
            SuperMolecule::alkene(size),
            SuperMolecule::alkyne(size),
            SuperMolecule::alcohol(size),
            SuperMolecule::methylether(size),
            SuperMolecule::ketone(size),
            SuperMolecule::amine(size),
        ];
        Self { supermolecules }
    }

    pub fn smiles(&self, y: Vec<usize>) -> String {
        let (y, c) = y.split_at(y.len() - self.supermolecules.len());
        for (&c, supermolecule) in c.iter().zip(self.supermolecules.iter()) {
            if c == 1 {
                return supermolecule.smiles(y.to_vec());
            }
        }
        unreachable!();
    }
}

impl MolecularRepresentation for SuperMoleculeDisjunct {
    type ChemicalRecord = SegmentAndBondCount;

    fn build(&self, y: Vec<f64>) -> Self::ChemicalRecord {
        let (y, c) = y.split_at(y.len() - self.supermolecules.len());
        // println!("{c:?}\n{y:?}");
        let mut cr = SegmentAndBondCount::new(HashMap::new(), HashMap::new());
        for (&c, supermolecule) in c.iter().zip(self.supermolecules.iter()) {
            if c == 0.0 {
                continue;
            }
            let SegmentAndBondCount { segments, bonds } = supermolecule.build(y.to_vec());
            for (segment, count) in segments.into_iter() {
                *cr.segments.entry(segment).or_insert(0.0) += count * c;
            }
            for (bond, count) in bonds.into_iter() {
                *cr.bonds.entry(bond).or_insert(0.0) += count * c;
            }
        }
        cr
    }

    fn variables(&self) -> usize {
        self.supermolecules
            .iter()
            .map(|s| s.variables())
            .max()
            .unwrap()
            + self.supermolecules.len()
    }
}

#[cfg(feature = "knitro_rs")]
mod knitro {
    use super::*;
    use crate::knitro::MolecularRepresentationKnitro;
    use knitro_rs::{Knitro, KnitroError, KN_VARTYPE_BINARY};

    impl MolecularRepresentationKnitro for SuperMoleculeDisjunct {
        fn setup_knitro(&self, kc: &Knitro) -> Result<(), KnitroError> {
            // structure variables
            let index_vars = kc.add_vars(self.variables(), Some(KN_VARTYPE_BINARY))?;
            let (y, c) = index_vars.split_at(index_vars.len() - self.supermolecules.len());

            // disjunction
            let coefs = vec![1.0; c.len()];
            let index_con = kc.add_con()?;
            kc.add_con_linear_struct_one(index_con, c, &coefs)?;
            kc.set_con_eqbnd(index_con, 1.0)?;

            for (&c, supermolecule) in c.iter().zip(self.supermolecules.iter()) {
                // binary variables
                let n_y = supermolecule.variables();
                let (index_vars, unused_vars) = y.split_at(n_y);

                // unused variables
                if !unused_vars.is_empty() {
                    let index_con = kc.add_con()?;
                    let coefs = vec![1.0; unused_vars.len()];
                    let big_m = coefs.len() as f64;
                    kc.add_con_linear_struct_one(index_con, unused_vars, &coefs)?;
                    kc.add_con_linear_term(index_con, c, big_m)?;
                    kc.set_con_upbnd(index_con, big_m)?;
                }

                // maximum size
                let (coefs, size) = supermolecule.size_constraint();
                let big_m = coefs.iter().sum::<f64>() - size;
                let index_con = kc.add_con()?;
                kc.add_con_linear_struct_one(index_con, index_vars, &coefs)?;
                kc.add_con_linear_term(index_con, c, big_m)?;
                kc.set_con_upbnd(index_con, size + big_m)?;

                // minimum size
                let index_con = kc.add_con()?;
                kc.add_con_linear_term(index_con, index_vars[0], 1.0)?;
                kc.add_con_linear_term(index_con, c, -1.0)?;
                kc.set_con_lobnd(index_con, 0.0)?;

                // bond constraints
                for bond in supermolecule.bond_constraints(index_vars[0]) {
                    let index_con = kc.add_con()?;
                    kc.add_con_linear_struct_one(index_con, &bond, &[1.0, -1.0])?;
                    kc.add_con_linear_term(index_con, c, -1.0)?;
                    kc.set_con_lobnd(index_con, -1.0)?;
                }

                // symmetry constraints
                for (index_vars, coefs) in supermolecule.symmetry_constraints(index_vars[0]) {
                    let big_m = coefs.iter().map(|&a| a.min(0.0)).sum::<f64>();
                    let index_con = kc.add_con()?;
                    kc.add_con_linear_struct_one(index_con, &index_vars, &coefs)?;
                    kc.add_con_linear_term(index_con, c, big_m)?;
                    kc.set_con_lobnd(index_con, big_m)?;
                }
            }
            Ok(())
        }
    }
}
