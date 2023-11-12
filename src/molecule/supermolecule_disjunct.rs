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

    pub fn non_associating(size: usize) -> Self {
        let supermolecules = vec![
            SuperMolecule::alkane(size),
            SuperMolecule::alkene(size),
            SuperMolecule::alkyne(size),
            SuperMolecule::methylether(size),
            SuperMolecule::ketone(size),
        ];
        Self { supermolecules }
    }
}

impl MolecularRepresentation<1> for SuperMoleculeDisjunct {
    type ChemicalRecord = SegmentAndBondCount;

    fn variables(&self) -> usize {
        self.supermolecules
            .iter()
            .map(|s| s.variables())
            .max()
            .unwrap()
            + self.supermolecules.len()
    }

    fn build(&self, y: &[f64]) -> [Self::ChemicalRecord; 1] {
        let (y, c) = y.split_at(y.len() - self.supermolecules.len());
        // println!("{c:?}\n{y:?}");
        let mut cr = SegmentAndBondCount::new(HashMap::new(), HashMap::new());
        for (&c, supermolecule) in c.iter().zip(self.supermolecules.iter()) {
            if c == 0.0 {
                continue;
            }
            let [SegmentAndBondCount { segments, bonds }] = supermolecule.build(y);
            for (segment, count) in segments.into_iter() {
                *cr.segments.entry(segment).or_insert(0.0) += count * c;
            }
            for (bond, count) in bonds.into_iter() {
                *cr.bonds.entry(bond).or_insert(0.0) += count * c;
            }
        }
        [cr]
    }

    fn smiles(&self, y: &[usize]) -> [String; 1] {
        let (y, c) = y.split_at(y.len() - self.supermolecules.len());
        for (&c, supermolecule) in c.iter().zip(self.supermolecules.iter()) {
            if c == 1 {
                return supermolecule.smiles(y);
            }
        }
        unreachable!();
    }
}

#[cfg(feature = "knitro_rs")]
mod knitro {
    use super::*;
    use crate::knitro::MolecularRepresentationKnitro;
    use knitro_rs::{Knitro, KnitroError, KN_VARTYPE_BINARY};

    impl MolecularRepresentationKnitro<1> for SuperMoleculeDisjunct {
        fn setup_knitro(&self, kc: &Knitro) -> Result<(), KnitroError> {
            // structure variables
            let index_vars = kc.add_vars(self.variables(), Some(KN_VARTYPE_BINARY))?;
            let (y, c) = index_vars.split_at(index_vars.len() - self.supermolecules.len());
            kc.set_var_primal_initial_values(y, &vec![0.5; y.len()])?;
            kc.set_var_primal_initial_values(c, &vec![1.0 / c.len() as f64; c.len()])?;

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

#[cfg(test)]
mod test {
    use super::*;
    use approx::assert_relative_eq;

    fn molecule() -> SuperMoleculeDisjunct {
        SuperMoleculeDisjunct {
            supermolecules: vec![
                SuperMolecule::alcohol(5),
                SuperMolecule::methylether(5),
                SuperMolecule::ketone(5),
            ],
        }
    }

    #[test]
    fn test_build() {
        println!("{}", molecule().variables());
        let [SegmentAndBondCount { segments, bonds }] =
            molecule().build(&[1.0, 0.8, 0.6, 0.1, 0.3, 0.2, 0.1, 0.5, 0.3, 0.2]);
        println!("{segments:?}\n{bonds:?}");
        assert_relative_eq!(segments["CH3"], 1.464);
        assert_relative_eq!(segments["CH2"], 1.17);
        assert_relative_eq!(segments[">CH"], 0.218);
        assert_relative_eq!(segments[">C<"], 0.008);
        assert_relative_eq!(segments["OH"], 0.5);
        assert_relative_eq!(segments["OCH3"], 0.3);
        assert_relative_eq!(segments[">C=O"], 0.06);
        assert_relative_eq!(segments["CH=O"], 0.14);
    }

    #[test]
    fn test_ethanol() {
        let [SegmentAndBondCount { segments, bonds }] =
            molecule().build(&[1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]);
        println!("{segments:?}\n{bonds:?}");
        assert_eq!(segments["CH3"], 1.0);
        assert_eq!(segments["CH2"], 1.0);
        assert_eq!(segments["OH"], 1.0);
        assert_relative_eq!(bonds[&["CH3".to_string(), "CH2".to_string()]], 1.0);
        assert_relative_eq!(bonds[&["OH".to_string(), "CH2".to_string()]], 1.0);
    }

    #[test]
    fn test_2_methoxybutane() {
        let [SegmentAndBondCount { segments, bonds }] =
            molecule().build(&[1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]);
        println!("{segments:?}\n{bonds:?}");
        assert_eq!(segments["CH3"], 2.0);
        assert_eq!(segments[">CH"], 1.0);
        assert_eq!(segments["OCH3"], 1.0);
        assert_relative_eq!(bonds[&["CH3".to_string(), ">CH".to_string()]], 2.0);
        assert_relative_eq!(bonds[&["OCH3".to_string(), ">CH".to_string()]], 1.0);
    }

    #[test]
    fn test_acetone() {
        let [SegmentAndBondCount { segments, bonds }] =
            molecule().build(&[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0]);
        println!("{segments:?}\n{bonds:?}");
        assert_eq!(segments["CH3"], 2.0);
        assert_eq!(segments[">C=O"], 1.0);
        assert_relative_eq!(bonds[&["CH3".to_string(), ">C=O".to_string()]], 2.0);
    }
}
