use super::polynomial::{Polynomial, Polynomial2};
use super::MolecularRepresentation;
use feos_core::parameter::{ChemicalRecord, Identifier};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::iter;

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct FunctionalGroup {
    groups: Vec<String>,
    smiles: Vec<String>,
    atoms: usize,
}

impl FunctionalGroup {
    pub fn ch3() -> Self {
        Self {
            groups: vec!["CH3".into()],
            smiles: vec!["C".into()],
            atoms: 1,
        }
    }

    pub fn ctch() -> Self {
        Self {
            groups: vec!["C≡CH".into()],
            smiles: vec!["C#C".into()],
            atoms: 2,
        }
    }

    pub fn oh() -> Self {
        Self {
            groups: vec!["OH".into()],
            smiles: vec!["O".into()],
            atoms: 1,
        }
    }

    pub fn nh2() -> Self {
        Self {
            groups: vec!["NH2".into()],
            smiles: vec!["N".into()],
            atoms: 1,
        }
    }

    pub fn och3() -> Self {
        Self {
            groups: vec!["OCH3".into()],
            smiles: vec!["CO".into()],
            atoms: 2,
        }
    }

    pub fn cdo() -> Self {
        Self {
            groups: vec![">C=O".into(), "CH=O".into()],
            smiles: vec!["C(=O)".into(), "".into()],
            atoms: 2,
        }
    }

    pub fn cdc() -> Self {
        Self {
            groups: vec!["=C<".into(), "=CH".into(), "=CH2".into()],
            smiles: vec!["C".into(), "".into(), "=C".into(), "".into()],
            atoms: 2,
        }
    }
}

pub struct SuperAlkyl;

impl SuperAlkyl {
    fn variables(size: usize) -> usize {
        1 + (1..=3)
            .map(|k| (size - 1) / k)
            .filter(|&s| s > 0)
            .map(Self::variables)
            .sum::<usize>()
    }

    fn bond_constraints(size: usize, index: &mut i32, constraints: &mut Vec<[i32; 2]>) {
        let this = *index;
        (1..=3)
            .map(|k| (size - 1) / k)
            .filter(|s| *s > 0)
            .for_each(|s| {
                constraints.push([this as i32, *index + 1]);
                *index += 1;
                Self::bond_constraints(s, index, constraints);
            });
    }

    fn symmetry_constraints(
        size: usize,
        index: &mut i32,
        level: f64,
        constraints: &mut Vec<(Vec<i32>, Vec<f64>)>,
    ) -> (Vec<i32>, Vec<f64>) {
        let mut indices = vec![*index];
        let mut levels = vec![level];
        let children: Vec<_> = (1..=3)
            .map(|k| (size - 1) / k)
            .filter(|s| *s > 0)
            .map(|s| {
                *index += 1;
                Self::symmetry_constraints(s, index, level + 1.0, constraints)
            })
            .collect();

        children
            .iter()
            .zip(children.iter().skip(1))
            .for_each(|(branch1, branch2)| {
                constraints.push(merge_symmetry_constraints(branch1, branch2));
            });

        children.into_iter().for_each(|(i, l)| {
            indices.extend(i);
            levels.extend(l);
        });
        (indices, levels)
    }

    fn build<I: Iterator<Item = f64>>(
        size: usize,
        y: &mut I,
        segments: &mut Polynomial<4>,
        bonds: &mut Polynomial2<4, 4>,
    ) -> Polynomial<4> {
        let value = y.next().unwrap();

        let children: Vec<_> = (1..=3)
            .map(|k| (size - 1) / k)
            .filter(|s| *s > 0)
            .map(|s| Self::build(s, y, segments, bonds))
            .collect();
        let res = calculate_bonds(children, bonds) * value;
        *segments += res;

        res
    }

    fn smiles<I: Iterator<Item = usize>>(size: usize, y: &mut I) -> Vec<char> {
        if size == 0 {
            return vec![];
        }
        let this = y.next().unwrap();
        let mut smiles = vec!['(', 'C'];
        (1..=3)
            .map(|k| (size - 1) / k)
            .filter(|&s| s > 0)
            .for_each(|s| smiles.extend(Self::smiles(s, y)));
        smiles.push(')');
        if this == 1 {
            smiles
        } else {
            vec![]
        }
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub struct SuperMolecule {
    pub size: usize,
    functional_groups: Vec<FunctionalGroup>,
    alkyls: Vec<usize>,
    symmetries: Vec<[usize; 2]>,
}

impl SuperMolecule {
    pub fn new(size: usize, functional_groups: Option<Vec<FunctionalGroup>>) -> Self {
        let functional_groups = functional_groups.unwrap_or_else(|| {
            vec![
                FunctionalGroup::ch3(),
                FunctionalGroup::ctch(),
                FunctionalGroup::oh(),
                FunctionalGroup::och3(),
            ]
        });
        Self {
            size,
            functional_groups,
            alkyls: vec![1],
            symmetries: vec![],
        }
    }

    pub fn alkane(size: usize) -> Self {
        Self {
            size,
            functional_groups: vec![FunctionalGroup::ch3()],
            alkyls: vec![1],
            symmetries: vec![],
        }
    }

    pub fn alkene(size: usize) -> Self {
        Self {
            size,
            functional_groups: vec![FunctionalGroup::cdc()],
            alkyls: vec![1, 2, 2, 4],
            symmetries: vec![[0, 1], [0, 2], [2, 3]],
        }
    }

    pub fn alkyne(size: usize) -> Self {
        Self {
            size,
            functional_groups: vec![FunctionalGroup::ctch()],
            alkyls: vec![1],
            symmetries: vec![],
        }
    }

    pub fn alcohol(size: usize) -> Self {
        Self {
            size,
            functional_groups: vec![FunctionalGroup::oh()],
            alkyls: vec![1],
            symmetries: vec![],
        }
    }

    pub fn methylether(size: usize) -> Self {
        Self {
            size,
            functional_groups: vec![FunctionalGroup::och3()],
            alkyls: vec![1],
            symmetries: vec![],
        }
    }

    pub fn ketone(size: usize) -> Self {
        Self {
            size,
            functional_groups: vec![FunctionalGroup::cdo()],
            alkyls: vec![1, 2],
            symmetries: vec![[0, 1]],
        }
    }

    pub fn amine(size: usize) -> Self {
        Self {
            size,
            functional_groups: vec![FunctionalGroup::nh2()],
            alkyls: vec![1],
            symmetries: vec![],
        }
    }

    pub fn all(size: usize) -> HashMap<String, Self> {
        HashMap::from([
            ("alkane".into(), Self::alkane(size)),
            ("alkene".into(), Self::alkene(size)),
            ("alkyne".into(), Self::alkyne(size)),
            ("alcohol".into(), Self::alcohol(size)),
            ("methylether".into(), Self::methylether(size)),
            ("ketone".into(), Self::ketone(size)),
            ("amine".into(), Self::amine(size)),
        ])
    }
}

impl SuperMolecule {
    fn functional_group_atoms(&self) -> usize {
        self.functional_groups
            .iter()
            .map(|f| f.atoms)
            .reduce(|a, b| a.min(b))
            .unwrap()
    }

    fn alkyl_tails(&self) -> Vec<usize> {
        self.alkyls
            .iter()
            .map(|a| (self.size - self.functional_group_atoms()) / a)
            .collect()
    }

    pub fn size_constraint(&self) -> (Vec<f64>, usize) {
        (
            self.functional_groups
                .iter()
                .map(|f| f.atoms as f64)
                .chain(iter::repeat(1.0).take(self.variables() - self.functional_groups.len()))
                .collect(),
            self.size,
        )
    }

    pub fn functional_group_constraint(&self) -> Vec<i32> {
        (0..self.functional_groups.len())
            .into_iter()
            .map(|f| f as i32)
            .collect()
    }

    pub fn bond_constraints(&self) -> Vec<[i32; 2]> {
        let mut index = self.functional_groups.len() as i32;
        let mut constraints = Vec::new();
        self.alkyl_tails()
            .iter()
            .filter(|&&s| s > 0)
            .for_each(|&s| {
                SuperAlkyl::bond_constraints(s, &mut index, &mut constraints);
                index += 1;
            });
        constraints
    }

    pub fn symmetry_constraints(&self) -> Vec<(Vec<i32>, Vec<f64>)> {
        let mut index = self.functional_groups.len() as i32;
        let mut constraints = Vec::new();
        let children: Vec<_> = self
            .alkyl_tails()
            .iter()
            .filter(|&&s| s > 0)
            .map(|&s| {
                let child = SuperAlkyl::symmetry_constraints(s, &mut index, 1.0, &mut constraints);
                index += 1;
                child
            })
            .collect();
        self.symmetries.iter().for_each(|&[a, b]| {
            if children.len() > b {
                constraints.push(merge_symmetry_constraints(&children[a], &children[b]))
            }
        });
        constraints
    }

    pub fn smiles(&self, y: Vec<usize>) -> String {
        let y_iter = &mut y.into_iter();
        let functional_group = self
            .functional_groups
            .iter()
            .filter(|_| y_iter.next().unwrap() == 1)
            .last()
            .unwrap();
        self.alkyl_tails()
            .into_iter()
            .zip(functional_group.smiles.iter())
            .fold(String::new(), |mut acc, (s, sm)| {
                acc.push_str(sm);
                acc.push_str(
                    &SuperAlkyl::smiles(s, y_iter)
                        .into_iter()
                        .collect::<String>(),
                );
                acc
            })
    }
}

impl MolecularRepresentation for SuperMolecule {
    fn variables(&self) -> usize {
        self.functional_groups.len()
            + self
                .alkyl_tails()
                .iter()
                .filter(|&&s| s > 0)
                .map(|&s| SuperAlkyl::variables(s))
                .sum::<usize>()
    }

    fn build(&self, y: Vec<f64>) -> ChemicalRecord {
        match self.alkyls.len() {
            1 => (),
            2 => return self.build_ketone(y),
            4 => return self.build_alkene(y),
            _ => unreachable!(),
        }
        let y_iter = &mut y.into_iter();
        let y_func: Vec<_> = self
            .functional_groups
            .iter()
            .map(|_| y_iter.next().unwrap())
            .collect();

        let mut c_segments = Polynomial::zero();
        let mut c_bonds = Polynomial2::zero();
        let s = SuperAlkyl::build(
            self.size - self.functional_group_atoms(),
            y_iter,
            &mut c_segments,
            &mut c_bonds,
        );

        let identifier = Identifier::new("", Some("SuperMolecule"), None, None, None, None);
        let cs = ["CH3", "CH2", ">CH", ">C<"];

        let mut segments: HashMap<_, _> = cs
            .iter()
            .map(|&s| s.to_string())
            .zip(c_segments.iter().cloned())
            .filter(|(_, g)| *g > 0.0)
            .collect();
        self.functional_groups
            .iter()
            .zip(y_func.iter())
            .filter(|(_, &y)| y > 0.0)
            .for_each(|(s, y)| *segments.entry(s.groups[0].clone()).or_insert(0.0) += y);

        let mut bonds = HashMap::new();
        fill_bond_map(cs, cs, &c_bonds, &mut bonds);
        s.iter().zip(cs.iter()).for_each(|(g, c)| {
            if *g > 0.0 {
                self.functional_groups
                    .iter()
                    .zip(y_func.iter())
                    .for_each(|(s, y)| {
                        if *y > 0.0 {
                            *bonds
                                .entry([s.groups[0].clone(), c.to_string()])
                                .or_insert(0.0) += y * g
                        }
                    });
            }
        });

        ChemicalRecord::new_count(identifier, segments, Some(bonds))
    }
}

impl SuperMolecule {
    fn build_ketone(&self, y: Vec<f64>) -> ChemicalRecord {
        let y_iter = &mut y.into_iter();
        y_iter.next();

        let cs = ["CH3", "CH2", ">CH", ">C<"];
        let mut c_segments = Polynomial::zero();
        let mut c_bonds = Polynomial2::zero();

        // formaldehyde is approximated by CH=O for the relaxation.
        let cds = ["CH=O", "CH=O", ">C=O"];
        let mut cd_bonds: Polynomial2<4, 3> = Polynomial2::zero();

        let alkyls = self
            .alkyl_tails()
            .iter()
            .filter(|&&s| s > 0)
            .map(|&s| SuperAlkyl::build(s, y_iter, &mut c_segments, &mut c_bonds))
            .collect();
        let cd_segments = calculate_bonds(alkyls, &mut cd_bonds);

        let identifier = Identifier::new("", Some("SuperKetone"), None, None, None, None);

        let segments: HashMap<_, _> = cs
            .iter()
            .chain(cds.iter())
            .map(|&s| s.to_string())
            .zip(c_segments.iter().chain(cd_segments.iter()).cloned())
            .filter(|(_, g)| *g > 0.0)
            .collect();

        let mut bonds = HashMap::new();
        fill_bond_map(cs, cs, &c_bonds, &mut bonds);
        fill_bond_map(cs, cds, &cd_bonds, &mut bonds);

        ChemicalRecord::new_count(identifier, segments, Some(bonds))
    }

    fn build_alkene(&self, y: Vec<f64>) -> ChemicalRecord {
        let y_iter = &mut y.into_iter();
        y_iter.next();

        let cs = ["CH3", "CH2", ">CH", ">C<"];
        let mut c_segments = Polynomial::zero();
        let mut c_bonds = Polynomial2::zero();

        let cds = ["=CH2", "=CH", "=C<"];
        let mut cd_segments = Vec::new();
        let mut cd_bonds: Polynomial2<4, 3> = Polynomial2::zero();

        let s = self.alkyl_tails();
        for tails in [[0, 1], [2, 3]] {
            let alkyls = tails
                .iter()
                .filter(|&&i| s[i] > 0)
                .map(|&i| SuperAlkyl::build(s[i], y_iter, &mut c_segments, &mut c_bonds))
                .collect();
            cd_segments.push(calculate_bonds(alkyls, &mut cd_bonds));
        }

        let identifier = Identifier::new("", Some("SuperAlkene"), None, None, None, None);

        let segments: HashMap<_, _> = cs
            .iter()
            .chain(cds.iter())
            .map(|&s| s.to_string())
            .zip(
                c_segments
                    .iter()
                    .chain((cd_segments[0] + cd_segments[1]).iter())
                    .cloned(),
            )
            .filter_map(|(s, g)| (g > 0.0).then(|| (s, g)))
            .collect();

        let mut bonds = HashMap::new();
        fill_bond_map(cs, cs, &c_bonds, &mut bonds);
        fill_bond_map(cs, cds, &cd_bonds, &mut bonds);
        let dd_bonds = cd_segments[0].outer_product(cd_segments[1]);
        fill_bond_map(cds, cds, &dd_bonds, &mut bonds);

        ChemicalRecord::new_count(identifier, segments, Some(bonds))
    }
}

fn merge_symmetry_constraints(
    branch1: &(Vec<i32>, Vec<f64>),
    branch2: &(Vec<i32>, Vec<f64>),
) -> (Vec<i32>, Vec<f64>) {
    let (indices1, levels1) = branch1;
    let (indices2, levels2) = branch2;
    (
        indices1
            .iter()
            .cloned()
            .chain(indices2.iter().cloned())
            .collect(),
        levels1
            .iter()
            .cloned()
            .chain(levels2.iter().map(|&l| -l))
            .collect(),
    )
}

fn calculate_bonds<const M: usize, const N: usize>(
    children: Vec<Polynomial<N>>,
    bonds: &mut Polynomial2<N, M>,
) -> Polynomial<M> {
    for (i, ci) in children.iter().enumerate() {
        *bonds += ci.outer_product(
            children
                .iter()
                .enumerate()
                .filter_map(|(j, cj)| (j != i).then(|| Polynomial::new(cj.sum())))
                .product::<Polynomial<M>>()
                * Polynomial::new(1.0),
        );
    }
    children.iter().map(|p| Polynomial::new(p.sum())).product()
}

fn fill_bond_map<const M: usize, const N: usize>(
    seg1: [&str; M],
    seg2: [&str; N],
    bond_poly: &Polynomial2<M, N>,
    bonds: &mut HashMap<[String; 2], f64>,
) {
    for (i, c1) in seg1.iter().enumerate() {
        for (j, c2) in seg2.iter().enumerate() {
            if bond_poly[(i, j)] > 0.0 {
                *bonds.entry([c1.to_string(), c2.to_string()]).or_insert(0.0) += bond_poly[(i, j)];
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use approx::assert_relative_eq;
    use itertools::Itertools;

    fn isomers(molecule: SuperMolecule) -> usize {
        let mut res = 0;
        let (size_constraint, size) = molecule.size_constraint();
        let functional_group_constraint = molecule.functional_group_constraint();
        let bond_constraints = molecule.bond_constraints();
        let symmetry_constraints = molecule.symmetry_constraints();
        for y in vec![0..=1; molecule.variables()]
            .into_iter()
            .multi_cartesian_product()
        {
            if y[functional_group_constraint.len()] == 0 {
                continue;
            }
            if y.iter()
                .zip(size_constraint.iter())
                .map(|(y, &s)| y * s as usize)
                .sum::<usize>()
                > size
            {
                continue;
            }
            if functional_group_constraint
                .iter()
                .map(|&i| y[i as usize])
                .sum::<usize>()
                != 1
            {
                continue;
            }
            if bond_constraints
                .iter()
                .any(|&[a, b]| y[a as usize] < y[b as usize])
            {
                continue;
            }
            if symmetry_constraints.iter().any(|(i, l)| {
                let x = i
                    .iter()
                    .zip(l.iter())
                    .map(|(&i, &l)| y[i as usize] as isize * l as isize)
                    .sum::<isize>();
                x < 0
            }) {
                continue;
            }
            println!("{:?}", y);
            res += 1;
        }
        res
    }

    // https://pubs.acs.org/doi/pdf/10.1021/ja01359a027
    #[test]
    fn test_isomers_alcohols() {
        assert_eq!(isomers(SuperMolecule::alcohol(4)), 4);
        assert_eq!(isomers(SuperMolecule::alcohol(5)), 8);
        assert_eq!(isomers(SuperMolecule::alcohol(6)), 16);
        assert_eq!(isomers(SuperMolecule::alcohol(7)), 33);
        assert_eq!(isomers(SuperMolecule::alcohol(8)), 72);
        // assert_eq!(isomers(SuperMolecule::alcohol(9)), 161);
    }

    #[test]
    fn test_isomers_alkenes() {
        assert_eq!(isomers(SuperMolecule::alkene(3)), 1);
        assert_eq!(isomers(SuperMolecule::alkene(4)), 4);
        assert_eq!(isomers(SuperMolecule::alkene(5)), 9);
        assert_eq!(isomers(SuperMolecule::alkene(6)), 23);
    }

    #[test]
    fn test_isomers_ketones() {
        assert_eq!(isomers(SuperMolecule::ketone(3)), 1);
        assert_eq!(isomers(SuperMolecule::ketone(4)), 3);
        assert_eq!(isomers(SuperMolecule::ketone(5)), 6);
        assert_eq!(isomers(SuperMolecule::ketone(6)), 13);
    }

    #[test]
    fn test_build() {
        let cr = SuperMolecule::alcohol(4).build(vec![1.0, 0.9, 0.8, 0.5, 0.3]);
        let (segments, bonds) = cr.segment_and_bond_count();
        println!("{segments:?}\n{bonds:?}");
        assert_relative_eq!(segments["CH3"], 1.326);
        assert_relative_eq!(segments["CH2"], 0.958);
        assert_relative_eq!(segments[">CH"], 0.216);
        assert_relative_eq!(bonds[&["CH2".to_string(), ">CH".to_string()]], 0.12);
        assert_relative_eq!(bonds[&["CH3".to_string(), ">CH".to_string()]], 0.36);
        assert_relative_eq!(bonds[&["CH2".to_string(), "CH2".to_string()]], 0.28);
        assert_relative_eq!(bonds[&["CH3".to_string(), "CH2".to_string()]], 0.84);
    }

    #[test]
    fn test_alcohol() {
        let cr = SuperMolecule::alcohol(5).build(vec![1.0; 8]);
        let (segments, bonds) = cr.segment_and_bond_count();
        println!("{cr}");
        assert_eq!(segments["CH3"], 4.0);
        assert_eq!(segments["CH2"], 1.0);
        assert_eq!(segments[">CH"], 1.0);
        assert_eq!(segments[">C<"], 1.0);
        assert_eq!(segments["OH"], 1.0);
        assert_relative_eq!(bonds[&["CH3".to_string(), ">C<".to_string()]], 2.0);
        assert_relative_eq!(bonds[&[">CH".to_string(), ">C<".to_string()]], 1.0);
        assert_relative_eq!(bonds[&["CH3".to_string(), ">CH".to_string()]], 1.0);
        assert_relative_eq!(bonds[&["CH2".to_string(), ">CH".to_string()]], 1.0);
        assert_relative_eq!(bonds[&["CH3".to_string(), "CH2".to_string()]], 1.0);
        assert_relative_eq!(bonds[&["OH".to_string(), ">C<".to_string()]], 1.0);
    }

    #[test]
    fn test_ketone() {
        let cr = SuperMolecule::ketone(5).build(vec![1.0; 6]);
        let (segments, bonds) = cr.segment_and_bond_count();
        println!("{cr}");
        assert_eq!(segments["CH3"], 3.0);
        assert_eq!(segments["CH2"], 1.0);
        assert_eq!(segments[">CH"], 1.0);
        assert_eq!(segments[">C=O"], 1.0);
        assert_relative_eq!(bonds[&["CH3".to_string(), ">CH".to_string()]], 1.0);
        assert_relative_eq!(bonds[&["CH2".to_string(), ">CH".to_string()]], 1.0);
        assert_relative_eq!(bonds[&["CH3".to_string(), "CH2".to_string()]], 1.0);
        assert_relative_eq!(bonds[&[">CH".to_string(), ">C=O".to_string()]], 1.0);
        assert_relative_eq!(bonds[&["CH3".to_string(), ">C=O".to_string()]], 1.0);
    }

    #[test]
    fn test_alkene() {
        let cr = SuperMolecule::alkene(5).build(vec![1.0; 7]);
        let (segments, bonds) = cr.segment_and_bond_count();
        println!("{cr}");
        assert_eq!(segments["CH3"], 4.0);
        assert_eq!(segments["CH2"], 1.0);
        assert_eq!(segments[">CH"], 1.0);
        assert_eq!(segments["=C<"], 1.0);
        assert_eq!(segments["=CH"], 1.0);
        assert_relative_eq!(bonds[&["CH3".to_string(), ">CH".to_string()]], 1.0);
        assert_relative_eq!(bonds[&["CH2".to_string(), ">CH".to_string()]], 1.0);
        assert_relative_eq!(bonds[&["CH3".to_string(), "CH2".to_string()]], 1.0);
        assert_relative_eq!(bonds[&[">CH".to_string(), "=C<".to_string()]], 1.0);
        assert_relative_eq!(bonds[&["CH3".to_string(), "=C<".to_string()]], 1.0);
        assert_relative_eq!(bonds[&["CH3".to_string(), "=CH".to_string()]], 1.0);
        assert_relative_eq!(bonds[&["=C<".to_string(), "=CH".to_string()]], 1.0);
    }
}
