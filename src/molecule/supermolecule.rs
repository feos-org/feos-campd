use super::polynomial::{Polynomial, Polynomial2};
use super::{
    ContinuousVariables, DiscreteVariables, LinearConstraint, MolecularRepresentation, Variable,
};
use feos::core::parameter::{Identifier, SegmentCount};
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::collections::HashMap;
use std::ops::{Add, Mul};

/// The molecular features used in the [GcPcSaftPropertyModel](crate::GcPcSaftPropertyModel).
#[derive(Clone)]
pub struct SegmentAndBondCount {
    pub segments: HashMap<String, f64>,
    pub bonds: HashMap<[String; 2], f64>,
}

impl SegmentAndBondCount {
    pub fn new(segments: HashMap<String, f64>, bonds: HashMap<[String; 2], f64>) -> Self {
        Self { segments, bonds }
    }
}

impl SegmentCount for SegmentAndBondCount {
    type Count = f64;

    fn identifier(&self) -> Cow<Identifier> {
        Cow::Owned(Identifier::default())
    }

    fn segment_count(&self) -> Cow<HashMap<String, Self::Count>> {
        Cow::Borrowed(&self.segments)
    }
}

impl Add for SegmentAndBondCount {
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self::Output {
        rhs.segments
            .into_iter()
            .for_each(|(s, n)| *self.segments.entry(s).or_insert(0.0) += n);
        rhs.bonds
            .into_iter()
            .for_each(|(b, n)| *self.bonds.entry(b).or_insert(0.0) += n);
        self
    }
}

impl Mul<f64> for SegmentAndBondCount {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        let segments = self
            .segments
            .into_iter()
            .map(|(s, n)| (s, n * rhs))
            .collect();
        let bonds = self.bonds.into_iter().map(|(b, n)| (b, n * rhs)).collect();
        SegmentAndBondCount::new(segments, bonds)
    }
}

#[derive(Clone, Serialize, Deserialize, Debug)]
struct FunctionalGroup {
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

struct SuperAlkyl;

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
                constraints.push([this, *index + 1]);
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

/// The molecule superstructure representation.
#[derive(Clone, Serialize, Deserialize)]
pub struct SuperMolecule {
    pub name: String,
    pub size: usize,
    functional_group: FunctionalGroup,
    alkyls: Vec<usize>,
    symmetries: Vec<[usize; 2]>,
}

impl SuperMolecule {
    pub fn alkane(size: usize) -> Self {
        Self {
            name: "alkane".into(),
            size,
            functional_group: FunctionalGroup::ch3(),
            alkyls: vec![1],
            symmetries: vec![],
        }
    }

    pub fn alkene(size: usize) -> Self {
        Self {
            name: "alkene".into(),
            size,
            functional_group: FunctionalGroup::cdc(),
            alkyls: vec![1, 2, 2, 4],
            symmetries: vec![[0, 1], [0, 2], [2, 3]],
        }
    }

    pub fn alkyne(size: usize) -> Self {
        Self {
            name: "alkyne".into(),
            size,
            functional_group: FunctionalGroup::ctch(),
            alkyls: vec![1],
            symmetries: vec![],
        }
    }

    pub fn alcohol(size: usize) -> Self {
        Self {
            name: "alcohol".into(),
            size,
            functional_group: FunctionalGroup::oh(),
            alkyls: vec![1],
            symmetries: vec![],
        }
    }

    pub fn methylether(size: usize) -> Self {
        Self {
            name: "methylether".into(),
            size,
            functional_group: FunctionalGroup::och3(),
            alkyls: vec![1],
            symmetries: vec![],
        }
    }

    pub fn ketone(size: usize) -> Self {
        Self {
            name: "ketone".into(),
            size,
            functional_group: FunctionalGroup::cdo(),
            alkyls: vec![1, 2],
            symmetries: vec![[0, 1]],
        }
    }

    pub fn amine(size: usize) -> Self {
        Self {
            name: "amine".into(),
            size,
            functional_group: FunctionalGroup::nh2(),
            alkyls: vec![1],
            symmetries: vec![],
        }
    }

    pub fn all(size: usize) -> [Self; 7] {
        [
            SuperMolecule::alkane(size),
            SuperMolecule::alkene(size),
            SuperMolecule::alkyne(size),
            SuperMolecule::alcohol(size),
            SuperMolecule::methylether(size),
            SuperMolecule::ketone(size),
            SuperMolecule::amine(size),
        ]
    }

    pub fn non_associating(size: usize) -> [Self; 5] {
        [
            SuperMolecule::alkane(size),
            SuperMolecule::alkene(size),
            SuperMolecule::alkyne(size),
            SuperMolecule::methylether(size),
            SuperMolecule::ketone(size),
        ]
    }
}

impl SuperMolecule {
    fn alkyl_tails(&self) -> Vec<usize> {
        self.alkyls
            .iter()
            .map(|a| (self.size - self.functional_group.atoms) / a)
            .collect()
    }

    pub fn size_constraint(&self) -> (Vec<f64>, f64) {
        (
            vec![1.0; self.structure_variables().len()],
            (self.size - self.functional_group.atoms) as f64,
        )
    }

    pub fn bond_constraints(&self, mut index: i32) -> Vec<[i32; 2]> {
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

    pub fn symmetry_constraints(&self, mut index: i32) -> Vec<(Vec<i32>, Vec<f64>)> {
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
}

impl MolecularRepresentation<1> for SuperMolecule {
    type ChemicalRecord = SegmentAndBondCount;

    fn structure_variables(&self) -> DiscreteVariables {
        let n = self
            .alkyl_tails()
            .iter()
            .filter(|&&s| s > 0)
            .map(|&s| SuperAlkyl::variables(s))
            .sum::<usize>();
        vec![Variable::binary(); n].into()
    }

    fn parameter_variables(&self) -> ContinuousVariables {
        vec![].into()
    }

    fn determine_parameters(&self, _: &[f64]) -> Vec<f64> {
        vec![]
    }

    fn constraints(&self, index_vars: &[i32], _: Option<&[i32]>) -> Vec<LinearConstraint> {
        // maximum size
        let mut constraints = Vec::new();
        let (coefs, size) = self.size_constraint();
        constraints.push(LinearConstraint::new(index_vars.to_vec(), coefs).upbnd(size));

        // minimum size
        constraints.push(LinearConstraint::new(vec![index_vars[0]], vec![1.0]).eqbnd(1.0));

        // bond constraints
        for bond in self.bond_constraints(index_vars[0]) {
            constraints.push(LinearConstraint::new(bond.to_vec(), vec![1.0, -1.0]).lobnd(0.0));
        }

        // symmetry constraints
        for (index_vars, coefs) in self.symmetry_constraints(index_vars[0]) {
            constraints.push(LinearConstraint::new(index_vars.to_vec(), coefs.to_vec()).lobnd(0.0));
        }
        constraints
    }

    fn build(&self, y: &[f64], _: &[f64]) -> [SegmentAndBondCount; 1] {
        match self.alkyls.len() {
            1 => (),
            2 => return [self.build_ketone(y)],
            4 => return [self.build_alkene(y)],
            _ => unreachable!(),
        }
        let y_iter = &mut y.iter().copied();

        let mut c_segments = Polynomial::zero();
        let mut c_bonds = Polynomial2::zero();
        let s = SuperAlkyl::build(
            self.size - self.functional_group.atoms,
            y_iter,
            &mut c_segments,
            &mut c_bonds,
        );

        let cs = ["CH3", "CH2", ">CH", ">C<"];

        let mut segments: HashMap<_, _> = cs
            .iter()
            .map(|&s| s.to_string())
            .zip(c_segments.iter().cloned())
            .filter(|(_, g)| *g > 0.0)
            .collect();
        *segments
            .entry(self.functional_group.groups[0].clone())
            .or_insert(0.0) += 1.0;

        let mut bonds = HashMap::new();
        fill_bond_map(cs, cs, &c_bonds, &mut bonds);
        s.iter().zip(cs.iter()).for_each(|(g, c)| {
            if *g > 0.0 {
                *bonds
                    .entry([self.functional_group.groups[0].clone(), c.to_string()])
                    .or_insert(0.0) += *g;
            }
        });

        [SegmentAndBondCount::new(segments, bonds)]
    }

    fn smiles(&self, y: &[usize]) -> [String; 1] {
        let y_iter = &mut y.iter().copied();
        [self
            .alkyl_tails()
            .into_iter()
            .zip(self.functional_group.smiles.iter())
            .fold(String::new(), |mut acc, (s, sm)| {
                acc.push_str(sm);
                acc.push_str(
                    &SuperAlkyl::smiles(s, y_iter)
                        .into_iter()
                        .collect::<String>(),
                );
                acc
            })]
    }
}

impl SuperMolecule {
    fn build_ketone(&self, y: &[f64]) -> SegmentAndBondCount {
        let y_iter = &mut y.iter().copied();

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

        SegmentAndBondCount::new(segments, bonds)
    }

    fn build_alkene(&self, y: &[f64]) -> SegmentAndBondCount {
        let y_iter = &mut y.iter().copied();

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
            .filter(|(_, g)| *g > 0.0)
            .collect();

        let mut bonds = HashMap::new();
        fill_bond_map(cs, cs, &c_bonds, &mut bonds);
        fill_bond_map(cs, cds, &cd_bonds, &mut bonds);
        let dd_bonds = cd_segments[0].outer_product(cd_segments[1]);
        fill_bond_map(cds, cds, &dd_bonds, &mut bonds);

        SegmentAndBondCount::new(segments, bonds)
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
                .filter(|&(j, _)| (j != i))
                .map(|(_, cj)| Polynomial::new(cj.sum()))
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
    use crate::MolecularRepresentation;
    use approx::assert_relative_eq;

    // https://pubs.acs.org/doi/pdf/10.1021/ja01359a027
    #[test]
    fn test_isomers_alcohols() {
        assert_eq!(SuperMolecule::alcohol(4).generate_solutions().len(), 4);
        assert_eq!(SuperMolecule::alcohol(5).generate_solutions().len(), 8);
        assert_eq!(SuperMolecule::alcohol(6).generate_solutions().len(), 16);
        assert_eq!(SuperMolecule::alcohol(7).generate_solutions().len(), 33);
        assert_eq!(SuperMolecule::alcohol(8).generate_solutions().len(), 72);
        // assert_eq!(SuperMolecule::alcohol(9).generate_solutions().len(), 161);
    }

    #[test]
    fn test_isomers_alkenes() {
        assert_eq!(SuperMolecule::alkene(3).generate_solutions().len(), 1);
        assert_eq!(SuperMolecule::alkene(4).generate_solutions().len(), 4);
        assert_eq!(SuperMolecule::alkene(5).generate_solutions().len(), 9);
        assert_eq!(SuperMolecule::alkene(6).generate_solutions().len(), 23);
    }

    #[test]
    fn test_isomers_ketones() {
        assert_eq!(SuperMolecule::ketone(3).generate_solutions().len(), 1);
        assert_eq!(SuperMolecule::ketone(4).generate_solutions().len(), 3);
        assert_eq!(SuperMolecule::ketone(5).generate_solutions().len(), 6);
        assert_eq!(SuperMolecule::ketone(6).generate_solutions().len(), 13);
    }

    #[test]
    fn test_isomers_disjunction() {
        let molecules = SuperMolecule::all(5);
        let total_disjunct = molecules.generate_solutions().len();
        let mut total = 0;
        for molecule in molecules {
            let solutions = molecule.generate_solutions().len();
            println!("{} {}", molecule.name, solutions);
            total += solutions;
        }
        assert_eq!(total_disjunct, total);
    }

    #[test]
    fn test_build() {
        let [SegmentAndBondCount { segments, bonds }] =
            SuperMolecule::alcohol(4).build(&[0.9, 0.8, 0.5, 0.3], &[]);
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
        let supermolecule = SuperMolecule::alcohol(5);
        let [SegmentAndBondCount { segments, bonds }] =
            supermolecule.build(&vec![1.0; supermolecule.structure_variables().len()], &[]);
        println!("{segments:?}\n{bonds:?}");
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
        let supermolecule = SuperMolecule::ketone(5);
        let [SegmentAndBondCount { segments, bonds }] =
            supermolecule.build(&vec![1.0; supermolecule.structure_variables().len()], &[]);
        println!("{segments:?}\n{bonds:?}");
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
        let supermolecule = SuperMolecule::alkene(5);
        let [SegmentAndBondCount { segments, bonds }] =
            supermolecule.build(&vec![1.0; supermolecule.structure_variables().len()], &[]);
        println!("{segments:?}\n{bonds:?}");
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

    fn molecule_disjunct() -> [SuperMolecule; 3] {
        [
            SuperMolecule::alcohol(5),
            SuperMolecule::methylether(5),
            SuperMolecule::ketone(5),
        ]
    }

    #[test]
    fn test_build_disjunct() {
        println!("{}", molecule_disjunct().structure_variables().len());
        let [SegmentAndBondCount { segments, bonds }] =
            molecule_disjunct().build(&[1.0, 0.8, 0.6, 0.1, 0.3, 0.2, 0.1, 0.5, 0.3, 0.2], &[]);
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
            molecule_disjunct().build(&[1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0], &[]);
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
            molecule_disjunct().build(&[1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0], &[]);
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
            molecule_disjunct().build(&[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0], &[]);
        println!("{segments:?}\n{bonds:?}");
        assert_eq!(segments["CH3"], 2.0);
        assert_eq!(segments[">C=O"], 1.0);
        assert_relative_eq!(bonds[&["CH3".to_string(), ">C=O".to_string()]], 2.0);
    }
}
