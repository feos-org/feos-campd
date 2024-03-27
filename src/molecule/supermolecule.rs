use super::polynomial::{
    Polynomial, Polynomial2,
    Variable::{self as PolyVar, Const, Var},
};
use super::{Constraint, Disjunction, MolecularRepresentation, StructureVariables};
use crate::variables::{ExplicitVariable, Variable};
use feos::core::parameter::{Identifier, SegmentCount};
use indexmap::IndexMap;
use std::borrow::Cow;
use std::collections::HashMap;
use std::iter::once;
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

#[derive(Clone, Debug)]
struct FunctionalGroup {
    groups: Vec<&'static str>,
    smiles: Vec<&'static str>,
    atoms: usize,
}

impl FunctionalGroup {
    pub fn ch3() -> Self {
        Self {
            groups: vec!["CH4", "CH3"],
            smiles: vec!["C"],
            atoms: 1,
        }
    }

    pub fn ctch() -> Self {
        Self {
            groups: vec!["CH#CH", "C≡CH"],
            smiles: vec!["C≡C"],
            atoms: 2,
        }
    }

    pub fn oh() -> Self {
        Self {
            groups: vec!["OH2", "OH"],
            smiles: vec!["O"],
            atoms: 1,
        }
    }

    pub fn nh2() -> Self {
        Self {
            groups: vec!["NH3", "NH2"],
            smiles: vec!["N"],
            atoms: 1,
        }
    }

    pub fn och3() -> Self {
        Self {
            groups: vec!["OHCH3", "OCH3"],
            smiles: vec!["CO"],
            atoms: 2,
        }
    }

    pub fn cdo() -> Self {
        Self {
            groups: vec!["CH2=O", "CH=O", ">C=O"],
            smiles: vec!["C(=O)", ""],
            atoms: 2,
        }
    }

    pub fn cdc() -> Self {
        Self {
            groups: vec!["=CH2", "=CH", "=C<"],
            smiles: vec!["C", "", "=C", ""],
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

    fn feature_variables<I: Iterator<Item = PolyVar>>(
        size: usize,
        variables: usize,
        y: &mut I,
        groups: &mut Polynomial<4>,
        bonds: &mut Polynomial2<4, 4>,
    ) -> (PolyVar, Polynomial<4>) {
        let this = y.next().unwrap();
        let mut indices = [this, Const(0), Const(0), Const(0)];
        let children: Vec<_> = (1..=3)
            .map(|k| (k, (size - 1) / k))
            .filter(|(_, s)| *s > 0)
            .map(|(k, s)| {
                let (var, child) = Self::feature_variables(s, variables, y, groups, bonds);
                indices[k] = var;
                child
            })
            .collect();
        children.into_iter().enumerate().for_each(|(i, child)| {
            let others = match i {
                0 => Polynomial::new(variables, [indices[0], indices[2], indices[3]]).upgrade(),
                1 => Polynomial::new(variables, [indices[0], indices[3]]).upgrade(),
                2 => Polynomial::new(variables, [indices[0]]).upgrade(),
                _ => unreachable!(),
            };
            *bonds += child * others;
        });
        let c_groups = Polynomial::new(variables, indices);
        *groups += &c_groups;
        (this, c_groups)
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
#[derive(Clone)]
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

    pub fn all(size: usize) -> Disjunction<Self, 7> {
        Disjunction([
            SuperMolecule::alkane(size),
            SuperMolecule::alkene(size),
            SuperMolecule::alkyne(size),
            SuperMolecule::alcohol(size),
            SuperMolecule::methylether(size),
            SuperMolecule::ketone(size),
            SuperMolecule::amine(size),
        ])
    }

    pub fn non_associating(size: usize) -> Disjunction<Self, 5> {
        Disjunction([
            SuperMolecule::alkane(size),
            SuperMolecule::alkene(size),
            SuperMolecule::alkyne(size),
            SuperMolecule::methylether(size),
            SuperMolecule::ketone(size),
        ])
    }

    pub fn get_initial_values<const N: usize>(
        molecules: &Disjunction<Self, N>,
        structure: &str,
        y0: &[f64],
    ) -> Vec<f64> {
        let mut y = y0.to_vec();
        let mut c = Vec::new();
        let n_y = molecules
            .0
            .iter()
            .map(|m| m.structure_variables().len())
            .max()
            .unwrap();
        if y0.len() != n_y {
            panic!("Need {n_y} structure variables!");
        }
        for m in &molecules.0 {
            if m.name == structure {
                y.extend_from_slice(y0);
                c.push(1.0);
            } else {
                y.extend_from_slice(&vec![0.0; y0.len()]);
                c.push(0.0);
            }
        }
        [y, c].concat()
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

impl MolecularRepresentation for SuperMolecule {
    fn structure_variables(&self) -> StructureVariables {
        let n = self
            .alkyl_tails()
            .iter()
            .filter(|&&s| s > 0)
            .map(|&s| SuperAlkyl::variables(s))
            .sum::<usize>();
        (0..n).map(|_| Variable::binary()).collect()
    }

    fn feature_variables(
        &self,
        index_structure_vars: &[i32],
    ) -> IndexMap<String, ExplicitVariable> {
        let variables = index_structure_vars[0] as usize + index_structure_vars.len();
        let y = &mut once(Const(1)).chain(index_structure_vars.iter().skip(1).map(|&y| Var(y)));
        let mut c_groups = Polynomial::zero(variables);
        let mut c_bonds = Polynomial2::zero(variables);
        let mut f_group_indices = Vec::new();
        let mut children = Vec::new();

        self.alkyl_tails()
            .into_iter()
            .enumerate()
            .filter(|&(_, s)| s > 0)
            .for_each(|(k, s)| {
                if k % 2 == 0 {
                    f_group_indices.push([Const(1), Const(0), Const(0)]);
                    children.push(Vec::new());
                }
                let (var, child) =
                    SuperAlkyl::feature_variables(s, variables, y, &mut c_groups, &mut c_bonds);
                f_group_indices.last_mut().unwrap()[k % 2 + 1] = var;
                children.last_mut().unwrap().push(child);
            });
        let f_groups: Vec<_> = f_group_indices
            .iter()
            .map(|&f| Polynomial::new(variables, f))
            .collect();

        let ff_bonds = (f_groups.len() == 2).then(|| f_groups[0].clone() * f_groups[1].clone());
        let f_group: Polynomial<3> = f_groups.into_iter().sum();

        let mut coefs = IndexMap::new();
        ["CH3", "CH2", ">CH", ">C<"]
            .into_iter()
            .chain(self.functional_group.groups.iter().copied())
            .map(String::from)
            .zip(c_groups.0.into_iter().chain(f_group.0))
            .for_each(|(k, (l, c))| {
                let (lin, con) = coefs.entry(k).or_insert((vec![0; variables], 0));
                lin.iter_mut().zip(l).for_each(|(lin, l)| *lin += l);
                *con += c;
            });

        let mut features: IndexMap<_, _> = coefs
            .into_iter()
            .map(|(k, (lin, con))| {
                let mut vars = Vec::new();
                let mut coefs = Vec::new();
                lin.into_iter().enumerate().for_each(|(k, v)| {
                    if v != 0 {
                        vars.push(k as i32);
                        coefs.push(v as f64);
                    }
                });
                (
                    k.clone(),
                    ExplicitVariable::new(k)
                        .linear_struct(vars, coefs)
                        .cons(con as f64),
                )
            })
            .collect();

        let mut f_bonds: Polynomial2<4, 3> = Polynomial2::zero(variables);
        children
            .into_iter()
            .zip(f_group_indices)
            .for_each(|(children, indices)| {
                children.into_iter().enumerate().for_each(|(i, child)| {
                    let others = match i {
                        0 => Polynomial::new(variables, [indices[0], indices[2]]).upgrade(),
                        1 => Polynomial::new(variables, [indices[0]]).upgrade(),
                        _ => unreachable!(),
                    };
                    f_bonds += child * others;
                })
            });

        let mut bonds = IndexMap::new();
        let c_groups = ["CH3", "CH2", ">CH", ">C<"];
        let f_groups = &self.functional_group.groups;
        fill_bonds(variables, &c_groups, &c_groups, c_bonds, &mut bonds);
        fill_bonds(variables, &c_groups, f_groups, f_bonds, &mut bonds);
        if let Some(ff_bonds) = ff_bonds {
            fill_bonds(variables, f_groups, f_groups, ff_bonds, &mut bonds);
        }

        features.extend(bonds.into_iter().map(|(k, Polynomial2([[b]]))| {
            let (qvars, qcoefs): (Vec<_>, Vec<_>) =
                b.0.into_iter()
                    .enumerate()
                    .flat_map(|(i1, b)| {
                        b.into_iter().enumerate().filter_map(move |(i2, b)| {
                            (b != 0).then_some(((i1 as i32, i2 as i32), b as f64))
                        })
                    })
                    .unzip();
            let (qvars1, qvars2) = qvars.into_iter().unzip();

            let (lvars, lcoefs): (Vec<_>, Vec<_>) =
                b.1.into_iter()
                    .enumerate()
                    .filter_map(|(i, b)| (b != 0).then_some((i as i32, b as f64)))
                    .unzip();

            let bond = ExplicitVariable::new(k.clone())
                .quadratic_struct(qvars1, qvars2, qcoefs)
                .linear_struct(lvars, lcoefs)
                .cons(b.2 as f64);
            (k, bond)
        }));

        features.retain(|_, v| !(v.lvars.is_empty() && v.qvars1.is_empty() && v.cons == 0.0));

        features
    }

    fn constraints(&self, index_structure_vars: &[i32]) -> Vec<Constraint> {
        // maximum size
        let mut constraints = Vec::new();
        let (coefs, size) = self.size_constraint();
        constraints.push(
            Constraint::new()
                .linear_struct(index_structure_vars.to_vec(), coefs)
                .upbnd(size),
        );

        // minimum size
        constraints.push(
            Constraint::new()
                .linear_struct(vec![index_structure_vars[0]], vec![1.0])
                .eqbnd(1.0),
        );

        // bond constraints
        for bond in self.bond_constraints(index_structure_vars[0]) {
            constraints.push(
                Constraint::new()
                    .linear_struct(bond.to_vec(), vec![1.0, -1.0])
                    .lobnd(0.0),
            );
        }

        // symmetry constraints
        for (index_structure_vars, coefs) in self.symmetry_constraints(index_structure_vars[0]) {
            constraints.push(
                Constraint::new()
                    .linear_struct(index_structure_vars.to_vec(), coefs.to_vec())
                    .lobnd(0.0),
            );
        }
        constraints
    }

    fn smiles(&self, y: &[usize]) -> String {
        let y_iter = &mut y.iter().copied();
        self.alkyl_tails()
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
            })
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

fn fill_bonds<const M: usize, const N: usize>(
    variables: usize,
    seg1: &[&str],
    seg2: &[&str],
    poly: Polynomial2<M, N>,
    bonds: &mut IndexMap<String, Polynomial2<1, 1>>,
) {
    seg1.iter().zip(poly.0).for_each(|(g2, b)| {
        seg2.iter().zip(b).for_each(|(g1, b)| {
            let mut g = [g1, g2];
            g.sort();
            let [g1, g2] = g;
            let k = format!("{g1}-{g2}");
            *bonds.entry(k).or_insert(Polynomial2::zero(variables)) += Polynomial2([[b]]);
        })
    });
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
        let molecules = SuperMolecule::all(3);
        let total_disjunct = molecules.generate_solutions().len();
        let mut total = 0;
        for molecule in molecules.0 {
            let solutions = molecule.generate_solutions().len();
            println!("{} {}", molecule.name, solutions);
            total += solutions;
        }
        assert_eq!(total_disjunct, total);
    }

    #[test]
    fn test_build() {
        let supermolecule = SuperMolecule::alcohol(4);
        let features = supermolecule.evaluate_feature_variables(&[1.0, 0.8, 0.5, 0.3]);
        println!("{features:?}");

        assert_relative_eq!(features["CH3"], 1.3);
        assert_relative_eq!(features["CH2"], 1.0);
        assert_relative_eq!(features[">CH"], 0.3);
        assert_relative_eq!(features["OH"], 1.0);
        assert_relative_eq!(features[">CH-CH2"], 0.15);
        assert_relative_eq!(features[">CH-CH3"], 0.39);
        assert_relative_eq!(features["CH2-CH2"], 0.35);
        assert_relative_eq!(features["CH2-CH3"], 0.61);
        assert_relative_eq!(features["CH3-OH"], 0.2);
        assert_relative_eq!(features["CH2-OH"], 0.5);
        assert_relative_eq!(features[">CH-OH"], 0.3);
    }

    #[test]
    fn test_alkane() {
        let supermolecule = SuperMolecule::alkane(5);
        let variables = supermolecule.structure_variables().len();
        let features = supermolecule.evaluate_feature_variables(&vec![1.0; variables]);
        println!("{features:?}");

        assert_eq!(features["CH3"], 5.0);
        assert_eq!(features["CH2"], 1.0);
        assert_eq!(features[">CH"], 1.0);
        assert_eq!(features[">C<"], 1.0);
        assert_relative_eq!(features[">C<->CH"], 1.0);
        assert_relative_eq!(features[">CH-CH3"], 1.0);
        assert_relative_eq!(features["CH2-CH3"], 1.0);
        assert_relative_eq!(features[">C<-CH3"], 3.0);
    }

    #[test]
    fn test_alcohol() {
        let supermolecule = SuperMolecule::alcohol(5);
        let variables = supermolecule.structure_variables().len();
        let features = supermolecule.evaluate_feature_variables(&vec![1.0; variables]);
        println!("{features:?}");
        assert_eq!(features["CH3"], 4.0);
        assert_eq!(features["CH2"], 1.0);
        assert_eq!(features[">CH"], 1.0);
        assert_eq!(features[">C<"], 1.0);
        assert_eq!(features["OH"], 1.0);
        assert_relative_eq!(features[">C<-CH3"], 2.0);
        assert_relative_eq!(features[">C<->CH"], 1.0);
        assert_relative_eq!(features[">CH-CH3"], 1.0);
        assert_relative_eq!(features[">CH-CH2"], 1.0);
        assert_relative_eq!(features["CH2-CH3"], 1.0);
        assert_relative_eq!(features[">C<-OH"], 1.0);
    }

    #[test]
    fn test_ketone() {
        let supermolecule = SuperMolecule::ketone(5);
        let variables = supermolecule.structure_variables().len();
        let features = supermolecule.evaluate_feature_variables(&vec![1.0; variables]);
        println!("{features:?}");
        assert_eq!(features["CH3"], 3.0);
        assert_eq!(features["CH2"], 1.0);
        assert_eq!(features[">CH"], 1.0);
        assert_eq!(features[">C=O"], 1.0);
        assert_relative_eq!(features[">CH-CH3"], 1.0);
        assert_relative_eq!(features[">CH-CH2"], 1.0);
        assert_relative_eq!(features["CH2-CH3"], 1.0);
        assert_relative_eq!(features[">C=O->CH"], 1.0);
        assert_relative_eq!(features[">C=O-CH3"], 1.0);
    }

    #[test]
    fn test_alkene() {
        let supermolecule = SuperMolecule::alkene(5);
        let variables = supermolecule.structure_variables().len();
        let features = supermolecule.evaluate_feature_variables(&vec![1.0; variables]);
        println!("{features:?}");
        assert_eq!(features["CH3"], 4.0);
        assert_eq!(features["CH2"], 1.0);
        assert_eq!(features[">CH"], 1.0);
        assert_eq!(features["=C<"], 1.0);
        assert_eq!(features["=CH"], 1.0);
        assert_relative_eq!(features[">CH-CH3"], 1.0);
        assert_relative_eq!(features[">CH-CH2"], 1.0);
        assert_relative_eq!(features["CH2-CH3"], 1.0);
        assert_relative_eq!(features["=C<->CH"], 1.0);
        assert_relative_eq!(features["=C<-CH3"], 1.0);
        assert_relative_eq!(features["=CH-CH3"], 1.0);
        assert_relative_eq!(features["=C<-=CH"], 1.0);
    }

    fn molecule_disjunct() -> Disjunction<SuperMolecule, 3> {
        Disjunction([
            SuperMolecule::alcohol(5),
            SuperMolecule::methylether(5),
            SuperMolecule::ketone(5),
        ])
    }

    #[test]
    fn test_ethanol() {
        let molecules = molecule_disjunct();
        let y = SuperMolecule::get_initial_values(
            &molecules,
            "alcohol",
            &[1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        );
        let features = molecule_disjunct().evaluate_feature_variables(&y);
        println!("{features:?}");
        assert_eq!(features["CH3"], 1.0);
        assert_eq!(features["CH2"], 1.0);
        assert_eq!(features["OH"], 1.0);
        assert_relative_eq!(features["CH2-CH3"], 1.0);
        assert_relative_eq!(features["CH2-OH"], 1.0);
    }

    #[test]
    fn test_2_methoxybutane() {
        let molecules = molecule_disjunct();
        let y = SuperMolecule::get_initial_values(
            &molecules,
            "methylether",
            &[1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        );
        let features = molecule_disjunct().evaluate_feature_variables(&y);
        println!("{features:?}");
        assert_eq!(features["CH3"], 2.0);
        assert_eq!(features[">CH"], 1.0);
        assert_eq!(features["OCH3"], 1.0);
        assert_relative_eq!(features[">CH-CH3"], 2.0);
        assert_relative_eq!(features[">CH-OCH3"], 1.0);
    }

    #[test]
    fn test_acetone() {
        let molecules = molecule_disjunct();
        let y = SuperMolecule::get_initial_values(
            &molecules,
            "ketone",
            &[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
        );
        let features = molecule_disjunct().evaluate_feature_variables(&y);
        println!("{features:?}");
        assert_eq!(features["CH3"], 2.0);
        assert_eq!(features[">C=O"], 1.0);
        assert_relative_eq!(features[">C=O-CH3"], 2.0);
    }
}
