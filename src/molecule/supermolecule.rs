use super::polynomial::{Polynomial, Polynomial2};
use super::{Disjunction, MolecularRepresentation};
use crate::{ChemicalRecord, GeneralConstraint};
use good_lp::{constraint, Constraint, Expression, Variable};
use num_dual::DualNum;
use std::collections::HashMap;
use std::vec;

#[derive(Clone, Copy)]
enum FunctionalGroup {
    Simple(&'static str, &'static str),
    Ketone([&'static str; 2], [&'static str; 2]),
    Alkene([&'static str; 3], [&'static str; 4]),
}

impl FunctionalGroup {
    const fn ch3() -> Self {
        Self::Simple("CH3", "C")
    }

    const fn ctch() -> Self {
        Self::Simple("C≡CH", "C#C")
    }

    const fn oh() -> Self {
        Self::Simple("OH", "O")
    }

    const fn nh2() -> Self {
        Self::Simple("NH2", "N")
    }

    const fn och3() -> Self {
        Self::Simple("OCH3", "CO")
    }

    const fn cdo() -> Self {
        Self::Ketone([">C=O", "CH=O"], ["C(=O)", ""])
    }

    const fn cdc() -> Self {
        Self::Alkene(["=C<", "=CH", "=CH2"], ["C", "", "=C", ""])
    }
}

impl FunctionalGroup {
    fn smiles(&self) -> Vec<&'static str> {
        match self {
            FunctionalGroup::Simple(_, s) => vec![s],
            FunctionalGroup::Ketone(_, s) => s.to_vec(),
            FunctionalGroup::Alkene(_, s) => s.to_vec(),
        }
    }

    fn base_group(&self) -> &'static str {
        match self {
            FunctionalGroup::Simple(g, _) => g,
            FunctionalGroup::Ketone(g, _) => g[0],
            FunctionalGroup::Alkene(g, _) => g[0],
        }
    }

    fn symmetries(&self) -> Vec<[usize; 2]> {
        match self {
            FunctionalGroup::Simple(_, _) => vec![],
            FunctionalGroup::Ketone(_, _) => vec![[0, 1]],
            FunctionalGroup::Alkene(_, _) => vec![[0, 1], [0, 2], [2, 3]],
        }
    }

    fn alkyls(&self) -> Vec<usize> {
        match self {
            FunctionalGroup::Simple(_, _) => vec![1],
            FunctionalGroup::Ketone(_, _) => vec![1, 2],
            FunctionalGroup::Alkene(_, _) => vec![1, 2, 2, 4],
        }
    }
}

struct SuperAlkyl;

impl SuperAlkyl {
    const fn variables(size: usize) -> usize {
        let mut res = 1;
        if size - 1 > 0 {
            res += Self::variables(size - 1)
        }
        if (size - 1) / 2 > 0 {
            res += Self::variables((size - 1) / 2)
        }
        if (size - 1) / 3 > 0 {
            res += Self::variables((size - 1) / 3)
        }
        res
    }

    fn bond_constraints(size: usize, index: &mut usize, constraints: &mut Vec<[usize; 2]>) {
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
        index: &mut usize,
        level: f64,
        constraints: &mut Vec<(Vec<usize>, Vec<f64>)>,
    ) -> (Vec<usize>, Vec<f64>) {
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

    fn build<I: Iterator<Item = D>, D: DualNum<f64> + Copy>(
        size: usize,
        y: &mut I,
        groups: &mut Polynomial<D, 4>,
        bonds: &mut Polynomial2<D, 4, 4>,
    ) -> Polynomial<D, 4> {
        let value = y.next().unwrap();

        let children: Vec<_> = (1..=3)
            .map(|k| (size - 1) / k)
            .filter(|s| *s > 0)
            .map(|s| Self::build(s, y, groups, bonds))
            .collect();
        let res = calculate_bonds(children, bonds) * value;
        *groups += res;

        res
    }

    fn smiles<I: Iterator<Item = f64>>(size: usize, y: &mut I) -> Vec<char> {
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
        if this == 1.0 {
            smiles
        } else {
            vec![]
        }
    }
}

/// The molecule superstructure representation.
#[derive(Clone, Copy)]
pub struct SuperMolecule {
    pub name: &'static str,
    min_size: usize,
    pub size: usize,
    functional_group: FunctionalGroup,
    functional_group_atoms: usize,
}

impl SuperMolecule {
    pub const fn alkane(min_size: usize, size: usize) -> Self {
        Self {
            name: "alkane",
            min_size,
            size,
            functional_group: FunctionalGroup::ch3(),
            functional_group_atoms: 1,
        }
    }

    pub const fn alkene(min_size: usize, size: usize) -> Self {
        Self {
            name: "alkene",
            min_size,
            size,
            functional_group: FunctionalGroup::cdc(),
            functional_group_atoms: 2,
        }
    }

    pub const fn alkyne(min_size: usize, size: usize) -> Self {
        Self {
            name: "alkyne",
            min_size,
            size,
            functional_group: FunctionalGroup::ctch(),
            functional_group_atoms: 2,
        }
    }

    pub const fn alcohol(min_size: usize, size: usize) -> Self {
        Self {
            name: "alcohol",
            min_size,
            size,
            functional_group: FunctionalGroup::oh(),
            functional_group_atoms: 1,
        }
    }

    pub const fn methylether(min_size: usize, size: usize) -> Self {
        Self {
            name: "methylether",
            min_size,
            size,
            functional_group: FunctionalGroup::och3(),
            functional_group_atoms: 2,
        }
    }

    pub const fn ketone(min_size: usize, size: usize) -> Self {
        Self {
            name: "ketone",
            min_size,
            size,
            functional_group: FunctionalGroup::cdo(),
            functional_group_atoms: 2,
        }
    }

    pub const fn amine(min_size: usize, size: usize) -> Self {
        Self {
            name: "amine",
            min_size,
            size,
            functional_group: FunctionalGroup::nh2(),
            functional_group_atoms: 1,
        }
    }

    pub const fn all(min_size: usize, size: usize) -> Disjunction<7> {
        Disjunction([
            SuperMolecule::alkane(min_size, size),
            SuperMolecule::alkene(min_size, size),
            SuperMolecule::alkyne(min_size, size),
            SuperMolecule::alcohol(min_size, size),
            SuperMolecule::methylether(min_size, size),
            SuperMolecule::ketone(min_size, size),
            SuperMolecule::amine(min_size, size),
        ])
    }

    pub const fn non_associating(min_size: usize, size: usize) -> Disjunction<5> {
        Disjunction([
            SuperMolecule::alkane(min_size, size),
            SuperMolecule::alkene(min_size, size),
            SuperMolecule::alkyne(min_size, size),
            SuperMolecule::methylether(min_size, size),
            SuperMolecule::ketone(min_size, size),
        ])
    }

    pub fn get_initial_values<const N: usize>(
        molecules: &Disjunction<N>,
        structure: &str,
        y0: &[f64],
    ) -> Vec<f64> {
        let mut c = Vec::new();
        let n_y = molecules.0.iter().map(|m| m.variables()).max().unwrap();
        if y0.len() != n_y {
            panic!("Need {n_y} structure variables!");
        }
        for m in &molecules.0 {
            if m.name == structure {
                c.push(1.0);
            } else {
                c.push(0.0);
            }
        }
        [y0, &c].concat()
    }
}

impl SuperMolecule {
    pub const fn variables(&self) -> usize {
        match self.functional_group {
            FunctionalGroup::Simple(_, _) => {
                SuperAlkyl::variables(self.size - self.functional_group_atoms)
            }
            FunctionalGroup::Ketone(_, _) => {
                let mut variables = 0;
                if self.size - self.functional_group_atoms > 0 {
                    variables += SuperAlkyl::variables(self.size - self.functional_group_atoms);
                }
                if (self.size - self.functional_group_atoms) / 2 > 0 {
                    variables +=
                        SuperAlkyl::variables((self.size - self.functional_group_atoms) / 2);
                }
                variables
            }
            FunctionalGroup::Alkene(_, _) => {
                let mut variables = 0;
                if self.size - self.functional_group_atoms > 0 {
                    variables += SuperAlkyl::variables(self.size - self.functional_group_atoms);
                }
                if (self.size - self.functional_group_atoms) / 2 > 0 {
                    variables +=
                        2 * SuperAlkyl::variables((self.size - self.functional_group_atoms) / 2);
                }
                if (self.size - self.functional_group_atoms) / 4 > 0 {
                    variables +=
                        SuperAlkyl::variables((self.size - self.functional_group_atoms) / 4);
                }
                variables
            }
        }
    }

    fn alkyl_tails(&self) -> Vec<usize> {
        self.functional_group
            .alkyls()
            .iter()
            .map(|a| (self.size - self.functional_group_atoms) / a)
            .collect()
    }

    pub fn size_constraint(&self) -> (f64, f64) {
        (
            self.min_size as f64 - self.functional_group_atoms as f64,
            self.size as f64 - self.functional_group_atoms as f64,
        )
    }

    pub fn bond_constraints(&self, mut index: usize) -> Vec<[usize; 2]> {
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

    pub fn symmetry_constraints(&self, mut index: usize) -> Vec<(Vec<usize>, Vec<f64>)> {
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
        self.functional_group
            .symmetries()
            .iter()
            .for_each(|&[a, b]| {
                if children.len() > b {
                    constraints.push(merge_symmetry_constraints(&children[a], &children[b]))
                }
            });
        constraints
    }

    pub fn constraints(&self, s: &[Variable]) -> Vec<(Vec<Variable>, Vec<f64>, GeneralConstraint)> {
        // size
        let mut constraints = Vec::new();
        let (min_size, size) = self.size_constraint();
        constraints.push((
            s.to_vec(),
            vec![1.0; s.len()],
            GeneralConstraint::Inequality(Some(min_size), Some(size)),
        ));

        // minimum size
        constraints.push((vec![s[0]], vec![1.0], GeneralConstraint::Equality(1.0)));

        // bond constraints
        for [a, b] in self.bond_constraints(0) {
            constraints.push((
                vec![s[a], s[b]],
                vec![1.0, -1.0],
                GeneralConstraint::Inequality(Some(0.0), None),
            ));
        }

        // symmetry constraints
        for (index_structure_vars, coefs) in self.symmetry_constraints(0) {
            let vars = index_structure_vars.into_iter().map(|i| s[i]).collect();
            constraints.push((vars, coefs, GeneralConstraint::Inequality(Some(0.0), None)));
        }

        constraints
    }

    pub fn build<D: DualNum<f64> + Copy>(&self, y: Vec<D>) -> ChemicalRecord<D> {
        match self.functional_group {
            FunctionalGroup::Simple(_, _) => (),
            FunctionalGroup::Ketone(_, _) => return self.build_ketone(y),
            FunctionalGroup::Alkene(_, _) => return self.build_alkene(y),
        }
        let y_iter = &mut y.into_iter();

        let mut c_groups = Polynomial::zero();
        let mut c_bonds = Polynomial2::zero();
        let s = SuperAlkyl::build(
            self.size - self.functional_group_atoms,
            y_iter,
            &mut c_groups,
            &mut c_bonds,
        );

        let cs = ["CH3", "CH2", ">CH", ">C<"];

        let mut groups: HashMap<_, _> = cs.into_iter().zip(c_groups.iter().cloned()).collect();
        *groups
            .entry(self.functional_group.base_group())
            .or_insert(D::zero()) += D::one();

        let mut bonds = HashMap::new();
        fill_bond_map(cs, cs, &c_bonds, &mut bonds);
        s.iter().zip(cs.iter()).for_each(|(&g, c)| {
            *bonds
                .entry([self.functional_group.base_group(), c])
                .or_insert(D::zero()) += g
        });

        ChemicalRecord::new(groups, bonds)
    }

    pub fn smiles(&self, s: &[f64]) -> String {
        let s_iter = &mut s.iter().copied();
        self.alkyl_tails()
            .into_iter()
            .zip(self.functional_group.smiles().iter())
            .fold(String::new(), |mut acc, (s, sm)| {
                acc.push_str(sm);
                acc.push_str(
                    &SuperAlkyl::smiles(s, s_iter)
                        .into_iter()
                        .collect::<String>(),
                );
                acc
            })
    }
}

impl<const N_Y: usize> MolecularRepresentation<N_Y> for SuperMolecule {
    fn structure_variables(&self) -> [(i32, i32); N_Y] {
        [(0, 1); N_Y]
    }

    fn constraints(&self, y: [Variable; N_Y]) -> Vec<Constraint> {
        assert_eq!(N_Y, self.variables());
        self.constraints(&y)
            .into_iter()
            .flat_map(|(vars, coefs, con)| {
                let expr: Expression = vars.iter().zip(&coefs).map(|(&v, &c)| v * c).sum();
                let mut cons = Vec::new();
                match con {
                    GeneralConstraint::Equality(eq) => cons.push(constraint!(expr == eq)),
                    GeneralConstraint::Inequality(lhs, rhs) => {
                        if let Some(lhs) = lhs {
                            cons.push(constraint!(expr.clone() >= lhs));
                        }
                        if let Some(rhs) = rhs {
                            cons.push(constraint!(expr <= rhs));
                        }
                    }
                }

                cons
            })
            .collect()
    }

    fn smiles(&self, y: [f64; N_Y]) -> String {
        self.smiles(&y)
    }

    fn build_molecule<D: DualNum<f64> + Copy>(&self, y: [D; N_Y]) -> ChemicalRecord<D> {
        self.build(y.to_vec())
    }
}

impl SuperMolecule {
    fn build_ketone<D: DualNum<f64> + Copy>(&self, y: Vec<D>) -> ChemicalRecord<D> {
        let y_iter = &mut y.into_iter();

        let cs = ["CH3", "CH2", ">CH", ">C<"];
        let mut c_groups = Polynomial::zero();
        let mut c_bonds = Polynomial2::zero();

        // formaldehyde is approximated by CH=O for the relaxation.
        let cds = ["CH=O", "CH=O", ">C=O"];
        let mut cd_bonds: Polynomial2<D, 4, 3> = Polynomial2::zero();

        let alkyls = self
            .alkyl_tails()
            .iter()
            .filter(|&&s| s > 0)
            .map(|&s| SuperAlkyl::build(s, y_iter, &mut c_groups, &mut c_bonds))
            .collect();
        let cd_groups = calculate_bonds(alkyls, &mut cd_bonds);

        let groups: HashMap<_, _> = cs
            .into_iter()
            .chain(cds)
            .zip(c_groups.iter().chain(cd_groups.iter()).cloned())
            .collect();

        let mut bonds = HashMap::new();
        fill_bond_map(cs, cs, &c_bonds, &mut bonds);
        fill_bond_map(cs, cds, &cd_bonds, &mut bonds);

        ChemicalRecord::new(groups, bonds)
    }

    fn build_alkene<D: DualNum<f64> + Copy>(&self, y: Vec<D>) -> ChemicalRecord<D> {
        let y_iter = &mut y.into_iter();

        let cs = ["CH3", "CH2", ">CH", ">C<"];
        let mut c_groups = Polynomial::zero();
        let mut c_bonds = Polynomial2::zero();

        let cds = ["=CH2", "=CH", "=C<"];
        let mut cd_groups = Vec::new();
        let mut cd_bonds: Polynomial2<D, 4, 3> = Polynomial2::zero();

        let s = self.alkyl_tails();
        for tails in [[0, 1], [2, 3]] {
            let alkyls = tails
                .iter()
                .filter(|&&i| s[i] > 0)
                .map(|&i| SuperAlkyl::build(s[i], y_iter, &mut c_groups, &mut c_bonds))
                .collect();
            cd_groups.push(calculate_bonds(alkyls, &mut cd_bonds));
        }

        let groups: HashMap<_, _> = cs
            .into_iter()
            .chain(cds)
            .zip(
                c_groups
                    .iter()
                    .chain((cd_groups[0] + cd_groups[1]).iter())
                    .cloned(),
            )
            .collect();

        let mut bonds = HashMap::new();
        fill_bond_map(cs, cs, &c_bonds, &mut bonds);
        fill_bond_map(cs, cds, &cd_bonds, &mut bonds);
        let dd_bonds = cd_groups[0].outer_product(cd_groups[1]);
        fill_bond_map(cds, cds, &dd_bonds, &mut bonds);

        ChemicalRecord::new(groups, bonds)
    }
}

fn merge_symmetry_constraints(
    branch1: &(Vec<usize>, Vec<f64>),
    branch2: &(Vec<usize>, Vec<f64>),
) -> (Vec<usize>, Vec<f64>) {
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

fn calculate_bonds<D: DualNum<f64> + Copy, const M: usize, const N: usize>(
    children: Vec<Polynomial<D, N>>,
    bonds: &mut Polynomial2<D, N, M>,
) -> Polynomial<D, M> {
    for (i, ci) in children.iter().enumerate() {
        *bonds += ci.outer_product(
            children
                .iter()
                .enumerate()
                .filter(|&(j, _)| (j != i))
                .map(|(_, cj)| Polynomial::new(cj.sum()))
                .product::<Polynomial<D, M>>()
                * Polynomial::new(D::one()),
        );
    }
    children.iter().map(|p| Polynomial::new(p.sum())).product()
}

fn fill_bond_map<D: DualNum<f64> + Copy, const M: usize, const N: usize>(
    seg1: [&'static str; M],
    seg2: [&'static str; N],
    bond_poly: &Polynomial2<D, M, N>,
    bonds: &mut HashMap<[&'static str; 2], D>,
) {
    for (i, c1) in seg1.iter().enumerate() {
        for (j, c2) in seg2.iter().enumerate() {
            *bonds.entry([c1, c2]).or_insert(D::zero()) += bond_poly[(i, j)];
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use approx::assert_relative_eq;
    use itertools::Itertools;
    use std::iter;

    fn isomers(molecule: SuperMolecule) -> usize {
        let mut res = 0;
        let (_, size) = molecule.size_constraint();
        let bond_constraints = molecule.bond_constraints(0);
        let symmetry_constraints = molecule.symmetry_constraints(0);
        for y in vec![0..=1; molecule.variables()]
            .into_iter()
            .multi_cartesian_product()
        {
            if y[0] == 0 {
                continue;
            }
            if y.iter()
                .zip(iter::repeat(1.0))
                .map(|(&y, s)| y as f64 * s)
                .sum::<f64>()
                > size
            {
                continue;
            }
            if bond_constraints.iter().any(|&[a, b]| y[a] < y[b]) {
                continue;
            }
            if symmetry_constraints.iter().any(|(i, l)| {
                let x = i
                    .iter()
                    .zip(l.iter())
                    .map(|(&i, &l)| y[i] as isize * l as isize)
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
        assert_eq!(isomers(SuperMolecule::alcohol(0, 4)), 4);
        assert_eq!(isomers(SuperMolecule::alcohol(0, 5)), 8);
        assert_eq!(isomers(SuperMolecule::alcohol(0, 6)), 16);
        assert_eq!(isomers(SuperMolecule::alcohol(0, 7)), 33);
        assert_eq!(isomers(SuperMolecule::alcohol(0, 8)), 72);
        // assert_eq!(isomers(SuperMolecule::alcohol(9)), 161);
    }

    #[test]
    fn test_isomers_alkenes() {
        assert_eq!(isomers(SuperMolecule::alkene(0, 3)), 1);
        assert_eq!(isomers(SuperMolecule::alkene(0, 4)), 4);
        assert_eq!(isomers(SuperMolecule::alkene(0, 5)), 9);
        assert_eq!(isomers(SuperMolecule::alkene(0, 6)), 23);
    }

    #[test]
    fn test_isomers_ketones() {
        assert_eq!(isomers(SuperMolecule::ketone(0, 3)), 1);
        assert_eq!(isomers(SuperMolecule::ketone(0, 4)), 3);
        assert_eq!(isomers(SuperMolecule::ketone(0, 5)), 6);
        assert_eq!(isomers(SuperMolecule::ketone(0, 6)), 13);
    }

    #[test]
    fn test_build() {
        let ChemicalRecord { groups, bonds } =
            SuperMolecule::alcohol(0, 4).build(vec![0.9, 0.8, 0.5, 0.3]);
        println!("{groups:?}\n{bonds:?}");
        assert_relative_eq!(groups["CH3"], 1.326);
        assert_relative_eq!(groups["CH2"], 0.958);
        assert_relative_eq!(groups[">CH"], 0.216);
        assert_relative_eq!(bonds[&["CH2", ">CH"]], 0.12);
        assert_relative_eq!(bonds[&["CH3", ">CH"]], 0.36);
        assert_relative_eq!(bonds[&["CH2", "CH2"]], 0.28);
        assert_relative_eq!(bonds[&["CH3", "CH2"]], 0.84);
    }

    #[test]
    fn test_alcohol() {
        let ChemicalRecord { groups, bonds } = SuperMolecule::alcohol(0, 5).build(vec![1.0; 7]);
        println!("{groups:?}\n{bonds:?}");
        assert_eq!(groups["CH3"], 4.0);
        assert_eq!(groups["CH2"], 1.0);
        assert_eq!(groups[">CH"], 1.0);
        assert_eq!(groups[">C<"], 1.0);
        assert_eq!(groups["OH"], 1.0);
        assert_relative_eq!(bonds[&["CH3", ">C<"]], 2.0);
        assert_relative_eq!(bonds[&[">CH", ">C<"]], 1.0);
        assert_relative_eq!(bonds[&["CH3", ">CH"]], 1.0);
        assert_relative_eq!(bonds[&["CH2", ">CH"]], 1.0);
        assert_relative_eq!(bonds[&["CH3", "CH2"]], 1.0);
        assert_relative_eq!(bonds[&["OH", ">C<"]], 1.0);
    }

    #[test]
    fn test_ketone() {
        let ChemicalRecord { groups, bonds } = SuperMolecule::ketone(0, 5).build(vec![1.0; 5]);
        println!("{groups:?}\n{bonds:?}");
        assert_eq!(groups["CH3"], 3.0);
        assert_eq!(groups["CH2"], 1.0);
        assert_eq!(groups[">CH"], 1.0);
        assert_eq!(groups[">C=O"], 1.0);
        assert_relative_eq!(bonds[&["CH3", ">CH"]], 1.0);
        assert_relative_eq!(bonds[&["CH2", ">CH"]], 1.0);
        assert_relative_eq!(bonds[&["CH3", "CH2"]], 1.0);
        assert_relative_eq!(bonds[&[">CH", ">C=O"]], 1.0);
        assert_relative_eq!(bonds[&["CH3", ">C=O"]], 1.0);
    }

    #[test]
    fn test_alkene() {
        let ChemicalRecord { groups, bonds } = SuperMolecule::alkene(0, 5).build(vec![1.0; 6]);
        println!("{groups:?}\n{bonds:?}");
        assert_eq!(groups["CH3"], 4.0);
        assert_eq!(groups["CH2"], 1.0);
        assert_eq!(groups[">CH"], 1.0);
        assert_eq!(groups["=C<"], 1.0);
        assert_eq!(groups["=CH"], 1.0);
        assert_relative_eq!(bonds[&["CH3", ">CH"]], 1.0);
        assert_relative_eq!(bonds[&["CH2", ">CH"]], 1.0);
        assert_relative_eq!(bonds[&["CH3", "CH2"]], 1.0);
        assert_relative_eq!(bonds[&[">CH", "=C<"]], 1.0);
        assert_relative_eq!(bonds[&["CH3", "=C<"]], 1.0);
        assert_relative_eq!(bonds[&["CH3", "=CH"]], 1.0);
        assert_relative_eq!(bonds[&["=C<", "=CH"]], 1.0);
    }
}
