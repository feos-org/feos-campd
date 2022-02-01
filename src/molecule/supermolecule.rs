use super::FunctionalGroup;
use feos_core::parameter::{ChemicalRecord, Identifier};
use num_dual::{StaticMat, StaticVec};
use num_traits::{NumAssign, Zero};
use petgraph::dot::{Config, Dot};
use petgraph::graph::{NodeIndex, UnGraph};
use petgraph::Graph;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct SuperMolecule {
    pub size: usize,
    functional_groups: Vec<FunctionalGroup>,
}

impl SuperMolecule {
    pub fn new(size: usize, functional_groups: Vec<FunctionalGroup>) -> Self {
        Self {
            size,
            functional_groups,
        }
    }

    pub fn variables(&self) -> usize {
        Self::variables_recursive(self.size) + self.functional_groups.len()
    }

    fn variables_recursive(size: usize) -> usize {
        1 + (1..=3)
            .map(|k| (size - 1) / k)
            .filter(|s| *s > 0)
            .map(Self::variables_recursive)
            .sum::<usize>()
    }

    pub fn functional_group_constraint(&self) -> Vec<usize> {
        (0..self.functional_groups.len()).collect()
    }

    pub fn bond_constraints(&self) -> Vec<(usize, usize)> {
        let mut index = self.functional_groups.len();
        let mut constraints = Vec::new();
        Self::bond_constraints_recursive(self.size, &mut index, &mut constraints);
        constraints
    }

    fn bond_constraints_recursive(
        size: usize,
        index: &mut usize,
        constraints: &mut Vec<(usize, usize)>,
    ) {
        let this = *index;
        (1..=3)
            .map(|k| (size - 1) / k)
            .filter(|s| *s > 0)
            .for_each(|s| {
                constraints.push((this, *index + 1));
                *index += 1;
                Self::bond_constraints_recursive(s, index, constraints);
            });
    }

    pub fn symmetry_constraints(&self) -> Vec<(Vec<usize>, Vec<isize>)> {
        let mut index = self.functional_groups.len();
        let mut constraints = Vec::new();
        Self::symmetry_constraints_recursive(self.size, &mut index, 0, &mut constraints);
        constraints
    }

    fn symmetry_constraints_recursive(
        size: usize,
        index: &mut usize,
        level: isize,
        constraints: &mut Vec<(Vec<usize>, Vec<isize>)>,
    ) -> (Vec<usize>, Vec<isize>) {
        let mut indices = vec![*index];
        let mut levels = vec![level];
        let children: Vec<_> = (1..=3)
            .map(|k| (size - 1) / k)
            .filter(|s| *s > 0)
            .map(|s| {
                *index += 1;
                Self::symmetry_constraints_recursive(s, index, level + 1, constraints)
            })
            .collect();

        children
            .iter()
            .cloned()
            .zip(children.iter().cloned().skip(1))
            .for_each(|((indices1, levels1), (indices2, levels2))| {
                constraints.push((
                    indices1.into_iter().chain(indices2.into_iter()).collect(),
                    levels1
                        .into_iter()
                        .chain(levels2.iter().map(|&l| -l))
                        .collect(),
                ));
            });

        children.into_iter().for_each(|(i, l)| {
            indices.extend(i);
            levels.extend(l);
        });
        (indices, levels)
    }

    pub fn build(&self, y: Vec<f64>) -> ChemicalRecord {
        let mut y_iter = y.into_iter();
        let y_func: Vec<_> = self
            .functional_groups
            .iter()
            .map(|_| y_iter.next().unwrap())
            .collect();

        let mut c_segments = StaticVec::zero();
        let mut c_bonds = StaticMat::zero();
        let s = Self::build_recursive(self.size, 1.0, &mut y_iter, &mut c_segments, &mut c_bonds);

        let identifier = Identifier::new("", Some("SuperMolecule"), None, None, None, None);
        let cs = [">C<", ">CH", "CH2", "CH3"];

        let mut segments: HashMap<_, _> = cs
            .iter()
            .map(|&s| s.to_string())
            .zip(c_segments.iter().cloned())
            .filter_map(|(s, g)| (g > 0.0).then(|| (s, g)))
            .collect();
        self.functional_groups
            .iter()
            .zip(y_func.iter())
            .filter(|(_, &y)| y > 0.0)
            .for_each(|(s, y)| *segments.entry(s.group.clone()).or_insert(0.0) += y);

        let mut bonds = HashMap::new();
        for (i, c1) in cs.iter().enumerate() {
            for (j, c2) in cs.iter().enumerate() {
                if c_bonds[(i, j)] > 0.0 {
                    *bonds.entry([c1.to_string(), c2.to_string()]).or_insert(0.0) +=
                        c_bonds[(i, j)];
                }
            }
        }
        s.iter().zip(cs.iter()).for_each(|(g, c)| {
            if *g > 0.0 {
                self.functional_groups
                    .iter()
                    .zip(y_func.iter())
                    .for_each(|(s, y)| {
                        if *y > 0.0 {
                            *bonds.entry([s.group.clone(), c.to_string()]).or_insert(0.0) += y * g
                        }
                    });
            }
        });

        ChemicalRecord::new_count(identifier, segments, Some(bonds))
    }

    fn build_recursive<T: NumAssign + Copy, I: Iterator<Item = T>>(
        size: usize,
        parent: T,
        y: &mut I,
        segments: &mut StaticVec<T, 4>,
        bonds: &mut StaticMat<T, 4, 4>,
    ) -> StaticVec<T, 4> {
        let value = y.next().unwrap() * parent;

        let children: Vec<_> = (1..=3)
            .map(|k| (size - 1) / k)
            .filter(|s| *s > 0)
            .map(|s| Self::build_recursive(s, value, y, segments, bonds))
            .collect();

        let mut res = StaticVec::zero();

        match children.len() {
            0 => res[3] = T::one(),
            1 => {
                let c1 = children[0].sum();
                res[2] = c1;
                res[3] = T::one() - c1;
                for (i, &c) in children[0].iter().enumerate() {
                    bonds[(i, 2)] += c;
                }
            }
            2 => {
                let c1 = children[0].sum();
                let c2 = children[1].sum();
                res[1] = c1 * c2;
                res[2] = c1 * (T::one() - c2) + (T::one() - c1) * c2;
                res[3] = (T::one() - c1) * (T::one() - c2);
                for (i, &g1) in children[0].iter().enumerate() {
                    for (j, &g2) in children[1].iter().enumerate() {
                        bonds[(i, 1)] += g1 * g2;
                        bonds[(j, 1)] += g1 * g2;
                    }
                }
                for (i, &g1) in children[0].iter().enumerate() {
                    bonds[(i, 2)] += g1 * (T::one() - c2);
                }
                for (j, &g2) in children[1].iter().enumerate() {
                    bonds[(j, 2)] += (T::one() - c1) * g2;
                }
            }
            3 => {
                let c1 = children[0].sum();
                let c2 = children[1].sum();
                let c3 = children[2].sum();
                res[0] = c1 * c2 * c3;
                res[1] = (T::one() - c1) * c2 * c3
                    + c1 * (T::one() - c2) * c3
                    + c1 * c2 * (T::one() - c3);
                res[2] = (T::one() - c1) * (T::one() - c2) * c3
                    + (T::one() - c1) * c2 * (T::one() - c3)
                    + c1 * (T::one() - c2) * (T::one() - c3);
                res[3] = (T::one() - c1) * (T::one() - c2) * (T::one() - c3);
                for (i, &g1) in children[0].iter().enumerate() {
                    for (j, &g2) in children[1].iter().enumerate() {
                        for (k, &g3) in children[2].iter().enumerate() {
                            bonds[(i, 0)] += g1 * g2 * g3;
                            bonds[(j, 0)] += g1 * g2 * g3;
                            bonds[(k, 0)] += g1 * g2 * g3;
                        }
                    }
                }
                for (i, &g1) in children[0].iter().enumerate() {
                    for (j, &g2) in children[1].iter().enumerate() {
                        bonds[(i, 1)] += g1 * g2 * (T::one() - c3);
                        bonds[(j, 1)] += g1 * g2 * (T::one() - c3);
                    }
                }
                for (i, &g1) in children[0].iter().enumerate() {
                    for (k, &g3) in children[2].iter().enumerate() {
                        bonds[(i, 1)] += g1 * (T::one() - c2) * g3;
                        bonds[(k, 1)] += g1 * (T::one() - c2) * g3;
                    }
                }
                for (j, &g2) in children[1].iter().enumerate() {
                    for (k, &g3) in children[2].iter().enumerate() {
                        bonds[(j, 1)] += (T::one() - c1) * g2 * g3;
                        bonds[(k, 1)] += (T::one() - c1) * g2 * g3;
                    }
                }
                for (i, &g1) in children[0].iter().enumerate() {
                    bonds[(i, 2)] += g1 * (T::one() - c2) * (T::one() - c3);
                }
                for (j, &g2) in children[1].iter().enumerate() {
                    bonds[(j, 2)] += (T::one() - c1) * g2 * (T::one() - c3);
                }
                for (k, &g3) in children[2].iter().enumerate() {
                    bonds[(k, 2)] += (T::one() - c1) * (T::one() - c2) * g3;
                }
            }
            _ => unreachable!(),
        };
        res *= value;
        *segments += res;

        res
    }

    #[allow(clippy::needless_collect)]
    pub fn indices(&self) -> String {
        let mut graph = Graph::new();
        let functional_groups: Vec<_> = self
            .functional_groups
            .iter()
            .map(|_| graph.add_node(()))
            .collect();
        let index = Self::indices_recursive(self.size, &mut graph, None);
        functional_groups.into_iter().for_each(|i| {
            graph.add_edge(i, index, ());
        });

        let graph = graph.map(|i, _| i.index(), |_, _| ());
        format!("{:?}", Dot::with_config(&graph, &[Config::EdgeNoLabel]))
    }

    fn indices_recursive(
        size: usize,
        graph: &mut Graph<(), ()>,
        parent: Option<NodeIndex>,
    ) -> NodeIndex {
        let this = graph.add_node(());
        if let Some(parent) = parent {
            graph.add_edge(parent, this, ());
        }
        (1..=3)
            .map(|k| (size - 1) / k)
            .filter(|s| *s > 0)
            .for_each(|s| {
                Self::indices_recursive(s, graph, Some(this));
            });
        this
    }

    pub fn molecule(&self, y: Vec<usize>) -> Molecule {
        let mut y_iter = y.into_iter();
        let mut graph = Graph::new();
        let mut index = NodeIndex::new(0);
        let mut functional_group = None;
        self.functional_groups.iter().for_each(|f| {
            if y_iter.next().unwrap() == 1 {
                functional_group = Some(f.clone());
                index = graph.add_node((true, f.group.clone()));
            }
        });
        Self::molecule_recursive(self.size, index, &mut y_iter, &mut graph);

        graph.retain_nodes(|g, y| g[y].0);
        let graph = graph.map(|_, n| n.1.clone(), |_, _| ());
        Molecule {
            graph,
            functional_group: functional_group.unwrap(),
        }
    }

    fn molecule_recursive<I: Iterator<Item = usize>>(
        size: usize,
        parent: NodeIndex,
        y: &mut I,
        graph: &mut Graph<(bool, String), ()>,
    ) -> bool {
        let value = y.next().unwrap();
        let index = graph.add_node((value == 1, "".into()));
        graph.add_edge(parent, index, ());
        let children = (1..=3)
            .map(|k| (size - 1) / k)
            .filter(|&s| s > 0)
            .filter(|&s| Self::molecule_recursive(s, index, y, graph))
            .count();

        graph[index].1 = (match children {
            0 => "CH3",
            1 => "CH2",
            2 => "CH",
            3 => "C",
            _ => unreachable!(),
        })
        .into();
        value == 1
    }
}

pub struct Molecule {
    graph: Graph<String, ()>,
    functional_group: FunctionalGroup,
}

impl Molecule {
    pub fn dot(&self) -> String {
        let ungraph: UnGraph<_, _> = self.graph.clone().into_edge_type();
        format!("{:?}", Dot::with_config(&ungraph, &[Config::EdgeNoLabel]))
    }

    pub fn smiles(&self) -> String {
        self.functional_group
            .smiles
            .chars()
            .chain(self.smiles_recursive(NodeIndex::new(1)))
            .collect()
    }

    fn smiles_recursive(&self, index: NodeIndex) -> Vec<char> {
        let mut children: Vec<_> = self
            .graph
            .neighbors(index)
            .map(|i| self.smiles_recursive(i))
            .collect();
        // children.reverse();
        // children.sort_by_key(|c| c.len());
        let mut smiles = vec!['C'];
        if let Some(last) = children.pop() {
            children.into_iter().for_each(|c| {
                smiles.push('(');
                smiles.extend(c);
                smiles.push(')');
            });
            smiles.extend(last);
        }
        smiles
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use itertools::Itertools;
    use std::io::prelude::*;

    fn isomers(size: usize) -> usize {
        let mut file = std::fs::File::create("octyl.dat").unwrap();
        let mut res = 0;
        let oh = FunctionalGroup {
            group: "OH".into(),
            smiles: vec!['O'],
        };
        let molecule = SuperMolecule::new(size, vec![oh]);
        let bond_constraints = molecule.bond_constraints();
        let symmetry_constraints = molecule.symmetry_constraints();
        for y in vec![0..=1; molecule.variables()]
            .into_iter()
            .multi_cartesian_product()
        {
            if y[0] + y[1] != 2 {
                continue;
            }
            if y.iter().sum::<isize>() as usize > size + 1 {
                continue;
            }
            if bond_constraints.iter().any(|(a, b)| y[*a] < y[*b]) {
                continue;
            }
            if symmetry_constraints.iter().any(|(i, l)| {
                let x = i
                    .iter()
                    .zip(l.iter())
                    .map(|(i, l)| y[*i] as isize * *l)
                    .sum::<isize>();
                x < 0
            }) {
                continue;
            }
            write!(file, "{:?}", y).unwrap();
            res += 1;
        }
        res
    }

    #[test]
    fn test_isomers() {
        assert_eq!(isomers(4), 8);
        assert_eq!(isomers(5), 16);
        assert_eq!(isomers(6), 33);
        assert_eq!(isomers(7), 72);
        assert_eq!(isomers(8), 161);
    }
}
