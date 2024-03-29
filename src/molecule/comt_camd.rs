use super::MolecularRepresentation;
use crate::variables::{Constraint, ExplicitVariable, StructureVariables, Variable};
use feos::core::parameter::ParameterError;
use indexmap::IndexMap;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::BufReader;
use std::iter;
use std::path::Path;

#[derive(Clone, Serialize, Deserialize)]
struct Structure {
    name: String,
    groups: Option<Vec<String>>,
    count: usize,
    rings: usize,
}

#[derive(Clone, Serialize, Deserialize)]
struct Group {
    name: String,
    open_bonds: usize,
    n_max: usize,
}

#[derive(Clone, Serialize, Deserialize)]
pub struct CoMTCAMD {
    structures: Vec<Structure>,
    groups: Vec<Group>,
}

impl CoMTCAMD {
    pub fn from_json<FP: AsRef<Path>>(file: FP) -> Result<Self, ParameterError> {
        Ok(serde_json::from_reader(BufReader::new(File::open(file)?))?)
    }

    pub fn from_molecules(molecules: Vec<String>) -> Self {
        let groups = molecules
            .into_iter()
            .map(|name| Group {
                name,
                open_bonds: 0,
                n_max: 1,
            })
            .collect();
        let structures = vec![Structure {
            name: "molecule".into(),
            groups: None,
            count: 1,
            rings: 0,
        }];
        Self { structures, groups }
    }

    pub fn from_json_molecules<FP: AsRef<Path>>(file: FP) -> Result<Self, ParameterError> {
        Ok(Self::from_molecules(serde_json::from_reader(
            BufReader::new(File::open(file)?),
        )?))
    }

    pub fn get_initial_values(&self, structure: &str, groups: &IndexMap<&str, usize>) -> Vec<f64> {
        self.structures
            .iter()
            .map(|s| if s.name == structure { 1.0 } else { 0.0 })
            .chain(self.groups.iter().flat_map(|g| {
                let count = *groups.get(&g.name as &str).unwrap_or(&0);
                std::iter::repeat(1.0)
                    .take(count)
                    .chain(std::iter::repeat(0.0).take(g.n_max - count))
            }))
            .collect()
    }

    pub fn get_initial_values_molecules(&self, molecule: &str) -> Vec<f64> {
        let mut groups = IndexMap::new();
        groups.insert(molecule, 1);
        self.get_initial_values("molecule", &groups)
    }
}

impl MolecularRepresentation for CoMTCAMD {
    fn structure_variables(&self) -> StructureVariables {
        self.structures
            .iter()
            .map(|s| (s.name.clone(), Variable::binary()))
            .chain(self.groups.iter().flat_map(|g| {
                (0..g.n_max).map(|i| {
                    let name = if g.n_max == 1 {
                        g.name.to_string()
                    } else {
                        format!("{}_{i}", g.name)
                    };
                    (name, Variable::binary())
                })
            }))
            .collect()
    }

    fn feature_variables(
        &self,
        index_structure_vars: &[i32],
    ) -> IndexMap<String, ExplicitVariable> {
        let mut index = index_structure_vars[self.structures.len()];
        self.groups
            .iter()
            .map(|g| {
                let mut vars = Vec::new();
                let mut coefs = Vec::new();
                for _ in 0..g.n_max {
                    vars.push(index);
                    coefs.push(1.0);
                    index += 1;
                }
                (
                    g.name.clone(),
                    ExplicitVariable::new(g.name.clone()).linear_struct(vars, coefs),
                )
            })
            .collect()
    }

    fn constraints(&self, index_structure_vars: &[i32]) -> Vec<Constraint> {
        let mut constraints = Vec::new();

        let (y_vars, n_vars) = index_structure_vars.split_at(self.structures.len());

        let mut n: IndexMap<_, _> = IndexMap::new();
        self.groups
            .iter()
            .flat_map(|g| std::iter::repeat(&g.name as &str).take(g.n_max))
            .zip(n_vars.to_vec())
            .for_each(|(g, n_var)| n.entry(g).or_insert(vec![]).push(n_var));

        // Exactly one structure can be active
        constraints.push(
            Constraint::new()
                .linear_struct(y_vars.to_vec(), vec![1.0; y_vars.len()])
                .eqbnd(1.0),
        );

        // Connect molecular structures and segments
        for (structure, &y) in self.structures.iter().zip(y_vars) {
            if let Some(groups) = &structure.groups {
                if !groups.is_empty() {
                    let (vars, coefs): (Vec<_>, Vec<_>) = groups
                        .iter()
                        .flat_map(|g| n[g as &str].iter().copied())
                        .zip(std::iter::repeat(1.0))
                        .chain(std::iter::once((y, -(structure.count as f64))))
                        .unzip();
                    constraints.push(Constraint::new().linear_struct(vars, coefs).eqbnd(0.0));
                }
            } else {
                let coefs = vec![1.0; n_vars.len()];
                constraints.push(
                    Constraint::new()
                        .linear_struct(n_vars.to_vec(), coefs)
                        .eqbnd(1.0),
                );
            };
        }

        // octett rule
        let n_open_bonds = self
            .groups
            .iter()
            .flat_map(|g| iter::repeat(g.open_bonds).take(g.n_max));
        let (vars, coefs): (Vec<_>, Vec<_>) = n_vars
            .iter()
            .zip(n_open_bonds)
            .map(|(&v, o)| (v, 2.0 - o as f64))
            .chain(
                self.structures
                    .iter()
                    .zip(y_vars)
                    .map(|(s, &y)| (y, 2.0 * (s.rings as f64 - 1.0))),
            )
            .unzip();
        constraints.push(Constraint::new().linear_struct(vars, coefs).eqbnd(0.0));

        // // single molecule
        // for &n in n_vars {
        //     let (vars, coefs): (Vec<_>, Vec<_>) = n_vars
        //         .iter()
        //         .zip(&self.n_open_bonds)
        //         .map(|(&v, &o)| (v, if v == n { 2.0 - o as f64 } else { 1.0 }))
        //         .chain([y["alkane"], y["alkene"], y["P"], y["M"]].map(|i| (i, -2.0)))
        //         .unzip();
        //     constraints.push(Constraint::new(vars, coefs).lobnd(0.0));
        // }

        // order groups
        for vars in n.values() {
            for vars in vars.windows(2) {
                constraints.push(
                    Constraint::new()
                        .linear_struct(vars.to_vec(), vec![1.0, -1.0])
                        .lobnd(0.0),
                )
            }
        }

        constraints
    }

    fn smiles(&self, y: &[usize]) -> String {
        let (_, n) = y.split_at(self.structures.len());
        let groups: Vec<_> = self
            .groups
            .iter()
            .flat_map(|g| std::iter::repeat(&g.name as &str).take(g.n_max))
            .zip(n)
            .filter_map(|(g, &n)| (n == 1).then_some(g))
            .collect();
        groups.join(",")
    }
}
