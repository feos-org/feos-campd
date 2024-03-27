use super::MolecularRepresentation;
use crate::variables::{Constraint, ExplicitVariable, StructureVariables, Variable};
use indexmap::IndexMap;

impl<M: MolecularRepresentation, const N: usize> MolecularRepresentation for [M; N] {
    fn structure_variables(&self) -> StructureVariables {
        let mut variables = self
            .iter()
            .map(|s| s.structure_variables())
            .max_by_key(|v| v.len())
            .unwrap();
        let n_y = variables.len();
        (0..N).for_each(|i| {
            (0..n_y).for_each(|j| {
                variables.insert(format!("v{i}_{j}"), Variable::binary());
            });
        });
        // variables.push(Variable::binary("c0".into()).init(1.0));
        (0..N).for_each(|i| {
            variables.insert(format!("c{i}"), Variable::binary());
        });
        variables
    }

    fn feature_variables(
        &self,
        index_structure_vars: &[i32],
    ) -> IndexMap<String, ExplicitVariable> {
        // structure variables
        let n_y = (index_structure_vars.len() - N) / (N + 1);
        let (y, mut c) = index_structure_vars.split_at(n_y);
        (0..N).for_each(|_| {
            let (_, v2) = c.split_at(n_y);
            c = v2;
        });

        let mut feature_variables: IndexMap<_, Vec<_>> = IndexMap::new();
        for (i, m) in self.iter().enumerate() {
            m.feature_variables(y)
                .into_iter()
                .for_each(|(k, v)| feature_variables.entry(k).or_default().push((i, v)));
        }

        feature_variables
            .into_iter()
            .map(|(k, variables)| {
                let mut lvars = Vec::new();
                let mut lcoefs = Vec::new();
                let mut qvars1 = Vec::new();
                let mut qvars2 = Vec::new();
                let mut qcoefs = Vec::new();

                variables.into_iter().for_each(|(i, v)| {
                    lvars.extend(v.lvars.into_iter().map(|y| y + ((i + 1) * n_y) as i32));
                    lcoefs.extend(v.lcoefs);
                    qvars1.extend(v.qvars1.into_iter().map(|y| y + ((i + 1) * n_y) as i32));
                    qvars2.extend(v.qvars2.into_iter().map(|y| y + ((i + 1) * n_y) as i32));
                    qcoefs.extend(v.qcoefs);
                    if v.cons != 0.0 {
                        lvars.push(c[i]);
                        lcoefs.push(v.cons);
                    }
                });

                (
                    k.clone(),
                    ExplicitVariable::new(k)
                        .linear_struct(lvars, lcoefs)
                        .quadratic_struct(qvars1, qvars2, qcoefs),
                )
            })
            .collect()
    }

    fn constraints(&self, index_structure_vars: &[i32]) -> Vec<Constraint> {
        // structure variables
        let n_y = (index_structure_vars.len() - N) / (N + 1);
        let (y, mut c) = index_structure_vars.split_at(n_y);
        let v: Vec<_> = (0..N)
            .map(|_| {
                let (v1, v2) = c.split_at(n_y);
                c = v2;
                v1
            })
            .collect();

        // disjunction
        let mut constraints = Vec::new();
        let coefs = vec![1.0; c.len()];
        constraints.push(
            Constraint::new()
                .linear_struct(c.to_vec(), coefs)
                .eqbnd(1.0),
        );

        // sum(v) = y
        for (i, y) in y.iter().enumerate() {
            let mut vars = vec![*y];
            let mut coefs = vec![1.0];
            for &v in &v {
                vars.push(v[i]);
                coefs.push(-1.0);
            }
            constraints.push(Constraint::new().linear_struct(vars, coefs).eqbnd(0.0));
        }

        // v <= c
        for (&v, &c) in v.iter().zip(c) {
            for &v in v {
                constraints.push(
                    Constraint::new()
                        .linear_struct(vec![v, c], vec![1.0, -1.0])
                        .upbnd(0.0),
                )
            }
        }

        for (i, (m, &c)) in self.iter().zip(c).enumerate() {
            let (used_vars, unused_vars) = y.split_at(m.structure_variables().len());
            let mut constr = m.constraints(used_vars);

            // unused variable constraint
            if !unused_vars.is_empty() {
                let coefs = vec![1.0; unused_vars.len()];
                constr.push(
                    Constraint::new()
                        .linear_struct(unused_vars.to_vec(), coefs)
                        .upbnd(0.0),
                );
            }

            // modified constraints with convex hull
            for constr in constr {
                let mut vars = constr.lvars;
                vars.iter_mut().for_each(|y| *y += ((i + 1) * n_y) as i32);
                let mut coefs = constr.lcoefs;
                if let Some(lobnd) = constr.lobnd {
                    if lobnd != 0.0 {
                        vars.push(c);
                        coefs.push(-lobnd);
                    }
                    constraints.push(Constraint::new().linear_struct(vars, coefs).lobnd(0.0));
                } else if let Some(upbnd) = constr.upbnd {
                    if upbnd != 0.0 {
                        vars.push(c);
                        coefs.push(-upbnd);
                    }
                    constraints.push(Constraint::new().linear_struct(vars, coefs).upbnd(0.0));
                }
            }
        }

        constraints
    }

    fn smiles(&self, y: &[usize]) -> Vec<String> {
        let (y, c) = y.split_at(y.len() - N);
        for (&c, m) in c.iter().zip(self.iter()) {
            if c == 1 {
                return m.smiles(y);
            }
        }
        unreachable!();
    }
}
