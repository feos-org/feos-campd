use super::{ContinuousVariables, DiscreteVariables, LinearConstraint, MolecularRepresentation};

impl<
        C,
        M1: MolecularRepresentation<1, ChemicalRecord = C>,
        M2: MolecularRepresentation<1, ChemicalRecord = C>,
    > MolecularRepresentation<2> for (M1, M2)
{
    type ChemicalRecord = C;

    fn build(&self, y: &[f64], p: &[f64]) -> [Self::ChemicalRecord; 2] {
        let (y1, y2) = y.split_at(self.0.structure_variables().len());
        let (p1, p2) = p.split_at(self.0.parameter_variables().len());
        let [cr1] = self.0.build(y1, p1);
        let [cr2] = self.1.build(y2, p2);
        [cr1, cr2]
    }

    fn structure_variables(&self) -> DiscreteVariables {
        let mut variables = self.0.structure_variables();
        variables.append(&mut self.1.structure_variables());
        variables
    }

    fn parameter_variables(&self) -> ContinuousVariables {
        let mut variables = self.0.parameter_variables();
        variables.append(&mut self.1.parameter_variables());
        variables
    }

    fn determine_parameters(&self, y: &[f64]) -> Vec<f64> {
        let (y1, y2) = y.split_at(self.0.structure_variables().len());
        let p1 = self.0.determine_parameters(y1);
        let p2 = self.0.determine_parameters(y2);
        [p1, p2].concat()
    }

    fn constraints(
        &self,
        index_structure_vars: &[i32],
        index_parameter_vars: Option<&[i32]>,
    ) -> Vec<LinearConstraint> {
        let (index_structure_vars1, index_structure_vars2) =
            index_structure_vars.split_at(self.0.structure_variables().len());
        let (index_parameter_vars1, index_parameter_vars2) =
            if let Some(index_parameter_vars) = index_parameter_vars {
                let (v1, v2) = index_parameter_vars.split_at(self.0.parameter_variables().len());
                (Some(v1), Some(v2))
            } else {
                (None, None)
            };
        let mut constraints = self
            .0
            .constraints(index_structure_vars1, index_parameter_vars1);
        constraints.append(
            &mut self
                .1
                .constraints(index_structure_vars2, index_parameter_vars2),
        );
        constraints
    }

    fn smiles(&self, y: &[usize]) -> [String; 2] {
        let (y1, y2) = y.split_at(self.0.structure_variables().len());
        let [smiles1] = self.0.smiles(y1);
        let [smiles2] = self.1.smiles(y2);
        [smiles1, smiles2]
    }
}
