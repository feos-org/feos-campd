use feos_core::EosError;
use good_lp::{Constraint, Variable};
use nalgebra::{SMatrix, SVector};
use num_dual::DualNum;
use std::array;
use std::collections::{HashMap, HashSet};

mod molecule;
mod process;
mod property;
mod solver;
pub use molecule::{CoMTCAMD, Disjunction, MolecularRepresentation, SuperMolecule};
pub use process::{ContinuousVariable, ProcessModel};
pub use property::{GcPcSaftPropertyModel, PcSaftPropertyModel, PropertyModel};
pub use solver::{
    GeneralConstraint, MixedIntegerNonLinearProgram, OptimizationResult, OuterApproximation,
};

/// Input for group-contribution models that allows for derivatives.
pub struct ChemicalRecord<D> {
    pub groups: HashMap<&'static str, D>,
    pub bonds: HashMap<[&'static str; 2], D>,
}

impl<D> ChemicalRecord<D> {
    pub fn new(groups: HashMap<&'static str, D>, bonds: HashMap<[&'static str; 2], D>) -> Self {
        Self { groups, bonds }
    }
}

#[derive(Clone)]
pub struct IntegratedDesign<M, R, P> {
    molecule: M,
    property: R,
    process: P,
}

impl<M, R, P> IntegratedDesign<M, R, P> {
    pub fn new(molecule: M, property: R, process: P) -> Self {
        Self {
            molecule,
            property,
            process,
        }
    }
}

impl<
        M: MolecularRepresentation<N_Y>,
        R: PropertyModel<N>,
        P: ProcessModel<R::EquationOfState, N_X, N>,
        const N_X: usize,
        const N_Y: usize,
        const N: usize,
    > MixedIntegerNonLinearProgram<N_X, N_Y, N> for IntegratedDesign<M, R, P>
{
    type Error = EosError;

    fn x_variables(&self) -> SVector<(f64, f64, f64), N_X> {
        SVector::from(self.process.variables().map(|v| (v.lobnd, v.upbnd, v.init)))
    }

    fn y_variables(&self) -> SMatrix<(i32, i32), N_Y, N> {
        SMatrix::from(array::from_fn(|_| self.molecule.structure_variables()))
    }

    fn linear_constraints(&self, y: SMatrix<Variable, N_Y, N>) -> Vec<Constraint> {
        y.data
            .0
            .iter()
            .flat_map(|&y| self.molecule.constraints(y))
            .collect()
    }

    fn constraints(&self) -> Vec<GeneralConstraint> {
        self.process.constraints()
    }

    fn evaluate<D: DualNum<f64> + Copy>(
        &self,
        x: SVector<D, N_X>,
        y: SMatrix<D, N_Y, N>,
    ) -> Result<(D, Vec<D>), EosError> {
        let y_set: HashSet<_> = y.data.0.iter().map(|y| y.map(|y| y.re() as i32)).collect();
        if y_set.len() != N {
            Err(EosError::IncompatibleComponents(N, y_set.len()))
        } else {
            let cr = y.data.0.map(|y| self.molecule.build_molecule(y));
            let eos = self.property.build_eos(cr.each_ref());
            self.process.evaluate(&eos, cr.each_ref(), x.data.0[0])
        }
    }

    fn y_to_string(&self, y: &SMatrix<f64, N_Y, N>) -> String {
        y.data.0.map(|y| self.molecule.smiles(y)).join(";")
    }

    fn exclude_solutions(&self, s: &[f64]) -> Vec<Vec<f64>> {
        let mut exclude = vec![s.to_vec()];
        if N == 2 {
            let (s1, s2) = s.split_at(s.len() / 2);
            exclude.push([s2, s1].concat())
        }
        exclude
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use approx::assert_relative_eq;
    use good_lp::highs;
    use process::OrganicRankineCycle;

    #[test]
    fn test_solve() {
        // Heptane
        let mut y = [0.0; 28];
        y[0] = 1.0;
        y[6] = 2.0;
        y[7] = 5.0;

        let property =
            PcSaftPropertyModel::full("../feos/parameters/pcsaft/rehner2023_homo.json", None)
                .unwrap();

        let campd = IntegratedDesign::new(CoMTCAMD, property, OrganicRankineCycle::default());
        let solver = OuterApproximation::new(&campd);
        solver.solve_ranking(SVector::from(y), highs, 5, &Default::default());
    }

    #[test]
    fn test_solve_process() {
        // Isopentane
        let mut y = [0.0; 28];
        y[0] = 1.0;
        y[6] = 3.0;
        y[7] = 1.0;
        y[8] = 1.0;

        let property =
            PcSaftPropertyModel::full("../feos/parameters/pcsaft/rehner2023_homo.json", None)
                .unwrap();

        let campd = IntegratedDesign::new(CoMTCAMD, property, OrganicRankineCycle::default());
        let mut oa = OuterApproximation::new(&campd);
        let result = oa.solve_nlp(SVector::from(y), vec![]).unwrap();
        println!("{:8.5} {:.5?}", result.objective.0, result.x.data.0[0]);
    }

    #[test]
    fn test_supermolecule_fixed() {
        const MOLECULE: SuperMolecule = SuperMolecule::alkene(0, 5);
        const N_Y: usize = MOLECULE.variables();

        // Propene
        let mut y = [0.0; N_Y];
        y[0] = 1.0;
        let campd = IntegratedDesign::new(
            MOLECULE,
            GcPcSaftPropertyModel,
            OrganicRankineCycle::default(),
        );
        let mut oa = OuterApproximation::new(&campd);
        let result = oa.solve_nlp(SVector::from(y), vec![]).unwrap();
        println!(
            "{:8.5} {:.5?} {}",
            result.objective.0,
            result.x.data.0[0],
            MOLECULE.smiles(&result.y.data.0[0])
        );
        assert_relative_eq!(
            result.objective.0,
            -0.4378352970105434,
            max_relative = 1e-10
        );
    }

    #[test]
    fn test_supermolecule() {
        let process = OrganicRankineCycle::default();
        let molecule = SuperMolecule::alkene(0, 5);
        let campd = IntegratedDesign::new(molecule, GcPcSaftPropertyModel, process);
        let solver = OuterApproximation::new(&campd);

        // Propene
        let y = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0];
        solver.solve_ranking(SVector::from(y), highs, 1, &Default::default());
    }

    #[test]
    fn test_supermolecule_disjunction_fixed() {
        let process = OrganicRankineCycle::default();
        let molecule = SuperMolecule::non_associating(0, 5);
        let campd = IntegratedDesign::new(molecule.clone(), GcPcSaftPropertyModel, process);
        let mut oa = OuterApproximation::new(&campd);

        // Propene
        let y0 = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let s = SuperMolecule::get_initial_values(&molecule, "alkene", &y0);
        let y: [f64; 12] = array::from_fn(|i| s[i]);
        let result = oa.solve_nlp(SVector::from(y), vec![]).unwrap();
        println!(
            "{:8.5} {:.5?} {}",
            result.objective.0,
            result.x.data.0[0],
            molecule.smiles(&result.y.data.0[0])
        );
        assert_relative_eq!(
            result.objective.0,
            -0.4378352970105434,
            max_relative = 1e-10
        );
    }

    #[test]
    fn test_supermolecule_disjunction() {
        let process = OrganicRankineCycle::default();
        let molecule = SuperMolecule::non_associating(0, 5);
        let campd = IntegratedDesign::new(molecule.clone(), GcPcSaftPropertyModel, process);
        let solver = OuterApproximation::new(&campd);

        // Propene
        let y0 = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0];
        let s = SuperMolecule::get_initial_values(&molecule, "alkene", &y0);
        let y: [f64; 12] = array::from_fn(|i| s[i]);
        solver.solve_ranking(SVector::from(y), highs, 3, &Default::default());
    }
}
