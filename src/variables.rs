use std::ops::{Deref, DerefMut};

#[cfg(feature = "knitro_rs")]
use crate::OptimizationMode;
#[cfg(feature = "knitro_rs")]
use knitro_rs::{Knitro, KnitroError};

#[derive(Clone, Copy)]
#[allow(dead_code)]
pub struct Variable<const INTEGER: bool, const FIXABLE: bool> {
    lobnd: Option<f64>,
    upbnd: Option<f64>,
}

pub struct Variables<const INTEGER: bool, const FIXABLE: bool>(Vec<Variable<INTEGER, FIXABLE>>);

impl<const INTEGER: bool, const FIXABLE: bool> Deref for Variables<INTEGER, FIXABLE> {
    type Target = Vec<Variable<INTEGER, FIXABLE>>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<const INTEGER: bool, const FIXABLE: bool> DerefMut for Variables<INTEGER, FIXABLE> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<const INTEGER: bool, const FIXABLE: bool> FromIterator<Variable<INTEGER, FIXABLE>>
    for Variables<INTEGER, FIXABLE>
{
    fn from_iter<T: IntoIterator<Item = Variable<INTEGER, FIXABLE>>>(iter: T) -> Self {
        Self(Vec::from_iter(iter))
    }
}

impl<const INTEGER: bool, const FIXABLE: bool> From<Vec<Variable<INTEGER, FIXABLE>>>
    for Variables<INTEGER, FIXABLE>
{
    fn from(value: Vec<Variable<INTEGER, FIXABLE>>) -> Self {
        Self(value)
    }
}

pub type ProcessVariables = Variables<false, false>;
pub type ParameterVariables = Variables<false, true>;
pub type StructureVariables = Variables<true, true>;

impl<const INTEGER: bool, const FIXABLE: bool> Variable<INTEGER, FIXABLE> {
    fn new(lobnd: Option<f64>, upbnd: Option<f64>) -> Self {
        Self { lobnd, upbnd }
    }
}

impl Variable<true, true> {
    pub fn binary() -> Self {
        Self::integer(0, 1)
    }

    pub fn integer(lobnd: usize, upbnd: usize) -> Self {
        let lobnd = lobnd as f64;
        let upbnd = upbnd as f64;
        Self::new(Some(lobnd), Some(upbnd))
    }
}

impl Variable<false, true> {
    pub fn free() -> Self {
        Self::new(None, None)
    }
}

impl Variable<false, false> {
    pub fn continuous(lobnd: f64, upbnd: f64) -> Self {
        Self::new(Some(lobnd), Some(upbnd))
    }
}

impl<const INTEGER: bool, const FIXABLE: bool> Variables<INTEGER, FIXABLE> {
    #[cfg(feature = "knitro_rs")]
    pub fn setup_knitro(
        self,
        kc: &Knitro,
        values: Option<&[f64]>,
        mode: OptimizationMode,
    ) -> Result<Vec<i32>, KnitroError> {
        use knitro_rs::KN_VARTYPE_INTEGER;

        let index_vars = self
            .0
            .into_iter()
            .map(|v| v.setup_knitro(kc))
            .collect::<Result<Vec<_>, KnitroError>>()?;
        if index_vars.is_empty() {
            return Ok(index_vars);
        }

        if let (_, OptimizationMode::Gradients) | (true, OptimizationMode::FixedMolecule) =
            (FIXABLE, mode)
        {
            kc.set_var_fxbnds(&index_vars, values.unwrap())?;
            kc.set_var_primal_initial_values(&index_vars, values.unwrap())?;
        } else if let Some(values) = values {
            kc.set_var_primal_initial_values(&index_vars, values)?;
        }

        if let (true, OptimizationMode::MolecularDesign) = (INTEGER, mode) {
            kc.set_var_types(&index_vars, &vec![KN_VARTYPE_INTEGER; index_vars.len()])?;
        }

        Ok(index_vars)
    }
}

impl<const INTEGER: bool, const FIXABLE: bool> Variable<INTEGER, FIXABLE> {
    #[cfg(feature = "knitro_rs")]
    pub fn setup_knitro(&self, kc: &Knitro) -> Result<i32, KnitroError> {
        let y = kc.add_var()?;
        if let Some(lobnd) = self.lobnd {
            kc.set_var_lobnd(y, lobnd)?;
        }
        if let Some(upbnd) = self.upbnd {
            kc.set_var_upbnd(y, upbnd)?;
        }
        Ok(y)
    }
}

pub struct LinearConstraint {
    pub vars: Vec<i32>,
    pub coefs: Vec<f64>,
    pub lobnd: Option<f64>,
    pub upbnd: Option<f64>,
}

impl LinearConstraint {
    pub fn new(vars: Vec<i32>, coefs: Vec<f64>) -> Self {
        Self {
            vars,
            coefs,
            lobnd: None,
            upbnd: None,
        }
    }

    pub fn lobnd(mut self, lobnd: f64) -> Self {
        self.lobnd = Some(lobnd);
        self
    }

    pub fn upbnd(mut self, upbnd: f64) -> Self {
        self.upbnd = Some(upbnd);
        self
    }

    pub fn eqbnd(self, eqbnd: f64) -> Self {
        self.lobnd(eqbnd).upbnd(eqbnd)
    }
}
