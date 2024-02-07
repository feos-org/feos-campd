use std::ops::{Deref, DerefMut};

#[cfg(feature = "knitro_rs")]
use crate::OptimizationMode;
#[cfg(feature = "knitro_rs")]
use knitro_rs::{Knitro, KnitroError};

#[derive(Clone, Copy)]
#[allow(dead_code)]
pub struct Variable<const INTEGER: bool> {
    lobnd: Option<f64>,
    upbnd: Option<f64>,
}

pub struct Variables<const INTEGER: bool>(Vec<Variable<INTEGER>>);

impl<const INTEGER: bool> Deref for Variables<INTEGER> {
    type Target = Vec<Variable<INTEGER>>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<const INTEGER: bool> DerefMut for Variables<INTEGER> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<const INTEGER: bool> FromIterator<Variable<INTEGER>> for Variables<INTEGER> {
    fn from_iter<T: IntoIterator<Item = Variable<INTEGER>>>(iter: T) -> Self {
        Self(Vec::from_iter(iter))
    }
}

impl<const INTEGER: bool> From<Vec<Variable<INTEGER>>> for Variables<INTEGER> {
    fn from(value: Vec<Variable<INTEGER>>) -> Self {
        Self(value)
    }
}

pub type ContinuousVariables = Variables<false>;
pub type DiscreteVariables = Variables<true>;

impl<const INTEGER: bool> Variable<INTEGER> {
    fn new(lobnd: Option<f64>, upbnd: Option<f64>) -> Self {
        Self { lobnd, upbnd }
    }
}

impl Variable<true> {
    pub fn binary() -> Self {
        Self::integer(0, 1)
    }

    pub fn integer(lobnd: usize, upbnd: usize) -> Self {
        let lobnd = lobnd as f64;
        let upbnd = upbnd as f64;
        Self::new(Some(lobnd), Some(upbnd))
    }
}

impl Variable<false> {
    pub fn free() -> Self {
        Self::continuous(None, None)
    }

    pub fn continuous(lobnd: Option<f64>, upbnd: Option<f64>) -> Self {
        Self::new(lobnd, upbnd)
    }
}

impl<const INTEGER: bool> Variables<INTEGER> {
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
            (INTEGER, mode)
        {
            kc.set_var_fxbnds(&index_vars, values.unwrap())?;
        } else if let Some(values) = values {
            kc.set_var_primal_initial_values(&index_vars, values)?;
        }

        if let (true, OptimizationMode::MolecularDesign) = (INTEGER, mode) {
            kc.set_var_types(&index_vars, &vec![KN_VARTYPE_INTEGER; index_vars.len()])?;
        }

        Ok(index_vars)
    }
}

impl<const INTEGER: bool> Variable<INTEGER> {
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
