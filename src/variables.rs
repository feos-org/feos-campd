use std::ops::{Deref, DerefMut};

use indexmap::IndexMap;
#[cfg(feature = "knitro_rs")]
use knitro_rs::{Knitro, KnitroError};

#[derive(Clone, Copy, Debug)]
#[allow(dead_code)]
pub struct Variable<const INTEGER: bool> {
    lobnd: Option<f64>,
    upbnd: Option<f64>,
    fxbnd: Option<f64>,
    init: Option<f64>,
}

#[derive(Debug)]
pub struct Variables<const INTEGER: bool>(IndexMap<String, Variable<INTEGER>>);

impl<const INTEGER: bool> Deref for Variables<INTEGER> {
    type Target = IndexMap<String, Variable<INTEGER>>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<const INTEGER: bool> DerefMut for Variables<INTEGER> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<const INTEGER: bool> FromIterator<(String, Variable<INTEGER>)> for Variables<INTEGER> {
    fn from_iter<T: IntoIterator<Item = (String, Variable<INTEGER>)>>(iter: T) -> Self {
        Self(IndexMap::from_iter(iter))
    }
}

impl<const INTEGER: bool> FromIterator<Variable<INTEGER>> for Variables<INTEGER> {
    fn from_iter<T: IntoIterator<Item = Variable<INTEGER>>>(iter: T) -> Self {
        let varname = if INTEGER { "y" } else { "x" };
        Self(IndexMap::from_iter(
            iter.into_iter()
                .enumerate()
                .map(|(i, v)| (format!("{varname}{i}"), v)),
        ))
    }
}

impl<const INTEGER: bool> IntoIterator for Variables<INTEGER> {
    type Item = <IndexMap<String, Variable<INTEGER>> as IntoIterator>::Item;
    type IntoIter = <IndexMap<String, Variable<INTEGER>> as IntoIterator>::IntoIter;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<const INTEGER: bool> From<Vec<(&str, Variable<INTEGER>)>> for Variables<INTEGER> {
    fn from(value: Vec<(&str, Variable<INTEGER>)>) -> Self {
        Self(value.into_iter().map(|(k, v)| (k.into(), v)).collect())
    }
}

pub type ProcessVariables = Variables<false>;
pub type StructureVariables = Variables<true>;

impl<const INTEGER: bool> Variable<INTEGER> {
    fn new(lobnd: Option<f64>, upbnd: Option<f64>, fxbnd: Option<f64>, init: Option<f64>) -> Self {
        Self {
            lobnd,
            upbnd,
            fxbnd,
            init,
        }
    }

    pub fn init(mut self, init: f64) -> Self {
        self.init = Some(init);
        self
    }
}

impl Variable<true> {
    pub fn binary() -> Self {
        Self::new(Some(0.0), Some(1.0), None, None)
    }
}

impl Variable<false> {
    pub fn continuous(lobnd: f64, upbnd: f64, init: f64) -> Self {
        Self::new(Some(lobnd), Some(upbnd), None, Some(init))
    }

    pub fn fixed(fxbnd: f64) -> Self {
        Self::new(None, None, Some(fxbnd), None)
    }
}

impl StructureVariables {
    #[cfg(feature = "knitro_rs")]
    pub fn setup_knitro(
        self,
        kc: &Knitro,
        values: Option<&[f64]>,
        relax: bool,
    ) -> Result<Vec<i32>, KnitroError> {
        use knitro_rs::KN_VARTYPE_INTEGER;

        let index_vars = self
            .0
            .into_iter()
            .map(|(n, v)| v.setup_knitro(kc, &n))
            .collect::<Result<Vec<_>, KnitroError>>()?;

        if !index_vars.is_empty() {
            if let Some(values) = values {
                kc.set_var_primal_init_values(&index_vars, values)?;
            }

            if !relax {
                kc.set_var_types(&index_vars, &vec![KN_VARTYPE_INTEGER; index_vars.len()])?;
            }
        }

        Ok(index_vars)
    }
}

impl ProcessVariables {
    #[cfg(feature = "knitro_rs")]
    pub fn setup_knitro(
        self,
        kc: &Knitro,
        values: Option<&[f64]>,
    ) -> Result<Vec<i32>, KnitroError> {
        let index_vars = self
            .0
            .into_iter()
            .map(|(n, v)| v.setup_knitro(kc, &n))
            .collect::<Result<Vec<_>, KnitroError>>()?;

        if !index_vars.is_empty() {
            if let Some(values) = values {
                kc.set_var_primal_init_values(&index_vars, values)?;
            }
        }

        Ok(index_vars)
    }
}

impl<const INTEGER: bool> Variable<INTEGER> {
    #[cfg(feature = "knitro_rs")]
    pub fn setup_knitro(&self, kc: &Knitro, name: &str) -> Result<i32, KnitroError> {
        let y = kc.add_var()?;
        kc.set_var_name(y, name)?;
        if let Some(lobnd) = self.lobnd {
            kc.set_var_lobnd(y, lobnd)?;
        }
        if let Some(upbnd) = self.upbnd {
            kc.set_var_upbnd(y, upbnd)?;
        }
        if let Some(fxbnd) = self.fxbnd {
            kc.set_var_fxbnd(y, fxbnd)?;
        }
        if let Some(init) = self.init {
            kc.set_var_primal_init_value(y, init)?;
        }
        Ok(y)
    }
}

#[derive(Default, Debug)]
pub struct Constraint {
    pub lvars: Vec<i32>,
    pub lcoefs: Vec<f64>,
    pub qvars1: Vec<i32>,
    pub qvars2: Vec<i32>,
    pub qcoefs: Vec<f64>,
    pub lobnd: Option<f64>,
    pub upbnd: Option<f64>,
}

impl Constraint {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn linear_struct(mut self, vars: Vec<i32>, coefs: Vec<f64>) -> Self {
        self.lvars = vars;
        self.lcoefs = coefs;
        self
    }

    pub fn quadratic_struct(mut self, vars1: Vec<i32>, vars2: Vec<i32>, coefs: Vec<f64>) -> Self {
        self.qvars1 = vars1;
        self.qvars2 = vars2;
        self.qcoefs = coefs;
        self
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

    #[cfg(feature = "knitro_rs")]
    pub fn setup_knitro(&self, kc: &Knitro) -> Result<(), KnitroError> {
        let c = kc.add_con()?;
        if !self.lvars.is_empty() {
            kc.add_con_linear_struct_one(c, &self.lvars, &self.lcoefs)?;
        }
        if !self.qvars1.is_empty() {
            kc.add_con_quadratic_struct_one(c, &self.qvars1, &self.qvars2, &self.qcoefs)?;
        }
        if let Some(lobnd) = self.lobnd {
            kc.set_con_lobnd(c, lobnd)?;
        }
        if let Some(upbnd) = self.upbnd {
            kc.set_con_upbnd(c, upbnd)?;
        }
        Ok(())
    }
}

#[derive(Default, Debug)]
pub struct ExplicitVariable {
    pub name: String,
    pub lvars: Vec<i32>,
    pub lcoefs: Vec<f64>,
    pub qvars1: Vec<i32>,
    pub qvars2: Vec<i32>,
    pub qcoefs: Vec<f64>,
    pub cons: f64,
}

impl ExplicitVariable {
    pub fn new(name: String) -> Self {
        Self {
            name,
            ..Default::default()
        }
    }

    pub fn linear_struct(mut self, vars: Vec<i32>, coefs: Vec<f64>) -> Self {
        self.lvars = vars;
        self.lcoefs = coefs;
        self
    }

    pub fn quadratic_struct(mut self, vars1: Vec<i32>, vars2: Vec<i32>, coefs: Vec<f64>) -> Self {
        self.qvars1 = vars1;
        self.qvars2 = vars2;
        self.qcoefs = coefs;
        self
    }

    pub fn cons(mut self, cons: f64) -> Self {
        self.cons = cons;
        self
    }

    pub fn add_cons(&mut self, cons: f64) {
        self.cons += cons;
    }

    pub fn evaluate(&self, vars: &[f64]) -> f64 {
        let lvars = self.lvars.iter().map(|&v| vars[v as usize]);
        let qvars1 = self.qvars1.iter().map(|&v| vars[v as usize]);
        let qvars2 = self.qvars2.iter().map(|&v| vars[v as usize]);
        self.cons
            + lvars.zip(&self.lcoefs).map(|(v, c)| v * c).sum::<f64>()
            + qvars1
                .zip(qvars2)
                .zip(&self.qcoefs)
                .map(|((v1, v2), c)| v1 * v2 * c)
                .sum::<f64>()
    }

    #[cfg(feature = "knitro_rs")]
    pub fn setup_knitro(&self, kc: &Knitro) -> Result<i32, KnitroError> {
        let y = kc.add_var()?;
        let c = kc.add_con()?;
        if !self.lvars.is_empty() {
            kc.add_con_linear_struct_one(c, &self.lvars, &self.lcoefs)?;
        }
        if !self.qvars1.is_empty() {
            kc.add_con_quadratic_struct_one(c, &self.qvars1, &self.qvars2, &self.qcoefs)?;
        }
        kc.add_con_linear_term(c, y, -1.0)?;
        kc.set_con_eqbnd(c, -self.cons)?;

        let lvars = kc.get_var_primal_values(&self.lvars)?;
        let qvars1 = kc.get_var_primal_values(&self.qvars1)?;
        let qvars2 = kc.get_var_primal_values(&self.qvars2)?;
        let y0 = self.cons
            + lvars
                .into_iter()
                .zip(&self.lcoefs)
                .map(|(v, c)| v * c)
                .sum::<f64>()
            + qvars1
                .into_iter()
                .zip(qvars2)
                .zip(&self.qcoefs)
                .map(|((v1, v2), c)| v1 * v2 * c)
                .sum::<f64>();
        kc.set_var_primal_init_value(y, y0)?;
        kc.set_var_name(y, &self.name)?;

        Ok(y)
    }
}
