#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

use std::error::Error;
use std::fmt;
use std::slice;

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

pub struct Knitro(KN_context_ptr);

impl Drop for Knitro {
    fn drop(&mut self) {
        unsafe {
            KN_free(&mut self.0);
        }
    }
}

pub struct KnitroCallback(CB_context_ptr);

#[derive(Debug)]
pub struct KnitroError(String, i32);

impl Error for KnitroError {}

impl fmt::Display for KnitroError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Knitro function {} returned code {}.", self.0, self.1)
    }
}

impl Knitro {
    fn handle_error(from: &str, res: i32) -> Result<(), KnitroError> {
        match res {
            0 | -403 => Ok(()),
            -199..=-100 => {
                println!("Warning: a feasible approximate solution was found.");
                Ok(())
            }
            _ => Err(KnitroError(from.into(), res)),
        }
    }

    pub fn new() -> Result<Self, KnitroError> {
        let mut kc: KN_context_ptr = &mut KN_context { _unused: [0; 0] };
        unsafe {
            Self::handle_error("KN_new", KN_new(&mut kc))?;
        }
        Ok(Self(kc))
    }

    pub fn add_var(&self) -> Result<i32, KnitroError> {
        let mut indexVar = 0;
        unsafe {
            Self::handle_error("KN_add_var", KN_add_var(self.0, &mut indexVar))?;
        }
        Ok(indexVar)
    }

    pub fn add_vars(&self, nV: usize) -> Result<Vec<i32>, KnitroError> {
        let mut indexVars = vec![0; nV];
        unsafe {
            Self::handle_error(
                "KN_add_vars",
                KN_add_vars(self.0, nV as i32, indexVars.as_mut_ptr()),
            )?;
        }
        Ok(indexVars)
    }

    pub fn set_var_type(&self, indexVar: i32, xType: u32) -> Result<(), KnitroError> {
        unsafe {
            Self::handle_error(
                "KN_set_var_type",
                KN_set_var_type(self.0, indexVar, xType as i32),
            )?;
        }
        Ok(())
    }

    pub fn set_var_types(&self, indexVars: &[i32], xTypes: &[u32]) -> Result<(), KnitroError> {
        let xTypes: Vec<_> = xTypes.iter().map(|&xType| xType as i32).collect();
        unsafe {
            Self::handle_error(
                "KN_set_var_types",
                KN_set_var_types(
                    self.0,
                    indexVars.len() as i32,
                    indexVars.as_ptr(),
                    xTypes.as_ptr(),
                ),
            )?;
        }
        Ok(())
    }

    pub fn set_var_lobnd(&self, indexVar: i32, xLoBnd: f64) -> Result<(), KnitroError> {
        unsafe {
            Self::handle_error(
                "KN_set_var_lobnd",
                KN_set_var_lobnd(self.0, indexVar, xLoBnd),
            )?;
        }
        Ok(())
    }

    pub fn set_var_upbnd(&self, indexVar: i32, xUpBnd: f64) -> Result<(), KnitroError> {
        unsafe {
            Self::handle_error(
                "KN_set_var_upbnd",
                KN_set_var_upbnd(self.0, indexVar, xUpBnd),
            )?;
        }
        Ok(())
    }

    pub fn set_var_fxbnd(&self, indexVar: i32, xFxBnd: f64) -> Result<(), KnitroError> {
        unsafe {
            Self::handle_error(
                "KN_set_var_fxbnd",
                KN_set_var_fxbnd(self.0, indexVar, xFxBnd),
            )?;
        }
        Ok(())
    }

    pub fn set_var_fxbnds(&self, indexVars: &[i32], xFxBnds: &[f64]) -> Result<(), KnitroError> {
        unsafe {
            Self::handle_error(
                "KN_set_var_fxbnds",
                KN_set_var_fxbnds(
                    self.0,
                    indexVars.len() as i32,
                    indexVars.as_ptr(),
                    xFxBnds.as_ptr(),
                ),
            )?;
        }
        Ok(())
    }

    pub fn set_var_primal_init_value(
        &self,
        indexVar: i32,
        xInitVal: f64,
    ) -> Result<(), KnitroError> {
        unsafe {
            Self::handle_error(
                "KN_set_var_primal_init_values",
                KN_set_var_primal_init_value(self.0, indexVar, xInitVal),
            )?;
        }
        Ok(())
    }

    pub fn set_var_primal_init_values(
        &self,
        indexVars: &[i32],
        xInitVals: &[f64],
    ) -> Result<(), KnitroError> {
        unsafe {
            Self::handle_error(
                "KN_set_var_primal_init_values",
                KN_set_var_primal_init_values(
                    self.0,
                    indexVars.len() as i32,
                    indexVars.as_ptr(),
                    xInitVals.as_ptr(),
                ),
            )?;
        }
        Ok(())
    }

    pub fn set_var_name(&self, indexVar: i32, cName: &str) -> Result<(), KnitroError> {
        let cName = std::ffi::CString::new(cName).unwrap().into_bytes_with_nul();
        let cName = cName.as_ptr() as *mut i8;
        unsafe {
            Self::handle_error("KN_set_var_names", KN_set_var_name(self.0, indexVar, cName))?;
        }
        Ok(())
    }

    pub fn get_number_vars(&self) -> Result<usize, KnitroError> {
        let mut nV = 0;
        unsafe {
            Self::handle_error("KN_get_number_vars", KN_get_number_vars(self.0, &mut nV))?;
        }
        Ok(nV as usize)
    }

    pub fn add_cons(&self, nC: usize) -> Result<Vec<i32>, KnitroError> {
        let mut indexCon = vec![0; nC];
        unsafe {
            Self::handle_error(
                "KN_add_cons",
                KN_add_cons(self.0, nC as i32, indexCon.as_mut_ptr()),
            )?;
        }
        Ok(indexCon)
    }

    pub fn add_con(&self) -> Result<i32, KnitroError> {
        let mut indexCon = 0;
        unsafe {
            Self::handle_error("KN_add_con", KN_add_con(self.0, &mut indexCon))?;
        }
        Ok(indexCon)
    }

    pub fn add_con_linear_struct_one(
        &self,
        indexCon: i32,
        indexVars: &[i32],
        coefs: &[f64],
    ) -> Result<(), KnitroError> {
        unsafe {
            Self::handle_error(
                "KN_add_con_linear_struct_one",
                KN_add_con_linear_struct_one(
                    self.0,
                    indexVars.len() as i64,
                    indexCon,
                    indexVars.as_ptr(),
                    coefs.as_ptr(),
                ),
            )?;
        }
        Ok(())
    }

    pub fn add_con_quadratic_struct_one(
        &self,
        indexCon: i32,
        indexVars1: &[i32],
        indexVars2: &[i32],
        coefs: &[f64],
    ) -> Result<(), KnitroError> {
        unsafe {
            Self::handle_error(
                "KN_add_con_quadratic_struct_one",
                KN_add_con_quadratic_struct_one(
                    self.0,
                    indexVars1.len() as i64,
                    indexCon,
                    indexVars1.as_ptr(),
                    indexVars2.as_ptr(),
                    coefs.as_ptr(),
                ),
            )?;
        }
        Ok(())
    }

    #[cfg(feature = "knitro_13")]
    pub fn add_con_linear_term(
        &self,
        indexCon: i32,
        indexVar: i32,
        coef: f64,
    ) -> Result<(), KnitroError> {
        unsafe {
            Self::handle_error(
                "KN_add_con_linear_term",
                KN_add_con_linear_term(self.0, indexCon, indexVar, coef),
            )?;
        }
        Ok(())
    }

    #[cfg(feature = "knitro_12")]
    pub fn add_con_linear_term(
        &self,
        indexCon: i32,
        indexVar: i32,
        coef: f64,
    ) -> Result<(), KnitroError> {
        let indexCons = [indexCon];
        let indexVars = [indexVar];
        let coefs = [coef];
        unsafe {
            Self::handle_error(
                "KN_add_con_linear_struct",
                KN_add_con_linear_struct(
                    self.0,
                    1,
                    indexCons.as_ptr(),
                    indexVars.as_ptr(),
                    coefs.as_ptr(),
                ),
            )?;
        }
        Ok(())
    }

    pub fn set_con_lobnd(&self, indexCon: i32, cLoBnd: f64) -> Result<(), KnitroError> {
        unsafe {
            Self::handle_error(
                "KN_set_con_lobnd",
                KN_set_con_lobnd(self.0, indexCon, cLoBnd),
            )?;
        }
        Ok(())
    }

    pub fn set_con_upbnd(&self, indexCon: i32, cUpBnd: f64) -> Result<(), KnitroError> {
        unsafe {
            Self::handle_error(
                "KN_set_con_upbnd",
                KN_set_con_upbnd(self.0, indexCon, cUpBnd),
            )?;
        }
        Ok(())
    }

    pub fn set_con_eqbnd(&self, indexCon: i32, cEqBnd: f64) -> Result<(), KnitroError> {
        unsafe {
            Self::handle_error(
                "KN_set_con_eqbnd",
                KN_set_con_eqbnd(self.0, indexCon, cEqBnd),
            )?;
        }
        Ok(())
    }

    pub fn get_number_cons(&self) -> Result<usize, KnitroError> {
        let mut nC = 0;
        unsafe {
            Self::handle_error("KN_get_number_cons", KN_get_number_cons(self.0, &mut nC))?;
        }
        Ok(nC as usize)
    }

    pub fn set_compcons(
        &self,
        indexComps1: &[i32],
        indexComps2: &[i32],
    ) -> Result<(), KnitroError> {
        unsafe {
            Self::handle_error(
                "KN_set_compcons",
                KN_set_compcons(
                    self.0,
                    indexComps1.len() as i32,
                    vec![KN_CCTYPE_VARVAR as i32; indexComps1.len()].as_ptr(),
                    indexComps1.as_ptr(),
                    indexComps2.as_ptr(),
                ),
            )?;
        }
        Ok(())
    }

    pub fn add_obj_linear_struct(
        &self,
        indexVars: &[i32],
        coefs: &[f64],
    ) -> Result<(), KnitroError> {
        unsafe {
            Self::handle_error(
                "KN_add_obj_linear_struct",
                KN_add_obj_linear_struct(
                    self.0,
                    indexVars.len() as i32,
                    indexVars.as_ptr(),
                    coefs.as_ptr(),
                ),
            )?;
        }
        Ok(())
    }

    pub fn add_eval_callback<C: EvalCallback>(
        &self,
        evalObj: bool,
        indexCons: &[i32],
        callback: &mut C,
    ) -> Result<KnitroCallback, KnitroError> {
        let mut cb: CB_context_ptr = &mut CB_context { _unused: [0; 0] };
        unsafe {
            Self::handle_error(
                "KN_add_eval_callback",
                KN_add_eval_callback(
                    self.0,
                    evalObj as i32,
                    indexCons.len() as i32,
                    indexCons.as_ptr(),
                    get_eval_callback(callback),
                    &mut cb,
                ),
            )?;
        }

        let callback = callback as *mut C;
        unsafe {
            Self::handle_error(
                "KN_set_cb_user_params",
                KN_set_cb_user_params(self.0, cb, callback as *mut std::os::raw::c_void),
            )?;
        }
        Ok(KnitroCallback(cb))
    }

    pub fn set_cb_grad<C: EvalCallback>(
        &self,
        cb: &mut KnitroCallback,
        objGradIndexVars: &[i32],
        jacIndexCons: &[i32],
        jacIndexVars: &[i32],
        callback: Option<&mut C>,
    ) -> Result<(), KnitroError> {
        let callback = callback.map_or_else(|| None, |c| get_eval_callback(c));
        unsafe {
            Self::handle_error(
                "KN_set_cb_grad",
                KN_set_cb_grad(
                    self.0,
                    cb.0,
                    objGradIndexVars.len() as i32,
                    objGradIndexVars.as_ptr(),
                    jacIndexCons.len() as i64,
                    jacIndexCons.as_ptr(),
                    jacIndexVars.as_ptr(),
                    callback,
                ),
            )?;
        }
        Ok(())
    }

    pub fn load_param_file(&self, filename: &str) -> Result<(), KnitroError> {
        let filename = std::ffi::CString::new(filename)
            .unwrap()
            .into_bytes_with_nul();
        let filename = filename.as_ptr() as *const i8;
        unsafe {
            Self::handle_error("KN_load_param_file", KN_load_param_file(self.0, filename))?;
        }
        Ok(())
    }

    pub fn solve(&self) -> Result<(), KnitroError> {
        unsafe {
            Self::handle_error("KN_solve", KN_solve(self.0))?;
        }
        Ok(())
    }

    pub fn get_obj_value(&self) -> Result<f64, KnitroError> {
        let mut obj = 0.0;
        unsafe {
            Self::handle_error("KN_get_obj_value", KN_get_obj_value(self.0, &mut obj))?;
        }
        Ok(obj)
    }

    pub fn get_var_primal_values_all(&self) -> Result<Vec<f64>, KnitroError> {
        let mut x = vec![0.0; self.get_number_vars()?];
        unsafe {
            Self::handle_error(
                "KN_get_var_primal_values_all",
                KN_get_var_primal_values_all(self.0, x.as_mut_ptr()),
            )?;
        }
        Ok(x)
    }

    pub fn get_var_primal_values(&self, indexVars: &[i32]) -> Result<Vec<f64>, KnitroError> {
        let mut x = vec![0.0; indexVars.len()];
        unsafe {
            Self::handle_error(
                "KN_get_var_primal_values",
                KN_get_var_primal_values(
                    self.0,
                    indexVars.len() as i32,
                    indexVars.as_ptr(),
                    x.as_mut_ptr(),
                ),
            )?;
        }
        Ok(x)
    }

    pub fn get_con_values(&self, indexCons: &[i32]) -> Result<Vec<f64>, KnitroError> {
        let mut c = vec![0.0; indexCons.len()];
        unsafe {
            Self::handle_error(
                "KN_get_con_values",
                KN_get_con_values(
                    self.0,
                    indexCons.len() as i32,
                    indexCons.as_ptr(),
                    c.as_mut_ptr(),
                ),
            )?;
        }
        Ok(c)
    }

    pub fn get_con_dual_values(&self, indexCons: &[i32]) -> Result<Vec<f64>, KnitroError> {
        let mut lambda = vec![0.0; indexCons.len()];
        unsafe {
            Self::handle_error(
                "KN_get_con_dual_values",
                KN_get_con_dual_values(
                    self.0,
                    indexCons.len() as i32,
                    indexCons.as_ptr(),
                    lambda.as_mut_ptr(),
                ),
            )?;
        }
        Ok(lambda)
    }

    pub fn get_objgrad_values_all(&self) -> Result<Vec<f64>, KnitroError> {
        let mut objGrad = vec![0.0; self.get_number_vars()?];
        unsafe {
            Self::handle_error(
                "KN_get_objgrad_values_all",
                KN_get_objgrad_values_all(self.0, objGrad.as_mut_ptr()),
            )?;
        }
        Ok(objGrad)
    }

    pub fn get_objgrad_values(&self) -> Result<(Vec<i32>, Vec<f64>), KnitroError> {
        let mut nnz = 0;
        unsafe {
            Self::handle_error("KN_get_objgrad_nnz", KN_get_objgrad_nnz(self.0, &mut nnz))?;
        }
        let nnz = nnz as usize;
        let mut indexVars = vec![0; nnz];
        let mut objGrad = vec![0.0; nnz];
        unsafe {
            Self::handle_error(
                "KN_get_objgrad_values",
                KN_get_objgrad_values(self.0, indexVars.as_mut_ptr(), objGrad.as_mut_ptr()),
            )?;
        }
        Ok((indexVars, objGrad))
    }

    pub fn get_jacobian_values(&self) -> Result<(Vec<i32>, Vec<i32>, Vec<f64>), KnitroError> {
        let mut nnz = 0;
        unsafe {
            Self::handle_error("KN_get_jacobian_nnz", KN_get_jacobian_nnz(self.0, &mut nnz))?;
        }
        let nnz = nnz as usize;
        let mut indexCons = vec![0; nnz];
        let mut indexVars = vec![0; nnz];
        let mut jac = vec![0.0; nnz];
        unsafe {
            Self::handle_error(
                "KN_get_jacobian_values",
                KN_get_jacobian_values(
                    self.0,
                    indexCons.as_mut_ptr(),
                    indexVars.as_mut_ptr(),
                    jac.as_mut_ptr(),
                ),
            )?;
        }
        Ok((indexCons, indexVars, jac))
    }
}

pub trait EvalCallback {
    fn callback(&self, x: &[f64], obj: &mut f64, c: &mut [f64]) -> i32;
}

unsafe extern "C" fn eval_callback_wrapper<C: EvalCallback>(
    kc: KN_context_ptr,
    cb: CB_context_ptr,
    evalRequest: KN_eval_request_ptr,
    evalResult: KN_eval_result_ptr,
    userParams: *mut ::std::os::raw::c_void,
) -> i32 {
    let userParams = &mut *(userParams as *mut C);
    let mut nV = 0;
    let res = KN_get_number_vars(kc, &mut nV);
    if res != 0 {
        return res;
    }
    let mut nC = 0;
    let res = KN_get_cb_number_cons(kc, cb, &mut nC);
    if res != 0 {
        return res;
    }
    let x = slice::from_raw_parts((*evalRequest).x, nV as usize);
    let c = slice::from_raw_parts_mut((*evalResult).c, nC as usize);
    let obj = &mut *(*evalResult).obj;
    userParams.callback(x, obj, c)
}

pub fn get_eval_callback<C: EvalCallback>(_: &C) -> KN_eval_callback {
    Some(eval_callback_wrapper::<C>)
}
