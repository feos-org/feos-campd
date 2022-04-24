#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

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

pub struct KnitroError(String, i32);

impl fmt::Debug for KnitroError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Knitro function {} returned code {}.", self.0, self.1)
    }
}

impl Knitro {
    fn handle_error(from: &str, res: i32) -> Result<(), KnitroError> {
        match res {
            0 => Ok(()),
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

    pub fn add_vars(&self, nV: usize, xType: Option<u32>) -> Result<Vec<i32>, KnitroError> {
        let mut indexVars = vec![0; nV];
        unsafe {
            Self::handle_error(
                "KN_add_vars",
                KN_add_vars(self.0, nV as i32, indexVars.as_mut_ptr()),
            )?;
        }
        if let Some(xType) = xType {
            let xType = vec![xType as i32; nV];
            unsafe {
                Self::handle_error(
                    "KN_set_var_types",
                    KN_set_var_types(self.0, nV as i32, indexVars.as_ptr(), xType.as_ptr()),
                )?;
            }
        }
        Ok(indexVars)
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

    pub fn set_var_primal_initial_values(
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
                "KN_add_con_linear_term",
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

    pub fn add_eval_callback<C: Callback>(
        &self,
        evalObj: bool,
        indexCons: &[i32],
        callback: &C,
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
                    get_callback(callback),
                    &mut cb,
                ),
            )?;
        }
        Ok(KnitroCallback(cb))
    }

    pub fn set_cb_user_params<C: Callback>(
        &self,
        cb: KnitroCallback,
        userParams: &mut C,
    ) -> Result<(), KnitroError> {
        let userParams = userParams as *mut C;
        unsafe {
            Self::handle_error(
                "KN_set_cb_user_params",
                KN_set_cb_user_params(self.0, cb.0, userParams as *mut std::os::raw::c_void),
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
}

pub trait Callback {
    fn callback(&self, x: &[f64], c: &mut [f64]) -> f64;
}

unsafe extern "C" fn callback_wrapper<C: Callback>(
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
    let obj = userParams.callback(x, c);
    *(*evalResult).obj = obj;
    0
}

pub fn get_callback<C>(_: &C) -> KN_eval_callback
where
    C: Callback,
{
    Some(callback_wrapper::<C>)
}
