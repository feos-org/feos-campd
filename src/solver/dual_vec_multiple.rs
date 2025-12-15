use nalgebra::allocator::Allocator;
use nalgebra::*;
use num_dual::*;
use num_traits::{Float, FloatConst, FromPrimitive, Inv, Num, One, Signed, Zero};
use std::fmt;
use std::iter::{Product, Sum};
use std::marker::PhantomData;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, RemAssign, Sub, SubAssign,
};

/// A vector dual number for the calculations of gradients or Jacobians.
#[derive(Clone, Debug)]
pub struct DualVecMult<T: DualNum<F>, F, N1: Dim, N2: Dim, N: Dim>
where
    DefaultAllocator: Allocator<N1, N> + Allocator<N2>,
{
    /// Real part of the dual number
    pub re: T,
    /// Derivative part of the dual number
    pub eps1: Derivative<T, F, N1, N>,
    /// Derivative part of the dual number
    pub eps2: Derivative<T, F, N2, U1>,
    f: PhantomData<F>,
}

pub type DualSVecMult<T, F, const N1: usize, const N2: usize, const N: usize> =
    DualVecMult<T, F, Const<N1>, Const<N2>, Const<N>>;

impl<T: DualNum<F> + Copy, F: Copy, const N1: usize, const N2: usize, const N: usize> Copy
    for DualSVecMult<T, F, N1, N2, N>
{
}

impl<T: DualNum<F>, F, N1: Dim, N2: Dim, N: Dim> DualVecMult<T, F, N1, N2, N>
where
    DefaultAllocator: Allocator<N1, N> + Allocator<N2>,
{
    /// Create a new dual number from its fields.
    #[inline]
    pub fn new(re: T, eps1: Derivative<T, F, N1, N>, eps2: Derivative<T, F, N2, U1>) -> Self {
        Self {
            re,
            eps1,
            eps2,
            f: PhantomData,
        }
    }
}

impl<T: DualNum<F> + Zero, F, N1: Dim, N2: Dim, N: Dim> DualVecMult<T, F, N1, N2, N>
where
    DefaultAllocator: Allocator<N1, N> + Allocator<N2>,
{
    /// Create a new dual number from the real part.
    #[inline]
    pub fn from_re(re: T) -> Self {
        Self::new(re, Derivative::none(), Derivative::none())
    }
}

/* chain rule */
impl<T: DualNum<F>, F: Float, N1: Dim, N2: Dim, N: Dim> DualVecMult<T, F, N1, N2, N>
where
    DefaultAllocator: Allocator<N1, N> + Allocator<N2>,
{
    #[inline]
    fn chain_rule(&self, f0: T, f1: T) -> Self {
        Self::new(f0, &self.eps1 * f1.clone(), &self.eps2 * f1)
    }
}

/* product rule */
impl<T: DualNum<F>, F: Float, N1: Dim, N2: Dim, N: Dim> Mul<&DualVecMult<T, F, N1, N2, N>>
    for &DualVecMult<T, F, N1, N2, N>
where
    DefaultAllocator: Allocator<N1, N> + Allocator<N2>,
{
    type Output = DualVecMult<T, F, N1, N2, N>;
    #[inline]
    fn mul(self, other: &DualVecMult<T, F, N1, N2, N>) -> Self::Output {
        DualVecMult::new(
            self.re.clone() * other.re.clone(),
            &self.eps1 * other.re.clone() + &other.eps1 * self.re.clone(),
            &self.eps2 * other.re.clone() + &other.eps2 * self.re.clone(),
        )
    }
}

/* quotient rule */
impl<T: DualNum<F>, F: Float, N1: Dim, N2: Dim, N: Dim> Div<&DualVecMult<T, F, N1, N2, N>>
    for &DualVecMult<T, F, N1, N2, N>
where
    DefaultAllocator: Allocator<N1, N> + Allocator<N2>,
{
    type Output = DualVecMult<T, F, N1, N2, N>;
    #[inline]
    fn div(self, other: &DualVecMult<T, F, N1, N2, N>) -> DualVecMult<T, F, N1, N2, N> {
        let inv = other.re.recip();
        let inv2 = inv.clone() * inv.clone();
        DualVecMult::new(
            self.re.clone() * inv,
            (&self.eps1 * other.re.clone() - &other.eps1 * self.re.clone()) * inv2.clone(),
            (&self.eps2 * other.re.clone() - &other.eps2 * self.re.clone()) * inv2,
        )
    }
}

/* string conversions */
impl<T: DualNum<F>, F, N1: Dim, N2: Dim, N: Dim> fmt::Display for DualVecMult<T, F, N1, N2, N>
where
    DefaultAllocator: Allocator<N1, N> + Allocator<N2>,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.re)?;
        self.eps1.fmt(f, "ε1")?;
        self.eps2.fmt(f, "ε2")
    }
}

impl_first_derivatives!(DualVecMult, [eps1, eps2], [N1, N2, N], [N1, N], [N2]);
impl_dual!(DualVecMult, [eps1, eps2], [N1, N2, N], [N1, N], [N2]);

/// Comparisons are only made based on the real part. This allows the code to follow the
/// same execution path as real-valued code would.
impl<T: DualNum<F> + PartialEq, F: Float, N1: Dim, N2: Dim, N: Dim> PartialEq
    for DualVecMult<T, F, N1, N2, N>
where
    DefaultAllocator: Allocator<N1, N> + Allocator<N2>,
{
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.re.eq(&other.re)
    }
}
