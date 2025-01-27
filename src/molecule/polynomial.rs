use num_dual::DualNum;
use std::iter::Product;
use std::ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign};
use std::slice::Iter;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Polynomial2<D, const M: usize, const N: usize>([[D; N]; M]);
pub type Polynomial<D, const N: usize> = Polynomial2<D, 1, N>;

impl<D: DualNum<f64> + Copy> Polynomial<D, 2> {
    pub fn new(p: D) -> Self {
        Polynomial2([[D::from(1.0) - p, p]])
    }
}

impl<D: DualNum<f64> + Copy, const M: usize, const N: usize> Polynomial2<D, M, N> {
    pub fn zero() -> Self {
        Polynomial2([[D::zero(); N]; M])
    }
}

impl<D: DualNum<f64> + Copy, const N: usize> Polynomial<D, N> {
    fn one() -> Self {
        let mut one = Self::zero();
        one[0] = D::one();
        one
    }

    pub fn sum(&self) -> D {
        self.0[0].into_iter().sum()
    }

    pub fn iter(&self) -> Iter<D> {
        self.0[0].iter()
    }
}

impl<D: DualNum<f64> + Copy, const N: usize> Index<usize> for Polynomial<D, N> {
    type Output = D;
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[0][index]
    }
}

impl<D: DualNum<f64> + Copy, const N: usize> IndexMut<usize> for Polynomial<D, N> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[0][index]
    }
}

impl<D: DualNum<f64> + Copy, const M: usize, const N: usize> Index<(usize, usize)>
    for Polynomial2<D, M, N>
{
    type Output = D;
    fn index(&self, index: (usize, usize)) -> &Self::Output {
        let (i, j) = index;
        &self.0[i][j]
    }
}

impl<D: DualNum<f64> + Copy, const M: usize, const N: usize> Add for Polynomial2<D, M, N> {
    type Output = Self;
    fn add(mut self, other: Self) -> Self {
        self += other;
        self
    }
}

impl<D: DualNum<f64> + Copy, const M: usize, const N: usize> AddAssign for Polynomial2<D, M, N> {
    fn add_assign(&mut self, other: Self) {
        for i in 0..M {
            for j in 0..N {
                self.0[i][j] += other.0[i][j]
            }
        }
    }
}

impl<D: DualNum<f64> + Copy, const N: usize> Mul<Polynomial<D, 2>> for Polynomial<D, N> {
    type Output = Self;
    fn mul(self, other: Polynomial<D, 2>) -> Self {
        let [[b0, b1]] = other.0;
        let mut res = self * b0;
        for i in 1..N {
            res[i] += self[i - 1] * b1;
        }
        res
    }
}

impl<D: DualNum<f64> + Copy, const M: usize, const N: usize> Mul<D> for Polynomial2<D, M, N> {
    type Output = Self;
    fn mul(mut self, other: D) -> Self {
        self *= other;
        self
    }
}

impl<D: DualNum<f64> + Copy, const M: usize, const N: usize> MulAssign<D> for Polynomial2<D, M, N> {
    fn mul_assign(&mut self, other: D) {
        for i in 0..M {
            for j in 0..N {
                self.0[i][j] *= other
            }
        }
    }
}

impl<D: DualNum<f64> + Copy, const N: usize> Product<Polynomial<D, 2>> for Polynomial<D, N> {
    fn product<I: Iterator<Item = Polynomial<D, 2>>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, a| acc * a)
    }
}

impl<D: DualNum<f64> + Copy, const M: usize> Polynomial<D, M> {
    pub fn outer_product<const N: usize>(self, other: Polynomial<D, N>) -> Polynomial2<D, M, N> {
        let mut mat = [[D::zero(); N]; M];
        for i in 0..M {
            for j in 0..N {
                mat[i][j] = self[i] * other[j]
            }
        }
        Polynomial2(mat)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_mul() {
        let p1 = Polynomial2([[1.0, 2.0]]);
        let p2 = Polynomial2([[4.0, 5.0, 6.0]]);
        assert_eq!(p2 * p1, Polynomial2([[4.0, 13.0, 16.0]]));
    }
}
