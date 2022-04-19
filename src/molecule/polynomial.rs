use std::iter::Product;
use std::ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign};
use std::slice::Iter;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Polynomial2<const M: usize, const N: usize>([[f64; N]; M]);
pub type Polynomial<const N: usize> = Polynomial2<1, N>;

impl Polynomial<2> {
    pub fn new(p: f64) -> Self {
        Polynomial2([[1.0 - p, p]])
    }
}
impl<const M: usize, const N: usize> Polynomial2<M, N> {
    pub fn zero() -> Self {
        Polynomial2([[0.0; N]; M])
    }
}

impl<const N: usize> Polynomial<N> {
    fn one() -> Self {
        let mut one = Self::zero();
        one[0] = 1.0;
        one
    }

    pub fn sum(&self) -> f64 {
        self.0[0].iter().sum()
    }

    pub fn iter(&self) -> Iter<f64> {
        self.0[0].iter()
    }
}

impl<const N: usize> Index<usize> for Polynomial<N> {
    type Output = f64;
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[0][index]
    }
}

impl<const N: usize> IndexMut<usize> for Polynomial<N> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[0][index]
    }
}

impl<const M: usize, const N: usize> Index<(usize, usize)> for Polynomial2<M, N> {
    type Output = f64;
    fn index(&self, index: (usize, usize)) -> &Self::Output {
        let (i, j) = index;
        &self.0[i][j]
    }
}

impl<const M: usize, const N: usize> Add for Polynomial2<M, N> {
    type Output = Self;
    fn add(mut self, other: Self) -> Self {
        self += other;
        self
    }
}

impl<const M: usize, const N: usize> AddAssign for Polynomial2<M, N> {
    fn add_assign(&mut self, other: Self) {
        for i in 0..M {
            for j in 0..N {
                self.0[i][j] += other.0[i][j]
            }
        }
    }
}

impl<const N: usize> Mul<Polynomial<2>> for Polynomial<N> {
    type Output = Self;
    fn mul(self, other: Polynomial<2>) -> Self {
        let [[b0, b1]] = other.0;
        let mut res = self * b0;
        for i in 1..N {
            res[i] += self[i - 1] * b1;
        }
        res
    }
}

impl<const M: usize, const N: usize> Mul<f64> for Polynomial2<M, N> {
    type Output = Self;
    fn mul(mut self, other: f64) -> Self {
        self *= other;
        self
    }
}

impl<const M: usize, const N: usize> MulAssign<f64> for Polynomial2<M, N> {
    fn mul_assign(&mut self, other: f64) {
        for i in 0..M {
            for j in 0..N {
                self.0[i][j] *= other
            }
        }
    }
}

impl<const N: usize> Product<Polynomial<2>> for Polynomial<N> {
    fn product<I: Iterator<Item = Polynomial<2>>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, a| acc * a)
    }
}

impl<const M: usize> Polynomial<M> {
    pub fn outer_product<const N: usize>(self, other: Polynomial<N>) -> Polynomial2<M, N> {
        let mut mat = [[0.0; N]; M];
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
