use std::array;
use std::iter::Sum;
use std::ops::{Add, AddAssign, Mul};

#[derive(Clone, Debug, PartialEq)]
#[allow(clippy::type_complexity)]
pub struct Polynomial2<const M: usize, const N: usize>(
    pub [[(Vec<Vec<isize>>, Vec<isize>, isize); N]; M],
);

#[derive(Clone, Debug, PartialEq)]
pub struct Polynomial<const N: usize>(pub [(Vec<isize>, isize); N]);

#[derive(Clone, Copy, Debug)]
pub enum Variable {
    Const(isize),
    Var(i32),
}

impl<const N: usize> Polynomial<N> {
    pub fn zero(variables: usize) -> Self {
        Self(array::from_fn(|_| (vec![0; variables], 0)))
    }

    pub fn new(variables: usize, p: [Variable; N]) -> Self {
        let mut coefs = Self::zero(variables).0;
        for i in 0..N {
            match p[i] {
                Variable::Const(c) => coefs[i].1 += c,
                Variable::Var(p) => coefs[i].0[p as usize] += 1,
            }
            if i + 1 < N {
                match p[i + 1] {
                    Variable::Const(c) => coefs[i].1 -= c,
                    Variable::Var(p) => coefs[i].0[p as usize] -= 1,
                }
            }
        }
        Self(coefs)
    }

    pub fn upgrade<const X: usize>(self) -> Polynomial<X> {
        let variables = self.0[0].0.len();
        Polynomial(array::from_fn(|i| {
            if i + N >= X {
                self.0[i + N - X].clone()
            } else {
                (vec![0; variables], 0)
            }
        }))
    }
}

impl<const M: usize, const N: usize> Polynomial2<M, N> {
    pub fn zero(variables: usize) -> Self {
        Self(array::from_fn(|_| {
            array::from_fn(|_| (vec![vec![0; variables]; variables], vec![0; variables], 0))
        }))
    }
}

impl<const N: usize> Add for Polynomial<N> {
    type Output = Self;
    fn add(mut self, other: Self) -> Self {
        self += &other;
        self
    }
}

impl<const N: usize> AddAssign<&Polynomial<N>> for Polynomial<N> {
    fn add_assign(&mut self, other: &Self) {
        for i in 0..N {
            for (s, &o) in self.0[i].0.iter_mut().zip(&other.0[i].0) {
                *s += o;
            }
            self.0[i].1 += other.0[i].1;
        }
    }
}

impl<const M: usize, const N: usize> AddAssign for Polynomial2<M, N> {
    fn add_assign(&mut self, other: Self) {
        for i in 0..M {
            for j in 0..N {
                for (s, o) in self.0[i][j].0.iter_mut().zip(&other.0[i][j].0) {
                    for (s, &o) in s.iter_mut().zip(o) {
                        *s += o;
                    }
                }
                for (s, o) in self.0[i][j].1.iter_mut().zip(&other.0[i][j].1) {
                    *s += o;
                }
                self.0[i][j].2 += other.0[i][j].2
            }
        }
    }
}

impl<const N: usize> Sum for Polynomial<N> {
    fn sum<I: Iterator<Item = Polynomial<N>>>(iter: I) -> Self {
        iter.reduce(|a, b| a + b).unwrap()
    }
}

impl<const M: usize, const N: usize> Mul<Polynomial<N>> for Polynomial<M> {
    type Output = Polynomial2<M, N>;
    fn mul(self, other: Polynomial<N>) -> Polynomial2<M, N> {
        Polynomial2(array::from_fn(|i| {
            array::from_fn(|j| {
                (
                    self.0[i]
                        .0
                        .iter()
                        .map(|&a| other.0[j].0.iter().map(|&b| a * b).collect())
                        .collect(),
                    self.0[i]
                        .0
                        .iter()
                        .zip(&other.0[j].0)
                        .map(|(&a, &b)| a * other.0[j].1 + self.0[i].1 * b)
                        .collect(),
                    self.0[i].1 * other.0[j].1,
                )
            })
        }))
    }
}
