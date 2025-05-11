use std::array::from_fn;
use std::ops::{Add, Mul};

use super::F;

pub const fn sum_up_to<const N: usize>() -> usize {
    N * (N + 1) / 2
}

fn sqr(x: F) -> F {
    x * x
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SymMatrix<const N: usize>
where
    [(); sum_up_to::<N>()]:,
{
    pub(crate) data: [F; sum_up_to::<N>()],
}

pub type SymMatrix3 = SymMatrix<3>;

impl<const N: usize> SymMatrix<N>
where
    [(); sum_up_to::<N>()]:,
{
    #[inline]
    pub fn new(data: [F; sum_up_to::<N>()]) -> Self {
        Self { data }
    }
    pub fn zero() -> Self {
        Self { data: [0.; _] }
    }
}

impl<const N: usize> Default for SymMatrix<N>
where
    [(); sum_up_to::<N>()]:,
{
    fn default() -> Self {
        Self::zero()
    }
}

impl SymMatrix3 {
    #[inline]
    pub fn outer([x, y, z]: [F; 3]) -> Self {
        Self::new([x * x, x * y, x * z, y * y, y * z, z * z])
    }
    #[rustfmt::skip]
    pub const SYM_IDX: [[usize;3];3] = [
      [0,1,2],
      [1,3,4],
      [2,4,5],
    ];

    #[inline]
    pub fn v(&self, x: usize, y: usize) -> F {
        self.data[Self::SYM_IDX[x][y]]
    }

    /// Constant version of access
    #[inline]
    pub const fn v_const<const X: usize, const Y: usize>(&self) -> F {
        self.data[Self::SYM_IDX[X][Y]]
    }

    pub fn vec_mul(&self, v: [F; 3]) -> [F; 3] {
        from_fn(|i| (0..3).map(|j| self.v(j, i) * v[j]).sum())
    }
    /// returns [eigenvalues, eigenvectors] of this Symmetric 3x3 matrix.
    #[inline]
    pub fn eigen(&self) -> ([F; 3], [[F; 3]; 3]) {
        let p1 =
            sqr(self.v_const::<0, 1>()) + sqr(self.v_const::<0, 2>()) + sqr(self.v_const::<1, 2>());
        let diag = [
            self.v_const::<0, 0>(),
            self.v_const::<1, 1>(),
            self.v_const::<2, 2>(),
        ];
        const ZERO_EPS: F = 1e-8;
        // diagonal
        if p1.abs() < ZERO_EPS {
            return (diag, [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]);
        }

        let d = self.v_const::<0, 1>();
        let e = self.v_const::<1, 2>();
        let f = self.v_const::<0, 2>();

        assert!(d != 0. || e != 0. || f != 0.);
        let vs = super::svd::eigen_jacobi(*self);
        let es = vs.map(|v| super::dot(v, self.vec_mul(v)));
        (es, vs)
    }
}

impl<const N: usize> Add for SymMatrix<N>
where
    [(); sum_up_to::<N>()]:,
{
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self::new(from_fn(|i| self.data[i] + o.data[i]))
    }
}

impl<const N: usize> Mul<F> for SymMatrix<N>
where
    [(); sum_up_to::<N>()]:,
{
    type Output = Self;
    fn mul(self, o: F) -> Self {
        Self::new(from_fn(|i| self.data[i] * o))
    }
}
