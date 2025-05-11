use super::sym::SymMatrix3;
use super::{F, add, dot, kmul};

use core::ops::{Add, AddAssign, Mul, MulAssign};

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct QuadricAccumulator {
    a: SymMatrix3,
    b: [F; 3],
    area: F,

    // volume constraints
    nv: [F; 3],
    dv: F,
}
const EPS: F = 1e-6;

impl AddAssign<Quadric> for QuadricAccumulator {
    fn add_assign(&mut self, q: Quadric) {
        self.a = self.a + q.a;
        self.b = add(self.b, q.b);
        self.area += q.area;

        self.nv = add(self.nv, q.nv);
        self.dv += q.dv;
    }
}

impl QuadricAccumulator {
    pub fn point(&self) -> [F; 3] {
        if self.area == 0. {
            return [0.; 3];
        }
        let inv_a = self.area.recip();
        debug_assert!(inv_a.is_finite(), "{}", self.area);

        let ([e0, e1, e2], [v0, v1, v2]) = self.a.eigen();
        [(e0, v0), (e1, v1), (e2, v2)]
            .into_iter()
            .map(|(e, v)| {
                if e.abs() < 1e-8 {
                    return [0.; 3];
                }
                kmul(-dot(self.b, v) / e, v)
            })
            .fold([0.; 3], add)
    }
    pub fn point_with_volume_opt(&self) -> Option<[F; 3]> {
        if self.area < F::EPSILON {
            return None;
        }
        let inv_a = self.area.recip();
        debug_assert!(inv_a.is_finite(), "{}", self.area);
        invert_quadric_volume(self.a, self.b, self.nv, self.dv)
    }
}

fn invert_quadric_volume(
    a: SymMatrix3,
    [r0, r1, r2]: [F; 3],
    [vx, vy, vz]: [F; 3],
    dv: F,
) -> Option<[F; 3]> {
    let [mxx, mxy, mxz, myy, myz, mzz] = a.data;

    let det2_01_01 = mxx * myy - mxy * mxy;
    let det2_01_02 = mxx * myz - mxz * mxy;
    let det2_01_12 = mxy * myz - mxz * myy;
    let det2_01_03 = mxx * vy - vx * mxy;
    let det2_01_13 = mxy * vy - vx * myy;
    let det2_01_23 = mxz * vy - vx * myz;

    // 3x3 sub-determinants required to calculate 4x4 determinant
    let invx = mzz * det2_01_13 - myz * det2_01_23 - vz * det2_01_12;
    let invy = mxz * det2_01_23 - mzz * det2_01_03 + vz * det2_01_02;
    let invz = myz * det2_01_03 - mxz * det2_01_13 - vz * det2_01_01;

    let det = invx * vx + invy * vy + invz * vz;

    if det < EPS {
        return None;
    }

    let denom = det.recip();

    // remaining 2x2 sub-determinants
    let det2_03_02 = mxx * vz - mxz * vx;
    let det2_03_12 = mxy * vz - mxz * vy;
    let det2_13_12 = myy * vz - myz * vy;

    let det2_03_03 = -vx * vx;
    let det2_03_13 = -vx * vy;
    let det2_03_23 = -vx * vz;

    let det2_13_13 = -vy * vy;
    let det2_13_23 = -vy * vz;

    // remaining 3x3 sub-determinants
    let imxx = mzz * det2_13_13 - myz * det2_13_23 - vz * det2_13_12;
    let imxy = myz * det2_03_23 - mzz * det2_03_13 + vz * det2_03_12;
    let imyy = mzz * det2_03_03 - mxz * det2_03_23 - vz * det2_03_02;

    let imxz = vy * det2_01_23 - vz * det2_01_13;
    let imyz = vz * det2_01_03 - vx * det2_01_23;
    let imzz = vx * det2_01_13 - vy * det2_01_03;

    let numer = [
        r0 * imxx + r1 * imxy + r2 * imxz - invx * dv,
        r0 * imxy + r1 * imyy + r2 * imyz - invy * dv,
        r0 * imxz + r1 * imyz + r2 * imzz - invz * dv,
    ];
    Some(kmul(denom, numer))
}

#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct Quadric {
    pub a: SymMatrix3,
    pub b: [F; 3],
    c: F,

    pub area: F,

    nv: [F; 3],
    dv: F,
}

impl Quadric {
    pub fn cost(&self, p: [F; 3]) -> F {
        let quadratic = dot(p, self.a.vec_mul(p));
        quadratic + 2. * dot(self.b, p) + self.c
    }
    pub fn new_plane(v: [F; 3], n: [F; 3], area: F) -> Self {
        let a = SymMatrix3::outer(n);
        let dist = -dot(n, v);
        let b = kmul(dist, n);
        let c = dist * dist;

        Self {
            a,
            b,
            c,

            area: 1.,

            nv: kmul(area / 3., n),
            dv: area * c / 3.,
        }
    }
}

impl Add for Quadric {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self {
            a: self.a + o.a,
            b: add(self.b, o.b),
            c: self.c + o.c,

            area: self.area + o.area,

            nv: add(self.nv, o.nv),
            dv: self.dv + o.dv,
        }
    }
}

impl Mul<F> for Quadric {
    type Output = Self;
    fn mul(self, o: F) -> Self {
        Self {
            a: self.a * o,
            b: kmul(o, self.b),
            c: o * self.c,

            area: self.area * o,

            nv: self.nv,
            dv: o * self.dv,
        }
    }
}

impl AddAssign for Quadric {
    fn add_assign(&mut self, o: Self) {
        *self = *self + o;
    }
}

impl MulAssign<F> for Quadric {
    fn mul_assign(&mut self, o: F) {
        *self = *self * o;
    }
}
