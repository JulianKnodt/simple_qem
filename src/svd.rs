use super::sym::SymMatrix3;
use super::{F, kmul, quat_to_mat};

#[inline]
fn approx_givens_quat(a11: F, a12: F, a22: F) -> [F; 2] {
    const G: F = 3. + 2. * std::f64::consts::SQRT_2 as F;

    let ch = 2. * (a11 - a22);
    let sh = a12;

    if ch == 0. && sh == 0. {
        return [1., 0.];
    }

    if G * sh * sh < ch * ch {
        let w = 1. / (sh * sh + ch * ch).sqrt();
        [w * ch, w * sh]
    } else {
        const PI: F = std::f64::consts::PI as F;
        [(PI / 8.).cos(), (PI / 8.).sin()]
    }
}

#[inline]
fn jacobi_conjugation(
    x: usize,
    y: usize,
    z: usize,
    s @ [s11, s21, s22, _, _, _]: [F; 6],
    q: [F; 4],
) -> ([F; 6], [F; 4]) {
    let [ch, sh] = approx_givens_quat(s11, s21, s22);

    let scale = ch * ch + sh * sh;
    let a = (ch * ch - sh * sh) / scale;
    let b = (2. * ch * sh) / scale;

    let s = conj_sym(s, a, b);

    let tmp = kmul(sh, [q[0], q[1], q[2]]);
    let sh = sh * q[3];
    let mut q = kmul(ch, q);

    q[z] += sh;
    q[3] -= tmp[z];
    q[x] += tmp[y];
    q[y] -= tmp[x];

    let [s11, s21, s22, s31, s32, s33] = s;
    ([s22, s32, s33, s21, s31, s11], q)
}

fn conj_sym([s11, s21, s22, s31, s32, s33]: [F; 6], a: F, b: F) -> [F; 6] {
    [
        a * (a * s11 + b * s21) + b * (a * s21 + b * s22),
        a * (-b * s11 + a * s21) + b * (-b * s21 + a * s22),
        -b * (-b * s11 + a * s21) + a * (-b * s21 + a * s22),
        a * s31 + b * s32,
        -b * s31 + a * s32,
        s33,
    ]
}

fn jacobi_eigen(mut s: [F; 6]) -> [F; 4] {
    let mut q = [0., 0., 0., 1.].map(F::from);
    for _ in 0..15 {
        let (ns, nq) = jacobi_conjugation(0, 1, 2, s, q);
        let (ns, nq) = jacobi_conjugation(1, 2, 0, ns, nq);
        let (ns, nq) = jacobi_conjugation(2, 0, 1, ns, nq);
        s = ns;
        q = nq;
    }
    q
}

/// Computes the eigenvectors of a symmetric matrix using jacobi iterations.
pub fn eigen_jacobi(s: SymMatrix3) -> [[F; 3]; 3] {
    let [s00, s01, s02, s11, s12, s22] = s.data;
    let alt_order = [s00, s01, s11, s02, s12, s22];
    super::transpose3(quat_to_mat(jacobi_eigen(alt_order)))
}
