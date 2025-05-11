use super::F;

use std::array::from_fn;

#[inline]
pub fn cross([x, y, z]: [F; 3], [a, b, c]: [F; 3]) -> [F; 3] {
    [y * c - z * b, z * a - x * c, x * b - y * a]
}

#[inline]
pub fn add<const N: usize>(a: [F; N], b: [F; N]) -> [F; N] {
    from_fn(|i| a[i] + b[i])
}

#[inline]
pub fn sub<const N: usize>(a: [F; N], b: [F; N]) -> [F; N] {
    from_fn(|i| a[i] - b[i])
}

/// L-2 norm of a vector
pub fn length<const N: usize>(v: [F; N]) -> F {
    dot(v, v).sqrt()
}

#[inline]
pub fn normalize<const N: usize>(v: [F; N]) -> [F; N] {
    let sum: F = v.iter().map(|v| v * v).sum();
    if sum < 1e-20 {
        return [0.; N];
    }
    let s = sum.sqrt().recip();
    v.map(|v| v * s)
}

#[inline]
pub fn dot<const N: usize>(a: [F; N], b: [F; N]) -> F {
    (0..N).map(|i| a[i] * b[i]).sum()
}

#[inline]
pub fn kmul<const N: usize>(k: F, xyz: [F; N]) -> [F; N] {
    xyz.map(|v| v * k)
}

/// returns each row of the matrix representing a quaternion
pub fn quat_to_mat([x, y, z, w]: [F; 4]) -> [[F; 3]; 3] {
    let qxx = x * x;
    let qyy = y * y;
    let qzz = z * z;
    let qxz = x * z;
    let qxy = x * y;
    let qyz = y * z;
    let qwx = w * x;
    let qwy = w * y;
    let qwz = w * z;

    [
        [1. - 2. * (qyy + qzz), 2. * (qxy - qwz), 2. * (qxz + qwy)],
        [2. * (qxy + qwz), 1. - 2. * (qxx + qzz), 2. * (qyz - qwx)],
        [2. * (qxz - qwy), 2. * (qyz + qwx), 1. - 2. * (qxx + qyy)],
    ]
}

pub fn transpose3(s: [[F; 3]; 3]) -> [[F; 3]; 3] {
    let [[a, b, c], [d, e, f], [g, h, i]] = s;
    [[a, d, g], [b, e, h], [c, f, i]]
}
