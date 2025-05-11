#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(generic_arg_infer)]
#![feature(cmp_minmax)]
#![feature(let_chains)]

pub type F = pars3d::F;

pub mod quadric;
pub mod svd;
pub mod sym;

mod vec;
use vec::*;

mod manifold;

mod parameters;
pub use parameters::Args;

mod qem;
pub use qem::simplify;
