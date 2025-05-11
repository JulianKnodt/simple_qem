pub use super::F;
use clap::Parser;

/// Mesh geometry decimation.
#[derive(Parser, Default, Debug)]
#[clap(group(
            clap::ArgGroup::new("target")
                .required(true)
                .args(&["tri_ratio", "tri_number"]),
        ))]
pub struct Args {
    /// Input mesh file.
    #[arg(short, long, required = true)]
    pub input: String,

    /// Output mesh file.
    #[arg(short, long, required = true)]
    pub output: String,

    /// Approximate ratio of output/input tris.
    #[arg(short = 't', long, group = "target")]
    pub tri_ratio: Option<F>,

    /// Approximate number of output tris.
    #[arg(short = 'T', long, group = "target")]
    pub tri_number: Option<usize>,

    /// Scalar weight of quadric defined along edges
    #[arg(long, default_value_t = 1.)]
    pub edge_weight: F,

    /// The minimum weight per edge.
    #[arg(long, default_value_t = 0.01)]
    pub min_edge_weight: F,

    /// Absolute epsilon to use when considering two edges to be nearly equivalent
    #[arg(long, default_value_t = 5e-6)]
    pub abs_eps: F,
}
