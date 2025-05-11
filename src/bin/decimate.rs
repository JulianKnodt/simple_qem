use clap::Parser;

use quad_dom_qem::{Args, simplify};

pub fn main() {
    let args = Args::parse();
    let mut scene = pars3d::load(&args.input).expect(&format!("Failed to open {}", args.input));
    let mut m = scene.into_flattened_mesh();
    // just the vertices and faces, fuse together identical positions
    m.geometry_only();
    simplify(&mut m.v, &mut m.f, &args);
    m.delete_empty_faces();
    m.delete_unused_vertices();
    m.repopulate_scene(&mut scene);
    pars3d::save(&args.output, &scene).expect(&format!("Failed to save to {}", args.output));
}
