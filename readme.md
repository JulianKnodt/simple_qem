# Simple QEM

Simple geometry-only Quad-Dominant Mesh Reduction

- Simplified rewrite of the original which does not allow for preserving any attributes,
  symmetry, joints, or wedges.
- Implements edge-weighting, and re-ordering of the priority queue.

### Example Usage

```
cargo run --release -- -i dense_cube.obj -o tmp.obj --target-tri-ratio 0.25 --abs-eps 1e-3
```

### Important

This is not the implementation used in the original paper, but a rewritten simplified version. The
original version is owned by Lightspeed Studios. This version has extremely minimal
functionality compared to that version.

