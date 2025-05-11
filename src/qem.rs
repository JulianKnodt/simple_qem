use std::collections::{HashMap, hash_map::Entry};

use ordered_float::NotNan;
use pars3d::FaceKind;
use priority_queue::PriorityQueue;

use super::{
    F, add, cross, dot, kmul, length,
    manifold::{CollapsibleManifold, EdgeKind},
    normalize,
    quadric::{Quadric, QuadricAccumulator},
    sub,
};

use super::parameters::Args;

/// In-place simplification of planar faces of a mesh.
/// Returns how many faces were removed
pub fn simplify(v: &mut [[F; 3]], f: &mut [FaceKind], args: &Args) {
    let num_tris = f.iter().map(FaceKind::num_tris).sum::<usize>();
    let target_num_tris = (num_tris as F * args.tri_ratio.unwrap_or(0.)).floor() as usize;
    let target_num_tris = target_num_tris.max(args.tri_number.unwrap_or(0));
    if num_tris <= target_num_tris {
        return;
    }

    // normalize all vertices to [-1., 1]
    use std::array::from_fn;
    // Normalize the geometry of this mesh to lay in the unit box.
    let [l, h] = v
        .iter()
        .copied()
        .fold([[F::INFINITY; 3], [F::NEG_INFINITY; 3]], |[l, h], n| {
            [from_fn(|i| l[i].min(n[i])), from_fn(|i| h[i].max(n[i]))]
        });
    let center = kmul(0.5, add(l, h));
    for v in v.iter_mut() {
        *v = sub(*v, center);
    }
    let largest_val = v
        .iter()
        .copied()
        .fold(0. as F, |m, [v0, v1, v2]| m.max(v0).max(v1).max(v2));
    let pos_scale = if largest_val == 0. {
        1.
    } else {
        largest_val.recip()
    };
    for v in v.iter_mut() {
        *v = kmul(pos_scale, *v);
    }

    let mut m = CollapsibleManifold::new_with(v.len(), |vi| (Quadric::default(), v[vi]));

    let mut num_edges = 0;
    let mut avg_edge_len = 0.;

    let mut edge_face_adj: HashMap<[usize; 2], EdgeKind> = HashMap::new();

    // faces for each vertex
    let mut face_verts = vec![vec![]; v.len()];
    let mut face_normals = vec![[0.; 3]; f.len()];
    for (fi, f) in f.iter().enumerate() {
        face_normals[fi] = normalize(f.normal(&v));

        for e in f.edges_ord() {
            edge_face_adj
                .entry(e)
                .and_modify(|p| {
                    p.insert(fi);
                })
                .or_insert_with(|| EdgeKind::Boundary(fi));
            num_edges += 1;
            let [e0, e1] = e.map(|vi| v[vi]);
            avg_edge_len += length(sub(e1, e0));
        }

        for &vi in f.as_slice() {
            face_verts[vi].push(fi);
        }
    }
    for fv in &mut face_verts {
        fv.sort_unstable();
        fv.dedup();
    }

    avg_edge_len /= num_edges as F;

    for f in f.iter() {
        for [e0, e1] in f.edges_ord() {
            m.add_edge(e0, e1);
        }
    }

    for (fi, f) in f.iter().enumerate() {
        if f.is_empty() {
            continue;
        }
        let area = f.area(&v);
        assert!(area >= 0.);
        let area = area + 1e-2; // < important for stability
        let n = normalize(face_normals[fi]);

        if length(n) == 0. {
            // Handle this better (there will be many degenerate triangles)
            continue;
        }

        let f_slice = f.as_slice();
        for (i, &vi) in f_slice.iter().enumerate() {
            let curr = v[vi];
            let pi = f_slice[i.checked_sub(1).unwrap_or_else(|| f.len() - 1)];
            let prev = v[pi];
            let ni = f_slice[(i + 1) % f.len()];
            let e = std::cmp::minmax(vi, ni);
            let next = v[ni];

            let interior_angle = {
                let e0 = normalize(sub(prev, curr));
                let e1 = normalize(sub(next, curr));
                dot(e0, e1).clamp(-1., 1.).acos()
            };
            let mut q = Quadric::new_plane(curr, n, area) * interior_angle;
            q.area = area;
            m.data[vi].0 += q;

            const PI: F = std::f64::consts::PI as F;

            macro_rules! dihedral_angle {
                ($f0: expr, $f1: expr) => {{
                    let fn0 = normalize(face_normals[$f0]);
                    let fn1 = normalize(face_normals[$f1]);
                    let angle = dot(fn0, fn1);
                    assert!((-1.0001..=1.0001).contains(&angle), "{angle}");
                    let angle = angle.clamp(-1., 1.);

                    let v = angle.acos();
                    assert!((0.0..=PI).contains(&v), "{v} {angle}");
                    v
                }};
            }

            let e_w = match edge_face_adj[&e] {
                EdgeKind::Boundary(_) => 4.,
                EdgeKind::Manifold([a, b]) => dihedral_angle!(a, b) / PI,
                EdgeKind::NonManifold(_) => 4.,
            };
            let e_w = e_w.max(args.min_edge_weight);

            let edge_dir = sub(curr, next);
            let edge_len = length(edge_dir);
            let edge_len = edge_len / avg_edge_len;
            if edge_len == 0. {
                continue;
            }
            let edge_dir = normalize(edge_dir);
            let edge_quadric = Quadric::new_plane(curr, normalize(cross(n, edge_dir)), 0.);

            let total_e_w = e_w * edge_len;
            let mut edge_quadric = edge_quadric * total_e_w.max(1e-4);
            edge_quadric.area = 0.;

            m.data[vi].0 += edge_quadric;
            m.data[ni].0 += edge_quadric;
        }
    }

    let mut curr_costs = vec![0.; v.len()];

    macro_rules! update_cost_of_edge {
        ($e0:expr, $e1: expr) => {{
            let [e0, e1] = std::cmp::minmax($e1, $e0);
            let mut q_acc = QuadricAccumulator::default();
            let &(q0, _) = m.get(e0);
            q_acc += q0;
            let &(q1, _) = m.get(e1);
            q_acc += q1;
            let p = q_acc
                .point_with_volume_opt()
                .unwrap_or_else(|| q_acc.point());
            debug_assert!(p.iter().copied().all(F::is_finite));

            let total_cost = (q0 + q1).cost(p).max(0.) - curr_costs[e0] - curr_costs[e1];

            NotNan::new(-total_cost).unwrap()
        }};
    }

    let mut recencies = HashMap::new();
    macro_rules! update_edge_face_adj {
        ($e0: expr, $e_dst: expr, $e0_adj:expr) => {{
            let e0 = $e0;
            let e_dst = $e_dst;
            debug_assert!(e0 <= e_dst);
            let adj = $e0_adj;

            debug_assert!(!m.is_deleted(adj));
            let e = std::cmp::minmax(adj, e0);

            let Some(prev_ef) = edge_face_adj.remove(&e) else {
                continue;
            };
            let new_e = std::cmp::minmax(adj, e_dst);

            match edge_face_adj.entry(new_e) {
                Entry::Vacant(v) => {
                    v.insert(prev_ef);
                }
                Entry::Occupied(mut o) => {
                    // only keep faces which have more than 2 vertices
                    let faces_culled = prev_ef
                        .as_slice()
                        .into_iter()
                        .chain(o.get().as_slice())
                        .filter_map(|&fi| {
                            f[fi].remap(|v| match v == e0 {
                                true => e_dst,
                                false => m.get_new_vertex(v),
                            });
                            (!f[fi].is_degenerate()).then_some(fi)
                        });
                    match EdgeKind::from_iter(faces_culled) {
                        Some(mut ek) => {
                            ek.dedup_by_key(|fi| f[fi].as_slice());
                            o.insert(ek)
                        }
                        None => o.remove(),
                    };
                }
            }

            if let Some(prev) = recencies.remove(&e) {
                recencies
                    .entry(new_e)
                    .and_modify(|p| *p += prev)
                    .or_insert(prev);
            }
        }};
    }

    let mut pq = PriorityQueue::new();
    for [e0, e1] in m.ord_edges() {
        pq.push([e0, e1], update_cost_of_edge!(e0, e1));
    }

    let mut curr_tris = f.iter().map(|f| f.num_tris()).sum::<usize>();

    let mut snd_pq = PriorityQueue::new();
    let mut did_update = vec![];
    'outer: while let Some((e, q)) = pq.pop() {
        assert!(snd_pq.is_empty());
        snd_pq.push(e, (0, q));
        recencies.clear();
        while let Some(([e0, e1], (rec, q_err))) = snd_pq.pop() {
            debug_assert!(e0 < e1);
            if m.is_deleted(e0) || m.is_deleted(e1) {
                continue;
            }
            if curr_tris <= target_num_tris {
                break 'outer;
            }

            let mut q_acc = QuadricAccumulator::default();
            let q0 = m.get(e0).0;
            let q1 = m.get(e1).0;
            q_acc += q0;
            q_acc += q1;
            let pos = q_acc
                .point_with_volume_opt()
                .unwrap_or_else(|| q_acc.point());
            let q01 = q0 + q1;

            // -- Commit

            if let Some(adj_faces) = edge_face_adj.get(&[e0, e1]) {
                for &af in adj_faces.as_slice() {
                    let Some([q0, q1]) = f[af].quad_opp_edge(e0, e1) else {
                        continue;
                    };
                    let r = recencies.entry(std::cmp::minmax(q0, q1)).or_insert(rec);
                    *r += 1;
                }
            };

            for adj in m.vertex_adj(e0) {
                update_edge_face_adj!(e0, e1, adj);
            }

            m.merge(e0, e1, |_, _| {
                curr_costs[e1] = q01.cost(pos).max(0.);
                (q01, pos)
            });
            debug_assert!(m.is_deleted(e0));
            debug_assert!(!m.is_deleted(e1));

            let [ef0, ef1] = face_verts.get_disjoint_mut([e0, e1]).unwrap();
            let ef1_len = ef1.len();
            for f in std::mem::take(ef0) {
                if !ef1[0..ef1_len].contains(&f) {
                    ef1.push(f);
                }
            }
            let prev_tri = ef1.iter().map(|&fi| f[fi].num_tris()).sum::<usize>();
            ef1.retain(|&fi| {
                f[fi].remap(|vi| m.get_new_vertex(vi));
                // important to not rotate here otherwise the ordering may change
                // leading to non-manifold edge introduction
                let retain = !f[fi].canonicalize_no_rotate();
                if !retain {
                    // necessary for triangle counting
                    f[fi] = FaceKind::empty();
                }
                retain
            });
            let new_tri = ef1.iter().map(|&fi| f[fi].num_tris()).sum::<usize>();
            curr_tris -= prev_tri - new_tri;

            did_update.clear();
            let e_dst = m.get_new_vertex(e1);
            for adj in m.vertex_adj(e_dst) {
                let prio = update_cost_of_edge!(e_dst, adj);
                let adj_e = std::cmp::minmax(e_dst, adj);
                snd_pq.remove(&adj_e);
                pq.push(adj_e, prio);
                did_update.push(adj_e);
            }

            for adj in m.vertex_adj(e_dst) {
                for adj2 in m.vertex_adj(adj) {
                    let adj_e = std::cmp::minmax(adj, adj2);
                    if adj2 == e_dst || did_update.contains(&adj_e) {
                        continue;
                    }

                    did_update.push(adj_e);
                    let prio = update_cost_of_edge!(adj, adj2);
                    let recency = recencies.get(&adj_e).copied().unwrap_or(0);

                    if !approx_eq(*prio, *q_err, args.abs_eps) {
                        snd_pq.remove(&adj_e);
                        pq.push(adj_e, prio);
                        continue;
                    }
                    let changed = snd_pq.change_priority(&adj_e, (recency, prio)).is_some();
                    if !changed {
                        pq.push(adj_e, prio);
                    }
                }
            }

            while let Some((_e, nq_err)) = pq.peek()
                && approx_eq(**nq_err, *q_err, args.abs_eps)
            {
                let (e, nq_err) = pq.pop().unwrap();
                let recency = recencies.get(&e).copied().unwrap_or(0);
                snd_pq.push(e, (recency, nq_err));
            }
        }
    }

    for (vi, &(_, p)) in m.vertices() {
        v[vi] = p;
    }

    // denormalize all output vertices
    let inv_pos_scale = pos_scale.recip();
    for v in v.iter_mut() {
        *v = add(kmul(inv_pos_scale, *v), center);
    }
}

fn approx_eq(a: F, b: F, abs_eps: F) -> bool {
    if a == b {
        return true;
    }
    (a - b).abs() < abs_eps
}
