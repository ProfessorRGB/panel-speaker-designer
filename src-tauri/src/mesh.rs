// Triangle meshing of a panel shape: constrained Delaunay triangulation of the
// flattened outline and holes, refined to a target element size with a
// minimum-angle guarantee.

use crate::geometry::{Outline, Pt};
use spade::{AngleLimit, ConstrainedDelaunayTriangulation, Point2, RefinementParameters, Triangulation};
use std::collections::{HashMap, HashSet};

pub struct Mesh {
    pub nodes: Vec<Pt>,
    pub tris: Vec<[usize; 3]>,  // counter-clockwise in (x, y)
}

impl Mesh {
    pub fn triangle_area(&self, t: usize) -> f64 {
        let [a, b, c] = self.tris[t].map(|i| self.nodes[i]);
        ((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])) / 2.0
    }
}

// Hard cap so a tiny element size can't run away with memory or time.
const MAX_NODES: usize = 40_000;

/// Meshes `outline` with elements of roughly size `h` (edge length).
pub fn mesh(outline: &Outline, h: f64) -> Result<Mesh, String> {
    let mut vertices = Vec::new();
    let mut edges = Vec::new();
    for ring in outline.rings() {
        let base = vertices.len();
        for p in ring {
            vertices.push(Point2::new(p[0], p[1]));
        }
        for i in 0..ring.len() {
            edges.push([base + i, base + (i + 1) % ring.len()]);
        }
    }

    let mut cdt = ConstrainedDelaunayTriangulation::<Point2<f64>>::bulk_load_cdt(vertices, edges)
        .map_err(|e| format!("Couldn't triangulate the outline: {e:?}"))?;

    // Equilateral triangle of side h.
    let max_area = 3f64.sqrt() / 4.0 * h * h;
    let result = cdt.refine(
        RefinementParameters::<f64>::new()
            .exclude_outer_faces(true)
            .with_max_allowed_area(max_area)
            .with_angle_limit(AngleLimit::from_deg(25.0))
            .with_max_additional_vertices(MAX_NODES),
    );
    if !result.refinement_complete {
        return Err("Mesh too fine: lower the max frequency or simplify the shape".into());
    }
    let excluded: HashSet<_> = result.excluded_faces.into_iter().collect();

    // Keep inner faces only and renumber the vertices they use.
    let mut index = HashMap::new();
    let mut nodes = Vec::new();
    let mut tris = Vec::new();
    for face in cdt.inner_faces() {
        if excluded.contains(&face.fix()) {
            continue;
        }
        let vs = face.vertices();
        let mut tri = [0usize; 3];
        for (k, v) in vs.iter().enumerate() {
            let id = v.fix().index();
            tri[k] = *index.entry(id).or_insert_with(|| {
                let p = v.position();
                nodes.push([p.x, p.y]);
                nodes.len() - 1
            });
        }
        tris.push(tri);
    }
    let mut mesh = Mesh { nodes, tris };
    for t in 0..mesh.tris.len() {
        if mesh.triangle_area(t) < 0.0 {
            mesh.tris[t].swap(1, 2);
        }
    }
    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{Path, Shape};

    #[test]
    fn mesh_covers_shape_with_holes() {
        let shape = Shape {
            outline: Path::rectangle(0.0, 0.0, 0.3, 0.2),
            holes: vec![Path::ellipse(0.1, 0.05, 0.1, 0.1)],
        };
        let h = 0.01;
        let outline = Outline::new(&shape, h);
        let m = mesh(&outline, h).unwrap();
        let area: f64 = (0..m.tris.len()).map(|t| m.triangle_area(t)).sum();
        assert!((area - outline.area()).abs() / outline.area() < 1e-9, "{area} vs {}", outline.area());
        assert!((0..m.tris.len()).all(|t| m.triangle_area(t) > 0.0));
        // No triangle centroid inside the hole.
        for t in &m.tris {
            let c = t.iter().fold([0.0, 0.0], |a, &i| [a[0] + m.nodes[i][0] / 3.0, a[1] + m.nodes[i][1] / 3.0]);
            assert!(outline.contains(c));
        }
    }
}
