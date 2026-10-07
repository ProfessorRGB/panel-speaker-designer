// One interface over both modal solvers: the Rayleigh-Ritz / closed-form
// solver for rectangles (fast and exact enough to be the reference) and FEA
// for every other shape. Callers get mode frequencies plus evaluators that
// sample mode shapes at grid cells or arbitrary points; points outside the
// panel material evaluate to NaN.

use crate::fem;
use crate::geometry::{Cutout, Outline, Panel, Pt};
use crate::mesh::{self, Mesh};
use crate::plate::{self, Boundary, GridBasis, Plate, PointBasis, SolveKey};
use std::f64::consts::PI;

// FEA resolution: elements per half-wavelength at the top frequency. Gives
// about 1–2% frequency error at the top of the band (see fem tests).
const ELEMENTS_PER_HALF_WAVE: f64 = 5.0;
const MAX_FEA_NODES: f64 = 15_000.0;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Discretisation {
    Analytic(SolveKey),
    Fea { h: f64 },
}

/// Everything the modes depend on, so equal keys give identical results.
#[derive(Clone, Debug, PartialEq)]
pub struct ModelKey {
    pub plate: Plate,
    pub boundary: Boundary,
    pub panel: Panel,
    pub cutouts: Vec<Cutout>,
    pub disc: Discretisation,
}

impl ModelKey {
    pub fn for_freq(plate: Plate, boundary: Boundary, panel: Panel, cutouts: Vec<Cutout>, freq_max: f64) -> ModelKey {
        let disc = if panel == Panel::Rectangle && cutouts.is_empty() {
            Discretisation::Analytic(SolveKey::for_freq(plate, boundary, freq_max))
        } else {
            let d = plate.rigidities();
            let d_min = d.d11.min(d.d22);
            let rho_h = plate.rho * plate.h;
            let omega = 2.0 * PI * freq_max;
            let k = (omega * omega * rho_h / d_min).powf(0.25);
            let mut h = (PI / k / ELEMENTS_PER_HALF_WAVE).min(plate.lx.min(plate.ly) / 12.0);
            // Keep the node count bounded. Refined meshes come out at about
            // 1.4 nodes per h² of area (measured).
            let nodes = 1.4 * plate.lx * plate.ly / (h * h);
            if nodes > MAX_FEA_NODES {
                h *= (nodes / MAX_FEA_NODES).sqrt();
            }
            Discretisation::Fea { h }
        };
        ModelKey { plate, boundary, panel, cutouts, disc }
    }

    fn fea_freq_limit(&self, h: f64) -> f64 {
        let d = self.plate.rigidities();
        let k = PI / (ELEMENTS_PER_HALF_WAVE * h);
        k * k * (d.d11.min(d.d22) / (self.plate.rho * self.plate.h)).sqrt() / (2.0 * PI)
    }
}

pub struct ModeEntry {
    pub freq: f64,
    pub label: Option<(u32, u32)>,  // nodal-line counts, where meaningful
}

enum Data {
    Analytic(plate::Solution),
    Fea { mesh: Mesh, w: Vec<Vec<f64>>, locator: Locator },
}

pub struct Model {
    pub key: ModelKey,
    pub modes: Vec<ModeEntry>,  // ascending, all ≤ freq_limit
    pub freq_limit: f64,
    pub outline: Outline,       // panel shape in metres, for masks and drawing
    data: Data,
}

impl Model {
    pub fn solve(key: ModelKey) -> Result<Model, String> {
        let p = key.plate;
        let shape = key.panel.to_shape(p.lx, p.ly, &key.cutouts);
        shape.validate(&key.cutouts)?;
        match key.disc {
            Discretisation::Analytic(sk) => {
                let sol = plate::solve(sk);
                let modes = sol.modes.iter()
                    .map(|m| ModeEntry { freq: m.freq, label: Some((m.m, m.n)) })
                    .collect();
                // Fine flattening only matters for curved outlines.
                let outline = Outline::new(&shape, p.lx.max(p.ly));
                Ok(Model { key, modes, freq_limit: sol.freq_limit, outline, data: Data::Analytic(sol) })
            }
            Discretisation::Fea { h } => {
                let outline = Outline::new(&shape, h);
                let mesh = mesh::mesh(&outline, h)?;
                let pinned: Vec<usize> = match key.boundary {
                    Boundary::Free => Vec::new(),
                    Boundary::SimplySupported => outline_nodes(&mesh, &outline.rings()[0], h),
                };
                let freq_limit = key.fea_freq_limit(h);
                let fe = fem::solve(&mesh, &p.rigidities(), p.rho * p.h, p.rho * p.h.powi(3) / 12.0, freq_limit, &pinned)?;
                let modes = fe.iter().map(|m| ModeEntry { freq: m.freq, label: None }).collect();
                let w = fe.into_iter().map(|m| m.w).collect();
                let locator = Locator::new(&mesh);
                Ok(Model { key, modes, freq_limit, outline, data: Data::Fea { mesh, w, locator } })
            }
        }
    }

    /// Evaluator for mode shapes at the centres of an n×n grid over the
    /// bounding box (row-major, rows along y).
    pub fn grid_eval(&self, n: usize) -> Evaluator<'_> {
        let points: Vec<Pt> = (0..n * n)
            .map(|i| [((i % n) as f64 + 0.5) / n as f64, ((i / n) as f64 + 0.5) / n as f64])
            .collect();
        match &self.data {
            Data::Analytic(sol) => Evaluator::Grid(sol, GridBasis::new(sol, n)),
            Data::Fea { .. } => self.point_eval(points),
        }
    }

    /// Evaluator for mode shapes at points given in normalised coordinates
    /// (0..1 across the bounding box).
    pub fn point_eval(&self, points: Vec<Pt>) -> Evaluator<'_> {
        match &self.data {
            Data::Analytic(sol) => Evaluator::Points(sol, PointBasis::new(sol, points.iter().map(|p| (p[0], p[1])).collect())),
            Data::Fea { mesh, w, locator } => {
                let (lx, ly) = (self.key.plate.lx, self.key.plate.ly);
                let located = points.iter().map(|p| locator.locate(mesh, [p[0] * lx, p[1] * ly])).collect();
                Evaluator::Mesh(mesh, w, located)
            }
        }
    }
}

/// Mesh nodes lying on a ring of the outline.
fn outline_nodes(mesh: &Mesh, ring: &[Pt], h: f64) -> Vec<usize> {
    let ring_only = Outline::from_rings(vec![ring.to_vec()]);
    (0..mesh.nodes.len())
        .filter(|&i| ring_only.distance_to_edge(mesh.nodes[i]) < 1e-6 * h)
        .collect()
}

pub enum Evaluator<'a> {
    Grid(&'a plate::Solution, GridBasis),
    Points(&'a plate::Solution, PointBasis),
    // Linear interpolation of nodal displacement within the containing triangle.
    Mesh(&'a Mesh, &'a Vec<Vec<f64>>, Vec<Option<(usize, [f64; 3])>>),
}

impl Evaluator<'_> {
    /// Shape of mode `k` at each point; NaN outside the panel.
    pub fn eval(&self, k: usize) -> Vec<f64> {
        match self {
            Evaluator::Grid(sol, basis) => basis.eval(sol, &sol.modes[k]),
            Evaluator::Points(sol, basis) => basis.eval(sol, &sol.modes[k]),
            Evaluator::Mesh(mesh, w, located) => located.iter()
                .map(|loc| match loc {
                    Some((t, l)) => {
                        let tri = mesh.tris[*t];
                        (0..3).map(|q| l[q] * w[k][tri[q]]).sum()
                    }
                    None => f64::NAN,
                })
                .collect(),
        }
    }
}

// ── Point location ───────────────────────────────────────────────────────────

/// Uniform bucket grid of triangles for point-in-mesh queries.
struct Locator {
    origin: Pt,
    cell: f64,
    cols: usize,
    rows: usize,
    buckets: Vec<Vec<usize>>,
}

impl Locator {
    fn new(mesh: &Mesh) -> Locator {
        let (mut min, mut max) = ([f64::MAX; 2], [f64::MIN; 2]);
        for p in &mesh.nodes {
            for a in 0..2 {
                min[a] = min[a].min(p[a]);
                max[a] = max[a].max(p[a]);
            }
        }
        let side = ((mesh.tris.len() as f64).sqrt() / 2.0).max(1.0);
        let cell = ((max[0] - min[0]).max(max[1] - min[1]) / side).max(1e-12);
        let cols = ((max[0] - min[0]) / cell).ceil() as usize + 1;
        let rows = ((max[1] - min[1]) / cell).ceil() as usize + 1;
        let mut buckets = vec![Vec::new(); cols * rows];
        for (t, tri) in mesh.tris.iter().enumerate() {
            let ps = tri.map(|i| mesh.nodes[i]);
            let c0 = ((ps.iter().map(|p| p[0]).fold(f64::MAX, f64::min) - min[0]) / cell) as usize;
            let c1 = ((ps.iter().map(|p| p[0]).fold(f64::MIN, f64::max) - min[0]) / cell) as usize;
            let r0 = ((ps.iter().map(|p| p[1]).fold(f64::MAX, f64::min) - min[1]) / cell) as usize;
            let r1 = ((ps.iter().map(|p| p[1]).fold(f64::MIN, f64::max) - min[1]) / cell) as usize;
            for r in r0..=r1.min(rows - 1) {
                for c in c0..=c1.min(cols - 1) {
                    buckets[r * cols + c].push(t);
                }
            }
        }
        Locator { origin: min, cell, cols, rows, buckets }
    }

    /// Containing triangle and barycentric coordinates, if inside the mesh.
    fn locate(&self, mesh: &Mesh, p: Pt) -> Option<(usize, [f64; 3])> {
        let c = ((p[0] - self.origin[0]) / self.cell).floor();
        let r = ((p[1] - self.origin[1]) / self.cell).floor();
        if c < 0.0 || r < 0.0 || c as usize >= self.cols || r as usize >= self.rows {
            return None;
        }
        let tol = 1e-9;
        for &t in &self.buckets[r as usize * self.cols + c as usize] {
            let [a, b, cc] = mesh.tris[t].map(|i| mesh.nodes[i]);
            let det = (b[0] - a[0]) * (cc[1] - a[1]) - (cc[0] - a[0]) * (b[1] - a[1]);
            let l1 = ((b[0] - p[0]) * (cc[1] - p[1]) - (cc[0] - p[0]) * (b[1] - p[1])) / det;
            let l2 = ((cc[0] - p[0]) * (a[1] - p[1]) - (a[0] - p[0]) * (cc[1] - p[1])) / det;
            let l3 = 1.0 - l1 - l2;
            if l1 >= -tol && l2 >= -tol && l3 >= -tol {
                return Some((t, [l1, l2, l3]));
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn acrylic(lx: f64, ly: f64) -> Plate {
        Plate { lx, ly, h: 0.003, ex: 3.2e9, ey: 3.2e9, g: 3.2e9 / 2.74, nu_xy: 0.37, rho: 1190.0 }
    }

    #[test]
    fn mesh_evaluator_reproduces_nodal_values_and_masks_outside() {
        let key = ModelKey::for_freq(acrylic(0.3, 0.3), Boundary::Free, Panel::Ellipse, vec![], 800.0);
        let model = Model::solve(key).unwrap();
        assert!(model.modes.len() >= 3);
        // Corners of the bounding box are outside a circle; the centre is inside.
        let values = model.point_eval(vec![[0.01, 0.01], [0.5, 0.5]]).eval(0);
        assert!(values[0].is_nan());
        assert!(values[1].is_finite());
        // Grid RMS over the inside cells ≈ 0.5 (shapes are normalised).
        let grid = model.grid_eval(80).eval(0);
        let inside: Vec<f64> = grid.into_iter().filter(|v| v.is_finite()).collect();
        let rms = (inside.iter().map(|v| v * v).sum::<f64>() / inside.len() as f64).sqrt();
        assert!((rms - 0.5).abs() < 0.03, "rms {rms}");
    }

    #[test]
    fn rounded_rectangle_tends_to_rectangle() {
        let p = acrylic(0.3, 0.2);
        let rect = Model::solve(ModelKey::for_freq(p, Boundary::Free, Panel::Rectangle, vec![], 1000.0)).unwrap();
        let nearly = Model::solve(ModelKey::for_freq(p, Boundary::Free, Panel::RoundedRectangle { radius: 0.002 }, vec![], 1000.0)).unwrap();
        for (a, b) in rect.modes.iter().zip(&nearly.modes).take(10) {
            assert!((a.freq - b.freq).abs() / a.freq < 0.02, "{} vs {}", a.freq, b.freq);
        }
    }

    #[test]
    fn tiny_hole_barely_changes_modes() {
        let p = acrylic(0.3, 0.2);
        let plain = Model::solve(ModelKey::for_freq(p, Boundary::Free, Panel::Rectangle, vec![], 1000.0)).unwrap();
        let holed = Model::solve(ModelKey::for_freq(
            p, Boundary::Free, Panel::Rectangle, vec![Cutout::Hole { x: 0.11, y: 0.07, d: 0.004 }], 1000.0,
        )).unwrap();
        assert!(matches!(holed.key.disc, Discretisation::Fea { .. }));
        for (a, b) in plain.modes.iter().zip(&holed.modes).take(12) {
            assert!((a.freq - b.freq).abs() / a.freq < 0.02, "{} vs {}", a.freq, b.freq);
        }
    }

    #[test]
    fn slotted_panel_converges_with_mesh_refinement() {
        let p = acrylic(0.3, 0.2);
        let slot = vec![Cutout::Slot { x: 0.15, y: 0.1, length: 0.12, width: 0.01, angle: 90.0 }];
        let coarse = ModelKey::for_freq(p, Boundary::Free, Panel::Rectangle, slot.clone(), 1500.0);
        let Discretisation::Fea { h } = coarse.disc else { panic!("expected FEA") };
        let fine = ModelKey { disc: Discretisation::Fea { h: h / 2.0 }, ..coarse.clone() };
        let (a, b) = (Model::solve(coarse).unwrap(), Model::solve(fine).unwrap());
        let plain = Model::solve(ModelKey::for_freq(p, Boundary::Free, Panel::Rectangle, vec![], 1500.0)).unwrap();
        for (x, y) in a.modes.iter().zip(&b.modes).take(15) {
            assert!((x.freq - y.freq).abs() / y.freq < 0.02, "{} vs {}", x.freq, y.freq);
        }
        // A slot across the middle cuts the panel's bending path along x,
        // so the first (2,0)-type bending modes must drop.
        assert!(b.modes[1].freq < plain.modes[1].freq * 0.98, "{} vs {}", b.modes[1].freq, plain.modes[1].freq);
    }
}
