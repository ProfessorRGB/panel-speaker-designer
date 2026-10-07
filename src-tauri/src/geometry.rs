// Panel outlines as closed paths of straight and cubic Bézier segments.
//
// Built-in shapes are generated as paths, so custom Bézier outlines (and
// holes) go through exactly the same meshing and analysis. Coordinates are in
// metres, with the origin at the top-left of the bounding box and y pointing
// down, matching the heat-map grid.

use serde::Deserialize;
use std::f64::consts::PI;

pub type Pt = [f64; 2];

#[derive(Clone, Debug, PartialEq, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Segment {
    Line { to: Pt },
    Cubic { c1: Pt, c2: Pt, to: Pt },
}

/// A closed path: it starts at `start` and the last segment returns to it
/// (if it doesn't end there, a straight closing edge is implied).
#[derive(Clone, Debug, PartialEq, Deserialize)]
pub struct Path {
    pub start: Pt,
    pub segments: Vec<Segment>,
}

/// A panel: one outline plus any holes.
#[derive(Clone, Debug, PartialEq)]
pub struct Shape {
    pub outline: Path,
    pub holes: Vec<Path>,
}

/// The built-in panel shapes, each filling a width × height bounding box.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Panel {
    Rectangle,
    RoundedRectangle { radius: f64 },
    Ellipse,
    Polygon { sides: usize },
}

impl Panel {
    pub fn to_shape(self, w: f64, h: f64) -> Shape {
        let outline = match self {
            Panel::Rectangle => Path::rectangle(0.0, 0.0, w, h),
            Panel::RoundedRectangle { radius } => Path::rounded_rectangle(0.0, 0.0, w, h, radius),
            Panel::Ellipse => Path::ellipse(0.0, 0.0, w, h),
            Panel::Polygon { sides } => Path::regular_polygon(0.0, 0.0, w, h, sides),
        };
        Shape { outline, holes: Vec::new() }
    }
}

// Bézier control distance for a quarter circle of unit radius.
const KAPPA: f64 = 0.552_284_749_830_793_4;

impl Path {
    pub fn rectangle(x: f64, y: f64, w: f64, h: f64) -> Path {
        Path {
            start: [x, y],
            segments: vec![
                Segment::Line { to: [x + w, y] },
                Segment::Line { to: [x + w, y + h] },
                Segment::Line { to: [x, y + h] },
                Segment::Line { to: [x, y] },
            ],
        }
    }

    /// Rectangle with circular-arc corners of radius `r` (clamped to fit).
    pub fn rounded_rectangle(x: f64, y: f64, w: f64, h: f64, r: f64) -> Path {
        let r = r.clamp(0.0, w.min(h) / 2.0);
        if r <= 0.0 {
            return Path::rectangle(x, y, w, h);
        }
        let k = r * KAPPA;
        let (x1, y1) = (x + w, y + h);
        Path {
            start: [x + r, y],
            segments: vec![
                Segment::Line { to: [x1 - r, y] },
                Segment::Cubic { c1: [x1 - r + k, y], c2: [x1, y + r - k], to: [x1, y + r] },
                Segment::Line { to: [x1, y1 - r] },
                Segment::Cubic { c1: [x1, y1 - r + k], c2: [x1 - r + k, y1], to: [x1 - r, y1] },
                Segment::Line { to: [x + r, y1] },
                Segment::Cubic { c1: [x + r - k, y1], c2: [x, y1 - r + k], to: [x, y1 - r] },
                Segment::Line { to: [x, y + r] },
                Segment::Cubic { c1: [x, y + r - k], c2: [x + r - k, y], to: [x + r, y] },
            ],
        }
    }

    /// Ellipse filling the box (a circle when w = h).
    pub fn ellipse(x: f64, y: f64, w: f64, h: f64) -> Path {
        let (rx, ry) = (w / 2.0, h / 2.0);
        let (cx, cy) = (x + rx, y + ry);
        let (kx, ky) = (rx * KAPPA, ry * KAPPA);
        Path {
            start: [cx + rx, cy],
            segments: vec![
                Segment::Cubic { c1: [cx + rx, cy + ky], c2: [cx + kx, cy + ry], to: [cx, cy + ry] },
                Segment::Cubic { c1: [cx - kx, cy + ry], c2: [cx - rx, cy + ky], to: [cx - rx, cy] },
                Segment::Cubic { c1: [cx - rx, cy - ky], c2: [cx - kx, cy - ry], to: [cx, cy - ry] },
                Segment::Cubic { c1: [cx + kx, cy - ry], c2: [cx + rx, cy - ky], to: [cx + rx, cy] },
            ],
        }
    }

    /// Regular polygon with `sides` vertices, stretched to fill the box, with
    /// a flat bottom edge.
    pub fn regular_polygon(x: f64, y: f64, w: f64, h: f64, sides: usize) -> Path {
        let sides = sides.max(3);
        let angle = |i: usize| PI / 2.0 + PI / sides as f64 + 2.0 * PI * i as f64 / sides as f64;
        let raw: Vec<Pt> = (0..sides).map(|i| [angle(i).cos(), angle(i).sin()]).collect();
        // Scale the unit polygon's bounding box onto the requested box.
        let (min_x, max_x) = raw.iter().fold((f64::MAX, f64::MIN), |(a, b), p| (a.min(p[0]), b.max(p[0])));
        let (min_y, max_y) = raw.iter().fold((f64::MAX, f64::MIN), |(a, b), p| (a.min(p[1]), b.max(p[1])));
        let pts: Vec<Pt> = raw.iter()
            .map(|p| [x + (p[0] - min_x) / (max_x - min_x) * w, y + (p[1] - min_y) / (max_y - min_y) * h])
            .collect();
        Path {
            start: pts[0],
            segments: pts[1..].iter().chain([&pts[0]]).map(|&to| Segment::Line { to }).collect(),
        }
    }

    /// Polygon approximation with no edge longer than `max_len` and curves
    /// followed to within `max_len / 50`. The closing point is not repeated.
    pub fn flatten(&self, max_len: f64) -> Vec<Pt> {
        let mut pts = vec![self.start];
        let mut cur = self.start;
        let tol = max_len / 50.0;
        for seg in &self.segments {
            match *seg {
                Segment::Line { to } => {
                    push_line(&mut pts, cur, to, max_len);
                    cur = to;
                }
                Segment::Cubic { c1, c2, to } => {
                    // Enough pieces for both the length limit and the sagitta
                    // tolerance (chord error ≈ L²·curvature/8, bounded via the
                    // control polygon's second differences).
                    let poly_len = dist(cur, c1) + dist(c1, c2) + dist(c2, to);
                    let dd = norm(sub(add(cur, c2), scale(c1, 2.0))).max(norm(sub(add(c1, to), scale(c2, 2.0))));
                    let by_len = (poly_len / max_len).ceil();
                    let by_tol = (3.0 * dd / (4.0 * tol)).sqrt().ceil();
                    let n = by_len.max(by_tol).max(1.0) as usize;
                    for i in 1..=n {
                        pts.push(cubic_at(cur, c1, c2, to, i as f64 / n as f64));
                    }
                    cur = to;
                }
            }
        }
        if dist(cur, self.start) > 1e-12 {
            push_line(&mut pts, cur, self.start, max_len);
        } else {
            pts.pop();
        }
        pts.dedup_by(|a, b| dist(*a, *b) < 1e-12);
        pts
    }
}

fn push_line(pts: &mut Vec<Pt>, from: Pt, to: Pt, max_len: f64) {
    let n = (dist(from, to) / max_len).ceil().max(1.0) as usize;
    for i in 1..=n {
        let t = i as f64 / n as f64;
        pts.push([from[0] + t * (to[0] - from[0]), from[1] + t * (to[1] - from[1])]);
    }
}

fn cubic_at(p0: Pt, p1: Pt, p2: Pt, p3: Pt, t: f64) -> Pt {
    let u = 1.0 - t;
    let (a, b, c, d) = (u * u * u, 3.0 * u * u * t, 3.0 * u * t * t, t * t * t);
    [
        a * p0[0] + b * p1[0] + c * p2[0] + d * p3[0],
        a * p0[1] + b * p1[1] + c * p2[1] + d * p3[1],
    ]
}

fn add(a: Pt, b: Pt) -> Pt { [a[0] + b[0], a[1] + b[1]] }
fn sub(a: Pt, b: Pt) -> Pt { [a[0] - b[0], a[1] - b[1]] }
fn scale(a: Pt, s: f64) -> Pt { [a[0] * s, a[1] * s] }
fn norm(a: Pt) -> f64 { a[0].hypot(a[1]) }
pub fn dist(a: Pt, b: Pt) -> f64 { norm(sub(a, b)) }

/// A shape flattened once for fast inside/distance queries.
pub struct Outline {
    rings: Vec<Vec<Pt>>,  // outline first, then holes
}

impl Outline {
    pub fn new(shape: &Shape, max_len: f64) -> Outline {
        let mut rings = vec![shape.outline.flatten(max_len)];
        rings.extend(shape.holes.iter().map(|h| h.flatten(max_len)));
        Outline { rings }
    }

    pub fn from_rings(rings: Vec<Vec<Pt>>) -> Outline {
        Outline { rings }
    }

    pub fn rings(&self) -> &[Vec<Pt>] {
        &self.rings
    }

    /// Inside the panel material: odd crossing count over all rings, so
    /// holes are excluded.
    pub fn contains(&self, p: Pt) -> bool {
        let mut inside = false;
        for ring in &self.rings {
            for i in 0..ring.len() {
                let a = ring[i];
                let b = ring[(i + 1) % ring.len()];
                if (a[1] > p[1]) != (b[1] > p[1]) {
                    let x = a[0] + (p[1] - a[1]) / (b[1] - a[1]) * (b[0] - a[0]);
                    if p[0] < x {
                        inside = !inside;
                    }
                }
            }
        }
        inside
    }

    /// Distance to the nearest edge (outline or hole).
    pub fn distance_to_edge(&self, p: Pt) -> f64 {
        let mut best = f64::MAX;
        for ring in &self.rings {
            for i in 0..ring.len() {
                let a = ring[i];
                let b = ring[(i + 1) % ring.len()];
                let ab = sub(b, a);
                let len2 = ab[0] * ab[0] + ab[1] * ab[1];
                let t = if len2 > 0.0 { ((p[0] - a[0]) * ab[0] + (p[1] - a[1]) * ab[1]) / len2 } else { 0.0 };
                let q = add(a, scale(ab, t.clamp(0.0, 1.0)));
                best = best.min(dist(p, q));
            }
        }
        best
    }

    /// Material area (outline minus holes), by the shoelace formula.
    #[cfg(test)]
    pub fn area(&self) -> f64 {
        let ring_area = |r: &Vec<Pt>| -> f64 {
            (0..r.len()).map(|i| {
                let (a, b) = (r[i], r[(i + 1) % r.len()]);
                a[0] * b[1] - b[0] * a[1]
            }).sum::<f64>().abs() / 2.0
        };
        ring_area(&self.rings[0]) - self.rings[1..].iter().map(ring_area).sum::<f64>()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn area_of(path: &Path) -> f64 {
        Outline::new(&Shape { outline: path.clone(), holes: vec![] }, 0.001).area()
    }

    #[test]
    fn ellipse_area_matches_formula() {
        let a = area_of(&Path::ellipse(0.0, 0.0, 0.3, 0.2));
        let exact = PI * 0.15 * 0.1;
        assert!((a - exact).abs() / exact < 1e-3, "{a} vs {exact}");
    }

    #[test]
    fn rounded_rectangle_area_matches_formula() {
        let r = 0.03;
        let a = area_of(&Path::rounded_rectangle(0.0, 0.0, 0.3, 0.2, r));
        let exact = 0.3 * 0.2 - (4.0 - PI) * r * r;
        assert!((a - exact).abs() / exact < 1e-3, "{a} vs {exact}");
    }

    #[test]
    fn polygon_fills_its_box() {
        let pts = Path::regular_polygon(0.0, 0.0, 0.3, 0.2, 6).flatten(1.0);
        assert_eq!(pts.len(), 6);
        let max_x = pts.iter().map(|p| p[0]).fold(f64::MIN, f64::max);
        let max_y = pts.iter().map(|p| p[1]).fold(f64::MIN, f64::max);
        assert!((max_x - 0.3).abs() < 1e-12 && (max_y - 0.2).abs() < 1e-12);
    }

    #[test]
    fn holes_are_outside_the_material() {
        let shape = Shape {
            outline: Path::rectangle(0.0, 0.0, 0.3, 0.2),
            holes: vec![Path::ellipse(0.1, 0.05, 0.1, 0.1)],
        };
        let o = Outline::new(&shape, 0.002);
        assert!(o.contains([0.05, 0.1]));
        assert!(!o.contains([0.15, 0.1]));
        assert!((o.distance_to_edge([0.05, 0.1]) - 0.05).abs() < 1e-3);
        let exact = 0.06 - PI * 0.05 * 0.05;
        assert!((o.area() - exact).abs() / exact < 1e-3);
    }
}
