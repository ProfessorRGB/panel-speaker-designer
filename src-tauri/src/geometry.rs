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

/// Panel outlines, each filling a width × height bounding box.
#[derive(Clone, Debug, PartialEq)]
pub enum Panel {
    Rectangle,
    RoundedRectangle { radius: f64 },
    Ellipse,
    Polygon { sides: usize },
    /// A user-drawn outline in normalised coordinates (0..1 across the box),
    /// so it scales with the width and height.
    Custom(Path),
}

impl Panel {
    pub fn to_shape(&self, w: f64, h: f64, cutouts: &[Cutout]) -> Shape {
        let outline = match self {
            Panel::Rectangle => Path::rectangle(0.0, 0.0, w, h),
            Panel::RoundedRectangle { radius } => Path::rounded_rectangle(0.0, 0.0, w, h, *radius),
            Panel::Ellipse => Path::ellipse(0.0, 0.0, w, h),
            Panel::Polygon { sides } => Path::regular_polygon(0.0, 0.0, w, h, *sides),
            Panel::Custom(path) => path.scaled(w, h),
        };
        Shape { outline, holes: cutouts.iter().map(Cutout::to_path).collect() }
    }
}

/// A hole cut through the panel. Positions are the cutout's centre, in
/// metres from the top-left of the bounding box.
#[derive(Clone, Debug, PartialEq, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Cutout {
    /// Round hole of diameter `d`.
    Hole { x: f64, y: f64, d: f64 },
    /// Slot with semicircular ends: overall `length` end to end, `width`
    /// across, rotated `angle` degrees clockwise from horizontal.
    Slot { x: f64, y: f64, length: f64, width: f64, angle: f64 },
}

impl Cutout {
    pub fn to_path(&self) -> Path {
        match *self {
            Cutout::Hole { x, y, d } => Path::ellipse(x - d / 2.0, y - d / 2.0, d, d),
            Cutout::Slot { x, y, length, width, angle } => {
                let r = width / 2.0;
                let s = (length - width).max(0.0) / 2.0;  // half the straight part
                let k = r * KAPPA;
                let (sin, cos) = angle.to_radians().sin_cos();
                let t = |p: Pt| -> Pt { [x + p[0] * cos - p[1] * sin, y + p[0] * sin + p[1] * cos] };
                Path {
                    start: t([-s, -r]),
                    segments: vec![
                        Segment::Line { to: t([s, -r]) },
                        Segment::Cubic { c1: t([s + k, -r]), c2: t([s + r, -k]), to: t([s + r, 0.0]) },
                        Segment::Cubic { c1: t([s + r, k]), c2: t([s + k, r]), to: t([s, r]) },
                        Segment::Line { to: t([-s, r]) },
                        Segment::Cubic { c1: t([-s - k, r]), c2: t([-s - r, k]), to: t([-s - r, 0.0]) },
                        Segment::Cubic { c1: t([-s - r, -k]), c2: t([-s - k, -r]), to: t([-s, -r]) },
                    ],
                }
            }
        }
    }

    fn name(&self) -> &'static str {
        match self {
            Cutout::Hole { .. } => "Hole",
            Cutout::Slot { .. } => "Slot",
        }
    }

    fn size_ok(&self) -> bool {
        match *self {
            Cutout::Hole { d, .. } => d > 0.0,
            Cutout::Slot { length, width, .. } => width > 0.0 && length > 0.0,
        }
    }
}

/// A straight rib bonded to one face of the panel, from (x1, y1) to (x2, y2)
/// in metres, with a `width` × `height` rectangular cross-section.
#[derive(Clone, Debug, PartialEq, Deserialize)]
pub struct Stiffener {
    pub x1: f64,
    pub y1: f64,
    pub x2: f64,
    pub y2: f64,
    pub width: f64,
    pub height: f64,
    pub e: f64,    // Young's modulus along the rib [Pa]
    pub g: f64,    // shear modulus [Pa]
    pub rho: f64,  // density [kg/m³]
}

impl Stiffener {
    pub fn start(&self) -> Pt { [self.x1, self.y1] }
    pub fn end(&self) -> Pt { [self.x2, self.y2] }
    pub fn length(&self) -> f64 { dist(self.start(), self.end()) }
}

/// Smallest gap allowed between a cutout and the panel edge or another
/// cutout. Narrower ligaments can't be meshed sensibly.
pub const MIN_GAP: f64 = 0.001;

impl Shape {
    /// Checks that every cutout lies inside the outline and keeps at least
    /// `MIN_GAP` from the edge and from other cutouts. `cutouts` are the
    /// sources of `self.holes`, in order, for error messages.
    pub fn validate(&self, cutouts: &[Cutout]) -> Result<(), String> {
        let res = |path: &Path| -> Vec<Pt> {
            let o = Outline::from_rings(vec![path.flatten(f64::MAX)]);
            let size = o.rings[0].iter().zip(o.rings[0].iter().skip(1)).map(|(a, b)| dist(*a, *b)).sum::<f64>();
            path.flatten((size / 200.0).max(1e-5))
        };
        let outer = res(&self.outline);
        if outer.len() >= 3 && self_intersects(&outer) {
            return Err("The outline crosses itself".into());
        }
        if outer.len() < 3 || Outline::from_rings(vec![outer.clone()]).area().abs() < 1e-8 {
            return Err("The outline needs at least three points enclosing an area".into());
        }
        let outer_only = Outline::from_rings(vec![outer.clone()]);
        let rings: Vec<Vec<Pt>> = self.holes.iter().map(res).collect();
        for (i, ring) in rings.iter().enumerate() {
            let label = || format!("{} {}", cutouts[i].name(), i + 1);
            if !cutouts[i].size_ok() {
                return Err(format!("{} needs a positive size", label()));
            }
            if ring.iter().any(|p| !outer_only.contains(*p)) {
                return Err(format!("{} extends past the panel edge", label()));
            }
            if ring_gap(ring, &outer) < MIN_GAP {
                return Err(format!("{} is closer than 1 mm to the panel edge", label()));
            }
            for (j, other) in rings.iter().enumerate().take(i) {
                let other_only = Outline::from_rings(vec![other.clone()]);
                let this_only = Outline::from_rings(vec![ring.clone()]);
                if ring.iter().any(|p| other_only.contains(*p))
                    || other.iter().any(|p| this_only.contains(*p))
                    || ring_gap(ring, other) < MIN_GAP
                {
                    return Err(format!("{} overlaps or nearly touches {} {}", label(), cutouts[j].name(), j + 1));
                }
            }
        }
        Ok(())
    }
}

impl Shape {
    /// Checks that every stiffener is at least 2 mm long, lies inside the
    /// panel material at least `MIN_GAP` from all edges (outer and cutouts),
    /// and doesn't cross another stiffener.
    pub fn validate_stiffeners(&self, stiffeners: &[Stiffener]) -> Result<(), String> {
        let size = stiffeners.iter().map(|s| s.length()).fold(0.0, f64::max).max(1e-3);
        let outline = Outline::new(self, size / 400.0);
        for (i, s) in stiffeners.iter().enumerate() {
            let label = format!("Stiffener {}", i + 1);
            if s.length() < 0.002 {
                return Err(format!("{label} is shorter than 2 mm"));
            }
            if [s.width, s.height, s.e, s.g, s.rho].iter().any(|v| !v.is_finite() || *v <= 0.0) {
                return Err(format!("{label} needs a positive size and material values"));
            }
            let samples: Vec<Pt> = (0..=200)
                .map(|k| {
                    let t = k as f64 / 200.0;
                    [s.x1 + t * (s.x2 - s.x1), s.y1 + t * (s.y2 - s.y1)]
                })
                .collect();
            if samples.iter().any(|p| !outline.contains(*p)) {
                return Err(format!("{label} runs off the panel or across a cutout"));
            }
            if samples.iter().any(|p| outline.distance_to_edge(*p) < MIN_GAP) {
                return Err(format!("{label} is closer than 1 mm to an edge"));
            }
            for (j, o) in stiffeners.iter().enumerate().take(i) {
                if segments_cross(s.start(), s.end(), o.start(), o.end()) {
                    return Err(format!("{label} crosses stiffener {}", j + 1));
                }
            }
        }
        Ok(())
    }
}

/// True if any two non-adjacent edges of a closed ring touch or cross.
fn self_intersects(ring: &[Pt]) -> bool {
    let n = ring.len();
    let edge = |i: usize| (ring[i], ring[(i + 1) % n]);
    // Sort edges by min x so only overlapping x-ranges are compared.
    let mut order: Vec<usize> = (0..n).collect();
    let min_x = |i: usize| edge(i).0[0].min(edge(i).1[0]);
    let max_x = |i: usize| edge(i).0[0].max(edge(i).1[0]);
    order.sort_by(|&a, &b| min_x(a).total_cmp(&min_x(b)));
    for (k, &i) in order.iter().enumerate() {
        for &j in &order[k + 1..] {
            if min_x(j) > max_x(i) {
                break;
            }
            let adjacent = (i + 1) % n == j || (j + 1) % n == i;
            if !adjacent {
                let (a, b) = edge(i);
                let (c, d) = edge(j);
                if segments_cross(a, b, c, d) {
                    return true;
                }
            }
        }
    }
    false
}

/// True if two segments intersect or touch.
fn segments_cross(a: Pt, b: Pt, c: Pt, d: Pt) -> bool {
    let orient = |p: Pt, q: Pt, r: Pt| (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]);
    let on = |p: Pt, q: Pt, r: Pt| {
        r[0] >= p[0].min(q[0]) - 1e-12 && r[0] <= p[0].max(q[0]) + 1e-12
            && r[1] >= p[1].min(q[1]) - 1e-12 && r[1] <= p[1].max(q[1]) + 1e-12
    };
    let (d1, d2, d3, d4) = (orient(c, d, a), orient(c, d, b), orient(a, b, c), orient(a, b, d));
    if ((d1 > 0.0 && d2 < 0.0) || (d1 < 0.0 && d2 > 0.0)) && ((d3 > 0.0 && d4 < 0.0) || (d3 < 0.0 && d4 > 0.0)) {
        return true;
    }
    (d1 == 0.0 && on(c, d, a)) || (d2 == 0.0 && on(c, d, b)) || (d3 == 0.0 && on(a, b, c)) || (d4 == 0.0 && on(a, b, d))
}

/// Smallest distance between the vertices of one ring and the edges of another.
fn ring_gap(a: &[Pt], b: &[Pt]) -> f64 {
    let ob = Outline::from_rings(vec![b.to_vec()]);
    let oa = Outline::from_rings(vec![a.to_vec()]);
    a.iter().map(|p| ob.distance_to_edge(*p))
        .chain(b.iter().map(|p| oa.distance_to_edge(*p)))
        .fold(f64::MAX, f64::min)
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

    pub fn scaled(&self, sx: f64, sy: f64) -> Path {
        let s = |p: Pt| [p[0] * sx, p[1] * sy];
        Path {
            start: s(self.start),
            segments: self.segments.iter().map(|seg| match *seg {
                Segment::Line { to } => Segment::Line { to: s(to) },
                Segment::Cubic { c1, c2, to } => Segment::Cubic { c1: s(c1), c2: s(c2), to: s(to) },
            }).collect(),
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
        }
        // The ring is closed implicitly; don't repeat the start point.
        pts.pop();
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
        // Small holes get at least 24 edges, so a 10 mm hole in a coarse mesh
        // is still round rather than a square.
        let flatten = |path: &Path| {
            let ring = path.flatten(max_len);
            let perimeter: f64 = (0..ring.len()).map(|i| dist(ring[i], ring[(i + 1) % ring.len()])).sum();
            if perimeter / max_len < 24.0 { path.flatten(perimeter / 24.0) } else { ring }
        };
        let mut rings = vec![shape.outline.flatten(max_len)];
        rings.extend(shape.holes.iter().map(flatten));
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
        self.distance_to_rings(p, 0..self.rings.len())
    }

    /// Distances to the outer outline and to the nearest hole (MAX if none).
    pub fn distance_to_outline_and_holes(&self, p: Pt) -> (f64, f64) {
        (self.distance_to_rings(p, 0..1), self.distance_to_rings(p, 1..self.rings.len()))
    }

    fn distance_to_rings(&self, p: Pt, which: std::ops::Range<usize>) -> f64 {
        let mut best = f64::MAX;
        for ring in &self.rings[which] {
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
    fn slot_area_and_orientation() {
        let slot = Cutout::Slot { x: 0.15, y: 0.1, length: 0.1, width: 0.02, angle: 90.0 };
        let o = Outline::new(&Shape { outline: slot.to_path(), holes: vec![] }, 0.0005);
        let exact = 0.08 * 0.02 + PI * 0.01 * 0.01;
        assert!((o.area() - exact).abs() / exact < 1e-3, "{} vs {exact}", o.area());
        // Rotated 90°: long axis vertical.
        assert!(o.contains([0.15, 0.14]) && !o.contains([0.19, 0.1]));
    }

    #[test]
    fn cutouts_are_validated() {
        let ok = |c: Vec<Cutout>| Panel::Rectangle.to_shape(0.3, 0.2, &c).validate(&c);
        assert!(ok(vec![Cutout::Hole { x: 0.1, y: 0.1, d: 0.05 }]).is_ok());
        assert!(ok(vec![Cutout::Hole { x: 0.01, y: 0.1, d: 0.05 }]).unwrap_err().contains("past the panel edge"));
        assert!(ok(vec![Cutout::Hole { x: 0.0255, y: 0.1, d: 0.05 }]).unwrap_err().contains("closer than 1 mm"));
        assert!(ok(vec![Cutout::Hole { x: 0.1, y: 0.1, d: 0.05 }, Cutout::Hole { x: 0.12, y: 0.1, d: 0.05 }])
            .unwrap_err().contains("overlaps"));
        assert!(ok(vec![Cutout::Hole { x: 0.1, y: 0.1, d: 0.02 }, Cutout::Hole { x: 0.1, y: 0.1, d: 0.06 }])
            .unwrap_err().contains("overlaps"));
        // Ellipse: a hole in the bounding-box corner is off the panel.
        let c = vec![Cutout::Hole { x: 0.03, y: 0.03, d: 0.02 }];
        assert!(Panel::Ellipse.to_shape(0.3, 0.2, &c).validate(&c).is_err());
    }

    #[test]
    fn stiffeners_are_validated() {
        let rib = |x1: f64, y1: f64, x2: f64, y2: f64| Stiffener {
            x1, y1, x2, y2, width: 0.005, height: 0.01, e: 10e9, g: 0.6e9, rho: 450.0,
        };
        let c = vec![Cutout::Hole { x: 0.15, y: 0.1, d: 0.04 }];
        let shape = Panel::Rectangle.to_shape(0.3, 0.2, &c);
        assert!(shape.validate_stiffeners(&[rib(0.02, 0.05, 0.28, 0.05)]).is_ok());
        assert!(shape.validate_stiffeners(&[rib(0.0005, 0.05, 0.28, 0.05)]).unwrap_err().contains("closer than 1 mm"));
        assert!(shape.validate_stiffeners(&[rib(0.02, 0.1, 0.28, 0.1)]).unwrap_err().contains("across a cutout"));
        assert!(shape.validate_stiffeners(&[rib(0.02, 0.05, 0.28, 0.05), rib(0.05, 0.02, 0.05, 0.18)])
            .unwrap_err().contains("crosses"));
    }

    #[test]
    fn custom_outlines_scale_and_are_validated() {
        // A normalised diamond scales to the box.
        let diamond = Path {
            start: [0.5, 0.0],
            segments: vec![
                Segment::Line { to: [1.0, 0.5] },
                Segment::Line { to: [0.5, 1.0] },
                Segment::Line { to: [0.0, 0.5] },
            ],
        };
        let shape = Panel::Custom(diamond).to_shape(0.3, 0.2, &[]);
        assert!(shape.validate(&[]).is_ok());
        assert!((Outline::new(&shape, 0.01).area() - 0.03).abs() < 1e-12);
        // A bow-tie crosses itself.
        let bowtie = Path {
            start: [0.0, 0.0],
            segments: vec![
                Segment::Line { to: [1.0, 1.0] },
                Segment::Line { to: [1.0, 0.0] },
                Segment::Line { to: [0.0, 1.0] },
            ],
        };
        let err = Panel::Custom(bowtie).to_shape(0.3, 0.2, &[]).validate(&[]).unwrap_err();
        assert!(err.contains("crosses itself"), "{err}");
        // A curve looping back over a straight edge crosses it too.
        let loopy = Path {
            start: [0.0, 0.5],
            segments: vec![
                Segment::Line { to: [1.0, 0.5] },
                Segment::Cubic { c1: [1.0, 1.0], c2: [-0.5, -0.5], to: [0.0, 0.5] },
            ],
        };
        assert!(Panel::Custom(loopy).to_shape(0.3, 0.2, &[]).validate(&[]).is_err());
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
