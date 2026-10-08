// Finite-element modal analysis of a thin plate of arbitrary shape.
//
// Element : DKT (Discrete Kirchhoff Triangle). Slopes s_x = ∂w/∂x and
//           s_y = ∂w/∂y are interpolated quadratically over the triangle,
//           with the Kirchhoff condition imposed at the corners and edge
//           midpoints (cubic w along each edge, linear normal slope). Nodal
//           unknowns are (w, ∂w/∂x, ∂w/∂y). Supports orthotropic bending
//           stiffness with principal axes along x and y.
// Mass    : lumped translational mass, plus the plate's (small) rotary
//           inertia ρh³/12 on the slope unknowns so the mass matrix is
//           positive definite.
// Solver  : shift-invert block Lanczos with full reorthogonalisation. The
//           sparse shifted stiffness is factored once with a supernodal
//           Cholesky. The block size handles the repeated frequencies of
//           symmetric shapes (a circle's modes come in pairs).

use crate::geometry::Pt;
use crate::mesh::Mesh;
use crate::plate::Rigidities;
use faer::linalg::matmul::matmul;
use faer::linalg::solvers::Solve;
use faer::sparse::{SparseColMat, Triplet};
use faer::{Accum, Mat, Par, Side};
use std::f64::consts::PI;

/// A stiffener as a chain of mesh nodes along a straight line, modelled as
/// Euler-Bernoulli beam elements (bending along the rib) plus St Venant
/// torsion, sharing the plate's nodal unknowns.
pub struct Beam {
    pub nodes: Vec<usize>,    // ordered along the rib
    pub ei: f64,              // bending stiffness added to the panel [N·m²]
    pub gj: f64,              // torsional stiffness [N·m²]
    pub mass_per_len: f64,    // [kg/m]
}

/// Stiffness of one beam element between nodes at `a` and `b`, in the plate
/// unknowns [w, ∂w/∂x, ∂w/∂y] of each node. Along the beam, the bending
/// slope is the directional derivative of w and the twist is the slope
/// across it.
fn beam_stiffness(a: Pt, b: Pt, ei: f64, gj: f64) -> [[f64; 6]; 6] {
    let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
    let len = dx.hypot(dy);
    let (c, s) = (dx / len, dy / len);
    // Local unknowns: [w1, θs1, θn1, w2, θs2, θn2].
    let mut kl = [[0.0; 6]; 6];
    let k = ei / len.powi(3);
    let (l, l2) = (len, len * len);
    let bend = [
        [12.0, 6.0 * l, -12.0, 6.0 * l],
        [6.0 * l, 4.0 * l2, -6.0 * l, 2.0 * l2],
        [-12.0, -6.0 * l, 12.0, -6.0 * l],
        [6.0 * l, 2.0 * l2, -6.0 * l, 4.0 * l2],
    ];
    let bidx = [0, 1, 3, 4];
    for (r, &i) in bidx.iter().enumerate() {
        for (q, &j) in bidx.iter().enumerate() {
            kl[i][j] = k * bend[r][q];
        }
    }
    let kt = gj / len;
    kl[2][2] = kt;
    kl[5][5] = kt;
    kl[2][5] = -kt;
    kl[5][2] = -kt;
    // Local from plate unknowns: θs = c·∂w/∂x + s·∂w/∂y, θn = −s·∂w/∂x + c·∂w/∂y.
    let mut t = [[0.0; 6]; 6];
    for n in 0..2 {
        let o = 3 * n;
        t[o][o] = 1.0;
        t[o + 1][o + 1] = c;
        t[o + 1][o + 2] = s;
        t[o + 2][o + 1] = -s;
        t[o + 2][o + 2] = c;
    }
    let mut kg = [[0.0; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            kg[i][j] = (0..6).map(|p| (0..6).map(|q| t[p][i] * kl[p][q] * t[q][j]).sum::<f64>()).sum();
        }
    }
    kg
}

pub struct FemMode {
    pub freq: f64,
    pub w: Vec<f64>,  // transverse displacement at each mesh node
}

// ── Element ──────────────────────────────────────────────────────────────────

type Row9 = [f64; 9];

/// DKT stiffness matrix for one triangle (counter-clockwise nodes).
fn dkt_stiffness(p: [Pt; 3], d: &[[f64; 3]; 3]) -> [[f64; 9]; 9] {
    let [(x1, y1), (x2, y2), (x3, y3)] = p.map(|q| (q[0], q[1]));
    let two_a = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
    let area = two_a / 2.0;
    // ∂L_i/∂x and ∂L_i/∂y for the area coordinates.
    let dl_dx = [(y2 - y3) / two_a, (y3 - y1) / two_a, (y1 - y2) / two_a];
    let dl_dy = [(x3 - x2) / two_a, (x1 - x3) / two_a, (x2 - x1) / two_a];

    // Slopes at the six quadratic nodes as linear functions of the nine
    // element unknowns [w1, sx1, sy1, w2, ...].
    let mut tx = [[0.0; 9]; 6];
    let mut ty = [[0.0; 9]; 6];
    for i in 0..3 {
        tx[i][3 * i + 1] = 1.0;
        ty[i][3 * i + 2] = 1.0;
    }
    const EDGES: [(usize, usize); 3] = [(0, 1), (1, 2), (2, 0)];
    for (e, &(i, j)) in EDGES.iter().enumerate() {
        let (dx, dy) = (p[j][0] - p[i][0], p[j][1] - p[i][1]);
        let len = dx.hypot(dy);
        let (c, s) = (dx / len, dy / len);
        // Tangential slope at the midpoint from a cubic w along the edge;
        // normal slope varies linearly.
        let mut ss: Row9 = [0.0; 9];
        let mut sn: Row9 = [0.0; 9];
        ss[3 * i] = -1.5 / len;
        ss[3 * j] = 1.5 / len;
        for k in [i, j] {
            ss[3 * k + 1] = -0.25 * c;
            ss[3 * k + 2] = -0.25 * s;
            sn[3 * k + 1] = -0.5 * s;
            sn[3 * k + 2] = 0.5 * c;
        }
        for q in 0..9 {
            tx[3 + e][q] = c * ss[q] - s * sn[q];
            ty[3 + e][q] = s * ss[q] + c * sn[q];
        }
    }

    // Curvatures are linear, so the 3 edge-midpoint rule integrates BᵀDB exactly.
    let mut k = [[0.0; 9]; 9];
    for l in [[0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5]] {
        // Quadratic shape function derivatives: corners L(2L−1), midsides 4LᵢLⱼ.
        let mut dn_dx = [0.0; 6];
        let mut dn_dy = [0.0; 6];
        for i in 0..3 {
            dn_dx[i] = (4.0 * l[i] - 1.0) * dl_dx[i];
            dn_dy[i] = (4.0 * l[i] - 1.0) * dl_dy[i];
        }
        for (e, &(i, j)) in EDGES.iter().enumerate() {
            dn_dx[3 + e] = 4.0 * (l[i] * dl_dx[j] + l[j] * dl_dx[i]);
            dn_dy[3 + e] = 4.0 * (l[i] * dl_dy[j] + l[j] * dl_dy[i]);
        }
        // B rows: κ = [∂s_x/∂x, ∂s_y/∂y, ∂s_x/∂y + ∂s_y/∂x]
        let mut b = [[0.0; 9]; 3];
        for a in 0..6 {
            for q in 0..9 {
                b[0][q] += dn_dx[a] * tx[a][q];
                b[1][q] += dn_dy[a] * ty[a][q];
                b[2][q] += dn_dy[a] * tx[a][q] + dn_dx[a] * ty[a][q];
            }
        }
        let weight = area / 3.0;
        for r in 0..9 {
            let db: [f64; 3] = std::array::from_fn(|m| (0..3).map(|n| d[m][n] * b[n][r]).sum());
            for c in 0..9 {
                k[c][r] += weight * (0..3).map(|m| b[m][c] * db[m]).sum::<f64>();
            }
        }
    }
    k
}

// ── Assembly and solve ───────────────────────────────────────────────────────

/// Modes up to `freq_max`, sorted by frequency. Edges are free except at
/// `pinned` nodes, where w = 0 (simply supported). Rigid-body modes of a
/// free plate are excluded. Shapes are normalised to RMS 0.5 over the plate.
pub fn solve(
    mesh: &Mesh,
    d: &Rigidities,
    rho_h: f64,
    rho_h3_12: f64,
    freq_max: f64,
    pinned: &[usize],
    beams: &[Beam],
) -> Result<Vec<FemMode>, String> {
    let nn = mesh.nodes.len();
    let n = 3 * nn;
    let dmat = [[d.d11, d.d12, 0.0], [d.d12, d.d22, 0.0], [0.0, 0.0, d.d66]];

    // Lumped mass (diagonal).
    let mut mass = vec![0.0; n];
    let mut node_area = vec![0.0; nn];
    for (t, tri) in mesh.tris.iter().enumerate() {
        let a3 = mesh.triangle_area(t) / 3.0;
        for &i in tri {
            node_area[i] += a3;
            mass[3 * i] += rho_h * a3;
            mass[3 * i + 1] += rho_h3_12 * a3;
            mass[3 * i + 2] += rho_h3_12 * a3;
        }
    }
    let total_area: f64 = node_area.iter().sum();
    for beam in beams {
        for pair in beam.nodes.windows(2) {
            let half = beam.mass_per_len * crate::geometry::dist(mesh.nodes[pair[0]], mesh.nodes[pair[1]]) / 2.0;
            mass[3 * pair[0]] += half;
            mass[3 * pair[1]] += half;
        }
    }

    // Shift so that K + σM is positive definite despite the rigid-body modes.
    // σ near the expected fundamental keeps the shifted problem well scaled.
    let d_min = d.d11.min(d.d22).min(d.d66);
    let f_scale = (d_min / rho_h).sqrt() / total_area / (2.0 * PI);
    let sigma = (2.0 * PI * f_scale).powi(2);

    let mut triplets = Vec::with_capacity(mesh.tris.len() * 45 + n);
    let mut k_diag_max: f64 = 0.0;
    for tri in &mesh.tris {
        let ke = dkt_stiffness(tri.map(|i| mesh.nodes[i]), &dmat);
        k_diag_max = (0..9).fold(k_diag_max, |a, q| a.max(ke[q][q]));
        for a in 0..9 {
            for b in 0..9 {
                let (r, c) = (3 * tri[a / 3] + a % 3, 3 * tri[b / 3] + b % 3);
                if r >= c {
                    triplets.push(Triplet::new(r, c, ke[a][b]));
                }
            }
        }
    }
    for beam in beams {
        for pair in beam.nodes.windows(2) {
            let ke = beam_stiffness(mesh.nodes[pair[0]], mesh.nodes[pair[1]], beam.ei, beam.gj);
            k_diag_max = (0..6).fold(k_diag_max, |a, q| a.max(ke[q][q]));
            for a in 0..6 {
                for b in 0..6 {
                    let (r, c) = (3 * pair[a / 3] + a % 3, 3 * pair[b / 3] + b % 3);
                    if r >= c {
                        triplets.push(Triplet::new(r, c, ke[a][b]));
                    }
                }
            }
        }
    }
    for (i, m) in mass.iter().enumerate() {
        triplets.push(Triplet::new(i, i, sigma * m));
    }
    // Simply-supported nodes: a stiff spring on w. Large enough to act as a
    // constraint, small enough to keep the factorisation well conditioned.
    for &i in pinned {
        triplets.push(Triplet::new(3 * i, 3 * i, 1e6 * k_diag_max));
    }
    let shifted = SparseColMat::<usize, f64>::try_new_from_triplets(n, n, &triplets)
        .map_err(|e| format!("Assembly failed: {e:?}"))?;
    let llt = shifted.sp_cholesky(Side::Lower)
        .map_err(|e| format!("Stiffness factorisation failed: {e:?}"))?;

    // Estimated mode count below freq_max (asymptotic modal density of a
    // plate, with an allowance for edge effects), to size the Krylov space.
    let omega_max = 2.0 * PI * freq_max;
    let d_eff = (d.d11 * d.d22).sqrt();
    let estimate = (total_area * omega_max / (4.0 * PI) * (rho_h / d_eff).sqrt() * 1.3) as usize + 12;

    let lambda_max = omega_max * omega_max;
    let pairs = block_lanczos(
        |x: &mut Mat<f64>| {
            for j in 0..x.ncols() {
                for i in 0..n {
                    x[(i, j)] *= mass[i];
                }
            }
            llt.solve_in_place(x.as_mut());
        },
        &mass,
        sigma,
        lambda_max,
        estimate,
    )?;

    // Rigid-body modes have (numerically) zero frequency; elastic modes sit
    // well above this threshold.
    let f_rigid = 0.05 * (d_min / rho_h).sqrt() / total_area;
    let mut modes: Vec<FemMode> = pairs.into_iter()
        .map(|(lambda, v)| (lambda.max(0.0).sqrt() / (2.0 * PI), v))
        .filter(|(f, _)| *f > f_rigid && *f <= freq_max)
        .map(|(freq, v)| {
            let mut w: Vec<f64> = (0..nn).map(|i| v[3 * i]).collect();
            // Normalise to RMS 0.5 over the plate, and make the largest
            // displacement positive so signs are deterministic.
            let ms: f64 = w.iter().zip(&node_area).map(|(x, a)| x * x * a).sum::<f64>() / total_area;
            let peak = w.iter().copied().fold(0.0, |a: f64, x| if x.abs() > a.abs() { x } else { a });
            let s = 0.5 / ms.sqrt() * peak.signum();
            w.iter_mut().for_each(|x| *x *= s);
            FemMode { freq, w }
        })
        .collect();
    modes.sort_by(|a, b| a.freq.total_cmp(&b.freq));
    Ok(modes)
}

// ── Eigensolver ──────────────────────────────────────────────────────────────

const BLOCK: usize = 8;

/// Eigenpairs (λ, x) of K x = λ M x with λ ≤ `lambda_max`, given `apply_op`
/// that overwrites X with (K + σM)⁻¹ M X. Works in the M-inner product, where
/// the operator is self-adjoint with eigenvalues ν = 1/(λ + σ).
fn block_lanczos(
    apply_op: impl Fn(&mut Mat<f64>),
    mass: &[f64],
    sigma: f64,
    lambda_max: f64,
    estimate: usize,
) -> Result<Vec<(f64, Vec<f64>)>, String> {
    let n = mass.len();
    let nu_min = 1.0 / (lambda_max + sigma);
    let mut capacity = (2 * estimate + 8 * BLOCK).min(n);
    let max_basis = (6 * estimate + 64 * BLOCK).min(n);
    let mut basis = Mat::<f64>::zeros(n, capacity);
    let mut used = 0;
    // Block tridiagonal projection: diagonal blocks A_j, sub-diagonal B_j.
    let mut diag_blocks: Vec<Mat<f64>> = Vec::new();
    let mut sub_blocks: Vec<Mat<f64>> = Vec::new();
    let mut rng = Lcg(0x9E37_79B9_7F4A_7C15);

    // Starting block: random, M-orthonormalised.
    let mut start = Mat::<f64>::from_fn(n, BLOCK, |_, _| rng.next());
    orthonormalise(&mut start, &basis, used, mass, &mut rng);
    for c in 0..BLOCK {
        basis.col_mut(c).copy_from(start.col(c));
    }
    used = BLOCK;
    let mut next_check = (estimate + 4 * BLOCK).min(n);

    loop {
        if used + BLOCK > capacity {
            // Convergence normally needs about twice the wanted count; far
            // beyond that something is wrong, so stop rather than grow the
            // basis toward the full problem size.
            if capacity + BLOCK > n || capacity >= max_basis {
                return Err("The eigen-solver didn't converge for this panel; try a lower max frequency".into());
            }
            // Grow the Krylov space.
            let new_cap = (capacity + capacity / 2).min(n);
            let mut bigger = Mat::<f64>::zeros(n, new_cap);
            bigger.as_mut().subcols_mut(0, used).copy_from(basis.as_ref().subcols(0, used));
            basis = bigger;
            capacity = new_cap;
        }

        // Expand: W = Op V_j, A_j = V_jᵀ M W, then the next block from W
        // orthogonalised against everything so far. Full reorthogonalisation
        // removes the A_j and B_jᵀ terms and any loss of orthogonality.
        let vj = basis.as_ref().subcols(used - BLOCK, BLOCK).to_owned();
        let mut w = vj.clone();
        apply_op(&mut w);
        let mw = scale_rows(&w, mass);
        let mut a = Mat::<f64>::zeros(BLOCK, BLOCK);
        matmul(a.as_mut(), Accum::Replace, vj.transpose(), mw.as_ref(), 1.0, Par::Seq);
        diag_blocks.push(Mat::from_fn(BLOCK, BLOCK, |r, c| 0.5 * (a[(r, c)] + a[(c, r)])));
        let b = orthonormalise(&mut w, &basis, used, mass, &mut rng);
        basis.as_mut().subcols_mut(used, BLOCK).copy_from(w.as_ref());
        sub_blocks.push(b);
        used += BLOCK;

        if used < next_check {
            continue;
        }
        next_check = used + (used / 4).max(4 * BLOCK);

        // Rayleigh-Ritz on the block tridiagonal projection of V_0 … V_j.
        let m = diag_blocks.len() * BLOCK;
        let t = Mat::<f64>::from_fn(m, m, |r, c| {
            let (br, bc) = (r / BLOCK, c / BLOCK);
            if br == bc {
                diag_blocks[br][(r % BLOCK, c % BLOCK)]
            } else if br == bc + 1 {
                sub_blocks[bc][(r % BLOCK, c % BLOCK)]
            } else if bc == br + 1 {
                sub_blocks[br][(c % BLOCK, r % BLOCK)]
            } else {
                0.0
            }
        });
        let evd = t.self_adjoint_eigen(Side::Lower).map_err(|e| format!("Eigen solve failed: {e:?}"))?;
        let (s, u) = (evd.S(), evd.U());
        // Residual of Ritz pair k: ‖B_{j+1} · (last block of y_k)‖.
        let last_b = &sub_blocks[diag_blocks.len() - 1];
        let residual = |k: usize| -> f64 {
            (0..BLOCK).map(|r| {
                (0..BLOCK).map(|c| last_b[(r, c)] * u[(m - BLOCK + c, k)]).sum::<f64>().powi(2)
            }).sum::<f64>().sqrt()
        };
        // A residual r gives an eigenvalue error of order r²/gap, so 1e-6 is
        // far below anything that matters here.
        let converged_to = |k: usize, tol: f64| residual(k) <= tol * s[k].abs().max(nu_min);
        let converged = |k: usize| converged_to(k, 1e-6);
        // Eigenvalues ascend; the wanted ν (largest) are at the end.
        let wanted: Vec<usize> = (0..m).rev().take_while(|&k| s[k] >= nu_min).collect();
        // Also require the first unwanted Ritz value to have roughly settled,
        // so no wanted eigenvalue is still approaching from below the cut-off.
        let boundary_ok = wanted.len() < m && converged_to(m - 1 - wanted.len(), 1e-3);
        if !(boundary_ok && wanted.iter().all(|&k| converged(k))) {
            continue;
        }

        let mut y = Mat::<f64>::zeros(m, wanted.len());
        for (col, &k) in wanted.iter().enumerate() {
            for r in 0..m {
                y[(r, col)] = u[(r, k)];
            }
        }
        let mut x = Mat::<f64>::zeros(n, wanted.len());
        matmul(x.as_mut(), Accum::Replace, basis.as_ref().subcols(0, m), y.as_ref(), 1.0, Par::rayon(0));
        return Ok(wanted.iter().enumerate()
            .map(|(col, &k)| (1.0 / s[k] - sigma, (0..n).map(|i| x[(i, col)]).collect()))
            .collect());
    }
}

fn scale_rows(x: &Mat<f64>, mass: &[f64]) -> Mat<f64> {
    Mat::from_fn(x.nrows(), x.ncols(), |i, j| x[(i, j)] * mass[i])
}

/// M-orthonormalises the columns of `x` against the first `used` columns of
/// `basis` and each other (two passes of block Gram-Schmidt, then modified
/// Gram-Schmidt within the block). Returns R with x_in ≈ basis·(…) + x_out·R.
/// Columns that collapse (deflation) are replaced with fresh random vectors.
fn orthonormalise(x: &mut Mat<f64>, basis: &Mat<f64>, used: usize, mass: &[f64], rng: &mut Lcg) -> Mat<f64> {
    let (n, b) = (x.nrows(), x.ncols());
    let v = basis.as_ref().subcols(0, used);
    let project_out = |x: &mut Mat<f64>| {
        if used == 0 {
            return;
        }
        for _ in 0..2 {
            let mx = scale_rows(x, mass);
            let mut c = Mat::<f64>::zeros(used, x.ncols());
            matmul(c.as_mut(), Accum::Replace, v.transpose(), mx.as_ref(), 1.0, Par::rayon(0));
            matmul(x.as_mut(), Accum::Add, v, c.as_ref(), -1.0, Par::rayon(0));
        }
    };
    let m_dot = |a: faer::ColRef<f64>, b: faer::ColRef<f64>| -> f64 {
        (0..n).map(|i| a[i] * mass[i] * b[i]).sum()
    };

    let scale_ref = (0..b).map(|c| m_dot(x.col(c), x.col(c)).sqrt()).fold(0.0, f64::max).max(1e-300);
    project_out(x);
    let mut r = Mat::<f64>::zeros(b, b);
    for c in 0..b {
        for _ in 0..2 {
            for p in 0..c {
                let coef = m_dot(x.col(p), x.col(c));
                r[(p, c)] += coef;
                let xp = x.col(p).to_owned();
                for i in 0..n {
                    x[(i, c)] -= coef * xp[i];
                }
            }
        }
        let norm = m_dot(x.col(c), x.col(c)).sqrt();
        if norm > 1e-10 * scale_ref {
            r[(c, c)] = norm;
            for i in 0..n {
                x[(i, c)] /= norm;
            }
        } else {
            // Deflation: the Krylov space has (numerically) closed in this
            // direction; continue with a random direction.
            r[(c, c)] = 0.0;
            let mut fresh = Mat::<f64>::from_fn(n, 1, |_, _| rng.next());
            project_out(&mut fresh);
            for p in 0..c {
                let coef = m_dot(x.col(p), fresh.col(0));
                for i in 0..n {
                    fresh[(i, 0)] -= coef * x[(i, p)];
                }
            }
            let fnorm = m_dot(fresh.col(0), fresh.col(0)).sqrt();
            for i in 0..n {
                x[(i, c)] = fresh[(i, 0)] / fnorm;
            }
        }
    }
    r
}

/// Small deterministic generator for start vectors, so results are repeatable.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64 / (1u64 << 53) as f64) - 0.5
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{Outline, Path, Shape};
    use crate::mesh;
    use crate::plate::{self, Boundary, Plate, SolveKey};

    fn acrylic(lx: f64, ly: f64) -> Plate {
        Plate { lx, ly, h: 0.003, ex: 3.2e9, ey: 3.2e9, g: 3.2e9 / 2.74, nu_xy: 0.37, rho: 1190.0 }
    }

    fn fem_modes(p: &Plate, path: Path, h: f64, freq_max: f64) -> Vec<FemMode> {
        let outline = Outline::new(&Shape { outline: path, holes: vec![] }, h);
        let m = mesh::mesh(&outline, h, &[]).unwrap();
        solve(&m, &p.rigidities(), p.rho * p.h, p.rho * p.h.powi(3) / 12.0, freq_max, &[], &[]).unwrap()
    }

    #[test]
    fn rectangle_matches_ritz_solver() {
        let p = acrylic(0.3, 0.2);
        let fe = fem_modes(&p, Path::rectangle(0.0, 0.0, 0.3, 0.2), 0.006, 2000.0);
        let ritz = plate::solve(SolveKey::for_freq(p, Boundary::Free, 2500.0));
        assert!(fe.len() >= 30, "only {} modes", fe.len());
        for (k, (a, b)) in fe.iter().zip(&ritz.modes).enumerate() {
            let err = (a.freq - b.freq) / b.freq;
            assert!(err.abs() < 0.02, "mode {k}: FEM {:.1} Hz vs Ritz {:.1} Hz", a.freq, b.freq);
        }
    }

    #[test]
    fn rectangle_with_orthotropic_material_matches_ritz_solver() {
        let p = Plate { ex: 3.0e9, ey: 0.09e9, g: 0.12e9, nu_xy: 0.3, rho: 130.0, ..acrylic(0.3, 0.2) };
        let fe = fem_modes(&p, Path::rectangle(0.0, 0.0, 0.3, 0.2), 0.006, 1500.0);
        let ritz = plate::solve(SolveKey::for_freq(p, Boundary::Free, 2000.0));
        for (k, (a, b)) in fe.iter().zip(&ritz.modes).enumerate().take(25) {
            let err = (a.freq - b.freq) / b.freq;
            assert!(err.abs() < 0.02, "mode {k}: FEM {:.1} Hz vs Ritz {:.1} Hz", a.freq, b.freq);
        }
    }

    #[test]
    fn simply_supported_rectangle_matches_exact_solution() {
        let p = acrylic(0.3, 0.2);
        let h = 0.006;
        let outline = Outline::new(&Shape { outline: Path::rectangle(0.0, 0.0, 0.3, 0.2), holes: vec![] }, h);
        let m = mesh::mesh(&outline, h, &[]).unwrap();
        let pinned: Vec<usize> = (0..m.nodes.len()).filter(|&i| outline.distance_to_edge(m.nodes[i]) < 1e-9).collect();
        let fe = solve(&m, &p.rigidities(), p.rho * p.h, p.rho * p.h.powi(3) / 12.0, 1500.0, &pinned, &[]).unwrap();
        let exact = plate::solve(SolveKey::for_freq(p, Boundary::SimplySupported, 2000.0));
        for (k, (a, b)) in fe.iter().zip(&exact.modes).enumerate() {
            let err = (a.freq - b.freq) / b.freq;
            assert!(err.abs() < 0.02, "mode {k}: FEM {:.1} Hz vs exact {:.1} Hz", a.freq, b.freq);
        }
    }

    #[test]
    fn free_circular_plate_matches_leissa() {
        // Leissa, "Vibration of Plates", free circular plate, ν = 0.33:
        // λ = ω a² √(ρh/D) with a the radius. (2,0) and (3,0) are double modes.
        let mut p = acrylic(0.3, 0.3);
        p.nu_xy = 0.33;
        p.g = p.ex / (2.0 * 1.33);
        let fe = fem_modes(&p, Path::ellipse(0.0, 0.0, 0.3, 0.3), 0.006, 1200.0);
        let d = p.rigidities().d11;
        let lam = |f: f64| 2.0 * PI * f * 0.15f64.powi(2) * (p.rho * p.h / d).sqrt();
        let got: Vec<f64> = fe.iter().map(|m| lam(m.freq)).collect();
        let expected = [5.253, 5.253, 9.084, 12.23, 12.23];
        assert!(got.len() >= expected.len(), "{got:?}");
        for (g, e) in got.iter().zip(expected) {
            assert!((g - e).abs() / e < 0.02, "λ = {g:.3}, expected {e} (all: {:?})", &got[..6.min(got.len())]);
        }
    }
}
