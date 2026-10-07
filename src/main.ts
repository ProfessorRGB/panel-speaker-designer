import { invoke } from "@tauri-apps/api/core";
import { ModeDensityChart, ResponseChart, Series } from "./charts";

// ── Types ────────────────────────────────────────────────────────────────────

// Cutouts as sent to Rust: positions and sizes in metres, angle in degrees.
type Cutout =
  | { kind: "hole"; x: number; y: number; d: number }
  | { kind: "slot"; x: number; y: number; length: number; width: number; angle: number };

interface PanelParams {
  shape: string;
  corner_r: number;
  sides: number;
  cutouts: Cutout[];
  lx: number;
  ly: number;
  h: number;
  ex: number;
  ey: number;
  g: number;
  nu: number;
  rho: number;
  boundary: string;
  freq_max: number;
  grid_n: number;
  exciter_d: number;
  eta: number;
  score: string;
}

interface ModeInfo {
  m: number | null;  // nodal-line counts; null where not meaningful
  n: number | null;
  freq: number;
}

// Grid cells: candidate positions, edge margin, and outside the panel.
const REGION_MARGIN = 1;
const REGION_OUTSIDE = 2;

interface CalculationResult {
  grid: (number | null)[];      // null outside the panel
  grid_raw: (number | null)[];
  region: number[];
  outline: [number, number][][];  // outline and hole rings, normalised
  grid_n: number;
  modes: ModeInfo[];
  mode_count: number;
  optimal_x: number;
  optimal_y: number;
  optimal_score_raw: number;
  bands: number[];
  response_opt: number[];
  raggedness_opt: number;
  solver: "analytic" | "fea";
  truncated_above: number | null;
}

// ── Material presets ─────────────────────────────────────────────────────────
// Moduli in MPa. Isotropic presets derive E_y and G from E_x and ν.
// Wood values are typical; real sheets vary widely, so measure if you can.
// Loss factors are for bare panels; mounting and surrounds add damping.

interface Material {
  ex: number;
  rho: number;
  nu: number;
  eta: number;
  ey?: number;  // set for orthotropic materials
  g?: number;
}

const MATERIALS: Record<string, Material> = {
  xps:      { ex: 20,    rho: 32,   nu: 0.35, eta: 0.05 },                     // XPS foam
  eps:      { ex: 5,     rho: 20,   nu: 0.10, eta: 0.06 },                     // EPS foam
  balsa:    { ex: 3000,  rho: 130,  nu: 0.30, eta: 0.03, ey: 90,   g: 120 },   // Balsa, grain along x
  birch:    { ex: 10000, rho: 680,  nu: 0.07, eta: 0.03, ey: 5500, g: 620 },   // Birch plywood, face grain along x
  acrylic:  { ex: 3200,  rho: 1190, nu: 0.37, eta: 0.04 },                     // PMMA
  aluminum: { ex: 69000, rho: 2700, nu: 0.33, eta: 0.005 },                    // Aluminium
  carbon:   { ex: 70000, rho: 1600, nu: 0.10, eta: 0.01 },                     // CFRP, quasi-isotropic estimate
};

// Matches the "clicked position" series colour in the response chart.
const PROBE_COLOR = "#d95926";

// Resolution of the mode shape used to trace node lines.
const SHAPE_GRID_N = 120;

// ── State ────────────────────────────────────────────────────────────────────

let lastResult: CalculationResult | null = null;
let lastParams: PanelParams | null = null;
let selectedModeIdx = -1;
let selectedShape: (number | null)[] | null = null;
let calcTimer: ReturnType<typeof setTimeout> | null = null;
let requestId = 0;  // newer calculations supersede older in-flight ones

// A position the user clicked, to compare its response with the optimum.
interface Probe { x: number; y: number; db: number[]; raggedness: number }
let probe: Probe | null = null;

// Cutouts as edited in the sidebar, in mm (angle in degrees).
let cutouts: Cutout[] = [];

// ── DOM refs ─────────────────────────────────────────────────────────────────

const $ = (id: string) => document.getElementById(id)!;

const selectShape    = $("shape")      as HTMLSelectElement;
const inputCornerR   = $("corner-r")   as HTMLInputElement;
const inputSides     = $("sides")      as HTMLInputElement;
const inputLx        = $("lx")         as HTMLInputElement;
const inputLy        = $("ly")         as HTMLInputElement;
const inputH         = $("h")          as HTMLInputElement;
const checkIsotropic = $("isotropic")  as HTMLInputElement;
const inputEx        = $("ex")         as HTMLInputElement;
const inputEy        = $("ey")         as HTMLInputElement;
const inputG         = $("g")          as HTMLInputElement;
const inputRho       = $("rho")        as HTMLInputElement;
const inputNu        = $("nu")         as HTMLInputElement;
const selectMaterial = $("material-preset") as HTMLSelectElement;
const selectBoundary = $("boundary")   as HTMLSelectElement;
const inputFreqMax   = $("freq-max")   as HTMLInputElement;
const inputExciterD  = $("exciter-d")  as HTMLInputElement;
const inputEta       = $("eta")        as HTMLInputElement;
const selectScore    = $("score")      as HTMLSelectElement;
const inputGridN     = $("grid-n")     as HTMLInputElement;
const selectMode     = $("mode-select") as HTMLSelectElement;

const canvas         = $("heatmap")    as HTMLCanvasElement;
const canvasWrap     = $("canvas-wrap");
const tooltip        = $("hover-tooltip");
const statusText     = $("status-text");
const badgeModes     = $("mode-count-badge");
const valOptimal     = $("val-optimal");
const valCursor      = $("val-cursor");
const valModes       = $("val-modes");
const valF1          = $("val-f1");
const valRagged      = $("val-ragged");
const responseLegend = $("response-legend");
const densityGaps    = $("density-gaps");
const chartTooltip   = $("chart-tooltip");

const responseChart = new ResponseChart($("response-chart") as HTMLCanvasElement, chartTooltip);
const densityChart  = new ModeDensityChart($("density-chart") as HTMLCanvasElement, chartTooltip);

// ── Colormap ─────────────────────────────────────────────────────────────────
// Perceptual: dark-blue → blue → teal → green → yellow → red

const COLORMAP: [number, [number, number, number]][] = [
  [0.00, [13,  17,  45]],
  [0.20, [32,  80, 180]],
  [0.40, [30, 160, 140]],
  [0.60, [50, 200,  60]],
  [0.80, [230, 200,  20]],
  [1.00, [210,  40,  20]],
];

function sampleColormap(t: number): [number, number, number] {
  t = Math.max(0, Math.min(1, t));
  for (let i = 1; i < COLORMAP.length; i++) {
    const [t0, c0] = COLORMAP[i - 1];
    const [t1, c1] = COLORMAP[i];
    if (t <= t1) {
      const f = (t - t0) / (t1 - t0);
      return [
        Math.round(c0[0] + f * (c1[0] - c0[0])),
        Math.round(c0[1] + f * (c1[1] - c0[1])),
        Math.round(c0[2] + f * (c1[2] - c0[2])),
      ];
    }
  }
  return COLORMAP[COLORMAP.length - 1][1];
}

// ── Canvas rendering ─────────────────────────────────────────────────────────

function getCanvasSize(lx: number, ly: number): { cw: number; ch: number; scale: number } {
  const wrap = canvasWrap.getBoundingClientRect();
  const pad = 48;
  const availW = wrap.width - pad;
  const availH = wrap.height - pad;
  const scale = Math.min(availW / lx, availH / ly);
  return {
    cw: Math.round(lx * scale) + pad,
    ch: Math.round(ly * scale) + pad,
    scale,
  };
}

function render() {
  if (!lastResult) return;

  const lxMm = parseFloat(inputLx.value);
  const lyMm = parseFloat(inputLy.value);
  const { cw, ch, scale } = getCanvasSize(lxMm, lyMm);

  canvas.width  = cw;
  canvas.height = ch;

  const ctx = canvas.getContext("2d")!;
  ctx.clearRect(0, 0, cw, ch);

  const pw = Math.round(lxMm * scale);   // panel pixel width
  const ph = Math.round(lyMm * scale);   // panel pixel height
  const ox = Math.round((cw - pw) / 2);  // panel offset x
  const oy = Math.round((ch - ph) / 2);

  const { grid, grid_n, region, outline, optimal_x, optimal_y } = lastResult;

  // ── Heat map ───────────────────────────────────────────────────────────────
  const cellW = pw / grid_n;
  const cellH = ph / grid_n;

  // The panel's real outline (and holes), used to clip and to draw the border.
  const outlinePath = new Path2D();
  for (const ring of outline) {
    ring.forEach(([x, y], i) => {
      const px = ox + x * pw, py = oy + y * ph;
      if (i === 0) outlinePath.moveTo(px, py); else outlinePath.lineTo(px, py);
    });
    outlinePath.closePath();
  }

  // Cells straddling a curved edge have their centre outside the panel;
  // borrow a neighbour's value so the clipped edge has no gaps.
  const display = fillOutside(grid, grid_n);
  ctx.save();
  ctx.clip(outlinePath, "evenodd");
  for (let row = 0; row < grid_n; row++) {
    for (let col = 0; col < grid_n; col++) {
      const i = row * grid_n + col;
      const v = display[i];
      if (v === null) continue;
      const [r, g, b] = sampleColormap(v);
      ctx.fillStyle = `rgb(${r},${g},${b})`;
      const x = ox + col * cellW, y = oy + row * cellH;
      ctx.fillRect(x, y, Math.ceil(cellW), Math.ceil(cellH));
      // Edge margin: excluded from the search, and the colour scale is set
      // by the interior, so these cells are clipped. Dim them.
      if (region[i] !== 0) {
        ctx.fillStyle = "rgba(0,0,0,0.45)";
        ctx.fillRect(x, y, Math.ceil(cellW), Math.ceil(cellH));
      }
    }
  }
  ctx.restore();

  // ── Panel border ──────────────────────────────────────────────────────────
  ctx.strokeStyle = "rgba(255,255,255,0.5)";
  ctx.lineWidth = 1;
  ctx.stroke(outlinePath);

  // ── Mode node lines ───────────────────────────────────────────────────────
  if (selectedShape) {
    drawNodeLines(ctx, selectedShape, SHAPE_GRID_N, ox, oy, pw, ph);
  }

  // ── Optimal position crosshair ────────────────────────────────────────────
  const optPx = ox + optimal_x * pw;
  const optPy = oy + optimal_y * ph;
  const r = 8;

  ctx.strokeStyle = "#ffffff";
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  ctx.moveTo(optPx - r - 4, optPy);
  ctx.lineTo(optPx + r + 4, optPy);
  ctx.moveTo(optPx, optPy - r - 4);
  ctx.lineTo(optPx, optPy + r + 4);
  ctx.stroke();

  ctx.strokeStyle = "#ffffff";
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  ctx.arc(optPx, optPy, r, 0, Math.PI * 2);
  ctx.stroke();

  // Inner dot
  ctx.fillStyle = "#ffffff";
  ctx.beginPath();
  ctx.arc(optPx, optPy, 2.5, 0, Math.PI * 2);
  ctx.fill();

  // ── Clicked comparison position ───────────────────────────────────────────
  if (probe) {
    const px = ox + probe.x * pw;
    const py = oy + probe.y * ph;
    ctx.lineWidth = 4;
    ctx.strokeStyle = "rgba(0,0,0,0.6)";
    ctx.beginPath();
    ctx.arc(px, py, 7, 0, Math.PI * 2);
    ctx.stroke();
    ctx.lineWidth = 2;
    ctx.strokeStyle = PROBE_COLOR;
    ctx.stroke();
  }

  // ── Dimension labels ──────────────────────────────────────────────────────
  ctx.fillStyle = "rgba(180,180,180,0.7)";
  ctx.font = `10px ${getComputedStyle(document.documentElement).getPropertyValue("--mono").trim() || "monospace"}`;
  ctx.textAlign = "center";
  ctx.fillText(`${lxMm} mm`, ox + pw / 2, oy + ph + 16);
  ctx.save();
  ctx.translate(ox - 14, oy + ph / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textAlign = "center";
  ctx.fillText(`${lyMm} mm`, 0, 0);
  ctx.restore();

}

/** Copy of `grid` with outside (null) cells next to the panel filled from
 *  their neighbours, two cells deep. */
function fillOutside(grid: (number | null)[], n: number): (number | null)[] {
  let cur = grid.slice();
  for (let pass = 0; pass < 2; pass++) {
    const next = cur.slice();
    for (let i = 0; i < n * n; i++) {
      if (cur[i] !== null) continue;
      const row = Math.floor(i / n), col = i % n;
      let sum = 0, count = 0;
      for (const [dr, dc] of [[-1, 0], [1, 0], [0, -1], [0, 1]]) {
        const r = row + dr, c = col + dc;
        if (r < 0 || c < 0 || r >= n || c >= n) continue;
        const v = cur[r * n + c];
        if (v !== null) { sum += v; count++; }
      }
      if (count) next[i] = sum / count;
    }
    cur = next;
  }
  return cur;
}

// Traces the zero contour of a mode shape (its node lines) with marching
// squares. `shape` is n×n, row-major, sampled at cell centres; null cells
// (outside the panel) are skipped.
function drawNodeLines(
  ctx: CanvasRenderingContext2D,
  shape: (number | null)[], n: number,
  ox: number, oy: number,
  pw: number, ph: number,
) {
  ctx.save();
  ctx.strokeStyle = "rgba(255, 220, 60, 0.85)";
  ctx.lineWidth = 1.5;
  ctx.beginPath();

  const px = (col: number) => ox + ((col + 0.5) / n) * pw;
  const py = (row: number) => oy + ((row + 0.5) / n) * ph;
  const at = (row: number, col: number) => shape[row * n + col] as number;

  for (let row = 0; row < n - 1; row++) {
    for (let col = 0; col < n - 1; col++) {
      if ([shape[row * n + col], shape[row * n + col + 1], shape[(row + 1) * n + col], shape[(row + 1) * n + col + 1]]
        .some((v) => v === null)) continue;
      // Corners clockwise from top-left, and the edges between them.
      const corners: [number, number][] = [[row, col], [row, col + 1], [row + 1, col + 1], [row + 1, col]];
      const crossings: [number, number][] = [];
      for (let k = 0; k < 4; k++) {
        const [r0, c0] = corners[k];
        const [r1, c1] = corners[(k + 1) % 4];
        const v0 = at(r0, c0);
        const v1 = at(r1, c1);
        if ((v0 < 0) !== (v1 < 0)) {
          const t = v0 / (v0 - v1);
          crossings.push([px(c0 + t * (c1 - c0)), py(r0 + t * (r1 - r0))]);
        }
      }
      // Two crossings: one segment. Four (a saddle): pair them in edge order.
      for (let k = 0; k + 1 < crossings.length; k += 2) {
        ctx.moveTo(...crossings[k]);
        ctx.lineTo(...crossings[k + 1]);
      }
    }
  }

  ctx.stroke();
  ctx.restore();
}

// ── Calculation ───────────────────────────────────────────────────────────────

function getParams(): PanelParams {
  return {
    shape:     selectShape.value,
    corner_r:  Math.max(0, parseFloat(inputCornerR.value) || 0) / 1000,
    sides:     Math.min(64, Math.max(3, parseInt(inputSides.value) || 6)),
    cutouts:   cutouts.map((c) => c.kind === "hole"
      ? { ...c, x: c.x / 1000, y: c.y / 1000, d: c.d / 1000 }
      : { ...c, x: c.x / 1000, y: c.y / 1000, length: c.length / 1000, width: c.width / 1000 }),
    lx:        parseFloat(inputLx.value) / 1000,
    ly:        parseFloat(inputLy.value) / 1000,
    h:         parseFloat(inputH.value)  / 1000,
    ex:        parseFloat(inputEx.value) * 1e6,   // MPa → Pa
    ey:        parseFloat(inputEy.value) * 1e6,
    g:         parseFloat(inputG.value)  * 1e6,
    nu:        parseFloat(inputNu.value),
    rho:       parseFloat(inputRho.value),
    boundary:  selectBoundary.value,
    freq_max:  parseFloat(inputFreqMax.value),
    grid_n:    Math.min(100, Math.max(4, parseInt(inputGridN.value) || 60)),
    exciter_d: Math.max(0, parseFloat(inputExciterD.value) || 0) / 1000,
    eta:       parseFloat(inputEta.value),
    score:     selectScore.value,
  };
}

// With "Isotropic" ticked, E_y and G follow from E_x and ν.
function syncIsotropic() {
  const iso = checkIsotropic.checked;
  inputEy.disabled = iso;
  inputG.disabled = iso;
  if (iso) {
    const ex = parseFloat(inputEx.value);
    const nu = parseFloat(inputNu.value);
    if (!isNaN(ex) && !isNaN(nu)) {
      inputEy.value = String(ex);
      inputG.value = String(Math.round((ex / (2 * (1 + nu))) * 100) / 100);
    }
  }
}

function setStatus(state: "calculating" | "done" | "error", msg: string) {
  statusText.className = state;
  statusText.textContent = msg;
}

async function calculate() {
  const id = ++requestId;
  setStatus("calculating", "Calculating…");

  const params = getParams();

  if (
    isNaN(params.lx) || params.lx <= 0 ||
    isNaN(params.ly) || params.ly <= 0 ||
    isNaN(params.h)  || params.h  <= 0 ||
    isNaN(params.ex) || params.ex <= 0 ||
    isNaN(params.ey) || params.ey <= 0 ||
    isNaN(params.g)  || params.g  <= 0 ||
    isNaN(params.nu) || params.nu <  0 ||
    isNaN(params.eta)|| params.eta <= 0 ||
    isNaN(params.rho)|| params.rho <= 0
  ) {
    setStatus("error", "Invalid parameters");
    return;
  }

  try {
    const result: CalculationResult = await invoke("compute_heatmap", { params });
    if (id !== requestId) return;
    lastResult = result;
    lastParams = params;
    updateUI(result, params);
    await Promise.all([loadSelectedShape(), loadProbe()]);
    if (id !== requestId) return;
    render();
    updateCharts();
    const solver = result.solver === "fea" ? " · FEA" : "";
    if (result.truncated_above !== null) {
      setStatus("done", `Modes above ${Math.round(result.truncated_above)} Hz omitted (solver limit)${solver}`);
    } else {
      setStatus("done", `Ready${solver}`);
    }
  } catch (err) {
    if (id !== requestId) return;
    setStatus("error", `Error: ${err}`);
  }
}

// Fetches the shape of the selected mode for the node-line overlay.
async function loadSelectedShape() {
  if (selectedModeIdx < 0 || !lastParams) {
    selectedShape = null;
    return;
  }
  try {
    selectedShape = await invoke("mode_shape", {
      params: lastParams,
      index: selectedModeIdx,
      n: SHAPE_GRID_N,
    });
  } catch {
    selectedShape = null;
  }
}

// Recomputes the clicked position's response for the current parameters.
async function loadProbe() {
  if (!probe || !lastParams) return;
  try {
    const curve: { db: number[]; raggedness: number } =
      await invoke("response_at", { params: lastParams, x: probe.x, y: probe.y });
    probe.db = curve.db;
    probe.raggedness = curve.raggedness;
  } catch {
    probe = null;
  }
}

function updateCharts() {
  if (!lastResult || !lastParams) return;
  const { bands, response_opt, raggedness_opt, modes } = lastResult;

  // Levels are relative: 0 dB is the mean of the optimal position's response.
  const ref = response_opt.length
    ? response_opt.reduce((a, b) => a + b, 0) / response_opt.length
    : 0;
  const series: Series[] = [{
    label: "Best position",
    color: "series1",
    values: response_opt.map((v) => v - ref),
  }];
  if (probe && probe.db.length === bands.length) {
    series.push({ label: "Clicked position", color: "series2", values: probe.db.map((v) => v - ref) });
  }
  responseChart.set(bands, series);

  const key = (color: string, text: string) =>
    `<span><span class="chart-swatch" style="background:${color}"></span>${text}</span>`;
  responseLegend.innerHTML =
    key("#3987e5", `Best ±${raggedness_opt.toFixed(1)} dB`) +
    (probe
      ? key(PROBE_COLOR, `Clicked ±${probe.raggedness.toFixed(1)} dB`) + '<button id="clear-probe">Clear</button>'
      : "<span>Click the heat map to compare a position</span>");
  document.getElementById("clear-probe")?.addEventListener("click", () => {
    probe = null;
    render();
    updateCharts();
  });

  const gaps = densityChart.set(modes.map((m) => m.freq), lastParams.freq_max);
  densityGaps.textContent = gaps === 0
    ? "No empty bands"
    : `${gaps} empty band${gaps === 1 ? "" : "s"}`;
}

function scheduleCalculate() {
  if (calcTimer) clearTimeout(calcTimer);
  calcTimer = setTimeout(calculate, 250);
}

function updateUI(result: CalculationResult, params: PanelParams) {
  // Optimal position
  const optXmm = (result.optimal_x * params.lx * 1000).toFixed(1);
  const optYmm = (result.optimal_y * params.ly * 1000).toFixed(1);
  const optXpct = (result.optimal_x * 100).toFixed(1);
  const optYpct = (result.optimal_y * 100).toFixed(1);
  valOptimal.textContent = `${optXmm} × ${optYmm} mm  (${optXpct}% × ${optYpct}%)`;
  valOptimal.classList.add("highlight");

  // Mode count
  valModes.textContent = String(result.mode_count);
  valRagged.textContent = result.bands.length ? `±${result.raggedness_opt.toFixed(1)} dB` : "—";
  badgeModes.textContent = `${result.mode_count} modes`;

  // Free-plate labels count node lines and are approximate; many modes mix
  // several patterns. FEA modes (non-rectangular shapes) have no label.
  const approx = params.boundary === "free" ? "≈" : "";
  const label = (m: ModeInfo) => (m.m === null ? "" : `${approx}(${m.m},${m.n})`);

  // Lowest frequency
  if (result.modes.length > 0) {
    valF1.textContent = `${result.modes[0].freq.toFixed(1)} Hz  ${label(result.modes[0])}`.trim();
  } else {
    valF1.textContent = "—";
  }

  // Populate mode selector
  const prev = selectedModeIdx;
  selectMode.innerHTML = '<option value="-1">None</option>';
  for (let i = 0; i < result.modes.length; i++) {
    const m = result.modes[i];
    const opt = document.createElement("option");
    opt.value = String(i);
    opt.textContent = `${i + 1}. ${label(m)}  ${m.freq.toFixed(0)} Hz`.replace(/\s+/g, " ");
    selectMode.appendChild(opt);
  }
  selectedModeIdx = prev < result.modes.length ? prev : -1;
  selectMode.value = String(selectedModeIdx);
}

// ── Canvas interactions ───────────────────────────────────────────────────────

function canvasToPanel(
  clientX: number, clientY: number,
): { x: number; y: number; normX: number; normY: number } | null {
  if (!lastResult) return null;

  const lxMm = parseFloat(inputLx.value);
  const lyMm = parseFloat(inputLy.value);
  const { cw, ch, scale } = getCanvasSize(lxMm, lyMm);

  const rect = canvas.getBoundingClientRect();
  const px = clientX - rect.left;
  const py = clientY - rect.top;

  const pw = Math.round(lxMm * scale);
  const ph = Math.round(lyMm * scale);
  const ox = Math.round((cw - pw) / 2);
  const oy = Math.round((ch - ph) / 2);

  const normX = (px - ox) / pw;
  const normY = (py - oy) / ph;

  if (normX < 0 || normX > 1 || normY < 0 || normY > 1) return null;

  return {
    x: normX * lxMm,
    y: normY * lyMm,
    normX,
    normY,
  };
}

canvas.addEventListener("mousemove", (e) => {
  const pos = canvasToPanel(e.clientX, e.clientY);
  if (!pos || !lastResult) {
    tooltip.classList.remove("visible");
    valCursor.textContent = "—";
    return;
  }

  const { normX, normY, x, y } = pos;
  const { grid, grid_n, region } = lastResult;
  const col = Math.min(Math.floor(normX * grid_n), grid_n - 1);
  const row = Math.min(Math.floor(normY * grid_n), grid_n - 1);
  const cell = row * grid_n + col;
  if (region[cell] === REGION_OUTSIDE) {
    tooltip.classList.remove("visible");
    valCursor.textContent = "—";
    return;
  }
  const normScore = grid[cell] ?? 0;

  valCursor.textContent = `${x.toFixed(1)} × ${y.toFixed(1)} mm`;

  const inMargin = region[cell] === REGION_MARGIN;
  const raw = lastResult.grid_raw[cell] ?? 0;
  const scoreText = inMargin
    ? "edge margin"
    : lastParams?.score === "coupling"
      ? `score ${(normScore * 100).toFixed(0)}%`
      : `raggedness ±${raw.toFixed(1)} dB`;
  tooltip.textContent = `${x.toFixed(1)} × ${y.toFixed(1)} mm  ·  ${scoreText}`;
  tooltip.classList.add("visible");

  const rect = canvas.getBoundingClientRect();
  let tx = e.clientX - rect.left + 12;
  let ty = e.clientY - rect.top  - 28;
  if (tx + 200 > rect.width) tx = e.clientX - rect.left - 180;
  tooltip.style.left = `${tx}px`;
  tooltip.style.top  = `${ty}px`;
});

canvas.addEventListener("click", async (e) => {
  const pos = canvasToPanel(e.clientX, e.clientY);
  if (!pos || !lastParams || !lastResult) return;
  const n = lastResult.grid_n;
  const cell = Math.min(Math.floor(pos.normY * n), n - 1) * n + Math.min(Math.floor(pos.normX * n), n - 1);
  if (lastResult.region[cell] === REGION_OUTSIDE) return;
  probe = { x: pos.normX, y: pos.normY, db: [], raggedness: 0 };
  const id = requestId;
  await loadProbe();
  if (id !== requestId) return;
  render();
  updateCharts();
});

canvas.addEventListener("mouseleave", () => {
  tooltip.classList.remove("visible");
  valCursor.textContent = "—";
});

// ── Input listeners ───────────────────────────────────────────────────────────

const numericInputs = [
  inputCornerR, inputSides, inputLx, inputLy, inputH, inputEx, inputEy, inputG, inputRho, inputNu,
  inputFreqMax, inputExciterD, inputEta, inputGridN,
];
numericInputs.forEach((el) => el.addEventListener("input", scheduleCalculate));
[inputEx, inputNu].forEach((el) => el.addEventListener("input", syncIsotropic));

checkIsotropic.addEventListener("change", () => {
  syncIsotropic();
  scheduleCalculate();
});

selectBoundary.addEventListener("change", scheduleCalculate);

// Shape-specific fields.
function syncShapeFields() {
  const shape = selectShape.value;
  $("field-corner-r").hidden = shape !== "rounded_rectangle";
  $("field-sides").hidden = shape !== "polygon";
  $("shape-hint").hidden = shape === "rectangle";
}
// ── Cutouts ──────────────────────────────────────────────────────────────────

const cutoutList = $("cutout-list");

type CutoutField = { key: string; label: string; unit: string; step: number };
const CUTOUT_FIELDS: Record<Cutout["kind"], CutoutField[]> = {
  hole: [
    { key: "x", label: "X", unit: "mm", step: 1 },
    { key: "y", label: "Y", unit: "mm", step: 1 },
    { key: "d", label: "Diameter", unit: "mm", step: 1 },
  ],
  slot: [
    { key: "x", label: "X", unit: "mm", step: 1 },
    { key: "y", label: "Y", unit: "mm", step: 1 },
    { key: "length", label: "Length", unit: "mm", step: 1 },
    { key: "width", label: "Width", unit: "mm", step: 0.5 },
    { key: "angle", label: "Angle (clockwise)", unit: "°", step: 5 },
  ],
};

function renderCutouts() {
  cutoutList.innerHTML = "";
  cutouts.forEach((c, i) => {
    const card = document.createElement("div");
    card.className = "cutout-card";
    const header = document.createElement("header");
    header.textContent = `${c.kind === "hole" ? "Hole" : "Slot"} ${i + 1}`;
    const remove = document.createElement("button");
    remove.type = "button";
    remove.textContent = "×";
    remove.title = "Remove this cutout";
    remove.addEventListener("click", () => {
      cutouts.splice(i, 1);
      renderCutouts();
      scheduleCalculate();
    });
    header.appendChild(remove);
    card.appendChild(header);

    for (const f of CUTOUT_FIELDS[c.kind]) {
      const label = document.createElement("label");
      label.textContent = f.label;
      const wrap = document.createElement("div");
      wrap.className = "input-with-unit";
      const input = document.createElement("input");
      input.type = "number";
      input.step = String(f.step);
      input.value = String((c as Record<string, number | string>)[f.key]);
      input.addEventListener("input", () => {
        const v = parseFloat(input.value);
        if (!isNaN(v)) {
          (c as Record<string, number | string>)[f.key] = v;
          scheduleCalculate();
        }
      });
      const unit = document.createElement("span");
      unit.className = "unit";
      unit.textContent = f.unit;
      wrap.append(input, unit);
      label.appendChild(wrap);
      card.appendChild(label);
    }
    cutoutList.appendChild(card);
  });
}

function addCutout(kind: Cutout["kind"]) {
  const w = parseFloat(inputLx.value) || 300;
  const h = parseFloat(inputLy.value) || 200;
  // Start somewhere plausible; the user positions it from there.
  cutouts.push(kind === "hole"
    ? { kind, x: Math.round(w * 0.3), y: Math.round(h * 0.5), d: Math.round(Math.min(w, h) * 0.1) }
    : { kind, x: Math.round(w * 0.5), y: Math.round(h * 0.5), length: Math.round(h * 0.5), width: 8, angle: 90 });
  probe = null;
  renderCutouts();
  scheduleCalculate();
}

$("add-hole").addEventListener("click", () => addCutout("hole"));
$("add-slot").addEventListener("click", () => addCutout("slot"));

selectShape.addEventListener("change", () => {
  syncShapeFields();
  probe = null;  // a clicked position may not be on the new shape
  scheduleCalculate();
});
selectScore.addEventListener("change", scheduleCalculate);

selectMaterial.addEventListener("change", () => {
  const preset = MATERIALS[selectMaterial.value];
  if (preset) {
    inputEx.value  = String(preset.ex);
    inputRho.value = String(preset.rho);
    inputNu.value  = String(preset.nu);
    inputEta.value = String(preset.eta);
    checkIsotropic.checked = preset.ey === undefined;
    if (preset.ey !== undefined && preset.g !== undefined) {
      inputEy.value = String(preset.ey);
      inputG.value  = String(preset.g);
    }
    syncIsotropic();
    scheduleCalculate();
  }
});

selectMode.addEventListener("change", async () => {
  selectedModeIdx = parseInt(selectMode.value);
  const id = requestId;
  await loadSelectedShape();
  if (id === requestId) render();
});


// Resize: re-render without recalculating
const resizeObserver = new ResizeObserver(() => render());
resizeObserver.observe(canvasWrap);
const chartObserver = new ResizeObserver(() => {
  responseChart.draw();
  densityChart.draw();
});
chartObserver.observe($("chart-panel"));

// ── Boot ──────────────────────────────────────────────────────────────────────

window.addEventListener("DOMContentLoaded", () => {
  syncIsotropic();
  syncShapeFields();
  calculate();
});
