import { invoke } from "@tauri-apps/api/core";

// ── Types ────────────────────────────────────────────────────────────────────

interface PanelParams {
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
}

interface ModeInfo {
  m: number;
  n: number;
  freq: number;
}

interface CalculationResult {
  grid: number[];
  grid_n: number;
  modes: ModeInfo[];
  mode_count: number;
  optimal_x: number;
  optimal_y: number;
  optimal_score_raw: number;
  margin_x: number;
  margin_y: number;
  truncated_above: number | null;
}

// ── Material presets ─────────────────────────────────────────────────────────
// Moduli in MPa. Isotropic presets derive E_y and G from E_x and ν.
// Wood values are typical; real sheets vary widely, so measure if you can.

interface Material {
  ex: number;
  rho: number;
  nu: number;
  ey?: number;  // set for orthotropic materials
  g?: number;
}

const MATERIALS: Record<string, Material> = {
  xps:      { ex: 20,    rho: 32,   nu: 0.35 },                     // XPS foam
  eps:      { ex: 5,     rho: 20,   nu: 0.10 },                     // EPS foam
  balsa:    { ex: 3000,  rho: 130,  nu: 0.30, ey: 90,   g: 120 },   // Balsa, grain along x
  birch:    { ex: 10000, rho: 680,  nu: 0.07, ey: 5500, g: 620 },   // Birch plywood, face grain along x
  acrylic:  { ex: 3200,  rho: 1190, nu: 0.37 },                     // PMMA
  aluminum: { ex: 69000, rho: 2700, nu: 0.33 },                     // Aluminium
  carbon:   { ex: 70000, rho: 1600, nu: 0.10 },                     // CFRP, quasi-isotropic estimate
};

// Resolution of the mode shape used to trace node lines.
const SHAPE_GRID_N = 120;

// ── State ────────────────────────────────────────────────────────────────────

let lastResult: CalculationResult | null = null;
let lastParams: PanelParams | null = null;
let selectedModeIdx = -1;
let selectedShape: number[] | null = null;
let calcTimer: ReturnType<typeof setTimeout> | null = null;
let requestId = 0;  // newer calculations supersede older in-flight ones

// ── DOM refs ─────────────────────────────────────────────────────────────────

const $ = (id: string) => document.getElementById(id)!;

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

  const { grid, grid_n, optimal_x, optimal_y, margin_x, margin_y } = lastResult;

  // ── Heat map ───────────────────────────────────────────────────────────────
  const cellW = pw / grid_n;
  const cellH = ph / grid_n;

  for (let row = 0; row < grid_n; row++) {
    for (let col = 0; col < grid_n; col++) {
      const v = grid[row * grid_n + col];
      const [r, g, b] = sampleColormap(v);
      ctx.fillStyle = `rgb(${r},${g},${b})`;
      ctx.fillRect(
        ox + col * cellW,
        oy + row * cellH,
        Math.ceil(cellW),
        Math.ceil(cellH),
      );
    }
  }

  // ── Edge margin ───────────────────────────────────────────────────────────
  // Excluded from the search; the colour scale is set by the interior, so
  // these cells are clipped. Dim them to show they aren't candidates.
  const mx = margin_x * pw;
  const my = margin_y * ph;
  ctx.fillStyle = "rgba(0,0,0,0.45)";
  ctx.beginPath();
  ctx.rect(ox, oy, pw, ph);
  ctx.rect(ox + mx, oy + ph - my, pw - 2 * mx, -(ph - 2 * my));  // reverse winding cuts a hole
  ctx.fill("evenodd");
  ctx.strokeStyle = "rgba(255,255,255,0.25)";
  ctx.setLineDash([3, 3]);
  ctx.strokeRect(ox + mx, oy + my, pw - 2 * mx, ph - 2 * my);
  ctx.setLineDash([]);

  // ── Panel border ──────────────────────────────────────────────────────────
  ctx.strokeStyle = "rgba(255,255,255,0.5)";
  ctx.lineWidth = 1;
  ctx.strokeRect(ox, oy, pw, ph);

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

// Traces the zero contour of a mode shape (its node lines) with marching
// squares. `shape` is n×n, row-major, sampled at cell centres.
function drawNodeLines(
  ctx: CanvasRenderingContext2D,
  shape: number[], n: number,
  ox: number, oy: number,
  pw: number, ph: number,
) {
  ctx.save();
  ctx.strokeStyle = "rgba(255, 220, 60, 0.85)";
  ctx.lineWidth = 1.5;
  ctx.beginPath();

  const px = (col: number) => ox + ((col + 0.5) / n) * pw;
  const py = (row: number) => oy + ((row + 0.5) / n) * ph;
  const at = (row: number, col: number) => shape[row * n + col];

  for (let row = 0; row < n - 1; row++) {
    for (let col = 0; col < n - 1; col++) {
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
    await loadSelectedShape();
    if (id !== requestId) return;
    render();
    if (result.truncated_above !== null) {
      setStatus("done", `Modes above ${Math.round(result.truncated_above)} Hz omitted (solver limit)`);
    } else {
      setStatus("done", "Ready");
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
  badgeModes.textContent = `${result.mode_count} modes`;

  // Free-plate labels count node lines and are approximate; many modes mix
  // several patterns.
  const approx = params.boundary === "free" ? "≈" : "";
  const label = (m: ModeInfo) => `${approx}(${m.m},${m.n})`;

  // Lowest frequency
  if (result.modes.length > 0) {
    valF1.textContent = `${result.modes[0].freq.toFixed(1)} Hz  ${label(result.modes[0])}`;
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
    opt.textContent = `${i + 1}. ${label(m)}  ${m.freq.toFixed(0)} Hz`;
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
  const { grid, grid_n } = lastResult;
  const col = Math.min(Math.floor(normX * grid_n), grid_n - 1);
  const row = Math.min(Math.floor(normY * grid_n), grid_n - 1);
  const normScore = grid[row * grid_n + col];

  valCursor.textContent = `${x.toFixed(1)} × ${y.toFixed(1)} mm`;

  const inMargin =
    normX < lastResult.margin_x || normX > 1 - lastResult.margin_x ||
    normY < lastResult.margin_y || normY > 1 - lastResult.margin_y;
  const scoreText = inMargin ? "edge margin" : `score ${(normScore * 100).toFixed(0)}%`;
  tooltip.textContent = `${x.toFixed(1)} × ${y.toFixed(1)} mm  ·  ${scoreText}`;
  tooltip.classList.add("visible");

  const rect = canvas.getBoundingClientRect();
  let tx = e.clientX - rect.left + 12;
  let ty = e.clientY - rect.top  - 28;
  if (tx + 200 > rect.width) tx = e.clientX - rect.left - 180;
  tooltip.style.left = `${tx}px`;
  tooltip.style.top  = `${ty}px`;
});

canvas.addEventListener("mouseleave", () => {
  tooltip.classList.remove("visible");
  valCursor.textContent = "—";
});

// ── Input listeners ───────────────────────────────────────────────────────────

const numericInputs = [
  inputLx, inputLy, inputH, inputEx, inputEy, inputG, inputRho, inputNu,
  inputFreqMax, inputExciterD, inputGridN,
];
numericInputs.forEach((el) => el.addEventListener("input", scheduleCalculate));
[inputEx, inputNu].forEach((el) => el.addEventListener("input", syncIsotropic));

checkIsotropic.addEventListener("change", () => {
  syncIsotropic();
  scheduleCalculate();
});

selectBoundary.addEventListener("change", scheduleCalculate);

selectMaterial.addEventListener("change", () => {
  const preset = MATERIALS[selectMaterial.value];
  if (preset) {
    inputEx.value  = String(preset.ex);
    inputRho.value = String(preset.rho);
    inputNu.value  = String(preset.nu);
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

// ── Boot ──────────────────────────────────────────────────────────────────────

window.addEventListener("DOMContentLoaded", () => {
  syncIsotropic();
  calculate();
});
