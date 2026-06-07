import { invoke } from "@tauri-apps/api/core";

// ── Types ────────────────────────────────────────────────────────────────────

interface PanelParams {
  lx: number;
  ly: number;
  h: number;
  e: number;
  rho: number;
  nu: number;
  boundary: string;
  freq_max: number;
  grid_n: number;
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
}

// ── Material presets ─────────────────────────────────────────────────────────

const MATERIALS: Record<string, { e: number; rho: number; nu: number }> = {
  xps:      { e: 20,    rho: 32,   nu: 0.35 },  // XPS foam, E in MPa
  eps:      { e: 5,     rho: 20,   nu: 0.10 },  // EPS foam
  balsa:    { e: 3700,  rho: 130,  nu: 0.30 },  // Balsa wood
  birch:    { e: 9700,  rho: 680,  nu: 0.30 },  // Birch plywood
  acrylic:  { e: 3200,  rho: 1190, nu: 0.37 },  // PMMA
  aluminum: { e: 69000, rho: 2700, nu: 0.33 },  // Aluminium
  carbon:   { e: 70000, rho: 1600, nu: 0.10 },  // CFRP (isotropic estimate)
};

// ── State ────────────────────────────────────────────────────────────────────

let lastResult: CalculationResult | null = null;
let selectedModeIdx = -1;
let calcTimer: ReturnType<typeof setTimeout> | null = null;

// ── DOM refs ─────────────────────────────────────────────────────────────────

const $ = (id: string) => document.getElementById(id)!;

const inputLx        = $("lx")         as HTMLInputElement;
const inputLy        = $("ly")         as HTMLInputElement;
const inputH         = $("h")          as HTMLInputElement;
const inputE         = $("e")          as HTMLInputElement;
const inputRho       = $("rho")        as HTMLInputElement;
const inputNu        = $("nu")         as HTMLInputElement;
const selectMaterial = $("material-preset") as HTMLSelectElement;
const selectBoundary = $("boundary")   as HTMLSelectElement;
const inputFreqMax   = $("freq-max")   as HTMLInputElement;
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

  const { grid, grid_n, optimal_x, optimal_y, modes } = lastResult;

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

  // ── Panel border ──────────────────────────────────────────────────────────
  ctx.strokeStyle = "rgba(255,255,255,0.5)";
  ctx.lineWidth = 1;
  ctx.strokeRect(ox, oy, pw, ph);

  // ── Mode node lines ───────────────────────────────────────────────────────
  if (selectedModeIdx >= 0 && selectedModeIdx < modes.length) {
    const mode = modes[selectedModeIdx];
    const boundary = selectBoundary.value;
    drawNodeLines(ctx, mode.m, mode.n, boundary, ox, oy, pw, ph);
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

function drawNodeLines(
  ctx: CanvasRenderingContext2D,
  m: number, n: number,
  boundary: string,
  ox: number, oy: number,
  pw: number, ph: number,
) {
  ctx.save();
  ctx.strokeStyle = "rgba(255, 220, 60, 0.85)";
  ctx.lineWidth = 1;
  ctx.setLineDash([4, 3]);

  const xLines: number[] = [];
  const yLines: number[] = [];

  if (boundary === "simply_supported") {
    for (let k = 1; k < m; k++) xLines.push(k / m);
    for (let k = 1; k < n; k++) yLines.push(k / n);
  } else {
    // Free: zeros of cos(mπx/L) at x=(2k+1)/(2m) for k=0..m-1
    for (let k = 0; k < m; k++) xLines.push((2 * k + 1) / (2 * m));
    for (let k = 0; k < n; k++) yLines.push((2 * k + 1) / (2 * n));
  }

  for (const t of xLines) {
    const px = ox + t * pw;
    ctx.beginPath();
    ctx.moveTo(px, oy);
    ctx.lineTo(px, oy + ph);
    ctx.stroke();
  }
  for (const t of yLines) {
    const py = oy + t * ph;
    ctx.beginPath();
    ctx.moveTo(ox, py);
    ctx.lineTo(ox + pw, py);
    ctx.stroke();
  }

  ctx.restore();
}

// ── Calculation ───────────────────────────────────────────────────────────────

function getParams(): PanelParams {
  return {
    lx:       parseFloat(inputLx.value) / 1000,
    ly:       parseFloat(inputLy.value) / 1000,
    h:        parseFloat(inputH.value)  / 1000,
    e:        parseFloat(inputE.value)  * 1e6,   // MPa → Pa
    rho:      parseFloat(inputRho.value),
    nu:       parseFloat(inputNu.value),
    boundary: selectBoundary.value,
    freq_max: parseFloat(inputFreqMax.value),
    grid_n:   Math.min(100, Math.max(4, parseInt(inputGridN.value) || 60)),
  };
}

function setStatus(state: "calculating" | "done" | "error", msg: string) {
  statusText.className = state;
  statusText.textContent = msg;
}

async function calculate() {
  setStatus("calculating", "Calculating…");

  const params = getParams();

  if (
    isNaN(params.lx) || params.lx <= 0 ||
    isNaN(params.ly) || params.ly <= 0 ||
    isNaN(params.h)  || params.h  <= 0 ||
    isNaN(params.e)  || params.e  <= 0 ||
    isNaN(params.rho)|| params.rho <= 0
  ) {
    setStatus("error", "Invalid parameters");
    return;
  }

  try {
    const result: CalculationResult = await invoke("compute_heatmap", { params });
    lastResult = result;
    updateUI(result, params);
    render();
    setStatus("done", "Ready");
  } catch (err) {
    setStatus("error", `Error: ${err}`);
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

  // Lowest frequency
  if (result.modes.length > 0) {
    valF1.textContent = `${result.modes[0].freq.toFixed(1)} Hz  (${result.modes[0].m},${result.modes[0].n})`;
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
    opt.textContent = `(${m.m},${m.n})  ${m.freq.toFixed(0)} Hz`;
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

  tooltip.textContent = `${x.toFixed(1)} × ${y.toFixed(1)} mm  ·  score ${(normScore * 100).toFixed(0)}%`;
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

const numericInputs = [inputLx, inputLy, inputH, inputE, inputRho, inputNu, inputFreqMax, inputGridN];
numericInputs.forEach((el) => el.addEventListener("input", scheduleCalculate));

selectBoundary.addEventListener("change", scheduleCalculate);

selectMaterial.addEventListener("change", () => {
  const preset = MATERIALS[selectMaterial.value];
  if (preset) {
    inputE.value   = String(preset.e);
    inputRho.value = String(preset.rho);
    inputNu.value  = String(preset.nu);
    scheduleCalculate();
  }
});

selectMode.addEventListener("change", () => {
  selectedModeIdx = parseInt(selectMode.value);
  render();
});


// Resize: re-render without recalculating
const resizeObserver = new ResizeObserver(() => render());
resizeObserver.observe(canvasWrap);

// ── Boot ──────────────────────────────────────────────────────────────────────

window.addEventListener("DOMContentLoaded", () => {
  calculate();
});
