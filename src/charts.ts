// Response and modal-density charts, drawn on <canvas> with no dependencies.
//
// Colours follow the app's dark surface. Series colours were checked for
// colour-vision separation and contrast against it; text always uses text
// colours, never series colours.

const COLORS = {
  surface: "#1a1a1a",
  grid: "#2e2e2e",
  axis: "#444",
  text: "#e8e8e8",
  muted: "#888",
  series1: "#3987e5",  // best position
  series2: "#d95926",  // clicked position
};

const FONT = `10px "SF Mono", ui-monospace, monospace`;

export interface Series {
  label: string;
  color: "series1" | "series2";
  values: number[];  // one per x value; NaN = missing
}

// ── Shared helpers ───────────────────────────────────────────────────────────

/** Sizes the canvas backing store to its CSS box at device pixel ratio. */
function setup(canvas: HTMLCanvasElement): { ctx: CanvasRenderingContext2D; w: number; h: number } {
  const rect = canvas.getBoundingClientRect();
  const dpr = window.devicePixelRatio || 1;
  canvas.width = Math.max(1, Math.round(rect.width * dpr));
  canvas.height = Math.max(1, Math.round(rect.height * dpr));
  const ctx = canvas.getContext("2d")!;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, rect.width, rect.height);
  ctx.font = FONT;
  return { ctx, w: rect.width, h: rect.height };
}

export function formatHz(f: number): string {
  if (f >= 1000) {
    const k = f / 1000;
    return `${k >= 10 ? k.toFixed(0) : k.toFixed(1).replace(/\.0$/, "")}k`;
  }
  return f >= 100 ? f.toFixed(0) : f.toFixed(0);
}

const FREQ_TICKS = [20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000];

function showTooltip(tip: HTMLElement, canvas: HTMLCanvasElement, x: number, y: number, html: string) {
  tip.innerHTML = html;
  tip.classList.add("visible");
  const box = canvas.getBoundingClientRect();
  const parent = tip.offsetParent?.getBoundingClientRect() ?? box;
  let left = box.left - parent.left + x + 12;
  const top = box.top - parent.top + Math.max(4, y - 40);
  if (left + tip.offsetWidth > parent.width - 4) left = box.left - parent.left + x - tip.offsetWidth - 12;
  tip.style.left = `${left}px`;
  tip.style.top = `${top}px`;
}

function swatch(color: string): string {
  return `<span class="chart-swatch" style="background:${color}"></span>`;
}

// ── Response chart ───────────────────────────────────────────────────────────

const PAD = { left: 40, right: 12, top: 10, bottom: 22 };

export class ResponseChart {
  private freqs: number[] = [];
  private series: Series[] = [];
  private hoverIdx = -1;

  constructor(private canvas: HTMLCanvasElement, private tip: HTMLElement) {
    canvas.addEventListener("mousemove", (e) => this.onHover(e));
    canvas.addEventListener("mouseleave", () => {
      this.hoverIdx = -1;
      this.tip.classList.remove("visible");
      this.draw();
    });
  }

  set(freqs: number[], series: Series[]) {
    this.freqs = freqs;
    this.series = series;
    this.hoverIdx = -1;
    this.draw();
  }

  private geometry(w: number, h: number) {
    const x0 = PAD.left, x1 = w - PAD.right, y0 = PAD.top, y1 = h - PAD.bottom;
    const fMin = this.freqs[0], fMax = this.freqs[this.freqs.length - 1];
    const lMin = Math.log(fMin), lSpan = Math.log(fMax) - lMin || 1;
    const xOf = (f: number) => x0 + ((Math.log(f) - lMin) / lSpan) * (x1 - x0);

    const all = this.series.flatMap((s) => s.values).filter((v) => isFinite(v));
    const step = 10;
    const yMax = Math.ceil(Math.max(...all) / step) * step;
    const yMin = Math.floor(Math.min(...all) / step) * step;
    const yOf = (v: number) => y1 - ((v - yMin) / (yMax - yMin || 1)) * (y1 - y0);
    return { x0, x1, y0, y1, xOf, yOf, yMin, yMax, step };
  }

  draw() {
    const { ctx, w, h } = setup(this.canvas);
    if (this.freqs.length < 2 || this.series.length === 0) {
      ctx.fillStyle = COLORS.muted;
      ctx.textAlign = "center";
      ctx.fillText("No modes in range", w / 2, h / 2);
      return;
    }
    const g = this.geometry(w, h);

    // Gridlines and y labels (hairline, solid, recessive)
    ctx.lineWidth = 1;
    ctx.textAlign = "right";
    ctx.textBaseline = "middle";
    for (let v = g.yMin; v <= g.yMax; v += g.step) {
      const y = Math.round(g.yOf(v)) + 0.5;
      ctx.strokeStyle = COLORS.grid;
      ctx.beginPath();
      ctx.moveTo(g.x0, y);
      ctx.lineTo(g.x1, y);
      ctx.stroke();
      ctx.fillStyle = COLORS.muted;
      ctx.fillText(`${v > 0 ? "+" : ""}${v}`, g.x0 - 6, y);
    }
    ctx.save();
    ctx.translate(10, (g.y0 + g.y1) / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = "center";
    ctx.fillText("dB", 0, 0);
    ctx.restore();

    // Frequency ticks
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    const fMin = this.freqs[0], fMax = this.freqs[this.freqs.length - 1];
    for (const f of FREQ_TICKS) {
      if (f < fMin || f > fMax) continue;
      const x = Math.round(g.xOf(f)) + 0.5;
      ctx.strokeStyle = COLORS.grid;
      ctx.beginPath();
      ctx.moveTo(x, g.y0);
      ctx.lineTo(x, g.y1);
      ctx.stroke();
      ctx.fillStyle = COLORS.muted;
      ctx.fillText(`${formatHz(f)}`, x, g.y1 + 6);
    }
    ctx.fillText("Hz", g.x1 - 6, g.y1 + 6);

    // Series: 2px lines, round joins. Later series drawn on top.
    for (const s of this.series) {
      ctx.strokeStyle = COLORS[s.color];
      ctx.lineWidth = 2;
      ctx.lineJoin = "round";
      ctx.lineCap = "round";
      ctx.beginPath();
      let started = false;
      s.values.forEach((v, i) => {
        if (!isFinite(v)) { started = false; return; }
        const x = g.xOf(this.freqs[i]), y = g.yOf(v);
        if (started) ctx.lineTo(x, y); else ctx.moveTo(x, y);
        started = true;
      });
      ctx.stroke();
    }

    // Hover crosshair with ringed dots
    if (this.hoverIdx >= 0) {
      const x = Math.round(g.xOf(this.freqs[this.hoverIdx])) + 0.5;
      ctx.strokeStyle = COLORS.axis;
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(x, g.y0);
      ctx.lineTo(x, g.y1);
      ctx.stroke();
      for (const s of this.series) {
        const v = s.values[this.hoverIdx];
        if (!isFinite(v)) continue;
        ctx.beginPath();
        ctx.arc(x, g.yOf(v), 4, 0, Math.PI * 2);
        ctx.fillStyle = COLORS[s.color];
        ctx.fill();
        ctx.lineWidth = 2;
        ctx.strokeStyle = COLORS.surface;
        ctx.stroke();
      }
    }
  }

  private onHover(e: MouseEvent) {
    if (this.freqs.length < 2) return;
    const rect = this.canvas.getBoundingClientRect();
    const mx = e.clientX - rect.left;
    const g = this.geometry(rect.width, rect.height);
    if (mx < g.x0 - 4 || mx > g.x1 + 4) {
      this.hoverIdx = -1;
      this.tip.classList.remove("visible");
      this.draw();
      return;
    }
    // Nearest band in log-frequency
    let best = 0, bestDist = Infinity;
    this.freqs.forEach((f, i) => {
      const d = Math.abs(g.xOf(f) - mx);
      if (d < bestDist) { bestDist = d; best = i; }
    });
    this.hoverIdx = best;
    this.draw();
    const rows = this.series
      .filter((s) => isFinite(s.values[best]))
      .map((s) => `${swatch(COLORS[s.color])}${s.label} <b>${s.values[best].toFixed(1)} dB</b>`)
      .join("<br>");
    showTooltip(this.tip, this.canvas, g.xOf(this.freqs[best]), e.clientY - rect.top,
      `<div class="chart-tip-title">${formatHz(this.freqs[best])} Hz</div>${rows}`);
  }
}

// ── Modal density chart ──────────────────────────────────────────────────────

// Nominal third-octave centres (IEC 61260), extended as needed.
function thirdOctaveBands(fLo: number, fHi: number): { centre: number; lo: number; hi: number }[] {
  const bands = [];
  for (let i = -20; i <= 14; i++) {
    const centre = 1000 * Math.pow(2, i / 3);
    const lo = centre / Math.pow(2, 1 / 6), hi = centre * Math.pow(2, 1 / 6);
    if (hi < fLo || lo > fHi) continue;
    bands.push({ centre, lo, hi });
  }
  return bands;
}

const NOMINAL = [16, 20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800,
  1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000];

function nominal(f: number): number {
  return NOMINAL.reduce((a, b) => (Math.abs(Math.log(b / f)) < Math.abs(Math.log(a / f)) ? b : a));
}

export class ModeDensityChart {
  private bands: { centre: number; lo: number; hi: number; count: number }[] = [];
  private hoverIdx = -1;

  constructor(private canvas: HTMLCanvasElement, private tip: HTMLElement) {
    canvas.addEventListener("mousemove", (e) => this.onHover(e));
    canvas.addEventListener("mouseleave", () => {
      this.hoverIdx = -1;
      this.tip.classList.remove("visible");
      this.draw();
    });
  }

  /** Returns how many bands above the first mode contain no modes. */
  set(modeFreqs: number[], fMax: number): number {
    if (modeFreqs.length === 0) {
      this.bands = [];
      this.draw();
      return 0;
    }
    this.bands = thirdOctaveBands(modeFreqs[0], fMax).map((b) => ({
      ...b,
      count: modeFreqs.filter((f) => f >= b.lo && f < b.hi).length,
    }));
    this.hoverIdx = -1;
    this.draw();
    return this.bands.filter((b) => b.count === 0).length;
  }

  private geometry(w: number, h: number) {
    const x0 = 30, x1 = w - 8, y0 = 10, y1 = h - PAD.bottom;
    const slot = (x1 - x0) / Math.max(1, this.bands.length);
    const barW = Math.min(24, Math.max(2, slot - 2));  // ≤24px, ≥2px gap
    const maxCount = Math.max(1, ...this.bands.map((b) => b.count));
    const step = maxCount <= 5 ? 1 : maxCount <= 10 ? 2 : maxCount <= 25 ? 5 : 10;
    const yMax = Math.ceil(maxCount / step) * step;
    const yOf = (v: number) => y1 - (v / yMax) * (y1 - y0);
    return { x0, x1, y0, y1, slot, barW, yOf, yMax, step };
  }

  draw() {
    const { ctx, w, h } = setup(this.canvas);
    if (this.bands.length === 0) {
      ctx.fillStyle = COLORS.muted;
      ctx.textAlign = "center";
      ctx.fillText("No modes in range", w / 2, h / 2);
      return;
    }
    const g = this.geometry(w, h);

    ctx.lineWidth = 1;
    ctx.textAlign = "right";
    ctx.textBaseline = "middle";
    for (let v = 0; v <= g.yMax; v += g.step) {
      const y = Math.round(g.yOf(v)) + 0.5;
      ctx.strokeStyle = v === 0 ? COLORS.axis : COLORS.grid;
      ctx.beginPath();
      ctx.moveTo(g.x0, y);
      ctx.lineTo(g.x1, y);
      ctx.stroke();
      ctx.fillStyle = COLORS.muted;
      ctx.fillText(String(v), g.x0 - 6, y);
    }

    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    this.bands.forEach((b, i) => {
      const cx = g.x0 + g.slot * (i + 0.5);
      if (b.count > 0) {
        // Bar: square at the baseline, 4px rounded data end.
        const top = g.yOf(b.count), x = cx - g.barW / 2;
        const r = Math.min(4, g.barW / 2, g.y1 - top);
        ctx.fillStyle = COLORS.series1;
        ctx.globalAlpha = this.hoverIdx === -1 || this.hoverIdx === i ? 1 : 0.55;
        ctx.beginPath();
        ctx.moveTo(x, g.y1);
        ctx.lineTo(x, top + r);
        ctx.arcTo(x, top, x + r, top, r);
        ctx.lineTo(x + g.barW - r, top);
        ctx.arcTo(x + g.barW, top, x + g.barW, top + r, r);
        ctx.lineTo(x + g.barW, g.y1);
        ctx.closePath();
        ctx.fill();
        ctx.globalAlpha = 1;
      } else {
        // Gap: an empty band. Marked in text ink so it never relies on colour.
        ctx.fillStyle = COLORS.muted;
        ctx.textBaseline = "bottom";
        ctx.fillText("0", cx, g.y1 - 2);
        ctx.textBaseline = "top";
      }
      // Label octave bands only, to avoid collisions.
      if (i % 3 === 0 || this.bands.length <= 8) {
        ctx.fillStyle = COLORS.muted;
        ctx.fillText(formatHz(nominal(b.centre)), cx, g.y1 + 6);
      }
    });
  }

  private onHover(e: MouseEvent) {
    if (this.bands.length === 0) return;
    const rect = this.canvas.getBoundingClientRect();
    const g = this.geometry(rect.width, rect.height);
    const mx = e.clientX - rect.left;
    const i = Math.floor((mx - g.x0) / g.slot);
    if (i < 0 || i >= this.bands.length) {
      this.hoverIdx = -1;
      this.tip.classList.remove("visible");
      this.draw();
      return;
    }
    this.hoverIdx = i;
    this.draw();
    const b = this.bands[i];
    const what = b.count === 0 ? "no modes (gap)" : `${b.count} mode${b.count === 1 ? "" : "s"}`;
    showTooltip(this.tip, this.canvas, g.x0 + g.slot * (i + 0.5), e.clientY - rect.top,
      `<div class="chart-tip-title">${formatHz(nominal(b.centre))} Hz ⅓-octave</div>${what}`);
  }
}
