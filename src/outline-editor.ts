// Bézier outline editor drawn over the panel view.
//
// The outline is a closed loop of anchors. Each anchor may have an incoming
// and an outgoing handle; a segment is a cubic Bézier if either of its end
// anchors has a handle on that side, otherwise a straight line. Coordinates
// are normalised to the panel's width × height box (0..1), so the outline
// scales with the box and matches what the solver receives.

export type Pt = [number, number];

export interface Anchor {
  p: Pt;
  in: Pt | null;
  out: Pt | null;
  smooth: boolean;  // handles stay collinear when one is dragged
}

/** The panel's pixel frame on the canvas. */
export interface Frame { ox: number; oy: number; pw: number; ph: number }

type Segment =
  | { kind: "line"; to: Pt }
  | { kind: "cubic"; c1: Pt; c2: Pt; to: Pt };

// Bézier control distance for a quarter circle of unit radius.
const KAPPA = 0.5522847498307934;

const HIT_RADIUS = 7;  // px
const UNDO_LIMIT = 100;

// ── Seeding from the built-in shapes (mirrors geometry.rs) ───────────────────

const corner = (x: number, y: number): Anchor => ({ p: [x, y], in: null, out: null, smooth: false });

/** Anchors reproducing a built-in shape in the normalised box. `w` and `h`
 *  are the box size in mm, needed to keep round corners round. */
export function anchorsFromShape(shape: string, w: number, h: number, cornerR: number, sides: number): Anchor[] {
  if (shape === "ellipse") {
    const k = 0.5 * KAPPA;
    return [
      { p: [1, 0.5], in: [1, 0.5 - k], out: [1, 0.5 + k], smooth: true },
      { p: [0.5, 1], in: [0.5 + k, 1], out: [0.5 - k, 1], smooth: true },
      { p: [0, 0.5], in: [0, 0.5 + k], out: [0, 0.5 - k], smooth: true },
      { p: [0.5, 0], in: [0.5 - k, 0], out: [0.5 + k, 0], smooth: true },
    ];
  }
  if (shape === "rounded_rectangle" && cornerR > 0) {
    const r = Math.min(cornerR, w / 2, h / 2);
    const rx = r / w, ry = r / h, kx = rx * KAPPA, ky = ry * KAPPA;
    return [
      { p: [rx, 0], in: [rx - kx, 0], out: null, smooth: true },
      { p: [1 - rx, 0], in: null, out: [1 - rx + kx, 0], smooth: true },
      { p: [1, ry], in: [1, ry - ky], out: null, smooth: true },
      { p: [1, 1 - ry], in: null, out: [1, 1 - ry + ky], smooth: true },
      { p: [1 - rx, 1], in: [1 - rx + kx, 1], out: null, smooth: true },
      { p: [rx, 1], in: null, out: [rx - kx, 1], smooth: true },
      { p: [0, 1 - ry], in: [0, 1 - ry + ky], out: null, smooth: true },
      { p: [0, ry], in: null, out: [0, ry - ky], smooth: true },
    ];
  }
  if (shape === "polygon") {
    const n = Math.max(3, sides);
    const raw: Pt[] = Array.from({ length: n }, (_, i) => {
      const a = Math.PI / 2 + Math.PI / n + (2 * Math.PI * i) / n;
      return [Math.cos(a), Math.sin(a)];
    });
    const xs = raw.map((p) => p[0]), ys = raw.map((p) => p[1]);
    const [x0, x1, y0, y1] = [Math.min(...xs), Math.max(...xs), Math.min(...ys), Math.max(...ys)];
    return raw.map(([x, y]) => corner((x - x0) / (x1 - x0), (y - y0) / (y1 - y0)));
  }
  return [corner(0, 0), corner(1, 0), corner(1, 1), corner(0, 1)];
}

/** Segments of the closed loop, in normalised coordinates. */
function segments(anchors: Anchor[]): { from: Pt; seg: Segment }[] {
  return anchors.map((a, i) => {
    const b = anchors[(i + 1) % anchors.length];
    const seg: Segment = a.out || b.in
      ? { kind: "cubic", c1: a.out ?? a.p, c2: b.in ?? b.p, to: b.p }
      : { kind: "line", to: b.p };
    return { from: a.p, seg };
  });
}

/** The outline as the solver's path format (normalised). */
export function toPath(anchors: Anchor[]): { start: Pt; segments: Segment[] } {
  return { start: anchors[0].p, segments: segments(anchors).map((s) => s.seg) };
}

function cubicAt(p0: Pt, p1: Pt, p2: Pt, p3: Pt, t: number): Pt {
  const u = 1 - t;
  const a = u * u * u, b = 3 * u * u * t, c = 3 * u * t * t, d = t * t * t;
  return [a * p0[0] + b * p1[0] + c * p2[0] + d * p3[0], a * p0[1] + b * p1[1] + c * p2[1] + d * p3[1]];
}

const lerp = (a: Pt, b: Pt, t: number): Pt => [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t];
const clamp01 = (p: Pt): Pt => [Math.min(1, Math.max(0, p[0])), Math.min(1, Math.max(0, p[1]))];

// ── Editor ───────────────────────────────────────────────────────────────────

type Drag =
  | { kind: "anchor"; index: number; start: Pt; orig: Anchor }
  | { kind: "handle"; index: number; side: "in" | "out" };

export class OutlineEditor {
  anchors: Anchor[] = [];
  active = false;
  private selected: number | null = null;
  private drag: Drag | null = null;
  private undo: string[] = [];

  constructor(
    private canvas: HTMLCanvasElement,
    private frame: () => Frame | null,
    private onChange: () => void,
    private onRedraw: () => void,
  ) {
    canvas.addEventListener("mousedown", (e) => this.onDown(e));
    window.addEventListener("mousemove", (e) => this.onMove(e));
    window.addEventListener("mouseup", () => this.onUp());
    canvas.addEventListener("dblclick", (e) => this.onDoubleClick(e));
    window.addEventListener("keydown", (e) => this.onKey(e));
  }

  setAnchors(anchors: Anchor[]) {
    this.anchors = anchors;
    this.selected = null;
    this.undo = [];
  }

  setActive(active: boolean) {
    this.active = active;
    this.selected = null;
    this.drag = null;
    this.canvas.style.cursor = active ? "default" : "";
    this.onRedraw();
  }

  // ── Drawing ────────────────────────────────────────────────────────────────

  draw(ctx: CanvasRenderingContext2D, f: Frame) {
    if (this.anchors.length < 2) return;
    const px = (p: Pt): Pt => [f.ox + p[0] * f.pw, f.oy + p[1] * f.ph];

    // The outline itself.
    ctx.save();
    ctx.beginPath();
    const [sx, sy] = px(this.anchors[0].p);
    ctx.moveTo(sx, sy);
    for (const { seg } of segments(this.anchors)) {
      if (seg.kind === "line") ctx.lineTo(...px(seg.to));
      else ctx.bezierCurveTo(...px(seg.c1), ...px(seg.c2), ...px(seg.to));
    }
    ctx.closePath();
    ctx.lineWidth = 2;
    ctx.strokeStyle = "#ffffff";
    ctx.stroke();

    if (!this.active) {
      ctx.restore();
      return;
    }

    // Handles, then anchors on top.
    ctx.lineWidth = 1;
    for (const a of this.anchors) {
      for (const h of [a.in, a.out]) {
        if (!h) continue;
        const [ax, ay] = px(a.p), [hx, hy] = px(h);
        ctx.strokeStyle = "rgba(255,255,255,0.6)";
        ctx.beginPath();
        ctx.moveTo(ax, ay);
        ctx.lineTo(hx, hy);
        ctx.stroke();
        ctx.beginPath();
        ctx.arc(hx, hy, 4, 0, Math.PI * 2);
        ctx.fillStyle = "#1e1e1e";
        ctx.fill();
        ctx.strokeStyle = "#ffffff";
        ctx.stroke();
      }
    }
    this.anchors.forEach((a, i) => {
      const [x, y] = px(a.p);
      const s = i === this.selected ? 9 : 7;
      ctx.fillStyle = i === this.selected ? "#0a84ff" : "#ffffff";
      ctx.strokeStyle = "#1e1e1e";
      ctx.lineWidth = 2;
      if (a.smooth) {
        ctx.beginPath();
        ctx.arc(x, y, s / 2 + 0.5, 0, Math.PI * 2);
        ctx.fill();
        ctx.stroke();
      } else {
        ctx.fillRect(x - s / 2, y - s / 2, s, s);
        ctx.strokeRect(x - s / 2, y - s / 2, s, s);
      }
    });
    ctx.restore();
  }

  // ── Interaction ────────────────────────────────────────────────────────────

  private toNorm(e: MouseEvent, f: Frame): Pt {
    const rect = this.canvas.getBoundingClientRect();
    return [(e.clientX - rect.left - f.ox) / f.pw, (e.clientY - rect.top - f.oy) / f.ph];
  }

  private hit(e: MouseEvent, f: Frame): Drag | { kind: "none" } {
    const rect = this.canvas.getBoundingClientRect();
    const mx = e.clientX - rect.left, my = e.clientY - rect.top;
    const near = (p: Pt) => Math.hypot(f.ox + p[0] * f.pw - mx, f.oy + p[1] * f.ph - my) <= HIT_RADIUS;
    // Handles first: they sit close to their anchors.
    for (let i = 0; i < this.anchors.length; i++) {
      const a = this.anchors[i];
      if (a.out && near(a.out)) return { kind: "handle", index: i, side: "out" };
      if (a.in && near(a.in)) return { kind: "handle", index: i, side: "in" };
    }
    for (let i = 0; i < this.anchors.length; i++) {
      if (near(this.anchors[i].p)) {
        return { kind: "anchor", index: i, start: this.toNorm(e, f), orig: structuredClone(this.anchors[i]) };
      }
    }
    return { kind: "none" };
  }

  private snapshot() {
    this.undo.push(JSON.stringify(this.anchors));
    if (this.undo.length > UNDO_LIMIT) this.undo.shift();
  }

  private onDown(e: MouseEvent) {
    const f = this.frame();
    if (!this.active || !f || e.button !== 0) return;
    const hit = this.hit(e, f);
    if (hit.kind === "none") {
      this.selected = null;
      this.onRedraw();
      return;
    }
    e.preventDefault();
    if (hit.kind === "anchor" && e.altKey) {
      this.snapshot();
      this.toggleSmooth(hit.index);
      this.selected = hit.index;
      this.changed();
      return;
    }
    this.snapshot();
    this.drag = hit;
    if (hit.kind === "anchor") this.selected = hit.index;
    this.onRedraw();
  }

  private onMove(e: MouseEvent) {
    const f = this.frame();
    if (!this.drag || !f) return;
    const p = this.toNorm(e, f);
    if (this.drag.kind === "anchor") {
      const { index, start, orig } = this.drag;
      // Clamp the move so the anchor and both handles stay in the box (a
      // Bézier stays inside its control points' hull, so the curve does too).
      let dx = p[0] - start[0], dy = p[1] - start[1];
      for (const q of [orig.p, orig.in, orig.out]) {
        if (!q) continue;
        dx = Math.min(Math.max(dx, -q[0]), 1 - q[0]);
        dy = Math.min(Math.max(dy, -q[1]), 1 - q[1]);
      }
      const move = (q: Pt | null): Pt | null => (q ? [q[0] + dx, q[1] + dy] : null);
      this.anchors[index] = { ...orig, p: move(orig.p)!, in: move(orig.in), out: move(orig.out) };
    } else {
      const { index, side } = this.drag;
      const a = this.anchors[index];
      const h = clamp01(p);
      a[side] = h;
      // Smooth anchors keep the opposite handle collinear, at its own length.
      const other = side === "in" ? "out" : "in";
      const o = a[other];
      if (a.smooth && o) {
        const len = Math.hypot(o[0] - a.p[0], o[1] - a.p[1]);
        const dx = a.p[0] - h[0], dy = a.p[1] - h[1];
        const d = Math.hypot(dx, dy);
        if (d > 1e-9) a[other] = clamp01([a.p[0] + (dx / d) * len, a.p[1] + (dy / d) * len]);
      }
    }
    this.changed();
  }

  private onUp() {
    this.drag = null;
  }

  /** Double-click on the outline inserts an anchor without changing its shape. */
  private onDoubleClick(e: MouseEvent) {
    const f = this.frame();
    if (!this.active || !f) return;
    const rect = this.canvas.getBoundingClientRect();
    const mx = e.clientX - rect.left, my = e.clientY - rect.top;
    let best = { d: Infinity, i: -1, t: 0 };
    segments(this.anchors).forEach(({ from, seg }, i) => {
      for (let k = 1; k < 64; k++) {
        const t = k / 64;
        const q = seg.kind === "line" ? lerp(from, seg.to, t) : cubicAt(from, seg.c1, seg.c2, seg.to, t);
        const d = Math.hypot(f.ox + q[0] * f.pw - mx, f.oy + q[1] * f.ph - my);
        if (d < best.d) best = { d, i, t };
      }
    });
    if (best.d > HIT_RADIUS * 1.5) return;
    this.snapshot();
    const i = best.i, j = (i + 1) % this.anchors.length;
    const a = this.anchors[i], b = this.anchors[j];
    let anchor: Anchor;
    if (!a.out && !b.in) {
      anchor = corner(...lerp(a.p, b.p, best.t));
    } else {
      // de Casteljau split at t.
      const p0 = a.p, p1 = a.out ?? a.p, p2 = b.in ?? b.p, p3 = b.p, t = best.t;
      const q0 = lerp(p0, p1, t), q1 = lerp(p1, p2, t), q2 = lerp(p2, p3, t);
      const r0 = lerp(q0, q1, t), r1 = lerp(q1, q2, t);
      const s = lerp(r0, r1, t);
      a.out = q0;
      b.in = q2;
      anchor = { p: s, in: r0, out: r1, smooth: true };
    }
    this.anchors.splice(i + 1, 0, anchor);
    this.selected = i + 1;
    this.changed();
  }

  private onKey(e: KeyboardEvent) {
    if (!this.active) return;
    const target = e.target as HTMLElement;
    if (target && ["INPUT", "SELECT", "TEXTAREA"].includes(target.tagName)) return;
    if ((e.key === "Delete" || e.key === "Backspace") && this.selected !== null) {
      e.preventDefault();
      if (this.anchors.length <= 3) return;
      this.snapshot();
      this.anchors.splice(this.selected, 1);
      this.selected = null;
      this.changed();
    } else if (e.key === "z" && (e.metaKey || e.ctrlKey) && !e.shiftKey) {
      const prev = this.undo.pop();
      if (prev) {
        e.preventDefault();
        this.anchors = JSON.parse(prev);
        this.selected = null;
        this.changed();
      }
    } else if (e.key === "Escape") {
      this.selected = null;
      this.onRedraw();
    }
  }

  /** Corner ↔ smooth. Smooth handles follow the neighbours' direction, a
   *  third of the way to each. */
  private toggleSmooth(i: number) {
    const a = this.anchors[i];
    if (a.in || a.out) {
      this.anchors[i] = { ...a, in: null, out: null, smooth: false };
      return;
    }
    const n = this.anchors.length;
    const prev = this.anchors[(i - 1 + n) % n].p, next = this.anchors[(i + 1) % n].p;
    const tx = next[0] - prev[0], ty = next[1] - prev[1];
    const tl = Math.hypot(tx, ty) || 1;
    const lin = Math.hypot(a.p[0] - prev[0], a.p[1] - prev[1]) / 3;
    const lout = Math.hypot(next[0] - a.p[0], next[1] - a.p[1]) / 3;
    this.anchors[i] = {
      ...a,
      in: clamp01([a.p[0] - (tx / tl) * lin, a.p[1] - (ty / tl) * lin]),
      out: clamp01([a.p[0] + (tx / tl) * lout, a.p[1] + (ty / tl) * lout]),
      smooth: true,
    };
  }

  private changed() {
    this.onChange();
    this.onRedraw();
  }
}
