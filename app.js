"use strict";
/* app.js — DOM/UI/charts for the Panel Data Multiverse Lab. Numeric core: engine.js (window.PanelEngine). */
var E = window.PanelEngine;

// ============== STATE ==============
var DATA = null;           // array of person-year records
var META = null;           // { industryCats, regionCats, years }
var S = {
  focal: "union", outcome: "lwage", seType: "cluster", viz: "data",
  ss: { estimator: "within", controls: { experience: true, hours: false, industry: false, region: false, health: false, race: false }, yearFE: false, sample: "full" },
  axes: { estimators: ["pooled", "between", "random", "within", "fd"], controlsVary: ["experience", "industry", "region", "health"], yearFEVary: true, sampleVary: ["full"] },
  multiResults: null, multiSummary: null, lessonStep: -1
};

var FOCAL_LABEL = { union: "union wage premium", married: "marriage wage premium", educ: "return to schooling" };
var FOCAL_SHORT = { union: "union", married: "married", educ: "educ" };
var FOCAL_HELP = {
  union: "Do unionized workers earn more? Watch how much of the gap is selection — higher-paid men sorting into union jobs.",
  married: "The “marriage premium”: married men out-earn single men. How much survives once we remove stable person traits?",
  educ: "Returns to a year of schooling. Note: education is fixed over 1980–87, so Fixed Effects and First Differences cannot identify it."
};
var ESTIMATOR_LABEL = { pooled: "Pooled OLS", between: "Between", random: "Random Effects", within: "Fixed Effects", fd: "First Diff.", cre: "Correlated RE", twfe: "Two-way FE" };
var CONTROL_LABEL = { experience: "experience (+²)", hours: "hours", industry: "industry", region: "region", health: "health", race: "race" };
var SAMPLE_LABEL = { full: "full sample", nohealth: "drop poor health", trim: "trim 1/99% outliers" };

// ============== FORMAT / SMALL HELPERS ==============
function fmt(x, d) { d = d == null ? 3 : d; return (x == null || Number.isNaN(x) || !isFinite(x)) ? "—" : Number(x).toFixed(d); }
function ciOf(b, se) { return [b - 1.96 * se, b + 1.96 * se]; }
function isSig(b, se) { if (!isFinite(se)) return false; var c = ciOf(b, se); return c[0] > 0 || c[1] < 0; }
function outcomeKey() { return S.outcome === "wage" ? "wage" : "lwage"; }
function outcomeLabel() { return S.outcome === "wage" ? "wage ($)" : "log wage"; }
function $(id) { return document.getElementById(id); }

// ============== SVG PRIMITIVES ==============
var SVGNS = "http://www.w3.org/2000/svg";
function svgEl(w, h) { var s = document.createElementNS(SVGNS, "svg"); s.setAttribute("viewBox", "0 0 " + w + " " + h); s.setAttribute("width", "100%"); s.setAttribute("height", "100%"); return s; }
function line(x1, y1, x2, y2, stroke, sw) { var e = document.createElementNS(SVGNS, "line"); e.setAttribute("x1", x1); e.setAttribute("y1", y1); e.setAttribute("x2", x2); e.setAttribute("y2", y2); e.setAttribute("stroke", stroke || "#6b7280"); e.setAttribute("stroke-width", sw || 1); return e; }
function rect(x, y, w, h, fill) { var e = document.createElementNS(SVGNS, "rect"); e.setAttribute("x", x); e.setAttribute("y", y); e.setAttribute("width", w); e.setAttribute("height", h); e.setAttribute("fill", fill || "#60a5fa"); return e; }
function circle(x, y, r, fill, op) { var e = document.createElementNS(SVGNS, "circle"); e.setAttribute("cx", x); e.setAttribute("cy", y); e.setAttribute("r", r); e.setAttribute("fill", fill || "#cbd5e1"); e.setAttribute("fill-opacity", op == null ? 0.85 : op); return e; }
function polyline(points, stroke, sw, op) { var e = document.createElementNS(SVGNS, "polyline"); e.setAttribute("points", points.map(function (p) { return p[0] + "," + p[1]; }).join(" ")); e.setAttribute("fill", "none"); e.setAttribute("stroke", stroke || "#93c5fd"); e.setAttribute("stroke-width", sw || 1); e.setAttribute("stroke-opacity", op == null ? 1 : op); return e; }
function txt(x, y, t, fill, fs, anchor) { var e = document.createElementNS(SVGNS, "text"); e.setAttribute("x", x); e.setAttribute("y", y); e.setAttribute("fill", fill || "#9aa0a6"); e.setAttribute("font-size", fs || 10); e.setAttribute("text-anchor", anchor || "middle"); e.textContent = t; return e; }
function scaleLin(d0, d1, r0, r1) { var d = (d1 - d0) || 1e-9, m = (r1 - r0) / d; return function (v) { return r0 + (v - d0) * m; }; }
function gridY(svg, m, w, h, n) { for (var i = 0; i <= n; i++) { var gy = m.t + i * (h - m.t - m.b) / n; svg.appendChild(line(m.l, gy, w - m.r, gy, "#1a1a1a", 1)); } }
var CHART_TITLES = {
  "chart-spaghetti": "Wage trajectories over time", "chart-dist": "Distribution of the outcome",
  "chart-bygroup": "Mean outcome by focal status and year", "chart-variance": "Within- vs between-person variation",
  "chart-speccurve": "Specification curve", "chart-estbars": "All estimators compared",
  "chart-rvf": "Residuals vs fitted", "chart-hist": "Residual distribution", "chart-qq": "Normal Q–Q plot"
};
function setChart(id, svg) {
  var c = $(id); c.innerHTML = ""; c.appendChild(svg);
  var b = document.createElement("button"); b.className = "chart-expand"; b.type = "button";
  b.setAttribute("aria-label", "Enlarge figure"); b.title = "Enlarge"; b.textContent = "⤢";
  b.addEventListener("click", function (e) { e.stopPropagation(); openFig(id); });
  c.appendChild(b);
}
function openFig(id) {
  var c = $(id); if (!c) return; var svg = c.querySelector("svg"); if (!svg) return;
  var clone = svg.cloneNode(true);
  $("figModalBody").innerHTML = ""; $("figModalBody").appendChild(clone);
  $("figModalTitle").textContent = CHART_TITLES[id] || "Figure";
  $("figModalCap").textContent = c.getAttribute("aria-label") || "";
  $("figModal").classList.add("open");
}
function closeFig() { $("figModal").classList.remove("open"); $("figModalBody").innerHTML = ""; }
function emptyChart(id, msg) { $(id).innerHTML = "<div style='padding:24px;color:#6b7280'>" + (msg || "No data") + "</div>"; }
function ariaChart(id, label) { var c = $(id); if (c) c.setAttribute("aria-label", label); }

// "nice" round tick values spanning [lo,hi]
function niceNum(x, round) {
  if (x <= 0 || !isFinite(x)) return 1;
  var e = Math.floor(Math.log(x) / Math.LN10), f = x / Math.pow(10, e), nf;
  if (round) nf = f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10;
  else nf = f <= 1 ? 1 : f <= 2 ? 2 : f <= 5 ? 5 : 10;
  return nf * Math.pow(10, e);
}
function niceTicks(lo, hi, n) {
  n = n || 5; if (!(hi > lo)) { hi = lo + 1; }
  var step = niceNum(niceNum(hi - lo, false) / (n - 1), true);
  var glo = Math.floor(lo / step) * step, ghi = Math.ceil(hi / step) * step, t = [];
  for (var v = glo; v <= ghi + step * 0.5; v += step) t.push(Math.abs(v) < step / 1e6 ? 0 : +v.toFixed(10));
  return { ticks: t, lo: glo, hi: ghi, step: step };
}
function tickFmt(step) { var d = step >= 1 ? 0 : Math.min(4, Math.ceil(-Math.log(step) / Math.LN10) + 1); return function (v) { return v.toFixed(d); }; }
// draw a full axis frame: gridlines, ticks, numeric labels, axis titles
function plotFrame(svg, m, w, h, xs, ys, xt, yt, xTitle, yTitle) {
  var x0 = m.l, x1 = w - m.r, yb = h - m.b, yt0 = m.t, xf = tickFmt(xt.step), yf = tickFmt(yt.step);
  yt.ticks.forEach(function (v) { var y = ys(v); if (y < yt0 - 0.5 || y > yb + 0.5) return; svg.appendChild(line(x0, y, x1, y, "#1b1b1b", 1)); svg.appendChild(line(x0 - 4, y, x0, y, "#5b6472", 1)); svg.appendChild(txt(x0 - 7, y + 3.5, yf(v), "#9aa0a6", 11, "end")); });
  xt.ticks.forEach(function (v) { var x = xs(v); if (x < x0 - 0.5 || x > x1 + 0.5) return; svg.appendChild(line(x, yt0, x, yb, "#1b1b1b", 1)); svg.appendChild(line(x, yb, x, yb + 4, "#5b6472", 1)); svg.appendChild(txt(x, yb + 16, xf(v), "#9aa0a6", 11, "middle")); });
  svg.appendChild(line(x0, yb, x1, yb, "#6b7280", 1.2));
  svg.appendChild(line(x0, yt0, x0, yb, "#6b7280", 1.2));
  svg.appendChild(txt((x0 + x1) / 2, h - 6, xTitle, "#cbd5e1", 12, "middle"));
  var yl = txt(13, (yt0 + yb) / 2, yTitle, "#cbd5e1", 12, "middle"); yl.setAttribute("transform", "rotate(-90 13 " + ((yt0 + yb) / 2) + ")"); svg.appendChild(yl);
}

// ============== META & GROUPING ==============
function computeMeta(rows) {
  function cats(key) { var c = {}; rows.forEach(function (r) { c[r[key]] = (c[r[key]] || 0) + 1; }); return Object.keys(c).sort(function (a, b) { return c[b] - c[a]; }); }
  var years = Array.from(new Set(rows.map(function (r) { return r.year; }))).sort(function (a, b) { return a - b; });
  return { industryCats: cats("industry"), regionCats: cats("region"), years: years };
}
function byPerson(rows) { var by = {}; rows.forEach(function (r) { (by[r.nr] || (by[r.nr] = [])).push(r); }); return by; }

// variance decomposition: share of variation that is within-person vs between-person
function varianceDecomp(rows, key) {
  var by = byPerson(rows), grand = E.mean(rows.map(function (r) { return r[key]; }));
  var ssb = 0, ssw = 0;
  Object.keys(by).forEach(function (id) {
    var g = by[id], mi = E.mean(g.map(function (r) { return r[key]; }));
    ssb += g.length * (mi - grand) * (mi - grand);
    g.forEach(function (r) { var d = r[key] - mi; ssw += d * d; });
  });
  var tot = ssb + ssw || 1;
  return { within: ssw / tot, between: ssb / tot };
}

// ============== DATA OVERVIEW CHARTS ==============
function drawSpaghetti() {
  var key = outcomeKey(), by = byPerson(DATA), ids = Object.keys(by);
  var pick = ids.slice().sort(function (a, b) { return (a * 7919) % 101 - (b * 7919) % 101; }).slice(0, 40);
  var w = 480, h = 300, m = { l: 48, r: 12, t: 12, b: 34 };
  var ys0 = DATA.map(function (r) { return r[key]; });
  var yMin = Math.min.apply(null, ys0), yMax = Math.max.apply(null, ys0);
  var pad = (yMax - yMin) * 0.05; yMin -= pad; yMax += pad;
  var yrs = META.years, xs = scaleLin(yrs[0], yrs[yrs.length - 1], m.l, w - m.r), ys = scaleLin(yMin, yMax, h - m.b, m.t);
  var svg = svgEl(w, h); gridY(svg, m, w, h, 6);
  var means = [], p25 = [], p75 = [];
  yrs.forEach(function (yr) { var v = DATA.filter(function (r) { return r.year === yr; }).map(function (r) { return r[key]; }); means.push([xs(yr), ys(E.mean(v))]); p25.push([xs(yr), ys(E.quantile(v, 0.25))]); p75.push([xs(yr), ys(E.quantile(v, 0.75))]); });
  pick.forEach(function (id) { var pts = by[id].slice().sort(function (a, b) { return a.year - b.year; }).map(function (r) { return [xs(r.year), ys(r[key])]; }); svg.appendChild(polyline(pts, "#4b5563", 1, 0.4)); });
  svg.appendChild(polyline(p25, "#a78bfa", 1.4, 0.9)); svg.appendChild(polyline(p75, "#a78bfa", 1.4, 0.9));
  svg.appendChild(polyline(means, "#93c5fd", 2.6, 1));
  svg.appendChild(line(m.l, h - m.b, w - m.r, h - m.b, "#374151", 1));
  yrs.forEach(function (yr) { svg.appendChild(txt(xs(yr), h - m.b + 14, String(yr).slice(2), "#9aa0a6", 9)); });
  [0, 0.5, 1].forEach(function (p) { svg.appendChild(txt(m.l - 6, ys(yMin + p * (yMax - yMin)) + 3, fmt(yMin + p * (yMax - yMin), 1), "#9aa0a6", 9, "end")); });
  svg.appendChild(txt(w / 2, h - 4, "year", "#6b7280", 10));
  setChart("chart-spaghetti", svg);
}

function drawDist() {
  var key = outcomeKey(), vals = DATA.map(function (r) { return r[key]; });
  var w = 480, h = 300, m = { l: 40, r: 12, t: 12, b: 34 };
  var mn = Math.min.apply(null, vals), mx = Math.max.apply(null, vals), bins = 30, bw = (mx - mn) / bins || 1, cnt = new Array(bins).fill(0);
  vals.forEach(function (v) { var i = Math.max(0, Math.min(bins - 1, Math.floor((v - mn) / bw))); cnt[i]++; });
  var maxC = Math.max.apply(null, cnt), xs = scaleLin(0, bins, m.l, w - m.r), ys = scaleLin(0, maxC, h - m.b, m.t);
  var svg = svgEl(w, h); gridY(svg, m, w, h, 5);
  var bwid = (w - m.l - m.r) / bins - 1;
  for (var i = 0; i < bins; i++) svg.appendChild(rect(xs(i) + 0.5, ys(cnt[i]), bwid, ys(0) - ys(cnt[i]), "#60a5fa"));
  svg.appendChild(line(m.l, h - m.b, w - m.r, h - m.b, "#374151", 1));
  [0, 0.5, 1].forEach(function (p) { svg.appendChild(txt(m.l + p * (w - m.l - m.r), h - m.b + 14, fmt(mn + p * (mx - mn), 1), "#9aa0a6", 9)); });
  svg.appendChild(txt(w / 2, h - 4, outcomeLabel(), "#6b7280", 10));
  setChart("chart-dist", svg);
  $("distLabel").textContent = outcomeLabel();
}

function focalGroup(r) { if (S.focal === "educ") return r.educ >= 12 ? 1 : 0; return r[S.focal]; }
function drawByGroup() {
  var key = outcomeKey(), yrs = META.years;
  var w = 480, h = 300, m = { l: 48, r: 12, t: 12, b: 34 };
  var s0 = [], s1 = [];
  yrs.forEach(function (yr) {
    var g0 = [], g1 = [];
    DATA.forEach(function (r) { if (r.year !== yr) return; (focalGroup(r) === 1 ? g1 : g0).push(r[key]); });
    s0.push({ yr: yr, v: g0.length ? E.mean(g0) : NaN }); s1.push({ yr: yr, v: g1.length ? E.mean(g1) : NaN });
  });
  var all = s0.concat(s1).map(function (p) { return p.v; }).filter(isFinite);
  var yMin = Math.min.apply(null, all), yMax = Math.max.apply(null, all), pad = (yMax - yMin) * 0.1; yMin -= pad; yMax += pad;
  var xs = scaleLin(yrs[0], yrs[yrs.length - 1], m.l, w - m.r), ys = scaleLin(yMin, yMax, h - m.b, m.t);
  var svg = svgEl(w, h); gridY(svg, m, w, h, 6);
  function draw(ser, col) { var pts = ser.filter(function (p) { return isFinite(p.v); }).map(function (p) { return [xs(p.yr), ys(p.v)]; }); svg.appendChild(polyline(pts, col, 2.4, 0.95)); pts.forEach(function (p) { svg.appendChild(circle(p[0], p[1], 2.6, col, 0.95)); }); }
  draw(s0, "#e5e7eb"); draw(s1, "#60a5fa");
  svg.appendChild(line(m.l, h - m.b, w - m.r, h - m.b, "#374151", 1));
  yrs.forEach(function (yr) { svg.appendChild(txt(xs(yr), h - m.b + 14, String(yr).slice(2), "#9aa0a6", 9)); });
  [0, 0.5, 1].forEach(function (p) { svg.appendChild(txt(m.l - 6, ys(yMin + p * (yMax - yMin)) + 3, fmt(yMin + p * (yMax - yMin), 1), "#9aa0a6", 9, "end")); });
  svg.appendChild(txt(w / 2, h - 4, S.focal === "educ" ? "year (focal: educ ≥ 12)" : "year", "#6b7280", 10));
  setChart("chart-bygroup", svg);
}

function drawVariance() {
  var keys = [S.focal, outcomeKey()], labels = [FOCAL_SHORT[S.focal], outcomeLabel()];
  var w = 480, h = 300, m = { l: 70, r: 80, t: 20, b: 30 };
  var svg = svgEl(w, h);
  var barH = 46, gap = 40, x0 = m.l, xW = w - m.l - m.r;
  keys.forEach(function (k, i) {
    var d = varianceDecomp(DATA, k), y = m.t + i * (barH + gap);
    svg.appendChild(rect(x0, y, xW * d.within, barH, "#22c55e"));
    svg.appendChild(rect(x0 + xW * d.within, y, xW * d.between, barH, "#6b7280"));
    svg.appendChild(txt(x0 - 8, y + barH / 2 + 4, labels[i], "#e5e7eb", 12, "end"));
    if (d.within > 0.06) svg.appendChild(txt(x0 + xW * d.within / 2, y + barH / 2 + 4, Math.round(d.within * 100) + "%", "#04210f", 12));
    if (d.between > 0.06) svg.appendChild(txt(x0 + xW * d.within + xW * d.between / 2, y + barH / 2 + 4, Math.round(d.between * 100) + "%", "#e5e7eb", 12));
  });
  svg.appendChild(txt(x0, m.t - 6, "within-person", "#22c55e", 11, "start"));
  svg.appendChild(txt(w - m.r, m.t - 6, "between-person", "#9aa0a6", 11, "end"));
  setChart("chart-variance", svg);
  var dv = varianceDecomp(DATA, S.focal);
  var note = Math.round(dv.within * 100) + "% of the variation in <strong>" + FOCAL_SHORT[S.focal] + "</strong> is within-person. ";
  note += dv.within < 0.02
    ? "Almost none — so Fixed Effects and First Differences <strong>cannot identify</strong> this effect."
    : "Fixed Effects uses only this within-person part; the " + Math.round(dv.between * 100) + "% between-person variation is discarded.";
  $("varianceNote").innerHTML = note;
}

// ============== MULTIVERSE ==============
function buildAxes() {
  return {
    estimators: S.axes.estimators.slice(),
    controlsVary: S.axes.controlsVary.slice(),
    controlsFixed: { experience: false, hours: false, industry: false, region: false, health: false, race: false },
    yearFEVary: S.axes.yearFEVary, yearFEFixed: false,
    sampleVary: S.axes.sampleVary.slice(),
    focal: S.focal, outcome: S.outcome, seType: S.seType
  };
}
function comboCount() {
  var a = S.axes;
  return a.estimators.length * (1 << a.controlsVary.length) * (a.yearFEVary ? 2 : 1) * Math.max(1, a.sampleVary.length);
}
function runMultiverse() {
  $("comboCount").innerHTML = "Computing " + comboCount() + " specifications…";
  setTimeout(function () {
    var t0 = performance.now();
    S.multiResults = E.enumerateMultiverse(DATA, buildAxes(), META);
    S.multiSummary = E.summarize(S.multiResults);
    var ms = Math.round(performance.now() - t0);
    $("comboCount").innerHTML = comboCount() + " specifications · " + ms + " ms · " + S.multiSummary.unidentified + " not identifiable";
    renderMulti();
  }, 20);
}

function drawSpecCurve() {
  var res = S.multiResults;
  if (!res) { emptyChart("chart-speccurve", "Click “Build the multiverse”"); return; }
  var ok = res.filter(function (r) { return r.identified && isFinite(r.b); }).slice().sort(function (a, b) { return a.b - b.b; });
  if (!ok.length) { emptyChart("chart-speccurve", "No identifiable specifications for this focal effect under these choices."); return; }
  var w = 920, h = 430, m = { l: 70, r: 16, t: 14, b: 8 };
  var rowsDef = [];
  S.axes.estimators.forEach(function (e) { rowsDef.push({ kind: "est", key: e, label: ESTIMATOR_LABEL[e] }); });
  S.axes.controlsVary.forEach(function (c) { rowsDef.push({ kind: "ctrl", key: c, label: CONTROL_LABEL[c] }); });
  if (S.axes.yearFEVary) rowsDef.push({ kind: "yfe", key: "yearFE", label: "year FE" });
  if (S.axes.sampleVary.length > 1) S.axes.sampleVary.forEach(function (s) { rowsDef.push({ kind: "samp", key: s, label: SAMPLE_LABEL[s] }); });
  var matrixH = rowsDef.length * 13 + 8, curveH = h - m.t - matrixH - 26;
  var n = ok.length, xAt = function (i) { return m.l + (n === 1 ? 0.5 : i / (n - 1)) * (w - m.l - m.r); };
  var lo = [], hi = [];
  ok.forEach(function (r) { var c = ciOf(r.b, r.se); lo.push(isFinite(c[0]) ? c[0] : r.b); hi.push(isFinite(c[1]) ? c[1] : r.b); });
  var yMin = Math.min(0, Math.min.apply(null, lo)), yMax = Math.max(0, Math.max.apply(null, hi));
  var pad = (yMax - yMin) * 0.08; yMin -= pad; yMax += pad;
  var ys = scaleLin(yMin, yMax, m.t + curveH, m.t);
  var svg = svgEl(w, h);
  for (var g = 0; g <= 4; g++) { var yy = m.t + g * curveH / 4; svg.appendChild(line(m.l, yy, w - m.r, yy, "#161616", 1)); svg.appendChild(txt(m.l - 6, yy + 3, fmt(yMax - g * (yMax - yMin) / 4, 2), "#9aa0a6", 9, "end")); }
  svg.appendChild(line(m.l, ys(0), w - m.r, ys(0), "#374151", 1.2));
  var med = S.multiSummary.median;
  svg.appendChild(line(m.l, ys(med), w - m.r, ys(med), "#f87171", 1.4));
  svg.appendChild(txt(w - m.r, ys(med) - 4, "median " + fmt(med, 3), "#f87171", 10, "end"));
  ok.forEach(function (r, i) {
    var x = xAt(i), c = ciOf(r.b, r.se), sig = isSig(r.b, r.se);
    if (isFinite(c[0])) svg.appendChild(line(x, ys(c[0]), x, ys(c[1]), sig ? "#22c55e" : "#6b7280", 0.8));
    svg.appendChild(circle(x, ys(r.b), 2.2, sig ? "#34d399" : "#9aa0a6", 0.95));
  });
  svg.appendChild(txt(m.l, m.t + curveH + 16, "← specifications sorted by estimate →", "#6b7280", 10, "start"));
  var my0 = m.t + curveH + 26;
  rowsDef.forEach(function (rd, ri) {
    var ry = my0 + ri * 13 + 6;
    svg.appendChild(txt(m.l - 6, ry + 3, rd.label, "#9aa0a6", 9, "end"));
    svg.appendChild(line(m.l, ry, w - m.r, ry, "#141414", 1));
    ok.forEach(function (r, i) {
      var active = rd.kind === "est" ? r.spec.estimator === rd.key
        : rd.kind === "ctrl" ? !!r.spec.controls[rd.key]
        : rd.kind === "yfe" ? !!r.spec.yearFE
        : r.spec.sample === rd.key;
      if (active) svg.appendChild(circle(xAt(i), ry, 1.7, isSig(r.b, r.se) ? "#34d399" : "#6b7280", 0.9));
    });
  });
  setChart("chart-speccurve", svg);
}

function renderMulti() {
  if (!S.multiResults) $("comboCount").innerHTML = "Will compute " + comboCount() + " specifications when you build it.";
  var sum = S.multiSummary, cards = $("multiCards"); cards.innerHTML = "";
  function card(title, val, foot, hi) { var d = document.createElement("div"); d.className = "card" + (hi ? " highlight" : ""); d.innerHTML = "<div class='title'>" + title + "</div><div class='value'>" + val + "</div><div class='foot'>" + (foot || "") + "</div>"; cards.appendChild(d); }
  if (sum) {
    card("Median estimate", fmt(sum.median, 3), FOCAL_LABEL[S.focal] + " (" + outcomeLabel() + ")", true);
    card("Range", fmt(sum.min, 3) + " … " + fmt(sum.max, 3), "min to max across the universe");
    card("Interquartile", fmt(sum.q25, 3) + " … " + fmt(sum.q75, 3), "middle 50% of specs");
    card("Significant", Math.round(sum.shareSig * 100) + "%", "of " + sum.identified + " identifiable specs");
    card("Sign flips", Math.round(sum.shareFlip * 100) + "%", "opposite sign to median");
    card("Not identifiable", sum.unidentified, "FE/FD on a fixed regressor");
  }
  drawSpecCurve();
  var box = $("multiInterpretation");
  if (sum) {
    var spread = sum.max - sum.min;
    var verdict = (sum.shareFlip > 0.02 || spread > Math.abs(sum.median))
      ? "<span class='badge badge-danger'>fragile</span> The estimate swings widely — the headline number depends heavily on analytic choices."
      : (sum.shareSig > 0.9 ? "<span class='badge badge-success'>robust</span> The effect keeps the same sign and is mostly significant across the universe."
        : "<span class='badge badge-warning'>mixed</span> The sign is stable but significance is sensitive to choices.");
    var driver = "";
    if (S.focal !== "educ") {
      var byEst = {};
      S.multiResults.forEach(function (r) { if (r.identified && isFinite(r.b)) (byEst[r.estimator] || (byEst[r.estimator] = [])).push(r.b); });
      var pooledM = byEst.pooled ? E.mean(byEst.pooled) : NaN, feM = byEst.within ? E.mean(byEst.within) : NaN;
      if (isFinite(pooledM) && isFinite(feM)) driver = " Within the universe, Pooled OLS averages " + fmt(pooledM, 3) + " but Fixed Effects averages " + fmt(feM, 3) + " — the gap is the part that was selection on stable traits.";
    } else {
      driver = " Every Fixed-Effects and First-Difference specification is dropped: schooling does not change within a person here, so those designs are silent about it.";
    }
    box.innerHTML = "<h4>The universe of the " + FOCAL_LABEL[S.focal] + "</h4><p>" + verdict + driver + "</p>";
  } else box.innerHTML = "<h4>Build the multiverse</h4><p>Choose which forks to open on the left, then click the button. Each combination is a complete, defensible analysis.</p>";
}

// ============== SINGLE SPEC ==============
function ssSpec(estimator) {
  return { focal: S.focal, outcome: S.outcome, controls: Object.assign({}, S.ss.controls), yearFE: S.ss.yearFE, sample: S.ss.sample, estimator: estimator, seType: S.seType };
}
function drawEstBars(items) {
  var w = 480, h = 300, m = { l: 50, r: 14, t: 16, b: 52 };
  var vals = [0];
  items.forEach(function (it) { if (it.identified && isFinite(it.b)) { vals.push(it.b); var c = ciOf(it.b, it.se); vals.push(c[0], c[1]); } });
  var yMin = Math.min.apply(null, vals), yMax = Math.max.apply(null, vals), pad = (yMax - yMin) * 0.15 || 0.1; yMin -= pad; yMax += pad;
  var bw = (w - m.l - m.r) / items.length, ys = scaleLin(yMin, yMax, h - m.b, m.t);
  var svg = svgEl(w, h); gridY(svg, m, w, h, 6);
  svg.appendChild(line(m.l, ys(0), w - m.r, ys(0), "#374151", 1.2));
  [0, 0.5, 1].forEach(function (p) { svg.appendChild(txt(m.l - 6, ys(yMin + p * (yMax - yMin)) + 3, fmt(yMin + p * (yMax - yMin), 2), "#9aa0a6", 9, "end")); });
  items.forEach(function (it, i) {
    var cx = m.l + i * bw + bw / 2;
    if (!it.identified || !isFinite(it.b)) {
      svg.appendChild(rect(cx - bw * 0.25, ys(0) - 1, bw * 0.5, 2, "#3a3a3a"));
      svg.appendChild(txt(cx, ys(0) - 8, "n.i.", "#6b7280", 10));
    } else {
      var y0 = ys(0), bar = ys(it.b), col = it.estimator === "within" ? "#a78bfa" : "#60a5fa";
      svg.appendChild(rect(cx - bw * 0.28, Math.min(bar, y0), bw * 0.56, Math.abs(bar - y0), col));
      var c = ciOf(it.b, it.se), ylo = ys(c[0]), yhi = ys(c[1]);
      svg.appendChild(line(cx, ylo, cx, yhi, "#e5e7eb", 1.6));
      svg.appendChild(line(cx - 6, ylo, cx + 6, ylo, "#e5e7eb", 1.6));
      svg.appendChild(line(cx - 6, yhi, cx + 6, yhi, "#e5e7eb", 1.6));
      svg.appendChild(txt(cx, Math.min(bar, y0) - 5, fmt(it.b, 3), "#cbd5e1", 9));
    }
    svg.appendChild(txt(cx, h - m.b + 14, ESTIMATOR_LABEL[it.estimator].split(" ")[0], "#9aa0a6", 9));
    svg.appendChild(txt(cx, h - m.b + 26, ESTIMATOR_LABEL[it.estimator].split(" ").slice(1).join(" "), "#6b7280", 8));
  });
  setChart("chart-estbars", svg);
}

function regressionTable(model, focal) {
  if (!model || !model.ok) return "<p class='muted'>This estimator does not identify a model for the current choices.</p>";
  var rows = model.coefNames.map(function (nm, i) {
    var b = model.beta[i], se = model.se[i], t = model.tstat[i], p = model.pval[i];
    var star = p < 0.001 ? "***" : p < 0.01 ? "**" : p < 0.05 ? "*" : p < 0.1 ? "." : "";
    var hl = nm === focal ? " style='background:#0c1929'" : "";
    return "<tr" + hl + "><td>" + nm + (nm === focal ? " ◀ focal" : "") + "</td><td class='num'>" + fmt(b, 4) + "</td><td class='num'>" + fmt(se, 4) + "</td><td class='num'>" + fmt(t, 2) + "</td><td class='num'>" + (p < 0.0001 ? "&lt;0.0001" : fmt(p, 4)) + " " + star + "</td></tr>";
  }).join("");
  var foot = "n = " + model.n + (model.nGroups ? " · persons = " + model.nGroups : "") + " · R² = " + fmt(model.r2, 3) + " · SE: " + model.seType + (model.theta != null ? " · θ = " + fmt(model.theta, 3) : "");
  return "<table class='stat-table'><thead><tr><th>Coefficient</th><th class='num'>Estimate</th><th class='num'>Std. Err.</th><th class='num'>t</th><th class='num'>p</th></tr></thead><tbody>" + rows + "</tbody></table><p class='help-text mt-2'>" + foot + "</p>";
}

var DIAG_W = 470, DIAG_H = 350, DIAG_M = { l: 58, r: 16, t: 16, b: 46 };
function diagSvg() { var s = svgEl(DIAG_W, DIAG_H); s.setAttribute("role", "img"); return s; }

function drawResiduals(model) {
  if (!model || !model.resid || !model.resid.length) { emptyChart("chart-rvf", "Not identified for this specification."); ariaChart("chart-rvf", "Residuals-versus-fitted plot, not available."); return; }
  var pts = model.fitted.map(function (f, i) { return { fit: f, r: model.resid[i] }; });
  var w = DIAG_W, h = DIAG_H, m = DIAG_M;
  var fmn = Math.min.apply(null, pts.map(function (p) { return p.fit; })), fmx = Math.max.apply(null, pts.map(function (p) { return p.fit; }));
  var rmx = Math.max.apply(null, pts.map(function (p) { return Math.abs(p.r); }));
  var xt = niceTicks(fmn, fmx, 5), yt = niceTicks(-rmx, rmx, 5);
  var xs = scaleLin(xt.lo, xt.hi, m.l, w - m.r), ys = scaleLin(yt.lo, yt.hi, h - m.b, m.t);
  var svg = diagSvg();
  plotFrame(svg, m, w, h, xs, ys, xt, yt, "Fitted value of the outcome", "Residual");
  var z = line(m.l, ys(0), w - m.r, ys(0), "#f87171", 1.4); z.setAttribute("stroke-dasharray", "5,4"); svg.appendChild(z);
  svg.appendChild(txt(w - m.r - 2, ys(0) - 5, "residual = 0", "#f87171", 10, "end"));
  pts.forEach(function (p) { svg.appendChild(circle(xs(p.fit), ys(p.r), 2.2, "#60a5fa", 0.4)); });
  setChart("chart-rvf", svg);
  ariaChart("chart-rvf", "Scatter plot of residuals against fitted values for " + ESTIMATOR_LABEL[S.ss.estimator] + "; a dashed line marks zero. Look for a flat, patternless band.");
}

function normalPdf(x, mu, sd) { return Math.exp(-(x - mu) * (x - mu) / (2 * sd * sd)) / (sd * Math.sqrt(2 * Math.PI)); }
function drawHist(model) {
  if (!model || !model.resid || !model.resid.length) { emptyChart("chart-hist", "Not identified for this specification."); ariaChart("chart-hist", "Residual histogram, not available."); return; }
  var r = model.resid, w = DIAG_W, h = DIAG_H, m = DIAG_M;
  var mn = Math.min.apply(null, r), mx = Math.max.apply(null, r), bins = 24, bw = (mx - mn) / bins || 1, cnt = new Array(bins).fill(0);
  r.forEach(function (v) { cnt[Math.max(0, Math.min(bins - 1, Math.floor((v - mn) / bw)))]++; });
  var mu = E.mean(r), sd = Math.sqrt(E.variance(r));
  var densAsCount = function (x) { return normalPdf(x, mu, sd) * r.length * bw; }; // scale density to count axis
  var peak = densAsCount(mu);
  var maxC = Math.max(Math.max.apply(null, cnt), peak);
  var xt = niceTicks(mn, mx, 5), yt = niceTicks(0, maxC, 5);
  var xs = scaleLin(xt.lo, xt.hi, m.l, w - m.r), ys = scaleLin(yt.lo, yt.hi, h - m.b, m.t);
  var svg = diagSvg();
  plotFrame(svg, m, w, h, xs, ys, xt, yt, "Residual", "Count");
  for (var i = 0; i < bins; i++) { var x0 = xs(mn + i * bw), x1 = xs(mn + (i + 1) * bw); svg.appendChild(rect(x0 + 0.5, ys(cnt[i]), Math.max(1, x1 - x0 - 1), ys(0) - ys(cnt[i]), "#3b82f6")); }
  var dpts = []; for (var s = 0; s <= 60; s++) { var xv = mn + (mx - mn) * s / 60; dpts.push([xs(xv), ys(densAsCount(xv))]); }
  svg.appendChild(polyline(dpts, "#f59e0b", 2, 0.95));
  svg.appendChild(txt(w - m.r - 2, m.t + 10, "normal fit", "#f59e0b", 10, "end"));
  setChart("chart-hist", svg);
  ariaChart("chart-hist", "Histogram of residuals for " + ESTIMATOR_LABEL[S.ss.estimator] + " with a fitted normal density overlaid. Look for an approximately symmetric, bell-shaped distribution.");
}

function drawQQ(model) {
  if (!model || !model.resid || model.resid.length < 5) { emptyChart("chart-qq", "Not identified for this specification."); ariaChart("chart-qq", "Normal quantile-quantile plot, not available."); return; }
  var r = model.resid.slice().sort(function (a, b) { return a - b; });
  function erfinv(x) { var a = 0.147, ln = Math.log(1 - x * x), s = 2 / (Math.PI * a) + ln / 2; return Math.sign(x) * Math.sqrt(Math.sqrt(s * s - ln / a) - s); }
  var n = r.length, w = DIAG_W, h = DIAG_H, m = DIAG_M, pts = [];
  for (var i = 0; i < n; i++) { var p = (i + 0.5) / n; pts.push({ t: Math.SQRT2 * erfinv(2 * p - 1), s: r[i] }); }
  var tmin = pts[0].t, tmax = pts[n - 1].t;
  var xt = niceTicks(tmin, tmax, 5), yt = niceTicks(r[0], r[n - 1], 5);
  var xs = scaleLin(xt.lo, xt.hi, m.l, w - m.r), ys = scaleLin(yt.lo, yt.hi, h - m.b, m.t);
  var svg = diagSvg();
  plotFrame(svg, m, w, h, xs, ys, xt, yt, "Theoretical quantile (standard normal)", "Sample quantile (residual)");
  // reference line through the data's robust slope: mean ± sd
  var mu = E.mean(r), sd = Math.sqrt(E.variance(r));
  var rx0 = xt.lo, rx1 = xt.hi;
  svg.appendChild(polyline([[xs(rx0), ys(mu + sd * rx0)], [xs(rx1), ys(mu + sd * rx1)]], "#f87171", 1.6));
  svg.appendChild(txt(w - m.r - 2, m.t + 10, "normal reference", "#f87171", 10, "end"));
  pts.forEach(function (p) { svg.appendChild(circle(xs(p.t), ys(p.s), 2.2, "#60a5fa", 0.45)); });
  setChart("chart-qq", svg);
  ariaChart("chart-qq", "Normal quantile-quantile plot of residuals for " + ESTIMATOR_LABEL[S.ss.estimator] + " with a reference line. Points on the line indicate normality.");
}

// ============== LIVE EQUATION ("under the hood") ==============
var SUB_IT = "<sub>it</sub>", SUB_I = "<sub>i</sub>", SUB_T = "<sub>t</sub>";
function vsym(nm, focal) {
  if (nm === focal) return "<span class='eq-focal'>" + FOCAL_SHORT[focal] + "</span>";
  if (nm === "expersq") return "exper²";
  if (nm === "(Intercept)") return "1";
  if (nm.indexOf("mean(") === 0) return "<span class='eq-note'>" + nm.replace("mean(" + focal + ")", "mean(" + FOCAL_SHORT[focal] + ")") + "</span>";
  return nm.replace(/^ind:/, "ind=").replace(/^reg:/, "reg=").replace(/^yr:/, "year=");
}
function fittedEquation(model, focal, est) {
  if (!model || !model.ok) return "<span class='eq-note'>not identified</span>";
  var lhs = est === "fd" ? "Δ " + yName() + SUB_IT : (est === "within" || est === "twfe") ? "(" + yName() + SUB_IT + " − " + yName() + "̄" + SUB_I + (est === "twfe" ? " − …" : "") + ")" : yName() + "̂" + SUB_IT;
  var shownKeys = { "(Intercept)": 1, "(drift)": 1 }; shownKeys[focal] = 1; shownKeys["exper"] = 1; shownKeys["expersq"] = 1; shownKeys["mean(" + focal + ")"] = 1;
  var parts = [], hidden = 0, first = true;
  model.coefNames.forEach(function (nm, i) {
    var b = model.beta[i]; if (!isFinite(b)) return;
    if (!shownKeys[nm]) { hidden++; return; }
    if (nm === "(Intercept)" || nm === "(drift)") { parts.push("<span class='eq-num'>" + fmt(b, 3) + "</span>"); first = false; return; }
    var sign = b >= 0 ? (first ? "" : " + ") : " − ";
    parts.push(sign + "<span class='eq-num'>" + fmt(Math.abs(b), 3) + "</span>·" + vsym(nm, focal) + (nm.indexOf("mean(") === 0 ? SUB_I : (est === "fd" ? "Δ" : SUB_IT)));
    first = false;
  });
  var tail = hidden ? " <span class='eq-note'>+ " + hidden + " more control term" + (hidden > 1 ? "s" : "") + "</span>" : "";
  return lhs + " = " + parts.join("") + tail;
}
function seFormula(seType, nGroups) {
  if (seType === "cluster") return "Var(β̂) = (X′X)⁻¹ (Σ<sub>g</sub> X<sub>g</sub>′û<sub>g</sub>û<sub>g</sub>′X<sub>g</sub>) (X′X)⁻¹ · CR1 — clustered by person (" + nGroups + " clusters)";
  if (seType === "robust") return "Var(β̂) = (X′X)⁻¹ (Σ<sub>i</sub> û<sub>i</sub>² x<sub>i</sub>x<sub>i</sub>′) (X′X)⁻¹ · HC1 — heteroskedasticity-robust";
  return "Var(β̂) = σ̂² (X′X)⁻¹ — classical, assuming i.i.d. homoskedastic errors";
}
function yName() { return S.outcome === "wage" ? "wage" : "lwage"; }
var EST_TRANSFORM = {
  pooled: { lbl: "Pooled OLS — fit on the levels, all variation", eq: function (f) { return yName() + SUB_IT + " = α + β·<span class='eq-focal'>" + f + "</span>" + SUB_IT + " + γ′z" + SUB_IT + " + (u" + SUB_I + " + ε" + SUB_IT + ")  <span class='cap'>← person effect left in the error</span>"; } },
  between: { lbl: "Between — collapse to one mean per person", eq: function (f) { return yName() + "̄" + SUB_I + " = α + β·<span class='eq-focal'>" + f + "</span>̄" + SUB_I + " + γ′z̄" + SUB_I + " + (u" + SUB_I + " + ε̄" + SUB_I + ")  <span class='cap'>← cross-person comparison only</span>"; } },
  within: { lbl: "Fixed effects — demean within person; uᵢ cancels", eq: function (f) { return "(" + yName() + SUB_IT + " − " + yName() + "̄" + SUB_I + ") = β·(<span class='eq-focal'>" + f + "</span>" + SUB_IT + " − <span class='eq-focal'>" + f + "</span>̄" + SUB_I + ") + γ′(z" + SUB_IT + " − z̄" + SUB_I + ") + (ε" + SUB_IT + " − ε̄" + SUB_I + ")"; } },
  fd: { lbl: "First differences — year-to-year change; uᵢ cancels", eq: function (f) { return "Δ" + yName() + SUB_IT + " = δ + β·Δ<span class='eq-focal'>" + f + "</span>" + SUB_IT + " + γ′Δz" + SUB_IT + " + Δε" + SUB_IT + "  <span class='cap'>← δ is the common drift</span>"; } },
  random: { lbl: "Random effects — GLS quasi-demeaning by θ", eq: function (f, m) { var th = m && isFinite(m.theta) ? fmt(m.theta, 3) : "θ"; return "(" + yName() + SUB_IT + " − <span class='eq-num'>" + th + "</span>·" + yName() + "̄" + SUB_I + ") = α(1−θ) + β·(<span class='eq-focal'>" + f + "</span>" + SUB_IT + " − θ·<span class='eq-focal'>" + f + "</span>̄" + SUB_I + ") + …  <span class='cap'>θ=0 ⇒ pooled, θ→1 ⇒ fixed effects</span>"; } },
  cre: { lbl: "Correlated RE (Mundlak) — add person-means to RE", eq: function (f) { return yName() + SUB_IT + " = α + β·<span class='eq-focal'>" + f + "</span>" + SUB_IT + " + δ·<span class='eq-focal'>" + f + "</span>̄" + SUB_I + " + γ′z" + SUB_IT + " + u" + SUB_I + " + ε" + SUB_IT + "  <span class='cap'>β equals the fixed-effects estimate</span>"; } },
  twfe: { lbl: "Two-way fixed effects — within person and within year", eq: function (f) { return "(" + yName() + SUB_IT + " − " + yName() + "̄" + SUB_I + " − " + yName() + "̄" + SUB_T + " + " + yName() + "̄) = β·(<span class='eq-focal'>" + f + "</span> demeaned both ways) + …  <span class='cap'>removes person and year effects</span>"; } }
};
function renderEquation(spec, model, est) {
  var box = $("eqContent");
  if (!$("showEq").checked) { box.innerHTML = "<span class='muted'>Equations hidden — tick “show equations”.</span>"; return; }
  var f = FOCAL_SHORT[spec.focal], T = EST_TRANSFORM[est] || EST_TRANSFORM.pooled;
  var ctrlText = controlsSummary(spec) + (spec.yearFE && est !== "twfe" ? " + year effects" : "");
  var html = "";
  html += "<div class='eqblock'><div class='lbl'>What the symbols mean</div><div class='eqmath' style='font-size:13px;line-height:1.9'>" +
    "<b>i</b> = a person (unit) &nbsp;·&nbsp; <b>t</b> = a year (time) &nbsp;·&nbsp; subscript <b>" + SUB_IT.replace(/<\/?sub>/g, "") + "</b> means “person i in year t”<br>" +
    "<b>" + yName() + "</b>" + SUB_IT + " = the outcome (" + outcomeLabel() + ") &nbsp;·&nbsp; <b class='eq-focal'>" + f + "</b>" + SUB_IT + " = the focal variable (the " + FOCAL_LABEL[spec.focal] + ")<br>" +
    "<b>α</b> = intercept / baseline &nbsp;·&nbsp; <b>β</b> = the effect being estimated (the coefficient on " + f + ") &nbsp;·&nbsp; <b>γ</b> = coefficients on the control variables <b>z</b><br>" +
    "<b>u</b>" + SUB_I + " = the <em>person effect</em> — stable traits of person i that don't change over time (unobserved) &nbsp;·&nbsp; <b>ε</b>" + SUB_IT + " = the <em>idiosyncratic error</em>, random year-to-year variation<br>" +
    "an <b>overbar</b> (e.g. " + yName() + "̄" + SUB_I + ") = that person's average across their years &nbsp;·&nbsp; <b>Δ</b> = year-to-year change &nbsp;·&nbsp; <b>θ</b> = random-effects weight (0 → pooled, 1 → fixed effects) &nbsp;·&nbsp; a <b>hat</b> (β̂) = an estimated value</div></div>" +
    "<div class='eqblock'><div class='lbl'>1 · Structural model</div><div class='eqmath'>" +
    yName() + SUB_IT + " = α + β·<span class='eq-focal'>" + f + "</span>" + SUB_IT + " + γ′z" + SUB_IT + " + u" + SUB_I + " + ε" + SUB_IT +
    "<br><span class='cap'>z = { " + ctrlText + " } · u" + SUB_I + " = stable person trait (unobserved) · the target is β, the " + FOCAL_LABEL[spec.focal] + "</span></div></div>";
  html += "<div class='eqblock'><div class='lbl'>2 · What " + ESTIMATOR_LABEL[est] + " actually fits</div><div class='eqmath'>" + T.eq(f, model) + "</div></div>";
  html += "<div class='eqblock'><div class='lbl'>3 · Fitted with your data</div><div class='eqmath'>" + fittedEquation(model, spec.focal, est) +
    ((est === "within" || est === "twfe") ? "  <span class='cap'>(+ " + (model && model.nGroups ? model.nGroups : "") + " person intercepts" + (est === "twfe" ? " + year effects" : "") + ", absorbed)</span>" : "") + "</div></div>";
  html += "<div class='eqblock'><div class='lbl'>4 · Standard errors (" + spec.seType + ")</div><div class='eqmath' style='font-size:13px'>" + seFormula(spec.seType, model && model.nGroups ? model.nGroups : "G") + "</div></div>";
  box.innerHTML = html;
}

function renderSingle() {
  var focal = S.focal;
  var bars = E.ALL_ESTIMATORS.map(function (e) { return E.runSpec(DATA, ssSpec(e), META); });
  drawEstBars(bars);
  // The classic Hausman test is valid under the conventional (non-robust) covariance, so it is
  // computed from classical-vcov FE and RE models regardless of the SE shown elsewhere.
  function classicalSpec(est) { var s = ssSpec(est); s.seType = "classical"; return s; }
  var feM = E.runSpec(DATA, classicalSpec("within"), META).model, reM = E.runSpec(DATA, classicalSpec("random"), META).model;
  var h = E.hausmanFull(feM, reM), hb = $("hausmanResult");
  if (h.ok) {
    var concl = h.pval < 0.05 ? "<span class='badge badge-success'>Reject H₀</span> Fixed effects preferred — the FE and RE coefficients differ systematically, signalling correlation between the person effect and the regressors." : "<span class='badge badge-info'>Cannot reject H₀</span> Random effects may be the more efficient choice here.";
    hb.innerHTML = "H = " + fmt(h.stat, 2) + " (" + h.df + " df), p = " + (h.pval < 0.0001 ? "&lt;0.0001" : fmt(h.pval, 4)) + "<div class='mt-2'>" + concl + "</div><div class='help-text mt-2'>Joint test over the " + h.df + " time-varying coefficient(s) common to FE and RE.</div>";
  } else hb.innerHTML = "<span class='muted'>Hausman test unavailable: " + (h.msg || "not identified") + ".</span>";
  var sel = E.runSpec(DATA, ssSpec(S.ss.estimator), META);
  renderEquation(ssSpec(S.ss.estimator), sel.model, S.ss.estimator);
  $("ssTitle").textContent = ESTIMATOR_LABEL[S.ss.estimator] + " · " + FOCAL_LABEL[focal];
  $("ssTable").innerHTML = sel.identified ? regressionTable(sel.model, focal) : "<div class='alert alert-warning'>" + ESTIMATOR_LABEL[S.ss.estimator] + " cannot identify <strong>" + FOCAL_SHORT[focal] + "</strong>: it does not vary within persons, so the within/difference transformation removes it. This is the central lesson about time-invariant regressors.</div>";
  drawResiduals(sel.model); drawHist(sel.model); drawQQ(sel.model);
  var pooled = bars.filter(function (b) { return b.estimator === "pooled"; })[0], within = bars.filter(function (b) { return b.estimator === "within"; })[0];
  var box = $("ssInterpretation"), msg = "";
  if (pooled && pooled.identified && within && within.identified) {
    var pct = Math.round(100 * (1 - within.b / pooled.b));
    msg = "Pooled OLS puts the " + FOCAL_LABEL[focal] + " at <strong>" + fmt(pooled.b, 3) + "</strong>; Fixed Effects cuts it to <strong>" + fmt(within.b, 3) + "</strong> — about " + pct + "% smaller. That difference is confounding by stable person traits that FE removes.";
  } else if (within && !within.identified) {
    msg = "Fixed Effects and First Differences return nothing for <strong>" + FOCAL_SHORT[focal] + "</strong>: it is constant within each person over 1980–87. Only Pooled OLS, Between and Random Effects can speak to it — and they lean on the assumption that schooling is uncorrelated with unobserved ability.";
  }
  box.innerHTML = "<h4>What the estimators disagree about</h4><p>" + msg + "</p>";
}

// ============== SIDEBAR / SHARED ==============
function renderSidebar() {
  $("focalHelp").textContent = FOCAL_HELP[S.focal];
  var persons = Object.keys(byPerson(DATA)).length;
  $("dataSummary").innerHTML = "<strong>" + persons + "</strong> men · <strong>" + DATA.length + "</strong> person-years · " + META.years[0] + "–" + META.years[META.years.length - 1] +
    "<br>Outcome: <span class='pill'>" + outcomeLabel() + "</span> · focal: <span class='pill'>" + FOCAL_SHORT[S.focal] + "</span>";
  $("learningObjective").innerHTML = "One dataset, one question — the <strong>" + FOCAL_LABEL[S.focal] + "</strong> — and many defensible analyses. The spread of results across reasonable analytic choices is the “garden of forking paths.” Build that garden by hand on real wage data and see which conclusions hold.";
}
// ---- Data Overview: descriptives & switchers ----
function varStats(key) {
  var by = byPerson(DATA), vals = DATA.map(function (r) { return r[key]; });
  var m = E.mean(vals), sd = Math.sqrt(E.variance(vals)), ids = Object.keys(by);
  var means = [], wss = 0;
  ids.forEach(function (id) { var g = by[id].map(function (r) { return r[key]; }), mi = E.mean(g); means.push(mi); g.forEach(function (v) { wss += (v - mi) * (v - mi); }); });
  return { mean: m, sd: sd, withinSD: Math.sqrt(wss / Math.max(1, DATA.length - ids.length)), betweenSD: Math.sqrt(E.variance(means)) };
}
function switchersInfo() {
  var by = byPerson(DATA), f = S.focal, a0 = 0, a1 = 0, sw = 0;
  Object.keys(by).forEach(function (id) {
    var v = by[id].map(function (r) { return r[f]; }), mn = Math.min.apply(null, v), mx = Math.max.apply(null, v);
    if (f === "educ") { if (mx > mn) sw++; else a0++; }
    else { if (mx === 0) a0++; else if (mn === 1) a1++; else sw++; }
  });
  return { a0: a0, a1: a1, sw: sw, N: Object.keys(by).length };
}
function renderDataExtras() {
  var persons = Object.keys(byPerson(DATA)).length, key = outcomeKey();
  // cards
  var cards = $("dataCards"); cards.innerHTML = "";
  function card(t, v, f) { var d = document.createElement("div"); d.className = "card"; d.innerHTML = "<div class='title'>" + t + "</div><div class='value'>" + v + "</div><div class='foot'>" + (f || "") + "</div>"; cards.appendChild(d); }
  card("Workers", persons, "panel units (i)");
  card("Person-years", DATA.length, META.years[0] + "–" + META.years[META.years.length - 1] + " (" + META.years.length + " waves)");
  card("Mean " + outcomeLabel(), fmt(E.mean(DATA.map(function (r) { return r[key]; })), 3), "");
  if (S.focal === "educ") {
    card("Mean schooling", fmt(E.mean(DATA.map(function (r) { return r.educ; })), 2) + " yrs", "fixed within person");
  } else {
    var g1 = DATA.filter(function (r) { return r[S.focal] === 1; }), g0 = DATA.filter(function (r) { return r[S.focal] === 0; });
    card("Share " + FOCAL_SHORT[S.focal] + " = 1", Math.round(100 * g1.length / DATA.length) + "%", "of person-years");
    card("Raw gap", fmt(E.mean(g1.map(function (r) { return r[key]; })) - E.mean(g0.map(function (r) { return r[key]; })), 3), "mean(1) − mean(0), no controls");
  }
  // descriptive table
  var rows = [["lwage", "log hourly wage"], ["wage", "hourly wage ($)"], [S.focal, FOCAL_SHORT[S.focal]], ["exper", "experience"], ["hours", "hours"], ["union", "union"], ["married", "married"], ["educ", "schooling"]]
    .filter(function (r, i, arr) { return arr.findIndex(function (x) { return x[0] === r[0]; }) === i; });
  var html = "<table class='stat-table'><thead><tr><th>Variable</th><th class='num'>Mean</th><th class='num'>SD</th><th class='num'>Within SD</th><th class='num'>Between SD</th></tr></thead><tbody>";
  rows.forEach(function (r) { var s = varStats(r[0]); html += "<tr><td>" + r[1] + "</td><td class='num'>" + fmt(s.mean, 2) + "</td><td class='num'>" + fmt(s.sd, 2) + "</td><td class='num'>" + fmt(s.withinSD, 2) + "</td><td class='num'>" + fmt(s.betweenSD, 2) + "</td></tr>"; });
  html += "</tbody></table><p class='help-text mt-2'>“Within SD” is variation over time for the same person; “Between SD” is variation across people. A variable with near-zero within SD (e.g. schooling) cannot be identified by Fixed Effects.</p>";
  $("descTable").innerHTML = html;
  // switchers
  var sw = switchersInfo(), pct = function (x) { return Math.round(100 * x / sw.N); };
  if (S.focal === "educ") {
    $("switchTable").innerHTML = "<table class='stat-table'><tbody><tr><td>Schooling constant within person</td><td class='num'>" + sw.a0 + " (" + pct(sw.a0) + "%)</td></tr><tr><td>Schooling changes within person</td><td class='num'>" + sw.sw + " (" + pct(sw.sw) + "%)</td></tr></tbody></table>";
    $("switchNote").innerHTML = "Essentially nobody changes years of schooling during 1980–87, so there are <strong>no switchers</strong> — Fixed Effects and First Differences have nothing to use and cannot identify the return to schooling.";
  } else {
    $("switchTable").innerHTML = "<table class='stat-table'><thead><tr><th>Group</th><th class='num'>Persons</th><th class='num'>Share</th></tr></thead><tbody>" +
      "<tr><td>Never " + FOCAL_SHORT[S.focal] + " (always 0)</td><td class='num'>" + sw.a0 + "</td><td class='num'>" + pct(sw.a0) + "%</td></tr>" +
      "<tr><td>Always " + FOCAL_SHORT[S.focal] + " (always 1)</td><td class='num'>" + sw.a1 + "</td><td class='num'>" + pct(sw.a1) + "%</td></tr>" +
      "<tr style='background:#0c1929'><td><strong>Switchers</strong></td><td class='num'><strong>" + sw.sw + "</strong></td><td class='num'><strong>" + pct(sw.sw) + "%</strong></td></tr></tbody></table>";
    $("switchNote").innerHTML = "Fixed Effects and First Differences identify the " + FOCAL_LABEL[S.focal] + " <strong>only from the " + sw.sw + " switchers</strong> (" + pct(sw.sw) + "%). The " + (sw.a0 + sw.a1) + " never/always persons contribute nothing to within estimates — which is why FE can be less precise.";
  }
}

function renderAll() {
  renderSidebar();
  if (S.viz === "data") { renderDataExtras(); drawSpaghetti(); drawDist(); drawByGroup(); drawVariance(); }
  else if (S.viz === "multi") renderMulti();
  else if (S.viz === "single") renderSingle();
}

// ============== GUIDED LEARNING ==============
function setFocal(f) {
  S.focal = f; S.multiResults = null; S.multiSummary = null;
  Array.prototype.forEach.call(document.getElementsByName("focal"), function (el) { el.checked = el.value === f; });
  renderAll(); syncHash();
}

var LESSON = [
  { n: "Step 1 of 6 · The question", t: "We'll estimate the <b>union wage premium</b> and learn why the estimate moves. We begin on Data Overview. Click Next.", act: function () { setFocal("union"); setViz("data"); } },
  { n: "Step 2 of 6 · Where the variation lives", t: "Find the <b>“Who identifies the effect? (switchers)”</b> table. Fixed Effects uses only the switchers — note how many there are versus the always/never groups.", act: function () { setViz("data"); } },
  { n: "Step 3 of 6 · Build the universe", t: "Open <b>Specification Curve</b> and click <b>Build the multiverse →</b>. Read the cards: median, range, % significant.", act: function () { setViz("multi"); } },
  { n: "Step 4 of 6 · Read the drivers", t: "In the dashboard below the curve, notice Fixed-Effects / First-Difference specs cluster at the <b>low</b> end and Pooled / Between at the <b>high</b> end. That gap is selection.", act: function () { setViz("multi"); } },
  { n: "Step 5 of 6 · Test the whole curve", t: "Click <b>Test the whole curve</b> to run randomization inference. Is the curve distinguishable from “no effect”?", act: function () { setViz("multi"); } },
  { n: "Step 6 of 6 · The identification limit", t: "Finally, set the focal effect to <b>Returns to schooling</b> in the left sidebar and rebuild. Half the universe vanishes — Fixed Effects cannot identify a regressor that never changes within a person. Lesson complete! ✕ to close.", act: function () { setViz("multi"); } }
];
function showLesson(i) {
  if (i < 0 || i >= LESSON.length) { closeLesson(); return; }
  S.lessonStep = i; var step = LESSON[i];
  $("lessonBar").style.display = "block"; document.body.style.paddingBottom = "96px";
  $("lessonStepN").textContent = step.n; $("lessonStepText").innerHTML = step.t;
  $("lessonPrev").disabled = (i === 0); $("lessonNext").textContent = (i === LESSON.length - 1) ? "Finish" : "Next →";
  if (step.act) step.act();
}
function closeLesson() { S.lessonStep = -1; $("lessonBar").style.display = "none"; document.body.style.paddingBottom = ""; }

// ============== SHAREABLE PERMALINKS ==============
function stateToHash() {
  var on = E.CONTROL_KEYS.concat(["race"]).filter(function (k) { return S.ss.controls[k]; });
  var p = { focal: S.focal, outcome: S.outcome, se: S.seType, tab: S.viz, est: S.ss.estimator, ctrl: on.join(","), yfe: S.ss.yearFE ? 1 : 0, samp: S.ss.sample };
  return "#" + Object.keys(p).map(function (k) { return k + "=" + encodeURIComponent(p[k]); }).join("&");
}
function syncHash() { try { history.replaceState(null, "", stateToHash()); } catch (e) {} }
function applyHash() {
  var h = (location.hash || "").replace(/^#/, ""); if (!h) return;
  var q = {}; h.split("&").forEach(function (s) { var kv = s.split("="); q[kv[0]] = decodeURIComponent(kv[1] || ""); });
  if (q.focal && FOCAL_LABEL[q.focal]) S.focal = q.focal;
  if (q.outcome === "wage" || q.outcome === "lwage") S.outcome = q.outcome;
  if (["classical", "robust", "cluster"].indexOf(q.se) >= 0) S.seType = q.se;
  if (q.est && E.ESTIMATORS_ALL.indexOf(q.est) >= 0) S.ss.estimator = q.est;
  if (q.samp && ["full", "nohealth", "trim"].indexOf(q.samp) >= 0) S.ss.sample = q.samp;
  if (q.yfe != null) S.ss.yearFE = q.yfe === "1";
  if (q.ctrl != null) { var arr = q.ctrl ? q.ctrl.split(",") : []; Object.keys(S.ss.controls).forEach(function (k) { S.ss.controls[k] = arr.indexOf(k) >= 0; }); }
  if (q.tab && ["data", "multi", "single", "export", "learn"].indexOf(q.tab) >= 0) S.viz = q.tab;
}
function syncDom() {
  Array.prototype.forEach.call(document.getElementsByName("focal"), function (el) { el.checked = el.value === S.focal; });
  $("outcome").value = S.outcome; $("seType").value = S.seType;
  $("ssEstimator").value = S.ss.estimator; $("ssSample").value = S.ss.sample; $("ssYearFE").checked = S.ss.yearFE;
}

// ============== UI BINDING ==============
function buildAxisControls() {
  var ec = $("axisEstimators"); ec.innerHTML = "";
  E.ESTIMATORS_ALL.forEach(function (e) {
    var l = document.createElement("label"); l.className = "opt";
    l.innerHTML = "<input type='checkbox' value='" + e + "' " + (S.axes.estimators.indexOf(e) >= 0 ? "checked" : "") + "> " + ESTIMATOR_LABEL[e];
    l.querySelector("input").addEventListener("change", function (ev) {
      var v = ev.target.value; if (ev.target.checked) { if (S.axes.estimators.indexOf(v) < 0) S.axes.estimators.push(v); } else S.axes.estimators = S.axes.estimators.filter(function (x) { return x !== v; });
      updateComboCount();
    });
    ec.appendChild(l);
  });
  var cc = $("axisControls"); cc.innerHTML = "";
  ["experience", "hours", "industry", "region", "health", "race"].forEach(function (c) {
    var l = document.createElement("label"); l.className = "opt";
    l.innerHTML = "<input type='checkbox' value='" + c + "' " + (S.axes.controlsVary.indexOf(c) >= 0 ? "checked" : "") + "> " + CONTROL_LABEL[c];
    l.querySelector("input").addEventListener("change", function (ev) {
      var v = ev.target.value; if (ev.target.checked) { if (S.axes.controlsVary.indexOf(v) < 0) S.axes.controlsVary.push(v); } else S.axes.controlsVary = S.axes.controlsVary.filter(function (x) { return x !== v; });
      updateComboCount();
    });
    cc.appendChild(l);
  });
}
function buildSSControls() {
  var box = $("ssControls"); box.innerHTML = "";
  ["experience", "hours", "industry", "region", "health", "race"].forEach(function (c) {
    var l = document.createElement("label"); l.className = "opt";
    l.innerHTML = "<input type='checkbox' value='" + c + "' " + (S.ss.controls[c] ? "checked" : "") + "> " + CONTROL_LABEL[c];
    l.querySelector("input").addEventListener("change", function (ev) { S.ss.controls[ev.target.value] = ev.target.checked; renderSingle(); syncHash(); });
    box.appendChild(l);
  });
}
function updateComboCount() { $("comboCount").innerHTML = (S.multiResults ? "Re-build to apply: " : "") + comboCount() + " specifications" + (comboCount() > 1500 ? " ⚠ large" : ""); }

function setViz(v) {
  S.viz = v;
  var map = { data: "dataArea", multi: "multiArea", single: "singleArea", export: "exportArea", learn: "learnArea" };
  Object.keys(map).forEach(function (k) { $(map[k]).style.display = k === v ? "block" : "none"; });
  ["data", "multi", "single", "export", "learn"].forEach(function (k) { $("viz-" + k).classList.toggle("active", k === v); });
  // the research-question sidebar is irrelevant on Methodology; give that tab the full width
  var full = (v === "learn");
  $("sidebarCol").style.display = full ? "none" : "";
  $("mainRow").style.gridTemplateColumns = full ? "1fr" : "";
  if (v === "export") { renderSidebar(); renderExport(); }
  else if (v !== "learn") renderAll();
  else renderSidebar();
  syncHash();
}

function controlsSummary(spec) {
  var on = Object.keys(spec.controls).filter(function (k) { return spec.controls[k] && !(k === "race" && spec.focal === "educ"); });
  return on.length ? on.map(function (k) { return CONTROL_LABEL[k]; }).join(", ") : "no controls";
}
function escapeHtml(s) { return s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;"); }

function renderExport() {
  var spec = ssSpec(S.ss.estimator);
  $("code-r").textContent = Exports.genR(spec);
  $("code-stata").textContent = Exports.genStata(spec);
  $("code-python").textContent = Exports.genPython(spec);
}

function buildLabReport() {
  var spec = ssSpec(S.ss.estimator), sel = E.runSpec(DATA, spec, META), sum = S.multiSummary;
  var table = sel.identified ? regressionTable(sel.model, S.focal) : "<p><em>This estimator does not identify the focal coefficient for the chosen options.</em></p>";
  var curve = Exports.chartSVGString("chart-speccurve") || "<p><em>(Build the specification curve to include it here.)</em></p>";
  var persons = Object.keys(byPerson(DATA)).length;
  var html = "<!doctype html><html lang='en'><head><meta charset='utf-8'><title>Panel Data Multiverse — Lab Report</title>" +
    "<style>body{font-family:system-ui,-apple-system,sans-serif;max-width:880px;margin:2rem auto;padding:0 1rem;color:#111;line-height:1.55}" +
    "table{border-collapse:collapse;width:100%;font-size:14px;margin:.5rem 0}th,td{border-bottom:1px solid #ddd;padding:6px 9px;text-align:left}" +
    ".num{text-align:right;font-variant-numeric:tabular-nums}pre{background:#f5f5f5;padding:10px;border-radius:6px;overflow:auto;font-size:12px}" +
    "h1{font-size:22px}h2{font-size:16px;margin-top:1.5rem;border-bottom:1px solid #eee;padding-bottom:4px}.muted{color:#666;font-size:12px}svg{max-width:100%;height:auto}</style></head><body>" +
    "<h1>Panel Data Multiverse — Lab Report</h1>" +
    "<p class='muted'>Generated " + new Date().toISOString().slice(0, 10) + " · NLSY wage panel: " + persons + " men, " + DATA.length + " person-years, 1980–1987.</p>" +
    "<h2>Specification</h2><p>" + Exports.label(spec) + "<br>Outcome: " + outcomeLabel() + " · Controls: " + controlsSummary(spec) +
    (spec.yearFE ? " · year fixed effects" : "") + " · Sample: " + SAMPLE_LABEL[spec.sample] + " · Standard errors: " + spec.seType + "</p>" +
    "<h2>Regression results</h2>" + table +
    (sum ? "<h2>Multiverse summary</h2><p>Across the most recent build: median estimate <strong>" + fmt(sum.median, 3) + "</strong>, range [" + fmt(sum.min, 3) + ", " + fmt(sum.max, 3) + "], " +
      Math.round(sum.shareSig * 100) + "% statistically significant, " + Math.round(sum.shareFlip * 100) + "% sign flips, " + sum.unidentified + " not identifiable.</p>" +
      "<h2>Specification curve</h2>" + curve : "") +
    "<h2>Reproducible code (R · plm)</h2><pre>" + escapeHtml(Exports.genR(spec)) + "</pre>" +
    "<hr><p class='muted'>Panel Data Multiverse Lab · estimates verified against R’s plm and Python’s linearmodels.</p></body></html>";
  Exports.download("panel-multiverse-report.html", html, "text/html");
}

function bindUI() {
  ["data", "multi", "single", "export", "learn"].forEach(function (k) { $("viz-" + k).addEventListener("click", function () { setViz(k); }); });
  // export tab: copy + download buttons
  Array.prototype.forEach.call(document.querySelectorAll("[data-copy]"), function (b) {
    b.addEventListener("click", function () { var t = $(b.getAttribute("data-copy")).textContent; navigator.clipboard.writeText(t).then(function () { var o = b.textContent; b.textContent = "Copied ✓"; setTimeout(function () { b.textContent = o; }, 1200); }); });
  });
  $("dl-multiverse-csv").addEventListener("click", function () {
    if (!S.multiResults) { alert("Build the multiverse first (Specification Curve tab)."); return; }
    Exports.download("multiverse_" + S.focal + "_" + outcomeKey() + ".csv", Exports.multiverseCSV(S.multiResults), "text/csv");
  });
  $("dl-curve-svg").addEventListener("click", function () { if (!document.querySelector("#chart-speccurve svg")) { alert("Build the specification curve first."); return; } Exports.exportSVG("chart-speccurve", "spec-curve-" + S.focal + ".svg"); });
  $("dl-curve-png").addEventListener("click", function () { if (!document.querySelector("#chart-speccurve svg")) { alert("Build the specification curve first."); return; } Exports.exportPNG("chart-speccurve", "spec-curve-" + S.focal + ".png", 2); });
  $("dl-report").addEventListener("click", buildLabReport);
  Array.prototype.forEach.call(document.getElementsByName("focal"), function (el) { el.addEventListener("change", function (e) { if (e.target.checked) setFocal(e.target.value); }); });
  $("outcome").addEventListener("change", function (e) { S.outcome = e.target.value; S.multiResults = null; renderAll(); syncHash(); });
  $("seType").addEventListener("change", function (e) { S.seType = e.target.value; S.multiResults = null; renderAll(); syncHash(); });
  $("btn-share").addEventListener("click", function () {
    var url = location.origin + location.pathname + stateToHash(); var b = this;
    navigator.clipboard.writeText(url).then(function () { var o = b.textContent; b.textContent = "🔗 Link copied ✓"; setTimeout(function () { b.textContent = o; }, 1500); });
  });
  // click anywhere on a chart (or its ⤢ button) to enlarge; close on backdrop, ✕, or Esc
  document.addEventListener("click", function (e) { var ch = e.target.closest ? e.target.closest(".chart") : null; if (ch && ch.id && ch.querySelector("svg")) openFig(ch.id); });
  $("figClose").addEventListener("click", closeFig);
  $("figModal").addEventListener("click", function (e) { if (e.target === $("figModal")) closeFig(); });
  document.addEventListener("keydown", function (e) { if (e.key === "Escape") closeFig(); });
  $("btn-lesson").addEventListener("click", function () { showLesson(0); });
  $("lessonNext").addEventListener("click", function () { showLesson(S.lessonStep + 1); });
  $("lessonPrev").addEventListener("click", function () { showLesson(S.lessonStep - 1); });
  $("lessonClose").addEventListener("click", closeLesson);
  $("axisYearFE").addEventListener("change", function (e) { S.axes.yearFEVary = e.target.checked; updateComboCount(); });
  Array.prototype.forEach.call(document.getElementsByClassName("axisSample"), function (el) {
    el.addEventListener("change", function () {
      S.axes.sampleVary = Array.prototype.filter.call(document.getElementsByClassName("axisSample"), function (x) { return x.checked; }).map(function (x) { return x.value; });
      if (!S.axes.sampleVary.length) { el.checked = true; S.axes.sampleVary = [el.value]; }
      updateComboCount();
    });
  });
  $("btn-runmulti").addEventListener("click", runMultiverse);
  $("showEq").addEventListener("change", function () { renderEquation(ssSpec(S.ss.estimator), E.runSpec(DATA, ssSpec(S.ss.estimator), META).model, S.ss.estimator); });
  $("btn-bootstrap").addEventListener("click", function () {
    var btn = this, lab = btn.textContent; btn.disabled = true; btn.textContent = "Bootstrapping… (B=500)";
    var box = $("bootstrapResult"); box.className = "alert alert-info mt-2"; box.textContent = "Resampling 545 persons with replacement, 500 times…";
    setTimeout(function () {
      var r = E.bootstrapFocal(DATA, ssSpec(S.ss.estimator), META, 500);
      var pt = E.runSpec(DATA, ssSpec(S.ss.estimator), META);
      if (r.B > 0) {
        box.className = "alert alert-success mt-2";
        box.innerHTML = "<strong>Cluster bootstrap — " + ESTIMATOR_LABEL[S.ss.estimator] + ", " + FOCAL_SHORT[S.focal] + "</strong><br>point estimate <strong>" + fmt(pt.b, 3) + "</strong> · 95% bootstrap CI <strong>[" + fmt(r.lo, 3) + ", " + fmt(r.hi, 3) + "]</strong> · bootstrap SE " + fmt(r.se, 3) + "<br><span class='help-text'>From " + r.B + " draws that resample whole persons (preserving within-person correlation). No distributional assumption.</span>";
      } else { box.className = "alert alert-warning mt-2"; box.textContent = "Not identified for this estimator/spec."; }
      btn.disabled = false; btn.textContent = lab;
    }, 30);
  });
  $("btn-inference").addEventListener("click", function () {
    var btn = this, lab = btn.textContent; btn.disabled = true; btn.textContent = "Running…";
    var axes = buildAxes(), obs = E.curveStats(E.enumerateMultiverse(DATA, axes, META));
    var B = 100, medGE = 0, sigGE = 0, b = 0;
    var box = $("inferenceBox"); box.style.display = "block"; box.className = "alert alert-info";
    setViz("multi"); box.style.display = "block"; // ensure visible on the multiverse tab
    function render(done) {
      var pMed = (medGE + 1) / (b + 1), pSig = (sigGE + 1) / (b + 1);
      box.innerHTML = "<strong>🎲 Randomization inference on the whole curve</strong>" +
        "<div class='progress-bar'><div class='progress-bar-fill' style='width:" + Math.round(100 * b / B) + "%'></div></div>" +
        (done ? "" : "Permuting the focal variable under “no effect” and rebuilding the curve… " + b + " / " + B + " draws") +
        (done ? ("Observed curve: median estimate <strong>" + fmt(obs.median, 3) + "</strong>, " + obs.nSigDom + " of " + obs.n + " specifications significant in the dominant direction.<br>" +
          "p(|median| ≥ observed) = <strong>" + pMed.toFixed(3) + "</strong> &nbsp;·&nbsp; p(# significant ≥ observed) = <strong>" + pSig.toFixed(3) + "</strong><br>" +
          "<span class='help-text'>" + (Math.min(pMed, pSig) < 0.05 ? "The observed curve is hard to reconcile with “the focal variable has no effect.”" : "The curve is not clearly distinguishable from what “no effect” would produce.") + "</span>") : "");
      box.className = "alert " + (done ? (Math.min(pMed, pSig) < 0.05 ? "alert-success" : "alert-warning") : "alert-info");
    }
    function step() {
      var chunk = Math.min(4, B - b);
      for (var c = 0; c < chunk; c++) { var s = E.nullReplication(DATA, axes, META); if (Math.abs(s.median) >= Math.abs(obs.median)) medGE++; if (s.nSigDom >= obs.nSigDom) sigGE++; b++; }
      render(b >= B);
      if (b < B) setTimeout(step, 0); else { btn.disabled = false; btn.textContent = lab; box.scrollIntoView({ behavior: "smooth", block: "center" }); }
    }
    render(false); setTimeout(step, 20);
  });
  $("ssEstimator").addEventListener("change", function (e) { S.ss.estimator = e.target.value; renderSingle(); syncHash(); });
  $("ssYearFE").addEventListener("change", function (e) { S.ss.yearFE = e.target.checked; renderSingle(); syncHash(); });
  $("ssSample").addEventListener("change", function (e) { S.ss.sample = e.target.value; renderSingle(); syncHash(); });
  Array.prototype.forEach.call(document.getElementsByClassName("collapsible"), function (el) {
    el.addEventListener("click", function () { el.classList.toggle("open"); var c = el.nextElementSibling; if (c) c.classList.toggle("show"); });
  });
}

// ============== INIT ==============
function init() {
  // guard against a stale cached engine.js running with a newer app.js
  if (!E || typeof E.curveStats !== "function" || !E.ESTIMATORS_ALL) {
    $("loading").innerHTML = "<div class='alert alert-warning'>An older cached version of the app is loaded. Please hard-refresh this page (⌘⇧R / Ctrl-Shift-R) to update. If it persists, empty the cache and reload, or open in a private window.</div>";
    return;
  }
  fetch("data/wagepan.json").then(function (r) { return r.json(); }).then(function (rows) {
    DATA = rows; META = computeMeta(rows);
    $("loading").style.display = "none"; $("appBody").style.display = "block";
    applyHash();
    buildAxisControls(); buildSSControls(); bindUI();
    syncDom(); updateComboCount();
    setViz(S.viz || "data");
    $("statusBar").innerHTML = "Loaded " + DATA.length + " person-years · estimators verified against R’s plm and Python’s linearmodels (see Methodology → Verification).";
  }).catch(function (err) {
    $("loading").innerHTML = "<div class='alert alert-warning'>Could not load data/wagepan.json — " + err + "</div>";
  });
}
init();
