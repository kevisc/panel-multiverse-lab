"use strict";
/*
 * engine.js — Panel Data Multiverse Lab numeric core.
 *
 * Pure computation, no DOM. Loadable in the browser (attaches to window.PanelEngine)
 * and in Node (module.exports) so the same code that runs the app also runs the
 * verification suite against canonical plm estimates (see verify.cjs).
 *
 * Estimators implemented: Pooled OLS, Between, Random Effects (Swamy–Arora style),
 * Fixed Effects (within), First Differences. Standard errors: classical, robust
 * (HC1), cluster-robust by person (CR1).
 */
(function (root, factory) {
  var mod = factory();
  if (typeof module === "object" && module.exports) module.exports = mod;
  if (root) root.PanelEngine = mod;
})(typeof self !== "undefined" ? self : this, function () {

  // ---------- basic stats ----------
  function mean(a) { return a.length ? a.reduce(function (x, y) { return x + y; }, 0) / a.length : NaN; }
  function variance(a) {
    var n = a.length; if (n <= 1) return 0;
    var m = mean(a), s = 0;
    for (var i = 0; i < n; i++) { var v = a[i] - m; s += v * v; }
    return s / (n - 1);
  }
  function quantile(a, q) {
    if (!a.length) return NaN;
    var s = a.slice().sort(function (x, y) { return x - y; });
    var p = (s.length - 1) * q, lo = Math.floor(p), hi = Math.ceil(p);
    if (lo === hi) return s[lo];
    return s[lo] * (hi - p) + s[hi] * (p - lo);
  }
  // mulberry32 — small, deterministic, seedable PRNG used by the cluster bootstrap
  // and the cluster-level permutation test so that p-values and bootstrap CIs are
  // bit-reproducible from a published seed.
  function mulberry32(seed) {
    var a = (seed | 0) >>> 0;
    return function () {
      a = (a + 0x6d2b79f5) >>> 0;
      var t = a;
      t = Math.imul(t ^ (t >>> 15), t | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }
  function pnorm(z) {
    var t = 1 / (1 + 0.2316419 * Math.abs(z));
    var d = 0.3989423 * Math.exp(-z * z / 2);
    var p = 1 - d * t * (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274))));
    return z > 0 ? p : 1 - p;
  }
  function pchisq(x, df) {
    // upper-tail p-value via Wilson–Hilferty
    if (x <= 0 || df <= 0) return 1;
    var z = Math.pow(x / df, 1 / 3) - (1 - 2 / (9 * df));
    z /= Math.sqrt(2 / (9 * df));
    return 1 - pnorm(z);
  }

  // ---------- linear algebra ----------
  function transpose(A) {
    var n = A.length, k = A[0].length, T = new Array(k);
    for (var j = 0; j < k; j++) { var row = new Array(n); for (var i = 0; i < n; i++) row[i] = A[i][j]; T[j] = row; }
    return T;
  }
  function matMul(A, B) {
    var n = A.length, k = A[0].length, m = B[0].length, out = new Array(n);
    for (var i = 0; i < n; i++) {
      var row = new Array(m).fill(0);
      for (var t = 0; t < k; t++) { var a = A[i][t]; if (a === 0) continue; for (var j = 0; j < m; j++) row[j] += a * B[t][j]; }
      out[i] = row;
    }
    return out;
  }
  function matVec(A, v) {
    var n = A.length, k = A[0].length, out = new Array(n);
    for (var i = 0; i < n; i++) { var s = 0; for (var j = 0; j < k; j++) s += A[i][j] * v[j]; out[i] = s; }
    return out;
  }
  function scale(A, s) { return A.map(function (r) { return r.map(function (v) { return v * s; }); }); }
  function diag(A) { var d = new Array(A.length); for (var i = 0; i < A.length; i++) d[i] = A[i][i]; return d; }
  function inverse(A) {
    var n = A.length, M = A.map(function (r) { return r.slice(); }), I = new Array(n);
    for (var i = 0; i < n; i++) { var r = new Array(n).fill(0); r[i] = 1; I[i] = r; }
    for (var col = 0; col < n; col++) {
      var piv = col, best = Math.abs(M[col][col]);
      for (var rr = col + 1; rr < n; rr++) { var v = Math.abs(M[rr][col]); if (v > best) { best = v; piv = rr; } }
      if (!isFinite(best) || best < 1e-12) return null;
      if (piv !== col) { var tm = M[col]; M[col] = M[piv]; M[piv] = tm; tm = I[col]; I[col] = I[piv]; I[piv] = tm; }
      var pv = M[col][col];
      for (var j = 0; j < n; j++) { M[col][j] /= pv; I[col][j] /= pv; }
      for (rr = 0; rr < n; rr++) {
        if (rr === col) continue;
        var f = M[rr][col]; if (f === 0) continue;
        for (j = 0; j < n; j++) { M[rr][j] -= f * M[col][j]; I[rr][j] -= f * I[col][j]; }
      }
    }
    return I;
  }
  function dotCol(col, y) { var s = 0; for (var i = 0; i < col.length; i++) s += col[i] * y[i]; return s; }

  // ---------- core OLS with classical / robust / cluster SEs ----------
  function ols(X, y, opts) {
    opts = opts || {};
    var names = opts.coefNames || [];
    var n = y.length;
    var fail = function () {
      var k0 = X && X[0] ? X[0].length : 0;
      return { beta: new Array(k0).fill(NaN), se: new Array(k0).fill(NaN), tstat: new Array(k0).fill(NaN),
        pval: new Array(k0).fill(NaN), fitted: [], resid: [], r2: NaN, rss: NaN, sigma2: NaN, vcov: [],
        n: n, k: k0, df: 0, seType: opts.seType || "classical", coefNames: names, ok: false };
    };
    if (!n || !X[0] || !X[0].length) return fail();
    var k = X[0].length;
    if (n <= k) return fail();
    var Xt = transpose(X), XtX = matMul(Xt, X), XtXinv = inverse(XtX);
    if (!XtXinv) return fail();
    var Xty = new Array(k);
    for (var j = 0; j < k; j++) Xty[j] = dotCol(Xt[j], y);
    var beta = matVec(XtXinv, Xty);
    var fitted = matVec(X, beta), resid = new Array(n), rss = 0, my = mean(y), tss = 0;
    for (var i = 0; i < n; i++) { var r = y[i] - fitted[i]; resid[i] = r; rss += r * r; var dy = y[i] - my; tss += dy * dy; }
    var df = Math.max(1, (opts.df != null ? opts.df : n - k));
    var sigma2 = rss / df, seType = opts.seType || "classical", vcov;
    if (seType === "robust") {
      var meat = []; for (i = 0; i < k; i++) meat.push(new Array(k).fill(0));
      for (i = 0; i < n; i++) { var u2 = resid[i] * resid[i]; for (var a = 0; a < k; a++) { var xa = X[i][a]; for (var b = 0; b < k; b++) meat[a][b] += u2 * xa * X[i][b]; } }
      vcov = scale(matMul(matMul(XtXinv, meat), XtXinv), n / df);
    } else if (seType === "cluster") {
      var ids = opts.clusterIds || [], by = {};
      for (i = 0; i < n; i++) { var c = String(ids[i]); (by[c] || (by[c] = [])).push(i); }
      var keys = Object.keys(by), G = keys.length;
      if (G < 2) { vcov = []; for (i = 0; i < k; i++) vcov.push(new Array(k).fill(NaN)); }
      else {
        var meatC = []; for (i = 0; i < k; i++) meatC.push(new Array(k).fill(0));
        for (var g = 0; g < G; g++) {
          var idx = by[keys[g]], sv = new Array(k).fill(0);
          for (var ii = 0; ii < idx.length; ii++) { var ri = idx[ii], ui = resid[ri]; for (j = 0; j < k; j++) sv[j] += X[ri][j] * ui; }
          for (var a2 = 0; a2 < k; a2++) for (var b2 = 0; b2 < k; b2++) meatC[a2][b2] += sv[a2] * sv[b2];
        }
        var adj = (G / (G - 1)) * ((n - 1) / df);
        vcov = scale(matMul(matMul(XtXinv, meatC), XtXinv), adj);
      }
    } else { seType = "classical"; vcov = scale(XtXinv, sigma2); }
    var se = diag(vcov).map(function (v) { return v >= 0 ? Math.sqrt(v) : NaN; });
    var tstat = new Array(k), pval = new Array(k);
    for (j = 0; j < k; j++) { tstat[j] = se[j] > 0 ? beta[j] / se[j] : NaN; pval[j] = isFinite(tstat[j]) ? 2 * (1 - pnorm(Math.abs(tstat[j]))) : NaN; }
    return { beta: beta, se: se, tstat: tstat, pval: pval, fitted: fitted, resid: resid,
      r2: tss > 0 ? 1 - rss / tss : NaN, rss: rss, sigma2: sigma2, vcov: vcov, n: n, k: k, df: df,
      seType: seType, coefNames: names, ok: true };
  }

  // ---------- grouping helpers ----------
  function groupIndex(ids) {
    var by = {}, order = [];
    for (var i = 0; i < ids.length; i++) { var id = String(ids[i]); if (!by[id]) { by[id] = []; order.push(id); } by[id].push(i); }
    return { by: by, order: order };
  }
  // per-person means broadcast back to each row (the Mundlak device for correlated random effects)
  function groupMeans(ids, M) {
    var grp = groupIndex(ids), n = M.length, k = M[0].length, out = new Array(n);
    grp.order.forEach(function (id) {
      var idx = grp.by[id], m = new Array(k).fill(0);
      idx.forEach(function (ii) { for (var j = 0; j < k; j++) m[j] += M[ii][j]; });
      for (var j = 0; j < k; j++) m[j] /= idx.length;
      idx.forEach(function (ii) { out[ii] = m.slice(); });
    });
    return out;
  }
  // Reduce a design to linearly independent columns, processed left-to-right via modified
  // Gram–Schmidt. Column order matters: earlier columns (intercept, then focal) are kept,
  // and any later column that is a linear combination of earlier ones is dropped — exactly
  // how plm/lm silently drop collinear regressors (e.g. year dummies that duplicate a linear
  // experience term within person, or constant columns in a between regression).
  function fullRank(X, names) {
    var n = X.length, k = X[0].length, keep = [], basis = [];
    for (var c = 0; c < k; c++) {
      var v = new Array(n), norm0 = 0;
      for (var i = 0; i < n; i++) { v[i] = X[i][c]; norm0 += v[i] * v[i]; }
      norm0 = Math.sqrt(norm0);
      var w = v.slice();
      for (var b = 0; b < basis.length; b++) { var d = 0; for (i = 0; i < n; i++) d += w[i] * basis[b][i]; for (i = 0; i < n; i++) w[i] -= d * basis[b][i]; }
      var norm = 0; for (i = 0; i < n; i++) norm += w[i] * w[i]; norm = Math.sqrt(norm);
      if (norm > 1e-7 * Math.max(1e-12, norm0)) { for (i = 0; i < n; i++) w[i] /= norm; basis.push(w); keep.push(c); }
    }
    if (keep.length === k) return { X: X, names: names, reduced: false };
    return { X: X.map(function (row) { return keep.map(function (c) { return row[c]; }); }), names: keep.map(function (c) { return names[c]; }), reduced: true };
  }
  // OLS after dropping collinear columns (order-preserving)
  function fitReduced(X, names, y, opts) {
    var rr = fullRank(X, names);
    var o = Object.assign({}, opts || {}); o.coefNames = rr.names;
    return ols(rr.X, y, o);
  }

  // ---------- Fixed Effects (within) ----------
  function feWithin(ids, X, y, opts) {
    opts = opts || {};
    var n = y.length, k = X[0].length, grp = groupIndex(ids), nG = grp.order.length;
    var Xw = new Array(n), Yw = new Array(n);
    grp.order.forEach(function (id) {
      var idx = grp.by[id], my = 0, mx = new Array(k).fill(0);
      for (var i = 0; i < idx.length; i++) { var ii = idx[i]; my += y[ii]; for (var j = 0; j < k; j++) mx[j] += X[ii][j]; }
      my /= idx.length; for (var j2 = 0; j2 < k; j2++) mx[j2] /= idx.length;
      for (i = 0; i < idx.length; i++) { ii = idx[i]; Yw[ii] = y[ii] - my; var row = new Array(k); for (j = 0; j < k; j++) row[j] = X[ii][j] - mx[j]; Xw[ii] = row; }
    });
    var rr = fullRank(Xw, opts.coefNames || []);
    var df = n - nG - rr.names.length;
    var fit = ols(rr.X, Yw, { seType: opts.seType, clusterIds: ids, df: df, coefNames: rr.names });
    fit.nGroups = nG;
    return fit;
  }

  // ---------- Between (person means) ----------
  function between(ids, X, y, opts) {
    opts = opts || {};
    var k = X[0].length, grp = groupIndex(ids), N = grp.order.length;
    var Xb = new Array(N), yb = new Array(N);
    grp.order.forEach(function (id, gi) {
      var idx = grp.by[id], my = 0, mx = new Array(k).fill(0);
      for (var i = 0; i < idx.length; i++) { var ii = idx[i]; my += y[ii]; for (var j = 0; j < k; j++) mx[j] += X[ii][j]; }
      my /= idx.length; for (var j2 = 0; j2 < k; j2++) mx[j2] /= idx.length;
      Xb[gi] = mx; yb[gi] = my;
    });
    var fit = fitReduced(Xb, opts.coefNames || [], yb, { seType: "classical" });
    fit.nGroups = N;
    return fit;
  }

  // ---------- Random Effects (Swamy–Arora style) ----------
  function randomEffects(ids, X, y, opts) {
    opts = opts || {};
    var n = y.length, k = X[0].length, names = opts.coefNames || [];
    var grp = groupIndex(ids), N = grp.order.length;
    var Ti = {}, ybar = {}, xbar = {};
    grp.order.forEach(function (id) {
      var idx = grp.by[id]; Ti[id] = idx.length;
      var my = 0, mx = new Array(k).fill(0);
      for (var i = 0; i < idx.length; i++) { var ii = idx[i]; my += y[ii]; for (var j = 0; j < k; j++) mx[j] += X[ii][j]; }
      my /= idx.length; for (var j2 = 0; j2 < k; j2++) mx[j2] /= idx.length;
      ybar[id] = my; xbar[id] = mx;
    });
    // sigma_e^2 from within regression on time-varying regressors (drop intercept col 0)
    var XnoInt = X.map(function (r) { return r.slice(1); });
    var fe = feWithin(ids, XnoInt, y, { seType: "classical", coefNames: names.slice(1) });
    var dfw = Math.max(1, fe.df || (n - N - (k - 1)));
    var sigma_e2 = (fe.rss || 0) / dfw; if (!isFinite(sigma_e2) || sigma_e2 <= 0) sigma_e2 = 1e-6;
    // sigma_u^2 from between regression residual variance
    var Xb = new Array(N), yb = new Array(N);
    grp.order.forEach(function (id, gi) { Xb[gi] = xbar[id].slice(); yb[gi] = ybar[id]; });
    // Between regression for the variance components: keep the intercept (a constant column
    // has zero variance but is not collinear) and let fullRank drop only redundant columns.
    var btw = fitReduced(Xb, names, yb, { seType: "classical" });
    var Tbar = mean(Object.keys(Ti).map(function (key) { return Ti[key]; }));
    var s2b = (btw.rss || 0) / Math.max(1, btw.df);
    var sigma_u2 = Math.max(0, s2b - sigma_e2 / Math.max(1e-9, Tbar || 1));
    // quasi-demeaning
    var Xq = new Array(n), Yq = new Array(n);
    var thetaBy = {};
    grp.order.forEach(function (id) { thetaBy[id] = 1 - Math.sqrt(sigma_e2 / Math.max(1e-9, sigma_e2 + Ti[id] * sigma_u2)); });
    for (var i = 0; i < n; i++) {
      var id = String(ids[i]), th = thetaBy[id] || 0, xb = xbar[id], rowQ = new Array(k);
      for (var j = 0; j < k; j++) rowQ[j] = X[i][j] - th * xb[j];
      Xq[i] = rowQ; Yq[i] = y[i] - th * ybar[id];
    }
    var fit = fitReduced(Xq, names, Yq, { seType: opts.seType, clusterIds: ids });
    fit.theta = 1 - Math.sqrt(sigma_e2 / Math.max(1e-9, sigma_e2 + (Tbar || 1) * sigma_u2));
    fit.sigma_e2 = sigma_e2; fit.sigma_u2 = sigma_u2; fit.nGroups = N;
    return fit;
  }

  // ---------- First Differences ----------
  function firstDifferences(ids, year, X, y, opts) {
    opts = opts || {};
    var k = X[0].length, grp = groupIndex(ids);
    var Xd = [], yd = [], idd = [];
    grp.order.forEach(function (id) {
      var idx = grp.by[id].slice().sort(function (a, b) { return year[a] - year[b]; });
      for (var t = 1; t < idx.length; t++) {
        var cur = idx[t], prev = idx[t - 1], row = new Array(k);
        for (var j = 0; j < k; j++) row[j] = X[cur][j] - X[prev][j];
        Xd.push(row); yd.push(y[cur] - y[prev]); idd.push(id);
      }
    });
    // include intercept (common drift) as first column; drop collinear differenced columns
    // (e.g. Δexper is ~constant within person and collinear with the drift term).
    var Xfull = Xd.map(function (r) { return [1].concat(r); });
    var names = ["(drift)"].concat(opts.coefNames || []);
    var fit = fitReduced(Xfull, names, yd, { seType: opts.seType, clusterIds: idd });
    fit.nGroups = grp.order.length;
    return fit;
  }

  // ---------- generalized Hausman test on a named coefficient ----------
  function coefAt(model, name) {
    var i = model && model.coefNames ? model.coefNames.indexOf(name) : -1;
    if (i < 0) return { b: NaN, se: NaN, varc: NaN, idx: -1 };
    var v = model.vcov && model.vcov[i] ? model.vcov[i][i] : NaN;
    return { b: model.beta[i], se: model.se[i], varc: v, idx: i };
  }
  // ============================================================
  //  Design construction for the wage panel
  // ============================================================
  var FOCAL = { union: "union", married: "married", educ: "educ" };

  // sample filter
  function applySample(rows, sample, outcomeKey) {
    if (sample === "nohealth") return rows.filter(function (r) { return r.poorhlth !== 1; });
    if (sample === "trim") {
      var v = rows.map(function (r) { return r[outcomeKey]; });
      var lo = quantile(v, 0.01), hi = quantile(v, 0.99);
      return rows.filter(function (r) { return r[outcomeKey] >= lo && r[outcomeKey] <= hi; });
    }
    return rows;
  }

  // build the raw (no-intercept) design columns from a spec
  function buildRaw(rows, spec, meta) {
    var focal = spec.focal, ctrl = spec.controls || {};
    var names = [], cols = [];
    function addCol(name, fn) { names.push(name); cols.push(rows.map(fn)); }
    addCol(focal, function (r) { return r[focal]; });
    if (ctrl.experience) { addCol("exper", function (r) { return r.exper; }); addCol("expersq", function (r) { return r.expersq; }); }
    if (ctrl.hours) addCol("hours", function (r) { return r.hours; });
    if (ctrl.health) addCol("poorhlth", function (r) { return r.poorhlth; });
    if (ctrl.race && focal !== "educ") { addCol("black", function (r) { return r.black; }); addCol("hisp", function (r) { return r.hisp; }); }
    if (ctrl.industry) meta.industryCats.slice(1).forEach(function (cat) { addCol("ind:" + cat, function (r) { return r.industry === cat ? 1 : 0; }); });
    if (ctrl.region) meta.regionCats.slice(1).forEach(function (cat) { addCol("reg:" + cat, function (r) { return r.region === cat ? 1 : 0; }); });
    if (spec.yearFE) meta.years.slice(1).forEach(function (yr) { addCol("yr:" + yr, function (r) { return r.year === yr ? 1 : 0; }); });
    // transpose cols -> M (n × p)
    var n = rows.length, p = names.length, M = new Array(n);
    for (var i = 0; i < n; i++) { var row = new Array(p); for (var j = 0; j < p; j++) row[j] = cols[j][i]; M[i] = row; }
    return { names: names, M: M, ids: rows.map(function (r) { return r.nr; }),
      year: rows.map(function (r) { return r.year; }), y: rows.map(function (r) { return r[spec.outcomeKey]; }) };
  }

  // run one specification → focal estimate
  function runSpec(allRows, spec, meta) {
    spec.outcomeKey = spec.outcome === "wage" ? "wage" : "lwage";
    var rows = applySample(allRows, spec.sample, spec.outcomeKey);
    var raw = buildRaw(rows, spec, meta);
    var focal = spec.focal, est = spec.estimator, model = null, focalName = focal;
    var withInt = raw.M.map(function (r) { return [1].concat(r); });
    var namesInt = ["(Intercept)"].concat(raw.names);

    if (est === "pooled") {
      model = fitReduced(withInt, namesInt, raw.y, { seType: spec.seType, clusterIds: raw.ids });
    } else if (est === "between") {
      model = between(raw.ids, withInt, raw.y, { coefNames: namesInt });
    } else if (est === "random") {
      model = randomEffects(raw.ids, withInt, raw.y, { seType: spec.seType, coefNames: namesInt });
    } else if (est === "within") {
      model = feWithin(raw.ids, raw.M, raw.y, { seType: spec.seType, coefNames: raw.names });
    } else if (est === "fd") {
      model = firstDifferences(raw.ids, raw.year, raw.M, raw.y, { seType: spec.seType, coefNames: raw.names });
    } else if (est === "cre") {
      // Correlated random effects (Mundlak): RE augmented with person-means of the regressors.
      // The coefficient on the original (time-varying) regressor then equals the fixed-effects one.
      var gm = groupMeans(raw.ids, raw.M);
      var aug = raw.M.map(function (r, i) { return [1].concat(r).concat(gm[i]); });
      var augNames = namesInt.concat(raw.names.map(function (nm) { return "mean(" + nm + ")"; }));
      model = randomEffects(raw.ids, aug, raw.y, { seType: spec.seType, coefNames: augNames });
    } else if (est === "twfe") {
      // two-way (entity + time) fixed effects = within with year dummies forced in
      var rawT = spec.yearFE ? raw : buildRaw(rows, Object.assign({}, spec, { yearFE: true, outcomeKey: spec.outcomeKey }), meta);
      model = feWithin(rawT.ids, rawT.M, rawT.y, { seType: spec.seType, coefNames: rawT.names });
    }
    var c = coefAt(model, focalName);
    var identified = isFinite(c.b);
    return { spec: spec, model: model, focalName: focalName, b: c.b, se: c.se, identified: identified,
      estimator: est, n: rows.length, nGroups: model ? model.nGroups : 0 };
  }

  // enumerate the multiverse over selected axes
  var ALL_ESTIMATORS = ["pooled", "between", "random", "within", "fd"]; // default multiverse axis
  var ESTIMATORS_ALL = ["pooled", "between", "random", "within", "fd", "cre", "twfe"]; // all selectable
  var CONTROL_KEYS = ["experience", "hours", "industry", "region", "health"]; // race handled separately
  function enumerateMultiverse(allRows, axes, meta) {
    // axes: { estimators:[...], controlsVary:[keys that toggle], controlsFixed:{}, yearFEVary:bool, yearFEFixed:bool,
    //         sampleVary:[...], focal, outcome, seType }
    var results = [];
    var ctrlKeys = axes.controlsVary || [];
    var nCombo = 1 << ctrlKeys.length;
    var yearOpts = axes.yearFEVary ? [false, true] : [axes.yearFEFixed];
    var sampleOpts = axes.sampleVary && axes.sampleVary.length ? axes.sampleVary : ["full"];
    axes.estimators.forEach(function (est) {
      for (var mask = 0; mask < nCombo; mask++) {
        var controls = Object.assign({}, axes.controlsFixed || {});
        ctrlKeys.forEach(function (key, bit) { controls[key] = !!(mask & (1 << bit)); });
        yearOpts.forEach(function (yfe) {
          sampleOpts.forEach(function (samp) {
            var spec = { focal: axes.focal, outcome: axes.outcome, controls: controls, yearFE: yfe,
              sample: samp, estimator: est, seType: axes.seType };
            results.push(runSpec(allRows, spec, meta));
          });
        });
      }
    });
    return results;
  }

  // ---------- full (multi-coefficient) Hausman test ----------
  // Classic Hausman over the slope coefficients common to FE and RE.
  function hausmanFull(fe, re) {
    var names = (fe.coefNames || []).filter(function (nm) {
      return nm !== "(Intercept)" && nm.indexOf("mean(") !== 0 && (re.coefNames || []).indexOf(nm) >= 0;
    });
    var k = names.length;
    if (!k) return { stat: NaN, pval: NaN, df: 0, ok: false, msg: "no common coefficients" };
    var iFe = names.map(function (nm) { return fe.coefNames.indexOf(nm); });
    var iRe = names.map(function (nm) { return re.coefNames.indexOf(nm); });
    var diff = names.map(function (nm, i) { return fe.beta[iFe[i]] - re.beta[iRe[i]]; });
    var Vd = [];
    for (var a = 0; a < k; a++) { Vd.push([]); for (var b = 0; b < k; b++) Vd[a].push(fe.vcov[iFe[a]][iFe[b]] - re.vcov[iRe[a]][iRe[b]]); }
    var Vinv = inverse(Vd);
    if (!Vinv) return { stat: NaN, pval: NaN, df: k, ok: false, msg: "Var(FE)−Var(RE) is not positive definite for this spec" };
    var tmp = matVec(Vinv, diff), H = 0;
    for (var i2 = 0; i2 < k; i2++) H += diff[i2] * tmp[i2];
    return { stat: H, pval: pchisq(H, k), df: k, ok: isFinite(H) && H >= 0, msg: H < 0 ? "negative statistic (Var difference indefinite)" : null };
  }

  // ---------- cluster (block) bootstrap of the focal coefficient ----------
  // Resamples persons with replacement; each draw is a fresh cluster so blocks stay distinct.
  // opts.rng: optional () => Number in [0,1) — when provided (e.g. seeded mulberry32),
  // the bootstrap is bit-reproducible; defaults to Math.random for backwards compatibility.
  function bootstrapFocal(allRows, spec, meta, B, opts) {
    B = B || 500;
    opts = opts || {};
    var rand = opts.rng || Math.random;
    var byId = {}, ids = [];
    allRows.forEach(function (r) { if (!byId[r.nr]) { byId[r.nr] = []; ids.push(r.nr); } byId[r.nr].push(r); });
    var draws = [];
    for (var b = 0; b < B; b++) {
      var samp = [];
      for (var i = 0; i < ids.length; i++) {
        var id = ids[(rand() * ids.length) | 0], cl = "b" + i;
        byId[id].forEach(function (r) { var c = Object.assign({}, r); c.nr = cl; samp.push(c); });
      }
      var res = runSpec(samp, spec, meta);
      if (res.identified && isFinite(res.b)) draws.push(res.b);
    }
    draws.sort(function (a, b) { return a - b; });
    return { lo: quantile(draws, 0.025), hi: quantile(draws, 0.975), se: Math.sqrt(variance(draws)), B: draws.length };
  }

  // ---------- specification-curve (joint) inference ----------
  // Two-sided 95% normal critical value (used for plotted CIs and the "significant in the
  // dominant direction" count; matches plm's default summary intervals).
  var Z_95 = 1.96;
  // Test statistics on a curve: the median estimate and the count significant in the dominant direction.
  function curveStats(results) {
    var ok = results.filter(function (r) { return r.identified && isFinite(r.b) && isFinite(r.se); });
    if (!ok.length) return { median: NaN, nSigDom: 0, n: 0 };
    var med = quantile(ok.map(function (r) { return r.b; }), 0.5), dom = med >= 0 ? 1 : -1, nSig = 0;
    ok.forEach(function (r) { var lo = r.b - Z_95 * r.se, hi = r.b + Z_95 * r.se; var sig = lo > 0 || hi < 0; if (sig && (r.b >= 0 ? 1 : -1) === dom) nSig++; });
    return { median: med, nSigDom: nSig, n: ok.length, dom: dom };
  }
  // Cluster-level (person-level) permutation of the focal regressor under the sharp null.
  // We assemble each person's focal sequence (ordered by year for determinism), shuffle the
  // mapping of persons to those sequences, then write each donor's sequence into the recipient's
  // rows in the recipient's year order. This preserves within-person dependence (autocorrelation
  // of the treatment trajectory) — the principled null for panel data, rather than a flat
  // row-level permutation which would destroy that structure.
  // opts.rng: optional seeded PRNG; defaults to Math.random.
  function permuteFocal(rows, focal, opts) {
    opts = opts || {};
    var rand = opts.rng || Math.random;
    var grp = groupIndex(rows.map(function (r) { return r.nr; }));
    var nPersons = grp.order.length;
    // sort each person's row indices by year, so seq positions align across donors/recipients
    var personOrder = grp.order.map(function (id) {
      return grp.by[id].slice().sort(function (a, b) {
        return (rows[a].year || 0) - (rows[b].year || 0);
      });
    });
    var personSeq = personOrder.map(function (idxs) {
      return idxs.map(function (i) { return rows[i][focal]; });
    });
    // Fisher-Yates over person indices (seeded if opts.rng provided)
    var perm = new Array(nPersons); for (var k = 0; k < nPersons; k++) perm[k] = k;
    for (var i = nPersons - 1; i > 0; i--) {
      var j = (rand() * (i + 1)) | 0, t = perm[i]; perm[i] = perm[j]; perm[j] = t;
    }
    // write donor's k-th value into recipient's k-th row (modulo length for unbalanced panels)
    var newVals = new Array(rows.length);
    for (var r = 0; r < nPersons; r++) {
      var rIdx = personOrder[r], donor = personSeq[perm[r]];
      for (var p = 0; p < rIdx.length; p++) newVals[rIdx[p]] = donor[p % donor.length];
    }
    return rows.map(function (row, i) { var c = Object.assign({}, row); c[focal] = newVals[i]; return c; });
  }
  // one randomization replication under the sharp null (focal permuted at the person level)
  function nullReplication(allRows, axes, meta, opts) {
    return curveStats(enumerateMultiverse(permuteFocal(allRows, axes.focal, opts), axes, meta));
  }
  // summary of a multiverse run (over identified specs only)
  function summarize(results) {
    var ok = results.filter(function (r) { return r.identified && isFinite(r.b); });
    var bs = ok.map(function (r) { return r.b; });
    var sig = ok.filter(function (r) { var lo = r.b - Z_95 * r.se, hi = r.b + Z_95 * r.se; return lo > 0 || hi < 0; });
    var med = ok.length ? quantile(bs, 0.5) : NaN;
    var refSign = med >= 0 ? 1 : -1;
    var flips = ok.filter(function (r) { return (r.b >= 0 ? 1 : -1) !== refSign; });
    return {
      total: results.length, identified: ok.length, unidentified: results.length - ok.length,
      median: med, min: ok.length ? Math.min.apply(null, bs) : NaN, max: ok.length ? Math.max.apply(null, bs) : NaN,
      q25: ok.length ? quantile(bs, 0.25) : NaN, q75: ok.length ? quantile(bs, 0.75) : NaN,
      shareSig: ok.length ? sig.length / ok.length : NaN, shareFlip: ok.length ? flips.length / ok.length : NaN
    };
  }

  return {
    mean: mean, variance: variance, quantile: quantile, pnorm: pnorm, pchisq: pchisq, mulberry32: mulberry32,
    ols: ols, feWithin: feWithin, between: between, randomEffects: randomEffects,
    firstDifferences: firstDifferences, hausmanFull: hausmanFull, coefAt: coefAt,
    groupMeans: groupMeans, bootstrapFocal: bootstrapFocal,
    curveStats: curveStats, permuteFocal: permuteFocal, nullReplication: nullReplication,
    buildRaw: buildRaw, runSpec: runSpec, enumerateMultiverse: enumerateMultiverse,
    summarize: summarize, applySample: applySample,
    ALL_ESTIMATORS: ALL_ESTIMATORS, ESTIMATORS_ALL: ESTIMATORS_ALL, CONTROL_KEYS: CONTROL_KEYS, FOCAL: FOCAL
  };
});
