"use strict";

// ============== UTILITY FUNCTIONS ==============
function mulberry32(a) {
  return function() {
    var t = a += 0x6d2b79f5;
    t = Math.imul(t ^ t >>> 15, t | 1);
    t ^= t + Math.imul(t ^ t >>> 7, t | 61);
    return ((t ^ t >>> 14) >>> 0) / 4294967296;
  };
}

function logistic(z) { return 1 / (1 + Math.exp(-z)); }
function clamp(x, a, b) { return Math.max(a, Math.min(b, x)); }
function mean(a) { return a.length ? a.reduce((x, y) => x + y, 0) / a.length : NaN; }
function variance(a) {
  var n = a.length;
  if (n <= 1) return 0;
  var m = mean(a), s = 0;
  for (var i = 0; i < n; i++) { var v = a[i] - m; s += v * v; }
  return s / (n - 1);
}
function quantile(a, q) {
  if (!a.length) return NaN;
  var s = a.slice().sort((x, y) => x - y);
  var p = (s.length - 1) * q;
  var lo = Math.floor(p), hi = Math.ceil(p);
  if (lo === hi) return s[lo];
  return s[lo] * (hi - p) + s[hi] * (p - lo);
}
function rnorm(rng, m = 0, s = 1) {
  var u = 0, v = 0;
  while (u === 0) u = rng();
  while (v === 0) v = rng();
  var z = Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  return m + s * z;
}
function fmt(x, d = 3) {
  return (Number.isNaN(x) || x === void 0 || !isFinite(x)) ? "—" : Number(x).toFixed(d);
}
function ci(b, se) { return [b - 1.96 * (se || NaN), b + 1.96 * (se || NaN)]; }
function pchisq(x, df) {
  // Approximate chi-squared p-value using Wilson-Hilferty transformation
  if (x <= 0 || df <= 0) return 1;
  var z = Math.pow(x / df, 1/3) - (1 - 2 / (9 * df));
  z /= Math.sqrt(2 / (9 * df));
  // Standard normal CDF approximation
  var t = 1 / (1 + 0.2316419 * Math.abs(z));
  var d = 0.3989423 * Math.exp(-z * z / 2);
  var p = d * t * (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274))));
  return z > 0 ? p : 1 - p;
}

// ============== LINEAR ALGEBRA & REGRESSION ==============
function pnorm(z) {
  var t = 1 / (1 + 0.2316419 * Math.abs(z));
  var d = 0.3989423 * Math.exp(-z * z / 2);
  var p = 1 - d * t * (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274))));
  return z > 0 ? p : 1 - p;
}

function vecDot(a, b) {
  var s = 0;
  for (var i = 0; i < a.length; i++) s += a[i] * b[i];
  return s;
}
function matIdentity(n) {
  var I = new Array(n);
  for (var i = 0; i < n; i++) {
    var row = new Array(n).fill(0);
    row[i] = 1;
    I[i] = row;
  }
  return I;
}
function matClone(A) { return A.map(r => r.slice()); }
function matTranspose(A) {
  var n = A.length, k = A[0].length;
  var T = new Array(k);
  for (var j = 0; j < k; j++) {
    var row = new Array(n);
    for (var i = 0; i < n; i++) row[i] = A[i][j];
    T[j] = row;
  }
  return T;
}
function matMul(A, B) {
  var n = A.length, k = A[0].length, m = B[0].length;
  var out = new Array(n);
  for (var i = 0; i < n; i++) {
    var row = new Array(m).fill(0);
    for (var t = 0; t < k; t++) {
      var a = A[i][t];
      if (a === 0) continue;
      for (var j = 0; j < m; j++) row[j] += a * B[t][j];
    }
    out[i] = row;
  }
  return out;
}
function matVecMul(A, v) {
  var n = A.length, k = A[0].length;
  var out = new Array(n);
  for (var i = 0; i < n; i++) {
    var s = 0;
    for (var j = 0; j < k; j++) s += A[i][j] * v[j];
    out[i] = s;
  }
  return out;
}
function matScale(A, s) { return A.map(r => r.map(v => v * s)); }
function matDiag(A) {
  var d = new Array(A.length);
  for (var i = 0; i < A.length; i++) d[i] = A[i][i];
  return d;
}

function matInverse(A) {
  var n = A.length;
  var M = matClone(A);
  var I = matIdentity(n);
  for (var col = 0; col < n; col++) {
    var pivotRow = col;
    var best = Math.abs(M[col][col]);
    for (var r = col + 1; r < n; r++) {
      var v = Math.abs(M[r][col]);
      if (v > best) { best = v; pivotRow = r; }
    }
    if (!isFinite(best) || best < 1e-12) return null;
    if (pivotRow !== col) {
      var tmp = M[col]; M[col] = M[pivotRow]; M[pivotRow] = tmp;
      tmp = I[col]; I[col] = I[pivotRow]; I[pivotRow] = tmp;
    }
    var pivot = M[col][col];
    for (var j = 0; j < n; j++) { M[col][j] /= pivot; I[col][j] /= pivot; }
    for (r = 0; r < n; r++) {
      if (r === col) continue;
      var factor = M[r][col];
      if (factor === 0) continue;
      for (j = 0; j < n; j++) {
        M[r][j] -= factor * M[col][j];
        I[r][j] -= factor * I[col][j];
      }
    }
  }
  return I;
}

function ols(X, y, opts = {}) {
  var n = y.length;
  if (!n) return { beta: [], se: [], tstat: [], pval: [], fitted: [], resid: [], r2: NaN, rss: NaN, sigma2: NaN, vcov: [], n: 0, k: 0, df: 0, seType: opts.seType || "classical", coefNames: opts.coefNames || [] };
  var k = X[0].length;
  if (n <= k) return { beta: new Array(k).fill(NaN), se: new Array(k).fill(NaN), tstat: new Array(k).fill(NaN), pval: new Array(k).fill(NaN), fitted: [], resid: [], r2: NaN, rss: NaN, sigma2: NaN, vcov: [], n, k, df: 0, seType: opts.seType || "classical", coefNames: opts.coefNames || [] };

  var Xt = matTranspose(X);
  var XtX = matMul(Xt, X);
  var XtXinv = matInverse(XtX);
  if (!XtXinv) return { beta: new Array(k).fill(NaN), se: new Array(k).fill(NaN), tstat: new Array(k).fill(NaN), pval: new Array(k).fill(NaN), fitted: [], resid: [], r2: NaN, rss: NaN, sigma2: NaN, vcov: [], n, k, df: 0, seType: opts.seType || "classical", coefNames: opts.coefNames || [] };

  var Xty = new Array(k).fill(0);
  for (var j = 0; j < k; j++) Xty[j] = vecDot(Xt[j], y);
  var beta = matVecMul(XtXinv, Xty);

  var fitted = matVecMul(X, beta);
  var resid = new Array(n);
  var rss = 0;
  var my = mean(y);
  var tss = 0;
  for (var i = 0; i < n; i++) {
    var r = y[i] - fitted[i];
    resid[i] = r;
    rss += r * r;
    var dy = y[i] - my;
    tss += dy * dy;
  }
  var df = Math.max(1, (opts.df != null ? opts.df : (n - k)));
  var sigma2 = rss / df;

  var seType = opts.seType || "classical";
  var vcov;
  if (seType === "classical") {
    vcov = matScale(XtXinv, sigma2);
  } else if (seType === "robust") {
    var meat = new Array(k);
    for (i = 0; i < k; i++) meat[i] = new Array(k).fill(0);
    for (i = 0; i < n; i++) {
      var u = resid[i];
      var wi = u * u;
      for (var a = 0; a < k; a++) {
        var xa = X[i][a];
        for (var b = 0; b < k; b++) meat[a][b] += wi * xa * X[i][b];
      }
    }
    vcov = matMul(matMul(XtXinv, meat), XtXinv);
    vcov = matScale(vcov, n / df); // HC1 adjustment
  } else if (seType === "cluster") {
    var clusterIds = opts.clusterIds || [];
    var by = {};
    for (i = 0; i < n; i++) {
      var cid = String(clusterIds[i]);
      if (!by[cid]) by[cid] = [];
      by[cid].push(i);
    }
    var groups = Object.keys(by);
    var G = groups.length;
    if (G < 2) {
      vcov = new Array(k);
      for (i = 0; i < k; i++) vcov[i] = new Array(k).fill(NaN);
    } else {
      var meatC = new Array(k);
      for (i = 0; i < k; i++) meatC[i] = new Array(k).fill(0);
      for (var g = 0; g < G; g++) {
        var idxs = by[groups[g]];
        var svec = new Array(k).fill(0);
        for (var ii = 0; ii < idxs.length; ii++) {
          var rowIdx = idxs[ii];
          var ui = resid[rowIdx];
          for (j = 0; j < k; j++) svec[j] += X[rowIdx][j] * ui;
        }
        for (var a2 = 0; a2 < k; a2++) for (var b2 = 0; b2 < k; b2++) meatC[a2][b2] += svec[a2] * svec[b2];
      }
      vcov = matMul(matMul(XtXinv, meatC), XtXinv);
      var adj = (G / (G - 1)) * ((n - 1) / df);
      vcov = matScale(vcov, adj);
    }
  } else {
    vcov = matScale(XtXinv, sigma2);
    seType = "classical";
  }

  var se = matDiag(vcov).map(v => (v >= 0 ? Math.sqrt(v) : NaN));
  var tstat = new Array(k), pval = new Array(k);
  for (j = 0; j < k; j++) {
    tstat[j] = se[j] > 0 ? beta[j] / se[j] : NaN;
    pval[j] = isFinite(tstat[j]) ? 2 * (1 - pnorm(Math.abs(tstat[j]))) : NaN;
  }

  var r2 = tss > 0 ? 1 - rss / tss : NaN;
  return { beta, se, tstat, pval, fitted, resid, r2, rss, sigma2, vcov, n, k, df, seType, coefNames: opts.coefNames || [] };
}

function feWithin(ids, X, y, opts = {}) {
  var n = y.length;
  if (!n) return { beta: [], se: [], tstat: [], pval: [], fitted: [], resid: [], r2: NaN, rss: NaN, sigma2: NaN, vcov: [], n: 0, k: 0, df: 0, seType: opts.seType || "classical", coefNames: opts.coefNames || [], nGroups: 0, Xw: [], Yw: [] };

  var k = X[0].length;
  var idxBy = {};
  for (var i = 0; i < n; i++) {
    var id = String(ids[i]);
    if (!idxBy[id]) idxBy[id] = [];
    idxBy[id].push(i);
  }
  var groups = Object.keys(idxBy);
  var nGroups = groups.length;

  var Xw = new Array(n);
  var Yw = new Array(n);
  for (var g = 0; g < nGroups; g++) {
    var idxs = idxBy[groups[g]];
    var my = 0;
    var mx = new Array(k).fill(0);
    for (i = 0; i < idxs.length; i++) {
      var ii = idxs[i];
      my += y[ii];
      for (var j = 0; j < k; j++) mx[j] += X[ii][j];
    }
    my /= idxs.length;
    for (j = 0; j < k; j++) mx[j] /= idxs.length;
    for (i = 0; i < idxs.length; i++) {
      ii = idxs[i];
      Yw[ii] = y[ii] - my;
      var row = new Array(k);
      for (j = 0; j < k; j++) row[j] = X[ii][j] - mx[j];
      Xw[ii] = row;
    }
  }

  var df = n - nGroups - k;
  var fit = ols(Xw, Yw, { seType: opts.seType, clusterIds: ids, df, coefNames: opts.coefNames || [] });
  fit.nGroups = nGroups;
  fit.Xw = Xw;
  fit.Yw = Yw;
  return fit;
}

function columnVariance(X, col) {
  var n = X.length;
  var m = 0;
  for (var i = 0; i < n; i++) m += X[i][col];
  m /= n;
  var s = 0;
  for (i = 0; i < n; i++) { var d = X[i][col] - m; s += d * d; }
  return s / Math.max(1, n - 1);
}

function randomEffects(ids, X, y, opts = {}) {
  var n = y.length;
  if (!n) return { beta: [], se: [], tstat: [], pval: [], fitted: [], resid: [], r2: NaN, rss: NaN, sigma2: NaN, vcov: [], n: 0, k: 0, df: 0, seType: opts.seType || "classical", coefNames: opts.coefNames || [], nGroups: 0, theta: NaN, sigma_e2: NaN, sigma_u2: NaN, Xq: [], Yq: [] };

  var k = X[0].length;
  var idxBy = {};
  for (var i = 0; i < n; i++) {
    var id = String(ids[i]);
    if (!idxBy[id]) idxBy[id] = [];
    idxBy[id].push(i);
  }
  var groups = Object.keys(idxBy);
  var N = groups.length;

  var Ti = {};
  var ybar = {};
  var xbar = {};
  for (var g = 0; g < N; g++) {
    var gid = groups[g];
    var idxs = idxBy[gid];
    Ti[gid] = idxs.length;
    var my = 0;
    var mx = new Array(k).fill(0);
    for (i = 0; i < idxs.length; i++) {
      var ii = idxs[i];
      my += y[ii];
      for (var j = 0; j < k; j++) mx[j] += X[ii][j];
    }
    my /= idxs.length;
    for (j = 0; j < k; j++) mx[j] /= idxs.length;
    ybar[gid] = my;
    xbar[gid] = mx;
  }

  // Estimate sigma_e^2 from within regression on time-varying regressors (drop intercept)
  var XnoInt = X.map(r => r.slice(1));
  var coefNamesNoInt = (opts.coefNames || []).slice(1);
  var feFit = feWithin(ids, XnoInt, y, { seType: "classical", coefNames: coefNamesNoInt });
  var dfw = Math.max(1, feFit.df || (n - N - (k - 1)));
  var sigma_e2 = (feFit.rss || 0) / dfw;
  if (!isFinite(sigma_e2) || sigma_e2 <= 0) sigma_e2 = 1e-6;

  // Between regression to estimate sigma_u^2 (drop columns with ~0 variance across groups)
  var Xb = new Array(N);
  var yb = new Array(N);
  for (g = 0; g < N; g++) {
    gid = groups[g];
    Xb[g] = xbar[gid].slice();
    yb[g] = ybar[gid];
  }
  var keep = [];
  for (var col = 0; col < k; col++) {
    var vcol = columnVariance(Xb, col);
    if (vcol > 1e-10) keep.push(col);
  }
  if (keep.length === 0) keep = [0];
  var XbR = Xb.map(row => keep.map(c => row[c]));
  var coefNamesB = keep.map(c => (opts.coefNames || [])[c] || ("x" + c));
  var between = ols(XbR, yb, { seType: "classical", coefNames: coefNamesB });
  var s2_between = (between.rss || 0) / Math.max(1, between.df || (N - XbR[0].length));
  var Tbar = mean(Object.values(Ti));
  var sigma_u2 = Math.max(0, s2_between - sigma_e2 / Math.max(1e-9, Tbar || 1));

  // Quasi-demeaning
  var Xq = new Array(n);
  var Yq = new Array(n);
  var thetaBy = {};
  for (g = 0; g < N; g++) {
    gid = groups[g];
    var th = 1 - Math.sqrt(sigma_e2 / Math.max(1e-9, sigma_e2 + Ti[gid] * sigma_u2));
    thetaBy[gid] = th;
  }
  for (i = 0; i < n; i++) {
    gid = String(ids[i]);
    var th2 = thetaBy[gid] || 0;
    var xbRow = xbar[gid];
    var rowQ = new Array(k);
    for (j = 0; j < k; j++) rowQ[j] = X[i][j] - th2 * xbRow[j];
    Xq[i] = rowQ;
    Yq[i] = y[i] - th2 * ybar[gid];
  }

  var fit = ols(Xq, Yq, { seType: opts.seType, clusterIds: ids, coefNames: opts.coefNames || [] });
  fit.theta = 1 - Math.sqrt(sigma_e2 / Math.max(1e-9, sigma_e2 + (Tbar || 1) * sigma_u2));
  fit.sigma_e2 = sigma_e2;
  fit.sigma_u2 = sigma_u2;
  fit.nGroups = N;
  fit.Xq = Xq;
  fit.Yq = Yq;
  return fit;
}

// ============== HAUSMAN TEST ==============
function coefIndex(model, name) {
  if (!model || !model.coefNames) return -1;
  return model.coefNames.indexOf(name);
}
function coefAt(model, name) {
  var idx = coefIndex(model, name);
  if (idx < 0) return { b: NaN, se: NaN, var: NaN, idx };
  var v = model.vcov && model.vcov[idx] ? model.vcov[idx][idx] : NaN;
  return { b: model.beta[idx], se: model.se[idx], var: v, idx };
}
function hausmanTest(fe, re) {
  var feD = coefAt(fe, "D");
  var reD = coefAt(re, "D");
  if (!isFinite(feD.b) || !isFinite(reD.b) || !isFinite(feD.var) || !isFinite(reD.var)) {
    return { stat: NaN, pval: NaN, conclusion: "Cannot compute" };
  }
  var diff = feD.b - reD.b;
  var varDiff = feD.var - reD.var;
  if (varDiff <= 0) {
    return { stat: NaN, pval: NaN, conclusion: "Variance difference non-positive; test unreliable" };
  }
  var H = diff * diff / varDiff;
  var pval = pchisq(H, 1);
  var conclusion = pval < 0.05
    ? "Reject H₀ (p < 0.05): FE preferred (evidence of correlation between uᵢ and D)"
    : "Cannot reject H₀: RE may be more efficient (no strong evidence of correlation)";
  return { stat: H, pval, conclusion };
}

// ============== STATE & SIMULATION ==============
var S = {
  seed: Math.floor(Math.random() * 100000),
  N: 500,
  T: 6,
  betaTrue: 0.5,
  sigmaB: 0.0,
  nGroups: 1,
  groupDeltaBeta: 0.0,
  sigmaU: 1.0,
  sigmaE: 1.0,
  sigmaME: 0.2,
  pBase: 0.5,
  rhoSel: 0.8,
  switchRate: 0.3,
  timeTrendOn: false,
  timeTrendSlope: 0.05,
  scatterWave: 1,
  spaghettiN: 20,
  colorByGroup: false,
  seType: "classical",
  specTrend: false,
  specQuad: false,
  specGroupFE: false,
  groupEstimator: "pooled",
  viz: "descr",
  modelDetail: "cs",
  diagModel: "pooled"
};

function simulate() {
  var rng = mulberry32(S.seed);
  var persons = new Array(S.N);
  for (var i = 0; i < S.N; i++) persons[i] = i;
  var u = persons.map(() => rnorm(rng, 0, S.sigmaU));
  var g = persons.map(() => (S.nGroups > 1 ? (1 + Math.floor(rng() * S.nGroups)) : 1));
  function groupBetaShift(gi) {
    if (S.nGroups <= 1) return 0;
    var center = (S.nGroups + 1) / 2;
    var denom = (S.nGroups - 1) / 2;
    return denom > 0 ? S.groupDeltaBeta * ((gi - center) / denom) : 0;
  }
  var betaI = persons.map((_, idx) => S.betaTrue + groupBetaShift(g[idx]) + rnorm(rng, 0, S.sigmaB));
  var alpha = Math.log(S.pBase / (1 - S.pBase));
  var rows = [];

  for (i = 0; i < S.N; i++) {
    var p_i = clamp(logistic(alpha + S.rhoSel * (u[i] / (S.sigmaU || 1))), 0.01, 0.99);
    var D = rng() < p_i ? 1 : 0;
    for (var t = 0; t < S.T; t++) {
      if (t > 0 && rng() < S.switchRate) D = rng() < p_i ? 1 : 0;
      var trend = S.timeTrendOn ? S.timeTrendSlope * (t + 1) : 0;
      var e_it = rnorm(rng, 0, S.sigmaE);
      var yStar = betaI[i] * D + u[i] + e_it + trend;
      var y = yStar + rnorm(rng, 0, S.sigmaME);
      rows.push({ pid: i, gid: g[i], beta_i: betaI[i], wave: t + 1, D: D, u: u[i], e: e_it, y: y, yStar: yStar, trend: trend });
    }
  }
  return rows;
}

// ============== LABELS & PRESETS ==============
function updateLabels() {
  document.getElementById("lblN").textContent = S.N;
  document.getElementById("lblT").textContent = S.T;
  document.getElementById("lblB").textContent = S.betaTrue.toFixed(2);
  document.getElementById("lblSb").textContent = S.sigmaB.toFixed(2);
  document.getElementById("lblG").textContent = String(S.nGroups);
  document.getElementById("lblGdb").textContent = S.groupDeltaBeta.toFixed(2);
  document.getElementById("lblSu").textContent = S.sigmaU.toFixed(2);
  document.getElementById("lblSe").textContent = S.sigmaE.toFixed(2);
  document.getElementById("lblSme").textContent = S.sigmaME.toFixed(2);
  document.getElementById("lblP").textContent = S.pBase.toFixed(2);
  document.getElementById("lblRho").textContent = S.rhoSel.toFixed(2);
  document.getElementById("lblSw").textContent = S.switchRate.toFixed(2);
  document.getElementById("lblTrend").textContent = S.timeTrendSlope.toFixed(2);

  var scatterEl = document.getElementById("scatterWave");
  scatterEl.max = String(S.T);
  if (S.scatterWave > S.T) S.scatterWave = S.T;
  scatterEl.value = String(S.scatterWave);
  document.getElementById("lblScatterWave").textContent = String(S.scatterWave);
  document.getElementById("lblSpaghetti").textContent = String(S.spaghettiN);
}

function applyPreset(name) {
  // Keep presets focused on selection/trends/noise; reset heterogeneity add-ons.
  S.sigmaB = 0;
  S.nGroups = 1;
  S.groupDeltaBeta = 0;

  if (name === "Baseline (No Bias)") {
    S.rhoSel = 0; S.pBase = 0.5; S.sigmaME = 0.2; S.timeTrendOn = false; S.switchRate = 0.4; S.sigmaU = 1.0;
  } else if (name === "Selection Bias") {
    S.rhoSel = 1.5; S.pBase = 0.5; S.sigmaME = 0.2; S.timeTrendOn = false; S.switchRate = 0.3; S.sigmaU = 1.2;
  } else if (name === "RCT-like") {
    S.rhoSel = 0.0; S.pBase = 0.5; S.switchRate = 0.5; S.timeTrendOn = false; S.sigmaME = 0.1;
  } else if (name === "Noisy Survey") {
    S.sigmaME = 1.0; S.sigmaE = 1.5; S.rhoSel = 0.3; S.timeTrendOn = false;
  } else if (name === "Time Trend") {
    S.timeTrendOn = true; S.timeTrendSlope = 0.15; S.rhoSel = 0.5;
  } else if (name === "Low Within Variation") {
    S.switchRate = 0.05; S.rhoSel = 0.8; S.timeTrendOn = false;
  }
  S.seed = Math.floor(Math.random() * 100000);
  syncInputs();
  computeAndRender();
}

function syncInputs() {
  document.getElementById("N").value = S.N;
  document.getElementById("T").value = S.T;
  document.getElementById("betaTrue").value = S.betaTrue;
  document.getElementById("sigmaB").value = S.sigmaB;
  document.getElementById("nGroups").value = S.nGroups;
  document.getElementById("groupDeltaBeta").value = S.groupDeltaBeta;
  document.getElementById("sigmaU").value = S.sigmaU;
  document.getElementById("sigmaE").value = S.sigmaE;
  document.getElementById("sigmaME").value = S.sigmaME;
  document.getElementById("pBase").value = S.pBase;
  document.getElementById("rhoSel").value = S.rhoSel;
  document.getElementById("switchRate").value = S.switchRate;
  document.getElementById("timeTrendSlope").value = S.timeTrendSlope;
  document.getElementById("timeTrendOn").checked = S.timeTrendOn;
  document.getElementById("scatterWave").value = S.scatterWave;
  document.getElementById("spaghettiN").value = S.spaghettiN;
  document.getElementById("colorByGroup").checked = S.colorByGroup;
  document.getElementById("seType").value = S.seType;
  document.getElementById("specTrend").checked = S.specTrend;
  document.getElementById("specQuad").checked = S.specQuad;
  document.getElementById("specGroupFE").checked = S.specGroupFE;
  document.getElementById("groupEstimator").value = S.groupEstimator;
  updateLabels();
}

// ============== SVG DRAWING ==============
function createSVG(w, h) {
  var svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
  svg.setAttribute("viewBox", "0 0 " + w + " " + h);
  svg.setAttribute("width", "100%");
  svg.setAttribute("height", "100%");
  return svg;
}

function line(x1, y1, x2, y2, stroke = "#6b7280", sw = 1) {
  var el = document.createElementNS("http://www.w3.org/2000/svg", "line");
  el.setAttribute("x1", x1); el.setAttribute("y1", y1);
  el.setAttribute("x2", x2); el.setAttribute("y2", y2);
  el.setAttribute("stroke", stroke); el.setAttribute("stroke-width", sw);
  return el;
}

function rect(x, y, w, h, fill = "#60a5fa") {
  var el = document.createElementNS("http://www.w3.org/2000/svg", "rect");
  el.setAttribute("x", x); el.setAttribute("y", y);
  el.setAttribute("width", w); el.setAttribute("height", h);
  el.setAttribute("fill", fill);
  return el;
}

function circle(x, y, r, fill = "#cbd5e1", op = 0.8) {
  var el = document.createElementNS("http://www.w3.org/2000/svg", "circle");
  el.setAttribute("cx", x); el.setAttribute("cy", y); el.setAttribute("r", r);
  el.setAttribute("fill", fill); el.setAttribute("fill-opacity", op);
  return el;
}

function polyline(points, stroke = "#93c5fd", sw = 1, op = 1) {
  var el = document.createElementNS("http://www.w3.org/2000/svg", "polyline");
  el.setAttribute("points", points.map(p => p[0] + "," + p[1]).join(" "));
  el.setAttribute("fill", "none");
  el.setAttribute("stroke", stroke);
  el.setAttribute("stroke-width", sw);
  el.setAttribute("stroke-opacity", op);
  return el;
}

function textSVG(x, y, txt, fill = "#9aa0a6", fs = 10, anchor = "middle") {
  var el = document.createElementNS("http://www.w3.org/2000/svg", "text");
  el.setAttribute("x", x); el.setAttribute("y", y);
  el.setAttribute("fill", fill); el.setAttribute("font-size", fs);
  el.setAttribute("text-anchor", anchor);
  el.textContent = txt;
  return el;
}

function axis(svg, x, y, len, horizontal, ticks) {
  if (horizontal) {
    svg.appendChild(line(x, y, x + len, y, "#374151", 1));
    for (var i = 0; i < ticks.length; i++) {
      var t = ticks[i], tx = x + t.pos * len;
      svg.appendChild(line(tx, y, tx, y + 4, "#4b5563", 1));
      svg.appendChild(textSVG(tx, y + 14, t.label, "#9aa0a6", 10, "middle"));
    }
  } else {
    svg.appendChild(line(x, y - len, x, y, "#374151", 1));
    for (var j = 0; j < ticks.length; j++) {
      var tt = ticks[j], ty = y - tt.pos * len;
      svg.appendChild(line(x - 4, ty, x, ty, "#4b5563", 1));
      svg.appendChild(textSVG(x - 6, ty + 3, tt.label, "#9aa0a6", 10, "end"));
    }
  }
}

function scaleLinear(d0, d1, r0, r1) {
  var d = d1 - d0 || 1e-9;
  var m = (r1 - r0) / d;
  return v => r0 + (v - d0) * m;
}

// ============== CHART DRAWING FUNCTIONS ==============
function groupColor(gid) {
  var palette = ["#60a5fa", "#34d399", "#fbbf24", "#a78bfa", "#f472b6"];
  var idx = Math.max(1, gid || 1) - 1;
  return palette[idx % palette.length];
}

function drawScatter(container, data, model, opts = {}) {
  container.innerHTML = "";
  var w = 480, h = 300, m = { l: 45, r: 15, t: 15, b: 35 };
  var yMin = Math.min.apply(null, data.map(d => d.y));
  var yMax = Math.max.apply(null, data.map(d => d.y));
  var pad = (yMax - yMin) * 0.1;
  yMin -= pad; yMax += pad;
  var xs = scaleLinear(-0.1, 1.1, m.l, w - m.r);
  var ys = scaleLinear(yMin, yMax, h - m.b, m.t);
  var svg = createSVG(w, h);

  // Grid
  for (var i = 0; i <= 10; i++) {
    var gx = m.l + i * (w - m.l - m.r) / 10;
    svg.appendChild(line(gx, m.t, gx, h - m.b, "#1a1a1a", 1));
  }
  for (i = 0; i <= 8; i++) {
    var gy = m.t + i * (h - m.t - m.b) / 8;
    svg.appendChild(line(m.l, gy, w - m.r, gy, "#1a1a1a", 1));
  }

  // Axes
  axis(svg, m.l, h - m.b, w - m.l - m.r, true, [
    { pos: (0 - (-0.1)) / 1.2, label: "D=0" },
    { pos: (1 - (-0.1)) / 1.2, label: "D=1" }
  ]);
  var yTicks = [0, 0.25, 0.5, 0.75, 1].map(p => ({ pos: p, label: fmt(yMin + p * (yMax - yMin), 1) }));
  axis(svg, m.l, h - m.b, h - m.t - m.b, false, yTicks);

  // Axis labels
  svg.appendChild(textSVG(w / 2, h - 5, "Treatment (D)", "#6b7280", 11, "middle"));

  // Points with jitter
  var rng = mulberry32(42);
  for (i = 0; i < data.length; i++) {
    var jx = (rng() - 0.5) * 0.15;
    var col = opts.colorByGroup ? groupColor(data[i].gid) : "#64748b";
    svg.appendChild(circle(xs(data[i].x + jx), ys(data[i].y), 2.5, col, opts.colorByGroup ? 0.75 : 0.6));
  }

  // Regression line
  if (model) {
    var b0 = coefAt(model, "Intercept").b;
    if (!isFinite(b0)) b0 = 0;
    var bD = coefAt(model, "D").b;
    if (isFinite(bD)) {
      var p1 = [xs(0), ys(b0)];
      var p2 = [xs(1), ys(b0 + bD)];
      svg.appendChild(polyline([p1, p2], "#93c5fd", 2.5, 1));
    }
  }

  container.appendChild(svg);
}

function drawLines(container, spaghetti, summary, T, opts = {}) {
  container.innerHTML = "";
  var w = 480, h = 300, m = { l: 45, r: 15, t: 15, b: 35 };
  var allY = [];
  for (var s = 0; s < spaghetti.length; s++) {
    for (var q = 0; q < spaghetti[s].series.length; q++) {
      allY.push(spaghetti[s].series[q].y);
    }
  }
  for (var r = 0; r < summary.length; r++) {
    allY.push(summary[r].mean, summary[r].median, summary[r].p25, summary[r].p75);
  }
  var yMin = Math.min.apply(null, allY), yMax = Math.max.apply(null, allY);
  var pad = (yMax - yMin) * 0.1;
  yMin -= pad; yMax += pad;
  var xs = scaleLinear(0.5, T + 0.5, m.l, w - m.r);
  var ys = scaleLinear(yMin, yMax, h - m.b, m.t);
  var svg = createSVG(w, h);

  // Grid
  for (var i = 1; i <= T; i++) {
    var gx = xs(i);
    svg.appendChild(line(gx, m.t, gx, h - m.b, "#1a1a1a", 1));
  }
  for (i = 0; i <= 8; i++) {
    var gy = m.t + i * (h - m.t - m.b) / 8;
    svg.appendChild(line(m.l, gy, w - m.r, gy, "#1a1a1a", 1));
  }

  // Axes
  var xTicks = [];
  for (i = 1; i <= T; i++) xTicks.push({ pos: (i - 0.5) / T, label: String(i) });
  axis(svg, m.l, h - m.b, w - m.l - m.r, true, xTicks);
  var yTicks = [0, 0.25, 0.5, 0.75, 1].map(p => ({ pos: p, label: fmt(yMin + p * (yMax - yMin), 1) }));
  axis(svg, m.l, h - m.b, h - m.t - m.b, false, yTicks);

  // Axis labels
  svg.appendChild(textSVG(w / 2, h - 5, "Wave (t)", "#6b7280", 11, "middle"));

  // Individual trajectories
  for (s = 0; s < spaghetti.length; s++) {
    var pts = [];
    for (q = 0; q < spaghetti[s].series.length; q++) {
      pts.push([xs(spaghetti[s].series[q].t), ys(spaghetti[s].series[q].y)]);
    }
    var col = opts.colorByGroup ? groupColor(spaghetti[s].gid) : "#4b5563";
    svg.appendChild(polyline(pts, col, 1, 0.35));
  }

  // Summary lines
  function pathFor(key, color, width) {
    var pts = [];
    for (i = 0; i < summary.length; i++) {
      pts.push([xs(summary[i].t), ys(summary[i][key])]);
    }
    svg.appendChild(polyline(pts, color, width, 1));
  }
  pathFor("p25", "#a78bfa", 1.5);
  pathFor("p75", "#a78bfa", 1.5);
  pathFor("median", "#e5e7eb", 2);
  pathFor("mean", "#93c5fd", 2.5);

  container.appendChild(svg);
}

function drawWaveMeans(container, rows, T) {
  container.innerHTML = "";
  var w = 480, h = 300, m = { l: 45, r: 15, t: 15, b: 35 };

  var series0 = [], series1 = [];
  for (var t = 1; t <= T; t++) {
    var y0 = [], y1 = [];
    for (var i = 0; i < rows.length; i++) {
      if (rows[i].wave !== t) continue;
      if (rows[i].D === 1) y1.push(rows[i].y);
      else y0.push(rows[i].y);
    }
    series0.push({ t, y: y0.length ? mean(y0) : NaN });
    series1.push({ t, y: y1.length ? mean(y1) : NaN });
  }

  var allY = series0.concat(series1).map(p => p.y).filter(v => isFinite(v));
  if (!allY.length) {
    container.innerHTML = "<div style='padding:20px;color:#6b7280'>No data</div>";
    return;
  }

  var yMin = Math.min.apply(null, allY), yMax = Math.max.apply(null, allY);
  var pad = (yMax - yMin) * 0.1;
  yMin -= pad; yMax += pad;
  var xs = scaleLinear(0.5, T + 0.5, m.l, w - m.r);
  var ys = scaleLinear(yMin, yMax, h - m.b, m.t);
  var svg = createSVG(w, h);

  for (var gx = 1; gx <= T; gx++) svg.appendChild(line(xs(gx), m.t, xs(gx), h - m.b, "#1a1a1a", 1));
  for (var gy = 0; gy <= 8; gy++) {
    var yy = m.t + gy * (h - m.t - m.b) / 8;
    svg.appendChild(line(m.l, yy, w - m.r, yy, "#1a1a1a", 1));
  }

  var xTicks = [];
  for (var tt = 1; tt <= T; tt++) xTicks.push({ pos: (tt - 0.5) / T, label: String(tt) });
  axis(svg, m.l, h - m.b, w - m.l - m.r, true, xTicks);
  var yTicks = [0, 0.25, 0.5, 0.75, 1].map(p => ({ pos: p, label: fmt(yMin + p * (yMax - yMin), 1) }));
  axis(svg, m.l, h - m.b, h - m.t - m.b, false, yTicks);

  svg.appendChild(textSVG(w / 2, h - 5, "Wave (t)", "#6b7280", 11, "middle"));

  function drawSeries(ser, stroke) {
    var pts = [];
    for (var i = 0; i < ser.length; i++) if (isFinite(ser[i].y)) pts.push([xs(ser[i].t), ys(ser[i].y)]);
    if (pts.length >= 2) svg.appendChild(polyline(pts, stroke, 2.5, 0.95));
    for (i = 0; i < ser.length; i++) if (isFinite(ser[i].y)) svg.appendChild(circle(xs(ser[i].t), ys(ser[i].y), 2.8, stroke, 0.9));
  }
  drawSeries(series0, "#e5e7eb");
  drawSeries(series1, "#60a5fa");

  container.appendChild(svg);
}

function drawBetaDist(container, rows, opts = {}) {
  container.innerHTML = "";
  var w = 480, h = 300, m = { l: 45, r: 15, t: 15, b: 35 };

  var byPid = {};
  for (var i = 0; i < rows.length; i++) {
    var pid = rows[i].pid;
    if (byPid[pid] == null) byPid[pid] = rows[i].beta_i;
  }
  var vals = Object.values(byPid).filter(v => isFinite(v));
  if (!vals.length) {
    container.innerHTML = "<div style='padding:20px;color:#6b7280'>No data</div>";
    return;
  }
  var mn = Math.min.apply(null, vals), mx = Math.max.apply(null, vals);
  if (mn === mx) { mn -= 0.5; mx += 0.5; }
  var bins = 25;
  var width = (mx - mn) / bins || 1;
  var cnt = new Array(bins).fill(0);
  for (i = 0; i < vals.length; i++) {
    var v = vals[i];
    var idx = Math.max(0, Math.min(bins - 1, Math.floor((v - mn) / width)));
    cnt[idx]++;
  }
  var maxCnt = Math.max.apply(null, cnt);
  var xs = scaleLinear(mn, mx, m.l, w - m.r);
  var ys = scaleLinear(0, maxCnt, h - m.b, m.t);
  var svg = createSVG(w, h);

  for (var gy = 0; gy <= 6; gy++) {
    var yy = m.t + gy * (h - m.t - m.b) / 6;
    svg.appendChild(line(m.l, yy, w - m.r, yy, "#1a1a1a", 1));
  }

  var bw = (w - m.l - m.r) / bins - 1;
  for (i = 0; i < bins; i++) {
    var x0 = m.l + i * ((w - m.l - m.r) / bins) + 0.5;
    var barH = ys(0) - ys(cnt[i]);
    svg.appendChild(rect(x0, ys(cnt[i]), bw, barH, "#a3a3a3"));
  }

  var meanLine = opts.betaMean;
  if (isFinite(meanLine)) {
    var xm = xs(meanLine);
    svg.appendChild(line(xm, m.t, xm, h - m.b, "#f87171", 2));
  }

  axis(svg, m.l, h - m.b, w - m.l - m.r, true, [
    { pos: 0, label: fmt(mn, 2) },
    { pos: 0.5, label: fmt((mn + mx) / 2, 2) },
    { pos: 1, label: fmt(mx, 2) }
  ]);
  axis(svg, m.l, h - m.b, h - m.t - m.b, false, [
    { pos: 0, label: "0" },
    { pos: 1, label: String(maxCnt) }
  ]);

  svg.appendChild(textSVG(w / 2, h - 5, "True individual effect (βᵢ)", "#6b7280", 11, "middle"));
  container.appendChild(svg);
}

function drawBarsWithCI(container, rows, betaTrue) {
  container.innerHTML = "";
  var w = 480, h = 300, m = { l: 55, r: 15, t: 20, b: 70 };
  var allY = [betaTrue, 0];
  for (var i = 0; i < rows.length; i++) {
    if (isFinite(rows[i].b)) allY.push(rows[i].b);
    if (rows[i].se && isFinite(rows[i].se)) {
      var c = ci(rows[i].b, rows[i].se);
      allY.push(c[0], c[1]);
    }
  }
  var yMin = Math.min.apply(null, allY), yMax = Math.max.apply(null, allY);
  var pad = (yMax - yMin) * 0.15;
  yMin -= pad; yMax += pad;

  var barWidth = (w - m.l - m.r) / rows.length;
  var xs = i => m.l + i * barWidth + barWidth / 2;
  var ys = scaleLinear(yMin, yMax, h - m.b, m.t);
  var svg = createSVG(w, h);

  // Grid
  for (i = 0; i <= 8; i++) {
    var gy = m.t + i * (h - m.t - m.b) / 8;
    svg.appendChild(line(m.l, gy, w - m.r, gy, "#1a1a1a", 1));
  }

  // Zero line
  var y0 = ys(0);
  svg.appendChild(line(m.l, y0, w - m.r, y0, "#374151", 1));

  // True beta line
  var ytrue = ys(betaTrue);
  svg.appendChild(line(m.l, ytrue, w - m.r, ytrue, "#f87171", 2));
  svg.appendChild(textSVG(w - m.r + 5, ytrue + 4, "True β", "#f87171", 10, "start"));

  // Y axis
  var yTicks = [0, 0.25, 0.5, 0.75, 1].map(p => ({ pos: p, label: fmt(yMin + p * (yMax - yMin), 2) }));
  axis(svg, m.l, h - m.b, h - m.t - m.b, false, yTicks);

  // Bars and CIs
  var bw = barWidth * 0.6;
  for (i = 0; i < rows.length; i++) {
    var cx = xs(i);
    var barY = ys(rows[i].b);
    var barH = Math.abs(barY - y0);
    var barTop = Math.min(barY, y0);

    // Bar
    svg.appendChild(rect(cx - bw / 2, barTop, bw, barH, rows[i].color || "#60a5fa"));

    // CI whiskers
    if (rows[i].se && isFinite(rows[i].se)) {
      var C = ci(rows[i].b, rows[i].se);
      var ylo = ys(C[0]), yhi = ys(C[1]);
      svg.appendChild(line(cx, ylo, cx, yhi, "#e5e7eb", 2));
      svg.appendChild(line(cx - 8, ylo, cx + 8, ylo, "#e5e7eb", 2));
      svg.appendChild(line(cx - 8, yhi, cx + 8, yhi, "#e5e7eb", 2));
    }

    // Label
    svg.appendChild(textSVG(cx, h - m.b + 15, rows[i].name, "#9aa0a6", 10, "middle"));
    svg.appendChild(textSVG(cx, h - m.b + 28, "β̂=" + fmt(rows[i].b, 2), "#6b7280", 9, "middle"));
  }

  container.appendChild(svg);
}

function drawResiduals(container, points) {
  container.innerHTML = "";
  var w = 480, h = 280, m = { l: 45, r: 15, t: 15, b: 35 };
  if (!points || points.length === 0) {
    container.innerHTML = "<div style='padding:20px;color:#6b7280'>No data available</div>";
    return;
  }

  var fitMin = points[0].fit, fitMax = points[0].fit;
  var rMin = points[0].resid, rMax = points[0].resid;
  for (var i = 1; i < points.length; i++) {
    if (points[i].fit < fitMin) fitMin = points[i].fit;
    if (points[i].fit > fitMax) fitMax = points[i].fit;
    if (points[i].resid < rMin) rMin = points[i].resid;
    if (points[i].resid > rMax) rMax = points[i].resid;
  }
  var padX = (fitMax - fitMin) * 0.1;
  var padY = (rMax - rMin) * 0.1;
  fitMin -= padX; fitMax += padX;
  rMin -= padY; rMax += padY;

  var xs = scaleLinear(fitMin, fitMax, m.l, w - m.r);
  var ys = scaleLinear(rMin, rMax, h - m.b, m.t);
  var svg = createSVG(w, h);

  // Grid
  for (i = 0; i <= 8; i++) {
    var gx = m.l + i * (w - m.l - m.r) / 8;
    var gy = m.t + i * (h - m.t - m.b) / 8;
    svg.appendChild(line(gx, m.t, gx, h - m.b, "#1a1a1a", 1));
    svg.appendChild(line(m.l, gy, w - m.r, gy, "#1a1a1a", 1));
  }

  // Zero line
  var y0 = ys(0);
  svg.appendChild(line(m.l, y0, w - m.r, y0, "#94a3b8", 1.5));

  // Points
  for (i = 0; i < points.length; i++) {
    svg.appendChild(circle(xs(points[i].fit), ys(points[i].resid), 2, "#64748b", 0.5));
  }

  // Axes
  var xTicks = [0, 0.5, 1].map(p => ({ pos: p, label: fmt(fitMin + p * (fitMax - fitMin), 1) }));
  var yTicks = [0, 0.5, 1].map(p => ({ pos: p, label: fmt(rMin + p * (rMax - rMin), 1) }));
  axis(svg, m.l, h - m.b, w - m.l - m.r, true, xTicks);
  axis(svg, m.l, h - m.b, h - m.t - m.b, false, yTicks);

  svg.appendChild(textSVG(w / 2, h - 5, "Fitted values", "#6b7280", 10, "middle"));

  container.appendChild(svg);
}

function drawHist(container, resid) {
  container.innerHTML = "";
  var w = 480, h = 280, m = { l: 45, r: 15, t: 15, b: 35 };
  if (resid.length === 0) {
    container.innerHTML = "<div style='padding:20px;color:#6b7280'>No data</div>";
    return;
  }

  var mn = Math.min.apply(null, resid), mx = Math.max.apply(null, resid);
  var bins = 25;
  var width = (mx - mn) / bins || 1;
  var cnt = new Array(bins).fill(0);
  for (var i = 0; i < resid.length; i++) {
    var v = resid[i];
    var idx = Math.max(0, Math.min(bins - 1, Math.floor((v - mn) / width)));
    cnt[idx]++;
  }

  var xs = scaleLinear(0, bins, m.l, w - m.r);
  var maxCnt = Math.max.apply(null, cnt);
  var ys = scaleLinear(0, maxCnt, h - m.b, m.t);
  var svg = createSVG(w, h);

  // Grid
  for (i = 0; i <= 5; i++) {
    var gy = m.t + i * (h - m.t - m.b) / 5;
    svg.appendChild(line(m.l, gy, w - m.r, gy, "#1a1a1a", 1));
  }

  // Bars
  var bw = (w - m.l - m.r) / bins - 1;
  for (i = 0; i < bins; i++) {
    var x = xs(i) + 0.5;
    var barH = ys(0) - ys(cnt[i]);
    svg.appendChild(rect(x, ys(cnt[i]), bw, barH, "#60a5fa"));
  }

  // Axes
  axis(svg, m.l, h - m.b, w - m.l - m.r, true, [
    { pos: 0, label: fmt(mn, 1) },
    { pos: 0.5, label: fmt((mn + mx) / 2, 1) },
    { pos: 1, label: fmt(mx, 1) }
  ]);
  axis(svg, m.l, h - m.b, h - m.t - m.b, false, [
    { pos: 0, label: "0" },
    { pos: 1, label: String(maxCnt) }
  ]);

  svg.appendChild(textSVG(w / 2, h - 5, "Residual value", "#6b7280", 10, "middle"));

  container.appendChild(svg);
}

function drawQQ(container, resid) {
  container.innerHTML = "";
  var w = 480, h = 280, m = { l: 45, r: 15, t: 15, b: 35 };
  var r = resid.slice().sort((a, b) => a - b);
  if (r.length < 5) {
    container.innerHTML = "<div style='padding:20px;color:#6b7280'>Not enough data</div>";
    return;
  }

  function erfinv(x) {
    var a = 0.147;
    var ln = Math.log(1 - x * x);
    var s = 2 / (Math.PI * a) + ln / 2;
    return Math.sign(x) * Math.sqrt(Math.sqrt(s * s - ln / a) - s);
  }

  var n = r.length, pts = [];
  for (var i = 0; i < n; i++) {
    var p = (i + 0.5) / n;
    var z = Math.SQRT2 * erfinv(2 * p - 1);
    pts.push({ theor: z, sample: r[i] });
  }

  var xs = scaleLinear(pts[0].theor, pts[pts.length - 1].theor, m.l, w - m.r);
  var ys = scaleLinear(r[0], r[r.length - 1], h - m.b, m.t);
  var svg = createSVG(w, h);

  // Grid
  for (i = 0; i <= 8; i++) {
    var gx = m.l + i * (w - m.l - m.r) / 8;
    var gy = m.t + i * (h - m.t - m.b) / 8;
    svg.appendChild(line(gx, m.t, gx, h - m.b, "#1a1a1a", 1));
    svg.appendChild(line(m.l, gy, w - m.r, gy, "#1a1a1a", 1));
  }

  // Reference line
  svg.appendChild(polyline([
    [xs(pts[0].theor), ys(r[0])],
    [xs(pts[pts.length - 1].theor), ys(r[r.length - 1])]
  ], "#94a3b8", 1.5));

  // Points
  for (i = 0; i < pts.length; i++) {
    svg.appendChild(circle(xs(pts[i].theor), ys(pts[i].sample), 2, "#64748b", 0.5));
  }

  // Axes
  axis(svg, m.l, h - m.b, w - m.l - m.r, true, [
    { pos: 0, label: fmt(pts[0].theor, 1) },
    { pos: 1, label: fmt(pts[pts.length - 1].theor, 1) }
  ]);
  axis(svg, m.l, h - m.b, h - m.t - m.b, false, [
    { pos: 0, label: fmt(r[0], 1) },
    { pos: 1, label: fmt(r[r.length - 1], 1) }
  ]);

  svg.appendChild(textSVG(w / 2, h - 5, "Theoretical quantiles", "#6b7280", 10, "middle"));

  container.appendChild(svg);
}

function drawHorizontalCIPlot(container, items, opts = {}) {
  container.innerHTML = "";
  if (!items || !items.length) {
    container.innerHTML = "<div style='padding:20px;color:#6b7280'>No data</div>";
    return;
  }

  var w = 480, h = 300, m = { l: 110, r: 15, t: 15, b: 25 };
  var xsMin = Infinity, xsMax = -Infinity;
  for (var i = 0; i < items.length; i++) {
    var it = items[i];
    if (!isFinite(it.b) || !isFinite(it.se)) continue;
    var C = ci(it.b, it.se);
    xsMin = Math.min(xsMin, C[0]);
    xsMax = Math.max(xsMax, C[1]);
  }
  if (!isFinite(xsMin) || !isFinite(xsMax)) {
    container.innerHTML = "<div style='padding:20px;color:#6b7280'>Not enough data</div>";
    return;
  }
  if (xsMin === xsMax) { xsMin -= 0.5; xsMax += 0.5; }
  var pad = (xsMax - xsMin) * 0.1;
  xsMin -= pad; xsMax += pad;

  var xs = scaleLinear(xsMin, xsMax, m.l, w - m.r);
  var yPos = idx => m.t + idx * ((h - m.t - m.b) / Math.max(1, items.length - 1));
  var svg = createSVG(w, h);

  // Grid
  for (var gx = 0; gx <= 8; gx++) {
    var xx = m.l + gx * (w - m.l - m.r) / 8;
    svg.appendChild(line(xx, m.t, xx, h - m.b, "#1a1a1a", 1));
  }

  // Reference lines
  if (opts.vline != null && isFinite(opts.vline)) {
    var xv = xs(opts.vline);
    svg.appendChild(line(xv, m.t, xv, h - m.b, "#f87171", 2));
  }
  if (opts.zeroLine) {
    var x0 = xs(0);
    svg.appendChild(line(x0, m.t, x0, h - m.b, "#374151", 1.5));
  }

  // X axis
  axis(svg, m.l, h - m.b, w - m.l - m.r, true, [
    { pos: 0, label: fmt(xsMin, 2) },
    { pos: 0.5, label: fmt((xsMin + xsMax) / 2, 2) },
    { pos: 1, label: fmt(xsMax, 2) }
  ]);

  for (i = 0; i < items.length; i++) {
    var it2 = items[i];
    var y = yPos(i);
    svg.appendChild(textSVG(m.l - 8, y + 4, it2.label, "#9aa0a6", 11, "end"));
    if (!isFinite(it2.b) || !isFinite(it2.se)) continue;
    var C2 = ci(it2.b, it2.se);
    svg.appendChild(line(xs(C2[0]), y, xs(C2[1]), y, "#e5e7eb", 2));
    svg.appendChild(circle(xs(it2.b), y, 4, it2.color || "#60a5fa", 0.95));
    if (it2.trueVal != null && isFinite(it2.trueVal)) {
      var xt = xs(it2.trueVal);
      svg.appendChild(line(xt, y - 6, xt, y + 6, "#f87171", 2));
    }
  }

  if (opts.xLabel) svg.appendChild(textSVG((m.l + (w - m.r)) / 2, h - 4, opts.xLabel, "#6b7280", 11, "middle"));
  container.appendChild(svg);
}

function drawSeSensitivity(container, which, cross, pooled, models) {
  var labels = {
    classical: { label: "Classical", color: "#60a5fa" },
    robust: { label: "Robust (HC1)", color: "#a78bfa" },
    cluster: { label: "Cluster (CR1)", color: "#34d399" }
  };

  function estimateModel(kind, seType) {
    if (kind === "cs") {
      var dsgn = buildDesign(cross, { includeIntercept: true, includeTrend: false, includeQuad: false, includeGroupFE: S.specGroupFE });
      return ols(dsgn.X, dsgn.y, { seType, clusterIds: dsgn.ids, coefNames: dsgn.coefNames });
    }
    if (kind === "pooled") {
      var dsgn2 = buildDesign(pooled, { includeIntercept: true, includeTrend: S.specTrend, includeQuad: S.specQuad, includeGroupFE: S.specGroupFE });
      return ols(dsgn2.X, dsgn2.y, { seType, clusterIds: dsgn2.ids, coefNames: dsgn2.coefNames });
    }
    if (kind === "fe") {
      var dsgn3 = buildDesign(pooled, { includeIntercept: false, includeTrend: S.specTrend, includeQuad: S.specQuad, includeGroupFE: false });
      return feWithin(dsgn3.ids, dsgn3.X, dsgn3.y, { seType, coefNames: dsgn3.coefNames });
    }
    if (kind === "re") {
      var dsgn4 = buildDesign(pooled, { includeIntercept: true, includeTrend: S.specTrend, includeQuad: S.specQuad, includeGroupFE: S.specGroupFE });
      return randomEffects(dsgn4.ids, dsgn4.X, dsgn4.y, { seType, coefNames: dsgn4.coefNames });
    }
    return null;
  }

  var base = models && models[which === "pooled" ? "pooledOLS" : which];
  var b0 = coefAt(base, "D").b;
  var items = [];
  ["classical", "robust", "cluster"].forEach(seType => {
    var fit = estimateModel(which, seType);
    var d = coefAt(fit, "D");
    items.push({ label: labels[seType].label, b: d.b, se: d.se, color: labels[seType].color, trueVal: b0 });
  });

  drawHorizontalCIPlot(container, items, { xLabel: "β̂ (for D)", zeroLine: true, vline: S.betaTrue });
}

function drawGroupSplit(container, pooled, opts = {}) {
  if (S.nGroups <= 1) {
    container.innerHTML = "<div style='padding:20px;color:#6b7280'>Set Groups (G) &gt; 1 to split regressions.</div>";
    return;
  }

  function groupBetaShift(gi) {
    if (S.nGroups <= 1) return 0;
    var center = (S.nGroups + 1) / 2;
    var denom = (S.nGroups - 1) / 2;
    return denom > 0 ? S.groupDeltaBeta * ((gi - center) / denom) : 0;
  }

  var items = [];
  for (var g = 1; g <= S.nGroups; g++) {
    var subset = pooled.filter(r => r.gid === g);
    var seType = opts.seType || "classical";
    var estimator = opts.estimator || "pooled";
    var fit;
    if (estimator === "fe") {
      var d1 = buildDesign(subset, { includeIntercept: false, includeTrend: !!opts.includeTrend, includeQuad: !!opts.includeQuad, includeGroupFE: false });
      fit = feWithin(d1.ids, d1.X, d1.y, { seType, coefNames: d1.coefNames });
    } else if (estimator === "re") {
      var d2 = buildDesign(subset, { includeIntercept: true, includeTrend: !!opts.includeTrend, includeQuad: !!opts.includeQuad, includeGroupFE: false });
      fit = randomEffects(d2.ids, d2.X, d2.y, { seType, coefNames: d2.coefNames });
    } else {
      var d3 = buildDesign(subset, { includeIntercept: true, includeTrend: !!opts.includeTrend, includeQuad: !!opts.includeQuad, includeGroupFE: false });
      fit = ols(d3.X, d3.y, { seType, clusterIds: d3.ids, coefNames: d3.coefNames });
    }
    var d = coefAt(fit, "D");
    items.push({ label: "G" + g, b: d.b, se: d.se, color: groupColor(g), trueVal: S.betaTrue + groupBetaShift(g) });
  }

  drawHorizontalCIPlot(container, items, { xLabel: "β̂ (for D)", zeroLine: true, vline: S.betaTrue });
}

// ============== MODEL DETAILS ==============
function fillModelDetails(which, cs, pooled, fe, re) {
  var content = "";
  var model, name, description, assumptions;

  if (which === "cs") {
    model = cs; name = "Cross-Section OLS (Selected Wave)";
    description = "Uses one wave of data, treating it as a simple cross-section.";
    assumptions = "Assumes no omitted variable bias: E[uᵢ|Dᵢ] = 0";
  } else if (which === "pooled") {
    model = pooled; name = "Pooled OLS (All Waves Stacked)";
    description = "Stacks all waves and ignores the panel structure. More data, but ignores within-person correlation.";
    assumptions = "Assumes strict exogeneity: E[uᵢ + εᵢₜ|Dᵢₜ] = 0 for all t";
  } else if (which === "fe") {
    model = fe; name = "Fixed Effects (Within Estimator)";
    description = "Demeans each variable by individual, eliminating all time-invariant factors including uᵢ.";
    assumptions = "Assumes strict exogeneity of time-varying shocks: E[εᵢₜ|Dᵢₛ, uᵢ] = 0 for all s,t";
  } else if (which === "re") {
    model = re; name = "Random Effects (GLS)";
    description = "Quasi-demeans with θ to account for serial correlation while retaining between-variation.";
    assumptions = "Assumes RE exogeneity: E[uᵢ|Dᵢₜ] = 0 (stronger than FE)";
  }

  var d = coefAt(model, "D");
  var sig = d.idx >= 0 ? (model.pval[d.idx] < 0.001 ? "***" : model.pval[d.idx] < 0.01 ? "**" : model.pval[d.idx] < 0.05 ? "*" : "") : "";
  var sigClass = d.idx >= 0 && model.pval[d.idx] < 0.05 ? "sig" : "not-sig";

  content += `<h3 style="margin-bottom:12px">${name}</h3>`;
  content += `<p class="help-text muted" style="margin-bottom:12px">${description}</p>`;
  content += `<div class="help-text muted" style="margin:-6px 0 10px">SE type: <strong>${model.seType || "classical"}</strong></div>`;

  content += `<table class="stat-table">`;
  content += `<thead><tr><th>Term</th><th class="num">Coef</th><th class="num">SE</th><th class="num">t</th><th class="num">p</th></tr></thead><tbody>`;
  for (var i = 0; i < model.coefNames.length; i++) {
    var nm = model.coefNames[i] || ("x" + i);
    var cls = nm === "D" ? sigClass : "";
    var stars = nm === "D" ? sig : "";
    content += `<tr><td>${nm === "t2" ? "t²" : nm}</td><td class="num"><span class="${cls}">${fmt(model.beta[i])} ${stars}</span></td><td class="num">${fmt(model.se[i])}</td><td class="num">${fmt(model.tstat[i], 2)}</td><td class="num">${fmt(model.pval[i], 4)}</td></tr>`;
  }
  content += `</tbody></table>`;

  content += `<div class="mt-3"><table class="stat-table"><tbody>`;
  content += `<tr><td>R²</td><td class="num">${fmt(model.r2)}</td></tr>`;
  content += `<tr><td>N observations</td><td class="num">${model.n}</td></tr>`;
  content += `<tr><td>Degrees of freedom</td><td class="num">${model.df}</td></tr>`;
  if (which === "fe" || which === "re") content += `<tr><td>N groups</td><td class="num">${model.nGroups}</td></tr>`;
  if (which === "re" && isFinite(model.theta)) content += `<tr><td>θ (quasi-demean factor)</td><td class="num">${fmt(model.theta)}</td></tr>`;
  content += `</tbody></table></div>`;

  // 95% CI
  if (d.idx >= 0) {
    var C = ci(d.b, d.se);
    content += `<div class="mt-3"><strong>95% CI for D:</strong> [${fmt(C[0])}, ${fmt(C[1])}]</div>`;
  }

  // Interpretation
  var bias = d.idx >= 0 ? (d.b - S.betaTrue) : NaN;
  var biasPercent = S.betaTrue !== 0 ? (bias / S.betaTrue * 100) : 0;
  content += `<div class="interpretation mt-3">`;
  content += `<h4>Interpretation</h4>`;
  content += `<p><strong>Bias (for D):</strong> β̂ - β = ${fmt(bias)} (${fmt(biasPercent, 1)}% of true mean effect)</p>`;
  if (isFinite(bias) && isFinite(d.se) && Math.abs(bias) < d.se) {
    content += `<p class="mt-2" style="color:#86efac">✓ Estimate is within 1 SE of true value</p>`;
  } else if (isFinite(bias) && isFinite(d.se) && Math.abs(bias) < 2 * d.se) {
    content += `<p class="mt-2" style="color:#fde047">⚠ Estimate is within 2 SE of true value (acceptable)</p>`;
  } else {
    content += `<p class="mt-2" style="color:#fca5a5">✗ Estimate is more than 2 SE from true value (biased)</p>`;
  }
  content += `</div>`;

  content += `<div class="alert alert-info mt-3"><strong>Key assumption:</strong> ${assumptions}</div>`;

  document.getElementById("modelDetailContent").innerHTML = content;
}

// ============== SUMMARY CONTENT ==============
function buildSummary(sim, cs, pooled, fe, re) {
  var n = sim.length;
  var yVals = sim.map(r => r.y);
  var yMean = mean(yVals), ySd = Math.sqrt(variance(yVals));
  var Dmean = mean(sim.map(r => r.D));

  // Variance decomposition
  var by = {}, means = [], wVars = [];
  for (var i = 0; i < sim.length; i++) {
    var r = sim[i];
    if (!by[r.pid]) by[r.pid] = [];
    by[r.pid].push(r.y);
  }
  for (var k in by) {
    var arr = by[k];
    var m = mean(arr);
    means.push(m);
    var vv = 0;
    for (var j = 0; j < arr.length; j++) { vv += (arr[j] - m) * (arr[j] - m); }
    vv /= Math.max(1, arr.length - 1);
    wVars.push(vv);
  }
  var vb = variance(means);
  var vw = wVars.length ? wVars.reduce((a, b) => a + b, 0) / wVars.length : 0;
  var tot = vb + vw;
  var iccEmp = tot === 0 ? NaN : vb / tot;

  // Within-person variation in D
  var dBy = {};
  for (i = 0; i < sim.length; i++) {
    var r2 = sim[i];
    if (!dBy[r2.pid]) dBy[r2.pid] = [];
    dBy[r2.pid].push(r2.D);
  }
  var switchCount = 0, totalPairs = 0;
  for (k in dBy) {
    var dArr = dBy[k];
    for (j = 1; j < dArr.length; j++) {
      if (dArr[j] !== dArr[j - 1]) switchCount++;
      totalPairs++;
    }
  }
  var actualSwitchRate = totalPairs > 0 ? switchCount / totalPairs : 0;

  var html = `
    <div class="grid-2">
      <div>
        <h3>Sample Characteristics</h3>
        <table class="stat-table">
          <tr><td>Total observations</td><td class="num">${S.N} × ${S.T} = ${n}</td></tr>
          <tr><td>Outcome mean (SD)</td><td class="num">${fmt(yMean, 2)} (${fmt(ySd, 2)})</td></tr>
          <tr><td>Treatment prevalence</td><td class="num">${fmt(Dmean * 100, 1)}%</td></tr>
          <tr><td>Actual switching rate</td><td class="num">${fmt(actualSwitchRate * 100, 1)}%</td></tr>
          <tr><td>Empirical ICC</td><td class="num">${fmt(iccEmp)}</td></tr>
        </table>
      </div>
      <div>
        <h3>Current Parameters</h3>
        <table class="stat-table">
          <tr><td>True β</td><td class="num">${fmt(S.betaTrue)}</td></tr>
          <tr><td>Random slope SD (σβ)</td><td class="num">${fmt(S.sigmaB)}</td></tr>
          <tr><td>Groups (G), Δβ</td><td class="num">${S.nGroups} , ${fmt(S.groupDeltaBeta)}</td></tr>
          <tr><td>Selection (ρ)</td><td class="num">${fmt(S.rhoSel)} ${S.rhoSel > 0.5 ? '<span class="badge badge-warning">Bias likely</span>' : '<span class="badge badge-success">Low</span>'}</td></tr>
          <tr><td>Time trend</td><td class="num">${S.timeTrendOn ? fmt(S.timeTrendSlope) + '/wave' : 'Off'}</td></tr>
          <tr><td>Measurement error</td><td class="num">${fmt(S.sigmaME)}</td></tr>
        </table>
      </div>
    </div>
  `;

  document.getElementById("summaryContent").innerHTML = html;
}

function buildDesign(rows, opts) {
  var includeIntercept = !!opts.includeIntercept;
  var includeTrend = !!opts.includeTrend;
  var includeQuad = !!opts.includeQuad;
  var includeGroupFE = !!opts.includeGroupFE;
  var X = new Array(rows.length);
  var y = new Array(rows.length);
  var ids = new Array(rows.length);
  var coefNames = [];
  if (includeIntercept) coefNames.push("Intercept");
  coefNames.push("D");
  if (includeTrend) coefNames.push("t");
  if (includeQuad) coefNames.push("t2");
  if (includeGroupFE && S.nGroups > 1) {
    for (var g = 2; g <= S.nGroups; g++) coefNames.push("G" + g);
  }

  for (var i = 0; i < rows.length; i++) {
    var r = rows[i];
    y[i] = r.y;
    ids[i] = r.pid;
    var row = [];
    if (includeIntercept) row.push(1);
    row.push(r.D);
    if (includeTrend) row.push(r.wave);
    if (includeQuad) row.push(r.wave * r.wave);
    if (includeGroupFE && S.nGroups > 1) {
      for (g = 2; g <= S.nGroups; g++) row.push(r.gid === g ? 1 : 0);
    }
    X[i] = row;
  }
  return { X, y, ids, coefNames };
}

// ============== MAIN RENDER ==============
function computeAndRender() {
  updateLabels();
  var sim = simulate();

  var cross = sim.filter(r => r.wave === S.scatterWave);
  var pooled = sim;

  var seType = S.seType || "classical";
  var csDsgn = buildDesign(cross, { includeIntercept: true, includeTrend: false, includeQuad: false, includeGroupFE: S.specGroupFE });
  var pooledDsgn = buildDesign(pooled, { includeIntercept: true, includeTrend: S.specTrend, includeQuad: S.specQuad, includeGroupFE: S.specGroupFE });
  var feDsgn = buildDesign(pooled, { includeIntercept: false, includeTrend: S.specTrend, includeQuad: S.specQuad, includeGroupFE: false });

  var cs = ols(csDsgn.X, csDsgn.y, { seType, clusterIds: csDsgn.ids, coefNames: csDsgn.coefNames });
  var pooledOLS = ols(pooledDsgn.X, pooledDsgn.y, { seType, clusterIds: pooledDsgn.ids, coefNames: pooledDsgn.coefNames });
  var fe = feWithin(feDsgn.ids, feDsgn.X, feDsgn.y, { seType, coefNames: feDsgn.coefNames });
  var re = randomEffects(pooledDsgn.ids, pooledDsgn.X, pooledDsgn.y, { seType, coefNames: pooledDsgn.coefNames });

  // Store for diagnostics
  window._lastResults = { sim, cross, pooled, cs, pooledOLS, fe, re };

  // Build summary
  buildSummary(sim, cs, pooledOLS, fe, re);

  // Tab rendering
  if (S.viz === "descr") {
    document.getElementById("descrArea").style.display = "block";
    document.getElementById("estArea").style.display = "none";
    document.getElementById("diagArea").style.display = "none";
    document.getElementById("learnArea").style.display = "none";

    drawScatter(document.getElementById("chart-scatter"), cross.map(r => ({ x: r.D, y: r.y, gid: r.gid })), cs, { colorByGroup: S.colorByGroup, nGroups: S.nGroups });

    // Spaghetti plot
    var ids = [], seen = {};
    for (var i = 0; i < pooled.length; i++) {
      var pid = pooled[i].pid;
      if (!seen[pid]) { seen[pid] = 1; ids.push(pid); }
    }
    var k = Math.min(Math.max(5, Math.floor(S.spaghettiN || 20)), ids.length);
    var rand = mulberry32(1234 + S.seed);
    var chosen = {}, arr = [];
    while (Object.keys(chosen).length < k && ids.length > 0) {
      var idx = Math.floor(rand() * ids.length);
      if (!chosen[idx]) { chosen[idx] = 1; arr.push(ids[idx]); }
    }
    var spaghetti = [], summary = [];
    for (i = 0; i < arr.length; i++) {
      var id = arr[i];
      var rows = pooled.filter(r => r.pid === id);
      var ser = rows.map(r => ({ t: r.wave, y: r.y }));
      spaghetti.push({ id, gid: rows.length ? rows[0].gid : 1, series: ser });
    }
    for (var t = 1; t <= S.T; t++) {
      var vals = pooled.filter(r => r.wave === t).map(r => r.y);
      if (vals.length) {
        summary.push({
          t,
          mean: mean(vals),
          median: quantile(vals, 0.5),
          p25: quantile(vals, 0.25),
          p75: quantile(vals, 0.75)
        });
      }
    }
    drawLines(document.getElementById("chart-lines"), spaghetti, summary, S.T, { colorByGroup: S.colorByGroup, nGroups: S.nGroups });

    drawWaveMeans(document.getElementById("chart-wave-means"), pooled, S.T);
    drawBetaDist(document.getElementById("chart-beta-dist"), pooled, { betaMean: S.betaTrue });
  }

  if (S.viz === "est") {
    document.getElementById("descrArea").style.display = "none";
    document.getElementById("estArea").style.display = "block";
    document.getElementById("diagArea").style.display = "none";
    document.getElementById("learnArea").style.display = "none";

    var csD = coefAt(cs, "D");
    var pooledD = coefAt(pooledOLS, "D");
    var feD = coefAt(fe, "D");
    var reD = coefAt(re, "D");
    var estRows = [
      { name: "Cross-Section", b: csD.b, se: csD.se, color: "#60a5fa" },
      { name: "Pooled OLS", b: pooledD.b, se: pooledD.se, color: "#818cf8" },
      { name: "Fixed Effects", b: feD.b, se: feD.se, color: "#34d399" },
      { name: "Random Effects", b: reD.b, se: reD.se, color: "#fbbf24" }
    ];
    drawBarsWithCI(document.getElementById("chart-bars"), estRows, S.betaTrue);

    // Variance decomposition table
    var iccTrue = (S.sigmaU * S.sigmaU) / (S.sigmaU * S.sigmaU + S.sigmaE * S.sigmaE + S.sigmaME * S.sigmaME);
    var vhtml = `
      <tr><td>Between-person variance (σ²<sub>u</sub>)</td><td class="num">${fmt(S.sigmaU * S.sigmaU)}</td><td class="muted">Individual heterogeneity</td></tr>
      <tr><td>Within-person variance (σ²<sub>ε</sub>)</td><td class="num">${fmt(S.sigmaE * S.sigmaE)}</td><td class="muted">Idiosyncratic shocks</td></tr>
      <tr><td>Measurement error (σ²<sub>me</sub>)</td><td class="num">${fmt(S.sigmaME * S.sigmaME)}</td><td class="muted">Adds noise, not bias</td></tr>
      <tr><td><strong>ICC (intraclass correlation)</strong></td><td class="num"><strong>${fmt(iccTrue)}</strong></td><td class="muted">σ²<sub>u</sub> / (σ²<sub>u</sub> + σ²<sub>ε</sub> + σ²<sub>me</sub>)</td></tr>
    `;
    document.getElementById("varianceTable").innerHTML = vhtml;

    // Hausman test
    var haus = hausmanTest(fe, re);
    var hausHtml = `
      <table class="stat-table">
        <tr><td>H statistic</td><td class="num">${fmt(haus.stat, 2)}</td></tr>
        <tr><td>p-value</td><td class="num">${fmt(haus.pval, 4)}</td></tr>
      </table>
      <p class="help-text mt-2">${haus.conclusion}</p>
    `;
    document.getElementById("hausmanResult").innerHTML = hausHtml;

    // Model details
    fillModelDetails(S.modelDetail, cs, pooledOLS, fe, re);

    // Bias alert
    var csBias = Math.abs(csD.b - S.betaTrue);
    var feBias = Math.abs(feD.b - S.betaTrue);
    var alertHtml = "";
    if (S.rhoSel > 0.3 && csBias > 0.1 && feBias < csBias * 0.5) {
      alertHtml = `<div class="alert alert-warning">
        <strong>Selection bias detected!</strong> Cross-section and Pooled OLS show substantial bias
        (${fmt(csBias, 2)} from true β) while Fixed Effects is closer (${fmt(feBias, 2)} from true β).
        This is because ρ = ${fmt(S.rhoSel)} creates correlation between treatment and unobserved individual effects.
      </div>`;
    } else if (S.rhoSel < 0.1 && csBias < 0.15) {
      alertHtml = `<div class="alert alert-success">
        <strong>No substantial bias.</strong> With low selection (ρ = ${fmt(S.rhoSel)}), all estimators
        perform similarly. This is similar to a randomized experiment.
      </div>`;
    }
    document.getElementById("biasAlert").innerHTML = alertHtml;

    drawSeSensitivity(document.getElementById("chart-se-sensitivity"), S.modelDetail, cross, pooled, { cs, pooledOLS, fe, re });
    drawGroupSplit(document.getElementById("chart-group-split"), pooled, { estimator: S.groupEstimator, seType, includeTrend: S.specTrend, includeQuad: S.specQuad, includeGroupFE: S.specGroupFE });
  }

  if (S.viz === "diag") {
    document.getElementById("descrArea").style.display = "none";
    document.getElementById("estArea").style.display = "none";
    document.getElementById("diagArea").style.display = "block";
    document.getElementById("learnArea").style.display = "none";

    var which = S.diagModel || "pooled";
    var points = [];
    if (which === "pooled") {
      points = pooledOLS.fitted.map((f, i) => ({ fit: f, resid: pooledOLS.resid[i] }));
    } else if (which === "cs") {
      points = cs.fitted.map((f, i) => ({ fit: f, resid: cs.resid[i] }));
    } else if (which === "fe") {
      points = fe.fitted.map((f, i) => ({ fit: f, resid: fe.resid[i] }));
    } else if (which === "re") {
      points = re.fitted.map((f, i) => ({ fit: f, resid: re.resid[i] }));
    }

    drawResiduals(document.getElementById("chart-rvf"), points);
    drawHist(document.getElementById("chart-hist"), points.map(p => p.resid));
    drawQQ(document.getElementById("chart-qq"), points.map(p => p.resid));

    // Diagnostic interpretation
    var residMean = mean(points.map(p => p.resid));
    var residVar = variance(points.map(p => p.resid));
    var interpHtml = `
      <p><strong>Residual mean:</strong> ${fmt(residMean, 4)} (should be ≈ 0)</p>
      <p><strong>Residual variance:</strong> ${fmt(residVar, 3)}</p>
      <p class="mt-2">The residuals appear ${Math.abs(residMean) < 0.01 ? "well-centered" : "slightly off-center"}.
      Check the Q-Q plot for normality—deviations in the tails suggest heavy-tailed or skewed errors.</p>
    `;
    document.getElementById("diagInterpretation").innerHTML = interpHtml;
  }

  if (S.viz === "learn") {
    document.getElementById("descrArea").style.display = "none";
    document.getElementById("estArea").style.display = "none";
    document.getElementById("diagArea").style.display = "none";
    document.getElementById("learnArea").style.display = "block";
  }

  // Cards
  var cards = [];
  function addCard(title, value, ciVals, foot, highlight) {
    cards.push({ title, value, ciVals, foot, highlight });
  }
  var csD2 = coefAt(cs, "D");
  var pooledD2 = coefAt(pooledOLS, "D");
  var feD2 = coefAt(fe, "D");
  var reD2 = coefAt(re, "D");
  var csCI = ci(csD2.b, csD2.se);
  var pooledCI = ci(pooledD2.b, pooledD2.se);
  var feCI = ci(feD2.b, feD2.se);
  var reCI = ci(reD2.b, reD2.se);

  addCard("True β", fmt(S.betaTrue), "", "Ground truth", true);
  addCard("Cross-Section", fmt(csD2.b), `[${fmt(csCI[0], 2)}, ${fmt(csCI[1], 2)}]`, `R²=${fmt(cs.r2, 2)}`);
  addCard("Pooled OLS", fmt(pooledD2.b), `[${fmt(pooledCI[0], 2)}, ${fmt(pooledCI[1], 2)}]`, `R²=${fmt(pooledOLS.r2, 2)}`);
  addCard("Fixed Effects", fmt(feD2.b), `[${fmt(feCI[0], 2)}, ${fmt(feCI[1], 2)}]`, `within R²=${fmt(fe.r2, 2)}`);
  addCard("Random Effects", fmt(reD2.b), `[${fmt(reCI[0], 2)}, ${fmt(reCI[1], 2)}]`, `θ=${fmt(re.theta, 2)}`);

  var cardsDiv = document.getElementById("cards");
  cardsDiv.innerHTML = "";
  for (i = 0; i < cards.length; i++) {
    var c = cards[i];
    var el = document.createElement("div");
    el.className = "card" + (c.highlight ? " highlight" : "");
    var inner = `<div class="title">${c.title}</div><div class="value">${c.value}</div>`;
    if (c.ciVals) inner += `<div class="muted" style="font-size:11px;margin-top:4px">95% CI: ${c.ciVals}</div>`;
    if (c.foot) inner += `<div class="foot">${c.foot}</div>`;
    el.innerHTML = inner;
    cardsDiv.appendChild(el);
  }
}

// ============== UI BINDINGS ==============
function bindRange(id, key) {
  var el = document.getElementById(id);
  el.value = S[key];
  el.addEventListener("input", function(e) {
    var v = parseFloat(e.target.value);
    if (key === "N" || key === "T" || key === "nGroups" || key === "scatterWave" || key === "spaghettiN") v = Math.round(v);
    S[key] = v;
    computeAndRender();
  });
}

function bindCheck(id, key) {
  var el = document.getElementById(id);
  el.checked = S[key];
  el.addEventListener("change", function(e) {
    S[key] = !!e.target.checked;
    computeAndRender();
  });
}

function bindSelect(id, key) {
  var el = document.getElementById(id);
  el.value = S[key];
  el.addEventListener("change", function(e) {
    S[key] = String(e.target.value);
    computeAndRender();
  });
}

function setViz(v) {
  S.viz = v;
  ["viz-descr", "viz-estimates", "viz-diagnostics", "viz-learn"].forEach(id => {
    document.getElementById(id).classList.remove("active");
  });
  if (v === "descr") document.getElementById("viz-descr").classList.add("active");
  if (v === "est") document.getElementById("viz-estimates").classList.add("active");
  if (v === "diag") document.getElementById("viz-diagnostics").classList.add("active");
  if (v === "learn") document.getElementById("viz-learn").classList.add("active");
  computeAndRender();
}

function markMD(which) {
  ["md-cs", "md-pooled", "md-fe", "md-re"].forEach(id => {
    document.getElementById(id).classList.remove("active");
  });
  document.getElementById("md-" + which).classList.add("active");
  S.modelDetail = which;
  var r = window._lastResults;
  if (r) fillModelDetails(which, r.cs, r.pooledOLS, r.fe, r.re);
}

function markDG(which) {
  ["dg-cs", "dg-pooled", "dg-fe", "dg-re"].forEach(id => {
    document.getElementById(id).classList.remove("active");
  });
  document.getElementById("dg-" + which).classList.add("active");
  S.diagModel = which;
  computeAndRender();
}

function bindUI() {
  // Tab navigation
  document.getElementById("viz-descr").addEventListener("click", () => setViz("descr"));
  document.getElementById("viz-estimates").addEventListener("click", () => setViz("est"));
  document.getElementById("viz-diagnostics").addEventListener("click", () => setViz("diag"));
  document.getElementById("viz-learn").addEventListener("click", () => setViz("learn"));

  // Parameter sliders
  bindRange("N", "N");
  bindRange("T", "T");
  bindRange("betaTrue", "betaTrue");
  bindRange("sigmaB", "sigmaB");
  bindRange("nGroups", "nGroups");
  bindRange("groupDeltaBeta", "groupDeltaBeta");
  bindRange("sigmaU", "sigmaU");
  bindRange("sigmaE", "sigmaE");
  bindRange("sigmaME", "sigmaME");
  bindRange("pBase", "pBase");
  bindRange("rhoSel", "rhoSel");
  bindRange("switchRate", "switchRate");
  bindRange("timeTrendSlope", "timeTrendSlope");
  bindCheck("timeTrendOn", "timeTrendOn");
  bindRange("scatterWave", "scatterWave");
  bindRange("spaghettiN", "spaghettiN");
  bindCheck("colorByGroup", "colorByGroup");
  bindSelect("seType", "seType");
  bindCheck("specTrend", "specTrend");
  bindCheck("specQuad", "specQuad");
  bindCheck("specGroupFE", "specGroupFE");
  bindSelect("groupEstimator", "groupEstimator");

  // Presets
  document.getElementById("preset-baseline").addEventListener("click", () => applyPreset("Baseline (No Bias)"));
  document.getElementById("preset-selection").addEventListener("click", () => applyPreset("Selection Bias"));
  document.getElementById("preset-rct").addEventListener("click", () => applyPreset("RCT-like"));
  document.getElementById("preset-noisy").addEventListener("click", () => applyPreset("Noisy Survey"));
  document.getElementById("preset-trend").addEventListener("click", () => applyPreset("Time Trend"));
  document.getElementById("preset-lowvar").addEventListener("click", () => applyPreset("Low Within Variation"));

  // New sample
  document.getElementById("btn-newsample").addEventListener("click", function() {
    S.seed = Math.floor(Math.random() * 100000);
    computeAndRender();
  });

  // Download CSV
  document.getElementById("btn-download").addEventListener("click", function() {
    var sim = simulate();
    var header = ["pid", "gid", "beta_i", "wave", "D", "u", "e", "yStar", "y", "trend"];
    function esc(s) { return '"' + String(s).replace(/"/g, '""') + '"'; }
    var csv = [header.join(",")].concat(
      sim.map(r => header.map(h => esc(r[h])).join(","))
    ).join("\n");
    var blob = new Blob([csv], { type: 'text/csv;charset=utf-8;' });
    var url = URL.createObjectURL(blob);
    var a = document.createElement("a");
    a.href = url;
    a.download = "panel_simulation.csv";
    a.click();
    URL.revokeObjectURL(url);
  });

  // Model detail chips
  document.getElementById("md-cs").addEventListener("click", () => markMD("cs"));
  document.getElementById("md-pooled").addEventListener("click", () => markMD("pooled"));
  document.getElementById("md-fe").addEventListener("click", () => markMD("fe"));
  document.getElementById("md-re").addEventListener("click", () => markMD("re"));

  // Diagnostic model chips
  document.getElementById("dg-cs").addEventListener("click", () => markDG("cs"));
  document.getElementById("dg-pooled").addEventListener("click", () => markDG("pooled"));
  document.getElementById("dg-fe").addEventListener("click", () => markDG("fe"));
  document.getElementById("dg-re").addEventListener("click", () => markDG("re"));

  // Collapsible sections
  document.querySelectorAll(".collapsible").forEach(el => {
    el.addEventListener("click", function() {
      this.classList.toggle("open");
      var contentId = this.id.replace("coll-", "content-");
      document.getElementById(contentId).classList.toggle("show");
    });
  });
}

function init() {
  syncInputs();
  bindUI();
  computeAndRender();
}

(function(root) {
  try {
    root.StatsLab = Object.assign(root.StatsLab || {}, {
      ols,
      feWithin,
      randomEffects,
      hausmanTest,
      coefAt,
      coefIndex,
      ci,
      pnorm,
      pchisq,
      matInverse,
      matMul,
      matTranspose,
      matVecMul
    });
  } catch (_) {
    // no-op
  }
})(typeof globalThis !== "undefined" ? globalThis : this);

if (typeof window !== "undefined" && window.addEventListener) {
  window.addEventListener("DOMContentLoaded", init);
}
