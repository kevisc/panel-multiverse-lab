ObjC.import("Foundation");
ObjC.import("stdlib");

function readText(path) {
  var p = $(path).stringByStandardizingPath;
  var data = $.NSFileManager.defaultManager.contentsAtPath(p);
  if (!data) throw new Error("Could not read file: " + path);
  return ObjC.unwrap($.NSString.alloc.initWithDataEncoding(data, $.NSUTF8StringEncoding));
}

function approxEqual(a, b, tol) {
  return Math.abs(a - b) <= tol;
}

function assert(cond, msg) {
  if (!cond) throw new Error(msg || "Assertion failed");
}

function olsClosedFormInterceptX(x, y) {
  var n = x.length;
  var sumx = 0,
    sumy = 0,
    sumx2 = 0,
    sumxy = 0;
  for (var i = 0; i < n; i++) {
    sumx += x[i];
    sumy += y[i];
    sumx2 += x[i] * x[i];
    sumxy += x[i] * y[i];
  }
  var det = n * sumx2 - sumx * sumx;
  assert(Math.abs(det) > 1e-12, "Singular XtX in closed-form OLS");
  var b0 = (sumx2 * sumy - sumx * sumxy) / det;
  var b1 = (-sumx * sumy + n * sumxy) / det;
  return { b0: b0, b1: b1, n: n, sumx: sumx, sumx2: sumx2 };
}

function invXtXInterceptX(n, sumx, sumx2) {
  var det = n * sumx2 - sumx * sumx;
  assert(Math.abs(det) > 1e-12, "Singular XtX inverse");
  var a00 = sumx2 / det;
  var a01 = -sumx / det;
  var a11 = n / det;
  return [
    [a00, a01],
    [a01, a11]
  ];
}

function robustVcovInterceptX(x, resid, invXtX, df) {
  var meat00 = 0,
    meat01 = 0,
    meat11 = 0;
  for (var i = 0; i < x.length; i++) {
    var u2 = resid[i] * resid[i];
    meat00 += u2;
    meat01 += u2 * x[i];
    meat11 += u2 * x[i] * x[i];
  }
  // inv * meat * inv (2x2)
  var A = invXtX;
  var B00 = meat00,
    B01 = meat01,
    B11 = meat11;
  var tmp00 = A[0][0] * B00 + A[0][1] * B01;
  var tmp01 = A[0][0] * B01 + A[0][1] * B11;
  var tmp10 = A[1][0] * B00 + A[1][1] * B01;
  var tmp11 = A[1][0] * B01 + A[1][1] * B11;
  var V00 = tmp00 * A[0][0] + tmp01 * A[0][1];
  var V01 = tmp00 * A[1][0] + tmp01 * A[1][1];
  var V11 = tmp10 * A[1][0] + tmp11 * A[1][1];
  var n = x.length;
  var hc1 = n / df;
  return [
    [V00 * hc1, V01 * hc1],
    [V01 * hc1, V11 * hc1]
  ];
}

function clusterVcovInterceptX(x, resid, clusters, invXtX, df) {
  var by = {};
  for (var i = 0; i < x.length; i++) {
    var c = String(clusters[i]);
    if (!by[c]) by[c] = [];
    by[c].push(i);
  }
  var keys = Object.keys(by);
  var G = keys.length;
  assert(G >= 2, "Need >=2 clusters");

  var meat00 = 0,
    meat01 = 0,
    meat11 = 0;
  for (var g = 0; g < keys.length; g++) {
    var idxs = by[keys[g]];
    var s0 = 0,
      s1 = 0;
    for (var j = 0; j < idxs.length; j++) {
      var ii = idxs[j];
      s0 += resid[ii];
      s1 += x[ii] * resid[ii];
    }
    meat00 += s0 * s0;
    meat01 += s0 * s1;
    meat11 += s1 * s1;
  }
  var A = invXtX;
  var B00 = meat00,
    B01 = meat01,
    B11 = meat11;
  var tmp00 = A[0][0] * B00 + A[0][1] * B01;
  var tmp01 = A[0][0] * B01 + A[0][1] * B11;
  var tmp10 = A[1][0] * B00 + A[1][1] * B01;
  var tmp11 = A[1][0] * B01 + A[1][1] * B11;
  var V00 = tmp00 * A[0][0] + tmp01 * A[0][1];
  var V01 = tmp00 * A[1][0] + tmp01 * A[1][1];
  var V11 = tmp10 * A[1][0] + tmp11 * A[1][1];

  var n = x.length;
  var adj = (G / (G - 1)) * ((n - 1) / df);
  return [
    [V00 * adj, V01 * adj],
    [V01 * adj, V11 * adj]
  ];
}

function main() {
  // Evaluate app.js in this context to get StatsLab exports.
  eval(readText("static/stats/app.js"));
  assert(typeof StatsLab === "object" && typeof StatsLab.ols === "function", "StatsLab exports missing");

  // Test 1: OLS betas (intercept + x) match closed form
  var x = [0, 1, 0, 1, 0, 1];
  var y = [1.0, 2.0, 1.2, 2.5, 0.7, 2.1];
  var X = x.map(v => [1, v]);

  var fit = StatsLab.ols(X, y, { seType: "classical", coefNames: ["Intercept", "x"] });
  var cf = olsClosedFormInterceptX(x, y);
  assert(approxEqual(fit.beta[0], cf.b0, 1e-10), "OLS intercept mismatch");
  assert(approxEqual(fit.beta[1], cf.b1, 1e-10), "OLS slope mismatch");

  // Test 2: HC1 robust vcov matches closed-form sandwich
  var df = fit.df;
  var invXtX = invXtXInterceptX(cf.n, cf.sumx, cf.sumx2);
  var resid = y.map((yi, i) => yi - (cf.b0 + cf.b1 * x[i]));
  var Vh = robustVcovInterceptX(x, resid, invXtX, df);

  var fitR = StatsLab.ols(X, y, { seType: "robust", coefNames: ["Intercept", "x"] });
  assert(approxEqual(fitR.vcov[0][0], Vh[0][0], 1e-10), "HC1 var(intercept) mismatch");
  assert(approxEqual(fitR.vcov[1][1], Vh[1][1], 1e-10), "HC1 var(slope) mismatch");

  // Test 3: Cluster vcov matches closed-form (2 clusters)
  var clusters = [1, 1, 1, 2, 2, 2];
  var Vc = clusterVcovInterceptX(x, resid, clusters, invXtX, df);
  var fitC = StatsLab.ols(X, y, { seType: "cluster", clusterIds: clusters, coefNames: ["Intercept", "x"] });
  assert(approxEqual(fitC.vcov[0][0], Vc[0][0], 1e-10), "Cluster var(intercept) mismatch");
  assert(approxEqual(fitC.vcov[1][1], Vc[1][1], 1e-10), "Cluster var(slope) mismatch");

  // Test 4: FE within recovers known slope
  var ids = [1, 1, 2, 2];
  var D = [0, 1, 0, 1];
  var y2 = [10, 15, -3, 2]; // y = a_i + 5*D
  var X2 = D.map(v => [v]);
  var fe = StatsLab.feWithin(ids, X2, y2, { seType: "classical", coefNames: ["D"] });
  assert(approxEqual(fe.beta[0], 5, 1e-12), "FE slope mismatch");

  console.log("OK: selftest passed");
}

try {
  main();
} catch (e) {
  console.log("FAIL:", String(e && e.message ? e.message : e));
  $.exit(1);
}

