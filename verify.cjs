/*
 * verify.cjs — checks the JS engine against canonical plm estimates (R).
 * Run: node verify.cjs
 * Canonical values produced by export_wagepan.R using plm (Croissant & Millo).
 */
const fs = require("fs");
const path = require("path");
const E = require("./engine.js");

const rows = JSON.parse(fs.readFileSync(path.join(__dirname, "data", "wagepan.json"), "utf8"));

// meta (reference categories = most frequent)
function cats(key) {
  const c = {};
  rows.forEach(r => { c[r[key]] = (c[r[key]] || 0) + 1; });
  return Object.keys(c).sort((a, b) => c[b] - c[a]);
}
const meta = { industryCats: cats("industry"), regionCats: cats("region"),
  years: Array.from(new Set(rows.map(r => r.year))).sort((a, b) => a - b) };

// canonical plm values: [beta, se] (se null where not checked); "DROP" = unidentified
const CANON = {
  union:   { pooled: [0.1673, 0.0181], between: [0.2567, 0.0503], random: [0.1017, null], within: [0.0833, 0.0193], fd: [0.0429, null] },
  married: { pooled: [0.1635, 0.0164], between: [0.2198, 0.0438], random: [0.0766, null], within: [0.0473, 0.0183], fd: [0.0384, null] },
  educ:    { pooled: [0.1021, 0.0047], between: [0.0986, 0.0112], random: [0.1029, null], within: "DROP", fd: "DROP" }
};

const baseControls = { experience: true, hours: false, industry: false, region: false, health: false, race: false };
let pass = 0, fail = 0;
function check(label, got, want, tol) {
  const ok = Math.abs(got - want) <= tol;
  console.log(`${ok ? "  ok  " : " FAIL "} ${label}: got ${got.toFixed(4)}  want ${want.toFixed(4)}  (tol ${tol})`);
  ok ? pass++ : fail++;
}

["union", "married", "educ"].forEach(focal => {
  console.log(`\n=== focal: ${focal} (+ exper + expersq, full sample, classical SE) ===`);
  ["pooled", "between", "random", "within", "fd"].forEach(est => {
    const spec = { focal, outcome: "lwage", controls: baseControls, yearFE: false, sample: "full", estimator: est, seType: "classical" };
    const res = E.runSpec(rows, spec, meta);
    const want = CANON[focal][est];
    if (want === "DROP") {
      const ok = !res.identified;
      console.log(`${ok ? "  ok  " : " FAIL "} ${est}: expected UNIDENTIFIED, got ${res.identified ? "b=" + res.b.toFixed(4) : "unidentified"}`);
      ok ? pass++ : fail++;
      return;
    }
    const betaTol = est === "random" ? 0.006 : 0.0015;
    check(`${est} beta`, res.b, want[0], betaTol);
    if (want[1] != null) check(`${est} se`, res.se, want[1], est === "between" ? 0.004 : 0.0015);
  });
});

// ---------- additional canonical-numbers tests ----------
// These benchmarks cover the parts of the engine that the original 23 tests left
// unexercised: the Mundlak equivalence (CRE focal coefficient = FE focal
// coefficient), two-way fixed effects (against plm effect="twoways"), HC1 and
// cluster-robust (CR1) standard errors (against sandwich::vcovHC), and the
// joint Hausman statistic (against plm::phtest).
function checkPass(label, cond, detail) {
  console.log(`${cond ? "  ok  " : " FAIL "} ${label}${detail ? "  — " + detail : ""}`);
  cond ? pass++ : fail++;
}

const specBase = (focal, est, seType, yearFE) => ({
  focal, outcome: "lwage", controls: baseControls,
  yearFE: !!yearFE, sample: "full", estimator: est,
  seType: seType || "classical"
});

console.log("\n=== Mundlak equivalence: CRE focal coefficient == FE focal coefficient ===");
for (const focal of ["union", "married"]) {
  const fe  = E.runSpec(rows, specBase(focal, "within"), meta);
  const cre = E.runSpec(rows, specBase(focal, "cre"),    meta);
  const diff = Math.abs(fe.b - cre.b);
  checkPass(
    `${focal} CRE = FE`,
    diff < 1e-4,
    `FE=${fe.b.toFixed(5)}  CRE=${cre.b.toFixed(5)}  |Δ|=${diff.toExponential(2)}`
  );
}

console.log("\n=== Two-way FE (entity + time): plm effect=\"twoways\" ===");
const twfe   = E.runSpec(rows, specBase("union", "twfe"),   meta);
const wYearFE = E.runSpec(rows, specBase("union", "within", "classical", true), meta);
// plm two-way FE union = 0.08130
check("TWFE union beta vs plm canonical", twfe.b, 0.08130, 0.0015);
checkPass(
  "TWFE = within + yearFE (engine self-equivalence)",
  Math.abs(twfe.b - wYearFE.b) < 1e-9,
  `TWFE=${twfe.b.toFixed(6)}  within+yearFE=${wYearFE.b.toFixed(6)}`
);

console.log("\n=== Heteroskedasticity-robust (HC1) and cluster-robust (CR1) SEs ===");
// Canonical values from R's sandwich package on the matching base-lm model
// (lm(lwage ~ union + exper + expersq) with sandwich::vcovHC(type='HC1') and
// sandwich::vcovCL(cluster=~nr, type='HC1')). Engine matches sandwich to 1e-5
// on pooled OLS.
const pooledHC1 = E.runSpec(rows, specBase("union", "pooled", "robust"), meta);
check("pooled union HC1 SE  (vs sandwich::vcovHC)",  pooledHC1.se, 0.01697, 0.0001);
const pooledCR1 = E.runSpec(rows, specBase("union", "pooled", "cluster"), meta);
check("pooled union CR1 SE  (vs sandwich::vcovCL)",  pooledCR1.se, 0.02997, 0.0001);
// Within-transformed regressions: the engine uses the Stata xtreg cluster
// small-sample correction (df = N(T-1) - k), while sandwich's vcovHC/vcovCL
// applied to the demeaned regression uses the lm convention (df = n - k_slopes).
// Both are defensible; we benchmark the engine's value with a documented
// tolerance so the test fails only on a substantive change.
const withinCR1 = E.runSpec(rows, specBase("union", "within", "cluster"), meta);
check("within union CR1 SE  (engine Stata-xtreg convention)", withinCR1.se, 0.02455, 0.0010);

console.log("\n=== Joint Hausman test (df = 3, against plm::phtest) ===");
const mkClassical = (est) => E.runSpec(rows, specBase("union", est, "classical"), meta).model;
const hjoint = E.hausmanFull(mkClassical("within"), mkClassical("random"));
checkPass(
  "Hausman df = 3",
  hjoint.df === 3,
  `df=${hjoint.df}`
);
check("Hausman chisq vs plm phtest", hjoint.stat, 139.94, 0.05);
checkPass(
  "Hausman p-value < 1e-20",
  hjoint.pval < 1e-20,
  `p = ${hjoint.pval.toExponential(2)}`
);

console.log(`\n${fail === 0 ? "ALL PASS" : "SOME FAILED"} — ${pass} passed, ${fail} failed`);
process.exit(fail === 0 ? 0 : 1);
