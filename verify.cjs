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

// Hausman (union + exper + expersq)
const mk = est => E.runSpec(rows, { focal: "union", outcome: "lwage", controls: baseControls, yearFE: false, sample: "full", estimator: est, seType: "classical" }, meta).model;
const h = E.hausman(mk("within"), mk("random"), "union");
console.log(`\n=== Hausman (union spec) ===`);
console.log(`  H=${h.stat.toFixed(2)} p=${h.pval.toExponential(2)} (plm: H=139.94 over 3 df; expect strong rejection on the union coef)`);

console.log(`\n${fail === 0 ? "ALL PASS" : "SOME FAILED"} — ${pass} passed, ${fail} failed`);
process.exit(fail === 0 ? 0 : 1);
