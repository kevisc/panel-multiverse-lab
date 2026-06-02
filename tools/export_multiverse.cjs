#!/usr/bin/env node
/* Headless multiverse export.
 *
 * Loads the bundled wage-panel extract and drives the JavaScript estimation
 * engine in Node to enumerate the full specification curve of each focal
 * estimand, then writes two CSV artefacts the downstream R / Python scripts
 * consume:
 *
 *   verification/js_multiverse.csv   the full grid (3 focals × 5 estimators
 *                                    × 16 control combinations × 2 year-FE
 *                                    options × 3 sample definitions = 1,440
 *                                    specifications). Schema:
 *                                      focal, estimator, experience, industry,
 *                                      region, health, yearFE, sample,
 *                                      identified, b, se
 *
 *   figures/multiverse_union.csv     the union-premium subset sorted by
 *                                    estimate, used by figures/make_figure.R.
 *                                    Schema:
 *                                      rank, b, lo, hi, sig, estimator,
 *                                      experience, industry, region, health,
 *                                      yearFE
 *
 * Run from anywhere; paths are resolved relative to this script's location,
 * so the script remains correct after the directory is moved or symlinked.
 *
 * Standard errors used here are cluster-robust by person (CR1), matching the
 * intervals reported in Figure 2 of the paper. The selection of estimators,
 * controls, year-FE and sample axes matches the engine's default multiverse
 * (which is also what the live application enumerates when the user clicks
 * "Build the multiverse" with its initial selections); it is the same grid
 * that variance_decomp.R decomposes into Table 1.
 */
'use strict';

const fs   = require('fs');
const path = require('path');

const ROOT = path.resolve(__dirname, '..');
const E    = require(path.join(ROOT, 'engine.js'));

const DATA_PATH       = path.join(ROOT, 'data', 'wagepan.json');
const GRID_CSV        = path.join(ROOT, 'verification', 'js_multiverse.csv');
const FIGURE_CSV      = path.join(ROOT, 'figures', 'multiverse_union.csv');

const FOCALS          = ['union', 'married', 'educ'];
const ESTIMATORS      = ['pooled', 'between', 'random', 'within', 'fd'];
const CONTROLS_VARY   = ['experience', 'industry', 'region', 'health'];
const SAMPLES         = ['full', 'nohealth', 'trim'];

function buildMeta(rows) {
  const counts = (key) => {
    const c = {};
    for (const r of rows) c[r[key]] = (c[r[key]] || 0) + 1;
    return Object.keys(c).sort((a, b) => c[b] - c[a]);
  };
  return {
    industryCats: counts('industry'),
    regionCats:   counts('region'),
    years:        [...new Set(rows.map((r) => r.year))].sort((a, b) => a - b)
  };
}

function ensureDir(file) {
  const dir = path.dirname(file);
  if (!fs.existsSync(dir)) fs.mkdirSync(dir, { recursive: true });
}

function fmt(x, digits) {
  return Number.isFinite(x) ? x.toFixed(digits) : 'NA';
}

function main() {
  if (!fs.existsSync(DATA_PATH)) {
    console.error('error: ' + DATA_PATH + ' is missing. Run Rscript data-prep.R first.');
    process.exit(1);
  }
  const rows = JSON.parse(fs.readFileSync(DATA_PATH, 'utf8'));
  const meta = buildMeta(rows);

  // ----- full grid → verification/js_multiverse.csv (3 estimands × 480 specs) -----
  const gridHeader = ['focal','estimator','experience','industry','region','health','yearFE','sample','identified','b','se'];
  const gridLines  = [gridHeader.join(',')];
  let nIdentified = 0;
  let nTotal      = 0;

  for (const focal of FOCALS) {
    const axes = {
      estimators:     ESTIMATORS.slice(),
      controlsVary:   CONTROLS_VARY.slice(),
      controlsFixed:  {},
      yearFEVary:     true,
      yearFEFixed:    false,
      sampleVary:     SAMPLES.slice(),
      focal:          focal,
      outcome:        'lwage',
      seType:         'cluster'
    };
    const results = E.enumerateMultiverse(rows, axes, meta);
    nTotal += results.length;
    for (const r of results) {
      if (r.identified) nIdentified++;
      const c = r.spec.controls || {};
      gridLines.push([
        focal,
        r.estimator,
        c.experience ? 1 : 0,
        c.industry   ? 1 : 0,
        c.region     ? 1 : 0,
        c.health     ? 1 : 0,
        r.spec.yearFE ? 1 : 0,
        r.spec.sample,
        r.identified ? 1 : 0,
        r.identified ? fmt(r.b, 6) : 'NA',
        r.identified && Number.isFinite(r.se) ? fmt(r.se, 6) : 'NA'
      ].join(','));
    }
  }
  ensureDir(GRID_CSV);
  fs.writeFileSync(GRID_CSV, gridLines.join('\n') + '\n');
  console.log('wrote ' + path.relative(ROOT, GRID_CSV) + '  (' + (gridLines.length - 1) + ' specs; ' + nIdentified + ' identified, ' + (nTotal - nIdentified) + ' not identified)');

  // ----- union subset → figures/multiverse_union.csv (sorted by estimate) -----
  const unionAxes = {
    estimators:    ESTIMATORS.slice(),
    controlsVary:  CONTROLS_VARY.slice(),
    controlsFixed: {},
    yearFEVary:    true,
    yearFEFixed:   false,
    sampleVary:    ['full'],
    focal:         'union',
    outcome:       'lwage',
    seType:        'cluster'
  };
  const unionResults = E.enumerateMultiverse(rows, unionAxes, meta)
    .filter((r) => r.identified && Number.isFinite(r.b))
    .sort((a, b) => a.b - b.b);

  const figHeader = ['rank','b','lo','hi','sig','estimator','experience','industry','region','health','yearFE'];
  const figLines  = [figHeader.join(',')];
  unionResults.forEach((r, idx) => {
    const lo = r.b - 1.96 * r.se;
    const hi = r.b + 1.96 * r.se;
    const sig = (lo > 0 || hi < 0) ? 1 : 0;
    const c = r.spec.controls || {};
    figLines.push([
      idx + 1,
      fmt(r.b, 5),
      fmt(lo, 5),
      fmt(hi, 5),
      sig,
      r.estimator,
      c.experience ? 1 : 0,
      c.industry   ? 1 : 0,
      c.region     ? 1 : 0,
      c.health     ? 1 : 0,
      r.spec.yearFE ? 1 : 0
    ].join(','));
  });
  ensureDir(FIGURE_CSV);
  fs.writeFileSync(FIGURE_CSV, figLines.join('\n') + '\n');
  console.log('wrote ' + path.relative(ROOT, FIGURE_CSV) + '  (' + unionResults.length + ' identified specifications, sorted)');
}

main();
