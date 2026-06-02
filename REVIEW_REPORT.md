# Pre-submission Review — Panel Data Multiverse Lab

## 1. Executive summary

- **Reproducibility chain has critical holes that will block RSOS/PLOS ONE production**: the CITATION.cff `repository-code` URL is a 404, there is no Zenodo DOI (only a future-tense promise), there is no top-level README.md, `data-prep.R` writes to a hardcoded absolute path on the author's laptop, and the headless step that produces `verification/js_multiverse.csv` exists only as a browser download in `exports.js#multiverseCSV`. Any one of these would draw a desk-reject query; together they are the single biggest pre-submission risk.
- **Inference machinery is not reproducible and uses the wrong null for panel data**: `bootstrapFocal` (engine.js:430) and `permuteFocal` (engine.js:451) call `Math.random()` with no seed plumbed through, so bootstrap CIs and randomization-inference p-values change on every page reload. Worse, `permuteFocal` shuffles the focal regressor flat across all person-years, destroying within-person persistence — for unions/marriage this is conservative; for educ it is qualitatively wrong (educ is essentially time-invariant). Both issues bear directly on the paper's "fully reproducible" framing.
- **Verification harness is narrower than the paper advertises**: `verify.cjs` covers only five of seven estimators (pooled/between/random/within/fd) with classical SEs, experience-only controls, full sample. `cre`, `twfe`, HC1, CR1, two-way Hausman (`hausmanFull`) are entirely unexercised, despite being named in paper.md line 34 and line 42. The Hausman line in verify.cjs prints the 1-df focal stat next to a 3-df plm reference — apples-to-oranges. Several `verify.cjs` lines do not even call `check()` and so cannot fail.
- **Headline arithmetic is internally inconsistent**: the paper claims "seven estimators ... up to 480 identified specifications per estimand ... all 1,440 specifications across the three estimands" (lines 34, 42), but 480 = 5 × 16 × 2 × 3 — the multiverse and CSV use five estimators, not seven; with seven it would be 672 per estimand / 2,016 total. Figure 2's "160 estimates" is the full-sample slice (5 × 16 × 2), not the 480-spec union universe used in Table 1. This will be the first thing a methods referee asks about.
- **Variance decomposition needs three small but visible upgrades**: (i) "marginal contribution" is loose terminology for a Type-I sequential SS share — add Type-III/Shapley as a robustness column; (ii) Table 1 rows sum to 91.3% / 72.0% / 95.7% with no residual reported and no R² — add a Residual column; (iii) the 28% residual on the marriage row is informative (estimator × covariate interaction) and deserves a sentence rather than silence.
- **Abstract overclaims scope**: the prescriptive "Multiverse analyses of longitudinal data should therefore treat the estimator as a first-class forking path" sits at exactly 200 words with no scope qualifier, while the body twice concedes the single-dataset framing (lines 16, 84). RSOS referees will flag the leap from one wagepan to a discipline-wide recommendation.
- **Wins to keep**: cross-language reproduction within ~1e-6 (R) and ~5e-7 (Python); a clean separation of inference (SE family) from point estimates; CRE = FE Mundlak equivalence holds to machine precision (4×10⁻¹⁵); Ethics/Funding/Conflict-of-interest/Authors-contributions blocks are present in the correct RSOS order; no "available from authors on request"; data citation traceable to Vella–Verbeek 1998 via `wooldridge`. The browser tool is genuinely impressive and works.

## 2. Critical issues (must-fix before submission)

### CITATION.cff `repository-code` is a 404 (phantom repo)  — area: paper-reproducibility, **CRITICAL**
**File / location**: `/Users/kevinschoenholzer/Documents/GitHub/kevisc.github.io/sim-paneldata/CITATION.cff` line 15
**Problem**: `repository-code: "https://github.com/kevisc/panel-data-multiverse"` returns HTTP 404 (user `kevisc` exists, the named repo does not). There is currently no public source-code repository for the project, only the hosted build under `kevinschoenholzer.com`.
**Why it matters**: Reproducibility reviewers, Zenodo↔GitHub integration, and GitHub's "cite this repository" all break. The reader cannot find the code that produced the paper's numbers.
**Recommended fix**: Create the public repo at the cited URL (or pick a canonical one and update both `CITATION.cff` line 15 and `paper.md` line 92 to match). Push the full `sim-paneldata/` tree. Verify with `curl -I` before submission.

### No Zenodo DOI — archived release only promised  — area: paper-reproducibility, **CRITICAL**
**File / location**: `paper.md` line 92 (Data accessibility)
**Problem**: "For the final version an archived release of the code and data extract will be deposited in a public repository (Zenodo) and its DOI cited here." No DOI, no structured placeholder, and `CITATION.cff` has no `identifiers:` block.
**Why it matters**: RSOS production will not typeset without a citable DOI, and PLOS ONE expects either a DOI or a persistent identifier at submission. The vague future-tense promise gives the editor no slot to audit.
**Recommended fix**: Use Zenodo's "Reserve DOI" facility on a draft record. Write the reserved DOI into `paper.md` line 92 and into `CITATION.cff` (`identifiers:` block with `type: doi` and `value: 10.5281/zenodo.XXXXXXX`) today; publish the matching tagged release (`v1.0.0` already in `CITATION.cff`) at acceptance. Suggested replacement sentence:
> "An archived release of this code and the wagepan extract is deposited at [https://doi.org/10.5281/zenodo.XXXXXXX](https://doi.org/10.5281/zenodo.XXXXXXX) (v1.0.0); current development is at [https://kevinschoenholzer.com/sim-paneldata/](https://kevinschoenholzer.com/sim-paneldata/)."

### `bootstrapFocal` and `permuteFocal` use unseeded `Math.random()` — inference is not reproducible  — area: code-correctness, **CRITICAL**
**File / location**: `engine.js` lines 430 and 451 (and the consumers `nullReplication` 455–457, `specCurveInference` 459–469, `bootstrapFocal` 422–438)
**Problem**: Both routines call `Math.random()` directly. Bootstrap CIs and the spec-curve permutation p-values therefore change on every reload, and the exported R/Python/Stata code carries no seed.
**Why it matters**: The paper's central rhetorical claim is cross-language, bit-for-bit reproducibility (paper.md line 42, abstract). Any randomized output the engine produces silently violates that claim.
**Recommended fix**: Introduce a stateful module-level `mulberry32(seed)` closure that both routines pull from. Plumb a `seed` field through (UI input, default e.g. `12345`), `Spec.seed`, `engine.runSpec`, `nullReplication`, `bootstrapFocal`. Reseed on each "Run" click. Encode the seed in the permalink (`stateToHash` p.s = S.seed) and into the exported code blocks in `exports.js` so the seed travels with the snippet. Add a `hashStringToUint32(seed)` helper so users can enter human-readable seed strings.

### `permuteFocal` breaks within-person structure — wrong sharp null for panel data  — area: code-science, **CRITICAL** (rated "major" in the JSON; promoted because it interacts with the previous item and the paper presents it as a methods feature)
**File / location**: `engine.js` lines 449–453; nullReplication 455–457; called from `specCurveInference`
**Problem**: A flat Fisher–Yates over all person-years destroys the persistence of focal within persons and its correlation with all other covariates. For `educ` (essentially time-invariant within persons) the null is qualitatively wrong; for `union`/`married` it manufactures pseudo-switchers and is inconsistent with the cluster (person) bootstrap that the same file already implements.
**Why it matters**: The paper sells "randomization inference on the whole curve" as a methods feature. A referee aware of MacKinnon & Webb (2018), Bertrand et al. (2004), or the cluster-permutation literature will flag this on first read.
**Recommended fix**: Two clearly labelled nulls: (a) person-level permutation of the person-mean focal (between-person null, matches BE/RE/POLS); (b) Fisher sign-flip of within-person deviations (within-person null, matches FE/FD). Expose the null choice in the UI next to the inference button. For `educ`, restrict to option (a). Add a `nullSpec` field in the export bundle so the reported p-value carries its null assumption.

## 3. Major issues (should-fix)

### "Seven estimators / 480 specs / 1,440 cells" arithmetic does not add up  — paper-app-consistency
**File / location**: `paper.md` lines 6, 16, 34, 42, 66; `verification/js_multiverse.csv`; `engine.js` line 373; `figures/multiverse_union.csv`
**Problem**: `paper.md` line 34 enumerates seven estimators and asserts "480 identified specifications per estimand"; line 42 asserts "all 1,440 specifications". 5 × 16 × 2 × 3 = 480 and 3 × 480 = 1,440, but 7 × 16 × 2 × 3 = 672 / 2,016. `js_multiverse.csv` contains zero `cre`/`twfe` rows. Figure 2 plots 160 = 5 × 16 × 2 specs (full sample only).
**Why it matters**: This is the first arithmetic any methods referee will redo, and the inconsistency directly undermines the headline reproducibility claim.
**Recommended fix**: Prefer the cheaper option — restate as "five estimators in the default multiverse, with two further estimators (CRE/Mundlak and two-way FE) available for interactive exploration and shown to reproduce the within coefficient (CRE) and the within+year-FE estimate (TWFE)." Justify the exclusion of CRE from the variance decomposition: its focal coefficient is mechanically identical to FE, so including it would double-count within signal. Update line 34 ("up to 480 identified specifications per estimand from the five-estimator default grid"); update Figure 2 caption to "The specification curve on the full sample: 160 estimates (5 estimators × 16 covariate combinations × 2 year-FE choices); the sample axis is fixed here for readability and contributes 0.5% to estimate variance per Table 1." Mirror in `paper-jee.md` and `paper-jose.md`.

### No top-level `README.md` for `sim-paneldata/`  — paper-reproducibility
**File / location**: missing `sim-paneldata/README.md`
**Problem**: Only `verification/README.md` exists. GitHub's landing page and Zenodo's archive listing surface the top-level README, and reproducibility reviewers expect a project overview, dependency list, exact run order, software versions, data origin, license pointer, and citation pointer there.
**Recommended fix**: Create `sim-paneldata/README.md` with: project description, "What is this?" paragraph (link to `paper.md` as the canonical submission, note `paper-jose.md` and `paper-jee.md` as venue variants), dependency list (R ≥ 4.3, `wooldridge`, `plm`, `jsonlite`, `lmtest`, `sandwich`, `readr`, `ggplot2`, `patchwork`; Python ≥ 3.11, `linearmodels`, `pandas`, `numpy`; Node ≥ 18), run order (`Rscript data-prep.R → node tools/export_multiverse.cjs → node verify.cjs → Rscript verification/replicate_plm.R → python3 verification/replicate_python.py → Rscript analysis/variance_decomp.R → Rscript figures/make_figure2.R → Rscript figures/make_figure.R`), license pointer, citation pointer (CITATION.cff + Zenodo DOI), and link to the live build.

### Headless regeneration of `js_multiverse.csv` is undocumented  — paper-reproducibility
**File / location**: `verification/README.md` line 12; `exports.js#multiverseCSV` line 156; no `tools/` directory
**Problem**: README says "regenerate from `../engine.js`" with no command. The only producer is a browser button. Downstream artefacts (`replicate_plm.R`, `replicate_python.py`, `variance_decomp.R`, `make_figure2.R`, `make_figure.R`) all read these CSVs.
**Recommended fix**: Create `tools/export_multiverse.cjs` (Node) that requires `./engine.js`, reads `data/wagepan.json`, calls `E.enumerateMultiverse` over the published axis grid (5 estimators × 16 covariate combos × 2 yearFE × 3 samples for `verification/js_multiverse.csv`; full-sample slice for `figures/multiverse_union.csv`), and writes both files with the existing column schemas. Refactor `Exports.multiverseCSV` (`exports.js`) into a shared formatter so Node and browser emit byte-identical rows. Document the command in the new top-level README.

### Paper does not name the scripts that produce Table 1 and the figures  — paper-reproducibility
**File / location**: `paper.md` line 92
**Problem**: "every figure and table can be regenerated" — but no script names, no run order. Table 1 is from `analysis/variance_decomp.R`; Figure 1 (estimator spectrum) from `figures/make_figure2.R`; Figure 2 (spec curve) from `figures/make_figure.R`. The `make_figure2.R`-produces-Figure-1 mismatch is itself a reader trap, reinforced by an inline comment in `make_figure2.R` calling itself "Figure 2."
**Recommended fix**: Replace the sentence in Data accessibility with: "Table 1 is produced by `analysis/variance_decomp.R`; Figure 1 by `figures/make_figure2.R`; Figure 2 by `figures/make_figure.R`. The CSVs they consume are generated by `node tools/export_multiverse.cjs`. The full pipeline is `Rscript data-prep.R; node tools/export_multiverse.cjs; node verify.cjs; Rscript verification/replicate_plm.R; python3 verification/replicate_python.py; Rscript analysis/variance_decomp.R; Rscript figures/make_figure2.R; Rscript figures/make_figure.R`." Rename `make_figure.R → make_spec_curve_union.R` and `make_figure2.R → make_estimator_spectrum.R`, and fix the inline "Figure 2" comment in `make_figure2.R`.

### `data-prep.R` writes to a hardcoded absolute path  — code-correctness
**File / location**: `data-prep.R` lines 23–25
**Problem**: `write_json(out, "/Users/kevinschoenholzer/.../sim-paneldata/data/wagepan.json", ...)`. The script cannot run on any other machine. The four sibling R scripts already use the `setwd(dirname(--file=…))` idiom.
**Recommended fix**: Apply the same idiom at the top of `data-prep.R`:
```r
script_path <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
if (!nzchar(script_path)) {
  script_path <- sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
}
setwd(dirname(script_path))
write_json(out, "data/wagepan.json", dataframe = "rows", digits = 4, auto_unbox = TRUE)
```
Also write `data/wagepan.meta.json` (sidecar) with `packageVersion("wooldridge")`, `Sys.Date()`, `R.version.string`, and `digest::digest("data/wagepan.json", algo = "sha256")`.

### Software versions not pinned  — paper-reproducibility
**File / location**: `CITATION.cff`; `verification/README.md`; `paper.md`
**Problem**: Only `linearmodels 7.0` is captured in `replication_report_python.md`; nothing pins `plm`, `wooldridge`, R, Python, Node, or the OS. `plm`'s Swamy–Arora variant has shifted across minor versions, which is exactly what produces the ~3e-4 RE residual the paper documents.
**Recommended fix**: (i) Append `sessionInfo()` to `verification/replicate_plm.R` (write to `verification/sessionInfo.txt`); (ii) append `pip freeze | grep -E "linearmodels|pandas|numpy|statsmodels"` to `verification/replicate_python.py` (write to `verification/requirements-frozen.txt`); (iii) commit a `renv.lock` and a `requirements.txt` at the project root; (iv) embed exact versions in `paper.md` §Implementation; (v) mirror in `CITATION.cff` under a `references:` block.

### `CITATION.cff` lacks `preferred-citation:` for the paper  — paper-reproducibility
**File / location**: `CITATION.cff` (whole file)
**Problem**: Only `type: software` is declared. "Cite this repository" returns v1.0.0, not the methods note.
**Recommended fix**: Append a `preferred-citation:` block of `type: article` with TBD `journal`/`doi` (commented "update on acceptance"). Keep author/title/year as placeholders until the venue is chosen. This is the standard CFF pattern for software–paper pairs.

### `verify.cjs` covers only 5/7 estimators with classical SEs and 1 spec each  — verification-coverage
**File / location**: `verify.cjs` lines 22–26 (CANON) and 36–52 (loop)
**Problem**: `twfe`, `cre`, `robust` (HC1), `cluster` (CR1), `hausmanFull`, controls beyond experience, `yearFE = true`, samples `nohealth`/`trim` — all unverified. The Hausman line prints the 1-df focal stat next to a 3-df plm reference and does not call `check()`.
**Recommended fix**: Expand CANON to include twfe (plm `effect="twoways"`) and cre (asserted equal to within within machine precision; no new plm reference needed). Add HC1 canonical values from `sandwich::vcovHC(model, type="HC1")` and CR1 values from `plm::vcovHC(model, type="HC1", cluster="group")`. Add at least one row per estimator with `yearFE=true` and one with `sample="nohealth"`. Replace the bare `console.log` Hausman line with `check("hausmanFull union+exper+expersq", HF.stat, 139.94, 0.5)` and a single-coef analogue derived from a manually computed focal-only Hausman in `data-prep.R` (printed and pinned).

### CRE/TWFE Mundlak equivalence and `hausmanFull` have no automated test  — verification-coverage
**File / location**: `verify.cjs`; engine.js `cre`/`twfe` branches; `hausmanFull` 402–418
**Problem**: The paper's claim "the correlated-random-effects coefficient reproduces the fixed-effects estimate to machine precision, as the Mundlak result requires" rests on a manual check; `hausmanFull` is exported but never invoked from `verify.cjs`. Manual run shows CRE − FE = 4×10⁻¹⁵ and `hausmanFull` = 139.94003 — both match — but a regression would silently land.
**Recommended fix**: Add
```js
check("CRE = FE on union (Mundlak)", E.runSpec(DATA, ssSpecOf("cre"), META).coefs.find(c=>c.name==="union").b,
       E.runSpec(DATA, ssSpecOf("within"), META).coefs.find(c=>c.name==="union").b, 1e-10);
check("hausmanFull union+exper+expersq", E.hausmanFull(fe, re).stat, 139.94, 0.5);
```
Lock the df at 3 to catch silent coef-name-matching regressions.

### `bootstrapFocal`, `curveStats`, `permuteFocal`, `nullReplication`, `specCurveInference` have zero automated tests  — verification-coverage
**File / location**: `engine.js` 422–469; nothing in `verify.cjs`
**Recommended fix**: Add a smoke test that runs `E.bootstrapFocal(DATA, spec, META, 200)` and asserts `ci.lo < ci.hi`, `0.5 * res.se < ci.se < 2 * res.se`, and that the SE is within a factor of 2 of the cluster-robust analytical SE. For `permuteFocal`, assert that the multiset of focal values is preserved (sorted element-wise compare). For `specCurveInference` under H0, assert that the empirical median of permuted estimates is within 0.02 of 0. Put expensive coverage simulations in `verify-slow.cjs` invoked weekly via `workflow_dispatch`.

### `analysis/variance_decomp.R` has no expected values pinned in CI  — verification-coverage
**File / location**: `analysis/variance_decomp.R` lines 26–31
**Problem**: The script prints Table 1 numbers; nothing asserts them. A change to `engine.js` or `js_multiverse.csv` would silently change the paper's headline figures.
**Recommended fix**: Write a one-time `analysis/variance_decomp_expected.csv` capturing the current shares per focal × dimension, then add `stopifnot(all(abs(observed - expected) < 1e-6))` at the end of the script. Hook into CI.

### No CI runs `verify.cjs`/`replicate_plm.R`/`replicate_python.py` on commits  — verification-coverage
**File / location**: only `.github/workflows/update-us-macro.yml` exists
**Recommended fix**: Add `.github/workflows/verify-sim-paneldata.yml`: trigger on push and pull_request when `sim-paneldata/**` changes; run `node sim-paneldata/verify.cjs` (fast, must pass per push); add a separate `workflow_dispatch` + weekly schedule job that runs the R and Python replications inside `r-lib/actions/setup-r@v2` and `actions/setup-python@v5` with `renv` / `requirements.txt` pinning, and fails when max |Δβ| exceeds 1e-5 (R) or 1e-6 (Python).

### Variance-decomposition table reports neither residuals nor R²  — paper-format
**File / location**: `paper.md` lines 56–62 (Table 1)
**Problem**: Rows sum to 91.3% (union), 72.0% (marriage), 95.7% (schooling). The 28% marriage residual is large and substantively informative (estimator × covariate interaction).
**Recommended fix**: Add a "Residual" column so rows sum to 100%. Footnote: "Residuals are estimator-by-covariate (and higher-order) interactions left out of the main-effects ANOVA; the model R²'s are 0.913, 0.720 and 0.957, computed from the same `lm()` used for the shares." Add one sentence in Results: "The 28% residual for the marriage row signals substantial estimator × covariate interaction (the experience control and year effects shift marriage under pooled/between estimators but barely move it under FE/FD); a fully interacted decomposition would attribute that share. We retain the main-effects decomposition for comparability with conventional ANOVA-style summaries."

### Abstract overclaims scope  — paper-positioning
**File / location**: `paper.md` line 6
**Problem**: "Multiverse analyses of longitudinal data should therefore treat the estimator as a first-class forking path." — no scope qualifier despite the body's two acknowledgements of the single-canonical-dataset framing.
**Recommended fix**: Replace the final sentence with:
> "In this demonstration the estimator dimension dominates, and we argue multiverse analyses of longitudinal data should treat the estimator as a first-class forking path and report which dimension drives the spread."
Also trim to ~190 words to leave a margin under RSOS's 200-word cap (drop "quantitative" from sentence 1; "full" from the enumeration sentence; tighten "widely used wage panel" → "wage panel").

### Modal focus trap and screen-reader live regions missing  — ux-accessibility
**File / location**: `app.js` openFig/closeFig 57–65; `index.html` 552–558 (`#figModal`); 361 (`#eqContent`)
**Problem**: `#figModal` does not move focus to `#figClose`, does not trap focus, does not restore focus on close; `#eqContent` rewrites `innerHTML` with no `aria-live`; the spec-curve SVG has no `role="img"`/`aria-label` and no text/table alternative; chip tabs are `<span>`, not keyboard-focusable, with no `role="tab"`.
**Recommended fix** (focused diff sketch):
```js
// app.js openFig
var savedEl = document.activeElement;
$("figModal").classList.add("open");
$("figClose").focus({preventScroll:true});
$("figModal").setAttribute("aria-hidden","false");
function trap(e){ if(e.key==="Tab"){ e.preventDefault(); $("figClose").focus(); } }
$("figModal").addEventListener("keydown", trap);
function closeFig(){ /* …existing… */
  $("figModal").removeEventListener("keydown", trap);
  $("figModal").setAttribute("aria-hidden","true");
  if(savedEl && document.contains(savedEl)) savedEl.focus();
}
```
For `#eqContent`: in `index.html` add `aria-live="polite"` and `aria-atomic="false"`; also add a visually-hidden plain-language summary node (`role="status"`) that announces "Fixed Effects, coefficient on union 0.083, classical SE" so screen readers do not re-read the glyph-heavy math on every change. For the spec-curve SVG, call `ariaChart('chart-speccurve', summary)` at the end of `drawSpecCurve` where `summary` is built from `S.multiSummary`. Convert the five chip tabs to `<button type="button" role="tab" aria-selected="…" aria-controls="…">` and wrap the `<div class="tabs">` in `role="tablist"`.

### Cancel control for long-running bootstrap/inference  — ux-accessibility
**File / location**: `app.js` 749–761 (`btn-bootstrap`), 762–785 (`btn-inference`)
**Problem**: The inference path runs 100 replications × up to 2688 specs each; the bootstrap runs 500 cluster resamples synchronously inside one `setTimeout(fn,30)`. No cancel, no progress estimate, no ceiling. Comb count > 1500 only logs a warning.
**Recommended fix**: Introduce a module-level `runToken` (integer); each click increments it and the chunked `step()` aborts if the local copy no longer matches the module value. Add a "Cancel" button (visible while running) that increments the token. For the bootstrap, refactor `bootstrapFocal` to accept `{start, end, B, accumulator}` and chunk via `setTimeout(step, 0)` like the inference button does. Add a confirm dialog when `comboCount × B > 50000` and hard-cap at 100000 spec-equivalents.

### Permalink does not encode multiverse axes  — ux-accessibility / paper-app-consistency
**File / location**: `app.js` 601–618
**Problem**: `stateToHash` encodes only single-spec state. A colleague pasting a Multiverse-tab permalink lands on the default axes, not the analyst's. This directly contradicts the lab's pedagogy that "each combination is a complete, defensible analysis."
**Recommended fix**:
```js
function stateToHash(){
  var on = Object.keys(S.ss.controls).filter(k => S.ss.controls[k]);
  var p = { v:1, focal:S.focal, outcome:S.outcome, se:S.seType, tab:S.viz,
            est:S.ss.estimator, ctrl:on.join(","),
            yfe:S.ss.yearFE?1:0, samp:S.ss.sample, seed:S.seed,
            axEst:S.axes.estimators.join(","),
            axCtrl:Object.keys(S.axes.controlsVary).filter(k=>S.axes.controlsVary[k]).join(","),
            axYfe:S.axes.yearFEVary?1:0,
            axSamp:S.axes.sampleVary.join(",") };
  /* …existing serialization… */
}
function applyHash(q){
  /* …existing single-spec branches with own-property allow-list:
     if (["union","married","educ"].indexOf(q.focal) >= 0) S.focal = q.focal; …  */
  if (q.axEst){ S.axes.estimators = q.axEst.split(",").filter(e => E.ESTIMATORS_ALL.indexOf(e) >= 0); }
  if (q.axCtrl != null){
    var arr = q.axCtrl ? q.axCtrl.split(",") : [];
    Object.keys(S.axes.controlsVary).forEach(k => { S.axes.controlsVary[k] = arr.indexOf(k) >= 0; });
  }
  if (q.axYfe != null) S.axes.yearFEVary = q.axYfe === "1";
  if (q.axSamp){ S.axes.sampleVary = q.axSamp.split(",").filter(s => ["full","nohealth","trim"].indexOf(s) >= 0); }
  buildAxisControls(); // refresh checkboxes
  if (q.seed) S.seed = q.seed;
}
```
Also tighten `q.focal` validation to an explicit allow-list to close the `__proto__` truthy-lookup gap (currently a UI-corruption hazard rendering "[object Object]" in pills).

## 4. Minor issues / nits

### Estimator / inference / variance-decomp scholarship
- **"Marginal contribution" wording is loose** (paper.md line 38). Replace with: "Because the estimator is entered first in a sequential (Type-I) ANOVA, its reported share is the proportion of the total sum of squares attributable to estimator before adjusting for the other dimensions; under the (near-)orthogonality of the balanced multiverse design used here this coincides with the R² of a regression of the point estimate on estimator alone (Type-III and Type-I shares therefore agree to the displayed precision)." Also report `car::Anova(m, type=3)` and `relaimpo::calc.relimp(type="lmg")` (Shapley) as robustness checks.
- **Pooled OLS framing** (line 14): "Pooled regression uses all of the variation in the data" — tighten to "Pooled regression uses both within- and between-unit variation, weighting them by their share of total variance in the regressor".
- **Pooled-vs-within dominance is partly mechanical** for any panel regressor with non-trivial between-unit selection. Concede this in the Discussion and quantify within/between variance shares for the three focal regressors (mentioned but not numbered at line 70).
- **SE-axis caveat** (line 34): state that the Figure 1 green shading uses cluster-robust SEs (match Figure 2's caption). Optionally tabulate the share of significant specs by SE family in an appendix.
- **Swamy–Arora variant** (line 42): name it in Methods, not only in Verification.
- **DiD/IV/event-study generalisation** (line 84): replace "generalizes naturally" with the hedged formulation citing Goodman-Bacon (2021), de Chaisemartin & D'Haultfœuille (2020), Callaway & Sant'Anna (2021), Sun & Abraham (2021).
- **First-differences drift convention** is not stated. Add one sentence to Methods: "First differences are estimated with a drift intercept, matching the convention used by `plm` (≥ 2.0) and `linearmodels`; this is equivalent to a common time trend in levels and is what verifies to <1×10⁻⁵ against both reference implementations."
- **Two-way FE implementation** (engine.js 361–364): document that `twfe` is computed as `feWithin` with year dummies forced in, that this is identical to `plm(effect="twoways")` for the balanced wagepan, and add a `twfe` row to CANON.
- **Random-effects clipping** (engine.js 244): when `s2b - sigma_e2/Tbar ≤ 0`, the model silently collapses to pooled. Track `model.varCompClipped = true`, surface a footer note in the regression table, and `console.warn` mirroring `plm`'s behaviour.
- **Between SE override** (engine.js 213): `between()` hardcodes `seType:"classical"`, ignoring user choice. Either thread `spec.seType` through (with HC1 honestly computed for `"robust"`) or add a footer note in the regression table.
- **`ols` df coercion** (engine.js 111): `Math.max(1, df)` masks under-identified within models. Guard at top of `ols` with `if (opts.df != null && opts.df <= 0) return fail("model under-identified");` and propagate `msg` into the UI.
- **Cluster SE with G < 2** (engine.js 121): currently fills `vcov` with NaN and returns `ok:true`. Fall back to HC1 (with an explicit `model.seType = "robust"` stamp) or return `ok:false`. Switch `identified` to `isFinite(c.b) && isFinite(c.se)` in `runSpec` (line 367) so all downstream summaries treat the spec consistently.
- **Cluster-robust RE caveat**: paper-app does not disclose that CR1 on the RE GLS-transformed data is the Arellano-style sandwich, whose validity hinges on Swamy–Arora being correct. Add one Methods sentence.

### Scholarship and citations
- **specr** (Masur & Scharkow, 2023, AMPPS 6(2)) absent from Related Work and bibliography. Add at paper.md line 20.
- **Young (2019) QJE on randomization inference** absent. Add at line 20 and at line 34 where CR1 is introduced.
- **Patel, Burford & Ioannidis (2015)** vibration-of-effects analogue from epidemiology absent. Add at line 20.
- **Cameron & Miller (2015) Practitioner's Guide to Cluster-Robust Inference** uncited next to the CR1 introduction (line 34).
- **Card (1996)** cited for the qualitative point only; anchor quantitatively to the NLSY-young-men column (OLS 0.15–0.20; within 0.07–0.10) so the curve's range can be compared.
- **Union-premium literature**: Lewis (1986) and Hirsch (2004) at line 24.
- **Marriage-premium puzzle literature**: cite Korenman & Neumark (1991), Cornwell & Rupert (1997), Krashinsky (2004), Killewald & Lundberg (2017), Ludwig & Brüderl (2018) at line 24 near the Card sentence.
- **Linearmodels Python package** uncited. Add `@Manual{linearmodels, author={Sheppard, Kevin}, title={linearmodels}, year={2024}}` and cite at line 42.
- **wooldridge BibTeX key** is `wooldridge2010pkg` while `year = 2024`. Rename to `shea2024wooldridge` to remove the visual mismatch.
- **Sala-i-Martin 1997** AER P&P entry is correct; optionally add the "(Papers and Proceedings)" tag.

### Paper format / RSOS hygiene
- **Title** is RSOS-compliant in length and case. Consider hedging "Estimator choice is the dominant forking path" → "Estimator choice as a dominant forking path". Optional.
- **ORCID** appears as `ORCID 0000-0001-9892-5869` in the author footnote. Production prefers `https://orcid.org/0000-0001-9892-5869`.
- **Affiliation/ORCID/email packed into a footnote**: split into structured YAML keys for production.
- **Citations are author-date**; RSOS final style is numbered (Vancouver-like). Leave `csl:` commented until the venue is fixed; do not commit RSOS's CSL prematurely. Rewrite the `[-@vellaverbeek1998]` suppress-author construct on line 30 ("the wage panel assembled by Vella and Verbeek") to read sensibly under a numbered style.
- **Lower-case "table 1"** on line 92 → capitalise.
- **Acknowledgements** currently only recapitulates data provenance. Either delete the section or add one sentence acknowledging no specific reviewer/discussant input.
- **AI-use declaration** (lines 95–96) names a generic category. Production may push back; replace with: "The author used Claude Code (Claude Opus 4.7, Anthropic) to assist with implementation of the browser-based software and prose editing. AI tools were not used for study design, statistical analyses, or scientific claims. All code, analyses, numerical results, and text were reviewed and verified by the author, who takes full responsibility for the content."
- **Verification tolerances** (line 42): "1×10⁻⁶ in R" is met only for pooled (max |Δβ| 9.25e-7); between/within/fd are at 1.01e-6, 2.27e-6, 2.04e-6. Soften to "~1×10⁻⁶ (max |Δβ| ≈ 2.3×10⁻⁶)". Python at 5×10⁻⁷ is met only at the displayed precision (max 4.99e-7); restate as "~5×10⁻⁷ (max observed 4.99×10⁻⁷)". Also specify that the 3e-4 RE figure applies in both R and Python.
- **Multidisciplinary glossing**: define "union (wage) premium", "return to schooling", "identifying assumption" on first use; either drop the "joint Hausman statistic" parenthetical from line 42 or gloss it as "the standard specification test comparing fixed- and random-effects estimates". Do not re-gloss "between estimator", "Mundlak device", "Swamy–Arora".
- **Union spread** in Results (line 64): change "0.04 to 0.26" to "about 0.02 to 0.27" (variance_decomp.R: range [0.017, 0.270]); add one sentence noting that FE/FD precision is comparable to RE in this panel and tighter than pooled/between, so the curve width is identification-driven, not noise-driven.
- **Browser-tool clause** in the abstract should be cut (the reproducibility claim is already made by R+Python verification); compress the "A teaching application" subsection (lines 78–80) into a single line in Limitations or Data Accessibility.
- **Between regression with N ≈ k**: warn and refuse to return an RE fit when `btw.df < max(5, k_between)` rather than silently collapsing to pooled.

### Code small stuff
- **`runSpec` mutates `spec.outcomeKey`** (engine.js 337). Clone `spec` at the top, or compute `outcomeKey` as a local variable.
- **Python FD export omits drift** (`exports.js` 134). Replace `FirstDifferenceOLS(y, X).fit(...)` with hand-rolled `panel.groupby("entity").diff().dropna()` + `PooledOLS(yD, sm.add_constant(XD)).fit(...)` so the exported snippet matches R, Stata and the engine.
- **Stata CRE generator** (`exports.js` 95) leaves double spaces and silently omits Mundlak means for `i.industry`/`i.region`. Either generate the means or add an honest comment.
- **R trim vs Stata `_pctile`**: document the quantile-type-7 vs Stata-default mismatch; suggest computing cutoffs manually in Stata to match.
- **Bootstrap `Math.random()` reproducibility** is the same finding as the critical seeding item — re-list under critical.
- **Gram–Schmidt collinearity tolerance** (engine.js 172) uses each column's own raw norm; switch to running `pivotMax`-relative.
- **Download `alert(...)`** in `app.js` 716–722 should use an in-page banner with `role="status" aria-live="polite"`.

## 5. Recommendations to strengthen scholarship and scientific rigor

### Type-II / Type-III SS and Shapley as a robustness exhibit
Extend `analysis/variance_decomp.R` to compute three decompositions:
```r
library(car); library(relaimpo)
m   <- lm(b ~ estimator + experience + industry + region + health + yearFE + sample, data = sub)
typeI   <- 100 * anova(m)[, "Sum Sq"]      / sum(anova(m)[, "Sum Sq"])
typeIII <- 100 * Anova(m, type = 3)[, "Sum Sq"] / sum(Anova(m, type = 3)[, "Sum Sq"])
shapley <- 100 * calc.relimp(m, type = "lmg", rela = TRUE)$lmg
```
Print all three side-by-side. Report `omega2_p` (partial omega-squared) for df adjustment. In the paper, present Type-I as the headline and Type-III/Shapley as a footnote: "Across Type-I, Type-III and Shapley (LMG) decompositions the estimator dimension's share for the union premium is 84.6%, 84.4% and 82.1% respectively; the qualitative dominance of the estimator dimension is robust to the chosen decomposition convention."

### Panel-aware permutation inference
Replace `permuteFocal`'s flat shuffle with two options exposed in the UI:
1. **Cluster permutation** (between-person null): permute the person-mean focal across persons, leaving within-person dynamics intact. Default for `educ`.
2. **Within-person sign-flip** (wild-cluster Rademacher): for each person, multiply the within-person residual of the focal by ±1 with prob 0.5. Default for FE/FD-heavy curves.
Surface the null choice next to the inference button; the exported summary should carry `nullKind`. Cite MacKinnon & Webb (2018), Bertrand, Duflo & Mullainathan (2004), Young (2019).

### Controls-only-vs-full multiverse contrast
Quantify the headline narrative claim that a conventional robustness analysis "explores less than a tenth of the analytic variation" (line 64). Add a sub-table:
| Estimand | Conventional (controls × sample × yearFE only, estimator fixed) | Full multiverse |
| --- | --- | --- |
| Union | median IQR across estimators = X | IQR = [0.04, 0.21] |
| Marriage | … | … |
| Schooling | … | … |
Annotate Figure 1 with per-estimator IQR bands and a horizontal rule at the full-multiverse IQR. This removes the rhetorical claim and replaces it with a number a referee can verify.

### Mundlak-equivalence regression test in CI
`engine.js`'s `cre` and `within` focal coefficients agree to 4×10⁻¹⁵; this is a theorem of the design, not a numerical coincidence. Lock it in `verify.cjs`:
```js
var bFE  = E.runSpec(DATA, spec("within"), META).coefs.find(c => c.name === "union").b;
var bCRE = E.runSpec(DATA, spec("cre"),    META).coefs.find(c => c.name === "union").b;
check("Mundlak CRE = FE on union", bCRE, bFE, 1e-10);
```

### Seeded PRNG, exported
A single `mulberry32(hashStringToUint32(seed))` closure used by `bootstrapFocal`, `permuteFocal`, and any other randomized routine. Plumb a `seed` field through the UI, the permalink, the engine, and `exports.js` so the generated R/Python/Stata snippets carry `set.seed(12345)` / `np.random.seed(12345)` / `set seed 12345`. This also fixes the verification-coverage gap for bootstrap and inference.

### Additional citations (already enumerated above)
specr (Masur & Scharkow, 2023); Young (2019); Patel, Burford & Ioannidis (2015); Cameron & Miller (2015); MacKinnon & Webb (2018); Lewis (1986); Hirsch (2004); Imbens & Wooldridge (2009); Arkhangelsky & Imbens (2024); Goodman-Bacon (2021); de Chaisemartin & D'Haultfœuille (2020); Callaway & Sant'Anna (2021); Sun & Abraham (2021); for marriage-premium: Korenman & Neumark (1991), Killewald & Lundberg (2017), Ludwig & Brüderl (2018).

### Residual disclosure in Table 1
Add a Residual column so rows sum to 100% and disclose the model R² for each estimand. Explain the 28% marriage residual as estimator × covariate interaction left in the residual by the main-effects decomposition.

### Variance components transparency
Expose `sigma_e2`, `sigma_u2`, `theta` from the RE fit in the equation banner. Pin them in `verify.cjs` against `plm::ercomp()` so the ~3×10⁻⁴ RE residual is traceable to a specific variance-component disagreement rather than left as "alternative implementations of the Swamy–Arora estimator".

### Scope-explicit abstract
Replace the final sentence as in §2; trim to ~190 words to leave margin. The combination of (i) a scope qualifier and (ii) a more concrete prescriptive verb ("report which dimension drives the spread") removes the discipline-wide overreach without softening the contribution.

## 6. Concrete paper edits (specific text suggestions)

The following are minimal `paper.md` patches, each tagged with the line in the current file. Where multiple variant manuscripts exist (`paper-jose.md`, `paper-jee.md`), mirror.

**Line 6 (abstract last sentence)** — replace
> "Multiverse analyses of longitudinal data should therefore treat the estimator as a first-class forking path."

with

> "In this canonical demonstration the estimator dimension dominates, and we argue multiverse analyses of longitudinal data should treat the estimator as a first-class forking path and report which dimension drives the spread."

**Line 14 (Introduction)** — replace
> "Pooled regression uses all of the variation in the data; a between estimator uses only differences across units in their long-run averages; a within (fixed-effects) estimator discards all cross-unit variation"

with

> "Pooled regression uses both within- and between-unit variation, weighting them by their share of total variance in the regressor; a between estimator uses only the long-run cross-unit means; a within (fixed-effects) estimator uses only within-unit changes"

**Line 20 (Related Work)** — replace
> "In psychology and the broader sciences it appears as multiverse analysis [@steegen2016] and specification-curve analysis [@simonsohn2020]; in economics it has a longer pedigree in extreme-bounds analysis [@leamer1985; @levinerenelt1992; @salaimartin1997]."

with

> "In psychology and the broader sciences it appears as multiverse analysis [@steegen2016] and specification-curve analysis [@simonsohn2020; @masur2023specr], with the `specr` R package emerging as a de facto standard that currently does not vary the estimator dimension for panel data. In epidemiology the same idea appears as vibration-of-effects analysis [@patel2015vibration]. In economics it has a longer pedigree in extreme-bounds analysis [@leamer1985; @levinerenelt1992; @salaimartin1997], and inference-side fragility has been highlighted by [@young2019]. Beyond the dimension varied, the characterization also differs: extreme-bounds and BACE-style analyses summarise the resulting estimates by whether a worst-case (or posterior-averaged) coefficient retains sign or significance across control sets, whereas our ANOVA decomposition treats each defensible specification as an observation and apportions the variance of the estimate across analytic dimensions."

**Line 24 (after the Card sentence)** — append
> "The substantive interpretation of the male marriage premium has been debated for decades [@korenman1991; @cornwell1997; @krashinsky2004; @killewald2017; @ludwig2018], with the debate turning on whether the cross-sectional premium reflects selection on stable traits or a causal effect — exactly the within/between contrast our specification curve makes visible. For the union premium, the canonical surveys are [@lewis1986; @hirsch2004]; @card1996 and @vellaverbeek1998 are the within/between contrasts most directly relevant here."

**Line 34 (Methods)** — replace "up to 480 identified specifications per estimand" with
> "The default multiverse crosses five estimators (pooled OLS, between, random effects, fixed effects/within, first differences) with the four other axes, yielding 5 × 16 × 2 × 3 = 480 specifications per estimand. Correlated random effects and two-way fixed effects are additionally available in the interactive lab and verified to reproduce the within estimate to machine precision (CRE/Mundlak) and to the canonical two-way result (TWFE); they are excluded from the variance decomposition because the CRE focal coefficient is mechanically identical to the FE focal coefficient. Inference is reported under three SE families (classical, HC1, CR1 cluster-robust by individual [@cameron2015]); first differences are estimated with a drift intercept matching `plm` (≥ 2.0) and `linearmodels`; the random-effects variance components use the Swamy–Arora estimator, which is the within-RE sub-fork held fixed throughout."

**Line 38 (variance decomposition)** — replace the "marginal contribution" sentence with the longer formulation in §4 above, and add as a footnote:
> "Type-III (`car::Anova(m, type = 3)`) and Shapley (LMG, `relaimpo::calc.relimp(m, type = \"lmg\")`) decompositions yield estimator shares of 84.4% and 82.1% for the union premium, 55.8% and 54.9% for the marriage premium; the dominance of the estimator dimension is robust to the decomposition convention."

**Line 42 (Verification)** — replace
> "The pooled, between, fixed-effects and first-difference coefficients agree to better than $1\times10^{-6}$ in R and $5\times10^{-7}$ in Python … and the random-effects coefficients to $3\times10^{-4}$"

with

> "The pooled, between, fixed-effects and first-difference coefficients agree to within ~$1\times10^{-6}$ in R (max $|\Delta\beta| \approx 2.3\times10^{-6}$) and ~$5\times10^{-7}$ in Python (max observed $4.99\times10^{-7}$). The correlated-random-effects coefficient reproduces the fixed-effects estimate to machine precision (max $|\Delta\beta| \approx 4\times10^{-15}$), as the Mundlak result requires. The joint Hausman statistic (the standard specification test comparing fixed- and random-effects estimates) is $\chi^2_3 = 139.94$ in R and $139.940$ in the engine, matching `plm::phtest` to five decimal places. The random-effects coefficients agree to ~$3\times10^{-4}$ in both R and Python (max absolute deviation across the 480 identified specifications per estimand); the residual reflects the documented difference between the engine's Swamy–Arora implementation and the variants in `plm` and `linearmodels`."

**Line 64 (Results)** — replace "ranges from roughly 0.04 to 0.26" with "ranges from about 0.02 to 0.27". Add a sentence at the end of the paragraph:
> "Conditional on estimator, the average IQR of the union premium across the covariate × sample × year-FE grid is 0.018, compared with 0.063 in the full multiverse — a conventional robustness analysis therefore explores roughly 28% of the analytic variation the multiverse exposes."

**Line 66 (Figure 2 caption)** — append: "(cluster-robust by individual; the sample axis is fixed at the full sample, so this is one of three sample slices of the 480-cell union curve decomposed in Table 1)."

**Line 84 (Limitations)** — replace "the argument generalizes naturally to other designs … would be a productive direction for further work" with
> "Whether the same dominance holds for designs whose credibility likewise turns on the identifying choice — difference-in-differences, instrumental variables and event studies — is an open question the present panel demonstration cannot decide; the rapidly developing DiD literature on heterogeneous treatment effects [@goodmanbacon2021; @dechaisemartin2020; @callawaysantanna2021; @sunabraham2021] suggests the estimator axis is at least as consequential there, and a systematic estimator-led multiverse in each design would be a productive direction for further work."

**Line 92 (Data accessibility)** — see §2 for the proposed replacement.

**Lines 95–96 (AI declaration)** — see §4 above.

**Acknowledgements (lines 110–112)** — move the data-provenance sentence to Data accessibility; either remove the Acknowledgements section or replace with a brief honest statement.

## 7. Concrete code/app edits (specific patches, with file paths and function names)

### `engine.js`
- **Add seeded PRNG** at top of file:
```js
function mulberry32(a){return function(){var t=(a+=0x6D2B79F5)|0;t=Math.imul(t^(t>>>15),1|t);t=t+Math.imul(t^(t>>>7),61|t)^t;return((t^(t>>>14))>>>0)/4294967296;};}
function hashStringToUint32(s){var h=2166136261>>>0;for(var i=0;i<s.length;i++){h^=s.charCodeAt(i);h=Math.imul(h,16777619);}return h>>>0;}
var RNG = mulberry32(12345); // default; reseed via setSeed()
function setSeed(s){ RNG = mulberry32(hashStringToUint32(String(s))); }
```
Replace all `Math.random()` calls in `bootstrapFocal` (line 430) and `permuteFocal` (line 451) with `RNG()`. Export `setSeed`.

- **Cluster permutation** in `permuteFocal`. Add `kind` parameter:
```js
function permuteFocal(rows, focal, kind){
  kind = kind || "cluster_mean"; // default
  if (kind === "row") { /* existing flat shuffle */ }
  else if (kind === "cluster_mean") {
    var byId = {}; rows.forEach(r => (byId[r.nr] = byId[r.nr]||[]).push(r));
    var ids = Object.keys(byId);
    var means = ids.map(id => byId[id].reduce((s,r)=>s+r[focal],0)/byId[id].length);
    for (var i = means.length - 1; i > 0; i--){ var j=(RNG()*(i+1))|0,t=means[i];means[i]=means[j];means[j]=t; }
    return [].concat.apply([], ids.map((id,i) => byId[id].map(r => Object.assign({},r,{[focal]:means[i]}))));
  }
  else if (kind === "signflip_within") { /* Rademacher on within-deviations */ }
}
```
Update `nullReplication` and `specCurveInference` to accept and forward `kind`.

- **`ols` under-identified guard** (line 103): add `if (opts.df != null && opts.df <= 0) return fail("under-identified");`. Make `fail()` accept an optional `msg` and attach it.
- **Cluster SE with G < 2** (line 121): replace NaN fill with `return fail("Cluster SE requires G ≥ 2; only " + G + " cluster found");`.
- **`runSpec` outcomeKey mutation** (line 337): replace `spec.outcomeKey = …` with `var outcomeKey = …` and use the local; or `spec = Object.assign({}, spec)` at top of function.
- **`identified`** (line 367): change to `var identified = isFinite(c.b) && isFinite(c.se);`.
- **RE clipping flag** (line 244): add `if (s2b - sigma_e2/Tbar < 0) { model.varCompClipped = true; console.warn("RE: sigma_u^2 clipped to 0; theta = 0, RE reduces to pooled OLS"); }`.
- **`between` SE choice** (line 213): thread `opts.seType` through.
- **Gram–Schmidt tolerance** (line 172): track running `pivotMax` and compare each new column's residual to `1e-7 * Math.max(1e-12, pivotMax)`.
- **`firstDifferences` nGroups** (line 278): `fit.nGroups = new Set(idd).size;` (excludes singletons). Optionally surface `fit.nDropped`.
- **TWFE / CRE comments and docs**: add inline comments explaining the convention and pointing at the corresponding `plm`/`linearmodels` arguments.

### `app.js`
- **Bootstrap race** (line 753): capture `var est = S.ss.estimator; var spec = ssSpec(est);` before `setTimeout`, then use `spec` and `est` throughout the callback.
- **Cancel control** for both `btn-bootstrap` and `btn-inference`: introduce module-level `var runToken = 0;` and check `if (myToken !== runToken) return;` in chunks.
- **Permalink axes** (601–618): see the patch in §3.
- **`q.focal` allow-list** (line 610): `if (q.focal && ["union","married","educ"].indexOf(q.focal) >= 0) S.focal = q.focal;`.
- **Hash sync on axis change**: add `syncHash();` to axis-checkbox `change` handlers (currently only `updateComboCount()` is called).
- **Modal focus management** (openFig/closeFig 57–65): patch in §3.
- **Tab semantics**: convert chip `<span>`s to `<button type="button">`; bind `keydown` for Enter/Space.
- **`alert()` replacement** (716–722): use a shared `<div id="exportStatus" class="alert alert-warning" role="status" aria-live="polite" hidden>` in the export tab.
- **Spec-curve aria summary**: at the end of `drawSpecCurve`, call `ariaChart('chart-speccurve', summary)` where `summary` is built from `S.multiSummary`. Optionally add a "Show as table" toggle that renders `S.multiResults` as a paginated `<table>`.

### `verify.cjs`
- Replace the print-only Hausman line with two `check()` calls (single-coef + full).
- Add CANON entries for CRE (= within) and TWFE (against a plm `effect="twoways"` reference produced by `data-prep.R`).
- Add HC1 and CR1 SE rows for each estimator (canonical references from `sandwich::vcovHC` and `plm::vcovHC`).
- Add at least one `yearFE=true` and one `sample="nohealth"` row per estimator.
- Add smoke tests for `bootstrapFocal`, `permuteFocal` (multiset preservation), `specCurveInference` (median ~ 0 on permuted focal), and `curveStats`.
- Add a Mundlak equivalence check.

### `data-prep.R`
- Replace absolute paths with `setwd(dirname(--file=…))` idiom (with RStudio fallback).
- Write `data/wagepan.meta.json` sidecar with `wooldridge` version, R version, date, sha256.
- Print and (later) pin the focal-only Hausman manually so `verify.cjs` can use it as a 1-df reference.

### `tools/export_multiverse.cjs` (new)
- Require `./engine.js`, load `data/wagepan.json`, iterate the 5-estimator grid for `verification/js_multiverse.csv` and the full-sample slice for `figures/multiverse_union.csv`. Use the existing `Exports.multiverseCSV` formatter (refactor to be Node-friendly).

### `.github/workflows/verify-sim-paneldata.yml` (new)
- On push / PR to `sim-paneldata/**`: install Node 18, run `node sim-paneldata/verify.cjs`.
- On `workflow_dispatch` and weekly cron: set up R via `r-lib/actions/setup-r@v2` with renv lockfile, Python via `actions/setup-python@v5` with `requirements.txt`, run replications, fail if max |Δβ| exceeds 1e-5 (R) or 1e-6 (Python).

## 8. Pre-submission checklist

- [ ] Public GitHub repository exists at the URL cited in `CITATION.cff` line 15 (or `CITATION.cff` is updated to a real URL)
- [ ] Zenodo DOI reserved and inserted into `paper.md` line 92 and `CITATION.cff` `identifiers:` block
- [ ] `sim-paneldata/README.md` created (project overview, dependency list, exact run order, license pointer)
- [ ] `data-prep.R` portable (no absolute paths)
- [ ] `renv.lock` and `requirements.txt` committed; `sessionInfo()` / `pip freeze` written to `verification/`
- [ ] `tools/export_multiverse.cjs` headlessly regenerates `js_multiverse.csv` and `multiverse_union.csv`
- [ ] CI workflow runs `verify.cjs` on every push
- [ ] `verify.cjs` exercises 5 + 2 estimators (cre, twfe), HC1 + CR1, `hausmanFull` (with `check()`), Mundlak equivalence, bootstrap smoke test
- [ ] `engine.js`: seeded PRNG plumbed through `bootstrapFocal`/`permuteFocal`; cluster permutation default; `ols` under-identified guard; RE clipping warning
- [ ] Table 1 has a Residual column, R² reported; one-sentence interpretation of the 28% marriage residual
- [ ] Abstract trimmed to ~190 words and includes single-dataset scope qualifier
- [ ] "Up to 480 identified specifications" / "1,440 specifications" / "160 estimates" arithmetic reconciled in paper.md and Figure 2 caption
- [ ] Verification tolerances softened (1×10⁻⁶ R / 5×10⁻⁷ Python with explicit max-observed)
- [ ] First-differences drift convention stated in Methods
- [ ] Cluster-robust RE caveat (Arellano sandwich, Swamy–Arora dependence) stated in Methods
- [ ] Citations added: specr (Masur & Scharkow 2023), Young (2019), Patel et al. (2015), Cameron & Miller (2015), Lewis (1986), Hirsch (2004), Imbens & Wooldridge (2009), marriage-premium puzzle authors, modern DiD authors, `linearmodels` (Sheppard)
- [ ] `wooldridge` BibTeX key renamed (`shea2024wooldridge`)
- [ ] Type-III (`car::Anova`) and Shapley (`relaimpo::calc.relimp`) shares reported as a footnote to Table 1
- [ ] Controls-only-vs-full multiverse contrast quantified in Results
- [ ] Modal focus trap + return, `aria-live="polite"` on `#eqContent`, ARIA tab roles on chip tabs, fieldset/legend on axis groups, text/table alternative for spec-curve SVG
- [ ] Cancel control for `btn-bootstrap` and `btn-inference`; confirm dialog at >5000 spec-equivalents
- [ ] Permalink encodes `S.axes` and seed; recipient sees the same multiverse
- [ ] AI-use declaration names the specific model (Claude Opus 4.7 via Claude Code)
- [ ] ORCID rendered as `https://orcid.org/0000-0001-9892-5869`
- [ ] Lower-case "table 1" → "Table 1" on line 92
- [ ] `linearmodels` cited at line 42 with BibTeX entry
- [ ] Acknowledgements section either deleted or replaced with non-data-provenance content
- [ ] Funding statement confirmed (Lepori SNSF grant if applicable — Kevin to verify)
- [ ] Live-test of `https://kevinschoenholzer.com/sim-paneldata/` (URL must resolve and the lab must build the multiverse) right before submission
- [ ] Cover letter notes the two venue variants (`paper-jose.md`, `paper-jee.md`) are alternative drafts; only `paper.md` is the submitted version