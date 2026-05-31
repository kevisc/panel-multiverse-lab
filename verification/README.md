# Verification & independent replication

The estimators behind the Panel Data Multiverse Lab are verified at three levels. All scripts
read the same bundled data (`../data/wagepan.json`) the application ships.

| Script | Language | What it checks |
|---|---|---|
| `../verify.cjs` | Node.js | Engine self-test against canonical `plm` values (23 benchmarks: pooled/between/RE/FE/FD coefficients and SEs, plus the Hausman test and the FE/FD non-identification of time-invariant regressors). |
| `replicate_plm.R` | R (`plm`) | Re-fits **every** specification in the default multiverse (3 estimands × 480 = 1,440) and compares the focal coefficient. |
| `replicate_python.py` | Python (`linearmodels`) | Independent re-fit of the same 1,440 specifications in a second language/ecosystem. |

`js_multiverse.csv` is the engine's output for all 1,440 specifications (regenerate from
`../engine.js`); the two replication scripts compare against it and write
`replication_report.md` and `replication_report_python.md`.

## Result

Across all 1,440 specifications, both independent implementations reproduce the JavaScript
engine: the pooled-OLS, between, fixed-effects, and first-difference coefficients agree to
better than 1×10⁻⁶ (R) and 5×10⁻⁷ (Python), and random effects to 3×10⁻⁴ (a Swamy–Arora
variance-component estimator). Identification agreement is 1,440/1,440 with `linearmodels`
and 1,436/1,440 with `plm`; the four exceptions are random-effects specifications whose only
regressor is time-invariant, which `plm` declines to fit (no within degrees of freedom) and
the tool estimates by feasible GLS.

## Reproduce

```sh
node ../verify.cjs              # engine self-test
Rscript replicate_plm.R         # R / plm  (installs wooldridge, plm if needed)
python3 replicate_python.py     # Python / linearmodels
```
