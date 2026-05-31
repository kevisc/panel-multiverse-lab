#!/usr/bin/env python3
"""Independent replication of the Panel Data Multiverse Lab in Python.

Re-fits every specification in the default multiverse (verification/js_multiverse.csv)
with linearmodels (PooledOLS, BetweenOLS, RandomEffects, PanelOLS entity-effects) and a
matching manual first-difference estimator, and compares the focal coefficient to the
browser tool's JavaScript engine. Writes replication_report_python.md.
"""
import json, os, sys
import numpy as np
import pandas as pd
from linearmodels.panel import PooledOLS, BetweenOLS, RandomEffects, PanelOLS

HERE = os.path.dirname(os.path.abspath(__file__))

# ---- load the same bundled data the app uses ----
rows = json.load(open(os.path.join(HERE, "..", "data", "wagepan.json")))
df = pd.DataFrame(rows)
# reference category = most frequent (matches the JS meta)
for col in ("industry", "region"):
    order = df[col].value_counts().index.tolist()
    df[col] = pd.Categorical(df[col], categories=order)
df["entity"] = df["nr"]; df["time"] = df["year"]

js = pd.read_csv(os.path.join(HERE, "js_multiverse.csv"))

def sample(samp):
    if samp == "nohealth": return df[df["poorhlth"] != 1].copy()
    if samp == "trim":
        lo, hi = np.quantile(df["lwage"], [0.01, 0.99])  # numpy default = linear (type 7)
        return df[(df["lwage"] >= lo) & (df["lwage"] <= hi)].copy()
    return df.copy()

def design(dd, focal, experience, industry, region, health, yearFE):
    """Return exog DataFrame (no constant), dropping controls constant in this sample."""
    cols = {focal: dd[focal].astype(float)}
    if experience: cols["exper"] = dd["exper"].astype(float); cols["expersq"] = dd["expersq"].astype(float)
    if health and dd["poorhlth"].nunique() > 1: cols["poorhlth"] = dd["poorhlth"].astype(float)
    X = pd.DataFrame(cols, index=dd.index)
    if industry and dd["industry"].nunique() > 1:
        X = pd.concat([X, pd.get_dummies(dd["industry"], prefix="ind", drop_first=True).astype(float)], axis=1)
    if region and dd["region"].nunique() > 1:
        X = pd.concat([X, pd.get_dummies(dd["region"], prefix="reg", drop_first=True).astype(float)], axis=1)
    if yearFE:
        X = pd.concat([X, pd.get_dummies(dd["time"], prefix="yr", drop_first=True).astype(float)], axis=1)
    return X

def add_const(X):
    return pd.concat([pd.Series(1.0, index=X.index, name="const"), X], axis=1)

def fit_focal(est, dd, focal, X):
    """Return focal coefficient or None if not identified."""
    panel = dd.set_index(["entity", "time"])
    y = panel["lwage"]
    Xp = X.copy(); Xp.index = panel.index
    try:
        if est == "pooled":
            r = PooledOLS(y, add_const(Xp)).fit()
        elif est == "between":
            r = BetweenOLS(y, add_const(Xp)).fit()
        elif est == "random":
            r = RandomEffects(y, add_const(Xp)).fit()
        elif est == "within":
            r = PanelOLS(y, Xp, entity_effects=True, drop_absorbed=True).fit()
        elif est == "fd":
            return fd_focal(dd, focal, X)
        else:
            return None
        return float(r.params[focal]) if focal in r.params.index else None
    except Exception:
        return None

def fd_focal(dd, focal, X):
    """First differences with a drift intercept, matching the engine and plm."""
    d = dd.sort_values(["entity", "time"]).copy()
    cols = list(X.columns)
    sub = pd.concat([d[["entity", "time", "lwage"]].reset_index(drop=True),
                     X.reset_index(drop=True)], axis=1)
    diffs = []
    for _, g in sub.groupby("entity"):
        g = g.sort_values("time")
        dd_ = g[["lwage"] + cols].diff().iloc[1:]
        diffs.append(dd_)
    D = pd.concat(diffs, ignore_index=True)
    yD = D["lwage"].values
    XD = D[cols].values
    # drift intercept + drop (near) constant differenced columns, as the engine does
    XD = np.column_stack([np.ones(len(yD)), XD]); names = ["(drift)"] + cols
    keep = [0] + [j for j in range(1, XD.shape[1]) if np.var(XD[:, j]) > 1e-9]
    XD = XD[:, keep]; names = [names[j] for j in keep]
    if focal not in names: return None
    beta, *_ = np.linalg.lstsq(XD, yD, rcond=None)
    return float(beta[names.index(focal)])

js["py_b"] = np.nan; js["py_ident"] = 0
for i, row in js.iterrows():
    dd = sample(row["sample"])
    X = design(dd, row["focal"], row["experience"], row["industry"], row["region"], row["health"], row["yearFE"])
    b = fit_focal(row["estimator"], dd, row["focal"], X)
    if b is not None and np.isfinite(b):
        js.at[i, "py_b"] = b; js.at[i, "py_ident"] = 1

ident_agree = (js["identified"] == js["py_ident"]).mean()
both = js[(js["identified"] == 1) & (js["py_ident"] == 1)].copy()
both["absdiff"] = (both["b"] - both["py_b"]).abs()

lines = ["# Python (linearmodels) replication of the Panel Data Multiverse Lab", "",
         f"Re-fit of every specification in the default multiverse (3 estimands x 480 = {len(js)}) "
         "using Python's `linearmodels` (and a matching manual first-difference estimator), compared to "
         "the browser tool's JavaScript engine.", "",
         f"- **Identification agreement:** {(js['identified']==js['py_ident']).sum()} / {len(js)} "
         f"({100*ident_agree:.1f}%).",
         f"- **Coefficient agreement** on the {len(both)} specifications identified by both, by estimator:", "",
         "| Estimator | n | max |Δβ| | mean |Δβ| |", "|---|---|---|---|"]
for e in ["pooled", "between", "within", "fd", "random"]:
    s = both.loc[both["estimator"] == e, "absdiff"]
    if len(s): lines.append(f"| {e} | {len(s)} | {s.max():.2e} | {s.mean():.2e} |")
lines += ["", "Generated by `verification/replicate_python.py` with linearmodels "
          f"{__import__('linearmodels').__version__}."]
open(os.path.join(HERE, "replication_report_python.md"), "w").write("\n".join(lines))

print(f"identification agreement: {100*ident_agree:.1f}%")
print(both.groupby("estimator")["absdiff"].agg(["size", "max", "mean"]))
print("wrote verification/replication_report_python.md")
