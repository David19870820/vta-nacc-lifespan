#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
06_preproc_sensitivity.py
=========================
Section 3.8 (recomputation on fully preprocessed data) and Section 3.9
(ROI definition sensitivity) of the manuscript.

The primary analysis uses raw BIDS images, which carry no nuisance
regression. This script repeats the age analysis on coupling recomputed
from the C-PAC preprocessed derivatives, which are motion-corrected and
nuisance-regressed, using the identical ROI definition and filtering
pipeline. It also tests twelve alternative ROI definitions.

Inputs
------
    data_derived/coupling_preprocessed.csv   C-PAC coupling, all ROI variants
    data_derived/full_all_scans.csv          session/age/sex/mean_fd plan
    data_derived/primary_onescan_per_subject.csv   raw-BIDS coupling (agreement test)
    results/fd_by_scan.csv                   per-scan FD, used to fill gaps

Outputs
-------
    results/preproc_analysis.txt
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from _stats import pearsonr, ols                     # noqa: E402

DER = ROOT / "data_derived"
RES = ROOT / "results"
RES.mkdir(exist_ok=True)

L = []


def P(s=""):
    print(s)
    L.append(str(s))


cp = pd.read_csv(DER / "coupling_preprocessed.csv")
plan = pd.read_csv(DER / "full_all_scans.csv")[["subject", "session", "age", "sex_M", "mean_fd"]]
fd = pd.read_csv(RES / "fd_by_scan.csv")

m = plan.merge(cp, on=["subject", "session"], how="inner")

# fill missing mean_fd from the per-scan FD table
fdmap = {}
for _, r in fd.iterrows():
    ses = str(r.get("raw_sessions_1400", "")).split(";") if pd.notna(r.get("raw_sessions_1400")) else []
    for s in ses:
        c = s + "_mean_fd"
        if c in fd.columns and pd.notna(r[c]):
            fdmap[(r["subject"], s)] = float(r[c])
m["fd_scan"] = [fdmap.get((s, ss), np.nan) for s, ss in zip(m["subject"], m["session"])]
m["mean_fd"] = m["mean_fd"].fillna(m["fd_scan"])

P("=" * 76)
P("PREPROCESSED-DATA ANALYSIS (C-PAC, full confound regression)")
P("=" * 76)
P("scans merged: %d  from %d subjects" % (len(m), m["subject"].nunique()))
P("columns available: %s" % ", ".join([c for c in cp.columns if c.startswith("c_")]))
P()

prim = m.drop_duplicates(subset="subject", keep="first").copy()
P("primary (one scan per subject): n = %d,  with mean FD: %d"
  % (len(prim), int(prim["mean_fd"].notna().sum())))
P()


def rep(label, d, xcols, ycol):
    y = d[ycol].values.astype(float)
    X = d[xcols].values.astype(float)
    b, se, t, p, dof = ols(y, X)
    i = xcols.index("age") if "age" in xcols else 0
    rl, pl = pearsonr(d["age"].values, y)
    P("  %-30s n=%3d | age b=%+.5f t=%+.3f p=%.4f | linear r=%+.4f p=%.4f"
      % (label, len(d), b[i + 1], t[i + 1], p[i + 1], rl, pl))


P("-" * 76)
P("A. MAIN RESULT, preprocessed data (c_primary = 4mm spheres)")
P("-" * 76)
d = prim.dropna(subset=["c_primary", "age"]).copy()
d["inv_age"] = 1.0 / d["age"]
rep("linear (age only)", d, ["age"], "c_primary")
rep("age + sex", d, ["age", "sex_M"], "c_primary")
rep("inverse-age + sex", d, ["inv_age", "sex_M"], "c_primary")
dd = d.dropna(subset=["mean_fd"])
rep("age + sex + mean FD", dd, ["age", "sex_M", "mean_fd"], "c_primary")
P()

P("-" * 76)
P("B. AGREEMENT WITH THE RAW-BIDS ANALYSIS")
P("-" * 76)
raw = pd.read_csv(DER / "primary_onescan_per_subject.csv")[["subject", "z"]].rename(columns={"z": "z_raw"})
cmpr = prim.merge(raw, on="subject", how="inner").dropna(subset=["c_primary", "z_raw"])
r_agree, p_agree = pearsonr(cmpr["z_raw"], cmpr["c_primary"])
P("  r(raw-BIDS coupling, preprocessed coupling) = %+.4f, p = %.3g, n = %d"
  % (r_agree, p_agree, len(cmpr)))
P("  mean raw = %+.4f (sd %.4f);  mean preprocessed = %+.4f (sd %.4f)"
  % (cmpr.z_raw.mean(), cmpr.z_raw.std(), cmpr.c_primary.mean(), cmpr.c_primary.std()))
P()

P("-" * 76)
P("C. ROI SENSITIVITY (age regression on each coupling definition)")
P("-" * 76)
for col in [c for c in cp.columns if c.startswith("c_")]:
    dd = prim.dropna(subset=[col, "age"]).copy()
    dd["inv_age"] = 1.0 / dd["age"]
    y = dd[col].values.astype(float)
    X = dd[["age", "sex_M"]].values.astype(float)
    b, se, t, p, dof = ols(y, X)
    rl, pl = pearsonr(dd["age"].values, y)
    P("  %-12s n=%3d | age b=%+.5f t=%+.3f p=%.4f | r=%+.4f p=%.4f"
      % (col, len(dd), b[1], t[1], p[1], rl, pl))
P()
P("  -> if all p > .05 across radii 2-6 mm, seed shifts +/-2 mm, and an")
P("     alternative NAcc target, the null is not an artefact of ROI choice.")

(RES / "preproc_analysis.txt").write_text("\n".join(L), encoding="utf-8")
print("\n[saved] results/preproc_analysis.txt")
