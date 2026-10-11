#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
00_quick_reproduce.py
=====================
Reproduce every statistic reported in the paper from the derived,
subject-level data shipped in data_derived/. Runs in a few seconds and
requires no imaging data.

Usage
-----
    python code/00_quick_reproduce.py

Inputs
------
    data_derived/primary_onescan_per_subject.csv   n = 307 (one scan per subject)
    data_derived/full_all_scans.csv                n = 417 (all scans)

Outputs
-------
    Console table; also written to results/quick_reproduce.txt
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from _stats import pearsonr, ols, welch_ttest          # noqa: E402

OUT = []
def P(s=""):
    print(s)
    OUT.append(str(s))


# ----------------------------------------------------------------------
# Load
# ----------------------------------------------------------------------
prim = pd.read_csv(ROOT / "data_derived" / "primary_onescan_per_subject.csv")
full = pd.read_csv(ROOT / "data_derived" / "full_all_scans.csv")

P("=" * 74)
P("VTA-NAcc resting-state coupling vs age - reproduction")
P("=" * 74)
P()
P("primary (one scan per subject) : %d subjects, %d with mean FD"
  % (len(prim), prim["mean_fd"].notna().sum()))
P("full (all scans)               : %d rows, %d with mean FD"
  % (len(full), full["mean_fd"].notna().sum()))
P()


def report(label, d, xcols, ycol="z"):
    y = d[ycol].values.astype(float)
    X = d[xcols].values.astype(float)
    b, se, t, p, dof = ols(y, X)
    parts = " | ".join("%s b=%.5f t=%.3f p=%.4f" % (c, b[i + 1], t[i + 1], p[i + 1])
                       for i, c in enumerate(xcols))
    rlin, plin = pearsonr(d["age"].values, y)
    P("  %-30s n=%3d | %s" % (label, len(d), parts))
    P("  %-30s        linear r=%.4f p=%.4f" % ("", rlin, plin))


# ----------------------------------------------------------------------
# 1. Primary analysis, no motion control (paper Table/Section 3.2)
# ----------------------------------------------------------------------
P("-" * 74)
P("1. PRIMARY ANALYSIS  (n = %d, one scan per subject)" % len(prim))
P("-" * 74)
p0 = prim.dropna(subset=["age", "z"]).copy()
p0["inv_age"] = 1.0 / p0["age"]
report("age", p0, ["age"])
report("age + sex", p0, ["age", "sex_M"])
report("inv_age + sex (Murty model)", p0, ["inv_age", "sex_M"])

r, p = pearsonr(p0["age"], p0["z"])
P()
P("  linear Pearson  r = %.4f,  p = %.4f,  n = %d" % (r, p, len(p0)))
P()

# ----------------------------------------------------------------------
# 2. Sensitivity: trimming the extremes
# ----------------------------------------------------------------------
P("-" * 74)
P("2. SENSITIVITY: excluding youngest/oldest participants")
P("-" * 74)
for pct in [0.05, 0.10]:
    lo, hi = p0["age"].quantile(pct), p0["age"].quantile(1 - pct)
    d = p0[(p0["age"] >= lo) & (p0["age"] <= hi)].copy()
    d["inv_age"] = 1.0 / d["age"]
    rl, pl = pearsonr(d["age"], d["z"])
    b, se, t, pp, dof = ols(d["z"].values, d[["inv_age", "sex_M"]].values)
    P("  exclude %2d%% each tail : n=%3d  age %2.0f-%2.0f | linear p=%.3f | inv-age p=%.3f"
      % (int(pct * 100), len(d), d["age"].min(), d["age"].max(), pl, pp[1]))
# coupling outliers
mu, sd = p0["z"].mean(), p0["z"].std()
d = p0[(p0["z"] - mu).abs() <= 3 * sd].copy()
d["inv_age"] = 1.0 / d["age"]
rl, pl = pearsonr(d["age"], d["z"])
b, se, t, pp, dof = ols(d["z"].values, d[["inv_age", "sex_M"]].values)
P("  exclude |z| > 3 SD     : n=%3d                | linear p=%.3f | inv-age p=%.3f"
  % (len(d), pl, pp[1]))
P()

# ----------------------------------------------------------------------
# 3. Head-motion diagnostics
# ----------------------------------------------------------------------
P("-" * 74)
P("3. HEAD-MOTION DIAGNOSTICS")
P("-" * 74)
md = p0.dropna(subset=["mean_fd"]).copy()
r, p = pearsonr(md["age"], md["mean_fd"])
P("  D1  r(age, mean FD)      = %+.4f,  p = %.4f,  n = %d" % (r, p, len(md)))
P("      mean FD: mean=%.3f  sd=%.3f  min=%.3f  max=%.3f"
  % (md["mean_fd"].mean(), md["mean_fd"].std(), md["mean_fd"].min(), md["mean_fd"].max()))
r, p = pearsonr(md["mean_fd"], md["z"])
P("  D2  r(mean FD, coupling) = %+.4f,  p = %.4f,  n = %d" % (r, p, len(md)))
inc = p0["mean_fd"].notna()
for col in ["age", "z"]:
    a = p0.loc[inc, col].dropna().values
    b = p0.loc[~inc, col].dropna().values
    t, pp = welch_ttest(a, b)
    P("  D3  %-3s  included %.2f (sd %.2f) | excluded %.2f (sd %.2f) | Welch t=%.2f p=%.3f"
      % (col, a.mean(), a.std(ddof=1), b.mean(), b.std(ddof=1), t, pp))
P()

# ----------------------------------------------------------------------
# 4. Motion as a covariate
# ----------------------------------------------------------------------
P("-" * 74)
P("4. AGE EFFECT WITH HEAD MOTION AS A COVARIATE")
P("-" * 74)
md["inv_age"] = 1.0 / md["age"]
report("age + sex", md, ["age", "sex_M"])
report("age + sex + mean FD", md, ["age", "sex_M", "mean_fd"])
report("inv_age + sex + mean FD", md, ["inv_age", "sex_M", "mean_fd"])
P()

# ----------------------------------------------------------------------
# 5. Excluding high-motion participants
# ----------------------------------------------------------------------
P("-" * 74)
P("5. EXCLUDING HIGH-MOTION PARTICIPANTS")
P("-" * 74)
for thr in [0.2, 0.3, 0.5, 0.7, 1.0]:
    d = md[md["mean_fd"] <= thr].copy()
    if len(d) < 30:
        continue
    rl, pl = pearsonr(d["age"], d["z"])
    b, se, t, pp, dof = ols(d["z"].values, d[["age", "sex_M"]].values)
    P("  mean FD <= %.1f mm : n=%3d | linear p=%.3f | age|sex p=%.3f"
      % (thr, len(d), pl, pp[1]))
P()

# ----------------------------------------------------------------------
# 6. Full-scan analysis, all scans (paper Section 3.2)
# ----------------------------------------------------------------------
P("-" * 74)
P("6. FULL-SCAN ANALYSIS  (all %d scans)" % len(full))
P("-" * 74)
fa = full.dropna(subset=["age", "z"]).copy()
fa["inv_age"] = 1.0 / fa["age"]
rlin, plin = pearsonr(fa["age"], fa["z"])
P("  n = %d rows from %d subjects | linear r=%.4f p=%.4f"
  % (len(fa), fa["subject"].nunique(), rlin, plin))
fas = fa.dropna(subset=["sex_M"])
b, se, t, pp, dof = ols(fas["z"].values, fas[["age", "sex_M"]].values)
P("  %-30s n=%3d | age b=%.5f t=%.3f p=%.4f  95%% CI=[%.5f, %.5f]"
  % ("age + sex (multiple regression)", len(fas), b[1], t[1], pp[1],
     b[1] - 1.959963985 * se[1], b[1] + 1.959963985 * se[1]))
report("inv_age + sex (Murty model)", fa, ["inv_age", "sex_M"])
P()

# ----------------------------------------------------------------------
# 7. Full-scan analysis restricted to scans with mappable FD (Section 3.4)
# ----------------------------------------------------------------------
P("-" * 74)
P("7. FULL-SCAN ANALYSIS with motion control")
P("-" * 74)
fu = full.dropna(subset=["mean_fd", "age", "z"]).copy()
rlin, plin = pearsonr(fu["age"], fu["z"])
P("  n = %d rows from %d subjects | linear r=%.4f p=%.4f"
  % (len(fu), fu["subject"].nunique(), rlin, plin))
report("age + sex + mean FD", fu, ["age", "sex_M", "mean_fd"])
P()

# ----------------------------------------------------------------------
P("=" * 74)
P("Coverage of this script")
P("=" * 74)
P("Reproduced here : Sections 3.1 (counts), 3.2, 3.3, 3.4, 3.5.")
P("Elsewhere       :")
P("  3.6 equivalence testing and Bayes factor  -> code/07_equivalence_lmm.py")
P("  3.7 linear mixed model                    -> code/07_equivalence_lmm.py")
P("  3.8 preprocessed recomputation            -> code/06_preproc_sensitivity.py")
P("  3.8 test-retest reliability               -> code/08_reliability.py")
P("  3.9 ROI definition sensitivity            -> code/06_preproc_sensitivity.py")
P("  Figures 1-4 and forest_values.csv         -> code/09_make_figures.py")
P("Every number printed above was regenerated from data_derived/ in this run.")
P("=" * 74)

(ROOT / "results").mkdir(exist_ok=True)
(ROOT / "results" / "quick_reproduce.txt").write_text("\n".join(OUT), encoding="utf-8")
print("\n[saved] results/quick_reproduce.txt")
