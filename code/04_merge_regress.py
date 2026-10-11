#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
04_merge_regress.py
===================
Merge coupling with phenotype and head motion, then run every regression
reported in the manuscript. Produces the derived tables shipped in
data_derived/ and the log in results/motion_analysis.txt.

Inputs (from earlier steps)
    data_derived/coupling_with_age.csv   (02)
    data_derived/fd_by_scan.csv          (03)
    data_raw/aws_links.csv               (01, phenotype)

Outputs
    data_derived/primary_onescan_per_subject.csv
    data_derived/full_all_scans.csv
    results/motion_analysis.txt

Usage
-----
    python code/04_merge_regress.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from _stats import pearsonr, ols, welch_ttest      # noqa: E402

DER = ROOT / "data_derived"
RES = ROOT / "results"
RES.mkdir(exist_ok=True)

L = []
def P(s=""):
    print(s)
    L.append(str(s))


df = pd.read_csv(DER / "coupling_with_age.csv")
fd = pd.read_csv(DER / "fd_by_scan.csv")

df["age"] = pd.to_numeric(df["age"], errors="coerce")
df["z"]   = pd.to_numeric(df["z"],   errors="coerce")
df["sex_M"] = (df["gender"].astype(str).str.upper().str[0] == "M").astype(int)

# ---------------- primary table: one scan per subject ----------------
first = df.drop_duplicates(subset="subject", keep="first").copy()
cols = [c for c in ["subject", "first_ses", "first_mean_fd",
                    "n_raw_1400", "n_coupling_rows"] if c in fd.columns]
prim = first.merge(fd[cols], on="subject", how="left")
prim = (prim[["subject", "age", "sex_M", "z", "r", "first_ses", "first_mean_fd"]]
        .rename(columns={"first_ses": "session", "first_mean_fd": "mean_fd"})
        .sort_values("subject").reset_index(drop=True))
prim.to_csv(DER / "primary_onescan_per_subject.csv", index=False, float_format="%.6f")

# ---------------- full table: every scan with a mappable FD ----------------
rows = []
for sub, g in df.groupby("subject", sort=False):
    g = g.copy()
    f = fd[fd["subject"] == sub]
    nraw = ncou = np.nan
    seslist = []
    if len(f):
        r0 = f.iloc[0]
        try:
            nraw, ncou = float(r0["n_raw_1400"]), float(r0["n_coupling_rows"])
        except Exception:
            pass
        if pd.notna(r0.get("raw_sessions_1400")):
            seslist = str(r0["raw_sessions_1400"]).split(";")
    if len(f) and nraw == ncou and len(seslist) == len(g):
        g["session"] = seslist
        g["mean_fd"] = [float(f.iloc[0][s + "_mean_fd"])
                        if (s + "_mean_fd") in f.columns and pd.notna(f.iloc[0][s + "_mean_fd"])
                        else np.nan for s in seslist]
    else:
        g["session"] = ""
        g["mean_fd"] = np.nan
    rows.append(g)
full = (pd.concat(rows, ignore_index=True)
        [["subject", "age", "sex_M", "z", "r", "session", "mean_fd"]]
        .sort_values(["subject", "session"]).reset_index(drop=True))
full.to_csv(DER / "full_all_scans.csv", index=False, float_format="%.6f")

P("primary rows: %d  (with FD: %d)" % (len(prim), int(prim["mean_fd"].notna().sum())))
P("full rows   : %d  (with FD: %d)" % (len(full), int(full["mean_fd"].notna().sum())))


def report(label, d, xcols, ycol="z"):
    b, se, t, p, dof = ols(d[ycol].values.astype(float), d[xcols].values.astype(float))
    parts = " | ".join("%s b=%.5f t=%.3f p=%.4f" % (c, b[i + 1], t[i + 1], p[i + 1])
                       for i, c in enumerate(xcols))
    rlin, plin = pearsonr(d["age"].values, d[ycol].values)
    P("  [%-26s] n=%3d | %s | linear r=%.4f p=%.4f"
      % (label, len(d), parts, rlin, plin))


p0 = prim.dropna(subset=["age", "z"]).copy()
p0["inv_age"] = 1.0 / p0["age"]
md = p0.dropna(subset=["mean_fd"]).copy()
md["inv_age"] = 1.0 / md["age"]

P()
P("=== PRIMARY (no motion control) ===")
report("age + sex", p0, ["age", "sex_M"])
report("inv_age + sex", p0, ["inv_age", "sex_M"])

P()
P("=== HEAD-MOTION DIAGNOSTICS ===")
r, p = pearsonr(md["age"], md["mean_fd"])
P("  D1  r(age, meanFD)      = %+.4f, p = %.4g, n = %d" % (r, p, len(md)))
r, p = pearsonr(md["mean_fd"], md["z"])
P("  D2  r(meanFD, coupling) = %+.4f, p = %.4g, n = %d" % (r, p, len(md)))
inc = p0["mean_fd"].notna()
for col in ["age", "z"]:
    a = p0.loc[inc, col].dropna().values
    b = p0.loc[~inc, col].dropna().values
    t, pp = welch_ttest(a, b)
    P("  D3  %-3s incl %.2f (sd %.2f) | excl %.2f (sd %.2f) | Welch t=%.2f p=%.3f"
      % (col, a.mean(), a.std(ddof=1), b.mean(), b.std(ddof=1), t, pp))

P()
P("=== MOTION AS COVARIATE ===")
report("age + sex + meanFD", md, ["age", "sex_M", "mean_fd"])
report("inv_age + sex + meanFD", md, ["inv_age", "sex_M", "mean_fd"])

P()
P("=== EXCLUDING HIGH MOTION ===")
for thr in [0.2, 0.3, 0.5, 0.7, 1.0]:
    d = md[md["mean_fd"] <= thr]
    if len(d) < 30:
        continue
    report("meanFD <= %.1f" % thr, d, ["age", "sex_M"])

P()
P("=== FULL-SCAN ===")
fu = full.dropna(subset=["mean_fd", "age", "z"])
report("age + sex + meanFD", fu, ["age", "sex_M", "mean_fd"])

(RES / "motion_analysis.txt").write_text("\n".join(L), encoding="utf-8")
print("\n[saved] results/motion_analysis.txt")
