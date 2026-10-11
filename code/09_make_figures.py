#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
09_make_figures.py
==================
Generate the four manuscript figures and the data table behind Figure 3,
from the derived subject-level data. No imaging data required.

  Fig1_age_vs_coupling.png    age vs coupling, with linear and inverse-age fits
  Fig2_motion_diagnostics.png (A) age vs mean FD  (B) mean FD vs coupling
  Fig3_forest.png             partial correlation of age with coupling
                              across all 15 model specifications
  Fig4_robustness.png         (A) pipeline agreement (B) BAS1-BAS2 test-retest
                              (C) 12 alternative ROI definitions

  forest_values.csv           the 15 rows plotted in Figure 3

Partial correlations are converted from the regression t statistic with
r = t / sqrt(t^2 + df), and 95% CIs are obtained on the Fisher-z scale with
se = 1/sqrt(df - 1); this reproduces the values in forest_values.csv.

Usage
-----
    python code/09_make_figures.py

Outputs
-------
    results/Fig1_age_vs_coupling.png
    results/Fig2_motion_diagnostics.png
    results/Fig3_forest.png
    results/Fig4_robustness.png
    results/forest_values.csv
"""
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from _stats import pearsonr, ols, welch_ttest          # noqa: E402

RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": "white", "savefig.dpi": 300,
                     "font.size": 9})

prim = pd.read_csv(ROOT / "data_derived" / "primary_onescan_per_subject.csv")
full = pd.read_csv(ROOT / "data_derived" / "full_all_scans.csv")
prep = pd.read_csv(ROOT / "data_derived" / "coupling_preprocessed.csv")

p0 = prim.dropna(subset=["age", "z"]).copy()
p0["inv_age"] = 1.0 / p0["age"]
md = p0.dropna(subset=["mean_fd"]).copy()
md["inv_age"] = 1.0 / md["age"]


def part_r(d, xcols, xname, ycol="z"):
    """Partial correlation of xname with y, plus 95% CI, from the regression t."""
    b, se, t, p, dof = ols(d[ycol].values.astype(float), d[xcols].values.astype(float))
    i = xcols.index(xname)
    tt = float(t[i + 1])
    r = tt / math.sqrt(tt * tt + dof)
    z = math.atanh(max(-0.999999, min(0.999999, r)))
    sez = 1.0 / math.sqrt(max(dof - 1, 2))
    return r, float(np.tanh(z - 1.959963985 * sez)), float(np.tanh(z + 1.959963985 * sez)), float(p[i + 1])


# ---------------------------------------------------------------- Figure 1
fig, ax = plt.subplots(figsize=(6.4, 4.2))
ax.scatter(p0["age"], p0["z"], s=12, alpha=0.55, edgecolors="none", color="#3b6ea5")
r, p = pearsonr(p0["age"], p0["z"])
xs = np.linspace(p0["age"].min(), p0["age"].max(), 200)
b, _, _, _, _ = ols(p0["z"].values, p0[["age"]].values)
ax.plot(xs, b[0] + b[1] * xs, color="#c0392b", lw=1.5,
        label="linear: r = %.3f, p = %.2f" % (r, p))
bi, _, _, pi_, _ = ols(p0["z"].values, p0[["inv_age"]].values)
ax.plot(xs, bi[0] + bi[1] / xs, color="#e67e22", lw=1.5, ls="--",
        label="inverse-age: p = %.2f" % pi_[1])
ax.axhline(0, color="grey", lw=0.7, ls=":")
ax.set_xlabel("Age (years)")
ax.set_ylabel("VTA\u2013NAcc coupling (Fisher z)")
ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(RESULTS / "Fig1_age_vs_coupling.png")
plt.close(fig)

# ---------------------------------------------------------------- Figure 2
fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.6))
r1, p1 = pearsonr(md["age"], md["mean_fd"])
ax = axes[0]
ax.scatter(md["age"], md["mean_fd"], s=12, alpha=0.5, edgecolors="none", color="#8e44ad")
b, _, _, _, _ = ols(md["mean_fd"].values, md[["age"]].values)
xs = np.linspace(md["age"].min(), md["age"].max(), 200)
ax.plot(xs, b[0] + b[1] * xs, color="#c0392b", lw=1.5)
ax.text(0.04, 0.95, "r = %+.3f\np = %.3f\nn = %d" % (r1, p1, len(md)),
        transform=ax.transAxes, va="top", fontsize=8)
ax.set_xlabel("Age (years)")
ax.set_ylabel("Mean FD (mm)")
ax.set_title("A", loc="left", fontweight="bold")
r2, p2 = pearsonr(md["mean_fd"], md["z"])
ax = axes[1]
ax.scatter(md["mean_fd"], md["z"], s=12, alpha=0.5, edgecolors="none", color="#16a085")
b, _, _, _, _ = ols(md["z"].values, md[["mean_fd"]].values)
xs = np.linspace(md["mean_fd"].min(), md["mean_fd"].max(), 200)
ax.plot(xs, b[0] + b[1] * xs, color="#c0392b", lw=1.5)
ax.text(0.04, 0.95, "r = %+.3f\np = %.3f\nn = %d" % (r2, p2, len(md)),
        transform=ax.transAxes, va="top", fontsize=8)
ax.set_xlabel("Mean FD (mm)")
ax.set_ylabel("VTA\u2013NAcc coupling (Fisher z)")
ax.set_title("B", loc="left", fontweight="bold")
fig.tight_layout()
fig.savefig(RESULTS / "Fig2_motion_diagnostics.png")
plt.close(fig)

# ---------------------------------------------------------------- Figure 3 + CSV
rows = []
rows.append(("Primary: linear", part_r(p0, ["age"], "age")))
rows.append(("Primary: age + sex", part_r(p0, ["age", "sex_M"], "age")))
rows.append(("Primary: inverse-age + sex", part_r(p0, ["inv_age", "sex_M"], "inv_age")))

d = p0[(p0["age"] >= p0["age"].quantile(0.05)) & (p0["age"] <= p0["age"].quantile(0.95))].copy()
d["inv_age"] = 1.0 / d["age"]
rows.append(("Exclude 5 pct tails (n={})".format(len(d)), part_r(d, ["age", "sex_M"], "age")))
d = p0[(p0["age"] >= p0["age"].quantile(0.10)) & (p0["age"] <= p0["age"].quantile(0.90))].copy()
d["inv_age"] = 1.0 / d["age"]
rows.append(("Exclude 10 pct tails (n={})".format(len(d)), part_r(d, ["age", "sex_M"], "age")))
mu, sd = p0["z"].mean(), p0["z"].std()
d = p0[(p0["z"] - mu).abs() <= 3 * sd].copy()
d["inv_age"] = 1.0 / d["age"]
rows.append(("Exclude coupling outliers (n={})".format(len(d)), part_r(d, ["age", "sex_M"], "age")))
rows.append(("+ mean FD covariate (n={})".format(len(md)), part_r(md, ["age", "sex_M", "mean_fd"], "age")))
rows.append(("Inverse-age + mean FD (n={})".format(len(md)), part_r(md, ["inv_age", "sex_M", "mean_fd"], "inv_age")))
for thr in [1.0, 0.7, 0.5, 0.3, 0.2]:
    d = md[md["mean_fd"] <= thr].copy()
    if len(d) < 30:
        continue
    rows.append(("Mean FD <= %.1f mm (n=%d)" % (thr, len(d)), part_r(d, ["age", "sex_M"], "age")))
fu = full.dropna(subset=["mean_fd", "z", "age"]).copy()
rows.append(("Full scan: linear (n={})".format(len(fu)), part_r(fu, ["age"], "age")))
rows.append(("Full scan: + sex + mean FD (n={})".format(len(fu)),
             part_r(fu, ["age", "sex_M", "mean_fd"], "age")))

pd.DataFrame([{"model": lab, "r_partial": v[0], "lo": v[1], "hi": v[2], "p": v[3]}
              for lab, v in rows]).to_csv(RESULTS / "forest_values.csv", index=False)

lab = [r[0] for r in rows]
est = np.array([r[1][0] for r in rows])
lo = np.array([r[1][1] for r in rows])
hi = np.array([r[1][2] for r in rows])
fig, ax = plt.subplots(figsize=(6.6, 0.30 * len(rows) + 1.2))
y = np.arange(len(rows))[::-1]
ax.axvspan(-0.10, 0.10, color="#f1c40f", alpha=0.13, zorder=0)
ax.errorbar(est, y, xerr=[est - lo, hi - est], fmt="o", color="#2c3e50",
            ecolor="#7f8c8d", elinewidth=0.9, capsize=2, ms=4, zorder=3)
ax.axvline(0, color="#c0392b", lw=1.0, ls="--", zorder=2)
ax.set_yticks(y)
ax.set_yticklabels(lab, fontsize=7.5)
ax.set_xlabel("Partial correlation of age with coupling (95% CI)")
ax.set_xlim(-0.42, 0.42)
fig.tight_layout()
fig.savefig(RESULTS / "Fig3_forest.png")
plt.close(fig)

# ---------------------------------------------------------------- Figure 4
ok = prep[prep["c_primary"].notna()]
K = p0[["subject", "session", "age", "sex_M", "z"]].merge(
    ok, on=["subject", "session"], how="inner").dropna(subset=["z", "c_primary"])
fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.4))
lim = [-0.9, 0.9]
ax = axes[0]
r, p = pearsonr(K["z"].values, K["c_primary"].values)
ax.scatter(K["z"], K["c_primary"], s=9, alpha=0.55, edgecolors="none", color="#2c7fb8")
ax.plot(lim, lim, color="grey", lw=0.7, ls=":")
b, _, _, _, _ = ols(K["c_primary"].values, K[["z"]].values)
xs = np.linspace(K["z"].min(), K["z"].max(), 100)
ax.plot(xs, b[0] + b[1] * xs, color="#c0392b", lw=1.2)
ax.text(0.04, 0.95, "r = %+.3f\np = %.4f\nn = %d" % (r, p, len(K)),
        transform=ax.transAxes, va="top", fontsize=8)
ax.set_xlabel("Raw BIDS coupling (z)")
ax.set_ylabel("C-PAC preprocessed coupling")
ax.set_title("A", loc="left", fontweight="bold")
ax.set_xlim(*lim); ax.set_ylim(*lim)

ax = axes[1]
bas = full[full["session"].isin(["ses-BAS1", "ses-BAS2"])].dropna(subset=["z"])
piv = bas.pivot_table(index="subject", columns="session", values="z", aggfunc="first").dropna()
r2, p2 = pearsonr(piv["ses-BAS1"].values, piv["ses-BAS2"].values)
ax.scatter(piv["ses-BAS1"], piv["ses-BAS2"], s=18, alpha=0.75, edgecolors="none", color="#7b3294")
ax.plot(lim, lim, color="grey", lw=0.7, ls=":")
b, _, _, _, _ = ols(piv["ses-BAS2"].values, piv[["ses-BAS1"]].values)
xs = np.linspace(piv["ses-BAS1"].min(), piv["ses-BAS1"].max(), 50)
ax.plot(xs, b[0] + b[1] * xs, color="#c0392b", lw=1.2)
ax.text(0.04, 0.95, "r = %+.3f\np = %.3f\nn = %d" % (r2, p2, len(piv)),
        transform=ax.transAxes, va="top", fontsize=8)
ax.set_xlabel("BAS1 coupling (z)")
ax.set_ylabel("BAS2 coupling (z)")
ax.set_title("B", loc="left", fontweight="bold")
ax.set_xlim(*lim); ax.set_ylim(*lim)

ax = axes[2]
LBL = {"c_rad2": "r=2 mm", "c_rad3": "r=3 mm", "c_primary": "r=4 mm",
       "c_rad5": "r=5 mm", "c_rad6": "r=6 mm", "c_vta_xp": "VTA +2x",
       "c_vta_xm": "VTA \u22122x", "c_vta_yp": "VTA +2y", "c_vta_ym": "VTA \u22122y",
       "c_vta_zp": "VTA +2z", "c_vta_zm": "VTA \u22122z", "c_nacc_w": "NAcc \u00b112"}
ORDER = ["c_rad2", "c_rad3", "c_primary", "c_rad5", "c_rad6", "c_vta_xp", "c_vta_xm",
         "c_vta_yp", "c_vta_ym", "c_vta_zp", "c_vta_zm", "c_nacc_w"]
ys = np.arange(len(ORDER))[::-1]
for yv, col in zip(ys, ORDER):
    if col not in K.columns:
        continue
    d2 = K.dropna(subset=[col, "age"])
    r3, p3 = pearsonr(d2["age"].values, d2[col].values)
    z = math.atanh(max(-0.999999, min(0.999999, r3)))
    se = 1.0 / math.sqrt(len(d2) - 3)
    l, h = float(np.tanh(z - 1.959963985 * se)), float(np.tanh(z + 1.959963985 * se))
    ax.errorbar([r3], [yv], xerr=[[r3 - l], [h - r3]], fmt="o", color="#2c3e50",
                ecolor="#7f8c8d", capsize=2, ms=3.5, elinewidth=0.9)
ax.axvline(0, color="#c0392b", lw=1.0, ls="--")
ax.axvspan(-0.10, 0.10, color="#f1c40f", alpha=0.13, zorder=0)
ax.set_yticks(ys)
ax.set_yticklabels([LBL[c] for c in ORDER], fontsize=7.5)
ax.set_xlabel("Age\u2013coupling correlation (95% CI)")
ax.set_title("C", loc="left", fontweight="bold")
ax.set_xlim(-0.35, 0.35)
fig.tight_layout()
fig.savefig(RESULTS / "Fig4_robustness.png")
plt.close(fig)

print("[saved] results/Fig1_age_vs_coupling.png")
print("[saved] results/Fig2_motion_diagnostics.png")
print("[saved] results/Fig3_forest.png")
print("[saved] results/Fig4_robustness.png")
print("[saved] results/forest_values.csv")
