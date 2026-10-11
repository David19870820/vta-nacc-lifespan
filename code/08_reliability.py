#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
08_reliability.py
=================
Reproduce the test-retest reliability estimates reported in Section 3.8 of
the manuscript (and used for the attenuation correction in Section 4.6).

Reliability is quantified with the standard one-way random-effects
intraclass correlation, ICC(1,1) (Shrout & Fleiss, 1979, Case 1):

    ICC(1,1) = (MSB - MSW) / (MSB + (k - 1) * MSW)

where MSB and MSW are the between- and within-subject mean squares from a
one-way ANOVA of coupling on subject, and k is the number of measurements
per subject. For the unbalanced case here we use k = S / G (total
measurements divided by number of subjects), which is the usual average-k
approximation; for the balanced BAS1/BAS2 case k is exactly 2.

Four estimates are reported, because they differ substantially and the
manuscript reports all of them:

  A. all sessions, subjects with >= 2 sessions, any session type
  B. BAS1-BAS2 pairs only (same session type, the type most primary scans
     come from)
  C. the motion-matched subset used by the linear mixed model
  D. the C-PAC preprocessed coupling, where both sessions were available

The script also reports the attenuation correction used in Section 4.6:
an observed correlation r is attenuated by sqrt(reliability), so the
equivalence bound of r = 0.136 corresponds to a true bound of
0.136 / sqrt(ICC).

Usage
-----
    python code/08_reliability.py

Outputs
-------
    results/reliability.json
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from _stats import pearsonr          # noqa: E402

RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)


def icc_1_1(groups):
    """One-way random-effects ICC(1,1) from a list of per-subject arrays."""
    groups = [np.asarray(g, float) for g in groups if len(g) >= 2]
    G = len(groups)
    if G < 2:
        return None
    S = sum(len(g) for g in groups)
    if S - G < 1:
        return None
    grand = float(np.concatenate(groups).mean())
    ssb = sum(len(g) * (float(g.mean()) - grand) ** 2 for g in groups)
    ssw = sum(float(((g - g.mean()) ** 2).sum()) for g in groups)
    msb = ssb / (G - 1)
    msw = ssw / (S - G)
    k = S / G
    denom = msb + (k - 1.0) * msw
    icc = (msb - msw) / denom if denom > 0 else float("nan")
    return {"n_subj": G, "n_meas": S, "k": k,
            "msb": msb, "msw": msw, "icc": icc}


prim = pd.read_csv(ROOT / "data_derived" / "primary_onescan_per_subject.csv")
full = pd.read_csv(ROOT / "data_derived" / "full_all_scans.csv")
prep = pd.read_csv(ROOT / "data_derived" / "coupling_preprocessed.csv")

out = []

# ---------------------------------------------------------------- A. all sessions
d = full[full["session"].notna() & full["z"].notna()].copy()
groups = [g["z"].values for _, g in d.groupby("subject") if len(g) >= 2]
a = icc_1_1(groups)
n_pairs = sum(len(g) * (len(g) - 1) // 2 for g in groups)
a.update({"label": "Raw BIDS, all session pairs (any session type)",
          "n_pairs": n_pairs,
          "note": "ICC uses a one-way ANOVA over subjects, k = S/G; "
                  "n_pairs counts all pairwise session combinations"})
out.append(a)

# ---------------------------------------------------------------- B. BAS1-BAS2
bas = d[d["session"].isin(["ses-BAS1", "ses-BAS2"])]
piv = bas.pivot_table(index="subject", columns="session", values="z",
                      aggfunc="first").dropna()
r_ab, p_ab = pearsonr(piv["ses-BAS1"].values, piv["ses-BAS2"].values)
b = icc_1_1([[x, y] for x, y in zip(piv["ses-BAS1"].values, piv["ses-BAS2"].values)])
b.update({"label": "Raw BIDS, BAS1-BAS2 pairs (same session type)",
          "n_pairs": len(piv), "r": r_ab, "p": p_ab})
out.append(b)

# ---------------------------------------------------------------- C. LMM subset
dm = full.dropna(subset=["mean_fd", "z"]).copy()
gm = [g["z"].values for _, g in dm.groupby("subject") if len(g) >= 2]
c = icc_1_1(gm)
c.update({"label": "LMM subset (motion-matched, >= 2 scans)",
          "n_scans": sum(len(g) for g in gm)})
out.append(c)

# ---------------------------------------------------------------- D. preprocessed
ok = prep[prep["c_primary"].notna()].copy()
# preprocessed coupling is only available for sessions C-PAC processed;
# require >= 2 such sessions per participant
gp = [g["c_primary"].values for _, g in ok.groupby("subject") if len(g) >= 2]
dd = icc_1_1(gp)
if dd is not None:
    dd.update({"label": "C-PAC preprocessed coupling, >= 2 sessions",
               "n_scans": sum(len(g) for g in gp)})
    out.append(dd)

# ---------------------------------------------------------------- attenuation
ATT = []
for rec in out:
    if not rec or not np.isfinite(rec.get("icc", float("nan"))):
        continue
    icc = rec["icc"]
    if icc <= 0:
        ATT.append({"label": rec["label"], "icc": icc,
                    "attenuation_factor": None,
                    "true_bound_from_0.136": None})
        continue
    factor = 1.0 / math.sqrt(icc)
    ATT.append({"label": rec["label"], "icc": icc,
                "attenuation_factor": factor,
                "true_bound_from_0.136": 0.1357 * factor})

payload = {"estimates": out, "attenuation": ATT,
           "note": "observed r is attenuated by sqrt(reliability); the true "
                   "bound is the observed equivalence bound divided by "
                   "sqrt(ICC)"}
(RESULTS / "reliability.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

for rec in out:
    print("%-52s n_subj=%-4s icc=%+.4f" % (rec["label"][:52], rec["n_subj"], rec["icc"]))
print()
for rec in ATT:
    if rec["attenuation_factor"]:
        print("  %-52s ICC=%.3f  -> true bound = %.3f"
              % (rec["label"][:52], rec["icc"], rec["true_bound_from_0.136"]))
print("\n[saved] results/reliability.json")
