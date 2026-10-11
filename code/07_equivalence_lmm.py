#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
07_equivalence_lmm.py
=====================
Reproduce Section 3.6 (equivalence testing and Bayes factor) and
Section 3.7 (linear mixed model) from the derived, subject-level data.

1. Two one-sided tests (TOST) for equivalence on the primary age-coupling
   correlation, performed on the Fisher-z scale. For an equivalence bound
   b the two one-sided p-values are

       p_plus  = Phi( (z - atanh(b)) / se )
       p_minus = Phi( -(z + atanh(b)) / se )

   with z = atanh(r), se = 1/sqrt(n-3); the TOST p-value is their maximum.
   The smallest supported bound is found by bisection.

2. Bayes factor for H0: rho = 0 against a two-sided alternative, using the
   Jeffreys-Zellner-Siow default prior as implemented for correlations
   (Rouder et al. 2009; Ly et al. 2016). The correlation maps to
   t = r*sqrt(n-2)/sqrt(1-r^2) with df = n-2 and effective N = n-2, and
   the BF integrates over the inverse-gamma prior on g with scale rscale:

       BF10 = int_0^inf (1+N g)^(-1/2)
                        * [1 + t^2/((1+N g) df)]^(-(df+1)/2)
                        / [1 + t^2/df]^(-(df+1)/2)
                        * rscale/(sqrt(2 pi)) g^(-3/2) exp(-rscale^2/(2g)) dg

   Two prior scales are reported (0.5 and 1.0), as in the manuscript.

3. Linear mixed model with a random intercept per participant, fitted by
   REML over all scans with usable motion estimates. The variance ratio
   lambda = sigma_u^2/sigma_e^2 is profiled on a grid; fixed-effect SEs use
   the REML estimate of sigma_e^2 and denominator df = G-1 (participants
   minus one), the conservative choice for a random-intercept model.

Only NumPy is required (SciPy is deliberately not a dependency of this
repository); the normal CDF uses math.erf.

Usage
-----
    python code/07_equivalence_lmm.py

Outputs
-------
    results/equiv_and_lmm.txt    human-readable log
    results/equiv_and_lmm.json   machine-readable side-car
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
from _stats import t_two_sided_p          # noqa: E402

OUT = []


def P(s=""):
    print(s)
    OUT.append(str(s))


def phi(x):
    """Standard normal CDF."""
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


# ----------------------------------------------------------------------
# Load
# ----------------------------------------------------------------------
prim = pd.read_csv(ROOT / "data_derived" / "primary_onescan_per_subject.csv")
full = pd.read_csv(ROOT / "data_derived" / "full_all_scans.csv")

prim = prim.dropna(subset=["age", "z"]).copy()
r_obs = float(np.corrcoef(prim["age"].values, prim["z"].values)[0, 1])
n_obs = len(prim)
z_obs = math.atanh(max(-0.999999999, min(0.999999999, r_obs)))
se_z = 1.0 / math.sqrt(n_obs - 3)


def ci_r(z, se, level):
    crit = {0.95: 1.959963985, 0.90: 1.644853627}[level]
    return (math.tanh(z - crit * se), math.tanh(z + crit * se))


# ----------------------------------------------------------------------
# 1. Equivalence testing (TOST)
# ----------------------------------------------------------------------
P("=" * 74)
P("1.  EQUIVALENCE TESTING (TOST), primary age-coupling correlation")
P("=" * 74)
P("observed r = %.4f,  n = %d" % (r_obs, n_obs))
lo95, hi95 = ci_r(z_obs, se_z, 0.95)
lo90, hi90 = ci_r(z_obs, se_z, 0.90)
P("95%% CI = [%+.3f, %+.3f]" % (lo95, hi95))
P("90%% CI = [%+.3f, %+.3f]   (this is what TOST at alpha=.05 uses)" % (lo90, hi90))
P()
P("  %-7s %-16s %-16s %-13s %s"
  % ("bound", "one-sided p(+)", "one-sided p(-)", "TOST p", "equivalent?"))
BOUNDS = [0.05, 0.08, 0.10, 0.12, 0.13, 0.14, 0.15, 0.20]
tost_rows = []
for b in BOUNDS:
    ab = math.atanh(b)
    p_plus = phi((z_obs - ab) / se_z)
    p_minus = phi(-(z_obs + ab) / se_z)
    tost = max(p_plus, p_minus)
    ok = tost < 0.05
    tost_rows.append((b, p_plus, p_minus, tost, ok))
    P("  %-7.2f %-16.5f %-16.5f %-13.5f %s"
      % (b, p_plus, p_minus, tost, "YES" if ok else "no"))


def tost_p(b):
    ab = math.atanh(b)
    return max(phi((z_obs - ab) / se_z), phi(-(z_obs + ab) / se_z))


lo_b, hi_b = 0.0, 0.9
for _ in range(200):
    mid = 0.5 * (lo_b + hi_b)
    if tost_p(mid) < 0.05:
        hi_b = mid
    else:
        lo_b = mid
bound_star = hi_b
P()
P("  -> smallest equivalence bound that is supported:  r = %.4f" % bound_star)
P("     i.e. we can exclude |rho| >= %.3f, but NOT small effects below that."
  % bound_star)
P()

# ----------------------------------------------------------------------
# 2. Bayes factor (JZS default prior for the correlation)
# ----------------------------------------------------------------------
P("=" * 74)
P("2.  BAYES FACTOR  H0: rho = 0   vs   H1: rho != 0")
P("=" * 74)


def jzs_bf10(r, n, rscale, n_grid=400000):
    """JZS Bayes factor for a correlation (Rouder et al. 2009; Ly et al. 2016)."""
    if n <= 2:
        return float("nan")
    df = n - 2
    neff = n - 2
    r = max(-0.999999999, min(0.999999999, r))
    t = r * math.sqrt(df) / math.sqrt(1.0 - r * r)
    g = np.exp(np.linspace(math.log(1e-10), math.log(1e7), n_grid))
    ng = neff * g
    num = (1.0 + ng) ** -0.5 * (1.0 + t * t / ((1.0 + ng) * df)) ** (-(df + 1) / 2.0)
    den = (1.0 + t * t / df) ** (-(df + 1) / 2.0)
    prior = (rscale / math.sqrt(2.0 * math.pi)) * g ** -1.5 \
        * np.exp(-rscale * rscale / (2.0 * g))
    return float(np.trapezoid((num / den) * prior, g))


bf = {}
for gamma in (0.5, 1.0):
    bf10 = jzs_bf10(r_obs, n_obs, gamma)
    bf[gamma] = (bf10, 1.0 / bf10)
    P("  Cauchy prior scale %.1f :  BF01 = %6.2f   (BF10 = %.3f)"
      % (gamma, 1.0 / bf10, bf10))
P("  Jeffreys: BF01 3-10 = moderate evidence for H0; >10 = strong")
P()

# ----------------------------------------------------------------------
# 3. Linear mixed model, random intercept per participant
# ----------------------------------------------------------------------
P("=" * 74)
P("3.  LINEAR MIXED MODEL, random intercept per participant (all scans)")
P("=" * 74)

d = full.dropna(subset=["age", "sex_M", "mean_fd", "z"]).copy()
y = d["z"].values.astype(float)
X = np.column_stack([np.ones(len(d)), d["age"].values, d["sex_M"].values,
                     d["mean_fd"].values])
gidx = pd.factorize(d["subject"].values)[0]
G = int(gidx.max()) + 1
N = len(d)
Z = np.zeros((N, G))
Z[np.arange(N), gidx] = 1.0
I = np.eye(N)


def reml_at(lam):
    """Profile REML log-likelihood at variance ratio lam = sigma_u^2/sigma_e^2."""
    V = lam * (Z @ Z.T) + I
    Vi = np.linalg.inv(V)
    XtViX = X.T @ Vi @ X
    sign_v, ld_v = np.linalg.slogdet(V)
    sign_x, ld_x = np.linalg.slogdet(XtViX)
    if sign_v <= 0 or sign_x <= 0:
        return -np.inf, None, None
    beta = np.linalg.solve(XtViX, X.T @ Vi @ y)
    resid = y - X @ beta
    quad = float(resid.T @ Vi @ resid)
    return -0.5 * (ld_v + ld_x + quad), beta, Vi


grid = np.concatenate([np.linspace(0.0, 1.0, 2001), np.linspace(1.0, 30.0, 581)])
best = None
for lam in grid:
    ll, beta, Vi = reml_at(lam)
    if best is None or ll > best[0]:
        best = (ll, lam, beta, Vi)
ll, lam_hat, beta, Vi = best
icc = lam_hat / (1.0 + lam_hat) if lam_hat > 0 else 0.0

resid = y - X @ beta
sigma_e2 = float(resid.T @ Vi @ resid) / (N - X.shape[1])
se_beta = np.sqrt(np.diag(np.linalg.inv(X.T @ Vi @ X)) * sigma_e2)
df_lmm = G - 1

P("  N = %d scans,  G = %d participants" % (N, G))
P("  sigma_u^2/sigma_e^2 (REML) = %.4f   -> ICC = %.3f" % (lam_hat, icc))
P("  sigma_e = %.4f" % math.sqrt(sigma_e2))
P()
P("  fixed effects (df = G-1 = %d):" % df_lmm)
names = ["(intercept)", "age", "sex_M", "mean_fd"]
for i, nm in enumerate(names):
    t = beta[i] / se_beta[i]
    p = t_two_sided_p(t, df_lmm)
    P("    %-12s b = %+.5f   se = %.5f   t = %+.3f   p = %.4f"
      % (nm, beta[i], se_beta[i], t, p))
P()

t_age = float(beta[1] / se_beta[1])
r_part = t_age / math.sqrt(t_age * t_age + df_lmm)
P("  age effect as partial r (for comparison with Fig. 3):  r = %+.4f" % r_part)
P("  note: ICC = 0 means between-person variance in coupling is fully")
P("        absorbed by the fixed effects; the LMM then reduces to OLS")
P("        with cluster-robust df. Reported here for completeness.")
P()

(ROOT / "results").mkdir(exist_ok=True)
(ROOT / "results" / "equiv_and_lmm.txt").write_text("\n".join(OUT), encoding="utf-8")
(ROOT / "results" / "equiv_and_lmm.json").write_text(
    json.dumps({
        "primary_r": r_obs,
        "n": n_obs,
        "ci95": [lo95, hi95],
        "ci90": [lo90, hi90],
        "tost": [{"bound": b, "p_plus": pp, "p_minus": pm, "tost_p": tp,
                  "equivalent": bool(ok)} for b, pp, pm, tp, ok in tost_rows],
        "smallest_supported_bound": bound_star,
        "bayes": {str(k): {"bf01": v[1], "bf10": v[0]} for k, v in bf.items()},
        "lmm": {
            "n_scans": N, "n_participants": G,
            "variance_ratio": float(lam_hat), "icc": float(icc),
            "sigma_e": math.sqrt(sigma_e2),
            "fixed_effects": {nm: {"b": float(beta[i]), "se": float(se_beta[i]),
                                   "t": float(beta[i] / se_beta[i]),
                                   "p": float(t_two_sided_p(beta[i] / se_beta[i], df_lmm))}
                              for i, nm in enumerate(names)},
        },
    }, indent=2), encoding="utf-8")

print("\n[saved] results/equiv_and_lmm.txt")
print("[saved] results/equiv_and_lmm.json")

