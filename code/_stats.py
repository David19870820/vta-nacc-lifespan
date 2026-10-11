# -*- coding: utf-8 -*-
"""
Dependency-free statistics helpers.

The analysis environment for this project had no SciPy available, so the
Student-t distribution is implemented here via the regularised incomplete
beta function (continued-fraction expansion, Numerical Recipes style).
Validated against known reference values:

    t = 2.0,  df = 10    -> two-sided p = 0.073388
    t = 1.96, df = 1e6   -> two-sided p = 0.049996

Only NumPy is required.
"""
import math
import numpy as np


def _betacf(a, b, x, maxit=300, eps=3e-16):
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < 1e-300:
        d = 1e-300
    d = 1.0 / d
    h = d
    for m in range(1, maxit + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-300:
            d = 1e-300
        c = 1.0 + aa / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-300:
            d = 1e-300
        c = 1.0 + aa / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < eps:
            break
    return h


def betai(a, b, x):
    """Regularised incomplete beta function I_x(a, b)."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbeta = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    bt = math.exp(lbeta + a * math.log(x) + b * math.log(1.0 - x))
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def t_two_sided_p(t, df):
    """Two-sided p-value for Student's t."""
    if df <= 0:
        return float("nan")
    x = df / (df + float(t) ** 2)
    return float(betai(df / 2.0, 0.5, x))


def pearsonr(x, y):
    """Pearson r and its two-sided p-value."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = len(x)
    xc, yc = x - x.mean(), y - y.mean()
    den = math.sqrt(float((xc ** 2).sum()) * float((yc ** 2).sum()))
    if den == 0:
        return float("nan"), float("nan")
    r = float((xc * yc).sum() / den)
    r = max(-0.999999999, min(0.999999999, r))
    if n <= 2:
        return r, float("nan")
    t = r * math.sqrt((n - 2) / (1 - r * r))
    return r, t_two_sided_p(t, n - 2)


def ols(y, X):
    """OLS with intercept. Returns beta, se, t, p, dof."""
    y = np.asarray(y, float)
    X = np.asarray(X, float)
    if X.ndim == 1:
        X = X[:, None]
    n = len(y)
    Xd = np.column_stack([np.ones(n), X])
    beta, *_ = np.linalg.lstsq(Xd, y, rcond=None)
    resid = y - Xd @ beta
    dof = n - Xd.shape[1]
    sigma2 = float(resid @ resid) / dof
    se = np.sqrt(sigma2 * np.diag(np.linalg.pinv(Xd.T @ Xd)))
    with np.errstate(divide="ignore", invalid="ignore"):
        t = beta / se
    p = np.array([t_two_sided_p(tt, dof) for tt in t])
    return beta, se, t, p, dof


def welch_ttest(a, b):
    """Welch's unequal-variance t-test."""
    a = np.asarray(a, float); a = a[~np.isnan(a)]
    b = np.asarray(b, float); b = b[~np.isnan(b)]
    na, nb = len(a), len(b)
    va, vb = a.var(ddof=1), b.var(ddof=1)
    se2 = va / na + vb / nb
    if se2 <= 0:
        return float("nan"), float("nan")
    t = (a.mean() - b.mean()) / math.sqrt(se2)
    df = se2 ** 2 / ((va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1))
    return float(t), t_two_sided_p(t, df)


if __name__ == "__main__":
    print("t=2.0,  df=10   p =", round(t_two_sided_p(2.0, 10), 6), "(expect 0.073388)")
    print("t=1.96, df=1e6  p =", round(t_two_sided_p(1.96, 1000000), 6), "(expect ~0.05)")
