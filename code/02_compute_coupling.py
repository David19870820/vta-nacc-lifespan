#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
02_compute_coupling.py
======================
Compute resting-state functional coupling between the VTA and the bilateral
NAcc for every downloaded scan.

ROI definition (Murty et al., 2014, to match the original study)
    VTA   ( 0, -16, -14)   4-mm sphere
    NAccL (-10, 12,  -8)   4-mm sphere
    NAccR ( 10, 12,  -8)   4-mm sphere

Pipeline
    1. extract mean time series per sphere          (nilearn NiftiSpheresMasker)
    2. detrend + z-standardise + band-pass 0.01-0.08 Hz
    3. Pearson r between VTA and each NAcc hemisphere, averaged
    4. Fisher z-transform; scans with < 30 volumes are dropped

Requires: nilearn, nibabel

Usage
-----
    python code/02_compute_coupling.py
"""
import glob, json, sys, warnings
from pathlib import Path
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RAW  = ROOT / "data_raw" / "NKI_Rockland"
OUT  = ROOT / "data_derived"
OUT.mkdir(parents=True, exist_ok=True)

ROIS   = [(0, -16, -14), (-10, 12, -8), (10, 12, -8)]   # VTA, NAcc L, NAcc R
RADIUS = 4.0
TR     = 1.4


def process_one(bold_path):
    from nilearn.maskers import NiftiSpheresMasker
    from nilearn.connectome import ConnectivityMeasure
    try:
        masker = NiftiSpheresMasker(seeds=ROIS, radius=RADIUS,
                                    detrend=True, standardize=True,
                                    low_pass=0.08, high_pass=0.01, t_r=TR)
        ts = masker.fit_transform(str(bold_path))
        if ts.shape[0] < 30:
            return dict(ok=False, err="fewer than 30 volumes")
        corr = ConnectivityMeasure(kind="correlation").fit_transform([ts])[0]
        r = float(np.mean([corr[0, 1], corr[0, 2]]))
        return dict(ok=True, r=r, z=float(np.arctanh(np.clip(r, -0.999999, 0.999999))),
                    n_vol=int(ts.shape[0]))
    except Exception as e:
        return dict(ok=False, err=str(e)[:120])


def main():
    pattern = str(RAW / "sub-*" / "ses-*" / "func" / "*task-rest_acq-1400_bold.nii.gz")
    files = sorted(glob.glob(pattern))
    print("found %d acq-1400 resting scans" % len(files))
    if not files:
        print("nothing to do - run 01_download_bids.py first")
        return

    try:
        from joblib import Parallel, delayed
        res = Parallel(n_jobs=-1, verbose=10)(delayed(process_one)(f) for f in files)
    except ImportError:
        res = [process_one(f) for f in files]

    rows = []
    for f, r in zip(files, res):
        p = Path(f)
        sub = p.parts[-4]
        ses = p.parts[-3]
        rows.append(dict(subject=sub, session=ses, **r))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "coupling_raw.csv", index=False)
    n_ok = int(df["ok"].sum())
    print("coupling computed for %d / %d scans -> %s" % (n_ok, len(df), OUT / "coupling_raw.csv"))

    # merge with age
    ph = RAW.parent / "aws_links.csv"
    if ph.exists():
        aws = pd.read_csv(ph)
        aws["age"] = pd.to_numeric(aws["age"], errors="coerce")
        aws = aws.dropna(subset=["age"])
        pheno = aws.groupby("subject").agg(age=("age", "first"),
                                           gender=("gender", "first")).reset_index()
        m = df[df["ok"]].copy()
        m["subject_clean"] = m["subject"].str.replace("sub-", "", regex=False)
        m = m.merge(pheno, left_on="subject_clean", right_on="subject", how="inner",
                    suffixes=("", "_pheno"))
        m.to_csv(OUT / "coupling_with_age.csv", index=False)
        print("merged with phenotype: %d rows -> %s" % (len(m), OUT / "coupling_with_age.csv"))
    else:
        print("aws_links.csv not found; skipped age merge")


if __name__ == "__main__":
    main()
