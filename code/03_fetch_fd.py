#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
03_fetch_fd.py
==============
Retrieve head-motion estimates (framewise displacement, FD) for the
TR = 1400 ms resting-state scans from the C-PAC derivatives of the
NKI-Rockland sample.

Why a separate step
-------------------
The raw BIDS release used in 02_compute_coupling.py does not ship motion
parameters. C-PAC derivatives DO provide FD.1D, including for the
acq-1400 resting scan:

    Outputs/C-PAC/<sub>/output/pipeline_analysis_nuisance/
        <sub>_<ses>/frame_wise_displacement_power/_scan_rest_acq-1400/FD.1D

Coverage note
-------------
C-PAC did not process every session. In this dataset FD for the
acq-1400 resting scan was retrievable for 307 subjects, of which 250 had
FD for their primary (first) scan. Participants whose session is missing
from C-PAC are simply left without a mean_fd value; they are retained in
the primary analysis but excluded from the motion-controlled analyses.

Usage
-----
    python code/03_fetch_fd.py
"""
import re, sys, time, urllib.request, urllib.parse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT  = ROOT / "data_derived"
OUT.mkdir(parents=True, exist_ok=True)

BUCKET = "https://fcp-indi.s3.amazonaws.com/"
CP     = "data/Projects/RocklandSample/Outputs/C-PAC/"
RAW    = "data/Projects/RocklandSample/RawDataBIDSLatest/"


def http_get(url, timeout=60, tries=3):
    for a in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            time.sleep(0.8 * (a + 1))
        except Exception:
            time.sleep(0.8 * (a + 1))
    return None


def list_dirs(prefix, want=None):
    q = urllib.parse.urlencode({"list-type": "2", "prefix": prefix,
                                "delimiter": "/", "max-keys": "1000"})
    b = http_get(BUCKET + "?" + q)
    if b is None:
        return []
    xml = b.decode("utf-8", "replace")
    out = [d[len(prefix):].rstrip("/") for d in re.findall(r"<Prefix>(.*?)</Prefix>", xml)]
    return [x for x in out if x.startswith(want)] if want else out


def list_keys(prefix):
    q = urllib.parse.urlencode({"list-type": "2", "prefix": prefix,
                                "delimiter": "/", "max-keys": "1000"})
    b = http_get(BUCKET + "?" + q)
    if b is None:
        return []
    return [k[len(prefix):] for k in re.findall(r"<Key>(.*?)</Key>", b.decode("utf-8", "replace"))]


def sessions_with_rest1400(sub):
    """Sessions of this subject that contain a rest_acq-1400 BOLD scan.

    Sorted alphabetically, which is the order used by 02_compute_coupling.py
    (glob + sorted), so rows and sessions line up one-to-one.
    """
    out = []
    for s in sorted(list_dirs(RAW + sub + "/", want="ses-")):
        keys = list_keys(RAW + sub + "/" + s + "/func/")
        if any(k.endswith("task-rest_acq-1400_bold.nii.gz") for k in keys):
            out.append(s)
    return out


def fd_for(sub, ses):
    key = ("%s%s/output/pipeline_analysis_nuisance/%s_%s/"
           "frame_wise_displacement_power/_scan_rest_acq-1400/FD.1D" % (CP, sub, sub, ses))
    b = http_get(BUCKET + key)
    if b is None:
        return None
    vals = [float(v) for v in b.decode("utf-8", "replace").split() if v.strip()]
    if not vals:
        return None
    a = np.asarray(vals, float)
    return dict(mean_fd=float(a.mean()), max_fd=float(a.max()),
                n_vol=int(a.size), pct_over_0p5=float((a > 0.5).mean()))


def one(sub):
    rec = {"subject": sub}
    ses = sessions_with_rest1400(sub)
    rec["raw_sessions_1400"] = ";".join(ses)
    rec["n_raw_1400"] = len(ses)
    fds = {}
    for s in ses:
        d = fd_for(sub, s)
        if d:
            fds[s] = d
    rec["n_fd_found"] = len(fds)
    for s in ses:
        if s in fds:
            f = fds[s]
            rec[s + "_mean_fd"] = round(f["mean_fd"], 4)
            rec[s + "_max_fd"]  = round(f["max_fd"], 4)
            rec[s + "_nvol"]    = f["n_vol"]
            rec[s + "_pctO5"]   = round(f["pct_over_0p5"], 3)
    rec["first_ses"] = ses[0] if ses else ""
    rec["first_mean_fd"] = (round(fds[ses[0]]["mean_fd"], 4)
                            if ses and ses[0] in fds else np.nan)
    allfd = [v["mean_fd"] for v in fds.values()]
    rec["mean_fd_all_ses"] = round(float(np.mean(allfd)), 4) if allfd else np.nan
    return rec


def main():
    src = OUT / "coupling_with_age.csv"
    if src.exists():
        df = pd.read_csv(src)
        subs = sorted(df["subject"].dropna().unique())
    else:
        print("coupling_with_age.csv not found; falling back to a bucket listing")
        subs = sorted(list_dirs(CP, want="sub-"))
    print("fetching FD for %d subjects ..." % len(subs), flush=True)

    t0 = time.time()
    rows = []
    with ThreadPoolExecutor(max_workers=12) as ex:
        for i, r in enumerate(ex.map(one, subs)):
            rows.append(r)
            if (i + 1) % 50 == 0:
                print("  %d/%d  %.0fs" % (i + 1, len(subs), time.time() - t0), flush=True)

    res = pd.DataFrame(rows)
    res.to_csv(OUT / "fd_by_scan.csv", index=False)
    print("saved fd_by_scan.csv (%d subjects, %d with first-scan FD)"
          % (len(res), int(res["first_mean_fd"].notna().sum())))


if __name__ == "__main__":
    main()
