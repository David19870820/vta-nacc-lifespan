#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
01_download_bids.py
===================
Download the raw NKI-Rockland resting-state scans (TR = 1400 ms) that this
study analysed, plus the phenotypic file needed for age.

WARNING - SIZE
    Each resting-state BOLD volume is roughly 80-120 MB. The bucket holds
    428 task-rest_acq-1400 scans for 307 participants; 417 of them yielded
    usable ROI time series and entered the analyses. A full download is
    therefore on the order of 40-50 GB. Use --limit to fetch a subset.

The data are already hosted publicly; this script simply fetches them.
Please do NOT redistribute the downloaded imaging data - see DATA_NOTES.md.

Usage
-----
    python code/01_download_bids.py --list-only          # just enumerate
    python code/01_download_bids.py --limit 5            # 5 subjects
    python code/01_download_bids.py                      # everything

Data source
-----------
    s3://fcp-indi/data/Projects/RocklandSample/RawDataBIDSLatest
    (open bucket, no credentials required)
"""
import argparse, re, sys, urllib.request, urllib.parse, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEST = ROOT / "data_raw" / "NKI_Rockland"
DEST.mkdir(parents=True, exist_ok=True)

BUCKET = "https://fcp-indi.s3.amazonaws.com/"
RAW    = "data/Projects/RocklandSample/RawDataBIDSLatest/"
PHENO  = "data/Projects/RocklandSample/aws_links.csv"


def http_get(url, timeout=120, tries=3, binary=True):
    for a in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            time.sleep(1.5 * (a + 1))
        except Exception:
            time.sleep(1.5 * (a + 1))
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


def find_targets():
    """Return [(subject, session, relpath)] for every rest_acq-1400 BOLD file."""
    subs = sorted(list_dirs(RAW, want="sub-"))
    print("subjects on the bucket: %d" % len(subs), flush=True)
    targets = []

    def one(sub):
        out = []
        for s in sorted(list_dirs(RAW + sub + "/", want="ses-")):
            keys = list_keys(RAW + sub + "/" + s + "/func/")
            for k in keys:
                if k.endswith("task-rest_acq-1400_bold.nii.gz"):
                    out.append((sub, s, sub + "/" + s + "/func/" + k))
        return out

    with ThreadPoolExecutor(max_workers=12) as ex:
        for i, r in enumerate(ex.map(one, subs)):
            targets.extend(r)
            if (i + 1) % 50 == 0:
                print("  scanned %d/%d subjects, %d scans found" % (i + 1, len(subs), len(targets)),
                      flush=True)
    return targets


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="download only the first N scans")
    ap.add_argument("--list-only", action="store_true")
    ap.add_argument("--pheno-only", action="store_true")
    a = ap.parse_args()

    # phenotypic file (small, needed for age)
    ph = DEST.parent / "aws_links.csv"
    if not ph.exists() or a.pheno_only:
        print("downloading aws_links.csv ...", flush=True)
        b = http_get(BUCKET + PHENO)
        if b:
            ph.write_bytes(b)
            print("  saved %s (%.1f MB)" % (ph, len(b) / 1024 / 1024))
        else:
            print("  FAILED to download aws_links.csv")
    else:
        print("aws_links.csv already present")
    if a.pheno_only:
        return

    targets = find_targets()
    print("rest_acq-1400 scans available: %d" % len(targets))
    if a.list_only:
        for t in targets[:20]:
            print("   ", t[2])
        return
    if a.limit:
        targets = targets[:a.limit]
        print("limiting to %d scans" % len(targets))

    done = 0
    for i, (sub, ses, rel) in enumerate(targets, 1):
        out = DEST / rel
        if out.exists() and out.stat().st_size > 1_000_000:
            done += 1
            continue
        out.parent.mkdir(parents=True, exist_ok=True)
        b = http_get(BUCKET + RAW + rel, timeout=900)
        if b:
            out.write_bytes(b)
            done += 1
            print("  [%d/%d] %s  %.1f MB" % (i, len(targets), rel, len(b) / 1024 / 1024), flush=True)
        else:
            print("  [%d/%d] %s  FAILED" % (i, len(targets), rel), flush=True)
    print("downloaded %d / %d scans -> %s" % (done, len(targets), DEST))


if __name__ == "__main__":
    main()
