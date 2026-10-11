#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""批量并行下载 + 提取（可续跑）"""
import os, sys, time, socket, traceback
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np, pandas as pd, nibabel as nib
from scipy.signal import butter, filtfilt, detrend as sp_detrend
import urllib.request

socket.setdefaulttimeout(45)

ROOT = ROOT = Path(__file__).resolve().parent.parent / "data_derived"
OUT  = ROOT / "preproc"; OUT.mkdir(parents=True, exist_ok=True)
TMP  = ROOT / "_tmp_nii"; TMP.mkdir(parents=True, exist_ok=True)
BUCKET = "https://fcp-indi.s3.amazonaws.com/"
CP     = "data/Projects/RocklandSample/Outputs/C-PAC/"
TR     = 1.4
NWORK  = 6

ROIS = {
    "VTA_r2": ((0,-16,-14),2.0), "VTA_r3": ((0,-16,-14),3.0), "VTA_r4": ((0,-16,-14),4.0),
    "VTA_r5": ((0,-16,-14),5.0), "VTA_r6": ((0,-16,-14),6.0),
    "VTA_xp": ((2,-16,-14),4.0), "VTA_xm": ((-2,-16,-14),4.0),
    "VTA_yp": ((0,-14,-14),4.0), "VTA_ym": ((0,-18,-14),4.0),
    "VTA_zp": ((0,-16,-12),4.0), "VTA_zm": ((0,-16,-16),4.0),
    "NAccL_r2": ((-10,12,-8),2.0), "NAccR_r2": ((10,12,-8),2.0),
    "NAccL_r3": ((-10,12,-8),3.0), "NAccR_r3": ((10,12,-8),3.0),
    "NAccL_r4": ((-10,12,-8),4.0), "NAccR_r4": ((10,12,-8),4.0),
    "NAccL_r5": ((-10,12,-8),5.0), "NAccR_r5": ((10,12,-8),5.0),
    "NAccL_r6": ((-10,12,-8),6.0), "NAccR_r6": ((10,12,-8),6.0),
    "NAccL_w":  ((-12,12,-8),4.0), "NAccR_w":  ((12,12,-8),4.0),
}
COUPLINGS = {
    "c_primary": ("VTA_r4", ["NAccL_r4","NAccR_r4"]),
    "c_rad2": ("VTA_r2", ["NAccL_r2","NAccR_r2"]),
    "c_rad3": ("VTA_r3", ["NAccL_r3","NAccR_r3"]),
    "c_rad5": ("VTA_r5", ["NAccL_r5","NAccR_r5"]),
    "c_rad6": ("VTA_r6", ["NAccL_r6","NAccR_r6"]),
    "c_vta_xp": ("VTA_xp", ["NAccL_r4","NAccR_r4"]),
    "c_vta_xm": ("VTA_xm", ["NAccL_r4","NAccR_r4"]),
    "c_vta_yp": ("VTA_yp", ["NAccL_r4","NAccR_r4"]),
    "c_vta_ym": ("VTA_ym", ["NAccL_r4","NAccR_r4"]),
    "c_vta_zp": ("VTA_zp", ["NAccL_r4","NAccR_r4"]),
    "c_vta_zm": ("VTA_zm", ["NAccL_r4","NAccR_r4"]),
    "c_nacc_w": ("VTA_r4", ["NAccL_w","NAccR_w"]),
}

def url_for(sub, ses):
    return ("%s%s%s/output/pipeline_analysis_nuisance/%s_%s/"
            "functional_preprocessed/_scan_rest_acq-1400/"
            "%s_%s_task-rest_acq-1400_bold_calc_tshift_resample_volreg_calc_maths.nii.gz"
            % (BUCKET, CP, sub, sub, ses, sub, ses))

def download(sub, ses, dest):
    last = -1; stall = 0
    for attempt in range(3):
        try:
            req = urllib.request.Request(url_for(sub, ses), headers={"User-Agent":"Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=45) as r, open(dest, "wb") as fh:
                while True:
                    b = r.read(1024*1024)
                    if not b: break
                    fh.write(b)
            if dest.exists() and dest.stat().st_size > 5_000_000:
                return True
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return "404"
            time.sleep(2*(attempt+1))
        except Exception:
            time.sleep(2*(attempt+1))
    return False

def roi_idx(img, mni, radius):
    v = np.round((np.linalg.inv(img.affine) @ np.array([mni[0],mni[1],mni[2],1.0]))[:3]).astype(int)
    zx,zy,zz = img.header.get_zooms()[:3]; sh = img.shape[:3]
    rx,ry,rz = (int(np.ceil(radius/z)) for z in (zx,zy,zz))
    xs=np.arange(max(0,v[0]-rx),min(sh[0],v[0]+rx+1))
    ys=np.arange(max(0,v[1]-ry),min(sh[1],v[1]+ry+1))
    zs=np.arange(max(0,v[2]-rz),min(sh[2],v[2]+rz+1))
    X,Y,Z=np.meshgrid(xs,ys,zs,indexing="ij")
    m=np.sqrt(((X-v[0])*zx)**2+((Y-v[1])*zy)**2+((Z-v[2])*zz)**2)<=radius
    return X[m],Y[m],Z[m]

def bandpass(x, lo=0.01, hi=0.08, tr=TR, order=3):
    nyq=0.5/tr
    b,a=butter(order,[lo/nyq,hi/nyq],btype="band")
    return filtfilt(b,a,x,axis=0)

def extract(path):
    img = nib.load(path)
    vol = np.asanyarray(img.dataobj, dtype=np.float32)
    n = vol.shape[3]
    if n < 30: return None, None
    ser = {}
    for name,(mni,rad) in ROIS.items():
        xi,yi,zi = roi_idx(img,mni,rad)
        s = vol[xi,yi,zi,:].mean(axis=0).astype(np.float64) if len(xi) else np.full(n,np.nan)
        if np.isnan(s).any():
            ser[name]=s
        else:
            x = bandpass(sp_detrend(s))
            ser[name]=(x-x.mean())/(x.std()+1e-12)
    out={}
    for name,(seed,tgts) in COUPLINGS.items():
        a=ser[seed]; rs=[]
        for t in tgts:
            b=ser[t]
            if np.isnan(a).any() or np.isnan(b).any(): continue
            rs.append(np.corrcoef(a,b)[0,1])
        out[name]=float(np.mean(rs)) if rs else np.nan
    return out, n

def main():
    plan = pd.read_csv(ROOT.parent / "data_derived" / "full_all_scans.csv")
    plan = plan[plan["session"].notna() & (plan["session"].astype(str).str.len()>0)]
    plan = plan[["subject","session"]].drop_duplicates().reset_index(drop=True)
    csv = OUT/"coupling_preproc.csv"
    done = set()
    prev = None
    if csv.exists():
        prev = pd.read_csv(csv); done = set(zip(prev.subject, prev.session))
    todo = [ (r.subject, r.session) for _,r in plan.iterrows() if (r.subject,r.session) not in done ]
    print("todo: %d  (done %d / %d)" % (len(todo), len(done), len(plan)), flush=True)
    if not todo: return

    ts_path = OUT/"roi_timeseries.npz"
    ts_store = dict(np.load(ts_path, allow_pickle=True)) if ts_path.exists() else {}
    new_rows = []
    t0 = time.time(); nfail = 0

    for i in range(0, len(todo), NWORK):
        batch = todo[i:i+NWORK]
        dests = [TMP/("b%d.nii.gz" % k) for k in range(len(batch))]
        with ThreadPoolExecutor(max_workers=NWORK) as ex:
            res = list(ex.map(lambda a: download(a[0], a[1], a[2]),
                              [(s, ss, d) for (s, ss), d in zip(batch, dests)]))
        for (sub,ses), d, rv in zip(batch, dests, res):
            rec = {"subject": sub, "session": ses}
            if rv == "404":
                rec["err"] = "HTTP 404"
            elif rv is not True:
                rec["err"] = "download failed"; nfail += 1
            else:
                try:
                    c, n = extract(d)
                    if c is None:
                        rec["err"] = "too few volumes"; nfail += 1
                    else:
                        rec.update(c); rec["n_vol"] = n
                except Exception as e:
                    rec["err"] = str(e)[:90]; nfail += 1
            new_rows.append(rec)
            try:
                if d.exists(): d.unlink()
            except Exception: pass
        allrows = (prev.to_dict("records") if prev is not None else []) + new_rows
        pd.DataFrame(allrows).to_csv(csv, index=False)
        el = time.time()-t0
        rate = len(new_rows)/el if el>0 else 0
        print("  %d/%d  %.0fs  %.1f/min  ~%.0f min left  fail=%d"
              % (len(new_rows), len(todo), el, rate*60,
                 (len(todo)-len(new_rows))/rate/60 if rate>0 else 0, nfail), flush=True)
    print("DONE new=%d fail=%d  %.0f min" % (len(new_rows), nfail, (time.time()-t0)/60), flush=True)

if __name__ == "__main__":
    main()

