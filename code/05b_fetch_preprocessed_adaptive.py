#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
自适应限速规避下载器
- 快的时候连续下载
- 单文件耗时超过 HARD_CAP 秒 或 速率低于 MIN_RATE 时，暂停 COOL 秒
- 可断点续跑
"""
import os, sys, time, socket
from pathlib import Path
import numpy as np, pandas as pd, nibabel as nib
from scipy.signal import butter, filtfilt, detrend as sp_detrend
import urllib.request

socket.setdefaulttimeout(25)
ROOT = Path(__file__).resolve().parent.parent
OUT  = ROOT / "data_derived" / "preproc_work"; OUT.mkdir(parents=True, exist_ok=True)
TMP  = ROOT / "data_derived" / "preproc_work" / "_tmp"; TMP.mkdir(parents=True, exist_ok=True)
BUCKET = "https://fcp-indi.s3.amazonaws.com/"
CP     = "data/Projects/RocklandSample/Outputs/C-PAC/"
TR     = 1.4
HARD_CAP = 75.0     # 单文件最长下载秒数
MIN_RATE = 1.8      # MB/s 低于此值则冷却
COOL     = 45       # 冷却秒数

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
    "c_primary": ("VTA_r4",["NAccL_r4","NAccR_r4"]), "c_rad2": ("VTA_r2",["NAccL_r2","NAccR_r2"]),
    "c_rad3": ("VTA_r3",["NAccL_r3","NAccR_r3"]), "c_rad5": ("VTA_r5",["NAccL_r5","NAccR_r5"]),
    "c_rad6": ("VTA_r6",["NAccL_r6","NAccR_r6"]),
    "c_vta_xp": ("VTA_xp",["NAccL_r4","NAccR_r4"]), "c_vta_xm": ("VTA_xm",["NAccL_r4","NAccR_r4"]),
    "c_vta_yp": ("VTA_yp",["NAccL_r4","NAccR_r4"]), "c_vta_ym": ("VTA_ym",["NAccL_r4","NAccR_r4"]),
    "c_vta_zp": ("VTA_zp",["NAccL_r4","NAccR_r4"]), "c_vta_zm": ("VTA_zm",["NAccL_r4","NAccR_r4"]),
    "c_nacc_w": ("VTA_r4",["NAccL_w","NAccR_w"]),
}

def url_for(s, ss):
    return ("%s%s%s/output/pipeline_analysis_nuisance/%s_%s/functional_preprocessed/"
            "_scan_rest_acq-1400/%s_%s_task-rest_acq-1400_bold_calc_tshift_resample_volreg_calc_maths.nii.gz"
            % (BUCKET, CP, s, s, ss, s, ss))

def download(s, ss, dest, cap=HARD_CAP):
    t0 = time.time(); got = 0
    try:
        req = urllib.request.Request(url_for(s, ss), headers={"User-Agent":"Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=25) as r, open(dest, "wb") as fh:
            while True:
                if time.time()-t0 > cap:
                    return ("slow", got)
                b = r.read(1024*256)
                if not b: break
                fh.write(b); got += len(b)
        return ("ok", got)
    except urllib.error.HTTPError as e:
        return ("404" if e.code == 404 else "err", got)
    except Exception:
        return ("err", got)

def roi_idx(img, mni, radius):
    v = np.round((np.linalg.inv(img.affine) @ np.array([mni[0],mni[1],mni[2],1.0]))[:3]).astype(int)
    zx,zy,zz = img.header.get_zooms()[:3]; sh = img.shape[:3]
    rx,ry,rz = (int(np.ceil(radius/z)) for z in (zx,zy,zz))
    xs=np.arange(max(0,v[0]-rx),min(sh[0],v[0]+rx+1)); ys=np.arange(max(0,v[1]-ry),min(sh[1],v[1]+ry+1))
    zs=np.arange(max(0,v[2]-rz),min(sh[2],v[2]+rz+1))
    X,Y,Z=np.meshgrid(xs,ys,zs,indexing="ij")
    m=np.sqrt(((X-v[0])*zx)**2+((Y-v[1])*zy)**2+((Z-v[2])*zz)**2)<=radius
    return X[m],Y[m],Z[m]

def extract(path):
    img = nib.load(path)
    vol = np.asanyarray(img.dataobj, dtype=np.float32)
    n = vol.shape[3]
    if n < 30: return None, None
    ser = {}
    for name,(mni,rad) in ROIS.items():
        xi,yi,zi = roi_idx(img,mni,rad)
        s = vol[xi,yi,zi,:].mean(axis=0).astype(np.float64) if len(xi) else np.full(n,np.nan)
        if np.isnan(s).any(): ser[name]=s
        else:
            b,a = butter(3, [0.01/(0.5/TR), 0.08/(0.5/TR)], btype="band")
            x = filtfilt(b,a,sp_detrend(s))
            ser[name]=(x-x.mean())/(x.std()+1e-12)
    out={}
    for name,(seed,tgts) in COUPLINGS.items():
        a=ser[seed]; rs=[]
        for t in tgts:
            bb=ser[t]
            if np.isnan(a).any() or np.isnan(bb).any(): continue
            rs.append(np.corrcoef(a,bb)[0,1])
        out[name]=float(np.mean(rs)) if rs else np.nan
    return out, n

def main():
    plan = pd.read_csv(ROOT / "data_derived" / "full_all_scans.csv")
    plan = plan[plan["session"].notna() & (plan["session"].astype(str).str.len()>0)]
    plan = plan[["subject","session"]].drop_duplicates().reset_index(drop=True)
    csv = OUT/"coupling_preproc.csv"
    prev = pd.read_csv(csv) if csv.exists() else pd.DataFrame(columns=["subject","session"])
    done = set(zip(prev.subject, prev.session))
    todo = [(r.subject, r.session) for _,r in plan.iterrows() if (r.subject,r.session) not in done]
    print("todo %d (done %d / %d)" % (len(todo), len(done), len(plan)), flush=True)
    rows = prev.to_dict("records")
    t_start = time.time(); nfail = 0; ncold = 0
    for i,(s,ss) in enumerate(todo, 1):
        dest = TMP/"cur.nii.gz"
        if dest.exists():
            try: dest.unlink()
            except Exception: pass
        t0 = time.time()
        st, got = download(s, ss, dest)
        dt = max(time.time()-t0, 1e-6)
        rec = {"subject": s, "session": ss}
        if st == "404": rec["err"] = "HTTP 404"
        elif st != "ok": rec["err"] = "download failed"; nfail += 1
        else:
            try:
                c, n = extract(dest)
                if c is None: rec["err"] = "too few volumes"; nfail += 1
                else: rec.update(c); rec["n_vol"] = n
            except Exception as e:
                rec["err"] = str(e)[:80]; nfail += 1
        rows.append(rec)
        try:
            if dest.exists(): dest.unlink()
        except Exception: pass
        rate = (got/1048576)/dt
        cooled = ""
        if st != "404" and rate < MIN_RATE:
            ncold += 1; time.sleep(COOL); cooled = "  [cool %ds]" % COOL
        pd.DataFrame(rows).to_csv(csv, index=False)
        el = (time.time()-t_start)/60
        print("  %3d/%d  %s %s  %5.1f MB %5.2f MB/s  elapsed %.1f min  fail=%d%s"
              % (i, len(todo), s, ss, got/1048576, rate, el, nfail, cooled), flush=True)
    print("DONE todo=%d fail=%d cold=%d  %.1f min" % (len(todo), nfail, ncold, (time.time()-t_start)/60), flush=True)

if __name__ == "__main__":
    main()
