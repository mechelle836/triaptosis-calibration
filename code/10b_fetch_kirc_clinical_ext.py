#!/usr/bin/env python3
"""Paper2-Triaptosis: KIRC 临床扩展补抓 (v2 — PanCanAtlas 2018 开放属性无 AGE/OS/分期)

路径 1: sample 级 clinical-data (kirc_tcga_pan_can_atlas_2018)
路径 2: legacy study kirc_tcga (TCGA 2013, 属性更全)
输出: data/external_validation/KIRC_clin_ext.csv (patientId 宽表)
"""
import json, os, urllib.request, csv, time

API = "https://www.cbioportal.org/api"
PROJ = "/Users/xiaolin/CodeBuddy/Claw/paper2_triaptosis"
OUT = os.path.join(PROJ, "data", "external_validation")
os.makedirs(OUT, exist_ok=True)

WANT_PAT = ["AGE", "PATIENT_AGE", "SEX", "GENDER", "OS_STATUS", "OS_MONTHS",
            "AJCC_PATHOLOGIC_TUMOR_STAGE", "AJCC_TUMOR_STAGE", "TUMOR_STAGE",
            "AJCC_PATHOLOGIC_M_STAGE", "AJCC_PATHOLOGIC_N_STAGE", "GRADE"]
WANT_SAMP = ["AJCC_PATHOLOGIC_TUMOR_STAGE", "TUMOR_STAGE", "SAMPLE_TYPE"]


def get(url, retries=4, timeout=180):
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={
                "Accept": "application/json", "User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            print(f"  GET fail {i+1}: {e}")
            time.sleep(2)
    return None


def fetch_all(study, dtype):
    rows, page, size = [], 0, 10000
    while True:
        c = get(f"{API}/studies/{study}/clinical-data"
                f"?clinicalDataType={dtype}&pageNumber={page}&pageSize={size}"
                f"&projection=DETAILED")
        if not isinstance(c, list) or not c:
            break
        rows.extend(c)
        if len(c) < size:
            break
        page += 1
    return rows


def widefy(recs, want):
    wide = {}
    for rec in recs:
        pid, sid = rec.get("patientId"), rec.get("sampleId")
        aid, val = rec.get("clinicalAttributeId"), rec.get("value")
        key = pid if sid is None else sid
        wide.setdefault(key, {"patientId": pid, "sampleId": sid})
        if aid in want:
            wide[key][aid] = val
    return list(wide.values())


for study in ["kirc_tcga", "kirc_tcga_pan_can_atlas_2018"]:
    print("=" * 70)
    print(f"Study: {study}")
    pat = fetch_all(study, "PATIENT")
    sam = fetch_all(study, "SAMPLE")
    print(f"  patient recs={len(pat)} sample recs={len(sam)}")
    if pat:
        print("  patient attrs:", sorted({r["clinicalAttributeId"] for r in pat}))
    if sam:
        print("  sample attrs:", sorted({r["clinicalAttributeId"] for r in sam}))

    patw = widefy(pat, WANT_PAT) if pat else []
    samw = widefy(sam, WANT_SAMP) if sam else []
    print(f"  patient rows={len(patw)} sample rows={len(samw)}")

    if patw:
        cols = ["patientId"] + [a for a in WANT_PAT
                                if any(r.get(a) not in (None, "") for r in patw)]
        write = os.path.join(OUT, f"{study}_clin_patient.csv")
        with open(write, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
            w.writeheader(); w.writerows(patw)
        print(f"  wrote {write}")
        for a in cols[1:]:
            n = sum(1 for r in patw if r.get(a) not in (None, "", "NA"))
            print(f"    attr {a}: {n}/{len(patw)}")
    if samw:
        cols = ["patientId", "sampleId"] + [a for a in WANT_SAMP
                                            if any(r.get(a) not in (None, "") for r in samw)]
        write = os.path.join(OUT, f"{study}_clin_sample.csv")
        with open(write, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
            w.writeheader(); w.writerows(samw)
        print(f"  wrote {write}")
print("DONE")
