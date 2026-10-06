#!/usr/bin/env python3
"""Paper2-Triaptosis: 外部验证数据补抓 (2026-09-02, v2 修复)

v2 修复:
  - /api/genes/{symbol} 返回字段是 entrezGeneId (不是 entrezId)
  - clinical-data/fetch 必须带 ids (患者列表), 否则返回不完整 (01 号脚本既有经验)
  - 先抓该研究全部临床记录, 本地按白名单过滤 + 打印实际属性名 (避免猜属性名)

任务 A: CPTAC-3 ccRCC (rcc_cptac_gdc) TAS-21 基因 mRNA 表达
任务 B: TCGA-KIRC (kirc_tcga_pan_can_atlas_2018) 患者级临床扩展 (AGE/SEX/分期)
"""
import json, os, time, urllib.request, csv

API = "https://www.cbioportal.org/api"
PROJ = "/Users/xiaolin/CodeBuddy/Claw/paper2_triaptosis"
OUT = os.path.join(PROJ, "data", "external_validation")
os.makedirs(OUT, exist_ok=True)

TAS21 = ["PIK3C3","PIK3R4","MTM1","KEAP1","NFE2L2","SH3GL3","ELMO2","NCKAP1","ACTN2",
         "KXD1","WDR91","FYCO1","RAB9A","SCARB2","ANXA8","CLCN4","ATP13A2",
         "UVRAG","ATG14","NRBF2","BECN1"]

WANT = ["AGE","SEX","AJCC_PATHOLOGIC_TUMOR_STAGE","AJCC_TUMOR_STAGE","TUMOR_STAGE",
        "AJCC_PATHOLOGIC_M_STAGE","AJCC_PATHOLOGIC_N_STAGE","GRADE","FUHRMAN_GRADE",
        "OS_STATUS","OS_MONTHS"]


def get(url, retries=4, timeout=120):
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


def post(url, body, retries=4, timeout=240):
    for i in range(retries):
        try:
            req = urllib.request.Request(
                url, data=json.dumps(body).encode(),
                headers={"Accept": "application/json", "Content-Type": "application/json",
                         "User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            print(f"  POST fail {i+1}: {e}")
            time.sleep(2)
    return None


def write_csv(path, rows, fieldnames=None):
    if not rows:
        print("  EMPTY ->", path)
        return
    if fieldnames is None:
        keys = []
        for r in rows:
            for k in r:
                if k not in keys:
                    keys.append(k)
        fieldnames = keys
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {path} n={len(rows)}")


# ================================================================
print("=" * 70)
print("Task A: CPTAC-3 (rcc_cptac_gdc) TAS-21 expression")
emap = {}
for g in TAS21:
    d = get(f"{API}/genes/{g}")
    if isinstance(d, dict) and d.get("entrezGeneId"):
        emap[g] = int(d["entrezGeneId"])
    else:
        print(f"  WARN: no entrez for {g}")
    time.sleep(0.3)
print(f"  entrez map: {len(emap)}/21")
ID2G = {v: k for k, v in emap.items()}
with open(os.path.join(OUT, "tas21_entrez_map.json"), "w") as f:
    json.dump(emap, f, indent=1)

expr = post(f"{API}/molecular-profiles/rcc_cptac_gdc_mrna_seq_tpm/molecular-data/fetch",
            {"entrezGeneIds": sorted(set(emap.values())),
             "sampleListId": "rcc_cptac_gdc_all"})
if not isinstance(expr, list):
    expr = post(f"{API}/molecular-profiles/rcc_cptac_gdc_rna_seq_mrna/molecular-data/fetch",
                {"entrezGeneIds": sorted(set(emap.values())),
                 "sampleListId": "rcc_cptac_gdc_all"})
    if isinstance(expr, list):
        print("  (used rna_seq_mrna fallback)")
mat, genes_seen = {}, set()
if isinstance(expr, list):
    for rec in expr:
        sid, val, gid = rec.get("sampleId"), rec.get("value"), rec.get("entrezGeneId")
        g = ID2G.get(int(gid) if gid is not None else -1)
        if not sid or val is None or not g:
            continue
        mat.setdefault(sid, {"sampleId": sid})
        mat[sid][g] = val
        genes_seen.add(g)
print(f"  samples={len(mat)} genes={sorted(genes_seen)}")
write_csv(os.path.join(OUT, "cptac_tas_expr21.csv"), list(mat.values()),
          fieldnames=["sampleId"] + sorted(emap))
print()

# ================================================================
print("=" * 70)
print("Task B: TCGA-KIRC clinical extension (kirc_tcga_pan_can_atlas_2018)")
STUDY = "kirc_tcga_pan_can_atlas_2018"
patients = get(f"{API}/studies/{STUDY}/patients", timeout=180)
pids = [p["patientId"] for p in patients] if isinstance(patients, list) else []
print(f"  patients={len(pids)}")

clin = post(f"{API}/studies/{STUDY}/clinical-data/fetch",
            {"attributeIds": WANT, "ids": pids})
if not isinstance(clin, list) or not clin:
    print("  POST fetch (with ids) empty -> GET all paged")
    clin, page, page_size = [], 0, 10000
    while True:
        chunk = get(f"{API}/studies/{STUDY}/clinical-data"
                    f"?pageNumber={page}&pageSize={page_size}&projection=DETAILED")
        if not isinstance(chunk, list) or not chunk:
            break
        clin.extend(chunk)
        if len(chunk) < page_size:
            break
        page += 1
print("  clinical records:", len(clin) if isinstance(clin, list) else None)

if isinstance(clin, list) and clin:
    seen_attrs = sorted({r.get("clinicalAttributeId") for r in clin
                         if r.get("clinicalAttributeId")})
    print("  attrs present:", seen_attrs)
    wide = {}
    for rec in clin:
        aid = rec.get("clinicalAttributeId")
        pid = rec.get("patientId")
        if not pid:
            continue
        wide.setdefault(pid, {"patientId": pid})
        if aid in WANT:
            wide[pid][aid] = rec.get("value")
    rows = list(wide.values())
    keep_cols = ["patientId"] + [a for a in WANT
                                 if any(r.get(a) not in (None, "") for r in rows)]
    write_csv(os.path.join(OUT, "KIRC_clin_ext.csv"), rows, fieldnames=keep_cols)
    for a in keep_cols[1:]:
        n = sum(1 for r in rows if r.get(a) not in (None, "", "NA"))
        print(f"  attr {a}: {n}/{len(rows)}")
print("DONE", OUT)
