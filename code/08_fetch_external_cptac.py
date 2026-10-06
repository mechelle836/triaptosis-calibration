#!/usr/bin/env python3
"""Paper2: 抓取外部验证队列 CPTAC (rcc_cptac_gdc, n=354) 21 个 triaptosis 基因表达 + OS 临床
独立于 TCGA 的蛋白质基因组学队列 (GDC CPTAC, 公开)
"""
import json, time, urllib.request, os

API = "https://www.cbioportal.org/api"
OUT = "/Users/xiaolin/CodeBuddy/Claw/paper2_triaptosis/data/external"
STUDY = "rcc_cptac_gdc"
PROFILE = f"{STUDY}_rna_seq_mrna"

def get(url, retries=3):
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=90) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            print(f"  GET retry {i+1}: {e}"); time.sleep(2)
    return None

def post(url, body, retries=3):
    for i in range(retries):
        try:
            req = urllib.request.Request(url, data=json.dumps(body).encode(),
                headers={"Accept": "application/json", "Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=240) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            print(f"  POST retry {i+1}: {e}"); time.sleep(2)
    return None

gene_map = json.load(open("/Users/xiaolin/CodeBuddy/Claw/paper2_triaptosis/data/gene_entrez_map.json"))
entrez_ids = list(gene_map.values())
id2gene = {v: k for k, v in gene_map.items()}

print("CPTAC 表达...", flush=True)
data = post(f"{API}/molecular-profiles/{PROFILE}/molecular-data/fetch",
            {"entrezGeneIds": entrez_ids, "sampleListId": f"{STUDY}_all"})
mat = {}
for rec in (data or []):
    s, v, gid = rec.get("sampleId"), rec.get("value"), rec.get("entrezGeneId")
    if s is None or v is None or gid not in id2gene: continue
    mat.setdefault(s, {})[id2gene[gid]] = v
print(f"  表达: {len(mat)} 样本", flush=True)
json.dump(mat, open(f"{OUT}/CPTAC_expr.json", "w"))

print("CPTAC 临床 (OS)...", flush=True)
patients = get(f"{API}/studies/{STUDY}/patients")
print(f"  patients: {len(patients) if patients else 0}", flush=True)
if patients:
    clin = post(f"{API}/studies/{STUDY}/clinical-data/fetch",
                {"attributeIds": ["OS_STATUS", "OS_MONTHS", "VITAL_STATUS"],
                 "ids": [p["patientId"] for p in patients]})
    cdict = {}
    for rec in (clin or []):
        pid, attr, val = rec.get("patientId"), rec.get("clinicalAttributeId"), rec.get("value")
        if pid and attr: cdict.setdefault(pid, {})[attr] = val
    json.dump(cdict, open(f"{OUT}/CPTAC_clin.json", "w"))
    n_os = sum(1 for v in cdict.values() if v.get("OS_MONTHS") and v.get("OS_STATUS"))
    n_dead = sum(1 for v in cdict.values() if str(v.get("OS_STATUS","")).upper() in ("1:DECEASED","DECEASED","1"))
    print(f"  临床: {len(cdict)} pts, OS完整={n_os}, 死亡={n_dead}", flush=True)
print("DONE", flush=True)
