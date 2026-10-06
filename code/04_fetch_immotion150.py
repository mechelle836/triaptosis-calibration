#!/usr/bin/env python3
"""Paper2: 抓取 IMmotion150 (rcc_iatlas_immotion150_2018) 21 个 triaptosis 基因表达 + RESPONSE 临床
数据源: cBioPortal 公开研究 (n=263 肾癌, Atezolizumab anti-PD-L1)
"""
import json, time, urllib.request, os

API = "https://www.cbioportal.org/api"
OUT = "/Users/xiaolin/CodeBuddy/Claw/paper2_triaptosis/data"
STUDY = "rcc_iatlas_immotion150_2018"
PROFILE = f"{STUDY}_rna_seq_mrna"

def get(url, retries=3):
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            print(f"  GET retry {i+1}: {e}"); time.sleep(2)
    return None

def post(url, body, retries=3):
    for i in range(retries):
        try:
            req = urllib.request.Request(url, data=json.dumps(body).encode(),
                headers={"Accept": "application/json", "Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=180) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            print(f"  POST retry {i+1}: {e}"); time.sleep(2)
    return None

gene_map = json.load(open(f"{OUT}/gene_entrez_map.json"))
entrez_ids = list(gene_map.values())
id2gene = {v: k for k, v in gene_map.items()}

print("IMmotion150 表达...", flush=True)
data = post(f"{API}/molecular-profiles/{PROFILE}/molecular-data/fetch",
            {"entrezGeneIds": entrez_ids, "sampleListId": f"{STUDY}_all"})
mat = {}
for rec in (data or []):
    s, v, gid = rec.get("sampleId"), rec.get("value"), rec.get("entrezGeneId")
    if s is None or v is None or gid not in id2gene: continue
    mat.setdefault(s, {})[id2gene[gid]] = v
print(f"  表达: {len(mat)} 样本", flush=True)
json.dump(mat, open(f"{OUT}/IMmotion150_expr.json", "w"))

print("IMmotion150 RESPONSE/临床...", flush=True)
samples = get(f"{API}/studies/{STUDY}/samples")
print(f"  samples: {len(samples) if samples else 0}", flush=True)
cdict = {}
for s in (samples or []):
    sid = s["sampleId"]
    d = get(f"{API}/studies/{STUDY}/samples/{sid}/clinical-data")
    if not d: continue
    rec = {r["clinicalAttributeId"]: r.get("value") for r in d
           if r["clinicalAttributeId"] in
           ("RESPONSE","RESPONDER","PFS_MONTHS","PFS_STATUS","CLINICAL_BENEFIT",
            "ICI_RX","ICI_PATHWAY","CLINICAL_STAGE","SAMPLE_COUNT")}
    if rec: cdict[sid] = rec
print(f"  临床: {len(cdict)}", flush=True)
resp = {}
for v in cdict.values():
    r = v.get("RESPONDER") or v.get("RESPONSE")
    if r: resp[r] = resp.get(r, 0) + 1
print("  RESPONDER 分布:", json.dumps(resp), flush=True)
json.dump(cdict, open(f"{OUT}/IMmotion150_clin.json", "w"))
print("DONE", flush=True)
