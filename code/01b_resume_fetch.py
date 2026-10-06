#!/usr/bin/env python3
"""断点续传版: 跳过 data/ 已有的 {cancer}_expr.json 与 {cancer}_clin.json"""
import json, time, urllib.request, os, sys

API = "https://www.cbioportal.org/api"
OUT = "/Users/xiaolin/CodeBuddy/Claw/paper2_triaptosis/data"

GENES = ["PIK3C3","PIK3R4","MTM1","KEAP1","NFE2L2","SH3GL3","ELMO2","NCKAP1","ACTN2",
         "KXD1","WDR91","FYCO1","RAB9A","SCARB2","ANXA8","CLCN4","ATP13A2",
         "UVRAG","ATG14","NRBF2","BECN1"]

CANCERS = {
    "ACC":"acc_tcga_pan_can_atlas_2018","BLCA":"blca_tcga_pan_can_atlas_2018",
    "BRCA":"brca_tcga_pan_can_atlas_2018","CESC":"cesc_tcga_pan_can_atlas_2018",
    "CHOL":"chol_tcga_pan_can_atlas_2018","COADREAD":"coadread_tcga_pan_can_atlas_2018",
    "DLBC":"dlbc_tcga_pan_can_atlas_2018","ESCA":"esca_tcga_pan_can_atlas_2018",
    "GBM":"gbm_tcga_pan_can_atlas_2018","HNSC":"hnsc_tcga_pan_can_atlas_2018",
    "KICH":"kich_tcga_pan_can_atlas_2018","KIRC":"kirc_tcga_pan_can_atlas_2018",
    "KIRP":"kirp_tcga_pan_can_atlas_2018","LAML":"laml_tcga_pan_can_atlas_2018",
    "LGG":"lgg_tcga_pan_can_atlas_2018","LIHC":"lihc_tcga_pan_can_atlas_2018",
    "LUAD":"luad_tcga_pan_can_atlas_2018","LUSC":"lusc_tcga_pan_can_atlas_2018",
    "MESO":"meso_tcga_pan_can_atlas_2018","OV":"ov_tcga_pan_can_atlas_2018",
    "PAAD":"paad_tcga_pan_can_atlas_2018","PCPG":"pcpg_tcga_pan_can_atlas_2018",
    "PRAD":"prad_tcga_pan_can_atlas_2018","SARC":"sarc_tcga_pan_can_atlas_2018",
    "SKCM":"skcm_tcga_pan_can_atlas_2018","STAD":"stad_tcga_pan_can_atlas_2018",
    "TGCT":"tgct_tcga_pan_can_atlas_2018","THCA":"thca_tcga_pan_can_atlas_2018",
    "THYM":"thym_tcga_pan_can_atlas_2018","UCEC":"ucec_tcga_pan_can_atlas_2018",
    "UCS":"ucs_tcga_pan_can_atlas_2018","UVM":"uvm_tcga_pan_can_atlas_2018",
}

def get(url, retries=4):
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=90) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            print(f"  GET retry {i+1}: {e}"); time.sleep(2)
    return None

def post(url, body, retries=4):
    for i in range(retries):
        try:
            req = urllib.request.Request(url,
                data=json.dumps(body).encode(),
                headers={"Accept": "application/json", "Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=240) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            print(f"  POST retry {i+1}: {e}"); time.sleep(2)
    return None

# 基因映射(已缓存则读取)
map_path = f"{OUT}/gene_entrez_map.json"
if os.path.exists(map_path):
    gene_map = json.load(open(map_path))
else:
    gene_map = {}
    for g in GENES:
        d = get(f"{API}/genes/{g}")
        if d and "entrezGeneId" in d: gene_map[g] = d["entrezGeneId"]
    json.dump(gene_map, open(map_path, "w"), indent=1)
entrez_ids = list(gene_map.values())
id2gene = {v: k for k, v in gene_map.items()}
print(f"基因映射: {len(gene_map)}", flush=True)

# profiles 全量(仅首次需要,缓存)
prof_path = f"{OUT}/_profiles_cache.json"
if os.path.exists(prof_path):
    all_profiles = json.load(open(prof_path))
else:
    all_profiles = get(f"{API}/molecular-profiles?projection=SUMMARY")
    json.dump(all_profiles, open(prof_path, "w"))
by_study = {}
for p in all_profiles:
    by_study.setdefault(p.get("studyId", ""), []).append(p)

summary = {}
todo_expr = [c for c in CANCERS if not os.path.exists(f"{OUT}/{c}_expr.json")]
todo_clin = [c for c in CANCERS if not os.path.exists(f"{OUT}/{c}_clin.json")]
print(f"待抓表达: {len(todo_expr)} | 待抓临床: {len(todo_clin)}", flush=True)

for cancer in todo_expr:
    study = CANCERS[cancer]
    profs = by_study.get(study, [])
    mrna = [p for p in profs
            if p["molecularAlterationType"] == "MRNA_EXPRESSION"
            and p["datatype"] == "CONTINUOUS"
            and p["molecularProfileId"] == f"{study}_rna_seq_v2_mrna"]
    if not mrna:
        mrna = [p for p in profs
                if p["molecularAlterationType"] == "MRNA_EXPRESSION"
                and p["datatype"] == "CONTINUOUS" and "rna_seq" in p["molecularProfileId"]]
    if not mrna:
        print(f"  {cancer}: 无 profile"); continue
    profile_id = mrna[0]["molecularProfileId"]
    data = post(f"{API}/molecular-profiles/{profile_id}/molecular-data/fetch",
                {"entrezGeneIds": entrez_ids, "sampleListId": f"{study}_all"})
    if not isinstance(data, list) or len(data) == 0:
        print(f"  {cancer}: 空"); continue
    mat = {}
    for rec in data:
        s, v, gid = rec.get("sampleId"), rec.get("value"), rec.get("entrezGeneId")
        if s is None or v is None or gid not in id2gene: continue
        mat.setdefault(s, {})[id2gene[gid]] = v
    json.dump(mat, open(f"{OUT}/{cancer}_expr.json", "w"))
    print(f"  expr {cancer}: {len(mat)} samples", flush=True)

for cancer in todo_clin:
    study = CANCERS[cancer]
    patients = get(f"{API}/studies/{study}/patients")
    if not patients:
        print(f"  {cancer}: patients失败"); continue
    clin = post(f"{API}/studies/{study}/clinical-data/fetch",
                {"attributeIds": ["OS_STATUS", "OS_MONTHS"], "ids": [p["patientId"] for p in patients]})
    if not isinstance(clin, list):
        print(f"  {cancer}: 临床失败"); continue
    cdict = {}
    for rec in clin:
        pid, attr, val = rec.get("patientId"), rec.get("clinicalAttributeId"), rec.get("value")
        if pid and attr: cdict.setdefault(pid, {})[attr] = val
    json.dump(cdict, open(f"{OUT}/{cancer}_clin.json", "w"))
    nd = sum(1 for v in cdict.values() if str(v.get("OS_STATUS","")).upper() in ("1:DECEASED","DECEASED","1"))
    print(f"  clin {cancer}: {len(cdict)} pts, {nd} dead", flush=True)

print("DONE", flush=True)
