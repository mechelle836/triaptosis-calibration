#!/usr/bin/env python3
"""
Paper2-Triaptosis: 通过 cBioPortal API 获取泛癌 21 个 triaptosis 基因表达 + 临床生存数据
机制: Triaptosis (Science 2024, VPS34/PI(3)P 内体依赖性氧化性细胞死亡)
数据源: TCGA PanCanAtlas 2018 (公开数据)
API 陷阱 (经验复用自第一篇 fetch_data.py):
  - ?studyId= 查询参数被忽略 -> 全量拉取后本地过滤
  - molecular-data/fetch: sampleListId 放 body
  - genes: 必须用路径参数 /api/genes/{symbol}
"""
import json, time, urllib.request, os, sys

API = "https://www.cbioportal.org/api"
OUT = "/Users/xiaolin/CodeBuddy/Claw/paper2_triaptosis/data"
os.makedirs(OUT, exist_ok=True)

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

# ---------- Step 1: 基因 Entrez ID ----------
print("=" * 60); print("Step 1: 基因 Entrez ID")
gene_map = {}
for g in GENES:
    d = get(f"{API}/genes/{g}")
    if d and "entrezGeneId" in d:
        gene_map[g] = d["entrezGeneId"]
    else:
        print(f"  {g}: NOT FOUND!")
print(f"  成功: {len(gene_map)}/{len(GENES)}")
with open(f"{OUT}/gene_entrez_map.json", "w") as f:
    json.dump(gene_map, f, indent=1)
entrez_ids = list(gene_map.values())
id2gene = {v: k for k, v in gene_map.items()}
if len(entrez_ids) < 15:
    print("!! 基因映射过少, 终止"); sys.exit(1)

# ---------- Step 2: 全量 profiles, 本地过滤 ----------
print("=" * 60); print("Step 2: 分子 profiles")
all_profiles = get(f"{API}/molecular-profiles?projection=SUMMARY")
by_study = {}
for p in all_profiles:
    by_study.setdefault(p.get("studyId", ""), []).append(p)

summary = {}
for cancer, study in CANCERS.items():
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
        print(f"  {cancer}: 无 RNA-seq profile"); continue
    profile_id = mrna[0]["molecularProfileId"]

    data = post(f"{API}/molecular-profiles/{profile_id}/molecular-data/fetch",
                {"entrezGeneIds": entrez_ids, "sampleListId": f"{study}_all"})
    if not isinstance(data, list) or len(data) == 0:
        print(f"  {cancer}: 表达数据为空"); continue

    mat = {}
    for rec in data:
        s, v, gid = rec.get("sampleId"), rec.get("value"), rec.get("entrezGeneId")
        if s is None or v is None or gid not in id2gene: continue
        mat.setdefault(s, {})[id2gene[gid]] = v

    n_samples = len(mat)
    n_genes = len({g for sv in mat.values() for g in sv})
    print(f"  {cancer}: {profile_id} | samples={n_samples} genes={n_genes}")
    with open(f"{OUT}/{cancer}_expr.json", "w") as f:
        json.dump(mat, f)
    summary[cancer] = {"profile": profile_id, "samples": n_samples, "genes": n_genes}

# ---------- Step 3: 临床生存数据 ----------
print("=" * 60); print("Step 3: 临床生存数据")
for cancer, study in CANCERS.items():
    if cancer not in summary: continue
    patients = get(f"{API}/studies/{study}/patients")
    if not patients:
        print(f"  {cancer}: 患者列表失败"); continue
    patient_ids = [p["patientId"] for p in patients]

    clin = post(f"{API}/studies/{study}/clinical-data/fetch",
                {"attributeIds": ["OS_STATUS", "OS_MONTHS"], "ids": patient_ids})
    if not isinstance(clin, list):
        print(f"  {cancer}: 临床数据失败"); continue

    cdict = {}
    for rec in clin:
        pid, attr, val = rec.get("patientId"), rec.get("clinicalAttributeId"), rec.get("value")
        if pid and attr:
            cdict.setdefault(pid, {})[attr] = val

    with open(f"{OUT}/{cancer}_clin.json", "w") as f:
        json.dump(cdict, f)

    n_os = sum(1 for v in cdict.values() if "OS_STATUS" in v and "OS_MONTHS" in v)
    n_dead = sum(1 for v in cdict.values()
                 if str(v.get("OS_STATUS", "")).upper() in ("1:DECEASED", "DECEASED", "1"))
    print(f"  {cancer}: patients={len(cdict)} OS完整={n_os} 死亡事件={n_dead}")

print("=" * 60)
print("汇总:")
print(json.dumps(summary, indent=1))
print(f"数据已保存至 {OUT}")
