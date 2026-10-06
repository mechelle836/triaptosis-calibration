#!/usr/bin/env python3
"""Fetch VHL, PBRM1, SETD2, BAP1 mutations for TCGA-KIRC PanCancer Atlas from cBioPortal.

Writes one row per sequenced patient. A gene is 1 only for a non-silent mutation;
patients absent from the sequenced list are not written, so they are not treated as wild type.
"""
import csv
import json
import os
import urllib.request

API = "https://www.cbioportal.org/api"
STUDY = "kirc_tcga_pan_can_atlas_2018"
PROFILE = f"{STUDY}_mutations"
SEQ_LIST = f"{STUDY}_sequenced"
GENES = {"VHL": 7428, "PBRM1": 55193, "SETD2": 29072, "BAP1": 8314}
SILENT = {"Silent", "Intron", "3'UTR", "5'UTR", "3'Flank", "5'Flank", "IGR", "RNA"}

P2 = "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
OUT_DIR = os.path.join(P2, "data", "external_validation")


def get(url):
    with urllib.request.urlopen(url, timeout=120) as r:
        return json.load(r)


def post(url, body):
    req = urllib.request.Request(
        url, data=json.dumps(body).encode(), headers={"Content-Type": "application/json"}
    )
    with urllib.request.urlopen(req, timeout=180) as r:
        return json.load(r)


sample_ids = get(f"{API}/sample-lists/{SEQ_LIST}/sample-ids")
patients = sorted({s[:12] for s in sample_ids})

muts = post(
    f"{API}/molecular-profiles/{PROFILE}/mutations/fetch?projection=DETAILED",
    {"sampleListId": SEQ_LIST, "entrezGeneIds": list(GENES.values())},
)

by_entrez = {v: k for k, v in GENES.items()}
hit = {g: set() for g in GENES}
raw_rows = []
for m in muts:
    gene = by_entrez.get(m["entrezGeneId"])
    mtype = m.get("mutationType", "")
    pid = m["patientId"][:12]
    raw_rows.append([pid, m["sampleId"], gene, mtype, m.get("proteinChange", "")])
    if gene and mtype not in SILENT:
        hit[gene].add(pid)

with open(os.path.join(OUT_DIR, "kirc_tcga_driver_mutations.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["patientId"] + list(GENES))
    for p in patients:
        w.writerow([p] + [int(p in hit[g]) for g in GENES])

with open(os.path.join(OUT_DIR, "kirc_tcga_driver_mutations_raw.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["patientId", "sampleId", "gene", "mutationType", "proteinChange"])
    w.writerows(raw_rows)

print("sequenced samples", len(sample_ids), "patients", len(patients))
print("mutation records", len(muts))
for g in GENES:
    print(g, len(hit[g]))
