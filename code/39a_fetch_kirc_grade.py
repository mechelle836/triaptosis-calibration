#!/usr/bin/env python3
"""Fetch histologic grade (sample attribute) and pathologic T stage for TCGA-KIRC PanCanAtlas.

Writes one row per patient. When a patient has several samples, the primary tumour (-01) is used.
"""
import csv
import json
import os
import urllib.request

API = "https://www.cbioportal.org/api"
STUDY = "kirc_tcga_pan_can_atlas_2018"
P2 = "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
OUT = os.path.join(P2, "data", "external_validation", "kirc_tcga_grade.csv")


def get(url):
    with urllib.request.urlopen(url, timeout=120) as r:
        return json.load(r)


samples = get(f"{API}/studies/{STUDY}/clinical-data?clinicalDataType=SAMPLE&attributeId=GRADE&projection=SUMMARY")
patients = get(f"{API}/studies/{STUDY}/clinical-data?clinicalDataType=PATIENT&attributeId=PATH_T_STAGE&projection=SUMMARY")

grade = {}
for rec in sorted(samples, key=lambda r: r["sampleId"]):
    pid = rec["patientId"]
    if pid not in grade or rec["sampleId"].endswith("-01"):
        grade[pid] = rec["value"]
tstage = {rec["patientId"]: rec["value"] for rec in patients}

pids = sorted(set(grade) | set(tstage))
with open(OUT, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["patientId", "GRADE", "PATH_T_STAGE"])
    for p in pids:
        w.writerow([p, grade.get(p, ""), tstage.get(p, "")])

vals = {}
for v in grade.values():
    vals[v] = vals.get(v, 0) + 1
print("patients", len(pids), "grade values", vals)
