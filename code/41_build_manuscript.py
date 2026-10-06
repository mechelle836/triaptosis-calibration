#!/usr/bin/env python3
"""Build manuscript_frontiers_oncology.md from manuscript_keys.md.

Replaces [@key; @key] citations with numbers in order of first appearance, writes the
reference list in that order, and fills Table S1 from results/01_screen.
Fails if a key has no reference entry or a reference is never cited.
"""
import csv
import os
import re
import statistics

P2 = "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
SRC = os.path.join(P2, "manuscript_keys.md")
DST = os.path.join(P2, "manuscript_frontiers_oncology.md")

REFS = {
    "jonasch2014": "Jonasch E, Gao J, Rathmell WK. Renal cell carcinoma. *BMJ*. 2014;349:g4797.",
    "choueiri2017": "Choueiri TK, Motzer RJ. Systemic therapy for metastatic renal-cell carcinoma. *N Engl J Med*. 2017;376:354–366.",
    "rini2005": "Rini BI, Small EJ. Biology and clinical development of vascular endothelial growth factor-targeted therapy in renal cell carcinoma. *J Clin Oncol*. 2005;23:1028–1043.",
    "tang2019": "Tang D, Kang R, Berghe TV, Vandenabeele P, Kroemer G. The molecular machinery of regulated cell death. *Cell Res*. 2019;29:347–364.",
    "swamynathan2024": "Swamynathan MM, Kuang S, Watrud KE, et al. Dietary pro-oxidant therapy by a vitamin K precursor targets PI 3-kinase VPS34 function. *Science*. 2024;386:eadk9167. PMID: 39446948.",
    "tang2025": "Tang D, Kang R, Kroemer G. Triaptosis: an endosome-dependent cell death modality. *Cell Res*. 2025;35:237–238. PMID: 39638924.",
    "xie2025": "Xie J, Zhang M, Qi M. Integrating machine learning algorithms to construct a triaptosis-related prognostic model in melanoma. *Cancer Manag Res*. 2025;17:1127–1141. PMID: 40535575.",
    "liu2025": "Liu X, Zhuang Z, Cheng J, et al. From single-cell and bulk transcriptomic integration to functional verification: triaptosis-associated lncRNA signature predicts survival and guides therapy in hepatocellular carcinoma. *Pharmaceuticals (Basel)*. 2025;18:1691. doi: 10.3390/ph18111691. PMID: 41304936.",
    "liu2026": "Liu J, Guan X, Dou M, et al. A novel triaptosis-related prognostic signature to assess the clinical value in HER2-low breast cancer: evidence from clinical cohorts and experimental validation. *Apoptosis*. 2026;31:89. doi: 10.1007/s10495-026-02313-2. PMID: 41793564.",
    "dai2026": "Dai H, Lu R, Zhang M, et al. Integrating single-cell multi-omics and machine learning to reveal triaptosis heterogeneity in clear cell renal cell carcinoma. *Hum Genomics*. 2026;20. doi: 10.1186/s40246-026-00982-3. PMID: 42098870.",
    "venet2011": "Venet D, Dumont JE, Detours V. Most random gene expression signatures are significantly associated with breast cancer outcome. *PLoS Comput Biol*. 2011;7:e1002240.",
    "hoadley2018": "Hoadley KA, Yau C, Hinoue T, et al. Cell-of-origin patterns dominate the molecular classification of 10,000 tumors from 33 types of cancer. *Cell*. 2018;173:291–304.e6. PMID: 29625048.",
    "gao2013": "Gao J, Aksoy BA, Dogrusoz U, et al. Integrative analysis of complex cancer genomics and clinical profiles using the cBioPortal. *Sci Signal*. 2013;6:pl1.",
    "sato2013": "Sato Y, Yoshizato T, Shiraishi Y, et al. Integrated molecular analysis of clear-cell renal cell carcinoma. *Nat Genet*. 2013;45:860–867. PMID: 23797736.",
    "clark2019": "Clark DJ, Dhanasekaran SM, Petralia F, et al. Integrated proteogenomic characterization of clear cell renal cell carcinoma. *Cell*. 2019;179:964–983.e31. PMID: 31675502.",
    "braun2020": "Braun DA, Hou Y, Bakouny Z, et al. Interplay of somatic alterations and immune infiltration modulates response to PD-1 blockade in advanced clear cell renal cell carcinoma. *Nat Med*. 2020;26:909–918. PMID: 32472114.",
    "mcdermott2018": "McDermott DF, Huseni MA, Atkins MB, et al. Clinical activity and molecular correlates of response to atezolizumab alone or in combination with bevacizumab versus sunitinib in renal cell carcinoma. *Nat Med*. 2018;24:749–757. PMID: 29867230.",
    "goldman2020": "Goldman MJ, Craft B, Hastie M, et al. Visualizing and interpreting cancer genomics data via the Xena platform. *Nat Biotechnol*. 2020;38:675–678.",
    "friedman2010": "Friedman J, Hastie T, Tibshirani R. Regularization paths for generalized linear models via coordinate descent. *J Stat Softw*. 2010;33:1–22.",
    "wu2021": "Wu T, Hu E, Xu S, et al. clusterProfiler 4.0: a universal enrichment tool for interpreting omics data. *Innovation (Camb)*. 2021;2:100141.",
    "therneau2000": "Therneau TM, Grambsch PM. *Modeling Survival Data: Extending the Cox Model*. New York: Springer; 2000.",
    "heagerty2000": "Heagerty PJ, Lumley T, Pepe MS. Time-dependent ROC curves for censored survival data and a diagnostic marker. *Biometrics*. 2000;56:337–344.",
    "uno2011": "Uno H, Cai T, Pencina MJ, D’Agostino RB, Wei LJ. On the C-statistics for evaluating overall adequacy of risk prediction procedures with censored survival data. *Stat Med*. 2011;30:1105–1117.",
    "gerds2013": "Gerds TA, Kattan MW, Schumacher M, Yu C. Estimating a time-dependent concordance index for survival prediction models with covariates dependent censoring. *Stat Med*. 2013;32:2173–2184.",
    "brooks2014": "Brooks SA, Brannon AR, Parker JS, et al. ClearCode34: a prognostic risk predictor for localized clear cell renal cell carcinoma. *Eur Urol*. 2014;66:77–84.",
    "tibshirani2002": "Tibshirani R, Hastie T, Narasimhan B, Chu G. Diagnosis of multiple cancer types by shrunken centroids of gene expression. *Proc Natl Acad Sci U S A*. 2002;99:6567–6572.",
    "develasco2017": "de Velasco G, Culhane AC, Fay AP, et al. Molecular subtypes improve prognostic value of International Metastatic Renal Cell Carcinoma Database Consortium prognostic model. *Oncologist*. 2017;22:286–292.",
    "rini2015": "Rini B, Goddard A, Knezevic D, et al. A 16-gene assay to predict recurrence after surgery in localised renal cell carcinoma: development and validation studies. *Lancet Oncol*. 2015;16:676–685.",
    "rini2018": "Rini BI, Escudier B, Martini JF, et al. Validation of the 16-gene recurrence score in patients with locoregional, high-risk renal cell carcinoma from a phase III trial of adjuvant sunitinib. *Clin Cancer Res*. 2018;24:4407–4415.",
    "hanzelmann2013": "Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for microarray and RNA-seq data. *BMC Bioinformatics*. 2013;14:7.",
    "charoentong2017": "Charoentong P, Finotello F, Angelova M, et al. Pan-cancer immunogenomic analyses reveal genotype-immunophenotype relationships and predictors of response to checkpoint blockade. *Cell Rep*. 2017;18:248–262. PMID: 28052254.",
    "yang2013": "Yang W, Soares J, Greninger P, et al. Genomics of Drug Sensitivity in Cancer (GDSC): a resource for therapeutic biomarker discovery in cancer cells. *Nucleic Acids Res*. 2013;41:D955–D961.",
    "johnson2007": "Johnson WE, Li C, Rabinovic A. Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics*. 2007;8:118–127.",
    "zhang2021": "Zhang Y, Narayanan SP, Mannan R, et al. Single-cell analyses of renal cell cancers reveal insights into tumor microenvironment, cell of origin, and therapy response. *Proc Natl Acad Sci U S A*. 2021;118:e2103240118. PMID: 34099557.",
    "hao2024": "Hao Y, Stuart T, Kowalski MH, et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis. *Nat Biotechnol*. 2024;42:293–304.",
    "aibar2017": "Aibar S, González-Blas CB, Moerman T, et al. SCENIC: single-cell regulatory network inference and clustering. *Nat Methods*. 2017;14:1083–1086.",
    "jin2021": "Jin S, Guerrero-Juarez CF, Zhang L, et al. Inference and analysis of cell-cell communication using CellChat. *Nat Commun*. 2021;12:1088.",
    "davidson2023": "Davidson G, Helleux A, Vano YA, et al. Mesenchymal-like tumor cells and myofibroblastic cancer-associated fibroblasts are associated with progression and immunotherapy response of clear cell renal cell carcinoma. *Cancer Res*. 2023;83:2952–2969. PMID: 37335139.",
    "keith2011": "Keith B, Johnson RS, Simon MC. HIF1α and HIF2α: sibling rivalry in hypoxic tumour growth and progression. *Nat Rev Cancer*. 2011;12:9–22.",
    "sauerbrei2014": "Sauerbrei W, Abrahamowicz M, Altman DG, le Cessie S, Carpenter J; STRATOS initiative. STRengthening analytical thinking for observational studies: the STRATOS initiative. *Stat Med*. 2014;33:5413–5432.",
}

s = open(SRC, encoding="utf-8").read()
body, refs_tail = s.split("@@REFERENCES@@")
order = []
for m in re.finditer(r"\[(@[^\]]+)\]", body):
    for k in re.findall(r"@([A-Za-z0-9]+)", m.group(1)):
        if k not in order:
            order.append(k)
missing = [k for k in order if k not in REFS]
unused = [k for k in REFS if k not in order]
if missing or unused:
    raise SystemExit(f"missing refs: {missing}; unused refs: {unused}")
num = {k: i + 1 for i, k in enumerate(order)}


def fmt(nums):
    nums = sorted(nums)
    out, i = [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        out.append(f"{nums[i]}–{nums[j]}" if j - i >= 2 else ", ".join(str(n) for n in nums[i:j + 1]))
        i = j + 1
    return "(" + ", ".join(out) + ")"


body = re.sub(r"\[(@[^\]]+)\]",
              lambda m: fmt([num[k] for k in re.findall(r"@([A-Za-z0-9]+)", m.group(1))]), body)
reflist = "\n".join(f"{num[k]}. {REFS[k]}" for k in order)
out = body + reflist + refs_tail

# Table S1 rows
R = os.path.join(P2, "results/01_screen")
rows = {r["cancer"]: r for r in csv.DictReader(open(os.path.join(R, "screen_null_nsig_by_cancer.csv")))}
byc = {r["cancer"]: r for r in csv.DictReader(open(os.path.join(R, "screen_null_by_cancer.csv")))}
orig = {r["cancer"]: r for r in csv.DictReader(open(os.path.join(R, "pancancer_screen_matrix.csv")))}
sets = {}
for r in csv.DictReader(open(os.path.join(R, "screen_null_sets.csv"))):
    sets.setdefault(r["cancer"], []).append(int(r["n_sig"]))


def q(v, p):
    v = sorted(v)
    h = (len(v) - 1) * p
    lo = int(h)
    return v[lo] + (v[min(lo + 1, len(v) - 1)] - v[lo]) * (h - lo)


lines = []
for c in sorted(rows, key=lambda c: (float(rows[c]["frac_ge_obs"]), -int(rows[c]["obs_nsig"]))):
    v = sets[c]
    med = statistics.median(v)
    lines.append("| {} | {} | {} | {:g} ({:g}–{:g}) | {:.3f} | {:.1f} | {:.1f} | {} |".format(
        c, orig[c]["n_sig_cox"], rows[c]["obs_nsig"], med, q(v, 0.05), q(v, 0.95),
        float(rows[c]["frac_ge_obs"]), 100 * float(byc[c]["frac"]), 100 * float(byc[c]["pass_soft_frac"]),
        "yes" if orig[c]["pass_soft"] == "TRUE" else "no"))
out = out.replace("@@TABLE_S1@@", "\n".join(lines))
if "@@" in out or re.search(r"\[@", out):
    raise SystemExit("unresolved placeholder")
open(DST, "w", encoding="utf-8").write(out)

abstract = out[out.index("## Abstract"):out.index("**Keywords**")]
abstract = re.sub(r"\*\*[A-Za-z]+\*\*:", "", abstract.replace("## Abstract", ""))
print("references:", len(order))
print("abstract words:", len(re.findall(r"\S+", abstract)))
