# Random-gene calibration of pan-cancer cell-death gene-set screens

Analysis code for the study *"Calibrating pan-cancer gene-set screening against
random gene sets localizes triaptosis prognostic signal to clear cell renal cell
carcinoma"*.

## Citation
Chu X, Chen C, Hu J, Xu L, Zhang X, Zhu X. (2026).
Putuo Hospital, Shanghai University of Traditional Chinese Medicine, Shanghai, China.

## What this repository provides
A reusable protocol for calibrating any mechanism-defined gene set against random
gene sets, and the results of applying it to a 21-gene triaptosis set across 32
TCGA cancer types.

**Headline methodological result**: an uncalibrated pan-cancer screen ranks
clear cell renal cell carcinoma (ccRCC) first for 53% of random 21-gene sets.
Only the per-gene prognostic yield discriminates true signal: ccRCC is the sole
cancer among 20 hard-pass cancers where triaptosis genes exceed the random
distribution (18/21 Cox-significant genes; no random set above 17; empirical
P = 0/200).

## Repository contents

| Path | Description |
|---|---|
| `code/01–02b` | Data acquisition and pan-cancer feasibility screening |
| `code/38–38b` | **Random-set calibration of the screen (200 sets)** |
| `code/07_triaptosis_null.R` | **Signature-level null distribution (195 sets)** |
| `code/34_pcd_benchmark.R` | **Ferroptosis/apoptosis/disulfidptosis comparison** |
| `code/36_enet_path.R` | Elastic-net TAS construction (α=0.5, 10-fold CV) |
| `code/07_model_validation.R` | Out-of-fold IPCW validation |
| `code/09–12, 33` | External validation (E-MTAB-1980, CPTAC-3, CheckMate) |
| `code/13, 32, 39` | Stage/grade/driver-mutation adjustment |
| `code/06, 22, 05` | Immune profiling and IMmotion150 ICI analysis |
| `code/19–30` | Single-cell and spatial transcriptomics (Supplementary) |
| `results/01_screen/` | Calibration outputs, per-cancer null distributions |
| `results/05_null/` | Signature-level null results and PCD benchmark |

## Reproducing the calibration

The calibration requires only public data (cBioPortal Pan-Cancer Atlas and UCSC
Xena) and no wet-lab input.

```r
# 1. Pan-cancer screen with the 21-gene triaptosis set
source("code/02_pancancer_screen.R")

# 2. Repeat the identical frozen screen with 200 random 21-gene sets
source("code/38_screen_null.R")       # seed 2026; outputs results/01_screen/

# 3. Signature-level null: 195 random sets through nested CV
source("code/07_triaptosis_null.R")   # 10 seeds, same outer folds
```

Expected outputs (all under `results/01_screen/`):
- `screen_null_summary.csv` — per-cancer random-set distribution
- `screen_null_reproduction.csv` — soft-pass reproduction rate
- `screen_null_source_agreement.csv` — cBioPortal vs Xena concordance

The same pipeline applies to any other mechanism-defined gene set by replacing the
gene list at the top of `02_pancancer_screen.R`.

## Software environment
- R 4.6.1 — survival 3.8.6, glmnet 5.0, survivalROC 1.0.3.1, timeROC 0.4.1,
  pec 2025.6.24, riskRegression 2026.3.11, GSVA 2.6.6, clusterProfiler 4.20.0,
  sva 3.60.0, pamr 1.57, rms 8.1.1, Seurat 5.5.1, AUCell 1.34.0, CellChat 2.2.0.9001
- Python 3.9.6 — data retrieval from public repositories
- Random seeds fixed throughout (screen calibration seed 2026; bootstrap seed 20260923)

## Data availability
All data are public: TCGA Pan-Cancer Atlas (cBioPortal), UCSC Xena PanCanAtlas
matrix, E-MTAB-1980 (ArrayExpress), CPTAC-3, IMmotion150 (GEO), CheckMate
009/010/025 (Braun et al. supplementary tables), GSE159115 (single cell),
GSE210041 (Visium).

## Licence
MIT (see LICENSE)
