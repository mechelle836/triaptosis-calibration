# =====================================================================
# 30_virtual_knockout.py —— scGen 虚拟敲除 TAS 基因
# 训练 scGen VAE, 敲除 TAS 基因, 预测对内体/PI3P/细胞死亡基因的影响
# =====================================================================
import scgen
import anndata as ad
import pandas as pd
import numpy as np

ROOT = "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
RES = f"{ROOT}/results/18_virtual_knockout"
import os
os.makedirs(RES, exist_ok=True)

# ---- 读 h5ad ------------------------------------------------------------
adata = ad.read(f"{ROOT}/results/11_scrna/scrna_gse159115.h5ad")
print(f"读入: {adata.n_obs} cells x {adata.n_vars} genes")

# ---- scGen 训练 (VAE) ---------------------------------------------------
scgen.SCVI.setup_anndata(adata)
model = scgen.SCVI(adata)
print("scGen 训练中 (VAE, ~30-60 分钟)...")
model.train()

# ---- TAS 10 基因 --------------------------------------------------------
tas10 = ["ATP13A2", "ELMO2", "FYCO1", "KXD1", "MTM1", "NRBF2", "PIK3R4",
         "RAB9A", "SH3GL3", "WDR91"]

# ---- 虚拟敲除: 敲除每个 TAS 基因, 预测对内体/PI3P/细胞死亡基因的影响 ------
# 内体/PI3P/细胞死亡相关基因 (从 triaptosis 机制)
target_genes = ["PIK3C3", "PIK3R4", "ATG14", "NRBF2", "BECN1", "UVRAG",
                "MTM1", "KXD1", "WDR91", "FYCO1", "RAB9A", "SCARB2",
                "ATP13A2", "CLCN4", "SH3GL3", "ELMO2", "NCKAP1", "ACTN2",
                "ANXA8", "KEAP1", "NFE2L2"]

results = []
for ko_gene in tas10:
    if ko_gene not in adata.var_names:
        continue
    try:
        # 虚拟敲除: 预测敲除 ko_gene 后的表达变化
        pred = scgen.tools.knockout(model, adata, ko_gene)
        # 计算对 target_genes 的影响 (敲除后表达变化)
        for tg in target_genes:
            if tg in adata.var_names:
                idx = adata.var_names.get_loc(tg)
                delta = pred[:, idx].mean() - adata.X[:, idx].mean()
                results.append({"knockout": ko_gene, "target": tg,
                                "delta_expr": float(delta)})
        print(f"敲除 {ko_gene}: 完成")
    except Exception as e:
        print(f"敲除 {ko_gene}: 失败 {e}")

# ---- 保存结果 -----------------------------------------------------------
df = pd.DataFrame(results)
df.to_csv(f"{RES}/virtual_knockout_effects.csv", index=False)
print(f"虚拟敲除完成: {len(df)} 条影响记录")

# ---- 汇总: 每个 TAS 基因敲除对 triaptosis 机制基因的总影响 ----------------
summary = df.groupby("knockout")["delta_expr"].agg(["mean", "std", "count"]).reset_index()
summary.columns = ["knockout", "mean_delta", "std_delta", "n_targets"]
summary = summary.sort_values("mean_delta")
summary.to_csv(f"{RES}/virtual_knockout_summary.csv", index=False)
print("每个 TAS 基因敲除的总影响:")
print(summary)