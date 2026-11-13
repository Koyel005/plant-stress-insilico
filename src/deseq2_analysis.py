"""
Step 1: Differential expression (drought vs watered) in wild-type Arabidopsis roots.
Dataset : GEO GSE119382. Raw counts come from a PUBLIC THIRD-PARTY REANALYSIS
          (github.com/YuntongMeng/arabidopsis-drought-rnaseq), not from my own read alignment.
Method  : PyDESeq2 (Python implementation of DESeq2), design ~condition.
Run     : python src/deseq2_analysis.py
"""
import numpy as np, pandas as pd
from pathlib import Path
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

ALPHA, LFC, MIN_TOTAL = 0.05, 1.0, 10
OUT = Path("results"); OUT.mkdir(exist_ok=True)

counts = pd.read_csv("data/counts_all_11samples_raw.tsv", sep="\t", index_col=0)
wt = [c for c in counts.columns if c.startswith("WT_")]          # 3 drought + 3 watered, wild type only
counts = counts[wt]
meta = pd.DataFrame({"condition": ["Drought" if "Drought" in c else "Watered" for c in wt]}, index=wt)
print(meta.condition.value_counts().to_string())

keep = counts.sum(axis=1) >= MIN_TOTAL                            # remove near-zero genes
print(f"{keep.sum()} of {len(keep)} genes kept (total count >= {MIN_TOTAL})")
dds = DeseqDataSet(counts=counts.loc[keep].T.astype(int), metadata=meta, design="~condition", quiet=True)
dds.deseq2()
stats = DeseqStats(dds, contrast=["condition", "Drought", "Watered"], quiet=True)
stats.summary()

res = stats.results_df.dropna(subset=["padj"]).copy()
res["regulation"] = "ns"
res.loc[(res.padj < ALPHA) & (res.log2FoldChange >= LFC), "regulation"] = "up"
res.loc[(res.padj < ALPHA) & (res.log2FoldChange <= -LFC), "regulation"] = "down"
res.index.name = "gene_id"
res.reset_index().to_csv(OUT / "deg_results.csv", index=False)

norm = pd.DataFrame(dds.layers["normed_counts"], index=dds.obs_names, columns=dds.var_names).T
norm.to_csv(OUT / "normalized_counts.csv")
print("genes tested (non-NA padj):", len(res))
print("RESULT:", (res.regulation == "up").sum(), "up,", (res.regulation == "down").sum(), "down")
