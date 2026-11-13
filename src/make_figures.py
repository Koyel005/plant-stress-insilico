"""Figures from results/deg_results.csv and results/normalized_counts.csv. Run: python src/make_figures.py"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from pathlib import Path
Path("figures").mkdir(exist_ok=True)
r = pd.read_csv("results/deg_results.csv", index_col=0)
norm = pd.read_csv("results/normalized_counts.csv", index_col=0)
col = {"up": "#d62728", "down": "#1f77b4", "ns": "#bbbbbb"}
marks = {"AT5G52310": "RD29A", "AT5G52300": "RD29B", "AT3G50970": "LTI30", "AT1G20440": "COR47",
         "AT2G39800": "P5CS1", "AT5G25610": "RD22"}

# 1 volcano
fig, ax = plt.subplots(figsize=(7, 6))
ax.scatter(r.log2FoldChange, -np.log10(r.padj.clip(lower=1e-300)), c=r.regulation.map(col), s=5, alpha=.7)
for x in (-1, 1): ax.axvline(x, ls="--", c="k", lw=.7)
ax.axhline(-np.log10(0.05), ls="--", c="k", lw=.7)
for g, n in marks.items():
    off = {"COR47": (-34, 10), "P5CS1": (6, -12)}.get(n, (4, 4))
    ax.annotate(n, (r.loc[g, "log2FoldChange"], -np.log10(r.loc[g, "padj"])), fontsize=8, fontweight="bold",
                xytext=off, textcoords="offset points", arrowprops=dict(arrowstyle="-", lw=.5))
ax.set(xlabel="log2 fold change (Drought / Watered)", ylabel="-log10 adjusted p",
       title=f"Volcano: {(r.regulation=='up').sum()} up, {(r.regulation=='down').sum()} down")
plt.tight_layout(); plt.savefig("figures/fig1_volcano.png", dpi=200); plt.close()

# 2 MA plot
fig, ax = plt.subplots(figsize=(7, 5))
ax.scatter(np.log10(r.baseMean + 1), r.log2FoldChange, c=r.regulation.map(col), s=5, alpha=.7)
ax.axhline(0, c="k", lw=.7)
ax.set(xlabel="log10 mean normalised count", ylabel="log2 fold change", title="MA plot")
plt.tight_layout(); plt.savefig("figures/fig2_MA.png", dpi=200); plt.close()

# 3 PCA on log2 normalised counts (top 500 most variable genes)
lg = np.log2(norm + 1); top = lg.var(axis=1).sort_values(ascending=False).head(500).index
p = PCA(n_components=2).fit(lg.loc[top].T); xy = p.transform(lg.loc[top].T)
fig, ax = plt.subplots(figsize=(6, 5))
for i, s in enumerate(lg.columns):
    c = "#d62728" if "Drought" in s else "#1f77b4"
    ax.scatter(*xy[i], c=c, s=70); ax.annotate(s.replace("WT_", ""), xy[i], fontsize=8, xytext=(5, 5), textcoords="offset points")
ax.set(xlabel=f"PC1 ({p.explained_variance_ratio_[0]*100:.1f}%)", ylabel=f"PC2 ({p.explained_variance_ratio_[1]*100:.1f}%)",
       title="PCA of samples (top 500 variable genes)")
plt.tight_layout(); plt.savefig("figures/fig3_PCA.png", dpi=200); plt.close()
print("PC1 %.1f%% PC2 %.1f%%" % tuple(p.explained_variance_ratio_[:2] * 100))

# 4 heatmap of top 40 up + top 40 down (by adjusted p), z-scored
sel = list(r[r.regulation == "up"].sort_values("padj").head(40).index) + list(r[r.regulation == "down"].sort_values("padj").head(40).index)
z = lg.loc[sel]; z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1), axis=0)
fig, ax = plt.subplots(figsize=(6, 10))
im = ax.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
ax.set_xticks(range(6)); ax.set_xticklabels([c.replace("WT_", "") for c in z.columns], rotation=45, ha="right", fontsize=8)
ax.set_yticks([]); ax.axhline(39.5, c="k", lw=1)
ax.set_title("Top 40 up (top) and top 40 down (bottom)\nz-scored log2 normalised counts", fontsize=10)
plt.colorbar(im, shrink=.4); plt.tight_layout(); plt.savefig("figures/fig4_heatmap.png", dpi=200); plt.close()

# 5 known drought markers
m = r.loc[list(marks)]; m.index = list(marks.values())
fig, ax = plt.subplots(figsize=(6, 4))
ax.bar(m.index, m.log2FoldChange, color="#d62728")
for i, (n, row) in enumerate(m.iterrows()): ax.text(i, row.log2FoldChange + .05, f"{row.log2FoldChange:.2f}", ha="center", fontsize=8)
ax.set(ylabel="log2 fold change (Drought / Watered)", title="Known drought-responsive genes (sanity check)")
plt.xticks(rotation=30); plt.tight_layout(); plt.savefig("figures/fig5_known_markers.png", dpi=200); plt.close()
