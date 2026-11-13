import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
steps = ["Public data\nGEO GSE119382\n(WT roots, 3 drought\nvs 3 watered)",
         "Raw counts\n(third-party\nreanalysis,\n32,833 genes)",
         "Filter\ntotal count >= 10\n(22,399 kept)",
         "PyDESeq2\n~condition\nWald test + BH",
         "Call DEGs\npadj < 0.05\n|log2FC| >= 1",
         "STRING network\ntop 500 up genes\nscore >= 0.4",
         "NetworkX\ndegree and\nbetweenness"]
fig, ax = plt.subplots(figsize=(14, 2.6)); ax.axis("off"); ax.set_xlim(0, 14); ax.set_ylim(0, 2.6)
for i, t in enumerate(steps):
    x = 0.15 + i * 1.98
    ax.add_patch(FancyBboxPatch((x, 0.35), 1.7, 1.9, boxstyle="round,pad=0.03", fc="#e8f1fb" if i < 5 else "#fdf0e0", ec="#345"))
    ax.text(x + 0.85, 1.3, t, ha="center", va="center", fontsize=8.5)
    if i < len(steps) - 1: ax.annotate("", xy=(x + 1.95, 1.3), xytext=(x + 1.73, 1.3), arrowprops=dict(arrowstyle="->"))
ax.text(0.15, 0.1, "Blue = re-run and verifiable from this repo.  Orange = STRING website export (data/string_interactions.tsv), then NetworkX.", fontsize=8)
plt.savefig("figures/fig0_workflow.png", dpi=200, bbox_inches="tight")
