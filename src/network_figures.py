"""Network figures + evidence audit from data/string_interactions.tsv and results/ppi_hubs.csv"""
import pandas as pd, numpy as np, networkx as nx, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
d = pd.read_csv("data/string_interactions.tsv", sep="\t"); d.columns = [c.lstrip("#") for c in d.columns]
d["pair"] = [tuple(sorted(p)) for p in zip(d.node1, d.node2)]; d = d.drop_duplicates("pair"); hubs = pd.read_csv("results/ppi_hubs.csv")
n = len(d)
ev = {"Text-mining": d.automated_textmining > 0, "Co-expression": d.coexpression > 0,
      "Database annotated": d.database_annotated > 0, "Experiments": d.experimentally_determined_interaction > 0,
      "Homology": d.homology > 0, "Gene neighbourhood": d.neighborhood_on_chromosome > 0}
ev = {k: round(v.sum() / n * 100, 1) for k, v in ev.items()}
print("unique edges", n, "| score min/median/max", d.combined_score.min(), d.combined_score.median(), d.combined_score.max())
print("evidence % of edges:", ev)
only_text = ((d.automated_textmining > 0) & (d.coexpression == 0) & (d.experimentally_determined_interaction == 0) & (d.database_annotated == 0)).mean() * 100
print("edges with text-mining as only listed evidence: %.1f%%" % only_text)

G = nx.Graph(); [G.add_edge(a, b) for a, b in zip(d.node1, d.node2)]
big = G.subgraph(max(nx.connected_components(G), key=len))
pos = nx.spring_layout(big, seed=7, k=0.35)
top = set(hubs.head(10).gene) | {"P5CSA"}
fig, ax = plt.subplots(figsize=(10, 10))
nx.draw_networkx_edges(big, pos, alpha=.15, width=.5, ax=ax)
sz = [20 + 12 * big.degree(x) for x in big]; col = ["#d62728" if x in top else "#7aa6d6" for x in big]
nx.draw_networkx_nodes(big, pos, node_size=sz, node_color=col, linewidths=.3, edgecolors="k", ax=ax)
nx.draw_networkx_labels(big, pos, labels={x: x for x in big if x in top}, font_size=8, font_weight="bold", ax=ax)
ax.set_title(f"STRING network of top 500 upregulated genes\nlargest component: {big.number_of_nodes()} genes, {big.number_of_edges()} edges (red = top hubs, size = degree)")
ax.axis("off"); plt.tight_layout(); plt.savefig("figures/fig6_network.png", dpi=200); plt.close()

fig, axs = plt.subplots(1, 2, figsize=(11, 4))
deg = [x for _, x in G.degree()]; axs[0].hist(deg, bins=range(1, max(deg) + 2), color="#7aa6d6", edgecolor="k")
axs[0].set(xlabel="degree", ylabel="number of genes", title="Degree distribution (most genes have few links)")
axs[1].barh(list(ev)[::-1], list(ev.values())[::-1], color="#d6a27a"); axs[1].set(xlabel="% of edges with this evidence", title="Evidence behind the 774 edges")
for i, v in enumerate(list(ev.values())[::-1]): axs[1].text(v + 1, i, f"{v}%", va="center", fontsize=8)
plt.tight_layout(); plt.savefig("figures/fig7_degree_evidence.png", dpi=200); plt.close()
