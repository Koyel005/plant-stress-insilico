"""
Step 4: hub ranking of the STRING network of the top 500 upregulated genes (NetworkX).
Workflow used: paste results/up_gene_ids_top500_for_STRING.txt into string-db.org (Multiple proteins,
Arabidopsis thaliana, confidence 0.4), export the network as TSV, then:
    python src/network_hubs.py --edges data/string_interactions.tsv
"""
import argparse, pandas as pd, networkx as nx
from pathlib import Path
p = argparse.ArgumentParser(); p.add_argument("--edges", required=True); p.add_argument("--outdir", default="results")
a = p.parse_args(); Path(a.outdir).mkdir(exist_ok=True)

d = pd.read_csv(a.edges, sep="\t")
d.columns = [c.lstrip("#") for c in d.columns]
n1, n2 = ("node1", "node2") if "node1" in d.columns else (d.columns[0], d.columns[1])
d = d.drop_duplicates([n1, n2]); d.to_csv(f"{a.outdir}/ppi_edges.csv", index=False)

G = nx.Graph()
for _, e in d.iterrows(): G.add_edge(e[n1], e[n2])
bet = nx.betweenness_centrality(G)
hubs = pd.DataFrame({"gene": list(G.nodes), "degree": [G.degree(n) for n in G.nodes],
                     "betweenness": [bet[n] for n in G.nodes]}).sort_values(["degree", "betweenness"], ascending=False)
hubs.to_csv(f"{a.outdir}/ppi_hubs.csv", index=False)
big = max(nx.connected_components(G), key=len)
print(f"Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges, largest component {len(big)}")
print(hubs.head(10).to_string(index=False))
print("Highest betweenness:", hubs.sort_values("betweenness", ascending=False).head(3).to_string(index=False))
