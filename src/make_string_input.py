"""Writes the top 500 upregulated genes (ranked by adjusted p-value) as a list to paste into string-db.org."""
import pandas as pd
deg = pd.read_csv("results/deg_results.csv")
up = deg[deg.regulation == "up"].sort_values("padj").head(500)["gene_id"]
up.to_csv("results/up_gene_ids_top500_for_STRING.txt", index=False, header=False)
print(len(up), "genes written")
