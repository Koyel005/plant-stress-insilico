# Drought response in *Arabidopsis thaliana* roots: PyDESeq2 + STRING network (GEO: GSE119382)

Computational reanalysis of public RNA-seq data (wild-type roots, 3 drought vs 3 watered).

- **Result:** 908 genes up, 578 down (padj < 0.05, |log2FC| >= 1); public reanalysis reports 904 / 576
- **Known drought genes recovered:** RD29A, RD29B, LTI30, COR47, RD22
- **Network:** top 500 upregulated genes -> STRING (web export included) -> 774 interactions among 333 genes; P5CSA has the highest betweenness; only 9.0% of edges are experimental
- **Full write-up:** [REPORT.md](REPORT.md)

## Important note on data
Raw counts come from a **public third-party reanalysis** (github.com/YuntongMeng/arabidopsis-drought-rnaseq) of GEO GSE119382.
I did not generate the sequencing data or the count table. My work starts from the count table.

## Run
```bash
pip install -r requirements.txt
python src/deseq2_analysis.py
python src/make_figures.py
python src/make_string_input.py   # top-500 list for string-db.org (web export -> data/string_interactions.tsv)
python src/network_hubs.py --edges data/string_interactions.tsv
python src/network_figures.py
```
