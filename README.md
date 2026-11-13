# Transcriptomic reanalysis of drought response in *Arabidopsis thaliana* roots (GEO: GSE119382)

Differential expression (PyDESeq2) and protein-interaction network analysis (STRING, NetworkX) of public RNA-seq data from wild-type roots, 3 drought vs 3 watered samples.

## Key results
- 908 genes upregulated and 578 downregulated (adjusted p < 0.05, |log2FC| >= 1); a public reanalysis reports 904 / 576
- Canonical drought genes recovered: RD29A, RD29B, LTI30, COR47, RD22
- STRING v12.5 network of the top 500 upregulated genes: 774 interactions among 333 genes; P5CSA has the highest betweenness; 9.0% of edges have experimental evidence

Full report: [REPORT.md](REPORT.md) (also [REPORT.pdf](REPORT.pdf))

## Data
Raw counts are from a public third-party reanalysis (github.com/YuntongMeng/arabidopsis-drought-rnaseq) of GEO GSE119382. Sequencing and count generation were not performed in this project.

## Run
```bash
pip install -r requirements.txt
python src/deseq2_analysis.py
python src/make_figures.py
python src/make_string_input.py    # gene list for string-db.org; export saved as data/string_interactions.tsv
python src/network_hubs.py --edges data/string_interactions.tsv
python src/network_figures.py
```
