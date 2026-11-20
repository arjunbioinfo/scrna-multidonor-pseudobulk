# scrna-multidonor-pseudobulk

[![CI](https://github.com/arjunbioinfo/scrna-multidonor-pseudobulk/actions/workflows/ci.yml/badge.svg)](https://github.com/arjunbioinfo/scrna-multidonor-pseudobulk/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

A reproducible, config-driven **Snakemake + Python** workflow for multi-donor single-cell RNA-seq **Harmony** batch integration, donor-level **pseudobulk PyDESeq2** differential expression, naive single-cell vs. pseudobulk false-positive benchmarking, and cell-type immune pathway enrichment.

Every push runs offline unit tests and executes the full Snakemake pipeline on real public 8-donor human PBMC data (**Kang et al., 2018**, IFN-beta stimulated vs. control, `GSE96583`) in GitHub Actions CI.

## Why Pseudobulk DE?

Treating individual cells from the same patient as independent replicates inflates false positives in case-control scRNA-seq studies ([Single-Cell Best Practices](https://www.sc-best-practices.org/conditions/differential_gene_expression.html)). This workflow aggregates raw counts per `(donor, condition, cell_type)`, models biological replicate variability with negative-binomial **PyDESeq2**, and directly compares the resulting DEGs against a naive single-cell Wilcoxon test.

## Pipeline Overview

1. **Download & Standardize (`download_data.py`)**: Fetches the Kang et al. (2018) 8-donor PBMC `.h5ad` from [Figshare (pertpy/scverse)](https://pertpy.scverse.org/en/1.0.2/_modules/pertpy/data/_datasets.html), standardizes `donor`, `condition`, and `cell_type` columns, preserves integer raw counts in `layers["counts"]`, and records dataset provenance in `data/raw/download_manifest.json`.
2. **QC & Harmony Integration (`qc_and_harmony.py`)**: Filters low-quality cells/genes, runs Scanpy HVG selection and PCA, integrates donors with `harmonypy`, computes pre- vs. post-Harmony UMAPs, and quantifies donor mixing (kNN Shannon entropy) and cell-type conservation (silhouette score).
3. **Pseudobulk PyDESeq2 & Benchmarking (`pseudobulk_deseq2.py`)**: Aggregates raw counts per `(donor, condition, cell_type)`, computes pseudobulk sample PCA, runs `PyDESeq2` (`stim` vs. `ctrl`) per cell type, benchmarks pseudobulk vs. naive single-cell Wilcoxon DEG counts on the top responsive population, and scores curated immune/interferon pathways.
4. **Visualization (`plot_summary.py`)**: Generates four 300 DPI summary figures in `figures/`.

## Quick Start

```bash
# 1. Create environment (Conda)
conda env create -f environment.yml
conda activate scrna-multidonor-pb

# ...or pip:
pip install -r requirements.txt

# 2. Run unit tests (offline)
make test

# 3. Fast end-to-end run on real Kang 2018 8-donor PBMC data (CI profile)
make ci

# 4. Full cohort run (all cells)
make run

# Or run via Docker
docker build -t scrna-multidonor-pb .
docker run --rm -v "$(pwd)/results:/app/results" -v "$(pwd)/figures:/app/figures" scrna-multidonor-pb
```

## Outputs

* `results/batch_metrics.tsv` — Donor kNN mixing entropy and silhouette scores before vs. after Harmony integration.
* `results/pseudobulk_pca.tsv` — Sample-level pseudobulk PCA coordinates across donors, conditions, and cell types.
* `results/pseudobulk_deseq2_all.tsv` — Full per-cell-type PyDESeq2 differential expression statistics (`baseMean`, `log2FoldChange`, `stat`, `pvalue`, `padj`).
* `results/celltype_deg_summary.tsv` — Up/down DEG counts per cell type (`padj < 0.05`, `|LFC| >= 1`).
* `results/pseudobulk_vs_wilcoxon.tsv` — DEG count and overlap comparison between donor-level PyDESeq2 and naive cell-level Wilcoxon.
* `results/pathway_enrichment.tsv` — Cell-type × immune pathway enrichment z-scores and BH-FDR.
* `figures/` — `harmony_integration_umap.png`, `pseudobulk_pca_and_degs.png`, `volcano_and_pb_vs_sc.png`, and `pathway_enrichment_heatmap.png`.

## Using Your Own Multi-Donor Cohort

Place your `.h5ad` file in `data/raw/kang_2018.h5ad` (or update `dataset.url` in `config/config.yaml`) with integer raw counts in `X` or `layers["counts"]` and donor, condition, and cell-type columns in `adata.obs`. Update `comparison.control` and `comparison.treated` in `config/config.yaml` and run `make run`.

## Repository Layout

* `Snakefile` — End-to-end Snakemake workflow
* `Makefile` — `setup`, `test`, `ci`, `run`, and `clean` shortcuts
* `config/config.yaml` / `config/ci.yaml` — Full-cohort and fast CI configurations
* `workflow/scripts/` — Modular Python scripts (`download_data.py`, `qc_and_harmony.py`, `pseudobulk_deseq2.py`, `plot_summary.py`)
* `notebooks/walkthrough.ipynb` — Interactive exploration notebook
* `tests/test_pipeline.py` — Offline 6-donor synthetic unit + smoke test suite
* `environment.yml` / `requirements.txt` — Conda and pip environments
* `Dockerfile` / `.dockerignore` — Container image definition
* `CITATION.cff` — Software citation metadata for GitHub and Zenodo
* `.github/workflows/ci.yml` — GitHub Actions CI (pytest + real Kang 2018 PBMC run)

## License

MIT — see [LICENSE](LICENSE).

## Author

**Mallikarjuna Thippana, PhD** — Computational Biologist, Karolinska Institutet  
[Portfolio](https://arjunbioinfo.github.io/) · [ORCID](https://orcid.org/0000-0001-8485-8626) · [Google Scholar](https://scholar.google.com/citations?user=NpLUj68AAAAJ) · [GitHub](https://github.com/arjunbioinfo)
