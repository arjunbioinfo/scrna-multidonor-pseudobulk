#!/usr/bin/env python3
"""Donor-level pseudobulk aggregation, PyDESeq2 DE, single-cell Wilcoxon benchmark, and pathway enrichment."""

from __future__ import annotations

import argparse
import inspect
import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy import sparse
from scipy.stats import norm
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# Curated immune & interferon response gene sets for cell-type pathway enrichment
IMMUNE_GENE_SETS: dict[str, list[str]] = {
    "Type_I_Interferon_Response": [
        "ISG15", "IFIT1", "IFIT2", "IFIT3", "IFITM1", "IFITM2", "IFITM3",
        "IFI6", "IFI27", "IFI35", "IFI44", "IFI44L", "MX1", "MX2",
        "OAS1", "OAS2", "OAS3", "OASL", "RSAD2", "BST2", "STAT1", "STAT2", "IRF7", "IRF9",
    ],
    "Chemokine_Cytokine_Signaling": [
        "CXCL10", "CXCL11", "CCL2", "CCL7", "CCL8", "CCL3", "CCL4",
        "IL1B", "IL6", "TNF", "TNFSF10", "TNFSF13B", "IL1RN", "CXCL8",
    ],
    "Antigen_Presentation_MHC": [
        "HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9",
        "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "CIITA",
    ],
    "NFkB_Inflammatory_Response": [
        "NFKBIA", "NFKB1", "NFKB2", "RELB", "TNFAIP3", "SOD2", "BCL2A1",
        "ICAM1", "TRAF1", "JUNB", "FOS", "IER3", "PTGS2", "CD83",
    ],
    "Cytotoxic_Effector_Program": [
        "NKG7", "GNLY", "GZMA", "GZMB", "GZMH", "GZMK", "PRF1",
        "CST7", "CCL5", "KLRD1", "KLRB1", "CTSW", "FGFBP2",
    ],
}


def build_pseudobulk(
    adata: ad.AnnData,
    min_cells: int = 10,
    min_counts: int = 1000,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Sum raw counts by (donor, condition, cell_type) into pseudobulk sample matrix."""
    raw_mat = adata.layers["counts"] if "counts" in adata.layers else adata.X
    genes = adata.var_names.to_numpy()
    obs = adata.obs[["donor", "condition", "cell_type"]].copy().reset_index(drop=True)

    counts_rows: list[np.ndarray] = []
    meta_rows: list[dict[str, object]] = []

    for (donor, cond, ct), grp in obs.groupby(["donor", "condition", "cell_type"], observed=True):
        idx = grp.index.to_numpy()
        n_cells = len(idx)
        if n_cells < min_cells:
            continue

        sub = raw_mat[idx]
        vec = np.asarray(sub.sum(axis=0)).ravel() if sparse.issparse(sub) else sub.sum(axis=0)
        vec = np.rint(np.clip(vec, 0, None)).astype(np.int64)
        total_counts = int(vec.sum())
        if total_counts < min_counts:
            continue

        sample_id = f"{donor}__{cond}__{ct}".replace(" ", "_")
        counts_rows.append(vec)
        meta_rows.append(
            {
                "sample_id": sample_id,
                "donor": str(donor),
                "condition": str(cond),
                "cell_type": str(ct),
                "n_cells": n_cells,
                "total_counts": total_counts,
            }
        )

    if not meta_rows:
        raise ValueError("No pseudobulk samples passed min_cells / min_counts thresholds.")

    meta_df = pd.DataFrame(meta_rows).set_index("sample_id")
    counts_df = pd.DataFrame(np.vstack(counts_rows), index=meta_df.index, columns=genes)
    return counts_df, meta_df


def compute_pseudobulk_pca(counts_df: pd.DataFrame, meta_df: pd.DataFrame) -> pd.DataFrame:
    """Compute log1p-CPM PCA across all pseudobulk samples."""
    lib_sizes = counts_df.sum(axis=1).replace(0, 1)
    cpm = counts_df.div(lib_sizes, axis=0) * 1e6
    log_cpm = np.log1p(cpm)
    # Keep top 2000 most variable genes across pseudobulks
    gene_var = log_cpm.var(axis=0)
    top_genes = gene_var.nlargest(min(2000, log_cpm.shape[1])).index
    X = log_cpm[top_genes].to_numpy()
    X = X - X.mean(axis=0, keepdims=True)

    n_comps = min(5, X.shape[0] - 1, X.shape[1] - 1)
    pca = PCA(n_components=max(2, n_comps), random_state=42)
    coords = pca.fit_transform(X)

    out = meta_df.reset_index().copy()
    out["PC1"] = coords[:, 0]
    out["PC2"] = coords[:, 1]
    out["pc1_var_pct"] = round(float(pca.explained_variance_ratio_[0] * 100.0), 2)
    out["pc2_var_pct"] = round(float(pca.explained_variance_ratio_[1] * 100.0), 2)
    return out


def make_deseq_dataset(counts_sub: pd.DataFrame, meta_sub: pd.DataFrame) -> DeseqDataSet:
    """Instantiate DeseqDataSet compatibly across PyDESeq2 0.4.x and 0.5.x APIs."""
    sig = inspect.signature(DeseqDataSet.__init__)
    kwargs: dict[str, object] = {
        "counts": counts_sub,
        "metadata": meta_sub,
        "refit_cooks": True,
        "n_cpus": 1,
        "quiet": True,
    }
    if "design" in sig.parameters:
        kwargs["design"] = "~condition"
    else:
        kwargs["design_factors"] = "condition"
    return DeseqDataSet(**kwargs)


def run_pydeseq2_per_celltype(
    counts_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    cond_ctrl: str,
    cond_treat: str,
    min_donors: int = 3,
    min_gene_counts: int = 20,
    padj_thresh: float = 0.05,
    lfc_thresh: float = 1.0,
    allowed_cell_types: set[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run PyDESeq2 differential expression separately for each qualifying cell type."""
    all_res: list[pd.DataFrame] = []
    summary_rows: list[dict[str, object]] = []

    for ct, m_ct in meta_df.groupby("cell_type"):
        if allowed_cell_types and ct not in allowed_cell_types:
            continue
        m_sub = m_ct[m_ct["condition"].isin([cond_ctrl, cond_treat])].copy()
        n_ctrl = int((m_sub["condition"] == cond_ctrl).sum())
        n_treat = int((m_sub["condition"] == cond_treat).sum())
        if n_ctrl < min_donors or n_treat < min_donors:
            logger.info("Skipping %s (ctrl=%d, treat=%d < min_donors=%d).", ct, n_ctrl, n_treat, min_donors)
            continue

        c_sub = counts_df.loc[m_sub.index]
        # Filter lowly expressed genes per cell type
        keep_genes = (c_sub.sum(axis=0) >= min_gene_counts) & ((c_sub > 0).sum(axis=0) >= min_donors)
        c_sub = c_sub.loc[:, keep_genes]
        if c_sub.shape[1] < 10:
            continue

        logger.info("Running PyDESeq2 for %s: %d samples (%d ctrl, %d treat), %d genes.", ct, len(m_sub), n_ctrl, n_treat, c_sub.shape[1])
        dds = make_deseq_dataset(c_sub, m_sub)
        dds.deseq2()
        stat_res = DeseqStats(dds, contrast=["condition", cond_treat, cond_ctrl], n_cpus=1, quiet=True)
        stat_res.summary()

        res = stat_res.results_df.reset_index().rename(columns={"index": "gene"})
        res["cell_type"] = ct
        res["significant"] = (res["padj"] < padj_thresh) & (res["log2FoldChange"].abs() >= lfc_thresh)
        all_res.append(res)

        n_up = int(((res["padj"] < padj_thresh) & (res["log2FoldChange"] >= lfc_thresh)).sum())
        n_down = int(((res["padj"] < padj_thresh) & (res["log2FoldChange"] <= -lfc_thresh)).sum())
        summary_rows.append(
            {
                "cell_type": ct,
                "n_ctrl_donors": n_ctrl,
                "n_treat_donors": n_treat,
                "n_genes_tested": int(len(res)),
                "n_deg_up": n_up,
                "n_deg_down": n_down,
                "n_deg_total": n_up + n_down,
            }
        )

    if not all_res:
        raise ValueError("PyDESeq2 did not produce results for any cell type.")

    deseq_all = pd.concat(all_res, ignore_index=True)
    summary_df = pd.DataFrame(summary_rows).sort_values("n_deg_total", ascending=False).reset_index(drop=True)
    return deseq_all, summary_df


def benchmark_pseudobulk_vs_wilcoxon(
    adata: ad.AnnData,
    deseq_all: pd.DataFrame,
    top_celltype: str,
    cond_ctrl: str,
    cond_treat: str,
    padj_thresh: float = 0.05,
    lfc_thresh: float = 1.0,
) -> pd.DataFrame:
    """Compare pseudobulk PyDESeq2 DEGs vs. naive cell-level Wilcoxon DEGs on top cell type."""
    sub = adata[(adata.obs["cell_type"] == top_celltype) & (adata.obs["condition"].isin([cond_ctrl, cond_treat]))].copy()
    sc.tl.rank_genes_groups(sub, groupby="condition", groups=[cond_treat], reference=cond_ctrl, method="wilcoxon")
    sc_df = sc.get.rank_genes_groups_df(sub, group=cond_treat)
    sc_sig = sc_df[(sc_df["pvals_adj"] < padj_thresh) & (sc_df["logfoldchanges"].abs() >= lfc_thresh)]

    pb_df = deseq_all[deseq_all["cell_type"] == top_celltype]
    pb_sig = pb_df[pb_df["significant"]]

    pb_genes = set(pb_sig["gene"])
    sc_genes = set(sc_sig["names"])
    overlap = len(pb_genes & sc_genes)

    return pd.DataFrame(
        [
            {
                "cell_type": top_celltype,
                "method": "Pseudobulk_PyDESeq2",
                "n_units": int(sub.obs["donor"].nunique() * 2),
                "n_degs": len(pb_genes),
                "n_overlap_with_other": overlap,
            },
            {
                "cell_type": top_celltype,
                "method": "SingleCell_Wilcoxon",
                "n_units": int(sub.n_obs),
                "n_degs": len(sc_genes),
                "n_overlap_with_other": overlap,
            },
        ]
    )


def compute_pathway_enrichment(deseq_all: pd.DataFrame) -> pd.DataFrame:
    """Compute Stouffer-like mean Wald statistic enrichment for immune pathways per cell type."""
    rows: list[dict[str, object]] = []
    for ct, grp in deseq_all.groupby("cell_type"):
        valid = grp.dropna(subset=["stat"]).set_index("gene")
        bg_mean = float(valid["stat"].mean())
        bg_std = float(valid["stat"].std(ddof=1)) or 1.0

        for pathway, gene_list in IMMUNE_GENE_SETS.items():
            present = [g for g in gene_list if g in valid.index]
            if len(present) < 3:
                continue
            vals = valid.loc[present, "stat"].to_numpy(dtype=float)
            mean_stat = float(np.mean(vals))
            z_score = (mean_stat - bg_mean) / (bg_std / np.sqrt(len(present)))
            pval = float(2.0 * norm.sf(abs(z_score)))
            rows.append(
                {
                    "cell_type": ct,
                    "pathway": pathway,
                    "n_genes_present": len(present),
                    "mean_wald_stat": mean_stat,
                    "enrichment_z": float(z_score),
                    "pvalue": pval,
                }
            )

    out = pd.DataFrame(rows)
    if not out.empty:
        _, fdr, _, _ = multipletests(out["pvalue"].to_numpy(), method="fdr_bh")
        out["fdr_bh"] = fdr
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Pseudobulk PyDESeq2 & pathway enrichment.")
    parser.add_argument("--in-h5ad", required=True, type=Path)
    parser.add_argument("--cond-ctrl", default="ctrl")
    parser.add_argument("--cond-treat", default="stim")
    parser.add_argument("--min-pb-cells", type=int, default=10)
    parser.add_argument("--min-pb-counts", type=int, default=1000)
    parser.add_argument("--min-donors", type=int, default=3)
    parser.add_argument("--min-gene-counts", type=int, default=20)
    parser.add_argument("--padj-thresh", type=float, default=0.05)
    parser.add_argument("--lfc-thresh", type=float, default=1.0)
    parser.add_argument("--cell-types", default="")
    parser.add_argument("--out-pca", required=True, type=Path)
    parser.add_argument("--out-deseq", required=True, type=Path)
    parser.add_argument("--out-summary", required=True, type=Path)
    parser.add_argument("--out-compare", required=True, type=Path)
    parser.add_argument("--out-pathways", required=True, type=Path)
    args = parser.parse_args()

    adata = ad.read_h5ad(args.in_h5ad)
    counts_df, meta_df = build_pseudobulk(adata, min_cells=args.min_pb_cells, min_counts=args.min_pb_counts)
    pb_pca = compute_pseudobulk_pca(counts_df, meta_df)

    allowed_cts = {c.strip() for c in args.cell_types.split(",") if c.strip()} or None
    deseq_all, summary_df = run_pydeseq2_per_celltype(
        counts_df,
        meta_df,
        cond_ctrl=args.cond_ctrl,
        cond_treat=args.cond_treat,
        min_donors=args.min_donors,
        min_gene_counts=args.min_gene_counts,
        padj_thresh=args.padj_thresh,
        lfc_thresh=args.lfc_thresh,
        allowed_cell_types=allowed_cts,
    )

    top_ct = str(summary_df.iloc[0]["cell_type"])
    compare_df = benchmark_pseudobulk_vs_wilcoxon(
        adata, deseq_all, top_ct, args.cond_ctrl, args.cond_treat, args.padj_thresh, args.lfc_thresh
    )
    pathways_df = compute_pathway_enrichment(deseq_all)

    for p in (args.out_pca, args.out_deseq, args.out_summary, args.out_compare, args.out_pathways):
        p.parent.mkdir(parents=True, exist_ok=True)
    pb_pca.to_csv(args.out_pca, sep="\t", index=False)
    deseq_all.to_csv(args.out_deseq, sep="\t", index=False)
    summary_df.to_csv(args.out_summary, sep="\t", index=False)
    compare_df.to_csv(args.out_compare, sep="\t", index=False)
    pathways_df.to_csv(args.out_pathways, sep="\t", index=False)
    logger.info("Completed pseudobulk PyDESeq2 across %d cell types.", len(summary_df))


if __name__ == "__main__":
    main()
