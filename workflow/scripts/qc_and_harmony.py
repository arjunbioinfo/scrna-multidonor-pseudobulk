#!/usr/bin/env python3
"""QC filtering, Harmony multi-donor batch integration, UMAP embedding, and mixing metrics."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import anndata as ad
import harmonypy as hm
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def stratified_subsample(
    adata: ad.AnnData, max_per_group: int, seed: int = 42
) -> ad.AnnData:
    """Deterministically cap cells per (donor, condition) while preserving all cell types."""
    if max_per_group <= 0:
        return adata

    rng = np.random.default_rng(seed)
    keep_idx: list[int] = []
    obs = adata.obs.reset_index(drop=True)
    for _, grp in obs.groupby(["donor", "condition"], observed=True):
        idx = grp.index.to_numpy()
        if len(idx) > max_per_group:
            idx = rng.choice(idx, size=max_per_group, replace=False)
        keep_idx.extend(sorted(idx.tolist()))

    logger.info("Stratified subsampling: %d -> %d cells (cap=%d/donor/condition).", adata.n_obs, len(keep_idx), max_per_group)
    return adata[keep_idx].copy()


def knn_mixing_entropy(embedding: np.ndarray, labels: np.ndarray, k: int = 30) -> float:
    """Compute normalized kNN Shannon entropy (0 = unmixed, 1 = perfectly mixed across batches)."""
    n_samples = len(labels)
    k_eff = min(k, max(2, n_samples - 1))
    nn = NearestNeighbors(n_neighbors=k_eff + 1, metric="euclidean").fit(embedding)
    _, indices = nn.kneighbors(embedding)
    nbr_indices = indices[:, 1:]

    uniq, encoded = np.unique(labels, return_inverse=True)
    n_classes = len(uniq)
    if n_classes <= 1:
        return 0.0

    max_ent = np.log(n_classes)
    entropies = np.empty(n_samples, dtype=float)
    for i in range(n_samples):
        counts = np.bincount(encoded[nbr_indices[i]], minlength=n_classes)
        probs = counts[counts > 0] / float(k_eff)
        entropies[i] = -np.sum(probs * np.log(probs)) / max_ent
    return float(np.mean(entropies))


def run_qc_and_harmony(
    adata: ad.AnnData,
    min_genes: int = 200,
    min_cells: int = 3,
    max_cells_per_group: int = 0,
    n_hvgs: int = 2000,
    n_pcs: int = 30,
    seed: int = 42,
) -> tuple[ad.AnnData, pd.DataFrame]:
    """Run QC, normalization, PCA, Harmony batch correction, UMAPs, and batch metrics."""
    np.random.seed(seed)
    adata = adata.copy()

    # Work from raw counts layer so normalization is always consistent
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()

    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata = stratified_subsample(adata, max_cells_per_group, seed=seed)

    # Refresh counts layer after filtering
    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    n_top = min(n_hvgs, adata.n_vars)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor="seurat")
    hvg_mask = adata.var["highly_variable"].to_numpy()
    if hvg_mask.sum() < 5:
        hvg_mask = np.ones(adata.n_vars, dtype=bool)

    adata_hvg = adata[:, hvg_mask].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    n_comps = min(n_pcs, adata_hvg.n_obs - 1, adata_hvg.n_vars - 1)
    sc.pp.pca(adata_hvg, n_comps=n_comps, random_state=seed)
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"].copy()

    # Pre-Harmony UMAP
    sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=15, random_state=seed)
    sc.tl.umap(adata, random_state=seed)
    adata.obsm["X_umap_unintegrated"] = adata.obsm["X_umap"].copy()

    # Harmony across donors
    ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, ["donor"], random_state=seed)
    z_corr = np.asarray(ho.Z_corr)
    if z_corr.shape[0] != adata.n_obs:
        z_corr = z_corr.T
    adata.obsm["X_pca_harmony"] = z_corr

    # Post-Harmony UMAP
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=15, random_state=seed)
    sc.tl.umap(adata, random_state=seed)
    adata.obsm["X_umap_harmony"] = adata.obsm["X_umap"].copy()

    # Quantify batch mixing (donor entropy) and biological conservation (cell-type silhouette)
    donors = adata.obs["donor"].to_numpy()
    cell_types = adata.obs["cell_type"].to_numpy()
    sub_n = min(3000, adata.n_obs)
    rng = np.random.default_rng(seed)
    s_idx = rng.choice(adata.n_obs, size=sub_n, replace=False) if adata.n_obs > sub_n else np.arange(adata.n_obs)

    rows = []
    for stage, rep in (("pre_harmony", "X_pca"), ("post_harmony", "X_pca_harmony")):
        emb = adata.obsm[rep][s_idx]
        donor_ent = knn_mixing_entropy(emb, donors[s_idx], k=30)
        ct_sil = float(silhouette_score(emb, cell_types[s_idx])) if len(np.unique(cell_types[s_idx])) > 1 else 0.0
        donor_sil = float(silhouette_score(emb, donors[s_idx])) if len(np.unique(donors[s_idx])) > 1 else 0.0
        rows.append(
            {
                "stage": stage,
                "donor_knn_entropy": donor_ent,
                "donor_silhouette": donor_sil,
                "celltype_silhouette": ct_sil,
                "n_cells": int(adata.n_obs),
                "n_genes": int(adata.n_vars),
            }
        )

    metrics_df = pd.DataFrame(rows)
    return adata, metrics_df


def main() -> None:
    parser = argparse.ArgumentParser(description="QC and Harmony integration for multi-donor scRNA-seq.")
    parser.add_argument("--in-h5ad", required=True, type=Path)
    parser.add_argument("--out-h5ad", required=True, type=Path)
    parser.add_argument("--out-metrics", required=True, type=Path)
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--min-cells", type=int, default=3)
    parser.add_argument("--max-cells-per-group", type=int, default=0)
    parser.add_argument("--n-hvgs", type=int, default=2000)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    adata = ad.read_h5ad(args.in_h5ad)
    integrated, metrics_df = run_qc_and_harmony(
        adata,
        min_genes=args.min_genes,
        min_cells=args.min_cells,
        max_cells_per_group=args.max_cells_per_group,
        n_hvgs=args.n_hvgs,
        n_pcs=args.n_pcs,
        seed=args.seed,
    )

    args.out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    args.out_metrics.parent.mkdir(parents=True, exist_ok=True)
    integrated.write_h5ad(args.out_h5ad)
    metrics_df.to_csv(args.out_metrics, sep="\t", index=False)
    logger.info("Saved Harmony-integrated object and batch metrics:\n%s", metrics_df.to_string(index=False))


if __name__ == "__main__":
    main()
