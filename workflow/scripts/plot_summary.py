#!/usr/bin/env python3
"""Generate publication-quality summary figures for multi-donor scRNA-seq integration & pseudobulk DE."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

sns.set_theme(style="whitegrid", context="paper", font_scale=1.1)


def plot_harmony_umaps(adata: ad.AnnData, metrics_df: pd.DataFrame, out_path: Path) -> None:
    """2x2 UMAP comparison before vs. after Harmony colored by donor and cell type."""
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 9.2))
    donors = adata.obs["donor"].astype(str).to_numpy()
    cts = adata.obs["cell_type"].astype(str).to_numpy()

    u_pre = adata.obsm["X_umap_unintegrated"]
    u_post = adata.obsm["X_umap_harmony"]

    pre_ent = float(metrics_df.loc[metrics_df["stage"] == "pre_harmony", "donor_knn_entropy"].iloc[0])
    post_ent = float(metrics_df.loc[metrics_df["stage"] == "post_harmony", "donor_knn_entropy"].iloc[0])

    def scatter_cat(ax, coords, labels, title, palette_name):
        uniq = sorted(np.unique(labels))
        colors = sns.color_palette(palette_name, n_colors=len(uniq))
        for label, col in zip(uniq, colors):
            mask = labels == label
            ax.scatter(coords[mask, 0], coords[mask, 1], s=5, alpha=0.65, color=col, label=label, edgecolors="none")
        ax.set_title(title, fontsize=10.5, fontweight="bold")
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_xticks([])
        ax.set_yticks([])

    scatter_cat(axes[0, 0], u_pre, donors, f"Pre-Harmony by Donor (kNN entropy={pre_ent:.2f})", "tab10")
    scatter_cat(axes[0, 1], u_post, donors, f"Post-Harmony by Donor (kNN entropy={post_ent:.2f})", "tab10")
    axes[0, 1].legend(title="Donor", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7.5, markerscale=2.2, frameon=True)

    scatter_cat(axes[1, 0], u_pre, cts, "Pre-Harmony by Cell Type", "Set2")
    scatter_cat(axes[1, 1], u_post, cts, "Post-Harmony by Cell Type", "Set2")
    axes[1, 1].legend(title="Cell Type", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7.5, markerscale=2.2, frameon=True)

    fig.suptitle("Multi-Donor Batch Integration with Harmony", fontsize=13, fontweight="bold", y=0.99)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_pseudobulk_pca_and_degs(pb_pca: pd.DataFrame, deg_summary: pd.DataFrame, out_path: Path) -> None:
    """2-panel figure: Pseudobulk PCA and per-cell-type up/down DEG counts."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.2, 4.9))

    sns.scatterplot(
        data=pb_pca,
        x="PC1",
        y="PC2",
        hue="cell_type",
        style="condition",
        s=65,
        alpha=0.85,
        ax=ax1,
    )
    pc1_v = pb_pca["pc1_var_pct"].iloc[0]
    pc2_v = pb_pca["pc2_var_pct"].iloc[0]
    ax1.set_xlabel(f"PC1 ({pc1_v:.1f}% variance)")
    ax1.set_ylabel(f"PC2 ({pc2_v:.1f}% variance)")
    ax1.set_title("Pseudobulk Sample PCA (Donor × Cell Type)", fontweight="bold")
    ax1.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7.5, frameon=True)

    df = deg_summary.sort_values("n_deg_total", ascending=True)
    y = np.arange(len(df))
    ax2.barh(y, df["n_deg_up"], color="#D95F02", label="Up (stim > ctrl)", height=0.6)
    ax2.barh(y, -df["n_deg_down"], color="#1B9E77", label="Down (stim < ctrl)", height=0.6)
    ax2.axvline(0, color="#333333", lw=0.9)
    ax2.set_yticks(y)
    ax2.set_yticklabels(df["cell_type"], fontsize=9)
    ax2.set_xlabel("Significant DEGs (PyDESeq2, padj < 0.05, |LFC| ≥ 1)")
    ax2.set_title("Cell-Type-Specific Pseudobulk DEGs", fontweight="bold")
    ax2.legend(frameon=True, fontsize=8.5, loc="lower right")

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_volcano_and_comparison(
    deseq_all: pd.DataFrame,
    pb_vs_sc: pd.DataFrame,
    padj_thresh: float,
    lfc_thresh: float,
    highlight_genes: list[str],
    out_path: Path,
) -> None:
    """2-panel figure: Top cell-type PyDESeq2 volcano + Pseudobulk vs. Single-Cell Wilcoxon DEG count."""
    top_ct = str(pb_vs_sc.iloc[0]["cell_type"])
    df = deseq_all[deseq_all["cell_type"] == top_ct].dropna(subset=["padj", "log2FoldChange"]).copy()
    df["neg_log10_padj"] = -np.log10(np.clip(df["padj"].to_numpy(dtype=float), 1e-300, 1.0))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.8, 4.9), gridspec_kw={"width_ratios": [1.55, 1.0]})

    sig_up = (df["padj"] < padj_thresh) & (df["log2FoldChange"] >= lfc_thresh)
    sig_dn = (df["padj"] < padj_thresh) & (df["log2FoldChange"] <= -lfc_thresh)
    ns = ~(sig_up | sig_dn)

    ax1.scatter(df.loc[ns, "log2FoldChange"], df.loc[ns, "neg_log10_padj"], s=10, color="#BDBDBD", alpha=0.5, label="NS")
    ax1.scatter(df.loc[sig_up, "log2FoldChange"], df.loc[sig_up, "neg_log10_padj"], s=16, color="#D95F02", alpha=0.8, label=f"Up (n={sig_up.sum()})")
    ax1.scatter(df.loc[sig_dn, "log2FoldChange"], df.loc[sig_dn, "neg_log10_padj"], s=16, color="#1B9E77", alpha=0.8, label=f"Down (n={sig_dn.sum()})")
    ax1.axhline(-np.log10(padj_thresh), color="#555555", ls="--", lw=0.8)
    ax1.axvline(lfc_thresh, color="#555555", ls="--", lw=0.8)
    ax1.axvline(-lfc_thresh, color="#555555", ls="--", lw=0.8)

    # Label canonical ISGs or top significant genes
    hl_df = df[df["gene"].isin(highlight_genes) & (sig_up | sig_dn)]
    if hl_df.empty:
        hl_df = df[sig_up | sig_dn].nlargest(6, "neg_log10_padj")
    for row in hl_df.head(8).itertuples(index=False):
        ax1.annotate(
            row.gene,
            xy=(row.log2FoldChange, row.neg_log10_padj),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=8.2,
            fontweight="bold",
        )

    ax1.set_xlabel("log2 Fold Change (stim / ctrl)")
    ax1.set_ylabel("-log10 adjusted p-value")
    ax1.set_title(f"Pseudobulk PyDESeq2 Volcano: {top_ct}", fontweight="bold")
    ax1.legend(frameon=True, fontsize=8.5, loc="upper left")

    # Right panel: Pseudobulk vs naive single-cell Wilcoxon DEG count
    colors = ["#386CB0", "#E7298A"]
    bars = ax2.bar(pb_vs_sc["method"], pb_vs_sc["n_degs"], color=colors[: len(pb_vs_sc)], width=0.55)
    for bar, row in zip(bars, pb_vs_sc.itertuples(index=False)):
        ax2.text(
            bar.get_x() + bar.get_width() / 2.0,
            bar.get_height() + max(pb_vs_sc["n_degs"].max() * 0.02, 1.0),
            f"{row.n_degs} DEGs\n(n={row.n_units} units)",
            ha="center",
            va="bottom",
            fontsize=8.5,
        )
    ax2.set_ylabel("Significant DEGs called")
    ax2.set_title(f"Pseudobulk vs. Naive Single-Cell ({top_ct})", fontweight="bold")
    ax2.set_ylim(0, max(pb_vs_sc["n_degs"].max() * 1.25, 10))
    ax2.tick_params(axis="x", labelsize=9)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_pathway_heatmap(pathways_df: pd.DataFrame, out_path: Path) -> None:
    """Heatmap of cell-type × immune pathway enrichment z-scores."""
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    if pathways_df.empty:
        ax.text(0.5, 0.5, "No pathways with ≥3 genes tested.", ha="center", va="center")
        fig.savefig(out_path, dpi=300)
        plt.close(fig)
        return

    piv = pathways_df.pivot(index="pathway", columns="cell_type", values="enrichment_z").fillna(0.0)
    piv.index = [p.replace("_", " ") for p in piv.index]
    vmax = max(float(np.abs(piv.to_numpy()).max()), 2.0)

    sns.heatmap(
        piv,
        cmap="RdBu_r",
        center=0.0,
        vmin=-vmax,
        vmax=vmax,
        annot=True,
        fmt=".1f",
        linewidths=0.5,
        cbar_kws={"label": "Pathway Enrichment Z-score (stim vs. ctrl)"},
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Cell-Type Immune Pathway Enrichment (Pseudobulk Wald Statistics)", fontweight="bold")
    plt.xticks(rotation=30, ha="right", fontsize=8.5)
    plt.yticks(fontsize=9)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot multi-donor scRNA-seq integration & pseudobulk DE outputs.")
    parser.add_argument("--in-h5ad", required=True, type=Path)
    parser.add_argument("--metrics", required=True, type=Path)
    parser.add_argument("--pb-pca", required=True, type=Path)
    parser.add_argument("--deseq-all", required=True, type=Path)
    parser.add_argument("--deg-summary", required=True, type=Path)
    parser.add_argument("--pb-vs-sc", required=True, type=Path)
    parser.add_argument("--pathways", required=True, type=Path)
    parser.add_argument("--padj-thresh", type=float, default=0.05)
    parser.add_argument("--lfc-thresh", type=float, default=1.0)
    parser.add_argument("--highlight-genes", default="ISG15,IFIT1,IFIT3,MX1,IFI6,CXCL10,OAS1")
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(args.in_h5ad)
    metrics_df = pd.read_csv(args.metrics, sep="\t")
    pb_pca = pd.read_csv(args.pb_pca, sep="\t")
    deseq_all = pd.read_csv(args.deseq_all, sep="\t")
    deg_summary = pd.read_csv(args.deg_summary, sep="\t")
    pb_vs_sc = pd.read_csv(args.pb_vs_sc, sep="\t")
    pathways_df = pd.read_csv(args.pathways, sep="\t")
    hl = [g.strip() for g in args.highlight_genes.split(",") if g.strip()]

    plot_harmony_umaps(adata, metrics_df, args.out_dir / "harmony_integration_umap.png")
    plot_pseudobulk_pca_and_degs(pb_pca, deg_summary, args.out_dir / "pseudobulk_pca_and_degs.png")
    plot_volcano_and_comparison(deseq_all, pb_vs_sc, args.padj_thresh, args.lfc_thresh, hl, args.out_dir / "volcano_and_pb_vs_sc.png")
    plot_pathway_heatmap(pathways_df, args.out_dir / "pathway_enrichment_heatmap.png")
    logger.info("Saved 4 summary figures to %s", args.out_dir)


if __name__ == "__main__":
    main()
