"""Unit and synthetic end-to-end smoke tests for scrna-multidonor-pseudobulk."""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

from download_data import ensure_raw_counts_layer, standardize_obs  # noqa: E402
from plot_summary import (  # noqa: E402
    plot_harmony_umaps,
    plot_pathway_heatmap,
    plot_pseudobulk_pca_and_degs,
    plot_volcano_and_comparison,
)
from pseudobulk_deseq2 import (  # noqa: E402
    IMMUNE_GENE_SETS,
    benchmark_pseudobulk_vs_wilcoxon,
    build_pseudobulk,
    compute_pathway_enrichment,
    compute_pseudobulk_pca,
    run_pydeseq2_per_celltype,
)
from qc_and_harmony import run_qc_and_harmony  # noqa: E402


def make_synthetic_multidonor_adata(seed: int = 42) -> ad.AnnData:
    """Create a 6-donor x 2-condition x 2-cell-type AnnData with injected IFN-beta response."""
    rng = np.random.default_rng(seed)
    isg_genes = IMMUNE_GENE_SETS["Type_I_Interferon_Response"][:12]
    chem_genes = IMMUNE_GENE_SETS["Chemokine_Cytokine_Signaling"][:8]
    bg_genes = [f"GENE_{i:03d}" for i in range(60)]
    genes = isg_genes + chem_genes + bg_genes
    n_genes = len(genes)

    rows = []
    obs_records = []
    donors = [f"D{i}" for i in range(1, 7)]
    conditions = ["ctrl", "stim"]
    cell_types = ["CD14_Mono", "CD4_T"]

    for d_idx, donor in enumerate(donors):
        donor_shift = rng.uniform(0.8, 1.3, size=n_genes)
        for cond in conditions:
            for ct in cell_types:
                n_cells = 18
                base_rate = np.full(n_genes, 3.0)
                # Cell-type identity signal
                if ct == "CD14_Mono":
                    base_rate[20:35] = 12.0
                else:
                    base_rate[35:50] = 12.0
                # Strong IFN-beta upregulation in stim across donors
                if cond == "stim":
                    base_rate[:12] = 28.0 if ct == "CD14_Mono" else 16.0
                    base_rate[12:20] = 18.0 if ct == "CD14_Mono" else 6.0

                lam = base_rate * donor_shift
                counts = rng.poisson(lam=lam, size=(n_cells, n_genes))
                rows.append(counts)
                for c_i in range(n_cells):
                    obs_records.append(
                        {
                            "replicate": donor,
                            "label": cond,
                            "cell_type": ct,
                            "cell_id": f"{donor}_{cond}_{ct}_{c_i}",
                        }
                    )

    X = sparse.csr_matrix(np.vstack(rows).astype(np.int64))
    obs = pd.DataFrame(obs_records).set_index("cell_id")
    var = pd.DataFrame(index=pd.Index(genes, name="gene"))
    return ad.AnnData(X=X, obs=obs, var=var)


def test_standardize_and_end_to_end_pipeline(tmp_path: Path) -> None:
    adata = make_synthetic_multidonor_adata(seed=42)
    adata = standardize_obs(adata)
    adata = ensure_raw_counts_layer(adata)
    assert {"donor", "condition", "cell_type"}.issubset(adata.obs.columns)
    assert "counts" in adata.layers

    integrated, metrics_df = run_qc_and_harmony(
        adata,
        min_genes=10,
        min_cells=3,
        max_cells_per_group=0,
        n_hvgs=60,
        n_pcs=15,
        seed=42,
    )
    assert len(metrics_df) == 2
    assert "X_pca_harmony" in integrated.obsm
    assert "X_umap_harmony" in integrated.obsm

    counts_df, meta_df = build_pseudobulk(integrated, min_cells=10, min_counts=100)
    assert len(meta_df) == 24  # 6 donors x 2 conditions x 2 cell types

    pb_pca = compute_pseudobulk_pca(counts_df, meta_df)
    assert {"PC1", "PC2"}.issubset(pb_pca.columns)

    deseq_all, summary_df = run_pydeseq2_per_celltype(
        counts_df,
        meta_df,
        cond_ctrl="ctrl",
        cond_treat="stim",
        min_donors=3,
        min_gene_counts=10,
        padj_thresh=0.05,
        lfc_thresh=1.0,
    )
    assert len(summary_df) == 2
    assert summary_df["n_deg_up"].sum() >= 10

    top_ct = str(summary_df.iloc[0]["cell_type"])
    pb_vs_sc = benchmark_pseudobulk_vs_wilcoxon(
        integrated, deseq_all, top_ct, "ctrl", "stim", padj_thresh=0.05, lfc_thresh=1.0
    )
    assert len(pb_vs_sc) == 2

    pathways_df = compute_pathway_enrichment(deseq_all)
    assert not pathways_df.empty
    ifn_rows = pathways_df[pathways_df["pathway"] == "Type_I_Interferon_Response"]
    assert (ifn_rows["enrichment_z"] > 2.0).all()

    # Verify all 4 figures render cleanly
    plot_harmony_umaps(integrated, metrics_df, tmp_path / "umap.png")
    plot_pseudobulk_pca_and_degs(pb_pca, summary_df, tmp_path / "pca_deg.png")
    plot_volcano_and_comparison(deseq_all, pb_vs_sc, 0.05, 1.0, ["ISG15", "IFIT1"], tmp_path / "volc.png")
    plot_pathway_heatmap(pathways_df, tmp_path / "path.png")
    for fname in ("umap.png", "pca_deg.png", "volc.png", "path.png"):
        assert (tmp_path / fname).stat().st_size > 1000
