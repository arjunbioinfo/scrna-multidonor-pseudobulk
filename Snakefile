"""Reproducible Snakemake workflow for multi-donor scRNA-seq Harmony integration and pseudobulk PyDESeq2."""

configfile: "config/config.yaml"

DATA_DIR = config.get("output_dirs", {}).get("data", "data/raw")
RESULTS_DIR = config.get("output_dirs", {}).get("results", "results")
FIGURES_DIR = config.get("output_dirs", {}).get("figures", "figures")


rule all:
    input:
        f"{RESULTS_DIR}/batch_metrics.tsv",
        f"{RESULTS_DIR}/pseudobulk_pca.tsv",
        f"{RESULTS_DIR}/pseudobulk_deseq2_all.tsv",
        f"{RESULTS_DIR}/celltype_deg_summary.tsv",
        f"{RESULTS_DIR}/pseudobulk_vs_wilcoxon.tsv",
        f"{RESULTS_DIR}/pathway_enrichment.tsv",
        f"{FIGURES_DIR}/harmony_integration_umap.png",
        f"{FIGURES_DIR}/pseudobulk_pca_and_degs.png",
        f"{FIGURES_DIR}/volcano_and_pb_vs_sc.png",
        f"{FIGURES_DIR}/pathway_enrichment_heatmap.png",


rule download_data:
    output:
        h5ad=f"{DATA_DIR}/kang_2018.h5ad",
        manifest=f"{DATA_DIR}/download_manifest.json",
    params:
        url=config["dataset"]["url"],
    shell:
        """
        python workflow/scripts/download_data.py \
            --url "{params.url}" \
            --out {output.h5ad} \
            --manifest {output.manifest}
        """


rule qc_and_harmony:
    input:
        h5ad=rules.download_data.output.h5ad,
    output:
        integrated=f"{RESULTS_DIR}/harmony_integrated.h5ad",
        metrics=f"{RESULTS_DIR}/batch_metrics.tsv",
    params:
        min_genes=config["qc"]["min_genes"],
        min_cells=config["qc"]["min_cells"],
        max_cells_group=config["qc"]["max_cells_per_donor_condition"],
        n_hvgs=config["integration"]["n_top_genes"],
        n_pcs=config["integration"]["n_pcs"],
        seed=config["random_seed"],
    shell:
        """
        python workflow/scripts/qc_and_harmony.py \
            --in-h5ad {input.h5ad} \
            --out-h5ad {output.integrated} \
            --out-metrics {output.metrics} \
            --min-genes {params.min_genes} \
            --min-cells {params.min_cells} \
            --max-cells-per-group {params.max_cells_group} \
            --n-hvgs {params.n_hvgs} \
            --n-pcs {params.n_pcs} \
            --seed {params.seed}
        """


rule pseudobulk_deseq2:
    input:
        integrated=rules.qc_and_harmony.output.integrated,
    output:
        pb_pca=f"{RESULTS_DIR}/pseudobulk_pca.tsv",
        deseq_all=f"{RESULTS_DIR}/pseudobulk_deseq2_all.tsv",
        deg_summary=f"{RESULTS_DIR}/celltype_deg_summary.tsv",
        pb_vs_sc=f"{RESULTS_DIR}/pseudobulk_vs_wilcoxon.tsv",
        pathways=f"{RESULTS_DIR}/pathway_enrichment.tsv",
    params:
        cond_ctrl=config["comparison"]["control"],
        cond_treat=config["comparison"]["treated"],
        min_pb_cells=config["pseudobulk"]["min_cells_per_sample"],
        min_pb_counts=config["pseudobulk"]["min_counts_per_sample"],
        min_donors=config["pseudobulk"]["min_donors_per_condition"],
        min_gene_counts=config["pseudobulk"]["min_total_gene_counts"],
        padj_thresh=config["pseudobulk"]["padj_threshold"],
        lfc_thresh=config["pseudobulk"]["lfc_threshold"],
        cell_types=",".join(config["pseudobulk"].get("cell_types", [])),
    shell:
        """
        python workflow/scripts/pseudobulk_deseq2.py \
            --in-h5ad {input.integrated} \
            --cond-ctrl "{params.cond_ctrl}" \
            --cond-treat "{params.cond_treat}" \
            --min-pb-cells {params.min_pb_cells} \
            --min-pb-counts {params.min_pb_counts} \
            --min-donors {params.min_donors} \
            --min-gene-counts {params.min_gene_counts} \
            --padj-thresh {params.padj_thresh} \
            --lfc-thresh {params.lfc_thresh} \
            --cell-types "{params.cell_types}" \
            --out-pca {output.pb_pca} \
            --out-deseq {output.deseq_all} \
            --out-summary {output.deg_summary} \
            --out-compare {output.pb_vs_sc} \
            --out-pathways {output.pathways}
        """


rule plot_summary:
    input:
        integrated=rules.qc_and_harmony.output.integrated,
        metrics=rules.qc_and_harmony.output.metrics,
        pb_pca=rules.pseudobulk_deseq2.output.pb_pca,
        deseq_all=rules.pseudobulk_deseq2.output.deseq_all,
        deg_summary=rules.pseudobulk_deseq2.output.deg_summary,
        pb_vs_sc=rules.pseudobulk_deseq2.output.pb_vs_sc,
        pathways=rules.pseudobulk_deseq2.output.pathways,
    output:
        umap_fig=f"{FIGURES_DIR}/harmony_integration_umap.png",
        pca_deg_fig=f"{FIGURES_DIR}/pseudobulk_pca_and_degs.png",
        volc_fig=f"{FIGURES_DIR}/volcano_and_pb_vs_sc.png",
        path_fig=f"{FIGURES_DIR}/pathway_enrichment_heatmap.png",
    params:
        padj_thresh=config["pseudobulk"]["padj_threshold"],
        lfc_thresh=config["pseudobulk"]["lfc_threshold"],
        highlight_genes=",".join(config.get("plotting", {}).get("highlight_genes", ["ISG15", "IFIT1", "IFIT3", "MX1", "IFI6", "CXCL10", "OAS1"])),
    shell:
        """
        python workflow/scripts/plot_summary.py \
            --in-h5ad {input.integrated} \
            --metrics {input.metrics} \
            --pb-pca {input.pb_pca} \
            --deseq-all {input.deseq_all} \
            --deg-summary {input.deg_summary} \
            --pb-vs-sc {input.pb_vs_sc} \
            --pathways {input.pathways} \
            --padj-thresh {params.padj_thresh} \
            --lfc-thresh {params.lfc_thresh} \
            --highlight-genes "{params.highlight_genes}" \
            --out-dir {FIGURES_DIR}
        """
