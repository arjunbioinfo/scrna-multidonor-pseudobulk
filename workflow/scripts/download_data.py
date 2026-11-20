#!/usr/bin/env python3
"""Download and standardize the public Kang et al. (2018) 8-donor PBMC scRNA-seq dataset."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import anndata as ad
import numpy as np
import requests
from scipy import sparse

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

DONOR_CANDIDATES = ("donor", "replicate", "patient", "ind", "Donor", "sample", "batch")
CONDITION_CANDIDATES = ("condition", "label", "stim", "perturbation", "treatment")
CELLTYPE_CANDIDATES = ("cell_type", "cell", "seurat_annotations", "celltype", "cluster")


def find_column(columns: list[str], candidates: tuple[str, ...], role: str) -> str:
    """Locate the first matching metadata column name."""
    for cand in candidates:
        if cand in columns:
            return cand
    raise ValueError(f"Could not find a {role} column in adata.obs (available: {columns}).")


def standardize_obs(adata: ad.AnnData) -> ad.AnnData:
    """Standardize obs metadata columns to donor, condition, and cell_type."""
    cols = list(adata.obs.columns)
    donor_col = find_column(cols, DONOR_CANDIDATES, "donor")
    cond_col = find_column(cols, CONDITION_CANDIDATES, "condition")
    ct_col = find_column(cols, CELLTYPE_CANDIDATES, "cell_type")

    adata.obs["donor"] = adata.obs[donor_col].astype(str)
    adata.obs["condition"] = adata.obs[cond_col].astype(str)
    adata.obs["cell_type"] = adata.obs[ct_col].astype(str)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    return adata


def ensure_raw_counts_layer(adata: ad.AnnData) -> ad.AnnData:
    """Ensure non-negative integer counts are stored in adata.layers['counts']."""
    if "counts" in adata.layers:
        mat = adata.layers["counts"]
    elif adata.raw is not None and adata.raw.shape == adata.shape:
        mat = adata.raw.X
    else:
        mat = adata.X

    if sparse.issparse(mat):
        mat = mat.tocsr().copy()
        mat.data = np.rint(np.clip(mat.data, 0, None)).astype(np.int64)
    else:
        mat = np.rint(np.clip(np.asarray(mat), 0, None)).astype(np.int64)

    adata.layers["counts"] = mat
    return adata


def download_file(url: str, dest: Path, timeout: int = 180) -> int:
    """Stream-download a remote URL if not already cached."""
    if dest.exists() and dest.stat().st_size > 0:
        logger.info("Using cached dataset: %s (%d bytes)", dest, dest.stat().st_size)
        return dest.stat().st_size

    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    logger.info("Downloading %s -> %s", url, dest)
    with requests.get(url, stream=True, timeout=timeout, headers={"User-Agent": "scrna-multidonor-pseudobulk/1.0"}) as r:
        r.raise_for_status()
        with open(tmp, "wb") as fh:
            for chunk in r.iter_content(chunk_size=1024 * 512):
                if chunk:
                    fh.write(chunk)
    with open(tmp, "rb") as fh:
        head = fh.read(520)
    # HDF5 files (and therefore .h5ad) start with this signature at offset 0 or 512.
    signature = b"\x89HDF\r\n\x1a\n"
    if not (head[:8] == signature or head[512:520] == signature):
        n = tmp.stat().st_size
        tmp.unlink()
        raise RuntimeError(
            f"Download from {url} is not an HDF5/.h5ad file ({n} bytes). "
            "The host may be blocking the request; download the file manually into "
            f"{dest} and re-run."
        )
    tmp.rename(dest)
    size = dest.stat().st_size
    logger.info("Saved %s (%d bytes)", dest, size)
    return size


def main() -> None:
    parser = argparse.ArgumentParser(description="Download and standardize multi-donor scRNA-seq h5ad.")
    parser.add_argument("--url", required=True, help="Direct HTTPS URL to .h5ad file.")
    parser.add_argument("--out", required=True, type=Path, help="Output standardized .h5ad path.")
    parser.add_argument("--manifest", required=True, type=Path, help="Output JSON provenance manifest.")
    args = parser.parse_args()

    raw_cache = args.out.with_suffix(".raw.h5ad")
    n_bytes = download_file(args.url, raw_cache)

    adata = ad.read_h5ad(raw_cache)
    adata = standardize_obs(adata)
    adata = ensure_raw_counts_layer(adata)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out)

    manifest = {
        "source_url": args.url,
        "raw_bytes": n_bytes,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_donors": int(adata.obs["donor"].nunique()),
        "conditions": sorted(adata.obs["condition"].unique().tolist()),
        "cell_types": sorted(adata.obs["cell_type"].unique().tolist()),
    }
    with open(args.manifest, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)
    logger.info(
        "Standardized dataset: %d cells, %d genes, %d donors, %d cell types.",
        manifest["n_cells"],
        manifest["n_genes"],
        manifest["n_donors"],
        len(manifest["cell_types"]),
    )


if __name__ == "__main__":
    main()
