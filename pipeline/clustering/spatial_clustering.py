"""
pipeline/clustering/spatial_clustering.py
Spatially-aware clustering: Leiden + spatial graph integration
Author: Shivani Patel
"""

from pathlib import Path
import logging
import numpy as np
import scanpy as sc
import squidpy as sq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

log = logging.getLogger(__name__)


def run_clustering(samples, cfg, outdir, threads=8):
    """Spatially-aware Leiden clustering across all samples."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    clust_cfg = cfg["clustering"]
    qc_dir    = Path("results/qc")

    for sample in samples:
        sid = sample["id"]
        log.info(f"  Clustering sample: {sid}")

        h5_path = qc_dir / f"{sid}_preprocessed.h5ad"
        if not h5_path.exists():
            log.warning(f"  Preprocessed file not found: {h5_path} — skipping.")
            continue

        adata = sc.read_h5ad(h5_path)

        # ── Spatial neighbor graph (Squidpy) ──────────────────────────────────
        if clust_cfg["use_spatial"]:
            log.info("    Building spatial neighbor graph...")
            sq.gr.spatial_neighbors(
                adata,
                coord_type="visium",
                n_rings=2,
                key_added="spatial_connectivities",
            )
            # Blend transcriptomic + spatial connectivity
            import scipy.sparse as sp
            alpha = 0.3   # spatial weight
            adata.obsp["blended_connectivities"] = (
                (1 - alpha) * adata.obsp["connectivities"] +
                alpha * adata.obsp["spatial_connectivities"]
            )
            conn_key = "blended_connectivities"
        else:
            conn_key = "connectivities"

        # ── Leiden at multiple resolutions ────────────────────────────────────
        best_res    = None
        best_adata  = None
        best_silhouette = -1

        for res in clust_cfg["resolution"]:
            sc.tl.leiden(
                adata,
                resolution   = res,
                obsp         = conn_key,
                key_added    = f"leiden_{res}",
                random_state = 42,
            )
            n_clusters = adata.obs[f"leiden_{res}"].nunique()
            log.info(f"    Resolution {res}: {n_clusters} clusters")

        # Use middle resolution as default
        default_res = clust_cfg["resolution"][len(clust_cfg["resolution"]) // 2]
        adata.obs["leiden"] = adata.obs[f"leiden_{default_res}"]
        log.info(f"    Default clustering: leiden_{default_res} "
                 f"({adata.obs['leiden'].nunique()} clusters)")

        # ── Marker genes per cluster ──────────────────────────────────────────
        log.info("    Finding cluster marker genes...")
        sc.tl.rank_genes_groups(
            adata,
            groupby  = "leiden",
            method   = "wilcoxon",
            use_raw  = True,
            pts      = True,
        )

        # ── Plots ─────────────────────────────────────────────────────────────
        plot_clusters(adata, sid, outdir, clust_cfg["resolution"])

        # ── Save ──────────────────────────────────────────────────────────────
        out_h5 = outdir / f"{sid}_clustered.h5ad"
        adata.write_h5ad(out_h5)
        log.info(f"    Saved: {out_h5}")


def plot_clusters(adata, sample_id, outdir, resolutions):
    """Spatial + UMAP cluster plots."""
    n_res = len(resolutions)
    fig, axes = plt.subplots(2, n_res + 1, figsize=(5 * (n_res + 1), 10))

    # UMAP row
    for i, res in enumerate(resolutions):
        key = f"leiden_{res}"
        sc.pl.umap(adata, color=key, ax=axes[0, i],
                   title=f"UMAP — leiden {res}", show=False, frameon=False)

    # Default leiden UMAP
    sc.pl.umap(adata, color="leiden", ax=axes[0, -1],
               title="UMAP — leiden (default)", show=False, frameon=False)

    # Spatial row
    for i, res in enumerate(resolutions):
        key = f"leiden_{res}"
        sc.pl.spatial(adata, color=key, ax=axes[1, i],
                      title=f"Spatial — leiden {res}",
                      spot_size=1.5, show=False)

    sc.pl.spatial(adata, color="leiden", ax=axes[1, -1],
                  title="Spatial — leiden (default)",
                  spot_size=1.5, show=False)

    plt.suptitle(f"Spatially-aware Clustering — {sample_id}", fontsize=13)
    plt.tight_layout()
    plt.savefig(outdir / f"{sample_id}_clusters.pdf", bbox_inches="tight")
    plt.savefig(outdir / f"{sample_id}_clusters.png", dpi=150, bbox_inches="tight")
    plt.close()

    # ── Marker gene dotplot ───────────────────────────────────────────────────
    try:
        fig2, ax2 = plt.subplots(figsize=(14, 6))
        sc.pl.rank_genes_groups_dotplot(
            adata, n_genes=5, ax=ax2,
            title=f"Top Marker Genes per Cluster — {sample_id}",
            show=False
        )
        plt.tight_layout()
        plt.savefig(outdir / f"{sample_id}_markers_dotplot.pdf", bbox_inches="tight")
        plt.close()
    except Exception as e:
        log.warning(f"  Marker dotplot failed: {e}")
