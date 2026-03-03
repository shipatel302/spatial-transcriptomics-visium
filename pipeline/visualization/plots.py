"""
pipeline/visualization/plots.py
Publication-ready spatial transcriptomics visualizations
Author: Shivani Patel
"""

from pathlib import Path
import logging
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap

log = logging.getLogger(__name__)

# ── Custom colormaps ──────────────────────────────────────────────────────────
TME_CMAP = LinearSegmentedColormap.from_list(
    "tme", ["#f8f9fa", "#2a9d8f", "#e76f51", "#264653"]
)


def run_visualization(samples, cfg, outdir, threads=8):
    """Generate all publication-ready figures."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    viz_cfg   = cfg["visualization"]
    clust_dir = Path("results/cluster")

    for sample in samples:
        sid = sample["id"]
        log.info(f"  Generating figures for: {sid}")

        h5_path = clust_dir / f"{sid}_clustered.h5ad"
        if not h5_path.exists():
            h5_path = Path("results/qc") / f"{sid}_preprocessed.h5ad"
        if not h5_path.exists():
            log.warning(f"  No AnnData — skipping {sid}")
            continue

        adata = sc.read_h5ad(h5_path)

        # Load deconvolution if available
        deconv_path = Path("results/deconv") / f"{sid}_cell_proportions.tsv"
        if deconv_path.exists():
            props = pd.read_csv(deconv_path, sep="\t", index_col=0)
            shared = props.index.intersection(adata.obs_names)
            for ct in props.columns:
                adata.obs.loc[shared, f"prop_{ct}"] = props.loc[shared, ct]

        # Load SVGs
        svg_path = Path("results/spatial_de") / f"{sid}_spatially_variable_genes.tsv"
        top_svgs = []
        if svg_path.exists():
            svg_df    = pd.read_csv(svg_path, sep="\t")
            gene_col  = "gene" if "gene" in svg_df.columns else svg_df.columns[0]
            top_svgs  = svg_df[gene_col].head(6).tolist()
            top_svgs  = [g for g in top_svgs if g in adata.var_names]

        # ── Figure 1: Overview panel ──────────────────────────────────────
        fig_overview(adata, sid, viz_cfg, outdir)

        # ── Figure 2: TME spatial map ─────────────────────────────────────
        if "tme_compartment" in adata.obs.columns:
            fig_tme(adata, sid, viz_cfg, outdir)

        # ── Figure 3: Top SVG expression ─────────────────────────────────
        if top_svgs:
            fig_svgs(adata, top_svgs, sid, viz_cfg, outdir)

        # ── Figure 4: QC on tissue ────────────────────────────────────────
        fig_qc_spatial(adata, sid, viz_cfg, outdir)

        log.info(f"    Figures saved for: {sid}")


def fig_overview(adata, sample_id, viz_cfg, outdir):
    """Overview: UMAP + spatial clusters + gene count."""
    fig = plt.figure(figsize=(18, 6))
    gs  = gridspec.GridSpec(1, 3, figure=fig, wspace=0.35)

    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    ax3 = fig.add_subplot(gs[2])

    # UMAP
    if "X_umap" in adata.obsm:
        sc.pl.umap(adata, color="leiden", ax=ax1,
                   title="UMAP — Leiden Clusters",
                   show=False, frameon=False, legend_loc="on data")

    # Spatial clusters
    sc.pl.spatial(adata, color="leiden", ax=ax2,
                  title="Spatial — Leiden Clusters",
                  spot_size=viz_cfg["spot_size"], show=False)

    # Gene count spatial
    sc.pl.spatial(adata, color="n_genes_by_counts", ax=ax3,
                  title="Genes Detected per Spot",
                  spot_size=viz_cfg["spot_size"],
                  cmap="viridis", show=False)

    plt.suptitle(f"Overview — {sample_id}", fontsize=14, y=1.02)
    _save(fig, outdir / f"{sample_id}_overview", viz_cfg)


def fig_tme(adata, sample_id, viz_cfg, outdir):
    """TME compartment spatial map."""
    prop_cols = [c for c in adata.obs.columns if c.startswith("prop_")]
    n = len(prop_cols) + 1
    n_cols = min(4, n)
    n_rows = int(np.ceil(n / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(5 * n_cols, 4.5 * n_rows))
    axes = np.array(axes).flatten()

    # TME compartment overview
    sc.pl.spatial(adata, color="tme_compartment", ax=axes[0],
                  title="TME Compartment", spot_size=viz_cfg["spot_size"],
                  show=False)

    for i, col in enumerate(prop_cols[:len(axes) - 1], start=1):
        ct_name = col.replace("prop_", "").replace("_", " ")
        sc.pl.spatial(adata, color=col, ax=axes[i],
                      title=ct_name, cmap="RdYlBu_r",
                      spot_size=viz_cfg["spot_size"],
                      vmin=0, show=False)

    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    plt.suptitle(f"Tumor Microenvironment — {sample_id}", fontsize=14)
    plt.tight_layout()
    _save(fig, outdir / f"{sample_id}_tme", viz_cfg)


def fig_svgs(adata, top_svgs, sample_id, viz_cfg, outdir):
    """Spatial expression of top spatially variable genes."""
    n = len(top_svgs)
    n_cols = min(3, n)
    n_rows = int(np.ceil(n / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(5.5 * n_cols, 4.5 * n_rows))
    axes = np.array(axes).flatten()

    for i, gene in enumerate(top_svgs):
        if gene in adata.var_names:
            sc.pl.spatial(adata, color=gene, ax=axes[i],
                          title=gene, cmap="magma",
                          spot_size=viz_cfg["spot_size"],
                          use_raw=True if adata.raw else False,
                          show=False)
        else:
            axes[i].axis("off")

    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    plt.suptitle(f"Top Spatially Variable Genes — {sample_id}", fontsize=14)
    plt.tight_layout()
    _save(fig, outdir / f"{sample_id}_top_svgs_spatial", viz_cfg)


def fig_qc_spatial(adata, sample_id, viz_cfg, outdir):
    """QC metrics overlaid on tissue."""
    metrics = [c for c in
               ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
               if c in adata.obs.columns]
    if not metrics:
        return

    fig, axes = plt.subplots(1, len(metrics),
                              figsize=(6 * len(metrics), 5))
    axes = np.array(axes).flatten()

    labels = {"total_counts": "Total UMI",
              "n_genes_by_counts": "Genes Detected",
              "pct_counts_mt": "% Mito"}

    for i, metric in enumerate(metrics):
        sc.pl.spatial(adata, color=metric, ax=axes[i],
                      title=labels.get(metric, metric),
                      spot_size=viz_cfg["spot_size"],
                      cmap="YlOrRd", show=False)

    plt.suptitle(f"QC Metrics on Tissue — {sample_id}", fontsize=14)
    plt.tight_layout()
    _save(fig, outdir / f"{sample_id}_qc_spatial", viz_cfg)


def _save(fig, base_path: Path, viz_cfg: dict):
    """Save figure in all requested formats."""
    for fmt in viz_cfg.get("save_formats", ["pdf", "png"]):
        path = base_path.with_suffix(f".{fmt}")
        dpi  = viz_cfg.get("dpi", 300) if fmt == "png" else None
        kwargs = {"bbox_inches": "tight"}
        if dpi:
            kwargs["dpi"] = dpi
        fig.savefig(path, **kwargs)
    plt.close(fig)
