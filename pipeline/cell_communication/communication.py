"""
pipeline/cell_communication/communication.py
Spatial cell-cell communication analysis using Squidpy
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

log = logging.getLogger(__name__)


def run_cell_communication(samples, cfg, outdir, threads=8):
    """Run spatial cell-cell communication analysis."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    comm_cfg  = cfg["communication"]
    deconv_dir = Path("results/deconv")
    clust_dir  = Path("results/cluster")

    for sample in samples:
        sid = sample["id"]
        log.info(f"  Cell-cell communication: {sid}")

        h5_path = clust_dir / f"{sid}_clustered.h5ad"
        if not h5_path.exists():
            h5_path = Path("results/qc") / f"{sid}_preprocessed.h5ad"

        if not h5_path.exists():
            log.warning(f"  No AnnData found — skipping {sid}")
            continue

        adata = sc.read_h5ad(h5_path)

        # Load TME compartments if available
        tme_path = deconv_dir / f"{sid}_tme_map.tsv"
        if tme_path.exists():
            tme_df = pd.read_csv(tme_path, sep="\t", index_col=0)
            shared = tme_df.index.intersection(adata.obs_names)
            adata.obs.loc[shared, "tme_compartment"] = tme_df.loc[shared, "tme_compartment"]
            cell_type_key = "tme_compartment"
        else:
            cell_type_key = "leiden"

        # ── Spatial neighbor graph ──────────────────────────────────────────
        if "spatial_connectivities" not in adata.obsp:
            sq.gr.spatial_neighbors(
                adata,
                coord_type = "visium",
                n_rings    = comm_cfg.get("n_rings", 2),
            )

        # ── Neighborhood enrichment ────────────────────────────────────────
        log.info("    Running neighborhood enrichment...")
        sq.gr.nhood_enrichment(
            adata,
            cluster_key  = cell_type_key,
            n_perms      = comm_cfg.get("n_permutations", 1000),
            seed         = 42,
            show_progress_bar = False,
        )

        # ── Ligand-receptor interaction (Squidpy) ───────────────────────────
        log.info("    Running ligand-receptor analysis...")
        try:
            sq.gr.ligrec(
                adata,
                n_perms       = comm_cfg.get("n_permutations", 1000),
                cluster_key   = cell_type_key,
                use_raw       = False,
                show_progress_bar = False,
                copy          = False,
            )
            lr_results = adata.uns["ligrec"]
            save_ligrec_results(lr_results, sid, outdir)
        except Exception as e:
            log.warning(f"    Ligand-receptor analysis failed: {e}")
            lr_results = None

        # ── Moran's I for cell type co-occurrence ───────────────────────────
        log.info("    Computing spatial autocorrelation...")
        try:
            sq.gr.spatial_autocorr(
                adata,
                mode   = "moran",
                n_perms = 100,
                n_jobs  = threads,
                show_progress_bar = False,
            )
        except Exception as e:
            log.warning(f"    Moran's I failed: {e}")

        # ── Plots ───────────────────────────────────────────────────────────
        plot_communication(adata, cell_type_key, sid, outdir)

        # ── Save ────────────────────────────────────────────────────────────
        adata.write_h5ad(outdir / f"{sid}_communication.h5ad")
        log.info(f"    Saved communication results: {sid}")


def save_ligrec_results(lr_results: dict, sample_id: str, outdir: Path):
    """Save ligand-receptor interaction results."""
    if lr_results is None:
        return

    # Extract significant interactions
    if "pvalues" in lr_results:
        pvals = lr_results["pvalues"]
        means = lr_results.get("means", None)

        sig_mask = pvals < 0.05
        sig_pairs = np.where(sig_mask)

        if means is not None and len(sig_pairs[0]) > 0:
            rows = []
            for i, j in zip(*sig_pairs):
                source = pvals.index[i] if hasattr(pvals, "index") else i
                target = pvals.columns[j] if hasattr(pvals, "columns") else j
                rows.append({
                    "source":       str(source),
                    "target":       str(target),
                    "pvalue":       float(pvals.iloc[i, j]),
                    "mean_expr":    float(means.iloc[i, j])
                    if means is not None else np.nan,
                })
            pd.DataFrame(rows).to_csv(
                outdir / f"{sample_id}_significant_interactions.tsv",
                sep="\t", index=False
            )


def plot_communication(adata, cell_type_key, sample_id, outdir):
    """Neighborhood enrichment + co-occurrence plots."""
    try:
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))

        # Neighborhood enrichment heatmap
        sq.pl.nhood_enrichment(
            adata,
            cluster_key = cell_type_key,
            ax          = axes[0],
            title       = f"Neighborhood Enrichment — {sample_id}",
            figsize     = (8, 6),
        )

        # Spatial scatter coloured by cell type
        sc.pl.spatial(
            adata,
            color    = cell_type_key,
            ax       = axes[1],
            title    = f"Spatial — {cell_type_key}",
            spot_size = 1.5,
            show     = False,
        )

        plt.suptitle(f"Cell-Cell Communication — {sample_id}", fontsize=13)
        plt.tight_layout()
        plt.savefig(outdir / f"{sample_id}_communication.pdf", bbox_inches="tight")
        plt.savefig(outdir / f"{sample_id}_communication.png", dpi=150, bbox_inches="tight")
        plt.close()

    except Exception as e:
        log.warning(f"  Communication plot failed: {e}")
