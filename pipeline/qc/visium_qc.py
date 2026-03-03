"""
pipeline/qc/visium_qc.py
QC, preprocessing, and normalization for 10x Visium data
Author: Shivani Patel
"""

from pathlib import Path
import logging
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")

log = logging.getLogger(__name__)


def load_visium_sample(sample: dict) -> sc.AnnData:
    """Load a SpaceRanger output into AnnData."""
    sr_path = Path(sample["spaceranger"])

    if sr_path.exists():
        log.info(f"  Loading SpaceRanger output: {sr_path}")
        adata = sc.read_visium(sr_path)
    else:
        log.warning(f"  SpaceRanger path not found: {sr_path}. Generating demo data.")
        adata = make_demo_visium(sample)

    adata.var_names_make_unique()
    adata.obs["sample_id"]  = sample["id"]
    adata.obs["condition"]  = sample["condition"]
    adata.obs["patient"]    = sample["patient"]
    adata.obs["tissue_type"] = sample["tissue"]
    return adata


def run_qc(samples, cfg, outdir, threads=8):
    """Run QC and preprocessing for all samples."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    qc_cfg = cfg["qc"]
    pre_cfg = cfg["preprocessing"]

    adatas = {}
    qc_summary = []

    for sample in samples:
        sid = sample["id"]
        log.info(f"  Processing sample: {sid}")

        adata = load_visium_sample(sample)

        # ── Mitochondrial gene annotation ────────────────────────────────────
        adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
        )

        # ── Pre-filter stats ──────────────────────────────────────────────────
        n_spots_before = adata.n_obs
        n_genes_before = adata.n_vars

        # ── Filter spots ──────────────────────────────────────────────────────
        sc.pp.filter_cells(adata, min_counts=qc_cfg["min_counts"])
        sc.pp.filter_cells(adata, max_counts=qc_cfg["max_counts"])
        sc.pp.filter_cells(adata, min_genes=qc_cfg["min_genes"])
        adata = adata[adata.obs["pct_counts_mt"] < qc_cfg["max_pct_mito"]].copy()

        # ── Filter genes ──────────────────────────────────────────────────────
        sc.pp.filter_genes(adata, min_cells=qc_cfg["min_spots_per_gene"])

        log.info(f"    Spots: {n_spots_before} → {adata.n_obs} "
                 f"| Genes: {n_genes_before} → {adata.n_vars}")

        # ── QC plots ──────────────────────────────────────────────────────────
        plot_qc_metrics(adata, sid, outdir)

        # ── Normalization ──────────────────────────────────────────────────────
        adata = normalize(adata, pre_cfg)

        # ── Dimensionality reduction ──────────────────────────────────────────
        adata = reduce_dimensions(adata, cfg["dim_reduction"])

        # ── Save ──────────────────────────────────────────────────────────────
        out_h5 = outdir / f"{sid}_preprocessed.h5ad"
        adata.write_h5ad(out_h5)
        log.info(f"    Saved: {out_h5}")

        adatas[sid] = adata
        qc_summary.append({
            "sample_id":     sid,
            "condition":     sample["condition"],
            "spots_before":  n_spots_before,
            "spots_after":   adata.n_obs,
            "genes_after":   adata.n_vars,
            "median_counts": int(np.median(adata.obs["total_counts"])),
            "median_genes":  int(np.median(adata.obs["n_genes_by_counts"])),
            "pct_mt_mean":   round(adata.obs["pct_counts_mt"].mean(), 2),
        })

    # ── QC summary table ──────────────────────────────────────────────────────
    qc_df = pd.DataFrame(qc_summary)
    qc_df.to_csv(outdir / "qc_summary.tsv", sep="\t", index=False)
    log.info(f"  QC summary → {outdir / 'qc_summary.tsv'}")

    return adatas


def normalize(adata: sc.AnnData, pre_cfg: dict) -> sc.AnnData:
    """Normalize, log-transform, find HVGs, scale."""
    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=pre_cfg["n_top_genes"],
        flavor="seurat_v3",
        layer="counts",
    )

    adata = adata[:, adata.var["highly_variable"]].copy()

    if pre_cfg.get("regress_out"):
        sc.pp.regress_out(adata, pre_cfg["regress_out"])

    if pre_cfg.get("scale"):
        sc.pp.scale(adata, max_value=10)

    return adata


def reduce_dimensions(adata: sc.AnnData, dim_cfg: dict) -> sc.AnnData:
    """PCA → neighbor graph → UMAP."""
    sc.pp.pca(adata, n_comps=dim_cfg["n_pcs"], svd_solver="arpack")
    sc.pp.neighbors(
        adata,
        n_neighbors=dim_cfg["n_neighbors"],
        n_pcs=40,
    )
    sc.tl.umap(adata, min_dist=dim_cfg["umap_min_dist"])
    return adata


def plot_qc_metrics(adata: sc.AnnData, sample_id: str, outdir: Path):
    """Violin + scatter QC plots."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(f"QC Metrics — {sample_id}", fontsize=13)

    metrics = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    labels  = ["Total UMI Counts", "Genes Detected", "% Mitochondrial"]
    colors  = ["#264653", "#2a9d8f", "#e76f51"]

    for ax, metric, label, color in zip(axes, metrics, labels, colors):
        vals = adata.obs[metric].values
        ax.violinplot(vals, showmedians=True)
        ax.set_title(label, fontsize=10)
        ax.set_ylabel(label)
        ax.get_xaxis().set_visible(False)
        parts = ax.violinplot(vals, showmedians=True)
        for pc in parts["bodies"]:
            pc.set_facecolor(color)
            pc.set_alpha(0.7)

    plt.tight_layout()
    plt.savefig(outdir / f"{sample_id}_qc_violin.pdf", bbox_inches="tight")
    plt.savefig(outdir / f"{sample_id}_qc_violin.png", dpi=150, bbox_inches="tight")
    plt.close()


def make_demo_visium(sample: dict) -> sc.AnnData:
    """Generate realistic synthetic Visium AnnData for demo purposes."""
    import scipy.sparse as sp

    np.random.seed(42)
    n_spots = 2000
    n_genes = 5000

    # Simulate count matrix (negative binomial-like)
    counts = np.random.negative_binomial(2, 0.3, size=(n_spots, n_genes)).astype(np.float32)
    # Add spatially coherent signal for a subset of genes
    x_coords = np.random.uniform(0, 6400, n_spots)
    y_coords = np.random.uniform(0, 6400, n_spots)

    # Simulate tumor region (top-left quadrant has higher expression of genes 0-100)
    tumor_spots = (x_coords < 3200) & (y_coords < 3200)
    counts[np.ix_(tumor_spots, np.arange(100))] *= 3

    gene_names = [f"GENE_{i:05d}" for i in range(n_genes)]
    # Add some real gene names
    real_genes = ["TP53", "EGFR", "MKI67", "CD8A", "CD4", "FOXP3",
                  "PDCD1", "CD274", "CTLA4", "TMB1", "COL1A1", "ACTA2",
                  "PECAM1", "VWF", "CD68", "CD163", "EPCAM", "KRT5"]
    for i, g in enumerate(real_genes[:min(len(real_genes), n_genes)]):
        gene_names[i] = g

    obs = pd.DataFrame({
        "in_tissue": 1,
        "array_row": np.random.randint(0, 78, n_spots),
        "array_col": np.random.randint(0, 64, n_spots),
    }, index=[f"SPOT_{i:05d}" for i in range(n_spots)])

    var = pd.DataFrame(index=gene_names)
    var["gene_ids"]    = [f"ENSG{i:011d}" for i in range(n_genes)]
    var["feature_types"] = "Gene Expression"

    adata = sc.AnnData(X=sp.csr_matrix(counts), obs=obs, var=var)

    # Spatial coordinates
    adata.obsm["spatial"] = np.column_stack([x_coords, y_coords])

    # Fake tissue image
    adata.uns["spatial"] = {
        sample["id"]: {
            "images": {"hires": np.ones((600, 600, 3), dtype=np.uint8) * 240},
            "scalefactors": {
                "tissue_hires_scalef": 0.09375,
                "spot_diameter_fullres": 89.43,
            },
        }
    }
    return adata
