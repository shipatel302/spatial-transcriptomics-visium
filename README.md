# spatial-transcriptomics-visium

**End-to-end spatial transcriptomics pipeline for 10x Visium — Tumor Microenvironment Analysis**

[![CI](https://github.com/shipatel302/spatial-transcriptomics-visium/actions/workflows/ci.yml/badge.svg)](https://github.com/shipatel302/spatial-transcriptomics-visium/actions)
![Python](https://img.shields.io/badge/Python-3.10-3776AB?logo=python)
![R](https://img.shields.io/badge/R-4.2-276DC3?logo=r)
![Scanpy](https://img.shields.io/badge/Scanpy-1.9-orange)
![License](https://img.shields.io/badge/license-MIT-lightgrey)

An end-to-end pipeline for analyzing 10x Visium spatial transcriptomics data, with a focus on tumor microenvironment (TME) characterization. Built in **Python (Scanpy/Squidpy)** and **R (RCTD/spacexr)** — covering QC, spatially-aware clustering, spatially variable gene detection, cell type deconvolution, cell-cell communication, and publication-ready visualizations.

---

## Pipeline Overview

```
SpaceRanger Output
    ↓
QC & Preprocessing (Scanpy)
    ↓ Filter spots · Normalize · HVG selection · PCA · UMAP
Spatially-aware Clustering (Leiden + Squidpy spatial graph)
    ↓ Multiple resolutions · Marker genes
Spatially Variable Gene Detection (SpatialDE / Moran's I)
    ↓ Top SVGs · Expression pattern discovery
Cell Type Deconvolution (RCTD / cell2location)
    ↓ Cell type proportions · TME compartment mapping
Cell-Cell Communication (Squidpy ligand-receptor)
    ↓ Neighborhood enrichment · Significant interactions
Publication-ready Visualizations
    ↓ Spatial maps · UMAP · TME overlays · SVG expression
```

---

## Key Features

| Feature | Tools | Output |
|---------|-------|--------|
| Spot QC | Scanpy | Z-score filtering, mito QC, violin plots |
| Spatially-aware clustering | Scanpy + Squidpy | Leiden clusters blended with spatial graph |
| Spatially variable genes | SpatialDE, Moran's I | Ranked SVG table + spatial expression plots |
| Cell type deconvolution | RCTD (R), cell2location | Proportion matrix + spatial maps |
| TME mapping | Custom | Tumor / Immune / Stromal compartment calls |
| Cell-cell communication | Squidpy | Neighborhood enrichment + ligand-receptor |
| Visualizations | Matplotlib + Scanpy | PDF + PNG publication figures |

---

## Repository Structure

```
spatial-transcriptomics-visium/
├── run_pipeline.py                  # Main entry point
├── config/config.yaml               # All parameters
├── pipeline/
│   ├── qc/visium_qc.py             # QC, normalization, PCA, UMAP
│   ├── clustering/spatial_clustering.py  # Leiden + spatial graph
│   ├── spatial_de/spatial_de.py    # SpatialDE / Moran's I
│   ├── deconvolution/deconvolution.py    # RCTD / cell2location
│   ├── cell_communication/communication.py  # Squidpy LR + nhood
│   └── visualization/plots.py      # Publication figures
├── R/
│   └── run_rctd.R                  # RCTD deconvolution (spacexr)
├── scripts/
│   └── run_demo.py                 # Demo run (simulated data)
├── envs/environment.yaml           # Conda environment
└── .github/workflows/ci.yml        # CI: lint + demo
```

---

## Quick Start

### 1. Clone

```bash
git clone https://github.com/shipatel302/spatial-transcriptomics-visium.git
cd spatial-transcriptomics-visium
```

### 2. Set up environment

```bash
conda env create -f envs/environment.yaml
conda activate spatial-tx
```

### 3. Run demo (no data needed)

```bash
python scripts/run_demo.py
```

Generates simulated 10x Visium data (HNSCC tumor microenvironment) and runs all pipeline steps.

### 4. Run on real data

Edit `config/config.yaml`:

```yaml
samples:
  - id:          "MY_SAMPLE"
    spaceranger: "path/to/spaceranger/output"
    condition:   "tumor"
    patient:     "P01"
    tissue:      "HNSCC"
```

Then run:

```bash
# Full pipeline
python run_pipeline.py --config config/config.yaml --outdir results/

# Specific steps
python run_pipeline.py --steps qc,cluster,spatial_de

# Single sample
python run_pipeline.py --sample MY_SAMPLE
```

---

## Input

SpaceRanger output directory containing:
```
<sample>/
├── filtered_feature_bc_matrix/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
└── spatial/
    ├── tissue_positions_list.csv
    ├── scalefactors_json.json
    └── tissue_hires_image.png
```

---

## Outputs

```
results/
├── qc/
│   ├── <sample>_preprocessed.h5ad
│   ├── <sample>_qc_violin.pdf
│   └── qc_summary.tsv
├── cluster/
│   ├── <sample>_clustered.h5ad
│   ├── <sample>_clusters.pdf
│   └── <sample>_markers_dotplot.pdf
├── spatial_de/
│   ├── <sample>_spatially_variable_genes.tsv
│   └── <sample>_top_svgs.pdf
├── deconv/
│   ├── <sample>_cell_proportions.tsv
│   ├── <sample>_tme_map.tsv
│   └── <sample>_deconv_spatial.pdf
├── comm/
│   ├── <sample>_significant_interactions.tsv
│   ├── <sample>_communication.h5ad
│   └── <sample>_communication.pdf
└── viz/
    ├── <sample>_overview.pdf
    ├── <sample>_tme.pdf
    ├── <sample>_top_svgs_spatial.pdf
    └── <sample>_qc_spatial.pdf
```

---

## Biological Context

This pipeline was designed for **HPV-associated head and neck squamous cell carcinoma (HNSCC)** TME analysis but is fully generalizable to any solid tumor spatial dataset. Key biological questions addressed:

- Where are tumor, immune, and stromal compartments spatially located?
- Which genes show spatially restricted expression patterns?
- What is the composition of the immune infiltrate at single-spot resolution?
- Which cell types are in physical proximity and potentially communicating?
- Are there spatially distinct immune niches (e.g., excluded vs. inflamed)?

---

## Author

**Shivani Patel**  
M.S. Bioinformatics Data Science, University of Delaware  
[linkedin.com/in/shivanip99](https://linkedin.com/in/shivanip99) · [github.com/shipatel302](https://github.com/shipatel302)

---

## License

MIT
