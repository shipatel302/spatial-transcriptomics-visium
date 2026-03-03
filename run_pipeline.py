#!/usr/bin/env python3
"""
spatial-transcriptomics-visium
==============================
End-to-end spatial transcriptomics pipeline for 10x Visium data

Pipeline:
    SpaceRanger QC → Scanpy preprocessing → Spatially-aware clustering
    → SpatialDE variable gene detection → Cell type deconvolution (STdeconvolve/RCTD)
    → Tumor microenvironment mapping → Cell-cell communication (Squidpy/CellChat)
    → Spatially resolved differential expression → Publication-ready visualizations

Author: Shivani Patel
"""

import argparse
import logging
import sys
from pathlib import Path

import yaml

from pipeline.qc.visium_qc import run_qc
from pipeline.clustering.spatial_clustering import run_clustering
from pipeline.spatial_de.spatial_de import run_spatial_de
from pipeline.deconvolution.deconvolution import run_deconvolution
from pipeline.cell_communication.communication import run_cell_communication
from pipeline.visualization.plots import run_visualization

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


def parse_args():
    p = argparse.ArgumentParser(
        description="Spatial Transcriptomics Pipeline — 10x Visium"
    )
    p.add_argument("--config",  default="config/config.yaml",
                   help="Path to config YAML")
    p.add_argument("--outdir",  default="results",
                   help="Output directory")
    p.add_argument("--steps",   default="all",
                   help="Comma-separated steps: qc,cluster,spatial_de,deconv,comm,viz")
    p.add_argument("--sample",  default=None,
                   help="Run on a single sample ID (default: all)")
    p.add_argument("--threads", type=int, default=8)
    return p.parse_args()


STEPS_ORDER = ["qc", "cluster", "spatial_de", "deconv", "comm", "viz"]

STEP_FNS = {
    "qc":         run_qc,
    "cluster":    run_clustering,
    "spatial_de": run_spatial_de,
    "deconv":     run_deconvolution,
    "comm":       run_cell_communication,
    "viz":        run_visualization,
}

STEP_LABELS = {
    "qc":         "QC & Preprocessing",
    "cluster":    "Spatially-aware Clustering",
    "spatial_de": "SpatialDE Variable Gene Detection",
    "deconv":     "Cell Type Deconvolution",
    "comm":       "Cell-Cell Communication",
    "viz":        "Publication-ready Visualizations",
}


def main():
    args = parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    steps = STEPS_ORDER if args.steps == "all" else args.steps.split(",")

    samples = cfg["samples"]
    if args.sample:
        samples = [s for s in samples if s["id"] == args.sample]
        if not samples:
            log.error(f"Sample '{args.sample}' not found in config.")
            sys.exit(1)

    print("\n╔══════════════════════════════════════════════════════╗")
    print("║   Spatial Transcriptomics Pipeline — 10x Visium     ║")
    print("║   Tumor Microenvironment Analysis                   ║")
    print("╚══════════════════════════════════════════════════════╝\n")
    log.info(f"Samples: {[s['id'] for s in samples]}")
    log.info(f"Steps:   {steps}")

    for step in steps:
        if step not in STEP_FNS:
            log.warning(f"Unknown step '{step}' — skipping.")
            continue

        log.info(f"── STEP: {STEP_LABELS[step]} {'─'*20}")
        STEP_FNS[step](
            samples  = samples,
            cfg      = cfg,
            outdir   = outdir / step,
            threads  = args.threads,
        )

    log.info(f"✓ Pipeline complete. Results in: {outdir}/")


if __name__ == "__main__":
    main()
