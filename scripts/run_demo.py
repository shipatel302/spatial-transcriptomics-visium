#!/usr/bin/env python3
"""
scripts/run_demo.py
Run the full spatial transcriptomics pipeline on simulated Visium data.
No real SpaceRanger output required.

Usage:
    python scripts/run_demo.py
"""

import sys
import logging
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

sys.path.insert(0, str(Path(__file__).parent.parent))

import yaml
from pipeline.qc.visium_qc import run_qc
from pipeline.clustering.spatial_clustering import run_clustering
from pipeline.spatial_de.spatial_de import run_spatial_de
from pipeline.deconvolution.deconvolution import run_deconvolution
from pipeline.cell_communication.communication import run_cell_communication
from pipeline.visualization.plots import run_visualization

print("""
╔══════════════════════════════════════════════════════╗
║   Spatial Transcriptomics Pipeline — Demo Run       ║
║   Simulated 10x Visium / HNSCC TME                  ║
╚══════════════════════════════════════════════════════╝
""")

with open("config/config.yaml") as f:
    cfg = yaml.safe_load(f)

OUTDIR = Path("results/demo")
SAMPLES = cfg["samples"]

steps = [
    ("QC & Preprocessing",            run_qc,               "qc"),
    ("Spatially-aware Clustering",     run_clustering,       "cluster"),
    ("Spatially Variable Genes",       run_spatial_de,       "spatial_de"),
    ("Cell Type Deconvolution",        run_deconvolution,    "deconv"),
    ("Cell-Cell Communication",        run_cell_communication, "comm"),
    ("Publication-ready Figures",      run_visualization,    "viz"),
]

for label, fn, subdir in steps:
    log.info(f"── {label} {'─' * (45 - len(label))}")
    try:
        fn(
            samples = SAMPLES,
            cfg     = cfg,
            outdir  = OUTDIR / subdir,
            threads = 4,
        )
    except Exception as e:
        log.error(f"  Step failed: {e}")
        log.info("  Continuing to next step...")

log.info(f"✓ Demo complete. Results in: {OUTDIR}/")
log.info("")
log.info("Output structure:")
for p in sorted(OUTDIR.rglob("*.tsv"))[:10]:
    log.info(f"  {p.relative_to(OUTDIR)}")
