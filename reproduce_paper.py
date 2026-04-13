#!/usr/bin/env python3
"""
reproduce_paper.py
==================
Independent Python reproduction of:

    Froese, T. et al. (in prep.). The healthy heart converges on a
    thermodynamic optimum (1/e) that predicts survival. Nature Aging.

One command reproduces all 7 statistical analyses and all 9 manuscript
figures from the preprocessed Zenodo beat-level CSVs.

Usage
-----
    pip install -r requirements.txt
    python reproduce_paper.py

Outputs
-------
    results/reproduction_log.txt   — full analysis log
    results/figures/*.png, *.pdf   — all main + supplementary figures
    results/hierarchical_subject_data.csv — pooled subject-level data
"""

import sys
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# Ensure the project root is on sys.path
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from code.config import DATA_DIR, RESULTS_DIR, FIGURES_DIR, ONE_OVER_E
from code.data_loader import verify_data_integrity, download_missing_csvs

# Analysis modules
from code.analyses import gold_standard, large_scale, hierarchical, cdc_vs_hr, code15, mortality, survival

# Figure modules
from code.figures import (
    fig1_aging, fig2_mortality,
    si_fig1_hr, si_fig2_systole,
    si_fig4_gold, si_fig5_large, si_fig6_code15,
    si_fig7_km, si_fig8_sex,
)

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
LOG_FILE = RESULTS_DIR / "reproduction_log.txt"
_log_fh = None


def log(msg):
    print(msg, flush=True)
    if _log_fh:
        _log_fh.write(msg + "\n")
        _log_fh.flush()


# ---------------------------------------------------------------------------
# Pooled dataset builder (used by multiple figures)
# ---------------------------------------------------------------------------

def _build_pooled_for_figures(pooled):
    """Tag pooled subjects with 'hgroup' for figure plotting.

    Uses the same 5-database pool from the hierarchical analysis.
    Mapping: clinical_group -> hgroup (identical names, just renamed
    column for backward compat with figure scripts).
    """
    pooled_fig = pooled.copy()
    pooled_fig["hgroup"] = pooled_fig["clinical_group"]
    pooled_fig = pooled_fig.dropna(subset=["age"])
    pooled_fig = pooled_fig[pooled_fig["age"] < 120]
    return pooled_fig


def _build_code15_merged():
    """Load CODE-15% subjects merged with mortality data for figures."""
    import pandas as pd
    from code.data_loader import load_beats, apply_quality_filters, subject_level_aggregate
    import numpy as np

    df = apply_quality_filters(load_beats("code15_beats.csv"), verbose=False)
    subj = subject_level_aggregate(df)
    subj = subj[subj["age"] < 120]

    exams = pd.read_csv(DATA_DIR / "code15_exams.csv", low_memory=False)
    subj_c = subj.copy()
    subj_c["patient_id"] = subj_c["unique_subject_id"].str.replace("CODE15_", "", regex=False).astype(str)
    exams["patient_id"] = exams["patient_id"].astype(str)
    merged = subj_c.merge(
        exams[["patient_id", "death", "timey", "is_male"]].drop_duplicates("patient_id"),
        on="patient_id", how="inner",
    )
    merged = merged.dropna(subset=["death", "timey"])
    merged["CDC_dev"] = np.abs(merged["CDC_median"] - ONE_OVER_E)
    return subj, merged


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    global _log_fh

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    _log_fh = open(LOG_FILE, "w", encoding="utf-8")

    t0 = time.time()

    log("=" * 60)
    log("CARDIAC DUTY CYCLE — INDEPENDENT PYTHON REPRODUCTION")
    log(f"Date: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"1/e = {ONE_OVER_E:.10f}")
    log("=" * 60)

    # ------------------------------------------------------------------
    # Step 0: Ensure data is available
    # ------------------------------------------------------------------
    log("\n" + "=" * 60)
    log("STEP 0: DATA INTEGRITY")
    log("=" * 60 + "\n")

    download_missing_csvs(log=log)
    data_ok = verify_data_integrity(log=log)
    if not data_ok:
        log("\nWARNING: Some files missing or have MD5 mismatches.")
        log("Download from: https://doi.org/10.5281/zenodo.19246123")

    # ------------------------------------------------------------------
    # Step 1: Statistical analyses (1–7)
    # ------------------------------------------------------------------
    log("\n" + "=" * 60)
    log("STEP 1: STATISTICAL ANALYSES")
    log("=" * 60)

    log("\n[1/7] Gold-standard analysis (LUDB, QTDB)...")
    gold_standard.run(log=log)

    log("\n[2/7] Large-scale analysis (Fantasia, AA, PTB, PTB-XL)...")
    large_scale.run(log=log)

    log("\n[3/7] Hierarchical linear model (all databases pooled)...")
    pooled = hierarchical.run(log=log)

    log("\n[4/7] CDC vs heart rate independence...")
    cdc_vs_hr.run(pooled, log=log)

    log("\n[5/7] CODE-15% standalone...")
    code15_subj = code15.run(log=log)

    log("\n[6/7] Age-stratified mortality...")
    merged, decade_results = mortality.run(code15_subj, log=log)

    log("\n[7/7] Survival curves + sex-stratified...")
    survival.run(merged, log=log)

    # ------------------------------------------------------------------
    # Step 2: Generate all figures
    # ------------------------------------------------------------------
    log("\n" + "=" * 60)
    log("STEP 2: FIGURES")
    log("=" * 60 + "\n")

    log("  Building pooled dataset for figures...")
    pooled_fig = _build_pooled_for_figures(pooled)
    code15_subj_fig, merged_fig = _build_code15_merged()

    fig1_aging.plot(pooled_fig, log=log)
    fig2_mortality.plot(merged_fig, decade_results, log=log)
    si_fig1_hr.plot(pooled_fig, merged_fig, log=log)
    si_fig2_systole.plot(pooled_fig, log=log)
    si_fig4_gold.plot(log=log)
    si_fig5_large.plot(log=log)
    si_fig6_code15.plot(code15_subj_fig, log=log)
    si_fig7_km.plot(merged_fig, log=log)
    si_fig8_sex.plot(pooled_fig, log=log)

    # ------------------------------------------------------------------
    # Done
    # ------------------------------------------------------------------
    elapsed = time.time() - t0
    log(f"\n{'=' * 60}")
    log(f"REPRODUCTION COMPLETE in {elapsed:.0f} seconds")
    log(f"Results:  {RESULTS_DIR}")
    log(f"Figures:  {FIGURES_DIR}")
    log(f"Log:      {LOG_FILE}")
    log(f"{'=' * 60}")

    _log_fh.close()


if __name__ == "__main__":
    main()
