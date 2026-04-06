"""Central configuration for the cardiac-duty-cycle reproduction."""

from pathlib import Path
import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data" / "preprocessed"
RESULTS_DIR = ROOT / "results"
FIGURES_DIR = RESULTS_DIR / "figures"

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
ONE_OVER_E = 1.0 / np.e  # 0.36787944...

# ---------------------------------------------------------------------------
# Quality-filter defaults (match MATLAB apply_quality_filters.m exactly)
# ---------------------------------------------------------------------------
RR_MIN = 0.3    # seconds  (200 bpm)
RR_MAX = 3.0    # seconds  (20 bpm)
RT_MIN = 0.15   # seconds
RT_MAX = 0.6    # seconds
CDC_MIN = 0.2
CDC_MAX = 0.6
MIN_BEATS_PER_SUBJECT = 2

# ---------------------------------------------------------------------------
# MD5 checksums from Zenodo (doi:10.5281/zenodo.19246123)
# ---------------------------------------------------------------------------
EXPECTED_MD5 = {
    "ludb_beats.csv":              "d50dae990e5a0ff311c124db15658f88",
    "qtdb_beats.csv":              "86e02d7eaaea5e102da056b7502a77cd",
    "ptb_beats.csv":               "079763962f299196ac478d853b9cdf4b",
    "ptbxl_beats.csv":             "6966af4646266d7c1951c925a23f2823",
    "fantasia_beats.csv":          "c3096c7fc4524bb322ee3b330b4cb304",
    "autonomic_aging_beats.csv":   "b80795d3adbb329f36a6773a5b316290",
    "code15_beats.csv":            "20b17df8d87b5512c0257d37a202fc97",
    "code15_exams.csv":            "0107516d3f63864498fb77d15799cc95",
}

# Files small enough to commit to git (< 1 MB each)
COMMITTED_CSVS = {"ludb_beats.csv", "qtdb_beats.csv", "ptb_beats.csv"}

# Files that must be downloaded from Zenodo
ZENODO_CSVS = set(EXPECTED_MD5) - COMMITTED_CSVS

ZENODO_DOI = "10.5281/zenodo.19246123"
ZENODO_RECORD_ID = "19246123"
ZENODO_BASE_URL = f"https://zenodo.org/records/{ZENODO_RECORD_ID}/files"

# ---------------------------------------------------------------------------
# Colour palette — matches the MATLAB Nature Aging figures
# ---------------------------------------------------------------------------
C_HC   = "#2166ac"   # blue    — Healthy Control (verified volunteers)
C_CN   = "#4dac26"   # green   — Clinically Normal (hospital patients, normal ECG)
C_PATH = "#d6604d"   # red     — Pathological
C_SD   = "#b2182b"   # crimson — Sudden Death (QTDB only)
C_T1   = "#4dac26"   # green   — Near 1/e
C_T2   = "#f4a582"   # amber   — Moderate deviation
C_T3   = "#d6604d"   # red     — Far from 1/e
