"""Data loading, quality filtering, and subject-level aggregation.

Replicates the MATLAB utilities:
  - apply_quality_filters.m
  - bootstrap_mode.m
  - build_beats_table.m / subject-level median aggregation
"""

import hashlib
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

from .config import (
    DATA_DIR, EXPECTED_MD5, ONE_OVER_E,
    RR_MIN, RR_MAX, RT_MIN, RT_MAX, CDC_MIN, CDC_MAX,
    MIN_BEATS_PER_SUBJECT,
    ZENODO_BASE_URL, ZENODO_CSVS,
)


# ------------------------------------------------------------------
# Data integrity
# ------------------------------------------------------------------

def md5_file(path):
    """Compute MD5 hex digest of a file."""
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_data_integrity(log=print):
    """Check MD5 checksums of all preprocessed CSVs. Returns True if all pass."""
    all_ok = True
    for fname, expected in EXPECTED_MD5.items():
        path = DATA_DIR / fname
        if not path.exists():
            log(f"  MISSING: {fname}")
            all_ok = False
            continue
        actual = md5_file(path)
        ok = actual == expected
        log(f"  [{'PASS' if ok else 'FAIL'}] {fname}")
        if not ok:
            log(f"    expected {expected}, got {actual}")
            all_ok = False
    return all_ok


def download_missing_csvs(log=print):
    """Download large CSVs from Zenodo if missing locally."""
    import urllib.request

    missing = [f for f in ZENODO_CSVS if not (DATA_DIR / f).exists()]
    if not missing:
        log("  All Zenodo CSVs already present.")
        return True

    log(f"  {len(missing)} file(s) missing — downloading from Zenodo...")
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    for fname in missing:
        url = f"{ZENODO_BASE_URL}/{fname}?download=1"
        dest = DATA_DIR / fname
        log(f"    Downloading {fname}...")
        try:
            urllib.request.urlretrieve(url, dest)
            actual_md5 = md5_file(dest)
            if actual_md5 != EXPECTED_MD5[fname]:
                log(f"    WARNING: MD5 mismatch for {fname}")
            else:
                log(f"    OK ({dest.stat().st_size / 1e6:.1f} MB)")
        except Exception as e:
            log(f"    FAILED: {e}")
            return False
    return True


# ------------------------------------------------------------------
# Loading
# ------------------------------------------------------------------

def load_beats(filename):
    """Load a beats CSV and coerce numeric columns."""
    path = DATA_DIR / filename
    df = pd.read_csv(path, low_memory=False)
    for col in ("r_sample", "t_end_sample", "next_r_sample", "fs", "age", "beat_idx"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    # PTB-XL encodes ages >= 90 as 300; recode to 90 (matches manuscript Methods)
    if "age" in df.columns:
        df.loc[df["age"] >= 300, "age"] = 90
    return df


# ------------------------------------------------------------------
# Quality filtering (replicates apply_quality_filters.m)
# ------------------------------------------------------------------

def apply_quality_filters(df, *, verbose=True, log=print):
    """Apply beat-level and subject-level quality filters.

    Returns a filtered DataFrame with added RR, RT, CDC, HR columns.
    """
    fs = df["fs"].values
    rr = (df["next_r_sample"].values - df["r_sample"].values) / fs
    rt = (df["t_end_sample"].values - df["r_sample"].values) / fs
    cdc = rt / rr

    mask = (
        (rr > 0) & (rt > 0) & (rt < rr)
        & (rr > RR_MIN) & (rr < RR_MAX)
        & (rt > RT_MIN) & (rt < RT_MAX)
        & (cdc > CDC_MIN) & (cdc < CDC_MAX)
    )

    n_total = len(df)
    n_valid = int(mask.sum())

    if verbose:
        log(f"  Quality filter: {n_valid}/{n_total} beats pass ({100 * n_valid / n_total:.1f}%)")

    out = df[mask].copy()
    out["RR"] = rr[mask]
    out["RT"] = rt[mask]
    out["CDC"] = cdc[mask]
    out["HR"] = 60.0 / rr[mask]

    # Subject-level filter
    if MIN_BEATS_PER_SUBJECT > 0 and "unique_subject_id" in out.columns:
        counts = out.groupby("unique_subject_id").size()
        valid_subj = counts[counts >= MIN_BEATS_PER_SUBJECT].index
        n_excl = len(counts) - len(valid_subj)
        if verbose and n_excl > 0:
            log(f"    Subjects with <{MIN_BEATS_PER_SUBJECT} valid beats: "
                f"{n_excl}/{len(counts)} excluded")
        out = out[out["unique_subject_id"].isin(valid_subj)]

    return out


# ------------------------------------------------------------------
# Subject-level aggregation
# ------------------------------------------------------------------

def subject_level_aggregate(df):
    """Aggregate beat-level data to one row per subject (median CDC, etc.)."""
    agg = df.groupby("unique_subject_id").agg(
        CDC_median=("CDC", "median"),
        HR_median=("RR", lambda x: 60.0 / np.median(x)),
        RR_median=("RR", "median"),
        RT_median=("RT", "median"),
        n_beats=("CDC", "count"),
        age=("age", "first"),
        sex=("sex", "first"),
        group=("group", "first"),
        source_database=("source_database", "first"),
        record_id=("record_id", "first"),
    ).reset_index()
    agg["DCDC"] = agg["CDC_median"] - ONE_OVER_E
    agg["CDC_dev"] = np.abs(agg["DCDC"])
    agg["RT_ms"] = agg["RT_median"] * 1000
    agg["RR_ms"] = agg["RR_median"] * 1000
    agg["diastole_ms"] = (agg["RR_median"] - agg["RT_median"]) * 1000
    return agg


# ------------------------------------------------------------------
# Mode estimation (replicates bootstrap_mode.m)
# ------------------------------------------------------------------

def bootstrap_mode(data, n_boot=5000, seed=42):
    """KDE mode with bootstrap 95 % CI. Matches MATLAB bootstrap_mode.m."""
    data = np.asarray(data, dtype=float)
    data = data[~np.isnan(data)]
    if len(data) < 3:
        return np.nan, (np.nan, np.nan)

    # Point estimate
    kde = gaussian_kde(data)
    x_grid = np.linspace(data.min(), data.max(), 300)
    mode_est = x_grid[np.argmax(kde(x_grid))]

    # Bootstrap
    rng = np.random.RandomState(seed)
    n = min(len(data), 5000)
    d = data if len(data) <= n else rng.choice(data, n, replace=False)

    boot_modes = np.zeros(n_boot)
    for b in range(n_boot):
        bs = rng.choice(d, len(d), replace=True)
        try:
            bkde = gaussian_kde(bs)
            bx = np.linspace(bs.min(), bs.max(), 200)
            boot_modes[b] = bx[np.argmax(bkde(bx))]
        except Exception:
            boot_modes[b] = np.nan

    ci = (float(np.nanpercentile(boot_modes, 2.5)),
          float(np.nanpercentile(boot_modes, 97.5)))
    return mode_est, ci


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

def stars(p):
    """Significance stars."""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"
