"""Analysis 6: Age-stratified mortality (CODE-15%)."""

import numpy as np
import pandas as pd

from ..data_loader import load_beats, apply_quality_filters, subject_level_aggregate
from ..config import DATA_DIR, ONE_OVER_E


def run(code15_subj=None, log=print):
    if code15_subj is None:
        log("  Loading CODE-15% beats...")
        beats = apply_quality_filters(load_beats("code15_beats.csv"), log=log)
        code15_subj = subject_level_aggregate(beats)

    log("  Loading CODE-15% exams...")
    exams = pd.read_csv(DATA_DIR / "code15_exams.csv", low_memory=False)
    exams["patient_id"] = pd.to_numeric(exams["patient_id"], errors="coerce").astype("Int64")
    exams_first = exams.sort_values("exam_id").groupby("patient_id").first().reset_index()

    subj = code15_subj.copy()
    subj["patient_id"] = pd.to_numeric(subj["record_id"], errors="coerce").astype("Int64")

    merged = subj.merge(
        exams_first[["patient_id", "death", "timey", "age", "is_male"]],
        on="patient_id", how="inner", suffixes=("", "_exam"),
    )
    merged["death"] = pd.to_numeric(merged["death"], errors="coerce")
    merged["timey"] = pd.to_numeric(merged["timey"], errors="coerce")
    merged = merged.dropna(subset=["death", "timey", "CDC_median"])

    log(f"  Merged: {len(merged)} patients with mortality data")

    merged["cdc_dev"] = np.abs(merged["CDC_median"] - ONE_OVER_E)
    merged["age_decade"] = (merged["age"] // 10).astype(int) * 10

    log(f"\n  {'Decade':<8} {'N':>6} {'T1_mort%':>10} {'T2_mort%':>10} {'T3_mort%':>10} {'Fold(T3/T1)':>12}")
    log(f"  {'-' * 60}")

    decade_results = {}
    overall_t1_deaths, overall_t1_n = 0, 0
    overall_t3_deaths, overall_t3_n = 0, 0

    for decade in sorted(merged["age_decade"].unique()):
        stratum = merged[merged["age_decade"] == decade]
        if len(stratum) < 30:
            continue

        t1_cut = stratum["cdc_dev"].quantile(1 / 3)
        t2_cut = stratum["cdc_dev"].quantile(2 / 3)

        t1 = stratum[stratum["cdc_dev"] <= t1_cut]
        t2 = stratum[(stratum["cdc_dev"] > t1_cut) & (stratum["cdc_dev"] <= t2_cut)]
        t3 = stratum[stratum["cdc_dev"] > t2_cut]

        m1 = t1["death"].mean() * 100 if len(t1) > 0 else np.nan
        m2 = t2["death"].mean() * 100 if len(t2) > 0 else np.nan
        m3 = t3["death"].mean() * 100 if len(t3) > 0 else np.nan
        fold = m3 / m1 if m1 > 0 else np.nan

        overall_t1_deaths += t1["death"].sum()
        overall_t1_n += len(t1)
        overall_t3_deaths += t3["death"].sum()
        overall_t3_n += len(t3)

        log(f"  {int(decade):<8} {len(stratum):>6} {m1:>10.1f} {m2:>10.1f} {m3:>10.1f} {fold:>12.2f}x")

        decade_results[int(decade)] = dict(
            n=len(stratum), mort_t1=m1, mort_t2=m2, mort_t3=m3,
            n_t1=len(t1), n_t2=len(t2), n_t3=len(t3), fold=fold,
        )

    if overall_t1_n > 0 and overall_t3_n > 0:
        overall_fold = (overall_t3_deaths / overall_t3_n) / (overall_t1_deaths / overall_t1_n)
        log(f"\n  Overall T3/T1 mortality fold-change: {overall_fold:.2f}x")

    return merged, decade_results
