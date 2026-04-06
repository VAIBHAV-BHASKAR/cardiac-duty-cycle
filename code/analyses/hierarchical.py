"""Analysis 3: Hierarchical linear model — all databases pooled."""

import pandas as pd
from scipy import stats

from ..data_loader import load_beats, apply_quality_filters, subject_level_aggregate, bootstrap_mode, stars
from ..config import ONE_OVER_E, RESULTS_DIR


def _map_clinical_group(row):
    """Map (source_database, group) → hierarchical clinical group."""
    db = row["source_database"]
    g = row["group"]
    if db in ("LUDB", "PTB", "PTB Diagnostic"):
        return "HealthyControl" if g == "healthy" else "Pathological"
    if db in ("Fantasia", "Autonomic Aging", "AutonomicAging"):
        return "HealthyControl"
    if db in ("PTB-XL", "PTBXL"):
        return "ClinicallyNormal" if g == "healthy" else "Pathological"
    if db in ("CODE-15%", "CODE15", "CODE-15"):
        return "ClinicallyNormal" if g == "healthy" else "Pathological"
    if db == "QTDB":
        return "ClinicallyNormal" if g == "healthy" else "Pathological"
    return "Pathological"


def run(log=print):
    datasets = [
        ("ludb_beats.csv", "LUDB"),
        ("qtdb_beats.csv", "QTDB"),
        ("ptb_beats.csv", "PTB"),
        ("ptbxl_beats.csv", "PTB-XL"),
        ("fantasia_beats.csv", "Fantasia"),
        ("autonomic_aging_beats.csv", "AA"),
        ("code15_beats.csv", "CODE15"),
    ]

    all_subjects = []
    for fname, label in datasets:
        log(f"  Loading {label}...")
        df = apply_quality_filters(load_beats(fname), log=log)
        subj = subject_level_aggregate(df)
        subj["dataset"] = label
        all_subjects.append(subj)

    pooled = pd.concat(all_subjects, ignore_index=True)
    pooled = pooled.dropna(subset=["age"])
    pooled["clinical_group"] = pooled.apply(_map_clinical_group, axis=1)

    log(f"\n  Pooled dataset: {len(pooled)} subjects")
    for cg in ("HealthyControl", "ClinicallyNormal", "Pathological"):
        subset = pooled[pooled["clinical_group"] == cg]
        if len(subset) == 0:
            continue
        mode, ci = bootstrap_mode(subset["CDC_median"].values)
        t_stat, p_val = stats.ttest_1samp(subset["CDC_median"], ONE_OVER_E)
        log(f"  {cg}: n={len(subset)}, mode={mode:.4f} [{ci[0]:.4f}, {ci[1]:.4f}], "
            f"vs 1/e: p={p_val:.2e} {stars(p_val)}")

    log("\n  OLS CDC ~ Age by clinical group:")
    for cg in ("HealthyControl", "ClinicallyNormal", "Pathological"):
        subset = pooled[pooled["clinical_group"] == cg]
        if len(subset) > 20:
            slope, intercept, r, p, se = stats.linregress(subset["age"], subset["CDC_median"])
            log(f"    {cg}: slope={slope:.6f}/yr, r={r:.4f}, p={p:.2e}")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    pooled.to_csv(RESULTS_DIR / "hierarchical_subject_data.csv", index=False)
    return pooled
