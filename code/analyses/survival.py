"""Analysis 7: Survival curves and sex-stratified analysis (CODE-15%)."""

import numpy as np
import pandas as pd
from scipy import stats

from ..data_loader import bootstrap_mode, stars
from ..config import ONE_OVER_E


def run(merged, log=print):
    merged = merged.dropna(subset=["death", "timey"]).copy()
    merged["cdc_dev"] = np.abs(merged["CDC_median"] - ONE_OVER_E)

    # Global tertiles
    t1_cut = merged["cdc_dev"].quantile(1 / 3)
    t2_cut = merged["cdc_dev"].quantile(2 / 3)

    t1 = merged[merged["cdc_dev"] <= t1_cut]
    t2 = merged[(merged["cdc_dev"] > t1_cut) & (merged["cdc_dev"] <= t2_cut)]
    t3 = merged[merged["cdc_dev"] > t2_cut]

    for label, grp in [("T1 (Near 1/e)", t1), ("T2 (Moderate)", t2), ("T3 (Far)", t3)]:
        n = len(grp)
        n_deaths = int(grp["death"].sum())
        mort = n_deaths / n * 100
        median_fu = grp["timey"].median()
        log(f"  {label}: n={n}, deaths={n_deaths} ({mort:.1f}%), median FU={median_fu:.1f} yr")

    contingency = np.array([
        [t1["death"].sum(), len(t1) - t1["death"].sum()],
        [t3["death"].sum(), len(t3) - t3["death"].sum()],
    ])
    chi2, p_lr = stats.chi2_contingency(contingency)[:2]
    log(f"\n  Chi-square (T1 vs T3): chi2={chi2:.2f}, p={p_lr:.2e}")

    # Sex-stratified
    if "is_male" in merged.columns:
        log("\n  Sex-stratified:")
        for sex_label, sex_val in [("Female", 0), ("Male", 1)]:
            sex_data = merged[merged["is_male"] == sex_val]
            if len(sex_data) < 100:
                continue
            t1s = sex_data[sex_data["cdc_dev"] <= t1_cut]
            t3s = sex_data[sex_data["cdc_dev"] > t2_cut]
            if len(t1s) > 0 and len(t3s) > 0:
                fold = (t3s["death"].mean()) / (t1s["death"].mean()) if t1s["death"].mean() > 0 else np.nan
                log(f"    {sex_label}: T1 mort={t1s['death'].mean()*100:.2f}%, "
                    f"T3 mort={t3s['death'].mean()*100:.2f}%, fold={fold:.2f}x")

    return dict(t1=t1, t2=t2, t3=t3, t1_cut=t1_cut, t2_cut=t2_cut, merged=merged)
