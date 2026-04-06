"""Analysis 5: CODE-15% standalone."""

from scipy import stats

from ..data_loader import load_beats, apply_quality_filters, subject_level_aggregate, bootstrap_mode, stars
from ..config import ONE_OVER_E


def run(log=print):
    log("  Loading CODE-15% beats...")
    beats = apply_quality_filters(load_beats("code15_beats.csv"), log=log)
    subj = subject_level_aggregate(beats)

    cn = subj[subj["group"] == "healthy"]
    pa = subj[subj["group"] == "pathological"]
    log(f"  Clinically Normal: n={len(cn)}")
    log(f"  Pathological: n={len(pa)}")

    cn_mode, cn_ci = bootstrap_mode(cn["CDC_median"].values)
    pa_mode, pa_ci = bootstrap_mode(pa["CDC_median"].values)
    log(f"  CN mode = {cn_mode:.4f} [{cn_ci[0]:.4f}, {cn_ci[1]:.4f}]")
    log(f"  Pathological mode = {pa_mode:.4f} [{pa_ci[0]:.4f}, {pa_ci[1]:.4f}]")
    log(f"  1/e = {ONE_OVER_E:.4f}")

    u, p = stats.mannwhitneyu(cn["CDC_median"], pa["CDC_median"], alternative="two-sided")
    log(f"  Wilcoxon: U={u:.0f}, p={p:.2e}")

    for label, data in [("CN", cn), ("Pathological", pa)]:
        valid = data.dropna(subset=["age"])
        if len(valid) > 20:
            slope, _, r, p_age, _ = stats.linregress(valid["age"], valid["CDC_median"])
            log(f"  {label} CDC~Age: slope={slope:.6f}/yr, r={r:.4f}, p={p_age:.2e}")

    return subj
