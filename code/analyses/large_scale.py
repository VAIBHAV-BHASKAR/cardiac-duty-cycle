"""Analysis 2: Large-scale databases (Fantasia, Autonomic Aging, PTB, PTB-XL)."""

from scipy import stats

from ..data_loader import load_beats, apply_quality_filters, subject_level_aggregate, bootstrap_mode, stars
from ..config import ONE_OVER_E


def run(log=print):
    results = {}

    for name, fname in [
        ("Fantasia", "fantasia_beats.csv"),
        ("Autonomic Aging", "autonomic_aging_beats.csv"),
        ("PTB", "ptb_beats.csv"),
        ("PTB-XL", "ptbxl_beats.csv"),
    ]:
        log(f"  Loading {name}...")
        df = apply_quality_filters(load_beats(fname), log=log)
        subj = subject_level_aggregate(df)

        log(f"  {name}: {len(subj)} subjects, groups: {list(subj['group'].unique())}")

        for g in subj["group"].unique():
            gdata = subj[subj["group"] == g]
            mode, ci = bootstrap_mode(gdata["CDC_median"].values)
            t_stat, p_val = stats.ttest_1samp(gdata["CDC_median"], ONE_OVER_E)
            log(f"    {g}: n={len(gdata)}, mode={mode:.4f} [{ci[0]:.4f}, {ci[1]:.4f}], "
                f"vs 1/e: t={t_stat:.2f}, p={p_val:.4e} {stars(p_val)}")

        if len(subj["group"].unique()) >= 2:
            hc = subj[subj["group"] == "healthy"]
            pa = subj[subj["group"] == "pathological"]
            if len(hc) > 0 and len(pa) > 0:
                u, p = stats.mannwhitneyu(hc["CDC_median"], pa["CDC_median"], alternative="two-sided")
                r_bs = 1 - 2 * u / (len(hc) * len(pa))
                log(f"    Wilcoxon: U={u:.0f}, p={p:.4e}, rank-biserial r={r_bs:.3f}")

        healthy = subj[subj["group"] == "healthy"]
        if len(healthy) > 20:
            valid_age = healthy.dropna(subset=["age"])
            if len(valid_age) > 20:
                r_age, p_age = stats.pearsonr(valid_age["age"], valid_age["CDC_median"])
                log(f"    CDC-Age correlation: r={r_age:.4f}, p={p_age:.4e} {stars(p_age)}")

        results[name] = dict(subjects=subj, n=len(subj))

    return results
