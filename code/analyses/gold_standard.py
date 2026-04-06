"""Analysis 1: Gold-standard databases (LUDB + QTDB) — full manual annotation."""

from scipy import stats

from ..data_loader import load_beats, apply_quality_filters, subject_level_aggregate, bootstrap_mode, stars
from ..config import ONE_OVER_E


def run(log=print):
    results = {}

    # --- LUDB ---
    log("  Loading LUDB...")
    ludb = apply_quality_filters(load_beats("ludb_beats.csv"), log=log)
    ludb_subj = subject_level_aggregate(ludb)

    hc = ludb_subj[ludb_subj["group"] == "healthy"]
    pa = ludb_subj[ludb_subj["group"] == "pathological"]
    log(f"\n  LUDB: {len(hc)} healthy, {len(pa)} pathological subjects")

    hc_mode, hc_ci = bootstrap_mode(hc["CDC_median"].values)
    pa_mode, pa_ci = bootstrap_mode(pa["CDC_median"].values)
    log(f"  Healthy mode = {hc_mode:.4f} [{hc_ci[0]:.4f}, {hc_ci[1]:.4f}]")
    log(f"  Pathological mode = {pa_mode:.4f} [{pa_ci[0]:.4f}, {pa_ci[1]:.4f}]")
    log(f"  1/e = {ONE_OVER_E:.4f}, delta(healthy) = {hc_mode - ONE_OVER_E:+.4f}")

    u_stat, p_val = stats.mannwhitneyu(hc["CDC_median"], pa["CDC_median"], alternative="two-sided")
    log(f"  Wilcoxon rank-sum: U={u_stat:.0f}, p={p_val:.4e} {stars(p_val)}")

    t_hc, p_hc = stats.ttest_1samp(hc["CDC_median"], ONE_OVER_E)
    t_pa, p_pa = stats.ttest_1samp(pa["CDC_median"], ONE_OVER_E)
    log(f"  Healthy vs 1/e: t={t_hc:.2f}, p={p_hc:.4e} {stars(p_hc)}")
    log(f"  Pathological vs 1/e: t={t_pa:.2f}, p={p_pa:.4e} {stars(p_pa)}")

    results["ludb"] = dict(
        subjects=ludb_subj, n_healthy=len(hc), n_path=len(pa),
        mode_hc=hc_mode, ci_hc=hc_ci, mode_pa=pa_mode, ci_pa=pa_ci,
    )

    # --- QTDB ---
    log("\n  Loading QTDB...")
    qtdb = apply_quality_filters(load_beats("qtdb_beats.csv"), log=log)
    qtdb_subj = subject_level_aggregate(qtdb)

    groups = {}
    for g in qtdb_subj["group"].unique():
        gdf = qtdb_subj[qtdb_subj["group"] == g]
        mode, ci = bootstrap_mode(gdf["CDC_median"].values)
        log(f"  QTDB {g}: n={len(gdf)}, mode={mode:.4f} [{ci[0]:.4f}, {ci[1]:.4f}]")
        groups[g] = gdf

    if len(groups) >= 2:
        kw_stat, kw_p = stats.kruskal(*[groups[g]["CDC_median"].values for g in groups])
        log(f"  Kruskal-Wallis: H={kw_stat:.2f}, p={kw_p:.4e} {stars(kw_p)}")

    results["qtdb"] = dict(subjects=qtdb_subj, groups={g: len(groups[g]) for g in groups})
    return results
