"""Analysis 4: CDC vs heart rate independence."""

from scipy import stats

from ..data_loader import stars
from ..config import ONE_OVER_E


def run(pooled, log=print):
    hc = pooled[pooled["clinical_group"] == "HealthyControl"].copy()
    log(f"  Healthy controls: n={len(hc)}")

    valid = hc.dropna(subset=["HR_median", "CDC_median"])
    slope, intercept, r, p, se = stats.linregress(valid["CDC_median"], valid["HR_median"])
    log(f"  HR ~ CDC regression: slope={slope:.2f}, intercept={intercept:.2f}, "
        f"R²={r**2:.4f}, p={p:.2e}")

    predicted_hr = slope * ONE_OVER_E + intercept
    log(f"  Predicted healthy HR at CDC=1/e: {predicted_hr:.1f} bpm")

    r_corr, p_corr = stats.pearsonr(valid["CDC_median"], valid["HR_median"])
    log(f"  CDC-HR correlation: r={r_corr:.4f}, p={p_corr:.2e}")

    return dict(slope=slope, intercept=intercept, r=r, predicted_hr=predicted_hr,
                healthy=valid)
