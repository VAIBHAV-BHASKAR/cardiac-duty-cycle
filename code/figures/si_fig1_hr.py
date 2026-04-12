"""SI Figure 1: CDC vs heart rate + hidden mortality risk.

Matches MATLAB plot_SI_Fig1.m: title updated to "4 databases", HR rounded.
"""

import numpy as np
import pandas as pd
from scipy import stats

from ..config import C_HC, C_T1, C_T2, C_T3, ONE_OVER_E
from ._helpers import save_fig, new_figure


def plot(pooled, merged, log=print):
    log("  Plotting SI Fig 1 — CDC vs HR...")
    fig, axes = new_figure(1, 2, figsize=(7.2, 3.5))

    # Panel a: Healthy CDC-HR regression
    ax = axes[0]
    healthy = pooled[pooled["hgroup"] == "HealthyControl"].dropna(
        subset=["HR_median", "CDC_median"])
    ax.scatter(healthy["CDC_median"], healthy["HR_median"], s=2, alpha=0.15,
               color=C_HC, rasterized=True)

    slope, intercept, r, p, se = stats.linregress(
        healthy["CDC_median"], healthy["HR_median"])
    x_line = np.linspace(0.25, 0.55, 100)
    ax.plot(x_line, slope * x_line + intercept, color=C_HC, linewidth=2)

    predicted_hr = slope * ONE_OVER_E + intercept
    ax.axvline(ONE_OVER_E, color="gray", linestyle="--", linewidth=0.7)
    ax.plot(ONE_OVER_E, predicted_hr, "ko", markersize=6, zorder=5)
    ax.annotate(f"1/$e$ → HR = {predicted_hr:.0f} bpm",
                xy=(ONE_OVER_E, predicted_hr),
                xytext=(0.42, predicted_hr + 15), fontsize=7,
                arrowprops=dict(arrowstyle="->", color="black", lw=0.8))

    ax.set_xlabel("Median CDC (RT/RR)", fontsize=9)
    ax.set_ylabel("Heart rate (bpm)", fontsize=9)
    ax.set_title("(a) 1/$e$ predicts healthy resting HR (4 databases)",
                 fontsize=9, fontweight="bold", loc="left")
    ax.text(0.45, 50, f"R = {r:.3f}\np < {max(p, 1e-300):.0e}",
            fontsize=7, color="#555")
    ax.tick_params(labelsize=8)

    # Panel b: Hidden risk in CN by HR stratum
    ax = axes[1]
    cn = merged[merged["group"] == "healthy"].copy()
    cn["HR_bin"] = pd.cut(cn["HR_median"], bins=np.arange(40, 110, 10))

    bar_data = []
    for hr_bin in sorted(cn["HR_bin"].dropna().unique()):
        d = cn[cn["HR_bin"] == hr_bin]
        if len(d) < 100:
            continue
        cuts = d["CDC_dev"].quantile([1 / 3, 2 / 3]).values
        for tname, mask, color in [
            ("Near 1/$e$", d["CDC_dev"] <= cuts[0], C_T1),
            ("Moderate", (d["CDC_dev"] > cuts[0]) & (d["CDC_dev"] <= cuts[1]), C_T2),
            ("Far from 1/$e$", d["CDC_dev"] > cuts[1], C_T3),
        ]:
            t = d[mask]
            if t["death"].sum() >= 5:
                bar_data.append(dict(
                    hr_bin=str(hr_bin), tertile=tname,
                    mort=100 * t["death"].mean(),
                    se=100 * np.sqrt(t["death"].mean() *
                                     (1 - t["death"].mean()) / len(t)),
                    color=color))

    if bar_data:
        bd = pd.DataFrame(bar_data)
        hr_bins = bd["hr_bin"].unique()
        x = np.arange(len(hr_bins))
        width = 0.25
        for i, (tert, color) in enumerate([
            ("Near 1/$e$", C_T1), ("Moderate", C_T2),
            ("Far from 1/$e$", C_T3),
        ]):
            vals, errs = [], []
            for hb in hr_bins:
                row = bd[(bd["hr_bin"] == hb) & (bd["tertile"] == tert)]
                vals.append(row["mort"].values[0] if len(row) > 0 else 0)
                errs.append(row["se"].values[0] if len(row) > 0 else 0)
            ax.bar(x + i * width, vals, width, yerr=errs, label=tert,
                   color=color, capsize=2, edgecolor="white", linewidth=0.5,
                   error_kw={"linewidth": 0.8})
        ax.set_xticks(x + width)
        ax.set_xticklabels([str(hb) for hb in hr_bins], fontsize=6,
                           rotation=45)

    ax.set_xlabel("Heart rate stratum (bpm)", fontsize=9)
    ax.set_ylabel("All-cause mortality (%)", fontsize=9)
    ax.set_title("(b) CDC risk beyond HR screening", fontsize=9,
                 fontweight="bold", loc="left")
    ax.legend(fontsize=6)
    ax.tick_params(labelsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    save_fig(fig, "SI_Fig1_CDC_vs_HR", log=log)
