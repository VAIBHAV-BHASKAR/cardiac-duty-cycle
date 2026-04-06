"""SI Figure 7: Kaplan-Meier survival curves."""

import numpy as np
import pandas as pd

from ..config import C_T1, C_T2, C_T3, ONE_OVER_E
from ._helpers import save_fig, new_figure


def _kaplan_meier(time, event):
    """Simple Kaplan-Meier estimator."""
    df = pd.DataFrame({"time": time, "event": event}).sort_values("time")
    times = sorted(df["time"].unique())

    km_times = [0.0]
    km_surv = [1.0]
    n_at_risk = len(df)
    surv = 1.0

    for t in times:
        d = int(((df["time"] == t) & (df["event"] == 1)).sum())
        c = int(((df["time"] == t) & (df["event"] == 0)).sum())
        if n_at_risk > 0 and d > 0:
            surv *= 1 - d / n_at_risk
        km_times.append(t)
        km_surv.append(surv)
        n_at_risk -= d + c

    return np.array(km_times), np.array(km_surv)


def plot(merged, log=print):
    log("  Plotting SI Fig 7 — Kaplan-Meier survival...")

    merged = merged.copy()
    merged["cdc_dev"] = np.abs(merged["CDC_median"] - ONE_OVER_E)
    t1_cut = merged["cdc_dev"].quantile(1 / 3)
    t2_cut = merged["cdc_dev"].quantile(2 / 3)
    merged["tertile"] = pd.cut(merged["cdc_dev"], bins=[-np.inf, t1_cut, t2_cut, np.inf],
                                labels=["Near 1/e", "Moderate", "Far from 1/e"])

    fig, axes = new_figure(3, 1, figsize=(5, 10))

    sex_groups = [("Overall", merged)]
    if "is_male" in merged.columns:
        sex_groups += [
            ("Female", merged[merged["is_male"] == 0]),
            ("Male", merged[merged["is_male"] == 1]),
        ]

    panel_labels = "abc"
    for ax_idx, (sex_label, data) in enumerate(sex_groups):
        ax = axes[ax_idx]
        for tert, color in [("Near 1/e", C_T1), ("Moderate", C_T2), ("Far from 1/e", C_T3)]:
            t = data[data["tertile"] == tert]
            if len(t) == 0:
                continue
            km_t, km_s = _kaplan_meier(t["timey"].values, t["death"].values)
            ax.step(km_t, km_s * 100, where="post", color=color, linewidth=1.5,
                    label=f"{tert} (N={len(t):,})")

        ax.set_xlabel("Follow-up (years)", fontsize=9)
        ax.set_ylabel("Survival (%)", fontsize=9)
        ax.set_title(f"({panel_labels[ax_idx]}) {sex_label}", fontsize=9, fontweight="bold", loc="left")
        ax.legend(fontsize=6, loc="lower left")
        ax.set_ylim(88, 100.5)
        ax.tick_params(labelsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    save_fig(fig, "SI_Fig7_KM_survival", log=log)
