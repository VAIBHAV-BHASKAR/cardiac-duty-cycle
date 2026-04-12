"""SI Figure 7: Kaplan-Meier survival curves.

Matches MATLAB plot_SI_Fig7.m: taller figure, manual panel positions,
larger fonts, wider y-range, legend boxes on, KM line width 2.0.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ..config import C_T1, C_T2, C_T3, ONE_OVER_E
from ._helpers import save_fig

_AX = 10
_TICK = 9
_TITLE = 11
_PANEL = 13
_LEG = 8
_KM_LW = 2.0
_LR_FS = 9

_LEG_FRAME = dict(facecolor="white", edgecolor=(0.7, 0.7, 0.7))


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


def _format_comma(n):
    return f"{n:,}"


def plot(merged, log=print):
    log("  Plotting SI Fig 7 — Kaplan-Meier survival...")

    merged = merged.copy()
    merged["cdc_dev"] = np.abs(merged["CDC_median"] - ONE_OVER_E)
    t1_cut = merged["cdc_dev"].quantile(1 / 3)
    t2_cut = merged["cdc_dev"].quantile(2 / 3)
    merged["tertile"] = pd.cut(
        merged["cdc_dev"], bins=[-np.inf, t1_cut, t2_cut, np.inf],
        labels=["Near 1/e", "Moderate", "Far from 1/e"])

    fig = plt.figure(figsize=(5, 11.0))  # taller: 28 cm

    positions = [
        [0.12, 0.72, 0.82, 0.23],
        [0.12, 0.40, 0.82, 0.23],
        [0.12, 0.08, 0.82, 0.23],
    ]

    sex_groups = [("Overall", merged)]
    if "is_male" in merged.columns:
        sex_groups += [
            ("Female", merged[merged["is_male"] == 0]),
            ("Male", merged[merged["is_male"] == 1]),
        ]

    panel_labels = "abc"
    for ax_idx, (sex_label, data) in enumerate(sex_groups):
        ax = fig.add_axes(positions[ax_idx])

        n_deaths_total = int(data["death"].sum())
        for tert, color in [("Near 1/e", C_T1), ("Moderate", C_T2),
                             ("Far from 1/e", C_T3)]:
            t = data[data["tertile"] == tert]
            if len(t) == 0:
                continue
            km_t, km_s = _kaplan_meier(t["timey"].values, t["death"].values)
            n_deaths = int(t["death"].sum())
            ax.step(km_t, km_s * 100, where="post", color=color,
                    linewidth=_KM_LW,
                    label=f"{tert} (N={_format_comma(len(t))}, "
                          f"{_format_comma(n_deaths)} deaths)")

        ax.set_xlabel("Follow-up (years)", fontsize=_AX)
        ax.set_ylabel("Survival (%)", fontsize=_AX)
        ax.set_title(
            f"({panel_labels[ax_idx]}) {sex_label}",
            fontsize=_TITLE, fontweight="bold", loc="left")
        leg = ax.legend(fontsize=_LEG, loc="lower left")
        leg.get_frame().set(**_LEG_FRAME)
        ax.set_ylim(88, 100.5)
        ax.tick_params(labelsize=_TICK)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    save_fig(fig, "SI_Fig7_KM_survival", log=log)
