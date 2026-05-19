"""SI Figure 4 (v1.3): Kaplan-Meier survival curves.

(Was SI Fig 7 in manuscript v1.0; renumbered to SI Fig 4 in v1.3.)

Matches MATLAB plot_SI_Fig4.m: taller figure, manual panel positions,
larger fonts, wider y-range, legend boxes on, KM line width 2.0.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2

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


def _logrank_p(time, event, group):
    """Multivariate log-rank test (k groups, chi-square with k-1 df).

    Standard formulation: at each event time, expected deaths per group
    follow the hypergeometric mean; variance the hypergeometric variance.
    """
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=int)
    group = np.asarray(group)

    keep = ~(np.isnan(time) | (event < 0))
    time, event, group = time[keep], event[keep], group[keep]

    labels = sorted(pd.unique(group).tolist())
    k = len(labels)
    if k < 2:
        return float("nan")

    event_times = np.unique(time[event == 1])
    oe = np.zeros(k)
    var = np.zeros((k, k))

    for t in event_times:
        at_risk = time >= t
        n_total = int(at_risk.sum())
        d_total = int(((time == t) & (event == 1)).sum())
        if n_total <= 1 or d_total == 0:
            continue
        n_per = np.array([int((at_risk & (group == g)).sum()) for g in labels])
        d_per = np.array([int(((time == t) & (event == 1) & (group == g)).sum())
                          for g in labels])
        e_per = n_per * d_total / n_total
        oe += d_per - e_per
        denom = (n_total ** 2) * (n_total - 1)
        for i in range(k):
            for j in range(k):
                if i == j:
                    var[i, i] += (n_per[i] * (n_total - n_per[i])
                                  * d_total * (n_total - d_total)) / denom
                else:
                    var[i, j] -= (n_per[i] * n_per[j]
                                  * d_total * (n_total - d_total)) / denom

    oe_red = oe[: k - 1]
    var_red = var[: k - 1, : k - 1]
    try:
        chi_sq = float(oe_red @ np.linalg.solve(var_red, oe_red))
    except np.linalg.LinAlgError:
        return float("nan")
    if chi_sq < 0 or not np.isfinite(chi_sq):
        return float("nan")
    return float(1.0 - chi2.cdf(chi_sq, k - 1))


def _format_comma(n):
    return f"{n:,}"


def plot(merged, log=print):
    log("  Plotting SI Fig 4 — Kaplan-Meier survival...")

    merged = merged.copy()
    merged["cdc_dev"] = np.abs(merged["CDC_median"] - ONE_OVER_E)
    t1_cut = merged["cdc_dev"].quantile(1 / 3)
    t2_cut = merged["cdc_dev"].quantile(2 / 3)
    merged["tertile"] = pd.cut(
        merged["cdc_dev"], bins=[-np.inf, t1_cut, t2_cut, np.inf],
        labels=["Near 1/e", "Moderate", "Far from 1/e"])

    # Wider figure to accommodate outside-right legends + log-rank annotation
    fig = plt.figure(figsize=(7.5, 11.0))  # taller: 28 cm

    # Narrower panel width leaves room for the right-side legend
    positions = [
        [0.10, 0.72, 0.62, 0.23],
        [0.10, 0.40, 0.62, 0.23],
        [0.10, 0.08, 0.62, 0.23],
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

        for tert, color in [("Near 1/e", C_T1), ("Moderate", C_T2),
                             ("Far from 1/e", C_T3)]:
            t = data[data["tertile"] == tert]
            if len(t) == 0:
                continue
            km_t, km_s = _kaplan_meier(t["timey"].values, t["death"].values)
            n_deaths = int(t["death"].sum())
            ax.step(km_t, km_s * 100, where="post", color=color,
                    linewidth=_KM_LW,
                    label=f"{tert}\nN = {_format_comma(len(t))}, "
                          f"{_format_comma(n_deaths)} deaths")

        # Panel title with comma-formatted aggregate n and deaths
        n_panel = len(data)
        d_panel = int(data["death"].sum())
        ax.set_xlabel("Follow-up (years)", fontsize=_AX)
        ax.set_ylabel("Survival (%)", fontsize=_AX)
        ax.set_title(
            f"({panel_labels[ax_idx]}) {sex_label} "
            f"(n = {_format_comma(n_panel)}; "
            f"{_format_comma(d_panel)} deaths)",
            fontsize=_TITLE, fontweight="bold", loc="left", pad=8)

        # Legend: outside the axes on the right, vertically centred
        leg = ax.legend(fontsize=_LEG,
                        bbox_to_anchor=(1.02, 0.5), loc="center left",
                        frameon=True, handlelength=1.5, borderaxespad=0.0)
        leg.get_frame().set(**_LEG_FRAME)

        # Y-axis: 88–100 captures the full T3 descent (Male panel needs this)
        ax.set_ylim(88, 100.5)
        ax.tick_params(labelsize=_TICK)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # Log-rank p-value across tertiles
        d_valid = data.dropna(subset=["tertile", "timey", "death"])
        p_val = _logrank_p(d_valid["timey"].values,
                           d_valid["death"].values.astype(int),
                           d_valid["tertile"].astype(str).values)
        if np.isfinite(p_val):
            p_str = ("Log-rank $p$ < 0.001"
                     if p_val < 1e-3
                     else f"Log-rank $p$ = {p_val:.3f}")
            ax.text(0.02, 0.06, p_str, transform=ax.transAxes,
                    fontsize=_LR_FS, fontweight="bold",
                    ha="left", va="bottom",
                    bbox=dict(facecolor="white", edgecolor=(0.5, 0.5, 0.5),
                              boxstyle="round,pad=0.3", alpha=0.85))

    save_fig(fig, "SI_Fig4_KM_survival", log=log)
