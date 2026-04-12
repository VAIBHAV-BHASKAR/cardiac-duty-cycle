"""SI Figure 5: Large-scale CDC distributions (4 databases).

Matches MATLAB plot_SI_Fig5.m: legend boxes on, annotation text wrapped to
two lines, panel titles show age range, explicit axes ordering.
"""

import numpy as np
from scipy.stats import gaussian_kde

from ..config import C_HC, C_PATH
from ..data_loader import (load_beats, apply_quality_filters,
                           subject_level_aggregate, bootstrap_mode)
from ._helpers import save_fig, new_figure

_LEG_FRAME = dict(facecolor="white", edgecolor=(0.7, 0.7, 0.7))


def plot(log=print):
    log("  Plotting SI Fig 5 — Large-scale distributions...")
    fig, axes = new_figure(2, 2, figsize=(7.2, 6))
    axes_flat = [axes[0][0], axes[0][1], axes[1][0], axes[1][1]]

    datasets = [
        ("Fantasia", "fantasia_beats.csv",
         "Database R-peaks,\ntangent T-end"),
        ("Autonomic Aging", "autonomic_aging_beats.csv",
         "Fully automatic\n(Pan-Tompkins + tangent)"),
        ("PTB: ages 17\u201387", "ptb_beats.csv",
         "Manual T-end (5 referees),\nPan-Tompkins R-peak"),
        ("PTB-XL: ages 2\u201390", "ptbxl_beats.csv",
         "ECGDeli automatic\n(R-peak + T-end)"),
    ]

    for idx, (title, fname, ann_note) in enumerate(datasets):
        ax = axes_flat[idx]
        df = apply_quality_filters(load_beats(fname), verbose=False)
        subj = subject_level_aggregate(df)

        for grp, color, label in [("healthy", C_HC, "Healthy"),
                                   ("pathological", C_PATH, "Pathological")]:
            g = subj[subj["group"] == grp]
            if len(g) == 0:
                continue
            mode, ci = bootstrap_mode(g["CDC_median"].values)
            ax.hist(g["DCDC"].values,
                    bins=min(40, max(10, len(g) // 5)), alpha=0.45,
                    color=color, density=True,
                    label=f"{label} (N={len(g):,})")
            if len(g) > 5:
                kde = gaussian_kde(g["DCDC"].values)
                x = np.linspace(g["DCDC"].min(), g["DCDC"].max(), 200)
                ax.plot(x, kde(x), color=color, linewidth=1.5)

        ax.axvline(0, color="#888888", linestyle="--", linewidth=0.8)
        ax.set_title(title, fontsize=9, fontweight="bold")
        # Only bottom row gets x-label
        if idx >= 2:
            ax.set_xlabel("$\\Delta$CDC", fontsize=8)
        else:
            ax.set_xlabel("")
        leg = ax.legend(fontsize=6, loc="best")
        leg.get_frame().set(**_LEG_FRAME)
        ax.tick_params(labelsize=7)
        ax.text(0.98, 0.98, ann_note, transform=ax.transAxes, fontsize=5.5,
                ha="right", va="top", style="italic", color="#666")

    save_fig(fig, "SI_Fig5_large_scale", log=log)
