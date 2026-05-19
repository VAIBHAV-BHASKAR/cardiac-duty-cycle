"""SI Figure 7 (v1.3): Large-scale CDC distributions (4 databases).

(Was SI Fig 5 in manuscript v1.0; renumbered to SI Fig 7 in v1.3.)

Note: The PTB-XL panel shows NPE vs PE; the others (Fantasia, AA) show only
healthy controls. PTB shows healthy vs pathological.
"""

import numpy as np
from scipy.stats import gaussian_kde

from ..config import C_HC, C_NPE, C_PE
from ..data_loader import (load_beats, apply_quality_filters,
                           subject_level_aggregate, bootstrap_mode)
from ._helpers import save_fig, new_figure

_LEG_FRAME = dict(facecolor="white", edgecolor=(0.7, 0.7, 0.7))


def plot(log=print):
    log("  Plotting SI Fig 7 — Large-scale distributions...")
    fig, axes = new_figure(2, 2, figsize=(7.2, 6))
    axes_flat = [axes[0][0], axes[0][1], axes[1][0], axes[1][1]]

    # (title, file, annotation, healthy_label, pathological_label, healthy_color)
    datasets = [
        ("Fantasia", "fantasia_beats.csv",
         "Database R-peaks,\ntangent T-end",
         "Healthy controls", "Pathological", C_HC),
        ("Autonomic Aging", "autonomic_aging_beats.csv",
         "Fully automatic\n(Pan-Tompkins + tangent)",
         "Healthy controls", "Pathological", C_HC),
        ("PTB: ages 17–87", "ptb_beats.csv",
         "Manual T-end (5 referees),\nPan-Tompkins R-peak",
         "Healthy controls", "Patients (pathological ECG findings)", C_HC),
        ("PTB-XL: ages 2–90", "ptbxl_beats.csv",
         "ECGDeli automatic\n(R-peak + T-end)",
         "Patients (non-pathological ECG)",
         "Patients (pathological ECG findings)", C_NPE),
    ]

    for idx, (title, fname, ann_note, h_label, p_label, h_color) in enumerate(datasets):
        ax = axes_flat[idx]
        df = apply_quality_filters(load_beats(fname), verbose=False)
        subj = subject_level_aggregate(df)

        for grp, color, label in [("healthy", h_color, h_label),
                                   ("pathological", C_PE, p_label)]:
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

    save_fig(fig, "SI_Fig7_large_scale", log=log)
