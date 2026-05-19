"""SI Figure 7 (v1.3): Large-scale CDC distributions (4 databases).

(Was SI Fig 5 in manuscript v1.0; renumbered to SI Fig 7 in v1.3.)

Note: The PTB-XL panel shows NPE vs PE; the others (Fantasia, AA) show only
healthy controls. PTB shows healthy vs pathological.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde, ranksums

from ..config import C_HC, C_NPE, C_PE
from ..data_loader import (load_beats, apply_quality_filters,
                           subject_level_aggregate, bootstrap_mode)
from ._helpers import save_fig

_LEG_FRAME = dict(facecolor="white", edgecolor=(0.7, 0.7, 0.7))
_LEG_FS = 7
_TICK_FS = 8


def plot(log=print):
    log("  Plotting SI Fig 7 — Large-scale distributions...")
    fig = plt.figure(figsize=(9.0, 6.5))
    plt.subplots_adjust(left=0.07, right=0.97, top=0.92, bottom=0.10,
                        wspace=0.55, hspace=0.50)
    axes_flat = [fig.add_subplot(2, 2, i + 1) for i in range(4)]

    # (title, file, annotation-bullets, healthy_label, pathological_label, healthy_color)
    # Bullet order is uniform across panels: R-peaks first, T-end second.
    datasets = [
        ("Fantasia", "fantasia_beats.csv",
         "• Database R-peaks\n• Tangent T-end",
         "Healthy controls", "Patients (pathological ECG findings)", C_HC),
        ("Autonomic Aging", "autonomic_aging_beats.csv",
         "• Pan-Tompkins R-peak\n• Tangent T-end",
         "Healthy controls", "Patients (pathological ECG findings)", C_HC),
        ("PTB: ages 17–87", "ptb_beats.csv",
         "• Pan-Tompkins R-peak\n• Manual T-end (5 referees)",
         "Healthy controls", "Patients (pathological ECG findings)", C_HC),
        ("PTB-XL: ages 2–90", "ptbxl_beats.csv",
         "• ECGDeli R-peak\n• ECGDeli T-end",
         "Patients (non-pathological ECG)",
         "Patients (pathological ECG findings)", C_NPE),
    ]

    for idx, (title, fname, ann_note, h_label, p_label, h_color) in enumerate(datasets):
        ax = axes_flat[idx]
        df = apply_quality_filters(load_beats(fname), verbose=False)
        subj = subject_level_aggregate(df)

        group_dcdc = {}
        for grp, color, label in [("healthy", h_color, h_label),
                                   ("pathological", C_PE, p_label)]:
            g = subj[subj["group"] == grp]
            if len(g) == 0:
                continue
            mode, ci = bootstrap_mode(g["CDC_median"].values)
            group_dcdc[grp] = g["DCDC"].values
            ax.hist(g["DCDC"].values,
                    bins=min(40, max(10, len(g) // 5)), alpha=0.45,
                    color=color, density=True,
                    label=f"{label}\n(N={len(g):,}, mode={mode:.3f})")
            if len(g) > 5:
                kde = gaussian_kde(g["DCDC"].values)
                x = np.linspace(g["DCDC"].min(), g["DCDC"].max(), 200)
                ax.plot(x, kde(x), color=color, linewidth=1.5)

        ax.axvline(0, color="#888888", linestyle="--", linewidth=0.8)
        ax.set_title(title, fontsize=10, fontweight="bold", pad=10)
        if idx >= 2:
            ax.set_xlabel("$\\Delta$CDC", fontsize=9)
        else:
            ax.set_xlabel("")
        ax.set_ylabel("Density", fontsize=9)
        ax.tick_params(labelsize=_TICK_FS)

        # Legend: outside the axes on the right, vertically centred
        leg = ax.legend(fontsize=_LEG_FS,
                        bbox_to_anchor=(1.02, 0.5), loc="center left",
                        frameon=True, handlelength=1.5, borderaxespad=0.0)
        leg.get_frame().set(**_LEG_FRAME)

        # Wilcoxon rank-sum p-value annotation (only when both groups present)
        if "healthy" in group_dcdc and "pathological" in group_dcdc:
            _, p_val = ranksums(group_dcdc["healthy"], group_dcdc["pathological"])
            p_str = ("Wilcoxon rank-sum $p$ < 0.001"
                     if p_val < 1e-3
                     else f"Wilcoxon rank-sum $p$ = {p_val:.3f}")
            # Placed just below the legend, in figure-fraction relative to axes
            ax.annotate(p_str,
                        xy=(1.02, 0.05), xycoords="axes fraction",
                        fontsize=_LEG_FS, ha="left", va="bottom")

        # Annotation method bullets — top-left inside the axes
        ax.text(0.03, 0.97, ann_note, transform=ax.transAxes,
                fontsize=_LEG_FS - 0.5, ha="left", va="top",
                style="italic", color="#555")

    # Skip tight_layout (incompatible with external legends); save_fig uses
    # bbox_inches="tight" which will capture the outside legends correctly.
    save_fig(fig, "SI_Fig7_large_scale", log=log)
