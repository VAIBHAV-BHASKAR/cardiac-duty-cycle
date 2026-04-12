"""SI Figure 4: Gold-standard CDC distributions (LUDB + QTDB)."""

import numpy as np
from scipy.stats import gaussian_kde

from ..config import C_HC, C_PATH, C_SD, ONE_OVER_E
from ..data_loader import load_beats, apply_quality_filters, subject_level_aggregate, bootstrap_mode
from ._helpers import save_fig, new_figure


def plot(log=print):
    log("  Plotting SI Fig 4 — Gold-standard distributions...")
    fig, axes = new_figure(2, 1, figsize=(4.7, 6))

    for ax_idx, (fname, title) in enumerate([
        ("ludb_beats.csv", "LUDB — Full manual annotation"),
        ("qtdb_beats.csv", "QTDB — Full manual annotation"),
    ]):
        ax = axes[ax_idx]
        df = apply_quality_filters(load_beats(fname), verbose=False)
        subj = subject_level_aggregate(df)

        color_map = {"healthy": C_HC, "pathological": C_PATH, "sudden_death": C_SD}
        label_map = {"healthy": "Healthy", "pathological": "Pathological", "sudden_death": "Sudden death"}

        for grp in ["healthy", "pathological", "sudden_death"]:
            g = subj[subj["group"] == grp]
            if len(g) == 0:
                continue
            mode, ci = bootstrap_mode(g["CDC_median"].values)
            ax.hist(g["DCDC"].values, bins=25, alpha=0.45, color=color_map[grp], density=True,
                    label=f"{label_map[grp]} (N={len(g)}, mode={mode:.3f})")
            if len(g) > 5:
                kde = gaussian_kde(g["DCDC"].values)
                x = np.linspace(g["DCDC"].min(), g["DCDC"].max(), 200)
                ax.plot(x, kde(x), color=color_map[grp], linewidth=1.5)

        ax.axvline(0, color="#888888", linestyle="--", linewidth=0.8, label="1/$e$")
        ax.set_xlabel("$\\Delta$CDC", fontsize=9)
        ax.set_ylabel("Density", fontsize=9)
        ax.set_title(title, fontsize=9, fontweight="bold")
        ax.legend(fontsize=6, loc="upper right")
        ax.tick_params(labelsize=8)
        # Place p-value text at mid-right (matches MATLAB position)
        if ax_idx == 1:
            ax.text(0.97, 0.5, "$\\it{p}$ = n.s. (K-W)",
                    transform=ax.transAxes, fontsize=6, ha="right",
                    va="center", color="#555")

    save_fig(fig, "SI_Fig4_gold_standard", log=log)
