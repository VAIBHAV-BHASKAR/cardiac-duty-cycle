"""Figure 1: CDC aging trajectories (2 panels).

Matches MATLAB plot_Fig1.m (commit 58e9817): 2x scale strategy for Nature
Aging print, multi-line legends, hardcoded panel positions.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

from ..config import C_HC, C_CN, C_PATH, FIGURES_DIR
from ..data_loader import bootstrap_mode
from ._helpers import save_fig


GROUPS = [
    ("HealthyControl", C_HC, "Healthy control"),
    ("ClinicallyNormal", C_CN, "Clinically normal"),
    ("Pathological", C_PATH, "Pathological"),
]

# 2x scale: draw at double Nature Aging column width so 50% print scaling
# yields correct font sizes.
_S = 2
_W = 18.3 * _S / 2.54   # cm -> inches
_H = 9.0  * _S / 2.54
_TICK  = 8 * _S
_AX    = 9 * _S
_PANEL = 11 * _S
_LEG   = 8 * _S
_DOT   = 4 * _S
_LW    = 1.5 * _S
_REF   = _LW * 0.7


def plot(pooled, log=print):
    log("  Plotting Fig 1 — CDC aging trajectories...")
    fig = plt.figure(figsize=(_W, _H))

    ax1 = fig.add_axes([0.10, 0.15, 0.38, 0.75])
    ax2 = fig.add_axes([0.58, 0.15, 0.38, 0.75])

    # Panel a: ΔCDC vs Age
    for grp, color, label in GROUPS:
        g = pooled[pooled["hgroup"] == grp]
        mode, ci = bootstrap_mode(g["CDC_median"].values)
        ax1.scatter(g["age"], g["DCDC"], s=_DOT, alpha=0.08, color=color,
                    rasterized=True)
        slope, intercept, r, p, se = stats.linregress(g["age"], g["DCDC"])
        x_line = np.linspace(max(10, g["age"].min()),
                             min(100, g["age"].max()), 100)
        ax1.plot(x_line, slope * x_line + intercept, color=color,
                 linewidth=_LW,
                 label=f"{label} (N={len(g):,})\n"
                       f"mode: {mode:.3f} [{ci[0]:.3f}, {ci[1]:.3f}]")

    ax1.axhline(0, color="#888888", linestyle="--", linewidth=_REF,
                label="Optimal (1/$e$)")
    ax1.set_xlabel("Age (years)", fontsize=_AX)
    ax1.set_ylabel("$\\Delta$CDC (CDC $-$ 1/$e$)", fontsize=_AX)
    ax1.set_xlim(0, 105)
    ax1.set_ylim(-0.15, 0.15)
    ax1.legend(fontsize=_LEG * 0.65, loc="lower right", framealpha=0.9)
    ax1.tick_params(labelsize=_TICK, length=4)
    ax1.text(-0.12, 1.04, "(a)", transform=ax1.transAxes, fontsize=_PANEL,
             fontweight="bold", va="bottom")

    # Panel b: Diastole vs Age
    for grp, color, label in GROUPS:
        g = pooled[pooled["hgroup"] == grp]
        ax2.scatter(g["age"], g["diastole_ms"], s=_DOT, alpha=0.08,
                    color=color, rasterized=True)
        slope, intercept, r, p, se = stats.linregress(g["age"],
                                                       g["diastole_ms"])
        x_line = np.linspace(max(10, g["age"].min()),
                             min(100, g["age"].max()), 100)
        ax2.plot(x_line, slope * x_line + intercept, color=color,
                 linewidth=_LW, label=label)

    ax2.set_xlabel("Age (years)", fontsize=_AX)
    ax2.set_ylabel("Thermodynamic diastole (ms)", fontsize=_AX)
    ax2.set_xlim(0, 105)
    ax2.legend(fontsize=_LEG, loc="upper right", framealpha=0.9)
    ax2.tick_params(labelsize=_TICK, length=4)
    ax2.text(-0.12, 1.04, "(b)", transform=ax2.transAxes, fontsize=_PANEL,
             fontweight="bold", va="bottom")

    for ax in [ax1, ax2]:
        for spine in ax.spines.values():
            spine.set_linewidth(0.5 * _S)

    save_fig(fig, "Fig1_CDC_aging", log=log)
