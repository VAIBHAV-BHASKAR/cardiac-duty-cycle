"""Figure 1: CDC aging trajectories (2 panels)."""

import numpy as np
from scipy import stats

from ..config import C_HC, C_CN, C_PATH, FIGURES_DIR
from ..data_loader import bootstrap_mode
from ._helpers import save_fig, new_figure


GROUPS = [
    ("HealthyControl", C_HC, "Healthy control"),
    ("ClinicallyNormal", C_CN, "Clinically normal"),
    ("Pathological", C_PATH, "Pathological"),
]


def plot(pooled, log=print):
    log("  Plotting Fig 1 — CDC aging trajectories...")
    fig, axes = new_figure(1, 2, figsize=(7.2, 3.5))

    # Panel a: ΔCDC vs Age
    ax = axes[0]
    for grp, color, label in GROUPS:
        g = pooled[pooled["hgroup"] == grp]
        mode, ci = bootstrap_mode(g["CDC_median"].values)
        ax.scatter(g["age"], g["DCDC"], s=0.3, alpha=0.08, color=color, rasterized=True)
        slope, intercept, r, p, se = stats.linregress(g["age"], g["DCDC"])
        x_line = np.linspace(max(10, g["age"].min()), min(100, g["age"].max()), 100)
        ax.plot(x_line, slope * x_line + intercept, color=color, linewidth=1.8,
                label=f"{label}\nmode={mode:.3f} [{ci[0]:.3f}, {ci[1]:.3f}]")

    ax.axhline(0, color="#888888", linestyle="--", linewidth=0.7)
    ax.set_xlabel("Age (years)", fontsize=9)
    ax.set_ylabel("$\\Delta$CDC (CDC $-$ 1/$e$)", fontsize=9)
    ax.set_title("(a)", fontsize=10, fontweight="bold", loc="left")
    ax.set_xlim(0, 105)
    ax.set_ylim(-0.15, 0.15)
    ax.legend(fontsize=5.5, loc="upper left", framealpha=0.9)
    ax.tick_params(labelsize=8)

    # Panel b: Diastole vs Age
    ax = axes[1]
    for grp, color, label in GROUPS:
        g = pooled[pooled["hgroup"] == grp]
        ax.scatter(g["age"], g["diastole_ms"], s=0.3, alpha=0.08, color=color, rasterized=True)
        slope, intercept, r, p, se = stats.linregress(g["age"], g["diastole_ms"])
        x_line = np.linspace(max(10, g["age"].min()), min(100, g["age"].max()), 100)
        ax.plot(x_line, slope * x_line + intercept, color=color, linewidth=1.8, label=label)

    ax.set_xlabel("Age (years)", fontsize=9)
    ax.set_ylabel("Thermodynamic diastole (ms)", fontsize=9)
    ax.set_title("(b)", fontsize=10, fontweight="bold", loc="left")
    ax.set_xlim(0, 105)
    ax.legend(fontsize=6, loc="upper right", framealpha=0.9)
    ax.tick_params(labelsize=8)

    save_fig(fig, "Fig1_CDC_aging", log=log)
