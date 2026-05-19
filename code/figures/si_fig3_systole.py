"""SI Figure 3 (v1.3): Systole & RR interval across the lifespan.

(Was SI Fig 2 in manuscript v1.0; renumbered to SI Fig 3 in v1.3.)
"""

import numpy as np
from scipy import stats

from ..config import C_HC, C_NPE, C_PE
from ._helpers import save_fig, new_figure

# v1.3 group taxonomy
GROUPS = [
    ("HC",  C_HC,  "Healthy controls"),
    ("NPE", C_NPE, "Patients (non-pathological ECG)"),
    ("PE",  C_PE,  "Patients (pathological ECG findings)"),
]


def plot(pooled, log=print):
    log("  Plotting SI Fig 3 — Systole & RR across lifespan...")
    fig, axes = new_figure(1, 2, figsize=(7.2, 3.5))

    # v1.3 SI Fig 3 caption: "Mechanical systole duration (ms, R-peak to T-end)"
    # (Main-text Results dropped the "Mechanical" qualifier; SI caption kept it.)
    for panel_idx, (y_col, y_label, title) in enumerate([
        ("RT_ms", "Mechanical systole duration (ms)", "(a) Systole duration"),
        ("RR_ms", "RR interval (ms)", "(b) RR interval"),
    ]):
        ax = axes[panel_idx]
        for grp, color, label in GROUPS:
            g = pooled[pooled["hgroup"] == grp].dropna(subset=["age", y_col])
            ax.scatter(g["age"], g[y_col], s=0.3, alpha=0.08, color=color, rasterized=True)
            slope, intercept, r, p, se = stats.linregress(g["age"], g[y_col])
            x_line = np.linspace(max(10, g["age"].min()), min(100, g["age"].max()), 100)
            ax.plot(x_line, slope * x_line + intercept, color=color, linewidth=1.8, label=label)

        ax.set_xlabel("Age (years)", fontsize=9)
        ax.set_ylabel(y_label, fontsize=9)
        ax.set_title(title, fontsize=9, fontweight="bold", loc="left")
        ax.set_xlim(0, 105)
        ax.legend(fontsize=6, loc="upper right" if panel_idx == 0 else "best")
        ax.tick_params(labelsize=8)

    save_fig(fig, "SI_Fig3_systole_RR", log=log)
