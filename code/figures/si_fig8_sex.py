"""SI Figure 8: Sex-stratified CDC aging.

Matches MATLAB plot_SI_Fig8.m: taller figure (18 cm), manual panel
positions, larger fonts, legend boxes on, slope annotations repositioned.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

from ..config import C_HC, C_CN, C_PATH
from ._helpers import save_fig

GROUPS = [
    ("HealthyControl", C_HC, "HC"),
    ("ClinicallyNormal", C_CN, "CN"),
    ("Pathological", C_PATH, "Path"),
]

_AX = 10
_TICK = 9
_TITLE = 11
_PANEL = 13
_LEG = 7.5
_SLOPE = 7.5
_LW = 2.0

_LEG_FRAME = dict(facecolor="white", edgecolor=(0.7, 0.7, 0.7))


def plot(pooled, log=print):
    log("  Plotting SI Fig 8 — Sex-stratified CDC aging...")

    pooled = pooled.copy()
    pooled["sex_clean"] = pooled["sex"].astype(str).str.upper().str.strip()
    pooled.loc[pooled["sex_clean"].isin(["F", "FEMALE", "0", "0.0"]),
               "sex_std"] = "F"
    pooled.loc[pooled["sex_clean"].isin(["M", "MALE", "1", "1.0"]),
               "sex_std"] = "M"

    fig = plt.figure(figsize=(7.2, 7.1))  # taller: 18 cm

    positions = [
        [0.10, 0.58, 0.38, 0.35],   # (a) top-left
        [0.57, 0.58, 0.38, 0.35],   # (b) top-right
        [0.10, 0.08, 0.38, 0.35],   # (c) bottom-left
        [0.57, 0.08, 0.38, 0.35],   # (d) bottom-right
    ]

    panel_labels = ["(a)", "(b)", "(c)", "(d)"]

    for col_idx, sex in enumerate(["F", "M"]):
        sex_data = pooled[pooled["sex_std"] == sex]

        # Row 0: ΔCDC vs Age
        ax = fig.add_axes(positions[col_idx])
        for grp, color, label in GROUPS:
            g = sex_data[sex_data["hgroup"] == grp].dropna(subset=["age"])
            if len(g) < 5:
                continue
            ax.scatter(g["age"], g["DCDC"], s=0.5, alpha=0.1, color=color,
                       rasterized=True)
            slope, intercept, r, p, se = stats.linregress(g["age"], g["DCDC"])
            x_line = np.linspace(max(10, g["age"].min()),
                                 min(100, g["age"].max()), 100)
            ax.plot(x_line, slope * x_line + intercept, color=color,
                    linewidth=_LW, label=f"{label} (r={r:.3f})")

        ax.axhline(0, color="#888", linestyle="--", linewidth=0.7)
        ax.set_title(f"{panel_labels[col_idx]} $\\Delta$CDC — {sex}",
                     fontsize=_TITLE, fontweight="bold", loc="left")
        ax.set_xlabel("Age (years)", fontsize=_AX)
        ax.set_ylabel("$\\Delta$CDC", fontsize=_AX)
        ax.set_xlim(0, 105)
        ax.set_ylim(-0.15, 0.15)
        leg = ax.legend(fontsize=_LEG,
                        loc="lower left" if col_idx == 0 else "best")
        leg.get_frame().set(**_LEG_FRAME)
        ax.tick_params(labelsize=_TICK)

        # Row 1: Diastole vs Age
        ax = fig.add_axes(positions[2 + col_idx])
        for grp, color, label in GROUPS:
            g = sex_data[sex_data["hgroup"] == grp].dropna(
                subset=["age", "diastole_ms"])
            if len(g) < 5:
                continue
            ax.scatter(g["age"], g["diastole_ms"], s=0.5, alpha=0.1,
                       color=color, rasterized=True)
            slope, intercept, r, p, se = stats.linregress(
                g["age"], g["diastole_ms"])
            x_line = np.linspace(max(10, g["age"].min()),
                                 min(100, g["age"].max()), 100)
            ax.plot(x_line, slope * x_line + intercept, color=color,
                    linewidth=_LW, label=label)

        ax.set_title(f"{panel_labels[2 + col_idx]} Diastole — {sex}",
                     fontsize=_TITLE, fontweight="bold", loc="left")
        ax.set_xlabel("Age (years)", fontsize=_AX)
        ax.set_ylabel("Diastole (ms)", fontsize=_AX)
        ax.set_xlim(0, 105)
        leg = ax.legend(fontsize=_LEG, loc="upper right")
        leg.get_frame().set(**_LEG_FRAME)
        ax.tick_params(labelsize=_TICK)

    save_fig(fig, "SI_Fig8_sex_stratified", log=log)
