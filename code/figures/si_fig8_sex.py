"""SI Figure 8: Sex-stratified CDC aging."""

import numpy as np
from scipy import stats

from ..config import C_HC, C_CN, C_PATH
from ._helpers import save_fig, new_figure

GROUPS = [
    ("HealthyControl", C_HC, "HC"),
    ("ClinicallyNormal", C_CN, "CN"),
    ("Pathological", C_PATH, "Path"),
]


def plot(pooled, log=print):
    log("  Plotting SI Fig 8 — Sex-stratified CDC aging...")

    pooled = pooled.copy()
    pooled["sex_clean"] = pooled["sex"].astype(str).str.upper().str.strip()
    pooled.loc[pooled["sex_clean"].isin(["F", "FEMALE", "0", "0.0"]), "sex_std"] = "F"
    pooled.loc[pooled["sex_clean"].isin(["M", "MALE", "1", "1.0"]), "sex_std"] = "M"

    fig, axes = new_figure(2, 2, figsize=(7.2, 6))

    for col_idx, sex in enumerate(["F", "M"]):
        sex_data = pooled[pooled["sex_std"] == sex]

        # Row 0: ΔCDC vs Age
        ax = axes[0][col_idx]
        for grp, color, label in GROUPS:
            g = sex_data[sex_data["hgroup"] == grp].dropna(subset=["age"])
            if len(g) < 5:
                continue
            ax.scatter(g["age"], g["DCDC"], s=0.5, alpha=0.1, color=color, rasterized=True)
            slope, intercept, r, p, se = stats.linregress(g["age"], g["DCDC"])
            x_line = np.linspace(max(10, g["age"].min()), min(100, g["age"].max()), 100)
            ax.plot(x_line, slope * x_line + intercept, color=color, linewidth=1.8,
                    label=f"{label} (r={r:.3f})")

        ax.axhline(0, color="#888", linestyle="--", linewidth=0.7)
        lbl = "(a)" if col_idx == 0 else "(b)"
        ax.set_title(f"{lbl} $\\Delta$CDC — {sex}", fontsize=9, fontweight="bold", loc="left")
        ax.set_xlabel("Age (years)", fontsize=8)
        ax.set_ylabel("$\\Delta$CDC", fontsize=8)
        ax.set_xlim(0, 105)
        ax.set_ylim(-0.15, 0.15)
        ax.legend(fontsize=5.5, loc="upper left")
        ax.tick_params(labelsize=7)

        # Row 1: Diastole vs Age
        ax = axes[1][col_idx]
        for grp, color, label in GROUPS:
            g = sex_data[sex_data["hgroup"] == grp].dropna(subset=["age", "diastole_ms"])
            if len(g) < 5:
                continue
            ax.scatter(g["age"], g["diastole_ms"], s=0.5, alpha=0.1, color=color, rasterized=True)
            slope, intercept, r, p, se = stats.linregress(g["age"], g["diastole_ms"])
            x_line = np.linspace(max(10, g["age"].min()), min(100, g["age"].max()), 100)
            ax.plot(x_line, slope * x_line + intercept, color=color, linewidth=1.8, label=label)

        lbl = "(c)" if col_idx == 0 else "(d)"
        ax.set_title(f"{lbl} Diastole — {sex}", fontsize=9, fontweight="bold", loc="left")
        ax.set_xlabel("Age (years)", fontsize=8)
        ax.set_ylabel("Diastole (ms)", fontsize=8)
        ax.set_xlim(0, 105)
        ax.legend(fontsize=5.5, loc="upper right")
        ax.tick_params(labelsize=7)

    save_fig(fig, "SI_Fig8_sex_stratified", log=log)
