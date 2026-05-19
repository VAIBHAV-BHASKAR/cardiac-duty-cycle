"""SI Figure 8 (v1.3): CODE-15% distributions and age trends.

(Was SI Fig 6 in manuscript v1.0; renumbered to SI Fig 8 in v1.3.)

Group naming: in CODE-15% the 'healthy' group corresponds to the
manuscript's "Patients (non-pathological ECG)" (NPE) — these are tele-ECG
referrals classified as normal_ecg, NOT screened healthy controls. Figure
labels follow the v1.3 patient-centric naming convention.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import gaussian_kde

from ..config import C_NPE, C_PE
from ..data_loader import bootstrap_mode
from ._helpers import save_fig

_AX = 11
_TICK = 10
_TITLE = 12
_PANEL = 14
_LEG = 9
_ANN = 10
_LW = 1.5

_LEG_FRAME = dict(facecolor="white", edgecolor=(0.7, 0.7, 0.7))


_NPE_LABEL = "Patients (non-pathological ECG)"
_PE_LABEL = "Patients (pathological ECG findings)"


def plot(code15_subj, log=print):
    log("  Plotting SI Fig 8 — CODE-15% distributions...")
    # Wider figure to host panel-a outside-right legend cleanly
    fig = plt.figure(figsize=(9.0, 5.5))

    # Panel a narrower (was 0.86 wide) to leave room for the outside-right
    # legend; panels b and c kept on the bottom row.
    ax_dist = fig.add_axes([0.07, 0.57, 0.62, 0.36])
    ax_age  = fig.add_axes([0.07, 0.08, 0.30, 0.36])
    ax_rr   = fig.add_axes([0.43, 0.08, 0.30, 0.36])

    cn   = code15_subj[code15_subj["group"] == "healthy"]
    path = code15_subj[code15_subj["group"] == "pathological"]
    n_total = len(cn) + len(path)

    cn_mode,   cn_ci   = bootstrap_mode(cn["CDC_median"].values)
    path_mode, path_ci = bootstrap_mode(path["CDC_median"].values)

    # Panel a: ΔCDC distribution
    ax_dist.hist(cn["DCDC"].values, bins=80, alpha=0.45, color=C_NPE,
                 density=True,
                 label=f"{_NPE_LABEL}\nN = {len(cn):,}, mode = {cn_mode:.3f}")
    ax_dist.hist(path["DCDC"].values, bins=80, alpha=0.45, color=C_PE,
                 density=True,
                 label=f"{_PE_LABEL}\nN = {len(path):,}, mode = {path_mode:.3f}")

    for grp_data, color in [(cn, C_NPE), (path, C_PE)]:
        kde = gaussian_kde(grp_data["DCDC"].values)
        x = np.linspace(-0.20, 0.22, 300)
        ax_dist.plot(x, kde(x), color=color, linewidth=_LW)

    ax_dist.axvline(0, color="#888888", linestyle="--", linewidth=0.8,
                    label="1/$e$")
    ax_dist.set_xlabel("$\\Delta$CDC", fontsize=_AX)
    ax_dist.set_ylabel("Density", fontsize=_AX)
    ax_dist.set_title(
        f"(a) CODE-15% (N = {n_total:,}; ages 17–100)",
        fontsize=_TITLE, fontweight="bold", loc="left", pad=8)
    ax_dist.set_xlim(-0.20, 0.22)
    ax_dist.tick_params(labelsize=_TICK)

    # Legend: outside the axes on the right, vertically centred, two-line entries
    leg = ax_dist.legend(fontsize=_LEG,
                         bbox_to_anchor=(1.02, 0.5), loc="center left",
                         frameon=True, handlelength=1.5, borderaxespad=0.0)
    leg.get_frame().set(**_LEG_FRAME)

    # Wilcoxon rank-sum p-value, placed just below the legend
    ax_dist.annotate("Wilcoxon rank-sum $p$ < 0.001",
                     xy=(1.02, 0.08), xycoords="axes fraction",
                     fontsize=_ANN, ha="left", va="bottom",
                     fontweight="bold")

    # Annotation method bullets — top-left inside the axes
    ax_dist.text(0.02, 0.97,
                 "• Pan-Tompkins R-peak\n• Tangent T-end",
                 transform=ax_dist.transAxes,
                 fontsize=_LEG, ha="left", va="top",
                 style="italic", color="#555")

    # Panel b: ΔCDC vs Age
    for grp, color, label in [(cn, C_NPE, _NPE_LABEL),
                               (path, C_PE, _PE_LABEL)]:
        g = grp.dropna(subset=["age"])
        g_plot = g.sample(min(5000, len(g)), random_state=42)
        ax_age.scatter(g_plot["age"], g_plot["DCDC"], s=0.3, alpha=0.05,
                       color=color, rasterized=True)
        slope, intercept, r, p, se = stats.linregress(g["age"], g["DCDC"])
        x_line = np.linspace(g["age"].min(), g["age"].max(), 100)
        ax_age.plot(x_line, slope * x_line + intercept, color=color,
                    linewidth=_LW, label=f"{label} (r = {r:.3f})")

    ax_age.axhline(0, color="#888888", linestyle="--", linewidth=0.7)
    # "Optimal (1/e)" label below the dashed line with white background box
    ax_age.text(0.98, 0.02, "Optimal (1/$e$)",
                transform=ax_age.transAxes,
                fontsize=_LEG, ha="right", va="bottom", color="#444",
                bbox=dict(facecolor="white", edgecolor="none",
                          alpha=0.8, pad=2))

    ax_age.set_xlabel("Age (years)", fontsize=_AX)
    ax_age.set_ylabel("$\\Delta$CDC", fontsize=_AX)
    ax_age.set_title("(b) $\\Delta$CDC vs Age",
                     fontsize=_TITLE, fontweight="bold", loc="left", pad=8)
    leg = ax_age.legend(fontsize=_LEG, loc="upper left")
    leg.get_frame().set(**_LEG_FRAME)
    ax_age.tick_params(labelsize=_TICK)

    # Panel c: RR vs Age (capped at 1800 ms to remove extreme outliers)
    for grp, color, label in [(cn, C_NPE, _NPE_LABEL),
                               (path, C_PE, _PE_LABEL)]:
        g = grp.dropna(subset=["age"])
        g_plot = g.sample(min(5000, len(g)), random_state=42)
        ax_rr.scatter(g_plot["age"], g_plot["RR_ms"], s=0.3, alpha=0.05,
                      color=color, rasterized=True)
        slope, intercept, r, p, se = stats.linregress(g["age"], g["RR_ms"])
        x_line = np.linspace(g["age"].min(), g["age"].max(), 100)
        ax_rr.plot(x_line, slope * x_line + intercept, color=color,
                   linewidth=_LW, label=f"{label} (r = {r:.3f})")

    ax_rr.set_xlabel("Age (years)", fontsize=_AX)
    ax_rr.set_ylabel("RR interval (ms)", fontsize=_AX)
    ax_rr.set_title("(c) RR vs Age",
                    fontsize=_TITLE, fontweight="bold", loc="left", pad=8)
    ax_rr.set_ylim(top=1800)  # cap removes ~20 subjects > 99.5th percentile
    leg = ax_rr.legend(fontsize=_LEG, loc="upper left")
    leg.get_frame().set(**_LEG_FRAME)
    ax_rr.tick_params(labelsize=_TICK)

    save_fig(fig, "SI_Fig8_code15", log=log)
