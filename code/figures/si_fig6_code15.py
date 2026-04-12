"""SI Figure 6: CODE-15% distributions + age trends.

Matches MATLAB plot_SI_Fig6.m: manual panel positions, larger fonts,
legend boxes on, age range in title, italic p-value.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import gaussian_kde

from ..config import C_CN, C_PATH
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


def plot(code15_subj, log=print):
    log("  Plotting SI Fig 6 — CODE-15% distributions...")
    fig = plt.figure(figsize=(7.2, 5.5))

    ax_dist = fig.add_axes([0.09, 0.57, 0.86, 0.37])
    ax_age  = fig.add_axes([0.09, 0.08, 0.40, 0.38])
    ax_rr   = fig.add_axes([0.57, 0.08, 0.40, 0.38])

    cn   = code15_subj[code15_subj["group"] == "healthy"]
    path = code15_subj[code15_subj["group"] == "pathological"]

    cn_mode,   cn_ci   = bootstrap_mode(cn["CDC_median"].values)
    path_mode, path_ci = bootstrap_mode(path["CDC_median"].values)

    ax_dist.hist(cn["DCDC"].values, bins=80, alpha=0.45, color=C_CN,
                 density=True,
                 label=f"Normal ECG (N={len(cn):,}, mode={cn_mode:.3f})")
    ax_dist.hist(path["DCDC"].values, bins=80, alpha=0.45, color=C_PATH,
                 density=True,
                 label=f"Pathological (N={len(path):,}, mode={path_mode:.3f})")

    for grp_data, color in [(cn, C_CN), (path, C_PATH)]:
        kde = gaussian_kde(grp_data["DCDC"].values)
        x = np.linspace(-0.15, 0.22, 300)
        ax_dist.plot(x, kde(x), color=color, linewidth=_LW)

    ax_dist.axvline(0, color="#888888", linestyle="--", linewidth=0.8,
                    label="1/$e$")
    ax_dist.set_xlabel("$\\Delta$CDC", fontsize=_AX)
    ax_dist.set_ylabel("Density", fontsize=_AX)
    ax_dist.set_title(
        "(a) CODE-15%: Clinically Normal vs Pathological (ages 17\u2013100)",
        fontsize=_TITLE, fontweight="bold", loc="left")
    leg = ax_dist.legend(fontsize=_LEG)
    leg.get_frame().set(**_LEG_FRAME)
    ax_dist.set_xlim(-0.15, 0.22)
    ax_dist.tick_params(labelsize=_TICK)
    # italic p-value annotation
    xlims = ax_dist.get_xlim()
    ylims = ax_dist.get_ylim()
    ax_dist.text(xlims[1] - 0.01 * (xlims[1] - xlims[0]),
                 ylims[0] + 0.5 * (ylims[1] - ylims[0]),
                 "$\\it{p}$ < 0.001", fontsize=_ANN, ha="right")

    # Panel b: CDC vs Age
    for grp, color, label in [(cn, C_CN, "Normal"),
                               (path, C_PATH, "Pathological")]:
        g = grp.dropna(subset=["age"])
        g_plot = g.sample(min(5000, len(g)), random_state=42)
        ax_age.scatter(g_plot["age"], g_plot["DCDC"], s=0.3, alpha=0.05,
                       color=color, rasterized=True)
        slope, intercept, r, p, se = stats.linregress(g["age"], g["DCDC"])
        x_line = np.linspace(g["age"].min(), g["age"].max(), 100)
        ax_age.plot(x_line, slope * x_line + intercept, color=color,
                    linewidth=_LW, label=f"{label} (r={r:.3f})")

    ax_age.axhline(0, color="#888888", linestyle="--", linewidth=0.7)
    ax_age.set_xlabel("Age (years)", fontsize=_AX)
    ax_age.set_ylabel("$\\Delta$CDC", fontsize=_AX)
    ax_age.set_title("(b) CDC vs Age", fontsize=_TITLE, fontweight="bold",
                     loc="left")
    leg = ax_age.legend(fontsize=_LEG)
    leg.get_frame().set(**_LEG_FRAME)
    ax_age.tick_params(labelsize=_TICK)

    # Panel c: RR vs Age
    for grp, color, label in [(cn, C_CN, "Normal"),
                               (path, C_PATH, "Pathological")]:
        g = grp.dropna(subset=["age"])
        g_plot = g.sample(min(5000, len(g)), random_state=42)
        ax_rr.scatter(g_plot["age"], g_plot["RR_ms"], s=0.3, alpha=0.05,
                      color=color, rasterized=True)
        slope, intercept, r, p, se = stats.linregress(g["age"], g["RR_ms"])
        x_line = np.linspace(g["age"].min(), g["age"].max(), 100)
        ax_rr.plot(x_line, slope * x_line + intercept, color=color,
                   linewidth=_LW, label=f"{label} (r={r:.3f})")

    ax_rr.set_xlabel("Age (years)", fontsize=_AX)
    ax_rr.set_ylabel("RR interval (ms)", fontsize=_AX)
    ax_rr.set_title("(c) RR vs Age", fontsize=_TITLE, fontweight="bold",
                    loc="left")
    leg = ax_rr.legend(fontsize=_LEG)
    leg.get_frame().set(**_LEG_FRAME)
    ax_rr.tick_params(labelsize=_TICK)

    save_fig(fig, "SI_Fig6_code15", log=log)
