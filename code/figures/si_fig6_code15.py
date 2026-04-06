"""SI Figure 6: CODE-15% distributions + age trends."""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import gaussian_kde

from ..config import C_CN, C_PATH
from ..data_loader import bootstrap_mode
from ._helpers import save_fig


def plot(code15_subj, log=print):
    log("  Plotting SI Fig 6 — CODE-15% distributions...")
    fig = plt.figure(figsize=(7.2, 5.5))

    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3)
    ax_dist = fig.add_subplot(gs[0, :])
    ax_age = fig.add_subplot(gs[1, 0])
    ax_rr = fig.add_subplot(gs[1, 1])

    cn = code15_subj[code15_subj["group"] == "healthy"]
    path = code15_subj[code15_subj["group"] == "pathological"]

    cn_mode, cn_ci = bootstrap_mode(cn["CDC_median"].values)
    path_mode, path_ci = bootstrap_mode(path["CDC_median"].values)

    ax_dist.hist(cn["DCDC"].values, bins=80, alpha=0.45, color=C_CN, density=True,
                 label=f"Normal ECG (N={len(cn):,}, mode={cn_mode:.3f})")
    ax_dist.hist(path["DCDC"].values, bins=80, alpha=0.45, color=C_PATH, density=True,
                 label=f"Pathological (N={len(path):,}, mode={path_mode:.3f})")

    for grp_data, color in [(cn, C_CN), (path, C_PATH)]:
        kde = gaussian_kde(grp_data["DCDC"].values)
        x = np.linspace(-0.15, 0.22, 300)
        ax_dist.plot(x, kde(x), color=color, linewidth=1.5)

    ax_dist.axvline(0, color="#888888", linestyle="--", linewidth=0.8, label="1/$e$")
    ax_dist.set_xlabel("$\\Delta$CDC", fontsize=9)
    ax_dist.set_ylabel("Density", fontsize=9)
    ax_dist.set_title("(a) CODE-15% CDC distributions", fontsize=9, fontweight="bold", loc="left")
    ax_dist.legend(fontsize=6)
    ax_dist.set_xlim(-0.15, 0.22)
    ax_dist.tick_params(labelsize=8)

    for grp, color, label in [(cn, C_CN, "Normal"), (path, C_PATH, "Pathological")]:
        g = grp.dropna(subset=["age"])
        g_plot = g.sample(min(5000, len(g)), random_state=42)
        ax_age.scatter(g_plot["age"], g_plot["DCDC"], s=0.3, alpha=0.05, color=color, rasterized=True)
        slope, intercept, r, p, se = stats.linregress(g["age"], g["DCDC"])
        x_line = np.linspace(g["age"].min(), g["age"].max(), 100)
        ax_age.plot(x_line, slope * x_line + intercept, color=color, linewidth=1.8,
                    label=f"{label} (r={r:.3f})")

    ax_age.axhline(0, color="#888888", linestyle="--", linewidth=0.7)
    ax_age.set_xlabel("Age (years)", fontsize=9)
    ax_age.set_ylabel("$\\Delta$CDC", fontsize=9)
    ax_age.set_title("(b) CDC vs Age", fontsize=9, fontweight="bold", loc="left")
    ax_age.legend(fontsize=6)
    ax_age.tick_params(labelsize=8)

    for grp, color, label in [(cn, C_CN, "Normal"), (path, C_PATH, "Pathological")]:
        g = grp.dropna(subset=["age"])
        g_plot = g.sample(min(5000, len(g)), random_state=42)
        ax_rr.scatter(g_plot["age"], g_plot["RR_ms"], s=0.3, alpha=0.05, color=color, rasterized=True)
        slope, intercept, r, p, se = stats.linregress(g["age"], g["RR_ms"])
        x_line = np.linspace(g["age"].min(), g["age"].max(), 100)
        ax_rr.plot(x_line, slope * x_line + intercept, color=color, linewidth=1.8,
                   label=f"{label} (r={r:.3f})")

    ax_rr.set_xlabel("Age (years)", fontsize=9)
    ax_rr.set_ylabel("RR interval (ms)", fontsize=9)
    ax_rr.set_title("(c) RR vs Age", fontsize=9, fontweight="bold", loc="left")
    ax_rr.legend(fontsize=6)
    ax_rr.tick_params(labelsize=8)

    save_fig(fig, "SI_Fig6_code15", log=log)
