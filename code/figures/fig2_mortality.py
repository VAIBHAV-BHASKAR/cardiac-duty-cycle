"""Figure 2: Age-stratified mortality by CDC-deviation tertile."""

import numpy as np
import pandas as pd

from ..config import C_T1, C_T2, C_T3, FIGURES_DIR
from ._helpers import save_fig, new_figure


def plot(merged, decade_results, log=print):
    log("  Plotting Fig 2 — Age-stratified mortality...")

    merged = merged.copy()
    merged["age_decade"] = (merged["age"] // 10).astype(int) * 10
    decades = sorted(d for d in decade_results if decade_results[d]["n"] >= 50 and 20 <= d <= 80)

    fig, ax = new_figure(figsize=(4.7, 3.5))
    ax = ax[0]

    x = np.arange(len(decades))
    width = 0.25

    bar_data = []
    for decade in decades:
        d = merged[merged["age_decade"] == decade]
        cuts = d["CDC_dev"].quantile([1 / 3, 2 / 3]).values
        for tname, mask, color in [
            ("Near 1/$e$", d["CDC_dev"] <= cuts[0], C_T1),
            ("Moderate", (d["CDC_dev"] > cuts[0]) & (d["CDC_dev"] <= cuts[1]), C_T2),
            ("Far from 1/$e$", d["CDC_dev"] > cuts[1], C_T3),
        ]:
            t = d[mask]
            mort = 100 * t["death"].mean()
            se = 100 * np.sqrt(t["death"].mean() * (1 - t["death"].mean()) / len(t)) if len(t) > 0 else 0
            bar_data.append(dict(decade=decade, tertile=tname, mort=mort, se=se, color=color))

    bd = pd.DataFrame(bar_data)

    for i, (tert, color) in enumerate([("Near 1/$e$", C_T1), ("Moderate", C_T2), ("Far from 1/$e$", C_T3)]):
        vals = [bd[(bd["decade"] == d) & (bd["tertile"] == tert)]["mort"].values[0] for d in decades]
        errs = [bd[(bd["decade"] == d) & (bd["tertile"] == tert)]["se"].values[0] for d in decades]
        ax.bar(x + i * width, vals, width, yerr=errs, label=tert, color=color,
               capsize=2, edgecolor="white", linewidth=0.5, error_kw={"linewidth": 0.8})

    for j, decade in enumerate(decades):
        near = bd[(bd["decade"] == decade) & (bd["tertile"] == "Near 1/$e$")]["mort"].values[0]
        far = bd[(bd["decade"] == decade) & (bd["tertile"] == "Far from 1/$e$")]["mort"].values[0]
        if near > 0:
            fold = far / near
            bar_top = far + bd[(bd["decade"] == decade) & (bd["tertile"] == "Far from 1/$e$")]["se"].values[0]
            ax.text(j + width, bar_top + 0.3, f"{fold:.1f}x", ha="center", va="bottom",
                    fontsize=5.5, color="#555")

    ax.set_xlabel("Age decade", fontsize=9)
    ax.set_ylabel("All-cause mortality (%)", fontsize=9)
    ax.set_xticks(x + width)
    ax.set_xticklabels([f"{d}s" for d in decades], fontsize=8)
    ax.legend(fontsize=7, loc="upper left")
    ax.tick_params(labelsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    save_fig(fig, "Fig2_mortality", log=log)
