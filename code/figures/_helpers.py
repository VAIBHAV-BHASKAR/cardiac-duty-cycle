"""Shared figure helpers."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ..config import FIGURES_DIR


def new_figure(nrows=1, ncols=1, **kwargs):
    fig, axes = plt.subplots(nrows, ncols, **kwargs)
    return fig, axes if hasattr(axes, "__len__") else [axes]


def save_fig(fig, name, log=print):
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    try:
        fig.tight_layout()
    except Exception:
        pass
    fig.savefig(FIGURES_DIR / f"{name}.png", dpi=200, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    log(f"    Saved {name}.png + .pdf")
