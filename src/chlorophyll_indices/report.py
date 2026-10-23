"""Figures for the lab + satellite evaluation."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .indices import INDICES

INDEX_ORDER = list(INDICES)
SERIES = (
    ("Lab, within species (mean)", "#9bb7d4"),
    ("Lab, all species pooled", "#2f5f8f"),
    ("Sentinel-2, leaf (LCC)", "#a8c98a"),
    ("Sentinel-2, canopy (CCC)", "#3f7f3a"),
)


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": "#55606a",
                         "axes.labelcolor": "#17212b", "xtick.color": "#3c4650", "ytick.color": "#3c4650"})
    return plt


def r2_comparison(out: Path, within: pd.Series, pooled: pd.DataFrame, sat_lcc: pd.DataFrame, sat_ccc: pd.DataFrame) -> None:
    plt = _plt()
    values = [
        within.reindex(INDEX_ORDER).to_numpy(),
        pooled.set_index("index").r2.reindex(INDEX_ORDER).to_numpy(),
        sat_lcc.set_index("index").r2.reindex(INDEX_ORDER).to_numpy(),
        sat_ccc.set_index("index").r2.reindex(INDEX_ORDER).to_numpy(),
    ]
    fig, ax = plt.subplots(figsize=(9, 4.4), dpi=150)
    x = np.arange(len(INDEX_ORDER))
    width = 0.2
    for i, ((label, colour), v) in enumerate(zip(SERIES, values)):
        bars = ax.bar(x + (i - 1.5) * width, v, width, color=colour, label=label)
        for bar, value in zip(bars, v):
            ax.text(bar.get_x() + bar.get_width() / 2, value + 0.012, f"{value:.2f}", ha="center", fontsize=7)
    ax.set_xticks(x, INDEX_ORDER)
    ax.set_ylim(0, 1)
    ax.set_ylabel("R² of a linear fit")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(fontsize=8, ncol=2, loc="upper left")
    ax.set_title("Red-edge indices (MTCI, IRECI) track chlorophyll best", loc="left", fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out / "r2_comparison.png")
    plt.close(fig)


def lab_scatter(out: Path, lab: pd.DataFrame, pooled: pd.DataFrame) -> None:
    plt = _plt()
    species = sorted(lab.species.unique(), key=lambda s: int(s.split()[-1]))
    colours = plt.get_cmap("tab10")(np.arange(len(species)))
    fig, axes = plt.subplots(1, 5, figsize=(15, 3.4), dpi=150, sharey=True)
    stats = pooled.set_index("index")
    for ax, name in zip(axes, INDEX_ORDER):
        for s, c in zip(species, colours):
            sub = lab[lab.species == s]
            ax.scatter(sub[name], sub.lcc_ug_cm2, s=11, color=c, alpha=0.8, label=s)
        xs = np.linspace(lab[name].min(), lab[name].max(), 20)
        ax.plot(xs, stats.loc[name, "slope"] * xs + stats.loc[name, "intercept"], color="#111", lw=1.4)
        ax.set_title(f"{name} · R² {stats.loc[name, 'r2']:.2f}", fontsize=10, fontweight="bold")
        ax.set_xlabel(name)
        ax.grid(alpha=0.2)
    axes[0].set_ylabel("Leaf chlorophyll (µg/cm²)")
    axes[-1].legend(fontsize=7, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    fig.suptitle("Lab spectra: 160 leaves, 8 species (pooled linear fit)", x=0.01, ha="left", fontsize=12,
                 fontweight="bold")
    fig.tight_layout()
    fig.savefig(out / "lab_scatter.png", bbox_inches="tight")
    plt.close(fig)


def species_heatmap(out: Path, table: pd.DataFrame) -> None:
    plt = _plt()
    table = table[INDEX_ORDER].sort_index(key=lambda s: s.str.split().str[-1].astype(int))
    fig, ax = plt.subplots(figsize=(6.4, 4.6), dpi=150)
    im = ax.imshow(table.to_numpy(), cmap="YlGnBu", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(INDEX_ORDER)), INDEX_ORDER)
    ax.set_yticks(range(len(table)), table.index)
    for (i, j), v in np.ndenumerate(table.to_numpy()):
        ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=8, color="white" if v > 0.6 else "#17212b")
    fig.colorbar(im, ax=ax, shrink=0.85).set_label("R² within species")
    ax.set_title("Index performance varies by species", loc="left", fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out / "species_r2_heatmap.png")
    plt.close(fig)


def satellite_scatter(out: Path, sat: pd.DataFrame, sat_ccc: pd.DataFrame) -> None:
    plt = _plt()
    best = sat_ccc.iloc[0]
    fig, ax = plt.subplots(figsize=(6.4, 4.6), dpi=150)
    sc = ax.scatter(sat[best["index"]], sat.ccc_g_m2, c=sat.lai, cmap="viridis", s=26, edgecolor="white", lw=0.4)
    xs = np.linspace(sat[best["index"]].min(), sat[best["index"]].max(), 20)
    ax.plot(xs, best.slope * xs + best.intercept, color="#111", lw=1.5)
    fig.colorbar(sc, ax=ax, shrink=0.85).set_label("Leaf area index")
    ax.set_xlabel(f"Sentinel-2 {best['index']}")
    ax.set_ylabel("Canopy chlorophyll, CCC (g/m²)")
    ax.grid(alpha=0.25)
    ax.set_title(f"{best['index']} vs canopy chlorophyll · R² {best.r2:.2f} · LOO RMSE {best.loo_rmse:.2f} g/m² · n {int(best.n)}",
                 loc="left", fontsize=10, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out / "satellite_ccc.png")
    plt.close(fig)
