"""Plotting utilities for the RNA Forge 5' Capping Efficiency desktop tool."""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import pandas as pd


def plot_eic(ax, eic_df: pd.DataFrame, species_name: str, rt_window: Optional[Tuple[float, float]] = None) -> None:
    """Draw one extracted ion chromatogram on an existing matplotlib Axes."""
    ax.clear()
    if eic_df is None or len(eic_df) == 0:
        ax.set_title(f"EIC for {species_name}")
        ax.set_xlabel("Retention Time (seconds)")
        ax.set_ylabel("Intensity")
        ax.text(0.5, 0.5, "No EIC points found", transform=ax.transAxes, ha="center", va="center")
        return

    ax.plot(eic_df["RT"], eic_df["Intensity"], marker="o", linewidth=1, markersize=3, label=species_name)
    if rt_window is not None:
        lo, hi = sorted(rt_window)
        ax.axvspan(lo, hi, alpha=0.15)
    ax.set_title(f"EIC for {species_name}")
    ax.set_xlabel("Retention Time (seconds)")
    ax.set_ylabel("Intensity")
    ax.legend(loc="best")
    ax.figure.tight_layout()


def save_eic_tiff(eic_df: pd.DataFrame, species_name: str, output_path: str | Path, rt_window=None, dpi: int = 300) -> None:
    fig, ax = plt.subplots(figsize=(7, 4.5))
    plot_eic(ax, eic_df, species_name, rt_window=rt_window)
    fig.savefig(output_path, dpi=dpi, format="tiff", bbox_inches="tight")
    plt.close(fig)
