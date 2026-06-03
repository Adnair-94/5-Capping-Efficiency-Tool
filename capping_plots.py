"""Plotting utilities for the RNA Forge 5' Capping Efficiency desktop tool."""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import pandas as pd


def _rt_for_display(eic_df: pd.DataFrame, rt_unit: str) -> pd.Series:
    rt = pd.to_numeric(eic_df["RT"], errors="coerce")
    if str(rt_unit).lower().startswith("min"):
        return rt / 60.0
    return rt


def _window_for_display(rt_window: Optional[Tuple[float, float]], rt_unit: str):
    if rt_window is None:
        return None
    lo, hi = sorted([float(rt_window[0]), float(rt_window[1])])
    if str(rt_unit).lower().startswith("min"):
        return lo / 60.0, hi / 60.0
    return lo, hi


def _rt_label(rt_unit: str) -> str:
    return "Retention Time (minutes)" if str(rt_unit).lower().startswith("min") else "Retention Time (seconds)"


def plot_eic(ax, eic_df: pd.DataFrame, species_name: str, rt_window: Optional[Tuple[float, float]] = None, rt_unit: str = "seconds") -> None:
    """Draw one extracted ion chromatogram on an existing matplotlib Axes.

    The stored EIC RT values are always seconds. ``rt_unit`` only controls display.
    """
    ax.clear()
    if eic_df is None or len(eic_df) == 0:
        ax.set_title(f"EIC for {species_name}")
        ax.set_xlabel(_rt_label(rt_unit))
        ax.set_ylabel("Intensity")
        ax.text(0.5, 0.5, "No EIC points found", transform=ax.transAxes, ha="center", va="center")
        return

    x = _rt_for_display(eic_df, rt_unit)
    ax.plot(x, eic_df["Intensity"], marker="o", linewidth=1, markersize=3, label=species_name)
    disp_window = _window_for_display(rt_window, rt_unit)
    if disp_window is not None:
        lo, hi = sorted(disp_window)
        ax.axvspan(lo, hi, alpha=0.15)
    ax.set_title(f"EIC for {species_name}")
    ax.set_xlabel(_rt_label(rt_unit))
    ax.set_ylabel("Intensity")
    ax.legend(loc="best")
    ax.figure.tight_layout()


def save_eic_tiff(eic_df: pd.DataFrame, species_name: str, output_path: str | Path, rt_window=None, dpi: int = 300, rt_unit: str = "seconds") -> None:
    fig, ax = plt.subplots(figsize=(7, 4.5))
    plot_eic(ax, eic_df, species_name, rt_window=rt_window, rt_unit=rt_unit)
    fig.savefig(output_path, dpi=dpi, format="tiff", bbox_inches="tight")
    plt.close(fig)
