"""
Core backend for the RNA Forge 5' Capping Efficiency desktop tool.

This module is a direct Python translation of the working Shiny app.R logic:
- species table presets, custom species generation, CSV/XLSX loading, and validation
- mzML/mzXML parsing through Pyteomics
- MS1 extracted-ion chromatogram generation
- optional MS2 confirmation near MS1 apex
- per-species RT-window AUC calculation by trapezoidal integration
- capping-efficiency calculation

Scientific rule intentionally preserved:
    capping efficiency = sum(AUC_capped) / (sum(AUC_capped) + sum(AUC_uncapped))

Implementation correction relative to the pasted app.R:
    The R app creates a `use_for_efficiency` flag but the final capEff block does not use it.
    This Python version honors that flag by default, because the UI explicitly defines it as
    the denominator/inclusion control. Set honor_use_for_efficiency=False in calculate_efficiency()
    to reproduce the literal R capEff behavior.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Tuple
import math
import os
import re
import sys

import numpy as np
import pandas as pd


PROTON_MASS = 1.007276466812
PHOSPHATE_STEP_MASS = 79.966331
REQUIRED_SPECIES_COLUMNS = ["name", "mz", "type"]
OPTIONAL_SPECIES_COLUMNS = [
    "ms2_mz",
    "mz_tol",
    "ms2_window_sec",
    "charge",
    "neutral_mass",
    "source",
    "use_for_efficiency",
    "notes",
]
SPECIES_COLUMNS = REQUIRED_SPECIES_COLUMNS + OPTIONAL_SPECIES_COLUMNS


@dataclass
class ModeState:
    species_table: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=SPECIES_COLUMNS))
    eic_df: pd.DataFrame = field(default_factory=lambda: empty_eic_df())
    eic_list: Dict[str, pd.DataFrame] = field(default_factory=dict)
    rt_windows: pd.DataFrame = field(default_factory=lambda: blank_windows(None))
    aucs: pd.Series = field(default_factory=lambda: pd.Series(dtype="float64"))


@dataclass
class AnalysisSession:
    ms_file: Optional[str] = None
    all: ModeState = field(default_factory=ModeState)
    filtered: ModeState = field(default_factory=ModeState)
    filter_message: str = ""
    loaded: bool = False

    def current_mode(self, requested_mode: str) -> str:
        if requested_mode == "filtered" and self.filtered.species_table is not None and len(self.filtered.species_table) > 0:
            return "filtered"
        return "all"

    def state(self, requested_mode: str) -> ModeState:
        mode = self.current_mode(requested_mode)
        return self.filtered if mode == "filtered" else self.all


def resource_path(relative_path: str) -> str:
    """Resolve files both during development and inside a PyInstaller onefile bundle."""
    if hasattr(sys, "_MEIPASS"):
        return os.path.join(sys._MEIPASS, relative_path)
    return os.path.join(os.path.abspath("."), relative_path)


# -----------------------------------------------------------------------------
# Species table defaults and generation
# -----------------------------------------------------------------------------

def default_species() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "name": ["m7GmpppAG_1144.143", "m7GmpppAG_1224.113", "m7GmpppAG_611.552", "pppGp_601.9497"],
            "mz": [1144.143, 1224.113, 611.552, 601.9497],
            "type": ["Capped", "Capped", "Capped", "Uncapped"],
            "ms2_mz": [np.nan, np.nan, np.nan, np.nan],
            "mz_tol": [np.nan, np.nan, np.nan, np.nan],
            "ms2_window_sec": [np.nan, np.nan, np.nan, np.nan],
            "charge": [np.nan, np.nan, np.nan, np.nan],
            "neutral_mass": [np.nan, np.nan, np.nan, np.nan],
            "source": ["manual_default"] * 4,
            "use_for_efficiency": [True, True, True, True],
            "notes": ["", "", "", ""],
        }
    )


def mz_from_neutral(neutral_mass: float, charge: int | float = 1) -> float:
    z = float(charge)
    return (float(neutral_mass) - z * PROTON_MASS) / z


def neutral_from_mz(mz: float, charge: int | float = 1) -> float:
    z = float(charge)
    return z * float(mz) + z * PROTON_MASS


def fragment_ppp_neutral_mass(fragment: str) -> float:
    masses = {"AG": 932.0095, "AUG": 1238.0348, "G": 602.9570}
    key = str(fragment).strip().upper()
    return masses.get(key, np.nan)


def fragment_ladder_neutral_masses(ppp_neutral_mass: float) -> Dict[str, float]:
    base = float(ppp_neutral_mass)
    return {
        "ppp": base,
        "pp": base - PHOSPHATE_STEP_MASS,
        "p": base - 2 * PHOSPHATE_STEP_MASS,
    }


def cap_presets() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "preset_id": ["cleancap_ag_gen1_ag", "cleancap_ag_gen2_ag", "cleancap_aug"],
            "label": [
                "CleanCap AG / AG T1 fragment (M capped = 1225.1219)",
                "CleanCap AG 3'-OMe Gen 2 / AG T1 fragment (M capped = 1239.1376)",
                "CleanCap AU/AUG / AUG T1 fragment (M capped = 1531.1472)",
            ],
            "cap_label": ["CleanCap_AG", "CleanCap_AG_3OMe_Gen2", "CleanCap_AUG"],
            "fragment": ["AG", "AG", "AUG"],
            "capped_neutral_mass": [1225.1219, 1239.1376, 1531.1472],
        }
    )


def make_species_row(
    name: str,
    mz: float,
    type_: str,
    charge=np.nan,
    neutral_mass=np.nan,
    source: str = "generated",
    use_for_efficiency: bool = True,
    notes: str = "",
) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "name": [str(name)],
            "mz": [float(mz)],
            "type": [str(type_)],
            "ms2_mz": [np.nan],
            "mz_tol": [np.nan],
            "ms2_window_sec": [np.nan],
            "charge": [charge],
            "neutral_mass": [neutral_mass],
            "source": [str(source)],
            "use_for_efficiency": [bool(use_for_efficiency)],
            "notes": [str(notes)],
        }
    )


def generate_ladder_species(
    fragment: str,
    ppp_neutral_mass: Optional[float] = None,
    include_z2: bool = False,
    include_g_ladder: bool = False,
    g_use_for_efficiency: bool = False,
    source: str = "generated",
) -> pd.DataFrame:
    fragment_clean = str(fragment).strip().upper()
    if ppp_neutral_mass is None or pd.isna(ppp_neutral_mass):
        ppp_neutral_mass = fragment_ppp_neutral_mass(fragment_clean)
    if pd.isna(ppp_neutral_mass):
        raise ValueError(
            "No triphosphate neutral mass is available for this fragment. "
            "For a custom fragment, provide the z=1 ppp-fragment m/z or neutral mass."
        )

    rows: List[pd.DataFrame] = []
    ladder = fragment_ladder_neutral_masses(ppp_neutral_mass)
    for state, neutral_mass in ladder.items():
        rows.append(
            make_species_row(
                name=f"5' {state}{fragment_clean}_z1",
                mz=mz_from_neutral(neutral_mass, 1),
                type_="Uncapped",
                charge=1,
                neutral_mass=neutral_mass,
                source=source,
                use_for_efficiency=True,
                notes=f"Generated uncapped {state} ladder for {fragment_clean} T1 fragment",
            )
        )
        if include_z2:
            rows.append(
                make_species_row(
                    name=f"5' {state}{fragment_clean}_z2",
                    mz=mz_from_neutral(neutral_mass, 2),
                    type_="Uncapped",
                    charge=2,
                    neutral_mass=neutral_mass,
                    source=source,
                    use_for_efficiency=True,
                    notes=f"Generated z=2 uncapped {state} ladder for {fragment_clean} T1 fragment",
                )
            )

    if include_g_ladder and fragment_clean != "G":
        g_ladder = fragment_ladder_neutral_masses(fragment_ppp_neutral_mass("G"))
        for state, neutral_mass in g_ladder.items():
            rows.append(
                make_species_row(
                    name=f"5' {state}G_z1",
                    mz=mz_from_neutral(neutral_mass, 1),
                    type_="Uncapped",
                    charge=1,
                    neutral_mass=neutral_mass,
                    source=source,
                    use_for_efficiency=bool(g_use_for_efficiency),
                    notes=(
                        "Generated G-only ladder; included in capping-efficiency denominator by user choice"
                        if g_use_for_efficiency
                        else "Generated G-only ladder; diagnostic only by default, excluded from denominator"
                    ),
                )
            )
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=SPECIES_COLUMNS)


def generate_preset_species(
    preset_id: str,
    include_capped_z2: bool = True,
    include_uncapped_z2: bool = False,
    include_ppp_pp_p: bool = True,
    include_g_ladder: bool = True,
    g_use_for_efficiency: bool = False,
) -> pd.DataFrame:
    presets = cap_presets()
    match = presets[presets["preset_id"] == preset_id]
    if len(match) != 1:
        raise ValueError("Unknown preset selected.")
    row = match.iloc[0]
    src = f"preset:{row['preset_id']}"
    out: List[pd.DataFrame] = []

    out.append(
        make_species_row(
            name=f"{row['cap_label']}_{row['fragment']}_capped_z1",
            mz=mz_from_neutral(row["capped_neutral_mass"], 1),
            type_="Capped",
            charge=1,
            neutral_mass=float(row["capped_neutral_mass"]),
            source=src,
            use_for_efficiency=True,
            notes=f"Preset capped species; cap remains attached to {row['fragment']} T1 fragment",
        )
    )
    if include_capped_z2:
        out.append(
            make_species_row(
                name=f"{row['cap_label']}_{row['fragment']}_capped_z2",
                mz=mz_from_neutral(row["capped_neutral_mass"], 2),
                type_="Capped",
                charge=2,
                neutral_mass=float(row["capped_neutral_mass"]),
                source=src,
                use_for_efficiency=True,
                notes=f"Preset capped z=2 species; cap remains attached to {row['fragment']} T1 fragment",
            )
        )
    if include_ppp_pp_p:
        out.append(
            generate_ladder_species(
                fragment=row["fragment"],
                ppp_neutral_mass=fragment_ppp_neutral_mass(row["fragment"]),
                include_z2=include_uncapped_z2,
                include_g_ladder=include_g_ladder,
                g_use_for_efficiency=g_use_for_efficiency,
                source=src,
            )
        )
    return pd.concat(out, ignore_index=True)


def generate_custom_species(
    custom_cap_name: str,
    capped_mz: float,
    capped_charge: int = 1,
    fragment: str = "AG",
    custom_ppp_mz: Optional[float] = None,
    custom_ppp_charge: int = 1,
    include_capped_z2: bool = True,
    include_uncapped_z2: bool = False,
    include_g_ladder: bool = False,
    g_use_for_efficiency: bool = False,
) -> pd.DataFrame:
    if pd.isna(capped_mz) or float(capped_mz) <= 0:
        raise ValueError("Enter a valid capped m/z for the custom cap species.")
    cap_z = int(capped_charge)
    capped_neutral = neutral_from_mz(capped_mz, cap_z)
    label = re.sub(r"[^A-Za-z0-9_]+", "_", str(custom_cap_name).strip()) or "CustomCap"
    frag = str(fragment).strip().upper()

    if frag == "CUSTOM":
        if custom_ppp_mz is None or pd.isna(custom_ppp_mz) or float(custom_ppp_mz) <= 0:
            raise ValueError("For a custom T1 fragment, enter the ppp-fragment m/z.")
        ppp_neutral = neutral_from_mz(custom_ppp_mz, custom_ppp_charge)
        frag_label = "Custom"
    else:
        ppp_neutral = fragment_ppp_neutral_mass(frag)
        frag_label = frag
        if pd.isna(ppp_neutral):
            raise ValueError("Unknown fragment. Use AG, AUG, G, or Custom.")

    out: List[pd.DataFrame] = [
        make_species_row(
            name=f"{label}_{frag_label}_capped_z{cap_z}",
            mz=float(capped_mz),
            type_="Capped",
            charge=cap_z,
            neutral_mass=capped_neutral,
            source="custom_generator",
            use_for_efficiency=True,
            notes="Custom capped species. Capped neutral mass calculated from m/z and z. The cap m7G is not treated as the RNase T1 cleavage G.",
        )
    ]
    if include_capped_z2 and cap_z != 2:
        out.append(
            make_species_row(
                name=f"{label}_{frag_label}_capped_z2",
                mz=mz_from_neutral(capped_neutral, 2),
                type_="Capped",
                charge=2,
                neutral_mass=capped_neutral,
                source="custom_generator",
                use_for_efficiency=True,
                notes="Calculated z=2 capped ion from custom capped neutral mass",
            )
        )
    out.append(
        generate_ladder_species(
            fragment=frag_label,
            ppp_neutral_mass=ppp_neutral,
            include_z2=include_uncapped_z2,
            include_g_ladder=include_g_ladder,
            g_use_for_efficiency=g_use_for_efficiency,
            source="custom_generator",
        )
    )
    return pd.concat(out, ignore_index=True)


# -----------------------------------------------------------------------------
# Species table validation
# -----------------------------------------------------------------------------

def normalize_type(values: Iterable) -> List[Optional[str]]:
    out: List[Optional[str]] = []
    for val in values:
        z = str(val).strip().lower()
        if z in {"capped", "cap", "clean cap", "cleancap"}:
            out.append("Capped")
        elif z in {"uncapped", "un-capped", "un capped", "un_cap", "uncap"}:
            out.append("Uncapped")
        else:
            out.append(None)
    return out


def as_num_safely(series) -> pd.Series:
    return pd.to_numeric(pd.Series(series), errors="coerce")


def as_logical_safely(series, default: bool = True) -> pd.Series:
    if series is None:
        return pd.Series(dtype=bool)
    out = []
    for val in series:
        z = str(val).strip().lower()
        if z in {"true", "t", "yes", "y", "1", "include", "included"}:
            out.append(True)
        elif z in {"false", "f", "no", "n", "0", "exclude", "excluded", "diagnostic"}:
            out.append(False)
        elif z in {"", "nan", "none", "na", "n/a"}:
            out.append(default)
        else:
            out.append(default)
    return pd.Series(out, dtype=bool)


def _first_col(df: pd.DataFrame, col: str):
    cols = [c for c in df.columns if c == col]
    if not cols:
        return None
    return df[cols[0]]


def validate_species_table(tbl: pd.DataFrame) -> Tuple[bool, pd.DataFrame, str]:
    """Validate and canonicalize a species table.

    Returns (ok, canonical_table, messages).
    """
    if tbl is None or len(tbl) == 0:
        return False, pd.DataFrame(columns=SPECIES_COLUMNS), "Species table is empty."

    df = pd.DataFrame(tbl).copy()
    names_raw = list(df.columns)
    names_lc = [str(c).strip().lower() for c in names_raw]

    canonical_map = {
        "name": "name",
        "species": "name",
        "species_name": "name",
        "mz": "mz",
        "m/z": "mz",
        "mass_to_charge": "mz",
        "type": "type",
        "class": "type",
        "category": "type",
        "ms2_mz": "ms2_mz",
        "ms2mz": "ms2_mz",
        "confirmation_mz": "ms2_mz",
        "diagnostic_mz": "ms2_mz",
        "mz_tol": "mz_tol",
        "mztol": "mz_tol",
        "tolerance": "mz_tol",
        "ms2_window_sec": "ms2_window_sec",
        "ms2windowsec": "ms2_window_sec",
        "ms2_window": "ms2_window_sec",
        "charge": "charge",
        "z": "charge",
        "neutral_mass": "neutral_mass",
        "neutralmass": "neutral_mass",
        "m": "neutral_mass",
        "source": "source",
        "preset": "source",
        "use_for_efficiency": "use_for_efficiency",
        "include_in_efficiency": "use_for_efficiency",
        "denominator": "use_for_efficiency",
        "notes": "notes",
        "note": "notes",
    }

    canonical_names = [canonical_map.get(c, names_raw[i]) for i, c in enumerate(names_lc)]
    if len(set(canonical_names)) != len(canonical_names):
        counts = {}
        unique_names = []
        for c in canonical_names:
            counts[c] = counts.get(c, 0) + 1
            unique_names.append(c if counts[c] == 1 else f"{c}_{counts[c]}")
        canonical_names = unique_names
    df.columns = canonical_names

    missing = [c for c in REQUIRED_SPECIES_COLUMNS if c not in df.columns]
    if missing:
        return (
            False,
            pd.DataFrame(columns=SPECIES_COLUMNS),
            "Missing required column(s): " + ", ".join(missing) + ". Required columns are name, mz, and type.",
        )

    out = pd.DataFrame(
        {
            "name": _first_col(df, "name").astype(str).str.strip(),
            "mz": as_num_safely(_first_col(df, "mz")),
            "type": normalize_type(_first_col(df, "type")),
        }
    )

    optional_defaults = {
        "ms2_mz": np.nan,
        "mz_tol": np.nan,
        "ms2_window_sec": np.nan,
        "charge": np.nan,
        "neutral_mass": np.nan,
        "source": "manual_or_uploaded",
        "use_for_efficiency": True,
        "notes": "",
    }
    for col, default in optional_defaults.items():
        val = _first_col(df, col)
        if val is None:
            out[col] = default
        elif col in {"ms2_mz", "mz_tol", "ms2_window_sec", "neutral_mass"}:
            out[col] = as_num_safely(val)
        elif col == "use_for_efficiency":
            out[col] = as_logical_safely(val, default=True)
        else:
            out[col] = val.astype(str)

    messages: List[str] = []
    bad_name = out.index[out["name"].isna() | (out["name"].str.strip() == "")].tolist()
    bad_mz = out.index[out["mz"].isna() | ~np.isfinite(out["mz"])].tolist()
    bad_type = out.index[out["type"].isna()].tolist()
    bad_tol = out.index[out["mz_tol"].notna() & ((~np.isfinite(out["mz_tol"])) | (out["mz_tol"] <= 0))].tolist()
    bad_window = out.index[
        out["ms2_window_sec"].notna() & ((~np.isfinite(out["ms2_window_sec"])) | (out["ms2_window_sec"] <= 0))
    ].tolist()

    if bad_name:
        messages.append("Rows with missing species names: " + ", ".join(str(i + 1) for i in bad_name))
    if bad_mz:
        messages.append("Rows with missing/non-numeric m/z: " + ", ".join(str(i + 1) for i in bad_mz))
    if bad_type:
        messages.append("Rows with invalid type; use Capped or Uncapped: " + ", ".join(str(i + 1) for i in bad_type))
    if bad_tol:
        messages.append("Rows with invalid mz_tol; use a positive number or leave blank: " + ", ".join(str(i + 1) for i in bad_tol))
    if bad_window:
        messages.append(
            "Rows with invalid ms2_window_sec; use a positive number or leave blank: " + ", ".join(str(i + 1) for i in bad_window)
        )

    if messages:
        return False, out[SPECIES_COLUMNS], "\n".join(messages)

    if out["name"].duplicated().any():
        original = out["name"].copy()
        seen = {}
        unique = []
        for name in original:
            seen[name] = seen.get(name, 0) + 1
            unique.append(name if seen[name] == 1 else f"{name}_{seen[name]}")
        out["name"] = unique
        messages.append("Duplicate species names were made unique automatically.")

    if not (out["type"] == "Capped").any():
        messages.append("Warning: no Capped species are present.")
    if not (out["type"] == "Uncapped").any():
        messages.append("Warning: no Uncapped species are present.")

    return True, out[SPECIES_COLUMNS], "\n".join(messages)


def read_species_file(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    ext = path.suffix.lower().lstrip(".")
    if ext == "csv":
        return pd.read_csv(path)
    if ext in {"xlsx", "xls"}:
        return pd.read_excel(path)
    raise ValueError("Unsupported species-table file type. Upload .csv, .xlsx, or .xls.")


def species_template() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "name": [
                "CleanCap_AG_AG_capped_z1",
                "CleanCap_AG_AG_capped_z2",
                "5' pppAG_z1",
                "5' ppAG_z1",
                "5' pAG_z1",
                "5' pppG_z1",
            ],
            "mz": [1224.1146, 611.5537, 931.0022, 851.0359, 771.0696, 601.9497],
            "type": ["Capped", "Capped", "Uncapped", "Uncapped", "Uncapped", "Uncapped"],
            "ms2_mz": [np.nan] * 6,
            "mz_tol": [np.nan] * 6,
            "ms2_window_sec": [np.nan] * 6,
            "charge": [1, 2, 1, 1, 1, 1],
            "neutral_mass": [1225.1219, 1225.1219, 932.0095, 852.0432, 772.0768, 602.9570],
            "source": ["template"] * 6,
            "use_for_efficiency": [True, True, True, True, True, False],
            "notes": [
                "Capped AG T1 fragment; cap remains attached after RNase T1 digestion",
                "Capped AG T1 fragment z=2",
                "Uncapped triphosphate AG ladder",
                "Uncapped diphosphate AG ladder",
                "Uncapped monophosphate AG ladder",
                "G-only ladder; diagnostic by default unless intentionally included in denominator",
            ],
        }
    )


def write_species_template(path: str | Path) -> None:
    species_template().to_csv(path, index=False)


def blank_windows(species_tbl: Optional[pd.DataFrame]) -> pd.DataFrame:
    if species_tbl is None or len(species_tbl) == 0:
        return pd.DataFrame({"Species": pd.Series(dtype=str), "RTmin": pd.Series(dtype=float), "RTmax": pd.Series(dtype=float)})
    return pd.DataFrame({"Species": species_tbl["name"].values, "RTmin": np.nan, "RTmax": np.nan})


def empty_eic_df() -> pd.DataFrame:
    return pd.DataFrame({"RT": pd.Series(dtype=float), "Intensity": pd.Series(dtype=float), "Species": pd.Series(dtype=str)})


# -----------------------------------------------------------------------------
# mzML/mzXML reading and EIC extraction
# -----------------------------------------------------------------------------

def _reader_for_file(ms_file: str | Path):
    ext = Path(ms_file).suffix.lower()
    try:
        if ext == ".mzml":
            from pyteomics import mzml

            return mzml.MzML(str(ms_file))
        if ext == ".mzxml":
            from pyteomics import mzxml

            return mzxml.MzXML(str(ms_file))
    except ImportError as exc:
        raise ImportError(
            "Pyteomics is required for mzML/mzXML parsing. Install with `pip install pyteomics`, "
            "or use the provided requirements.txt before building the executable."
        ) from exc
    raise ValueError("Unsupported file type. Please select a .mzML or .mzXML file.")


def _get_ms_level(spectrum: dict) -> Optional[int]:
    for key in ("ms level", "msLevel"):
        if key in spectrum:
            try:
                return int(spectrum[key])
            except Exception:
                pass
    # mzXML often stores msLevel under params-like keys; try a permissive fallback.
    for key, value in spectrum.items():
        if str(key).lower().replace(" ", "") == "mslevel":
            try:
                return int(value)
            except Exception:
                return None
    return None


def _get_rt_seconds(spectrum: dict) -> float:
    candidates = []
    scan_list = spectrum.get("scanList") or spectrum.get("scan list")
    if isinstance(scan_list, dict):
        scans = scan_list.get("scan") or []
        if isinstance(scans, list) and scans:
            candidates.append(scans[0])
        elif isinstance(scans, dict):
            candidates.append(scans)
    candidates.append(spectrum)

    for container in candidates:
        if not isinstance(container, dict):
            continue
        for key in ("scan start time", "retentionTime", "retention time", "startTime"):
            if key in container:
                val = container[key]
                unit = str(getattr(val, "unit_info", "") or container.get("unitName", "") or container.get("unit", "")).lower()
                try:
                    rt = float(val)
                except Exception:
                    continue
                if "minute" in unit or unit == "min":
                    return rt * 60.0
                return rt
    return np.nan


def _sum_intensity_in_window(mz_array, intensity_array, target_mz: float, tol: float) -> float:
    if mz_array is None or intensity_array is None:
        return 0.0
    mzs = np.asarray(mz_array, dtype=float)
    ints = np.asarray(intensity_array, dtype=float)
    if mzs.size == 0 or ints.size == 0:
        return 0.0
    mask = (mzs >= target_mz - tol) & (mzs <= target_mz + tol)
    if not np.any(mask):
        return 0.0
    return float(np.nansum(ints[mask]))


def _get_tol_for_row(sdt: pd.DataFrame, i: int, default_mz_tol: float) -> float:
    try:
        row_tol = float(sdt.iloc[i].get("mz_tol", np.nan))
    except Exception:
        row_tol = np.nan
    return row_tol if np.isfinite(row_tol) and row_tol > 0 else float(default_mz_tol)


def _get_ms2_window_for_row(sdt: pd.DataFrame, i: int, default_ms2_window: float) -> float:
    try:
        row_window = float(sdt.iloc[i].get("ms2_window_sec", np.nan))
    except Exception:
        row_window = np.nan
    return row_window if np.isfinite(row_window) and row_window > 0 else float(default_ms2_window)


def _get_ms2_mz_for_row(sdt: pd.DataFrame, i: int) -> float:
    try:
        row_ms2 = float(sdt.iloc[i].get("ms2_mz", np.nan))
    except Exception:
        row_ms2 = np.nan
    return row_ms2 if np.isfinite(row_ms2) and row_ms2 > 0 else float(sdt.iloc[i]["mz"])


def _get_apex_rt(eic: pd.DataFrame) -> float:
    if eic is None or len(eic) == 0:
        return np.nan
    ints = pd.to_numeric(eic["Intensity"], errors="coerce").to_numpy(dtype=float)
    if ints.size == 0 or np.all(np.isnan(ints)) or np.nanmax(ints) <= 0:
        return np.nan
    idx = int(np.nanargmax(ints))
    return float(eic.iloc[idx]["RT"])


def _has_ms2_signal_near_apex(ms2_eic: pd.DataFrame, apex_rt: float, window_sec: float) -> bool:
    if ms2_eic is None or len(ms2_eic) == 0 or not np.isfinite(apex_rt):
        return False
    rt = pd.to_numeric(ms2_eic["RT"], errors="coerce").to_numpy(dtype=float)
    ints = pd.to_numeric(ms2_eic["Intensity"], errors="coerce").to_numpy(dtype=float)
    return bool(np.any((np.abs(rt - apex_rt) <= float(window_sec)) & (ints > 0)))


def build_analysis(
    ms_file: str | Path,
    species_tbl: pd.DataFrame,
    default_mz_tol: float = 0.01,
    default_ms2_window: float = 30.0,
    progress_cb: Optional[Callable[[float, str], None]] = None,
) -> Tuple[AnalysisSession, str]:
    """Load mzML/mzXML, extract EICs, and build all/filtered analysis states.

    progress_cb receives (fraction_0_to_1, status_text).
    """
    ok, sdt, msg = validate_species_table(species_tbl)
    if not ok:
        raise ValueError(msg)

    ms_file = str(ms_file)
    n_species = len(sdt)
    ms1_rts: List[float] = []
    ms1_intensities = [list() for _ in range(n_species)]
    ms2_rts: List[float] = []
    ms2_intensities = [list() for _ in range(n_species)]

    if progress_cb:
        progress_cb(0.0, "Opening MS file...")

    with _reader_for_file(ms_file) as reader:
        # Pyteomics readers are iterators; total scan count is not always cheap/available.
        # Progress is therefore approximate while spectra are being streamed.
        for scan_idx, spectrum in enumerate(reader, start=1):
            ms_level = _get_ms_level(spectrum)
            if ms_level not in {1, 2}:
                continue
            rt = _get_rt_seconds(spectrum)
            if not np.isfinite(rt):
                continue
            mz_array = spectrum.get("m/z array")
            intensity_array = spectrum.get("intensity array")

            if ms_level == 1:
                ms1_rts.append(rt)
                for i in range(n_species):
                    tol_i = _get_tol_for_row(sdt, i, default_mz_tol)
                    ms1_intensities[i].append(_sum_intensity_in_window(mz_array, intensity_array, float(sdt.iloc[i]["mz"]), tol_i))
            elif ms_level == 2:
                ms2_rts.append(rt)
                for i in range(n_species):
                    tol_i = _get_tol_for_row(sdt, i, default_mz_tol)
                    ms2_target = _get_ms2_mz_for_row(sdt, i)
                    ms2_intensities[i].append(_sum_intensity_in_window(mz_array, intensity_array, ms2_target, tol_i))

            if progress_cb and scan_idx % 250 == 0:
                progress_cb(0.10, f"Scanned {scan_idx:,} spectra...")

    eic_list_all: Dict[str, pd.DataFrame] = {}
    eic_df_parts: List[pd.DataFrame] = []
    ms2_eic_list: Dict[str, pd.DataFrame] = {}

    for i in range(n_species):
        name = str(sdt.iloc[i]["name"])
        eic = pd.DataFrame({"RT": ms1_rts, "Intensity": ms1_intensities[i], "Species": name})
        eic = eic.sort_values("RT", kind="mergesort").reset_index(drop=True)
        eic_list_all[name] = eic
        eic_df_parts.append(eic)

        ms2_eic = pd.DataFrame({"RT": ms2_rts, "Intensity": ms2_intensities[i], "Species": name})
        ms2_eic = ms2_eic.sort_values("RT", kind="mergesort").reset_index(drop=True)
        ms2_eic_list[name] = ms2_eic

        if progress_cb:
            progress_cb(0.45 + 0.25 * ((i + 1) / max(n_species, 1)), f"Built EIC {i + 1}/{n_species}")

    eic_df_all = pd.concat(eic_df_parts, ignore_index=True) if eic_df_parts else empty_eic_df()

    confirmed_idx: List[int] = []
    for i in range(n_species):
        name = str(sdt.iloc[i]["name"])
        apex_rt = _get_apex_rt(eic_list_all[name])
        window_i = _get_ms2_window_for_row(sdt, i, default_ms2_window)
        if _has_ms2_signal_near_apex(ms2_eic_list[name], apex_rt, window_i):
            confirmed_idx.append(i)
        if progress_cb:
            progress_cb(0.70 + 0.30 * ((i + 1) / max(n_species, 1)), f"MS2 check {i + 1}/{n_species}")

    all_state = ModeState(
        species_table=sdt.reset_index(drop=True),
        eic_df=eic_df_all,
        eic_list=eic_list_all,
        rt_windows=blank_windows(sdt),
        aucs=pd.Series([np.nan] * n_species, dtype="float64"),
    )

    if len(confirmed_idx) == 0:
        filtered_table = pd.DataFrame(columns=SPECIES_COLUMNS)
        filtered_state = ModeState(
            species_table=filtered_table,
            eic_df=empty_eic_df(),
            eic_list={},
            rt_windows=blank_windows(filtered_table),
            aucs=pd.Series(dtype="float64"),
        )
        filter_message = "No MS2-confirmed species found. Filtered mode will fall back to All Species mode."
    else:
        filtered_table = sdt.iloc[confirmed_idx].reset_index(drop=True)
        filtered_eic_list = {str(sdt.iloc[i]["name"]): eic_list_all[str(sdt.iloc[i]["name"])] for i in confirmed_idx}
        filtered_eic_df = pd.concat(filtered_eic_list.values(), ignore_index=True) if filtered_eic_list else empty_eic_df()
        filtered_state = ModeState(
            species_table=filtered_table,
            eic_df=filtered_eic_df,
            eic_list=filtered_eic_list,
            rt_windows=blank_windows(filtered_table),
            aucs=pd.Series([np.nan] * len(filtered_table), dtype="float64"),
        )
        names = ", ".join(filtered_table["name"].astype(str).tolist())
        filter_message = f"MS2-confirmed: {len(confirmed_idx)} of {n_species} species: {names}"

    session = AnalysisSession(ms_file=ms_file, all=all_state, filtered=filtered_state, filter_message=filter_message, loaded=True)
    return session, msg


def integrate_eic_window(eic: pd.DataFrame, rt_min: float, rt_max: float) -> float:
    if eic is None or len(eic) == 0:
        return np.nan
    lo, hi = sorted([float(rt_min), float(rt_max)])
    rt = pd.to_numeric(eic["RT"], errors="coerce").to_numpy(dtype=float)
    ints = pd.to_numeric(eic["Intensity"], errors="coerce").to_numpy(dtype=float)
    idx = np.where((rt >= lo) & (rt <= hi) & (ints > 0))[0]
    if len(idx) < 2:
        return np.nan
    # np.trapezoid is used where available; np.trapz fallback keeps compatibility.
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(ints[idx], rt[idx]))
    return float(np.trapz(ints[idx], rt[idx]))


def save_rt_window(session: AnalysisSession, requested_mode: str, species_name: str, rt_min: float, rt_max: float) -> float:
    state = session.state(requested_mode)
    if state.species_table is None or len(state.species_table) == 0:
        raise ValueError("No species available in this mode.")
    matches = state.species_table.index[state.species_table["name"].astype(str) == str(species_name)].tolist()
    if len(matches) != 1:
        raise ValueError("Selected species was not found in the current table.")
    idx = matches[0]
    eic = state.eic_list.get(str(species_name))
    auc = integrate_eic_window(eic, rt_min, rt_max)
    state.rt_windows.loc[idx, "RTmin"] = min(float(rt_min), float(rt_max))
    state.rt_windows.loc[idx, "RTmax"] = max(float(rt_min), float(rt_max))
    state.aucs.iloc[idx] = auc
    return auc


def auc_table(session: AnalysisSession, requested_mode: str) -> pd.DataFrame:
    state = session.state(requested_mode)
    ct = state.species_table
    if ct is None or len(ct) == 0:
        return pd.DataFrame()
    rt_tbl = state.rt_windows.reset_index(drop=True)
    out = pd.DataFrame(
        {
            "Species": ct["name"].values,
            "mz": pd.to_numeric(ct["mz"], errors="coerce").values,
            "charge": ct["charge"].values if "charge" in ct.columns else np.nan,
            "neutral_mass": pd.to_numeric(ct["neutral_mass"], errors="coerce").values if "neutral_mass" in ct.columns else np.nan,
            "type": ct["type"].values,
            "use_for_efficiency": ct["use_for_efficiency"].values if "use_for_efficiency" in ct.columns else True,
            "source": ct["source"].values if "source" in ct.columns else "",
            "RTmin": rt_tbl["RTmin"].values if "RTmin" in rt_tbl.columns else np.nan,
            "RTmax": rt_tbl["RTmax"].values if "RTmax" in rt_tbl.columns else np.nan,
            "AUC": state.aucs.values,
            "notes": ct["notes"].values if "notes" in ct.columns else "",
        }
    )
    return out


def calculate_efficiency(auc_tbl: pd.DataFrame, honor_use_for_efficiency: bool = True) -> Dict[str, float | str | bool]:
    if auc_tbl is None or len(auc_tbl) == 0:
        return {"ok": False, "message": "No species available in this mode."}

    df = auc_tbl.copy()
    if honor_use_for_efficiency and "use_for_efficiency" in df.columns:
        include = as_logical_safely(df["use_for_efficiency"], default=True).values
        df = df.loc[include].copy()

    if not (df["type"] == "Capped").any() or not (df["type"] == "Uncapped").any():
        return {"ok": False, "message": "Define at least one included Capped and one included Uncapped species."}

    capped_aucs = pd.to_numeric(df.loc[df["type"] == "Capped", "AUC"], errors="coerce")
    uncapped_aucs = pd.to_numeric(df.loc[df["type"] == "Uncapped", "AUC"], errors="coerce")

    if capped_aucs.isna().any() or uncapped_aucs.isna().any():
        return {
            "ok": False,
            "message": "Capping efficiency could not be calculated. Save valid RT windows/AUCs for all included Capped and Uncapped species in the current mode.",
        }
    capped_sum = float(capped_aucs.sum())
    uncapped_sum = float(uncapped_aucs.sum())
    denom = capped_sum + uncapped_sum
    if denom <= 0:
        return {"ok": False, "message": "Capping efficiency could not be calculated because total included AUC is zero."}
    eff = capped_sum / denom
    return {
        "ok": True,
        "capping_efficiency_fraction": eff,
        "capping_efficiency_percent": 100.0 * eff,
        "capped_auc_sum": capped_sum,
        "uncapped_auc_sum": uncapped_sum,
        "included_species_count": int(len(df)),
        "message": f"Capping efficiency: {100.0 * eff:.2f}%\nCapped AUC sum: {capped_sum:.4g}\nUncapped AUC sum: {uncapped_sum:.4g}",
    }
