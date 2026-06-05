"""
Core backend for the RNA Forge 5' Capping Efficiency desktop tool.

This module is a direct Python translation of the working Shiny app.R logic:
- species table presets, custom species generation, CSV/XLSX loading, and validation
- mzML/mzXML parsing through a bundled pure-Python XML/base64 parser
- MS1 extracted-ion chromatogram generation
- optional MS1 and MS2 confirmation near MS1 apex
- per-species intensity-maxima reporting
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
import base64
import math
import os
import re
import sys
import xml.etree.ElementTree as ET
import zlib

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
    intensity_summary: pd.DataFrame = field(default_factory=lambda: pd.DataFrame())


@dataclass
class AnalysisSession:
    ms_file: Optional[str] = None
    all: ModeState = field(default_factory=ModeState)
    filtered: ModeState = field(default_factory=ModeState)
    filter_message: str = ""
    ms1_message: str = ""
    uncapped_intensity_message: str = ""
    labile_loss_message: str = ""
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
    cols = {
        "Species": pd.Series(dtype=str),
        "RTmin": pd.Series(dtype=float),
        "RTmax": pd.Series(dtype=float),
        "RTmin_min": pd.Series(dtype=float),
        "RTmax_min": pd.Series(dtype=float),
        "AUC_min_intensity": pd.Series(dtype=float),
        "AUC_min_consecutive_points": pd.Series(dtype=float),
        "AUC_points_used": pd.Series(dtype=float),
        "AUC_segments_used": pd.Series(dtype=float),
        "AUC_max_consecutive_points": pd.Series(dtype=float),
        "AUC_peak_passed_min_points": pd.Series(dtype=object),
    }
    if species_tbl is None or len(species_tbl) == 0:
        return pd.DataFrame(cols)
    n = len(species_tbl)
    return pd.DataFrame({
        "Species": species_tbl["name"].values,
        "RTmin": np.nan,
        "RTmax": np.nan,
        "RTmin_min": np.nan,
        "RTmax_min": np.nan,
        "AUC_min_intensity": np.nan,
        "AUC_min_consecutive_points": np.nan,
        "AUC_points_used": np.nan,
        "AUC_segments_used": np.nan,
        "AUC_max_consecutive_points": np.nan,
        "AUC_peak_passed_min_points": pd.Series([pd.NA] * n, dtype="object"),
    })


def empty_eic_df() -> pd.DataFrame:
    return pd.DataFrame({"RT": pd.Series(dtype=float), "Intensity": pd.Series(dtype=float), "Species": pd.Series(dtype=str)})


# -----------------------------------------------------------------------------
# mzML/mzXML reading and EIC extraction
# -----------------------------------------------------------------------------

# This implementation intentionally avoids pyteomics. The previous build depended on
# PyInstaller correctly bundling pyteomics and failed on some Windows builds. For this
# app, we only need a narrow subset of mzML/mzXML functionality: scan RT, MS level,
# m/z array, and intensity array. The parser below streams XML and decodes standard
# base64 binary arrays directly, which makes the executable more transferable.


def _local_name(tag: str) -> str:
    """Return the XML local tag name without namespace."""
    if "}" in tag:
        return tag.rsplit("}", 1)[1]
    return tag


def _cv_params(elem) -> List[Tuple[str, str, str, str]]:
    """Collect cvParam tuples: accession, name, value, unitName."""
    out: List[Tuple[str, str, str, str]] = []
    for child in elem.iter():
        if _local_name(child.tag) == "cvParam":
            out.append(
                (
                    child.attrib.get("accession", ""),
                    child.attrib.get("name", ""),
                    child.attrib.get("value", ""),
                    child.attrib.get("unitName", "") or child.attrib.get("unitAccession", ""),
                )
            )
    return out


def _has_cv(elem, accession: str | None = None, name: str | None = None) -> bool:
    name_l = name.lower() if name else None
    for acc, nm, _val, _unit in _cv_params(elem):
        if accession and acc == accession:
            return True
        if name_l and nm.lower() == name_l:
            return True
    return False


def _cv_value(elem, accessions: Iterable[str] = (), names: Iterable[str] = ()) -> Tuple[Optional[str], Optional[str]]:
    accessions = set(accessions)
    names_l = {n.lower() for n in names}
    for acc, nm, val, unit in _cv_params(elem):
        if (accessions and acc in accessions) or (names_l and nm.lower() in names_l):
            return val, unit
    return None, None


def _parse_iso_duration_to_seconds(value: str | None) -> float:
    """Parse mzXML-style retentionTime values such as PT123.4S or PT5.2M."""
    if value is None:
        return np.nan
    s = str(value).strip()
    if not s:
        return np.nan
    try:
        return float(s)
    except Exception:
        pass
    m = re.match(r"^PT(?:(?P<h>[0-9.]+)H)?(?:(?P<m>[0-9.]+)M)?(?:(?P<s>[0-9.]+)S)?$", s, flags=re.I)
    if not m:
        return np.nan
    total = 0.0
    if m.group("h"):
        total += float(m.group("h")) * 3600.0
    if m.group("m"):
        total += float(m.group("m")) * 60.0
    if m.group("s"):
        total += float(m.group("s"))
    return total


def _decode_base64_array(binary_text: str | None, dtype: np.dtype, compressed: bool, byte_order: str = "little") -> np.ndarray:
    if binary_text is None:
        return np.array([], dtype=float)
    btxt = "".join(str(binary_text).split())
    if not btxt:
        return np.array([], dtype=float)
    raw = base64.b64decode(btxt)
    if compressed:
        raw = zlib.decompress(raw)
    endian = ">" if str(byte_order).lower() in {"network", "big", "bigendian", "big-endian"} else "<"
    dt = np.dtype(dtype).newbyteorder(endian)
    if len(raw) == 0:
        return np.array([], dtype=float)
    return np.frombuffer(raw, dtype=dt).astype(float, copy=False)


def _mzml_ms_level(spectrum_elem) -> Optional[int]:
    val, _unit = _cv_value(spectrum_elem, accessions={"MS:1000511"}, names={"ms level"})
    if val is None:
        return None
    try:
        return int(float(val))
    except Exception:
        return None


def _mzml_rt_seconds(spectrum_elem) -> float:
    # Restrict to scan elements when possible to avoid unrelated timestamps.
    for elem in spectrum_elem.iter():
        if _local_name(elem.tag) != "scan":
            continue
        val, unit = _cv_value(elem, accessions={"MS:1000016"}, names={"scan start time"})
        if val is None:
            continue
        try:
            rt = float(val)
        except Exception:
            continue
        unit_l = str(unit or "").lower()
        if "minute" in unit_l or unit_l.endswith("0031"):
            return rt * 60.0
        return rt
    return np.nan


def _mzml_binary_arrays(spectrum_elem) -> Tuple[np.ndarray, np.ndarray]:
    mz_array = np.array([], dtype=float)
    intensity_array = np.array([], dtype=float)

    for bda in spectrum_elem.iter():
        if _local_name(bda.tag) != "binaryDataArray":
            continue

        is_mz = _has_cv(bda, accession="MS:1000514", name="m/z array")
        is_intensity = _has_cv(bda, accession="MS:1000515", name="intensity array")
        if not (is_mz or is_intensity):
            continue

        if _has_cv(bda, accession="MS:1000521", name="32-bit float"):
            dtype = np.float32
        elif _has_cv(bda, accession="MS:1000523", name="64-bit float"):
            dtype = np.float64
        else:
            # mzML profile/centroid spectra should normally use 32/64-bit floats.
            dtype = np.float64

        compressed = _has_cv(bda, accession="MS:1000574", name="zlib compression")
        binary_text = None
        for child in bda.iter():
            if _local_name(child.tag) == "binary":
                binary_text = child.text
                break

        arr = _decode_base64_array(binary_text, dtype=dtype, compressed=compressed, byte_order="little")
        if is_mz:
            mz_array = arr
        elif is_intensity:
            intensity_array = arr

    return mz_array, intensity_array


def _iter_mzml_spectra(ms_file: str | Path):
    for event, elem in ET.iterparse(str(ms_file), events=("end",)):
        if _local_name(elem.tag) != "spectrum":
            continue
        ms_level = _mzml_ms_level(elem)
        if ms_level in {1, 2}:
            rt = _mzml_rt_seconds(elem)
            mz_array, intensity_array = _mzml_binary_arrays(elem)
            yield ms_level, rt, mz_array, intensity_array
        elem.clear()


def _mzxml_scan_arrays(scan_elem) -> Tuple[np.ndarray, np.ndarray]:
    peaks_elem = None
    for child in scan_elem:
        if _local_name(child.tag) == "peaks":
            peaks_elem = child
            break
    if peaks_elem is None:
        return np.array([], dtype=float), np.array([], dtype=float)

    precision = str(peaks_elem.attrib.get("precision", "32"))
    dtype = np.float64 if precision == "64" else np.float32
    byte_order = peaks_elem.attrib.get("byteOrder", "network")
    compression_type = str(peaks_elem.attrib.get("compressionType", "none")).lower()
    compressed = "zlib" in compression_type

    arr = _decode_base64_array(peaks_elem.text, dtype=dtype, compressed=compressed, byte_order=byte_order)
    if arr.size < 2:
        return np.array([], dtype=float), np.array([], dtype=float)
    if arr.size % 2 != 0:
        arr = arr[:-1]
    pairs = arr.reshape((-1, 2))
    return pairs[:, 0], pairs[:, 1]


def _iter_mzxml_spectra(ms_file: str | Path):
    for event, elem in ET.iterparse(str(ms_file), events=("end",)):
        if _local_name(elem.tag) != "scan":
            continue
        try:
            ms_level = int(float(elem.attrib.get("msLevel", "0")))
        except Exception:
            ms_level = None
        if ms_level in {1, 2}:
            rt = _parse_iso_duration_to_seconds(elem.attrib.get("retentionTime"))
            mz_array, intensity_array = _mzxml_scan_arrays(elem)
            yield ms_level, rt, mz_array, intensity_array
        elem.clear()


def _iter_spectra(ms_file: str | Path):
    ext = Path(ms_file).suffix.lower()
    if ext == ".mzml":
        yield from _iter_mzml_spectra(ms_file)
    elif ext == ".mzxml":
        yield from _iter_mzxml_spectra(ms_file)
    else:
        raise ValueError("Unsupported file type. Please select a .mzML or .mzXML file.")


def _sum_intensity_in_window(mz_array, intensity_array, target_mz: float, tol: float) -> float:
    if mz_array is None or intensity_array is None:
        return 0.0
    mzs = np.asarray(mz_array, dtype=float)
    ints = np.asarray(intensity_array, dtype=float)
    if mzs.size == 0 or ints.size == 0:
        return 0.0
    if mzs.size != ints.size:
        n = min(mzs.size, ints.size)
        mzs = mzs[:n]
        ints = ints[:n]
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


def _max_intensity_and_apex(eic: pd.DataFrame) -> Tuple[float, float]:
    if eic is None or len(eic) == 0:
        return np.nan, np.nan
    rt = pd.to_numeric(eic["RT"], errors="coerce").to_numpy(dtype=float)
    ints = pd.to_numeric(eic["Intensity"], errors="coerce").to_numpy(dtype=float)
    valid = np.isfinite(rt) & np.isfinite(ints)
    if not np.any(valid):
        return np.nan, np.nan
    rt_v = rt[valid]
    int_v = ints[valid]
    if int_v.size == 0 or np.nanmax(int_v) <= 0:
        return 0.0, np.nan
    idx = int(np.nanargmax(int_v))
    return float(int_v[idx]), float(rt_v[idx])


def _ms2_near_apex_max(ms2_eic: pd.DataFrame, apex_rt: float, window_sec: float) -> float:
    if ms2_eic is None or len(ms2_eic) == 0 or not np.isfinite(apex_rt):
        return 0.0
    rt = pd.to_numeric(ms2_eic["RT"], errors="coerce").to_numpy(dtype=float)
    ints = pd.to_numeric(ms2_eic["Intensity"], errors="coerce").to_numpy(dtype=float)
    mask = np.isfinite(rt) & np.isfinite(ints) & (np.abs(rt - apex_rt) <= float(window_sec))
    if not np.any(mask):
        return 0.0
    return float(np.nanmax(ints[mask]))


def intensity_summary_from_state(state: ModeState, ms1_intensity_threshold: float = 5000.0, uncapped_percent_threshold: float = 1.0) -> pd.DataFrame:
    """Build a species-level intensity maxima table for a mode state.

    The uncapped screen compares each uncapped species MS1 maximum intensity against
    a percentage of the cumulative capped MS1 maximum intensity. This is a diagnostic
    rule, not the capping-efficiency calculation.
    """
    ct = state.species_table
    if ct is None or len(ct) == 0:
        return pd.DataFrame()
    rows = []
    for _, row in ct.reset_index(drop=True).iterrows():
        name = str(row["name"])
        eic = state.eic_list.get(name)
        max_i, apex_rt = _max_intensity_and_apex(eic)
        rows.append({
            "Species": name,
            "mz": pd.to_numeric(row.get("mz"), errors="coerce"),
            "type": row.get("type", ""),
            "charge": row.get("charge", np.nan),
            "use_for_efficiency": row.get("use_for_efficiency", True),
            "MS1_max_intensity": max_i,
            "MS1_apex_RT_sec": apex_rt,
            "MS1_apex_RT_min": (apex_rt / 60.0) if np.isfinite(apex_rt) else np.nan,
            "RT_at_MS1_max_sec": apex_rt,
            "MS1_threshold": float(ms1_intensity_threshold),
            "MS1_confirmed": bool(np.isfinite(max_i) and max_i >= float(ms1_intensity_threshold)),
        })
    out = pd.DataFrame(rows)
    capped_sum = float(pd.to_numeric(out.loc[out["type"] == "Capped", "MS1_max_intensity"], errors="coerce").fillna(0).sum())
    uncapped_limit = capped_sum * float(uncapped_percent_threshold) / 100.0
    out["cumulative_capped_MS1_max_intensity"] = capped_sum
    out["uncapped_percent_of_capped_sum"] = np.nan
    out["uncapped_less_than_threshold"] = pd.Series([pd.NA] * len(out), dtype="object")
    if capped_sum > 0:
        uncapped_mask = out["type"] == "Uncapped"
        out.loc[uncapped_mask, "uncapped_percent_of_capped_sum"] = (
            pd.to_numeric(out.loc[uncapped_mask, "MS1_max_intensity"], errors="coerce").fillna(0) / capped_sum * 100.0
        )
        out.loc[uncapped_mask, "uncapped_less_than_threshold"] = (
            pd.to_numeric(out.loc[uncapped_mask, "MS1_max_intensity"], errors="coerce").fillna(0) < uncapped_limit
        )
    out["uncapped_threshold_percent"] = float(uncapped_percent_threshold)
    out["uncapped_threshold_intensity"] = uncapped_limit
    return out


def make_ms1_message(summary: pd.DataFrame, threshold: float) -> str:
    if summary is None or len(summary) == 0:
        return "No MS1 intensity summary available."
    confirmed = summary.loc[summary["MS1_confirmed"] == True, "Species"].astype(str).tolist()
    return f"MS1-confirmed above {float(threshold):.4g}: {len(confirmed)} of {len(summary)} species: {', '.join(confirmed) if confirmed else 'none'}"


def make_uncapped_intensity_message(summary: pd.DataFrame, percent_threshold: float = 1.0) -> str:
    if summary is None or len(summary) == 0:
        return "No uncapped intensity screen available."
    capped_sum = float(pd.to_numeric(summary.get("cumulative_capped_MS1_max_intensity", pd.Series([0])), errors="coerce").fillna(0).iloc[0])
    if capped_sum <= 0:
        return "Uncapped intensity screen could not be evaluated because cumulative capped MS1 maximum intensity is zero."
    uncapped = summary.loc[summary["type"] == "Uncapped"].copy()
    if len(uncapped) == 0:
        return "Uncapped intensity screen: no uncapped species present."
    below = uncapped.loc[uncapped["uncapped_less_than_threshold"] == True, "Species"].astype(str).tolist()
    above = uncapped.loc[uncapped["uncapped_less_than_threshold"] == False, "Species"].astype(str).tolist()
    limit = capped_sum * float(percent_threshold) / 100.0
    return (
        f"Uncapped intensity screen: threshold = {float(percent_threshold):.3g}% of cumulative capped MS1 maxima "
        f"({capped_sum:.4g}) = {limit:.4g}. "
        f"Below threshold: {', '.join(below) if below else 'none'}. "
        f"At/above threshold: {', '.join(above) if above else 'none'}."
    )





def _row_neutral_mass(row: pd.Series) -> float:
    """Return neutral mass using neutral_mass where available, otherwise m/z and charge."""
    nm = pd.to_numeric(pd.Series([row.get("neutral_mass", np.nan)]), errors="coerce").iloc[0]
    if np.isfinite(nm) and nm > 0:
        return float(nm)
    mz = pd.to_numeric(pd.Series([row.get("mz", np.nan)]), errors="coerce").iloc[0]
    z = pd.to_numeric(pd.Series([row.get("charge", np.nan)]), errors="coerce").iloc[0]
    if not np.isfinite(z) or z <= 0:
        z = 1.0
    if not np.isfinite(mz) or mz <= 0:
        return np.nan
    return neutral_from_mz(float(mz), float(z))


def annotate_labile_phosphate_losses(summary: pd.DataFrame, species_tbl: pd.DataFrame, rt_tolerance_sec: float = 5.0, enabled: bool = True) -> pd.DataFrame:
    """Flag co-eluted lower-phosphate uncapped species as likely labile phosphate-loss ions.

    A lower-phosphate species is flagged when there is another uncapped species in the
    same charge state whose neutral mass is approximately +79.966331 Da and whose MS1
    apex RT is within ``rt_tolerance_sec``. This targets cases such as pppG and ppG
    co-eluting, where the lower phosphate state can represent in-source/labile phosphate
    loss rather than an independent uncapped molecule.
    """
    if summary is None or len(summary) == 0:
        return pd.DataFrame() if summary is None else summary
    out = summary.copy()
    out["labile_parent_species"] = ""
    out["labile_parent_mz"] = np.nan
    out["labile_parent_RT_sec"] = np.nan
    out["labile_RT_delta_sec"] = np.nan
    out["labile_mass_delta_Da"] = np.nan
    out["ignored_as_labile_phosphate_loss"] = False
    out["ignore_for_efficiency_reason"] = ""
    if not enabled or species_tbl is None or len(species_tbl) == 0:
        return out

    sdt = species_tbl.reset_index(drop=True).copy()
    sdt["_neutral_calc"] = [_row_neutral_mass(row) for _, row in sdt.iterrows()]
    sdt["_charge_calc"] = pd.to_numeric(sdt.get("charge", pd.Series([np.nan] * len(sdt))), errors="coerce")
    sdt.loc[~np.isfinite(sdt["_charge_calc"]) | (sdt["_charge_calc"] <= 0), "_charge_calc"] = 1.0

    # Use a conservative neutral-mass tolerance. Uploaded/generated masses should be exact;
    # this tolerance only protects against rounding in displayed species tables.
    mass_tol = 0.10
    rt_tol = float(rt_tolerance_sec) if np.isfinite(float(rt_tolerance_sec)) else 5.0

    for i, lower in sdt.iterrows():
        if str(lower.get("type", "")) != "Uncapped":
            continue
        lower_name = str(lower.get("name", ""))
        lower_summary = out.loc[out["Species"].astype(str) == lower_name]
        if len(lower_summary) != 1:
            continue
        lower_rt = pd.to_numeric(lower_summary.iloc[0].get("MS1_apex_RT_sec", np.nan), errors="coerce")
        lower_nm = float(lower.get("_neutral_calc", np.nan))
        lower_z = float(lower.get("_charge_calc", 1.0))
        if not np.isfinite(lower_nm) or not np.isfinite(lower_rt):
            continue

        best = None
        for j, parent in sdt.iterrows():
            if i == j or str(parent.get("type", "")) != "Uncapped":
                continue
            parent_z = float(parent.get("_charge_calc", 1.0))
            if abs(parent_z - lower_z) > 0.01:
                continue
            parent_nm = float(parent.get("_neutral_calc", np.nan))
            if not np.isfinite(parent_nm):
                continue
            delta = parent_nm - lower_nm
            if abs(delta - PHOSPHATE_STEP_MASS) > mass_tol:
                continue
            parent_name = str(parent.get("name", ""))
            parent_summary = out.loc[out["Species"].astype(str) == parent_name]
            if len(parent_summary) != 1:
                continue
            parent_rt = pd.to_numeric(parent_summary.iloc[0].get("MS1_apex_RT_sec", np.nan), errors="coerce")
            if not np.isfinite(parent_rt):
                continue
            rt_delta = abs(float(parent_rt) - float(lower_rt))
            if rt_delta <= rt_tol:
                score = (rt_delta, abs(delta - PHOSPHATE_STEP_MASS))
                if best is None or score < best[0]:
                    best = (score, parent_name, float(parent.get("mz", np.nan)), float(parent_rt), rt_delta, delta)
        if best is not None:
            _score, parent_name, parent_mz, parent_rt, rt_delta, delta = best
            idx = out.index[out["Species"].astype(str) == lower_name]
            out.loc[idx, "labile_parent_species"] = parent_name
            out.loc[idx, "labile_parent_mz"] = parent_mz
            out.loc[idx, "labile_parent_RT_sec"] = parent_rt
            out.loc[idx, "labile_RT_delta_sec"] = rt_delta
            out.loc[idx, "labile_mass_delta_Da"] = delta
            out.loc[idx, "ignored_as_labile_phosphate_loss"] = True
            out.loc[idx, "ignore_for_efficiency_reason"] = "co-eluted labile phosphate loss"
    return out


def make_labile_loss_message(summary: pd.DataFrame) -> str:
    if summary is None or len(summary) == 0 or "ignored_as_labile_phosphate_loss" not in summary.columns:
        return "Labile phosphate-loss screen unavailable."
    flagged = summary.loc[summary["ignored_as_labile_phosphate_loss"] == True]
    if len(flagged) == 0:
        return "Labile phosphate-loss screen: no co-eluted lower-phosphate species flagged."
    parts = []
    for _, row in flagged.iterrows():
        parts.append(f"{row['Species']} ignored as labile loss of {row.get('labile_parent_species', '')}")
    return "Labile phosphate-loss screen: " + "; ".join(parts) + "."


def build_analysis(
    ms_file: str | Path,
    species_tbl: pd.DataFrame,
    default_mz_tol: float = 0.01,
    default_ms2_window: float = 30.0,
    ms1_intensity_threshold: float = 5000.0,
    uncapped_percent_threshold: float = 1.0,
    ignore_labile_phosphate_loss: bool = True,
    labile_rt_tolerance_sec: float = 5.0,
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

    # Stream spectra once and build EICs for all user-defined species.
    for scan_idx, (ms_level, rt, mz_array, intensity_array) in enumerate(_iter_spectra(ms_file), start=1):
        if ms_level not in {1, 2}:
            continue
        if not np.isfinite(rt):
            continue

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

    if len(ms1_rts) == 0:
        raise ValueError("No MS1 spectra were found in the selected file. Check that the file is a valid centroid/profile mzML or mzXML export containing MS1 scans.")

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
    all_state.intensity_summary = intensity_summary_from_state(
        all_state,
        ms1_intensity_threshold=ms1_intensity_threshold,
        uncapped_percent_threshold=uncapped_percent_threshold,
    )
    all_state.intensity_summary = annotate_labile_phosphate_losses(
        all_state.intensity_summary,
        all_state.species_table,
        rt_tolerance_sec=labile_rt_tolerance_sec,
        enabled=ignore_labile_phosphate_loss,
    )

    # Add MS2-near-apex maxima to the intensity summary.
    ms2_maxima = []
    for i in range(n_species):
        name = str(sdt.iloc[i]["name"])
        apex_rt = float(all_state.intensity_summary.iloc[i]["MS1_apex_RT_sec"]) if len(all_state.intensity_summary) > i else np.nan
        window_i = _get_ms2_window_for_row(sdt, i, default_ms2_window)
        ms2_maxima.append(_ms2_near_apex_max(ms2_eic_list.get(name), apex_rt, window_i))
    all_state.intensity_summary["MS2_near_apex_max_intensity"] = ms2_maxima

    if len(confirmed_idx) == 0:
        filtered_table = pd.DataFrame(columns=SPECIES_COLUMNS)
        filtered_state = ModeState(
            species_table=filtered_table,
            eic_df=empty_eic_df(),
            eic_list={},
            rt_windows=blank_windows(filtered_table),
            aucs=pd.Series(dtype="float64"),
            intensity_summary=pd.DataFrame(),
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
        filtered_state.intensity_summary = all_state.intensity_summary[
            all_state.intensity_summary["Species"].astype(str).isin(filtered_table["name"].astype(str))
        ].reset_index(drop=True)
        names = ", ".join(filtered_table["name"].astype(str).tolist())
        filter_message = f"MS2-confirmed: {len(confirmed_idx)} of {n_species} species: {names}"

    ms1_message = make_ms1_message(all_state.intensity_summary, ms1_intensity_threshold)
    uncapped_message = make_uncapped_intensity_message(all_state.intensity_summary, uncapped_percent_threshold)
    labile_message = make_labile_loss_message(all_state.intensity_summary)
    session = AnalysisSession(
        ms_file=ms_file,
        all=all_state,
        filtered=filtered_state,
        filter_message=filter_message,
        ms1_message=ms1_message,
        uncapped_intensity_message=uncapped_message,
        labile_loss_message=labile_message,
        loaded=True,
    )
    return session, msg


def _contiguous_runs(indices: np.ndarray) -> List[np.ndarray]:
    if indices is None or len(indices) == 0:
        return []
    idx = np.asarray(indices, dtype=int)
    splits = np.where(np.diff(idx) > 1)[0] + 1
    return [run for run in np.split(idx, splits) if len(run) > 0]


def integrate_eic_window_with_qc(
    eic: pd.DataFrame,
    rt_min: float,
    rt_max: float,
    min_intensity: float = 0.0,
    min_consecutive_points: int = 2,
) -> Dict[str, float | int | bool]:
    """Integrate selected EIC points while enforcing a minimum contiguous peak width.

    A point is eligible when it falls inside the selected RT window and its intensity is
    greater than ``min_intensity``. Eligible points are split into contiguous runs in scan
    order. Only runs with at least ``min_consecutive_points`` are integrated. This avoids
    treating isolated 2-3 point spikes as valid peaks.
    """
    result = {
        "AUC": np.nan,
        "AUC_points_used": 0,
        "AUC_segments_used": 0,
        "AUC_max_consecutive_points": 0,
        "AUC_peak_passed_min_points": False,
    }
    if eic is None or len(eic) == 0:
        return result
    lo, hi = sorted([float(rt_min), float(rt_max)])
    rt = pd.to_numeric(eic["RT"], errors="coerce").to_numpy(dtype=float)
    ints = pd.to_numeric(eic["Intensity"], errors="coerce").to_numpy(dtype=float)
    min_i = float(min_intensity) if np.isfinite(float(min_intensity)) else 0.0
    min_pts = max(2, int(min_consecutive_points))
    eligible = np.where(np.isfinite(rt) & np.isfinite(ints) & (rt >= lo) & (rt <= hi) & (ints > min_i))[0]
    runs = _contiguous_runs(eligible)
    result["AUC_max_consecutive_points"] = int(max((len(run) for run in runs), default=0))
    valid_runs = [run for run in runs if len(run) >= min_pts]
    result["AUC_segments_used"] = int(len(valid_runs))
    result["AUC_points_used"] = int(sum(len(run) for run in valid_runs))
    result["AUC_peak_passed_min_points"] = bool(len(valid_runs) > 0)
    if not valid_runs:
        return result
    auc = 0.0
    for run in valid_runs:
        if len(run) >= 2:
            if hasattr(np, "trapezoid"):
                auc += float(np.trapezoid(ints[run], rt[run]))
            else:
                auc += float(np.trapz(ints[run], rt[run]))
    result["AUC"] = float(auc) if auc > 0 else np.nan
    return result


def integrate_eic_window(eic: pd.DataFrame, rt_min: float, rt_max: float, min_intensity: float = 0.0, min_consecutive_points: int = 2) -> float:
    return float(integrate_eic_window_with_qc(eic, rt_min, rt_max, min_intensity, min_consecutive_points)["AUC"])


def save_rt_window(
    session: AnalysisSession,
    requested_mode: str,
    species_name: str,
    rt_min: float,
    rt_max: float,
    min_intensity: float = 0.0,
    min_consecutive_points: int = 2,
) -> float:
    state = session.state(requested_mode)
    if state.species_table is None or len(state.species_table) == 0:
        raise ValueError("No species available in this mode.")
    matches = state.species_table.index[state.species_table["name"].astype(str) == str(species_name)].tolist()
    if len(matches) != 1:
        raise ValueError("Selected species was not found in the current table.")
    idx = matches[0]
    eic = state.eic_list.get(str(species_name))
    qc = integrate_eic_window_with_qc(eic, rt_min, rt_max, min_intensity=min_intensity, min_consecutive_points=min_consecutive_points)
    auc = float(qc["AUC"])
    lo = min(float(rt_min), float(rt_max))
    hi = max(float(rt_min), float(rt_max))
    for col in ["RTmin", "RTmax", "RTmin_min", "RTmax_min", "AUC_min_intensity", "AUC_min_consecutive_points", "AUC_points_used", "AUC_segments_used", "AUC_max_consecutive_points", "AUC_peak_passed_min_points"]:
        if col not in state.rt_windows.columns:
            state.rt_windows[col] = pd.NA if col == "AUC_peak_passed_min_points" else np.nan
    state.rt_windows.loc[idx, "RTmin"] = lo
    state.rt_windows.loc[idx, "RTmax"] = hi
    state.rt_windows.loc[idx, "RTmin_min"] = lo / 60.0
    state.rt_windows.loc[idx, "RTmax_min"] = hi / 60.0
    state.rt_windows.loc[idx, "AUC_min_intensity"] = float(min_intensity)
    state.rt_windows.loc[idx, "AUC_min_consecutive_points"] = int(min_consecutive_points)
    state.rt_windows.loc[idx, "AUC_points_used"] = int(qc["AUC_points_used"])
    state.rt_windows.loc[idx, "AUC_segments_used"] = int(qc["AUC_segments_used"])
    state.rt_windows.loc[idx, "AUC_max_consecutive_points"] = int(qc["AUC_max_consecutive_points"])
    state.rt_windows.loc[idx, "AUC_peak_passed_min_points"] = bool(qc["AUC_peak_passed_min_points"])
    state.aucs.iloc[idx] = auc
    # Do not raise when the selected region fails the peak criteria. The failed
    # species remains visible in the RT/AUC table as "not found beyond the
    # selected limits" and can be ignored automatically in the efficiency
    # calculation. This avoids one missing/weak species blocking the whole assay.
    return auc


def intensity_table(session: AnalysisSession, requested_mode: str) -> pd.DataFrame:
    state = session.state(requested_mode)
    if state.intensity_summary is None:
        return pd.DataFrame()
    return state.intensity_summary.copy()


def _species_detection_status(row: pd.Series) -> Tuple[bool, str]:
    """Classify whether a species has been found beyond the user-set limits.

    A species is treated as found for efficiency only when it has passed MS1
    confirmation, has a saved RT window, passes the contiguous-point AUC rule,
    and has a finite positive AUC. Species that fail these checks are kept in
    exported tables but can be automatically ignored in the capping-efficiency
    calculation.
    """
    ms1_confirmed = str(row.get("MS1_confirmed", "")).strip().lower() in {"true", "1", "yes"}
    rtmin = pd.to_numeric(pd.Series([row.get("RTmin", np.nan)]), errors="coerce").iloc[0]
    rtmax = pd.to_numeric(pd.Series([row.get("RTmax", np.nan)]), errors="coerce").iloc[0]
    auc = pd.to_numeric(pd.Series([row.get("AUC", np.nan)]), errors="coerce").iloc[0]
    peak_passed_raw = row.get("AUC_peak_passed_min_points", pd.NA)
    peak_passed = str(peak_passed_raw).strip().lower() in {"true", "1", "yes"}

    if not ms1_confirmed:
        return False, "species not found beyond MS1 intensity threshold"
    if not np.isfinite(rtmin) or not np.isfinite(rtmax):
        return False, "RT window not saved"
    if not peak_passed:
        return False, "species not found beyond AUC intensity/contiguous-point limit"
    if not np.isfinite(auc) or auc <= 0:
        return False, "species not found beyond AUC limit"
    return True, "found beyond limits"


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
            "RTmin_min": rt_tbl["RTmin_min"].values if "RTmin_min" in rt_tbl.columns else (pd.to_numeric(rt_tbl["RTmin"], errors="coerce") / 60.0).values if "RTmin" in rt_tbl.columns else np.nan,
            "RTmax_min": rt_tbl["RTmax_min"].values if "RTmax_min" in rt_tbl.columns else (pd.to_numeric(rt_tbl["RTmax"], errors="coerce") / 60.0).values if "RTmax" in rt_tbl.columns else np.nan,
            "AUC_min_intensity": rt_tbl["AUC_min_intensity"].values if "AUC_min_intensity" in rt_tbl.columns else np.nan,
            "AUC_min_consecutive_points": rt_tbl["AUC_min_consecutive_points"].values if "AUC_min_consecutive_points" in rt_tbl.columns else np.nan,
            "AUC_points_used": rt_tbl["AUC_points_used"].values if "AUC_points_used" in rt_tbl.columns else np.nan,
            "AUC_segments_used": rt_tbl["AUC_segments_used"].values if "AUC_segments_used" in rt_tbl.columns else np.nan,
            "AUC_max_consecutive_points": rt_tbl["AUC_max_consecutive_points"].values if "AUC_max_consecutive_points" in rt_tbl.columns else np.nan,
            "AUC_peak_passed_min_points": rt_tbl["AUC_peak_passed_min_points"].values if "AUC_peak_passed_min_points" in rt_tbl.columns else pd.NA,
            "AUC": state.aucs.values,
            "notes": ct["notes"].values if "notes" in ct.columns else "",
        }
    )
    intensity = intensity_table(session, requested_mode)
    if intensity is not None and len(intensity) > 0 and "Species" in intensity.columns:
        keep_cols = [
            "Species",
            "MS1_max_intensity",
            "MS1_apex_RT_sec",
            "MS1_apex_RT_min",
            "RT_at_MS1_max_sec",
            "MS1_confirmed",
            "MS1_threshold",
            "MS2_near_apex_max_intensity",
            "uncapped_percent_of_capped_sum",
            "uncapped_less_than_threshold",
            "uncapped_threshold_percent",
            "uncapped_threshold_intensity",
            "labile_parent_species",
            "labile_parent_mz",
            "labile_parent_RT_sec",
            "labile_RT_delta_sec",
            "labile_mass_delta_Da",
            "ignored_as_labile_phosphate_loss",
            "ignore_for_efficiency_reason",
        ]
        keep_cols = [c for c in keep_cols if c in intensity.columns]
        out = out.merge(intensity[keep_cols], on="Species", how="left")
    include_base = as_logical_safely(out["use_for_efficiency"], default=True).values if "use_for_efficiency" in out.columns else np.array([True] * len(out))
    labile_ignored = out["ignored_as_labile_phosphate_loss"].fillna(False).astype(bool).values if "ignored_as_labile_phosphate_loss" in out.columns else np.array([False] * len(out))
    out["use_for_efficiency_effective"] = include_base & (~labile_ignored)

    detection = out.apply(_species_detection_status, axis=1)
    out["species_found_beyond_limit"] = [bool(x[0]) for x in detection]
    out["species_detection_status"] = [str(x[1]) for x in detection]
    out["included_for_efficiency_final"] = out["use_for_efficiency_effective"].fillna(False).astype(bool) & out["species_found_beyond_limit"].fillna(False).astype(bool)
    if "ignored_as_labile_phosphate_loss" in out.columns:
        labile_mask = out["ignored_as_labile_phosphate_loss"].fillna(False).astype(bool)
        out.loc[labile_mask, "species_detection_status"] = (
            out.loc[labile_mask, "species_detection_status"].astype(str)
            + "; ignored as co-eluted labile phosphate-loss species"
        )
    return out



# -----------------------------------------------------------------------------
# Composite EIC helpers
# -----------------------------------------------------------------------------

def _safe_float(value) -> float:
    try:
        v = float(value)
        return v if np.isfinite(v) else np.nan
    except Exception:
        return np.nan


def infer_species_charge(row) -> float:
    """Infer charge state from the explicit charge column or common name patterns."""
    charge = _safe_float(row.get("charge", np.nan))
    if np.isfinite(charge) and charge > 0:
        return charge
    name = str(row.get("Species", row.get("name", "")))
    # Match z1, z=1, _z1, -z2, charge2, etc.
    m = re.search(r"(?:^|[^A-Za-z0-9])z\s*=?\s*(\d+)(?:$|[^A-Za-z0-9])", name, flags=re.I)
    if not m:
        m = re.search(r"(?:charge|chg)\s*=?\s*(\d+)", name, flags=re.I)
    if m:
        try:
            z = float(m.group(1))
            return z if z > 0 else np.nan
        except Exception:
            return np.nan
    return np.nan


def candidate_neutral_masses_for_grouping(row, allowed_charges=(1, 2)) -> list[tuple[int, float]]:
    """Return plausible neutral masses for capped charge-state grouping.

    If charge and neutral_mass are supplied, they are used directly. If charge is
    missing, z=1 and z=2 candidates are generated from m/z so that species named
    only by their m/z values can still be paired.
    """
    mz = _safe_float(row.get("mz", np.nan))
    explicit_neutral = _safe_float(row.get("neutral_mass", np.nan))
    charge = infer_species_charge(row)
    out: list[tuple[int, float]] = []
    if np.isfinite(explicit_neutral) and explicit_neutral > 0:
        if np.isfinite(charge) and charge > 0:
            out.append((int(round(charge)), float(explicit_neutral)))
        else:
            for z in allowed_charges:
                out.append((int(z), float(explicit_neutral)))
        return out
    if not np.isfinite(mz) or mz <= 0:
        return out
    if np.isfinite(charge) and charge > 0:
        z = int(round(charge))
        return [(z, neutral_from_mz(mz, z))]
    for z in allowed_charges:
        out.append((int(z), neutral_from_mz(mz, z)))
    return out


def _rows_share_neutral_mass(row_a, row_b, tolerance_da: float = 0.05) -> tuple[bool, float]:
    cand_a = candidate_neutral_masses_for_grouping(row_a)
    cand_b = candidate_neutral_masses_for_grouping(row_b)
    best_delta = np.inf
    best_mass = np.nan
    for za, ma in cand_a:
        for zb, mb in cand_b:
            delta = abs(float(ma) - float(mb))
            if delta < best_delta:
                best_delta = delta
                best_mass = (float(ma) + float(mb)) / 2.0
            # Prefer grouping different charge states, but allow explicit same-neutral grouping.
            if delta <= float(tolerance_da) and (za != zb or np.isfinite(_safe_float(row_a.get("neutral_mass", np.nan)))):
                return True, (float(ma) + float(mb)) / 2.0
    return False, best_mass


def _combine_eic_frames(eics: list[pd.DataFrame]) -> pd.DataFrame:
    """Sum EIC intensities across traces, interpolating to a common RT grid if needed."""
    valid = []
    for eic in eics:
        if eic is None or len(eic) == 0 or "RT" not in eic.columns or "Intensity" not in eic.columns:
            continue
        rt = pd.to_numeric(eic["RT"], errors="coerce").to_numpy(dtype=float)
        intensity = pd.to_numeric(eic["Intensity"], errors="coerce").to_numpy(dtype=float)
        mask = np.isfinite(rt) & np.isfinite(intensity)
        if np.any(mask):
            order = np.argsort(rt[mask])
            valid.append((rt[mask][order], intensity[mask][order]))
    if not valid:
        return empty_eic_df()
    if len(valid) == 1:
        return pd.DataFrame({"RT": valid[0][0], "Intensity": valid[0][1], "Species": "combined"})
    same_grid = all(len(rt) == len(valid[0][0]) and np.allclose(rt, valid[0][0], rtol=0, atol=1e-9) for rt, _ in valid[1:])
    if same_grid:
        grid = valid[0][0]
        summed = np.zeros_like(grid, dtype=float)
        for _rt, intensity in valid:
            summed += intensity
    else:
        grid = np.unique(np.concatenate([rt for rt, _intensity in valid]))
        summed = np.zeros_like(grid, dtype=float)
        for rt, intensity in valid:
            # Outside each EIC's acquired RT span, contribution is zero.
            summed += np.interp(grid, rt, intensity, left=0.0, right=0.0)
    return pd.DataFrame({"RT": grid, "Intensity": summed, "Species": "combined"})



def combine_capped_charge_state_auc_table(
    auc_tbl: pd.DataFrame,
    neutral_mass_tolerance_da: float = 0.05,
) -> pd.DataFrame:
    """Return an AUC table where capped z=1/z=2 rows are combined.

    This is the analysis-level counterpart to the composite EIC charge-state
    summing. Uncapped rows are left unchanged. Capped rows that resolve to the
    same neutral mass are collapsed into one row and their AUC values are summed.

    The transform is intentionally conservative: it combines rows only when the
    explicit neutral_mass/charge metadata or inferred m/z/charge candidates match
    within ``neutral_mass_tolerance_da``. It does not combine unrelated capped
    species such as different capped fragments.
    """
    if auc_tbl is None or len(auc_tbl) == 0:
        return pd.DataFrame() if auc_tbl is None else auc_tbl.copy()

    tbl = auc_tbl.copy().reset_index(drop=True)
    if "type" not in tbl.columns or "Species" not in tbl.columns:
        return tbl

    non_capped = tbl.loc[tbl["type"].astype(str) != "Capped"].copy()
    capped = tbl.loc[tbl["type"].astype(str) == "Capped"].copy().reset_index(drop=True)
    if len(capped) <= 1:
        tbl["combined_capped_charge_states_for_efficiency"] = False
        tbl["source_species_for_combined_capped"] = ""
        return tbl

    rows_out: list[dict] = []
    assigned: set[int] = set()

    for i, row_i in capped.iterrows():
        if i in assigned:
            continue
        group_indices = [i]
        assigned.add(i)
        group_mass_candidates = []

        for j, row_j in capped.iterrows():
            if j in assigned:
                continue
            shares, mass = _rows_share_neutral_mass(row_i, row_j, tolerance_da=neutral_mass_tolerance_da)
            if shares:
                group_indices.append(j)
                assigned.add(j)
                if np.isfinite(mass):
                    group_mass_candidates.append(float(mass))

        group = capped.iloc[group_indices].copy()
        source_species = group["Species"].astype(str).tolist()
        out_row = group.iloc[0].to_dict()
        out_row["combined_capped_charge_states_for_efficiency"] = bool(len(group) > 1)
        out_row["source_species_for_combined_capped"] = ", ".join(source_species) if len(group) > 1 else ""

        inferred_masses = []
        charges = []
        for _, grow in group.iterrows():
            for _z, mass in candidate_neutral_masses_for_grouping(grow):
                if np.isfinite(mass):
                    inferred_masses.append(float(mass))
            ch = infer_species_charge(grow)
            if np.isfinite(ch):
                charges.append(str(int(round(ch))))

        med_mass = float(np.median(inferred_masses)) if inferred_masses else (float(np.median(group_mass_candidates)) if group_mass_candidates else np.nan)
        charge_label = "+".join(sorted(set(charges), key=lambda x: int(x) if str(x).isdigit() else 99)) if charges else "combined"
        if len(group) > 1:
            out_row["Species"] = f"Combined capped M~{med_mass:.2f} z{charge_label}" if np.isfinite(med_mass) else "Combined capped charge states"
            out_row["charge"] = charge_label
            out_row["neutral_mass"] = med_mass

        # Numerical fields where a combined capped species should represent the
        # summed charge-state contribution.
        for col in ["AUC", "AUC_points_used", "AUC_segments_used", "MS1_max_intensity", "MS2_near_apex_max_intensity"]:
            if col in group.columns:
                vals = pd.to_numeric(group[col], errors="coerce")
                out_row[col] = float(vals.sum()) if vals.notna().any() else np.nan

        for col in ["AUC_max_consecutive_points"]:
            if col in group.columns:
                vals = pd.to_numeric(group[col], errors="coerce")
                out_row[col] = float(vals.max()) if vals.notna().any() else np.nan

        for col in ["RTmin", "MS1_apex_RT_sec", "RT_at_MS1_max_sec"]:
            if col in group.columns:
                vals = pd.to_numeric(group[col], errors="coerce")
                out_row[col] = float(vals.min()) if vals.notna().any() else np.nan
        for col in ["RTmax"]:
            if col in group.columns:
                vals = pd.to_numeric(group[col], errors="coerce")
                out_row[col] = float(vals.max()) if vals.notna().any() else np.nan

        if "RTmin" in out_row:
            out_row["RTmin_min"] = out_row["RTmin"] / 60.0 if np.isfinite(_safe_float(out_row["RTmin"])) else np.nan
        if "RTmax" in out_row:
            out_row["RTmax_min"] = out_row["RTmax"] / 60.0 if np.isfinite(_safe_float(out_row["RTmax"])) else np.nan
        if "MS1_apex_RT_sec" in out_row:
            out_row["MS1_apex_RT_min"] = out_row["MS1_apex_RT_sec"] / 60.0 if np.isfinite(_safe_float(out_row["MS1_apex_RT_sec"])) else np.nan

        # Boolean status columns: combined species passes if at least one charge
        # state row passes. Rows not found beyond limits can still be ignored by
        # calculate_efficiency after this table is built.
        for col in ["AUC_peak_passed_min_points", "MS1_confirmed", "species_found_beyond_limit", "included_for_efficiency_final"]:
            if col in group.columns:
                out_row[col] = bool(group[col].astype(str).str.lower().isin(["true", "1", "yes"]).any())

        if "species_detection_status" in group.columns:
            statuses = [str(x) for x in group["species_detection_status"].dropna().tolist() if str(x).strip()]
            out_row["species_detection_status"] = "combined capped charge states; " + " | ".join(statuses) if len(group) > 1 else (statuses[0] if statuses else "")

        rows_out.append(out_row)

    capped_out = pd.DataFrame(rows_out)
    if len(non_capped) > 0:
        for col in ["combined_capped_charge_states_for_efficiency", "source_species_for_combined_capped"]:
            if col not in non_capped.columns:
                non_capped[col] = False if col.startswith("combined") else ""
        out = pd.concat([capped_out, non_capped], ignore_index=True, sort=False)
    else:
        out = capped_out

    if len(out) > 0 and "type" in out.columns:
        out["_sort"] = out["type"].map({"Capped": 0, "Uncapped": 1}).fillna(2)
        out = out.sort_values(["_sort", "Species"]).drop(columns=["_sort"]).reset_index(drop=True)
    return out

def composite_eic_data(
    session: AnalysisSession,
    requested_mode: str,
    combine_capped_charge_states: bool = False,
    neutral_mass_tolerance_da: float = 0.05,
) -> tuple[dict, pd.DataFrame]:
    """Return EIC map and species table for the final composite plot.

    When ``combine_capped_charge_states`` is True, capped z=1/z=2 traces that
    resolve to the same neutral mass are summed and shown as one capped trace.
    This changes composite-plot visualization only; capping efficiency remains
    based on the AUC table unless the user explicitly changes future logic.
    """
    state = session.state(requested_mode)
    tbl = auc_table(session, requested_mode)
    if tbl is None or len(tbl) == 0:
        return {}, pd.DataFrame()
    include_col = "included_for_efficiency_final" if "included_for_efficiency_final" in tbl.columns else ("use_for_efficiency_effective" if "use_for_efficiency_effective" in tbl.columns else "use_for_efficiency")
    include = tbl[include_col].astype(str).str.lower().isin(["true", "1", "yes"])
    plot_tbl = tbl.loc[include & tbl["type"].astype(str).isin(["Capped", "Uncapped"])].copy().reset_index(drop=True)
    if len(plot_tbl) == 0:
        return {}, plot_tbl
    if not combine_capped_charge_states:
        eic_map = {}
        for _, row in plot_tbl.iterrows():
            species = str(row.get("Species", ""))
            if species:
                eic_map[species] = state.eic_list.get(species, empty_eic_df())
        return eic_map, plot_tbl

    eic_map: dict[str, pd.DataFrame] = {}
    rows_out: list[dict] = []

    uncapped = plot_tbl.loc[plot_tbl["type"].astype(str) == "Uncapped"].copy()
    for _, row in uncapped.iterrows():
        species = str(row.get("Species", ""))
        if not species:
            continue
        eic_map[species] = state.eic_list.get(species, empty_eic_df())
        rows_out.append(row.to_dict())

    capped = plot_tbl.loc[plot_tbl["type"].astype(str) == "Capped"].copy().reset_index(drop=True)
    assigned = set()
    for i, row_i in capped.iterrows():
        if i in assigned:
            continue
        group_indices = [i]
        group_mass_candidates = []
        assigned.add(i)
        for j, row_j in capped.iterrows():
            if j in assigned:
                continue
            shares, mass = _rows_share_neutral_mass(row_i, row_j, tolerance_da=neutral_mass_tolerance_da)
            if shares:
                group_indices.append(j)
                assigned.add(j)
                if np.isfinite(mass):
                    group_mass_candidates.append(float(mass))
        group = capped.iloc[group_indices].copy()
        source_species = group["Species"].astype(str).tolist()
        eics = [state.eic_list.get(sp, empty_eic_df()) for sp in source_species]
        combined_eic = _combine_eic_frames(eics)
        inferred_masses = []
        for _, grow in group.iterrows():
            for _z, mass in candidate_neutral_masses_for_grouping(grow):
                if np.isfinite(mass):
                    inferred_masses.append(float(mass))
        if inferred_masses:
            med_mass = float(np.median(inferred_masses))
        elif group_mass_candidates:
            med_mass = float(np.median(group_mass_candidates))
        else:
            med_mass = np.nan
        charges = []
        for _, grow in group.iterrows():
            ch = infer_species_charge(grow)
            if np.isfinite(ch):
                charges.append(str(int(round(ch))))
        if not charges and len(group) > 1:
            charges = ["1", "2"]
        charge_label = "+".join(sorted(set(charges), key=lambda x: int(x) if x.isdigit() else 99)) if charges else "combined"
        if len(group) > 1:
            label = f"Combined capped M~{med_mass:.2f} z{charge_label}" if np.isfinite(med_mass) else "Combined capped charge states"
        else:
            label = str(group.iloc[0].get("Species", "Capped"))
        combined_eic = combined_eic.copy()
        if len(combined_eic) > 0:
            combined_eic["Species"] = label
        eic_map[label] = combined_eic
        out_row = group.iloc[0].to_dict()
        out_row["Species"] = label
        out_row["type"] = "Capped"
        out_row["combined_capped_charge_states"] = bool(len(group) > 1)
        out_row["source_species"] = ", ".join(source_species)
        if "AUC" in group.columns:
            out_row["AUC"] = pd.to_numeric(group["AUC"], errors="coerce").sum()
        for col in ["RTmin", "RTmax"]:
            vals = pd.to_numeric(group[col], errors="coerce") if col in group.columns else pd.Series(dtype=float)
            out_row[col] = float(vals.min() if col == "RTmin" else vals.max()) if vals.notna().any() else np.nan
        if "RTmin" in out_row:
            out_row["RTmin_min"] = out_row["RTmin"] / 60.0 if np.isfinite(_safe_float(out_row["RTmin"])) else np.nan
        if "RTmax" in out_row:
            out_row["RTmax_min"] = out_row["RTmax"] / 60.0 if np.isfinite(_safe_float(out_row["RTmax"])) else np.nan
        rows_out.append(out_row)

    out_tbl = pd.DataFrame(rows_out)
    # Plot capped combined trace(s) first, then uncapped traces.
    if len(out_tbl) > 0 and "type" in out_tbl.columns:
        out_tbl["_sort"] = out_tbl["type"].map({"Capped": 0, "Uncapped": 1}).fillna(2)
        out_tbl = out_tbl.sort_values(["_sort", "Species"]).drop(columns=["_sort"]).reset_index(drop=True)
    return eic_map, out_tbl

def calculate_efficiency(
    auc_tbl: pd.DataFrame,
    honor_use_for_efficiency: bool = True,
    ignore_not_found: bool = True,
    allow_bounded_efficiency: bool = True,
    lower_bound_percent: float = 99.0,
    upper_bound_percent: float = 1.0,
    combine_capped_charge_states: bool = False,
) -> Dict[str, float | str | bool]:
    if auc_tbl is None or len(auc_tbl) == 0:
        return {"ok": False, "message": "No species available in this mode."}

    df0 = auc_tbl.copy()
    if honor_use_for_efficiency:
        if "use_for_efficiency_effective" in df0.columns:
            include = as_logical_safely(df0["use_for_efficiency_effective"], default=True).values
        elif "use_for_efficiency" in df0.columns:
            include = as_logical_safely(df0["use_for_efficiency"], default=True).values
        else:
            include = np.array([True] * len(df0))
        df0 = df0.loc[include].copy()

    initially_included_count = int(len(df0))
    ignored_not_found_count = 0
    ignored_not_found_species: List[str] = []

    if ignore_not_found:
        if "species_found_beyond_limit" in df0.columns:
            found = as_logical_safely(df0["species_found_beyond_limit"], default=False).values
        else:
            found = pd.to_numeric(df0.get("AUC", pd.Series([np.nan] * len(df0))), errors="coerce").fillna(0).values > 0
        ignored = df0.loc[~found].copy()
        ignored_not_found_count = int(len(ignored))
        if len(ignored) > 0 and "Species" in ignored.columns:
            ignored_not_found_species = ignored["Species"].astype(str).tolist()
        df = df0.loc[found].copy()
    else:
        df = df0.copy()

    if len(df) == 0:
        return {
            "ok": False,
            "message": "Capping efficiency could not be calculated because no included species were found beyond the selected limits.",
            "initially_included_species_count": initially_included_count,
            "ignored_not_found_species_count": ignored_not_found_count,
            "ignored_not_found_species": ", ".join(ignored_not_found_species),
        }

    combined_note = ""
    if combine_capped_charge_states:
        before_capped = int((df["type"] == "Capped").sum()) if "type" in df.columns else 0
        df = combine_capped_charge_state_auc_table(df)
        after_capped = int((df["type"] == "Capped").sum()) if "type" in df.columns else 0
        if before_capped != after_capped:
            combined_note = f"\nCapped charge states combined for efficiency: {before_capped} rows -> {after_capped} capped species"
        else:
            combined_note = "\nCapped charge-state combination enabled; no eligible capped z=1/z=2 pair was found"

    capped_aucs = pd.to_numeric(df.loc[df["type"] == "Capped", "AUC"], errors="coerce")
    uncapped_aucs = pd.to_numeric(df.loc[df["type"] == "Uncapped", "AUC"], errors="coerce")

    if not ignore_not_found and (capped_aucs.isna().any() or uncapped_aucs.isna().any()):
        return {
            "ok": False,
            "message": "Capping efficiency could not be calculated. Save valid RT windows/AUCs for all included Capped and Uncapped species in the current mode, or enable ignoring species not found beyond limits.",
        }

    capped_aucs = capped_aucs.dropna()
    uncapped_aucs = uncapped_aucs.dropna()
    capped_sum = float(capped_aucs.sum()) if len(capped_aucs) > 0 else 0.0
    uncapped_sum = float(uncapped_aucs.sum()) if len(uncapped_aucs) > 0 else 0.0

    ignored_line = ""
    if ignored_not_found_count > 0:
        ignored_line = f"\nIgnored species not found beyond limits: {ignored_not_found_count}"
        if ignored_not_found_species:
            ignored_line += f" ({', '.join(ignored_not_found_species)})"

    capped_present = len(capped_aucs) > 0 and capped_sum > 0
    uncapped_present = len(uncapped_aucs) > 0 and uncapped_sum > 0

    if capped_present and not uncapped_present:
        if allow_bounded_efficiency:
            return {
                "ok": True,
                "bounded_result": True,
                "bound_type": "greater_than",
                "capping_efficiency_percent": float(lower_bound_percent),
                "capped_auc_sum": capped_sum,
                "uncapped_auc_sum": 0.0,
                "included_species_count": int(len(df)),
                "initially_included_species_count": initially_included_count,
                "ignored_not_found_species_count": ignored_not_found_count,
                "ignored_not_found_species": ", ".join(ignored_not_found_species),
                "combine_capped_charge_states": bool(combine_capped_charge_states),
                "message": f"Capping efficiency: >{float(lower_bound_percent):.0f}%\nNo included uncapped species were found beyond the selected limits.\nCapped AUC sum: {capped_sum:.4g}\nUncapped AUC sum: 0{ignored_line}",
            }
        return {"ok": False, "message": "No included uncapped species were found beyond the selected limits." + ignored_line}

    if uncapped_present and not capped_present:
        if allow_bounded_efficiency:
            return {
                "ok": True,
                "bounded_result": True,
                "bound_type": "less_than",
                "capping_efficiency_percent": float(upper_bound_percent),
                "capped_auc_sum": 0.0,
                "uncapped_auc_sum": uncapped_sum,
                "included_species_count": int(len(df)),
                "initially_included_species_count": initially_included_count,
                "ignored_not_found_species_count": ignored_not_found_count,
                "ignored_not_found_species": ", ".join(ignored_not_found_species),
                "combine_capped_charge_states": bool(combine_capped_charge_states),
                "message": f"Capping efficiency: <{float(upper_bound_percent):.0f}%\nNo included capped species were found beyond the selected limits.\nCapped AUC sum: 0\nUncapped AUC sum: {uncapped_sum:.4g}{ignored_line}",
            }
        return {"ok": False, "message": "No included capped species were found beyond the selected limits." + ignored_line}

    denom = capped_sum + uncapped_sum
    if denom <= 0:
        return {"ok": False, "message": "Capping efficiency could not be calculated because total included AUC is zero." + ignored_line}
    eff = capped_sum / denom
    return {
        "ok": True,
        "bounded_result": False,
        "capping_efficiency_fraction": eff,
        "capping_efficiency_percent": 100.0 * eff,
        "capped_auc_sum": capped_sum,
        "uncapped_auc_sum": uncapped_sum,
        "included_species_count": int(len(df)),
        "initially_included_species_count": initially_included_count,
        "ignored_not_found_species_count": ignored_not_found_count,
        "ignored_not_found_species": ", ".join(ignored_not_found_species),
        "combine_capped_charge_states": bool(combine_capped_charge_states),
        "message": f"Capping efficiency: {100.0 * eff:.2f}%\nCapped AUC sum: {capped_sum:.4g}\nUncapped AUC sum: {uncapped_sum:.4g}{ignored_line}",
    }
