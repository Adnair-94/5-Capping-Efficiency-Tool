"""Export helpers for the RNA Forge 5' Capping Efficiency desktop tool."""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZipFile, ZIP_DEFLATED

import pandas as pd

from capping_core import AnalysisSession, auc_table, calculate_efficiency, intensity_table, composite_eic_data
from capping_plots import save_eic_tiff, save_composite_eic_tiff


def export_csv(df: pd.DataFrame, path: str | Path) -> None:
    df.to_csv(path, index=False)


def export_excel_workbook(
    session: AnalysisSession,
    requested_mode: str,
    path: str | Path,
    ignore_not_found: bool = True,
    allow_bounded_efficiency: bool = True,
) -> None:
    state = session.state(requested_mode)
    auc_tbl = auc_table(session, requested_mode)
    eff = calculate_efficiency(
        auc_tbl,
        ignore_not_found=ignore_not_found,
        allow_bounded_efficiency=allow_bounded_efficiency,
    )
    summary = pd.DataFrame([eff])
    with pd.ExcelWriter(path, engine="xlsxwriter") as writer:
        summary.to_excel(writer, index=False, sheet_name="summary")
        state.rt_windows.to_excel(writer, index=False, sheet_name="rt_windows")
        auc_tbl.to_excel(writer, index=False, sheet_name="auc_table")
        intensity_table(session, requested_mode).to_excel(writer, index=False, sheet_name="intensity_summary")
        state.species_table.to_excel(writer, index=False, sheet_name="species_table")
        state.eic_df.to_excel(writer, index=False, sheet_name="eic_points")


def export_current_eic_tiff(session: AnalysisSession, requested_mode: str, species_name: str, path: str | Path, dpi: int = 300, rt_unit: str = "seconds") -> None:
    state = session.state(requested_mode)
    eic = state.eic_list.get(species_name)
    rt_window = None
    if state.rt_windows is not None and len(state.rt_windows) > 0:
        row = state.rt_windows[state.rt_windows["Species"].astype(str) == str(species_name)]
        if len(row) == 1 and pd.notna(row.iloc[0]["RTmin"]) and pd.notna(row.iloc[0]["RTmax"]):
            rt_window = (float(row.iloc[0]["RTmin"]), float(row.iloc[0]["RTmax"]))
    save_eic_tiff(eic, species_name, path, rt_window=rt_window, dpi=dpi, rt_unit=rt_unit)


def export_all_eic_tiffs_zip(session: AnalysisSession, requested_mode: str, zip_path: str | Path, dpi: int = 300, rt_unit: str = "seconds") -> None:
    state = session.state(requested_mode)
    with TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        files = []
        for species_name, eic in state.eic_list.items():
            safe = "".join(c if c.isalnum() or c in "._-" else "_" for c in species_name)
            out = tmpdir_path / f"{safe}_EIC.tiff"
            rt_window = None
            row = state.rt_windows[state.rt_windows["Species"].astype(str) == str(species_name)] if state.rt_windows is not None else pd.DataFrame()
            if len(row) == 1 and pd.notna(row.iloc[0]["RTmin"]) and pd.notna(row.iloc[0]["RTmax"]):
                rt_window = (float(row.iloc[0]["RTmin"]), float(row.iloc[0]["RTmax"]))
            save_eic_tiff(eic, species_name, out, rt_window=rt_window, dpi=dpi, rt_unit=rt_unit)
            files.append(out)
        with ZipFile(zip_path, "w", compression=ZIP_DEFLATED) as zf:
            for file in files:
                zf.write(file, arcname=file.name)



def export_composite_eic_tiff(
    session: AnalysisSession,
    requested_mode: str,
    path: str | Path,
    dpi: int = 300,
    rt_unit: str = "seconds",
    combine_capped_charge_states: bool = False,
) -> None:
    """Export a final overlay plot of included capped and uncapped EICs.

    Species flagged as ignored labile phosphate-loss products are omitted through the
    use_for_efficiency_effective/final inclusion logic. If requested, capped z=1/z=2
    traces with matching neutral mass are summed for the composite plot only.
    """
    state = session.state(requested_mode)
    eic_map, species_table = composite_eic_data(
        session,
        requested_mode,
        combine_capped_charge_states=combine_capped_charge_states,
    )
    save_composite_eic_tiff(eic_map, species_table, path, rt_windows=state.rt_windows, dpi=dpi, rt_unit=rt_unit)
