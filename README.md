# RNA Forge 5′ Capping Efficiency Desktop Tool

This is a Python/PySide6 desktop translation of the working single-file Shiny `app.R` for LC-MS-based 5′ capping efficiency analysis.

## What is preserved from app.R

- `.mzML` and `.mzXML` input
- user-defined capped and uncapped species table
- CSV/XLSX species-table upload and validation
- preset and custom species generation
- per-row `mz_tol`, optional `ms2_mz`, and `ms2_window_sec`
- MS1 EIC extraction
- optional MS2 confirmation near the MS1 apex
- interactive RT-window selection
- trapezoidal MS1 AUC integration
- capping efficiency calculation from summed capped and uncapped AUCs
- RT-window CSV export and AUC-table CSV export

## Important implementation note

The pasted app.R creates a `use_for_efficiency` flag and labels the G-only ladder as diagnostic by default, but the R `capEff` block does not actually use the flag when calculating capping efficiency. This desktop version honors `use_for_efficiency` by default because that is the intended denominator-control behavior. The core function `calculate_efficiency(..., honor_use_for_efficiency=False)` can reproduce the literal R behavior if needed.

## Development install

```bat
python -m venv .venv
.venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python main.py
```

## Build onefile Windows executable

```bat
build_windows.bat
```

The executable will be created in `dist\RNA_Forge_Capping_Efficiency_Tool.exe`.

## Clean-machine testing order

1. Test `capping_core.py` alone using a known mzML/mzXML and species table.
2. Compare numerical AUC and capping efficiency against the Shiny app.R output.
3. Compare EIC plots visually against app.R.
4. Test CSV, Excel, TIFF, and ZIP exports.
5. Test the GUI workflow on the build machine.
6. Build the onefile `.exe`.
7. Test the `.exe` on the build machine.
8. Test the `.exe` on a clean Windows machine with no R/Python installation.
