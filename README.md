# RNA Forge 5' Capping Efficiency Tool — desktop version

This is a Python/PySide6 desktop translation of the Shiny app.R logic for LC-MS-based 5' capping efficiency analysis.

## Version 2 notes

This version fixes three practical problems from the first desktop port:

1. The GUI is reorganized into panes and tabs instead of a compressed single sidebar.
2. The EIC viewer uses embedded Plotly via Qt WebEngine, so box/lasso selection is closer to the original Shiny/plotly workflow.
3. PyInstaller build commands explicitly bundle Pyteomics, Plotly, and PySide6 WebEngine components.

## Install for local development

```bash
python -m pip install -r requirements.txt
python main.py
```

## Build on Windows

```bat
build_windows.bat
```

The `.exe` will be created in the `dist` folder.

## Recommended GitHub Actions build

The workflow file is included at:

```text
.github/workflows/build-windows-exe.yml
```

Push the files to GitHub, open the Actions tab, and run **Build Windows EXE**.

## Important scientific behavior

Capping efficiency is calculated as:

```text
sum(AUC_capped) / [sum(AUC_capped) + sum(AUC_uncapped)]
```

This Python implementation honors the `use_for_efficiency` flag. Species marked as diagnostic-only, such as the G-only ladder by default, are not included in the efficiency denominator unless `use_for_efficiency = TRUE`.

## mzML/mzXML support

mzML and mzXML parsing is handled by the bundled pure-Python reader in `capping_core.py`; no external MS parser is required at runtime.


## v3 correction

This version removes the Pyteomics runtime dependency and uses a bundled pure-Python mzML/mzXML reader. This avoids the previous Windows executable failure where the app opened but could not parse MS files because Pyteomics was not bundled.
