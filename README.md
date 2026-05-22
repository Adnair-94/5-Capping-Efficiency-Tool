# RNA Forge 5′ Capping Efficiency Tool — Desktop Version v4

Standalone Python/PySide6 desktop implementation of the RNA Forge 5′ capping efficiency workflow.

## v4 fix
The EIC viewer now loads Plotly HTML from a temporary local HTML file instead of `QWebEngineView.setHtml()`.
This avoids Qt WebEngine data-URL size limits that can prevent large inline Plotly figures from displaying after PyInstaller packaging.

## Build
Use GitHub Actions or run:

```bat
build_windows.bat
```

## Notes
- The scientific logic is in `capping_core.py`.
- The GUI is in `main.py`.
- Export functions are in `capping_exports.py`.
- Plot export functions are in `capping_plots.py`.
