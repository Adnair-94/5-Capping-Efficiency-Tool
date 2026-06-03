# RNA Forge 5' Capping Efficiency Tool — desktop version

Standalone Python/PySide6 desktop implementation of the RNA Forge 5' capping efficiency workflow.

## Supported MS input files

- `.mzML`
- `.mzXML`

MGF support is intentionally disabled in this version. MGF files are usually MS2 peak-list exports and are not suitable for MS1 EIC/AUC-based capping-efficiency quantitation unless produced in a non-standard MS1-containing form.

## Main analysis behavior

- Reads user-defined capped and uncapped species from the editable table or uploaded species table.
- Extracts MS1 EICs for each species.
- Optionally filters species using MS2 signal near the MS1 apex.
- Reports MS1 maximum intensity and apex retention time for each species.
- Reports an MS1 confirmation message based on the user-set MS1 intensity threshold.
- Saves user-selected RT windows from the Plotly EIC view.
- Integrates AUC using only points above the user-set AUC point minimum intensity.
- Calculates capping efficiency as:

  capped AUC sum / (capped AUC sum + uncapped AUC sum)

  while honoring `use_for_efficiency`.

## Reporting fields

The AUC and intensity-summary exports include:

- `MS1_max_intensity`
- `MS1_apex_RT_sec`
- `MS1_apex_RT_min`
- `RT_at_MS1_max_sec`
- `MS1_confirmed`
- `MS1_threshold`
- `RTmin`
- `RTmax`
- `RTmin_min`
- `RTmax_min`
- `AUC_min_intensity`
- `MS2_near_apex_max_intensity`
- `uncapped_percent_of_capped_sum`
- `uncapped_less_than_threshold`
- `uncapped_threshold_percent`
- `uncapped_threshold_intensity`

## v8 patch

- Added an RT display-unit selector: seconds or minutes.
- EIC x-axis can now be displayed in minutes.
- Plotly box/lasso selections made in minutes are converted back to seconds for stored RT windows and AUC integration.
- Exported EIC TIFFs use the selected RT display unit.
- RT-window tables include both second and minute columns.

## Build on Windows

Run:

```bat
build_windows.bat
```

or use the included GitHub Actions workflow.
