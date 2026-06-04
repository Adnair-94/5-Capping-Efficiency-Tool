# RNA Forge 5' Capping Efficiency Tool — desktop version v9

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
- Requires a user-set minimum number of contiguous above-threshold EIC points before a selected window is accepted as a valid peak.
- Flags co-eluted lower-phosphate uncapped species as likely labile phosphate-loss ions and excludes them from the effective efficiency denominator when enabled.
- Calculates capping efficiency as:

  capped AUC sum / (capped AUC sum + uncapped AUC sum)

  while honoring `use_for_efficiency` and labile-loss exclusion flags.

## v9 changes

- Added **Min contiguous points for peak** input. Default = 6.
- AUC is only accepted when the selected RT window contains at least the configured number of contiguous points above the AUC intensity threshold.
- Added Plotly x-axis range slider for RT navigation.
- Added **Ignore co-eluted labile phosphate-loss species** option.
- Added labile phosphate-loss RT tolerance input. Default = 5 seconds.
- Example: if pppG and ppG co-elute and ppG is one phosphate lower than pppG, ppG is flagged as a labile loss and excluded from effective capping-efficiency calculations.
- Added final composite EIC plot for effective included capped and uncapped species.
- Added export of final composite EIC TIFF.

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
- `AUC_min_consecutive_points`
- `AUC_points_used`
- `AUC_segments_used`
- `AUC_max_consecutive_points`
- `AUC_peak_passed_min_points`
- `MS2_near_apex_max_intensity`
- `uncapped_percent_of_capped_sum`
- `uncapped_less_than_threshold`
- `uncapped_threshold_percent`
- `uncapped_threshold_intensity`
- `labile_parent_species`
- `labile_parent_mz`
- `labile_parent_RT_sec`
- `labile_RT_delta_sec`
- `labile_mass_delta_Da`
- `ignored_as_labile_phosphate_loss`
- `ignore_for_efficiency_reason`
- `use_for_efficiency_effective`

## Build on Windows

Run:

```bat
build_windows.bat
```

or use the included GitHub Actions workflow.
