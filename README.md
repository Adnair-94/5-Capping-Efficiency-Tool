# RNA Forge 5' Capping Efficiency Tool — desktop version v13

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
- Saves RT windows either from Plotly box/lasso selection or from numeric RT min/max input.
- Integrates AUC using only points above the user-set AUC point minimum intensity.
- Requires a user-set minimum number of contiguous above-threshold EIC points before a selected window is accepted as a valid peak.
- Allows users to include/exclude individual species from final efficiency and composite plotting through the Analysis/plot selection table.
- Flags co-eluted lower-phosphate uncapped species as likely labile phosphate-loss ions and excludes them from the effective efficiency denominator when enabled.
- Calculates capping efficiency as:

  capped AUC sum / (capped AUC sum + uncapped AUC sum)

  while honoring `use_for_efficiency`, species selection, not-found exclusion, and labile-loss exclusion flags.

## v13 changes

- Added numeric RT-window entry as an alternative to Plotly box/lasso selection.
- Added an **Analysis/plot selection** table. Uncheck species/signals that should be ignored in the efficiency calculation and final composite plot.
- Added a **Composite plot** tab with editable x-axis and y-axis ranges.
- Composite x-axis values use the current RT display unit, seconds or minutes.
- Final composite Plotly and TIFF output now reflect only species selected/included for analysis.
- Axis controls allow users to hide unrelated RT regions or unwanted peaks in the final composite plot without changing the stored AUC calculation.

## Existing retained behavior

- RT display can switch between seconds and minutes.
- Plotly x-axis range slider is retained for RT navigation.
- Optional exclusion of co-eluted labile phosphate-loss species is retained.
- Species that do not pass MS1/AUC/contiguous-point limits are retained in tables but marked as not found beyond limits.
- Optional automatic exclusion of not-found species from the efficiency calculation is retained.
- Optional bounded calls are retained: `>99%` when capped species are detected but no included uncapped species are found, and `<1%` when uncapped species are detected but no included capped species are found.
- Optional capped z=1/z=2 combination is retained for final efficiency and composite plotting.

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
- `species_found_beyond_limit`
- `species_detection_status`
- `included_for_efficiency_final`

## Build on Windows

Run:

```bat
build_windows.bat
```

or use the included GitHub Actions workflow.


## v14 notes

This package is the GitHub-uploadable merged version: it keeps capped charge-state combination and retains the later controls for numeric RT-window entry, species inclusion/ignore handling, composite EIC axis editing, selected-species composite plotting, missing-species handling, MS1 intensity thresholds, contiguous-point AUC validation, RT unit switching, and labile phosphate-loss flagging.

Upload the contents of this folder to the repository root. Do not upload the ZIP as the only repo content. The workflow file must remain at `.github/workflows/build-windows-exe.yml`.
