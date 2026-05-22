"""PySide6 desktop GUI for the RNA Forge 5' Capping Efficiency Tool."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QFileDialog,
    QFormLayout,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QProgressDialog,
    QSpinBox,
    QDoubleSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.widgets import SpanSelector

from capping_core import (
    SPECIES_COLUMNS,
    AnalysisSession,
    auc_table,
    build_analysis,
    calculate_efficiency,
    default_species,
    generate_custom_species,
    generate_preset_species,
    cap_presets,
    read_species_file,
    save_rt_window,
    species_template,
    validate_species_table,
)
from capping_exports import (
    export_all_eic_tiffs_zip,
    export_csv,
    export_current_eic_tiff,
    export_excel_workbook,
)
from capping_plots import plot_eic


class CappingEfficiencyWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("RNA Forge 5' Capping Efficiency Tool")
        self.resize(1500, 900)

        self.ms_file: Optional[str] = None
        self.species_df: pd.DataFrame = default_species()
        self.session: Optional[AnalysisSession] = None
        self.selected_window: Optional[Tuple[float, float]] = None

        self._build_ui()
        self._load_species_table(self.species_df)
        self._set_status("Using default species table. You can edit it manually or upload/generate a species table.")

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------
    def _build_ui(self):
        root = QWidget()
        self.setCentralWidget(root)
        main_layout = QHBoxLayout(root)

        left = QWidget()
        left_layout = QVBoxLayout(left)
        left.setMaximumWidth(650)
        main_layout.addWidget(left)

        right = QWidget()
        right_layout = QVBoxLayout(right)
        main_layout.addWidget(right, stretch=1)

        # MS file controls
        file_group = QGroupBox("MS input")
        file_layout = QVBoxLayout(file_group)
        self.ms_label = QLabel("No MS file selected")
        choose_ms_btn = QPushButton("Choose .mzML or .mzXML file")
        choose_ms_btn.clicked.connect(self.choose_ms_file)
        file_layout.addWidget(choose_ms_btn)
        file_layout.addWidget(self.ms_label)
        left_layout.addWidget(file_group)

        # Analysis options
        opts_group = QGroupBox("Analysis options")
        opts_form = QFormLayout(opts_group)
        self.default_mz_tol = QDoubleSpinBox()
        self.default_mz_tol.setDecimals(5)
        self.default_mz_tol.setRange(0.00001, 10.0)
        self.default_mz_tol.setSingleStep(0.001)
        self.default_mz_tol.setValue(0.01)
        self.default_ms2_window = QDoubleSpinBox()
        self.default_ms2_window.setDecimals(1)
        self.default_ms2_window.setRange(0.1, 10000.0)
        self.default_ms2_window.setSingleStep(1.0)
        self.default_ms2_window.setValue(30.0)
        self.mode_combo = QComboBox()
        self.mode_combo.addItem("MS2-confirmed Only", "filtered")
        self.mode_combo.addItem("All Species", "all")
        self.mode_combo.currentIndexChanged.connect(self.refresh_after_mode_change)
        opts_form.addRow("Default m/z tolerance (Th)", self.default_mz_tol)
        opts_form.addRow("Default MS2 window (+/- seconds)", self.default_ms2_window)
        opts_form.addRow("Species mode", self.mode_combo)
        left_layout.addWidget(opts_group)

        # Preset generator
        gen_group = QGroupBox("Species generator")
        gen_layout = QGridLayout(gen_group)
        self.preset_combo = QComboBox()
        presets = cap_presets()
        for _, row in presets.iterrows():
            self.preset_combo.addItem(str(row["label"]), str(row["preset_id"]))
        self.include_capped_z2 = QCheckBox("Include capped z = 2 ion")
        self.include_capped_z2.setChecked(True)
        self.include_uncapped_z2 = QCheckBox("Include uncapped z = 2 ladder")
        self.include_g_ladder = QCheckBox("Include G-only pppG/ppG/pG ladder")
        self.include_g_ladder.setChecked(True)
        self.g_as_denominator = QCheckBox("Include G-only ladder in efficiency denominator")
        self.generated_mode = QComboBox()
        self.generated_mode.addItem("Replace table", "replace")
        self.generated_mode.addItem("Append to table", "append")
        load_preset_btn = QPushButton("Load preset species")
        load_preset_btn.clicked.connect(self.load_preset_species)
        gen_layout.addWidget(QLabel("Preset"), 0, 0)
        gen_layout.addWidget(self.preset_combo, 0, 1, 1, 2)
        gen_layout.addWidget(self.include_capped_z2, 1, 0, 1, 3)
        gen_layout.addWidget(self.include_uncapped_z2, 2, 0, 1, 3)
        gen_layout.addWidget(self.include_g_ladder, 3, 0, 1, 3)
        gen_layout.addWidget(self.g_as_denominator, 4, 0, 1, 3)
        gen_layout.addWidget(QLabel("Action"), 5, 0)
        gen_layout.addWidget(self.generated_mode, 5, 1)
        gen_layout.addWidget(load_preset_btn, 5, 2)
        left_layout.addWidget(gen_group)

        # Custom generator (compact)
        custom_group = QGroupBox("Custom cap species generator")
        custom_layout = QGridLayout(custom_group)
        self.custom_cap_name = QLineEdit("CustomCap")
        self.custom_capped_mz = QDoubleSpinBox(); self.custom_capped_mz.setDecimals(5); self.custom_capped_mz.setRange(0, 100000); self.custom_capped_mz.setValue(0)
        self.custom_capped_charge = QSpinBox(); self.custom_capped_charge.setRange(1, 20); self.custom_capped_charge.setValue(1)
        self.custom_fragment = QComboBox(); self.custom_fragment.addItems(["AG", "AUG", "G", "Custom"])
        self.custom_ppp_mz = QDoubleSpinBox(); self.custom_ppp_mz.setDecimals(5); self.custom_ppp_mz.setRange(0, 100000); self.custom_ppp_mz.setValue(0)
        self.custom_ppp_charge = QSpinBox(); self.custom_ppp_charge.setRange(1, 20); self.custom_ppp_charge.setValue(1)
        self.custom_include_capped_z2 = QCheckBox("Include calculated capped z = 2 ion"); self.custom_include_capped_z2.setChecked(True)
        self.custom_include_uncapped_z2 = QCheckBox("Include uncapped z = 2 ladder")
        self.custom_include_g_ladder = QCheckBox("Include G-only ladder")
        self.custom_g_as_denominator = QCheckBox("Include G-only ladder in denominator")
        gen_custom_btn = QPushButton("Generate custom species"); gen_custom_btn.clicked.connect(self.generate_custom_species)
        custom_layout.addWidget(QLabel("Cap label"), 0, 0); custom_layout.addWidget(self.custom_cap_name, 0, 1)
        custom_layout.addWidget(QLabel("Capped m/z"), 1, 0); custom_layout.addWidget(self.custom_capped_mz, 1, 1)
        custom_layout.addWidget(QLabel("Capped z"), 2, 0); custom_layout.addWidget(self.custom_capped_charge, 2, 1)
        custom_layout.addWidget(QLabel("5' T1 fragment"), 3, 0); custom_layout.addWidget(self.custom_fragment, 3, 1)
        custom_layout.addWidget(QLabel("Custom ppp m/z"), 4, 0); custom_layout.addWidget(self.custom_ppp_mz, 4, 1)
        custom_layout.addWidget(QLabel("Custom ppp z"), 5, 0); custom_layout.addWidget(self.custom_ppp_charge, 5, 1)
        custom_layout.addWidget(self.custom_include_capped_z2, 6, 0, 1, 2)
        custom_layout.addWidget(self.custom_include_uncapped_z2, 7, 0, 1, 2)
        custom_layout.addWidget(self.custom_include_g_ladder, 8, 0, 1, 2)
        custom_layout.addWidget(self.custom_g_as_denominator, 9, 0, 1, 2)
        custom_layout.addWidget(gen_custom_btn, 10, 0, 1, 2)
        left_layout.addWidget(custom_group)

        # Species table controls
        table_group = QGroupBox("Species table")
        table_layout = QVBoxLayout(table_group)
        table_buttons = QHBoxLayout()
        upload_btn = QPushButton("Upload species table")
        upload_btn.clicked.connect(self.upload_species_table)
        template_btn = QPushButton("Save species template")
        template_btn.clicked.connect(self.save_species_template)
        reset_btn = QPushButton("Reset default")
        reset_btn.clicked.connect(self.reset_species)
        add_btn = QPushButton("Add row")
        add_btn.clicked.connect(self.add_species_row)
        remove_btn = QPushButton("Remove selected row")
        remove_btn.clicked.connect(self.remove_selected_species_row)
        for btn in [upload_btn, template_btn, reset_btn, add_btn, remove_btn]:
            table_buttons.addWidget(btn)
        table_layout.addLayout(table_buttons)
        self.species_table = QTableWidget()
        self.species_table.setColumnCount(len(SPECIES_COLUMNS))
        self.species_table.setHorizontalHeaderLabels(SPECIES_COLUMNS)
        table_layout.addWidget(self.species_table, stretch=1)
        left_layout.addWidget(table_group, stretch=1)

        # Run/export controls
        run_group = QGroupBox("Run and export")
        run_layout = QGridLayout(run_group)
        load_btn = QPushButton("Load Data & Build EICs")
        load_btn.clicked.connect(self.load_ms_data)
        save_rt_btn = QPushButton("Export RT windows CSV")
        save_rt_btn.clicked.connect(self.export_rt_csv)
        save_auc_btn = QPushButton("Export AUC CSV")
        save_auc_btn.clicked.connect(self.export_auc_csv)
        save_excel_btn = QPushButton("Export Excel workbook")
        save_excel_btn.clicked.connect(self.export_excel)
        save_tiff_btn = QPushButton("Export current EIC TIFF")
        save_tiff_btn.clicked.connect(self.export_current_tiff)
        save_zip_btn = QPushButton("Export all EIC TIFFs ZIP")
        save_zip_btn.clicked.connect(self.export_tiff_zip)
        run_layout.addWidget(load_btn, 0, 0, 1, 2)
        run_layout.addWidget(save_rt_btn, 1, 0)
        run_layout.addWidget(save_auc_btn, 1, 1)
        run_layout.addWidget(save_excel_btn, 2, 0, 1, 2)
        run_layout.addWidget(save_tiff_btn, 3, 0)
        run_layout.addWidget(save_zip_btn, 3, 1)
        left_layout.addWidget(run_group)

        self.status_box = QTextEdit()
        self.status_box.setReadOnly(True)
        self.status_box.setMaximumHeight(110)
        left_layout.addWidget(self.status_box)

        # Right side: species selector, plot, tables, efficiency
        top_right = QHBoxLayout()
        self.selected_species_combo = QComboBox()
        self.selected_species_combo.currentIndexChanged.connect(self.plot_selected_species)
        self.save_window_btn = QPushButton("Save RT Window for Selected Species")
        self.save_window_btn.clicked.connect(self.save_selected_window)
        top_right.addWidget(QLabel("Selected species"))
        top_right.addWidget(self.selected_species_combo, stretch=1)
        top_right.addWidget(self.save_window_btn)
        right_layout.addLayout(top_right)

        self.figure = Figure(figsize=(7, 4.5))
        self.canvas = FigureCanvas(self.figure)
        self.ax = self.figure.add_subplot(111)
        self.toolbar = NavigationToolbar(self.canvas, self)
        right_layout.addWidget(self.toolbar)
        right_layout.addWidget(self.canvas, stretch=2)
        self.span = SpanSelector(self.ax, self._on_span_select, "horizontal", useblit=True, props=dict(alpha=0.2))

        self.rt_table = QTableWidget()
        self.auc_table_widget = QTableWidget()
        tables_layout = QHBoxLayout()
        rt_group = QGroupBox("RT windows")
        rt_layout = QVBoxLayout(rt_group); rt_layout.addWidget(self.rt_table)
        auc_group = QGroupBox("AUC table")
        auc_layout = QVBoxLayout(auc_group); auc_layout.addWidget(self.auc_table_widget)
        tables_layout.addWidget(rt_group)
        tables_layout.addWidget(auc_group)
        right_layout.addLayout(tables_layout, stretch=1)

        self.efficiency_box = QTextEdit()
        self.efficiency_box.setReadOnly(True)
        self.efficiency_box.setMaximumHeight(100)
        right_layout.addWidget(self.efficiency_box)

    # ------------------------------------------------------------------
    # Table helpers
    # ------------------------------------------------------------------
    def _set_status(self, text: str):
        self.status_box.setPlainText(str(text))

    def _message(self, title: str, text: str, icon=QMessageBox.Information):
        box = QMessageBox(self)
        box.setIcon(icon)
        box.setWindowTitle(title)
        box.setText(text)
        box.exec()

    def _df_to_table(self, table: QTableWidget, df: pd.DataFrame):
        if df is None:
            df = pd.DataFrame()
        table.blockSignals(True)
        table.setRowCount(len(df))
        table.setColumnCount(len(df.columns))
        table.setHorizontalHeaderLabels([str(c) for c in df.columns])
        for r in range(len(df)):
            for c, col in enumerate(df.columns):
                val = df.iloc[r, c]
                text = "" if pd.isna(val) else str(val)
                item = QTableWidgetItem(text)
                table.setItem(r, c, item)
        table.blockSignals(False)
        table.resizeColumnsToContents()

    def _species_table_to_df(self) -> pd.DataFrame:
        cols = [self.species_table.horizontalHeaderItem(c).text() for c in range(self.species_table.columnCount())]
        rows = []
        for r in range(self.species_table.rowCount()):
            row = {}
            for c, col in enumerate(cols):
                item = self.species_table.item(r, c)
                row[col] = "" if item is None else item.text()
            rows.append(row)
        return pd.DataFrame(rows, columns=cols)

    def _load_species_table(self, df: pd.DataFrame):
        self.species_df = df.copy()
        for col in SPECIES_COLUMNS:
            if col not in self.species_df.columns:
                self.species_df[col] = np.nan if col not in {"source", "notes"} else ""
        self.species_df = self.species_df[SPECIES_COLUMNS]
        self._df_to_table(self.species_table, self.species_df)
        self.session = None
        self.refresh_tables()

    def _apply_species_df(self, df: pd.DataFrame, mode: str, message: str):
        ok, valid, msg = validate_species_table(df)
        if not ok:
            self._set_status("Species validation failed:\n" + msg)
            self._message("Species validation failed", msg, QMessageBox.Critical)
            return
        if mode == "append":
            merged = pd.concat([self._species_table_to_df(), valid], ignore_index=True)
            ok, valid, msg2 = validate_species_table(merged)
            msg = "\n".join(x for x in [msg, msg2] if x)
        self._load_species_table(valid)
        self._set_status(message + ("\n" + msg if msg else ""))

    # ------------------------------------------------------------------
    # Species actions
    # ------------------------------------------------------------------
    def choose_ms_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "Choose MS file", "", "MS files (*.mzML *.mzXML)")
        if path:
            self.ms_file = path
            self.ms_label.setText(path)
            self.session = None
            self.refresh_tables()
            self._set_status("New MS file selected. Click 'Load Data & Build EICs'.")

    def upload_species_table(self):
        path, _ = QFileDialog.getOpenFileName(self, "Upload species table", "", "Species tables (*.csv *.xlsx *.xls)")
        if not path:
            return
        try:
            raw = read_species_file(path)
            ok, valid, msg = validate_species_table(raw)
            if not ok:
                raise ValueError(msg)
            self._load_species_table(valid)
            self._set_status("Species table uploaded and validated." + ("\n" + msg if msg else ""))
        except Exception as exc:
            self._message("Upload failed", str(exc), QMessageBox.Critical)
            self._set_status("Species table upload failed:\n" + str(exc))

    def save_species_template(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save species template", "capping_species_template.csv", "CSV files (*.csv)")
        if path:
            species_template().to_csv(path, index=False)
            self._set_status(f"Species template saved: {path}")

    def reset_species(self):
        self._load_species_table(default_species())
        self._set_status("Default species table restored. Click 'Load Data & Build EICs' if data were already loaded.")

    def add_species_row(self):
        df = self._species_table_to_df()
        new = {col: "" for col in SPECIES_COLUMNS}
        new.update({"name": f"species_{len(df) + 1}", "type": "Capped", "use_for_efficiency": "TRUE", "source": "manual"})
        df = pd.concat([df, pd.DataFrame([new])], ignore_index=True)
        self._load_species_table(df)
        self._set_status("Species row added. Edit the row, then click 'Load Data & Build EICs'.")

    def remove_selected_species_row(self):
        row = self.species_table.currentRow()
        if row < 0:
            return
        df = self._species_table_to_df().drop(index=row).reset_index(drop=True)
        self._load_species_table(df)
        self._set_status("Species row removed. Click 'Load Data & Build EICs' to rebuild EICs before analysis.")

    def load_preset_species(self):
        try:
            gen = generate_preset_species(
                preset_id=self.preset_combo.currentData(),
                include_capped_z2=self.include_capped_z2.isChecked(),
                include_uncapped_z2=self.include_uncapped_z2.isChecked(),
                include_g_ladder=self.include_g_ladder.isChecked(),
                g_use_for_efficiency=self.g_as_denominator.isChecked(),
            )
            self._apply_species_df(gen, self.generated_mode.currentData(), "Preset species table generated. Review denominator flags before analysis.")
        except Exception as exc:
            self._message("Preset generation failed", str(exc), QMessageBox.Critical)

    def generate_custom_species(self):
        try:
            custom_ppp = self.custom_ppp_mz.value() if self.custom_ppp_mz.value() > 0 else None
            gen = generate_custom_species(
                custom_cap_name=self.custom_cap_name.text(),
                capped_mz=self.custom_capped_mz.value(),
                capped_charge=self.custom_capped_charge.value(),
                fragment=self.custom_fragment.currentText(),
                custom_ppp_mz=custom_ppp,
                custom_ppp_charge=self.custom_ppp_charge.value(),
                include_capped_z2=self.custom_include_capped_z2.isChecked(),
                include_uncapped_z2=self.custom_include_uncapped_z2.isChecked(),
                include_g_ladder=self.custom_include_g_ladder.isChecked(),
                g_use_for_efficiency=self.custom_g_as_denominator.isChecked(),
            )
            self._apply_species_df(gen, self.generated_mode.currentData(), "Custom species table generated. Review the generated uncapped ladder before analysis.")
        except Exception as exc:
            self._message("Custom generation failed", str(exc), QMessageBox.Critical)

    # ------------------------------------------------------------------
    # Analysis
    # ------------------------------------------------------------------
    def load_ms_data(self):
        if not self.ms_file:
            self._message("No MS file", "Select a .mzML or .mzXML file first.", QMessageBox.Warning)
            return
        species = self._species_table_to_df()
        ok, valid, msg = validate_species_table(species)
        if not ok:
            self._message("Species validation failed", msg, QMessageBox.Critical)
            self._set_status("Species table validation failed:\n" + msg)
            return
        self._load_species_table(valid)

        progress = QProgressDialog("Loading MS file and extracting EICs...", "Cancel", 0, 100, self)
        progress.setWindowModality(Qt.WindowModal)
        progress.setMinimumDuration(0)

        def progress_cb(frac, text):
            progress.setValue(max(0, min(100, int(frac * 100))))
            progress.setLabelText(text)
            QApplication.processEvents()
            if progress.wasCanceled():
                raise RuntimeError("Analysis cancelled by user.")

        try:
            session, validation_msg = build_analysis(
                self.ms_file,
                self.species_df,
                default_mz_tol=self.default_mz_tol.value(),
                default_ms2_window=self.default_ms2_window.value(),
                progress_cb=progress_cb,
            )
            progress.setValue(100)
            self.session = session
            status_parts = ["MS file loaded and EICs built.", session.filter_message]
            if validation_msg:
                status_parts.append(validation_msg)
            self._set_status("\n".join(status_parts))
            self.refresh_after_mode_change()
        except Exception as exc:
            progress.close()
            self._message("Analysis failed", str(exc), QMessageBox.Critical)
            self._set_status("Analysis failed:\n" + str(exc))

    def requested_mode(self) -> str:
        return self.mode_combo.currentData() or "filtered"

    def refresh_after_mode_change(self):
        self.refresh_species_choices()
        self.refresh_tables()
        self.plot_selected_species()

    def refresh_species_choices(self):
        self.selected_species_combo.blockSignals(True)
        self.selected_species_combo.clear()
        if self.session and self.session.loaded:
            state = self.session.state(self.requested_mode())
            for name in state.species_table["name"].astype(str).tolist():
                self.selected_species_combo.addItem(name)
        self.selected_species_combo.blockSignals(False)

    def refresh_tables(self):
        if self.session and self.session.loaded:
            state = self.session.state(self.requested_mode())
            self._df_to_table(self.rt_table, state.rt_windows)
            tbl = auc_table(self.session, self.requested_mode())
            self._df_to_table(self.auc_table_widget, tbl)
            eff = calculate_efficiency(tbl, honor_use_for_efficiency=True)
            self.efficiency_box.setPlainText(str(eff.get("message", "")))
        else:
            self._df_to_table(self.rt_table, pd.DataFrame({"Message": ["No MS data loaded."]}))
            self._df_to_table(self.auc_table_widget, pd.DataFrame({"Message": ["No AUC table available."]}))
            self.efficiency_box.setPlainText("No MS data loaded.")

    def _on_span_select(self, xmin, xmax):
        self.selected_window = (float(xmin), float(xmax))
        self._set_status(f"Selected RT window: [{xmin:.2f}, {xmax:.2f}] seconds. Click 'Save RT Window for Selected Species'.")
        self.plot_selected_species()

    def plot_selected_species(self):
        self.ax.clear()
        if not (self.session and self.session.loaded) or self.selected_species_combo.count() == 0:
            self.ax.set_title("No MS data loaded")
            self.ax.set_xlabel("Retention Time (seconds)")
            self.ax.set_ylabel("Intensity")
            self.canvas.draw_idle()
            return
        species = self.selected_species_combo.currentText()
        state = self.session.state(self.requested_mode())
        eic = state.eic_list.get(species, pd.DataFrame())
        rt_window = self.selected_window
        # Use stored window if no new unsaved span has been selected.
        if rt_window is None and state.rt_windows is not None and len(state.rt_windows) > 0:
            row = state.rt_windows[state.rt_windows["Species"].astype(str) == species]
            if len(row) == 1 and pd.notna(row.iloc[0]["RTmin"]) and pd.notna(row.iloc[0]["RTmax"]):
                rt_window = (float(row.iloc[0]["RTmin"]), float(row.iloc[0]["RTmax"]))
        plot_eic(self.ax, eic, species, rt_window=rt_window)
        self.canvas.draw_idle()

    def save_selected_window(self):
        if not (self.session and self.session.loaded):
            self._message("No data", "Load MS data first.", QMessageBox.Warning)
            return
        if self.selected_window is None:
            self._message("No RT window", "Drag across the plot to select an RT window first.", QMessageBox.Warning)
            return
        species = self.selected_species_combo.currentText()
        try:
            auc = save_rt_window(self.session, self.requested_mode(), species, self.selected_window[0], self.selected_window[1])
            self._set_status(f"Window for {species} saved: [{self.selected_window[0]:.2f}, {self.selected_window[1]:.2f}] seconds; AUC = {auc:.4g}")
            self.selected_window = None
            self.refresh_tables()
            self.plot_selected_species()
        except Exception as exc:
            self._message("Could not save window", str(exc), QMessageBox.Critical)

    # ------------------------------------------------------------------
    # Exports
    # ------------------------------------------------------------------
    def _require_session(self) -> bool:
        if not (self.session and self.session.loaded):
            self._message("No data", "Load MS data first.", QMessageBox.Warning)
            return False
        return True

    def export_rt_csv(self):
        if not self._require_session():
            return
        mode = self.session.current_mode(self.requested_mode())
        stem = Path(self.ms_file).stem if self.ms_file else "rt_windows"
        path, _ = QFileDialog.getSaveFileName(self, "Export RT windows", f"{stem}_{mode}_rt_windows.csv", "CSV files (*.csv)")
        if path:
            export_csv(self.session.state(self.requested_mode()).rt_windows, path)
            self._set_status(f"RT windows exported: {path}")

    def export_auc_csv(self):
        if not self._require_session():
            return
        mode = self.session.current_mode(self.requested_mode())
        stem = Path(self.ms_file).stem if self.ms_file else "auc_table"
        path, _ = QFileDialog.getSaveFileName(self, "Export AUC table", f"{stem}_{mode}_auc_table.csv", "CSV files (*.csv)")
        if path:
            export_csv(auc_table(self.session, self.requested_mode()), path)
            self._set_status(f"AUC table exported: {path}")

    def export_excel(self):
        if not self._require_session():
            return
        mode = self.session.current_mode(self.requested_mode())
        stem = Path(self.ms_file).stem if self.ms_file else "capping_efficiency"
        path, _ = QFileDialog.getSaveFileName(self, "Export Excel workbook", f"{stem}_{mode}_capping_efficiency.xlsx", "Excel files (*.xlsx)")
        if path:
            export_excel_workbook(self.session, self.requested_mode(), path)
            self._set_status(f"Excel workbook exported: {path}")

    def export_current_tiff(self):
        if not self._require_session():
            return
        species = self.selected_species_combo.currentText()
        safe = "".join(c if c.isalnum() or c in "._-" else "_" for c in species)
        path, _ = QFileDialog.getSaveFileName(self, "Export current EIC TIFF", f"{safe}_EIC.tiff", "TIFF files (*.tiff *.tif)")
        if path:
            export_current_eic_tiff(self.session, self.requested_mode(), species, path)
            self._set_status(f"Current EIC TIFF exported: {path}")

    def export_tiff_zip(self):
        if not self._require_session():
            return
        mode = self.session.current_mode(self.requested_mode())
        stem = Path(self.ms_file).stem if self.ms_file else "eic_plots"
        path, _ = QFileDialog.getSaveFileName(self, "Export all EIC TIFFs ZIP", f"{stem}_{mode}_EIC_plots.zip", "ZIP files (*.zip)")
        if path:
            export_all_eic_tiffs_zip(self.session, self.requested_mode(), path)
            self._set_status(f"EIC TIFF ZIP exported: {path}")


def main():
    app = QApplication(sys.argv)
    win = CappingEfficiencyWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
