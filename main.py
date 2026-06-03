"""PySide6 desktop GUI for the RNA Forge 5' Capping Efficiency Tool.

Version 8 changes:
- Adds an RT display-unit selector for EIC plots: seconds or minutes.
- Plotly selections made in minutes are converted back to seconds before AUC integration.
- Exported TIFF EIC plots use the selected RT display unit.
- Stored analysis values remain in seconds, with minute columns reported in output tables.
"""

from __future__ import annotations

import json
import sys
import tempfile
import html as html_lib
from pathlib import Path
from typing import Optional

import pandas as pd

from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import QColor, QPalette
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
    QScrollArea,
    QSpinBox,
    QDoubleSpinBox,
    QSplitter,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

try:
    from PySide6.QtWebEngineWidgets import QWebEngineView
    WEBENGINE_AVAILABLE = True
except Exception:  # pragma: no cover - only used when QtWebEngine is unavailable
    QWebEngineView = None
    WEBENGINE_AVAILABLE = False

from capping_core import (
    SPECIES_COLUMNS,
    AnalysisSession,
    auc_table,
    intensity_table,
    build_analysis,
    calculate_efficiency,
    cap_presets,
    default_species,
    generate_custom_species,
    generate_preset_species,
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


APP_NAME = "RNA Forge 5' Capping Efficiency Tool"


def safe_filename(text: str) -> str:
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in str(text))


class CappingEfficiencyWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(APP_NAME)
        self.resize(1650, 950)

        self.ms_file: Optional[str] = None
        self.species_df: pd.DataFrame = default_species()
        self.session: Optional[AnalysisSession] = None
        self._plot_temp_dir = Path(tempfile.mkdtemp(prefix="rna_forge_eic_"))
        self._plot_counter = 0

        self._build_ui()
        self._load_species_table(self.species_df)
        self._set_status(
            "Using default species table. Upload/generate species, select an MS file, then click 'Load Data & Build EICs'."
        )
        self._show_blank_plot("No MS data loaded")
        self.refresh_tables()

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------
    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QHBoxLayout(central)

        splitter = QSplitter(Qt.Horizontal)
        root.addWidget(splitter)

        left = QWidget()
        left.setMinimumWidth(460)
        left.setMaximumWidth(620)
        left_layout = QVBoxLayout(left)
        splitter.addWidget(left)

        right = QWidget()
        right_layout = QVBoxLayout(right)
        splitter.addWidget(right)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)

        # Input and run controls kept visible at the top.
        input_group = QGroupBox("MS input and analysis")
        input_layout = QGridLayout(input_group)
        self.ms_label = QLabel("No MS file selected")
        self.ms_label.setWordWrap(True)
        choose_ms_btn = QPushButton("Choose .mzML or .mzXML file")
        choose_ms_btn.clicked.connect(self.choose_ms_file)
        load_btn = QPushButton("Load Data && Build EICs")
        load_btn.clicked.connect(self.load_ms_data)

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

        self.ms1_intensity_threshold = QDoubleSpinBox()
        self.ms1_intensity_threshold.setDecimals(1)
        self.ms1_intensity_threshold.setRange(0.0, 1e12)
        self.ms1_intensity_threshold.setSingleStep(1000.0)
        self.ms1_intensity_threshold.setValue(5000.0)

        self.auc_min_intensity = QDoubleSpinBox()
        self.auc_min_intensity.setDecimals(1)
        self.auc_min_intensity.setRange(0.0, 1e12)
        self.auc_min_intensity.setSingleStep(1000.0)
        self.auc_min_intensity.setValue(5000.0)

        self.uncapped_percent_threshold = QDoubleSpinBox()
        self.uncapped_percent_threshold.setDecimals(3)
        self.uncapped_percent_threshold.setRange(0.0, 100.0)
        self.uncapped_percent_threshold.setSingleStep(0.1)
        self.uncapped_percent_threshold.setValue(1.0)

        self.mode_combo = QComboBox()
        self.mode_combo.addItem("MS2-confirmed Only", "filtered")
        self.mode_combo.addItem("All Species", "all")
        self.mode_combo.currentIndexChanged.connect(self.refresh_after_mode_change)

        self.rt_unit_combo = QComboBox()
        self.rt_unit_combo.addItem("Seconds", "seconds")
        self.rt_unit_combo.addItem("Minutes", "minutes")
        self.rt_unit_combo.currentIndexChanged.connect(self.refresh_rt_unit_views)

        input_layout.addWidget(choose_ms_btn, 0, 0, 1, 2)
        input_layout.addWidget(load_btn, 0, 2, 1, 1)
        input_layout.addWidget(self.ms_label, 1, 0, 1, 3)
        input_layout.addWidget(QLabel("Default m/z tolerance (Th)"), 2, 0)
        input_layout.addWidget(self.default_mz_tol, 2, 1)
        input_layout.addWidget(QLabel("Default MS2 window (+/- sec)"), 3, 0)
        input_layout.addWidget(self.default_ms2_window, 3, 1)
        input_layout.addWidget(QLabel("MS1 confirmation threshold"), 4, 0)
        input_layout.addWidget(self.ms1_intensity_threshold, 4, 1)
        input_layout.addWidget(QLabel("AUC point min intensity"), 5, 0)
        input_layout.addWidget(self.auc_min_intensity, 5, 1)
        input_layout.addWidget(QLabel("Uncapped screen threshold (%)"), 6, 0)
        input_layout.addWidget(self.uncapped_percent_threshold, 6, 1)
        input_layout.addWidget(QLabel("Species mode"), 7, 0)
        input_layout.addWidget(self.mode_combo, 7, 1, 1, 2)
        input_layout.addWidget(QLabel("RT display unit"), 8, 0)
        input_layout.addWidget(self.rt_unit_combo, 8, 1, 1, 2)
        left_layout.addWidget(input_group)

        self.left_tabs = QTabWidget()
        left_layout.addWidget(self.left_tabs, stretch=1)
        self._build_species_tab()
        self._build_preset_tab()
        self._build_custom_tab()
        self._build_exports_tab()

        self.status_box = QTextEdit()
        self.status_box.setReadOnly(True)
        self.status_box.setMaximumHeight(120)
        left_layout.addWidget(self.status_box)

        # Right side: plot controls and Plotly view.
        top_right = QHBoxLayout()
        self.selected_species_combo = QComboBox()
        self.selected_species_combo.currentIndexChanged.connect(self.plot_selected_species)
        self.save_window_btn = QPushButton("Save selected RT window")
        self.save_window_btn.clicked.connect(self.save_selected_window)
        top_right.addWidget(QLabel("Selected species"))
        top_right.addWidget(self.selected_species_combo, stretch=1)
        top_right.addWidget(self.save_window_btn)
        right_layout.addLayout(top_right)

        self.plot_hint = QLabel(
            "Use Plotly box/lasso select on the EIC, then click 'Save selected RT window'. "
            "The selected x-range is converted to seconds for storage/AUC. AUC uses only points above the AUC point minimum intensity."
        )
        self.plot_hint.setWordWrap(True)
        right_layout.addWidget(self.plot_hint)

        if WEBENGINE_AVAILABLE:
            self.plot_view = QWebEngineView()
            right_layout.addWidget(self.plot_view, stretch=4)
        else:
            self.plot_view = QTextEdit()
            self.plot_view.setReadOnly(True)
            right_layout.addWidget(self.plot_view, stretch=4)
            self._set_status("Qt WebEngine is unavailable. Plotly box/lasso selection cannot be used in this build.")

        lower_tabs = QTabWidget()
        self.rt_table = QTableWidget()
        self.auc_table_widget = QTableWidget()
        self.intensity_table_widget = QTableWidget()
        lower_tabs.addTab(self.rt_table, "RT windows")
        lower_tabs.addTab(self.auc_table_widget, "AUC table")
        lower_tabs.addTab(self.intensity_table_widget, "Intensity summary")
        right_layout.addWidget(lower_tabs, stretch=2)

        self.efficiency_box = QTextEdit()
        self.efficiency_box.setReadOnly(True)
        self.efficiency_box.setMaximumHeight(105)
        right_layout.addWidget(self.efficiency_box)

    def _build_species_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)
        btn_row = QGridLayout()
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
        btn_row.addWidget(upload_btn, 0, 0)
        btn_row.addWidget(template_btn, 0, 1)
        btn_row.addWidget(reset_btn, 1, 0)
        btn_row.addWidget(add_btn, 1, 1)
        btn_row.addWidget(remove_btn, 2, 0, 1, 2)
        layout.addLayout(btn_row)
        self.species_table = QTableWidget()
        self.species_table.setColumnCount(len(SPECIES_COLUMNS))
        self.species_table.setHorizontalHeaderLabels(SPECIES_COLUMNS)
        self.species_table.itemChanged.connect(self._species_edited)
        layout.addWidget(self.species_table, stretch=1)
        self.left_tabs.addTab(tab, "Species table")

    def _build_preset_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)
        form = QFormLayout()
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
        form.addRow("Preset", self.preset_combo)
        form.addRow("", self.include_capped_z2)
        form.addRow("", self.include_uncapped_z2)
        form.addRow("", self.include_g_ladder)
        form.addRow("", self.g_as_denominator)
        form.addRow("Generated species action", self.generated_mode)
        layout.addLayout(form)
        layout.addWidget(load_preset_btn)
        layout.addStretch(1)
        self.left_tabs.addTab(tab, "Preset generator")

    def _build_custom_tab(self):
        outer = QWidget()
        outer_layout = QVBoxLayout(outer)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        inner = QWidget()
        form = QFormLayout(inner)

        self.custom_cap_name = QLineEdit("CustomCap")
        self.custom_capped_mz = QDoubleSpinBox()
        self.custom_capped_mz.setDecimals(5)
        self.custom_capped_mz.setRange(0, 100000)
        self.custom_capped_mz.setValue(0)
        self.custom_capped_charge = QSpinBox()
        self.custom_capped_charge.setRange(1, 20)
        self.custom_capped_charge.setValue(1)
        self.custom_fragment = QComboBox()
        self.custom_fragment.addItems(["AG", "AUG", "G", "Custom"])
        self.custom_ppp_mz = QDoubleSpinBox()
        self.custom_ppp_mz.setDecimals(5)
        self.custom_ppp_mz.setRange(0, 100000)
        self.custom_ppp_mz.setValue(0)
        self.custom_ppp_charge = QSpinBox()
        self.custom_ppp_charge.setRange(1, 20)
        self.custom_ppp_charge.setValue(1)
        self.custom_include_capped_z2 = QCheckBox("Include calculated capped z = 2 ion")
        self.custom_include_capped_z2.setChecked(True)
        self.custom_include_uncapped_z2 = QCheckBox("Include uncapped z = 2 ladder")
        self.custom_include_g_ladder = QCheckBox("Include G-only ladder")
        self.custom_g_as_denominator = QCheckBox("Include G-only ladder in denominator")
        gen_custom_btn = QPushButton("Generate custom species")
        gen_custom_btn.clicked.connect(self.generate_custom_species)

        form.addRow("Cap label", self.custom_cap_name)
        form.addRow("Capped m/z", self.custom_capped_mz)
        form.addRow("Capped charge z", self.custom_capped_charge)
        form.addRow("5' T1 fragment", self.custom_fragment)
        form.addRow("Custom ppp-fragment m/z", self.custom_ppp_mz)
        form.addRow("Custom ppp-fragment charge z", self.custom_ppp_charge)
        form.addRow("", self.custom_include_capped_z2)
        form.addRow("", self.custom_include_uncapped_z2)
        form.addRow("", self.custom_include_g_ladder)
        form.addRow("", self.custom_g_as_denominator)
        form.addRow("", gen_custom_btn)
        scroll.setWidget(inner)
        outer_layout.addWidget(scroll)
        self.left_tabs.addTab(outer, "Custom generator")

    def _build_exports_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)
        export_rt_btn = QPushButton("Export RT windows CSV")
        export_rt_btn.clicked.connect(self.export_rt_csv)
        export_auc_btn = QPushButton("Export AUC CSV")
        export_auc_btn.clicked.connect(self.export_auc_csv)
        export_intensity_btn = QPushButton("Export intensity summary CSV")
        export_intensity_btn.clicked.connect(self.export_intensity_csv)
        export_excel_btn = QPushButton("Export Excel workbook")
        export_excel_btn.clicked.connect(self.export_excel)
        export_tiff_btn = QPushButton("Export current EIC TIFF")
        export_tiff_btn.clicked.connect(self.export_current_tiff)
        export_zip_btn = QPushButton("Export all EIC TIFFs ZIP")
        export_zip_btn.clicked.connect(self.export_tiff_zip)
        for btn in [export_rt_btn, export_auc_btn, export_intensity_btn, export_excel_btn, export_tiff_btn, export_zip_btn]:
            layout.addWidget(btn)
        layout.addStretch(1)
        self.left_tabs.addTab(tab, "Exports")

    # ------------------------------------------------------------------
    # Generic helpers
    # ------------------------------------------------------------------
    def _set_status(self, text: str):
        if hasattr(self, "status_box"):
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
                table.setItem(r, c, QTableWidgetItem(text))
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
                self.species_df[col] = ""
        self.species_df = self.species_df[SPECIES_COLUMNS]
        self.species_table.blockSignals(True)
        self.species_table.setRowCount(len(self.species_df))
        self.species_table.setColumnCount(len(SPECIES_COLUMNS))
        self.species_table.setHorizontalHeaderLabels(SPECIES_COLUMNS)
        for r in range(len(self.species_df)):
            for c, col in enumerate(SPECIES_COLUMNS):
                val = self.species_df.iloc[r][col]
                text = "" if pd.isna(val) else str(val)
                self.species_table.setItem(r, c, QTableWidgetItem(text))
        self.species_table.blockSignals(False)
        self.species_table.resizeColumnsToContents()
        self.session = None
        self.refresh_species_choices()
        self.refresh_tables()
        self._show_blank_plot("No MS data loaded")

    def _species_edited(self, _item):
        self.session = None
        self.refresh_species_choices()
        self.refresh_tables()
        self._show_blank_plot("Species table edited. Reload MS data.")
        self._set_status("Species table edited. Click 'Load Data & Build EICs' to rebuild EICs.")

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
            if not ok:
                self._message("Species validation failed", msg, QMessageBox.Critical)
                return
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
            self._show_blank_plot("MS file selected. Load data to build EICs.")
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
        self._load_species_table(pd.concat([df, pd.DataFrame([new])], ignore_index=True))
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
        self.species_df = valid

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
                ms1_intensity_threshold=self.ms1_intensity_threshold.value(),
                uncapped_percent_threshold=self.uncapped_percent_threshold.value(),
                progress_cb=progress_cb,
            )
            progress.setValue(100)
            self.session = session
            status_parts = ["MS file loaded and EICs built.", session.ms1_message, session.filter_message, session.uncapped_intensity_message]
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

    def rt_display_unit(self) -> str:
        if not hasattr(self, "rt_unit_combo"):
            return "seconds"
        return self.rt_unit_combo.currentData() or "seconds"

    def rt_axis_label(self) -> str:
        return "Retention Time (minutes)" if self.rt_display_unit() == "minutes" else "Retention Time (seconds)"

    def rt_seconds_to_display(self, value: float) -> float:
        return float(value) / 60.0 if self.rt_display_unit() == "minutes" else float(value)

    def rt_display_to_seconds(self, value: float) -> float:
        return float(value) * 60.0 if self.rt_display_unit() == "minutes" else float(value)

    def refresh_rt_unit_views(self):
        self.plot_selected_species()

    def refresh_after_mode_change(self):
        self.refresh_species_choices()
        self.refresh_tables()
        self.plot_selected_species()

    def refresh_species_choices(self):
        if not hasattr(self, "selected_species_combo"):
            return
        self.selected_species_combo.blockSignals(True)
        self.selected_species_combo.clear()
        if self.session and self.session.loaded:
            state = self.session.state(self.requested_mode())
            for name in state.species_table["name"].astype(str).tolist():
                self.selected_species_combo.addItem(name)
        self.selected_species_combo.blockSignals(False)

    def refresh_tables(self):
        if not hasattr(self, "rt_table"):
            return
        if self.session and self.session.loaded:
            state = self.session.state(self.requested_mode())
            self._df_to_table(self.rt_table, state.rt_windows)
            tbl = auc_table(self.session, self.requested_mode())
            self._df_to_table(self.auc_table_widget, tbl)
            self._df_to_table(self.intensity_table_widget, intensity_table(self.session, self.requested_mode()))
            eff = calculate_efficiency(tbl, honor_use_for_efficiency=True)
            self.efficiency_box.setPlainText(str(eff.get("message", "")))
        else:
            self._df_to_table(self.rt_table, pd.DataFrame({"Message": ["No MS data loaded."]}))
            self._df_to_table(self.auc_table_widget, pd.DataFrame({"Message": ["No AUC table available."]}))
            self._df_to_table(self.intensity_table_widget, pd.DataFrame({"Message": ["No intensity summary available."]}))
            self.efficiency_box.setPlainText("No MS data loaded.")

    # ------------------------------------------------------------------
    # Plotly EIC view and RT-window saving
    # ------------------------------------------------------------------
    def _plotly_html(self, eic: pd.DataFrame, species: str, rt_window=None, rt_unit: str = "seconds") -> str:
        import plotly.graph_objects as go

        fig = go.Figure()
        if eic is not None and len(eic) > 0:
            fig.add_trace(
                go.Scatter(
                    x=(pd.to_numeric(eic["RT"], errors="coerce") / 60.0) if rt_unit == "minutes" else pd.to_numeric(eic["RT"], errors="coerce"),
                    y=pd.to_numeric(eic["Intensity"], errors="coerce"),
                    mode="lines+markers",
                    name=species,
                    marker={"size": 5},
                    line={"width": 1.5},
                )
            )
        if rt_window is not None:
            lo, hi = sorted([float(rt_window[0]), float(rt_window[1])])
            if rt_unit == "minutes":
                lo, hi = lo / 60.0, hi / 60.0
            fig.add_vrect(x0=lo, x1=hi, fillcolor="LightSkyBlue", opacity=0.25, line_width=0)
        fig.update_layout(
            template="plotly_white",
            title=f"EIC for {species}",
            xaxis_title="Retention Time (minutes)" if rt_unit == "minutes" else "Retention Time (seconds)",
            yaxis_title="Intensity",
            dragmode="select",
            hovermode="closest",
            margin={"l": 70, "r": 25, "t": 55, "b": 60},
            modebar_add=["select2d", "lasso2d"],
        )
        html = fig.to_html(include_plotlyjs=True, full_html=True, config={"displaylogo": False, "displayModeBar": True})
        js = r"""
<script>
window.__selectedRTWindow = null;
function attachSelectionHandler(){
  const gd = document.querySelector('.plotly-graph-div');
  if(!gd || !gd.on){ setTimeout(attachSelectionHandler, 100); return; }
  gd.on('plotly_selected', function(eventData){
    if(eventData && eventData.points && eventData.points.length >= 2){
      const xs = eventData.points.map(p => Number(p.x)).filter(x => Number.isFinite(x));
      if(xs.length >= 2){
        window.__selectedRTWindow = [Math.min.apply(null, xs), Math.max.apply(null, xs)];
      }
    }
  });
  gd.on('plotly_deselect', function(){ window.__selectedRTWindow = null; });
}
attachSelectionHandler();
</script>
"""
        return html.replace("</body>", js + "\n</body>")


    def _load_html_in_plot_view(self, html_text: str):
        """Load HTML through a local file instead of setHtml().

        QWebEngineView.setHtml() internally uses a data URL in many Qt builds.
        Large inline Plotly HTML can exceed that path and silently fail, leaving the
        old blank page visible. A temporary local file avoids that limitation and is
        more reliable after PyInstaller packaging.
        """
        if not WEBENGINE_AVAILABLE:
            return
        self._plot_counter += 1
        path = self._plot_temp_dir / f"eic_plot_{self._plot_counter:04d}.html"
        path.write_text(html_text, encoding="utf-8")
        self.plot_view.load(QUrl.fromLocalFile(str(path)))

    def _show_blank_plot(self, message: str):
        if not hasattr(self, "plot_view"):
            return
        if WEBENGINE_AVAILABLE:
            escaped_message = html_lib.escape(str(message))
            html = f"""
<!doctype html>
<html><body style="font-family: Arial, sans-serif; background: white; color: #222;">
<div style="height: 80vh; display: flex; align-items: center; justify-content: center; border: 1px solid #ddd;">
  <div style="text-align: center; font-size: 20px;">{escaped_message}</div>
</div>
<script>window.__selectedRTWindow = null;</script>
</body></html>
"""
            self._load_html_in_plot_view(html)
        else:
            self.plot_view.setPlainText(message)

    def plot_selected_species(self):
        if not hasattr(self, "plot_view"):
            return
        if not (self.session and self.session.loaded) or self.selected_species_combo.count() == 0:
            self._show_blank_plot("No MS data loaded")
            return
        species = self.selected_species_combo.currentText()
        state = self.session.state(self.requested_mode())
        eic = state.eic_list.get(species, pd.DataFrame())
        rt_window = None
        if state.rt_windows is not None and len(state.rt_windows) > 0:
            row = state.rt_windows[state.rt_windows["Species"].astype(str) == species]
            if len(row) == 1 and pd.notna(row.iloc[0]["RTmin"]) and pd.notna(row.iloc[0]["RTmax"]):
                rt_window = (float(row.iloc[0]["RTmin"]), float(row.iloc[0]["RTmax"]))
        if WEBENGINE_AVAILABLE:
            try:
                self._load_html_in_plot_view(self._plotly_html(eic, species, rt_window=rt_window, rt_unit=self.rt_display_unit()))
            except Exception as exc:
                self._show_blank_plot("EIC plot failed to render")
                self._set_status("EIC plot failed to render:\n" + str(exc))
        else:
            self.plot_view.setPlainText("Qt WebEngine is unavailable. Rebuild with PySide6 QtWebEngine support.")

    def save_selected_window(self):
        if not (self.session and self.session.loaded):
            self._message("No data", "Load MS data first.", QMessageBox.Warning)
            return
        if not WEBENGINE_AVAILABLE:
            self._message("Plot selection unavailable", "Qt WebEngine is unavailable in this build.", QMessageBox.Warning)
            return
        self.plot_view.page().runJavaScript("JSON.stringify(window.__selectedRTWindow || null)", self._save_selected_window_from_js)

    def _save_selected_window_from_js(self, result):
        try:
            if result in (None, "null", ""):
                self._message(
                    "No RT window selected",
                    "Use the Plotly box or lasso select tool on the EIC, then click 'Save selected RT window'.",
                    QMessageBox.Warning,
                )
                return
            selected = json.loads(result)
            if not isinstance(selected, list) or len(selected) != 2:
                raise ValueError("Invalid selection returned from Plotly.")
            x_min, x_max = float(selected[0]), float(selected[1])
            if x_min == x_max:
                raise ValueError("Selection has zero width.")
            rt_min, rt_max = self.rt_display_to_seconds(x_min), self.rt_display_to_seconds(x_max)
            species = self.selected_species_combo.currentText()
            min_i = self.auc_min_intensity.value()
            auc = save_rt_window(self.session, self.requested_mode(), species, rt_min, rt_max, min_intensity=min_i)
            unit_label = "minutes" if self.rt_display_unit() == "minutes" else "seconds"
            self._set_status(
                f"Window for {species} saved: [{x_min:.3f}, {x_max:.3f}] {unit_label} "
                f"([{min(rt_min, rt_max):.2f}, {max(rt_min, rt_max):.2f}] seconds); "
                f"AUC = {auc:.4g}; included points Intensity > {min_i:.4g}"
            )
            self.refresh_tables()
            self.plot_selected_species()
        except Exception as exc:
            self._message("Could not save RT window", str(exc), QMessageBox.Critical)

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

    def export_intensity_csv(self):
        if not self._require_session():
            return
        mode = self.session.current_mode(self.requested_mode())
        stem = Path(self.ms_file).stem if self.ms_file else "intensity_summary"
        path, _ = QFileDialog.getSaveFileName(self, "Export intensity summary", f"{stem}_{mode}_intensity_summary.csv", "CSV files (*.csv)")
        if path:
            export_csv(intensity_table(self.session, self.requested_mode()), path)
            self._set_status(f"Intensity summary exported: {path}")

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
        path, _ = QFileDialog.getSaveFileName(self, "Export current EIC TIFF", f"{safe_filename(species)}_EIC.tiff", "TIFF files (*.tiff *.tif)")
        if path:
            export_current_eic_tiff(self.session, self.requested_mode(), species, path, rt_unit=self.rt_display_unit())
            self._set_status(f"Current EIC TIFF exported: {path}")

    def export_tiff_zip(self):
        if not self._require_session():
            return
        mode = self.session.current_mode(self.requested_mode())
        stem = Path(self.ms_file).stem if self.ms_file else "eic_plots"
        path, _ = QFileDialog.getSaveFileName(self, "Export all EIC TIFFs ZIP", f"{stem}_{mode}_EIC_plots.zip", "ZIP files (*.zip)")
        if path:
            export_all_eic_tiffs_zip(self.session, self.requested_mode(), path, rt_unit=self.rt_display_unit())
            self._set_status(f"EIC TIFF ZIP exported: {path}")


def apply_light_palette(app: QApplication):
    app.setStyle("Fusion")
    palette = QPalette()
    palette.setColor(QPalette.Window, QColor(245, 245, 245))
    palette.setColor(QPalette.WindowText, QColor(0, 0, 0))
    palette.setColor(QPalette.Base, QColor(255, 255, 255))
    palette.setColor(QPalette.AlternateBase, QColor(240, 240, 240))
    palette.setColor(QPalette.ToolTipBase, QColor(255, 255, 220))
    palette.setColor(QPalette.ToolTipText, QColor(0, 0, 0))
    palette.setColor(QPalette.Text, QColor(0, 0, 0))
    palette.setColor(QPalette.Button, QColor(245, 245, 245))
    palette.setColor(QPalette.ButtonText, QColor(0, 0, 0))
    palette.setColor(QPalette.BrightText, QColor(255, 0, 0))
    palette.setColor(QPalette.Highlight, QColor(40, 110, 180))
    palette.setColor(QPalette.HighlightedText, QColor(255, 255, 255))
    app.setPalette(palette)


def main():
    app = QApplication(sys.argv)
    apply_light_palette(app)
    win = CappingEfficiencyWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
