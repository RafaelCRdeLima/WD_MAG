"""Run queue widget: one row per run (star construction or a field
seed), live status read from its log file by a periodic QTimer -- never
a blocking wait. Generic over run_id (string) so the same widget serves
both Passo A's single star-construction run and Passo B's N seed runs.
"""

from PyQt6.QtWidgets import (
    QHeaderView,
    QLabel,
    QProgressBar,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

COLUMNS = ["run", "resolução", "status", "t/t_dyn", "progresso"]


class RunQueue(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)

        self.empty_label = QLabel("Nenhuma run em andamento.")
        layout.addWidget(self.empty_label)

        self.table = QTableWidget(0, len(COLUMNS))
        self.table.setHorizontalHeaderLabels(COLUMNS)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.setVisible(False)
        layout.addWidget(self.table)

        self._row_by_id = {}

    def add_run(self, run_id: str, resolution: int, target_ttdyn: float):
        self.empty_label.setVisible(False)
        self.table.setVisible(True)

        row = self.table.rowCount()
        self.table.insertRow(row)
        self.table.setItem(row, 0, QTableWidgetItem(str(run_id)))
        self.table.setItem(row, 1, QTableWidgetItem(f"{resolution}³"))
        self.table.setItem(row, 2, QTableWidgetItem("na fila"))
        self.table.setItem(row, 3, QTableWidgetItem("0.000"))

        bar = QProgressBar()
        bar.setRange(0, 1000)
        bar.setValue(0)
        self.table.setCellWidget(row, 4, bar)

        self._row_by_id[run_id] = {"row": row, "target_ttdyn": target_ttdyn, "bar": bar}

    def update_run(self, run_id: str, status: str, ttdyn: float):
        info = self._row_by_id.get(run_id)
        if info is None:
            return
        row = info["row"]
        self.table.setItem(row, 2, QTableWidgetItem(status))
        self.table.setItem(row, 3, QTableWidgetItem(f"{ttdyn:.3f}"))
        target = info["target_ttdyn"] or 1.0
        frac = max(0.0, min(1.0, ttdyn / target))
        info["bar"].setValue(int(frac * 1000))
