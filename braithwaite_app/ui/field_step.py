"""Passo B: campo. Desabilitado até o Passo A produzir uma estrela com
janela de medição válida (valid=True) -- sem opção de forçar, mesma
política do portão de neutronização.
"""

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import (
    QCheckBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QRadioButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)


class FieldStepWidget(QWidget):
    run_requested = pyqtSignal(dict)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._star_info = None  # dict from StarStepWidget.star_ready, or None
        self._build_ui()
        self.set_star(None)

    def _build_ui(self):
        root = QVBoxLayout(self)
        self.box = QGroupBox("Passo B -- Campo (precisa de uma estrela válida no Passo A)")
        layout = QVBoxLayout(self.box)

        field_row = QHBoxLayout()
        self.range_checkbox = QCheckBox("faixa (em vez de valor único)")
        self.range_checkbox.toggled.connect(self._on_range_toggled)
        self.emag_min = QDoubleSpinBox()
        self.emag_min.setDecimals(3)
        self.emag_min.setRange(0.001, 1.0)
        self.emag_min.setSingleStep(0.01)
        self.emag_min.setValue(0.15)
        self.emag_max = QDoubleSpinBox()
        self.emag_max.setDecimals(3)
        self.emag_max.setRange(0.001, 1.0)
        self.emag_max.setSingleStep(0.01)
        self.emag_max.setValue(0.15)
        self.emag_max.setEnabled(False)
        field_row.addWidget(QLabel("E_mag/|W| valor/mínimo:"))
        field_row.addWidget(self.emag_min)
        field_row.addWidget(self.range_checkbox)
        field_row.addWidget(QLabel("máximo:"))
        field_row.addWidget(self.emag_max)
        field_row.addStretch(1)
        layout.addLayout(field_row)

        form = QFormLayout()
        self.n_seeds_spin = QSpinBox()
        self.n_seeds_spin.setRange(1, 100)
        self.n_seeds_spin.setValue(10)
        form.addRow("Número de sementes:", self.n_seeds_spin)

        res_row = QHBoxLayout()
        self.res_64 = QRadioButton("64³")
        self.res_128 = QRadioButton("128³")
        self.res_64.setChecked(True)
        res_row.addWidget(self.res_64)
        res_row.addWidget(self.res_128)
        res_row.addStretch(1)
        form.addRow("Resolução:", res_row)
        layout.addLayout(form)

        run_row = QHBoxLayout()
        self.run_button = QPushButton("Rodar estudo")
        self.run_button.clicked.connect(self._on_run_clicked)
        run_row.addStretch(1)
        run_row.addWidget(self.run_button)
        layout.addLayout(run_row)

        self.status_label = QLabel("")
        self.status_label.setWordWrap(True)
        layout.addWidget(self.status_label)

        root.addWidget(self.box)

    def _on_range_toggled(self, checked):
        self.emag_max.setEnabled(checked)
        if not checked:
            self.emag_max.setValue(self.emag_min.value())

    def set_star(self, star_info: dict | None):
        """Called by MainWindow on star_ready / star_invalidated. `star_info`
        is None (no valid star) or the dict emitted by StarStepWidget
        (must include at least: cache_key, resolution, window)."""
        self._star_info = star_info
        valid = star_info is not None
        self.box.setEnabled(valid)
        if valid:
            self.box.setTitle(
                f"Passo B -- Campo (estrela: rho_c={star_info.get('rho_c', float('nan')):.3e}, "
                f"janela válida [{star_info['window'][0]:.2f}, {star_info['window'][1]:.2f}] t/t_dyn)"
            )
            # keep the field study's resolution matched to the star's, since
            # the IC is cached per-resolution -- running field at a different
            # resolution than the star was built for is not meaningful.
            if star_info.get("resolution") == 128:
                self.res_128.setChecked(True)
            else:
                self.res_64.setChecked(True)
        else:
            self.box.setTitle("Passo B -- Campo (precisa de uma estrela válida no Passo A)")

    def _on_run_clicked(self):
        if self._star_info is None:
            self.status_label.setText("Bloqueado: nenhuma estrela válida do Passo A.")
            return
        config = {
            "star": self._star_info,
            "emag_over_w_min": self.emag_min.value(),
            "emag_over_w_max": self.emag_max.value() if self.range_checkbox.isChecked() else self.emag_min.value(),
            "n_seeds": self.n_seeds_spin.value(),
            "resolution": 128 if self.res_128.isChecked() else 64,
        }
        self.run_requested.emit(config)
