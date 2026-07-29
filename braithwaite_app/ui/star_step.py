"""Passo A: estrela de fundo. Não pré-relaxa -- constrói a IC determinística,
verifica o patch de gravidade, roda uma evolução curta (Fase 3 liga isto a
uma run real; aqui é só a estrutura) e mostra a janela de medição como
resultado. Nunca exibe um booleano de "estabilizou".
"""

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import (
    QDoubleSpinBox,
    QFormLayout,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from core.scf_store import neutronization_check
from ui.widgets import ScientificDoubleSpinBox


class StarStepWidget(QWidget):
    """Emits `star_ready(dict)` with the star's cache key and measurement
    window once construction finishes and the window is valid. Emits
    `star_invalidated()` whenever the current selection stops being usable
    (new rho_c typed, gate failed, etc.) so Passo B can re-disable itself.
    Emits `build_requested(dict)` when the user clicks the build button --
    MainWindow owns the actual launch + polling (same pattern as
    FieldStepWidget.run_requested) and calls back into
    `on_run_started`/`on_result` below.
    """

    star_ready = pyqtSignal(dict)
    star_invalidated = pyqtSignal()
    build_requested = pyqtSignal(dict)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._last_result = None  # set once construction succeeds (Fase 3+)
        self._build_ui()

    def _build_ui(self):
        root = QVBoxLayout(self)
        box = QGroupBox("Passo A -- Estrela de fundo")
        layout = QVBoxLayout(box)

        form = QFormLayout()
        self.rho_c_spin = ScientificDoubleSpinBox(decimals=3)
        self.rho_c_spin.setRange(1e5, 1e13)
        self.rho_c_spin.setValue(1.0e9)  # default: the star behind the 10-seed study
        form.addRow("rho_c (g/cm3):", self.rho_c_spin)

        self.mu_e_spin = QDoubleSpinBox()
        self.mu_e_spin.setDecimals(2)
        self.mu_e_spin.setRange(1.5, 2.5)
        self.mu_e_spin.setValue(2.0)
        form.addRow("mu_e:", self.mu_e_spin)

        layout.addLayout(form)

        scope_note = QLabel(
            "rho_c e mu_e são os DOIS únicos graus de liberdade físicos desta "
            "estrela, não uma limitação da tela: o método de Braithwaite exige "
            "uma estrela de fundo sem campo e sem rotação por construção (o campo "
            "é semeado e liberado depois, não imposto aqui -- ver o texto no topo "
            "da janela); e a EOS usada (ztwd, degenerada a T=0) não depende de "
            "temperatura, só de densidade e composição (mu_e). Fixados "
            "campo-free + sem rotação + T=0, a estrutura inteira (massa, raio, "
            "perfil) fica determinada por (rho_c, mu_e) -- não há um terceiro "
            "parâmetro físico para expor."
        )
        scope_note.setWordWrap(True)
        scope_note.setStyleSheet("color: #555; font-style: italic;")
        layout.addWidget(scope_note)

        self.rho_c_spin.valueChanged.connect(self._on_params_changed)
        self.mu_e_spin.valueChanged.connect(self._on_params_changed)

        self.neutronization_warning = QLabel("")
        self.neutronization_warning.setWordWrap(True)
        self.neutronization_warning.setStyleSheet("color: #b00020; font-weight: bold;")
        layout.addWidget(self.neutronization_warning)

        row = QHBoxLayout()
        self.build_button = QPushButton("Construir / verificar estrela")
        self.build_button.clicked.connect(self._on_build_clicked)
        row.addWidget(self.build_button)
        row.addStretch(1)
        layout.addLayout(row)

        grid = QGridLayout()
        self.result_labels = {}
        for i, key in enumerate(["rho_c medido", "VE", "janela de medição"]):
            grid.addWidget(QLabel(key), 0, i)
            val = QLabel("—")
            self.result_labels[key] = val
            grid.addWidget(val, 1, i)
        layout.addLayout(grid)

        self.status_label = QLabel("")
        self.status_label.setWordWrap(True)
        layout.addWidget(self.status_label)

        root.addWidget(box)
        self._on_params_changed()

    def _on_params_changed(self):
        self._last_result = None
        self.star_invalidated.emit()
        self.result_labels["rho_c medido"].setText("—")
        self.result_labels["VE"].setText("—")
        self.result_labels["janela de medição"].setText("—")
        self.status_label.setText("")

        rho_c = self.rho_c_spin.value()
        mu_e = self.mu_e_spin.value()
        check = neutronization_check(rho_c, mu_e)
        if check["ok"]:
            self.neutronization_warning.setText("")
            self.build_button.setEnabled(True)
        else:
            self.neutronization_warning.setText(
                f"BLOQUEADO: rho_c = {rho_c:.3e} g/cm³ >= limiar de neutronização "
                f"{check['threshold']:.3e} g/cm³ para mu_e={mu_e:.2f} (Boshkayev et al. "
                "2013, ApJ 762, 117). Não é possível construir esta estrela."
            )
            self.build_button.setEnabled(False)

    def current_params(self) -> dict:
        return {"rho_c": self.rho_c_spin.value(), "mu_e": self.mu_e_spin.value()}

    def _on_build_clicked(self):
        self.build_button.setEnabled(False)
        self.status_label.setText("Lançando construção da estrela em segundo plano...")
        self.build_requested.emit(self.current_params())

    def on_run_started(self, run_id: str):
        self.status_label.setText(f"Construindo estrela (run '{run_id}')... a interface segue responsiva.")

    def on_result(self, result: dict):
        """MainWindow calls this once the star-construction run ends and
        find_measurement_window() has been evaluated on its log. `result`
        must include: rho_c_measured, VE, window_result (the dict returned
        by find_measurement_window), cache_key, resolution.
        """
        self.build_button.setEnabled(True)
        self.result_labels["rho_c medido"].setText(f"{result['rho_c_measured']:.4e}")
        self.result_labels["VE"].setText(f"{result.get('VE', float('nan')):.3e}")

        window_result = result["window_result"]
        if window_result["valid"]:
            lo, hi = window_result["window"]
            self.result_labels["janela de medição"].setText(f"[{lo:.3f}, {hi:.3f}] t/t_dyn")
            self.status_label.setText("Estrela pronta -- janela de medição válida.")
            self.star_ready.emit({
                "rho_c": result["rho_c_measured"],
                "cache_key": result["cache_key"],
                "resolution": result["resolution"],
                "window": window_result["window"],
                "t_dyn_s": result["t_dyn_s"],
            })
        else:
            self.result_labels["janela de medição"].setText("sem janela válida")
            self.status_label.setText(window_result["reason"])
            self.star_invalidated.emit()

    def on_run_failed(self, message: str):
        self.build_button.setEnabled(True)
        self.status_label.setText(f"Falha na construção da estrela: {message}")
        self.star_invalidated.emit()
