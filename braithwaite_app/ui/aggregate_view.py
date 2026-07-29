"""Passo 8 -- vista agregada: a figura do resultado. Lê a fonte de
verdade persistida (core/persistence.py), nunca recolapsa em média +
barra de erro -- o rank-ordering por semente É o resultado (prova que a
espalha é real, não fase de uma oscilação comum; ver docs/teoria.md e a
revisão desta sessão). Sobreposição de resolução e exportação para
publicação incluídas.
"""

from pathlib import Path

from PyQt6.QtWidgets import QComboBox, QFileDialog, QHBoxLayout, QLabel, QPushButton, QVBoxLayout, QWidget
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from core import persistence
from core.extraction import extract_full_series

POLOIDAL_TOROIDAL_THRESHOLD = 0.5


class AggregateView(QWidget):
    def __init__(self, run_dir: Path, parent=None):
        super().__init__(parent)
        self.run_dir = Path(run_dir)

        layout = QVBoxLayout(self)

        top_row = QHBoxLayout()
        self.star_selector = QComboBox()
        self.star_selector.currentTextChanged.connect(self.refresh)
        refresh_btn = QPushButton("recarregar do disco")
        refresh_btn.clicked.connect(self.refresh)
        export_btn = QPushButton("exportar figura (PDF)")
        export_btn.clicked.connect(self._export)
        top_row.addWidget(QLabel("estrela:"))
        top_row.addWidget(self.star_selector, stretch=1)
        top_row.addWidget(refresh_btn)
        top_row.addWidget(export_btn)
        layout.addLayout(top_row)

        self.status_label = QLabel("")
        layout.addWidget(self.status_label)

        self.figure = Figure(figsize=(10, 5))
        self.canvas = FigureCanvasQTAgg(self.figure)
        layout.addWidget(self.canvas)
        self._ax_ranked, self._ax_overlay = self.figure.subplots(1, 2)

        self.reload_star_list()

    def reload_star_list(self):
        current = self.star_selector.currentText()
        df = persistence.load_all_results()
        stars = df[df["row_type"] == "star"]["star_cache_key"].dropna().unique().tolist()
        self.star_selector.blockSignals(True)
        self.star_selector.clear()
        for s in stars:
            self.star_selector.addItem(s)
        if current and self.star_selector.findText(current) >= 0:
            self.star_selector.setCurrentText(current)
        self.star_selector.blockSignals(False)
        self.refresh()

    def refresh(self):
        star_key = self.star_selector.currentText()
        for ax in (self._ax_ranked, self._ax_overlay):
            ax.clear()

        if not star_key:
            self.status_label.setText("Nenhuma estrela com resultados persistidos ainda.")
            self.canvas.draw()
            return

        df = persistence.load_all_results()
        star_row = df[(df["row_type"] == "star") & (df["star_cache_key"] == star_key)]
        seed_rows = df[(df["row_type"] == "seed") & (df["star_cache_key"] == star_key)]
        seed_rows = seed_rows.dropna(subset=["E_tor_over_Emag", "seed"])

        if seed_rows.empty:
            self.status_label.setText(f"Estrela '{star_key}': ainda sem sementes medidas.")
            self.canvas.draw()
            return

        self._plot_ranked(seed_rows)
        self._plot_overlay(star_row, seed_rows)

        n = len(seed_rows)
        vmin, vmax = seed_rows["E_tor_over_Emag"].min(), seed_rows["E_tor_over_Emag"].max()
        self.status_label.setText(
            f"Estrela '{star_key}': {n} sementes, E_tor/E_mag em [{vmin:.3f}, {vmax:.3f}] "
            f"-- todas poloidal-dominadas: {bool((seed_rows['E_tor_over_Emag'] < POLOIDAL_TOROIDAL_THRESHOLD).all())}"
        )
        self.figure.tight_layout()
        self.canvas.draw()

    def _plot_ranked(self, seed_rows):
        ax = self._ax_ranked
        ordered = seed_rows.sort_values("E_tor_over_Emag")
        x = range(len(ordered))
        y = ordered["E_tor_over_Emag"].tolist()
        labels = [str(int(s)) for s in ordered["seed"].tolist()]

        res_colors = {64: "tab:blue", 128: "tab:orange"}
        colors = [res_colors.get(int(r), "gray") for r in ordered["resolution"].tolist()]

        ax.scatter(x, y, c=colors, s=60, zorder=3)
        for xi, yi, seed_label in zip(x, y, labels):
            ax.annotate(seed_label, (xi, yi), textcoords="offset points", xytext=(0, 8), fontsize=8, ha="center")

        if len(y) > 1:
            ax.axhspan(min(y), max(y), alpha=0.08, color="tab:blue")
        ax.axhline(POLOIDAL_TOROIDAL_THRESHOLD, color="red", linestyle="--", linewidth=1,
                   label="limiar poloidal/toroidal (0.5)")
        ax.set_ylim(0, max(POLOIDAL_TOROIDAL_THRESHOLD * 1.2, max(y) * 1.2 if y else 1))
        ax.set_xticks([])
        ax.set_ylabel("E_tor/E_mag (janela válida)")
        ax.set_title("Por semente, ordenado -- NÃO é média+erro")
        ax.legend(fontsize=8, loc="upper left")

    def _plot_overlay(self, star_row, seed_rows):
        ax = self._ax_overlay
        if star_row.empty:
            return
        star = star_row.iloc[-1]
        if not bool(star.get("window_valid", False)):
            return
        lo, hi = float(star["window_lo"]), float(star["window_hi"])
        t_dyn_s_guess = None  # not persisted directly; recovered via log if available below

        for _, seed_row in seed_rows.iterrows():
            run_id_guess = Path(str(seed_row.get("plotfile_path", ""))).name
            # plotfile name is plt_<run_id><step>; recover run_id by stripping
            # the trailing zero-padded step and the "plt_" prefix
            if not run_id_guess.startswith("plt_"):
                continue
            run_id = run_id_guess[len("plt_"):-5]
            log_path = self.run_dir / f"run_{run_id}.log"
            if not log_path.exists():
                continue
            from core.run_launcher import read_progress
            progress = read_progress(log_path, run_id)
            if progress.t_dyn_s is None or progress.w_abs_erg is None:
                continue
            series = extract_full_series(self.run_dir, run_id, progress.w_abs_erg, progress.t_dyn_s)
            series = [s for s in series if lo <= s["t_ttdyn"] <= hi]
            if not series:
                continue
            t = [s["t_ttdyn"] for s in series]
            ratio = [s["E_tor_over_Emag"] for s in series]
            ax.plot(t, ratio, marker="o", markersize=3, label=f"seed {int(seed_row['seed'])}")

        ax.axvspan(lo, hi, alpha=0.08, color="green")
        ax.set_xlabel("t/t_dyn")
        ax.set_ylabel("E_tor/E_mag")
        ax.set_title("Sobrepostas na janela válida -- não cruzam = espalha real")
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=7, loc="best")

    def _export(self):
        path, _ = QFileDialog.getSaveFileName(self, "Exportar figura", "braithwaite_result.pdf", "PDF (*.pdf)")
        if path:
            self.figure.savefig(path, bbox_inches="tight")
            self.status_label.setText(f"Exportado para {path}")
