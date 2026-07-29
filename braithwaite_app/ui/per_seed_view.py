"""Passo 6 -- vista por-semente: quatro gráficos lidos dos plotfiles
reais de UMA run (não do resumo persistido -- esse é um único ponto por
semente; aqui é a trajetória completa). Janela de validade sempre
sombreada; fora dela, esmaecido/tracejado com rótulo "estrela de fundo
derivando" -- nunca com o mesmo peso visual da janela válida.
"""

from pathlib import Path

from PyQt6.QtWidgets import QComboBox, QLabel, QVBoxLayout, QWidget
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from core.extraction import extract_full_series
from core.star_builder import parse_rho_c_log

DIVB_HEALTHY_MAX = 1e-8  # order of magnitude above the ~1e-10 baseline seen all session


class PerSeedView(QWidget):
    def __init__(self, run_dir: Path, parent=None):
        super().__init__(parent)
        self.run_dir = Path(run_dir)
        self._runs: dict[str, dict] = {}  # run_id -> context (log_path, t_dyn_s, w_abs_erg, window, rho_c_ic)

        layout = QVBoxLayout(self)
        self.selector = QComboBox()
        self.selector.currentTextChanged.connect(self._on_selected)
        layout.addWidget(self.selector)

        self.status_label = QLabel("Nenhuma run carregada.")
        layout.addWidget(self.status_label)

        self.figure = Figure(figsize=(9, 7))
        self.canvas = FigureCanvasQTAgg(self.figure)
        layout.addWidget(self.canvas)
        self._axes = self.figure.subplots(2, 2)

    def register_run(self, run_id: str, log_path: Path, t_dyn_s: float,
                      w_abs_erg: float, window: tuple[float, float], rho_c_ic: float):
        self._runs[run_id] = {
            "log_path": log_path, "t_dyn_s": t_dyn_s, "w_abs_erg": w_abs_erg,
            "window": window, "rho_c_ic": rho_c_ic,
        }
        if self.selector.findText(run_id) < 0:
            self.selector.addItem(run_id)

    def _on_selected(self, run_id: str):
        if not run_id or run_id not in self._runs:
            return
        ctx = self._runs[run_id]
        try:
            self._plot(run_id, ctx)
            self.status_label.setText(f"'{run_id}' -- janela válida [{ctx['window'][0]:.3f}, {ctx['window'][1]:.3f}] t/t_dyn")
        except Exception as e:  # pragma: no cover -- surfaced to the user, not swallowed
            self.status_label.setText(f"Erro lendo plotfiles de '{run_id}': {e}")

    def _plot(self, run_id: str, ctx: dict):
        lo, hi = ctx["window"]
        rho_series = parse_rho_c_log(ctx["log_path"], ctx["t_dyn_s"])
        field_series = extract_full_series(self.run_dir, run_id, ctx["w_abs_erg"], ctx["t_dyn_s"])

        for ax in self._axes.flat:
            ax.clear()

        ax_ratio, ax_vs_rho, ax_rho, ax_divb = self._axes.flat

        # 1) E_tor/E_mag(t) -- valid window shaded, outside faded/dashed
        if field_series:
            t = [r["t_ttdyn"] for r in field_series]
            ratio = [r["E_tor_over_Emag"] for r in field_series]
            in_window = [lo <= ti <= hi for ti in t]
            self._plot_with_window(ax_ratio, t, ratio, in_window, lo, hi)
        ax_ratio.set_xlabel("t/t_dyn")
        ax_ratio.set_ylabel("E_tor/E_mag")
        ax_ratio.set_title("E_tor/E_mag(t)")

        # 2) E_tor/E_mag vs rho_c -- the validity-boundary evidence panel
        if field_series and rho_series:
            rho_at = _nearest_lookup(rho_series)
            pct = []
            ratio2 = []
            in_window2 = []
            for r in field_series:
                rho = rho_at(r["t_ttdyn"])  # both series are in t/t_dyn units
                if rho is None:
                    continue
                pct.append(100.0 * (rho - ctx["rho_c_ic"]) / ctx["rho_c_ic"])
                ratio2.append(r["E_tor_over_Emag"])
                in_window2.append(lo <= r["t_ttdyn"] <= hi)
            self._plot_with_window(ax_vs_rho, pct, ratio2, in_window2, None, None, x_is_time=False)
        ax_vs_rho.set_xlabel("rho_c desvio do IC (%)")
        ax_vs_rho.set_ylabel("E_tor/E_mag")
        ax_vs_rho.set_title("E_tor/E_mag vs rho_c (delimita a validade)")

        # 3) rho_c(t) with window marked
        if rho_series:
            t = [r[0] for r in rho_series]
            rho = [r[1] for r in rho_series]
            ax_rho.plot(t, rho, color="tab:blue", linewidth=1)
            ax_rho.axvspan(lo, hi, alpha=0.15, color="green")
        ax_rho.set_xlabel("t/t_dyn")
        ax_rho.set_ylabel("rho_c (g/cm3)")
        ax_rho.set_title("rho_c(t) -- janela sombreada")

        # 4) Div_B interior(t) -- health monitor, red if it leaves ~1e-10
        if field_series:
            t = [r["t_ttdyn"] for r in field_series]
            divb = [r["divB_interior_max"] for r in field_series]
            colors = ["red" if d > DIVB_HEALTHY_MAX else "tab:green" for d in divb]
            ax_divb.scatter(t, divb, c=colors, s=14)
            ax_divb.plot(t, divb, color="tab:green", linewidth=0.5, alpha=0.5)
            ax_divb.axhspan(0, DIVB_HEALTHY_MAX, alpha=0.08, color="green")
        ax_divb.set_yscale("log")
        ax_divb.set_xlabel("t/t_dyn")
        ax_divb.set_ylabel("|Div B| interior (max)")
        ax_divb.set_title("Div B interior -- mascarado, nunca a borda")

        self.figure.tight_layout()
        self.canvas.draw()

    @staticmethod
    def _plot_with_window(ax, x, y, in_window, lo, hi, x_is_time=True):
        if not x:
            return
        x_in = [xi for xi, w in zip(x, in_window) if w]
        y_in = [yi for yi, w in zip(y, in_window) if w]
        x_out = [xi for xi, w in zip(x, in_window) if not w]
        y_out = [yi for yi, w in zip(y, in_window) if not w]
        if x_out:
            ax.plot(x_out, y_out, "--", color="0.7", linewidth=1,
                    label="fora de validade -- estrela de fundo derivando")
        if x_in:
            ax.plot(x_in, y_in, "-", color="tab:blue", linewidth=1.5, label="janela válida")
        if x_is_time and lo is not None:
            ax.axvspan(lo, hi, alpha=0.12, color="green")
        ax.legend(fontsize=7, loc="best")


def _nearest_lookup(rho_series):
    """Closure over a sorted (t/t_dyn, rho_c) series -> nearest-rho_c(t_ttdyn) lookup."""
    def lookup(t_ttdyn):
        if not rho_series:
            return None
        best = min(rho_series, key=lambda r: abs(r[0] - t_ttdyn))
        return best[1]
    return lookup
