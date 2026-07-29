"""Braithwaite desktop app -- replaces the Streamlit Tab 5 skeleton.

Scope (do not expand without re-reading the spec this was built from):
this app owns ONLY the Braithwaite stability study (background star ->
seed field -> Castro relaxation -> Bt/Bp diagnostics). Tabs 1-4 of the
Streamlit dashboard (equilibrium, sweep, export, run registry) are
untouched and continue to be the interface for the SCF side of the
project; this app reuses their persistence layer (dashboard/store.py)
but does not reimplement or replace them.

Two sequential physics steps, not independent controls:
  Passo A (ui/star_step.py) -- background star. Does NOT pre-relax to a
    settled state (tested and refuted this session -- rho_c never
    settles under any IC/damping combination tried, because Castro's MHD
    reconstruction has no well-balancing). Builds a deterministic IC,
    verifies the gravity r=0 patch, runs a short evolution, and reports
    a MEASUREMENT WINDOW [t_field_relax, X] (or "no valid window") --
    never a settled/stabilized boolean.
  Passo B (ui/field_step.py) -- the actual field study. Disabled until
    Passo A produces a star with a valid window. No override.

The one rule that can't be violated: this app never blocks waiting on
Castro. Runs (including Passo A's short evolution -- it IS a Castro run)
are launched as detached subprocesses (setsid/nohup) with stdout/stderr
to a fixed log file; a QTimer polls the log periodically. Closing and
reopening the window must reconnect to any run still going (TODO: the
reconnect-after-reopen path itself is not yet implemented -- today a
fresh launch always starts a new run; see core/persistence.py once
Fase 4 lands the shared results store, which is what reconnect will read
from).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QSplitter, QTabWidget, QWidget, QVBoxLayout,
)
from PyQt6.QtCore import Qt, QTimer

from ui.star_step import StarStepWidget
from ui.field_step import FieldStepWidget
from ui.run_queue import RunQueue
from ui.per_seed_view import PerSeedView
from ui.aggregate_view import AggregateView
from ui.field_view import FieldView
from core.run_launcher import RunSpec, launch, read_progress, WD_BRAITHWAITE_DIR
from core.star_builder import find_measurement_window, parse_rho_c_log
from core import scf_store
from core import persistence
from core import extraction

POLL_INTERVAL_MS = 3000
STAR_RUN_ID = "star_build"
STAR_STOP_TIME_S = 0.4  # covers t/t_dyn~1.45 for the validated t_dyn=0.2758s star;
                          # TODO (Fase 5): compute analytically per rho_c instead of
                          # this fixed value once the app builds stars other than
                          # the one validated this session.


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Braithwaite — estudo de estabilidade (wd-magnetizada)")
        self.resize(1050, 900)

        control_tab = QWidget()
        layout = QVBoxLayout(control_tab)

        splitter = QSplitter(Qt.Orientation.Vertical)
        self.star_step = StarStepWidget()
        self.field_step = FieldStepWidget()
        self.run_queue = RunQueue()
        splitter.addWidget(self.star_step)
        splitter.addWidget(self.field_step)
        splitter.addWidget(self.run_queue)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 1)
        splitter.setStretchFactor(2, 1)
        layout.addWidget(splitter)

        self.per_seed_view = PerSeedView(WD_BRAITHWAITE_DIR)
        self.aggregate_view = AggregateView(WD_BRAITHWAITE_DIR)
        self.field_view = FieldView(WD_BRAITHWAITE_DIR)

        tabs = QTabWidget()
        tabs.addTab(control_tab, "Estudo")
        tabs.addTab(self.per_seed_view, "Semente")
        tabs.addTab(self.aggregate_view, "Agregado")
        tabs.addTab(self.field_view, "Campo")
        self.setCentralWidget(tabs)

        self.star_step.star_ready.connect(self.field_step.set_star)
        self.star_step.star_invalidated.connect(lambda: self.field_step.set_star(None))
        self.star_step.build_requested.connect(self._on_star_build_requested)
        self.field_step.run_requested.connect(self._on_run_requested)

        self._active_runs = {}  # run_id -> {"log_path": Path, "kind": "star"|"seed", "target_ttdyn": float}

        self._timer = QTimer(self)
        self._timer.timeout.connect(self._poll_runs)
        self._timer.start(POLL_INTERVAL_MS)

    def _on_star_build_requested(self, params: dict):
        rho_c, mu_e = params["rho_c"], params["mu_e"]

        # Real SCF solve, synchronous -- measured 16-25ms for the
        # field-free case this session, confirmed by test to run well
        # under any threshold that would need backgrounding. This is what
        # makes rho_c/mu_e an actual control instead of the silent
        # stale-file reuse the review caught.
        try:
            manifest = scf_store.build_model_dat(
                rho_c, mu_e, out_path=WD_BRAITHWAITE_DIR / "model.dat"
            )
        except RuntimeError as e:
            self.star_step.on_run_failed(f"SCF não convergiu: {e}")
            return

        # Neutronization + VE gates apply to the CONVERGED result, not
        # just the pre-flight target -- same policy as Tab 3's export (R5),
        # no override.
        neut = scf_store.neutronization_check(rho_c, mu_e)
        if not neut["ok"]:
            self.star_step.on_run_failed(
                f"Estrela convergida está acima do limiar de neutronização "
                f"({rho_c:.3e} >= {neut['threshold']:.3e}) -- não deveria ter "
                "chegado até aqui; o portão da tela já deveria ter bloqueado."
            )
            return
        if manifest["VE"] >= 1e-3:
            self.star_step.on_run_failed(
                f"VE = {manifest['VE']:.3e} >= 1e-3 (critério V3 do plano) -- "
                "equilíbrio não confiável o bastante para virar condição "
                "inicial do Castro. Sem opção de forçar."
            )
            return

        spec = RunSpec(
            run_id=STAR_RUN_ID,
            n_cell=64,
            rng_seed=0,          # unused: e_mag_over_w=0.0 means no field is generated
            e_mag_over_w=0.0,    # star construction: field-free, per the corrected
                                  # premise -- no pre-relaxation, just enough evolution
                                  # to find the measurement window
            stop_time_s=STAR_STOP_TIME_S,
            max_step=2000,
            plot_int=10,
            check_int=200,
            np=1,
        )
        log_path = launch(spec)
        self._active_runs[STAR_RUN_ID] = {
            "log_path": log_path, "kind": "star",
            "target_ttdyn": 1.45, "rho_c_ic": None, "t_dyn_s": None,
            "manifest": manifest, "mu_e": mu_e,
        }
        self.run_queue.add_run(STAR_RUN_ID, resolution=64, target_ttdyn=1.45)
        self.star_step.on_run_started(STAR_RUN_ID)

    def _on_run_requested(self, config: dict):
        star = config["star"]
        cache_key = star["cache_key"]
        resolution = config["resolution"]
        lo, hi = star["window"]
        t_dyn_s = star["t_dyn_s"]
        # small margin over the window's upper bound so the chosen
        # extraction plotfile isn't the very last (possibly still-flushing)
        # step of the run
        stop_time_s = (hi + 0.05) * t_dyn_s

        emag_min, emag_max = config["emag_over_w_min"], config["emag_over_w_max"]
        n_seeds = config["n_seeds"]
        if emag_min == emag_max:
            emag_values = [emag_min] * n_seeds
        else:
            step = (emag_max - emag_min) / max(1, n_seeds - 1)
            emag_values = [emag_min + i * step for i in range(n_seeds)]

        seeds = _pick_new_seeds(cache_key, emag_min, emag_max, n_seeds)

        for seed, e_mag in zip(seeds, emag_values):
            run_id = f"seed{seed}_{cache_key[:6]}"
            spec = RunSpec(
                run_id=run_id, n_cell=resolution, rng_seed=seed, e_mag_over_w=e_mag,
                stop_time_s=stop_time_s, max_step=2000, plot_int=10, check_int=200,
                np=4 if resolution == 128 else 1,
            )
            log_path = launch(spec)
            self._active_runs[run_id] = {
                "log_path": log_path, "kind": "seed",
                "rho_c_ic": None, "t_dyn_s": None, "w_abs_erg": None,
                "star_cache_key": cache_key, "resolution": resolution,
                "e_mag_over_w_target": e_mag, "seed": seed, "window": (lo, hi),
            }
            self.run_queue.add_run(run_id, resolution=resolution, target_ttdyn=hi)

        self.field_step.status_label.setText(f"{len(seeds)} runs de campo lançadas: sementes {seeds}")

    def _poll_runs(self):
        for run_id, info in list(self._active_runs.items()):
            progress = read_progress(info["log_path"], run_id)

            if info["rho_c_ic"] is None and progress.rho_c is not None:
                info["rho_c_ic"] = progress.rho_c
            if info["t_dyn_s"] is None and progress.t_dyn_s is not None:
                info["t_dyn_s"] = progress.t_dyn_s
            if info.get("w_abs_erg") is None and progress.w_abs_erg is not None:
                info["w_abs_erg"] = progress.w_abs_erg

            self.run_queue.update_run(run_id, progress.status, progress.t_ttdyn)

            if progress.status == "error":
                self._finish_run(run_id, info, error=progress.error_text)
            elif progress.status == "ended":
                self._finish_run(run_id, info, error=None)

    def _finish_run(self, run_id: str, info: dict, error: str | None):
        del self._active_runs[run_id]
        if info["kind"] == "seed":
            self._finish_seed_run(run_id, info, error)
            return

        if error:
            self.star_step.on_run_failed(error)
            return

        if info["t_dyn_s"] is None or info["rho_c_ic"] is None:
            self.star_step.on_run_failed("run terminou sem produzir série de rho_c utilizável")
            return
        rho_series = parse_rho_c_log(info["log_path"], info["t_dyn_s"])
        if not rho_series:
            self.star_step.on_run_failed("run terminou sem produzir série de rho_c utilizável")
            return

        window_result = find_measurement_window(rho_series, info["rho_c_ic"])
        manifest = info.get("manifest", {})
        cache_key = persistence.star_cache_key(
            rho_c=info["rho_c_ic"], mu_e=info.get("mu_e", 2.0), resolution=64,
            scf_params_hash=manifest.get("scf_params_hash", ""),
            git_commit_scf=manifest.get("git_commit_scf", ""),
            git_commit_wd_braithwaite=persistence.store.git_commit_hash(WD_BRAITHWAITE_DIR),
        )
        persistence.save_star_result(
            cache_key, rho_c=info["rho_c_ic"], mu_e=info.get("mu_e", 2.0), resolution=64,
            window_result=window_result, VE=manifest.get("VE", float("nan")),
            scf_params_hash=manifest.get("scf_params_hash", ""),
            git_commit_scf=manifest.get("git_commit_scf", ""),
        )
        self.star_step.on_result({
            "rho_c_measured": info["rho_c_ic"],
            "VE": manifest.get("VE", float("nan")),
            "window_result": window_result,
            "cache_key": cache_key,
            "resolution": 64,
            "t_dyn_s": info["t_dyn_s"],
        })
        self.aggregate_view.reload_star_list()

    def _finish_seed_run(self, run_id: str, info: dict, error: str | None):
        if error:
            self.field_step.status_label.setText(f"Run '{run_id}' falhou: {error}")
            return
        if info.get("t_dyn_s") is None or info.get("w_abs_erg") is None:
            self.field_step.status_label.setText(
                f"Run '{run_id}' terminou sem t_dyn/|W| utilizáveis no log."
            )
            return

        lo, hi = info["window"]
        plotfile = _find_plotfile_in_window(
            WD_BRAITHWAITE_DIR, run_id, info["t_dyn_s"], lo, hi
        )
        if plotfile is None:
            self.field_step.status_label.setText(
                f"Run '{run_id}': nenhum plotfile caiu dentro da janela [{lo:.3f}, {hi:.3f}]."
            )
            return

        measurement = extraction.extract_field_measurement(
            plotfile, info["w_abs_erg"], info["t_dyn_s"]
        )
        persistence.save_seed_result(
            info["star_cache_key"], seed=info["seed"], resolution=info["resolution"],
            e_mag_over_w_target=info["e_mag_over_w_target"], measurement=measurement,
            plotfile_path=str(plotfile),
        )
        self.field_step.status_label.setText(
            f"Run '{run_id}' (seed {info['seed']}): E_tor/E_mag = "
            f"{measurement['E_tor_over_Emag']:.4f} em t/t_dyn={measurement['t_ttdyn']:.3f}"
        )

        # Feed the two plot views: per-seed trajectory (selected explicitly
        # so the just-finished run is what the user sees first) and the
        # aggregate ranked scatter (re-reads the persisted store, which
        # save_seed_result above just updated).
        self.per_seed_view.register_run(
            run_id, info["log_path"], info["t_dyn_s"], info["w_abs_erg"],
            info["window"], info["rho_c_ic"],
        )
        self.per_seed_view.selector.setCurrentText(run_id)
        self.aggregate_view.reload_star_list()


def _pick_new_seeds(star_cache_key: str, emag_min: float, emag_max: float, n_seeds: int) -> list[int]:
    """Sequential integer seeds (1, 2, 3, ...), skipping any already
    persisted for this star at this field-energy target -- avoids
    silently re-running (and double-counting) a seed already measured.
    """
    existing = persistence.load_seeds_for_star(star_cache_key)
    used = set()
    if not existing.empty:
        same_field = existing[
            (existing["e_mag_over_w_target"] >= emag_min - 1e-9)
            & (existing["e_mag_over_w_target"] <= emag_max + 1e-9)
        ]
        used = {int(s) for s in same_field["seed"].dropna().tolist()}

    seeds = []
    candidate = 1
    while len(seeds) < n_seeds:
        if candidate not in used:
            seeds.append(candidate)
        candidate += 1
    return seeds


def _find_plotfile_in_window(run_dir: Path, run_id: str, t_dyn_s: float,
                              lo_ttdyn: float, hi_ttdyn: float) -> Path | None:
    """Picks the plotfile with the largest t/t_dyn that still falls
    inside [lo_ttdyn, hi_ttdyn] -- scanning from the newest plotfile
    backward (plot numbering correlates with time) rather than reading
    every plotfile's header.
    """
    candidates = sorted(run_dir.glob(f"plt_{run_id}0*"), reverse=True)
    for plotfile in candidates:
        info = extraction.run_finterior("density", plotfile)
        t_ttdyn = info["time"] / t_dyn_s
        if lo_ttdyn <= t_ttdyn <= hi_ttdyn:
            return plotfile
        if t_ttdyn < lo_ttdyn:
            break  # scanning backward past the window's start -- stop
    return None


def main():
    app = QApplication(sys.argv)
    win = MainWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
