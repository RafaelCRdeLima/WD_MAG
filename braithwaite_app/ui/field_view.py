"""Campo -- spatial visualization of the magnetic field inside one
plotfile: a 2D matplotlib slice (fast, always available) and a 3D
pyvista streamline view (heavier -- built on demand by a button, not on
every selector change, since it reads the full 3D grid via `yt`).

Purely visual/diagnostic -- reads plotfiles directly (core/field_reader.py,
yt-based) and never writes to core/persistence.py's results store; the
numbers behind the science stay on the finterior-verified path in
core/extraction.py (see that module's docstring for why the two paths
must not merge).
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QApplication, QCheckBox, QComboBox, QHBoxLayout, QLabel, QPushButton,
    QSizePolicy, QSpinBox, QSplitter, QVBoxLayout, QWidget,
)
from pyvistaqt import QtInteractor

from core import field_reader

_SLICE_AXES = ["z", "y", "x"]
_3D_COLOR_OPTIONS = ["|B|", "fração toroidal (E_tor/E_mag local)"]
_SEED_SHELLS = ((0.3, 0.2), (0.6, 0.3), (0.9, 0.5))  # (fraction of R_star, share of seeds) -- weighted toward the surface, since that's where a line has to start to show whether it closes outside the star (the exterior-dipole question) without first burning its length crossing the interior


def _fibonacci_sphere_directions(n: int) -> np.ndarray:
    """n unit vectors spread evenly over a sphere (golden-angle spiral) --
    deterministic, no RNG, so the same plotfile always seeds the same
    lines instead of a different tangle on every click of 'gerar'."""
    if n <= 1:
        return np.array([[0.0, 0.0, 1.0]])
    i = np.arange(n)
    y = 1 - 2 * i / (n - 1)
    r_xy = np.sqrt(np.clip(1 - y**2, 0, None))
    golden_angle = np.pi * (3 - np.sqrt(5))
    theta = golden_angle * i
    return np.stack([np.cos(theta) * r_xy, y, np.sin(theta) * r_xy], axis=1)


def _seed_points(center: tuple, r_star: float, n_seeds: int) -> np.ndarray:
    """Seeds spread across three depths (30/60/90% of R_star), weighted
    toward the surface shell, rather than one shell at a domain-scale
    radius -- that single-shell-outside-the-star placement was the bug
    behind the earlier short-fragment lines: those seeds landed in the
    near-vacuum exterior tail, where the field is too weak to integrate
    far in either direction."""
    center = np.asarray(center)
    shells = [
        center + _fibonacci_sphere_directions(max(1, round(n_seeds * share))) * (r_star * frac)
        for frac, share in _SEED_SHELLS
    ]
    return np.vstack(shells)


class FieldView(QWidget):
    def __init__(self, run_dir: Path, parent=None):
        super().__init__(parent)
        self.run_dir = Path(run_dir)
        self._plotfiles: list[Path] = []

        root = QVBoxLayout(self)
        root.setSpacing(3)

        # --- controls: two compact rows, everything above the panels ---
        row1 = QHBoxLayout()
        row1.addWidget(QLabel("run:"))
        self.run_selector = QComboBox()
        self.run_selector.currentTextChanged.connect(self._on_run_selected)
        row1.addWidget(self.run_selector, 2)
        row1.addWidget(QLabel("passo:"))
        self.step_selector = QComboBox()
        self.step_selector.currentTextChanged.connect(lambda _t: self._replot_slice())
        row1.addWidget(self.step_selector, 2)
        reload_btn = QPushButton("recarregar lista de runs")
        reload_btn.clicked.connect(self.reload_runs)
        row1.addWidget(reload_btn)
        root.addLayout(row1)

        row2 = QHBoxLayout()
        row2.addWidget(QLabel("campo:"))
        self.field_selector = QComboBox()
        self.field_selector.addItems(list(field_reader.FIELD_SPECS.keys()))
        self.field_selector.currentTextChanged.connect(lambda _t: self._replot_slice())
        row2.addWidget(self.field_selector)
        row2.addWidget(QLabel("corte:"))
        self.axis_selector = QComboBox()
        self.axis_selector.addItems(_SLICE_AXES)
        self.axis_selector.currentTextChanged.connect(lambda _t: self._replot_slice())
        row2.addWidget(self.axis_selector)
        self.interior_only_check = QCheckBox("só interior da estrela")
        self.interior_only_check.stateChanged.connect(lambda _s: self._replot_slice())
        row2.addWidget(self.interior_only_check)

        row2.addSpacing(20)
        gen_btn = QPushButton("gerar linhas de campo 3D")
        gen_btn.clicked.connect(self._plot_3d)
        row2.addWidget(gen_btn)
        row2.addWidget(QLabel("colorir por:"))
        self.color3d_selector = QComboBox()
        self.color3d_selector.addItems(_3D_COLOR_OPTIONS)
        row2.addWidget(self.color3d_selector)
        row2.addWidget(QLabel("linhas:"))
        self.n_lines_spin = QSpinBox()
        self.n_lines_spin.setRange(6, 90)
        self.n_lines_spin.setValue(12)
        row2.addWidget(self.n_lines_spin)
        row2.addStretch(1)
        root.addLayout(row2)

        # thin single-line header, not a floating word-wrapped block --
        # this used to eat a huge chunk of vertical space by accident
        self.status_label = QLabel("")
        self.status_label.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Fixed)
        self.status_label.setMaximumHeight(20)
        root.addWidget(self.status_label)

        # second line, dedicated to the shape/field-extreme numbers
        # (R_eq/R_pol, Bt_max/Bp_max) populated by "gerar linhas de campo
        # 3D" -- kept separate from status_label so the run/plotfile
        # message above doesn't get too long to read at a glance.
        self.shape_label = QLabel("")
        self.shape_label.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Fixed)
        self.shape_label.setMaximumHeight(20)
        root.addWidget(self.shape_label)

        # --- the two visualization panels, filling essentially everything else ---
        splitter = QSplitter(Qt.Orientation.Horizontal)

        self.figure = Figure()
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        splitter.addWidget(self.canvas)

        self.plotter = QtInteractor(splitter)
        self.plotter.interactor.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        splitter.addWidget(self.plotter.interactor)

        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 1)
        root.addWidget(splitter, 1)  # stretch=1: this row absorbs all leftover space, not the controls above it

        self.reload_runs()

    def reload_runs(self):
        current = self.run_selector.currentText()
        try:
            ids = field_reader.list_run_ids(self.run_dir)
        except OSError as e:
            self.status_label.setText(f"Erro lendo {self.run_dir}: {e}")
            return
        self.run_selector.blockSignals(True)
        self.run_selector.clear()
        self.run_selector.addItems(ids)
        if current and self.run_selector.findText(current) >= 0:
            self.run_selector.setCurrentText(current)
        self.run_selector.blockSignals(False)
        self._on_run_selected(self.run_selector.currentText())

    def _on_run_selected(self, run_id: str):
        self.step_selector.blockSignals(True)
        self.step_selector.clear()

        if not run_id:
            self._plotfiles = []
            self.step_selector.blockSignals(False)
            self._replot_slice()
            return

        self.status_label.setText(f"Filtrando passos válidos de '{run_id}'...")
        QApplication.processEvents()  # paint the status line before the per-plotfile time reads below

        all_plotfiles = field_reader.list_plotfiles_for_run(self.run_dir, run_id)
        t_dyn_s = field_reader.read_run_t_dyn(self.run_dir, run_id)
        window = (
            field_reader.find_validity_window_for_run(self.run_dir, run_id, t_dyn_s)
            if t_dyn_s else None
        )

        if window is not None:
            lo, hi = window
            in_window = []
            for pf in all_plotfiles:
                try:
                    t_ttdyn = field_reader.plotfile_time_s(pf) / t_dyn_s
                except Exception:
                    continue  # unreadable/partially-written plotfile -- skip, don't break the whole list
                if lo <= t_ttdyn <= hi:
                    in_window.append(pf)
            self._plotfiles = in_window
            if in_window:
                self.status_label.setText(
                    f"'{run_id}': {len(in_window)}/{len(all_plotfiles)} passos dentro da janela "
                    f"de validade [{lo:.3f}, {hi:.3f}] t/t_dyn"
                )
            else:
                self.status_label.setText(
                    f"'{run_id}': nenhum dos {len(all_plotfiles)} passos cai dentro da janela de "
                    f"validade [{lo:.3f}, {hi:.3f}] t/t_dyn (estrela ainda não chegou lá, ou já passou)"
                )
        else:
            # No field-free "star*" reference with a matching t_dyn was found on disk --
            # showing every step unfiltered is safer than silently guessing a window.
            self._plotfiles = all_plotfiles
            self.status_label.setText(
                f"'{run_id}': janela de validade não encontrada (sem run campo-livre 'star*' de "
                f"referência com o mesmo t_dyn em disco) -- mostrando todos os {len(all_plotfiles)} "
                "passos, sem filtro"
            )

        # in_window/all_plotfiles are already in chronological (step) order --
        # list_plotfiles_for_run sorts by step number -- so this preserves that order.
        self.step_selector.addItems([p.name for p in self._plotfiles])
        if self._plotfiles:
            self.step_selector.setCurrentIndex(len(self._plotfiles) - 1)  # latest step still WITHIN the window
        self.step_selector.blockSignals(False)
        self._replot_slice()

    def _current_plotfile(self) -> Path | None:
        idx = self.step_selector.currentIndex()
        if idx < 0 or idx >= len(self._plotfiles):
            return None
        return self._plotfiles[idx]

    def _replot_slice(self):
        plotfile = self._current_plotfile()
        self.figure.clear()
        if plotfile is None:
            self.canvas.draw()
            return
        try:
            s = field_reader.load_slice(
                plotfile, self.field_selector.currentText(), self.axis_selector.currentText(),
            )
        except Exception as e:  # pragma: no cover -- surfaced to the user, not swallowed
            self.status_label.setText(f"Erro no corte 2D de '{plotfile.name}': {e}")
            self.canvas.draw()
            return

        ax = self.figure.add_subplot(111)
        array = s["array"]
        if self.interior_only_check.isChecked():
            array = np.where(s["interior_mask"], array, np.nan)

        if s["diverging"]:
            cmap = plt.get_cmap("coolwarm").copy()
            norm = None
        else:
            cmap = plt.get_cmap("inferno").copy()
            norm = "log" if s["log"] else None
        vmin, vmax = s["interior_vrange"]
        cmap.set_bad(color="white", alpha=0.0)  # NaN (masked toroidal_frac exterior, or the "só interior" toggle) -> blank, not a solid color

        im = ax.imshow(
            array, origin="lower", extent=s["extent_cm"], cmap=cmap,
            norm=norm, vmin=vmin, vmax=vmax,
        )
        self.figure.colorbar(im, ax=ax, label=f"{s['field']} ({s['unit']}) — escala ajustada ao interior da estrela")

        # Stellar-surface reference: without this, a slice of any field
        # gives no way to tell what's inside the star vs. the field-free
        # exterior/vacuum floor around it -- a density contour at a small
        # fraction of the slice's own peak marks that boundary regardless
        # of which physical field is being displayed.
        density = s["density_array"]
        rho_peak = float(np.nanmax(density))
        if rho_peak > 0:
            ax.contour(
                density, levels=[rho_peak * 1e-3], extent=s["extent_cm"], origin="lower",
                colors="lime" if s["diverging"] else "white", linewidths=1.3, linestyles="--",
            )
            ax.plot([], [], "--", color="lime" if s["diverging"] else "white",
                    linewidth=1.3, label="borda da estrela (densidade)")
            ax.legend(fontsize=7, loc="upper right", facecolor="0.3", labelcolor="white")

        ax.set_xlabel("cm")
        ax.set_ylabel("cm")
        ax.set_title(f"corte em {s['axis']} — t={s['time_s']:.4g} s")
        self.figure.tight_layout()
        self.canvas.draw()
        self.status_label.setText(f"'{plotfile.name}': t={s['time_s']:.4g} s")

    def _plot_3d(self):
        plotfile = self._current_plotfile()
        if plotfile is None:
            self.status_label.setText("Nenhum plotfile selecionado.")
            return

        self.status_label.setText(f"Lendo grade 3D de '{plotfile.name}'...")
        QApplication.processEvents()  # paint the status line before the blocking yt read below

        try:
            g = field_reader.load_vector_grid(plotfile)
        except Exception as e:  # pragma: no cover -- surfaced to the user, not swallowed
            self.status_label.setText(f"Erro lendo grade 3D de '{plotfile.name}': {e}")
            return

        nx, ny, nz = g["dims"]
        ox, oy, oz = g["origin_cm"]
        sx, sy, sz = g["spacing_cm"]
        # first cell CENTER (covering_grid data is cell-centered; ImageData
        # points are placed at cell centers here so no interpolation is
        # needed between the two)
        point_origin = (ox + sx / 2, oy + sy / 2, oz + sz / 2)
        domain_center = (ox + sx * nx / 2, oy + sy * ny / 2, oz + sz * nz / 2)

        grid = pv.ImageData(dimensions=(nx, ny, nz), spacing=(sx, sy, sz), origin=point_origin)
        bx, by, bz = g["B_x"], g["B_y"], g["B_z"]
        grid["B"] = np.stack(
            [bx.ravel(order="F"), by.ravel(order="F"), bz.ravel(order="F")], axis=1
        )
        grid["Bmag"] = np.sqrt(bx**2 + by**2 + bz**2).ravel(order="F")
        grid["toroidal_frac"] = g["toroidal_frac"].ravel(order="F")
        grid["density"] = g["density"].ravel(order="F")

        self.plotter.clear()

        # Stellar surface -- the physical reference frame everything else
        # is drawn against; without it there's no way to tell which
        # tangle of lines is inside the star and which is exterior tail.
        rho_max = float(g["density"].max())
        surface_added = False
        if rho_max > 0:
            try:
                envelope = grid.contour(isosurfaces=[rho_max * 1e-3], scalars="density")
                if envelope.n_points:
                    self.plotter.add_mesh(
                        envelope, color="lightsteelblue", opacity=0.18,
                        specular=0.3, smooth_shading=True,
                    )
                    surface_added = True
            except Exception:
                pass  # falls back to axes/lines alone -- not fatal

        color_choice = self.color3d_selector.currentText()
        by_toroidal = color_choice.startswith("fração")
        color_kwargs = (
            dict(scalars="toroidal_frac", cmap="coolwarm", clim=(0.0, 1.0),
                 # short title -- the full explanation is the separate
                 # "azul=poloidal / vermelho=toroidal" text added below,
                 # a longer title here collided with the tick labels
                 scalar_bar_args={
                     "title": "E_tor/E_mag", "n_labels": 2,
                     "width": 0.5, "height": 0.08, "position_x": 0.25, "position_y": 0.03,
                 })
            if by_toroidal
            else dict(
                scalars="Bmag", cmap="plasma", log_scale=True,
                # log_scale: |B| spans many orders of magnitude (weak
                # exterior tail to strong interior core) -- same reason
                # the 2D slice above uses a log norm; few labels + a
                # compact fixed-width format keep the tick text from
                # overlapping into unreadable clutter.
                scalar_bar_args={
                    "title": "|B| (G)", "n_labels": 2, "fmt": "%.0e",
                    "width": 0.5, "height": 0.08, "position_x": 0.25, "position_y": 0.03,
                },
            )
        )

        r_star = field_reader.estimate_star_radius(g["density"], g["spacing_cm"])
        axis_radii = field_reader.estimate_axis_radii(
            g["density"], g["dims"], g["spacing_cm"], g["origin_cm"]
        )
        bt_bp = field_reader.estimate_bt_bp_max(
            g["B_x"], g["B_y"], g["B_z"], g["dims"], g["origin_cm"], g["spacing_cm"]
        )
        lines_ok = False
        try:
            seeds = _seed_points(domain_center, r_star * 0.95, self.n_lines_spin.value())
            source = pv.PolyData(seeds)
            # Real flow-line integration (VTK's adaptive RK45 stream
            # tracer), not a short local sample: 'cl' step units make
            # each step a fraction of one cell regardless of |B|, both
            # directions from every seed (integration_direction='both'
            # -- so a line that starts mid-arc still gets its other
            # half), and max_length set from the star's own size so a
            # line can complete a multi-R_star loop instead of being cut
            # off at an arbitrary domain-scale distance.
            streamlines = grid.streamlines_from_source(
                source, vectors="B", integration_direction="both",
                integrator_type=45, step_unit="cl",
                initial_step_length=0.2, min_step_length=0.01, max_step_length=0.5,
                max_steps=5000, max_length=max(r_star, 1.0) * 3.5,
                terminal_speed=1e-12,
            )
            if streamlines.n_points:
                tube_radius = max(r_star, 1.0) * 0.006
                self.plotter.add_mesh(streamlines.tube(radius=tube_radius), **color_kwargs)
                lines_ok = True
        except Exception as e:  # pragma: no cover -- 3D view degrades, doesn't crash the app
            self.status_label.setText(f"Linhas de campo falharam em '{plotfile.name}': {e}")

        if by_toroidal:
            self.plotter.add_text(
                "azul = poloidal (0)   vermelho = toroidal (1)",
                position="lower_left", font_size=9, color="black",
            )

        self.plotter.show_bounds(
            grid=True, location="outer", color="black",
            xtitle="x (cm)", ytitle="y (cm)", ztitle="z (cm)",
            n_xlabels=3, n_ylabels=3, n_zlabels=3, fmt="%.1e",
        )
        self.plotter.reset_camera()

        status = f"3D: '{plotfile.name}', t={g['time_s']:.4g} s"
        if not surface_added:
            status += " (superfície da estrela não pôde ser gerada)"
        if not lines_ok:
            status += " (linhas de campo falharam)"
        self.status_label.setText(status)

        shape_parts = []
        r_eq_km, r_pol_km = axis_radii["r_eq_cm"] / 1e5, axis_radii["r_pol_cm"] / 1e5
        if r_eq_km > 0 and np.isfinite(r_eq_km) and np.isfinite(r_pol_km):
            ratio = r_pol_km / r_eq_km
            if abs(ratio - 1.0) < 0.005:
                shape = "≈esférica"
            elif ratio < 1.0:
                shape = "oblata"
            else:
                shape = "prolata"
            shape_parts.append(
                f"R_eq={r_eq_km:.1f} km, R_pol={r_pol_km:.1f} km (R_pol/R_eq={ratio:.3f}, {shape})"
            )
        if np.isfinite(bt_bp["ratio"]):
            shape_parts.append(
                f"Bt_max={bt_bp['bt_max_G']:.3e} G, Bp_max={bt_bp['bp_max_G']:.3e} G "
                f"(Bt_max/Bp_max={bt_bp['ratio']:.3f})"
            )
        self.shape_label.setText("  |  ".join(shape_parts))
