"""Writes a magnetized, axisymmetric SCF equilibrium as a 2D meridional model
file carrying the VECTOR POTENTIAL, for a Castro or FLASH problem to read.

Why not castro_model_writer.py. That one asserts check_field_free_non_rotating
and check_spherical_symmetry and then averages over theta to make a 1D radial
profile. It cannot be extended by relaxing a gate: averaging over theta is
exactly the operation that destroys what has to be exported here. This is a
separate writer for a separate job.

Why the potential and not the field. Interpolating B onto faces leaves a
divergence that constrained transport then carries forever. Writing A and
taking the DISCRETE curl on the staggered mesh makes div B = 0 an identity of
the discretisation instead of an error to be cleaned. For an axisymmetric
field the potential is analytic, so nothing is fitted:

    poloidal   A_phi = u / varpi           (the flux function IS the potential)
    toroidal   A_z   = -int_0^varpi B_phi(varpi', z) dvarpi'

with A_varpi = 0. Then, in cylindrical coordinates,

    B_varpi = -(1/varpi) du/dz      B_z = (1/varpi) du/dvarpi     [from A_phi]
    B_phi   = -dA_z/dvarpi                                        [from A_z]

which returns the field it was built from. investigations/vector_potential_
export.py measured that construction end to end and got div B = 1.7e-16
normalised; this module is that construction packaged, with the verification
kept as a gate rather than as a one-off.

FILE FORMAT (plain text, structured, easy for both C++ and Fortran)

    # comment lines, provenance
    NVARPI NZ
    varpi_0 ... varpi_{NVARPI-1}       one value per line
    z_0     ... z_{NZ-1}               one value per line
    density A_phi A_z                  NVARPI*NZ lines, z index fastest

Units are CGS throughout: cm, g/cm^3, and A in gauss*cm.
"""

import json
import subprocess
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator as RGI

CGS_NOTE = "cgs: varpi,z in cm; density in g/cm^3; A_phi,A_z in G*cm"


def to_meridional(r, theta, fields, vp, zz):
    """Resample (r, theta) fields onto the meridional (varpi, z) grid."""
    VP, ZZ = np.meshgrid(vp, zz, indexing="ij")
    RR = np.sqrt(VP ** 2 + ZZ ** 2)
    TT = np.arccos(np.clip(ZZ / np.maximum(RR, 1e-30), -1.0, 1.0))
    return [RGI((r, theta), f, bounds_error=False, fill_value=0.0)((RR, TT))
            for f in fields]


def vector_potential(vp, u_m, bphi_m):
    """A_phi = u/varpi and A_z = -int_0^varpi B_phi dvarpi', on the grid.

    The A_z integral is trapezoidal and cumulative from the axis outward,
    which is the same direction the analytic definition integrates in; no
    constant of integration is free, because A_z(0, z) = 0 by construction.
    """
    VP = vp[:, None]
    A_phi = np.where(VP > 0.0, u_m / np.maximum(VP, 1e-30), 0.0)
    A_z = -np.concatenate(
        [np.zeros((1, u_m.shape[1])),
         np.cumsum(0.5 * (bphi_m[1:] + bphi_m[:-1]) * np.diff(vp)[:, None],
                   axis=0)], axis=0)
    return A_phi, A_z


def verify_curl_on_cartesian(vp, zz, A_phi, A_z, half, n_cart=64):
    """Sample A on the staggered Cartesian mesh, curl it, and measure.

    Returns (rel_divB, b_max, dx). rel_divB is |div B|max normalised by
    max|B|/dx -- the only scale that makes a divergence dimensionless here.
    The mesh uses the half-shift geometry (a cell CENTRE at r = 0) that the
    Castro problems in this repository require.
    """
    f_phi = RGI((vp, zz), A_phi, bounds_error=False, fill_value=0.0)
    f_z = RGI((vp, zz), A_z, bounds_error=False, fill_value=0.0)

    dx = 2.0 * half / n_cart
    lo = -((n_cart + 1) / 2.0) * dx
    ctr = lo + dx * (np.arange(n_cart) + 0.5)
    fac = lo + dx * np.arange(n_cart + 1)

    def A_at(X, Y, Z, comp):
        vpg = np.sqrt(X ** 2 + Y ** 2)
        pts = np.stack([vpg.ravel(), Z.ravel()], axis=-1)
        ap = f_phi(pts).reshape(X.shape)
        az = f_z(pts).reshape(X.shape)
        inv = np.where(vpg > 0, 1.0 / np.maximum(vpg, 1e-30), 0.0)
        return {"x": -ap * Y * inv, "y": ap * X * inv, "z": az}[comp]

    Ax = A_at(*np.meshgrid(ctr, fac, fac, indexing="ij"), "x")
    Ay = A_at(*np.meshgrid(fac, ctr, fac, indexing="ij"), "y")
    Az = A_at(*np.meshgrid(fac, fac, ctr, indexing="ij"), "z")

    Bx = (np.diff(Az, axis=1) - np.diff(Ay, axis=2)) / dx
    By = (np.diff(Ax, axis=2) - np.diff(Az, axis=0)) / dx
    Bz = (np.diff(Ay, axis=0) - np.diff(Ax, axis=1)) / dx

    div = (np.diff(Bx, axis=0) + np.diff(By, axis=1)
           + np.diff(Bz, axis=2)) / dx
    b_max = max(np.abs(Bx).max(), np.abs(By).max(), np.abs(Bz).max())
    rel = np.abs(div).max() / (b_max / dx)
    return float(rel), float(b_max), float(dx)


def verify_meridional_curl(vp, zz, A_phi, A_z, u_m, bphi_m):
    """Curl A back on the meridional grid and compare with the input field.

    Independent of verify_curl_on_cartesian: that one checks div B = 0, which
    the construction guarantees whether or not A is the RIGHT potential. This
    checks that it is. Returns relative sup-norm errors (poloidal, toroidal),
    measured over the interior only, since one-sided differences at the axis
    and at the outer edge are first order and would dominate the norm.
    """
    s = (slice(1, -1), slice(1, -1))
    dvp, dz = vp[1] - vp[0], zz[1] - zz[0]
    VP = vp[:, None]

    # from A_phi: B_varpi = -(1/varpi) du/dz, B_z = (1/varpi) du/dvarpi
    u_from_A = A_phi * VP
    b_vp = -np.gradient(u_from_A, dz, axis=1) / np.maximum(VP, 1e-30)
    b_z = np.gradient(u_from_A, dvp, axis=0) / np.maximum(VP, 1e-30)
    b_vp_ref = -np.gradient(u_m, dz, axis=1) / np.maximum(VP, 1e-30)
    b_z_ref = np.gradient(u_m, dvp, axis=0) / np.maximum(VP, 1e-30)
    scale_p = max(np.abs(b_vp_ref[s]).max(), np.abs(b_z_ref[s]).max(), 1e-300)
    err_pol = max(np.abs((b_vp - b_vp_ref)[s]).max(),
                  np.abs((b_z - b_z_ref)[s]).max()) / scale_p

    # from A_z: B_phi = -dA_z/dvarpi
    b_phi = -np.gradient(A_z, dvp, axis=0)
    scale_t = max(np.abs(bphi_m[s]).max(), 1e-300)
    err_tor = np.abs((b_phi - bphi_m)[s]).max() / scale_t
    return float(err_pol), float(err_tor)


def _git_hash(path):
    try:
        return subprocess.run(
            ["git", "-C", str(path), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=10).stdout.strip() or "no-git"
    except Exception:
        return "no-git"


def write_model(vp, zz, rho_m, A_phi, A_z, out_path, params,
                diagnostics=None, density_floor=0.0, v_phi=None):
    """Write the model file and a <path>.manifest.json sidecar.

    Returns the manifest dict. Nothing is validated here -- call the verify_*
    functions and put their results in `diagnostics`, so what was checked is
    recorded next to the data rather than asserted in a docstring.
    """
    out_path = Path(out_path)
    rho_out = np.maximum(rho_m, density_floor)

    with out_path.open("w") as f:
        f.write("# magnetized axisymmetric white dwarf model\n")
        f.write(f"# {CGS_NOTE}\n")
        if v_phi is None:
            f.write("# columns: density A_phi A_z, z index fastest\n")
        else:
            # v2 adds the azimuthal velocity. Rotation cannot be carried by
            # the vector potential, and Castro's rotation support is a
            # rotating frame at constant Omega, which a j-constant law is
            # not -- so differential rotation has to arrive as an initial
            # velocity field and be allowed to evolve.
            f.write("# columns: density A_phi A_z v_phi, z index fastest\n")
        for k, v in sorted(params.items()):
            f.write(f"# param {k} = {v}\n")
        if diagnostics:
            for k, v in sorted(diagnostics.items()):
                f.write(f"# check {k} = {v}\n")
        f.write(f"{len(vp)} {len(zz)}\n")
        for x in vp:
            f.write(f"{x:.10e}\n")
        for x in zz:
            f.write(f"{x:.10e}\n")
        for i in range(len(vp)):
            for j in range(len(zz)):
                row = (f"{rho_out[i, j]:.10e} {A_phi[i, j]:.10e} "
                       f"{A_z[i, j]:.10e}")
                if v_phi is not None:
                    row += f" {v_phi[i, j]:.10e}"
                f.write(row + "\n")

    manifest = {
        "file": out_path.name,
        "format": ("axisym-vector-potential-v1" if v_phi is None
                   else "axisym-vector-potential-v2"),
        "has_velocity": v_phi is not None,
        "units": CGS_NOTE,
        "n_varpi": int(len(vp)),
        "n_z": int(len(zz)),
        "varpi_max_cm": float(vp[-1]),
        "z_range_cm": [float(zz[0]), float(zz[-1])],
        "params": {k: (float(v) if isinstance(v, (int, float, np.floating))
                       else v) for k, v in params.items()},
        "checks": diagnostics or {},
        "git_commit": _git_hash(Path(__file__).resolve().parent.parent),
    }
    with out_path.with_suffix(out_path.suffix + ".manifest.json").open("w") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)
    return manifest
