"""Field-line tracing on the uniform covering grid of a plotfile.

Used by extract_fields.py. Kept in its own module because it is the one
piece of numerics in these figure scripts that is not a call into the
app: yt's own Streamlines object integrates on the AMR hierarchy, which
this single-level 64^3 data does not need, and it does not give the
per-vertex quantities the figure colors by. The integrator is a plain
RK4 on a trilinearly interpolated field -- no physics, only geometry.

Coordinate conventions match braithwaite_app/core/field_reader.py: the
grid is cell-centered, origin_cm is the domain's lower corner, and the
toroidal direction is azimuthal about z (a bookkeeping axis for this
non-rotating star, not a preferred one).
"""

import numpy as np


def _interp(field, pos, origin, spacing, dims):
    """Trilinear interpolation of a (nx, ny, nz, 3) field at pos (…, 3).
    Returns zeros outside the grid, which stops the tracer."""
    g = (pos - origin) / spacing - 0.5      # cell-centered
    i0 = np.floor(g).astype(int)
    frac = g - i0
    out = np.zeros(pos.shape[:-1] + (field.shape[-1],))
    ok = np.all((i0 >= 0) & (i0 + 1 < np.asarray(dims)), axis=-1)
    if not np.any(ok):
        return out
    i, j, k = i0[ok, 0], i0[ok, 1], i0[ok, 2]
    fx, fy, fz = frac[ok, 0:1], frac[ok, 1:2], frac[ok, 2:3]
    c = 0.0
    for di, wx in ((0, 1 - fx), (1, fx)):
        for dj, wy in ((0, 1 - fy), (1, fy)):
            for dk, wz in ((0, 1 - fz), (1, fz)):
                c = c + wx * wy * wz * field[i + di, j + dj, k + dk]
    out[ok] = c
    return out


def decompose(bx, by, bz, origin, spacing, dims):
    """Total, poloidal and toroidal parts of B as (nx, ny, nz, 3) arrays.

    B_t = (B . e_phi) e_phi with e_phi = (-y, x, 0)/varpi, and
    B_p = B - B_t.  On the axis (varpi = 0) the azimuthal direction is
    undefined and the toroidal part is taken as zero, the same
    convention field_reader.py uses for its Bt field.
    """
    nx, ny, nz = dims
    ax = (origin[0] + spacing[0] * (np.arange(nx) + 0.5))[:, None, None]
    ay = (origin[1] + spacing[1] * (np.arange(ny) + 0.5))[None, :, None]
    varpi = np.sqrt(ax**2 + ay**2) + np.zeros((nx, ny, nz))
    safe = np.where(varpi > 0, varpi, 1.0)
    ephi = np.stack([-(ay + np.zeros_like(varpi)) / safe,
                     (ax + np.zeros_like(varpi)) / safe,
                     np.zeros_like(varpi)], axis=-1)
    ephi[varpi <= 0] = 0.0

    total = np.stack([bx, by, bz], axis=-1)
    bphi = np.sum(total * ephi, axis=-1, keepdims=True)
    tor = bphi * ephi
    return total, total - tor, tor


def trace(field, seeds, origin, spacing, dims, inside, step_cells=0.35,
          max_steps=1200, mag_floor_frac=1e-3):
    """RK4 field lines through `field`, integrated both ways from each
    seed and joined into one polyline per seed.

    Stops a line when it leaves `inside` (the stellar-interior mask, so
    lines do not wander through the vacuum where the field is numerical
    floor), when the local magnitude drops below mag_floor_frac of the
    interior maximum, or at max_steps.

    Returns a list of (n, 3) arrays in cm.
    """
    step = step_cells * float(np.min(spacing))
    mag = np.linalg.norm(field, axis=-1)
    floor = mag[inside].max() * mag_floor_frac
    dims_arr = np.asarray(dims)

    def unit(pos):
        v = _interp(field, pos[None, :], origin, spacing, dims_arr)[0]
        n = np.linalg.norm(v)
        return (v / n, n) if n > 0 else (None, 0.0)

    def in_star(pos):
        g = np.floor((pos - origin) / spacing).astype(int)
        if np.any(g < 0) or np.any(g >= dims_arr):
            return False
        return bool(inside[tuple(g)])

    def march(seed, sign):
        pts = []
        pos = np.array(seed, dtype=float)
        for _ in range(max_steps):
            k1, m = unit(pos)
            if k1 is None or m < floor or not in_star(pos):
                break
            k2, _ = unit(pos + sign * 0.5 * step * k1)
            if k2 is None:
                break
            k3, _ = unit(pos + sign * 0.5 * step * k2)
            if k3 is None:
                break
            k4, _ = unit(pos + sign * step * k3)
            if k4 is None:
                break
            pos = pos + sign * step * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
            pts.append(pos.copy())
        return pts

    lines = []
    for seed in seeds:
        back = march(seed, -1.0)[::-1]
        fwd = march(seed, +1.0)
        line = np.array(back + [np.asarray(seed, dtype=float)] + fwd)
        if len(line) > 8:
            lines.append(line)
    return lines


def sample_along(field, line, origin, spacing, dims):
    """Magnitude of `field` at each vertex of a traced line."""
    v = _interp(field, line, origin, spacing, np.asarray(dims))
    return np.linalg.norm(v, axis=-1)
