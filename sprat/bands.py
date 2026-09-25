"""MPB band structures (mode = bands).

* ``bulk``:  TM and TE bands of the infinite lattice along Gamma-X-M-Gamma, the TM gap and
             the lattice-constant calibration a = f_mid * target_wavelength_nm;
* ``sweep``: the TM gap against the analyte index;
* ``w1``:    the W1 supercell projection (guided bands inside the gap).

MPB is imported inside the functions; ``SPRAT_FAKE_FDTD=1`` returns a synthetic result.
"""

from __future__ import annotations

import os
from typing import Any

import numpy as np


def _mpb():
    try:
        import meep as mp
        from meep import mpb
    except ImportError as e:                      # pragma: no cover
        raise ImportError("MPB (pymeep) is not installed; see README, Installation") from e
    return mp, mpb


def bulk_bands(n: float, rod_radius: float, rod_eps: float, polarisation: str = "tm", num_bands: int = 8,
               k_points: int = 24, resolution: int = 32) -> tuple[np.ndarray, dict[str, Any]]:
    mp, mpb = _mpb()
    geometry = [mp.Cylinder(radius=rod_radius, material=mp.Medium(epsilon=rod_eps), height=mp.inf,
                            axis=mp.Vector3(0, 0, 1))]
    kpts = mp.interpolate(k_points, [mp.Vector3(), mp.Vector3(0.5), mp.Vector3(0.5, 0.5), mp.Vector3()])
    ms = mpb.ModeSolver(geometry=geometry, geometry_lattice=mp.Lattice(size=mp.Vector3(1, 1)),
                        default_material=mp.Medium(index=n), k_points=kpts, resolution=resolution, num_bands=num_bands)
    if polarisation == "tm":
        ms.run_tm()
    else:
        ms.run_te()
    freqs = np.array(ms.all_freqs)
    return freqs, gap_between(freqs, 0, 1)


def gap_between(freqs: np.ndarray, lower_band: int, upper_band: int) -> dict[str, Any]:
    lo, hi = float(freqs[:, lower_band].max()), float(freqs[:, upper_band].min())
    if hi <= lo:
        return dict(exists=False, lower=None, upper=None, mid=None, relative_width_pct=0.0)
    mid = 0.5 * (lo + hi)
    return dict(exists=True, lower=lo, upper=hi, mid=mid, relative_width_pct=100.0 * (hi - lo) / mid)


def w1_projection(n: float, rod_radius: float, rod_eps: float, supercell: int = 15, num_bands: int = 24,
                  k_points: int = 32, resolution: int = 32) -> tuple[np.ndarray, np.ndarray]:
    mp, mpb = _mpb()
    half = (supercell - 1) // 2
    geometry = [mp.Cylinder(radius=rod_radius, material=mp.Medium(epsilon=rod_eps), height=mp.inf,
                            axis=mp.Vector3(0, 0, 1), center=mp.Vector3(0, j))
                for j in range(-half, half + 1) if j != 0]
    kpts = mp.interpolate(k_points, [mp.Vector3(0.0), mp.Vector3(0.5)])
    ms = mpb.ModeSolver(geometry=geometry, geometry_lattice=mp.Lattice(size=mp.Vector3(1, supercell)),
                        default_material=mp.Medium(index=n), k_points=kpts, resolution=resolution, num_bands=num_bands)
    ms.run_tm()
    return np.array(ms.all_freqs), np.array([k.x for k in kpts])


def run(record: dict, structure: dict, params: dict) -> dict:
    b, la = params["bands"], structure["lattice"]
    n = float(params["analyte"]["n"])
    r, eps = float(la["rod_radius"]), float(la["rod_eps"])
    target = float(la["target_wavelength_nm"])
    task = b["task"]
    if os.environ.get("SPRAT_FAKE_FDTD"):
        return _fake(record, task, n, target)
    res, kp, nb = int(b["resolution"]), int(b["k_points"]), int(b["num_bands"])
    result: dict[str, Any] = dict(task=task, resolution=res, k_points=kp, num_bands=nb, rod_radius=r, rod_eps=eps)
    if task == "calibration":
        _, gap = bulk_bands(n, r, eps, "tm", nb, kp, res)
        result.update(gap=gap, n=n, a_nm=(gap["mid"] * target if gap["exists"] else None), f_mid=gap.get("mid"),
                      target_wavelength_nm=target, rod_radius=r)
    elif task == "bulk":
        for pol in ("tm", "te"):
            freqs, gap = bulk_bands(n, r, eps, pol, nb, kp, res)
            result[pol] = dict(gap=gap, freqs=freqs.tolist())
        g = result["tm"]["gap"]
        result["gap"] = g
        result["n"] = n
        if g["exists"]:
            result["a_nm"] = g["mid"] * target
            result["target_wavelength_nm"] = target
            result["rod_radius_nm"] = r * result["a_nm"]
        else:
            result["a_nm"] = None
    elif task == "sweep":
        rows = []
        for na in np.arange(float(b["sweep_start"]), float(b["sweep_stop"]) + 1e-9, float(b["sweep_step"])):
            _, gap = bulk_bands(float(na), r, eps, "tm", nb, kp, res)
            rows.append(dict(n=round(float(na), 6), gap=gap))
        result["sweep"] = rows
    else:
        freqs, kx = w1_projection(n, r, eps, int(b["supercell"]), nb, max(kp, 32), res)
        _, gap = bulk_bands(n, r, eps, "tm", 8, kp, res)
        inside = []
        if gap["exists"]:
            inside = [dict(band=int(j), f_min=float(freqs[:, j].min()), f_max=float(freqs[:, j].max()))
                      for j in range(freqs.shape[1]) if gap["lower"] < freqs[:, j].mean() < gap["upper"]]
        result.update(kx_over_2pi=kx.tolist(), freqs=freqs.tolist(), gap=gap, guided_bands=inside, n=n,
                      supercell=int(b["supercell"]))
    record["result"] = result
    try:
        import meep as mp
        record["software"]["meep"] = getattr(mp, "__version__", "")
    except ImportError:
        pass
    return record


def _fake(record: dict, task: str, n: float, target: float) -> dict:
    from .records import gap_edges
    lo, hi = gap_edges(n)
    gap = dict(exists=True, lower=lo, upper=hi, mid=0.5 * (lo + hi), relative_width_pct=100 * (hi - lo) / (0.5 * (lo + hi)))
    result: dict[str, Any] = dict(task=task, n=n, gap=gap)
    if task in ("bulk", "calibration"):
        result.update(a_nm=gap["mid"] * target, target_wavelength_nm=target, tm=dict(gap=gap, freqs=[]),
                      te=dict(gap=dict(exists=False), freqs=[]))
    elif task == "sweep":
        result["sweep"] = [dict(n=x, gap=dict(exists=True, lower=gap_edges(x)[0], upper=gap_edges(x)[1],
                                              mid=0.5 * sum(gap_edges(x)),
                                              relative_width_pct=100 * (gap_edges(x)[1] - gap_edges(x)[0]) / (0.5 * sum(gap_edges(x)))))
                           for x in np.round(np.arange(1.0, 1.501, 0.05), 6)]
    else:
        result.update(kx_over_2pi=[], freqs=[], guided_bands=[])
    record["result"] = result
    record["software"]["meep"] = "fake"
    return record


__all__ = ["bulk_bands", "gap_between", "w1_projection", "run"]
