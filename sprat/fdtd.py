"""The four FDTD run modes (Meep), one record per call.

* ``reference``: cavity-less transmission, the normalisation of the spectra.
* ``spectrum``:  transmission P_out(f) with the input flux plane as well.
* ``harminv``:   complex eigenfrequency by harmonic inversion of E_z at the probe point.
* ``field``:     DFT map of E_z at ``f_res`` and the analyte energy fraction.

Physics and numerics follow the FDTD script of the 2026 campaign exactly: TM polarisation
(E_z), Gaussian source with cutoff 5, subpixel smoothing, mirror symmetry in x when the
defect is on the axis, ``stop_when_dft_decayed`` with hard ceilings tied to the expected Q,
and Q_lim = pi f_cen t, whose ratio to Q is the margin of the sampling rule (harmonic inversion resolves decays
longer than the signal, so a margin of at least one is conservative).

Meep is imported inside the functions, so this module can be imported without it.  Setting
the environment variable ``SPRAT_FAKE_FDTD=1`` replaces the engine by a synthetic result
(used by the runner tests on machines without Meep).
"""

from __future__ import annotations

import math
import os
import time
from typing import Any

import numpy as np

from . import geometry as geo
from .cost import estimate
from .params import q_limit


def _meep():
    try:
        import meep as mp
    except ImportError as e:                        # pragma: no cover - depends on the machine
        raise ImportError("Meep is not installed. Install pymeep from conda-forge (see README, Installation) "
                          "or set SPRAT_FAKE_FDTD=1 for a synthetic test run.") from e
    return mp


def meep_version() -> str:
    try:
        import meep as mp
        return getattr(mp, "__version__", "")
    except ImportError:
        return ""


# --------------------------------------------------------------------------- Meep objects
def meep_geometry(structure: dict, params: dict, g: geo.Geometry):
    mp = _meep()
    si = mp.Medium(epsilon=float(structure["lattice"]["rod_eps"]))
    return [mp.Cylinder(radius=r, material=si, height=mp.inf, axis=mp.Vector3(0, 0, 1), center=mp.Vector3(x, y))
            for x, y, r in g.rods]


def boundary_layers(params: dict, g: geo.Geometry):
    mp = _meep()
    if g.absorber > 0:
        return [mp.Absorber(g.absorber, direction=mp.X), mp.PML(g.pml, direction=mp.Y)]
    return [mp.PML(g.pml)]


def analyte_medium(params: dict):
    mp = _meep()
    n, k = float(params["analyte"]["n"]), float(params["analyte"]["k"])
    if k > 0:
        eps = n * n
        sigma_d = 2 * math.pi * float(params["source"]["fcen"]) * 2 * n * k / eps
        return mp.Medium(epsilon=eps, D_conductivity=sigma_d)
    return mp.Medium(index=n)


def symmetries(structure: dict, params: dict, g: geo.Geometry):
    mp = _meep()
    sector, phase = geo.symmetry_sector(structure, params, g)
    return ([] if sector == "none" else [mp.Mirror(mp.X, phase=phase)]), sector


def _simulation(structure: dict, params: dict, g: geo.Geometry, sources, sym):
    mp = _meep()
    return mp.Simulation(cell_size=mp.Vector3(g.sx, g.sy, 0), geometry=meep_geometry(structure, params, g),
                         sources=sources, boundary_layers=boundary_layers(params, g),
                         default_material=analyte_medium(params), resolution=int(params["numerics"]["resolution"]),
                         eps_averaging=bool(params["numerics"]["subpixel"]), symmetries=sym,
                         Courant=float(params["numerics"]["courant"]), force_complex_fields=False)


def _pixels(g: geo.Geometry, res: int) -> int:
    return int(round(g.sx * res)) * int(round(g.sy * res))


# --------------------------------------------------------------------------- modes
def run_spectrum(structure: dict, params: dict, prov: dict) -> dict[str, Any]:
    """Transmission (mode spectrum) or the cavity-less reference; no mirror symmetry (the source is at one end)."""
    mp = _meep()
    g = geo.build(structure, params)
    src_x, in_x, out_x = geo.source_and_flux_x(g)
    p = params
    fcen, df = float(p["source"]["fcen"]), float(p["source"]["fwidth"])
    src = [mp.Source(mp.GaussianSource(fcen, fwidth=df, cutoff=float(p["source"]["cutoff"])), component=mp.Ez,
                     center=mp.Vector3(src_x, 0), size=mp.Vector3(0, float(p["spectrum"]["source_width"])))]
    sim = _simulation(structure, p, g, src, [])
    w = float(p["spectrum"]["flux_width"])
    flux_in = sim.add_flux(fcen, df, int(p["spectrum"]["nfreq"]), mp.FluxRegion(center=mp.Vector3(in_x, 0), size=mp.Vector3(0, w)))
    flux_out = sim.add_flux(fcen, df, int(p["spectrum"]["nfreq"]), mp.FluxRegion(center=mp.Vector3(out_x, 0), size=mp.Vector3(0, w)))
    q = max(float(p["run"]["q_est"]), 50.0)
    t_max = min(float(p["spectrum"]["t_max"]), 10.0 * q / (math.pi * fcen))
    t_min = max(4.0 * q / (math.pi * fcen) * 0.25, 200.0)
    sim.run(until_after_sources=mp.stop_when_dft_decayed(tol=float(p["spectrum"]["dft_tol"]),
                                                         minimum_run_time=t_min, maximum_run_time=t_max))
    prov.update(steps=int(sim.timestep()), sim_time=float(sim.round_time()),
                pixels=_pixels(g, int(p["numerics"]["resolution"])), symmetry="none")
    return dict(f=[float(x) for x in mp.get_flux_freqs(flux_out)], flux_in=[float(x) for x in mp.get_fluxes(flux_in)],
                flux_out=[float(x) for x in mp.get_fluxes(flux_out)], t_end=float(sim.round_time()),
                source_x=src_x, flux_in_x=in_x, flux_out_x=out_x)


def _harminv_once(structure: dict, params: dict, t_run: float):
    mp = _meep()
    g = geo.build(structure, params)
    probe = mp.Vector3(*geo.probe_point(structure, params))
    sym, sector = symmetries(structure, params, g)
    p = params
    fcen, df = float(p["source"]["fcen"]), float(p["source"]["fwidth"])
    src = [mp.Source(mp.GaussianSource(fcen, fwidth=df, cutoff=float(p["source"]["cutoff"])), component=mp.Ez, center=probe)]
    sim = _simulation(structure, p, g, src, sym)
    h = mp.Harminv(mp.Ez, probe, fcen, df)
    sim.run(mp.after_sources(h), until_after_sources=float(t_run))
    modes = [dict(f=float(m.freq), Q=float(abs(m.Q)), amplitude=float(abs(m.amp)), error=float(abs(m.err)))
             for m in h.modes if m.freq > 0 and abs(m.Q) > 5]
    modes.sort(key=lambda m: -m["amplitude"])
    prov = dict(steps=int(sim.timestep()), sim_time=float(sim.round_time()),
                pixels=_pixels(g, int(p["numerics"]["resolution"])), symmetry=sector)
    return modes, prov


def run_harminv(structure: dict, params: dict, prov: dict) -> dict[str, Any]:
    p = params
    fcen = float(p["source"]["fcen"])
    t_run = float(p["harminv"]["t"])
    modes, o = _harminv_once(structure, p, t_run)
    history = [dict(t=t_run, Q=(modes[0]["Q"] if modes else None))]
    if p["harminv"]["auto_t"] and modes:
        q = modes[0]["Q"]
        t_need = float(p["harminv"]["auto_factor"]) * q / (math.pi * fcen)
        if t_need > 1.5 * t_run:
            t_run = float(min(t_need, float(p["harminv"]["t_max"])))
            modes, o = _harminv_once(structure, p, t_run)
            history.append(dict(t=t_run, Q=(modes[0]["Q"] if modes else None)))
    prov.update(o)
    return dict(modes=modes, t_used=t_run, auto_history=history, Q_limit=q_limit(fcen, t_run))


def run_field(structure: dict, params: dict, prov: dict, records_dir: str, label: str) -> dict[str, Any]:
    mp = _meep()
    g = geo.build(structure, params)
    probe = mp.Vector3(*geo.probe_point(structure, params))
    sym, sector = symmetries(structure, params, g)
    p = params
    fcen, df = float(p["source"]["fcen"]), float(p["source"]["fwidth"])
    f_res = float(p["field"]["f_res"])
    if f_res <= 0:
        raise ValueError("[field] f_res must be given (run harminv first)")
    src = [mp.Source(mp.GaussianSource(fcen, fwidth=df, cutoff=float(p["source"]["cutoff"])), component=mp.Ez, center=probe)]
    sim = _simulation(structure, p, g, src, sym)
    Lx, Ly = g.interior
    vol = mp.Volume(center=mp.Vector3(), size=mp.Vector3(Lx, Ly))
    try:
        dft = sim.add_dft_fields([mp.Ez], [f_res], where=vol)
    except TypeError:                                   # older Meep signature
        dft = sim.add_dft_fields([mp.Ez], f_res, f_res, 1, where=vol)
    q = max(float(p["run"]["q_est"]), 50.0)
    t_max = min(float(p["field"]["t_max"]), 8.0 * q / (math.pi * f_res))
    sim.run(until_after_sources=mp.stop_when_dft_decayed(tol=float(p["spectrum"]["dft_tol"]),
                                                         minimum_run_time=200.0, maximum_run_time=t_max))
    ez = sim.get_dft_array(dft, mp.Ez, 0)
    eps = sim.get_array(vol=vol, component=mp.Dielectric)
    if ez.shape != eps.shape:
        n0, n1 = min(ez.shape[0], eps.shape[0]), min(ez.shape[1], eps.shape[1])
        ez, eps = ez[:n0, :n1], eps[:n0, :n1]
    result = analyse_field(ez, eps, float(p["analyte"]["n"]), float(structure["lattice"]["rod_eps"]),
                           float(structure["lattice"]["a_nm"]), f_res)
    field_file = os.path.join(records_dir, label + ".npz")
    np.savez_compressed(field_file, Ez=ez.astype(np.complex64), eps=eps.astype(np.float32), sx=g.sx, sy=g.sy,
                        Lx=Lx, Ly=Ly, pml=g.pml, edge_x=g.edge_x, resolution=int(p["numerics"]["resolution"]),
                        f_res=f_res, wavelength_nm=result["wavelength_nm"], n=float(p["analyte"]["n"]),
                        rod_eps=float(structure["lattice"]["rod_eps"]))
    prov.update(steps=int(sim.timestep()), sim_time=float(sim.round_time()),
                pixels=_pixels(g, int(p["numerics"]["resolution"])), symmetry=sector)
    result["field_file"] = os.path.basename(field_file)
    return result


def analyse_field(ez: np.ndarray, eps: np.ndarray, n: float, rod_eps: float, a_nm: float, f_res: float) -> dict[str, Any]:
    """The masked analyte energy fraction of the campaign and the first-order sensitivity it implies."""
    E2 = np.abs(np.asarray(ez, dtype=np.complex128)) ** 2
    eps = np.asarray(eps, dtype=np.float64)
    threshold = 0.5 * (n * n + rod_eps)
    mask = eps < threshold
    num = float((eps * E2)[mask].sum())
    den = float((eps * E2).sum())
    frac = num / den if den > 0 else float("nan")
    lam = a_nm / f_res
    return dict(analyte_energy_fraction=frac, wavelength_nm=float(lam), S_first_order_nm_per_RIU=float(lam / n * frac))


# --------------------------------------------------------------------------- driver
def run(record: dict, structure: dict, params: dict, records_dir: str = ".") -> dict:
    mode = params["run"]["mode"]
    prov = record["provenance"]
    label = record["task"]["label"]
    if os.environ.get("SPRAT_FAKE_FDTD"):
        return _fake(record, structure, params, records_dir)
    mp = _meep()
    if params["run"].get("quiet", True):
        mp.verbosity(0)
    record["software"]["meep"] = meep_version()
    if mode in ("spectrum", "reference"):
        result = run_spectrum(structure, params, prov)
    elif mode == "harminv":
        result = run_harminv(structure, params, prov)
    elif mode == "field":
        result = run_field(structure, params, prov, records_dir, label)
    else:
        raise ValueError(mode)
    if "modes" in result:
        a_nm = float(structure["lattice"]["a_nm"])
        for m in result["modes"]:
            m["wavelength_nm"] = a_nm / m["f"]
            m["fwhm_nm"] = m["wavelength_nm"] / m["Q"]
    record["result"] = result
    return record


def _fake(record: dict, structure: dict, params: dict, records_dir: str) -> dict:
    """Synthetic result for tests: the cost model's steps and pixels, a mode at fcen with Q = q_est."""
    mode = params["run"]["mode"]
    est = estimate(structure, params)
    g = geo.build(structure, params)
    res = int(params["numerics"]["resolution"])
    sector, _ = geo.symmetry_sector(structure, params, g)
    steps = int(est["total_time"] * res / float(params["numerics"]["courant"]))
    time.sleep(float(os.environ.get("SPRAT_FAKE_SLEEP", "0.05")))
    record["software"]["meep"] = "fake"
    record["provenance"].update(steps=steps, sim_time=est["total_time"], pixels=_pixels(g, res),
                                symmetry=(sector if mode in ("harminv", "field") else "none"))
    fcen = float(params["source"]["fcen"])
    a_nm = float(structure["lattice"]["a_nm"])
    if mode == "harminv":
        q = float(params["run"]["q_est"])
        t = float(params["harminv"]["t"])
        record["result"] = dict(modes=[dict(f=fcen, Q=q, amplitude=1.0, error=1e-9, wavelength_nm=a_nm / fcen,
                                            fwhm_nm=a_nm / fcen / q)], t_used=t, auto_history=[dict(t=t, Q=q)],
                                Q_limit=q_limit(fcen, t))
    elif mode in ("spectrum", "reference"):
        nf = int(params["spectrum"]["nfreq"])
        f = np.linspace(fcen - params["source"]["fwidth"] / 2, fcen + params["source"]["fwidth"] / 2, nf)
        base = np.exp(-((f - fcen) / (params["source"]["fwidth"] / 2)) ** 2)
        if mode == "spectrum":
            q = float(params["run"]["q_est"])
            e = 2 * (f - fcen) * q / fcen
            base = base * (e ** 2 / (1 + e ** 2))
        record["result"] = dict(f=f.tolist(), flux_in=(base * 1.02).tolist(), flux_out=base.tolist(),
                                t_end=est["total_time"])
    elif mode == "field":
        Lx, Ly = g.interior
        nx, ny = int(round(Lx * res)), int(round(Ly * res))
        eps = np.full((nx, ny), float(params["analyte"]["n"]) ** 2, dtype=np.float32)
        ez = np.ones((nx, ny), dtype=np.complex64)
        f_res = float(params["field"]["f_res"])
        result = analyse_field(ez, eps, float(params["analyte"]["n"]), float(structure["lattice"]["rod_eps"]), a_nm, f_res)
        path = os.path.join(records_dir, record["task"]["label"] + ".npz")
        os.makedirs(records_dir, exist_ok=True)
        np.savez_compressed(path, Ez=ez, eps=eps, sx=g.sx, sy=g.sy, Lx=Lx, Ly=Ly, pml=g.pml, edge_x=g.edge_x,
                            resolution=res, f_res=f_res, wavelength_nm=result["wavelength_nm"],
                            n=float(params["analyte"]["n"]), rod_eps=float(structure["lattice"]["rod_eps"]))
        result["field_file"] = os.path.basename(path)
        record["result"] = result
    return record


__all__ = ["run", "run_spectrum", "run_harminv", "run_field", "analyse_field", "meep_version"]
