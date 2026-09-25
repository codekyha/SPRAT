"""Cost model: pixel-steps, memory, calibration and wall-time estimates.

With Courant factor C the time step is dt = C / r, so a run of T time units on an
s_x x s_y cell at resolution r costs

    W = N_x N_y N_t = (s_x r)(s_y r) T / dt = s_x s_y T r^3 / C   pixel-steps,

which is 2 s_x s_y T r^3 for C = 0.5.  Mirror symmetry halves it.  The reference
computation of the paper (30a x 29a, r = 24, t_harminv = 17 750) is 2.16e11 pixel-steps and
took 677 s on one core of an AMD EPYC 7742 node running 32 tasks; throughput on another
machine is measured with ``sprat calibrate``.
"""

from __future__ import annotations

import json
import math
import os
import platform
import time
from typing import Any

from . import geometry as geo
from .params import source_time

REFERENCE_PIXEL_STEPS = 2.16e11        # one reference-run unit (RRU)
REFERENCE_THROUGHPUT = 3.19e8          # pixel-steps per second per task on the Altay node (32 tasks per node)
BYTES_PER_PIXEL = 64.0                 # 2D TM: six field arrays plus two permittivity arrays


def pixel_steps(sx: float, sy: float, total_time: float, resolution: int, courant: float = 0.5,
                symmetric: bool = False) -> float:
    w = sx * sy * float(total_time) * float(resolution) ** 3 / float(courant)
    return 0.5 * w if symmetric else w


def memory_mb(sx: float, sy: float, resolution: int, nfreq: int = 0, flux_length: float = 3.0) -> float:
    fields = BYTES_PER_PIXEL * sx * sy * resolution ** 2
    dft = 2 * 2 * (flux_length * resolution) * nfreq * 16.0        # two planes, two components, complex doubles
    return (fields + dft) / 1048576.0


def run_time(structure: dict, p: dict) -> float:
    """Total simulated time (source + measurement) used by the cost estimate of one task."""
    mode = p["run"]["mode"]
    fcen = float(p["source"]["fcen"])
    ts = source_time(p["source"]["fwidth"], p["source"]["cutoff"])
    q = max(float(p["run"]["q_est"]), 50.0)
    if mode == "harminv":
        t = float(p["harminv"]["t"])
        if p["harminv"]["auto_t"]:
            t = max(t, min(p["harminv"]["auto_factor"] * q / (math.pi * fcen), p["harminv"]["t_max"]))
        return ts + t
    if mode == "spectrum":
        return ts + min(p["spectrum"]["t_max"], 6.0 * q / (math.pi * fcen))
    if mode == "reference":
        return ts + 600.0
    if mode == "field":
        f = float(p["field"]["f_res"]) or fcen
        return ts + min(p["field"]["t_max"], 4.0 * q / (math.pi * f))
    return 0.0


def estimate(structure: dict, p: dict) -> dict[str, Any]:
    """Pixel-steps, memory and RRU for one resolved task."""
    mode = p["run"]["mode"]
    if mode in ("bands", "pwe"):
        return dict(pixel_steps=0.0, memory_mb=300.0, rru=0.0, symmetric=False, sx=0.0, sy=0.0,
                    total_time=0.0, fixed_seconds=(600.0 if mode == "pwe" and not p["pwe"]["quick"] else 120.0))
    g = geo.build(structure, p)
    sector, _ = geo.symmetry_sector(structure, p, g)
    symmetric = sector != "none" and mode in ("harminv", "field")
    T = run_time(structure, p)
    res = int(p["numerics"]["resolution"])
    w = pixel_steps(g.sx, g.sy, T, res, p["numerics"]["courant"], symmetric)
    nfreq = int(p["spectrum"]["nfreq"]) if mode in ("spectrum", "reference") else 0
    return dict(pixel_steps=w, memory_mb=memory_mb(g.sx, g.sy, res, nfreq), rru=w / REFERENCE_PIXEL_STEPS,
                symmetric=symmetric, sx=g.sx, sy=g.sy, total_time=T, fixed_seconds=0.0)


def seconds(est: dict[str, Any], throughput: float) -> float:
    return est["fixed_seconds"] + est["pixel_steps"] / max(throughput, 1.0)


# --------------------------------------------------------------------------- calibration
CALIBRATION_FILE = "calibration.json"


def read_calibration(directory: str = ".") -> dict[str, Any]:
    for cand in (os.path.join(directory, CALIBRATION_FILE), CALIBRATION_FILE):
        if os.path.exists(cand):
            with open(cand) as fh:
                d = json.load(fh)
            d.setdefault("source", cand)
            return d
    return dict(throughput=REFERENCE_THROUGHPUT, source="built-in Altay reference (run `sprat calibrate`)",
                host="", cores=0, measured=False)


def write_calibration(directory: str, throughput: float, samples: list[dict[str, Any]]) -> dict[str, Any]:
    d = dict(throughput=float(throughput), measured=True, host=platform.node(), machine=platform.machine(),
             cores=os.cpu_count(), date=time.strftime("%Y-%m-%dT%H:%M:%S"), samples=samples,
             unit="pixel-steps per second per single-core task")
    os.makedirs(directory, exist_ok=True)
    with open(os.path.join(directory, CALIBRATION_FILE), "w") as fh:
        json.dump(d, fh, indent=1)
    return d


def makespan(durations: list[float], workers: int) -> float:
    """Longest-processing-time-first list scheduling: the wall time of the batch."""
    import heapq
    if not durations:
        return 0.0
    workers = max(1, int(workers))
    heap = [0.0] * workers
    heapq.heapify(heap)
    for t in sorted(durations, reverse=True):
        m = heapq.heappop(heap)
        heapq.heappush(heap, m + t)
    return max(heap)


__all__ = ["REFERENCE_PIXEL_STEPS", "REFERENCE_THROUGHPUT", "pixel_steps", "memory_mb", "run_time", "estimate",
           "seconds", "read_calibration", "write_calibration", "makespan", "CALIBRATION_FILE"]
