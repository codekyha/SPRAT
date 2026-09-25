"""Environment check, throughput calibration and the Meep self-test."""

from __future__ import annotations

import os
import platform
import shutil
import time

import psutil

from . import __version__, cost
from .params import parse_params_text, resolve
from .structure import default_structure

REFERENCE_Q = 6926.72             # the archived reference record (r_d = 0.060a, N_sep = 4, n_cl = 12, res 24, t = 17 750)
REFERENCE_F = 0.304614
REFERENCE_T = 17750.0
REGRESSION_TOL_Q = 1e-3           # |Q/Q_ref - 1|
REGRESSION_TOL_F = 1e-5           # |f_r - f_ref|


def check(records_dir: str = "records") -> int:
    ok = True
    print(f"sprat {__version__}  python {platform.python_version()}  {platform.system()} {platform.machine()}")
    for mod in ("numpy", "scipy", "psutil", "matplotlib", "pandas"):
        try:
            m = __import__(mod)
            print(f"  {mod:12s} {getattr(m, '__version__', '?')}")
        except ImportError:
            print(f"  {mod:12s} missing" + ("" if mod in ("matplotlib", "pandas") else "  <-- required"))
            ok &= mod in ("matplotlib", "pandas")
    try:
        import meep as mp
        print(f"  {'meep':12s} {getattr(mp, '__version__', '?')}  (FDTD and MPB available)")
        try:
            from meep import mpb  # noqa: F401  (probe only)
            del mpb
        except ImportError:
            print("  MPB module missing: band structures unavailable")
    except ImportError:
        print("  meep         missing: FDTD/MPB modes unavailable; analysis, figures and the plane-wave layer still work")
        if platform.system() == "Windows":
            print("  (native Windows has no pymeep; use WSL2, see README)")
    phys, logical = psutil.cpu_count(logical=False), psutil.cpu_count(logical=True)
    vm = psutil.virtual_memory()
    print(f"  cores        {phys} physical, {logical} logical")
    print(f"  memory       {vm.total / 1073741824:.1f} GB total, {vm.available / 1073741824:.1f} GB available")
    du = shutil.disk_usage(os.path.abspath(records_dir) if os.path.isdir(records_dir) else ".")
    print(f"  disk         {du.free / 1073741824:.1f} GB free")
    try:
        os.sched_getaffinity(0)
        print("  pinning      available (sched_setaffinity)")
    except AttributeError:
        print("  pinning      not available on this platform (tasks run unpinned)")
    cal = cost.read_calibration(records_dir)
    print(f"  calibration  {cal['throughput']:.3e} pixel-steps/s per task ({cal.get('source')})")
    print("environment OK" if ok else "environment INCOMPLETE")
    return 0 if ok else 1


def _reference_task(resolution: int, t: float, records_dir: str, label: str):
    s = default_structure()
    p = parse_params_text(f"[run]\nmode = harminv\nq_est = 6927\nlabel = {label}\n[numerics]\nresolution = {resolution}\n"
                          f"[harminv]\nt = {t}\n")
    rp = resolve(s, p)
    return s, rp


def calibrate(records_dir: str = "records", resolutions=(12, 24), t: float = 2000.0) -> dict:
    """Run the reference geometry at the given resolutions and measure pixel-steps per second."""
    from .execute import execute_task
    samples = []
    for res in resolutions:
        label = f"calibration_res{res}_t{t:g}"
        s, rp = _reference_task(res, t, records_dir, label)
        est = cost.estimate(s, rp)
        task = dict(id=0, label=label, mode="harminv", structure=s, params=rp, overrides={}, record=os.path.join(records_dir, label + ".json"),
                    estimate=est, est_seconds=0.0, created=time.strftime("%FT%T"))
        rec = execute_task(task, overwrite=True)
        o = rec["provenance"]
        samples.append(dict(resolution=res, t=t, wall_s=o["wall_s"], pixel_steps=o.get("pixel_steps"), throughput=o.get("throughput"),
                            steps=o.get("steps"), pixels=o.get("pixels")))
        print(f"  res {res}: {o['wall_s']:.1f} s, {o.get('throughput', 0):.3e} pixel-steps/s")
    thr = max(x["throughput"] for x in samples if x["throughput"]) if samples else cost.REFERENCE_THROUGHPUT
    # the highest resolution is the production-like sample; use it
    prod = [x for x in samples if x["resolution"] == max(resolutions)]
    if prod and prod[0]["throughput"]:
        thr = prod[0]["throughput"]
    d = cost.write_calibration(records_dir, thr, samples)
    print(f"calibration written: {thr:.3e} pixel-steps/s per single-core task -> {os.path.join(records_dir, cost.CALIBRATION_FILE)}")
    return d


def selftest(records_dir: str = "selftest_records", reference: bool = False, workers: int = 2) -> int:
    """Smoke test of the four FDTD modes and the bands at resolution 8; --reference adds the regression run."""
    from .execute import execute_task
    from .records import strongest_mode
    s = default_structure()
    ok = True
    cases = [
        ("harminv", "[run]\nmode = harminv\nq_est = 500\nlabel = selftest_harminv\n[numerics]\nresolution = 8\n[harminv]\nt = 600\n"),
        ("reference", "[run]\nmode = reference\nlabel = selftest_reference\n[numerics]\nresolution = 8\n[spectrum]\nnfreq = 51\nt_max = 400\n"),
        ("spectrum", "[run]\nmode = spectrum\nq_est = 500\nlabel = selftest_spectrum\n[numerics]\nresolution = 8\n[spectrum]\nnfreq = 51\nt_max = 400\n"),
        ("field", "[run]\nmode = field\nq_est = 500\nlabel = selftest_field\n[numerics]\nresolution = 8\n[field]\nf_res = 0.3046\nt_max = 300\n"),
        ("bands", "[run]\nmode = bands\nlabel = selftest_bands\n[bands]\nresolution = 16\nk_points = 6\n"),
    ]
    for mode, text in cases:
        rp = resolve(s, parse_params_text(text))
        task = dict(id=0, label=rp["run"]["label"], mode=mode, structure=s, params=rp, overrides={},
                    record=os.path.join(records_dir, rp["run"]["label"] + ".json"), estimate=cost.estimate(s, rp),
                    est_seconds=0.0, created=time.strftime("%FT%T"))
        t0 = time.time()
        try:
            rec = execute_task(task, overwrite=True)
            extra = ""
            m = strongest_mode(rec)
            if m:
                extra = f"Q = {m['Q']:.1f} at f = {m['f']:.5f}"
            elif mode == "bands":
                extra = f"TM gap {rec['result']['gap']}"
            print(f"  {mode:10s} OK in {time.time() - t0:.0f} s  {extra}")
        except Exception as e:                  # noqa: BLE001
            ok = False
            print(f"  {mode:10s} FAILED: {type(e).__name__}: {e}")
    if reference:
        print("reference regression (about 11 minutes on one core at the Altay throughput) ...")
        s2, rp = _reference_task(24, REFERENCE_T, records_dir, "selftest_reference_regression")
        task = dict(id=0, label=rp["run"]["label"], mode="harminv", structure=s2, params=rp, overrides={},
                    record=os.path.join(records_dir, rp["run"]["label"] + ".json"), estimate=cost.estimate(s2, rp),
                    est_seconds=0.0, created=time.strftime("%FT%T"))
        rec = execute_task(task, overwrite=True)
        m = strongest_mode(rec)
        dq = abs(m["Q"] / REFERENCE_Q - 1)
        df = abs(m["f"] - REFERENCE_F)
        passed = dq < REGRESSION_TOL_Q and df < REGRESSION_TOL_F
        ok &= passed
        print(f"  regression: Q = {m['Q']:.2f} (reference {REFERENCE_Q}, |dQ/Q| = {dq:.2e}), f = {m['f']:.6f} "
              f"(reference {REFERENCE_F}, |df| = {df:.2e}) -> {'PASS' if passed else 'FAIL'}")
    print("selftest PASS" if ok else "selftest FAIL")
    return 0 if ok else 1


__all__ = ["check", "calibrate", "selftest", "REFERENCE_Q", "REFERENCE_F", "REFERENCE_T"]
