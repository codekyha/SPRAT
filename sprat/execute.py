"""Execute one task (the worker side of the runner).

The heavy engines are imported lazily so that the planner, the runner and the analysis
layer work on a machine without Meep.
"""

from __future__ import annotations

import json
import os
import sys
import time
from typing import Any

from . import records
from .structure import structure_from_dict


def execute_task(task: dict[str, Any], overwrite: bool = False) -> dict[str, Any]:
    """Run the task and write its record; return the record."""
    structure = structure_from_dict(task["structure"])
    params = task["params"]
    mode = params["run"]["mode"]
    path = task["record"]
    if os.path.exists(path) and not (overwrite or params["run"].get("overwrite")):
        old = records.load(path)
        probe = records.new_record(structure, params)
        if records.same_task(old, probe):
            print(f"SKIPPED (record exists): {path}")
            return old
        raise FileExistsError(f"{path} exists with different parameters (label collision); "
                              "choose another label or pass --overwrite")
    rec = records.new_record(structure, params, task=dict(
        label=task["label"], tag=params["run"].get("tag", ""), created=task.get("created", ""),
        structure_file=task.get("structure_file", ""), params_file=task.get("params_file", ""),
        overrides=task.get("overrides", {}), source="sprat"))
    t0 = time.time()
    if mode in ("harminv", "spectrum", "reference", "field"):
        from . import fdtd
        rec = fdtd.run(rec, structure, params, records_dir=os.path.dirname(path) or ".")
    elif mode == "bands":
        from . import bands
        rec = bands.run(rec, structure, params)
    elif mode == "pwe":
        from . import pwe
        rec = pwe.run(rec, structure, params)
    else:
        raise ValueError(f"unknown mode {mode!r}")
    rec["provenance"]["wall_s"] = time.time() - t0
    o = rec["provenance"]
    if o.get("steps") and o.get("pixels"):
        effective = o["pixels"] * (0.5 if o.get("symmetry") in ("even", "odd") else 1.0)
        o["pixel_steps"] = effective * o["steps"]
        o["throughput"] = o["pixel_steps"] / max(o["wall_s"], 1e-9)
    records.save(rec, path)
    print(f"WRITTEN: {path}  [{o['wall_s']:.1f} s, {o.get('throughput', 0):.2e} pixel-steps/s]")
    m = records.strongest_mode(rec)
    if m:
        warn = "  <-- WARNING: Q above the harmonic-inversion ceiling" if m["Q"] > rec["result"].get("Q_limit", 1e99) else ""
        print(f"  wavelength = {m['wavelength_nm']:.3f} nm   Q = {m['Q']:.1f}   fwhm = {m['fwhm_nm']:.4f} nm{warn}")
    return rec


def task_from_record(record: dict[str, Any], out: str | None = None) -> dict[str, Any]:
    """A record carries its own structure and parameters, so it can be re-run from itself; the new
    record goes to ``out`` (default: ``<label>_rerun.json`` beside the original)."""
    label = record["task"]["label"]
    src = record.get("_file") or ""
    return dict(id=0, label=label, structure=record["structure"], params=record["params"],
                record=out or os.path.join(os.path.dirname(src) or ".", label + "_rerun.json"),
                structure_file=record["task"].get("structure_file", ""), params_file=record["task"].get("params_file", ""),
                overrides=record["task"].get("overrides", {}), created=time.strftime("%Y-%m-%dT%H:%M:%S"))


def run_task_file(path: str, overwrite: bool = False, out: str | None = None) -> int:
    """Run one task file (from ``sprat plan``) or re-run one record from itself (a ``sprat-record-*`` JSON)."""
    with open(path, encoding="utf-8") as fh:
        task = json.load(fh)
    if str(task.get("schema", "")).startswith("sprat-record"):
        task["_file"] = path
        task = task_from_record(task, out)
    elif out:
        task["record"] = out
    try:
        execute_task(task, overwrite=overwrite)
    except Exception as e:                      # noqa: BLE001 - the runner records the failure
        print(f"FAILED: {task.get('label')}: {type(e).__name__}: {e}", file=sys.stderr)
        raise
    return 0


__all__ = ["execute_task", "run_task_file", "task_from_record"]
