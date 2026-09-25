"""From the two text files to a task list with cost estimates.

``sprat plan structure.phc params.par -o tasks`` writes ``tasks.jsonl`` (one fully resolved
task per line, the input of ``sprat run``), ``tasks.tsv`` (a readable summary) and
``tasks.summary.md`` (counts, pixel-steps, hours at the calibrated throughput, the longest
task, memory per task, the recommended worker count).
"""

from __future__ import annotations

import json
import os
import time
from typing import Any

from . import cost
from .params import apply_override, expand_sweeps, load_params, parse_override, resolve
from .records import record_path
from .structure import load_structure, structure_to_json

TSV_COLUMNS = ["id", "label", "mode", "radius", "row", "dx", "dy", "n", "resolution", "cladding_rows",
               "guide_periods", "termination_x", "rows", "t_harminv", "q_est", "margin", "pixel_steps", "rru",
               "memory_mb", "est_hours"]


def build_tasks(structure_file: str, params_file: str, overrides: list[str] | None = None,
                records_dir: str = "records", throughput: float | None = None) -> list[dict[str, Any]]:
    structure = load_structure(structure_file)
    params = load_params(params_file)
    for text in overrides or []:
        k, v = parse_override(text)
        apply_override(structure, params, k, v)
    if structure["lattice"]["a_nm"] == "calibrate":
        structure["lattice"]["a_nm"] = calibrated_a_nm(records_dir)
    thr = throughput or cost.read_calibration(records_dir)["throughput"]
    tasks = []
    for i, (s, p, ov) in enumerate(expand_sweeps(structure, params), start=1):
        rp = resolve(s, p)
        est = cost.estimate(s, rp)
        label = rp["run"]["label"]
        tasks.append(dict(
            id=i, label=label, mode=rp["run"]["mode"], structure=structure_to_json(s), params=rp,
            overrides=ov, structure_file=os.path.abspath(structure_file), params_file=os.path.abspath(params_file),
            record=record_path(records_dir, label), estimate=est, est_seconds=cost.seconds(est, thr),
            created=time.strftime("%Y-%m-%dT%H:%M:%S")))
    _check_duplicates(tasks)
    return tasks


def _check_duplicates(tasks: list[dict]) -> None:
    seen: dict[str, dict] = {}
    for t in tasks:
        if t["label"] in seen:
            raise ValueError(f"two tasks share the label {t['label']!r} (ids {seen[t['label']]['id']} and {t['id']}); "
                             "add the varied key to [run] label or make the labels distinct with tag")
        seen[t["label"]] = t


def calibrated_a_nm(records_dir: str) -> float:
    """a_nm = f_mid * target from the most recent bulk band record in the records directory."""
    from .records import load_dir
    cands = [r for r in load_dir(records_dir, "bands_bulk_*.json")
             if r["result"].get("task") == "bulk" and r["result"].get("a_nm")]
    if not cands:
        raise FileNotFoundError(f"a_nm = calibrate: no bulk band record with a gap in {records_dir}; "
                                "run a [run] mode = bands task first or set a_nm explicitly")
    cands.sort(key=lambda r: r["provenance"].get("started", ""))
    return float(cands[-1]["result"]["a_nm"])


def write_tasks(tasks: list[dict], out_prefix: str, throughput: float, workers: int | None = None) -> dict[str, str]:
    os.makedirs(os.path.dirname(os.path.abspath(out_prefix)) or ".", exist_ok=True)
    jsonl = out_prefix + ".jsonl"
    with open(jsonl, "w", encoding="utf-8") as fh:
        for t in tasks:
            fh.write(json.dumps(t) + "\n")
    tsv = out_prefix + ".tsv"
    with open(tsv, "w", encoding="utf-8") as fh:
        fh.write("\t".join(TSV_COLUMNS) + "\n")
        for t in tasks:
            fh.write("\t".join(str(v) for v in _tsv_row(t)) + "\n")
    md = out_prefix + ".summary.md"
    with open(md, "w", encoding="utf-8") as fh:
        fh.write(summary(tasks, throughput, workers))
    return dict(jsonl=jsonl, tsv=tsv, summary=md)


def _tsv_row(t: dict) -> list:
    s, p, e = t["structure"], t["params"], t["estimate"]
    d, c = s["defect"], s["cell"]
    return [t["id"], t["label"], t["mode"], d["radius"], d["row"], d["dx"], d["dy"], p["analyte"]["n"],
            p["numerics"]["resolution"], c["cladding_rows"], c["guide_periods"], c["termination_x"],
            ";".join(f"{j}:{r}" for j, r in sorted(s["rows"].items())),
            p["harminv"]["t"] if t["mode"] == "harminv" else "", p["run"]["q_est"], p["harminv"]["margin"],
            f"{e['pixel_steps']:.3e}", f"{e['rru']:.2f}", f"{e['memory_mb']:.0f}", f"{t['est_seconds'] / 3600:.2f}"]


def summary(tasks: list[dict], throughput: float, workers: int | None = None) -> str:
    import psutil
    n = len(tasks)
    total_px = sum(t["estimate"]["pixel_steps"] for t in tasks)
    secs = [t["est_seconds"] for t in tasks]
    total_h = sum(secs) / 3600
    longest = max(tasks, key=lambda t: t["est_seconds"]) if tasks else None
    phys = psutil.cpu_count(logical=False) or os.cpu_count() or 1
    workers = workers or phys
    mem_task = max((t["estimate"]["memory_mb"] for t in tasks), default=0.0) + 300.0
    avail_gb = psutil.virtual_memory().available / 1073741824
    max_by_mem = max(1, int(avail_gb * 1024 / mem_task))
    lines = ["# Task plan", "", f"generated {time.strftime('%Y-%m-%d %H:%M:%S')}",
             f"throughput used for the estimate: {throughput:.3e} pixel-steps/s per task", "",
             f"- tasks: {n}", f"- total: {total_px:.3e} pixel-steps = {total_px / cost.REFERENCE_PIXEL_STEPS:.1f} RRU "
             f"= {total_h:.1f} single-core hours",
             (f"- longest task: {longest['label']} ({longest['est_seconds'] / 3600:.2f} h)" if longest else "- longest task: none"),
             f"- memory per task (largest): {mem_task:.0f} MB including the interpreter",
             f"- this machine: {phys} physical cores, {avail_gb:.1f} GB available; at most {max_by_mem} tasks fit in memory",
             ""]
    lines.append("| workers | wall time (longest-first) |")
    lines.append("|---|---|")
    for w in sorted({1, 2, 4, 8, 16, 32, workers, min(workers, max_by_mem)}):
        lines.append(f"| {w} | {cost.makespan(secs, w) / 3600:.1f} h |")
    lines += ["", "| mode | tasks | hours |", "|---|---|---|"]
    by_mode: dict[str, list[float]] = {}
    for t in tasks:
        by_mode.setdefault(t["mode"], []).append(t["est_seconds"])
    for m, v in by_mode.items():
        lines.append(f"| {m} | {len(v)} | {sum(v) / 3600:.1f} |")
    return "\n".join(lines) + "\n"


def read_tasks(path: str) -> list[dict[str, Any]]:
    tasks = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                tasks.append(json.loads(line))
    return tasks


__all__ = ["build_tasks", "write_tasks", "summary", "read_tasks", "calibrated_a_nm", "TSV_COLUMNS"]
