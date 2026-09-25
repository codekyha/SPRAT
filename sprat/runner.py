"""Workstation parallel runner.

Every task is executed by a separate interpreter (``python -m sprat.cli run-one``), so the
memory of a finished Meep run is released and a crash cannot take the batch down.  The
runner keeps ``workers`` slots busy, pins each slot to its own cores (Linux, including WSL2;
macOS has no affinity API), orders the tasks longest-first, skips tasks whose record exists,
retries failures, writes a job log and prints an ETA from the cost model.
"""

from __future__ import annotations

import json
import os
import signal
import statistics
import subprocess
import sys
import time
from dataclasses import dataclass, field
from typing import Any

import psutil

from . import records
from .plan import read_tasks


@dataclass
class Slot:
    index: int
    cores: list[int]
    proc: subprocess.Popen | None = None
    task: dict | None = None
    started: float = 0.0
    attempt: int = 0
    log: str = ""


@dataclass
class RunSummary:
    total: int = 0
    done: int = 0
    skipped: int = 0
    failed: int = 0
    collisions: int = 0
    wall_s: float = 0.0
    joblog: str = ""
    failed_labels: list = field(default_factory=list)


def physical_cores() -> int:
    return psutil.cpu_count(logical=False) or os.cpu_count() or 1


def plan_slots(workers: int, cores_per_task: int, pin: bool) -> list[Slot]:
    """Contiguous core blocks per slot; wraps around when the machine has fewer cores."""
    total = os.cpu_count() or 1
    try:
        available = sorted(os.sched_getaffinity(0))          # respects an outer taskset / cgroup
    except AttributeError:
        available = list(range(total))
    slots = []
    for k in range(workers):
        start = (k * cores_per_task) % max(len(available), 1)
        cores = [available[(start + i) % len(available)] for i in range(cores_per_task)] if pin else []
        slots.append(Slot(index=k, cores=cores))
    return slots


def pin(proc: subprocess.Popen, cores: list[int]) -> None:
    if not cores:
        return
    try:
        os.sched_setaffinity(proc.pid, cores)
    except (AttributeError, OSError):
        try:
            psutil.Process(proc.pid).cpu_affinity(cores)
        except (AttributeError, psutil.Error):
            pass


def memory_cap(tasks: list[dict], workers: int, reserve_mb: float = 300.0) -> int:
    need = max((t.get("estimate", {}).get("memory_mb", 0.0) for t in tasks), default=0.0) + reserve_mb
    avail = psutil.virtual_memory().available / 1048576.0
    return max(1, min(workers, int(avail / need))) if need > 0 else workers


def _record_ok(task: dict, overwrite: bool) -> str:
    """'skip' if the record exists for the same task, 'collision' if it exists for a different one, else 'run'."""
    path = task["record"]
    if overwrite or not os.path.exists(path):
        return "run"
    try:
        old = records.load(path)
    except (ValueError, json.JSONDecodeError):
        return "run"
    probe = dict(structure=task["structure"], params=task["params"])
    return "skip" if records.same_task(old, probe) else "collision"


def run(tasks_file: str, workers: int | None = None, pin_cores: bool = True, cores_per_task: int = 1,
        order: str = "longest", overwrite: bool = False, retry: int = 1, dry_run: bool = False,
        log_dir: str = "logs", memory_guard: bool = True, python: str = sys.executable,
        poll_s: float = 1.0, quiet: bool = False) -> RunSummary:
    tasks = read_tasks(tasks_file)
    summary = RunSummary(total=len(tasks))
    run_id = time.strftime("%Y%m%d_%H%M%S")
    campaign = os.path.splitext(os.path.basename(tasks_file))[0]
    ldir = os.path.join(log_dir, f"{campaign}_{run_id}")
    os.makedirs(os.path.join(ldir, "tasks"), exist_ok=True)

    pending, skipped = [], []
    for t in tasks:
        state = _record_ok(t, overwrite)
        if state == "skip":
            skipped.append(t)
        elif state == "collision":
            summary.collisions += 1
            print(f"COLLISION: {t['record']} exists with different parameters; task {t['id']} not run "
                  "(use a different label or --overwrite)")
        else:
            pending.append(t)
    summary.skipped = len(skipped)
    if order == "longest":
        pending.sort(key=lambda t: -t.get("est_seconds", 0.0))
    elif order == "shortest":
        pending.sort(key=lambda t: t.get("est_seconds", 0.0))

    workers = workers or physical_cores()
    if memory_guard:
        capped = memory_cap(pending, workers)
        if capped < workers:
            print(f"memory guard: {workers} workers requested, {capped} fit in the available memory")
            workers = capped
    est_total = sum(t.get("est_seconds", 0.0) for t in pending)
    print(f"[sprat run] {len(tasks)} tasks: {len(pending)} to run, {len(skipped)} already done, "
          f"{summary.collisions} collisions; workers = {workers}, cores/task = {cores_per_task}, "
          f"pinning = {'on' if pin_cores else 'off'}; estimated {est_total / 3600:.1f} single-core hours")
    if dry_run or not pending:
        return summary

    slots = plan_slots(workers, cores_per_task, pin_cores)
    joblog = os.path.join(ldir, "joblog.tsv")
    summary.joblog = joblog
    with open(joblog, "w", encoding="utf-8") as fh:
        fh.write("id\tlabel\trc\tattempt\twall_s\tslot\tcores\tstarted\tfinished\test_s\tthroughput\n")

    stop = {"flag": False, "kill": False}

    def on_sigint(signum, frame):            # noqa: ARG001
        if stop["flag"]:
            stop["kill"] = True
            print("\n[sprat run] second interrupt: terminating running tasks")
        else:
            stop["flag"] = True
            print("\n[sprat run] interrupt: no new tasks will start; running tasks continue (Ctrl-C again to kill)")

    old_handler = signal.signal(signal.SIGINT, on_sigint)
    t_start = time.time()
    speed_ratios: list[float] = []
    attempts: dict[int, int] = {}
    env = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", SPRAT_JOB_ID=run_id,
               PYTHONUNBUFFERED="1")
    try:
        while pending or any(s.proc is not None for s in slots):
            # launch
            for s in slots:
                if s.proc is None and pending and not stop["flag"]:
                    t = pending.pop(0)
                    attempts[t["id"]] = attempts.get(t["id"], 0) + 1
                    tpath = os.path.join(ldir, "tasks", f"{t['id']}.json")
                    with open(tpath, "w", encoding="utf-8") as fh:
                        json.dump(t, fh)
                    s.log = os.path.join(ldir, f"task_{t['id']}_{attempts[t['id']]}.out")
                    out = open(s.log, "w", encoding="utf-8")
                    cmd = [python, "-m", "sprat.cli", "run-one", "--task", tpath] + (["--overwrite"] if overwrite else [])
                    e = dict(env, SPRAT_SLOT=str(s.index))
                    s.proc = subprocess.Popen(cmd, stdout=out, stderr=subprocess.STDOUT, env=e)
                    s.task, s.started, s.attempt = t, time.time(), attempts[t["id"]]
                    pin(s.proc, s.cores)
                    if not quiet:
                        print(f"[start] {t['id']:>4} slot {s.index:>2} cores {s.cores or '-'}  {t['label']}  "
                              f"(est {t.get('est_seconds', 0) / 3600:.2f} h)")
            # poll
            for s in slots:
                if s.proc is None:
                    continue
                if stop["kill"]:
                    s.proc.terminate()
                rc = s.proc.poll()
                if rc is None:
                    continue
                wall = time.time() - s.started
                t = s.task
                thr = ""
                if rc == 0 and os.path.exists(t["record"]):
                    try:
                        thr = f"{records.load(t['record'])['provenance'].get('throughput', 0):.3e}"
                    except Exception:      # noqa: BLE001
                        thr = ""
                    est = t.get("est_seconds", 0.0)
                    if est > 0 and wall > 0:
                        speed_ratios.append(est / wall)
                    summary.done += 1
                else:
                    if s.attempt <= retry and not stop["flag"]:
                        print(f"[retry] {t['id']} rc={rc} (attempt {s.attempt}); see {s.log}")
                        pending.insert(0, t)
                    else:
                        summary.failed += 1
                        summary.failed_labels.append(t["label"])
                        print(f"[FAILED] {t['id']} {t['label']} rc={rc}; see {s.log}")
                with open(joblog, "a", encoding="utf-8") as fh:
                    fh.write(f"{t['id']}\t{t['label']}\t{rc}\t{s.attempt}\t{wall:.1f}\t{s.index}\t"
                             f"{'-'.join(map(str, s.cores))}\t{time.strftime('%FT%T', time.localtime(s.started))}\t"
                             f"{time.strftime('%FT%T')}\t{t.get('est_seconds', 0):.0f}\t{thr}\n")
                if rc == 0 and not quiet:
                    remaining = sum(x.get("est_seconds", 0.0) for x in pending) + \
                        sum(max(x.task.get("est_seconds", 0.0) - (time.time() - x.started), 0.0)
                            for x in slots if x.proc is not None and x is not s)
                    factor = statistics.median(speed_ratios) if speed_ratios else 1.0
                    eta = remaining / max(factor, 1e-9) / max(workers, 1)
                    print(f"[done ] {t['id']:>4} in {wall / 3600:.2f} h  ({summary.done} done, {len(pending)} pending; "
                          f"speed x{factor:.2f} of estimate; ETA about {eta / 3600:.1f} h)")
                s.proc, s.task = None, None
            time.sleep(poll_s)
    finally:
        signal.signal(signal.SIGINT, old_handler)
    summary.wall_s = time.time() - t_start
    print(f"[sprat run] finished: {summary.done} done, {summary.skipped} skipped, {summary.failed} failed, "
          f"{summary.collisions} collisions in {summary.wall_s / 3600:.2f} h; joblog {joblog}")
    return summary


def status(tasks_file: str, log_dir: str = "logs") -> dict[str, Any]:
    tasks = read_tasks(tasks_file)
    done = [t for t in tasks if os.path.exists(t["record"])]
    pending = [t for t in tasks if not os.path.exists(t["record"])]
    campaign = os.path.splitext(os.path.basename(tasks_file))[0]
    logs = sorted(d for d in (os.listdir(log_dir) if os.path.isdir(log_dir) else []) if d.startswith(campaign + "_"))
    failed = []
    thr = []
    if logs:
        jl = os.path.join(log_dir, logs[-1], "joblog.tsv")
        if os.path.exists(jl):
            with open(jl, encoding="utf-8") as fh:
                rows = [line.rstrip("\n").split("\t") for line in fh][1:]
            failed = [r[1] for r in rows if r[2] != "0"]
            thr = [float(r[10]) for r in rows if len(r) > 10 and r[10]]
    remaining_est = sum(t.get("est_seconds", 0.0) for t in pending)
    return dict(total=len(tasks), done=len(done), pending=len(pending), failed_in_last_run=len(failed),
                remaining_single_core_hours=remaining_est / 3600,
                median_throughput=(statistics.median(thr) if thr else None), last_log=(logs[-1] if logs else ""))


__all__ = ["run", "status", "physical_cores", "plan_slots", "memory_cap", "RunSummary"]
