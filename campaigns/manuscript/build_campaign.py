#!/usr/bin/env python3
"""Build tasks.jsonl, the task list of the data set of the paper, from the imported deposit.

Every record of the deposit (858 harmonic inversions, 7 spectra and their 7 references in the 25a cell, 7 spectra and
their 7 references with the guide continued into the absorber, 2 field maps, 4 band computations, the plane-wave
layer) becomes one fully resolved task with an English label, so that ``sprat reproduce campaigns/manuscript``
regenerates the whole data set point by point.  Exact duplicates (same structure and parameters) are merged; a label
that would collide with a different task gets a numeric suffix.

Usage: python build_campaign.py <records dir produced by `sprat import-legacy`> [task-list dir] [--append]

``--append`` keeps the existing tasks.jsonl line for line and adds a task for every imported deposited file that no
task names yet (SPRAT 1.2.0 added the 14 absorber-terminated spectrum runs of raw_h14b this way); records written by
SPRAT itself are not deposited files and are skipped.
"""
from __future__ import annotations

import json
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

from sprat import cost                              # noqa: E402
from sprat.analysis import collect                  # noqa: E402
from sprat.params import resolve                    # noqa: E402
from sprat.structure import structure_from_dict, structure_to_json   # noqa: E402

TIER = {"legacy-verification": "verification", "legacy-campaign": "campaign"}


def _key(s, rp) -> tuple[str, str]:
    return (json.dumps(structure_to_json(s), sort_keys=True),
            json.dumps({k: v for k, v in rp.items() if k != "run"} | {"run": {k: v for k, v in rp["run"].items()
                                                                            if k not in ("label", "tag", "quiet", "overwrite")}},
                       sort_keys=True))


def _task(r: dict, next_id: int, labels: dict) -> tuple[tuple[str, str], dict]:
    s = structure_from_dict(r["structure"])
    p = json.loads(json.dumps(r["params"]))
    mode = p["run"]["mode"]
    legacy_label = r["task"]["label"]
    tier = TIER.get(r["task"].get("source", ""), "campaign")
    if mode == "harminv" and r["result"].get("t_used"):
        p["harminv"]["t"] = float(r["result"]["t_used"])
    p["run"]["label"] = "auto"
    p["run"]["tag"] = "verification" if tier == "verification" else (p["run"].get("tag") or "campaign")
    if r["task"].get("legacy_wp") == "smoke":              # the smoke tests of 2026 (labels v4_duman_*)
        p["run"]["tag"], tier = "smoke", "smoke"
    if mode == "bands":
        tier = "bands"
    if mode == "pwe":
        tier = "pwe"
    rp = resolve(s, p)
    label = rp["run"]["label"]
    if label in labels:
        n = 2
        while f"{label}_v{n}" in labels:
            n += 1
        label = f"{label}_v{n}"
        rp["run"]["label"] = label
    est = cost.estimate(s, rp)
    t = dict(id=next_id, label=label, mode=mode, tier=tier, structure=structure_to_json(s), params=rp,
             overrides={}, structure_file="campaigns/manuscript/structure_reference.phc", params_file="campaigns/manuscript/stages",
             record=os.path.join("records", label + ".json"), estimate=est, est_seconds=cost.seconds(est, cost.REFERENCE_THROUGHPUT),
             created=time.strftime("%Y-%m-%dT%H:%M:%S"), legacy_labels=[legacy_label],
             quick=(mode in ("bands", "pwe") or tier == "smoke" or legacy_label == "harminv_tarama_rd0.060_dx+0.000_dy+0.000_s4_na1.3300_res24_nk12_simoto"))
    if t["quick"] and tier == "campaign":
        t["tier"] = "reference"
    return _key(s, rp), t


def _summary(out_dir: str, tasks: list[dict], n_records: int, dropped: int, note: str = "") -> None:
    hours = sum(t["est_seconds"] for t in tasks) / 3600
    by: dict[str, list] = {}
    for t in tasks:
        by.setdefault(t["tier"], [0, 0.0])
        by[t["tier"]][0] += 1
        by[t["tier"]][1] += t["est_seconds"] / 3600
    with open(os.path.join(out_dir, "tasks_summary.md"), "w", encoding="utf-8") as fh:
        fh.write(f"# Task list of the data set\n\n{note}{len(tasks)} tasks from {n_records} deposited records "
                 f"({dropped} exact duplicates merged); {hours:.0f} single-core hours at the Altay throughput\n\n"
                 "| tier | tasks | single-core hours |\n|---|---|---|\n")
        for k, (n, h) in sorted(by.items()):
            fh.write(f"| {k} | {n} | {h:.1f} |\n")
    print(f"{len(tasks)} tasks in {out_dir}/tasks.jsonl ({dropped} duplicates merged, {hours:.0f} single-core hours)")
    for k, (n, h) in sorted(by.items()):
        print(f"  {k:14s} {n:4d} tasks  {h:7.1f} h")


def main(records_dir: str, out_dir: str) -> None:
    recs = collect.load(records_dir)
    tasks, seen, labels = [], {}, {}
    dropped = 0
    for r in sorted(recs, key=lambda r: r["task"]["label"]):
        key, t = _task(r, len(tasks) + 1, labels)
        if key in seen:
            seen[key]["legacy_labels"].append(r["task"]["label"])
            dropped += 1
            continue
        labels[t["label"]] = True
        seen[key] = t
        tasks.append(t)
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "tasks.jsonl"), "w", encoding="utf-8") as fh:
        for t in tasks:
            fh.write(json.dumps(t) + "\n")
    _summary(out_dir, tasks, len(recs), dropped, f"built {time.strftime('%Y-%m-%d')}: ")


def append(records_dir: str, out_dir: str) -> None:
    path = os.path.join(out_dir, "tasks.jsonl")
    with open(path, encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    tasks = [json.loads(x) for x in lines]
    covered = {lab for t in tasks for lab in t["legacy_labels"]}
    labels = {t["label"]: True for t in tasks}
    seen = {_key(structure_from_dict(t["structure"]), t["params"]): t for t in tasks}
    recs = collect.load(records_dir)
    new = []
    for r in sorted(recs, key=lambda r: r["task"]["label"]):
        if r["task"]["label"] in covered or not r["task"].get("source", "").startswith("legacy"):
            continue                                        # a task exists, or the record is not a deposited file
        key, t = _task(r, len(tasks) + len(new) + 1, labels)
        if key in seen:
            raise SystemExit(f"{r['task']['label']}: the same task as {seen[key]['label']}; append does not merge")
        labels[t["label"]] = True
        seen[key] = t
        new.append(t)
    with open(path, "a", encoding="utf-8") as fh:
        for t in new:
            fh.write(json.dumps(t) + "\n")
    print(f"appended {len(new)} tasks")
    merged = sum(len(t["legacy_labels"]) - 1 for t in tasks + new)
    _summary(out_dir, tasks + new, sum(len(t["legacy_labels"]) for t in tasks + new), merged,
             f"built 2026-09-23; {len(new)} tasks appended {time.strftime('%Y-%m-%d')} (the absorber-terminated spectrum runs): ")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--append"]
    if not args or args[0] in ("-h", "--help") or not os.path.isdir(args[0]):
        sys.exit(__doc__)
    target = args[1] if len(args) > 1 else os.path.dirname(os.path.abspath(__file__))
    (append if "--append" in sys.argv else main)(args[0], target)
