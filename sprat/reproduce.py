"""``sprat reproduce``: plan and run a campaign of the ``campaigns/`` directory.

A campaign directory holds ``tasks.jsonl`` (every task of the data set, fully resolved, built
from the deposited records by ``campaigns/manuscript/build_campaign.py``), the structure and
parameter files of its stages for reading, and ``expected/`` with the deposited registry and
flat table for ``sprat audit``.  Tiers select subsets by tag:

* ``quick``: bands, plane-wave layer, the smoke tests and the reference record (minutes);
* ``core``:  every record of the systematic runs (the figures of the paper), without the verification runs;
* ``full``:  everything.
"""

from __future__ import annotations

import os
from typing import Any

from . import cost, plan, runner

QUICK_TAGS = {"bands", "pwe", "smoke", "reference"}


def load_campaign(campaign_dir: str) -> list[dict[str, Any]]:
    path = os.path.join(campaign_dir, "tasks.jsonl")
    if not os.path.exists(path):
        raise FileNotFoundError(f"{path} not found; build it with campaigns/manuscript/build_campaign.py")
    return plan.read_tasks(path)


def select_tier(tasks: list[dict], tier: str) -> list[dict]:
    if tier == "full":
        return tasks
    if tier == "core":
        return [t for t in tasks if t.get("tier") != "verification"]
    if tier == "quick":
        return [t for t in tasks if t.get("tier") in QUICK_TAGS or t.get("quick")]
    raise ValueError(tier)


def prepare(campaign_dir: str, records_dir: str, tier: str = "full", throughput: float | None = None) -> str:
    """Write <records_dir>/_campaign/tasks_<tier>.jsonl (the resolved task list of the tier, with the record
    paths pointing into records_dir).  The campaign directory itself is never written to."""
    tasks = select_tier(load_campaign(campaign_dir), tier)
    thr = throughput or cost.read_calibration(records_dir)["throughput"]
    out = []
    for i, t in enumerate(tasks, start=1):
        t = dict(t)
        t["id"] = i
        t["record"] = os.path.join(records_dir, t["label"] + ".json")
        t["est_seconds"] = cost.seconds(t["estimate"], thr)
        out.append(t)
    os.makedirs(os.path.join(records_dir, "_campaign"), exist_ok=True)
    prefix = os.path.join(records_dir, "_campaign", f"tasks_{tier}")
    files = plan.write_tasks(out, prefix, thr)
    return files["jsonl"]


def reproduce_cli(campaign_dir: str, records_dir: str = "records", tier: str = "full", workers: int | None = None,
                  dry_run: bool = False) -> None:
    jsonl = prepare(campaign_dir, records_dir, tier)
    with open(jsonl.replace(".jsonl", ".summary.md"), encoding="utf-8") as fh:
        print(fh.read())
    if dry_run:
        print(f"dry run: task list written to {jsonl}")
        return
    s = runner.run(jsonl, workers=workers)
    if s.failed:
        raise SystemExit(1)
    print(f"done; run `sprat audit {records_dir} --expected "
          f"{os.path.join(campaign_dir, 'expected', 'numbers_registry_manuscript_v10_2.json')} -o tables`")


__all__ = ["load_campaign", "select_tier", "prepare", "reproduce_cli"]
