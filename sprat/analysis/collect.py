"""Collect a records directory into the flat table, and select records by parameters."""

from __future__ import annotations

import json
import os
from typing import Any, Iterable

from .. import records


def load(records_dir: str) -> list[dict]:
    return records.load_dir(records_dir)


def table(recs: list[dict], gap_from_bands: bool = False) -> list[dict[str, Any]]:
    """Flatten the records.  The validity criterion of the paper uses the fixed MPB gap-edge table
    (``records.GAP_EDGES_MPB``); ``gap_from_bands=True`` takes the edges from a bands sweep record instead."""
    gap = records.gap_table_from_records(recs) if gap_from_bands else None
    return [records.flatten(r, gap) for r in recs]


def collect_cli(records_dir: str, out: str) -> None:
    recs = load(records_dir)
    rows = table(recs)
    os.makedirs(out, exist_ok=True)
    records.write_table(rows, os.path.join(out, "records.csv"), os.path.join(out, "records.jsonl"))
    c = records.counts(rows)
    with open(os.path.join(out, "counts.json"), "w", encoding="utf-8") as fh:
        json.dump(c, fh, indent=1)
    print(json.dumps(c, indent=1))
    print(f"wrote {out}/records.csv, records.jsonl, counts.json ({len(rows)} rows)")


# --------------------------------------------------------------------------- selection
def _close(a, b, tol=1e-9) -> bool:
    return a is not None and b is not None and abs(float(a) - float(b)) <= tol


def matches(r: dict, **crit) -> bool:
    """Match a record against flat criteria: radius, row, dx, dy, n, resolution, cladding_rows, guide_periods,
    termination_x, mode, rows (string 'j:r;...'), tag, source, min_margin."""
    s, p = r["structure"], r["params"]
    d, c = s["defect"], s["cell"]
    flat = dict(mode=p["run"]["mode"], radius=d["radius"], row=d["row"], dx=d["dx"], dy=d["dy"], n=p["analyte"]["n"],
                resolution=p["numerics"]["resolution"], cladding_rows=c["cladding_rows"], guide_periods=c["guide_periods"],
                termination_x=c["termination_x"], pml=p["numerics"]["pml"],
                rows=";".join(f"{int(j)}:{float(v):g}" for j, v in sorted(((int(k), v) for k, v in (s.get("rows") or {}).items()))),
                tag=p["run"].get("tag", ""), source=r.get("task", {}).get("source", ""))
    for k, v in crit.items():
        if k == "min_margin":
            m = margin(r)
            if m is None or m < v:
                return False
        elif k == "no_rows":
            if v and flat["rows"]:
                return False
        elif isinstance(v, float):
            if not _close(flat.get(k), v, 1e-6):
                return False
        elif flat.get(k) != v:
            return False
    return True


def margin(r: dict) -> float | None:
    m = records.strongest_mode(r)
    ql = r["result"].get("Q_limit")
    return (ql / m["Q"]) if (m and ql) else None


def select(recs: Iterable[dict], **crit) -> list[dict]:
    return [r for r in recs if matches(r, **crit)]


def best(recs: Iterable[dict], **crit) -> dict | None:
    """The record with the longest harminv signal (highest Q_lim) among those matching."""
    cands = select(recs, **crit)
    if not cands:
        return None
    return max(cands, key=lambda r: (r["result"].get("t_used") or r["params"]["harminv"]["t"] or 0.0))


def Q_of(r: dict | None) -> float | None:
    m = records.strongest_mode(r) if r else None
    return float(m["Q"]) if m else None


def f_of(r: dict | None) -> float | None:
    m = records.strongest_mode(r) if r else None
    return float(m["f"]) if m else None


def by_key(recs: Iterable[dict], key) -> dict:
    """Group records by a key function, keeping the best (longest signal) per group."""
    out: dict = {}
    for r in recs:
        k = key(r)
        if k not in out or (r["result"].get("t_used") or 0) > (out[k]["result"].get("t_used") or 0):
            out[k] = r
    return out


__all__ = ["load", "table", "collect_cli", "matches", "select", "best", "margin", "Q_of", "f_of", "by_key"]
