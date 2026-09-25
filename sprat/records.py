"""The record: one self-describing JSON file per run, and the flat table built from many.

Record layout (schema ``sprat-record-1.0``)::

    {
     "schema":     "sprat-record-1.0",
     "software":   {"name", "version", "meep", "python", "code_sha256"},
     "task":       {"label", "tag", "created", "structure_file", "params_file", "overrides"},
     "structure":  {... the resolved .phc content ...},
     "params":     {... the resolved .par content, sweeps applied ...},
     "provenance": {"host", "started", "wall_s", "steps", "pixels", "pixel_steps", "throughput",
                    "symmetry", "slot", "job_id", "legacy_file", "legacy_sha256"},
     "result":     {mode-specific, see below}
    }

``result`` by mode:

* harminv:   ``modes`` (list of {f, Q, amplitude, error, wavelength_nm, fwhm_nm}, strongest first),
             ``t_used``, ``auto_history``, ``Q_limit``;
* spectrum / reference: ``f``, ``flux_in``, ``flux_out``, ``t_end``;
* field:     ``analyte_energy_fraction``, ``wavelength_nm``, ``S_first_order_nm_per_RIU``,
             ``field_file`` (the .npz beside the record);
* bands:     ``task``, ``gap`` (lower, upper, mid, relative_width_pct, exists), ``freqs``, ``k``, ``a_nm``;
* pwe:       the plane-wave results (see ``sprat.pwe``).
"""

from __future__ import annotations

import csv
import glob
import hashlib
import json
import os
import platform
import sys
import time
from typing import Any, Iterable

from . import SCHEMA_VERSION, __version__
from .structure import structure_from_dict, structure_to_json

# TM gap edges (lower, upper) against the analyte index from the MPB sweep of the campaign
# (records 'bands_sweep' override this table when present).
GAP_EDGES_MPB = [(1.00, 0.28202, 0.41790), (1.05, 0.28112, 0.40654), (1.10, 0.28017, 0.39520),
                 (1.15, 0.27918, 0.38399), (1.20, 0.27814, 0.37300), (1.25, 0.27707, 0.36229),
                 (1.30, 0.27595, 0.35191), (1.35, 0.27480, 0.34188), (1.40, 0.27360, 0.33223),
                 (1.45, 0.27236, 0.32296), (1.50, 0.27109, 0.31407)]

GAP_WINDOW = (0.2, 0.8)          # admissible position of the resonance inside the gap
MARGIN_MIN = 1.0                 # Q_lim / Q >= 1: the mode is resolved
DEAD_ZONE_RADIUS = 0.14          # r_d >= 0.14a: the strongest mode is a mode of the finite waveguide at every radius
PRELIMINARY_SWEEP = (6, 20)      # (cladding rows, resolution) of the preliminary sweep of the 2026 campaign
GUIDE_MODE_RADIUS = 0.13         # in that sweep the admitted records at N_sep = 1 or r_d >= 0.13a hold waveguide modes


def waveguide_mode(row: dict) -> bool:
    """The rule for the admitted records of the preliminary sweep (n_cl = 6, resolution 20): at N_sep = 1 or
    r_d >= 0.13a the strongest mode that harmonic inversion returns is a low-Q Fabry-Perot mode of the finite waveguide
    (Q 50 to 77, median 53, wavelength near 1623 nm), whereas the admitted cavity records of the same sweep have Q 145
    to 660. Weaker modes listed beside it are not used. The validity criterion tests the numerics only and admits these
    records. The rule applies to admitted rows (the flag ``waveguide_mode`` of the flat table is set on admitted rows
    only): at r_d = 0.13a some excluded records of the sweep hold the cavity mode, close to the lower edge of the gap."""
    if row.get("radius") is None or row.get("row") is None:
        return False
    return (row["cladding_rows"], row["resolution"]) == PRELIMINARY_SWEEP and (
        row["row"] == 1 or row["radius"] >= GUIDE_MODE_RADIUS - 1e-9)

FLAT_COLUMNS = ["file", "source", "mode", "tag", "label", "radius", "dx", "dy", "row", "n", "k", "resolution",
                "cladding_rows", "guide_periods", "a_nm", "fcen", "fwidth", "pml", "pad_x", "pad_y", "termination_x",
                "absorber_periods", "rows", "t_harminv", "symmetry", "Q", "f_r", "wavelength_nm", "fwhm_nm",
                "harminv_error", "amplitude", "Q_limit", "margin", "gap_position", "resolved", "valid",
                "exclusion_reason", "waveguide_mode", "mode_count", "steps", "pixels", "throughput", "wall_s", "host",
                "job_id", "code_sha256", "started", "legacy_file"]


# --------------------------------------------------------------------------- gap edges
def gap_edges(n: float, table: list[tuple[float, float, float]] | None = None) -> tuple[float, float]:
    """Linear interpolation of the TM gap edges at analyte index n."""
    tab = table or GAP_EDGES_MPB
    if n <= tab[0][0]:
        return tab[0][1], tab[0][2]
    if n >= tab[-1][0]:
        return tab[-1][1], tab[-1][2]
    for (n0, l0, u0), (n1, l1, u1) in zip(tab, tab[1:]):
        if n0 <= n <= n1:
            w = (n - n0) / (n1 - n0)
            return l0 + w * (l1 - l0), u0 + w * (u1 - u0)
    raise ValueError(n)


def gap_table_from_records(records: Iterable[dict]) -> list[tuple[float, float, float]] | None:
    """Build the gap-edge table from a 'bands' sweep record if one is present."""
    for r in records:
        if r.get("params", {}).get("run", {}).get("mode") == "bands" and r["result"].get("task") == "sweep":
            rows = [(x["n"], x["gap"]["lower"], x["gap"]["upper"]) for x in r["result"]["sweep"] if x["gap"]["exists"]]
            if len(rows) >= 2:
                return sorted(rows)
    return None


# --------------------------------------------------------------------------- record i/o
def code_sha256(path: str | None = None) -> str:
    """SHA-256 of the FDTD module (the code that writes the records)."""
    try:
        p = path or os.path.join(os.path.dirname(os.path.abspath(__file__)), "fdtd.py")
        with open(p, "rb") as fh:
            return hashlib.sha256(fh.read()).hexdigest()
    except OSError:
        return ""


def software_block(meep_version: str = "") -> dict[str, Any]:
    return dict(name="sprat", version=__version__, meep=meep_version, python=platform.python_version(),
                code_sha256=code_sha256())


def new_record(structure: dict, params: dict, task: dict | None = None) -> dict[str, Any]:
    return dict(schema=SCHEMA_VERSION, software=software_block(), task=dict(task or {}),
                structure=structure_to_json(structure), params=params,
                provenance=dict(host=platform.node(), started=time.strftime("%Y-%m-%dT%H:%M:%S"),
                                job_id=os.environ.get("SPRAT_JOB_ID", ""), slot=os.environ.get("SPRAT_SLOT", "")),
                result={})


def record_path(records_dir: str, label: str) -> str:
    return os.path.join(records_dir, label + ".json")


def save(record: dict, path: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(record, fh, indent=1, default=_json_default)
    os.replace(tmp, path)


def _json_default(o):
    try:
        import numpy as np
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (np.bool_,)):
            return bool(o)
    except ImportError:
        pass
    raise TypeError(f"not JSON serialisable: {type(o)}")


def load(path: str) -> dict[str, Any]:
    with open(path, encoding="utf-8") as fh:
        r = json.load(fh)
    if r.get("schema") != SCHEMA_VERSION:
        raise ValueError(f"{path}: not a {SCHEMA_VERSION} record (schema = {r.get('schema')!r}); "
                         "legacy files must be converted with `sprat import-legacy`")
    r["_file"] = path
    return r


def load_dir(records_dir: str, pattern: str = "*.json") -> list[dict[str, Any]]:
    out = []
    for path in sorted(glob.glob(os.path.join(records_dir, pattern))):
        name = os.path.basename(path)
        if name.startswith("_") or name in ("calibration.json",):
            continue
        try:
            out.append(load(path))
        except (ValueError, json.JSONDecodeError) as e:
            print(f"skipped {path}: {e}", file=sys.stderr)
    return out


def structure_of(record: dict) -> dict:
    return structure_from_dict(record["structure"])


def mode_of(record: dict) -> str:
    return record["params"]["run"]["mode"]


def strongest_mode(record: dict) -> dict | None:
    modes = record.get("result", {}).get("modes") or []
    if not modes:
        return None
    return max(modes, key=lambda m: m.get("amplitude", 0.0))


def same_task(a: dict, b: dict) -> bool:
    """Two records describe the same task when their structure and parameters agree."""
    return _strip_labels(a) == _strip_labels(b)


def _strip_labels(r: dict) -> tuple:
    p = json.loads(json.dumps(r["params"]))
    p.get("run", {}).pop("label", None)
    p.get("run", {}).pop("tag", None)
    p.get("run", {}).pop("overwrite", None)
    p.get("run", {}).pop("quiet", None)
    return json.dumps(r["structure"], sort_keys=True), json.dumps(p, sort_keys=True)


# --------------------------------------------------------------------------- flattening
def flatten(record: dict, gap_table: list | None = None) -> dict[str, Any]:
    """One row of the flat table with the validity flags of the paper."""
    s, p, o, res = record["structure"], record["params"], record.get("provenance", {}), record.get("result", {})
    d, c, la = s["defect"], s["cell"], s["lattice"]
    m = strongest_mode(record)
    Q = m["Q"] if m else None
    f = m["f"] if m else None
    Qlim = res.get("Q_limit")
    margin = (Qlim / Q) if (Q and Qlim) else None
    n = float(p["analyte"]["n"])
    lo, hi = gap_edges(n, gap_table)
    position = ((f - lo) / (hi - lo)) if f else None
    resolved = margin is not None and margin >= MARGIN_MIN
    valid = bool(resolved and position is not None and GAP_WINDOW[0] <= position <= GAP_WINDOW[1])
    if valid or p["run"]["mode"] != "harminv":         # the validity flags belong to harmonic inversions only
        reason = ""
    elif Q is None:
        reason = "no mode"
    elif not resolved:
        reason = "Q > Q_lim"
    else:
        reason = "gap position outside 0.2-0.8"
    row = dict(
        file=os.path.basename(record.get("_file", "")), source=record.get("task", {}).get("source", "sprat"),
        mode=p["run"]["mode"], tag=p["run"].get("tag", ""), label=record["task"].get("label", p["run"].get("label")),
        radius=d["radius"] if d["type"] == "rod" else None, dx=d["dx"], dy=d["dy"], row=d["row"] if d["type"] == "rod" else None,
        n=n, k=p["analyte"]["k"], resolution=p["numerics"]["resolution"], cladding_rows=c["cladding_rows"],
        guide_periods=c["guide_periods"], a_nm=la["a_nm"], fcen=p["source"]["fcen"], fwidth=p["source"]["fwidth"],
        pml=p["numerics"]["pml"], pad_x=c["pad_x"], pad_y=c["pad_y"], termination_x=c["termination_x"],
        absorber_periods=(c["absorber_periods"] if c["termination_x"] == "absorber" else None),
        rows=";".join(f"{j}:{r}" for j, r in sorted(((int(k), v) for k, v in (s.get("rows") or {}).items()))),
        t_harminv=res.get("t_used", p["harminv"]["t"] if p["run"]["mode"] == "harminv" else None),
        symmetry=o.get("symmetry"), Q=Q, f_r=f, wavelength_nm=(m.get("wavelength_nm") if m else None),
        fwhm_nm=(m.get("fwhm_nm") if m else None), harminv_error=(m.get("error") if m else None),
        amplitude=(m.get("amplitude") if m else None), Q_limit=Qlim, margin=margin, gap_position=position,
        resolved=resolved, valid=valid, exclusion_reason=reason, mode_count=len(res.get("modes") or []),
        steps=o.get("steps"), pixels=o.get("pixels"), throughput=o.get("throughput"), wall_s=o.get("wall_s"),
        host=o.get("host"), job_id=o.get("job_id", ""), code_sha256=record.get("software", {}).get("code_sha256", ""),
        started=o.get("started"), legacy_file=o.get("legacy_file", ""))
    # admitted records of the preliminary sweep that hold a mode of the finite waveguide (see waveguide_mode)
    row["waveguide_mode"] = bool(valid and p["run"]["mode"] == "harminv" and waveguide_mode(row))
    return row


def write_table(rows: list[dict], csv_path: str, jsonl_path: str | None = None) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(csv_path)), exist_ok=True)
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FLAT_COLUMNS, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    if jsonl_path:
        with open(jsonl_path, "w", encoding="utf-8") as fh:
            for r in rows:
                fh.write(json.dumps(r, default=_json_default) + "\n")


def admitted(record: dict) -> bool:
    """True when a harmonic-inversion record satisfies the validity criterion of the paper (margin >= 1 and position
    inside the TM gap between 0.2 and 0.8).  Fits and statistics use admitted records only; excluded records may be
    shown in a figure, marked as such."""
    return bool(flatten(record)["valid"])


def counts(rows: list[dict]) -> dict[str, int]:
    h = [r for r in rows if r["mode"] == "harminv"]
    ver = [r for r in h if r["source"] == "legacy-verification" or r["tag"] == "verification"]
    ref = [r for r in h if not r["valid"] and r["radius"] == 0.06 and r["row"] == 4 and r["cladding_rows"] == 12 and r["resolution"] == 24
           and r["termination_x"] == "pml" and r["guide_periods"] == 25 and r["dx"] == 0 and r["dy"] == 0 and not r["rows"]]
    return dict(harminv=len(h), verification=len(ver), campaign=len(h) - len(ver), excluded_reference=len(ref),
                valid=sum(1 for r in h if r["valid"]), resolved=sum(1 for r in h if r["resolved"]),
                excluded_margin=sum(1 for r in h if r["exclusion_reason"] == "Q > Q_lim"),
                excluded_gap=sum(1 for r in h if r["exclusion_reason"] == "gap position outside 0.2-0.8"),
                dead_zone=sum(1 for r in h if r["radius"] >= DEAD_ZONE_RADIUS - 1e-9),
                dead_zone_admitted=sum(1 for r in h if r["valid"] and r["radius"] >= DEAD_ZONE_RADIUS - 1e-9),
                waveguide_mode_admitted=sum(1 for r in h if r.get("waveguide_mode") is True),
                spectra=sum(1 for r in rows if r["mode"] == "spectrum"),
                references=sum(1 for r in rows if r["mode"] == "reference"),
                spectra_absorber=sum(1 for r in rows if r["mode"] == "spectrum" and r["termination_x"] == "absorber"),
                references_absorber=sum(1 for r in rows if r["mode"] == "reference" and r["termination_x"] == "absorber"),
                fields=sum(1 for r in rows if r["mode"] == "field"),
                bands=sum(1 for r in rows if r["mode"] == "bands"), pwe=sum(1 for r in rows if r["mode"] == "pwe"))


__all__ = ["GAP_EDGES_MPB", "GAP_WINDOW", "MARGIN_MIN", "DEAD_ZONE_RADIUS", "PRELIMINARY_SWEEP", "GUIDE_MODE_RADIUS",
           "waveguide_mode", "FLAT_COLUMNS", "gap_edges", "gap_table_from_records",
           "code_sha256", "software_block", "new_record", "record_path", "save", "load", "load_dir", "structure_of",
           "mode_of", "strongest_mode", "same_task", "flatten", "admitted", "write_table", "counts"]
