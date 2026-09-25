"""Record-level statistics quoted in the paper and its supplementary material.

* ``census``: the admission criterion by source and layer (supplementary tables S1 and S2), the records above the
  harmonic-inversion ceiling, the margins achieved by layer and by verification series (table S3), the target
  margins of the verification runs, and the offset of the source centre frequency from the resonance (section S3).
* ``repeats``: points run more than once, and how closely the repeats at a higher margin reproduce the run with the
  lowest margin (section 2.4 and section S3).
* ``band_gap``: the band-structure records summarised as in section 3.1.
"""

from __future__ import annotations

import json
from typing import Any

import numpy as np

from ..records import GAP_WINDOW, MARGIN_MIN, flatten, strongest_mode
from .collect import margin as _margin
from .sensitivity import is_verification

GAP_FAIL = "gap position outside 0.2-0.8"
REPEAT_FACTOR = 1.5          # a repeat reaches at least 1.5 times the margin of the run it repeats


def _layer(r: dict) -> str:
    return r["task"].get("legacy_layer") or r["params"]["run"].get("tag") or ""


def _source(r: dict) -> str:
    return "verification" if is_verification(r) else "campaign"


def _stats(v: list[float]) -> dict[str, Any]:
    return dict(records=len(v), min=float(min(v)), max=float(max(v)), median=float(np.median(v))) if v else dict(records=0)


def census(recs: list[dict]) -> dict[str, Any]:
    h = [r for r in recs if r["params"]["run"]["mode"] == "harminv"]
    fl = [(r, flatten(r)) for r in h]
    table: dict[str, dict] = {}
    for r, f in fl:
        t = table.setdefault(f"{_source(r)}/{_layer(r)}", dict(records=0, admitted=0, margin_fail=0, gap_fail=0))
        t["records"] += 1
        t["admitted"] += int(bool(f["valid"]))
        t["margin_fail"] += int(f["exclusion_reason"] == "Q > Q_lim")
        t["gap_fail"] += int(f["exclusion_reason"] == GAP_FAIL)
    out: dict[str, Any] = dict(by_source_layer=table,
                               total={k: sum(v[k] for v in table.values()) for k in ("records", "admitted", "margin_fail", "gap_fail")})
    inside = lambda f: f["gap_position"] is not None and GAP_WINDOW[0] <= f["gap_position"] <= GAP_WINDOW[1]     # noqa: E731
    over = [f for _, f in fl if f["margin"] is not None and f["margin"] < MARGIN_MIN]
    if over:
        out["over_ceiling"] = dict(records=len(over), margin_range=[min(f["margin"] for f in over), max(f["margin"] for f in over)],
                                   radii=sorted({round(f["radius"], 4) for f in over}),
                                   also_outside_gap=sum(1 for f in over if not inside(f)),
                                   removed_by_margin_alone=sum(1 for f in over if inside(f)))
    gap = [f for _, f in fl if f["exclusion_reason"] == GAP_FAIL]
    if gap:
        byr: dict = {}
        for f in gap:
            byr.setdefault(round(f["radius"], 4), []).append(f["gap_position"])
        out["outside_gap_by_radius"] = {f"{rd:.3f}": dict(records=len(v), gap_position_range=[min(v), max(v)]) for rd, v in sorted(byr.items())}
        out["outside_gap_min_radius"] = min(byr)
        out["outside_gap_analyte_indices"] = sorted({round(f["n"], 4) for f in gap})
    # the ceiling Q_lim = pi f_cen t is built on the source centre, not on the resonance: the two differ by 1 - f_r / f_cen
    dev = [(100.0 * abs(1.0 - f["f_r"] / f["fcen"]), round(f["radius"], 4)) for _, f in fl if f["f_r"] and f["fcen"]]
    if dev:
        mx = max(dev)
        ref = [100.0 * abs(1.0 - f["f_r"] / f["fcen"]) for r, f in fl if f["f_r"] and f["radius"] == 0.06 and f["row"] == 4
               and f["cladding_rows"] == 12 and f["resolution"] == 24 and f["termination_x"] == "pml" and f["guide_periods"] == 25
               and f["n"] == 1.33 and f["dx"] == 0 and f["dy"] == 0 and not f["rows"] and _source(r) == "campaign"]
        out["fcen_offset_pct"] = dict(max=mx[0], at_radius=mx[1], reference=(ref[0] if ref else None),
                                      fcen=sorted({round(f["fcen"], 6) for _, f in fl if f["fcen"]}))
    # margins achieved, by source and layer
    out["margin_by_source_layer"] = {k: _stats([f["margin"] for r, f in fl if f"{_source(r)}/{_layer(r)}" == k and f["margin"] is not None])
                                     for k in table}
    # verification runs: margins by series (task.legacy_wp) and the target each was run to, recovered from the rule
    # t = 100 ceil(p Q_est / (100 pi f_cen)) of section S3 with the recorded t, Q_est and f_cen
    ver = [(r, f) for r, f in fl if _source(r) == "verification" and _layer(r) == "production"]
    series: dict[str, list] = {}
    targets = []
    for r, f in ver:
        wp = r["task"].get("legacy_wp") or ""
        series.setdefault(wp, []).append(f["margin"])
        q_est = r["params"]["run"].get("q_est")
        t = r["result"].get("t_used") or r["params"]["harminv"].get("t")
        # the ceiling in the rule makes t overshoot, so the target is t pi f_cen / Q_est truncated to 0.1
        p = (np.floor(float(t) * np.pi * float(f["fcen"]) / float(q_est) * 10 + 1e-6) / 10) if (q_est and t and f["fcen"]) else None
        p = float(p) if p is not None else None
        targets.append((wp, p, f["margin"]))
    out["verification_by_series"] = {wp: dict(records=len(v), margin_range=[min(v), max(v)]) for wp, v in sorted(series.items())}
    in_series = [(wp, p, m) for wp, p, m in targets if wp.startswith("H")]
    tg: dict = {}
    for _, p, _ in in_series:
        tg[f"{p:g}" if p is not None else "unknown"] = tg.get(f"{p:g}" if p is not None else "unknown", 0) + 1
    out["verification_targets"] = dict(production_records=len(ver), series_records=len(in_series), by_target=tg,
                                       below_target=sum(1 for _, p, m in targets if p is not None and m is not None and m < p - 1e-9),
                                       rule="p = t pi f_cen / Q_est truncated to 0.1 (the inverse of equation S2)")
    return out


def _point_key(r: dict) -> tuple:
    s, p = r["structure"], r["params"]
    d, c = s["defect"], s["cell"]
    return (round(float(d["radius"]), 4), int(d["row"]), round(float(d["dx"]), 4), round(float(d["dy"]), 4), round(float(p["analyte"]["n"]), 4),
            int(p["numerics"]["resolution"]), int(c["cladding_rows"]), int(c["guide_periods"]), c["termination_x"], float(p["numerics"]["pml"]),
            json.dumps(s.get("rows") or {}, sort_keys=True), (float(c.get("absorber_periods") or 0) if c["termination_x"] == "absorber" else None))


POINTS = {
    "reference_N5": dict(radius=0.06, row=5, cladding_rows=12, termination_x="pml"),
    "reference": dict(radius=0.06, row=4, cladding_rows=12, termination_x="pml"),
    "absorber_rd0.085": dict(radius=0.085, row=4, cladding_rows=12, termination_x="absorber"),
}


def repeats(recs: list[dict], factor: float = REPEAT_FACTOR) -> dict[str, Any]:
    """Points run more than once at identical parameters.  The base run is the campaign run of lowest margin (any run
    if the point has no campaign run); a repeat is a run of the same point with at least ``factor`` times that
    margin.  The reproduction is the largest |Q_repeat / Q_base - 1|."""
    groups: dict = {}
    for r in recs:
        if r["params"]["run"]["mode"] == "harminv" and strongest_mode(r) and _margin(r):
            groups.setdefault(_point_key(r), []).append(r)
    rows = []
    for k, lst in groups.items():
        if len(lst) < 2:
            continue
        camp = [r for r in lst if not is_verification(r)]
        base = min(camp or lst, key=_margin)
        rep = [r for r in lst if r is not base and _margin(r) >= factor * _margin(base)]
        if not rep:
            continue
        qb = strongest_mode(base)["Q"]
        rows.append(dict(point=dict(radius=k[0], row=k[1], dx=k[2], dy=k[3], n=k[4], resolution=k[5], cladding_rows=k[6], guide_periods=k[7],
                                    termination_x=k[8], pml=k[9], rows=k[10]),
                         base=dict(label=base["task"]["label"], Q=qb, margin=_margin(base)),
                         repeats=[dict(label=r["task"]["label"], Q=strongest_mode(r)["Q"], margin=_margin(r)) for r in rep],
                         max_deviation_pct=100.0 * max(abs(strongest_mode(r)["Q"] / qb - 1.0) for r in rep)))
    ref = dict(dx=0.0, dy=0.0, n=1.33, resolution=24, guide_periods=25, pml=1.0, rows="{}")
    named = {}
    for name, crit in POINTS.items():
        want = dict(ref, **crit)
        hit = [g for g in rows if all((abs(g["point"][k] - v) < 1e-9) if isinstance(v, float) else g["point"][k] == v for k, v in want.items())]
        named[name] = hit[0] if hit else None
    return dict(rule="base: lowest-margin run of a point among the systematic runs (any run if it has none); repeats: runs of the same point with at least %.1f times its margin" % factor,
                groups=rows, points=named)


def band_gap(recs: list[dict], pwe: dict | None = None, design_range=(1.30, 1.45)) -> dict[str, Any] | None:
    """The band-structure records: the TM gap at the analyte index of the bulk record, the lattice-constant calibration,
    the wavelength window of the gap, the gap against the analyte index, and the plane-wave cross-check."""
    b = {r["result"].get("task"): r for r in recs if r["params"]["run"]["mode"] == "bands"}
    if not b:
        return None
    out: dict[str, Any] = {}
    bulk, cal, sw = b.get("bulk"), b.get("calibration"), b.get("sweep")
    g = None
    if bulk:
        res = bulk["result"]
        g = (res.get("tm") or {}).get("gap") or res.get("gap")
        out["bulk"] = dict(n=res.get("n"), lower=g["lower"], upper=g["upper"], mid=g["mid"], relative_width_pct=g["relative_width_pct"],
                           te_gap_exists=bool(((res.get("te") or {}).get("gap") or {}).get("exists")))
    a_used = [float(r["structure"]["lattice"]["a_nm"]) for r in recs if r["params"]["run"]["mode"] == "harminv"]
    a_used = max(set(a_used), key=a_used.count) if a_used else None
    if cal:
        res = cal["result"]
        a_cal = float(res["a_nm"])
        out["calibration"] = dict(a_calibrated_nm=a_cal, f_mid=res.get("f_mid"), target_nm=res.get("target_wavelength_nm"), a_used_nm=a_used,
                                  rod_radius=res.get("rod_radius"),
                                  rod_radius_nm=(float(res.get("rod_radius") or 0.2) * a_used if a_used else None))
        if g:
            out["window_nm"] = dict(short=a_cal / g["upper"], long=a_cal / g["lower"], width=a_cal / g["lower"] - a_cal / g["upper"])
    if sw:
        pts = sorted([x for x in sw["result"]["sweep"] if x["gap"]["exists"]], key=lambda x: x["n"])
        if pts:
            a, z = pts[0], pts[-1]
            dr = [x for x in pts if design_range[0] - 1e-9 <= x["n"] <= design_range[1] + 1e-9]
            nar = min(dr, key=lambda x: x["gap"]["relative_width_pct"]) if dr else None
            out["sweep"] = dict(n=[a["n"], z["n"]], relative_width_pct=[a["gap"]["relative_width_pct"], z["gap"]["relative_width_pct"]],
                                lower=[a["gap"]["lower"], z["gap"]["lower"]], upper=[a["gap"]["upper"], z["gap"]["upper"]],
                                lower_change_pct=100.0 * (1.0 - z["gap"]["lower"] / a["gap"]["lower"]),
                                upper_change_pct=100.0 * (1.0 - z["gap"]["upper"] / a["gap"]["upper"]),
                                narrowest_in_design_range=(dict(n=nar["n"], relative_width_pct=nar["gap"]["relative_width_pct"]) if nar else None))
    if pwe and g and pwe.get("bulk_gap"):
        pb = pwe["bulk_gap"]
        out["plane_wave_check"] = dict(f_lo=pb["f_lo"], f_hi=pb["f_hi"],
                                       max_relative_difference_pct=100.0 * max(abs(pb["f_lo"] / g["lower"] - 1), abs(pb["f_hi"] / g["upper"] - 1)))
    return out


__all__ = ["census", "repeats", "band_gap", "POINTS", "REPEAT_FACTOR"]
