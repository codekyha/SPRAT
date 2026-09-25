"""Grading of the frozen prediction sets against the records.

Every criterion is evaluated from records selected by their parameters (never by name), so
the grader works on the imported 2026 deposit and on a regenerated data set alike.
Deviations are reported, never reconciled: a FAIL stays a FAIL and gets an explanation in
the report.
"""

from __future__ import annotations

import json
import math
import os
import time
from typing import Any

from .. import pwe as pwe_mod
from ..records import strongest_mode
from .collect import Q_of, best, f_of, select
from .sensitivity import _pick

CRITERIA_FILE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "predictions_criteria.json")


def load_criteria(path: str | None = None) -> dict:
    with open(path or CRITERIA_FILE, encoding="utf-8") as fh:
        return json.load(fh)


def _margin(r):
    m = strongest_mode(r)
    ql = r["result"].get("Q_limit")
    return (ql / m["Q"]) if (m and ql) else None


def _fmt(x, nd=2):
    return "-" if x is None else (f"{x:.{nd}f}" if isinstance(x, float) else str(x))


def _fit_lnQ(pairs):
    n = len(pairs)
    if n < 2:
        return None, None
    xs = [p[0] for p in pairs]
    ys = [math.log(p[1]) for p in pairs]
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    kap = sxy / sxx
    a = my - kap * mx
    if n > 2:
        s2 = sum((y - a - kap * x) ** 2 for x, y in zip(xs, ys)) / (n - 2)
        return kap, math.sqrt(s2 / sxx)
    return kap, None


class Grader:
    """Holds the records, the plane-wave result and the analysis outputs the criteria refer to."""

    REF = dict(mode="harminv", n=1.33, dx=0.0, dy=0.0, resolution=24, guide_periods=25, pml=1.0, no_rows=True)

    def __init__(self, recs: list[dict], pwe: dict | None, barrier: dict | None = None, fields: dict | None = None,
                 min_margin: float = 4.0):
        self.recs, self.pwe, self.barrier, self.fields, self.min_margin = recs, pwe, barrier, fields, min_margin
        self.rows: list[dict] = []

    # ---- helpers
    def rec(self, **crit):
        return best(self.recs, **dict(self.REF, **crit))

    def verified(self, **crit):
        """The longest-signal record with margin >= min_margin (the rule of the reruns, fixed before them)."""
        cands = [r for r in select(self.recs, **dict(self.REF, **crit)) if (_margin(r) or 0) >= self.min_margin]
        if not cands:
            return None
        return max(cands, key=lambda r: r["result"].get("t_used") or 0)

    def add(self, cid, quantity, expected, measured, verdict, note=""):
        self.rows.append(dict(id=cid, quantity=quantity, expected=expected, measured=measured, verdict=verdict, note=note))

    # ---- the criteria
    def grade(self, criteria: dict) -> list[dict]:
        A = criteria["sets"]["first"]["archive_constants"]
        for c in criteria["criteria"]:
            cid = c["id"]
            try:
                getattr(self, "c_" + cid.replace("-", "_"), None)(c, A) if hasattr(self, "c_" + cid) else self.c_default(c)
            except Exception as e:                       # noqa: BLE001 - a broken criterion must not stop the report
                self.add(cid, c["quantity"], json.dumps(c.get("expected")), f"error: {type(e).__name__}: {e}", "ERROR")
        return self.rows

    def c_default(self, c):
        if c.get("closed"):
            self.add(c["id"], c["quantity"], str(c["expected"]), "closed without the study", "CLOSED",
                     "see predictions/frozen, PB_closure; the closure there rests on a bound from the stored map that version 9 "
                     "withdraws (the stored map is not the permittivity E_z sees); the Hellmann-Feynman estimate on the same field, "
                     "consistent with the discretisation, reproduces the measured sensitivity to 1.5 %, and the resolution study "
                     "remains unrun")
        else:
            self.add(c["id"], c["quantity"], json.dumps(c["expected"]), "not evaluated", "INFO", c.get("note", ""))

    def c_P1(self, c, A):
        self._pwe_channels(c, "reference")

    def c_P2(self, c, A):
        self._pwe_channels(c, "second")

    def _pwe_channels(self, c, which):
        if not self.pwe:
            return self.add(c["id"], c["quantity"], str(c["expected"]), "no pwe record", "MISSING")
        k = pwe_mod.kappa_block(self.pwe, which)
        lo, hi = c["expected"]
        # The criterion registered the two values of the campaign's plane-wave layer.  A layer written by SPRAT 1.1 or
        # later keeps physical roots only and reports the slower campaign value as the slowest truncation artefact of
        # the square basis; that value is what the first half of the criterion reproduces.
        slow = k.get("kappa_boundary_artefact") or k["kappa_pred"]
        ok = abs(slow / lo - 1) < c["tolerance_rel"] and abs(k["kappa_pred_zone_edge"] / hi - 1) < c["tolerance_rel"]
        note = c.get("note", "")
        if k.get("kappa_boundary_artefact"):
            note = (note + "; " if note else "") + ("the slower value is a truncation artefact of the square plane-wave basis, not a Bloch "
                                                    "channel; the physical channel is the zone-edge root (supplementary section S6)")
        self.add(c["id"], c["quantity"], f"{lo} / {hi} (+- {100 * c['tolerance_rel']:.0f} %)",
                 f"{slow:.4f} / {k['kappa_pred_zone_edge']:.4f}", "PASS" if ok else "FAIL", note)

    def c_P3(self, c, A):
        k6 = self.verified(radius=0.06, row=6, cladding_rows=14, termination_x="pml")
        k5 = self.verified(radius=0.06, row=5, cladding_rows=13, termination_x="pml") or self.verified(radius=0.06, row=5, cladding_rows=12, termination_x="pml")
        if not (k6 and k5):
            return self.add(c["id"], c["quantity"], str(c["expected"]), "-", "MISSING")
        r = Q_of(k6) / Q_of(k5)
        lo, hi = c["expected"]
        self.add(c["id"], c["quantity"], f"{lo} to {hi}",
                 f"Q6 = {Q_of(k6):.0f} (margin {_margin(k6):.2f}), Q5 = {Q_of(k5):.0f} [{k5['task']['label']}], ratio {r:.3f}, ln {math.log(r):.3f}",
                 "PASS" if lo <= r <= hi else "FAIL",
                 "the clearance-8 leak carries a larger share of 1/Q at N = 6 than at N = 5; the leak-removed step is in the barrier analysis")

    def c_P4(self, c, A):
        k = self.verified(radius=0.06, row=5, cladding_rows=12, termination_x="pml")
        self._range(c, k, extra=lambda k: f"step 4->5: {math.log(Q_of(k) / A['Q4_cl12']):.3f}")

    def c_P13(self, c, A):
        k = self.verified(radius=0.06, row=5, cladding_rows=13, termination_x="pml")
        self._range(c, k, extra=lambda k: f"step 4->5: {math.log(Q_of(k) / A['Q4_cl12']):.3f}")

    def _range(self, c, k, extra=None):
        if not k:
            return self.add(c["id"], c["quantity"], str(c["expected"]), "-", "MISSING")
        lo, hi = c["expected"]
        q = Q_of(k)
        note = extra(k) if extra else ""
        self.add(c["id"], c["quantity"], f"{lo} to {hi}", f"{q:.0f}, margin {_margin(k):.2f} [{k['task']['label']}]", "PASS" if lo <= q <= hi else "FAIL", note)

    def c_P4b(self, c, A):
        k = self.verified(radius=0.06, row=4, cladding_rows=12, termination_x="pml")
        if not k:
            return self.add(c["id"], c["quantity"], f"{c['expected']} +- 1 %", "-", "MISSING")
        d = Q_of(k) / c["expected"] - 1
        self.add(c["id"], c["quantity"], f"{c['expected']} +- {100 * c['tolerance_rel']:.0f} %",
                 f"{Q_of(k):.1f} ({100 * d:+.2f} %), margin {_margin(k):.2f}", "PASS" if abs(d) <= c["tolerance_rel"] else "FAIL")

    def c_P5(self, c, A):
        if not self.pwe:
            return self.add(c["id"], c["quantity"], str(c["expected"]), "no pwe record", "MISSING")
        ng = pwe_mod.kappa_block(self.pwe, "reference")["n_g"]
        ngm = self.pwe["fabry_perot"]["ng_mean"]
        D = self.pwe["fabry_perot"]["D_over_a"]
        lo, hi = c["expected"]
        ok = abs(ng - lo) <= c["tolerance_abs"] and abs(ngm - hi) <= c["tolerance_abs"]
        self.add(c["id"], c["quantity"], f"{lo}, {hi} +- {c['tolerance_abs']}", f"n_g = {ng:.3f}, <n_g> = {ngm:.3f}, D = {D:.2f} a (W1 end 12.5 a)",
                 "PASS" if ok else "FAIL", "plane-wave layer")

    def c_P6(self, c, A):
        sel = select(self.recs, **dict(self.REF, guide_periods=49, row=4, cladding_rows=12, termination_x="pml"))
        pts = sorted((float(r["structure"]["defect"]["radius"]), f_of(r), Q_of(r)) for r in _pick(sel, lambda r: round(float(r["structure"]["defect"]["radius"]), 4)).values())
        if len(pts) < 5:
            return self.add(c["id"], c["quantity"], f"{c['expected']} +- 10 %", f"{len(pts)} points", "MISSING")
        peaks = [pts[i] for i in range(1, len(pts) - 1) if pts[i][2] > pts[i - 1][2] and pts[i][2] > pts[i + 1][2]]
        if len(peaks) < 2:
            return self.add(c["id"], c["quantity"], f"{c['expected']} +- 10 %", f"{len(peaks)} peak(s) over {len(pts)} points", "UNDETERMINED",
                            "one peak: the range is too short or there is no oscillation")
        peaks = sorted(peaks, key=lambda p: -p[2])[:2]
        df = abs(peaks[0][1] - peaks[1][1])
        d = df / c["expected"] - 1
        self.add(c["id"], c["quantity"], f"{c['expected']} +- {100 * c['tolerance_rel']:.0f} %",
                 f"peaks at r_d {peaks[0][0]:.4f}, {peaks[1][0]:.4f}; delta_f = {df:.5f} ({100 * d:+.1f} %)", "PASS" if abs(d) <= c["tolerance_rel"] else "FAIL")

    def c_P7(self, c, A):
        sel = select(self.recs, **dict(self.REF, row=4, cladding_rows=12, termination_x="absorber"))
        ab = _pick(sel, lambda r: round(float(r["structure"]["defect"]["radius"]), 4))
        if len(ab) < 4:
            self.add(c["id"], "absorber peak/trough", f"< {c['expected']['peak_trough_max']}", f"{len(ab)} points", "MISSING")
        else:
            qs = [Q_of(v) for v in ab.values()]
            ratio = max(qs) / min(qs)
            self.add(c["id"], f"absorber peak/trough over {len(ab)} r_d points", f"< {c['expected']['peak_trough_max']}",
                     f"{ratio:.2f} (max {max(qs):.0f}, min {min(qs):.0f}; PML {A['Q_peak_rd0.065'] / A['Q_trough_rd0.085']:.2f})",
                     "PASS" if ratio < c["expected"]["peak_trough_max"] else "FAIL")
        ser = [(n, self.verified(radius=0.06, row=n, cladding_rows=12, termination_x="absorber")) for n in (3, 4, 5)]
        ser = [(n, Q_of(k)) for n, k in ser if k]
        if len(ser) >= 2:
            kap, err = _fit_lnQ(ser)
            d = kap / A["kappa_paper"] - 1
            self.add(c["id"], "kappa with the absorber (N_sep " + ",".join(str(n) for n, _ in ser) + ")", f"{A['kappa_paper']} +- {100 * c['expected']['kappa_change_max']:.0f} %",
                     f"{kap:.4f}" + (f" +- {err:.4f}" if err else "") + f" ({100 * d:+.2f} %); steps " + ", ".join(f"{math.log(ser[i + 1][1] / ser[i][1]):.3f}" for i in range(len(ser) - 1)),
                     "PASS" if abs(d) <= c["expected"]["kappa_change_max"] else "FAIL",
                     "raw absorber values at mixed clearances; the leak-removed series is in the barrier analysis")
            k10 = self.verified(radius=0.06, row=5, cladding_rows=13, termination_x="absorber")
            if k10:
                ser2 = sorted([s for s in ser if s[0] != 5] + [(5, Q_of(k10))])
                kap2, err2 = _fit_lnQ(ser2)
                self.add(c["id"], "kappa with the absorber, N_sep = 5 at n_cl = 13", "1.90 to 1.96 (channels)", f"{kap2:.4f}" + (f" +- {err2:.4f}" if err2 else ""), "INFO")
        else:
            self.add(c["id"], "kappa with the absorber", "-", "-", "MISSING")

    def c_P8(self, c, A):
        if not self.fields or not self.fields.get("reference"):
            return self.add(c["id"], c["quantity"], str(c["expected"]), "no field map at n = 1.33", "MISSING")
        e = self.fields["reference"]
        b = c["expected"]
        ok = b["eta_a"][0] <= e["eta_lin"] <= b["eta_a"][1] and abs(e["S_lin"] / b["S_th_nm_per_RIU"] - 1) <= b["tolerance_rel"]
        self.add(c["id"], c["quantity"], f"eta_a {b['eta_a']}, S_th {b['S_th_nm_per_RIU']} +- {100 * b['tolerance_rel']:.0f} %",
                 f"eta_mask {e['eta_mask']:.4f}, eta_lin {e['eta_lin']:.4f}, eta_pure {e['eta_pure']:.4f}; S_lin {e['S_lin']:.1f} nm/RIU ({100 * (e['S_lin'] / b['S_th_nm_per_RIU'] - 1):+.1f} %)",
                 "PASS" if ok else "FAIL",
                 ("graded with the linear unmixing, as registered; the unmixing reads the stored permittivity map as the arithmetic "
                  "mean of the two permittivities, which it is not (Meep writes the harmonic mean of the eigenvalues of the "
                  "smoothed tensor); the estimator consistent with the discretisation gives eta_a = %.4f and S_th = %.1f nm/RIU"
                  % (e["eta_hf"], e["S_hf"])) if e.get("eta_hf") else "graded with the linear unmixing, as registered")

    def c_P9(self, c, A):
        for n, key in ((4, "N4_max_rel"), (5, "N5_max_rel")):
            k = self.rec(radius=0.06, row=n, cladding_rows=12, termination_x="pml", pml=2.0)
            r = self.verified(radius=0.06, row=n, cladding_rows=12, termination_x="pml")
            if not (k and r):
                self.add(c["id"], f"PML 2a against 1a (N_sep = {n})", f"< {100 * c['expected'][key]:.0f} %", "-", "MISSING")
                continue
            d = Q_of(k) / Q_of(r) - 1
            self.add(c["id"], f"PML 2a against 1a (N_sep = {n})", f"< {100 * c['expected'][key]:.0f} %", f"{Q_of(k):.0f} against {Q_of(r):.0f} ({100 * d:+.2f} %)",
                     "PASS" if abs(d) <= c["expected"][key] else "FAIL")

    def c_P10(self, c, A):
        ref = self.verified(radius=0.06, row=4, cladding_rows=12, termination_x="pml")
        if not ref:
            return self.add(c["id"], c["quantity"], "-", "no reference", "MISSING")
        for rb, band in c["expected"].items():
            cands = [r for r in select(self.recs, mode="harminv", n=1.33, dx=0.0, dy=0.0, resolution=24, guide_periods=25, pml=1.0,
                                       radius=0.06, row=4, cladding_rows=12, termination_x="pml")
                     if {int(k): float(v) for k, v in (r["structure"].get("rows") or {}).items()} == {2: float(rb)}]
            if not cands:
                self.add(c["id"], f"Q(r_b = {rb}a)/Q_ref", f"{band}", "-", "MISSING")
                continue
            k = max(cands, key=lambda r: r["result"].get("t_used") or 0)
            r = Q_of(k) / Q_of(ref)
            mid = math.log(math.sqrt(band[0] * band[1]))
            yb = (math.log(band[1]) - math.log(band[0])) / 2
            tol = c["tolerance_ln_rel"] * abs(mid) + yb
            self.add(c["id"], f"Q(r_b = {rb}a)/Q_ref", f"{band[0]} to {band[1]} (+- 30 % in ln)",
                     f"{r:.3f} (ln {math.log(r):+.3f}; expected {mid:+.3f} +- {tol:.3f})", "PASS" if abs(math.log(r) - mid) <= tol else "FAIL")

    def c_P12(self, c, A):
        k = self.verified(radius=0.06, row=4, cladding_rows=12, termination_x="absorber")
        if not k:
            return self.add(c["id"], c["quantity"], str(c["expected"]), "-", "MISSING")
        lo, hi = c["expected"]
        q = Q_of(k)
        self.add(c["id"], c["quantity"], f"{lo} to {hi} (model {c['model_range']})", f"{q:.0f} = {q / 6927:.3f} x 6927, margin {_margin(k):.2f}",
                 "PASS" if lo <= q <= hi else "FAIL", "PASS confirms the Fabry-Perot reading: the absolute Q of the 25a cell contains the cell factor")

    def c_H10(self, c, A):
        k = self.verified(radius=0.06, row=5, cladding_rows=13, termination_x="absorber")
        self._range(c, k)

    def c_H8(self, c, A):
        k = self.rec(radius=0.06, row=5, cladding_rows=12, termination_x="pml", resolution=32)
        r = self.verified(radius=0.06, row=5, cladding_rows=12, termination_x="pml")
        if not (k and r):
            return self.add(c["id"], c["quantity"], f"< {100 * c['expected']['max_rel']:.0f} %", "-", "MISSING")
        d = Q_of(k) / Q_of(r) - 1
        self.add(c["id"], c["quantity"], f"< {100 * c['expected']['max_rel']:.0f} %", f"res 32 {Q_of(k):.0f} against res 24 {Q_of(r):.0f} ({100 * d:+.2f} %)",
                 "PASS" if abs(d) <= c["expected"]["max_rel"] else "FAIL")

    # ---- second set (needs the barrier analysis)
    def _second(self):
        return (self.barrier or {}).get("second_radius")

    def c_PA1(self, c, A):
        s = self._second()
        if not s:
            return self.add(c["id"], c["quantity"], str(c["expected"]), "second-radius series missing", "MISSING")
        lo, hi = c["expected"]
        self.add(c["id"], c["quantity"], f"[{lo}, {hi}]", f"{s['kappa_absorber']:.4f}", "PASS" if lo <= s["kappa_absorber"] <= hi else "FAIL", c.get("note", ""))

    def c_PA2(self, c, A):
        s = self._second()
        if not s:
            return self.add(c["id"], c["quantity"], "< 0.02", "-", "MISSING")
        d = s["kappa_difference"]
        self.add(c["id"], c["quantity"], f"< {c['expected']['max']}", f"{d:.4f}", "PASS" if d < c["expected"]["max"] else "FAIL",
                 "the difference equals the drift of the cell factor G with N (d lnG per row %.4f and %.4f), which the criterion assumed absent; "
                 "the absorber value carries no such term" % (s["dlnG_per_row"]["3->4"], s["dlnG_per_row"]["4->5"]))

    def c_PA3(self, c, A):
        s = self._second()
        if not s:
            return self.add(c["id"], c["quantity"], f"{c['expected']} +- {c['tolerance_abs']}", "-", "MISSING")
        g = s["G"]["4"]
        self.add(c["id"], c["quantity"], f"{c['expected']} +- {c['tolerance_abs']}", f"{g:.4f}", "PASS" if abs(g - c["expected"]) <= c["tolerance_abs"] else "FAIL")

    def c_PA4(self, c, A):
        k = self.verified(radius=0.10, row=4, cladding_rows=12, termination_x="pml")
        self._rel(c, k)

    def c_PA5(self, c, A):
        k = self.verified(radius=0.10, row=4, cladding_rows=12, termination_x="absorber")
        self._rel(c, k)

    def _rel(self, c, k):
        if not k:
            return self.add(c["id"], c["quantity"], f"{c['expected']} +- {100 * c['tolerance_rel']:.0f} %", "-", "MISSING")
        d = Q_of(k) / c["expected"] - 1
        self.add(c["id"], c["quantity"], f"{c['expected']} +- {100 * c['tolerance_rel']:.0f} %", f"{Q_of(k):.2f} ({100 * d:+.3f} %)",
                 "PASS" if abs(d) <= c["tolerance_rel"] else "FAIL")

    def c_PA6(self, c, A):
        s = self._second()
        kref = (self.barrier or {}).get("series", {}).get("kappa")
        if not s or kref is None:
            return self.add(c["id"], c["quantity"], "true", "-", "MISSING")
        self.add(c["id"], c["quantity"], f"kappa(0.100a) < {kref:.4f}", f"{s['kappa_absorber']:.4f} < {kref:.4f}", "PASS" if s["kappa_absorber"] < kref else "FAIL")


def grade(recs: list[dict], pwe: dict | None, barrier: dict | None = None, fields: dict | None = None,
          criteria_path: str | None = None) -> dict[str, Any]:
    criteria = load_criteria(criteria_path)
    g = Grader(recs, pwe, barrier, fields)
    rows = g.grade(criteria)
    verdicts = {}
    for r in rows:
        verdicts[r["verdict"]] = verdicts.get(r["verdict"], 0) + 1
    if criteria_path:
        from .pipeline import portable_path
        crit_name = portable_path(criteria_path)
    else:
        crit_name = "sprat/data/predictions_criteria.json"      # the packaged criteria, named relative to the package
    return dict(generated=time.strftime("%Y-%m-%d %H:%M:%S"), criteria_file=crit_name,
                frozen=criteria["sets"], rows=rows, verdict_counts=verdicts,
                rule="Every FAIL is reported first and explained afterwards; no criterion is edited after the results are seen.")


def report_markdown(g: dict) -> str:
    lines = ["# Criteria fixed before the runs: grading report", "", f"generated {g['generated']}", "",
             "Frozen sets: " + "; ".join(f"{k}: {v['frozen_file']} (sha256 {v['sha256'][:12]}..., frozen {v['frozen_on']})" for k, v in g["frozen"].items()), "",
             "| id | quantity | expected | measured | verdict | note |", "|---|---|---|---|---|---|"]
    for r in g["rows"]:
        lines.append(f"| {r['id']} | {r['quantity']} | {r['expected']} | {r['measured']} | **{r['verdict']}** | {r['note']} |")
    lines += ["", "Verdicts: " + ", ".join(f"{k} {v}" for k, v in sorted(g["verdict_counts"].items())), "", g["rule"], ""]
    return "\n".join(lines)


def grade_cli(records_dir: str, out: str, criteria_path: str | None = None) -> None:
    from . import pipeline
    res = pipeline.analyze(records_dir, out, predictions_only=True, criteria_path=criteria_path)
    print(report_markdown(res["predictions"]))


__all__ = ["CRITERIA_FILE", "load_criteria", "Grader", "grade", "report_markdown", "grade_cli"]
