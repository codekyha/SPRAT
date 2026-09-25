"""Transmission spectra under the two guide terminations: the absorber comparison of the paper (H14b).

The spectra of the reference geometry in the 25a cell, whose guide ends in the perfectly matched layer (PML), are
compared with spectra of the same geometry whose guide is continued into an 8a adiabatic absorber; every other
parameter of each pair is identical (gate G1). Both are fitted with the Fano routine the criterion names,

    T(lambda) = T_bg + T0 (q + e)^2 / (1 + e^2),   e = 2 (lambda - lambda_r) / fwhm,

over the whole frequency window, each spectrum normalised by the cavity-less reference run of its own termination.

Criterion D-22, fixed on 24 September 2026 before the runs were submitted:
    factor = |q|_PML / |q|_absorber at each analyte index;
    factor >= 3 supports the attribution of the residual asymmetry to the reflections at the guide ends there,
    factor <= 1.5 refutes it there, anything in between is inconclusive;
    pass = supports at all seven indices, every gate passed and S1 satisfied.
Amendment A1, fixed on 24 September 2026 at 15:46:39 (+03:00), before the first absorber cavity spectrum was
written (15:53:07): the ratio is graded only at n_a = 1.325 to 1.425, where |q|_PML exceeds 0.2; at 1.300 and
1.450 the PML asymmetry is close to a zero of q, which follows the Fabry-Perot phase, and the factor is reported
without a grade. The statistic, the thresholds, the gates at all seven indices and S1 are unchanged.
Amendment A2, fixed with A1 (a measurement, not a test): the cell factor G(n_a) = Q_PML / Q_absorber, both from the
fit, and the reflectionless quality factor Q_w = Q_PML,harminv / G, at every index whose two rows pass every gate,
provided S1 holds. Q_PML,harminv is the q_est of the PML spectrum run, the harmonic-inversion Q of the analyte sweep at
that index, rounded. Reported with it: G at n_a = 1.33 from harmonic inversion (the reference record over the
absorber record), and the opening of the linewidth over the sweep in the 25a cell and without reflections.

Gates (a row that fails one is reported and never graded):
    G1  structure and parameters of the two spectra and of the two references identical apart from the label and
        the termination
    G2  identical frequency grids (relative tolerance 1e-12)
    G3  the fitted resonance at least five linewidths inside both window edges
    G4  fit rms <= 0.05, lambda_r inside the window, 0 < fwhm < window / 5
    G5  t_end / tau >= 3, tau = Q / (pi f_r)
    S1  (set) the absorber Q interpolated to n_a = 1.33 from 1.325 and 1.350 within 5 % of the harmonic-inversion Q
        of the same absorber cell (the record v4_H4b_s4, 4495.1)

The fit and the gates reproduce the scripts ``48_fano_karsilastir_S4b.py`` and ``49_H14b_A1A2.py`` that graded the
runs (both in ``raw_h14b.tar.gz`` of the deposit); the fit is kept line for line, including its starting values.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy.optimize import least_squares

from .collect import best, select

INDICES = (1.300, 1.325, 1.350, 1.375, 1.400, 1.425, 1.450)
A1_GRADED = (1.325, 1.350, 1.375, 1.400, 1.425)
SUPPORTS, REFUTES = 3.0, 1.5
S1_REFERENCE_Q, S1_TOLERANCE = 4495.1, 0.05
EXCLUDED_KEYS = {"params.run.label", "structure.cell.termination_x", "structure.cell.absorber_periods"}
CRITERION = ("factor = |q|_PML / |q|_absorber; >= 3 supports the guide-end attribution at that index, <= 1.5 refutes it "
             "there, otherwise inconclusive; pass = supports at all seven indices, every gate passed, S1 satisfied "
             "(D-22, fixed 2026-09-24 before the runs were submitted)")
AMENDMENT = ("A1: the ratio is graded at n_a = 1.325 to 1.425 only, where |q|_PML exceeds 0.2; 1.300 and 1.450 are reported, "
             "not graded; statistic, thresholds, gates and S1 unchanged. A2: G(n_a) = Q_PML / Q_absorber (fit) and "
             "Q_w = Q_PML,harminv / G at every index whose rows pass every gate, provided S1 holds. Fixed 2026-09-24 "
             "15:46:39 +03:00, before the first absorber cavity spectrum was written (15:53:07)")


# --------------------------------------------------------------------------- the fit the criterion names
def fano(lam, T_bg, T0, q, lam_r, dlam):
    eps = 2.0 * (lam - lam_r) / dlam
    return T_bg + T0 * (q + eps) ** 2 / (1.0 + eps ** 2)


def fit(lam: np.ndarray, T: np.ndarray) -> dict[str, float]:
    """The Fano fit of 48_fano_karsilastir_S4b.py, line for line (whole window, least squares)."""
    i = int(np.argmin(T))
    lam0 = lam[i]
    half = 0.5 * (np.median(T) + T[i])
    try:
        left = np.where(T[:i] > half)[0][-1]
        right = i + np.where(T[i:] > half)[0][0]
        d0 = abs(lam[right] - lam[left])
    except Exception:                                         # noqa: BLE001 (as in the script)
        d0 = abs(lam[-1] - lam[0]) / 50.0
    p0 = [float(np.median(T)), float(np.median(T)), 0.05, float(lam0), float(max(d0, 1e-4))]
    r = least_squares(lambda p: fano(lam, *p) - T, p0, x_scale=[1, 1, 1, abs(lam0), max(d0, 1e-4)], max_nfev=20000)
    T_bg, T0, q, lam_r, dlam = r.x
    return dict(T_bg=float(T_bg), T0=float(T0), q=float(q), lambda_r=float(lam_r), FWHM_nm=float(abs(dlam)),
                rms=float(np.sqrt(np.mean(r.fun ** 2))), T_min=float(T.min()))


# --------------------------------------------------------------------------- records
def _flat(d: dict, prefix: str = "") -> dict[str, Any]:
    out: dict[str, Any] = {}
    for k, v in d.items():
        key = f"{prefix}.{k}" if prefix else str(k)
        if isinstance(v, dict):
            out.update(_flat(v, key))
        else:
            out[key] = v
    return out


def identity(a: dict, b: dict) -> list[str]:
    """Differences between two records in structure and parameters, apart from the label and the termination (G1)."""
    fa = _flat(dict(structure=a["structure"], params=a["params"]))
    fb = _flat(dict(structure=b["structure"], params=b["params"]))
    diff = []
    for k in sorted(set(fa) | set(fb)):
        if k in EXCLUDED_KEYS:
            continue
        va, vb = fa.get(k, "<missing>"), fb.get(k, "<missing>")
        if isinstance(va, (int, float)) and isinstance(vb, (int, float)) and not isinstance(va, bool) and not isinstance(vb, bool):
            same = float(va) == float(vb)
        else:
            same = va == vb
        if not same:
            diff.append(f"{k}: PML {va!r}, absorber {vb!r}")
    if b["structure"]["cell"]["termination_x"] != "absorber":
        diff.append("termination of the second run is not the absorber")
    return diff


def _pair(recs: list[dict], n: float, termination: str) -> tuple[dict | None, dict | None]:
    crit = dict(n=n, radius=0.06, row=4, dx=0.0, dy=0.0, cladding_rows=12, resolution=24, guide_periods=25,
                termination_x=termination)
    sp = select(recs, mode="spectrum", **crit)
    ref = select(recs, mode="reference", n=n, cladding_rows=12, resolution=24, guide_periods=25, termination_x=termination)
    return (sp[0] if len(sp) == 1 else None), (ref[0] if len(ref) == 1 else None)


def _row(recs: list[dict], n: float) -> dict[str, Any]:
    r: dict[str, Any] = dict(n=n, gates={}, notes=[])
    a_s, a_n = _pair(recs, n, "pml")
    e_s, e_n = _pair(recs, n, "absorber")
    missing = [name for name, x in (("PML spectrum", a_s), ("PML reference", a_n), ("absorber spectrum", e_s),
                                    ("absorber reference", e_n)) if x is None]
    if missing:
        r["status"] = "missing: " + ", ".join(missing)
        return r
    r["labels"] = dict(pml=a_s["task"]["label"], pml_reference=a_n["task"]["label"], absorber=e_s["task"]["label"],
                       absorber_reference=e_n["task"]["label"])
    g1 = identity(a_s, e_s) + ["reference: " + x for x in identity(a_n, e_n)]
    r["gates"]["G1"] = not g1
    r["notes"] += g1
    fa = np.asarray(a_s["result"]["f"], float)
    grids = [np.asarray(x["result"]["f"], float) for x in (a_n, e_s, e_n)]
    g2 = all(g.shape == fa.shape and np.allclose(fa, g, rtol=1e-12, atol=0.0) for g in grids)
    r["gates"]["G2"] = g2
    if not g2:
        r["notes"].append("frequency grids differ")
    if not (r["gates"]["G1"] and g2):
        r["status"] = "failed G1/G2 (protocol)"
        return r
    a_nm = float(a_s["structure"]["lattice"]["a_nm"])
    lam = a_nm / fa
    lo, hi = float(lam.min()), float(lam.max())
    for side, sp, ref in (("pml", a_s, a_n), ("absorber", e_s, e_n)):
        P = np.asarray(sp["result"]["flux_out"], float)
        Pn = np.asarray(ref["result"]["flux_out"], float)
        T = P / np.where(Pn == 0, np.nan, Pn)
        ok = np.isfinite(T)
        x = fit(lam[ok], T[ok])
        Q = x["lambda_r"] / x["FWHM_nm"] if x["FWHM_nm"] > 0 else float("nan")
        tau = Q / (math.pi * a_nm / x["lambda_r"]) if Q == Q else float("nan")
        near = ok & (np.abs(lam - x["lambda_r"]) < 5 * x["FWHM_nm"])
        t_end = float(sp["result"]["t_end"])
        x.update(Q=Q, t_end=t_end, t_end_over_tau=float(t_end / tau) if tau == tau and tau > 0 else float("nan"),
                 T_max_near=float(np.nanmax(T[near])) if near.any() else float("nan"),
                 edge_margin_linewidths=float(min(x["lambda_r"] - lo, hi - x["lambda_r"]) / x["FWHM_nm"]) if x["FWHM_nm"] > 0 else float("nan"),
                 q_harminv_run=float(sp["params"]["run"]["q_est"]))
        g3 = x["edge_margin_linewidths"] >= 5.0
        g4 = (x["rms"] <= 0.05) and (lo < x["lambda_r"] < hi) and (0.0 < x["FWHM_nm"] < (hi - lo) / 5.0)
        g5 = x["t_end_over_tau"] >= 3.0
        x["gates"] = dict(G3=bool(g3), G4=bool(g4), G5=bool(g5))
        r[side] = x
        for g, v in x["gates"].items():
            r["gates"][f"{g}_{side}"] = v
    if not all(r["gates"].values()):
        r["status"] = "failed " + ", ".join(k for k, v in r["gates"].items() if not v)
        return r
    qa, qe = abs(r["pml"]["q"]), abs(r["absorber"]["q"])
    r["factor"] = (qa / qe) if qe > 0 else float("inf")
    r["grade"] = "supports" if r["factor"] >= SUPPORTS else ("refutes" if r["factor"] <= REFUTES else "inconclusive")
    r["status"] = "graded"
    return r


def _verdict(rows: dict[float, dict], graded_keys, s1) -> dict[str, Any]:
    complete = all(rows[n].get("status") == "graded" for n in rows) and s1 is True
    g = [n for n in graded_keys if rows.get(n, {}).get("status") == "graded"]
    ref = [n for n in g if rows[n]["grade"] == "refutes"]
    inc = [n for n in g if rows[n]["grade"] == "inconclusive"]
    if not g:
        text = "no graded index"
    elif ref:
        text = "refuted at n_a = " + ", ".join("%.3f" % n for n in ref)
    elif not inc and len(g) == len(graded_keys) and complete:
        text = "pass"
    elif not inc:
        text = "supports at every graded index; not a pass: set incomplete"
    else:
        text = "not a pass: inconclusive at n_a = " + ", ".join("%.3f" % n for n in inc)
    return dict(complete=complete, graded=g, refutes=ref, inconclusive=inc, verdict=text)


def analyse(recs: list[dict]) -> dict[str, Any] | None:
    """The comparison of the two terminations; None when the records hold no absorber spectrum."""
    if not select(recs, mode="spectrum", termination_x="absorber"):
        return None
    rows = {n: _row(recs, n) for n in INDICES}
    h4b = best(recs, mode="harminv", radius=0.06, row=4, cladding_rows=12, resolution=24, termination_x="absorber",
               guide_periods=25, n=1.33, dx=0.0, dy=0.0, no_rows=True)
    ref = best(recs, mode="harminv", radius=0.06, row=4, cladding_rows=12, resolution=24, termination_x="pml",
               guide_periods=25, n=1.33, dx=0.0, dy=0.0, no_rows=True)
    from ..records import strongest_mode
    q_h4b = q_ref = None
    if h4b is not None:
        m = strongest_mode(h4b)
        q_h4b = float(m["Q"]) if m else None
    if ref is not None:
        m = strongest_mode(ref)
        q_ref = float(m["Q"]) if m else None
    s1 = None
    s1_value = None
    g_rows = {n: r for n, r in rows.items() if r.get("status") == "graded"}
    if 1.325 in g_rows and 1.350 in g_rows:
        q133 = g_rows[1.325]["absorber"]["Q"] + (1.33 - 1.325) / 0.025 * (g_rows[1.350]["absorber"]["Q"] - g_rows[1.325]["absorber"]["Q"])
        s1_value = q133
        s1 = abs(q133 / S1_REFERENCE_Q - 1.0) <= S1_TOLERANCE
    d22 = _verdict(rows, INDICES, s1)
    a1 = _verdict(rows, A1_GRADED, s1)
    # A2: the cell factor and the reflectionless quality factor
    a2: dict[str, Any] = dict(evaluated=bool(s1))
    if s1:
        tab = {}
        for n, r in rows.items():
            if r.get("status") != "graded":
                continue
            G = r["pml"]["Q"] / r["absorber"]["Q"]
            qh = r["pml"]["q_harminv_run"]
            tab[n] = dict(G=G, Q_pml_fit=r["pml"]["Q"], Q_absorber_fit=r["absorber"]["Q"], Q_pml_harminv=qh, Q_w=qh / G,
                          fwhm_w_nm=r["pml"]["lambda_r"] / (qh / G), fwhm_cell_nm=r["pml"]["lambda_r"] / qh)
        a2["rows"] = tab
        if 1.325 in tab and 1.350 in tab:
            a2["G_1p33"] = tab[1.325]["G"] + (1.33 - 1.325) / 0.025 * (tab[1.350]["G"] - tab[1.325]["G"])
        if q_ref and q_h4b:
            a2.update(G_1p33_harminv=q_ref / q_h4b, Q_pml_harminv_1p33=q_ref, Q_absorber_harminv_1p33=q_h4b)
        if 1.300 in tab and 1.450 in tab:
            a2["Q_w_fall"] = tab[1.300]["Q_w"] / tab[1.450]["Q_w"]
            a2["fwhm_w_opening"] = tab[1.450]["fwhm_w_nm"] / tab[1.300]["fwhm_w_nm"]
            a2["Q_cell_fall_harminv"] = tab[1.300]["Q_pml_harminv"] / tab[1.450]["Q_pml_harminv"]
            a2["fwhm_cell_opening"] = tab[1.450]["fwhm_cell_nm"] / tab[1.300]["fwhm_cell_nm"]
    graded = [r for r in rows.values() if r.get("status") == "graded"]
    summary = {}
    if graded:
        inner = [rows[n] for n in A1_GRADED if rows[n].get("status") == "graded"]
        summary = dict(
            abs_q_absorber_max=max(abs(r["absorber"]["q"]) for r in graded),
            abs_q_absorber_min=min(abs(r["absorber"]["q"]) for r in graded),
            abs_q_pml_inner_min=min(abs(r["pml"]["q"]) for r in inner) if inner else None,
            factor_inner_min=min(r["factor"] for r in inner) if inner else None,
            factor_inner_max=max(r["factor"] for r in inner) if inner else None,
            T_max_near_absorber_max=max(r["absorber"]["T_max_near"] for r in graded),
            T_max_near_absorber_min=min(r["absorber"]["T_max_near"] for r in graded),
            T_max_near_pml_max=max(r["pml"]["T_max_near"] for r in graded),
            t_end_over_tau_absorber=[min(r["absorber"]["t_end_over_tau"] for r in graded),
                                     max(r["absorber"]["t_end_over_tau"] for r in graded)],
            rms_absorber_max=max(r["absorber"]["rms"] for r in graded),
            edge_margin_absorber_min=min(r["absorber"]["edge_margin_linewidths"] for r in graded))
    return dict(criterion=CRITERION, amendment=AMENDMENT, indices=list(INDICES), a1_graded=list(A1_GRADED),
                rows=[rows[n] for n in INDICES], S1=s1, S1_value=s1_value, S1_reference=S1_REFERENCE_Q,
                harminv_absorber_Q_1p33=q_h4b, harminv_pml_Q_1p33=q_ref, D22=d22, A1=a1, A2=a2, summary=summary)


def absorber_spectra(recs: list[dict]) -> list[dict[str, Any]]:
    """The absorber-terminated spectra normalised as the comparison normalises them (SPRAT 1.2.1, figure 6 inset).

    One entry per analyte index whose absorber spectrum and absorber reference are both present: the wavelength on
    the raw lattice constant (ascending) and T = flux_out / flux_out of the cavity-less absorber reference."""
    out = []
    for n in INDICES:
        sp, ref = _pair(recs, n, "absorber")
        if sp is None or ref is None:
            continue
        f = np.asarray(sp["result"]["f"], float)
        P = np.asarray(sp["result"]["flux_out"], float)
        Pn = np.asarray(ref["result"]["flux_out"], float)
        if Pn.shape != P.shape:
            continue
        lam = float(sp["structure"]["lattice"]["a_nm"]) / f
        order = np.argsort(lam)
        T = P / np.where(Pn == 0, np.nan, Pn)
        out.append(dict(n=n, label=sp["task"]["label"], reference_label=ref["task"]["label"],
                        lambda_nm=lam[order].tolist(), T=T[order].tolist()))
    return out


__all__ = ["fano", "fit", "identity", "analyse", "absorber_spectra", "INDICES", "A1_GRADED", "CRITERION", "AMENDMENT"]
