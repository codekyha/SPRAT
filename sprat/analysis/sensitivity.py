"""Campaign-level series: sensitivity against the analyte index, the defect-radius laws, the raw
barrier series, the cladding ladder and the displacement tolerance.

Selection policy: every series is built from the records that match its parameters; when the
same point was run more than once, ``policy='longest'`` keeps the run with the longest
harmonic-inversion signal (highest Q_lim); ``policy='campaign'`` first excludes the
verification reruns (tag ``verification``, or the imported v4 set) and then keeps the longest
signal, so that the campaign figures of the paper are reproduced from the campaign records.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from ..records import admitted, strongest_mode
from .collect import select

REFERENCE = dict(radius=0.060, row=4, dx=0.0, dy=0.0, resolution=24, cladding_rows=12, guide_periods=25,
                 termination_x="pml", pml=1.0, no_rows=True)



def _r2(y: np.ndarray, residual: np.ndarray) -> float | None:
    """Coefficient of determination; None when the data have no variance (a synthetic or degenerate set)."""
    ss = float(((y - y.mean()) ** 2).sum())
    return float(1 - (residual ** 2).sum() / ss) if ss > 0 else None

def is_verification(r: dict) -> bool:
    """Verification runs: the margin-rule reruns (tag 'verification', or the imported v4 set)."""
    return r["task"].get("source") == "legacy-verification" or r["params"]["run"].get("tag") == "verification"


def _t_used(r: dict) -> float:
    return r["result"].get("t_used") or r["params"]["harminv"]["t"] or 0.0


def _created(r: dict) -> str:
    return r["task"].get("created") or r.get("provenance", {}).get("started") or ""


def _pick(recs: list[dict], key, policy: str = "longest") -> dict:
    """Group by key(record); keep one record per group according to the policy.

    ``longest``: the longest harmonic-inversion signal; ``campaign``: the same after excluding the
    verification reruns; ``first``: the earliest campaign record (the run as first made, before any
    repeat with a longer signal); ``shortest``: the shortest signal (the original of a point that was
    later repeated at a higher margin)."""
    groups: dict = {}
    for r in recs:
        groups.setdefault(key(r), []).append(r)
    out = {}
    for k, lst in groups.items():
        if policy in ("campaign", "first"):
            camp = [r for r in lst if not is_verification(r)]
            lst = camp or lst
        if policy == "first":
            out[k] = min(lst, key=lambda r: (_created(r), _t_used(r)))
        elif policy == "shortest":
            out[k] = min(lst, key=_t_used)
        else:
            out[k] = max(lst, key=_t_used)
    return out


def _linfit(x, y) -> dict[str, Any]:
    x, y = np.asarray(x, float), np.asarray(y, float)
    c = np.polyfit(x, y, 1)
    res = y - np.polyval(c, x)
    return dict(slope=float(c[0]), intercept=float(c[1]), R2=_r2(y, res), rms=float(np.sqrt((res ** 2).mean())))


def _mode_values(r: dict) -> tuple[float, float, float, float]:
    m = strongest_mode(r)
    return float(m["f"]), float(m["Q"]), float(m["wavelength_nm"]), float(r["result"].get("Q_limit") or float("nan"))


# --------------------------------------------------------------------------- analyte sweep
def analyte_sweep(recs: list[dict], policy: str = "campaign", **geometry) -> dict[str, Any]:
    """lambda_r(n), Q(n), local S by central differences, FOM, DL; the reference geometry by default."""
    g = dict(REFERENCE)
    g.update(geometry)
    sel = select(recs, mode="harminv", **g)
    pts = _pick(sel, lambda r: round(float(r["params"]["analyte"]["n"]), 6), policy)
    if len(pts) < 5:
        raise ValueError(f"analyte sweep: only {len(pts)} points found")
    n = np.array(sorted(pts))
    lam = np.array([_mode_values(pts[x])[2] for x in n])
    Q = np.array([_mode_values(pts[x])[1] for x in n])
    Qlim = np.array([_mode_values(pts[x])[3] for x in n])
    margin = Qlim / Q
    n_c = n[1:-1]
    S_c = (lam[2:] - lam[:-2]) / (n[2:] - n[:-2])
    lam_c, Q_c = lam[1:-1], Q[1:-1]
    fwhm_c = lam_c / Q_c
    FOM_c = S_c / fwhm_c
    DL10_c = (fwhm_c / 10.0) / S_c
    DL20_c = (fwhm_c / 20.0) / S_c
    i = int(np.argmin(np.abs(n_c - 1.33)))
    chord = np.polyfit(n, lam, 1)
    resid = lam - np.polyval(chord, n)
    R2 = float(1 - (resid ** 2).sum() / ((lam - lam.mean()) ** 2).sum())
    dSdn = float((S_c[i + 1] - S_c[i - 1]) / (n_c[i + 1] - n_c[i - 1])) if 0 < i < len(n_c) - 1 else float("nan")
    cubic = np.polyfit(n, lam, 3)
    a_nm = float(pts[n[0]]["structure"]["lattice"]["a_nm"])
    out = dict(n=n.tolist(), lambda_nm=lam.tolist(), Q=Q.tolist(), margin=margin.tolist(), labels=[pts[x]["task"]["label"] for x in n],
               n_c=n_c.tolist(), S_c=S_c.tolist(), lambda_c=lam_c.tolist(), Q_c=Q_c.tolist(), fwhm_c=fwhm_c.tolist(),
               FOM_c=FOM_c.tolist(), DL10_c=DL10_c.tolist(), DL20_c=DL20_c.tolist(),
               at_1p33=dict(n=float(n_c[i]), S=float(S_c[i]), Q=float(Q_c[i]), lambda_nm=float(lam_c[i]), fwhm_nm=float(fwhm_c[i]),
                            FOM=float(FOM_c[i]), DL10=float(DL10_c[i]), DL20=float(DL20_c[i]), dS_dn=dSdn),
               chord=dict(S=float(chord[0]), intercept=float(chord[1]), R2=R2, rms_nm=float(np.sqrt((resid ** 2).mean()))),
               cubic_coefficients_descending=cubic.tolist(), Q_drop_factor=float(Q[0] / Q[-1]), a_nm=a_nm)
    # validity window of the sweep: position of each resonance inside the TM gap (MPB edge table) and its margin
    from ..records import gap_edges
    pos = []
    for x in n:
        lo, hi = gap_edges(float(x))
        pos.append((_mode_values(pts[x])[0] - lo) / (hi - lo))
    out["gap_position"] = [float(p) for p in pos]
    out["gap_position_range_pct"] = [100.0 * float(min(pos)), 100.0 * float(max(pos))]
    out["margin_range"] = [float(np.nanmin(margin)), float(np.nanmax(margin))]
    # change of the local sensitivity and of the linewidth across the central-difference range
    out["S_change_pct"] = float(100.0 * (S_c[-1] / S_c[0] - 1.0))
    out["fwhm_factor"] = float(fwhm_c[-1] / fwhm_c[0])
    # the end point of the sweep has no central difference: continue the last local value with the
    # curvature dS/dn_a quoted at n_a = 1.33 (the extrapolation the paper compares the n_a = 1.45 field map with)
    if np.isfinite(dSdn):
        out["S_extrapolated"] = dict(n=float(n[-1]), S=float(S_c[-1] + dSdn * (n[-1] - n_c[-1])), from_n=float(n_c[-1]),
                                     method="S(n_c,last) + (dS/dn_a at 1.33) (n_last - n_c,last)")
    # scaling to the target wavelength (section 2.2): lambda_r(1.33) -> lambda_t; Q, FOM and DL are invariant
    tgt = pts[n[0]]["structure"]["lattice"].get("target_wavelength_nm")
    if tgt:
        k = float(tgt) / float(lam_c[i])
        rod = float(pts[n[0]]["structure"]["lattice"].get("rod_radius", 0.2))
        out["scaled_to_target"] = dict(target_nm=float(tgt), factor=k, a_nm=a_nm * k, rod_radius_nm=rod * a_nm * k,
                                       S=float(S_c[i]) * k, fwhm_nm=float(fwhm_c[i]) * k)
    return out


# --------------------------------------------------------------------------- defect radius
def _split_admitted(pts: dict) -> tuple[dict, dict]:
    """Split a {key: record} selection into the admitted records and those the validity criterion excludes."""
    return ({k: r for k, r in pts.items() if admitted(r)}, {k: r for k, r in pts.items() if not admitted(r)})


def _excluded_block(ex: dict) -> dict[str, Any]:
    ks = sorted(ex)
    return dict(rd=[float(k) for k in ks], f=[_mode_values(ex[k])[0] for k in ks], lambda_nm=[_mode_values(ex[k])[2] for k in ks],
                Q=[_mode_values(ex[k])[1] for k in ks], labels=[ex[k]["task"]["label"] for k in ks])


def radius_sweep(recs: list[dict], policy: str = "campaign", r_min: float = 0.050, r_max: float = 0.115, include_excluded: bool = False,
                 **geometry) -> dict[str, Any]:
    """The fine sweep Q(r_d), lambda_r(r_d) in the reference geometry (PML, 25a guide).  The arrays hold the admitted
    records; a point excluded by the validity criterion (r_d = 0.110a, gap position 0.179) is listed under
    ``excluded`` and enters no fit or statistic."""
    g = dict(REFERENCE)
    g.pop("radius")
    g.update(geometry)
    sel = [r for r in select(recs, mode="harminv", n=1.33, **g) if r_min - 1e-9 <= float(r["structure"]["defect"]["radius"]) <= r_max + 1e-9]
    pts, ex = _split_admitted(_pick(sel, lambda r: round(float(r["structure"]["defect"]["radius"]), 4), policy))
    if include_excluded:                    # every record, as the analysis of manuscript version 8 used them
        pts, ex = {**pts, **ex}, {}
    rd = np.array(sorted(pts))
    if len(rd) < 4:
        raise ValueError(f"radius sweep: only {len(rd)} points found")
    lam = np.array([_mode_values(pts[x])[2] for x in rd])
    Q = np.array([_mode_values(pts[x])[1] for x in rd])
    f = np.array([_mode_values(pts[x])[0] for x in rd])
    margin = np.array([_mode_values(pts[x])[3] for x in rd]) / Q
    a_nm = float(pts[rd[0]]["structure"]["lattice"]["a_nm"])
    d = {round(float(x), 4): (float(l), float(q)) for x, l, q in zip(rd, lam, Q)}
    out: dict[str, Any] = dict(rd=rd.tolist(), f=f.tolist(), lambda_nm=lam.tolist(), Q=Q.tolist(), margin=margin.tolist(),
                               labels=[pts[x]["task"]["label"] for x in rd], a_nm=a_nm, excluded=_excluded_block(ex))
    lin = np.polyfit(rd, lam, 1)
    res = lam - np.polyval(lin, rd)
    out["chord"] = dict(dlambda_drd_nm_per_a=float(lin[0]), dlambda_drd_nm_per_nm=float(lin[0] / a_nm),
                        R2=_r2(lam, res), rms_nm=float(np.sqrt((res ** 2).mean())))
    h = 0.005
    if all(k in d for k in (0.050, 0.055, 0.065, 0.070)):
        dldx = (d[0.050][0] - 8 * d[0.055][0] + 8 * d[0.065][0] - d[0.070][0]) / (12 * h)
        out["local_dlambda_drd"] = dict(nm_per_a=float(dldx), nm_per_nm=float(dldx / a_nm), stencil="five-point, h = 0.005a, at 0.060a")
    if all(k in d for k in (0.055, 0.060, 0.065)):
        f1 = (d[0.065][1] - d[0.055][1]) / (2 * h)
        f2 = (d[0.055][1] - 2 * d[0.060][1] + d[0.065][1]) / h ** 2
        out["local_Q"] = dict(dQ_drd_per_a=float(f1), dQ_drd_per_nm=float(f1 / a_nm),
                              extremum_rd=(float(0.060 - f1 / f2) if f2 != 0 else None),
                              extremum_Q=(float(d[0.060][1] - 0.5 * f1 ** 2 / f2) if f2 != 0 else None))
    qlin = np.polyfit(rd, Q, 1)
    out["Q_linear_R2"] = _r2(Q, Q - np.polyval(qlin, rd))
    peaks = [float(rd[i]) for i in range(1, len(rd) - 1) if Q[i] > Q[i - 1] and Q[i] > Q[i + 1]]
    troughs = [float(rd[i]) for i in range(1, len(rd) - 1) if Q[i] < Q[i - 1] and Q[i] < Q[i + 1]]
    out["peaks"], out["troughs"] = peaks, troughs
    if peaks and [t for t in troughs if t > peaks[0]]:
        out["peak_trough_spacing_a"] = float(min(t for t in troughs if t > peaks[0]) - peaks[0])
    out["peak_trough_ratio"] = float(Q.max() / Q.min())
    if len(peaks) >= 2:
        p1, p2 = peaks[0], peaks[-1]
        out["oscillation"] = dict(period_a=float(p2 - p1), period_nm=float((p2 - p1) * a_nm),
                                  delta_f=float(a_nm / d[round(p1, 4)][0] - a_nm / d[round(p2, 4)][0]))
        out["oscillation"]["n_g_from_period_25a"] = float(1.0 / (2 * 25 * out["oscillation"]["delta_f"])) if out["oscillation"]["delta_f"] else None
    return out


def base_series(recs: list[dict], cladding_rows: int = 8, resolution: int = 20, rows=(2, 3, 4), policy: str = "campaign",
                include_excluded: bool = False) -> dict[str, Any]:
    """The n_cl = 8, res = 20 series of table 1 and figure 3(a,b): f_r(r_d) and Q(r_d) for N_sep = 2, 3, 4.  The arrays
    hold the admitted records; the points the validity criterion excludes (r_d = 0.120a, gap position 0.13-0.15) are
    listed under ``excluded`` of each series and enter no fit or mean."""
    out: dict[str, Any] = dict(cladding_rows=cladding_rows, resolution=resolution, series={})
    all_rd: set[float] = set()
    for N in rows:
        sel = select(recs, mode="harminv", n=1.33, dx=0.0, dy=0.0, row=N, cladding_rows=cladding_rows, resolution=resolution,
                     guide_periods=25, termination_x="pml", no_rows=True)
        pts, ex = _split_admitted(_pick(sel, lambda r: round(float(r["structure"]["defect"]["radius"]), 4), policy))
        if include_excluded:                # every record, as the analysis of manuscript version 8 used them
            pts, ex = {**pts, **ex}, {}
        rd = sorted(pts)
        all_rd.update(rd)
        out["series"][N] = dict(rd=[float(x) for x in rd], f=[_mode_values(pts[x])[0] for x in rd], Q=[_mode_values(pts[x])[1] for x in rd],
                                labels=[pts[x]["task"]["label"] for x in rd], excluded=_excluded_block(ex))
    # quadratic f_r(r_d) law on the N = 4 series
    s4 = out["series"].get(4)
    if s4 and len(s4["rd"]) >= 4:
        rd, f = np.array(s4["rd"]), np.array(s4["f"])
        c2, cov2 = np.polyfit(rd, f, 2, cov="unscaled")
        res2 = f - np.polyval(c2, rd)
        s2 = float((res2 ** 2).sum() / max(len(rd) - 3, 1))
        se2 = np.sqrt(np.diag(cov2) * s2)
        c1 = np.polyfit(rd, f, 1)
        res1 = f - np.polyval(c1, rd)
        out["quadratic_fit"] = dict(coefficients_ascending=c2[::-1].tolist(), standard_errors_ascending=se2[::-1].tolist(),
                                    R2=float(1 - (res2 ** 2).sum() / ((f - f.mean()) ** 2).sum()),
                                    rms=float(np.sqrt((res2 ** 2).mean())), n=int(len(rd)), radius_range=[float(rd.min()), float(rd.max())])
        out["linear_fit"] = dict(coefficients_ascending=c1[::-1].tolist(), R2=float(1 - (res1 ** 2).sum() / ((f - f.mean()) ** 2).sum()),
                                 rms=float(np.sqrt((res1 ** 2).mean())))
        if len(rd) >= 5:
            c3 = np.polyfit(rd, f, 3)
            res3 = f - np.polyval(c3, rd)
            out["cubic_fit"] = dict(coefficients_ascending=c3[::-1].tolist(), R2=_r2(f, res3), rms=float(np.sqrt((res3 ** 2).mean())))
    # the same laws for every row count: how much the intercept and the linear slope depend on N_sep
    per_row = {}
    for N, s in out["series"].items():
        if len(s["rd"]) >= 4:
            q = np.polyfit(s["rd"], s["f"], 2)
            per_row[N] = dict(linear=_linfit(s["rd"], s["f"]), quadratic_coefficients_ascending=q[::-1].tolist())
    if len(per_row) >= 2:
        ic = [v["linear"]["intercept"] for v in per_row.values()]
        sl = [abs(v["linear"]["slope"]) for v in per_row.values()]
        out["fits_by_row"] = per_row
        out["intercept_spread_pct"] = float(100.0 * (max(ic) - min(ic)) / min(ic))
        out["linear_slope_spread_pct"] = float(100.0 * (max(sl) - min(sl)) / min(sl))
        # the same spreads for the quadratic law, the form the paper uses
        qi = [v["quadratic_coefficients_ascending"][0] for v in per_row.values()]
        q2 = [abs(v["quadratic_coefficients_ascending"][2]) for v in per_row.values()]
        out["quadratic_intercept_spread_pct"] = float(100.0 * (max(qi) - min(qi)) / min(qi))
        out["quadratic_curvature_spread_pct"] = float(100.0 * (max(q2) - min(q2)) / min(q2))
    # raw kappa from the common radii of the three rows
    common = sorted(set.intersection(*[set(out["series"][N]["rd"]) for N in rows if N in out["series"]]) if out["series"] else [])
    if common and len(rows) >= 2:
        steps = {}
        lnQ = {N: {rd: math.log(q) for rd, q in zip(out["series"][N]["rd"], out["series"][N]["Q"])} for N in rows}
        for a, b in zip(rows, rows[1:]):
            vals = [lnQ[b][rd] - lnQ[a][rd] for rd in common]
            steps[f"{a}->{b}"] = dict(mean=float(np.mean(vals)), sem=float(np.std(vals, ddof=1) / math.sqrt(len(vals))) if len(vals) > 1 else None,
                                      values=[float(v) for v in vals])
        Ns = np.array([N for N in rows for _ in common], float)
        L = np.array([lnQ[N][rd] for N in rows for rd in common])
        k = float(np.polyfit(Ns, L, 1)[0])
        out["raw_kappa"] = dict(kappa=k, per_row_factor=float(math.exp(k)), steps=steps, common_rd=[float(x) for x in common],
                                mean_Q=[float(np.mean([out["series"][N]["Q"][out["series"][N]["rd"].index(rd)] for rd in common])) for N in rows])
        out["raw_kappa"]["mean_Q_ratio_last_over_first"] = out["raw_kappa"]["mean_Q"][-1] / out["raw_kappa"]["mean_Q"][0]
    return out


def kappa_series(recs: list[dict], radius: float = 0.060, policy: str = "campaign", **geometry) -> dict[str, Any]:
    """Q(N_sep) at fixed radius in the reference cell (the raw, PML series of figure 3c)."""
    g = dict(REFERENCE)
    g.pop("row")
    g["radius"] = radius
    g.update(geometry)
    sel = select(recs, mode="harminv", n=1.33, **g)
    pts = _pick(sel, lambda r: int(r["structure"]["defect"]["row"]), policy)
    N = np.array(sorted(pts), float)
    if len(N) < 2:
        raise ValueError("kappa series: fewer than two rows found")
    Q = np.array([_mode_values(pts[int(x)])[1] for x in N])
    lam = np.array([_mode_values(pts[int(x)])[2] for x in N])
    margin = np.array([_mode_values(pts[int(x)])[3] for x in N]) / Q
    c, cov = np.polyfit(N, np.log(Q), 1, cov=True) if len(N) > 2 else (np.polyfit(N, np.log(Q), 1), np.full((2, 2), np.nan))
    res = np.log(Q) - np.polyval(c, N)
    out = dict(N=[int(x) for x in N], Q=Q.tolist(), lambda_nm=lam.tolist(), margin=margin.tolist(),
               labels=[pts[int(x)]["task"]["label"] for x in N], kappa=float(c[0]), kappa_sd=float(np.sqrt(cov[0, 0])),
               per_row_factor=float(math.exp(c[0])), intercept=float(c[1]),
               R2=float(1 - (res ** 2).sum() / ((np.log(Q) - np.log(Q).mean()) ** 2).sum()) if len(N) > 2 else None,
               steps=np.diff(np.log(Q)).tolist(),
               dlambda_per_row_nm=[float(lam[i + 1] - lam[i]) for i in range(len(N) - 1)])
    # the series as first run (figure 3c of the paper): the earliest campaign record of each point, before any
    # repeat with a longer signal (the N_sep = 5 point was first limited by the ceiling on the sampling time)
    first = _pick(sel, lambda r: int(r["structure"]["defect"]["row"]), "first")
    if len(first) == len(pts):
        Qf = np.array([_mode_values(first[int(x)])[1] for x in N])
        mf = np.array([_mode_values(first[int(x)])[3] for x in N]) / Qf
        cf, covf = np.polyfit(N, np.log(Qf), 1, cov=True) if len(N) > 2 else (np.polyfit(N, np.log(Qf), 1), np.full((2, 2), np.nan))
        out["first_run"] = dict(N=[int(x) for x in N], Q=Qf.tolist(), margin=mf.tolist(), labels=[first[int(x)]["task"]["label"] for x in N],
                                kappa=float(cf[0]), kappa_sd=float(np.sqrt(covf[0, 0])), per_row_factor=float(math.exp(cf[0])),
                                clearance=[int(first[int(x)]["structure"]["cell"]["cladding_rows"]) - int(x) for x in N])
    return out


def resolved_pair_kappa(recs: list[dict], radius: float = 0.100) -> dict[str, Any] | None:
    """The converged pair of the campaign at r_d = 0.100a: N = 2 (n_cl 8, res 20) and N = 4 (n_cl 12, res 32)."""
    q = {}
    for r in select(recs, mode="harminv", n=1.33, dx=0.0, dy=0.0, radius=radius, termination_x="pml", guide_periods=25, no_rows=True):
        s, p = r["structure"], r["params"]
        c = s["cell"]["cladding_rows"] - s["defect"]["row"]
        m = strongest_mode(r)
        ql = r["result"].get("Q_limit")
        if c >= 6 and m and ql and m["Q"] <= ql:
            key = (s["defect"]["row"], s["cell"]["cladding_rows"], p["numerics"]["resolution"])
            q.setdefault(key, float(m["Q"]))
    if (4, 12, 32) in q and (2, 8, 20) in q:
        return dict(Q2=q[(2, 8, 20)], Q4=q[(4, 12, 32)], kappa=float(math.log(q[(4, 12, 32)] / q[(2, 8, 20)]) / 2),
                    Q5=q.get((5, 12, 32)))
    return None


def cladding_ladder(recs: list[dict], radius: float = 0.100, resolution: int = 32) -> dict[str, Any]:
    out: dict[str, Any] = dict(radius=radius, resolution=resolution, series={})
    for N in (4, 5):
        sel = select(recs, mode="harminv", n=1.33, dx=0.0, dy=0.0, radius=radius, row=N, resolution=resolution, termination_x="pml",
                     guide_periods=25, no_rows=True)
        pts = _pick(sel, lambda r: int(r["structure"]["cell"]["cladding_rows"]), "campaign")
        cl = sorted(pts)
        out["series"][N] = dict(cladding_rows=cl, Q=[_mode_values(pts[c])[1] for c in cl], Q_limit=[_mode_values(pts[c])[3] for c in cl],
                                labels=[pts[c]["task"]["label"] for c in cl])
        if cl:
            qc = _mode_values(pts[cl[-1]])[1]
            out["series"][N]["clearance"] = [int(c) - N for c in cl]
            out["series"][N]["deviation_from_thickest_pct"] = [100.0 * (_mode_values(pts[c])[1] / qc - 1.0) for c in cl]
    return out


def _vertex_fixed_quadratic(x, lam) -> dict[str, Any]:
    """lambda = lambda0 - c x^2 (vertex held at x = 0), least squares; and the free quadratic for the vertex position."""
    x, lam = np.asarray(x, float), np.asarray(lam, float)
    A = np.vstack([np.ones_like(x), -x ** 2]).T
    (l0, c), *_ = np.linalg.lstsq(A, lam, rcond=None)
    res = lam - (l0 - c * x ** 2)
    out = dict(c_nm_per_a2=float(c), lambda0_nm=float(l0), rms_nm=float(np.sqrt((res ** 2).mean())))
    if len(x) >= 3:
        q = np.polyfit(x, lam, 2)
        out["free_quadratic_vertex_a"] = float(-q[1] / (2 * q[0])) if q[0] != 0 else None
        out["free_quadratic_c_nm_per_a2"] = float(-q[0])
    return out


def displacement(recs: list[dict], policy: str = "campaign", gain_target: float = 1.30) -> dict[str, Any]:
    """Displacement of the defect rod: Q, lambda_r and f_r against dy (two radii) and dx (two values of dy).

    For dy: the linear fit of ln Q (the coefficient kappa_dy of section 3.4), the best gain and the largest loss,
    the displacement a gain ``gain_target`` would need on the linear law, and the vertex-fixed quadratic
    lambda_r = lambda_0 - c dy^2 with the shift and the frequency change at the largest |dy|.  For dx: the fall of
    Q over the range and the same quadratic in dx."""
    out: dict[str, Any] = dict(dy={}, dx={})
    for rd in (0.05, 0.06):
        sel = select(recs, mode="harminv", n=1.33, dx=0.0, radius=rd, row=4, cladding_rows=12, resolution=24, termination_x="pml",
                     guide_periods=25, no_rows=True)
        pts = _pick(sel, lambda r: round(float(r["structure"]["defect"]["dy"]), 4), policy)
        dys = sorted(pts)
        Q = [_mode_values(pts[d])[1] for d in dys]
        q0 = _mode_values(pts[0.0])[1] if 0.0 in pts else None
        out["dy"][rd] = dict(dy=dys, Q=Q, Q0=q0, ratio=[q / q0 for q in Q] if q0 else None, labels=[pts[d]["task"]["label"] for d in dys])
        e = out["dy"][rd]
        e["lambda_nm"] = [_mode_values(pts[d])[2] for d in dys]
        e["f"] = [_mode_values(pts[d])[0] for d in dys]
        if len(dys) >= 3 and q0:
            fit = _linfit(dys, np.log(Q))
            e["lnQ_fit"] = dict(kappa_dy_per_a=fit["slope"], intercept=fit["intercept"], R2=fit["R2"])
            e["best_gain"] = float(max(e["ratio"]))
            e["dy_at_best_gain"] = float(dys[int(np.argmax(Q))])
            e["worst_ratio"] = float(min(e["ratio"]))
            e["worst_loss_pct"] = float(100.0 * (1.0 - min(e["ratio"])))
            e["dy_for_gain_target"] = dict(gain=gain_target, dy_a=(math.log(gain_target) / fit["slope"]) if fit["slope"] > 0 else None)
            e["lambda_fit"] = _vertex_fixed_quadratic(dys, e["lambda_nm"])
            dmax = max(abs(d) for d in dys)
            ends = [d for d in dys if abs(abs(d) - dmax) < 1e-9]
            i0 = dys.index(0.0)
            e["at_max_abs_dy"] = dict(dy=dmax, lambda_shift_nm=float(np.mean([e["lambda_nm"][dys.index(d)] - e["lambda_nm"][i0] for d in ends])),
                                      f_change=float(np.mean([e["f"][dys.index(d)] - e["f"][i0] for d in ends])),
                                      lambda_shift_each_nm={f"{d:+.3f}": e["lambda_nm"][dys.index(d)] - e["lambda_nm"][i0] for d in ends})
    for dy in (0.0, 0.1):
        sel = select(recs, mode="harminv", n=1.33, dy=dy, radius=0.06, row=4, cladding_rows=12, resolution=24, termination_x="pml",
                     guide_periods=25, no_rows=True)
        pts = _pick(sel, lambda r: round(float(r["structure"]["defect"]["dx"]), 4), policy)
        dxs = sorted(pts)
        Q = [_mode_values(pts[d])[1] for d in dxs]
        q0 = _mode_values(pts[0.0])[1] if 0.0 in pts else None
        out["dx"][dy] = dict(dx=dxs, Q=Q, Q0=q0, ratio=[q / q0 for q in Q] if q0 else None, labels=[pts[d]["task"]["label"] for d in dxs])
        e = out["dx"][dy]
        e["lambda_nm"] = [_mode_values(pts[d])[2] for d in dxs]
        e["f"] = [_mode_values(pts[d])[0] for d in dxs]
        if len(dxs) >= 3 and q0:
            e["fall_over_range_pct"] = float(100.0 * (1.0 - Q[-1] / q0))
            e["lambda_fit"] = _vertex_fixed_quadratic(dxs, e["lambda_nm"])
            e["lambda_shift_at_max_nm"] = float(e["lambda_nm"][-1] - e["lambda_nm"][dxs.index(0.0)])
            # the fall of Q is even in dx (mirror symmetry) and quadratic: Q(dx)/Q(0) - 1 = c dx^2 through the origin,
            # and with a free offset b, which measures the step between the symmetric run at dx = 0 (mirror symmetry
            # exploited) and the displaced runs (no symmetry)
            from ..records import flatten
            xs = np.array([d for d in dxs if d > 0], float)
            rel = np.array([Q[dxs.index(d)] / q0 - 1.0 for d in xs])
            c0_ = float((xs ** 2 @ rel) / (xs ** 4).sum())
            A2 = np.vstack([np.ones_like(xs), xs ** 2]).T
            (b_, c_), *_ = np.linalg.lstsq(A2, rel, rcond=None)
            e["quadratic_fall"] = dict(c_per_a2=c0_, residual_max_pct=float(100.0 * np.abs(rel - c0_ * xs ** 2).max()),
                                       rel_change_pct=[float(100.0 * v) for v in rel], dx=[float(v) for v in xs],
                                       with_offset=dict(offset_pct=float(100.0 * b_), c_per_a2=float(c_),
                                                        residual_max_pct=float(100.0 * np.abs(rel - b_ - c_ * xs ** 2).max())),
                                       symmetry={f"{d:.2f}": flatten(pts[d]).get("symmetry") for d in dxs})
    return out


def preliminary_sweep(recs: list[dict], cladding_rows: int = 6, resolution: int = 20, dx: float = 0.08,
                      dead_zone_radius: float = 0.14) -> dict[str, Any] | None:
    """The preliminary sweep of section 3.4 in a different cell (n_cl = 6, resolution 20): the ratio Q(dx)/Q(0) over
    every (r_d, dy, N_sep) group below the dead zone in which both records satisfy the validity criterion of
    section 2.4, and the dead zone itself (r_d >= 0.14a), where harmonic inversion returns only modes of the finite
    waveguide. The admitted records at N_sep = 1 or r_d >= 0.13a hold such modes too (records.waveguide_mode); the ratio
    over all groups is kept for the superseded statement, and 'cavity_groups' gives it over the cavity groups alone."""
    from ..records import flatten
    sel = select(recs, mode="harminv", n=1.33, cladding_rows=cladding_rows, resolution=resolution, termination_x="pml",
                 guide_periods=25, no_rows=True)
    if not sel:
        return None
    cav = [r for r in sel if float(r["structure"]["defect"]["radius"]) < dead_zone_radius - 1e-9]
    dead = [r for r in sel if float(r["structure"]["defect"]["radius"]) >= dead_zone_radius - 1e-9]
    pts = _pick(cav, lambda r: (round(float(r["structure"]["defect"]["radius"]), 4), round(float(r["structure"]["defect"]["dy"]), 4),
                                int(r["structure"]["defect"]["row"]), round(float(r["structure"]["defect"]["dx"]), 4)), "longest")
    groups: dict = {}
    for (rd, dy, N, x), r in pts.items():
        groups.setdefault((rd, dy, N), {})[x] = (_mode_values(r)[1], flatten(r)["valid"])
    x1 = round(dx, 4)
    ratios = [g[x1][0] / g[0.0][0] for g in groups.values() if 0.0 in g and x1 in g and g[0.0][1] and g[x1][1]]
    out: dict[str, Any] = dict(cladding_rows=cladding_rows, resolution=resolution, records=len(sel), groups=len(ratios), dx=dx,
                               radius_range=[min(k[0] for k in groups), max(k[0] for k in groups)] if groups else None,
                               row_range=[min(k[2] for k in groups), max(k[2] for k in groups)] if groups else None)
    if ratios:
        out.update(median_ratio=float(np.median(ratios)), quartiles=[float(np.percentile(ratios, 25)), float(np.percentile(ratios, 75))])
    if dead:
        out["dead_zone"] = dict(radius_min=dead_zone_radius, records=len(dead), Q_median=float(np.median([_mode_values(r)[1] for r in dead])),
                                lambda_median_nm=float(np.median([_mode_values(r)[2] for r in dead])))
    # the waveguide modes (records.waveguide_mode): N_sep = 1 at every radius and r_d >= 0.13a; the admitted cavity
    # records of the sweep are the others; the dx ratios over the cavity groups alone
    from ..records import waveguide_mode
    flat = [(r, flatten(r)) for r in sel]
    wg = [(r, f) for r, f in flat if f["valid"] and waveguide_mode(f)]
    cv = [(r, f) for r, f in flat if f["valid"] and not waveguide_mode(f)]
    if wg:
        qw = [_mode_values(r)[1] for r, _ in wg]
        out["waveguide_modes"] = dict(admitted=len(wg), Q_min=float(min(qw)), Q_median=float(np.median(qw)), Q_max=float(max(qw)),
                                      lambda_median_nm=float(np.median([_mode_values(r)[2] for r, _ in wg])),
                                      rule="N_sep = 1 or r_d >= 0.13a in the preliminary sweep")
    if cv:
        qc = [_mode_values(r)[1] for r, _ in cv]
        out["cavity_modes"] = dict(admitted=len(cv), Q_min=float(min(qc)), Q_max=float(max(qc)),
                                   lambda_range_nm=[float(min(_mode_values(r)[2] for r, _ in cv)), float(max(_mode_values(r)[2] for r, _ in cv))])
    cav_ratios = [g[x1][0] / g[0.0][0] for (rd, dy, N), g in groups.items()
                  if 0.0 in g and x1 in g and g[0.0][1] and g[x1][1] and N != 1 and rd < 0.13 - 1e-9]
    if cav_ratios:
        out["cavity_groups"] = dict(groups=len(cav_ratios), median_ratio=float(np.median(cav_ratios)),
                                    quartiles=[float(np.percentile(cav_ratios, 25)), float(np.percentile(cav_ratios, 75))])
    return out


def convergence(recs: list[dict]) -> dict[str, Any]:
    """Resolution and cladding convergence tables of the campaign (tables 2 and 3)."""
    out: dict[str, Any] = {}
    sel = select(recs, mode="harminv", n=1.33, dx=0.0, dy=0.0, cladding_rows=8, termination_x="pml", guide_periods=25, no_rows=True)
    rows, modes = {}, {}
    for r in sel:
        s, p = r["structure"], r["params"]
        if p["numerics"]["resolution"] in (20, 32):
            key = (round(float(s["defect"]["radius"]), 4), int(s["defect"]["row"]))
            mv = _mode_values(r)
            rows.setdefault(key, {})[p["numerics"]["resolution"]] = mv[1]
            modes.setdefault(key, {})[p["numerics"]["resolution"]] = mv
    out["resolution_20_vs_32_at_ncl8"] = [dict(radius=k[0], row=k[1], Q20=v.get(20), Q32=v.get(32),
                                               change_pct=(100 * (v[32] / v[20] - 1) if 20 in v and 32 in v else None))
                                          for k, v in sorted(rows.items()) if len(v) == 2]
    ch = [abs(x["change_pct"]) for x in out["resolution_20_vs_32_at_ncl8"] if x["change_pct"] is not None]
    if ch:
        out["resolution_max_change_pct"] = float(max(ch))
    # SPRAT 1.2.1: the resonance frequency of the same pairs: relative change of f_r and the shift
    # of the raw resonance wavelength lambda_r = a / f_r between the two resolutions
    fr = [dict(radius=k[0], row=k[1], f20=v[20][0], f32=v[32][0], lambda20_nm=v[20][2], lambda32_nm=v[32][2],
               f_change_pct=100.0 * (v[32][0] / v[20][0] - 1.0), dlambda_nm=v[32][2] - v[20][2])
          for k, v in sorted(modes.items()) if 20 in v and 32 in v]
    if fr:
        out["resolution_20_vs_32_frequency"] = fr
        out["resolution_f_max_change_pct"] = float(max(abs(x["f_change_pct"]) for x in fr))
        out["resolution_dlambda_max_nm"] = float(max(abs(x["dlambda_nm"]) for x in fr))
    # cladding 10 against 12 at resolution 24 (N_sep = 4): the same rule at the production resolution
    sel = select(recs, mode="harminv", n=1.33, dx=0.0, dy=0.0, row=4, resolution=24, termination_x="pml", guide_periods=25, no_rows=True)
    pts = _pick([r for r in sel if r["structure"]["cell"]["cladding_rows"] in (10, 12)],
                lambda r: (round(float(r["structure"]["defect"]["radius"]), 4), int(r["structure"]["cell"]["cladding_rows"])), "campaign")
    pairs = []
    for rd in sorted({k[0] for k in pts}):
        if (rd, 10) in pts and (rd, 12) in pts:
            q10, q12 = _mode_values(pts[(rd, 10)])[1], _mode_values(pts[(rd, 12)])[1]
            pairs.append(dict(radius=rd, Q10=q10, Q12=q12, change_pct=100.0 * (q10 / q12 - 1.0)))
    out["cladding_10_vs_12_at_res24"] = pairs
    return out


__all__ = ["REFERENCE", "analyte_sweep", "radius_sweep", "base_series", "kappa_series", "resolved_pair_kappa",
           "cladding_ladder", "displacement", "preliminary_sweep", "convergence"]
