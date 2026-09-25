"""The barrier analysis of the verification runs: leak removal, cell factor, the reflectionless
leak-free series and the decay constant kappa, with the Fabry-Perot fits and the error budget.

Model.  For a defect N_sep rows from the guide with clearance c = n_cl - N_sep to the cladding
edge, in a cell of guide termination T (PML or absorber),

    1/Q = e^{-kappa N} / (A G_T) + 1/Q_top(c),      Q_top(c) = Q_top(8) e^{kappa_top (c - 8)},

where G_PML is the Fabry-Perot cell factor of the PML-terminated 25a guide and G_abs = 1.
Three PML points at clearance 8 (N = 4, 5, 6) determine kappa, Q_top(8) and A exactly (the
three-point solve); the n_cl = 12 / 13 pair at N = 5 gives Q_top(7) and kappa_top; every
record is then leak-corrected, the PML records are divided by the directly measured
G = Q_PML / Q_abs at the same N, and the four-point series Q_w(N), N = 3..6, gives kappa.
"""

from __future__ import annotations

import math
import random
from typing import Any

import numpy as np
from scipy.optimize import least_squares

from .. import pwe as pwe_mod
from ..records import admitted, strongest_mode
from .collect import Q_of, best, f_of, select
from .sensitivity import _pick, radius_sweep

MC_SIGMA_Q = 0.0023          # independent discretisation error per Q (0.23 %)
MC_SIGMA_G6 = 0.0011         # bound on the drift of G between N = 5 and 6 (0.11 %)
MC_DRAWS = 20000
MC_SEED = 7


def _q(rec):
    return Q_of(rec)


def _lab(rec):
    return rec["task"]["label"] if rec else None


def _margin(rec):
    m = strongest_mode(rec)
    ql = rec["result"].get("Q_limit")
    return (ql / m["Q"]) if (m and ql) else None


# --------------------------------------------------------------------------- inputs
def inputs(recs: list[dict], radius: float = 0.060) -> dict[str, Any]:
    """The records of the chain, selected by parameters (longest signal wins)."""
    base = dict(mode="harminv", n=1.33, dx=0.0, dy=0.0, resolution=24, guide_periods=25, pml=1.0, no_rows=True, radius=radius)
    P = dict(base, termination_x="pml")
    A = dict(base, termination_x="absorber")
    sel = {
        "Q3_pml_cl12": best(recs, row=3, cladding_rows=12, **P), "Q4_pml_cl12": best(recs, row=4, cladding_rows=12, **P),
        "Q5_pml_cl12": best(recs, row=5, cladding_rows=12, **P), "Q5_pml_cl13": best(recs, row=5, cladding_rows=13, **P),
        "Q6_pml_cl14": best(recs, row=6, cladding_rows=14, **P),
        "A3_abs_cl12": best(recs, row=3, cladding_rows=12, **A), "A4_abs_cl12": best(recs, row=4, cladding_rows=12, **A),
        "A5_abs_cl12": best(recs, row=5, cladding_rows=12, **A), "A5_abs_cl13": best(recs, row=5, cladding_rows=13, **A),
        "Q4_pml_nx49": best(recs, mode="harminv", n=1.33, dx=0.0, dy=0.0, resolution=24, guide_periods=49, pml=1.0, no_rows=True,
                            radius=radius, row=4, cladding_rows=12, termination_x="pml"),
        "Q5_pml_cl12_res32": best(recs, mode="harminv", n=1.33, dx=0.0, dy=0.0, resolution=32, guide_periods=25, pml=1.0,
                                  no_rows=True, radius=radius, row=5, cladding_rows=12, termination_x="pml"),
        "Q4_pml_cl12_pml2": best(recs, mode="harminv", n=1.33, dx=0.0, dy=0.0, resolution=24, guide_periods=25, pml=2.0,
                                 no_rows=True, radius=radius, row=4, cladding_rows=12, termination_x="pml"),
        "Q5_pml_cl12_pml2": best(recs, mode="harminv", n=1.33, dx=0.0, dy=0.0, resolution=24, guide_periods=25, pml=2.0,
                                 no_rows=True, radius=radius, row=5, cladding_rows=12, termination_x="pml"),
    }
    return {k: dict(label=_lab(v), Q=_q(v), f=f_of(v), margin=_margin(v), clearance=(v["structure"]["cell"]["cladding_rows"] - v["structure"]["defect"]["row"]) if v else None,
                    termination=(v["structure"]["cell"]["termination_x"] if v else None),
                    host=(v["provenance"].get("host", "") if v else ""), job_id=(v["provenance"].get("job_id", "") if v else ""),
                    code_sha256=(v["software"].get("code_sha256", "") if v else ""), record=v) for k, v in sel.items()}


def absorber_sweep(recs: list[dict], row: int = 4, cladding_rows: int = 12, policy: str = "longest",
                   excluded: bool | None = False) -> dict[float, dict]:
    """Absorber record at every radius, N_sep = 4, n_cl = 12, res 24: the longest signal by default, or with
    ``policy='shortest'`` the original of a point that was later repeated at a higher margin.  Only records admitted
    by the validity criterion are returned; ``excluded=True`` returns the excluded ones instead (r_d = 0.110a) and
    ``excluded=None`` every record."""
    sel = select(recs, mode="harminv", n=1.33, dx=0.0, dy=0.0, resolution=24, guide_periods=25, pml=1.0, no_rows=True,
                 termination_x="absorber", row=row, cladding_rows=cladding_rows)
    pts = _pick(sel, lambda r: round(float(r["structure"]["defect"]["radius"]), 4), policy)
    if excluded is None:                    # every record, as the analysis of manuscript version 8 used them
        return pts
    return {k: r for k, r in pts.items() if admitted(r) != excluded}


def series_records(recs: list[dict], points: list[dict]) -> list[dict | None]:
    """One record per point of a series: prefer the record from the job (provenance.job_id) that ran the most points of
    the series, so that a series is one set of runs; the longest signal decides within a job.  A point that was also
    run in another sweep is thus taken from the series it belongs to."""
    cands = [select(recs, **c) for c in points]
    cover: dict = {}
    for cs in cands:
        for j in {(r.get("provenance", {}).get("job_id") or "") for r in cs}:
            if j:
                cover[j] = cover.get(j, 0) + 1
    return [max(cs, key=lambda r: (cover.get(r.get("provenance", {}).get("job_id") or "", 0), r["result"].get("t_used") or 0.0)) if cs else None
            for cs in cands]


def guide49_sweep(recs: list[dict]) -> dict[float, dict]:
    sel = select(recs, mode="harminv", n=1.33, dx=0.0, dy=0.0, resolution=24, guide_periods=49, pml=1.0, no_rows=True,
                 termination_x="pml", row=4, cladding_rows=12)
    return _pick(sel, lambda r: round(float(r["structure"]["defect"]["radius"]), 4), "longest")


# --------------------------------------------------------------------------- the chain
def chain(q4, q5c8, q5c7, q6, a3, a4, a5c7, a5c8, dlnG6: float = 0.0) -> dict[str, Any] | None:
    """Three-point solve -> leak law -> G -> four-point kappa.  None if the inputs are inconsistent."""
    i4, i5, i6 = 1 / q4, 1 / q5c8, 1 / q6
    ek = (i4 - i5) / (i5 - i6)
    if ek <= 1:
        return None
    k3 = math.log(ek)
    it8 = i6 - (i5 - i6) / (ek - 1)
    if it8 <= 0:
        return None
    d = 1 / q5c7 - 1 / q5c8
    if d <= 0:
        return None
    kt = math.log((it8 + d) / it8)
    amp = 1 / ((i5 - i6) / (math.exp(-5 * k3) - math.exp(-6 * k3)))

    def dl(q, c):
        return 1 / (1 / q - it8 * math.exp(-kt * (c - 8)))

    G4 = dl(q4, 8) / dl(a4, 8)
    G5 = dl(q5c8, 8) / dl(a5c8, 8)
    G5b = dl(q5c7, 7) / dl(a5c7, 7)
    G6 = G5 * math.exp(dlnG6)
    Qw = {3: dl(a3, 9), 4: dl(a4, 8), 5: dl(a5c8, 8), 6: dl(q6, 8) / G6}
    N = np.array([3.0, 4.0, 5.0, 6.0])
    L = np.log([Qw[n] for n in (3, 4, 5, 6)])
    k, b = np.polyfit(N, L, 1)
    rr = L - (k * N + b)
    se_k = float(math.sqrt((rr @ rr) / (len(N) - 2) / ((N - N.mean()) ** 2).sum()))
    return dict(kappa=float(k), kappa_se_ols=se_k, A=float(math.exp(b)), lnA=float(b), A_three_point=float(amp),
                Q_top_8=1 / it8, inv_Q_top_8=it8, Q_top_7=1 / (it8 + d), inv_Q_top_7=it8 + d, kappa_top=kt, kappa_three_point=k3,
                G={"4": G4, "5": G5, "5_cl12": G5b}, Q_w=Qw, steps=[math.log(Qw[n + 1] / Qw[n]) for n in (3, 4, 5)],
                residuals_lnQ=[float(x) for x in rr], dl=dl,
                corrected_pml_steps={"4->5": math.log((1 / (i5 - it8)) / (1 / (i4 - it8))), "5->6": math.log((1 / (i6 - it8)) / (1 / (i5 - it8)))},
                raw_pml_steps={"4->5": math.log(q5c8 / q4), "5->6": math.log(q6 / q5c8)},
                leak_share_pct={"N4_c8": 100 * it8 / i4, "N5_c8": 100 * it8 / i5, "N6_c8": 100 * it8 / i6, "N5_c7": 100 * (it8 + d) * q5c7})


def monte_carlo(q4, q5c8, q5c7, q6, a3, a4, a5c7, a5c8, draws: int = MC_DRAWS, seed: int = MC_SEED) -> dict[str, Any]:
    """Error budget: independent 0.23 % per Q plus the 0.11 % bound on G(N = 6); the same draw order as the 2026 analysis."""
    rng = random.Random(seed)
    ks, qts, qws = [], [], []
    for _ in range(draws):
        p = lambda: 1 + rng.gauss(0, MC_SIGMA_Q)         # noqa: E731
        r = chain(q4 * p(), q5c8 * p(), q5c7 * p(), q6 * p(), a3 * p(), a4 * p(), a5c7 * p(), a5c8 * p(), rng.gauss(0, MC_SIGMA_G6))
        if r:
            ks.append(r["kappa"])
            qts.append(r["Q_top_8"])
            qws.append(r["Q_w"][4])
    ks, qts, qws = np.array(ks), np.array(qts), np.array(qws)
    return dict(draws=int(len(ks)), kappa_mean=float(ks.mean()), kappa_sd=float(ks.std()), Q_top_8_mean=float(qts.mean()),
                Q_top_8_sd=float(qts.std()), Q_w4_mean=float(qws.mean()), Q_w4_sd=float(qws.std()),
                sigma_Q=MC_SIGMA_Q, sigma_lnG6=MC_SIGMA_G6, seed=seed)


# --------------------------------------------------------------------------- Fabry-Perot
def airy(p, beta_arr):
    R, D, phi = p[0], p[1], p[2]
    return (1 + R - 2 * math.sqrt(max(R, 1e-12)) * np.cos(2 * beta_arr * D + phi)) / (1 - R)


def matched_ratio_fit(RD, GG, BET) -> dict[str, Any]:
    """Two-reflector fit of the measured cell factor G(r_d) = Q_PML / Q_abs against beta(f_r)."""
    RD, GG, BET = map(np.asarray, (RD, GG, BET))
    sol = least_squares(lambda p: np.log(airy(p, BET)) - np.log(GG), [0.078, 12.5, 4.4], bounds=([1e-4, 8, -20], [0.6, 20, 20]))
    res = sol.fun
    cov = np.linalg.pinv(sol.jac.T @ sol.jac) * (res @ res / max(len(res) - 3, 1))
    se = np.sqrt(np.diag(cov))
    scl = least_squares(lambda p: np.log(airy(p[:3], BET) * p[3]) - np.log(GG), [0.078, 12.5, 4.4, 1.0],
                        bounds=([1e-4, 8, -20, 0.5], [0.6, 20, 20, 2.0]))
    sesc = np.sqrt(np.diag(np.linalg.pinv(scl.jac.T @ scl.jac) * (scl.fun @ scl.fun / max(len(scl.fun) - 4, 1))))
    return dict(points=int(len(RD)), rd=RD.tolist(), G=GG.tolist(), beta_over_pi=(BET / math.pi).tolist(),
                R=float(sol.x[0]), R_se=float(se[0]), D=float(sol.x[1]), D_se=float(se[1]), phi=float(sol.x[2]), phi_se=float(se[2]),
                rms_lnG=float(math.sqrt((res @ res) / len(res))), residuals=[float(x) for x in res],
                scale_c=float(scl.x[3]), scale_c_se=float(sesc[3]), R_with_scale=float(scl.x[0]), D_with_scale=float(scl.x[1]),
                rms_lnG_with_scale=float(math.sqrt((scl.fun @ scl.fun) / len(scl.fun))))


def shape_fits(rds, fs, Qs, beta, D0: float = 12.5, f_ref: float = 0.30461) -> dict[str, Any]:
    """The four Fabry-Perot shape fits of the PML sweep Q(r_d): two-mirror / one-mirror x D fixed / D free.

    ln Q = a + c1 (f - f_ref) - ln F(phi),  phi = 2 beta(f) D + phi_r,
    two mirrors: F = (1 - R) / (1 + R - 2 sqrt(R) cos phi);  one mirror: F = 1 + r cos phi.
    """
    rds, fs, Qs = map(np.asarray, (rds, fs, Qs))
    lnQ = np.log(Qs)
    b = np.array([beta(f) for f in fs])

    def F2(R, phi):
        return (1 - R) / (1 + R - 2 * math.sqrt(max(R, 1e-12)) * np.cos(phi))

    def F1(r, phi):
        return 1 + r * np.cos(phi)

    def fit(model, D_free):
        bestfit = None
        for phi0 in (0.0, math.pi / 2, math.pi, 3 * math.pi / 2):
            if model == "two_mirror":
                def res(p):
                    a, c1, R, phr = p[:4]
                    D = p[4] if D_free else D0
                    return a + c1 * (fs - f_ref) - np.log(F2(R, 2 * b * D + phr)) - lnQ
                x0 = [math.log(Qs.mean()), 0.0, 0.1, phi0] + ([D0] if D_free else [])
                lo = [-np.inf, -np.inf, 0.0, -2 * math.pi] + ([6.0] if D_free else [])
                hi = [np.inf, np.inf, 0.7, 4 * math.pi] + ([20.0] if D_free else [])
            else:
                def res(p):
                    a, c1, r, phr = p[:4]
                    D = p[4] if D_free else D0
                    return a + c1 * (fs - f_ref) - np.log(F1(r, 2 * b * D + phr)) - lnQ
                x0 = [math.log(Qs.mean()), 0.0, 0.3, phi0] + ([D0] if D_free else [])
                lo = [-np.inf, -np.inf, 0.0, -2 * math.pi] + ([6.0] if D_free else [])
                hi = [np.inf, np.inf, 0.95, 4 * math.pi] + ([20.0] if D_free else [])
            try:
                s = least_squares(res, x0, bounds=(lo, hi))
            except Exception:           # noqa: BLE001
                continue
            rms = float(np.sqrt(np.mean(s.fun ** 2)))
            if bestfit is None or rms < bestfit["rms_lnQ"]:
                p = s.x
                out = dict(model=model, D_free=D_free, a=float(p[0]), c1=float(p[1]), rms_lnQ=rms, phi_r=float(p[3] % (2 * math.pi)),
                           D=float(p[4]) if D_free else D0, residuals=[float(v) for v in s.fun])
                if model == "two_mirror":
                    R = float(p[2])
                    out.update(R=R, peak_trough_predicted=((1 + math.sqrt(R)) / (1 - math.sqrt(R))) ** 2)
                else:
                    r = float(p[2])
                    out.update(r=r, R=r * r, peak_trough_predicted=(1 + r) / (1 - r))
                bestfit = out
        return bestfit

    return dict(points=int(len(rds)), rd=rds.tolist(), f=fs.tolist(), Q=Qs.tolist(), beta_over_pi=(b / math.pi).tolist(),
                peak_trough_measured=float(Qs.max() / Qs.min()), D_fixed=D0, f_ref=f_ref,
                fits={"two_mirror_D12.5": fit("two_mirror", False), "two_mirror_D_free": fit("two_mirror", True),
                      "one_mirror_D12.5": fit("one_mirror", False), "one_mirror_D_free": fit("one_mirror", True)})


# --------------------------------------------------------------------------- global fit
def global_fit(points: list[dict], G_of_f, k3: float, amp: float, inv_top8: float, kt: float) -> dict[str, Any]:
    """Four-parameter fit (ln A, kappa, ln 1/Q_top(8), kappa_top) to all nine records with the fitted G(f) for PML."""
    def resid(p):
        lnA, kap, li8, ktop = p
        out = []
        for pt in points:
            G = G_of_f(pt["f"]) if pt["termination"] == "pml" else 1.0
            inv = math.exp(-(lnA + kap * pt["N"])) / G + math.exp(li8) * math.exp(-ktop * (pt["clearance"] - 8))
            out.append(math.log(inv) - math.log(1 / pt["Q"]))
        return np.array(out)
    p0 = [math.log(amp), k3, math.log(inv_top8), kt]
    sol = least_squares(resid, p0)
    res = sol.fun
    n, npar = len(res), 4
    cov = np.linalg.pinv(sol.jac.T @ sol.jac) * (res @ res / max(n - npar, 1))
    se = np.sqrt(np.diag(cov))
    lnA, kap, li8, ktop = sol.x
    return dict(kappa=float(kap), kappa_se=float(se[1]), Q_top_8=float(math.exp(-li8)), Q_top_se_pct=float(100 * se[2]),
                kappa_top=float(ktop), kappa_top_se=float(se[3]), Q_w={"4": float(math.exp(lnA + 4 * kap)), "5": float(math.exp(lnA + 5 * kap)),
                                                                           "6": float(math.exp(lnA + 6 * kap))},
                rms_lnQ=float(np.sqrt((res ** 2).mean())), max_abs_residual=float(np.abs(res).max()),
                residuals={pt["label"]: float(x) for pt, x in zip(points, res)})


# --------------------------------------------------------------------------- everything
def analyse(recs: list[dict], pwe: dict, field_eta: dict | None = None, S_measured: float | None = None,
            radius: float = 0.060, second_radius: float = 0.100) -> dict[str, Any]:
    inp = inputs(recs, radius)
    need = ["Q4_pml_cl12", "Q5_pml_cl13", "Q5_pml_cl12", "Q6_pml_cl14", "A3_abs_cl12", "A4_abs_cl12", "A5_abs_cl12", "A5_abs_cl13"]
    missing = [k for k in need if inp[k]["record"] is None]
    if missing:
        raise ValueError("barrier analysis needs the verification records " + ", ".join(missing))
    Q4, Q5c8, Q5c7, Q6 = (inp[k]["Q"] for k in ("Q4_pml_cl12", "Q5_pml_cl13", "Q5_pml_cl12", "Q6_pml_cl14"))
    A3, A4, A5c7, A5c8 = (inp[k]["Q"] for k in ("A3_abs_cl12", "A4_abs_cl12", "A5_abs_cl12", "A5_abs_cl13"))
    c0 = chain(Q4, Q5c8, Q5c7, Q6, A3, A4, A5c7, A5c8)
    if c0 is None:
        raise ValueError("the three-point solve is inconsistent (no positive leak)")
    dl = c0.pop("dl")
    mc = monte_carlo(Q4, Q5c8, Q5c7, Q6, A3, A4, A5c7, A5c8)
    KP = pwe_mod.kappa_block(pwe, "reference")
    k_slow, k_edge = KP["kappa_pred"], KP["kappa_pred_zone_edge"]
    out: dict[str, Any] = dict(
        inputs={k: {kk: vv for kk, vv in v.items() if kk != "record"} for k, v in inp.items()},
        three_point_solve=dict(kappa=c0["kappa_three_point"], Q_top_8=c0["Q_top_8"], A=c0["A_three_point"],
                               raw_pml_steps=c0["raw_pml_steps"], corrected_pml_steps=c0["corrected_pml_steps"],
                               dev_zone_edge_pct=100 * (c0["kappa_three_point"] / k_edge - 1), dev_slowest_pct=100 * (c0["kappa_three_point"] / k_slow - 1)),
        leak=dict(Q_top_8=c0["Q_top_8"], Q_top_7=c0["Q_top_7"], kappa_top=c0["kappa_top"], leak_share_pct=c0["leak_share_pct"],
                  dinvQ_pml_N5=1 / Q5c7 - 1 / Q5c8, dinvQ_abs_N5=1 / A5c7 - 1 / A5c8,
                  ratio_abs_over_pml=(1 / A5c7 - 1 / A5c8) / (1 / Q5c7 - 1 / Q5c8)),
        G=c0["G"],
        series=dict(Q_w={str(k): v for k, v in c0["Q_w"].items()}, steps=c0["steps"], kappa=c0["kappa"], kappa_se_ols=c0["kappa_se_ols"],
                    kappa_sd_monte_carlo=mc["kappa_sd"], A=c0["A"], per_row_factor=math.exp(c0["kappa"]),
                    residuals_lnQ=c0["residuals_lnQ"], dev_zone_edge_pct=100 * (c0["kappa"] / k_edge - 1),
                    dev_slowest_pct=100 * (c0["kappa"] / k_slow - 1), kappa_absorber_only=0.5 * (c0["steps"][0] + c0["steps"][1])),
        monte_carlo=mc, channels=dict(slowest=k_slow, zone_edge=k_edge, n_g=KP.get("n_g"), next=KP.get("kappa_next"),
                                      re_ky_next_over_pi=KP.get("re_ky_next_over_pi"),
                                      boundary_artefact=KP.get("kappa_boundary_artefact"),
                                      kx0=KP.get("kappa_at_kx0"), kxpi=KP.get("kappa_at_kxpi"),
                                      boundary_artefact_kx0=KP.get("kappa_boundary_artefact_kx0"),
                                      boundary_artefact_kxpi=KP.get("kappa_boundary_artefact_kxpi"),
                                      truncation_check=pwe.get("truncation_check")))
    # ---- quality of the four-point series and agreement of the three routes (supplementary section S6)
    s_ = out["series"]
    s_["max_abs_residual_lnQ"] = float(max(abs(x) for x in c0["residuals_lnQ"]))
    s_["range_factor"] = c0["Q_w"][6] / c0["Q_w"][3]
    ratios = [math.exp(x) for x in c0["steps"]]
    s_["per_row_ratios"] = ratios
    s_["per_row_ratio_spread_pct"] = 100.0 * (max(ratios) - min(ratios)) / min(ratios)
    s_["routes"] = dict(four_point=c0["kappa"], three_point=c0["kappa_three_point"], absorber_only=s_["kappa_absorber_only"])
    s_["routes_spread_pct"] = 100.0 * (max(s_["routes"].values()) - min(s_["routes"].values())) / c0["kappa"]
    out["G_spread_pct"] = 100.0 * (max(c0["G"].values()) - min(c0["G"].values())) / min(c0["G"].values())
    # SPRAT 1.2.1: the channel is evaluated at the guided-mode wavevector k_x = beta because the
    # cavity couples to the guide only through the Fourier component of its tail that is phase-matched to the guided
    # mode; the same frequency at k_x = 0 (or pi/a) would predict another per-row factor
    ch_ = out["channels"]
    tc_ = ch_.get("truncation_check") or {}
    if tc_.get("kx_over_pi") is not None:
        ch_["beta_over_pi"] = float(tc_["kx_over_pi"])
    if ch_.get("kx0") is not None:
        ch_["per_row_factor_kx0"] = math.exp(ch_["kx0"])
    if ch_.get("kxpi") is not None:
        ch_["per_row_factor_kxpi"] = math.exp(ch_["kxpi"])
    ch_["per_row_factor_zone_edge"] = math.exp(k_edge)
    out["leak"]["Q_top_8_rel_sd_pct"] = 100.0 * mc["Q_top_8_sd"] / c0["Q_top_8"]
    out["reference_campaign"] = reference_campaign_record(recs, radius)
    # ---- Fabry-Perot: the PML fine sweep (campaign policy) and its shape fits
    fw = np.array(pwe["W1"]["f"])
    kx = np.array(pwe["W1"]["kx_over_pi"])
    beta = lambda f: math.pi * float(np.interp(f, fw, kx))           # noqa: E731 (linear, as in the 2026 fit)
    try:
        sweep = radius_sweep(recs, policy="campaign", r_min=0.045, r_max=0.115)
    except ValueError:
        sweep = None
    if sweep and len(sweep["rd"]) >= 6:
        out["fp_shape_fits"] = shape_fits(sweep["rd"], sweep["f"], sweep["Q"], beta)
        ab = absorber_sweep(recs)
        RD, GG, BET = [], [], []
        for rd, qp, f in zip(sweep["rd"], sweep["Q"], sweep["f"]):
            k = round(rd, 4)
            if k in ab:
                RD.append(rd)
                GG.append(dl(qp, 8) / dl(Q_of(ab[k]), 8))
                BET.append(beta(f))
        if len(RD) >= 4:
            out["matched_ratio_fit"] = matched_ratio_fit(RD, GG, BET)
            out["matched_ratio_fit"]["absorber_labels"] = [ab[round(r, 4)]["task"]["label"] for r in RD]
            # the same fit with the original absorber records (the earliest run at each radius) in place of the repeats
            ab0 = absorber_sweep(recs, policy="first")
            if all(round(r, 4) in ab0 for r in RD):
                G0 = [dl(qp, 8) / dl(Q_of(ab0[round(rd, 4)]), 8) for rd, qp in zip(sweep["rd"], sweep["Q"]) if round(rd, 4) in ab]
                m0 = matched_ratio_fit(RD, G0, BET)
                mr = out["matched_ratio_fit"]
                out["matched_ratio_fit_original_points"] = dict(
                    R=m0["R"], D=m0["D"], rms_lnG=m0["rms_lnG"], absorber_labels=[ab0[round(r, 4)]["task"]["label"] for r in RD],
                    replaced=[ab[round(r, 4)]["task"]["label"] for r in RD if ab[round(r, 4)] is not ab0[round(r, 4)]],
                    R_difference=mr["R"] - m0["R"], D_difference_a=mr["D"] - m0["D"])
        # the same fit and statistic over every record, as manuscript version 8 computed them (for the superseded list)
        try:
            sweep_all = radius_sweep(recs, policy="campaign", r_min=0.045, r_max=0.115, include_excluded=True)
            ab_all = absorber_sweep(recs, excluded=None)
            RDa, GGa, BETa = [], [], []
            for rd, qp, f in zip(sweep_all["rd"], sweep_all["Q"], sweep_all["f"]):
                if round(rd, 4) in ab_all:
                    RDa.append(rd); GGa.append(dl(qp, 8) / dl(Q_of(ab_all[round(rd, 4)]), 8)); BETa.append(beta(f))
            v8: dict[str, Any] = dict(matched_ratio_fit=matched_ratio_fit(RDa, GGa, BETa) if len(RDa) >= 4 else None)
            for name, rds_, qs_ in (("pml_25a", sweep_all["rd"], sweep_all["Q"]), ("absorber", sorted(ab_all), [Q_of(ab_all[k]) for k in sorted(ab_all)])):
                y = np.log(qs_)
                r = y - np.polyval(np.polyfit(rds_, y, 3), rds_)
                v8.setdefault("detrended_oscillation", {})[name] = dict(detrended_ptp=float(r.max() - r.min()), points=len(rds_))
            out["v8_all_points"] = v8
        except ValueError:
            pass
        abx = absorber_sweep(recs, excluded=True)
        out["absorber_sweep"] = dict(rd=sorted(ab), Q=[Q_of(ab[k]) for k in sorted(ab)], margin=[_margin(ab[k]) for k in sorted(ab)],
                                     labels=[ab[k]["task"]["label"] for k in sorted(ab)],
                                     peak_trough=(max(Q_of(v) for v in ab.values()) / min(Q_of(v) for v in ab.values())) if ab else None,
                                     excluded=dict(rd=sorted(abx), Q=[Q_of(abx[k]) for k in sorted(abx)],
                                                   labels=[abx[k]["task"]["label"] for k in sorted(abx)]))
        g49 = guide49_sweep(recs)
        out["guide49_sweep"] = dict(rd=sorted(g49), Q=[Q_of(g49[k]) for k in sorted(g49)], labels=[g49[k]["task"]["label"] for k in sorted(g49)])
        # detrended oscillation amplitude (post-hoc statistic of the improvement plan)
        for name, rds_, qs_ in (("pml_25a", sweep["rd"], sweep["Q"]), ("absorber", sorted(ab), [Q_of(ab[k]) for k in sorted(ab)])):
            if len(rds_) >= 5:
                y = np.log(qs_)
                r = y - np.polyval(np.polyfit(rds_, y, 3), rds_)
                out.setdefault("detrended_oscillation", {})[name] = dict(raw_max_over_min=float(np.exp(y).max() / np.exp(y).min()),
                                                                       detrended_ptp=float(r.max() - r.min()), factor=float(math.exp(r.max() - r.min())))
        # ---- global fit with the two-mirror D-free shape fit as G(f)
        fpp = out["fp_shape_fits"]["fits"]["two_mirror_D_free"]
        if fpp:
            R_FP, PHI_R, D_FP = fpp["R"], fpp["phi_r"], fpp["D"]
            beta_c = pwe_mod.beta_of_f(pwe)

            def G_of_f(f):
                phi = 2 * float(beta_c(f)) * D_FP + PHI_R
                return (1 + R_FP - 2 * math.sqrt(R_FP) * math.cos(phi)) / (1 - R_FP)
            pts = []
            for key, N, c in (("Q3_pml_cl12", 3, 9), ("Q4_pml_cl12", 4, 8), ("Q5_pml_cl12", 5, 7), ("Q5_pml_cl13", 5, 8), ("Q6_pml_cl14", 6, 8),
                              ("A3_abs_cl12", 3, 9), ("A4_abs_cl12", 4, 8), ("A5_abs_cl12", 5, 7), ("A5_abs_cl13", 5, 8)):
                v = inp[key]
                if v["record"] is not None:
                    pts.append(dict(label=key, N=N, clearance=c, Q=v["Q"], f=v["f"], termination=v["termination"]))
            out["global_fit"] = global_fit(pts, G_of_f, c0["kappa_three_point"], c0["A_three_point"], c0["inv_Q_top_8"], c0["kappa_top"])
            out["global_fit"]["G_fp_used"] = dict(R=R_FP, phi_r=PHI_R, D=D_FP, source="two_mirror_D_free shape fit")
            out["global_fit"]["records"] = [dict(pt, G_fp=(G_of_f(pt["f"]) if pt["termination"] == "pml" else 1.0)) for pt in pts]
    # ---- second radius series (both terminations), leak law transferred from the reference radius
    out["second_radius"] = second_radius_series(recs, pwe, dl, second_radius, c0["kappa"], inv_top8=c0["inv_Q_top_8"], kappa_top=c0["kappa_top"])
    # ---- clearance rule
    req = lambda qw: 8 + math.log(100 * qw / c0["Q_top_8"]) / c0["kappa_top"]      # noqa: E731
    out["clearance_rule"] = dict(rule="c >= 8 + ln(100 Q_w / Q_top(8)) / kappa_top",
                                 required={f"N{n}": req(c0["Q_w"][n]) for n in (4, 5, 6)})
    target_nm = float(inp["Q4_pml_cl12"]["record"]["structure"]["lattice"].get("target_wavelength_nm") or 1550.0)
    # ---- performance in water
    if field_eta and S_measured:
        out["performance"] = performance(c0, pwe, field_eta, S_measured, target_nm=target_nm,
                                         Q_reference_campaign=(out["reference_campaign"] or {}).get("Q"))
        out["clearance_rule"]["required"]["Q_w_opt"] = req(out["performance"]["Q_w_opt"])
        # the supplementary states the rule at the optimum rounded as printed (1.2 x 10^4)
        q2 = float("%.2g" % out["performance"]["Q_w_opt"])
        out["clearance_rule"]["required_at_rounded_optimum"] = dict(Q_w=q2, c=req(q2))
    out["clearance_rule"]["rows"] = {k: int(math.ceil(v - 1e-12)) for k, v in out["clearance_rule"]["required"].items()}
    # ---- tolerance: dQ/dr_d on both sweeps and the cross-term
    out["tolerance"] = tolerance(inp, out.get("absorber_sweep"), sweep, c0, dl, target_nm=target_nm)
    # ---- resolution and PML checks
    if inp["Q5_pml_cl12_res32"]["record"] is not None:
        out["resolution_check"] = dict(Q_res32=inp["Q5_pml_cl12_res32"]["Q"], Q_res24=Q5c7, change_pct=100 * (inp["Q5_pml_cl12_res32"]["Q"] / Q5c7 - 1))
        # SPRAT 1.2.1: the resonance of the same pair; the wavelength shift is scaled to the target
        # wavelength with the rule of section 2.2 (a -> a lambda_t / lambda_r of the reference point)
        f24, f32 = inp["Q5_pml_cl12"]["f"], inp["Q5_pml_cl12_res32"]["f"]
        a_nm = float(inp["Q4_pml_cl12"]["record"]["structure"]["lattice"]["a_nm"])
        lam_ref = a_nm / inp["Q4_pml_cl12"]["f"]
        dlam = a_nm / f32 - a_nm / f24
        out["resolution_check"].update(f_res24=f24, f_res32=f32, f_change_pct=100.0 * (f32 / f24 - 1.0), dlambda_nm=dlam,
                                       dlambda_nm_at_target=dlam * target_nm / lam_ref, target_nm=target_nm)
        ref_rec = inp["Q4_pml_cl12"]["record"]
        out["resolution_check"]["defect_radius_pixels"] = (float(ref_rec["structure"]["defect"]["radius"])
                                                           * float(ref_rec["params"]["numerics"]["resolution"]))
    pmlc = {}
    for key, ref in (("Q4_pml_cl12_pml2", "Q4_pml_cl12"), ("Q5_pml_cl12_pml2", "Q5_pml_cl12")):
        if inp[key]["record"] is not None:
            pmlc[key] = dict(Q_pml2=inp[key]["Q"], Q_pml1=inp[ref]["Q"], change_pct=100 * (inp[key]["Q"] / inp[ref]["Q"] - 1))
    if pmlc:
        out["pml_thickness_check"] = pmlc
    if inp["Q4_pml_nx49"]["record"] is not None:
        out["G_49a"] = dl(inp["Q4_pml_nx49"]["Q"], 8) / dl(A4, 8)
        out["spread_25a_over_49a"] = dl(Q4, 8) / dl(inp["Q4_pml_nx49"]["Q"], 8)
    out["barrier_row_knob"] = barrier_row_knob(recs, Q4, pwe)
    out["barrier_row_series"] = barrier_row_series(out["barrier_row_knob"], inp["Q4_pml_cl12"], c0["kappa"])
    out["ncl8_leak_correction"] = ncl8_leak_correction(recs, c0)
    n8a = ncl8_leak_correction(recs, c0, include_excluded=True)
    if n8a:
        out.setdefault("v8_all_points", {})["ncl8_leak_correction"] = n8a
    return out


def reference_campaign_record(recs: list[dict], radius: float = 0.060) -> dict | None:
    """The campaign record of the reference geometry in the production layer: the quality factor the earlier
    manuscript quoted (the 'v3 headline' of the registry), before the verification runs."""
    from .sensitivity import is_verification
    cands = [r for r in select(recs, mode="harminv", n=1.33, dx=0.0, dy=0.0, resolution=24, guide_periods=25, pml=1.0, no_rows=True,
                               radius=radius, row=4, cladding_rows=12, termination_x="pml")
             if not is_verification(r) and r["params"]["run"].get("tag") == "production"]
    if not cands:
        return None
    r = min(cands, key=lambda r: (r["task"].get("created") or "", r["result"].get("t_used") or 0.0))
    return dict(label=r["task"]["label"], Q=Q_of(r), f=f_of(r), margin=_margin(r))


def ncl8_leak_correction(recs, c0, include_excluded: bool = False) -> dict[str, Any] | None:
    """Leak correction of the n_cl = 8 table (table 1): ln(Q4/Q3) radius by radius, raw and after subtracting
    1/Q_top(c) at the clearances 4 (N = 4) and 5 (N = 3) with the measured leak law."""
    from .sensitivity import base_series
    bs = base_series(recs, cladding_rows=8, resolution=20, rows=(3, 4), include_excluded=include_excluded)
    s3, s4 = bs["series"].get(3), bs["series"].get(4)
    if not s3 or not s4:
        return None
    Qtop = lambda c: c0["Q_top_8"] * math.exp(c0["kappa_top"] * (c - 8))      # noqa: E731
    common = sorted(set(s3["rd"]) & set(s4["rd"]))
    raw, cor = [], []
    for rd in common:
        q3 = s3["Q"][s3["rd"].index(rd)]
        q4 = s4["Q"][s4["rd"].index(rd)]
        raw.append(math.log(q4 / q3))
        cor.append(math.log((1 / (1 / q4 - 1 / Qtop(4))) / (1 / (1 / q3 - 1 / Qtop(5)))))
    n = len(common)
    # the same correction from the quality factors rounded to whole numbers, as table 1 of the paper prints them
    cor_r = []
    for rd in common:
        q3 = float(round(s3["Q"][s3["rd"].index(rd)]))
        q4 = float(round(s4["Q"][s4["rd"].index(rd)]))
        cor_r.append(math.log((1 / (1 / q4 - 1 / Qtop(4))) / (1 / (1 / q3 - 1 / Qtop(5)))))
    return dict(rd=[float(x) for x in common], raw=raw, corrected=cor,
                raw_mean=float(np.mean(raw)), raw_sem=float(np.std(raw, ddof=1) / math.sqrt(n)) if n > 1 else None,
                corrected_mean=float(np.mean(cor)), corrected_sem=float(np.std(cor, ddof=1) / math.sqrt(n)) if n > 1 else None,
                Q_top_clearance4=Qtop(4), Q_top_clearance5=Qtop(5),
                corrected_from_rounded_Q=cor_r, corrected_mean_from_rounded_Q=float(np.mean(cor_r)),
                max_route_difference=float(max(abs(a - b) for a, b in zip(cor, cor_r))),
                corrected_mean_below_reference_pct=float(100.0 * (1.0 - np.mean(cor) / c0["kappa"])),
                scatter_ratio_raw_over_corrected=(float(np.std(raw, ddof=1) / np.std(cor, ddof=1)) if n > 1 else None))


def second_radius_series(recs, pwe, dl, radius: float, kappa_ref: float, inv_top8: float | None = None,
                         kappa_top: float | None = None, draws: int = MC_DRAWS, seed: int = MC_SEED,
                         sigma_leak: float = 0.5) -> dict[str, Any] | None:
    """The kappa series at a second defect radius under both terminations, with the leak law of the reference radius.

    Records: ``series_records`` (one job per series where the records allow it).  Uncertainty of the absorber
    kappa: the quadrature of (i) a Monte Carlo over the independent discretisation error per Q (0.23 %) and a
    +-50 % (one sigma) error of the transferred leak, and (ii) half the spread of the two per-row steps."""
    rows = ((3, 12), (4, 12), (5, 13))
    series = {}
    for term in ("absorber", "pml"):
        recs_t = series_records(recs, [dict(mode="harminv", n=1.33, dx=0.0, dy=0.0, resolution=24, guide_periods=25, pml=1.0, no_rows=True,
                                            radius=radius, row=N, cladding_rows=nk, termination_x=term) for N, nk in rows])
        if any(r is None for r in recs_t):
            return None
        pts = [dict(N=N, clearance=nk - N, Q=Q_of(r), f=f_of(r), margin=_margin(r), label=r["task"]["label"], Q_leak_removed=dl(Q_of(r), nk - N))
               for (N, nk), r in zip(rows, recs_t)]
        k = float(np.polyfit([p["N"] for p in pts], [math.log(p["Q_leak_removed"]) for p in pts], 1)[0])
        series[term] = dict(points=pts, kappa=k, steps=[math.log(pts[i + 1]["Q_leak_removed"] / pts[i]["Q_leak_removed"]) for i in range(2)])
    fbar = float(np.mean([p["f"] for p in series["absorber"]["points"]]))
    cs, ce = pwe_mod.channels_at(pwe, fbar)
    beta_bar = float(pwe_mod.beta_of_f(pwe)(fbar)) / math.pi          # SPRAT 1.2.1: the guided-mode wavevector there
    G = {str(pa["N"]): pp["Q_leak_removed"] / pa["Q_leak_removed"] for pa, pp in zip(series["absorber"]["points"], series["pml"]["points"])}
    df = {str(pa["N"]): pp["f"] - pa["f"] for pa, pp in zip(series["absorber"]["points"], series["pml"]["points"])}
    sens = {}
    for s in (0.0, 0.5, 1.0, 1.5, 2.0):
        pts = series["absorber"]["points"]
        vals = [1 / (1 / p["Q"] - s * (1 / p["Q"] - 1 / p["Q_leak_removed"])) for p in pts]
        sens[f"{s:.1f}"] = float(np.polyfit([p["N"] for p in pts], np.log(vals), 1)[0])
    ka = series["absorber"]["kappa"]
    apts = series["absorber"]["points"]
    NN = np.array([p["N"] for p in apts], float)
    leak = [1 / p["Q"] - 1 / p["Q_leak_removed"] for p in apts]
    rng = random.Random(seed)
    ks = []
    for _ in range(draws):
        s = 1 + rng.gauss(0, sigma_leak)
        inv = [1 / (p["Q"] * (1 + rng.gauss(0, MC_SIGMA_Q))) - s * lk for p, lk in zip(apts, leak)]
        if min(inv) > 0:
            ks.append(float(np.polyfit(NN, -np.log(inv), 1)[0]))
    mc_sd = float(np.std(ks)) if ks else None
    half = 0.5 * abs(series["absorber"]["steps"][0] - series["absorber"]["steps"][1])
    dlnG = {"3->4": math.log(G["4"] / G["3"]), "4->5": math.log(G["5"] / G["4"])}
    return dict(radius=radius, series=series, kappa_absorber=ka, kappa_pml=series["pml"]["kappa"],
                kappa_difference=abs(ka - series["pml"]["kappa"]), f_mean=fbar,
                channels_at_measured_f=dict(slowest=cs, zone_edge=ce, beta_over_pi=beta_bar), channels_at_f=dict(f=fbar),
                dev_zone_edge_pct=100 * (ka / ce - 1), dev_slowest_pct=100 * (ka / cs - 1),
                G=G, dlnG_per_row=dlnG, dlnG_mean=float(np.mean(list(dlnG.values()))),
                f_pull_pml_minus_absorber=df, leak_transfer_sensitivity=sens,
                leak_transfer_range_pct=100.0 * (max(sens.values()) - min(sens.values())) / ka,
                leak_share_max_pct=100.0 * max(lk * p["Q"] for lk, p in zip(leak, apts)),
                kappa_absorber_uncertainty=dict(monte_carlo_sd=mc_sd, half_step_spread=half,
                                                total=(math.sqrt(mc_sd ** 2 + half ** 2) if mc_sd is not None else None),
                                                draws=len(ks), sigma_Q=MC_SIGMA_Q, sigma_leak_scale=sigma_leak, seed=seed),
                kappa_lower_than_reference=ka < kappa_ref)


def vertical_channel(Q_w: float, Q_abs: float, reductions=(0.03, 0.10, 0.50), Q_perp_values=(1e4, 1e5)) -> dict[str, Any]:
    """A third loss channel (radiation out of the plane) in 1/Q_tot = 1/Q_w + 1/Q_perp + 1/Q_abs (section 4).

    With Q_0 = (1/Q_w + 1/Q_abs)^-1, the Q_perp that lowers Q_tot by a fraction r is Q_0 (1 - r) / r, and a given
    Q_perp lowers it by Q_0 / (Q_0 + Q_perp).  Q_perp overtakes the loss to the medium below Q_abs and becomes the
    largest of the three below Q_w."""
    q0 = 1.0 / (1.0 / Q_w + 1.0 / Q_abs)
    key = lambda q: "1e%d" % round(math.log10(q)) if abs(math.log10(q) - round(math.log10(q))) < 1e-12 else "%g" % q    # noqa: E731
    return dict(Q_tot_no_vertical=q0,
                Q_perp_for_reduction={"%g" % (100 * r): q0 * (1.0 - r) / r for r in reductions},
                reduction_pct_at_Q_perp={key(q): 100.0 * q0 / (q0 + q) for q in Q_perp_values},
                Q_tot_at_Q_perp={key(q): 1.0 / (1.0 / q0 + 1.0 / q) for q in Q_perp_values},
                overtakes_absorption_below=Q_abs, largest_channel_below=Q_w,
                relation="1/Q_tot = 1/Q_w + 1/Q_perp + 1/Q_abs")


def performance(c0: dict, pwe: dict, field_eta: dict, S_measured: float, target_nm: float = 1550.0,
                Q_reference_campaign: float | None = None) -> dict[str, Any]:
    """Figures of merit and the loss budget in water (section 4).

    For weak absorption the loss rate is the imaginary part of the same first-order shift that gives the sensitivity:
    with eps_a = (n_a + i k)^2, Q_abs = n_a / (2 k eta_a), and eta_a = n_a S / lambda_r by equation (4), so
    Q_abs = lambda_r / (2 k S) follows from the measured sensitivity alone (scale invariant: lambda_r / S is).  The
    energy fraction used is therefore the one the measured sensitivity requires; Q_abs for the field-map estimators is
    kept in ``absorption`` for comparison."""
    KW = pwe["loss"]["k_water"]
    lam = field_eta["wavelength_nm"]
    n_a = float(field_eta.get("n", 1.33))
    eta_S = S_measured * n_a / lam
    q_abs = lambda eta: n_a / (2 * KW * eta)        # noqa: E731
    Qabs = q_abs(eta_S)
    Mx = lambda x: x * (1 + 2 * x) / (1 + x) ** 3     # noqa: E731
    xo = (1 + math.sqrt(3)) / 2
    QW = c0["Q_w"]
    QT4 = 1 / (1 / QW[4] + 1 / Qabs)
    QT5 = 1 / (1 / QW[5] + 1 / Qabs)
    fom = lambda q: S_measured * q / lam            # noqa: E731
    # the earlier headline: the campaign quality factor of the reference record (production layer)
    q_ref = Q_reference_campaign
    out = dict(S_measured=S_measured, wavelength_nm=lam, k_water=KW, eta_used=eta_S, eta_basis="required by the measured sensitivity",
               Q_abs=Qabs, Q_abs_from_S=lam / (2 * KW * S_measured),
               FOM_lossless_N4=fom(QW[4]), FOM_water_N4=fom(QT4), FOM_water_N5=fom(QT5), FOM_ratio_N5_over_N4=fom(QT5) / fom(QT4),
               Q_tot={"N4": QT4, "N5": QT5}, T_min={"N4": (QT4 / Qabs) ** 2, "N5": (QT5 / Qabs) ** 2},
               metric_fraction={"N4": Mx(QW[4] / Qabs) / Mx(xo), "N5": Mx(QW[5] / Qabs) / Mx(xo)}, Q_w_opt=xo * Qabs, x_opt=xo,
               water_reduction_N4_pct=100.0 * (1.0 - QT4 / QW[4]),
               detection_limit_factor_vs_6927=(q_ref / QW[4]) if q_ref else None, Q_reference_campaign=q_ref)
    # absorption of water: the coefficient behind k, and Q_abs for the other energy fractions of section 3.7
    absn = dict(target_nm=target_nm, alpha_per_cm=4 * math.pi * KW / (target_nm * 1e-7), k_water=KW, n_a=n_a,
                eta_required=eta_S, Q_abs_eta_required=Qabs, Q_reduction_by_water_pct=100.0 * (1.0 - QT4 / QW[4]))
    for key, name in (("eta_hf", "Q_abs_eta_hf"), ("eta_pure", "Q_abs_eta_low"), ("eta_mask", "Q_abs_eta_mask"),
                      # SPRAT 1.2.1: the mixed points of the grid counted wholly as silicon or wholly as
                      # analyte, the two ways a conductivity that is not smoothed at interfaces can weight them
                      ("eta_hf_low", "Q_abs_mixed_as_silicon"), ("eta_hf_high", "Q_abs_mixed_as_analyte")):
        if field_eta.get(key):
            absn[name] = q_abs(field_eta[key])
    if field_eta.get("eta_pure"):
        # the v8 value: the pure-pixel estimate over the stored-map denominator, and its consequences
        Qo = q_abs(field_eta["eta_pure"])
        QT4o, QT5o = 1 / (1 / QW[4] + 1 / Qo), 1 / (1 / QW[5] + 1 / Qo)
        absn["v8_basis"] = dict(eta=field_eta["eta_pure"], Q_abs=Qo, Q_tot={"N4": QT4o, "N5": QT5o}, FOM_water_N4=fom(QT4o),
                                FOM_water_N5=fom(QT5o), T_min={"N4": (QT4o / Qo) ** 2, "N5": (QT5o / Qo) ** 2},
                                metric_fraction={"N4": Mx(QW[4] / Qo) / Mx(xo), "N5": Mx(QW[5] / Qo) / Mx(xo)}, Q_w_opt=xo * Qo)
    out["absorption"] = absn
    out["vertical"] = vertical_channel(QW[4], Qabs)
    return out


def tolerance(inp: dict, absorber_sweep: dict | None, sweep: dict | None, c0: dict, dl, target_nm: float = 1550.0) -> dict[str, Any]:
    """Fabrication tolerance at the reference radius (section 3.8): dQ/dr_d on both sweeps, the cross-term of restoring
    lambda_r after one added row, d lnG/dr_d, and the number of linewidths a 1 nm radius error moves the resonance
    by at the target wavelength (d lambda_r/d r_d is dimensionless, so the shift in nm is scale invariant)."""
    out: dict[str, Any] = {}
    a_nm = 481.4
    if inp["Q4_pml_cl12"]["record"] is not None:
        a_nm = float(inp["Q4_pml_cl12"]["record"]["structure"]["lattice"]["a_nm"])
    sl = lambda xs, ys, i: (ys[i + 1] - ys[i - 1]) / (xs[i + 1] - xs[i - 1]) / a_nm       # noqa: E731
    # d lambda_r / d r_d at the reference radius: the five-point stencil of the fine sweep
    dlam_drd = (sweep or {}).get("local_dlambda_drd", {}).get("nm_per_nm")
    if inp["Q5_pml_cl13"]["record"] is not None and inp["Q4_pml_cl12"]["record"] is not None:
        dlam = a_nm / inp["Q5_pml_cl13"]["f"] - a_nm / inp["Q4_pml_cl12"]["f"]
        out["dlambda_per_row_nm"] = dlam
        drd = dlam / dlam_drd if dlam_drd else None
        out["drd_to_restore_lambda_nm"] = drd
        for name, sw, qref in (("absorber", absorber_sweep, c0["Q_w"][4]), ("pml_25a", sweep, dl(inp["Q4_pml_cl12"]["Q"], 8))):
            if sw and 0.06 in [round(x, 4) for x in sw["rd"]]:
                xs, ys = [float(x) for x in sw["rd"]], [float(q) for q in sw["Q"]]
                i = [round(x, 4) for x in xs].index(0.06)
                if 0 < i < len(xs) - 1:
                    s = sl(xs, ys, i)
                    qhere = ys[i] if name == "pml_25a" else qref
                    out[name] = dict(dQ_drd_per_nm=s, fabrication_tolerance_pct_per_nm=100 * abs(s) / qhere,
                                     dlnQ_drd_per_a=(math.log(ys[i + 1]) - math.log(ys[i - 1])) / (xs[i + 1] - xs[i - 1]))
                    if drd is not None:
                        out[name].update(cross_term_pct_of_row_step=abs(s * drd / qref) / c0["kappa"] * 100,
                                         cross_term_max_pct=max(abs(sl(xs, ys, j) * drd / ys[j]) / c0["kappa"] * 100 for j in range(1, len(xs) - 1)))
                    if dlam_drd:
                        out[name]["Q_here"] = qhere
        if "absorber" in out and "pml_25a" in out and absorber_sweep and sweep:
            # d lnG / d r_d over the same three radii: G = Q_PML(leak removed) / Q_abs(leak removed)
            xa = [round(float(x), 4) for x in absorber_sweep["rd"]]
            xp = [round(float(x), 4) for x in sweep["rd"]]
            if all(k in xa and k in xp for k in (0.055, 0.065)):
                lnG = {k: math.log(dl(float(sweep["Q"][xp.index(k)]), 8) / dl(float(absorber_sweep["Q"][xa.index(k)]), 8)) for k in (0.055, 0.065)}
                out["dlnG_drd_per_a"] = (lnG[0.065] - lnG[0.055]) / 0.01
    if dlam_drd:
        lw = dict(dlambda_drd_nm_per_nm=dlam_drd, target_nm=target_nm)
        for name in ("pml_25a", "absorber"):
            if name in out:
                lw[name] = dlam_drd * 1.0 * out[name]["Q_here"] / target_nm
        out["linewidths_per_nm_radius_error"] = lw
    return out


def barrier_row_knob(recs, Q_ref: float, pwe: dict) -> dict[str, Any]:
    out = {}
    pred = {x["r_b"]: x for x in pwe.get("barrier_radius", [])}
    for r in select(recs, mode="harminv", n=1.33, dx=0.0, dy=0.0, resolution=24, guide_periods=25, pml=1.0, termination_x="pml",
                    radius=0.060, row=4, cladding_rows=12):
        rows = {int(k): float(v) for k, v in (r["structure"].get("rows") or {}).items()}
        if rows == {2: rows.get(2)} and 2 in rows:
            rb = rows[2]
            a_nm = float(r["structure"]["lattice"]["a_nm"])
            m = strongest_mode(r)
            out[f"{rb:g}"] = dict(Q=Q_of(r), ratio=Q_of(r) / Q_ref, label=r["task"]["label"],
                                  predicted_ratio=(pred[rb]["Q_ratio_one_row"] if rb in pred else None),
                                  f=f_of(r), lambda_nm=a_nm / f_of(r), lambda_nm_stored=(m.get("wavelength_nm") if m else None),
                                  margin=_margin(r), host=r["provenance"].get("host", ""), job_id=r["provenance"].get("job_id", ""),
                                  code_sha256=r["software"].get("code_sha256", ""), a_nm=a_nm)
    return out


def barrier_row_series(knob: dict, ref: dict, kappa: float, nominal: float = 0.20) -> dict[str, Any] | None:
    """The radius of the second barrier row (section 3.4): Q, f_r and lambda_r = a/f_r over the series with the nominal
    radius taken from the reference record, the resonance shift from the smallest to the largest radius, the change
    of Q and its share of a row step in ln Q."""
    if not knob or ref.get("record") is None:
        return None
    a_nm = float(ref["record"]["structure"]["lattice"]["a_nm"])
    pts = {float(k): dict(Q=v["Q"], f=v["f"], label=v["label"]) for k, v in knob.items()}
    pts.setdefault(nominal, dict(Q=ref["Q"], f=ref["f"], label=ref["label"]))
    rb = sorted(pts)
    lo, hi = pts[rb[0]], pts[rb[-1]]
    ln_ratio = math.log(hi["Q"] / lo["Q"])
    pred = [v for v in knob.values() if v.get("predicted_ratio")]
    closer = sum(1 for v in pred if abs(math.log(v["ratio"])) < abs(math.log(v["predicted_ratio"])))
    return dict(rb=rb, Q=[pts[x]["Q"] for x in rb], f=[pts[x]["f"] for x in rb], lambda_nm=[a_nm / pts[x]["f"] for x in rb],
                labels=[pts[x]["label"] for x in rb], a_nm=a_nm,
                lambda_shift_nm=a_nm / hi["f"] - a_nm / lo["f"], f_change=hi["f"] - lo["f"],
                Q_change_pct=100.0 * (hi["Q"] / lo["Q"] - 1.0), ln_Q_ratio=ln_ratio, row_step_fraction_pct=100.0 * ln_ratio / kappa,
                ratios_with_prediction=len(pred), ratios_closer_to_unity_than_predicted=closer)


__all__ = ["inputs", "chain", "monte_carlo", "airy", "matched_ratio_fit", "shape_fits", "global_fit", "analyse",
           "second_radius_series", "performance", "vertical_channel", "tolerance", "barrier_row_knob", "barrier_row_series",
           "absorber_sweep", "series_records", "reference_campaign_record", "guide49_sweep"]
