"""Plane-wave expansion (TM, square lattice of dielectric rods): the engine-independent checks.

Units a = c = 1, f = a / lambda.  The module ports the plane-wave layer of the 2026 campaign:

1. the bulk TM gap (validation against MPB);
2. the W1 guided band beta(f) and group index n_g(f) from a 1 x 13 supercell;
3. the complex band structure at k_x = beta(f): the decay constant kappa = 2 Im k_y a of the
   least-evanescent Bloch channel, which in this gap is the zone-edge channel (Re k_y = pi/a).
   Only physical roots are kept: the companion linearisation of a truncated basis also returns
   roots whose eigenvectors sit on the outermost G_y rings of the basis (truncation artefacts);
   they are rejected by their edge weight and by requiring the replica with |Re k_y| <= pi/a
   (``decay_roots``; the rejected roots are reported by ``boundary_roots`` and in the
   ``truncation_check`` block, square against circular basis);
3b. a reciprocity check: decay of the W1 guided-mode field at the defect sites;
4. the Fabry-Perot reading of the Q(r_d) oscillation, D = pi / (beta(f_1) - beta(f_2));
5. the barrier-row radius knob: kappa per row against the rod radius of one barrier row;
6. the water loss budget and the notch coupling optimum M(x) = x(1+2x)/(1+x)^3.

Only numpy and scipy are needed.  ``run`` is the task entry point (mode = pwe).
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.linalg import eig, eigh
from scipy.special import j1

# ------------------------------------------------------------------ the solver
def gvecs(M: int, Ny: int = 1) -> tuple[np.ndarray, np.ndarray]:
    m, n = np.meshgrid(np.arange(-M, M + 1), np.arange(-M * Ny, M * Ny + 1), indexing="ij")
    return 2 * np.pi * m.ravel(), 2 * np.pi * n.ravel() / Ny


def eps_matrix(Gx, Gy, eps_rod: float, eps_bg: float, r: float, Ny: int = 1, rods_y=(0,)) -> np.ndarray:
    """Fourier matrix eps(G - G') for rods of radius r at (0, y_j) in a 1 x Ny cell."""
    dGx, dGy = Gx[:, None] - Gx[None, :], Gy[:, None] - Gy[None, :]
    g = np.hypot(dGx, dGy)
    x = g * r
    safe = np.where(x > 1e-12, x, 1.0)
    form = np.where(x > 1e-12, 2 * j1(safe) / safe, 1.0)
    S = sum(np.exp(-1j * dGy * y) for y in rods_y)
    E = (eps_rod - eps_bg) * (np.pi * r ** 2 / Ny) * form * S
    E[g < 1e-12] += eps_bg
    return E


def bands(kx: float, ky: float, Gx, Gy, E, nb: int) -> np.ndarray:
    A = np.diag((kx + Gx) ** 2 + (ky + Gy) ** 2).astype(complex)
    return np.sqrt(np.abs(eigh(A, E, eigvals_only=True, subset_by_index=[0, nb - 1]))) / (2 * np.pi)


EDGE_TOL = 1e-3            # largest eigenvector weight on the two outermost G_y rings of a physical root
RE_TOL = 5e-3              # slack on |Re k_y| <= pi/a for the central replica (rad/a)


def _companion(f: float, kx: float, Gx, Gy, E):
    N = len(Gx)
    w2 = (2 * np.pi * f) ** 2
    D0 = np.diag((kx + Gx) ** 2 + Gy ** 2).astype(complex) - w2 * E
    return np.block([[np.zeros((N, N)), np.eye(N)], [-D0, -np.diag(2 * Gy).astype(complex)]])


def _roots_with_edge_weight(f: float, kx: float, Gx, Gy, E) -> tuple[np.ndarray, np.ndarray]:
    """All roots decaying in +y and the weight of each eigenvector on the two outermost G_y rings."""
    N = len(Gx)
    k, V = eig(_companion(f, kx, Gx, Gy, E))
    my = np.abs(Gy) / (2 * np.pi)
    edge = my >= my.max() - 1 - 1e-9
    P = np.abs(V[:N, :]) ** 2
    P = P / P.sum(axis=0, keepdims=True)
    w = P[edge, :].sum(axis=0)
    keep = k.imag > 1e-7
    return k[keep], w[keep]


def decay_roots(f: float, kx: float, Gx, Gy, E, physical: bool = True, edge_tol: float = EDGE_TOL) -> np.ndarray:
    """Complex k_y at fixed f and k_x decaying in +y (companion linearisation of the TM problem, bulk crystal).

    ``physical=True`` keeps the Bloch channels only: roots whose eigenvector puts more than ``edge_tol`` of its
    weight on the two outermost G_y rings of the truncated basis are truncation artefacts, and of each physical
    root only the replica with |Re k_y| <= pi/a is kept.  ``physical=False`` returns every root, as the campaign
    layer did (its slowest root was an artefact of the square truncation)."""
    if not physical:
        k = eig(_companion(f, kx, Gx, Gy, E), right=False)
        return k[k.imag > 1e-7]
    k, w = _roots_with_edge_weight(f, kx, Gx, Gy, E)
    return k[(w < edge_tol) & (np.abs(k.real) <= np.pi + RE_TOL)]


def boundary_roots(f: float, kx: float, Gx, Gy, E, min_weight: float = 0.3) -> np.ndarray:
    """The roots rejected as truncation artefacts: eigenvectors localised on the outermost G_y rings."""
    k, w = _roots_with_edge_weight(f, kx, Gx, Gy, E)
    return k[w > min_weight]


def slowest_decay(f: float, kx: float, Gx, Gy, E) -> float:
    """Im k_y of the least-evanescent physical channel."""
    return float(decay_roots(f, kx, Gx, Gy, E).imag.min())


def zone_edge_decay(f: float, kx: float, Gx, Gy, E, tol: float = 2e-3) -> float:
    """Im k_y of the physical channel at the zone edge (Re k_y = pi/a)."""
    k = decay_roots(f, kx, Gx, Gy, E)
    re = np.abs((k.real + np.pi) % (2 * np.pi) - np.pi)
    return float(k.imag[np.abs(re - np.pi) < tol].min())


def gvecs_circle(M: int) -> tuple[np.ndarray, np.ndarray]:
    """A circular truncation |m| <= M of the reciprocal lattice (for the basis-shape check)."""
    m, n = np.meshgrid(np.arange(-M, M + 1), np.arange(-M, M + 1), indexing="ij")
    sel = (m ** 2 + n ** 2) <= M * M
    return 2 * np.pi * m[sel].ravel().astype(float), 2 * np.pi * n[sel].ravel().astype(float)


def truncation_check(f: float, kx: float, eps_rod: float, eps_bg: float, r: float) -> dict[str, Any]:
    """The physical channel and the slowest truncation artefact for square bases M = 7, 10, 13 and a circular
    basis |m| <= 16: the physical value does not depend on the basis, the artefact does."""
    out: dict[str, Any] = dict(f=f, kx_over_pi=kx / np.pi, edge_tol=EDGE_TOL, bases={})
    for name, (Gx, Gy) in (("square M=7", gvecs(7)), ("square M=10", gvecs(10)), ("square M=13", gvecs(13)),
                           ("circle |m|<=16", gvecs_circle(16))):
        E = eps_matrix(Gx, Gy, eps_rod, eps_bg, r)
        kk, w = _roots_with_edge_weight(f, kx, Gx, Gy, E)
        phys = kk[(w < EDGE_TOL) & (np.abs(kk.real) <= np.pi + RE_TOL)]
        art = kk[w > 0.3]
        i = int(np.argmin(phys.imag))
        out["bases"][name] = dict(plane_waves=int(len(Gx)), kappa_physical=float(2 * phys.imag[i]),
                                  re_ky_over_pi=float(abs(phys.real[i]) / np.pi),
                                  kappa_slowest_artefact=(float(2 * art.imag.min()) if len(art) else None))
    return out


# ------------------------------------------------------------------ the feasibility set
def feasibility(eps_rod: float = 11.90, n_analyte: float = 1.33, r: float = 0.20, reference_f: float = 0.30461,
                second_f: float = 0.29230, kappa_reference: float = 1.9595, kappa_second: float = 1.7386,
                quick: bool = False, log=print) -> dict[str, Any]:
    eps_bg = n_analyte ** 2
    out: dict[str, Any] = {"model": dict(eps_rod=eps_rod, n_analyte=n_analyte, r_over_a=r, polarisation="TM (E_z)")}
    t0 = time.time()
    # 1. bulk gap
    M = 7 if quick else 10
    Gx, Gy = gvecs(M)
    Eb = eps_matrix(Gx, Gy, eps_rod, eps_bg, r)
    path = [(np.pi * t, 0.0) for t in np.linspace(0, 1, 21)] + [(np.pi, np.pi * t) for t in np.linspace(0, 1, 21)[1:]] + \
           [(np.pi * t, np.pi * t) for t in np.linspace(1, 0, 29)[1:]]
    b = np.array([bands(kx, ky, Gx, Gy, Eb, 3) for kx, ky in path])
    out["bulk_gap"] = dict(M=M, f_lo=float(b[:, 0].max()), f_hi=float(b[:, 1].min()))
    GAP = (out["bulk_gap"]["f_lo"], out["bulk_gap"]["f_hi"])
    log(f"pwe 1 gap {out['bulk_gap']} ({time.time() - t0:.0f} s)")
    # 2. W1 band (1 x 13 supercell, row y = 0 removed)
    p, Mx = 6, 5
    Gxs, Gys = gvecs(Mx, Ny=2 * p + 1)
    Es = eps_matrix(Gxs, Gys, eps_rod, eps_bg, r, Ny=2 * p + 1, rods_y=[j for j in range(-p, p + 1) if j != 0])
    svals = np.linspace(0.40, 0.66, 27)
    fw = []
    for sv in svals:
        g = [x for x in bands(sv * np.pi, 0.0, Gxs, Gys, Es, 2 * p + 6) if GAP[0] + 1e-4 < x < GAP[1] - 1e-4]
        fw.append(g[0])
    fw = np.array(fw)
    kf, fk = CubicSpline(fw, svals), CubicSpline(svals, fw)
    beta = lambda f: float(kf(f)) * np.pi                       # noqa: E731
    n_g = lambda f: 1.0 / (2.0 * float(fk(float(kf(f)), 1)))    # noqa: E731
    out["W1"] = dict(supercell=f"1x{2 * p + 1}", Mx=Mx, kx_over_pi=svals.tolist(), f=fw.tolist())
    log(f"pwe 2 W1 band done ({time.time() - t0:.0f} s)")
    # 3. kappa at the two resonance frequencies
    pts = {"reference": (reference_f, kappa_reference), "second": (second_f, kappa_second)}
    out["kappa"] = {}
    def channels(f, kx, with_next=False):
        """(slowest physical, zone-edge physical, slowest truncation artefact) as kappa = 2 Im k_y a; one eigensolve.
        With ``with_next`` the next physical channel (the slowest one distinct from the first) and its Re k_y / (pi/a)
        are appended."""
        kk, w = _roots_with_edge_weight(f, kx, Gx, Gy, Eb)
        k = kk[(w < EDGE_TOL) & (np.abs(kk.real) <= np.pi + RE_TOL)]
        art = kk[w > 0.3]
        re = np.abs((k.real + np.pi) % (2 * np.pi) - np.pi)
        res_ = (float(2 * k.imag.min()), float(2 * k.imag[np.abs(re - np.pi) < 2e-3].min()),
                float(2 * art.imag.min()) if len(art) else None)
        if with_next:
            im0 = k.imag.min()
            rest = k[k.imag > im0 * (1 + 1e-6) + 1e-9]
            j = int(np.argmin(rest.imag))
            res_ += (float(2 * rest.imag[j]), float(abs(rest.real[j]) / np.pi))
        return res_

    for name, (f, km) in pts.items():
        kp, kz, ka, kn, rn = channels(f, beta(f), with_next=True)
        k0, _, a0 = channels(f, 0.0)
        kpi, _, api = channels(f, np.pi)
        out["kappa"][name] = dict(f=f, beta_over_pi=beta(f) / np.pi, n_g=n_g(f), kappa_pred=kp,
                                  kappa_pred_zone_edge=kz, kappa_measured=km,
                                  deviation_pct=100 * (kp - km) / km, deviation_zone_edge_pct=100 * (kz - km) / km,
                                  kappa_at_kx0=k0, kappa_at_kxpi=kpi, kappa_next=kn, re_ky_next_over_pi=rn,
                                  kappa_boundary_artefact=ka, kappa_boundary_artefact_kx0=a0, kappa_boundary_artefact_kxpi=api)
        log(f"pwe 3 {name}: kappa physical {kp:.4f} (zone edge {kz:.4f}); next physical channel {kn:.4f} "
            f"(Re k_y = {rn:.3f} pi/a); slowest truncation artefact {ka}")
    fr0 = pts["reference"][0]
    out["truncation_check"] = truncation_check(fr0, beta(fr0), eps_rod, eps_bg, r) if not quick else None
    fgrid = np.linspace(0.2835, 0.3240, 10 if quick else 28)
    curve = []
    for f in fgrid:
        kp, kz, _ = channels(f, beta(f))
        curve.append(dict(f=float(f), beta_over_pi=beta(f) / np.pi, n_g=n_g(f), kappa_pred=kp, kappa_pred_zone_edge=kz))
    out["kappa_curve"] = curve
    # 3b. reciprocity: decay of the guided-mode field at the defect sites, 1 x 21 supercell
    pr = 10
    Nyr = 2 * pr + 1
    Gxr, Gyr = gvecs(5, Ny=Nyr)
    Er = eps_matrix(Gxr, Gyr, eps_rod, eps_bg, r, Ny=Nyr, rods_y=[j for j in range(-pr, pr + 1) if j != 0])
    out["reciprocity"] = {}
    for name, (f, km) in pts.items():
        bb = beta(f)
        A = np.diag((bb + Gxr) ** 2 + Gyr ** 2).astype(complex)
        wv, V = eigh(A, Er, subset_by_index=[0, 2 * pr + 5])
        fsr = np.sqrt(np.abs(wv)) / (2 * np.pi)
        i = int(np.argmin(np.abs(fsr - f)))
        v = V[:, i]
        site = [abs(np.sum(v * np.exp(1j * Gyr * j))) ** 2 for j in range(1, 9)]
        out["reciprocity"][name] = dict(f_mode=float(fsr[i]),
                                        kappa_rows=[float(-np.log(site[j + 1] / site[j])) for j in range(7)])
    # 4. Fabry-Perot reading of the Q(r_d) oscillation (the campaign's quadratic f_r(r_d) law)
    ffit = lambda rd: 0.31182 - 0.01751 * rd - 1.77681 * rd ** 2          # noqa: E731
    f1 = ffit(0.065) + (reference_f - ffit(0.060))
    df = 0.01365
    sig = float(np.hypot(abs(-0.01751 - 2 * 1.77681 * 0.065) * 0.0025, abs(-0.01751 - 2 * 1.77681 * 0.105) * 0.0025))
    D = lambda d: float(np.pi / (beta(f1) - beta(f1 - d)))                 # noqa: E731
    ng_mean = (beta(f1) - beta(f1 - df)) / (2 * np.pi * df)
    out["fabry_perot"] = dict(f_peak1=f1, delta_f=df, sigma_delta_f=sig, D_over_a=D(df), D_range=[D(df + sig), D(df - sig)],
                              ng_mean=ng_mean, W1_end_over_a=12.5, peak_trough_ratio=6949 / 2037,
                              R_two_mirror=((np.sqrt(6949 / 2037) - 1) / (np.sqrt(6949 / 2037) + 1)) ** 2,
                              R_one_mirror=((6949 / 2037 - 1) / (6949 / 2037 + 1)) ** 2,
                              phase_change_rad={"N2->N3": 2 * 12.5 * 2 * np.pi * n_g(reference_f) * 0.0024,
                                                "N3->N4": 2 * 12.5 * 2 * np.pi * n_g(reference_f) * 0.0003},
                              predicted_period_for_W1_49a=1 / (2 * ng_mean * 24.5))
    # 5. barrier-row radius
    fr, br = reference_f, beta(reference_f)
    out["barrier_radius"] = []
    for rb in (0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23):
        kr = 2 * slowest_decay(fr, br, Gx, Gy, eps_matrix(Gx, Gy, eps_rod, eps_bg, rb))
        out["barrier_radius"].append(dict(r_b=rb, kappa_row=kr))
    k20 = [x["kappa_row"] for x in out["barrier_radius"] if x["r_b"] == 0.20][0]
    for x in out["barrier_radius"]:
        x["Q_ratio_one_row"] = float(np.exp(x["kappa_row"] - k20))
    # 6. loss budget and coupling optimum
    k_w = 10.8 * 1550e-7 / (4 * np.pi)          # water extinction at 1550 nm from alpha = 10.8 /cm
    Mx_ = lambda x: x * (1 + 2 * x) / (1 + x) ** 3    # noqa: E731
    xo = (1 + np.sqrt(3)) / 2
    out["loss"] = dict(k_water=k_w, Q_abs={"eta=0.544": n_analyte / (2 * k_w * 0.544), "eta=0.618": n_analyte / (2 * k_w * 0.618)},
                       rows=[dict(Q_abs=q, Qtot_N4=1 / (1 / 6927 + 1 / q), Qtot_N5=1 / (1 / 46029 + 1 / q),
                                  Tmin_N4=(1 / (1 / 6927 + 1 / q) / q) ** 2, Tmin_N5=(1 / (1 / 46029 + 1 / q) / q) ** 2,
                                  metric_frac_N4=Mx_(6927 / q) / Mx_(xo), metric_frac_N5=Mx_(46029 / q) / Mx_(xo), Qw_opt=xo * q)
                             for q in (8.1e3, 8.65e3, 9.2e3)], x_opt=xo)
    kr = out["kappa"]["reference"]
    out["pred_Q"] = dict(Q5_from_Q4=[6927 * np.exp(kr["kappa_pred"]), 6927 * np.exp(kr["kappa_pred_zone_edge"])],
                         Q6_over_Q5=[np.exp(kr["kappa_pred"]), np.exp(kr["kappa_pred_zone_edge"])],
                         Q6_from_46029=[46029 * np.exp(kr["kappa_pred"]), 46029 * np.exp(kr["kappa_pred_zone_edge"])])
    out["runtime_s"] = round(time.time() - t0, 1)
    return _plain(out)


def _plain(o):
    if isinstance(o, dict):
        return {k: _plain(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_plain(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    return o


# ------------------------------------------------------------------ access helpers
def kappa_block(pwe: dict, which: str = "reference") -> dict:
    """The kappa entry of a pwe result, for the new ('reference'/'second') or the legacy keys."""
    k = pwe["kappa"]
    legacy = {"reference": "r_d=0.060a (reference)", "second": "r_d=0.100a"}
    if which in k:
        return k[which]
    return k[legacy[which]]


def beta_of_f(pwe: dict):
    """Interpolator f -> beta (rad/a) from the W1 band of a pwe result."""
    from scipy.interpolate import interp1d
    fw = np.array(pwe["W1"]["f"])
    kx = np.pi * np.array(pwe["W1"]["kx_over_pi"])
    return interp1d(fw, kx, kind="cubic", fill_value="extrapolate")


def channels_at(pwe: dict, f: float) -> tuple[float, float]:
    """(slowest, zone-edge) kappa interpolated to f on the kappa curve.  In a record written by SPRAT 1.1 both are the
    physical zone-edge channel; in the campaign layer (pwe_legacy) the first was the slowest truncation artefact."""
    cur = pwe["kappa_curve"]
    fs = [c["f"] for c in cur]
    return (float(np.interp(f, fs, [c["kappa_pred"] for c in cur])),
            float(np.interp(f, fs, [c["kappa_pred_zone_edge"] for c in cur])))


def run(record: dict, structure: dict, params: dict) -> dict:
    p = params["pwe"]
    la = structure["lattice"]
    quiet = params["run"].get("quiet", True)
    res = feasibility(eps_rod=float(la["rod_eps"]), n_analyte=float(params["analyte"]["n"]), r=float(la["rod_radius"]),
                      reference_f=float(p["reference_f"]), second_f=float(p["second_f"]),
                      kappa_reference=float(p["kappa_reference"]), kappa_second=float(p["kappa_second"]),
                      quick=bool(p["quick"]), log=(lambda *a: None) if quiet else print)
    record["result"] = res
    return record


__all__ = ["gvecs", "gvecs_circle", "eps_matrix", "bands", "decay_roots", "boundary_roots", "slowest_decay", "zone_edge_decay",
           "truncation_check", "EDGE_TOL", "feasibility",
           "kappa_block", "beta_of_f", "channels_at", "run"]
