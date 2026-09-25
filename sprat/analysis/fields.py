"""Field-map analyses: the analyte energy fraction and the row-by-row decay of the stored resonant field resolved in k_x.

The analyte energy fraction.  In TM polarisation E_z is tangential to every interface, so under Meep's anisotropic
subpixel smoothing (eps_averaging=True) the permittivity that E_z sees at its own Yee points is the arithmetic mean
<eps> = f eps_rod + (1 - f) eps_a over the voxel (f is the silicon fill fraction), and the Hellmann-Feynman theorem for
the discretised operator gives

    dlambda_r/dn_a = (lambda_r / n_a) eta_a,   eta_a = sum (1 - f) eps_a |E_z|^2 / sum <eps> |E_z|^2   (sums on the E_z grid).

The permittivity map stored with the field (``sim.get_array(component=mp.Dielectric)``) is not <eps>: Meep writes, at
the centred points of the grid, the harmonic mean of the eigenvalues of the smoothed tensor, three over the trace of
its inverse, with each diagonal component of the inverse averaged from its own Yee points
(12 / [sum of the 4-point sums of (eps^-1)_xx, (eps^-1)_yy, (eps^-1)_zz]).
``meep_smoothing`` reconstructs the smoothed tensor from the geometry (exact disc-voxel areas, the normal of the rod
surface), checks the reconstruction against the stored map point by point and returns <eps> and 1 - f on the E_z grid,
averaged to the centred points of the stored field.  From it:

    eta_hf       the Hellmann-Feynman estimate consistent with the discretisation,
    eta_hf_low   pure-analyte points only (every mixed point counted as silicon), a lower value,
    eta_hf_high  eta_hf_low plus the energy of the mixed points, an upper value.

The estimators of the original analysis read the stored map as if it were <eps> (``stored-map estimators``): the masked
estimator of the campaign, the linear unmixing, the pure-pixel value and the exact circle-pixel geometry, each over the
denominator sum eps_stored |E|^2.  Since eps_stored < <eps> in a mixed point, that denominator is too small and those
values are biased upwards; they are kept for the comparison with the original analysis.  The first-order sensitivity is
S = lambda_r eta / n_a.
"""

from __future__ import annotations

import math
import os
from typing import Any

import numpy as np

from .. import pwe as pwe_mod
from .collect import select


def load_map(record: dict) -> dict[str, Any]:
    path = os.path.join(os.path.dirname(record.get("_file", "")), record["result"]["field_file"])
    z = np.load(path)
    d = {k: z[k] for k in z.files}
    d["path"] = path
    return d


# --------------------------------------------------------------------------- Meep's smoothed permittivity
_GL_X, _GL_W = np.polynomial.legendre.leggauss(24)


def _disc_rect_area(cx: float, cy: float, r: float, x0: float, x1: float, y0: float, y1: float) -> float:
    """Area of the intersection of a disc with an axis-aligned rectangle: the chord length integrated over x with
    Gauss-Legendre on the pieces between the points where the circle crosses the rectangle (exact to rounding)."""
    a, b = max(x0, cx - r), min(x1, cx + r)
    if b <= a:
        return 0.0
    bps = [a, b]
    for yy in (y0, y1):
        dy = yy - cy
        if abs(dy) < r:
            s = math.sqrt(r * r - dy * dy)
            bps += [x for x in (cx - s, cx + s) if a < x < b]
    bps.sort()
    area = 0.0
    for p, q in zip(bps[:-1], bps[1:]):
        if q - p <= 1e-15:
            continue
        x = 0.5 * (p + q) + 0.5 * (q - p) * _GL_X
        half = np.sqrt(np.maximum(r * r - (x - cx) ** 2, 0.0))
        chord = np.clip(cy + half, y0, y1) - np.clip(cy - half, y0, y1)
        area += 0.5 * (q - p) * float(np.sum(_GL_W * np.maximum(chord, 0.0)))
    return area


def _rod_pattern(r: float, fx: float, fy: float, res: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, int, int]:
    """Silicon fill fraction and surface normal on the half-pixel grid around one rod whose centre lies (fx, fy)
    half-pixels from a half-grid point.  Returns (F, NX, NY, u0, v0): the arrays and the half-grid offset of their
    first element relative to that point."""
    h = 1.0 / res
    m = int(math.ceil(2 * res * (r + h))) + 1
    u = np.arange(-m, m + 1)
    F = np.zeros((u.size, u.size))
    NX = np.zeros_like(F)
    NY = np.zeros_like(F)
    for i, du in enumerate(u):
        x = (du - fx) * h / 2                       # location relative to the rod centre
        for j, dv in enumerate(u):
            y = (dv - fy) * h / 2
            if abs(x) > r + h or abs(y) > r + h:
                continue
            fr = _disc_rect_area(0.0, 0.0, r, x - h / 2, x + h / 2, y - h / 2, y + h / 2) / (h * h)
            if fr > 0:
                F[i, j] = fr
                n = math.hypot(x, y)
                NX[i, j], NY[i, j] = ((x / n, y / n) if n > 0 else (1.0, 0.0))
    return F, NX, NY, -m, -m


def meep_smoothing(rec: dict, d: dict) -> dict[str, Any] | None:
    """Reconstruct Meep's smoothed permittivity on the grid of a stored field map (see the module docstring).

    Returns the reconstruction of the stored map at its centred points (``stored``), <eps> and 1 - f averaged from the
    four surrounding E_z points (``eps_zz``, ``analyte``) and the same evaluated at the centred point itself
    (``eps_zz_c``, ``analyte_c``), or None when the map is not on a centred grid about the origin."""
    from .. import geometry as geo
    s, p = rec["structure"], rec["params"]
    eps = np.asarray(d["eps"], dtype=np.float64)
    nx, ny = eps.shape
    res = int(round(float(d.get("resolution", p["numerics"]["resolution"]))))
    if nx % 2 or ny % 2:
        return None
    e_rod = float(s["lattice"]["rod_eps"])
    n_a = float(p["analyte"]["n"])
    e_a = n_a * n_a
    g = geo.build(s, p)
    F = np.zeros((2 * nx + 1, 2 * ny + 1))
    NX = np.zeros_like(F)
    NY = np.zeros_like(F)
    cache: dict = {}
    for cx, cy, r in g.rods:
        uc, vc = 2 * res * cx + nx, 2 * res * cy + ny          # rod centre in half-grid index units
        iu, iv = int(math.floor(uc + 0.5)), int(math.floor(vc + 0.5))
        fx, fy = round(uc - iu, 9), round(vc - iv, 9)
        key = (round(r, 9), fx, fy)
        if key not in cache:
            cache[key] = _rod_pattern(r, fx, fy, res)
        P, PX, PY, u0, v0 = cache[key]
        a0, b0 = iu + u0, iv + v0
        a1, b1 = a0 + P.shape[0], b0 + P.shape[1]
        ca0, cb0, ca1, cb1 = max(a0, 0), max(b0, 0), min(a1, F.shape[0]), min(b1, F.shape[1])
        if ca0 >= ca1 or cb0 >= cb1:
            continue
        sub = (slice(ca0 - a0, ca1 - a0), slice(cb0 - b0, cb1 - b0))
        tgt = (slice(ca0, ca1), slice(cb0, cb1))
        hit = P[sub] > 0
        F[tgt] += P[sub]
        NX[tgt] = np.where(hit, PX[sub], NX[tgt])
        NY[tgt] = np.where(hit, PY[sub], NY[tgt])
    F = np.clip(F, 0.0, 1.0)
    avg = F * e_rod + (1.0 - F) * e_a
    avinv = F / e_rod + (1.0 - F) / e_a
    mixed = (F > 0) & (F < 1)
    ixx = np.where(mixed, NX ** 2 * avinv + (1 - NX ** 2) / avg, 1.0 / avg)
    iyy = np.where(mixed, NY ** 2 * avinv + (1 - NY ** 2) / avg, 1.0 / avg)
    izz = 1.0 / avg
    E0, E1, O0, O1 = slice(0, -2, 2), slice(2, None, 2), slice(1, -1, 2), slice(1, -1, 2)
    # centred point (i, j) = half-grid (2i+1, 2j+1); E_z at the four even neighbours, E_x at (odd, even), E_y at (even, odd)
    zz4 = izz[E0, E0] + izz[E1, E0] + izz[E0, E1] + izz[E1, E1]
    xx4 = 2 * (ixx[O0, E0] + ixx[O0, E1])
    yy4 = 2 * (iyy[E0, O1] + iyy[E1, O1])
    stored = 12.0 / (zz4 + xx4 + yy4)
    eps_zz = 0.25 * (avg[E0, E0] + avg[E1, E0] + avg[E0, E1] + avg[E1, E1])
    analyte = 1.0 - 0.25 * (F[E0, E0] + F[E1, E0] + F[E0, E1] + F[E1, E1])
    inner = (slice(2, -2), slice(2, -2))
    dev = np.abs(stored - eps)[inner]
    return dict(stored=stored, eps_zz=eps_zz, analyte=analyte, eps_zz_c=avg[O0, O1], analyte_c=1.0 - F[O0, O1],
                reconstruction_max_abs=float(dev.max()), reconstruction_rms=float(np.sqrt((dev ** 2).mean())),
                rod_patterns=len(cache), resolution=res)


def hellmann_feynman(rec: dict, d: dict, S_measured: float | None = None) -> dict[str, Any] | None:
    """eta_a and S_th from the reconstructed <eps> (module docstring), with the lower and upper values and the
    energy share of the mixed points, on the E_z grid (4-point average) and at the centred points."""
    sm = meep_smoothing(rec, d)
    if sm is None:
        return None
    n = float(rec["params"]["analyte"]["n"])
    e_a = n * n
    lam = float(d.get("wavelength_nm", rec["result"].get("wavelength_nm")))
    Ez = np.asarray(d["Ez"])[: sm["eps_zz"].shape[0], : sm["eps_zz"].shape[1]]
    E2 = np.abs(Ez.astype(np.complex128)) ** 2
    out: dict[str, Any] = dict(reconstruction_max_abs=sm["reconstruction_max_abs"], reconstruction_rms=sm["reconstruction_rms"],
                               rod_patterns=sm["rod_patterns"], wavelength_nm=lam, n=n)
    for tag, ezz, an in (("", sm["eps_zz"], sm["analyte"]), ("_centred", sm["eps_zz_c"], sm["analyte_c"])):
        tot = float((ezz * E2).sum())
        mixed = (an > 1e-3) & (an < 1 - 1e-3)
        pure = an >= 1 - 1e-12
        eta = float((an * e_a * E2).sum() / tot)
        low = float((e_a * E2)[pure].sum() / tot)
        mix = float((ezz * E2)[mixed].sum() / tot)
        out.update({f"eta_hf{tag}": eta, f"eta_hf_low{tag}": low, f"eta_hf_high{tag}": low + mix, f"mixed_energy_fraction_hf{tag}": mix,
                    f"mixed_point_fraction_hf{tag}": float(mixed.mean()), f"S_hf{tag}": lam * eta / n, f"S_hf_low{tag}": lam * low / n,
                    f"S_hf_high{tag}": lam * (low + mix) / n})
    if S_measured is not None:
        for tag in ("", "_centred"):
            out[f"dev_hf{tag}_pct"] = 100.0 * (out[f"S_hf{tag}"] / S_measured - 1.0)
        out["dev_hf_low_pct"] = 100.0 * (out["S_hf_low"] / S_measured - 1.0)
        out["S_measured_used"] = S_measured
    return out


def _interior(d: dict, rec: dict) -> tuple[float, float]:
    if "Lx" in d:
        return float(d["Lx"]), float(d["Ly"])
    pml = float(d.get("pml", d.get("dpml", 1.0)))
    return float(d["sx"]) - 2 * pml, float(d["sy"]) - 2 * pml


def eta_estimators(rec: dict, S_measured: float | None = None) -> dict[str, Any]:
    d = load_map(rec)
    s, p = rec["structure"], rec["params"]
    eps = np.asarray(d["eps"], dtype=np.float64)
    Ez = np.asarray(d["Ez"])
    n0, n1 = min(Ez.shape[0], eps.shape[0]), min(Ez.shape[1], eps.shape[1])
    Ez, eps = Ez[:n0, :n1], eps[:n0, :n1]
    n = float(p["analyte"]["n"])
    eps_rod = float(s["lattice"]["rod_eps"])
    eps_a = n * n
    lam = float(d.get("wavelength_nm", rec["result"].get("wavelength_nm")))
    E2 = np.abs(Ez.astype(np.complex128)) ** 2
    W = eps * E2
    tot = float(W.sum())
    threshold = 0.5 * (eps_a + eps_rod)
    eta_mask = float(W[eps < threshold].sum() / tot)
    w = np.clip((eps_rod - eps) / (eps_rod - eps_a), 0.0, 1.0)
    pure, solid = w >= 0.999, w <= 0.001
    mixed = ~(pure | solid)
    eta_low = float(W[pure].sum() / tot)
    frac_mixed = float(W[mixed].sum() / tot)
    eta_high = eta_low + frac_mixed
    eta_lin = float((w * eps_a * E2).sum() / tot)
    # exact circle-pixel geometry on a refined subgrid (the third estimator)
    Lx, Ly = _interior(d, rec)
    nx, ny = eps.shape
    x = np.linspace(-Lx / 2, Lx / 2, nx)
    y = np.linspace(-Ly / 2, Ly / 2, ny)
    dx, dy = x[1] - x[0], y[1] - y[0]
    SUB = 7
    ox = (np.arange(SUB) + 0.5) / SUB - 0.5
    XX, YY = np.meshgrid(x, y, indexing="ij")
    inside = np.zeros((nx, ny), float)
    from .. import geometry as geo
    g = geo.build(s, p)
    for cx, cy, rad in g.rods:
        sel = (np.abs(XX - cx) < rad + dx) & (np.abs(YY - cy) < rad + dy)
        if not sel.any():
            continue
        xs, ys = XX[sel], YY[sel]
        acc = np.zeros(xs.shape, float)
        for a_ in ox:
            for b_ in ox:
                acc += ((xs + a_ * dx - cx) ** 2 + (ys + b_ * dy - cy) ** 2) < rad ** 2
        inside[sel] += acc / SUB ** 2
    inside = np.clip(inside, 0, 1)
    w_geo = 1.0 - inside
    eta_geo = float((w_geo * eps_a * E2).sum() / tot)          # over the stored-map denominator, as in the original analysis
    eps_model = w_geo * eps_a + (1 - w_geo) * eps_rod
    rms_eps = float(np.sqrt(((eps_model - eps) ** 2).mean()))
    S = lambda e: lam * e / n      # noqa: E731
    out = dict(label=rec["task"]["label"], field_file=rec["result"]["field_file"], n=n, wavelength_nm=lam, shape=[int(n0), int(n1)],
               mixed_pixel_fraction_by_count=float(mixed.mean()), mixed_energy_fraction=frac_mixed,
               eta_mask=eta_mask, eta_lin=eta_lin, eta_pure=eta_low, eta_low=eta_low, eta_high=eta_high, eta_geo=eta_geo,
               rms_eps_model=rms_eps, S_mask=S(eta_mask), S_lin=S(eta_lin), S_pure=S(eta_low), S_low=S(eta_low), S_high=S(eta_high),
               S_geo=S(eta_geo), eta_recorded=rec["result"].get("analyte_energy_fraction"))
    s_ref = S_measured if (S_measured is not None and abs(n - 1.33) < 1e-6) else None
    if s_ref is not None:
        out["S_measured"] = S_measured
        out["eta_required_by_S"] = S_measured * n / lam
        out["dev_low_pct"] = 100 * (out["S_low"] / S_measured - 1)
        out["dev_high_pct"] = 100 * (out["S_high"] / S_measured - 1)
        out["dev_lin_pct"] = 100 * (out["S_lin"] / S_measured - 1)
        out["dev_mask_pct"] = 100 * (out["S_mask"] / S_measured - 1)
        out["dev_geo_pct"] = 100 * (out["S_geo"] / S_measured - 1)
    # the estimator consistent with the discretisation (module docstring)
    hf = hellmann_feynman(rec, d, s_ref)
    if hf:
        out.update({k: v for k, v in hf.items() if k not in ("wavelength_nm", "n")})
    return out


def row_decay(rec: dict, pwe: dict | None = None, window: float = 9.0) -> dict[str, Any]:
    """Per-row decay of the stored field: barrier rows at k_x = beta, cladding rows at k_x = 0, complex k_y fits."""
    d = load_map(rec)
    s = rec["structure"]
    eps = np.asarray(d["eps"], dtype=np.float64)
    Ez = np.asarray(d["Ez"])
    nx, ny = eps.shape
    Lx, Ly = _interior(d, rec)
    x = np.linspace(-Lx / 2, Lx / 2, nx)
    y = np.linspace(-Ly / 2, Ly / 2, ny)
    f_r = float(d.get("f_res", d.get("f_rez", rec["params"]["field"]["f_res"])))
    lam = float(d.get("wavelength_nm", rec["result"].get("wavelength_nm")))
    n = float(rec["params"]["analyte"]["n"])
    E2 = np.abs(Ez) ** 2
    ch = {}
    beta = math.pi * 0.5464741196463214
    if pwe:
        k = pwe_mod.kappa_block(pwe, "reference")
        ch = dict(beta_over_pi=k["beta_over_pi"], kappa_slow=k["kappa_pred"], kappa_edge=k["kappa_pred_zone_edge"],
                  kappa_kx0=k["kappa_at_kx0"], kappa_kxpi=k["kappa_at_kxpi"])
        beta = math.pi * k["beta_over_pi"]
    cav = int(s["defect"]["row"])
    iy = lambda v: int(np.argmin(abs(y - v)))      # noqa: E731
    m = np.abs(x) <= window
    wnd = np.hanning(int(m.sum()))
    xs = x[m]
    rows = {}
    jmax = int(np.floor(Ly / 2 - 0.6))
    for j in range(0, jmax + 1):
        sl = Ez[m, iy(j)]
        rows[j] = dict(U=float((eps[:, np.abs(y - j) < 0.5] * E2[:, np.abs(y - j) < 0.5]).sum()),
                       cbeta=complex(np.sum(sl * wnd * np.exp(-1j * beta * xs))), c0=complex(np.sum(sl * wnd)))

    def decay(js, key):
        v = [rows[j][key] if key == "U" else abs(rows[j][key]) ** 2 for j in js]
        return [math.log(v[i] / v[i + 1]) for i in range(len(v) - 1)]

    bar = list(range(cav - 1, 0, -1))
    cl = list(range(cav + 1, jmax))
    out: dict[str, Any] = dict(label=rec["task"]["label"], n=n, f_res=f_r, wavelength_nm=lam, beta_over_pi=beta / math.pi, window_a=window,
                               rows=len(rows), channels=ch,
                               barrier=dict(rows=bar, U=[rows[j]["U"] for j in bar], decay_U=decay(bar, "U"), decay_kbeta=decay(bar, "cbeta"),
                                            decay_k0=decay(bar, "c0"), phase_kbeta=[float(np.angle(rows[j]["cbeta"])) for j in bar]),
                               cladding=dict(rows=cl, U=[rows[j]["U"] for j in cl], decay_U=decay(cl, "U"), decay_kbeta=decay(cl, "cbeta"),
                                             decay_k0=decay(cl, "c0")))
    for name, js, key in (("cladding_kbeta", cl[2:-1], "cbeta"), ("cladding_k0", cl[2:-1], "c0"), ("barrier_kbeta", bar, "cbeta")):
        if len(js) < 3:
            continue
        c = np.array([rows[j][key] for j in js], complex)
        jj = np.array(js, float)
        A = np.vstack([jj, np.ones_like(jj)]).T
        lnc = np.log(c)
        lnc = lnc.real + 1j * np.unwrap(lnc.imag)
        sol, *_ = np.linalg.lstsq(A, lnc, rcond=None)
        ky = sol[0] / 1j
        out[name] = dict(rows=js, re_ky_over_pi=float(ky.real / math.pi), two_im_ky=float(-2 * ky.imag), rms=float(np.abs(lnc - A @ sol).std()))
    return out


def analyse(recs: list[dict], pwe: dict | None = None, S_measured: float | None = None) -> dict[str, Any]:
    fields = sorted(select(recs, mode="field"), key=lambda r: r["params"]["analyte"]["n"])
    out: dict[str, Any] = dict(eta=[], row_decay=[])
    for r in fields:
        try:
            out["eta"].append(eta_estimators(r, S_measured))
            out["row_decay"].append(row_decay(r, pwe))
        except (FileNotFoundError, KeyError) as e:
            out.setdefault("skipped", []).append(dict(label=r["task"]["label"], reason=str(e)))
    at133 = [e for e in out["eta"] if abs(e["n"] - 1.33) < 1e-6]
    out["reference"] = at133[0] if at133 else None
    return out


__all__ = ["load_map", "meep_smoothing", "hellmann_feynman", "eta_estimators", "row_decay", "analyse"]
