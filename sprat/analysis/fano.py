"""Fano line-shape fits of the transmission spectra and the coupling decomposition.

    T(lambda) = base + T0 (q + e)^2 / (1 + e^2),   e = 2 (lambda - lambda_r) / fwhm.

A spectrum is normalised by the cavity-less reference run at the same analyte index,
T = flux_out / flux_out_reference (the input-plane normalisation flux_out / flux_in is also
reported: it collapses at resonance because the input plane sees the reflected wave).
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy.optimize import curve_fit

from .collect import best, select


def fano(lam, base, T0, q, lam_r, fwhm):
    e = 2.0 * (lam - lam_r) / fwhm
    return base + T0 * (q + e) ** 2 / (1.0 + e ** 2)


def fit_fano(lam: np.ndarray, T: np.ndarray, lam0: float, fwhm0: float, window: float = 8.0) -> dict[str, Any]:
    """Fit within |lambda - lam0| < window * fwhm0 (the harmonic-inversion estimate seeds the fit)."""
    m = np.abs(lam - lam0) < window * fwhm0
    if m.sum() < 8:
        m = np.ones_like(lam, dtype=bool)
    p, cov = curve_fit(fano, lam[m], T[m], p0=[0.0, 1.0, 0.0, lam0, fwhm0], maxfev=200000)
    err = np.sqrt(np.abs(np.diag(cov)))
    resid = T[m] - fano(lam[m], *p)
    r2 = 1.0 - float(resid.var() / T[m].var()) if T[m].var() > 0 else float("nan")
    return dict(base=float(p[0]), T0=float(p[1]), q=float(p[2]), lambda_r=float(p[3]), fwhm=float(abs(p[4])),
                error_lambda_r=float(err[3]), error_fwhm=float(err[4]), R2=r2, Q=float(p[3] / abs(p[4])), points=int(m.sum()))


def coupling_branches(Q_loaded: float, T_min: float) -> dict[str, dict[str, float]]:
    """Q_coupling and Q_intrinsic from the loaded Q and the notch depth (both branches).

        1/Q_L = 1/Q_i + 1/Q_c,   T_min = ((Q_c - Q_i)/(Q_c + Q_i))^2.
    """
    t = float(np.sqrt(max(min(T_min, 1.0), 0.0)))
    out = {}
    for name, ratio in (("undercoupled", (1 - t) / (1 + t)), ("overcoupled", (1 + t) / (1 - t) if t < 1 else float("inf"))):
        if not np.isfinite(ratio) or ratio <= 0:
            continue
        Qc = Q_loaded * (1 + ratio) / ratio
        out[name] = dict(Q_coupling=float(Qc), Q_intrinsic=float(ratio * Qc), ratio_Qi_over_Qc=float(ratio))
    return out


def analyse_spectra(recs: list[dict], a_nm: float | None = None, termination: str | None = "pml") -> list[dict[str, Any]]:
    """Every spectrum record with its reference at the same analyte index and the harminv record of the same geometry.

    Only the spectra of the given guide termination are analysed (default: the PML-terminated spectra of the paper,
    section 3.5); the absorber-terminated spectra are compared with them in ``termination.analyse``."""
    out = []
    spectra = select(recs, mode="spectrum") if termination is None else select(recs, mode="spectrum", termination_x=termination)
    for sp in sorted(spectra, key=lambda r: r["params"]["analyte"]["n"]):
        s, p = sp["structure"], sp["params"]
        n = float(p["analyte"]["n"])
        a = a_nm or float(s["lattice"]["a_nm"])
        ref = best(recs, mode="reference", n=n, cladding_rows=s["cell"]["cladding_rows"], guide_periods=s["cell"]["guide_periods"],
                   resolution=p["numerics"]["resolution"], termination_x=s["cell"]["termination_x"])
        h = best(recs, mode="harminv", n=n, radius=float(s["defect"]["radius"]), row=s["defect"]["row"], dx=float(s["defect"]["dx"]),
                 dy=float(s["defect"]["dy"]), cladding_rows=s["cell"]["cladding_rows"], resolution=p["numerics"]["resolution"],
                 termination_x=s["cell"]["termination_x"], guide_periods=s["cell"]["guide_periods"])
        f = np.array(sp["result"]["f"], dtype=float)
        df_grid = float(np.median(np.abs(np.diff(np.sort(f))))) if f.size > 1 else float("nan")
        Po = np.array(sp["result"]["flux_out"], dtype=float)
        Pi = np.array(sp["result"]["flux_in"], dtype=float)
        lam = a / f
        order = np.argsort(lam)
        lam, Po, Pi = lam[order], Po[order], Pi[order]
        normalisation = "input plane"
        T = Po / Pi
        if ref is not None:
            Pn = np.array(ref["result"]["flux_out"], dtype=float)[order]
            if Pn.shape == Po.shape:
                T = Po / Pn
                normalisation = "cavity-less reference"
        T_input = Po / Pi
        if h is not None:
            from ..records import strongest_mode
            m = strongest_mode(h)
            lam_h, Q_h = float(m["wavelength_nm"]), float(m["Q"])
        else:
            i0 = int(np.argmin(T))
            lam_h, Q_h = float(lam[i0]), 2000.0
        fwhm_h = lam_h / Q_h
        fit = fit_fano(lam, T, lam_h, fwhm_h)
        fit_input = fit_fano(lam, T_input, lam_h, 0.5 * fwhm_h)
        far = np.abs(lam - lam_h) > 15 * fwhm_h
        near = ~far
        lam_fine = np.linspace(fit["lambda_r"] - 5 * fit["fwhm"], fit["lambda_r"] + 5 * fit["fwhm"], 4001)
        T_fit_min = float(np.clip(fano(lam_fine, fit["base"], fit["T0"], fit["q"], fit["lambda_r"], fit["fwhm"]).min(), 0.0, None))
        T_far = float(np.median(T[far])) if far.sum() else float("nan")
        T_far_eff = T_far if (np.isfinite(T_far) and T_far > 1e-6) else 1.0
        extinction_dB = -10.0 * math.log10(max(T_fit_min / T_far_eff, 1e-12))
        out.append(dict(
            label=sp["task"]["label"], n=n, normalisation=normalisation, reference_label=(ref["task"]["label"] if ref else ""),
            harminv_label=(h["task"]["label"] if h else ""), points=int(len(lam)),
            lambda_nm=lam.tolist(), T=T.tolist(), T_input=T_input.tolist(), fit=fit, fit_input_normalisation=fit_input,
            lambda_r_harminv=lam_h, Q_harminv=Q_h, fwhm_harminv=fwhm_h,
            fwhm_ratio_harminv_over_fit=float(fwhm_h / fit["fwhm"]),
            fwhm_ratio_harminv_over_input_fit=float(fwhm_h / fit_input["fwhm"]),
            lambda_r_difference_nm=float(abs(fit["lambda_r"] - lam_h)),
            T_min=float(T.min()), T_fit_min=T_fit_min, T_far=T_far, extinction_dB=float(extinction_dB),
            frequency_step=df_grid, T_min_neighbours=_neighbours(T), T_max_near=float(T[near].max()) if near.any() else None,
            far_window_points=int(far.sum()),
            coupling=coupling_branches(fit["Q"], T_fit_min / T_far_eff)))
    return out


def _neighbours(T: np.ndarray) -> list[float]:
    """Transmission at the two samples next to the smallest one: how far the frequency grid resolves the dip."""
    i = int(np.argmin(T))
    return [float(T[j]) for j in (i - 1, i + 1) if 0 <= j < T.size]


def summary(spectra: list[dict]) -> dict[str, Any]:
    if not spectra:
        return {}
    r = np.array([s["fwhm_ratio_harminv_over_fit"] for s in spectra])
    ri = np.array([s["fwhm_ratio_harminv_over_input_fit"] for s in spectra])
    tmin = np.array([s["T_min"] for s in spectra])
    # the inverse of equation (3) on the raw notch depth: Q_i = Q_L / sqrt(T_min), which scatters with the numerical floor
    qi = np.array([s["fit"]["Q"] / np.sqrt(s["T_min"]) for s in spectra if s["T_min"] > 0])
    agree = [100.0 * (s["fit"]["fwhm"] / s["fwhm_harminv"] - 1.0) for s in spectra]
    return dict(count=len(spectra), n=[s["n"] for s in spectra], q=[s["fit"]["q"] for s in spectra],
                abs_q_max=float(max(abs(s["fit"]["q"]) for s in spectra)), abs_q_min=float(min(abs(s["fit"]["q"]) for s in spectra)),
                fwhm_ratio_min=float(r.min()), fwhm_ratio_max=float(r.max()),
                fwhm_ratio_input_min=float(ri.min()), fwhm_ratio_input_max=float(ri.max()),
                T_far_min=float(np.nanmin([s["T_far"] for s in spectra])), T_far_max=float(np.nanmax([s["T_far"] for s in spectra])),
                lambda_r_difference_max_nm=float(max(s["lambda_r_difference_nm"] for s in spectra)),
                agreement_pct=agree, agreement_max_abs_pct=float(max(abs(a) for a in agree)),
                T_min_min=float(tmin.min()), T_min_max=float(tmin.max()),
                frequency_step=float(np.median([s["frequency_step"] for s in spectra])),
                T_min_neighbour_min=float(min(min(s["T_min_neighbours"]) for s in spectra)),
                T_min_neighbour_max=float(max(max(s["T_min_neighbours"]) for s in spectra)),
                T_max_near_max=float(max(s["T_max_near"] for s in spectra if s["T_max_near"] is not None)),
                spectra_with_far_window=int(sum(1 for s in spectra if s["far_window_points"] > 0)),
                Q_intrinsic_from_T_min_min=(float(qi.min()) if qi.size else None), Q_intrinsic_from_T_min_max=(float(qi.max()) if qi.size else None),
                points=sorted({int(s["points"]) for s in spectra}))


__all__ = ["fano", "fit_fano", "coupling_branches", "analyse_spectra", "summary"]
