"""Quantities that combine several analysis blocks; they are computed after all of them.

* the quadratic law of the n_cl = 8 series applied to the fine sweep of the reference cell (section 3.2), and its
  derivative at the reference radius against the five-point stencil (section 3.8);
* the cell-factor estimates of section 3.4: d lnG/df = (d lnG/d r_d) / (df/d r_d), with df/d r_d measured on the fine
  sweep of the reference cell at 0.060a, applied to the frequency change over the barrier-row series and at
  |dy| = 0.10a;
* the displacement coefficient over the row coefficient, the raw 25a-cell kappa at the second radius against the
  reflectionless one, the three figures of merit, and the analyte energy fraction of the second field map (the
  Hellmann-Feynman estimate, and the stored-map pure-pixel value that version 8 read as a lower bound) against the
  measured sensitivity extrapolated to its analyte index.
"""

from __future__ import annotations

import math
from typing import Any


def _get(d, *path):
    for p in path:
        if not isinstance(d, dict) or p not in d:
            return None
        d = d[p]
    return d


def combine(A: dict[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    bs, rs, ba = A.get("base_series") or {}, A.get("radius_sweep") or {}, A.get("barrier") or {}
    sw, fe, disp = A.get("analyte_sweep") or {}, A.get("fields") or {}, A.get("displacement") or {}
    qf = bs.get("quadratic_fit")
    if qf and rs.get("rd"):
        c = qf["coefficients_ascending"]
        law = lambda x: c[0] + c[1] * x + c[2] * x * x          # noqa: E731
        dev = [f - law(x) for x, f in zip(rs["rd"], rs["f"])]
        i = max(range(len(dev)), key=lambda j: abs(dev[j]))
        out["law_transfer_to_fine_sweep"] = dict(points=len(dev), max_abs_deviation=abs(dev[i]), at_radius=rs["rd"][i],
                                                 deviations=dev, radius_range=[rs["rd"][0], rs["rd"][-1]])
        if rs.get("local_dlambda_drd") and 0.06 in [round(x, 4) for x in rs["rd"]]:
            f0 = rs["f"][[round(x, 4) for x in rs["rd"]].index(0.06)]
            d_law = -(c[1] + 2 * c[2] * 0.06) / f0 ** 2
            st = rs["local_dlambda_drd"]["nm_per_nm"]
            out["dlambda_drd_at_reference"] = dict(law_nm_per_nm=d_law, stencil_nm_per_nm=st, agreement_pct=100.0 * (st / d_law - 1.0),
                                                   note="law: -F'(r_d) / f_r^2 with F the quadratic law of the n_cl = 8 series")
    bsa, rsa = A.get("base_series_all_points") or {}, A.get("radius_sweep_all_points") or {}
    if bsa.get("quadratic_fit") and rsa.get("rd"):
        c = bsa["quadratic_fit"]["coefficients_ascending"]
        dev = [f - (c[0] + c[1] * x + c[2] * x * x) for x, f in zip(rsa["rd"], rsa["f"])]
        i = max(range(len(dev)), key=lambda j: abs(dev[j]))
        out["law_transfer_all_points"] = dict(points=len(dev), max_abs_deviation=abs(dev[i]), at_radius=rsa["rd"][i])
        if 0.06 in [round(x, 4) for x in rsa["rd"]] and rsa.get("local_dlambda_drd"):
            f0 = rsa["f"][[round(x, 4) for x in rsa["rd"]].index(0.06)]
            out["dlambda_drd_all_points"] = dict(law_nm_per_nm=-(c[1] + 2 * c[2] * 0.06) / f0 ** 2,
                                                 stencil_nm_per_nm=rsa["local_dlambda_drd"]["nm_per_nm"])
    dlnG_drd = _get(ba, "tolerance", "dlnG_drd_per_a")
    # df/dr_d at the reference radius, measured on the fine sweep of the same cell (five-point stencil of lambda_r):
    # df/d(r_d/a) = -f_r^2 dlambda_r/dr_d; the linear law of the n_cl = 8 series, which the paper rejects, is kept only
    # as the fallback it was in SPRAT 1.0
    df_drd, df_src = None, None
    st = _get(rs, "local_dlambda_drd", "nm_per_nm")
    if st and rs.get("rd") and 0.06 in [round(x, 4) for x in rs["rd"]]:
        f0 = rs["f"][[round(x, 4) for x in rs["rd"]].index(0.06)]
        df_drd, df_src = -f0 * f0 * st, "fine sweep of the reference cell, five-point stencil at 0.060a"
    if df_drd is None:
        df_drd, df_src = _get(bs, "fits_by_row", 4, "linear", "slope"), "linear law of the n_cl = 8 series at N_sep = 4"
    if dlnG_drd is not None and df_drd:
        dlnG_df = dlnG_drd / df_drd
        ce: dict[str, Any] = dict(dlnG_drd_per_a=dlnG_drd, df_drd_per_a=df_drd, df_drd_source=df_src, dlnG_df=dlnG_df)
        brs = ba.get("barrier_row_series")
        if brs:
            ce["barrier_row"] = dict(f_change=brs["f_change"], dlnG=dlnG_df * brs["f_change"], change_pct=100.0 * dlnG_df * brs["f_change"])
        e = _get(disp, "dy", 0.06, "at_max_abs_dy")
        if e:
            ce["dy_max"] = dict(dy=e["dy"], f_change=e["f_change"], dlnG=dlnG_df * e["f_change"], change_pct=100.0 * dlnG_df * e["f_change"])
        # the axial displacement: the resonance shift at the largest dx moves the cell factor; the intrinsic change of Q is
        # the measured change of the 25a cell less that estimate
        for dyv in (0.0, 0.1):
            ex = _get(disp, "dx", dyv)
            if ex and ex.get("dx") and 0.0 in ex["dx"] and ex.get("Q0"):
                i0, i1 = ex["dx"].index(0.0), len(ex["dx"]) - 1
                fch = ex["f"][i1] - ex["f"][i0]
                meas = math.log(ex["Q"][i1] / ex["Q0"])
                ce.setdefault("dx_max", {})["%.2f" % dyv] = dict(dx=ex["dx"][i1], f_change=fch, dlnG=dlnG_df * fch,
                                                                  change_pct=100.0 * dlnG_df * fch, measured_pct=100.0 * meas,
                                                                  intrinsic_pct=100.0 * (meas - dlnG_df * fch))
        out["cell_factor_estimates"] = ce
    kdy = _get(disp, "dy", 0.06, "lnQ_fit", "kappa_dy_per_a")
    kap = _get(ba, "series", "kappa")
    if kdy is not None and kap:
        out["displacement_over_row_coefficient"] = dict(kappa_dy=kdy, kappa=kap, ratio=kdy / kap)
    rp, sr = A.get("resolved_pair"), ba.get("second_radius")
    if rp and sr:
        out["second_radius_raw_against_reflectionless"] = dict(raw=rp["kappa"], reflectionless=sr["kappa_absorber"],
                                                              deviation_pct=100.0 * (rp["kappa"] / sr["kappa_absorber"] - 1.0))
    pf = ba.get("performance")
    fom_cell = _get(sw, "at_1p33", "FOM")
    if pf and fom_cell:
        out["figures_of_merit"] = dict(cell=fom_cell, reflectionless=pf["FOM_lossless_N4"], water=pf["FOM_water_N4"],
                                       cell_over_water=fom_cell / pf["FOM_water_N4"])
    se = sw.get("S_extrapolated")
    maps = [e for e in (fe.get("eta") or []) if se and abs(e["n"] - se["n"]) < 1e-6]
    if maps:
        e = maps[0]
        out["second_map_bound"] = dict(n=e["n"], S_low=e["S_low"], S_measured_extrapolated=se["S"], deviation_pct=100.0 * (e["S_low"] / se["S"] - 1.0))
        if e.get("S_hf") is not None:
            out["second_map_hf"] = dict(n=e["n"], S_hf=e["S_hf"], S_hf_centred=e.get("S_hf_centred"), S_measured_extrapolated=se["S"],
                                        deviation_pct=100.0 * (e["S_hf"] / se["S"] - 1.0),
                                        deviation_centred_pct=(100.0 * (e["S_hf_centred"] / se["S"] - 1.0) if e.get("S_hf_centred") else None),
                                        S_hf_low=e.get("S_hf_low"), S_hf_high=e.get("S_hf_high"))
    cl = _get(A, "convergence", "cladding_10_vs_12_at_res24")
    if cl:
        v = [abs(x["change_pct"]) for x in cl]
        out["cladding_uncertainty_pct"] = [min(v), max(v)]
    if pf and kap:
        out["water_reduction_of_reflectionless_Q_pct"] = _get(pf, "absorption", "Q_reduction_by_water_pct")
    pv = physical_values(A)
    if pv:
        out["physical_values"] = pv
    return out


TABLE5_N = (1.305, 1.330, 1.350, 1.375, 1.400, 1.425, 1.445)


def physical_values(A: dict[str, Any]) -> dict[str, Any] | None:
    """The reflectionless quality factor and figure of merit across the analyte sweep (SPRAT 1.2.1).

    Q_w at an analyte index is taken, in this order of preference:
      * at n_a = 1.33, the harmonic-inversion value of the reference point with the guide continued into the absorber
        and the cladding leak removed (the series Q_w(4) of the barrier analysis, section 3.3);
      * at an index of the absorber spectra (amendment A2 of the criterion of the absorber spectra), Q_w = Q(harmonic
        inversion, 25a cell) / G with G from the spectra;
      * between two such indices, linear interpolation in n_a of that Q_w, which varies smoothly with the index (G does
        not: it follows the round-trip phase of the guide).
    FOM_w = S Q_w / lambda_r with the local sensitivity and the resonance of the analyte sweep. The curve of the figure
    uses the spectral values and their interpolation at every point of the sweep; the table uses the rule above.
    """
    sw, ba, te = A.get("analyte_sweep") or {}, A.get("barrier") or {}, A.get("termination") or {}
    a2 = _get(te, "A2", "rows")
    if not sw.get("n_c") or not a2:
        return None
    nodes = sorted((float(k), float(v["Q_w"])) for k, v in a2.items())
    xs = [x for x, _ in nodes]

    def interp(n: float) -> tuple[float, str] | None:
        for x, q in nodes:
            if abs(n - x) < 1e-9:
                return q, "spectra"
        for (x0, q0), (x1, q1) in zip(nodes, nodes[1:]):
            if x0 < n < x1:
                return q0 + (n - x0) / (x1 - x0) * (q1 - q0), "interpolated"
        return None

    nc, Sc, lc, Qc = (list(map(float, sw[k])) for k in ("n_c", "S_c", "lambda_c", "Q_c"))
    curve = dict(n=[], Q_w=[], FOM_w=[])
    for n, S, lam in zip(nc, Sc, lc):
        v = interp(n)
        if v is None:
            continue
        curve["n"].append(n)
        curve["Q_w"].append(v[0])
        curve["FOM_w"].append(S * v[0] / lam)
    qw_ref = _get(ba, "series", "Q_w", "4")
    if qw_ref is None:
        qw_ref = _get(ba, "series", "Q_w", 4)
    rows = []
    for n0 in TABLE5_N:
        idx = [i for i, n in enumerate(nc) if abs(n - n0) < 1e-9]
        if not idx:
            continue
        i = idx[0]
        if abs(n0 - 1.33) < 1e-9 and qw_ref is not None:
            q, how = float(qw_ref), "harmonic inversion, absorber, leak removed"
        else:
            v = interp(n0)
            if v is None:
                continue
            q, how = v
        rows.append(dict(n=n0, lambda_nm=lc[i], S=Sc[i], Q_cell=Qc[i], FOM_cell=Sc[i] * Qc[i] / lc[i], Q_w=q, FOM_w=Sc[i] * q / lc[i],
                         cell_factor=Qc[i] / q, source=how))
    out: dict[str, Any] = dict(rule="Q_w: harmonic inversion with the absorber and the leak removed at n_a = 1.33; Q(25a cell)/G of the "
                                    "absorber spectra at their indices; linear interpolation in n_a of that Q_w between them",
                               nodes=dict(n=xs, Q_w=[q for _, q in nodes]), curve=curve, table=rows)
    v133 = interp(1.33)
    if v133 and qw_ref:
        out["at_1p33"] = dict(Q_w_harminv=float(qw_ref), Q_w_spectra_interpolated=v133[0],
                              deviation_pct=100.0 * (v133[0] / float(qw_ref) - 1.0))
    return out


__all__ = ["combine"]
