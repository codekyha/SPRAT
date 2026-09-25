"""The numbers registry: the record-derived numbers of the paper with the record or the derivation behind each.

Rule of the paper: no number enters the manuscript without a row in the registry.  Measured rows carry the record
label, host, job id and the SHA-256 of the code that produced the record; derived rows carry
the derivation.  The quantity names are those of the deposited registry (version 2.0.2 of the
data set) and of the manuscript registry, so that ``sprat audit`` can compare them row by row.

SPRAT 1.1 (manuscript version 9) changes four things, and the rows follow them:

* the complex band structure keeps physical roots only: the slower value of the original plane-wave layer (1.9015) is a
  truncation artefact of the square plane-wave basis and is reported as such, the zone-edge channel is the one
  physical channel at the operating point, and the cladding channels at k_x = 0 and pi/a are the physical ones;
* fits and statistics use the records admitted by the validity criterion only (the r_d = 0.120a points of the
  n_cl = 8 series and the r_d = 0.110a point of the fine sweep are excluded);
* the analyte energy fraction is the Hellmann-Feynman estimate over the permittivity E_z sees (``fields``); the
  stored-map values of the original analysis keep names that say what they are;
* the absorption-limited quality factor follows from the measured sensitivity, Q_abs = lambda_r / (2 k S).

Statements of the manuscript that these changes supersede are listed by ``superseded``.

SPRAT 1.2 (manuscript version 10) adds the section ``termination``: the transmission spectra with the guide continued
into the absorber compared with those of the 25a cell (``termination.analyse``), graded by the criterion fixed before
the runs and by its amendment, and the cell factor across the analyte sweep. The dataset section gains the share of the
harmonic-inversion records that are admitted and hold a cavity mode. Two quantities are renamed: "records from the
systematic runs" (earlier "records from the parameter campaign") and "eta_a, threshold mask of the original analysis"
(earlier "eta_a, threshold mask of the campaign").

SPRAT 1.2.1 (manuscript version 10.2) adds the section ``reflectionless``: the reflectionless quality factor and figure
of merit at the rows of table 5 (``combined.physical_values``) and the spectral value interpolated to n_a = 1.33 against
the harmonic-inversion one. It also adds the convergence of the resonance frequency with resolution (sections
``systematic_runs`` and ``derived``), the guided-mode wavevector of the barrier channel with the per-row factors the same
frequency gives at k_x = 0 and pi/a (``derived``), and the absorption-limited quality factor with the mixed points of
the grid counted wholly as analyte or wholly as silicon (``absorption``).
"""

from __future__ import annotations

import math
import time
from typing import Any

from .. import __version__


def _M(recs_by_key: dict, key: str, quantity: str, note: str = "") -> dict[str, Any]:
    v = recs_by_key[key]
    return dict(quantity=quantity, value=round(v["Q"], 2) if v["Q"] is not None else None, unit="", kind="measured",
                label=v["label"], host=v.get("host", ""), job_id=v.get("job_id", ""), code_sha256=v.get("code_sha256", ""),
                f_r=(round(v["f"], 7) if v["f"] else None), margin=(round(v["margin"], 2) if v["margin"] else None), note=note)


def _D(quantity: str, value, unit: str, derivation: str, uncertainty=None, note: str = "") -> dict[str, Any]:
    d = dict(quantity=quantity, value=value, unit=unit, kind="derived", derivation=derivation, note=note)
    if uncertainty is not None:
        d["uncertainty"] = uncertainty
    return d


def build(analysis: dict[str, Any], counts: dict[str, Any] | None = None) -> dict[str, Any]:
    ba = analysis.get("barrier") or {}
    sw = analysis.get("analyte_sweep") or {}
    fe = analysis.get("fields") or {}
    inp = ba.get("inputs", {})
    reg: dict[str, Any] = dict(spec="sprat/numbers-registry-1.0", created=time.strftime("%Y-%m-%d"), software=f"sprat {__version__}",
                               rule="No number enters the manuscript without a row here. Measured rows carry the record label, host, "
                                    "job id and code SHA-256; derived rows carry the derivation. Rebuild with `sprat analyze`.",
                               measured=[], derived=[], second_radius=[], field_map=[], systematic_runs=[], dataset=[])
    if not ba:
        return reg
    # ---- measured
    for key, q, note in (("Q4_pml_cl12", "Q, N_sep = 4, PML, n_cl = 12 (clearance 8)", "reference geometry, 25a cell"),
                         ("Q5_pml_cl13", "Q, N_sep = 5, PML, n_cl = 13 (clearance 8)", "clearance-matched to N = 4 and N = 6"),
                         ("Q5_pml_cl12", "Q, N_sep = 5, PML, n_cl = 12 (clearance 7)", "carries the clearance-7 leak"),
                         ("Q6_pml_cl14", "Q, N_sep = 6, PML, n_cl = 14 (clearance 8)", ""),
                         ("A3_abs_cl12", "Q, N_sep = 3, absorber 8a, n_cl = 12 (clearance 9)", "reflectionless"),
                         ("A4_abs_cl12", "Q, N_sep = 4, absorber 8a, n_cl = 12 (clearance 8)", "reflectionless reference"),
                         ("A5_abs_cl13", "Q, N_sep = 5, absorber 8a, n_cl = 13 (clearance 8)", "reflectionless"),
                         ("A5_abs_cl12", "Q, N_sep = 5, absorber 8a, n_cl = 12 (clearance 7)", "gives the leak cross-check"),
                         ("Q4_pml_nx49", "Q at the reference point, PML, 49a guide", "second cell length; G < 1 here"),
                         ("Q5_pml_cl12_res32", "Q, N_sep = 5, resolution 32 (H8)", "resolution check")):
        if inp.get(key, {}).get("label"):
            reg["measured"].append(_M(inp, key, q, note))
    # ---- derived
    s, l, g, tp, mc = ba["series"], ba["leak"], ba["G"], ba["three_point_solve"], ba["monte_carlo"]
    ch = ba["channels"]
    D = reg["derived"]
    D += [
        _D("kappa, barrier decay constant", round(s["kappa"], 4), "", "four-point reflectionless leak-free fit, N_sep = 3-6", round(mc["kappa_sd"], 4),
           "Monte Carlo over 0.23 % independent discretisation error per Q plus the 0.11 % G(N=6) bound"),
        _D("kappa, three-point PML solve", round(tp["kappa"], 4), "", "e^kappa = (1/Q4 - 1/Q5)/(1/Q5 - 1/Q6) at clearance 8; independent of G"),
        _D("kappa, absorber-only series", round(s["kappa_absorber_only"], 4), "", "mean of the two leak-corrected absorber steps N=3->4->5"),
        _D("per-row factor e^kappa", round(s["per_row_factor"], 3), "", "exp of the fitted kappa"),
        _D("A, barrier prefactor", round(s["A"], 4), "", "intercept of the four-point fit, Q_w = A e^{kappa N}"),
        _D("Q_top, clearance 8", round(l["Q_top_8"]), "", "three-point solve, exact", round(mc["Q_top_8_sd"])),
        _D("Q_top, clearance 7", round(l["Q_top_7"]), "", "three-point solve plus the n_cl 12/13 pair at N=5"),
        _D("kappa_top, per cladding row", round(l["kappa_top"], 3), "", "ln[Q_top(8)/Q_top(7)]"),
        _D("leak share of 1/Q at clearance 8", [round(l["leak_share_pct"][k], 2) for k in ("N4_c8", "N5_c8", "N6_c8")], "%",
           "1/Q_top(8) divided by 1/Q at N_sep = 4, 5, 6"),
        _D("leak share of 1/Q at clearance 7, N_sep = 5", round(l["leak_share_pct"]["N5_c7"], 2), "%", "1/Q_top(7) divided by 1/Q at N_sep = 5"),
        _D("G, cell Fabry-Perot factor at the reference point", round(g["4"], 4), "", "Q_PML/Q_absorber after leak removal, N_sep = 4", None,
           "N_sep = 5 gives %.4f (n_cl 13) and %.4f (n_cl 12)" % (g["5"], g["5_cl12"])),
        _D("Q_w, reflectionless leak-free series", {k: round(v, 1) for k, v in s["Q_w"].items()}, "",
           "leak removed from every record, PML records divided by the measured G"),
        _D("per-row steps of the series", [round(x, 4) for x in s["steps"]], "", "ln(Q_{N+1}/Q_N) of the four-point series"),
        _D("kappa_pred, zone-edge channel", round(ch["zone_edge"], 4), "",
           "2 Im k_y a of the slowest physical channel (Re k_y = pi/a) at the guided-mode wavevector, plane-wave cut-off 10", None,
           "deviation of the measurement: %+.3f %%" % s["dev_zone_edge_pct"]),
        _D("deviation of the four-point kappa from the zone-edge channel", round(s["dev_zone_edge_pct"], 2), "%",
           "100 (kappa / kappa_pred - 1)"),
        _D("clearance rule", "c >= 8 + ln(100 Q_w / Q_top(8)) / kappa_top", "", "require Q_top(c) >= 100 Q_w", None,
           "required clearances: " + ", ".join(f"{k} {v:.2f}" for k, v in ba["clearance_rule"]["required"].items())),
    ]
    if ch.get("next") is not None:
        D.append(_D("kappa_pred, next physical channel", round(ch["next"], 4), "",
                    "the slowest physical root after the zone-edge one, at the same f and k_x", None,
                    "Re k_y = %.3f pi/a" % (ch.get("re_ky_next_over_pi") or 0.0)))
    if ch.get("boundary_artefact") is not None:
        D.append(_D("kappa of the slowest truncation artefact of the square plane-wave basis", round(ch["boundary_artefact"], 4), "",
                    "rejected root: eigenvector on the two outermost G_y rings of the basis, no replica with |Re k_y| <= pi/a", None,
                    "the 'least-evanescent complex pair' of the original plane-wave layer"))
    # SPRAT 1.2.1: the channel is evaluated at the guided-mode wavevector; the same frequency at k_x = 0 or pi/a
    if ch.get("beta_over_pi") is not None:
        D.append(_D("guided-mode wavevector beta at the reference resonance", round(ch["beta_over_pi"], 4), "pi/a",
                    "W1 band of the plane-wave record at f_r = 0.30461"))
    if ch.get("per_row_factor_kx0") is not None:
        D.append(_D("per-row factor of the slowest channel at k_x = 0, reference resonance", round(ch["per_row_factor_kx0"], 2), "",
                    "exp(kappa_pred at k_x = 0 = %.4f): the channel evaluated off the guided-mode wavevector" % ch["kx0"]))
    if ch.get("per_row_factor_kxpi") is not None:
        D.append(_D("per-row factor of the slowest channel at k_x = pi/a, reference resonance", round(ch["per_row_factor_kxpi"], 1), "",
                    "exp(kappa_pred at k_x = pi/a = %.4f)" % ch["kxpi"]))
    tc = ch.get("truncation_check") or {}
    if tc.get("bases"):
        B = tc["bases"]
        D.append(_D("zone-edge channel for four plane-wave bases", {k: round(v["kappa_physical"], 4) for k, v in B.items()}, "",
                    "square cut-offs M = 7, 10, 13 and a circular basis |m| <= 16, at the reference resonance"))
        D.append(_D("slowest truncation artefact for four plane-wave bases",
                    {k: (round(v["kappa_slowest_artefact"], 4) if v["kappa_slowest_artefact"] is not None else None) for k, v in B.items()}, "",
                    "same bases; the artefact moves with the shape of the basis, the physical root does not"))
        kp = [v["kappa_physical"] for v in B.values()]
        D.append(_D("spread of the zone-edge channel over the four bases", float("%.2g" % (max(kp) - min(kp))), "", "largest minus smallest"))
    if "G_49a" in ba:
        D.append(_D("G, 49a cell at the reference point", round(ba["G_49a"], 4), "", "Q_PML(49a)/Q_absorber after leak removal", None, "below unity"))
        D.append(_D("spread of the quality factor across the three terminations at the reference point", round(ba["spread_25a_over_49a"], 2), "",
                    "Q(25a cell) / Q(49a cell) at r_d = 0.060a, N_sep = 4, both leak-removed"))
    pc = ba.get("pml_thickness_check")
    if pc and "Q4_pml_cl12_pml2" in pc and "Q5_pml_cl12_pml2" in pc:
        D.append(_D("change of Q with the matched layer doubled from 1a to 2a, N_sep = 4 and 5",
                    [round(abs(pc[k]["change_pct"]), 2) for k in ("Q4_pml_cl12_pml2", "Q5_pml_cl12_pml2")], "%",
                    "100 |Q(PML 2a) / Q(PML 1a) - 1| at the reference geometry and at N_sep = 5, n_cl = 12"))
    rc = ba.get("resolution_check") or {}
    if rc.get("f_change_pct") is not None:                     # SPRAT 1.2.1: the resonance of the resolution pair
        D += [_D("change of f_r from resolution 24 to 32 at the reference radius, N_sep = 5", round(rc["f_change_pct"], 3), "%",
                 "100 (f32 / f24 - 1), records v4_H8_s5 and v4_H1b_s5 (n_cl = 12)"),
              _D("shift of lambda_r from resolution 24 to 32 at the reference radius, scaled to 1550 nm", round(abs(rc["dlambda_nm_at_target"]), 2),
                 "nm", "|a/f32 - a/f24| x lambda_t / lambda_r of the reference point"),
              _D("defect radius of the reference geometry in pixels at resolution 24", round(rc["defect_radius_pixels"], 2), "",
                 "r_d x resolution")]
    mr = ba.get("matched_ratio_fit")
    if mr:
        D += [_D("R, W1 end power reflectance", round(mr["R"], 3), "", "two-reflector fit to the %d matched PML/absorber ratios" % mr["points"], round(mr["R_se"], 3)),
              _D("D, effective reflector distance", round(mr["D"], 2), "a", "same fit", round(mr["D_se"], 2), "geometric W1 end at 12.5a"),
              _D("rms of the two-reflector model in ln G", round(mr["rms_lnG"], 3), "", "%d matched pairs" % mr["points"]),
              _D("absorber scale factor c", round(mr["scale_c"], 3), "", "free overall scale in the same fit", round(mr["scale_c_se"], 3),
                 "consistent with unity: the absorber reproduces an infinite guide")]
        if mr.get("G"):
            D.append(_D("range of the cell factor G over the fine sweep", [round(min(mr["G"]), 2), round(max(mr["G"]), 2)], "",
                        "smallest and largest of the %d matched ratios Q_PML/Q_absorber after leak removal, r_d = %g to %ga"
                        % (mr["points"], min(mr["rd"]), max(mr["rd"])), None,
                        "the 25a cell raises or lowers the quality factor by up to a factor of two"))
    rsw = analysis.get("radius_sweep") or {}
    if rsw.get("peak_trough_spacing_a") is not None:
        D.append(_D("spacing of the peak and the trough of Q over the admitted radii of the fine sweep, 25a cell",
                    round(rsw["peak_trough_spacing_a"], 3), "a", "peak at %s, trough at %s" % (rsw.get("peaks"), rsw.get("troughs")), None,
                    "half the period of the oscillation"))
    pf = ba.get("performance")
    if pf and fe.get("reference"):
        D += [_D("eta_a required by the measured sensitivity", round(pf["eta_used"], 4), "", "n_a S / lambda_r, equation (4)"),
              _D("FOM at the reference point, lossless 2D", round(pf["FOM_lossless_N4"]), "1/RIU", "S Q_w / lambda_r"),
              _D("FOM at the reference point, in water", round(pf["FOM_water_N4"]), "1/RIU", "S Q_tot / lambda_r with 1/Q_tot = 1/Q_w + 1/Q_abs"),
              _D("FOM at N_sep = 5, in water", round(pf["FOM_water_N5"]), "1/RIU", "same")]
        if pf.get("detection_limit_factor_vs_6927"):
            D.append(_D("detection-limit factor against the v3 headline", round(pf["detection_limit_factor_vs_6927"], 3), "",
                        "%.2f / Q_w(4), the 25a-cell Q of the reference record over the reflectionless one" % pf["Q_reference_campaign"]))
        D += [
              _D("Q_abs in water", round(pf["Q_abs"], 1), "", "lambda_r / (2 k_water S) = n_a / (2 k_water eta_a) with eta_a = n_a S / lambda_r = %.4f"
                 % pf["eta_used"]),
              _D("reduction of the reflectionless quality factor by water", round(pf["water_reduction_N4_pct"], 1), "%", "100 (1 - Q_tot / Q_w) at N_sep = 4"),
              _D("Q_tot in water", {k: round(v) for k, v in pf["Q_tot"].items()}, "", "1/Q_tot = 1/Q_w + 1/Q_abs"),
              _D("T_min in water", {k: round(v, 3) for k, v in pf["T_min"].items()}, "", "(Q_tot/Q_abs)^2"),
              _D("coupling metric fraction", {k: round(v, 3) for k, v in pf["metric_fraction"].items()}, "",
                 "M(x)/M(x_opt) with x = Q_w/Q_abs, x_opt = (1+sqrt 3)/2"),
              _D("Q_w at the coupling optimum", round(pf["Q_w_opt"]), "", "x_opt Q_abs")]
    tol = ba.get("tolerance", {})
    if "dlambda_per_row_nm" in tol:
        D.append(_D("d lambda_r / d N_sep at the operating point", round(tol["dlambda_per_row_nm"], 3), "nm",
                    "lambda_r(N=5) - lambda_r(N=4) at r_d = 0.060a"))
    for key, name in (("absorber", "reflectionless"), ("pml_25a", "25a cell")):
        if key in tol:
            t = tol[key]
            D += [_D(f"dQ/dr_d at the reference point, {name}", round(t["dQ_drd_per_nm"], 1), "1/nm",
                     "central difference of the %s sweep over r_d = 0.055-0.065a" % ("absorber" if key == "absorber" else "PML")),
                  _D(f"cross-term, {name}", round(t["cross_term_pct_of_row_step"], 3), "% of a row step",
                     "restoring lambda_r after one row costs d_r_d = %.4f nm; the resulting dQ as a fraction of a row step" % tol["drd_to_restore_lambda_nm"],
                     None, "maximum over the whole sweep: %.2f %%" % t["cross_term_max_pct"]),
                  _D(f"fabrication tolerance, {name}", round(t["fabrication_tolerance_pct_per_nm"], 1), "% per +-1 nm radius error",
                     "|dQ/dr_d| x 1 nm / Q")]
    if "absorber" in tol:
        D.append(_D("d lnQ_w/dr_d at the reference point, reflectionless", round(tol["absorber"]["dlnQ_drd_per_a"], 1), "1/a",
                    "central difference of ln Q over the absorber sweep, r_d = 0.055 to 0.065a"))
    if "dlnG_drd_per_a" in tol:
        D.append(_D("d lnG/dr_d at the reference point", round(tol["dlnG_drd_per_a"], 1), "1/a", "central difference of ln G over the same three radii"))
    do = ba.get("detrended_oscillation")
    if do and "pml_25a" in do and "absorber" in do:
        D.append(_D("suppression of the Q(r_d) oscillation by the absorber", round(do["pml_25a"]["detrended_ptp"] / do["absorber"]["detrended_ptp"], 1), "",
                    "ratio of the detrended peak-to-peak in ln Q, PML 25a cell over absorber"))
        D += [_D("detrended peak-to-peak of ln Q over the fine sweep, 25a cell", round(do["pml_25a"]["detrended_ptp"], 3), "",
                 "cubic trend removed, admitted radii"),
              _D("detrended peak-to-peak of ln Q over the fine sweep, absorber", round(do["absorber"]["detrended_ptp"], 3), "", "same")]
    knob_ = ba.get("barrier_row_knob") or {}
    lr = {k: math.log(v["ratio"]) / math.log(v["predicted_ratio"]) for k, v in sorted(knob_.items(), key=lambda t: float(t[0]))
          if v.get("predicted_ratio") and abs(math.log(v["predicted_ratio"])) > 1e-12}
    if lr:
        D.append(_D("measured over predicted ln Q ratio for the barrier-row radius", [round(v, 2) for v in lr.values()], "",
                    "ln(Q(r_b)/Q_ref) over the plane-wave prediction of the physical channel, r_b = " + ", ".join(lr)))
    if inp.get("Q4_pml_cl12", {}).get("Q") and l.get("Q_top_8"):
        D.append(_D("T_min expected from the cladding leak at the reference geometry", float("%.2g" % ((inp["Q4_pml_cl12"]["Q"] / l["Q_top_8"]) ** 2)), "",
                    "(Q_L / Q_i)^2 with Q_i = Q_top(8), the only intrinsic loss of the two-dimensional model"))
    pwc = (analysis.get("band_gap") or {}).get("plane_wave_check") or {}
    if pwc.get("max_relative_difference_pct") is not None:
        D.append(_D("largest relative difference of the gap edges, MPB against the plane-wave layer", round(pwc["max_relative_difference_pct"], 3), "%",
                    "at n_a = 1.33"))
    n8 = ba.get("ncl8_leak_correction")
    if n8:
        D += [_D("leak-corrected per-row step of the n_cl = 8 series, radius by radius", [round(x, 3) for x in n8["corrected"]], "",
                 "ln(Q4/Q3) after subtracting 1/Q_top(c) at clearances 5 and 4, for r_d = " + ", ".join(f"{x:g}" for x in n8["rd"])),
              _D("mean per-row step of the n_cl = 8 series, raw", round(n8["raw_mean"], 3), "", "mean of ln(Q4/Q3) over the defect radii of table 1",
                 round(n8["raw_sem"], 3) if n8["raw_sem"] else None),
              _D("mean per-row step of the n_cl = 8 series, leak removed", round(n8["corrected_mean"], 3), "",
                 "same, after subtracting 1/Q_top(c) with the measured leak law", round(n8["corrected_sem"], 3) if n8["corrected_sem"] else None),
              _D("leak-removed mean of the n_cl = 8 series below the reference kappa", round(n8["corrected_mean_below_reference_pct"], 1), "%",
                 "100 (1 - mean / kappa)"),
              _D("scatter of the n_cl = 8 steps, raw over leak removed", round(n8["scatter_ratio_raw_over_corrected"], 2), "",
                 "ratio of the standard deviations")]
    # ---- second radius
    sr = ba.get("second_radius")
    if sr:
        reg["second_radius"] += [
            _D("kappa at r_d = 0.100a, absorber series, leak-removed", round(sr["kappa_absorber"], 4), "",
               "N_sep = 3, 4, 5 at clearances 9, 8, 8; no Fabry-Perot correction enters", 0.007,
               "steps %.4f and %.4f; zone-edge channel at the measured f is %.4f (%+.2f %%)"
               % (sr["series"]["absorber"]["steps"][0], sr["series"]["absorber"]["steps"][1], sr["channels_at_measured_f"]["zone_edge"],
                  sr["dev_zone_edge_pct"])),
            _D("kappa_pred, zone-edge channel at the second radius", round(sr["channels_at_measured_f"]["zone_edge"], 4), "",
               "the zone-edge channel of the kappa curve at the mean resonance frequency %.7f of the absorber series" % sr["f_mean"]),
            _D("deviation of kappa at r_d = 0.100a from the zone-edge channel", round(sr["dev_zone_edge_pct"], 2), "%", "100 (kappa / kappa_pred - 1)"),
            _D("kappa at r_d = 0.100a, PML series, leak-removed", round(sr["kappa_pml"], 4), "", "same records with PML termination", None,
               "lower than the absorber value by %.4f, which equals the drift of G with N" % (sr["kappa_absorber"] - sr["kappa_pml"])),
            _D("mean resonance frequency of the r_d = 0.100a series", round(sr["f_mean"], 7), "", "mean f_r of the three absorber records"),
            _D("G at r_d = 0.100a", {k: round(v, 4) for k, v in sr["G"].items()}, "", "Q_PML/Q_absorber after leak removal"),
            _D("resonance frequency pulled by the 25a cell, r_d = 0.100a", [round(sr["f_pull_pml_minus_absorber"][k], 8) for k in ("3", "4", "5")], "",
               "f_r(PML) - f_r(absorber) at N_sep = 3, 4, 5"),
        ]
        fp_ = sr["f_pull_pml_minus_absorber"]
        if all(fp_.get(k) for k in ("3", "4", "5")):
            reg["second_radius"].append(_D("per-row decay of the frequency pull of the 25a cell, r_d = 0.100a",
                                           [round(math.log(fp_["3"] / fp_["4"]), 2), round(math.log(fp_["4"] / fp_["5"]), 2)], "",
                                           "ln of the ratio of successive pulls, N_sep = 3 to 4 and 4 to 5", None,
                                           "close to the barrier coefficient: the pull scales with the rate of coupling to the guide"))
    # ---- field map
    rd = fe.get("row_decay") or []
    et = fe.get("eta") or []
    if rd:
        r0 = rd[0]
        if "cladding_k0" in r0:
            reg["field_map"] += [_D("2 Im k_y at k_x = 0, cladding rows %d to %d" % (r0["cladding_k0"]["rows"][0], r0["cladding_k0"]["rows"][-1]),
                                    round(abs(r0["cladding_k0"]["two_im_ky"]), 4), "",
                                    "complex fit of the Hann-windowed k_x = 0 projection of the stored DFT field, row by row", None,
                                    ("plane-wave prediction %.4f" % r0["channels"]["kappa_kx0"]) if r0.get("channels") else ""),
                                 _D("Re k_y at k_x = 0, cladding", round(abs(r0["cladding_k0"]["re_ky_over_pi"]), 4), "pi/a", "same fit")]
            chn = r0.get("channels") or {}
            if chn.get("kappa_kx0"):
                reg["field_map"] += [
                    _D("kappa_pred at k_x = 0 (cladding channel)", round(chn["kappa_kx0"], 4), "", "physical plane-wave root at the reference resonance"),
                    _D("kappa_pred at k_x = pi/a", round(chn["kappa_kxpi"], 4), "", "same, at the edge of the zone in k_x"),
                    _D("deviation of the cladding decay of the field map from the k_x = 0 channel",
                       round(100.0 * (abs(r0["cladding_k0"]["two_im_ky"]) / chn["kappa_kx0"] - 1.0), 2), "%", "100 (2 Im k_y / kappa_pred(k_x = 0) - 1)")]
        if "cladding_kbeta" in r0 and (r0.get("channels") or {}).get("kappa_edge"):
            ck = r0["cladding_kbeta"]
            reg["field_map"] += [_D("2 Im k_y at k_x = beta, cladding rows %d to %d" % (ck["rows"][0], ck["rows"][-1]), round(abs(ck["two_im_ky"]), 4), "",
                                    "complex fit of the k_x = beta projection of the stored DFT field, row by row", None,
                                    "Re k_y = %.4f pi/a; zone-edge channel %.4f" % (abs(ck["re_ky_over_pi"]), r0["channels"]["kappa_edge"])),
                                 _D("deviation of the cladding decay at k_x = beta from the zone-edge channel",
                                    round(100.0 * (abs(ck["two_im_ky"]) / r0["channels"]["kappa_edge"] - 1.0), 2), "%",
                                    "100 (2 Im k_y / kappa_pred(zone edge) - 1)")]
        reg["field_map"] += [_D("per-row cladding decay of the row energy, rows %d to %d" % (r0["cladding"]["rows"][0], r0["cladding"]["rows"][-1]),
                                [round(v, 3) for v in r0["cladding"]["decay_U"]], "", "energy in the band |y - j| < 0.5 of the stored field"),
                             _D("per-row barrier decay at k_x = beta, rows 3 to 2", round(r0["barrier"]["decay_kbeta"][0], 4), "",
                                "the one barrier step free of guide and cavity near field")]
        if len(rd) > 1 and "cladding_k0" in rd[1]:
            reg["field_map"].append(_D("2 Im k_y at k_x = 0, cladding, n_a = %.2f" % rd[1]["n"], round(abs(rd[1]["cladding_k0"]["two_im_ky"]), 4), "",
                                       "same analysis on the second field map"))
    if et:
        e0 = et[0]
        if e0.get("eta_hf") is not None:
            reg["field_map"] += [
                _D("eta_a, Hellmann-Feynman estimate consistent with the discretisation", round(e0["eta_hf"], 4), "",
                   "sum (1 - f) eps_a |E_z|^2 / sum <eps> |E_z|^2 with <eps> reconstructed on the E_z grid from the geometry"),
                _D("S_th from the Hellmann-Feynman estimate", round(e0["S_hf"], 1), "nm/RIU", "lambda_r eta_a / n_a", None,
                   ("measured S = %.2f nm/RIU" % e0["S_measured"]) if "S_measured" in e0 else ""),
                _D("deviation of S_th from the measured sensitivity, Hellmann-Feynman estimate", round(e0.get("dev_hf_pct", float("nan")), 2), "%",
                   "100 (S_th / S - 1)"),
                _D("eta_a, Hellmann-Feynman estimate with the permittivity at the points of the stored field", round(e0["eta_hf_centred"], 4), "",
                   "same, <eps> and 1 - f evaluated at the centred points instead of averaged from the E_z points", None,
                   "S_th = %.1f nm/RIU (%+.2f %%)" % (e0["S_hf_centred"], e0.get("dev_hf_centred_pct", float("nan")))),
                _D("eta_a, pure-analyte points only, consistent permittivity", round(e0["eta_hf_low"], 4), "",
                   "every mixed point counted as silicon", None, "S_th = %.1f nm/RIU (%+.2f %%)" % (e0["S_hf_low"], e0.get("dev_hf_low_pct", float("nan")))),
                _D("eta_a, mixed points counted as analyte, consistent permittivity", round(e0["eta_hf_high"], 4), "",
                   "pure-analyte points plus the whole energy of the mixed points", None, "S_th = %.1f nm/RIU" % e0["S_hf_high"]),
                _D("energy fraction at mixed points, consistent permittivity", round(e0["mixed_energy_fraction_hf"], 4), "",
                   "points whose E_z neighbours mix the two media"),
                _D("fraction of the points that are mixed", round(e0["mixed_point_fraction_hf"], 4), "", "by count"),
                _D("largest deviation of the reconstructed permittivity map from the stored one", float("%.2g" % e0["reconstruction_max_abs"]), "",
                   "Meep writes 12 / (sum of the 4-point sums of the diagonal of eps^-1) at the centred points; interior points", None,
                   "rms %.1e" % e0["reconstruction_rms"])]
        reg["field_map"] += [_D("eta_a, pure pixels over the stored-map denominator", round(e0["eta_low"], 4), "",
                                "the pure-pixel value of the earlier analysis: the stored map read as <eps>"),
                             _D("S_th from the pure pixels over the stored-map denominator", round(e0["S_low"], 1), "nm/RIU", "lambda_r eta / n_a", None,
                                ("deviation %+.2f %%" % e0["dev_low_pct"]) if "dev_low_pct" in e0 else ""),
                             _D("eta_a, linear unmixing of the stored map", round(e0["eta_lin"], 4), "", "the estimator registered as P8", None,
                                "S_th = %.1f nm/RIU" % e0["S_lin"]),
                             _D("eta_a, circle-pixel geometry over the stored-map denominator", round(e0["eta_geo"], 4), "", "earlier analysis", None,
                                "S_th = %.1f nm/RIU" % e0["S_geo"]),
                             _D("eta_a, threshold mask of the original analysis", round(e0["eta_mask"], 4), "", "stored eps below the mean of the two permittivities", None,
                                "S_th = %.1f nm/RIU" % e0["S_mask"]),
                             _D("energy fraction in mixed pixels of the stored map", round(e0["mixed_energy_fraction"], 4), "",
                                "boundary pixels with 0.001 < w < 0.999, over the stored-map denominator")]
        if len(et) > 1:
            e1 = et[1]
            if e1.get("S_hf") is not None:
                reg["field_map"] += [_D("S_th from the Hellmann-Feynman estimate, n_a = %.2f" % e1["n"], round(e1["S_hf"], 1), "nm/RIU", "same on the second map"),
                                     _D("eta_a, Hellmann-Feynman estimate, n_a = %.2f" % e1["n"], round(e1["eta_hf"], 4), "", "same on the second map")]
            reg["field_map"] += [_D("S_th from the pure pixels over the stored-map denominator, n_a = %.2f" % e1["n"], round(e1["S_low"], 1), "nm/RIU",
                                    "same on the second map"),
                                 _D("lambda_r of the n_a = %.2f field map" % e1["n"], round(e1["wavelength_nm"], 2), "nm", "stored in the map")]
        cbm = (analysis.get("combined") or {}).get("second_map_hf")
        if cbm:
            reg["field_map"].append(_D("deviation of S_th from the measured sensitivity at n_a = %.2f, Hellmann-Feynman estimate" % cbm["n"],
                                       round(cbm["deviation_pct"], 2), "%", "against %.1f nm/RIU extrapolated from the analyte sweep" % cbm["S_measured_extrapolated"]))
    # ---- numbers of the systematic runs
    if sw:
        a = sw["at_1p33"]
        reg["systematic_runs"] += [_D("S at n_a = 1.33, local (central difference)", round(a["S"], 2), "nm/RIU", "31-point analyte sweep, central differences"),
                            _D("S, single linear slope over the sweep", round(sw["chord"]["S"], 2), "nm/RIU", "chord fit; R^2 = %.6f" % sw["chord"]["R2"]),
                            _D("Q at n_a = 1.33 (25a cell)", round(a["Q"], 2), "", "reference record of the systematic runs"),
                            _D("FOM at n_a = 1.33 (25a cell)", round(a["FOM"], 1), "1/RIU", "S / (lambda_r / Q)"),
                            _D("detection limit at n_a = 1.33 (25a cell, R = fwhm/10)", float("%.3e" % a["DL10"]), "RIU", "(fwhm/10)/S"),
                            _D("Q drop over the sweep", round(sw["Q_drop_factor"], 2), "", "Q(1.300)/Q(1.450)"),
                            _D("resonance shift over the analyte sweep", round(sw["lambda_nm"][-1] - sw["lambda_nm"][0], 1), "nm",
                               "lambda_r(%.3f) - lambda_r(%.3f) at the raw lattice constant" % (sw["n"][-1], sw["n"][0]))]
        sc = sw.get("scaled_to_target") or {}
        if sw.get("a_nm") and a.get("lambda_nm"):
            k_ = float(sc.get("target_nm", 1550.0)) / a["lambda_nm"]           # the scaling a -> a lambda_t / lambda_r of section 2.2
            a_s = sw["a_nm"] * k_
            reg["systematic_runs"] += [_D("lattice constant scaled to 1550 nm", round(a_s, 1), "nm", "a lambda_t / lambda_r at n_a = 1.33"),
                                _D("rod radius scaled to 1550 nm", round(0.20 * a_s, 1), "nm", "0.20 a"),
                                _D("defect radius of the reference geometry scaled to 1550 nm", round(0.060 * a_s, 1), "nm", "0.060 a"),
                                _D("sensitivity scaled to 1550 nm", round(a["S"] * k_, 1), "nm/RIU", "S lambda_t / lambda_r"),
                                _D("linewidth of the 25a cell scaled to 1550 nm", round(a["fwhm_nm"] * k_, 3), "nm", "fwhm lambda_t / lambda_r")]
    ks = analysis.get("kappa_series")
    if ks:
        reg["systematic_runs"].append(_D("kappa, raw PML series of the systematic runs", round(ks["kappa"], 4), "", "log-linear fit to Q(N_sep) in the 25a cell, N_sep = %s" % ks["N"],
                                  round(ks["kappa_sd"], 4), "superseded by the reflectionless leak-free series"))
    bs = analysis.get("base_series")
    if bs and "quadratic_fit" in bs:
        qf = bs["quadratic_fit"]
        rr = qf.get("radius_range") or [None, None]
        reg["systematic_runs"].append(_D("f_r(r_d) quadratic law, n_cl = 8", [round(c, 6) for c in qf["coefficients_ascending"]], "",
                                  "polynomial fit to the %d admitted radii %g to %ga, ascending coefficients; R^2 = %.4f"
                                  % (qf["n"], rr[0], rr[1], qf["R2"])))
        if qf.get("standard_errors_ascending"):
            reg["systematic_runs"].append(_D("standard errors of the quadratic law", [float("%.2g" % c) for c in qf["standard_errors_ascending"]], "",
                                      "from the residual variance with n - 3 degrees of freedom"))
        reg["systematic_runs"] += [_D("R^2 of the quadratic law", round(qf["R2"], 4), "", "same fit"),
                            _D("rms residual of the quadratic law", float("%.2g" % qf["rms"]), "a/lambda", "same fit")]
        if bs.get("linear_fit"):
            reg["systematic_runs"] += [_D("R^2 of the linear law", round(bs["linear_fit"]["R2"], 4), "", "straight line through the same points"),
                                _D("rms residual of the linear law", float("%.2g" % bs["linear_fit"]["rms"]), "a/lambda", "same")]
        if bs.get("cubic_fit"):
            reg["systematic_runs"].append(_D("R^2 of the cubic law", round(bs["cubic_fit"]["R2"], 4), "", "cubic through the same points"))
    if bs and bs.get("fits_by_row"):
        fr = bs["fits_by_row"]
        reg["systematic_runs"] += [_D("quadratic intercepts of f_r(r_d) for N_sep = 2, 3, 4", [round(fr[k]["quadratic_coefficients_ascending"][0], 4) for k in sorted(fr, key=int)],
                               "a/lambda", "quadratic law of each row count, admitted radii"),
                            _D("quadratic curvatures of f_r(r_d) for N_sep = 2, 3, 4", [round(fr[k]["quadratic_coefficients_ascending"][2], 3) for k in sorted(fr, key=int)],
                               "", "same fits"),
                            _D("spread of the quadratic intercepts over N_sep = 2 to 4", round(bs["quadratic_intercept_spread_pct"], 2), "%",
                               "100 (max - min) / min")]
    rk = (bs or {}).get("raw_kappa")
    if rk:
        st23 = rk["steps"].get("2->3")
        if st23:
            reg["systematic_runs"].append(_D("mean per-row step ln(Q3/Q2) of the n_cl = 8 series", round(st23["mean"], 3), "",
                                      "mean over the %d admitted radii" % len(rk["common_rd"]), round(st23["sem"], 3) if st23.get("sem") else None))
        reg["systematic_runs"].append(_D("pooled kappa of the n_cl = 8 series", round(rk["kappa"], 3), "",
                                  "log-linear fit of ln Q against N_sep = 2, 3, 4 over the %d admitted radii" % len(rk["common_rd"])))
    cv = analysis.get("convergence") or {}
    r2032 = cv.get("resolution_20_vs_32_at_ncl8") or []
    if r2032:
        reg["systematic_runs"] += [_D("change of Q from resolution 20 to 32, N_sep = 4, n_cl = 8, r_d = " + ", ".join("%.3f" % e["radius"] for e in r2032) + "a",
                               [round(e["change_pct"], 1) for e in r2032], "%", "100 (Q32 / Q20 - 1), table 2"),
                            _D("largest change of Q from resolution 20 to 32", round(max(abs(e["change_pct"]) for e in r2032), 1), "%",
                               "largest magnitude over the same radii")]
    fr = cv.get("resolution_20_vs_32_frequency") or []
    if fr:                                                      # SPRAT 1.2.1: the resonance frequency of the same pairs
        reg["systematic_runs"] += [_D("change of f_r from resolution 20 to 32, N_sep = 4, n_cl = 8, r_d = "
                                      + ", ".join("%.3f" % e["radius"] for e in fr) + "a",
                                      [round(e["f_change_pct"], 4) for e in fr], "%", "100 (f32 / f20 - 1), the pairs of table 2"),
                                   _D("largest change of f_r from resolution 20 to 32", round(max(abs(e["f_change_pct"]) for e in fr), 3), "%",
                                      "largest magnitude over the same radii"),
                                   _D("largest shift of lambda_r from resolution 20 to 32", round(max(abs(e["dlambda_nm"]) for e in fr), 2), "nm",
                                      "raw lattice constant, same radii")]
    c1012 = cv.get("cladding_10_vs_12_at_res24") or []
    if c1012:
        reg["systematic_runs"].append(_D("change of Q from n_cl = 10 to 12 at resolution 24, r_d = " + ", ".join("%.3f" % e["radius"] for e in c1012) + "a",
                                  [round(abs(e["change_pct"]), 1) for e in c1012], "%", "100 |Q10 / Q12 - 1|"))
    ps = analysis.get("preliminary_sweep") or {}
    wm, cm, cg = ps.get("waveguide_modes") or {}, ps.get("cavity_modes") or {}, ps.get("cavity_groups") or {}
    if wm:
        reg["systematic_runs"] += [_D("Q of the admitted records holding a mode of the waveguide: smallest, median, largest",
                               [round(wm["Q_min"], 1), round(wm["Q_median"], 1), round(wm["Q_max"], 1)], "",
                               "preliminary sweep, n_cl = 6, resolution 20, N_sep = 1 or r_d >= 0.13a (records.waveguide_mode)"),
                            _D("median wavelength of the admitted records holding a mode of the waveguide", round(wm["lambda_median_nm"]), "nm",
                               "same records, raw lattice constant")]
    if cm:
        reg["systematic_runs"] += [_D("admitted cavity records of the preliminary sweep", cm["admitted"], "", "admitted rows of the sweep less the waveguide modes"),
                            _D("Q of the admitted cavity records of the preliminary sweep: smallest, largest", [round(cm["Q_min"]), round(cm["Q_max"])], "",
                               "same records"),
                            _D("wavelength range of the admitted cavity records of the preliminary sweep", [round(x) for x in cm["lambda_range_nm"]], "nm",
                               "same records, raw lattice constant")]
    if cg:
        reg["systematic_runs"].append(_D("median ratio Q(dx = 0.08a)/Q(0) over the cavity groups of the preliminary sweep", round(cg["median_ratio"], 3), "",
                                  "%d (r_d, dy, N_sep) groups with N_sep >= 2 and r_d < 0.13a; interquartile range %.3f to %.3f"
                                  % (cg["groups"], cg["quartiles"][0], cg["quartiles"][1])))
    ab = ba.get("absorber_sweep")
    if ab and 0.085 in [round(x, 4) for x in ab["rd"]]:
        i = [round(x, 4) for x in ab["rd"]].index(0.085)
        if ab["margin"][i] and ab["margin"][i] >= 4:
            reg["second_radius"].append(_D("Q(absorber, r_d = 0.085a) repeated at margin 4.5", round(ab["Q"][i], 2), "",
                                           "the longest-signal absorber record at r_d = 0.085a (%s)" % ab["labels"][i]))
    _manuscript_rows(reg, analysis)
    if analysis.get("termination"):
        reg["termination"] = _termination_rows(analysis["termination"])
    rl = _reflectionless_rows(analysis)
    if rl:
        reg["reflectionless"] = rl
    if counts:
        reg["dataset"] += [_D("harmonic-inversion records in the deposit", counts["harminv"], "", "rows of records.csv with mode harminv"),
                           _D("records from the systematic runs", counts.get("campaign", counts["harminv"] - counts.get("verification", 0)), "",
                              "rows whose source is not a verification run"),
                           _D("records from the verification runs", counts.get("verification", 0), "", "rows with tag verification or source legacy-verification"),
                           _D("records admitted by the criterion of section 2.4", counts["valid"], "", "margin >= 1 and gap position 0.2-0.8"),
                           _D("records failing the margin condition", counts["excluded_margin"], "", "Q > Q_lim"),
                           _D("records outside the gap window that clear the margin", counts["excluded_gap"], "", ""),
                           _D("admitted records in the dead zone", counts.get("dead_zone_admitted", 0), "",
                              "admitted rows with r_d >= 0.14a (preliminary sweep), where the strongest mode is a mode of the "
                              "finite waveguide at every radius"),
                           _D("admitted records holding a mode of the waveguide", counts.get("waveguide_mode_admitted", 0), "",
                              "admitted rows of the preliminary sweep (n_cl = 6, resolution 20) at N_sep = 1 or r_d >= 0.13a; "
                              "column waveguide_mode of records.csv"),
                           _D("admitted records holding a cavity mode", counts["valid"] - counts.get("waveguide_mode_admitted", 0), "",
                              "admitted rows less those holding a mode of the waveguide"),
                           _D("share of the harmonic-inversion records admitted and holding a cavity mode",
                              round(100.0 * (counts["valid"] - counts.get("waveguide_mode_admitted", 0)) / counts["harminv"], 2), "%",
                              "admitted records holding a cavity mode over the harmonic-inversion records"),
                           _D("records clearing the margin alone", counts["resolved"], "", "Q <= Q_lim"),
                           _D("transmission spectra", counts["spectra"], "", "seven per guide termination when the absorber spectra are present"),
                           _D("transmission spectra with the guide continued into the absorber", counts.get("spectra_absorber", 0), "",
                              "spectrum records with termination absorber (raw_h14b)"),
                           _D("discrete-Fourier-transform field maps", counts["fields"], "", ""),
                           _D("band-structure computations", counts["bands"], "", "bands records, including the lattice-constant calibration"),
                           _D("excluded records lying in the reference geometry", counts.get("excluded_reference", 0), "",
                              "excluded rows with r_d = 0.060a, N_sep = 4, n_cl = 12, res 24, PML, 25a guide")]
    return reg


TABLE5_SPAN = "n_a = 1.305, 1.330, 1.350, 1.375, 1.400, 1.425, 1.445"


def _reflectionless_rows(analysis: dict[str, Any]) -> list[dict[str, Any]]:
    """SPRAT 1.2.1: the reflectionless quality factor and figure of merit at the rows of table 5
    (``combined.physical_values``), and the spectral value interpolated to n_a = 1.33 against the harmonic-inversion one."""
    R: list[dict[str, Any]] = []
    pv = ((analysis.get("combined") or {}).get("physical_values") or {})
    tab = pv.get("table") or []
    if len(tab) == 7:
        rule = ("harmonic inversion with the absorber and the leak removed at n_a = 1.330 (Q_w(4) of the series); Q(harmonic inversion, "
                "25a cell) / G of the absorber spectra at 1.350 to 1.425; linear interpolation in n_a of that Q_w at 1.305 and 1.445")
        R += [_D("reflectionless Q_w at the rows of table 5, " + TABLE5_SPAN, [round(r["Q_w"]) for r in tab], "", rule),
              _D("reflectionless FOM_w at the rows of table 5, " + TABLE5_SPAN, [round(r["FOM_w"]) for r in tab], "1/RIU",
                 "S Q_w / lambda_r with the local sensitivity and resonance of the analyte sweep"),
              _D("sources of Q_w at the rows of table 5", [r["source"] for r in tab], "", "see the reflectionless Q_w row")]
    a133 = pv.get("at_1p33") or {}
    if a133:
        R += [_D("spectral Q_w interpolated to n_a = 1.33", round(a133["Q_w_spectra_interpolated"]), "",
                 "Q_w of the absorber spectra at 1.325 and 1.350 interpolated linearly in n_a"),
              _D("deviation of the spectral Q_w interpolated to n_a = 1.33 from the harmonic-inversion value",
                 round(a133["deviation_pct"], 1), "%", "against the harmonic-inversion Q_w(4) with the leak removed")]
    return R


def _termination_rows(t: dict[str, Any]) -> list[dict[str, Any]]:
    """The rows of the comparison of the two guide terminations (analysis key ``termination``)."""
    rows = t["rows"]
    if not rows or any(r.get("status") != "graded" for r in rows):
        return [_D("comparison of the two terminations", "incomplete: " + "; ".join("%.3f %s" % (r["n"], r.get("status")) for r in rows
                                                                                       if r.get("status") != "graded"), "", "")]
    a1 = set(t.get("a1_graded") or [])
    a2 = (t.get("A2") or {}).get("rows") or {}
    span = "n_a = 1.300 to 1.450 in steps of 0.025"
    fit = "Fano fit over the whole window (the routine of 48_fano_karsilastir_S4b.py), each spectrum over the reference of its own termination"
    su = t.get("summary") or {}
    R = [_D("Fano q in the 25a cell, " + span, [round(r["pml"]["q"], 4) for r in rows], "", fit, None, "the spectra of section 3.5"),
         _D("Fano q with the guide continued into the absorber, " + span, [round(r["absorber"]["q"], 4) for r in rows], "", fit, None,
            "every other run parameter as the 25a-cell spectrum of the same index"),
         _D("factor |q| in the 25a cell over |q| with the absorber, " + span, [round(r["factor"], 2) for r in rows], "", "ratio of the two fits"),
         _D("grade by the criterion fixed before the runs (D-22), " + span, [r["grade"] for r in rows], "",
            ">= 3 supports, <= 1.5 refutes, otherwise inconclusive"),
         _D("grade by the amended criterion (A1), " + span, [r["grade"] if r["n"] in a1 else "not graded" for r in rows], "",
            "the ratio graded where |q| in the 25a cell exceeds 0.2"),
         _D("verdict by the criterion fixed before the runs (D-22)", t["D22"]["verdict"], "", "pass = supports at all seven indices"),
         _D("verdict by the amended criterion (A1)", t["A1"]["verdict"], "", "pass = supports at the five graded indices"),
         _D("A1 threshold of |q| in the 25a cell for grading the ratio", 0.2, "", "fixed 2026-09-24 15:46:39 +03:00 before the first absorber cavity spectrum"),
         _D("Q from the Fano fit with the absorber, " + span, [round(r["absorber"]["Q"]) for r in rows], "", "lambda_r / fwhm of the fit"),
         _D("Q from the Fano fit in the 25a cell, " + span, [round(r["pml"]["Q"]) for r in rows], "", "lambda_r / fwhm of the fit"),
         _D("t_end/tau of the absorber spectra, " + span, [round(r["absorber"]["t_end_over_tau"], 2) for r in rows], "",
            "tau = Q / (pi f_r); gate G5 requires 3"),
         _D("largest transmission near resonance with the absorber, " + span, [round(r["absorber"]["T_max_near"], 3) for r in rows], "",
            "within five linewidths of the fitted resonance")]
    if su:
        R += [_D("largest |q| with the absorber", round(su["abs_q_absorber_max"], 3), "", "over the seven spectra"),
              _D("factor at the five graded indices: smallest, largest", [round(su["factor_inner_min"], 2), round(su["factor_inner_max"], 2)], "",
                 "n_a = 1.325 to 1.425"),
              _D("largest transmission near resonance with the absorber, over the seven spectra", round(su["T_max_near_absorber_max"], 3), "",
                 "within five linewidths of the fitted resonance"),
              _D("largest fit rms: 25a cell, absorber", [round(max(r["pml"]["rms"] for r in rows), 3), round(su["rms_absorber_max"], 3)], "",
                 "gate G4 requires 0.05"),
              _D("smallest distance of an absorber resonance from the window edge", round(su["edge_margin_absorber_min"], 1), "linewidths",
                 "gate G3 requires 5")]
    if t.get("S1_value") is not None:
        R += [_D("S1: absorber Q interpolated to n_a = 1.33", round(t["S1_value"]), "", "linear interpolation between n_a = 1.325 and 1.350"),
              _D("S1: deviation of the interpolated absorber Q from the harmonic-inversion value", round(100.0 * (t["S1_value"] / 4495.1 - 1.0), 1),
                 "%", "against 4495.1, the absorber cell at n_a = 1.33 (record v4_H4b_s4); tolerance 5 %")]
    if a2:
        ns = [r["n"] for r in rows]
        G = [a2[n]["G"] for n in ns]
        Qw = [a2[n]["Q_w"] for n in ns]
        R += [_D("cell factor G = Q(25a cell) / Q(absorber) from the spectra, " + span, [round(x, 3) for x in G], "",
                 "ratio of the two Fano quality factors (amendment A2)"),
              _D("reflectionless Q_w = Q(harmonic inversion, 25a cell) / G, " + span, [round(x) for x in Qw], "",
                 "harmonic-inversion Q of the 25a cell (8327 to 1662) over the cell factor"),
              _D("Q of the 25a cell by harmonic inversion, " + span, [round(a2[n]["Q_pml_harminv"]) for n in ns], "",
                 "the harmonic-inversion Q of the analyte sweep at that index, rounded (the q_est of each spectrum run)"),
              _D("cell factor G at n_a = 1.300 and 1.450", [round(G[0], 2), round(G[-1], 2)], "", "from the spectra"),
              _D("reflectionless Q_w at n_a = 1.300 and 1.450", [round(Qw[0], -1), round(Qw[-1], -1)], "", "from the spectra")]
        if t["A2"].get("Q_w_fall"):
            R += [_D("fall of the reflectionless Q_w over the analyte sweep", round(t["A2"]["Q_w_fall"], 2), "", "Q_w(1.300) / Q_w(1.450)"),
                  _D("opening of the reflectionless linewidth over the analyte sweep", round(t["A2"]["fwhm_w_opening"], 2), "",
                     "(lambda_r / Q_w) at 1.450 over that at 1.300"),
                  _D("opening of the linewidth of the 25a cell over the analyte sweep", round(t["A2"]["fwhm_cell_opening"], 2), "",
                     "(lambda_r / Q) at 1.450 over that at 1.300, harmonic-inversion Q, lambda_r of the fit")]
        if t["A2"].get("G_1p33"):
            R.append(_D("cell factor G interpolated to n_a = 1.33 from the spectra", round(t["A2"]["G_1p33"], 3), "",
                        "linear interpolation between 1.325 and 1.350"))
        if t["A2"].get("G_1p33_harminv"):
            R.append(_D("cell factor at n_a = 1.33 from harmonic inversion, 25a cell over absorber", round(t["A2"]["G_1p33_harminv"], 4), "",
                        "%.2f / %.2f, the longest-signal records of the two terminations; 1.5423 with the cladding leak removed"
                        % (t["A2"]["Q_pml_harminv_1p33"], t["A2"]["Q_absorber_harminv_1p33"])))
    return R


def _g(d, *path):
    for p in path:
        if not isinstance(d, dict) or p not in d:
            return None
        d = d[p]
    return d


def _manuscript_rows(reg: dict, A: dict) -> None:
    """The rows the manuscript registry added in v6 and v7 (sections vertical, barrier_row_v7, absorption_v7 and
    derived_v7 there; vertical, barrier_row, absorption and derived here) and the superseded statements, with the
    quantity names of the manuscript registry so that ``sprat audit`` can compare them row by row."""
    ba = A.get("barrier") or {}
    pf, cb, rp = ba.get("performance") or {}, A.get("combined") or {}, (A.get("repeats") or {}).get("points") or {}
    v = pf.get("vertical")
    if v:
        q0 = v["Q_tot_no_vertical"]
        red, rp_ = v["Q_perp_for_reduction"], v["reduction_pct_at_Q_perp"]
        reg["vertical"] = [
            _D("Q_tot with no vertical channel", round(q0, 1), "", "1/Q_tot = 1/Q_w + 1/Q_abs at the reference point"),
            _D("Q_perp leaving Q_tot within 3 per cent", round(red["3"], 1), "",
               "1/Q_tot = 1/Q_w + 1/Q_perp + 1/Q_abs solved for a 3 per cent reduction: Q_perp = Q_0 (1 - r) / r"),
            _D("Q_perp costing 10 per cent of Q_tot", round(red["10"], 1), "", "same, for a 10 per cent reduction"),
            _D("Q_perp costing 50 per cent of Q_tot", round(red["50"], 1), "", "same, for a 50 per cent reduction"),
            _D("reduction of Q_tot at Q_perp = 1e4", round(rp_["1e4"], 3), "per cent",
               "Q_0 / (Q_0 + Q_perp): Q_tot falls from %.0f to %.0f" % (q0, v["Q_tot_at_Q_perp"]["1e4"])),
            _D("reduction of Q_tot at Q_perp = 1e5", round(rp_["1e5"], 3), "per cent",
               "Q_0 / (Q_0 + Q_perp): Q_tot falls from %.0f to %.0f" % (q0, v["Q_tot_at_Q_perp"]["1e5"]))]
    knob, brs = ba.get("barrier_row_knob") or {}, ba.get("barrier_row_series")
    if knob and brs:
        reg["barrier_row"] = []
        for rb in sorted(knob, key=float):
            k = knob[rb]
            reg["barrier_row"].append(dict(quantity=f"Q, N_sep = 4, PML 25a cell, second barrier row radius {rb}a", value=round(k["Q"], 2), unit="",
                                           kind="measured", label=k["label"], host=k["host"], job_id=k["job_id"], code_sha256=k["code_sha256"],
                                           f_r=round(k["f"], 7), margin=(round(k["margin"], 2) if k["margin"] else None),
                                           lambda_r_nm_raw=(round(k["lambda_nm_stored"], 3) if k.get("lambda_nm_stored") else None), note=""))
        reg["barrier_row"].append(_D("resonance shift over the barrier-row series 0.18a to 0.22a", round(brs["lambda_shift_nm"], 4), "nm",
                                     "a/f_r at the raw lattice constant a = %.1f nm, %ga minus %ga" % (brs["a_nm"], brs["rb"][-1], brs["rb"][0]), None,
                                     "quoted as %.1f nm in section 3.4" % brs["lambda_shift_nm"]))
    ab = pf.get("absorption")
    if ab:
        reg["absorption"] = [_D("Q_abs quoted in section 4", round(pf["Q_abs"], 1), "",
                                "lambda_r / (2 k S) = 2 pi n_a / (eta_a alpha lambda) with eta_a = %.4f (required by the measured sensitivity), "
                                "alpha = %.2f 1/cm, lambda = %.0f nm" % (pf["eta_used"], ab["alpha_per_cm"], ab["target_nm"]), None, "quoted as 9.2e3")]
        for key, q in (("Q_abs_eta_hf", "Q_abs at the Hellmann-Feynman eta_a of the field map"),
                       ("Q_abs_eta_low", "Q_abs at the stored-map pure-pixel eta_a (the v8 value)")):
            if ab.get(key):
                reg["absorption"].append(_D(q, round(ab[key], 1), "", "n_a / (2 k eta_a)"))
        if ab.get("Q_abs_mixed_as_silicon") and ab.get("Q_abs_mixed_as_analyte"):   # SPRAT 1.2.1
            reg["absorption"].append(_D("Q_abs with the mixed points counted wholly as analyte and wholly as silicon",
                                        [float("%.2g" % ab["Q_abs_mixed_as_analyte"]), float("%.2g" % ab["Q_abs_mixed_as_silicon"])], "",
                                        "n_a / (2 k eta_a) with the Hellmann-Feynman eta_a of the field map, mixed points as analyte "
                                        "(0.749) or as silicon (0.518)"))
    D = reg["derived"]
    # two statements of the manuscript corrected in v8 against the records (section corrected_v8 of its registry)
    tf_ = cb.get("law_transfer_to_fine_sweep") or {}
    if tf_.get("max_abs_deviation") is not None:
        D.append(_D("largest deviation of the quadratic law from the fine sweep in the reference geometry",
                    float("%.2g" % abs(tf_["max_abs_deviation"])), "a/lambda",
                    "law fitted to the admitted radii of the n_cl = 8, resolution 20 sweep, evaluated at the %d admitted radii of the fine sweep; "
                    "largest at r_d = %.3fa" % (tf_.get("points", 0), tf_.get("at_radius", float("nan")))))
        rs_ = A.get("radius_sweep") or {}
        if rs_.get("rd") and rs_.get("f"):
            f_at = rs_["f"][[round(x, 4) for x in rs_["rd"]].index(round(tf_["at_radius"], 4))]
            D.append(_D("largest deviation of the quadratic law from the fine sweep, in wavelength at 1550 nm",
                        round(1550.0 * abs(tf_["max_abs_deviation"]) / f_at, 1), "nm", "lambda delta f / f_r at the radius of the largest deviation"))
    dl_ = cb.get("dlambda_drd_at_reference")
    if dl_:
        D += [_D("d lambda_r / d r_d at the reference radius, five-point stencil of the fine sweep", round(dl_["stencil_nm_per_nm"], 2), "nm/nm",
                 "step 0.005a, scale invariant"),
              _D("d lambda_r / d r_d at the reference radius from the quadratic law", round(dl_["law_nm_per_nm"], 2), "nm/nm", dl_["note"]),
              _D("agreement of the stencil with the quadratic law", round(dl_["agreement_pct"], 2), "%", "100 (stencil / law - 1)")]
    for dy_key, lab in ((0.0, "dy = 0"), (0.1, "dy = 0.10a")):
        qd = _g(A, "displacement", "dx", dy_key, "quadratic_fall") or _g(A, "displacement", "dx", str(dy_key), "quadratic_fall")
        if qd:
            D += [_D(f"fall of Q at dx = 0.08a, {lab}", round(-qd["rel_change_pct"][-1], 2), "%", "100 (1 - Q(0.08a)/Q(0))"),
                  _D(f"quadratic coefficient of Q(dx)/Q(0) - 1, {lab}", round(qd["c_per_a2"], 2), "1/a^2",
                     "least squares through the origin over dx = 0.02a to 0.08a", None,
                     "largest residual %.3f %%; with a free offset %.3f %% the coefficient is %.2f" % (qd["residual_max_pct"], qd["with_offset"]["offset_pct"],
                                                                                                        qd["with_offset"]["c_per_a2"])),
                  _D(f"largest residual of the quadratic law in dx, {lab}", round(qd["residual_max_pct"], 3), "%", "same fit"),
                  _D(f"offset of the quadratic law in dx with a free constant, {lab}", round(qd["with_offset"]["offset_pct"], 3), "%",
                     "step between the run at dx = 0 (mirror symmetry exploited) and the displaced runs (none)")]
    sp_ = A.get("spectra_summary") or {}
    if sp_.get("frequency_step"):
        D += [_D("frequency step of the transmission spectra", float("%.2g" % sp_["frequency_step"]), "a/lambda", "median spacing of the 601 points"),
              _D("transmission at the samples next to the minimum", [float("%.2g" % sp_["T_min_neighbour_min"]), float("%.2g" % sp_["T_min_neighbour_max"])], "",
                 "smallest and largest over the seven spectra"),
              _D("smallest transmission sample over the seven spectra", [float("%.3g" % sp_["T_min_min"]), float("%.3g" % sp_["T_min_max"])], "",
                 "the range of the smallest sample of each spectrum"),
              _D("normalised transmission far from resonance", [round(sp_["T_far_min"], 3), round(sp_["T_far_max"], 3)], "",
                 "median of T more than fifteen linewidths from resonance; %d of the %d spectra extend that far"
                 % (sp_.get("spectra_with_far_window", 0), sp_.get("count", 0))),
              _D("largest normalised transmission near resonance", round(sp_["T_max_near_max"], 2), "", "within fifteen linewidths, over the seven spectra")]
    ss_ = A.get("spectra_summary") or {}
    if ss_.get("fwhm_ratio_input_min") is not None and ss_.get("fwhm_ratio_input_max") is not None:
        for q_, k_ in (("smallest linewidth narrowing by the input-plane normalisation", "fwhm_ratio_input_min"),
                       ("largest linewidth narrowing by the input-plane normalisation", "fwhm_ratio_input_max")):
            D.append(_D(q_, round(ss_[k_], 2), "", "harmonic-inversion linewidth over the Fano linewidth of the spectrum "
                                                   "normalised to the input plane of the same run, over the seven spectra"))
    dr = cb.get("displacement_over_row_coefficient")
    if dr:
        D.append(_D("ratio of the transverse-displacement coefficient to the row coefficient", round(dr["ratio"], 4), "",
                    "%.4f / %.4f (kappa_dy of the linear fit of ln Q against dy, over the four-point kappa)" % (dr["kappa_dy"], dr["kappa"])))
    dsp = A.get("displacement") or {}
    ch_ = [abs(r - 1.0) for e in list((dsp.get("dy") or {}).values()) + list((dsp.get("dx") or {}).values())
           for r in (e.get("ratio") or [])]
    if ch_:
        D.append(_D("largest change of Q over the displacement series of table 4", round(100.0 * max(ch_), 1), "%",
                    "largest |Q/Q(0) - 1| over the transverse series (r_d = 0.050a, 0.060a) and the axial series (dy = 0, 0.10a)",
                    None, "stated as the bound 'less than 7 per cent' in the abstract, the introduction and the conclusions"))
    for key, q in (("reference_N5", "reproduction of the N_sep = 5 point by repeats at higher margin"),
                   ("reference", "reproduction of the reference point by a repeat at higher margin")):
        g = rp.get(key)
        if g:
            D.append(_D(q, round(g["max_deviation_pct"], 5), "%", "systematic run %.2f (margin %.2f) against " % (g["base"]["Q"], g["base"]["margin"])
                        + ", ".join("%.2f (margin %.2f)" % (x["Q"], x["margin"]) for x in g["repeats"]), None,
                        "largest |Q_repeat / Q_base - 1|"))
    if brs:
        D.append(_D("change of Q over the barrier-row series", round(brs["Q_change_pct"], 2), "%",
                    "Q(%ga) / Q(%ga) - 1" % (brs["rb"][-1], brs["rb"][0]), None,
                    "ln ratio %.4f, i.e. %.1f %% of the row step %.4f" % (brs["ln_Q_ratio"], brs["row_step_fraction_pct"], ba["series"]["kappa"])))
    ce = cb.get("cell_factor_estimates")
    if ce and ce.get("barrier_row"):
        D.append(_D("cell-factor change over the barrier-row series (estimate)", round(ce["barrier_row"]["change_pct"], 2), "%",
                    "d ln G / d f = %.1f per a / (df/dr_d = %.4f) = %.1f; delta f = %.3g" % (ce["dlnG_drd_per_a"], ce["df_drd_per_a"], ce["dlnG_df"],
                                                                                          ce["barrier_row"]["f_change"]), None, "estimate"))
    if ce and ce.get("dx_max"):
        dm = ce["dx_max"]
        D.append(_D("cell-factor change at dx = 0.08a (estimate), dy = 0 and 0.10a", [round(dm[k]["change_pct"], 2) for k in sorted(dm)], "%",
                    "resonance change at the largest dx times d ln G/df = %.1f" % ce["dlnG_df"], None, "estimate"))
        D.append(_D("intrinsic change of Q at dx = 0.08a after the cell-factor estimate, dy = 0 and 0.10a",
                    [round(dm[k]["intrinsic_pct"], 2) for k in sorted(dm)], "%", "measured ln change of the 25a cell less the cell-factor estimate",
                    None, "estimate"))
    if ce and ce.get("dy_max"):
        D.append(_D("cell-factor change at |dy| = 0.10a (estimate)", round(abs(ce["dy_max"]["change_pct"]), 2), "%",
                    "mean resonance change %.3g at |dy| = %.2fa times d ln G/df = %.1f" % (ce["dy_max"]["f_change"], ce["dy_max"]["dy"], ce["dlnG_df"]),
                    None, "the cell factor falls (%.2f %%), equally for both signs" % ce["dy_max"]["change_pct"]))
    if pf:
        D.append(_D("absorption cap quoted in the abstract", round(pf["Q_abs"], 1), "", "Q_abs = lambda_r / (2 k S), quoted as 'near 9200'"))
        D.append(_D("coupling optimum", round(pf["Q_w_opt"], 1), "", "(1 + sqrt 3)/2 x Q_abs"))
    lw = _g(ba, "tolerance", "linewidths_per_nm_radius_error")
    if lw and "absorber" in lw:
        D.append(_D("reflectionless linewidth count for a 1 nm radius error", round(lw["absorber"], 3), "",
                    "%.4f nm / (%.0f nm / %.1f)" % (lw["dlambda_drd_nm_per_nm"], lw["target_nm"], ba["series"]["Q_w"]["4"])))
    reg["superseded"] = superseded(A)


def superseded(A: dict) -> list[dict[str, str]]:
    """The superseded statements of the manuscript registry, with every number recomputed from the records."""
    ba, ks, sw, fe = A.get("barrier") or {}, A.get("kappa_series") or {}, A.get("analyte_sweep") or {}, A.get("fields") or {}
    rp, rs, cb = A.get("resolved_pair") or {}, A.get("radius_sweep") or {}, A.get("combined") or {}
    s, mc, pf, l = ba.get("series") or {}, ba.get("monte_carlo") or {}, ba.get("performance") or {}, ba.get("leak") or {}
    sr, tol, e = ba.get("second_radius") or {}, ba.get("tolerance") or {}, fe.get("reference") or {}
    out = []
    try:
        qref = (ba.get("reference_campaign") or {}).get("Q")
        out.append(dict(v3_value="kappa = %.4f +- %.4f" % (ks["kappa"], ks["kappa_sd"]), v4_value="kappa = %.4f +- %.4f" % (s["kappa"], mc["kappa_sd"]),
                        reason="log-linear fit to raw PML values in a 25a cell with mixed clearances and a non-asymptotic N_sep = 2 point"))
        out.append(dict(v3_value="per-row factor %.1f" % ks["per_row_factor"], v4_value="%.2f" % s["per_row_factor"], reason="follows from kappa"))
        out.append(dict(v3_value="Q = %.0f at the reference point" % qref, v4_value="Q_w = %.0f reflectionless leak-free" % s["Q_w"]["4"],
                        reason="%.0f is the 25a cell value and contains G = %.4f" % (qref, ba["G"]["4"])))
        out.append(dict(v3_value="FOM = %.0f 1/RIU" % sw["at_1p33"]["FOM"],
                        v4_value="%d lossless 2D, %d in water" % (round(pf["FOM_lossless_N4"]), round(pf["FOM_water_N4"])),
                        reason="cell-inflated Q, and the water-corrected value was never stated"))
        dev2 = (cb.get("second_map_bound") or {}).get("deviation_pct")
        out.append(dict(v3_value="eta_a = %.3f, a %.1f %% gap" % (e["eta_mask"], e["dev_mask_pct"]),
                        v4_value="eta_a >= %.3f, S_th = %.1f nm/RIU against a measured %.1f (%+.2f %%)" % (e["eta_low"], e["S_low"], e["S_measured"], e["dev_low_pct"]),
                        reason="the masked estimator counts whole boundary pixels as analyte; version 4 replaced it by the pure-pixel value, then "
                               "read as a lower bound, which at n_a = 1.45 lay %.1f %% above the measured sensitivity (version 9 corrects the "
                               "reading, see below)" % dev2))
        out.append(dict(v3_value="clearance rule n_cl - N_sep >= 6", v4_value="c >= 8 + ln(100 Q_w/Q_top(8))/kappa_top",
                        reason="clearance 7 carries %.2f %% of 1/Q at N_sep = 5 and clearance 8 carries %.2f %% at N_sep = 6"
                               % (l["leak_share_pct"]["N5_c7"], l["leak_share_pct"]["N6_c8"])))
        out.append(dict(v3_value="dQ/dr_d = %+.0f 1/nm, Q rises on the flank of the peak at %.4fa" % (rs["local_Q"]["dQ_drd_per_nm"], rs["local_Q"]["extremum_rd"]),
                        v4_value="dQ/dr_d = %.0f 1/nm, Q falls smoothly" % tol["absorber"]["dQ_drd_per_nm"],
                        reason="the peak and the sign are Fabry-Perot features of the 25a cell"))
        do = ba["detrended_oscillation"]
        out.append(dict(v3_value="a linear coefficient is not meaningful for Q (section 3.8)",
                        v4_value="reflectionless Q(r_d) is smooth; a linear coefficient is meaningful",
                        reason="the oscillation that motivated the disclaimer is the cell, suppressed %.1fx by the absorber"
                               % (do["pml_25a"]["detrended_ptp"] / do["absorber"]["detrended_ptp"])))
        out.append(dict(v3_value="kappa = %.2f at r_d = 0.100a" % rp["kappa"],
                        v4_value="%.3f +- %.3f, reflectionless and leak-free" % (sr["kappa_absorber"], sr["kappa_absorber_uncertainty"]["total"]),
                        reason="the v3 value is a raw PML value in a 25a cell; job A measured it on an absorber series that needs no Fabry-Perot "
                               "correction, and it matches the zone-edge channel to %.1f %%" % sr["dev_zone_edge_pct"]))
    except (KeyError, TypeError):
        return out
    out += superseded_v9(A)
    return out


def superseded_v9(A: dict) -> list[dict[str, str]]:
    """The statements of manuscript version 8 that version 9 corrects, both sides recomputed from the records: the
    original plane-wave layer against the physical roots, every record against the admitted ones, the stored-map
    energy fraction against the consistent one, and the absorption basis."""
    ba, fe, cb = A.get("barrier") or {}, A.get("fields") or {}, A.get("combined") or {}
    lg, bs, bsa = A.get("pwe_legacy_channels") or {}, A.get("base_series") or {}, A.get("base_series_all_points") or {}
    out: list[dict[str, str]] = []

    def add(fn):
        try:
            out.append(fn())
        except (KeyError, TypeError, ValueError, IndexError, ZeroDivisionError):
            pass

    s, ch, sr = ba.get("series") or {}, ba.get("channels") or {}, ba.get("second_radius") or {}
    add(lambda: dict(
        v3_value="zone-edge root %.4f and %.4f, least-evanescent complex pair %.4f and %.4f; measured %+.2f and %+.2f %% from the first, %+.2f and %+.2f %% from the second"
                 % (lg["reference"]["zone_edge"], lg["second"]["zone_edge"], lg["reference"]["slowest"], lg["second"]["slowest"],
                    100 * (s["kappa"] / lg["reference"]["zone_edge"] - 1), 100 * (sr["kappa_absorber"] / lg["second"]["zone_edge"] - 1),
                    100 * (s["kappa"] / lg["reference"]["slowest"] - 1), 100 * (sr["kappa_absorber"] / lg["second"]["slowest"] - 1)),
        v4_value="one physical channel, at the zone edge: %.4f and %.4f; measured %+.2f and %+.2f %%; the slower root %.4f is a truncation artefact"
                 % (ch["zone_edge"], sr["channels_at_measured_f"]["zone_edge"], s["dev_zone_edge_pct"], sr["dev_zone_edge_pct"], ch["boundary_artefact"]),
        reason="the original plane-wave layer kept every eigenvalue of the companion problem of a truncated square basis; the slower roots have their "
               "eigenvectors on the outermost G_y rings and no replica with |Re k_y| <= pi/a, and with a circular basis the same root moves to "
               "%.4f while the physical one stays at %.4f; the rod permittivity is now 11.9025, as in the FDTD runs"
               % (ch["truncation_check"]["bases"]["circle |m|<=16"]["kappa_slowest_artefact"], ch["truncation_check"]["bases"]["circle |m|<=16"]["kappa_physical"])))
    rd0 = (fe.get("row_decay") or [{}])[0]
    add(lambda: dict(
        v3_value="cladding channels %.4f at k_x = 0 and %.4f at k_x = pi/a; field-map decay %.4f, %.2f %% from the first"
                 % (lg["reference"]["kx0"], lg["reference"]["kxpi"], abs(rd0["cladding_k0"]["two_im_ky"]),
                    100 * (abs(rd0["cladding_k0"]["two_im_ky"]) / lg["reference"]["kx0"] - 1)),
        v4_value="%.4f at k_x = 0 and %.4f at k_x = pi/a; the field-map decay lies %.2f %% above the first"
                 % (ch["kx0"], ch["kxpi"], 100 * (abs(rd0["cladding_k0"]["two_im_ky"]) / ch["kx0"] - 1)),
        reason="the same truncation artefacts (%.4f and %.4f) had been taken for the cladding channels" % (ch["boundary_artefact_kx0"], ch["boundary_artefact_kxpi"])))
    e0 = fe.get("reference") or {}
    e1 = ([e for e in (fe.get("eta") or []) if abs(e["n"] - 1.45) < 1e-6] or [{}])[0]
    smb, smh = cb.get("second_map_bound") or {}, cb.get("second_map_hf") or {}
    add(lambda: dict(
        v3_value="eta_a >= %.3f, S_th = %.1f nm/RIU (%+.2f %%); %.1f nm/RIU at n_a = 1.45 (%+.2f %%)"
                 % (e0["eta_low"], e0["S_low"], e0["dev_low_pct"], e1["S_low"], smb["deviation_pct"]),
        v4_value="eta_a = %.3f, S_th = %.1f nm/RIU (%+.2f %%); %.1f nm/RIU at n_a = 1.45 (%+.2f %%)"
                 % (e0["eta_hf"], e0["S_hf"], e0["dev_hf_pct"], e1["S_hf"], smh["deviation_pct"]),
        reason="the stored permittivity map is what Meep writes, 12 / (sum of the diagonal of the inverse smoothed tensor) at the centred points "
               "(reproduced to %.1e), not the <eps> that E_z sees; over the consistent denominator the pure-analyte value is %.3f "
               "(%.1f nm/RIU, %+.2f %%) and is no bound above the measurement; mixed points hold %.3f of the energy"
               % (e0["reconstruction_max_abs"], e0["eta_hf_low"], e0["S_hf_low"], e0["dev_hf_low_pct"], e0["mixed_energy_fraction_hf"])))
    pf = ba.get("performance") or {}
    v8b = (pf.get("absorption") or {}).get("v8_basis") or {}
    add(lambda: dict(
        v3_value="Q_abs = %.0f at eta_a = %.4f; Q_tot = %.0f and %.0f; T_min = %.3f and %.3f; metric %.3f and %.3f; optimum %.0f; FOM in water %.0f and %.0f"
                 % (v8b["Q_abs"], v8b["eta"], v8b["Q_tot"]["N4"], v8b["Q_tot"]["N5"], v8b["T_min"]["N4"], v8b["T_min"]["N5"],
                    v8b["metric_fraction"]["N4"], v8b["metric_fraction"]["N5"], v8b["Q_w_opt"], v8b["FOM_water_N4"], v8b["FOM_water_N5"]),
        v4_value="Q_abs = lambda_r / (2 k S) = %.0f; Q_tot = %.0f and %.0f; T_min = %.3f and %.3f; metric %.3f and %.3f; optimum %.0f; FOM in water %.0f and %.0f"
                 % (pf["Q_abs"], pf["Q_tot"]["N4"], pf["Q_tot"]["N5"], pf["T_min"]["N4"], pf["T_min"]["N5"], pf["metric_fraction"]["N4"],
                    pf["metric_fraction"]["N5"], pf["Q_w_opt"], pf["FOM_water_N4"], pf["FOM_water_N5"]),
        reason="the loss to the medium is the imaginary part of the same first-order shift that gives the sensitivity, so the energy "
               "fraction the measured sensitivity requires, %.4f, applies" % pf["eta_used"]))
    qa, q = bsa.get("quadratic_fit") or {}, bs.get("quadratic_fit") or {}
    add(lambda: dict(
        v3_value="a/lambda_r = %.5f %+.5f r_d/a %+.5f (r_d/a)^2 over %d radii, R^2 = %.4f; linear R^2 = %.3f, rms %.1e; quadratic rms %.1e"
                 % (qa["coefficients_ascending"][0], qa["coefficients_ascending"][1], qa["coefficients_ascending"][2], qa["n"], qa["R2"],
                    bsa["linear_fit"]["R2"], bsa["linear_fit"]["rms"], qa["rms"]),
        v4_value="a/lambda_r = %.5f %+.5f r_d/a %+.5f (r_d/a)^2 over %d radii, R^2 = %.4f; linear R^2 = %.3f, rms %.1e; quadratic rms %.1e"
                 % (q["coefficients_ascending"][0], q["coefficients_ascending"][1], q["coefficients_ascending"][2], q["n"], q["R2"],
                    bs["linear_fit"]["R2"], bs["linear_fit"]["rms"], q["rms"]),
        reason="the radius 0.12a lies outside the gap window at every row count and is excluded by the validity criterion"))
    tfa, tf9 = cb.get("law_transfer_all_points") or {}, cb.get("law_transfer_to_fine_sweep") or {}
    dla, dl9 = cb.get("dlambda_drd_all_points") or {}, cb.get("dlambda_drd_at_reference") or {}
    add(lambda: dict(
        v3_value="the law predicts the %d points of the fine sweep to within %.1e; its derivative at 0.060a is %.2f nm/nm against %.2f from the stencil"
                 % (tfa["points"], tfa["max_abs_deviation"], dla["law_nm_per_nm"], dla["stencil_nm_per_nm"]),
        v4_value="the law predicts the %d points of the fine sweep to within %.1e; its derivative at 0.060a is %.2f nm/nm against %.2f from the stencil"
                 % (tf9["points"], tf9["max_abs_deviation"], dl9["law_nm_per_nm"], dl9["stencil_nm_per_nm"]),
        reason="excluded points removed from the law and from the fine sweep (0.110a, gap position 0.18)"))
    fra, fr9 = bsa.get("fits_by_row") or {}, bs.get("fits_by_row") or {}
    add(lambda: dict(
        v3_value="linear intercepts %.4f, %.4f, %.4f (%.1f %%) and %.0f %% in the linear slope"
                 % tuple([fra[k]["linear"]["intercept"] for k in sorted(fra, key=int)] + [bsa["intercept_spread_pct"], bsa["linear_slope_spread_pct"]]),
        v4_value="quadratic intercepts %.4f, %.4f, %.4f (%.2f %%) and curvatures %.2f, %.2f, %.3f"
                 % tuple([fr9[k]["quadratic_coefficients_ascending"][0] for k in sorted(fr9, key=int)] + [bs["quadratic_intercept_spread_pct"]]
                         + [fr9[k]["quadratic_coefficients_ascending"][2] for k in sorted(fr9, key=int)]),
        reason="the paper uses the quadratic law, so the dependence on the row count is quoted for it; admitted radii only"))
    rka, rk9 = bsa.get("raw_kappa") or {}, bs.get("raw_kappa") or {}
    v8a = ba.get("v8_all_points") or {}
    n8a, n89 = v8a.get("ncl8_leak_correction") or {}, ba.get("ncl8_leak_correction") or {}
    add(lambda: dict(
        v3_value="%d radii: ln(Q3/Q2) = %.2f +- %.2f, ln(Q4/Q3) = %.3f +- %.3f, leak removed %.3f +- %.3f, pooled kappa %.2f"
                 % (len(n8a["rd"]), rka["steps"]["2->3"]["mean"], rka["steps"]["2->3"]["sem"], n8a["raw_mean"], n8a["raw_sem"],
                    n8a["corrected_mean"], n8a["corrected_sem"], rka["kappa"]),
        v4_value="%d radii: ln(Q3/Q2) = %.2f +- %.2f, ln(Q4/Q3) = %.3f +- %.3f, leak removed %.3f +- %.3f, pooled kappa %.2f"
                 % (len(n89["rd"]), rk9["steps"]["2->3"]["mean"], rk9["steps"]["2->3"]["sem"], n89["raw_mean"], n89["raw_sem"],
                    n89["corrected_mean"], n89["corrected_sem"], rk9["kappa"]),
        reason="the radius 0.12a is excluded by the gap criterion; the correction still reduces the scatter, by a factor %.2f"
               % n89["scatter_ratio_raw_over_corrected"]))
    m8, m9 = v8a.get("matched_ratio_fit") or {}, ba.get("matched_ratio_fit") or {}
    do8, do9 = v8a.get("detrended_oscillation") or {}, ba.get("detrended_oscillation") or {}
    add(lambda: dict(
        v3_value="%d matched pairs: R = %.3f +- %.3f, D = %.2f +- %.2f a, c = %.3f +- %.3f; suppression %.1f"
                 % (m8["points"], m8["R"], m8["R_se"], m8["D"], m8["D_se"], m8["scale_c"], m8["scale_c_se"],
                    do8["pml_25a"]["detrended_ptp"] / do8["absorber"]["detrended_ptp"]),
        v4_value="%d matched pairs: R = %.3f +- %.3f, D = %.2f +- %.2f a, c = %.3f +- %.3f; suppression %.1f"
                 % (m9["points"], m9["R"], m9["R_se"], m9["D"], m9["D_se"], m9["scale_c"], m9["scale_c_se"],
                    do9["pml_25a"]["detrended_ptp"] / do9["absorber"]["detrended_ptp"]),
        reason="the pair at 0.110a is excluded by the gap criterion (position 0.18)"))
    dx0 = ((A.get("displacement") or {}).get("dx") or {})
    q0_, q1_ = (dx0.get(0.0) or dx0.get("0.0") or {}).get("quadratic_fall"), (dx0.get(0.1) or dx0.get("0.1") or {}).get("quadratic_fall")
    dxm = ((cb.get("cell_factor_estimates") or {}).get("dx_max") or {})
    add(lambda: dict(
        v3_value="displacement along the guide has no measurable effect: %.1f %% and %.1f %% at 0.08a, below the discretisation error"
                 % (-q0_["rel_change_pct"][-1], -q1_["rel_change_pct"][-1]),
        v4_value="a quadratic decrease of the 25a-cell value, %.2f %% and %.2f %% at 0.08a, %.2f and %.2f per a^2; the cell factor, "
                 "estimated from the resonance shift, falls by %.2f %% and %.2f %%"
                 % (-q0_["rel_change_pct"][-1], -q1_["rel_change_pct"][-1], -q0_["c_per_a2"], -q1_["c_per_a2"],
                    -dxm["0.00"]["change_pct"], -dxm["0.10"]["change_pct"]),
        reason="a ratio at fixed resolution carries no common discretisation error, and Q is even in dx by the mirror symmetry of "
               "lattice and cell; the cell-factor estimate exceeds the measured change, so the intrinsic change is below about one "
               "per cent and its sign is not resolved"))
    ps = A.get("preliminary_sweep") or {}
    add(lambda: dict(
        v3_value="a preliminary sweep (n_cl = 6, resolution 20) gives a median Q(dx = 0.08a)/Q(0) of %.2f (%.2f to %.2f) over %d groups"
                 % (ps["median_ratio"], ps["quartiles"][0], ps["quartiles"][1], ps["groups"]),
        v4_value="%d of the %d groups hold modes of the waveguide; the %d cavity groups give %.3f (%.3f to %.3f); the statement is withdrawn"
                 % (ps["groups"] - ps["cavity_groups"]["groups"], ps["groups"], ps["cavity_groups"]["groups"], ps["cavity_groups"]["median_ratio"],
                    ps["cavity_groups"]["quartiles"][0], ps["cavity_groups"]["quartiles"][1]),
        reason="the admitted records of that sweep at N_sep = 1 or r_d >= 0.13a hold low-Q modes of the finite waveguide (Q %.1f to %.1f), "
               "which the validity criterion admits" % (ps["waveguide_modes"]["Q_min"], ps["waveguide_modes"]["Q_max"])))
    rsa, rs9 = A.get("radius_sweep_all_points") or {}, A.get("radius_sweep") or {}
    add(lambda: dict(
        v3_value="a second peak of %.0f at %.3fa; period %.3fa, that is %.0f nm; Delta f = %.5f between the peaks"
                 % (rsa["Q"][rsa["rd"].index(rsa["peaks"][-1])], rsa["peaks"][-1], rsa["oscillation"]["period_a"],
                    rsa["oscillation"]["period_nm"], rsa["oscillation"]["delta_f"]),
        v4_value="one peak (%.3fa) and one trough (%.3fa) among the admitted radii, %.3fa apart; Q rises again to %.0f at %.3fa, the last "
                 "admitted radius" % (rs9["peaks"][0], rs9["troughs"][0], rs9["peak_trough_spacing_a"], rs9["Q"][-1], rs9["rd"][-1]),
        reason="the second peak needs the radius 0.110a, which lies outside the gap window"))
    return out


SECTIONS = ("measured", "derived", "second_radius", "field_map", "vertical", "barrier_row", "absorption", "systematic_runs", "termination",
            "dataset")


def markdown(reg: dict) -> str:
    lines = ["# Numbers registry", "", f"{reg.get('software', '')}, created {reg.get('created', '')}", "", reg.get("rule", ""), ""]
    for sec in SECTIONS:
        rows = reg.get(sec) or []
        if not rows:
            continue
        lines += [f"## {sec.replace('_', ' ')}", "", "| quantity | value | unit | source |", "|---|---|---|---|"]
        for r in rows:
            src = r.get("label") or r.get("derivation", "")
            unc = f" +- {r['uncertainty']}" if r.get("uncertainty") is not None else ""
            lines.append(f"| {r['quantity']} | {r['value']}{unc} | {r.get('unit', '')} | {src} |")
        lines.append("")
    if reg.get("superseded"):
        lines += ["## superseded", "", "| earlier value | current value | reason |", "|---|---|---|"]
        lines += [f"| {r['v3_value']} | {r['v4_value']} | {r['reason']} |" for r in reg["superseded"]] + [""]
    return "\n".join(lines)


__all__ = ["build", "markdown", "superseded", "superseded_v9", "SECTIONS"]
