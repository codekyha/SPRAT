"""Figure 4: the barrier mechanism and the finite-cell correction.

(a) Q(r_d) at N_sep = 4 for three guide terminations: PML 25a cell, PML 49a cell, absorber.
(b) The measured cell factor G(r_d) = Q_PML / Q_absorber at matched r_d (leak removed), with the
    two-reflector fit.
(c) Barrier coefficient against resonant frequency: the zone-edge channel of the complex band
    structure (the slowest physical Bloch channel at the guided-mode wavevector; the slower root of
    a square plane-wave basis is a truncation artefact and is not drawn), the reflectionless leak-free values at r_d = 0.060a and 0.100a, and the
    log-linear fits to the raw 25a-cell values at the same two radii.
(d) Q^-1 exp(kappa N) against exp(kappa N) for the PML series at clearances 8 and 7; the
    intercept is the barrier prefactor 1/A and the slopes are 1/Q_top.

A port of the figure script of the paper (v7, ``sekil_mekanizma_v7.py``): same reduction, same
printed numbers, same layout, sizes and colours.  Colour carries one meaning throughout: blue =
PML-terminated 25a cell (raw), vermillion = absorber-terminated (reflectionless), orange = 49a
cell or clearance 7, green = the zone-edge channel of the complex band structure, open grey = a
record outside the admission window of the paper (shown, not fitted).

Inputs: the ``barrier`` block of ``analysis.json`` (three-point solve, leak law, cell factor,
four-point series, matched-ratio fit, the absorber and 49a sweeps, the second radius), the
``kappa_series``, ``resolved_pair`` and ``radius_sweep`` blocks, and the complex-band curve of the
plane-wave record.
"""

from __future__ import annotations

import math

import numpy as np

from . import style
from .inputs import Inputs, MissingData, key

CANVAS_CM = (16.0, 11.8)
KAPPA_ERR_SECOND = 0.0065       # fallback error bar at r_d = 0.100a when the analysis carries no uncertainty
CHANNEL_RANGE = (0.2860, 0.3120)


def airy(p, beta):
    """Two-reflector cell factor G(beta) = (1 + R - 2 sqrt(R) cos(2 beta D + phi)) / (1 - R)."""
    R, D, phi = p
    return (1 + R - 2 * math.sqrt(R) * np.cos(2 * beta * D + phi)) / (1 - R)


def _one_digit(x: float) -> float:
    """An uncertainty rounded to one significant digit, as drawn (Monte Carlo 0.00195 -> 0.002)."""
    return float(f"{x:.1g}") if x else 0.0


def numbers(inp: Inputs) -> dict:
    """The quantities drawn in figure 4, from the analysis outputs and the records."""
    ba = inp.block("barrier")
    for k in ("inputs", "three_point_solve", "leak", "G", "series", "matched_ratio_fit", "absorber_sweep",
              "guide49_sweep", "fp_shape_fits"):
        if not ba.get(k):
            raise MissingData(f"analysis block 'barrier.{k}' is missing")
    inp_q = {k: v["Q"] for k, v in ba["inputs"].items() if v.get("Q") is not None}
    need = ("Q4_pml_cl12", "Q5_pml_cl13", "Q5_pml_cl12", "Q6_pml_cl14")
    if any(k not in inp_q for k in need):
        raise MissingData("barrier inputs incomplete: " + ", ".join(k for k in need if k not in inp_q))
    mr = ba["matched_ratio_fit"]
    ref_sweep = inp.block("radius_sweep")
    f_ref = ref_sweep["f"][[round(x, 4) for x in ref_sweep["rd"]].index(0.06)]
    pwe = inp.pwe
    curve = pwe.get("kappa_curve") or []
    if not curve:
        raise MissingData("the plane-wave record has no kappa_curve")
    sr = ba.get("second_radius")
    out = dict(
        Q4=inp_q["Q4_pml_cl12"], Q5c8=inp_q["Q5_pml_cl13"], Q5c7=inp_q["Q5_pml_cl12"], Q6=inp_q["Q6_pml_cl14"],
        K3=ba["three_point_solve"]["kappa"], AMP=ba["three_point_solve"]["A"],
        INV8=1 / ba["leak"]["Q_top_8"], INV7=1 / ba["leak"]["Q_top_7"], KTOP=ba["leak"]["kappa_top"],
        G5=key(ba["G"], 5), KAP=ba["series"]["kappa"], KAP_ERR=_one_digit(ba["series"].get("kappa_sd_monte_carlo", 0.0)),
        RD=np.asarray(mr["rd"], float), GG=np.asarray(mr["G"], float), BET=math.pi * np.asarray(mr["beta_over_pi"], float),
        FIT=(mr["R"], mr["D"], mr["phi"]), R_SE=mr["R_se"], D_SE=mr["D_se"], RMS=mr["rms_lnG"],
        fp_rd=ba["fp_shape_fits"]["rd"], fp_Q=ba["fp_shape_fits"]["Q"],
        g49=list(zip(ba["guide49_sweep"]["rd"], ba["guide49_sweep"]["Q"])),
        pabs=list(zip(ba["absorber_sweep"]["rd"], ba["absorber_sweep"]["Q"])),
        pabs_ex=list(zip((ba["absorber_sweep"].get("excluded") or {}).get("rd", []), (ba["absorber_sweep"].get("excluded") or {}).get("Q", []))),
        fp_ex=list(zip((ref_sweep.get("excluded") or {}).get("rd", []), (ref_sweep.get("excluded") or {}).get("Q", []))),
        CF=np.array([r["f"] for r in curve]), CE=np.array([r["kappa_pred_zone_edge"] for r in curve]),
        CS=np.array([r["kappa_pred"] for r in curve]), F_REF=float(f_ref),
        KAP2=(sr["kappa_absorber"] if sr else None), F2=(sr["f_mean"] if sr else None),
        KAP2_ERR=(((sr or {}).get("kappa_absorber_uncertainty") or {}).get("total") or KAPPA_ERR_SECOND),
        RAW_REF=inp.block("kappa_series")["kappa"],
        guide_end=0.5 * float(inp.reference()["structure"]["cell"]["guide_periods"]))
    try:                                   # the raw fit at the second radius, placed at the plane-wave mode frequency
        rp = inp.block("resolved_pair")
        out["RAW_SECOND"] = (float(pwe["kappa"]["second"]["f"]), rp["kappa"])
    except (MissingData, KeyError, TypeError):
        out["RAW_SECOND"] = None
    return out


def report_lines(v: dict) -> list[str]:
    """The numbers printed by the figure script of the paper, in the same format."""
    L = ["  kappa (four-point)      = %.4f" % v["KAP"]]
    if v["KAP2"] is not None:
        L.append("  kappa (r_d = 0.100a)    = %.4f at f = %.7f" % (v["KAP2"], v["F2"]))
    L += ["  kappa (three-point)     = %.4f" % v["K3"],
          "  Q_top(8) = %.4g   Q_top(7) = %.4g   kappa_top = %.4f" % (1 / v["INV8"], 1 / v["INV7"], v["KTOP"]),
          "  G(N=5) = %.4f" % v["G5"],
          "  two-mirror fit: R = %.4f +- %.4f, D = %.3f +- %.3f a, rms(lnG) = %.4f"
          % (v["FIT"][0], v["R_SE"], v["FIT"][1], v["D_SE"], v["RMS"])]
    return L


def figure4(inp: Inputs, log=print):
    v = numbers(inp)
    fmt = style.formatter
    letter = style.LETTER_MECHANISM
    BLUE, ORANGE, VERM, GREEN, DARK = style.BLUE, style.ORANGE, style.VERMILLION, style.GREEN, style.DARK
    fig = style.new_figure(CANVAS_CM[0] * style.CM * 72, CANVAS_CM[1] * style.CM * 72)
    axs = fig.subplots(2, 2)
    a, b, c, d = axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1]
    for ax, s in zip((a, b, c, d), "abcd"):
        style.panel_letter(ax, f"({s})", **letter)

    # ---- (a) three terminations
    a.axvline(0.060, color="0.75", lw=0.7, ls="-.", zorder=0)
    a.annotate("reference\nradius", xy=(0.0617, 1.93e3), fontsize=7.0, color=DARK, va="top")
    a.plot(v["fp_rd"], v["fp_Q"], "o-", color=BLUE, mfc="white", mec=BLUE, mew=1.0, lw=1.1, ms=4.0,
           label="PML, $25a$ guide")
    a.plot([p[0] for p in v["g49"]], [p[1] for p in v["g49"]], "^--", color=ORANGE, mfc="white", mec=ORANGE,
           mew=1.0, lw=1.0, ms=4.0, label="PML, $49a$ guide")
    a.plot([p[0] for p in v["pabs"]], [p[1] for p in v["pabs"]], "s-", color=VERM, mfc=VERM, mec=VERM,
           lw=1.4, ms=3.8, label="absorber, reflectionless")
    if v["fp_ex"] or v["pabs_ex"]:
        # explained in the caption; no legend entry, so that the legend keeps its place below the curves
        a.plot([p[0] for p in v["fp_ex"]], [p[1] for p in v["fp_ex"]], "o", color=style.GREY, mfc="white", mec=style.GREY,
               mew=0.9, ms=4.0, ls="none")
        a.plot([p[0] for p in v["pabs_ex"]], [p[1] for p in v["pabs_ex"]], "s", color=style.GREY, mfc="white",
               mec=style.GREY, mew=0.9, ms=3.8, ls="none")
    a.set_xlabel("Defect radius, $r_d/a$")
    a.set_ylabel("Quality factor, $Q$")
    a.set_yscale("log")
    a.set_ylim(1.4e3, 1.15e4)
    a.set_xlim(0.045, 0.115)
    a.legend(loc="upper right", handlelength=2.2, labelspacing=0.25, borderaxespad=0.2)
    a.xaxis.set_major_formatter(fmt(3))

    # ---- (b) measured G with the two-reflector fit
    RD, GG, BET = v["RD"], v["GG"], v["BET"]
    bb = np.linspace(BET.min(), BET.max(), 400)
    rr = np.interp(bb, BET[::-1], RD[::-1])
    b.axhline(1.0, color="0.70", lw=0.7, ls=":", zorder=0)
    b.plot(rr, airy(v["FIT"], bb), "-", color=VERM, lw=1.3, label="two-reflector fit")
    b.plot(RD, GG, "o", color=BLUE, mfc=BLUE, mec=BLUE, ms=4.4, ls="none", label="measured")
    ex = {round(r, 4): q for r, q in v["fp_ex"]}
    exr = [(r, ex[round(r, 4)] / q) for r, q in v["pabs_ex"] if round(r, 4) in ex]
    if exr:
        b.plot([p[0] for p in exr], [p[1] for p in exr], "o", color=style.GREY, mfc="white", mec=style.GREY, mew=0.9,
               ms=4.4, ls="none")
    b.set_xlabel("Defect radius, $r_d/a$")
    b.set_ylabel("Cell factor, $G$")
    b.set_xlim(0.045, 0.115)
    b.set_ylim(0.35, 2.25)
    b.text(0.047, 2.20, "$R = %.3f \\pm %.3f$\n$D = %.1f \\pm %.1f\\,a$  (geometric $%sa$)"
           % (v["FIT"][0], v["R_SE"], v["FIT"][1], v["D_SE"], style.num(v["guide_end"], 1)),
           fontsize=7.2, va="top", color=DARK)
    b.legend(loc="lower left", handlelength=2.2, labelspacing=0.25, borderaxespad=0.2)
    b.xaxis.set_major_formatter(fmt(3))
    b.yaxis.set_major_formatter(fmt(1))

    # ---- (c) channel identification across the guided band, at two defect radii
    CF, CE = v["CF"], v["CE"]           # v["CS"] (the slower root of a square basis) is an artefact and is not drawn
    sel = (CF >= CHANNEL_RANGE[0]) & (CF <= CHANNEL_RANGE[1])
    c.plot(CF[sel], CE[sel], "--", color=GREEN, lw=1.4, label="zone-edge channel")
    c.errorbar([v["F_REF"]], [v["KAP"]], yerr=[v["KAP_ERR"]], fmt="o", color=VERM, mfc=VERM, mec=VERM,
               ms=5.0, capsize=2.4, elinewidth=1.0, label="reflectionless, leak-free")
    if v["KAP2"] is not None:
        c.errorbar([v["F2"]], [v["KAP2"]], yerr=[v["KAP2_ERR"]], fmt="o", color=VERM, mfc=VERM, mec=VERM, ms=5.0,
                   capsize=2.4, elinewidth=1.0)
    raw_f, raw_k = [v["F_REF"]], [v["RAW_REF"]]
    if v["RAW_SECOND"] is not None:
        raw_f.append(v["RAW_SECOND"][0])
        raw_k.append(v["RAW_SECOND"][1])
    c.plot(raw_f, raw_k, "s", color=BLUE, mfc="white", mec=BLUE, mew=1.1, ms=5.0, ls="none",
           label="raw fit, $25a$ cell")
    tag = dict(fontsize=7.0, color=DARK, ha="left", va="top", textcoords="offset points",
               bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none"))
    c.annotate("$r_d = 0.060a$", xy=(v["F_REF"], v["KAP"]), xytext=(7, -2), **tag)
    if v["KAP2"] is not None:
        c.annotate("$r_d = 0.100a$", xy=(v["F2"], v["KAP2"]), xytext=(7, -2), **tag)
    c.set_xlabel("Resonant frequency, $a/\\lambda_r$")
    c.set_ylabel("Barrier coefficient, $\\kappa$")
    c.set_xlim(0.2865, 0.3115)
    c.set_ylim(1.55, 2.02)
    c.set_xticks([0.290, 0.295, 0.300, 0.305, 0.310])
    c.legend(loc="lower right", handlelength=2.2, labelspacing=0.25, borderaxespad=0.2)
    c.xaxis.set_major_formatter(fmt(3))
    c.yaxis.set_major_formatter(fmt(2))

    # ---- (d) linearised leak: Q^{-1} e^{kN} = 1/A + (1/Q_top) e^{kN}
    K3, AMP, INV8, INV7 = v["K3"], v["AMP"], v["INV8"], v["INV7"]
    EX = np.array([math.exp(K3 * n) for n in (4, 5, 6)])
    YY = np.array([1 / v["Q4"], 1 / v["Q5c8"], 1 / v["Q6"]]) * EX
    y7 = math.exp(K3 * 5) / v["Q5c7"]
    xl = np.linspace(0, EX.max() * 1.10, 50)
    d.axhline(1 / AMP, color="0.70", lw=0.7, ls=":", zorder=0)
    d.plot(xl * 1e-5, 1 / AMP + INV8 * xl, "-", color=BLUE, lw=1.3,
           label="clearance 8, $Q_{\\mathrm{top}} = %.2f \\times 10^{6}$" % (1e-6 / INV8))
    xl7 = np.linspace(0, math.exp(K3 * 5) * 1.45, 20)          # clearance 7: the intercept and one point
    d.plot(xl7 * 1e-5, 1 / AMP + INV7 * xl7, "--", color=ORANGE, lw=1.2,
           label="clearance 7, $Q_{\\mathrm{top}} = %.2f \\times 10^{5}$" % (1e-5 / INV7))
    d.plot(EX * 1e-5, YY, "o", color=BLUE, mfc=BLUE, mec=BLUE, ms=5.0, ls="none")
    d.plot([math.exp(K3 * 5) * 1e-5], [y7], "s", color=ORANGE, mfc="white", mec=ORANGE, mew=1.1, ms=5.0, ls="none")
    place = {4: ("left", "top", 6, -3.5), 5: ("left", "top", 6, 0), 6: ("right", "baseline", -6, 6)}
    for n, x, y in zip((4, 5, 6), EX, YY):
        ha, va, dx, dy = place[n]
        d.annotate("$N_{\\mathrm{sep}} = %d$" % n, xy=(x * 1e-5, y), xytext=(dx, dy), textcoords="offset points",
                   fontsize=7.0, color=DARK, ha=ha, va=va)
    d.set_xlabel("$\\exp(\\kappa N_{\\mathrm{sep}})\\ \\times 10^{-5}$")
    d.set_ylabel("$Q^{-1}\\exp(\\kappa N_{\\mathrm{sep}})$")
    d.set_xlim(0, EX.max() * 1.10 * 1e-5)
    d.set_ylim(1 / AMP - 0.012, max(YY.max(), y7) + 0.020)
    d.annotate("intercept $= A^{-1} = %.4f$" % (1 / AMP), xy=(EX.max() * 1.05e-5, 1 / AMP), xytext=(-2, -3),
               textcoords="offset points", fontsize=7.0, color=DARK, ha="right", va="top")
    d.set_xticks([0.0, 0.3, 0.6, 0.9, 1.2])
    d.legend(loc="upper left", handlelength=2.2, labelspacing=0.25, borderaxespad=0.2)
    d.xaxis.set_major_formatter(fmt(1))
    d.yaxis.set_major_formatter(fmt(2))

    fig.tight_layout(pad=0.4, w_pad=2.0, h_pad=1.6)
    for line in report_lines(v):
        log(line)
    return fig


__all__ = ["figure4", "numbers", "report_lines", "airy", "KAPPA_ERR_SECOND"]
