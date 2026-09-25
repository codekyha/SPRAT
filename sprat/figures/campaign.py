"""Figures 2, 3, 5, 6 and 7: the band gap, the defect and coupling laws, the displacement of the
defect rod, the transmission spectra and the sensing performance.

Every series is read from the analysis outputs (``analysis.json``, ``spectra.json``), except
the band edges, which come from the band records, and the resonance of the reference geometry,
which comes from its record.  The layout (canvas and axes in points) reproduces the figures of
the paper; the data, the fits and the printed numbers follow the records.
"""

from __future__ import annotations

import math

import matplotlib
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

from . import style
from .inputs import Inputs, MissingData, key

N_LABEL = "Analyte refractive index, $n_a$"
RD_LABEL = "Defect radius, $r_d/a$"
Q_LABEL = "Quality factor, $Q$"
SEP = "$N_{\\mathrm{sep}}$"
CL = "$n_{\\mathrm{cl}}$"
OPEN = dict(mfc="white", mew=1.0)          # open markers of the campaign figures


def _sci(x, pos=None) -> str:
    """Tick label of a fixed log tick: 10^3, 2x10^3, 2.5x10^-5 (no operator spacing around the times sign)."""
    if x <= 0:
        return ""
    e = math.floor(math.log10(x) + 1e-9)
    c = round(x / 10 ** e, 6)
    return "$10^{%d}$" % e if abs(c - 1) < 1e-9 else "$%g{\\times}10^{%d}$" % (c, e)


SCI = FuncFormatter(_sci)

LAYOUT = {
    2: dict(canvas=(454.52, 174.06), axes=[(38.81, 18.13, 173.99, 124.16), (268.47, 18.13, 173.99, 124.16)]),
    3: dict(canvas=(455.89, 346.86), axes=[(38.81, 18.15, 172.5, 125.4), (269.96, 18.15, 172.5, 125.4),
                                           (38.81, 188.69, 172.5, 125.4), (269.96, 188.69, 172.5, 125.4)]),
    5: dict(canvas=(450.66, 175.06), axes=[(39.81, 18.13, 172.5, 124.16), (270.96, 18.13, 172.5, 124.16)]),
    6: dict(canvas=(450.98, 190.54), axes=[(35.60, 18.45, 177.47, 140.32), (266.31, 18.45, 177.47, 140.32)],
            colorbar=(46.0, 139.83, 50.0, 4.91), inset=(291.16, 68.97, 60.34, 30.87),
            inset_a=(161.0, 86.0, 47.0, 44.0)),
    7: dict(canvas=(461.33, 330.96), axes=[(40.76, 18.05, 155.67, 120.14), (252.47, 18.05, 155.67, 120.14),
                                           (40.76, 179.04, 155.67, 120.15), (252.47, 179.04, 155.67, 120.15)],
            inset=(131.05, 92.54, 59.15, 31.24)),
}


def _canvas(number: int):
    lay = LAYOUT[number]
    fig = style.new_figure(*lay["canvas"])
    axs = [style.axes_pt(fig, *box) for box in lay["axes"]]
    for ax, letter in zip(axs, "abcd"):
        style.panel_letter(ax, f"({letter})")
    return fig, axs


def _design_range(inp: Inputs) -> tuple[float, float]:
    """The analyte range of the design: the span of the analyte sweep of the reference geometry."""
    n = inp.block("analyte_sweep", "n")
    return float(min(n)), float(max(n))


def _loglocator(ax, ticks, minor=True):
    """Fixed major ticks on a log axis; the default minor ticks stay unless ``minor`` is False."""
    ax.yaxis.set_major_locator(FixedLocator(ticks))
    if not minor:
        ax.yaxis.set_minor_locator(NullLocator())


# ============================================================================ figure 2
def figure2(inp: Inputs):
    """TM gap edges and relative gap width against the analyte index (MPB sweep record)."""
    sweep = [x for x in inp.bands("sweep")["sweep"] if x["gap"]["exists"]]
    if len(sweep) < 2:
        raise MissingData("the bands sweep record holds fewer than two points with a gap")
    sweep.sort(key=lambda x: x["n"])
    n = np.array([x["n"] for x in sweep])
    lo = np.array([x["gap"]["lower"] for x in sweep])
    hi = np.array([x["gap"]["upper"] for x in sweep])
    width = np.array([x["gap"]["relative_width_pct"] for x in sweep])
    sw = inp.block("analyte_sweep")
    n_sw = np.asarray(sw["n"], float)
    f_sw = float(sw["a_nm"]) / np.asarray(sw["lambda_nm"], float)
    d0, d1 = _design_range(inp)

    fig, (a, b) = _canvas(2)
    a.fill_between(n, lo, hi, color=style.GAP, lw=0, label="TM band gap")
    a.plot(n, lo, color=style.BLUE, lw=1.1, label="lower edge")
    a.plot(n, hi, color=style.VERMILLION, lw=1.1, label="upper edge")
    a.axvspan(d0, d1, color=style.BAND, alpha=0.12, lw=0, label="design range")
    a.plot(n_sw, f_sw, color="k", lw=0.9, ls=":", label="$f_r$ of the reference design")
    a.set_xlim(n[0], n[-1])
    a.set_ylim(0.26, 0.43)
    a.set_yticks(np.arange(0.26, 0.4201, 0.02))
    a.set_xlabel(N_LABEL)
    a.set_ylabel("Normalised frequency, $a/\\lambda$")
    a.legend(loc="upper right")

    b.axvspan(d0, d1, color=style.BAND, alpha=0.12, lw=0)
    b.plot(n, width, "o-", color=style.BLUE, lw=1.0, ms=4.0, mec=style.BLUE, **OPEN)
    b.set_xlim(n[0], n[-1])
    b.set_ylim(10, 42)
    b.set_xlabel(N_LABEL)
    b.set_ylabel("Relative gap width (%)")
    return fig


# ============================================================================ figure 3
def figure3(inp: Inputs):
    """Defect-radius laws (a, b), coupling against the separating rows (c) and the cladding ladder (d)."""
    bs = inp.block("base_series")
    rs = inp.block("radius_sweep")
    ks = inp.block("kappa_series")
    rp = inp.block("resolved_pair")
    cl = inp.block("cladding_ladder")
    gap = inp.bands("bulk")["gap"]
    rows = sorted(int(k) for k in bs["series"])
    series = {N: key(bs["series"], N) for N in rows}
    ncl, res = bs["cladding_rows"], bs["resolution"]
    look = {2: ("^", style.ORANGE, ":"), 3: ("s", style.SKY, "--"), 4: ("o", style.GREEN, "-")}
    fine_label = f"{SEP} = 4, {CL} = 12, res 24"

    fig, (a, b, c, d) = _canvas(3)

    # ---- (a) resonance frequency against the defect radius
    a.axhspan(gap["lower"], gap["upper"], color=style.GAP, lw=0, zorder=0)
    a.text(0.975, 0.936, "TM band gap", transform=a.transAxes, ha="right", va="top", fontsize=6.8, color=style.NOTE)
    for N in rows:
        m, col, _ = look.get(N, ("o", style.GREY, "-"))
        a.plot(series[N]["rd"], series[N]["f"], m, color=col, mec=col, ms=4.2, ls="none", zorder=3,
               label=f"{SEP} = {N}, {CL} = {ncl}, res {res}", **OPEN)
    a.plot(rs["rd"], rs["f"], "D", color="k", ms=3.2, mew=1.0, ls="none", zorder=4, label=fine_label)
    ex_rd = [x for N in rows for x in (series[N].get("excluded") or {}).get("rd", [])] + list((rs.get("excluded") or {}).get("rd", []))
    ex_f = [x for N in rows for x in (series[N].get("excluded") or {}).get("f", [])] + list((rs.get("excluded") or {}).get("f", []))
    if ex_rd:
        a.plot(ex_rd, ex_f, "o", color=style.GREY, mfc="white", mec=style.GREY, mew=0.8, ms=4.2, ls="none", zorder=3,
               label="outside the admission window")
    x = np.linspace(min(series[rows[-1]]["rd"]), max(series[rows[-1]]["rd"]), 200)
    q = bs.get("quadratic_fit")
    if q:
        c0, c1, c2 = q["coefficients_ascending"]
        a.plot(x, c0 + c1 * x + c2 * x ** 2, "-", color=style.VERMILLION, lw=1.1, zorder=2,
               label="quadratic fit, equation (5)")
    lf = bs.get("linear_fit")
    if lf:
        c0, c1 = lf["coefficients_ascending"]
        a.plot(x, c0 + c1 * x, ":", color=style.GREY, lw=1.0, zorder=2, label="linear fit")
    a.set_xlim(-0.005, 0.125)
    a.set_ylim(0.255, 0.35)
    a.set_yticks([0.26, 0.28, 0.30, 0.32, 0.34])
    a.set_xlabel(RD_LABEL)
    a.set_ylabel("Resonant frequency, $a/\\lambda_r$")
    a.legend(loc="lower left", fontsize=6.3)

    # ---- (b) quality factor against the defect radius
    for N in rows:
        m, col, ls = look.get(N, ("o", style.GREY, "-"))
        b.semilogy(series[N]["rd"], series[N]["Q"], marker=m, ls=ls, color=col, mec=col, lw=0.9, ms=3.6, mfc="white",
                   mew=0.9, zorder=3, label=f"{SEP} = {N}, {CL} = {ncl}, res {res}")
    b.semilogy(rs["rd"], rs["Q"], "D-", color="k", lw=1.1, ms=3.4, mew=1.0, zorder=4, label=fine_label)
    ex_q = [x for N in rows for x in (series[N].get("excluded") or {}).get("Q", [])] + list((rs.get("excluded") or {}).get("Q", []))
    if ex_rd:
        b.semilogy(ex_rd, ex_q, "o", color=style.GREY, mfc="white", mec=style.GREY, mew=0.8, ms=3.6, ls="none", zorder=3,
                   label="outside the admission window")
    peaks = rs.get("peaks") or []
    if len(peaks) >= 2:
        p1, p2 = peaks[0], peaks[-1]
        y = 1.25e4
        b.annotate("", xy=(p2, y), xytext=(p1, y), arrowprops=dict(arrowstyle="<->", lw=0.7, color=style.NOTE,
                                                                   mutation_scale=8))
        b.text(0.5 * (p1 + p2), y * 1.12, style.num(p2 - p1, 3) + "$a$", ha="center", va="bottom", fontsize=6.8,
               color=style.NOTE)
    b.set_xlim(-0.005, 0.125)
    b.set_ylim(40, 1.2e5)
    b.set_yticks([1e2, 1e3, 1e4])
    b.set_xlabel(RD_LABEL)
    b.set_ylabel(Q_LABEL)
    b.legend(loc="upper left", ncol=2, fontsize=6.3, columnspacing=1.0)

    # ---- (c) quality factor against the separating rows
    common = bs.get("raw_kappa", {}).get("common_rd") or []
    for rd in common:
        qs = [series[N]["Q"][series[N]["rd"].index(rd)] for N in rows]
        c.semilogy(rows, qs, "-", color="0.82", lw=0.6, zorder=1)
    mean_q = bs.get("raw_kappa", {}).get("mean_Q")
    if mean_q:
        # filled, so that the means cannot be read as the open grey circles that mark excluded records in (a) and (b)
        c.semilogy(rows, mean_q, "o", color=style.GREY, mec=style.GREY, mfc=style.GREY, ms=4.4, ls="none", zorder=3,
                   label=f"{CL} = {ncl}, res {res}: individual $r_d$ (lines), mean (circles)")
    NN = np.linspace(1.75, 5.3, 50)
    if rp:
        c.semilogy([2, 4], [rp["Q2"], rp["Q4"]], "s", color=style.BLUE, mec=style.BLUE, ms=4.4, mew=1.0, ls="none",
                   zorder=4, label=f"$r_d$ = 0.100$a$, {CL} \u2212 {SEP} \u2265 6")
    c.semilogy(ks["N"], ks["Q"], "D", color="k", ms=4.6, mew=1.0, ls="none", zorder=5,
               label=f"$r_d$ = 0.060$a$, {CL} = 12, res 24")
    c.semilogy(NN, np.exp(ks["intercept"] + ks["kappa"] * NN), "-", color="k", lw=1.0, zorder=2,
               label=f"raw fit, \u03ba = {ks['kappa']:.2f} ($r_d$ = 0.060$a$)")
    if rp:
        c.semilogy(NN, rp["Q2"] * np.exp(rp["kappa"] * (NN - 2)), "--", color=style.BLUE, lw=1.0, zorder=2,
                   label=f"raw fit, \u03ba = {rp['kappa']:.2f} ($r_d$ = 0.100$a$)")
    c.set_xlim(1.6, 5.4)
    c.set_ylim(40, 3e5)
    c.set_xticks([2, 3, 4, 5])
    c.set_xlabel(f"Separating rows, {SEP}")
    c.set_ylabel(Q_LABEL)
    c.legend(loc="upper left", fontsize=6.3)

    # ---- (d) cladding ladder at r_d = 0.100a, resolution 32
    s4, s5 = key(cl["series"], 4), key(cl["series"], 5)
    d.semilogy(s4["cladding_rows"], s4["Q"], "o-", color=style.GREEN, mec=style.GREEN, lw=0.9, ms=4.4, zorder=3,
               label=f"{SEP} = 4", **OPEN)
    d.semilogy(s5["cladding_rows"], s5["Q"], "s-", color=style.VERMILLION, mec=style.VERMILLION, lw=0.9, ms=4.4,
               zorder=3, label=f"{SEP} = 5", **OPEN)
    qlim = s4["Q_limit"][-1]
    d.axhline(qlim, color=style.NOTE, lw=0.8, ls="--", label="$Q_{\\mathrm{lim}}$")
    q_top = s4["Q"][s4["cladding_rows"].index(max(s4["cladding_rows"]))]
    d.axhline(q_top, color=style.GREEN, lw=0.6, ls=":")
    for nc, q in zip(s4["cladding_rows"], s4["Q"]):
        if nc != max(s4["cladding_rows"]):
            # the first label sits below and to the left of its point, clear of the steep N_sep = 5 line that
            # starts at the same abscissa; the others sit centred below their points
            first = nc == min(s4["cladding_rows"])
            d.annotate(style.pct(100 * (q / q_top - 1)), xy=(nc, q), xytext=(-3.5, -4.0) if first else (0, -7.0),
                       textcoords="offset points", ha="right" if first else "center", va="top", fontsize=6.5,
                       color=style.NOTE)
    d.set_xlim(7.1, 12.6)
    d.set_ylim(1500, 4e4)
    d.set_xticks(s4["cladding_rows"])
    _loglocator(d, [2e3, 5e3, 1e4, 2e4])
    d.yaxis.set_major_formatter(SCI)
    d.set_xlabel(f"Cladding periods, {CL}")
    d.set_ylabel(Q_LABEL)
    d.legend(loc="upper left")
    return fig


# ============================================================================ figure 5


def figure5(inp: Inputs):
    """Quality factor against the displacement of the defect rod, across (a) and along (b) the guide."""
    disp = inp.block("displacement")
    fig, (a, b) = _canvas(5)
    # no discretisation band: a ratio at fixed resolution carries no common-mode discretisation error.  Panel (b) has
    # its own, expanded scale, on which the second-order fall with dx and its quadratic fit are visible.
    for ax in (a, b):
        ax.axhline(1.0, color="0.6", lw=0.6, zorder=1)
    a.set_ylim(0.92, 1.06)
    a.set_yticks([0.92, 0.96, 1.00, 1.04])
    b.set_ylim(0.9925, 1.0015)
    b.set_yticks([0.993, 0.995, 0.997, 0.999, 1.001])
    b.yaxis.set_major_formatter(style.formatter(3))

    for rd, m, col in ((0.05, "s", style.BLUE), (0.06, "o", style.VERMILLION)):
        s = key(disp["dy"], rd)
        a.plot(s["dy"], s["ratio"], marker=m, color=col, mec=col, lw=0.9, ms=4.4, zorder=3,
               label=f"$r_d$ = {rd:.3f}$a$ ($Q_0$ = {s['Q0']:.0f})", **OPEN)
    a.legend(loc="lower right")
    a.set_xlim(-0.11, 0.11)
    a.set_xticks([-0.10, -0.05, 0.0, 0.05, 0.10])
    a.set_xlabel("Transverse displacement, $\\delta_y/a$")
    a.set_ylabel("$Q(\\delta_y)\\,/\\,Q(0)$")

    xq = np.linspace(0.0, 0.082, 100)
    for dy, m, col in ((0.0, "o", style.VERMILLION), (0.1, "^", style.GREEN)):
        s = key(disp["dx"], dy)
        b.plot(s["dx"], s["ratio"], marker=m, color=col, mec=col, lw=0, ms=4.4, zorder=3,
               label=f"$\\delta_y$ = {dy:.2f}$a$ ($Q_0$ = {s['Q0']:.0f})", **OPEN)
        qf = s.get("quadratic_fall")
        if qf:
            b.plot(xq, 1.0 + qf["c_per_a2"] * xq ** 2, "--", color=col, lw=0.9, zorder=2)
    b.plot([], [], "--", color="0.35", lw=0.9, label="$1 - c\\,\\delta_x^2$ fits")
    b.legend(loc="lower left", title="$r_d$ = 0.060$a$")
    b.set_xlim(-0.005, 0.085)
    b.set_xticks([0.0, 0.02, 0.04, 0.06, 0.08])
    b.set_xlabel("Axial displacement, $\\delta_x/a$")
    b.set_ylabel("$Q(\\delta_x)\\,/\\,Q(0)$")
    return fig


# ============================================================================ figure 6
def figure6(inp: Inputs):
    """Normalised spectra collapsed on the linewidth (a) and the independent linewidth check (b)."""
    spectra = inp.spectra
    ns = np.array([s["n"] for s in spectra])
    norm = Normalize(vmin=float(ns.min()), vmax=float(ns.max()))
    cmap = matplotlib.colormaps["viridis"]
    lay = LAYOUT[6]
    fig, (a, b) = _canvas(6)

    # ---- (a) the seven spectra on the normalised detuning 2 (lambda - lambda_r) / fwhm
    x = np.linspace(-8, 8, 801)
    for s in spectra:
        f = s["fit"]
        col = cmap(norm(s["n"]))
        lam, T = np.asarray(s["lambda_nm"]), np.asarray(s["T"])
        e = 2 * (lam - f["lambda_r"]) / f["fwhm"]
        m = np.abs(e) <= 8.0
        a.plot(e[m], T[m], "o", ms=2.0, mfc="none", mec=col, mew=0.6, ls="none", zorder=2)
        a.plot(x, f["base"] + f["T0"] * (f["q"] + x) ** 2 / (1 + x ** 2), "-", color=col, lw=0.8, zorder=3)
    a.plot(x, x ** 2 / (1 + x ** 2), "--", color="k", lw=1.1, zorder=4)
    handles = [Line2D([], [], marker="o", ls="none", ms=3.0, mfc="none", mec="0.4", mew=1.0),
               Line2D([], [], color="0.4", lw=0.8), Line2D([], [], color="k", lw=1.1, ls="--")]
    a.legend(handles, ["FDTD", "Fano line-shape fit", "Lorentzian notch, $q$ = 0"], loc="upper left")
    a.set_xlim(-8, 8)
    a.set_ylim(-0.04, 1.45)
    a.set_yticks(np.arange(0.0, 1.41, 0.2))
    a.set_xlabel("Normalised detuning, $2(\\lambda-\\lambda_r)/\\Delta\\lambda$")
    a.set_ylabel("Normalised transmission, $T$")
    cax = style.axes_pt(fig, *lay["colorbar"])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
    cb.set_ticks([float(ns.min()), 0.5 * float(ns.min() + ns.max()), float(ns.max())])
    cb.ax.xaxis.set_major_formatter(style.formatter(3))
    cb.ax.tick_params(labelsize=6.0, length=2.0, width=0.6, pad=1.5)
    cb.outline.set_linewidth(0.5)
    cb.set_label("$n_a$", fontsize=7.0, labelpad=0.5)
    _absorber_inset(fig, inp, lay["inset_a"], cmap, norm)

    # ---- (b) linewidth: harmonic inversion against the two spectrum normalisations
    fwhm_h = np.array([s["fwhm_harminv"] for s in spectra])
    fwhm_fit = np.array([s["fit"]["fwhm"] for s in spectra])
    fwhm_in = np.array([s["fit_input_normalisation"]["fwhm"] for s in spectra])
    b.plot(ns, fwhm_h, "-", color="k", lw=1.0, zorder=2, label="$\\lambda_r/Q$, harmonic inversion")
    b.plot(ns, fwhm_fit, "o", color=style.BLUE, mec=style.BLUE, ms=4.4, ls="none", zorder=3,
           label="spectrum, cavity-free reference", **OPEN)
    b.plot(ns, fwhm_in, "x", color=style.VERMILLION, ms=4.4, mew=1.0, ls="none", zorder=3,
           label="spectrum, input-plane reference")
    b.set_xlim(1.285, 1.465)
    b.set_ylim(0, 1.3)
    b.set_xticks(np.round(np.arange(1.300, 1.4501, 0.025), 3))
    b.set_xlabel(N_LABEL)
    b.set_ylabel("Linewidth, $\\Delta\\lambda$ (nm)")
    b.legend(loc="upper left")
    ins = style.axes_pt(fig, *lay["inset"])
    ins.axhline(0, color="0.7", lw=0.6, zorder=1)
    ins.plot(ns, 100 * (fwhm_fit / fwhm_h - 1), "s-", color=style.BLUE, lw=0.8, ms=2.4, mew=0.6, zorder=2)
    ins.set_xlim(1.285, 1.465)
    ins.set_ylim(-4, 2)
    ins.set_xticks([1.30, 1.45])
    ins.set_yticks([-3, 0])
    ins.tick_params(labelsize=5.6, length=2.0, pad=1.5)
    ins.set_ylabel("deviation (%)", fontsize=6.0, labelpad=1)
    return fig


INSET_N = 1.375            # the analyte index of the inset: the largest asymmetry in the 25a cell


def _absorber_inset(fig, inp: Inputs, box, cmap, norm) -> None:
    """Inset of figure 6(a) (SPRAT 1.2.1): one analyte index under the two guide terminations, each
    spectrum on its own normalised detuning from the fit of the comparison (analysis.json, key termination)."""
    try:
        ab = inp.spectra_absorber
        rows = inp.block("termination", "rows")
        pml = inp.spectra
    except MissingData:
        return
    row = [r for r in rows if abs(float(r["n"]) - INSET_N) < 1e-9 and r.get("status") == "graded"]
    s_ab = [s for s in ab if abs(float(s["n"]) - INSET_N) < 1e-9]
    s_pm = [s for s in pml if abs(float(s["n"]) - INSET_N) < 1e-9]
    if not (row and s_ab and s_pm):
        return
    ins = style.axes_pt(fig, *box)
    ins.axhline(1.0, color="0.75", lw=0.5, ls=":", zorder=1)
    for side, s, col in (("pml", s_pm[0], cmap(norm(INSET_N))), ("absorber", s_ab[0], "k")):
        f = row[0][side]
        lam, T = np.asarray(s["lambda_nm"], float), np.asarray(s["T"], float)
        e = 2 * (lam - f["lambda_r"]) / f["FWHM_nm"]
        m = np.abs(e) <= 6.0
        ins.plot(e[m], T[m], "-", color=col, lw=0.9, zorder=3 if side == "absorber" else 2)
    ins.set_xlim(-6, 6)
    ins.set_ylim(-0.04, 1.36)
    ins.set_xticks([-5, 0, 5])
    ins.set_yticks([0, 1])
    ins.tick_params(labelsize=5.6, length=2.0, pad=1.5)
    ins.text(0.04, 0.95, "%.3f" % INSET_N, transform=ins.transAxes, ha="left", va="top", fontsize=5.6)


# ============================================================================ figure 7
def figure7(inp: Inputs):
    """Calibration (a), local sensitivity (b), quality factor (c), figure of merit and detection limit (d)."""
    sw = inp.block("analyte_sweep")
    n, lam, Q = (np.asarray(sw[k], float) for k in ("n", "lambda_nm", "Q"))
    nc, S, Qc, FOM = (np.asarray(sw[k], float) for k in ("n_c", "S_c", "Q_c", "FOM_c"))
    at = sw["at_1p33"]
    chord = sw["chord"]
    lin = np.polyval([chord["S"], chord["intercept"]], n)
    cubic = np.poly1d(sw["cubic_coefficients_descending"])
    lay = LAYOUT[7]
    fig, (a, b, c, d) = _canvas(7)
    xlim = (1.293, 1.457)
    xt = np.round(np.arange(1.300, 1.4501, 0.025), 3)

    # ---- (a) calibration curve
    nn = np.linspace(n[0], n[-1], 200)
    a.plot(n, lam, "o", color=style.BLUE, mec=style.BLUE, ms=3.4, mfc="white", mew=0.9, ls="none", zorder=3, label="FDTD")
    a.plot(nn, cubic(nn), "-", color="k", lw=1.0, zorder=4, label="cubic fit")
    a.plot(n[[0, -1]], lin[[0, -1]], ":", color=style.GREY, lw=1.0, zorder=2, label="linear fit")
    a.set_ylim(1556, 1660)
    a.set_yticks([1560, 1580, 1600, 1620, 1640])
    a.legend(loc="upper left")
    ins = style.axes_pt(fig, *lay["inset"])
    ins.axhline(0, color="0.7", lw=0.6, zorder=1)
    ins.plot(n, lam - lin, "-", color=style.BLUE, lw=0.9, zorder=2)
    ins.set_ylim(-1.15, 1.0)
    ins.set_yticks([-1, 0, 1])
    ins.set_xticks([1.30, 1.45])
    ins.tick_params(labelsize=5.6, length=2.0, pad=1.5)
    ins.set_ylabel("linear residual (nm)", fontsize=6.0, labelpad=1)

    # ---- (b) local sensitivity
    b.plot(nc, S, "-", color=style.BLUE, lw=1.2, zorder=2, label="$S$ = d$\\lambda_r$/d$n_a$")
    b.axhline(chord["S"], color=style.GREY, lw=1.0, ls=":", zorder=1,
              label=f"linear-fit slope, {chord['S']:.1f} nm/RIU")
    b.plot([at["n"]], [at["S"]], "o", color=style.BLUE, mec=style.BLUE, ms=4.4, mew=1.0, zorder=3)
    b.annotate(f"{at['S']:.1f} nm/RIU", xy=(at["n"], at["S"]), xytext=(7, 3.5), textcoords="offset points",
               ha="left", va="bottom", fontsize=7.0)
    b.set_ylim(570, 672)
    b.set_yticks([580, 600, 620, 640, 660])
    b.legend(loc="lower left")

    # ---- (c) quality factor: harmonic inversion and the transmission spectra
    c.semilogy(nc, Qc, "-", color=style.BLUE, lw=1.2, zorder=2, label="harmonic inversion")
    try:
        sp = inp.spectra
        c.semilogy([s["n"] for s in sp], [s["fit"]["Q"] for s in sp], "s", color=style.VERMILLION, mec=style.VERMILLION,
                   ms=4.4, ls="none", zorder=3, label="transmission spectrum", **OPEN)
    except MissingData:
        pass
    # SPRAT 1.2.1: the reflectionless quality factor of the absorber spectra (amendment A2)
    phys = _physical(inp)
    if phys:
        c.semilogy(phys["nodes"]["n"], phys["nodes"]["Q_w"], "D-", color="k", mec="k", ms=3.4, lw=0.8, zorder=4,
                   label="reflectionless, $Q_{\\mathrm{w}}$")
    c.set_ylim(1200, 12000)
    _loglocator(c, [2000, 4000, 8000], minor=False)
    c.yaxis.set_major_formatter(style.formatter(0))
    c.legend(loc="lower left")

    # ---- (d) figure of merit and detection limit DL = fwhm / (10 S) = 1 / (10 FOM)
    d.semilogy(nc, FOM, "-", color=style.BLUE, lw=1.2, zorder=2, label="25a cell")
    d.plot([at["n"]], [at["FOM"]], "o", color=style.BLUE, mec=style.BLUE, ms=4.4, mew=1.0, zorder=3)
    d.annotate(f"FOM = {at['FOM']:.0f} RIU$^{{-1}}$", xy=(at["n"], at["FOM"]), xytext=(7, 3.5),
               textcoords="offset points", ha="left", va="bottom", fontsize=7.0)
    if phys:
        d.semilogy(phys["curve"]["n"], phys["curve"]["FOM_w"], "--", color="k", lw=1.1, zorder=2, label="reflectionless")
        ref = [r for r in phys["table"] if abs(r["n"] - at["n"]) < 1e-9]
        if ref:
            d.plot([ref[0]["n"]], [ref[0]["FOM_w"]], "o", color="k", mec="k", ms=4.0, mew=1.0, zorder=3)
            d.annotate(f"{ref[0]['FOM_w']:.0f}", xy=(ref[0]["n"], ref[0]["FOM_w"]),
                       xytext=(-4, -2.5), textcoords="offset points", ha="right", va="top", fontsize=7.0)
        d.legend(loc="lower left")
    lo, hi = 450, 4500
    d.set_ylim(lo, hi)
    _loglocator(d, [500, 1000, 2000, 4000], minor=False)
    d.yaxis.set_major_formatter(style.formatter(0))
    d.set_ylabel("Figure of merit, FOM (RIU$^{-1}$)")
    r = d.twinx()
    for side in ("left", "top", "bottom"):
        r.spines[side].set_visible(False)
    r.spines["right"].set_visible(True)
    r.set_yscale("log")
    r.set_ylim(1 / (10 * lo), 1 / (10 * hi))
    r.yaxis.set_major_locator(FixedLocator([2.5e-5, 5e-5, 1e-4, 2e-4]))
    r.yaxis.set_minor_locator(NullLocator())
    r.yaxis.set_major_formatter(SCI)
    r.set_ylabel("Detection limit, DL = $\\Delta\\lambda/(10S)$ (RIU)")

    for ax in (a, b, c, d):
        ax.set_xlim(*xlim)
        ax.set_xticks(xt)
        ax.set_xlabel(N_LABEL)
    a.set_ylabel("Resonant wavelength, $\\lambda_r$ (nm)")
    b.set_ylabel("Sensitivity, $S$ (nm/RIU)")
    c.set_ylabel(Q_LABEL)
    return fig


def _physical(inp: Inputs) -> dict | None:
    """The reflectionless values across the analyte sweep (analysis.json, combined.physical_values), or None."""
    try:
        return inp.block("combined", "physical_values")
    except MissingData:
        return None


__all__ = ["figure2", "figure3", "figure5", "figure6", "figure7", "LAYOUT"]
