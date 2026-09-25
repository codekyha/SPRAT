"""Figure 1: geometry of the side-coupled cavity, drawn from the reference record.

(a) The whole computational cell: PML, analyte, rods, the W1 guide, the line source, the two
flux planes and the zoom window of (b).  (b) The neighbourhood of the cavity with the lattice
constant, the defect radius, the separating rows N_sep and the rod radius.  (c) The definition
of the displacement components of the defect rod, drawn for illustration with dx = 0.10a and
dy = 0.08a (the reference geometry has dx = dy = 0).

The cell, the rod list and the positions of the source and the flux planes come from
``sprat.geometry`` applied to the structure and parameter blocks of the reference record, so
the drawing is the geometry that the FDTD runs used.
"""

from __future__ import annotations

from matplotlib.patches import Circle, Rectangle

from .. import geometry
from ..structure import structure_from_dict
from . import style
from .inputs import Inputs

CANVAS = (391.09, 180.68)                                 # points; the page is then cropped to the content
AXES = {"a": (37.08, 18.27, 135.38, 130.96), "b": (215.59, 26.45, 78.57, 114.59), "c": (332.61, 59.85, 51.28, 47.79)}
ZOOM_HALF_WIDTH = 2.4                                     # (b): x in [-2.4, 2.4]a, y from -1 to N_sep + 2
ILLUSTRATION = dict(dx=0.10, dy=0.08)                     # (c): displacement drawn for illustration
WINDOW = 0.32                                             # half width of the window of (c) (0.64a x 0.64a)
LABEL_BOX = dict(facecolor="white", alpha=0.90, edgecolor="none", pad=0.9)


def _dimension(ax, p0, p1, extension, label, horizontal=True, side=-1, lw=0.6, fontsize=8.0, head=6, gap=1.6):
    """Engineering-drawing dimension line between p0 and p1 with extension lines of length ``extension``
    (units of a) and the label ``gap`` points beyond the dimension line.

    side = -1: the dimension line lies below (horizontal) / right of (vertical) the points; +1: above / left.
    """
    (x0, y0), (x1, y1) = p0, p1
    u = abs(extension)
    arrow = dict(arrowstyle="<|-|>", lw=lw, color="k", shrinkA=0, shrinkB=0, mutation_scale=head)
    if horizontal:
        yk = (min(y0, y1) - u) if side < 0 else (max(y0, y1) + u)
        for xx, yy in ((x0, y0), (x1, y1)):
            ax.plot([xx, xx], [yy, yk], "-", color="k", lw=lw * 0.7, zorder=7)
        ax.annotate("", xy=(x1, yk), xytext=(x0, yk), arrowprops=arrow, zorder=7)
        ax.annotate(label, xy=(0.5 * (x0 + x1), yk), xytext=(0, side * gap), textcoords="offset points", ha="center",
                    va="top" if side < 0 else "bottom", fontsize=fontsize, zorder=8)
    else:
        xk = (max(x0, x1) + u) if side < 0 else (min(x0, x1) - u)
        for xx, yy in ((x0, y0), (x1, y1)):
            ax.plot([xx, xk], [yy, yy], "-", color="k", lw=lw * 0.7, zorder=7)
        ax.annotate("", xy=(xk, y1), xytext=(xk, y0), arrowprops=arrow, zorder=7)
        ax.annotate(label, xy=(xk, 0.5 * (y0 + y1)), xytext=(-side * gap, 0), textcoords="offset points",
                    ha="left" if side < 0 else "right", va="center", fontsize=fontsize, zorder=8)


def _pml(ax, sx, sy, ex, ey):
    """Absorbing layers: ``ex`` thick at the guide ends (PML or absorber), ``ey`` at the top and bottom."""
    for (x, y, w, h) in ((-sx / 2, -sy / 2, sx, ey), (-sx / 2, sy / 2 - ey, sx, ey), (-sx / 2, -sy / 2, ex, sy),
                         (sx / 2 - ex, -sy / 2, ex, sy)):
        ax.add_patch(Rectangle((x, y), w, h, facecolor=style.PML_FILL, edgecolor="none", zorder=0.5))
        ax.add_patch(Rectangle((x, y), w, h, facecolor="none", hatch="////", edgecolor="0.62", lw=0.0, zorder=0.6))


def _rods(ax, rods, defect_center, lw_defect):
    for (x, y, r) in rods:
        is_defect = defect_center is not None and abs(x - defect_center[0]) < 1e-9 and abs(y - defect_center[1]) < 1e-9
        if is_defect:
            ax.add_patch(Circle((x, y), r, facecolor=style.VERMILLION, edgecolor="k", lw=lw_defect, zorder=4))
        else:
            ax.add_patch(Circle((x, y), r, facecolor=style.ROD, edgecolor="none", zorder=3))


def figure1(inp: Inputs):
    rec = inp.reference()
    s, p = structure_from_dict(rec["structure"]), rec["params"]
    g = geometry.build(s, p)
    sx, sy, d = g.sx, g.sy, g.pml
    r0 = float(s["lattice"]["rod_radius"])
    rd, jc = float(s["defect"]["radius"]), int(s["defect"]["row"])
    dx0, dy0 = float(s["defect"]["dx"]), float(s["defect"]["dy"])
    src_x, in_x, out_x = geometry.source_and_flux_x(g)

    fig = style.new_figure(*CANVAS)
    axa, axb, axc = (style.axes_pt(fig, *AXES[k]) for k in "abc")

    # ---------------------------------------------------------------- (a) the computational cell
    axa.add_patch(Rectangle((-sx / 2, -sy / 2), sx, sy, facecolor=style.ANALYTE, edgecolor="k", lw=0.6, zorder=0))
    axa.add_patch(Rectangle((-sx / 2, -0.5), sx, 1.0, facecolor=style.GUIDE, edgecolor="none", zorder=1))
    _pml(axa, sx, sy, g.edge_x, d)
    _rods(axa, g.rods, g.defect_center, 0.4)
    axa.plot([src_x, src_x], [-1.0, 1.0], "-", color=style.BLUE, lw=1.8, solid_capstyle="butt", zorder=5)
    axa.annotate("", xy=(src_x + 1.8, 0), xytext=(src_x + 0.35, 0), zorder=6,
                 arrowprops=dict(arrowstyle="-|>", lw=0.9, color=style.BLUE, mutation_scale=7))
    for x in (in_x, out_x):
        axa.plot([x, x], [-1.6, 1.6], color=style.BLUE, lw=0.9, dashes=(3, 2), zorder=5)
    for x, text, y in ((src_x, "source", -5.5), (in_x, "input flux", 5.5), (out_x, "output flux", 5.5)):
        axa.text(x + 0.30, y, text, rotation=90, ha="left", va="center", fontsize=6.8, zorder=8, bbox=LABEL_BOX)
    axa.text(8.0, sy / 2 - d / 2, "PML", ha="center", va="center", fontsize=6.5, color=style.DARK, zorder=8,
             bbox=dict(facecolor=style.PML_FILL, edgecolor="none", pad=0.6))
    zx0, zx1, zy0, zy1 = -ZOOM_HALF_WIDTH, ZOOM_HALF_WIDTH, -1.0, jc + 2.0
    axa.add_patch(Rectangle((zx0, zy0), zx1 - zx0, zy1 - zy0, facecolor="none", edgecolor="k", lw=0.9,
                            ls=(0, (3.5, 2)), zorder=7))
    axa.set_xlim(-sx / 2 - 0.3, sx / 2 + 0.3)
    axa.set_ylim(-sy / 2 - 0.3, sy / 2 + 0.3)
    axa.set_aspect("equal")
    axa.set_xticks([-10, 0, 10])
    axa.set_yticks([-10, -5, 0, 5, 10])
    axa.set_xlabel("$x/a$")
    axa.set_ylabel("$y/a$")
    style.panel_letter(axa, "(a)")

    # ---------------------------------------------------------------- (b) the cavity region
    axb.add_patch(Rectangle((zx0, zy0), zx1 - zx0, zy1 - zy0, facecolor=style.ANALYTE, edgecolor="k", lw=0.6, zorder=0))
    axb.add_patch(Rectangle((zx0, -0.5), zx1 - zx0, 1.0, facecolor=style.GUIDE, edgecolor="none", zorder=1))
    near = [c for c in g.rods if zx0 - 1 < c[0] < zx1 + 1 and zy0 - 1 < c[1] < zy1 + 1]
    _rods(axb, near, g.defect_center, 1.0)
    axb.add_patch(Circle((0.0, jc), r0, facecolor="none", edgecolor=style.NOTE, lw=0.9, ls=(0, (2.2, 1.6)), zorder=3.5))
    axb.add_patch(Rectangle((dx0 - WINDOW, jc + dy0 - WINDOW), 2 * WINDOW, 2 * WINDOW, facecolor="none", edgecolor="k",
                            lw=0.7, ls=(0, (2, 1.4)), zorder=7))
    _dimension(axb, (1.0, jc + 1.0), (2.0, jc + 1.0), 0.42, "$a$", horizontal=True, side=+1)
    axb.annotate("", xy=(-1.5, float(jc)), xytext=(-1.5, 0.0), zorder=7,
                 arrowprops=dict(arrowstyle="<|-|>", lw=0.6, color="k", shrinkA=0, shrinkB=0, mutation_scale=6))
    axb.text(-1.42, 1.61, "$N_{\\mathrm{sep}}$ = %d" % jc, ha="left", va="center", fontsize=7.0, zorder=8, bbox=LABEL_BOX)
    axb.annotate("$r_d$", xy=(dx0 - 0.72 * rd, jc + dy0 + 0.72 * rd), xytext=(-1.30, jc + 1.19), fontsize=8.0, ha="center",
                 va="center", zorder=8, bbox=LABEL_BOX, arrowprops=dict(arrowstyle="-", lw=0.5, color=style.NOTE))
    axb.annotate("$r$", xy=(2.0 - 0.72 * r0, 2.0 - 0.72 * r0), xytext=(1.45, 1.55), fontsize=8.0, ha="center",
                 va="center", zorder=8, bbox=LABEL_BOX, arrowprops=dict(arrowstyle="-", lw=0.5, color=style.NOTE))
    axb.text(zx1 - 0.12, 0.0, "W1", ha="right", va="center", fontsize=7.5, zorder=8, bbox=LABEL_BOX)
    axb.set_xlim(zx0, zx1)
    axb.set_ylim(zy0, zy1)
    axb.set_aspect("equal")
    axb.set_xticks([-2, 0, 2])
    axb.set_yticks(range(int(zy0), int(zy1) + 1))
    axb.set_xlabel("$x/a$")
    axb.set_ylabel("$y/a$")
    style.panel_letter(axb, "(b)")

    # ---------------------------------------------------------------- (c) the displacement of the defect rod
    ddx, ddy = ILLUSTRATION["dx"], ILLUSTRATION["dy"]
    axc.add_patch(Rectangle((-9, -9), 18, 18, facecolor=style.ANALYTE, edgecolor="none", zorder=0))
    axc.add_patch(Circle((0, 0), r0, facecolor="none", edgecolor=style.NOTE, lw=0.9, ls=(0, (2.2, 1.6)), zorder=3))
    axc.add_patch(Circle((ddx, ddy), rd, facecolor=style.VERMILLION, edgecolor="k", lw=1.0, zorder=4))
    axc.plot(0, 0, "+", color=style.NOTE, ms=5, mew=0.8, zorder=5)
    axc.plot(ddx, ddy, "+", color="k", ms=5, mew=0.8, zorder=5)
    _dimension(axc, (0.0, 0.0), (ddx, ddy), 0.20, "$\\delta_x$", horizontal=True, side=-1)
    _dimension(axc, (0.0, 0.0), (ddx, ddy), 0.16, "$\\delta_y$", horizontal=False, side=-1)
    axc.set_xlim(-WINDOW - 0.02, WINDOW + 0.22)
    axc.set_ylim(-WINDOW - 0.16, WINDOW + 0.02)
    axc.set_aspect("equal")
    axc.set_xticks([])
    axc.set_yticks([])
    for sp in axc.spines.values():
        sp.set_visible(True)
    axc.text(0.5, -0.045, "%s$a$ \u00d7 %s$a$" % (style.num(2 * WINDOW, 2), style.num(2 * WINDOW, 2)),
             transform=axc.transAxes, ha="center", va="top", fontsize=7.0)
    style.panel_letter(axc, "(c)")
    return fig


__all__ = ["figure1"]
