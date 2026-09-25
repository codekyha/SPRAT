"""House style of the seven figures of the paper.

Type is Arial where it is installed and otherwise Liberation Sans, which has the same metrics;
mathtext uses the same family.  Colours come from the Okabe-Ito palette.  Axes are open (no top
or right spine) with outward ticks, panel letters are bold and sit outside the axes at the top
left, fonts are embedded as TrueType (``pdf.fonttype = 42``) and the figures are 16 cm wide.

Two size schemes exist: ``campaign`` for figures 1-3 and 5-7 (7 pt tick labels, 8 pt axis labels,
9 pt panel letters) and ``mechanism`` for figure 4 (the sizes of its own script, which is
included at ``\\textwidth``).

Figures are built on :class:`matplotlib.figure.Figure` directly (no pyplot, no global backend
change) and saved by :func:`save` as PDF and PNG.  The PDF carries the title "Figure N" and the
author and no creator, producer or creation date, so that two builds from the same records give
identical files.
"""

from __future__ import annotations

import contextlib
import os
import warnings
from functools import lru_cache

import matplotlib
from matplotlib.figure import Figure
from matplotlib.ticker import FuncFormatter

CM = 1 / 2.54                   # inches per centimetre
PT = 1 / 72                     # inches per point
WIDTH_CM = 16.0                 # width of every figure
AUTHOR = "Hasan Oguz"

# Okabe-Ito palette
BLACK = "#000000"
ORANGE = "#E69F00"
SKY = "#56B4E9"
GREEN = "#009E73"
YELLOW = "#F0E442"
BLUE = "#0072B2"
VERMILLION = "#D55E00"
PURPLE = "#CC79A7"
OKABE_ITO = (BLACK, ORANGE, SKY, GREEN, YELLOW, BLUE, VERMILLION, PURPLE)

# neutral tones and the fills of the schematic and the gap plots
GREY = "0.55"                   # secondary lines (fits, reference lines)
DARK = "0.25"                   # annotation text of figure 4
NOTE = "0.3"                    # annotation text and lines of figures 1-3, 5-7
ROD = "#4A5568"                 # dielectric rods
ANALYTE = "#E6F2FA"             # analyte background
GUIDE = "#CFE6F5"               # W1 guide row
GAP = "#DCEBF5"                 # TM band gap
PML_FILL = "0.88"               # perfectly matched layer
BAND = "0.5"                    # shaded ranges, drawn with an alpha

LETTER = dict(x=-0.16, y=1.02, size=9.0)            # panel letters of the campaign scheme
LETTER_MECHANISM = dict(x=-0.13, y=1.06, size=10.5)  # panel letters of figure 4


@lru_cache(maxsize=1)
def font_family() -> str:
    """Arial if installed, otherwise Liberation Sans (metrically identical); DejaVu Sans only as a last resort."""
    from matplotlib import font_manager
    names = {f.name for f in font_manager.fontManager.ttflist}
    for fam in ("Arial", "Liberation Sans"):
        if fam in names:
            return fam
    warnings.warn("neither Arial nor Liberation Sans is installed; the figures fall back to DejaVu Sans "
                  "(install fonts-liberation for the house style)", RuntimeWarning, stacklevel=2)
    return "DejaVu Sans"


def rc(scheme: str = "campaign") -> dict:
    """rcParams of the house style; ``scheme`` is 'campaign' (figures 1-3, 5-7) or 'mechanism' (figure 4)."""
    fam = font_family()
    base = {
        "font.family": "sans-serif", "font.sans-serif": [fam, "Liberation Sans", "DejaVu Sans"],
        "mathtext.fontset": "custom", "mathtext.rm": fam, "mathtext.it": fam + ":italic",
        "mathtext.bf": fam + ":bold", "mathtext.sf": fam, "mathtext.cal": fam,
        "axes.unicode_minus": True, "legend.frameon": False,
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.direction": "out", "ytick.direction": "out",
        "pdf.fonttype": 42, "ps.fonttype": 42, "figure.dpi": 150, "savefig.bbox": "tight",
    }
    if scheme == "campaign":
        base.update({
            "font.size": 7.0, "axes.labelsize": 8.0, "xtick.labelsize": 7.0, "ytick.labelsize": 7.0,
            "legend.fontsize": 6.8, "legend.title_fontsize": 6.8, "legend.handlelength": 1.8,
            "legend.labelspacing": 0.25, "legend.borderaxespad": 0.3, "axes.linewidth": 0.6,
            "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.minor.size": 2.0, "ytick.minor.size": 2.0,
            "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.minor.width": 0.6, "ytick.minor.width": 0.6,
            "lines.linewidth": 1.0, "lines.markersize": 4.4, "lines.markeredgewidth": 1.0,
            "hatch.linewidth": 0.8, "savefig.dpi": 300, "savefig.pad_inches": 0.1,
        })
    elif scheme == "mechanism":
        base.update({
            "font.size": 8.6, "axes.labelsize": 8.8, "xtick.labelsize": 8.0, "ytick.labelsize": 8.0,
            "legend.fontsize": 7.4, "axes.linewidth": 0.7,
            "xtick.major.width": 0.7, "ytick.major.width": 0.7, "xtick.minor.width": 0.5, "ytick.minor.width": 0.5,
            "lines.linewidth": 1.2, "lines.markersize": 4.2, "savefig.dpi": 300, "savefig.pad_inches": 0.02,
        })
    else:
        raise ValueError(f"unknown style scheme {scheme!r}")
    return base


@contextlib.contextmanager
def context(scheme: str = "campaign"):
    """Apply the house style for the duration of a figure build (creation and saving)."""
    with matplotlib.rc_context(rc(scheme)):
        yield


# --------------------------------------------------------------------------- figure and axes
def new_figure(width_pt: float, height_pt: float) -> Figure:
    """A figure of the given canvas size in points (no pyplot state involved)."""
    return Figure(figsize=(width_pt * PT, height_pt * PT))


def axes_pt(fig: Figure, left: float, top: float, width: float, height: float, **kw):
    """Add axes placed in points from the top-left corner of the canvas.

    The layout of every campaign figure is given this way, so the axes keep their size whatever
    the length of the labels; the saved file is then cropped to its content (``bbox = tight``).
    """
    W, H = (v * 72 for v in fig.get_size_inches())
    return fig.add_axes((left / W, 1 - (top + height) / H, width / W, height / H), **kw)


def panel_letter(ax, text: str, x: float | None = None, y: float | None = None, size: float | None = None) -> None:
    """Bold panel letter outside the axes, at the top left (axes-fraction position)."""
    ax.text(LETTER["x"] if x is None else x, LETTER["y"] if y is None else y, text, transform=ax.transAxes,
            fontweight="bold", fontsize=LETTER["size"] if size is None else size, ha="left", va="bottom")


def open_axes(ax) -> None:
    """Hide the top and right spines (the rcParams do this for new axes; twin axes need it explicitly)."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# --------------------------------------------------------------------------- numbers
def num(x: float, d: int | None = None, minus: str = "\u2212") -> str:
    """Format a number for a label: ``d`` decimals (default 0 above 100, 1 above 10, else 3);
    a space groups thousands only from five digits on (12 345, but 6927); ``minus`` is the sign
    used for negative numbers (the typographic minus by default, as on the tick labels)."""
    if d is None:
        d = 0 if abs(x) >= 100 else (1 if abs(x) >= 10 else 3)
    s = f"{abs(x):,.{d}f}"
    whole, _, frac = s.partition(".")
    digits = whole.replace(",", "")
    whole = whole.replace(",", " ") if len(digits) >= 5 else digits
    out = f"{whole}.{frac}" if frac else whole
    return (minus + out) if x < 0 and float(out.replace(" ", "")) != 0 else out


def formatter(d: int | None = None, minus: str = "\u2212") -> FuncFormatter:
    return FuncFormatter(lambda x, pos: num(x, d, minus))


def pct(x: float, d: int = 0) -> str:
    """Signed percentage with the typographic minus, e.g. -31%."""
    return num(x, d) + "%"


# --------------------------------------------------------------------------- saving
def save(fig: Figure, out_dir: str, number: int, dpi: int = 300) -> list[str]:
    """Write ``figureN.pdf`` and ``figureN.png``; return both paths."""
    os.makedirs(out_dir, exist_ok=True)
    stem = os.path.join(out_dir, f"figure{number}")
    fig.savefig(stem + ".pdf", metadata={"Title": f"Figure {number}", "Author": AUTHOR, "Creator": None,
                                          "Producer": None, "CreationDate": None})
    fig.savefig(stem + ".png", dpi=dpi, metadata={"Software": None})
    return [stem + ".pdf", stem + ".png"]


__all__ = ["CM", "PT", "WIDTH_CM", "AUTHOR", "BLACK", "ORANGE", "SKY", "GREEN", "YELLOW", "BLUE", "VERMILLION",
           "PURPLE", "OKABE_ITO", "GREY", "DARK", "NOTE", "ROD", "ANALYTE", "GUIDE", "GAP", "PML_FILL", "BAND",
           "LETTER", "LETTER_MECHANISM", "font_family", "rc", "context", "new_figure", "axes_pt", "panel_letter",
           "open_axes", "num", "formatter", "pct", "save"]
