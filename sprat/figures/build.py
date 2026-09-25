"""``sprat figures``: build the seven figures of the paper from a records directory.

    build_all(records_dir, out_dir, which=None, tables=None, log=print) -> {"figure1": [pdf, png], ...}

``tables`` is the directory of the analysis outputs (``analysis.json``, ``spectra.json``); by
default ``<out_dir>/../tables``.  When it holds no ``analysis.json``, the analysis is run first
(``sprat.analysis.pipeline.analyze``) and writes its outputs there.  ``which`` selects figures:
``None`` or ``"all"``, a comma-separated string (``"1,4,7"`` or ``"figure4"``) or a sequence of
numbers.  A figure whose data the records do not hold is skipped with a message; the returned
mapping lists the figures that were written.
"""

from __future__ import annotations

import os
import time
from typing import Iterable

from . import campaign, mechanism, schematic, style
from .inputs import Inputs, MissingData

FIGURES = {
    1: ("geometry of the cavity and the guide", schematic.figure1, "campaign"),
    2: ("TM gap edges against the analyte index", campaign.figure2, "campaign"),
    3: ("resonance and Q against defect radius, rows and cladding", campaign.figure3, "campaign"),
    4: ("barrier mechanism and finite-cell correction", mechanism.figure4, "mechanism"),
    5: ("displacement of the defect rod", campaign.figure5, "campaign"),
    6: ("transmission spectra and linewidth check", campaign.figure6, "campaign"),
    7: ("refractive-index sensing performance", campaign.figure7, "campaign"),
}


def parse_which(which: str | int | Iterable | None) -> list[int]:
    """Figure numbers from None / 'all' / '1,3,7' / 'figure4' / an int / a sequence of those."""
    if which is None:
        return sorted(FIGURES)
    if isinstance(which, (int, str)):
        items = [which] if isinstance(which, int) else [w for w in which.replace(" ", "").split(",") if w]
    else:
        items = list(which)
    out: list[int] = []
    for w in items:
        s = str(w).strip().lower()
        if s == "all":
            return sorted(FIGURES)
        s = s.removeprefix("figure").removeprefix("fig")
        try:
            n = int(s)
        except ValueError:
            raise ValueError(f"--which: {w!r} is not a figure number (1-7, comma-separated, or all)") from None
        if n not in FIGURES:
            raise ValueError(f"--which: there is no figure {n} (the paper has figures 1-7)")
        if n not in out:
            out.append(n)
    if not out:
        raise ValueError("--which: no figure selected")
    return sorted(out)


def default_tables(out_dir: str) -> str:
    return os.path.join(os.path.dirname(os.path.abspath(out_dir)), "tables")


def ensure_tables(records_dir: str, tables_dir: str, log=print) -> str:
    """Run the analysis into ``tables_dir`` unless its ``analysis.json`` already exists."""
    if not os.path.isfile(os.path.join(tables_dir, "analysis.json")):
        from ..analysis import pipeline
        log(f"no analysis outputs in {tables_dir}: running the analysis first")
        pipeline.analyze(records_dir, tables_dir, log=log)
    return tables_dir


def build_one(number: int, inp: Inputs, out_dir: str, log=print, dpi: int = 300) -> list[str]:
    title, fn, scheme = FIGURES[number]
    with style.context(scheme):
        fig = fn(inp, log=log) if number == 4 else fn(inp)
        files = style.save(fig, out_dir, number, dpi=dpi)
    return files


def build_all(records_dir: str, out_dir: str, which=None, tables: str | None = None, log=print,
              dpi: int = 300) -> dict[str, list[str]]:
    """Build the selected figures; return ``{"figureN": [pdf, png]}`` for those written."""
    numbers = parse_which(which)
    if not os.path.isdir(records_dir):
        raise FileNotFoundError(records_dir)
    tables_dir = tables or default_tables(out_dir)
    if any(n != 1 for n in numbers):
        ensure_tables(records_dir, tables_dir, log=log)
    inp = Inputs(records_dir, tables_dir)
    log(f"figures from {records_dir} (analysis outputs {tables_dir}), type {style.font_family()}")
    done: dict[str, list[str]] = {}
    skipped: dict[str, str] = {}
    for n in numbers:
        t0 = time.time()
        try:
            files = build_one(n, inp, out_dir, log=log, dpi=dpi)
        except MissingData as e:
            skipped[f"figure{n}"] = str(e)
            log(f"figure {n} ({FIGURES[n][0]}): skipped, {e}")
            continue
        done[f"figure{n}"] = files
        log(f"figure {n} ({FIGURES[n][0]}): {', '.join(os.path.basename(f) for f in files)} ({time.time() - t0:.1f} s)")
    log(f"wrote {len(done)} of {len(numbers)} figures to {out_dir}" + (f"; skipped {', '.join(skipped)}" if skipped else ""))
    return done


__all__ = ["FIGURES", "parse_which", "default_tables", "ensure_tables", "build_one", "build_all"]
