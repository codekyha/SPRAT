"""What the figures read: the records directory and the analysis outputs (``tables``).

The analysis outputs are read from ``analysis.json`` and ``spectra.json`` as written by
``sprat analyze``; the few quantities that the analysis does not export (the geometry of the
reference cell, the band records, the complex-band curve of the plane-wave record) are taken
from the records themselves.  Every accessor raises :class:`MissingData` with a plain message
when the records do not support a figure, so that ``sprat figures`` can skip it and say why.
"""

from __future__ import annotations

import json
import os
from functools import cached_property
from typing import Any

from ..analysis import collect

# the reference geometry of the paper (same parameters as sprat.analysis.sensitivity.REFERENCE)
REFERENCE = dict(radius=0.060, row=4, dx=0.0, dy=0.0, resolution=24, cladding_rows=12, guide_periods=25,
                 termination_x="pml", pml=1.0, no_rows=True)
REFERENCE_N = 1.33


class MissingData(LookupError):
    """The records (or the analysis outputs) do not hold what a figure needs."""


def key(d: dict, k) -> Any:
    """``d[k]`` for a mapping read back from JSON, where numeric keys became strings ('4', '0.05', '0.0')."""
    if k in d:
        return d[k]
    try:
        kf = float(k)
    except (TypeError, ValueError):
        raise KeyError(k) from None
    for kk, v in d.items():
        try:
            if abs(float(kk) - kf) < 1e-9:
                return v
        except (TypeError, ValueError):
            continue
    raise KeyError(k)


class Inputs:
    """Lazy access to the records and to the analysis outputs of one records directory."""

    def __init__(self, records_dir: str, tables_dir: str | None):
        self.records_dir = records_dir
        self.tables_dir = tables_dir

    # ------------------------------------------------------------------ files
    @cached_property
    def records(self) -> list[dict]:
        recs = collect.load(self.records_dir)
        if not recs:
            raise MissingData(f"no records in {self.records_dir}")
        return recs

    def _table(self, name: str):
        if not self.tables_dir:
            raise MissingData("no analysis outputs (run `sprat analyze` or pass --tables)")
        path = os.path.join(self.tables_dir, name)
        if not os.path.isfile(path):
            raise MissingData(f"{path} is missing (run `sprat analyze`)")
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)

    @cached_property
    def analysis(self) -> dict:
        return self._table("analysis.json")

    @cached_property
    def spectra(self) -> list[dict]:
        """The analysed spectra (``spectra.json``); recomputed from the records if the file is absent."""
        try:
            sp = self._table("spectra.json")
        except MissingData:
            from ..analysis import fano
            sp = fano.analyse_spectra(self.records)
        if not sp:
            raise MissingData("no transmission spectra among the records")
        return sorted(sp, key=lambda s: s["n"])

    @cached_property
    def spectra_absorber(self) -> list[dict]:
        """The absorber-terminated spectra (``spectra_absorber.json``, SPRAT 1.2.1); recomputed from the records if the
        file is absent."""
        try:
            sp = self._table("spectra_absorber.json")
        except MissingData:
            from ..analysis import termination
            sp = termination.absorber_spectra(self.records)
        if not sp:
            raise MissingData("no absorber-terminated spectra among the records")
        return sorted(sp, key=lambda s: s["n"])

    # ------------------------------------------------------------------ analysis blocks
    def block(self, *path: str) -> Any:
        """A nested block of ``analysis.json``, e.g. ``block('barrier', 'leak')``."""
        cur: Any = self.analysis
        for i, p in enumerate(path):
            if not isinstance(cur, dict) or p not in cur:
                reason = (self.analysis.get("errors") or {}).get(path[0], "") if i == 0 else ""
                where = ".".join(path[: i + 1])
                raise MissingData(f"analysis block '{where}' is missing" + (f" ({reason})" if reason else ""))
            cur = cur[p]
        if cur is None:
            raise MissingData(f"analysis block '{'.'.join(path)}' is empty")
        return cur

    # ------------------------------------------------------------------ records
    def reference(self) -> dict:
        """The harmonic-inversion record of the reference geometry at n_a = 1.33 (longest signal)."""
        r = collect.best(self.records, mode="harminv", n=REFERENCE_N, **REFERENCE)
        if r is None:
            raise MissingData("no harmonic-inversion record of the reference geometry (r_d = 0.060a, N_sep = 4, "
                              "n_cl = 12, resolution 24, 25a guide, PML) at n_a = 1.33")
        return r

    @cached_property
    def pwe(self) -> dict:
        """Result block of the newest plane-wave record (the choice of ``sprat analyze``)."""
        cands = [r for r in self.records if r["params"]["run"]["mode"] == "pwe"]
        if not cands:
            raise MissingData("no plane-wave record (run a task with mode = pwe)")
        cands.sort(key=lambda r: (r["provenance"].get("started") or "", r["task"]["label"]))
        return cands[-1]["result"]

    def bands(self, task: str) -> dict:
        """Result block of the newest bands record of a task ('bulk', 'sweep', 'w1', 'calibration')."""
        cands = [r for r in self.records if r["params"]["run"]["mode"] == "bands" and r["result"].get("task") == task]
        if not cands:
            raise MissingData(f"no bands record with task = {task}")
        cands.sort(key=lambda r: (r["provenance"].get("started") or "", r["task"]["label"]))
        return cands[-1]["result"]


__all__ = ["Inputs", "MissingData", "REFERENCE", "REFERENCE_N", "key"]
