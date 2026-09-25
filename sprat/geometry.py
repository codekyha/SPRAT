"""From a structure and a parameter block to the computational cell and the rod list.

Everything here is plain Python (no Meep) so that the geometry can be tested and drawn
without the FDTD engine.  ``sprat.fdtd`` converts the result to Meep objects.

Cell size (units of a):

* PML termination:      s_x = n_x + 2 pad_x + 2 d_PML,
* absorber termination: s_x = n_x + 2 n_abs  (the lattice and the guide continue into the absorber),
* both:                 s_y = (2 n_cl + 1) + 2 pad_y + 2 d_PML.

The PML (or absorber) lies inside the cell, as in Meep.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class Geometry:
    sx: float
    sy: float
    rods: list[tuple[float, float, float]]          # (x, y, radius), units of a
    termination_x: str                                # 'pml' | 'absorber'
    pml: float                                        # PML thickness
    absorber: float                                   # absorber thickness (0 with pml termination)
    defect_center: tuple[float, float] | None         # (x, y) of the defect rod centre, None without a cavity
    mirror_x: bool                                    # rod list symmetric under x -> -x
    columns: int = 0
    notes: list[str] = field(default_factory=list)

    @property
    def edge_x(self) -> float:
        """Thickness of the absorbing layer along x (PML or absorber)."""
        return self.absorber if self.absorber > 0 else self.pml

    @property
    def interior(self) -> tuple[float, float]:
        """Size of the cell without the absorbing layers, (L_x, L_y)."""
        return self.sx - 2 * self.edge_x, self.sy - 2 * self.pml


def cell_size(structure: dict, pml: float) -> tuple[float, float]:
    c = structure["cell"]
    if c["termination_x"] == "absorber":
        sx = c["guide_periods"] + 2.0 * c["absorber_periods"]
    else:
        sx = c["guide_periods"] + 2.0 * c["pad_x"] + 2.0 * pml
    sy = (2 * c["cladding_rows"] + 1) + 2.0 * c["pad_y"] + 2.0 * pml
    return float(sx), float(sy)


def build(structure: dict, params: dict, with_defect: bool | None = None) -> Geometry:
    """Build the rod list.

    ``with_defect`` defaults to ``mode != 'reference'``: the cavity-less reference run keeps
    the full lattice.
    """
    c, d, la = structure["cell"], structure["defect"], structure["lattice"]
    pml = float(params["numerics"]["pml"])
    mode = params["run"]["mode"]
    if with_defect is None:
        with_defect = mode != "reference"
    sx, sy = cell_size(structure, pml)
    absorber = float(c["absorber_periods"]) if c["termination_x"] == "absorber" else 0.0
    n_col = int(round(sx)) if absorber > 0 else int(c["guide_periods"])
    half = (n_col - 1) / 2.0
    ncl = int(c["cladding_rows"])
    r0 = float(la["rod_radius"])
    rows = structure.get("rows", {})
    guide = structure["waveguide"]["type"] == "w1"
    rods: list[tuple[float, float, float]] = []
    for i in range(n_col):
        x = -half + i
        for j in range(-ncl, ncl + 1):
            if j == 0 and guide:
                continue
            rods.append((float(x), float(j), float(rows.get(j, r0))))
    center = None
    notes: list[str] = []
    if with_defect and d["type"] == "rod":
        jc = int(d["row"])
        before = len(rods)
        rods = [r for r in rods if not (abs(r[0]) < 1e-9 and abs(r[1] - jc) < 1e-9)]
        if len(rods) == before:
            notes.append(f"no lattice rod at (0, {jc}) to replace")
        if d["radius"] > 1e-6:
            rods.append((float(d["dx"]), float(jc + d["dy"]), float(d["radius"])))
        center = (float(d["dx"]), float(jc + d["dy"]))
    for e in structure.get("rods", []):
        if e["op"] == "add":
            rods.append((float(e["x"]), float(e["y"]), float(e["radius"])))
        else:
            before = len(rods)
            rods = [r for r in rods if not (abs(r[0] - e["x"]) < 1e-9 and abs(r[1] - e["y"]) < 1e-9)]
            if len(rods) == before:
                notes.append(f"remove {e['x']} {e['y']}: no rod there")
    mirror = _mirror_symmetric(rods)
    return Geometry(sx=sx, sy=sy, rods=rods, termination_x=c["termination_x"], pml=pml, absorber=absorber,
                    defect_center=center, mirror_x=mirror, columns=n_col, notes=notes)


def _mirror_symmetric(rods: list[tuple[float, float, float]], tol: float = 1e-9) -> bool:
    keyset = {(round(x, 6), round(y, 6), round(r, 6)) for x, y, r in rods}
    return all((round(-x, 6), round(y, 6), round(r, 6)) in keyset for x, y, r in rods)


def symmetry_sector(structure: dict, params: dict, geometry: Geometry) -> tuple[str, int]:
    """Return ('none' | 'even' | 'odd', phase) following the rule of the original code.

    'auto' gives the even mirror sector when the defect displacement dx is zero and the rod
    list is mirror-symmetric; an explicit 'even' or 'odd' is honoured only when dx is zero.
    """
    s = params["numerics"]["symmetry"]
    dx = float(structure["defect"]["dx"]) if structure["defect"]["type"] == "rod" else 0.0
    if s == "auto":
        s = "even" if (abs(dx) < 1e-12 and geometry.mirror_x) else "none"
    if s == "none" or abs(dx) > 1e-12 or not geometry.mirror_x:
        return "none", 0
    return s, (1 if s == "even" else -1)


def probe_point(structure: dict, params: dict) -> tuple[float, float]:
    """Harminv probe / source position: the defect centre plus the probe offset."""
    d = structure["defect"]
    cx = float(d["dx"]) if d["type"] == "rod" else 0.0
    cy = float(d["row"] + d["dy"]) if d["type"] == "rod" else float(d["row"])
    return cx + float(params["harminv"]["probe_dx"]), cy + float(params["harminv"]["probe_dy"])


def source_and_flux_x(geometry: Geometry) -> tuple[float, float, float]:
    """x positions of the line source, the input flux plane and the output flux plane."""
    edge = geometry.edge_x
    source_x = -geometry.sx / 2 + edge + 0.5
    return source_x, source_x + 1.0, geometry.sx / 2 - edge - 0.5


def describe(geometry: Geometry) -> dict[str, Any]:
    return dict(sx=geometry.sx, sy=geometry.sy, rods=len(geometry.rods), columns=geometry.columns,
                termination_x=geometry.termination_x, pml=geometry.pml, absorber=geometry.absorber,
                mirror_x=geometry.mirror_x, interior=geometry.interior, notes=geometry.notes)


__all__ = ["Geometry", "cell_size", "build", "symmetry_sector", "probe_point", "source_and_flux_x", "describe"]
