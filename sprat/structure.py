"""The structure file (``.phc``): lattice, cell, waveguide, defect, row overrides, free rods.

Units are the lattice constant ``a`` unless a key says otherwise.  Row indices are counted
from the guide row (``j = 0``); positive rows lie on the cavity side.

Example (the reference geometry of the paper)::

    [lattice]
    type = square
    a_nm = 481.4
    rod_radius = 0.20
    rod_eps = 11.9025
    background = analyte

    [cell]
    guide_periods = 25
    cladding_rows = 12
    pad_x = 1.5
    pad_y = 1.0
    termination_x = pml
    absorber_periods = 8

    [waveguide]
    type = W1

    [defect]
    type = rod
    radius = 0.060
    row = 4
    dx = 0.0
    dy = 0.0

    [rows]
    # 2 = 0.18          row 2 gets radius 0.18a (the barrier-row knob)

    [rods]
    # add 3 -2 0.20     extra rod at (x, y) with radius r
    # remove 0 4        delete the rod at (x, y)
"""

from __future__ import annotations

import copy
import json
import math
from typing import Any

from ._textfile import Field, Section, TextFileError, apply_schema, parse_text, schema_markdown

STRUCTURE_SCHEMA: dict[str, Section] = {
    "lattice": Section(
        "The periodic lattice and its materials.",
        {
            "type": Field("enum", "square", "lattice type (only square is implemented in 1.x)", choices=("square",)),
            "a_nm": Field("float_or_calibrate", 481.4,
                          "lattice constant in nanometres, or 'calibrate' to take f_mid * target_wavelength_nm "
                          "from the most recent bulk band record", unit="nm", minimum=1.0),
            "target_wavelength_nm": Field("float", 1550.0, "wavelength that the band-gap centre is calibrated to",
                                          unit="nm", minimum=1.0),
            "rod_radius": Field("float", 0.20, "rod radius r/a", unit="a", minimum=0.0, maximum=0.5),
            "rod_eps": Field("float", 11.9025, "rod permittivity (silicon: 3.45^2)", minimum=1.0),
            "background": Field("enum", "analyte", "medium between the rods: the analyte of the parameter file",
                                choices=("analyte",)),
        }),
    "cell": Section(
        "The finite computational cell.",
        {
            "guide_periods": Field("int", 25, "number of lattice columns along the guide (n_x)", minimum=3),
            "cladding_rows": Field("int", 12, "rows of rods on each side of the guide (n_cl)", minimum=1),
            "pad_x": Field("float", 1.5, "free analyte between the last column and the PML (PML termination only)",
                           unit="a", minimum=0.0),
            "pad_y": Field("float", 1.0, "free analyte between the last row and the PML", unit="a", minimum=0.0),
            "termination_x": Field("enum", "pml", "guide termination along x: 'pml' (padding + PML) or 'absorber' "
                                   "(lattice and guide continue into an adiabatic absorber)",
                                   choices=("pml", "absorber")),
            "absorber_periods": Field("float", 8.0, "thickness of the adiabatic absorber (absorber termination only)",
                                      unit="a", minimum=1.0),
        }),
    "waveguide": Section(
        "The line defect.",
        {"type": Field("enum", "w1", "W1: the row j = 0 is removed; none: full lattice", choices=("w1", "none"))}),
    "defect": Section(
        "The point defect (the cavity).",
        {
            "type": Field("enum", "rod", "rod: the lattice rod at (0, row) is replaced by a rod of the given radius; "
                          "none: no cavity", choices=("rod", "none")),
            "radius": Field("float", 0.060, "defect rod radius r_d/a (0 removes the rod entirely)", unit="a",
                            minimum=0.0, maximum=0.5),
            "row": Field("int", 4, "row of the defect counted from the guide (N_sep)", minimum=1),
            "dx": Field("float", 0.0, "displacement of the defect rod along the guide", unit="a"),
            "dy": Field("float", 0.0, "displacement of the defect rod across the guide", unit="a"),
        }),
    "rows": Section("Per-row radius overrides: `<row index> = <radius/a>`; rows are counted from the guide, "
                    "positive on the cavity side, negative on the other side.", free_keys=True, free_type="float"),
    "rods": Section("Free-form edits applied after everything else, one per line: `add x y r` or `remove x y` "
                    "(units of a).", bare_lines=True),
}


def parse_structure_text(text: str, path: str = "<structure>") -> dict[str, Any]:
    """Parse and validate a structure file; return the resolved nested dict."""
    raw = apply_schema(parse_text(text, path), STRUCTURE_SCHEMA, path)
    s: dict[str, Any] = {k: dict(v) for k, v in raw.items()}
    # rows: keys are row indices
    rows: dict[int, float] = {}
    for key, value in raw["rows"].items():
        if key == "_lines":
            continue
        try:
            j = int(key)
        except ValueError:
            raise TextFileError(f"{path}: [rows] key {key!r} must be an integer row index") from None
        if j == 0:
            raise TextFileError(f"{path}: [rows] row 0 is the guide row and has no rods")
        if not (0.0 <= float(value) <= 0.5):
            raise TextFileError(f"{path}: [rows] radius of row {j} must lie in [0, 0.5]")
        rows[j] = float(value)
    s["rows"] = rows
    # rods: bare lines
    edits: list[dict[str, Any]] = []
    for line in raw["rods"].get("_lines", []):
        parts = line.split()
        try:
            if parts[0].lower() == "add" and len(parts) == 4:
                edits.append({"op": "add", "x": float(parts[1]), "y": float(parts[2]), "radius": float(parts[3])})
            elif parts[0].lower() == "remove" and len(parts) == 3:
                edits.append({"op": "remove", "x": float(parts[1]), "y": float(parts[2])})
            else:
                raise ValueError
        except (ValueError, IndexError):
            raise TextFileError(f"{path}: [rods] line {line!r} must be 'add x y r' or 'remove x y'") from None
    s["rods"] = edits
    validate_structure(s, path)
    return s


def validate_structure(s: dict[str, Any], path: str = "<structure>") -> None:
    cell, defect = s["cell"], s["defect"]
    if defect["type"] == "rod" and defect["row"] > cell["cladding_rows"]:
        raise TextFileError(f"{path}: defect row {defect['row']} lies beyond cladding_rows = {cell['cladding_rows']}")
    for j in s["rows"]:
        if abs(j) > cell["cladding_rows"]:
            raise TextFileError(f"{path}: [rows] row {j} lies beyond cladding_rows = {cell['cladding_rows']}")
    if cell["termination_x"] == "absorber" and abs(cell["absorber_periods"] - round(cell["absorber_periods"])) > 1e-9:
        # allowed, but the column count is then rounded; say so
        pass


def load_structure(path: str) -> dict[str, Any]:
    with open(path, encoding="utf-8") as fh:
        return parse_structure_text(fh.read(), path)


def structure_from_dict(d: dict[str, Any]) -> dict[str, Any]:
    """Normalise a structure dict read back from JSON (rows keys become ints)."""
    s = copy.deepcopy(d)
    s["rows"] = {int(k): float(v) for k, v in (s.get("rows") or {}).items()}
    s.setdefault("rods", [])
    return s


def structure_to_json(s: dict[str, Any]) -> dict[str, Any]:
    d = copy.deepcopy(s)
    d["rows"] = {str(k): v for k, v in sorted(s["rows"].items())}
    return d


def structure_signature(s: dict[str, Any]) -> str:
    """A short deterministic text for labels: the knobs that change the geometry."""
    c, d = s["cell"], s["defect"]
    parts = []
    if d["type"] == "rod":
        parts.append(f"rd{d['radius']:.4f}_dx{d['dx']:+.3f}_dy{d['dy']:+.3f}_row{d['row']}")
    else:
        parts.append("nodefect")
    parts.append(f"cl{c['cladding_rows']}_nx{c['guide_periods']}")
    if c["termination_x"] == "absorber":
        parts.append(f"abs{c['absorber_periods']:g}")
    else:
        parts.append("pml")
    if s["rows"]:
        parts.append("rows" + "_".join(f"{j}-{r:g}" for j, r in sorted(s["rows"].items())))
    if s["rods"]:
        parts.append(f"edits{len(s['rods'])}")
    if s["waveguide"]["type"] == "none":
        parts.append("noguide")
    if abs(s["lattice"]["rod_radius"] - 0.20) > 1e-12:
        parts.append(f"r{s['lattice']['rod_radius']:g}")
    return "_".join(parts)


def structure_docs() -> str:
    return schema_markdown(STRUCTURE_SCHEMA, "Structure file (.phc) reference")


def dumps(s: dict[str, Any]) -> str:
    return json.dumps(structure_to_json(s), indent=1)


def structure_text(s: dict[str, Any], comment: str = "") -> str:
    """Write a structure dict back as a .phc file (used by the campaign builder)."""
    L = []
    if comment:
        L += [f"# {line}" for line in comment.splitlines()]
    la, c, d = s["lattice"], s["cell"], s["defect"]
    L += ["[lattice]", f"type = {la['type']}", f"a_nm = {la['a_nm']}",
          f"target_wavelength_nm = {la['target_wavelength_nm']}", f"rod_radius = {la['rod_radius']}",
          f"rod_eps = {la['rod_eps']}", f"background = {la['background']}", "",
          "[cell]", f"guide_periods = {c['guide_periods']}", f"cladding_rows = {c['cladding_rows']}",
          f"pad_x = {c['pad_x']}", f"pad_y = {c['pad_y']}", f"termination_x = {c['termination_x']}",
          f"absorber_periods = {c['absorber_periods']:g}", "",
          "[waveguide]", f"type = {s['waveguide']['type']}", "",
          "[defect]", f"type = {d['type']}", f"radius = {d['radius']}", f"row = {d['row']}",
          f"dx = {d['dx']}", f"dy = {d['dy']}", "", "[rows]"]
    for j, r in sorted(s["rows"].items()):
        L.append(f"{j} = {r}")
    L += ["", "[rods]"]
    for e in s["rods"]:
        if e["op"] == "add":
            L.append(f"add {e['x']} {e['y']} {e['radius']}")
        else:
            L.append(f"remove {e['x']} {e['y']}")
    return "\n".join(L) + "\n"


def default_structure() -> dict[str, Any]:
    return parse_structure_text("[lattice]\n", "<default>")


__all__ = ["STRUCTURE_SCHEMA", "parse_structure_text", "load_structure", "structure_from_dict",
           "structure_to_json", "structure_signature", "structure_docs", "structure_text",
           "default_structure", "validate_structure", "TextFileError", "math"]
