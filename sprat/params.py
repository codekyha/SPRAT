"""The parameter file (``.par``): run mode, analyte, numerics, source, stopping rules, sweeps.

Example::

    [run]
    mode = harminv              # harminv | spectrum | reference | field | bands | pwe
    label = auto
    tag = reference
    q_est = 6927                # expected Q: sets the harminv time (with margin) and the stopping ceilings

    [analyte]
    n = 1.33
    k = 0.0

    [numerics]
    resolution = 24
    courant = 0.5
    pml = 1.0
    subpixel = true
    symmetry = auto

    [source]
    fcen = auto                 # auto = a_nm / target_wavelength_nm
    fwidth = 0.06
    cutoff = 5

    [harminv]
    t = auto                    # auto = margin * q_est / (pi fcen), rounded up to 100
    margin = 4.5

    [sweep]
    structure.defect.radius = 0.050 : 0.110 : 0.005
    analyte.n = 1.300, 1.330, 1.450

Any key of the structure file can be swept or overridden through the dotted path
``structure.<section>.<key>``; any key of this file through ``<section>.<key>``.
"""

from __future__ import annotations

import copy
import itertools
import math
from typing import Any

from ._textfile import (Field, Section, TextFileError, apply_schema, dotted_get, dotted_set, parse_list,
                        parse_scalar, parse_text, schema_markdown)
from .structure import STRUCTURE_SCHEMA, structure_signature, validate_structure

MODES = ("harminv", "spectrum", "reference", "field", "bands", "pwe")

PARAMS_SCHEMA: dict[str, Section] = {
    "run": Section(
        "What to run and how to name it.",
        {
            "mode": Field("enum", "harminv", "harminv: complex eigenfrequency (Q, f_r) by harmonic inversion; "
                          "spectrum: transmission P_out(f); reference: cavity-less transmission (normalisation); "
                          "field: DFT field map at f_res and the analyte energy fraction; bands: MPB band structure; "
                          "pwe: plane-wave-expansion checks", choices=MODES),
            "label": Field("str", "auto", "record name; 'auto' builds it from the resolved values"),
            "tag": Field("str", "", "free text stored in every record (the name of a stage or a series)"),
            "overwrite": Field("bool", False, "rerun even if the record exists"),
            "q_est": Field("float", 2000.0, "expected quality factor: sets the harminv time through the margin rule "
                           "and the run-time ceilings of the spectrum and field modes", minimum=1.0),
            "quiet": Field("bool", True, "silence the Meep log"),
        }),
    "analyte": Section(
        "The medium between the rods.",
        {
            "n": Field("float", 1.33, "refractive index of the analyte", minimum=1.0),
            "k": Field("float", 0.0, "extinction coefficient; > 0 adds the corresponding conductivity", minimum=0.0),
        }),
    "numerics": Section(
        "Discretisation and boundary numerics.",
        {
            "resolution": Field("int", 24, "grid points per lattice constant", unit="1/a", minimum=4),
            "courant": Field("float", 0.5, "Courant factor S = c dt / dx", minimum=0.05, maximum=0.7),
            "pml": Field("float", 1.0, "PML thickness (all sides with pml termination; y sides with the absorber)",
                         unit="a", minimum=0.1),
            "subpixel": Field("bool", True, "subpixel smoothing of the permittivity"),
            "symmetry": Field("enum", "auto", "mirror symmetry in x: auto (even when dx = 0 and the structure is "
                              "mirror-symmetric), none, even, odd", choices=("auto", "none", "even", "odd")),
        }),
    "source": Section(
        "The Gaussian pulse.",
        {
            "fcen": Field("float_or_auto", "auto", "centre frequency a/lambda; auto = a_nm / target_wavelength_nm"),
            "fwidth": Field("float", 0.06, "frequency width a/lambda", minimum=1e-4),
            "cutoff": Field("float", 5.0, "Gaussian cutoff in widths (Meep default 5)", minimum=1.0),
        }),
    "harminv": Section(
        "Harmonic inversion of E_z at a probe point next to the defect.",
        {
            "t": Field("float_or_auto", "auto", "signal length after the source is off (Meep time units); "
                       "auto = margin * q_est / (pi fcen) rounded up to a multiple of 100", unit="a/c"),
            "margin": Field("float", 4.5, "target Q_lim / Q_est for the auto rule", minimum=1.0),
            "t_max": Field("float", 3.0e4, "ceiling of the second pass when auto_t is on", unit="a/c", minimum=1.0),
            "auto_t": Field("bool", False, "two-pass rule: rerun with t = auto_factor * Q / (pi fcen) when the first "
                            "pass finds a Q that needs more than 1.5 t"),
            "auto_factor": Field("float", 4.0, "factor of the two-pass rule", minimum=1.0),
            "probe_dx": Field("float", 0.13, "probe/source offset from the defect centre along x", unit="a"),
            "probe_dy": Field("float", 0.07, "probe/source offset from the defect centre along y", unit="a"),
        }),
    "spectrum": Section(
        "Transmission spectrum and cavity-less reference.",
        {
            "nfreq": Field("int", 201, "number of flux frequencies", minimum=3),
            "dft_tol": Field("float", 1e-6, "relative tolerance of the DFT-decay stopping rule", minimum=1e-12),
            "t_max": Field("float", 4.0e4, "hard ceiling of the run time after the source", unit="a/c", minimum=1.0),
            "flux_width": Field("float", 3.0, "height of the flux planes across the guide", unit="a", minimum=0.5),
            "source_width": Field("float", 2.0, "height of the line source across the guide", unit="a", minimum=0.5),
        }),
    "field": Section(
        "DFT field map at one frequency and the analyte energy fraction.",
        {
            "f_res": Field("float", 0.0, "frequency of the map (a/lambda); required, normally the resonance of the "
                           "matching harminv record", minimum=0.0),
            "t_max": Field("float", 2.0e4, "hard ceiling of the run time after the source", unit="a/c", minimum=1.0),
        }),
    "bands": Section(
        "MPB band structures (mode = bands).",
        {
            "task": Field("enum", "bulk", "bulk: TM and TE bands along Gamma-X-M-Gamma and the gap; sweep: TM gap "
                          "against the analyte index; w1: W1 supercell projection; calibration: the TM gap centre "
                          "and the lattice constant a = f_mid * target_wavelength_nm only",
                          choices=("bulk", "sweep", "w1", "calibration")),
            "resolution": Field("int", 32, "MPB grid points per lattice constant", minimum=8),
            "k_points": Field("int", 24, "interpolation points per segment of the k path", minimum=2),
            "num_bands": Field("int_or_auto", "auto", "bands to compute; auto = 8 (bulk) or 24 (w1)"),
            "supercell": Field("int", 15, "rows of the W1 supercell (odd)", minimum=3),
            "sweep_start": Field("float", 1.00, "first analyte index of the sweep", minimum=1.0),
            "sweep_stop": Field("float", 1.50, "last analyte index of the sweep", minimum=1.0),
            "sweep_step": Field("float", 0.05, "step of the sweep", minimum=1e-4),
        }),
    "pwe": Section(
        "Plane-wave-expansion checks (mode = pwe): bulk gap, W1 band, complex band structure, Fabry-Perot "
        "reading, barrier-row radius, loss budget.",
        {
            "quick": Field("bool", False, "smaller plane-wave cut-off (M = 7 instead of 10) for a fast check"),
            "reference_f": Field("float", 0.30461, "resonance frequency of the reference point (r_d = 0.060a)"),
            "second_f": Field("float", 0.29230, "resonance frequency of the second radius (r_d = 0.100a)"),
            "kappa_reference": Field("float", 1.9595, "measured kappa at the reference point, for the deviation column"),
            "kappa_second": Field("float", 1.7386, "measured kappa at the second radius, for the deviation column"),
        }),
    "sweep": Section("Sweeps: `<dotted key> = a:b:step` or `v1, v2, ...`; `zip = key1, key2` pairs two lists "
                     "element by element instead of taking their product.", free_keys=True, free_type="str"),
}

# dotted keys of the parameter file that may not be swept
_NOT_SWEEPABLE = {"run.label", "run.mode", "sweep"}


def parse_params_text(text: str, path: str = "<params>") -> dict[str, Any]:
    """Parse and validate a parameter file (sweeps are kept as lists of typed values)."""
    raw = apply_schema(parse_text(text, path), PARAMS_SCHEMA, path)
    p: dict[str, Any] = {k: dict(v) for k, v in raw.items()}
    sweeps: dict[str, list] = {}
    zips: list[list[str]] = []
    for key, value in raw["sweep"].items():
        if key == "_lines":
            continue
        if key == "zip":
            for group in str(value).split(";"):
                names = [x.strip().lower() for x in group.split(",") if x.strip()]
                if len(names) < 2:
                    raise TextFileError(f"{path}: [sweep] zip needs at least two keys, got {value!r}")
                zips.append(names)
            continue
        fld = _lookup_field(key, path)
        sweeps[key] = parse_list(str(value), fld.type, f"{path} [sweep] {key}", fld)
    for group in zips:
        for name in group:
            if name not in sweeps:
                raise TextFileError(f"{path}: [sweep] zip refers to {name!r}, which is not swept")
        lengths = {len(sweeps[n]) for n in group}
        if len(lengths) != 1:
            raise TextFileError(f"{path}: [sweep] zipped keys {group} have different lengths")
    p["sweep"] = {"values": sweeps, "zip": zips}
    if p["run"]["mode"] == "field" and p["field"]["f_res"] <= 0 and "field.f_res" not in sweeps:
        raise TextFileError(f"{path}: mode = field needs [field] f_res > 0")
    return p


def _lookup_field(dotted: str, path: str) -> Field:
    """The schema field behind a dotted key of either file."""
    parts = dotted.split(".")
    if parts[0] == "structure":
        if len(parts) == 3 and parts[1] == "rows":
            return Field("float", None, "row radius", minimum=0.0, maximum=0.5)
        if len(parts) != 3 or parts[1] not in STRUCTURE_SCHEMA or parts[2] not in STRUCTURE_SCHEMA[parts[1]].fields:
            raise TextFileError(f"{path}: unknown structure key {dotted!r}")
        return STRUCTURE_SCHEMA[parts[1]].fields[parts[2]]
    if dotted in _NOT_SWEEPABLE or len(parts) != 2 or parts[0] not in PARAMS_SCHEMA \
            or parts[1] not in PARAMS_SCHEMA[parts[0]].fields:
        raise TextFileError(f"{path}: unknown or non-sweepable parameter key {dotted!r}")
    return PARAMS_SCHEMA[parts[0]].fields[parts[1]]


def load_params(path: str) -> dict[str, Any]:
    with open(path, encoding="utf-8") as fh:
        return parse_params_text(fh.read(), path)


def parse_override(text: str, path: str = "--set") -> tuple[str, Any]:
    """``section.key=value`` (or ``structure.section.key=value``) -> (dotted, typed value)."""
    if "=" not in text:
        raise TextFileError(f"{path}: expected key=value, got {text!r}")
    key, _, value = text.partition("=")
    key = key.strip().lower()
    fld = _lookup_field(key, path)
    return key, parse_scalar(value, fld.type, f"{path} {key}", fld)


def apply_override(structure: dict, params: dict, dotted: str, value: Any) -> None:
    parts = dotted.split(".")
    if parts[0] == "structure":
        if parts[1] == "rows":
            structure["rows"][int(parts[2])] = float(value)
        else:
            structure[parts[1]][parts[2]] = value
    else:
        params[parts[0]][parts[1]] = value


def expand_sweeps(structure: dict, params: dict) -> list[tuple[dict, dict, dict]]:
    """Return the list of (structure, params, overrides) for every point of the sweep grid.

    Sweeps not mentioned in a ``zip`` group form a Cartesian product; each zip group advances
    together.  The returned dicts are deep copies with the sweep values applied and the
    ``sweep`` block removed from the parameters.
    """
    sw = params.get("sweep", {"values": {}, "zip": []})
    values, groups = sw["values"], sw["zip"]
    grouped: set[str] = {n for g in groups for n in g}
    axes: list[list[dict[str, Any]]] = []
    for name in values:
        if name in grouped:
            continue
        axes.append([{name: v} for v in values[name]])
    for g in groups:
        axes.append([{n: values[n][i] for n in g} for i in range(len(values[g[0]]))])
    points: list[tuple[dict, dict, dict]] = []
    for combo in itertools.product(*axes) if axes else [()]:
        s = copy.deepcopy(structure)
        p = copy.deepcopy(params)
        p.pop("sweep", None)
        overrides: dict[str, Any] = {}
        for part in combo:
            for k, v in part.items():
                apply_override(s, p, k, v)
                overrides[k] = v
        validate_structure(s)
        points.append((s, p, overrides))
    return points


# --------------------------------------------------------------------------- resolution
def source_time(fwidth: float, cutoff: float) -> float:
    """Total duration of the Gaussian source in Meep time units: 2 cutoff / fwidth."""
    return 2.0 * cutoff / float(fwidth)


def harminv_time(margin: float, q_est: float, fcen: float) -> float:
    """t = margin Q_est / (pi fcen), rounded up to a multiple of 100."""
    return 100.0 * math.ceil(margin * q_est / (math.pi * fcen) / 100.0)


def q_limit(fcen: float, t: float) -> float:
    """Q_lim = pi fcen t, the Q of a mode that decays by a factor e over the signal; the margin Q_lim / Q >= 1 is a
    conservative sampling rule, since harmonic inversion resolves decays longer than the signal."""
    return math.pi * fcen * t


def resolve(structure: dict, params: dict) -> dict[str, Any]:
    """Fill every 'auto' of the parameter block from the structure (a_nm must be numeric)."""
    p = copy.deepcopy(params)
    p.pop("sweep", None)
    a_nm = structure["lattice"]["a_nm"]
    if a_nm == "calibrate":
        raise TextFileError("a_nm = calibrate must be resolved from a bulk band record before planning "
                            "(run `sprat bands` first or give a_nm explicitly)")
    if p["source"]["fcen"] == "auto":
        p["source"]["fcen"] = float(a_nm) / float(structure["lattice"]["target_wavelength_nm"])
    if p["harminv"]["t"] == "auto":
        p["harminv"]["t"] = harminv_time(p["harminv"]["margin"], p["run"]["q_est"], p["source"]["fcen"])
    if p["bands"]["num_bands"] == "auto":
        p["bands"]["num_bands"] = 24 if p["bands"]["task"] == "w1" else 8
    if p["run"]["label"] == "auto":
        p["run"]["label"] = auto_label(structure, p)
    return p


def auto_label(structure: dict, p: dict) -> str:
    """Deterministic record name from the resolved values (English, file-system safe)."""
    mode = p["run"]["mode"]
    if mode == "bands":
        b = p["bands"]
        if b["task"] == "sweep":
            return f"bands_sweep_r{structure['lattice']['rod_radius']:g}_res{b['resolution']}"
        return f"bands_{b['task']}_na{p['analyte']['n']:.4f}_r{structure['lattice']['rod_radius']:g}_res{b['resolution']}"
    if mode == "pwe":
        return f"pwe_na{p['analyte']['n']:.4f}_r{structure['lattice']['rod_radius']:g}" + ("_quick" if p["pwe"]["quick"] else "")
    sig = structure_signature(structure)
    if mode == "reference":
        # the reference run has no defect; drop the defect part of the signature
        sig = "_".join(x for x in sig.split("_") if not (x.startswith("rd") or x.startswith("dx") or x.startswith("dy") or x.startswith("row")))
    core = f"{mode}_{sig}_na{p['analyte']['n']:.4f}_res{p['numerics']['resolution']}_pml{p['numerics']['pml']:g}"
    if mode == "harminv":
        core += f"_t{p['harminv']['t']:g}"
    elif mode in ("spectrum", "reference"):
        core += f"_nf{p['spectrum']['nfreq']}"
    elif mode == "field":
        core += f"_f{p['field']['f_res']:.6f}"
    if p["analyte"]["k"] > 0:
        core += f"_k{p['analyte']['k']:g}"
    if p["numerics"]["symmetry"] != "auto":
        core += f"_sym{p['numerics']['symmetry']}"
    return core


def params_docs() -> str:
    return schema_markdown(PARAMS_SCHEMA, "Parameter file (.par) reference")


def params_text(p: dict[str, Any], comment: str = "", sweeps: dict[str, str] | None = None) -> str:
    """Write a parameter dict as a .par file (used by the campaign builder)."""
    L = [f"# {line}" for line in comment.splitlines()] if comment else []
    for name, sec in PARAMS_SCHEMA.items():
        if name == "sweep":
            continue
        L.append(f"[{name}]")
        for key in sec.fields:
            v = p[name][key]
            if isinstance(v, bool):
                v = "true" if v else "false"
            L.append(f"{key} = {v}")
        L.append("")
    L.append("[sweep]")
    for k, v in (sweeps or {}).items():
        L.append(f"{k} = {v}")
    return "\n".join(L) + "\n"


def default_params() -> dict[str, Any]:
    return parse_params_text("[run]\n", "<default>")


__all__ = ["PARAMS_SCHEMA", "MODES", "parse_params_text", "load_params", "parse_override", "apply_override",
           "expand_sweeps", "resolve", "auto_label", "source_time", "harminv_time", "q_limit", "params_docs",
           "params_text", "default_params", "dotted_get", "dotted_set"]
