"""Importer for the 2026 campaign records (Turkish keys) into the English record schema.

The deposited data pack (Zenodo, ``zenodo_pack_v2``) holds the raw per-run JSON files
exactly as the campaign scripts wrote them.  ``import_deposit`` converts every one of them
into a ``sprat-record-1.0`` record without modifying the originals; each converted record
carries ``provenance.legacy_file`` and ``provenance.legacy_sha256`` so the chain back to the
deposit is explicit.  The key map is ``LEGACY_KEYS`` (rendered into docs/LEGACY_KEYS.md).
"""

from __future__ import annotations

import copy
import glob
import hashlib
import json
import os
import re
import tarfile
from typing import Any

import numpy as np

from . import SCHEMA_VERSION, __version__, records
from .params import default_params
from .structure import default_structure, structure_to_json

LEGACY_MEEP = "1.30.0"

MODE_MAP = {"harminv": "harminv", "spektrum": "spectrum", "norm": "reference", "alan": "field"}
LAYER_MAP = {"tarama": "sweep", "uretim": "production", "elle": "manual"}
SYMMETRY_MAP = {"oto": "auto", "yok": "none", "cift": "even", "tek": "odd"}

# Turkish key -> (English location, meaning); the documentation table of the importer.
LEGACY_KEYS = [
    ("parametreler", "structure + params", "the parameter block splits into the structure and the parameter file"),
    ("parametreler.mod", "params.run.mode", "harminv -> harminv, spektrum -> spectrum, norm -> reference, alan -> field"),
    ("parametreler.katman", "params.run.tag", "tarama -> sweep, uretim -> production, elle -> manual"),
    ("parametreler.etiket", "task.label", "the record name, kept verbatim"),
    ("parametreler.a_nm", "structure.lattice.a_nm", "lattice constant (nm)"),
    ("parametreler.r", "structure.lattice.rod_radius", "rod radius r/a"),
    ("parametreler.rd", "structure.defect.radius", "defect radius r_d/a"),
    ("parametreler.dx, dy", "structure.defect.dx, dy", "defect displacement (a)"),
    ("parametreler.nsep", "structure.defect.row", "separating rows N_sep"),
    ("parametreler.nx", "structure.cell.guide_periods", "guide periods"),
    ("parametreler.nkaplama", "structure.cell.cladding_rows", "cladding rows per side n_cl"),
    ("parametreler.pad, ypad", "structure.cell.pad_x, pad_y", "paddings (a)"),
    ("parametreler.sonlandirma", "structure.cell.termination_x", "pml -> pml, emici -> absorber"),
    ("parametreler.emici_kalinlik", "structure.cell.absorber_periods", "absorber thickness (a)"),
    ("parametreler.bariyer_satir_r / bariyer_satirlari", "structure.rows", "'j:r' -> {j: r}"),
    ("parametreler.kavite", "(mode)", "False only for the cavity-less reference run (mode reference)"),
    ("parametreler.na", "params.analyte.n", "analyte index"),
    ("parametreler.kappa", "params.analyte.k", "analyte extinction coefficient"),
    ("parametreler.res", "params.numerics.resolution", "grid points per a"),
    ("parametreler.dpml", "params.numerics.pml", "PML thickness (a)"),
    ("parametreler.simetri", "params.numerics.symmetry", "oto -> auto, yok -> none, cift -> even, tek -> odd"),
    ("parametreler.fcen, df", "params.source.fcen, fwidth", "source centre and width"),
    ("parametreler.t_harminv", "params.harminv.t", "requested harminv signal length"),
    ("parametreler.t_max_harminv", "params.harminv.t_max", "ceiling of the two-pass rule"),
    ("parametreler.otomatik_t, otomatik_kat", "params.harminv.auto_t, auto_factor", "two-pass rule"),
    ("parametreler.q_tahmin", "params.run.q_est", "expected Q"),
    ("parametreler.nfreq, dft_tol, t_max_spektrum", "params.spectrum.nfreq, dft_tol, t_max", "spectrum settings"),
    ("parametreler.f_rez, t_alan_max", "params.field.f_res, t_max", "field-map settings"),
    ("parametreler.sessiz", "params.run.quiet", "Meep log silenced"),
    ("olcum", "provenance", "the measurement block"),
    ("olcum.baslangic", "provenance.started", "start time"),
    ("olcum.dugum", "provenance.host", "node name"),
    ("olcum.is_no, dizi_no", "provenance.job_id, array_id", "SLURM job and array ids"),
    ("olcum.adim", "provenance.steps", "time steps"),
    ("olcum.zaman", "provenance.sim_time", "simulated time"),
    ("olcum.piksel", "provenance.pixels", "grid points"),
    ("olcum.piksel_adim", "provenance.pixel_steps", "work"),
    ("olcum.piksel_adim_s", "provenance.throughput", "pixel-steps per second"),
    ("olcum.duvar_s", "provenance.wall_s", "wall seconds"),
    ("olcum.simetri", "provenance.symmetry", "sector actually used"),
    ("olcum.kod_sha256", "software.code_sha256", "SHA-256 of the script that wrote the record (v4 only)"),
    ("olcum.paket", "provenance.legacy_package", "'v4' for the verification runs"),
    ("sonuc", "result", "the result block"),
    ("sonuc.modlar[i].f, Q, genlik, hata, lambda_nm, FWHM_nm", "result.modes[i].f, Q, amplitude, error, wavelength_nm, fwhm_nm",
     "harmonic-inversion modes, strongest first"),
    ("sonuc.t_harminv_kullanilan", "result.t_used", "signal length actually used"),
    ("sonuc.otomatik_gecmis", "result.auto_history", "the passes of the two-pass rule"),
    ("sonuc.Q_cozunurluk_siniri", "result.Q_limit", "pi fcen t"),
    ("sonuc.f, P_in, P_out, t_bitis", "result.f, flux_in, flux_out, t_end", "spectrum and reference runs"),
    ("sonuc.f_kesri, lambda_r_nm, S_teori_nm_RIU", "result.analyte_energy_fraction, wavelength_nm, S_first_order_nm_per_RIU",
     "field runs"),
    ("alan_<label>.npz: Ez, eps, sx, sy, dpml, res, f_rez, lambda_nm", "<label>.npz: Ez, eps, sx, sy, pml, resolution, f_res, wavelength_nm (+ Lx, Ly, edge_x, n, rod_eps)",
     "field maps; the importer adds the interior size and the medium"),
    ("bant_bulk_na*.json, bant_tarama.json, bant_w1_na*.json, a_nm.json", "bands records (task bulk / sweep / w1)",
     "gap: var -> exists, alt -> lower, ust -> upper, orta -> mid, yuzde -> relative_width_pct"),
    ("pwe_v4_feasibility_results.json", "pwe record", "kappa keys 'r_d=0.060a (reference)' -> reference, 'r_d=0.100a' -> second"),
]


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------- one record
def convert_record(d: dict[str, Any], filename: str, sha: str = "", source: str = "legacy") -> dict[str, Any]:
    p, o, s = d.get("parametreler", {}), d.get("olcum", {}), d.get("sonuc", {})
    mode = MODE_MAP[p.get("mod", "harminv")]
    st = default_structure()
    la, c, df = st["lattice"], st["cell"], st["defect"]
    la["a_nm"] = float(p.get("a_nm", 481.4))
    la["rod_radius"] = float(p.get("r", 0.20))
    c.update(guide_periods=int(p.get("nx", 25)), cladding_rows=int(p.get("nkaplama", 12)),
             pad_x=float(p.get("pad", 1.5)), pad_y=float(p.get("ypad", 1.0)),
             termination_x=("absorber" if p.get("sonlandirma", "pml") == "emici" else "pml"),
             absorber_periods=float(p.get("emici_kalinlik", 8.0)))
    df.update(type="rod", radius=float(p.get("rd", 0.0)), row=int(p.get("nsep", 4)), dx=float(p.get("dx", 0.0)),
              dy=float(p.get("dy", 0.0)))
    rows: dict[int, float] = {}
    if p.get("bariyer_satirlari"):
        rows = {int(k): float(v) for k, v in p["bariyer_satirlari"].items()}
    elif p.get("bariyer_satir_r"):
        for part in str(p["bariyer_satir_r"]).split(","):
            j, r = part.split(":")
            rows[int(j)] = float(r)
    st["rows"] = rows
    pr = default_params()
    pr.pop("sweep", None)
    pr["run"].update(mode=mode, label=p.get("etiket", os.path.splitext(filename)[0]),
                     tag=LAYER_MAP.get(p.get("katman", ""), p.get("katman", "")), overwrite=False,
                     q_est=float(p.get("q_tahmin", 2000.0)), quiet=bool(p.get("sessiz", True)))
    pr["analyte"].update(n=float(p.get("na", 1.33)), k=float(p.get("kappa", 0.0)))
    pr["numerics"].update(resolution=int(p.get("res", 24)), courant=0.5, pml=float(p.get("dpml", 1.0)), subpixel=True,
                          symmetry=SYMMETRY_MAP.get(p.get("simetri", "oto"), "auto"))
    fcen = float(p.get("fcen") or (la["a_nm"] / 1550.0))
    pr["source"].update(fcen=fcen, fwidth=float(p.get("df", 0.06)), cutoff=5.0)
    t_req = p.get("t_harminv")
    t_used = s.get("t_harminv_kullanilan", t_req)
    pr["harminv"].update(t=float(t_used if t_used is not None else (t_req or 0.0)), margin=4.5,
                         t_max=float(p.get("t_max_harminv", 3.0e4)), auto_t=bool(p.get("otomatik_t", False)),
                         auto_factor=float(p.get("otomatik_kat", 4.0)), probe_dx=0.13, probe_dy=0.07)
    pr["spectrum"].update(nfreq=int(p.get("nfreq", 201)), dft_tol=float(p.get("dft_tol", 1e-6)),
                          t_max=float(p.get("t_max_spektrum", 4.0e4)), flux_width=3.0, source_width=2.0)
    pr["field"].update(f_res=float(p.get("f_rez", 0.0)), t_max=float(p.get("t_alan_max", 2.0e4)))
    rec = dict(
        schema=SCHEMA_VERSION,
        software=dict(name="sprat-legacy-import", version=__version__, meep=LEGACY_MEEP, python="",
                      code_sha256=o.get("kod_sha256", ""), legacy_script="02_kavite.py"),
        task=dict(label=pr["run"]["label"], tag=pr["run"]["tag"], created=o.get("baslangic", ""), source=source,
                  legacy_layer=LAYER_MAP.get(p.get("katman", ""), p.get("katman", "")), legacy_wp=_wp(pr["run"]["label"])),
        structure=structure_to_json(st), params=pr,
        provenance=dict(host=o.get("dugum", ""), started=o.get("baslangic", ""), wall_s=o.get("duvar_s"),
                        steps=o.get("adim"), sim_time=o.get("zaman"), pixels=o.get("piksel"),
                        pixel_steps=o.get("piksel_adim"), throughput=o.get("piksel_adim_s"),
                        symmetry=SYMMETRY_MAP.get(o.get("simetri", ""), o.get("simetri", "")),
                        job_id=str(o.get("is_no", "")), array_id=str(o.get("dizi_no", "")),
                        legacy_package=o.get("paket", ""), legacy_file=filename, legacy_sha256=sha),
        result={})
    if mode == "harminv":
        modes = [dict(f=float(m["f"]), Q=float(m["Q"]), amplitude=float(m.get("genlik", 0.0)),
                      error=float(m.get("hata", 0.0)), wavelength_nm=float(m.get("lambda_nm", la["a_nm"] / m["f"])),
                      fwhm_nm=float(m.get("FWHM_nm", (la["a_nm"] / m["f"]) / m["Q"])))
                 for m in (s.get("modlar") or [])]
        modes.sort(key=lambda m: -m["amplitude"])
        rec["result"] = dict(modes=modes, t_used=float(t_used) if t_used is not None else None,
                             auto_history=[dict(t=h.get("t"), Q=h.get("Q")) for h in (s.get("otomatik_gecmis") or [])],
                             Q_limit=s.get("Q_cozunurluk_siniri"))
    elif mode in ("spectrum", "reference"):
        rec["result"] = dict(f=s.get("f", []), flux_in=s.get("P_in", []), flux_out=s.get("P_out", []),
                             t_end=s.get("t_bitis"))
    elif mode == "field":
        rec["result"] = dict(analyte_energy_fraction=s.get("f_kesri"), wavelength_nm=s.get("lambda_r_nm"),
                             S_first_order_nm_per_RIU=s.get("S_teori_nm_RIU"))
    return rec


_WP_RE = re.compile(r"^v4_(H\d+[a-z]?|regresyon|duman)_")
_WP_MAP = {"regresyon": "regression", "duman": "smoke"}


def _wp(label: str) -> str:
    """Work-package code of a 2026 verification label (``v4_H4b_...`` -> ``H4b``); the two
    non-numbered packages are translated (``regresyon`` -> ``regression``, ``duman`` -> ``smoke``)."""
    m = _WP_RE.match(label or "")
    return _WP_MAP.get(m.group(1), m.group(1)) if m else ""


_LABEL_RE = re.compile(r"rd(?P<rd>[0-9.]+)_dx(?P<dx>[+-][0-9.]+)_dy(?P<dy>[+-][0-9.]+)_s(?P<nsep>\d+)_na(?P<na>[0-9.]+)"
                       r"_res(?P<res>\d+)_nk(?P<nk>\d+)")


def field_record_from_npz(npz_path: str, out_dir: str, sha: str = "") -> dict[str, Any] | None:
    """Build a field record from a legacy ``alan_<label>.npz`` map (the campaign stored no JSON for them)."""
    name = os.path.basename(npz_path)
    label = name[5:-4] if name.startswith("alan_") else name[:-4]
    m = _LABEL_RE.search(label)
    if not m:
        return None
    z = np.load(npz_path)
    d = dict(parametreler=dict(mod="alan", katman=("uretim" if "uretim" in label else "tarama"), etiket=label,
                               a_nm=481.4, r=0.20, rd=float(m["rd"]), dx=float(m["dx"]), dy=float(m["dy"]),
                               nsep=int(m["nsep"]), na=float(m["na"]), res=int(m["res"]), nkaplama=int(m["nk"]),
                               nx=25, pad=1.5, ypad=1.0, dpml=float(z["dpml"]), f_rez=float(z["f_res"] if "f_res" in z else z["f_rez"]),
                               simetri="oto"),
             olcum={}, sonuc=dict(lambda_r_nm=float(z["lambda_nm"])))
    rec = convert_record(d, name, sha, source="legacy-campaign")
    eps = np.asarray(z["eps"], dtype=np.float64)
    ez = np.asarray(z["Ez"])
    n0, n1 = min(ez.shape[0], eps.shape[0]), min(ez.shape[1], eps.shape[1])
    ez, eps = ez[:n0, :n1], eps[:n0, :n1]
    from .fdtd import analyse_field
    res = analyse_field(ez, eps, float(m["na"]), 11.9025, 481.4, float(d["parametreler"]["f_rez"]))
    sx, sy, dpml = float(z["sx"]), float(z["sy"]), float(z["dpml"])
    new_npz = os.path.join(out_dir, label + ".npz")
    np.savez_compressed(new_npz, Ez=ez.astype(np.complex64), eps=eps.astype(np.float32), sx=sx, sy=sy,
                        Lx=sx - 2 * dpml, Ly=sy - 2 * dpml, pml=dpml, edge_x=dpml, resolution=int(z["res"]),
                        f_res=float(d["parametreler"]["f_rez"]), wavelength_nm=float(z["lambda_nm"]), n=float(m["na"]),
                        rod_eps=11.9025)
    res["field_file"] = os.path.basename(new_npz)
    res["legacy_field_file"] = name
    rec["result"] = res
    rec["provenance"]["symmetry"] = "even" if abs(float(m["dx"])) < 1e-12 else "none"
    return rec


# --------------------------------------------------------------------------- bands and pwe
def _gap(g: dict) -> dict:
    return dict(exists=bool(g.get("var")), lower=g.get("alt"), upper=g.get("ust"), mid=g.get("orta"),
                relative_width_pct=float(g.get("yuzde", 0.0)))


def bands_records(bands_dir: str) -> list[dict[str, Any]]:
    out = []
    a_nm_file = os.path.join(bands_dir, "a_nm.json")
    a_nm_info = json.load(open(a_nm_file)) if os.path.exists(a_nm_file) else {}
    for path in sorted(glob.glob(os.path.join(bands_dir, "bant_bulk_na*.json"))):
        d = json.load(open(path))
        n = float(re.search(r"na([0-9]+(?:\.[0-9]+)?)", os.path.basename(path)).group(1))
        rec = _bands_skeleton("bulk", n, os.path.basename(path), sha256_file(path), a_nm_info.get("res", 32))
        rec["result"].update(tm=dict(gap=_gap(d["tm"]["gap"]), freqs=d["tm"]["freqs"]),
                             te=dict(gap=_gap(d["te"]["gap"]), freqs=d["te"]["freqs"]), gap=_gap(d["tm"]["gap"]),
                             a_nm=d.get("a_nm"), target_wavelength_nm=a_nm_info.get("lambda_hedef", 1550.0),
                             a_nm_file=a_nm_info)
        out.append(rec)
    path = os.path.join(bands_dir, "a_nm_bant_ortasi.json")
    if os.path.exists(path):
        d = json.load(open(path))
        rec = _bands_skeleton("calibration", float(d.get("na", 1.33)), os.path.basename(path), sha256_file(path), d.get("res", 32))
        rec["params"]["run"]["label"] = "bands_calibration_band_centre_legacy"
        rec["task"]["label"] = rec["params"]["run"]["label"]
        rec["result"].update(a_nm=d.get("a_nm"), f_mid=d.get("f_orta"), target_wavelength_nm=d.get("lambda_hedef", 1550.0),
                             rod_radius=d.get("r_over_a", 0.20), note="lattice constant from the TM band-gap centre, a = f_mid * target")
        out.append(rec)
    path = os.path.join(bands_dir, "bant_tarama.json")
    if os.path.exists(path):
        d = json.load(open(path))
        rec = _bands_skeleton("sweep", 1.33, "bant_tarama.json", sha256_file(path), a_nm_info.get("res", 32))
        rec["result"]["sweep"] = [dict(n=float(x["na"]), gap=_gap(x)) for x in d]
        out.append(rec)
    for path in sorted(glob.glob(os.path.join(bands_dir, "bant_w1_na*.json"))):
        d = json.load(open(path))
        n = float(re.search(r"na([0-9]+(?:\.[0-9]+)?)", os.path.basename(path)).group(1))
        rec = _bands_skeleton("w1", n, os.path.basename(path), sha256_file(path), a_nm_info.get("res", 32))
        rec["result"].update(kx_over_2pi=d["kx"], freqs=d["freqs"], gap=_gap(d["gap"]), supercell=15)
        out.append(rec)
    return out


def _bands_skeleton(task: str, n: float, filename: str, sha: str, res: int) -> dict[str, Any]:
    st = default_structure()
    pr = default_params()
    pr.pop("sweep", None)
    pr["run"].update(mode="bands", label=f"bands_{task}_na{n:.4f}_legacy" if task != "sweep" else "bands_sweep_legacy",
                     tag="legacy")
    pr["analyte"]["n"] = n
    pr["bands"].update(task=task, resolution=int(res), num_bands=(24 if task == "w1" else 8))
    return dict(schema=SCHEMA_VERSION,
                software=dict(name="sprat-legacy-import", version=__version__, meep=LEGACY_MEEP, python="",
                              code_sha256="", legacy_script="01_bant_yapisi.py"),
                task=dict(label=pr["run"]["label"], tag="legacy", created="", source="legacy-campaign"),
                structure=structure_to_json(st), params=pr,
                provenance=dict(host="", started="", legacy_file=filename, legacy_sha256=sha),
                result=dict(task=task, n=n))


def pwe_record(path: str) -> dict[str, Any]:
    d = json.load(open(path))
    d = copy.deepcopy(d)
    ren = {"r_d=0.060a (reference)": "reference", "r_d=0.100a": "second"}
    for block in ("kappa", "reciprocity"):
        if block in d:
            d[block] = {ren.get(k, k): v for k, v in d[block].items()}
    for v in d.get("kappa", {}).values():
        if "kappa_manuscript" in v:
            v["kappa_measured"] = v.pop("kappa_manuscript")
    if "model" in d and "eps_si" in d["model"]:
        d["model"]["eps_rod"] = d["model"].pop("eps_si")
    st = default_structure()
    pr = default_params()
    pr.pop("sweep", None)
    pr["run"].update(mode="pwe", label="pwe_legacy", tag="legacy")
    return dict(schema=SCHEMA_VERSION,
                software=dict(name="sprat-legacy-import", version=__version__, meep="", python="", code_sha256="",
                              legacy_script="pwe_v4_feasibility.py"),
                task=dict(label="pwe_legacy", tag="legacy", created="", source="legacy-verification"),
                structure=structure_to_json(st), params=pr,
                provenance=dict(host="", started="", legacy_file=os.path.basename(path), legacy_sha256=sha256_file(path)),
                result=d)


# --------------------------------------------------------------------------- the deposit
def _extract(deposit: str, name: str, work: str, fallback: bool = True) -> str | None:
    """Unpack <name>.tar.gz of the deposit into work/<name>/ (or use an existing directory). With ``fallback`` the
    deposit directory itself stands in for a missing <name> (the layout of the earliest packs); without it a missing
    component is None."""
    tgz = os.path.join(deposit, name + ".tar.gz")
    if os.path.exists(tgz):
        dest = os.path.join(work, name)
        if not os.path.isdir(dest):
            os.makedirs(dest, exist_ok=True)
            with tarfile.open(tgz) as tf:
                tf.extractall(dest)
        inner = os.path.join(dest, name)
        return inner if os.path.isdir(inner) else dest
    for cand in ((os.path.join(deposit, name), deposit) if fallback else (os.path.join(deposit, name),)):
        if os.path.isdir(cand):
            return cand
    return None


def import_deposit(deposit: str, out: str, log=print) -> dict[str, int]:
    """Convert a Zenodo data pack (packed or unpacked) into English records under ``out``."""
    os.makedirs(out, exist_ok=True)
    work = os.path.join(out, "_legacy_extract")
    counts = dict(harminv=0, verification=0, spectra=0, references=0, fields=0, bands=0, pwe=0, skipped=0)

    def write(rec: dict, source: str) -> None:
        rec["task"]["source"] = source
        records.save(rec, records.record_path(out, rec["task"]["label"]))

    hj = _extract(deposit, "harminv_json", work)
    if hj:
        for path in sorted(glob.glob(os.path.join(hj, "*.json"))):
            try:
                rec = convert_record(json.load(open(path)), os.path.basename(path), sha256_file(path), "legacy-campaign")
            except Exception as e:                       # noqa: BLE001
                log(f"skipped {path}: {e}")
                counts["skipped"] += 1
                continue
            write(rec, "legacy-campaign")
            counts["harminv"] += 1
    v4 = _extract(deposit, "v4_runs", work)
    if v4:
        sdir = os.path.join(v4, "sonuclar")
        for path in sorted(glob.glob(os.path.join(sdir, "v4_*.json"))):
            try:
                rec = convert_record(json.load(open(path)), os.path.basename(path), sha256_file(path), "legacy-verification")
            except Exception as e:                       # noqa: BLE001
                log(f"skipped {path}: {e}")
                counts["skipped"] += 1
                continue
            write(rec, "legacy-verification")
            counts["verification"] += 1
        pwe_path = os.path.join(sdir, "pwe_v4_feasibility_results.json")
        if os.path.exists(pwe_path):
            write(pwe_record(pwe_path), "legacy-verification")
            counts["pwe"] += 1
    sp = _extract(deposit, "spectra", work)
    if sp:
        for path in sorted(glob.glob(os.path.join(sp, "*.json"))):
            rec = convert_record(json.load(open(path)), os.path.basename(path), sha256_file(path), "legacy-campaign")
            write(rec, "legacy-campaign")
            counts["spectra" if rec["params"]["run"]["mode"] == "spectrum" else "references"] += 1
    # spectra of the verification runs (raw_h14b: the same seven spectra with the guide continued into the absorber)
    vs = _extract(deposit, "v4_spectra", work, fallback=False)
    if vs:
        for path in sorted(glob.glob(os.path.join(vs, "*.json"))):
            rec = convert_record(json.load(open(path)), os.path.basename(path), sha256_file(path), "legacy-verification")
            rec["task"]["legacy_wp"] = rec["task"].get("legacy_wp") or ("H14b" if rec["task"]["label"].endswith("_em8b") else "")
            write(rec, "legacy-verification")
            counts["spectra" if rec["params"]["run"]["mode"] == "spectrum" else "references"] += 1
    fm = _extract(deposit, "field_maps", work)
    if fm:
        for path in sorted(glob.glob(os.path.join(fm, "*.npz"))):
            rec = field_record_from_npz(path, out, sha256_file(path))
            if rec:
                write(rec, "legacy-campaign")
                counts["fields"] += 1
    bd = _extract(deposit, "bands", work)
    if bd:
        for rec in bands_records(bd):
            write(rec, "legacy-campaign")
            counts["bands"] += 1
    man_path = os.path.join(out, "_IMPORT_MANIFEST.json")
    name = os.path.basename(os.path.abspath(deposit))      # the name of the deposit directory only: no path of the machine
    entry = dict(deposit=name, counts=counts)
    earlier = []
    if os.path.isfile(man_path):                            # a second deposit imported into the same directory
        try:
            with open(man_path, encoding="utf-8") as fh:
                old_man = json.load(fh)
            earlier = old_man.get("deposits") or [dict(deposit=old_man.get("deposit"), counts=old_man.get("counts"))]
        except (ValueError, OSError):
            earlier = []
    deposits = [d for d in earlier if d.get("deposit") != name] + [entry]
    with open(man_path, "w", encoding="utf-8") as fh:
        json.dump(dict(deposit=deposits[0]["deposit"] if len(deposits) == 1 else [d["deposit"] for d in deposits],
                       counts=counts if len(deposits) == 1 else {k: sum(int((d.get("counts") or {}).get(k, 0)) for d in deposits)
                                                                   for k in counts},
                       deposits=deposits, importer=f"sprat {__version__}",
                       key_map=[dict(legacy=a, english=b, meaning=c) for a, b, c in LEGACY_KEYS]), fh, indent=1)
    log(f"imported {counts}")
    return counts


def keys_markdown() -> str:
    lines = ["# Legacy key map", "",
             "The original scripts of 2026 wrote their records with Turkish keys.  `sprat import-legacy` converts them with this map; "
             "the original files are never modified and every converted record stores `provenance.legacy_file` and "
             "`provenance.legacy_sha256`.", "",
             "| legacy | English | meaning |", "|---|---|---|"]
    for a, b, c in LEGACY_KEYS:
        lines.append(f"| `{a}` | `{b}` | {c} |")
    return "\n".join(lines) + "\n"


def import_cli(deposit: str, out: str) -> None:
    import_deposit(deposit, out)


__all__ = ["LEGACY_KEYS", "convert_record", "field_record_from_npz", "bands_records", "pwe_record", "import_deposit",
           "keys_markdown", "import_cli", "sha256_file"]
