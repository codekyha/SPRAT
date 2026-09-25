"""``sprat deposit``: build the Zenodo data pack from a records directory.

Layout of the pack::

    <out>/
      README.md                     contents, origin of the records, data dictionary, reproduction
      DATA_AVAILABILITY.md          the statement for the paper, then the counts of this pack and their source
      LICENSE.txt                   CC BY 4.0 for the data; MIT for the scripts inside raw_legacy.tar.gz
      zenodo_metadata.json          title, version, creators, description, keywords, licences, notes and
                                    the related identifiers that are known
      records.csv / records.jsonl   the flat table with the validity flags
      counts.json
      harminv_records.tar.gz        the harmonic-inversion records (JSON, English keys)
      spectra_records.tar.gz        spectra and cavity-less references
      field_records.tar.gz          field-map records and their .npz maps
      bands_records.tar.gz          band computations
      pwe_records.tar.gz            the plane-wave layer
      raw_legacy.tar.gz             (optional) the original files of the 2026 computational runs
      raw_h14b.tar.gz               (optional) the spectra with the guide continued into the absorber, with the job
                                    files, the criterion fixed before the runs, its amendment and the grading
      <included files>              (optional) files passed with --include, e.g. the supplementary PDF
                                    and the source archive sprat-<version>.zip
      analysis.tar.gz               unpacks to analysis/: the analysis outputs (registry, predictions,
                                    report) and manuscript_registry.json, the numbers registry of the
                                    manuscript
      predictions.tar.gz            unpacks to predictions/: the criteria fixed before the runs
      parameters_manifest.json      every constant of the model, the numerics and the analysis
      MANIFEST_sha256.txt           checksums of every file above

The DOI of the data record (``--data-doi``) and of the software record (``--software-doi``) enter
the statement, the README and the metadata; when the software is archived in the data record the
two are the same DOI, and ``--include sprat-<version>.zip`` puts its source archive in the pack
(with such an archive included and no software DOI given, the data DOI stands for both).
Until they are known the placeholders NNNNNNN and SSSSSSS stand in the texts; the metadata never
carries a placeholder as an identifier.
"""

from __future__ import annotations

import glob
import hashlib
import json
import os
import re
import shutil
import tarfile
import time
from typing import Any

from . import __version__, records
from .analysis import collect, pipeline, predictions
from .legacy import LEGACY_MEEP

REPOSITORY_URL = "https://github.com/codekyha/sprat"
DATA_DOI_PLACEHOLDER = "10.5281/zenodo.NNNNNNN"
SOFTWARE_DOI_PLACEHOLDER = "10.5281/zenodo.SSSSSSS"
# The numbers registry of the manuscript, kept in the repository checkout; copied into the pack as
# analysis/manuscript_registry.json.
MANUSCRIPT_REGISTRY = os.path.join("campaigns", "manuscript", "expected", "numbers_registry_manuscript_v10_2.json")

PAPER_TITLE = ("Discrete quality-factor control in a side-coupled photonic crystal microcavity: "
               "evanescent Bloch tunnelling and the finite-cell correction")
DATA_TITLE = "Data, code and supplementary material for: " + PAPER_TITLE
AUTHOR = dict(name="Oguz, Hasan", affiliation="Istanbul Okan University; Pamukkale University", orcid="0000-0001-7484-4415")
NOTES_FUNDING = ("The computational runs were performed on the UHeM Altay cluster (grant 5027772026). "
                 "Funding: Istanbul Okan University BAP OBAP2026010006; Pamukkale University BAP 2025ALDEP037.")
# Versions of the 2026 deposit that were prepared but never published (3.1.0 is the first published version; 3.0.0 was
# built from the same records with SPRAT 1.0.0 and replaced before publication).
UNPUBLISHED_LEGACY_VERSIONS = ("1.0.0", "2.0.0", "2.0.1", "2.0.2", "2.1.0", "2.2.0")
# The documents of the 2026 deposit that the files of this pack supersede.
LEGACY_SUPERSEDED = ("supplementary.pdf", "DATA_AVAILABILITY.md", "zenodo_metadata.json")

DATA_DICTIONARY = [
    ("structure.lattice.a_nm, rod_radius, rod_eps", "lattice constant (nm), rod radius (a), rod permittivity"),
    ("structure.cell.guide_periods, cladding_rows, pad_x, pad_y, termination_x, absorber_periods",
     "guide periods n_x, cladding rows per side n_cl, paddings (a), guide termination (pml | absorber), absorber thickness (a)"),
    ("structure.defect.radius, row, dx, dy", "defect radius r_d (a), separating rows N_sep, displacement (a)"),
    ("structure.rows", "per-row radius overrides {row: radius}"),
    ("params.run.mode, tag, q_est", "run mode, stage of the computational runs, expected Q used by the stopping rules"),
    ("params.analyte.n, k", "analyte index and extinction coefficient"),
    ("params.numerics.resolution, pml, symmetry", "grid points per a, PML thickness (a), mirror sector"),
    ("params.source.fcen, fwidth", "Gaussian source centre and width (a/lambda)"),
    ("params.harminv.t, margin, auto_t", "harmonic-inversion signal length, target margin, two-pass rule"),
    ("result.modes[i].f, Q, amplitude, error, wavelength_nm, fwhm_nm", "modes, strongest first; modes[0] is the resonance"),
    ("result.Q_limit", "Q_lim = pi f_cen t of the sampling rule (margin = Q_lim / Q)"),
    ("result.f, flux_in, flux_out", "spectrum and reference runs"),
    ("result.analyte_energy_fraction, wavelength_nm, field_file", "field runs; the .npz holds Ez, eps, sx, sy, Lx, Ly, resolution, f_res"),
    ("task.label, tag, source", "record name, stage of the computational runs, origin: legacy-campaign (systematic runs) or "
     "legacy-verification (verification runs) for records converted "
     "from the 2026 files, sprat for records written by SPRAT"),
    ("provenance.host, started, wall_s, steps, pixels, throughput, job_id", "where and how long the run took"),
    ("provenance.legacy_file, legacy_sha256", "the original 2026 file the record was converted from, and its checksum"),
    ("software.name, version, meep, legacy_script", "the program that wrote the record and its version, the Meep version of the run; "
     "for converted records sprat-legacy-import and the 2026 script that wrote the original file"),
    ("software.code_sha256", "SHA-256 of the code that wrote the record: the 2026 script for the verification runs "
     "(the systematic runs did not record it), sprat/fdtd.py for records written by SPRAT"),
    ("records.csv: margin, gap_position, resolved, valid, exclusion_reason",
     "Q_lim/Q, position in the TM gap (0-1), margin >= 1, margin >= 1 and gap position 0.2-0.8, reason if excluded"),
    ("records.csv: waveguide_mode",
     "True for the admitted records of the preliminary sweep (n_cl = 6, resolution 20) at N_sep = 1 or r_d >= 0.13a, whose "
     "strongest mode is a low-Q mode of the finite waveguide (Q 50 to 77); cavity records are those with valid = True and "
     "waveguide_mode = False"),
]


def _sha(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _owner_free(ti: tarfile.TarInfo) -> tarfile.TarInfo:
    """Archive members carry no user or group of the build machine."""
    ti.uid = ti.gid = 0
    ti.uname = ti.gname = ""
    return ti


def _tar(dest: str, files: list[str], root_name: str) -> int:
    with tarfile.open(dest, "w:gz") as tf:
        for f in files:
            tf.add(f, arcname=os.path.join(root_name, os.path.basename(f)), filter=_owner_free)
    return len(files)


# --------------------------------------------------------------------------- identifiers and origin
_DOI_RE = re.compile(r"^10\.\d{4,9}/\S+$")


def _tar_dir(dest: str, src: str, root_name: str) -> int:
    """Archive the directory ``src`` as ``root_name/...`` (sorted, owner-free) and return the number of files."""
    files = sorted(p for p in glob.glob(os.path.join(src, "**", "*"), recursive=True) if os.path.isfile(p))
    with tarfile.open(dest, "w:gz") as tf:
        for f in files:
            tf.add(f, arcname=os.path.join(root_name, os.path.relpath(f, src)), filter=_owner_free)
    return len(files)


def normalise_doi(doi: str) -> str:
    """``10.5281/zenodo.123``, ``doi:10.5281/...`` or ``https://doi.org/10.5281/...`` -> ``10.5281/zenodo.123``."""
    d = (doi or "").strip()
    for prefix in ("https://doi.org/", "http://doi.org/", "https://dx.doi.org/", "http://dx.doi.org/", "doi:"):
        if d.lower().startswith(prefix):
            d = d[len(prefix):]
            break
    if not _DOI_RE.match(d):
        raise ValueError(f"not a DOI: {doi!r} (expected the form 10.5281/zenodo.1234567)")
    return d


def same_record(data_doi: str, software_doi: str) -> bool:
    """True when the software is archived in the data record itself (one DOI for both)."""
    d, w = normalise_doi(data_doi), normalise_doi(software_doi)
    return (not is_placeholder(d)) and d == w


def sprat_archive_of(included: list[str] | None) -> str | None:
    """The SPRAT source archive among the included files, if any (sprat-<version>.zip or .tar.gz)."""
    for f in included or []:
        if f.lower().startswith("sprat-") and f.lower().endswith((".zip", ".tar.gz")):
            return f
    return None


def is_placeholder(doi: str) -> bool:
    """The default placeholders, or any identifier with a run of one repeated capital (XXXX, NNNNNNN)."""
    return doi in (DATA_DOI_PLACEHOLDER, SOFTWARE_DOI_PLACEHOLDER) or re.search(r"([A-Z])\1{3,}", doi) is not None


def origin(recs: list[dict]) -> dict[str, int]:
    """How many records were converted from the 2026 originals, how many SPRAT produced with Meep, and how many SPRAT
    computed with its plane-wave layer (``sprat_pwe``, counted apart: they are band-structure computations, not runs)."""
    imported = sum(1 for r in recs if str(r.get("task", {}).get("source", "")).startswith("legacy"))
    own_pwe = [r for r in recs if not str(r.get("task", {}).get("source", "")).startswith("legacy") and records.mode_of(r) == "pwe"]
    versions = sorted({str((r.get("software") or {}).get("version") or "") for r in own_pwe} - {""})
    out = dict(imported=imported, sprat=len(recs) - imported - len(own_pwe), sprat_pwe=len(own_pwe))
    if versions:
        out["sprat_pwe_version"] = ", ".join(versions)          # the SPRAT version that computed them, from the records
    return out


_WORDS = ("No", "One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine")


def _plural(n: int, one: str, many: str, start: bool = False) -> str:
    """'n things' with the noun in the right number; at the start of a sentence small counts are written in words."""
    count = _WORDS[n] if start and 0 <= n < len(_WORDS) else str(n)
    return f"{count} {one if n == 1 else many}"


def _how(org: dict[str, int]) -> str:
    if not org.get("sprat"):
        return "converted and analysed"
    if not org.get("imported"):
        return "produced and analysed"
    return "converted (the 2026 records) or produced (the regenerated ones) and analysed"


def origin_text(org: dict[str, int], sprat_meep: str = "", raw_included: bool = False, h14b: bool = False) -> str:
    """The paragraph of the README that says who produced the records."""
    legacy = ("produced on the UHeM Altay cluster in 2026 by the scripts of the computational runs: "
              f"`02_kavite.py` (Meep {LEGACY_MEEP}; with the v4 patch for the verification runs) wrote the harmonic-inversion, "
              "spectrum, reference and field-map runs, `01_bant_yapisi.py` (MPB) the band structures and `pwe_v4_feasibility.py` "
              "the plane-wave layer")
    convert = ("`sprat import-legacy` converted them into the SPRAT record schema: it translates the keys and the coded values "
               "(the mode `spektrum` becomes `spectrum`, for example) and stores the name and the SHA-256 of the original file in "
               "each record (`provenance.legacy_file`, `provenance.legacy_sha256`)")
    meep = f" with Meep {sprat_meep}" if sprat_meep else " with Meep"
    if org.get("imported") and not org.get("sprat"):
        lead = f"{org['imported']} records were" if org.get("sprat_pwe") else "The records were"
        text = (f"{lead} {legacy}. {convert[0].upper() + convert[1:]}. SPRAT converts and analyses the records and "
                "can regenerate each of them with Meep.")
    elif org.get("sprat") and not org.get("imported"):
        text = (f"The records were produced by SPRAT{meep}. Each record stores the resolved structure and parameters, the "
                "software version and the SHA-256 of the FDTD module that wrote it (`software.code_sha256`).")
    else:
        text = (f"{org.get('imported', 0)} records were {legacy}; {convert}. The other "
                f"{_plural(org.get('sprat', 0), 'record was', 'records were')} produced by SPRAT{meep}. `task.source` tells the "
                "two apart.")
    if org.get("sprat_pwe"):
        text += (f" {_plural(org['sprat_pwe'], 'plane-wave record was', 'plane-wave records were', start=True)} computed by the plane-wave "
                 f"layer of SPRAT {org.get('sprat_pwe_version') or __version__} (the complex band structure with its physical roots; no Meep).")
    if raw_included and org.get("imported"):
        text += (" The original files and the scripts that wrote them are in `raw_legacy.tar.gz`"
                 + (" and `raw_h14b.tar.gz`." if h14b else "."))
    return text


def _join(items: list[str]) -> str:
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1]


def statement(counts: dict, data_doi: str = DATA_DOI_PLACEHOLDER, software_doi: str = SOFTWARE_DOI_PLACEHOLDER,
              org: dict | None = None, raw_included: bool = False, predictions_included: bool = True,
              supplementary: bool = False) -> str:
    """The data availability statement of the paper, matched to the content of this pack."""
    org = org or dict(imported=1, sprat=0)
    held = ["the harmonic-inversion records with their geometry parameters and validity flags"]
    if counts.get("spectra"):
        held.append("the transmission spectra with their reference runs")
    if counts.get("fields"):
        held.append("the stored field maps")
    if counts.get("bands"):
        held.append("the band-structure computations")
    if predictions_included:
        held.append("the criteria fixed before the verification runs with their grading")
    held.append("the analysis outputs")
    if raw_included and org.get("imported"):
        held.append("the scripts that produced the records")
    if supplementary:
        held.append("the supplementary document, which is also available with this article")
    head = ("The data that support the findings of this study are openly available at the following URL/DOI: "
            f"https://doi.org/{data_doi}. The deposit holds {_join(held)}. ")
    if same_record(data_doi, software_doi):
        return head + (f"The same record holds SPRAT version {__version__}, the software with which the records were {_how(org)}; "
                       f"SPRAT is maintained at {REPOSITORY_URL}.")
    return head + (f"The records were {_how(org)} with SPRAT version {__version__}, archived at https://doi.org/{software_doi} "
                   f"and maintained at {REPOSITORY_URL}.")


def _repository_root() -> str | None:
    """The SPRAT checkout that holds predictions/ and campaigns/: the parent of the package (editable
    install) or the current directory (the command run from the checkout)."""
    for cand in (os.path.dirname(os.path.dirname(os.path.abspath(__file__))), os.getcwd()):
        if os.path.isdir(os.path.join(cand, "predictions", "frozen")):
            return cand
    return None


def legacy_info(raw_legacy: str | None) -> dict[str, Any]:
    """What the directory given to --raw-legacy is: its deposit version (from its zenodo_metadata.json)
    and which of its documents are present."""
    info: dict[str, Any] = dict(included=False, version=None, files=0, readme=False, superseded=[], code_licence=None)
    if not (raw_legacy and os.path.isdir(raw_legacy)):
        return info
    meta = os.path.join(raw_legacy, "zenodo_metadata.json")
    if os.path.isfile(meta):
        try:
            with open(meta, encoding="utf-8") as fh:
                info["version"] = json.load(fh).get("version")
        except (ValueError, OSError):
            pass
    info["readme"] = os.path.isfile(os.path.join(raw_legacy, "README.md"))
    info["superseded"] = [f for f in LEGACY_SUPERSEDED if os.path.isfile(os.path.join(raw_legacy, f))]
    if os.path.isfile(os.path.join(raw_legacy, "LICENSE_code.txt")):
        info["code_licence"] = "raw_legacy/LICENSE_code.txt"
    changes = sorted(f for f in os.listdir(raw_legacy) if f.startswith("CHANGES_FROM_") and f.endswith(".md"))
    info["changes"] = changes[0] if changes else None
    return info


def _legacy_what(info: dict[str, Any]) -> str:
    v = info.get("version")
    exact = f"byte for byte apart from the changes listed in its {info['changes']}" if info.get("changes") else "byte for byte"
    if v and v in UNPUBLISHED_LEGACY_VERSIONS:
        return f"the unpublished version {v} of this deposit, {exact}"
    if v:
        return f"version {v} of this deposit, {exact}"
    return "the original files of the 2026 computational runs, byte for byte"


def _legacy_description(info: dict[str, Any], version: str) -> str:
    text = (f"{_legacy_what(info)}: the per-run files with Turkish keys as the scripts wrote them, and the scripts of the systematic runs, "
            "the verification runs and the analysis.")
    if info.get("readme"):
        text += " Its README.md gives the data dictionary of the Turkish keys"
        if info.get("superseded"):
            text += f"; its own {_join(info['superseded'])} are superseded by the files of version {version}"
        text += "."
    elif info.get("superseded"):
        text += f" Its own {_join(info['superseded'])} are superseded by the files of version {version}."
    return text


_ABSOLUTE_PATH_RE = re.compile(r"(?<![\w.:/])(?:/(?:home|tmp|Users|root|mnt)/|[A-Za-z]:\\{1,2}Users\\)")


def _check_build_paths(out: str, paths: list[str | None], log=print) -> list[str]:
    """Warn about text files of the pack (archives excluded) that name a directory of the build machine."""
    needles = sorted({os.path.abspath(p) for p in paths if p} | {os.path.expanduser("~")}, key=len, reverse=True)
    needles = [n for n in needles if n.count(os.sep) >= 2]
    found = []
    for path in sorted(glob.glob(os.path.join(out, "**", "*"), recursive=True)):
        if not os.path.isfile(path) or path.endswith((".gz", ".zip", ".pdf", ".npz", ".png")):
            continue
        with open(path, "rb") as fh:
            text = fh.read().decode("utf-8", "replace")
        hit = next((n for n in needles if n in text), None)
        if hit is None:
            m = _ABSOLUTE_PATH_RE.search(text)
            hit = text[m.start():m.start() + 60].split('"')[0].split("\n")[0] if m else None
        if hit:
            rel = os.path.relpath(path, out)
            found.append(rel)
            log(f"warning: {rel} contains the absolute path {hit}; rebuild it with relative paths before publishing the pack")
    return found


# --------------------------------------------------------------------------- the pack
def build(records_dir: str, out: str, raw_legacy: str | None = None, tables: str | None = None, version: str = "3.1.0",
          include: list[str] | None = None, data_doi: str = DATA_DOI_PLACEHOLDER,
          software_doi: str = SOFTWARE_DOI_PLACEHOLDER, raw_h14b: str | None = None, arxiv_id: str | None = None,
          log=print) -> dict[str, Any]:
    data_doi, software_doi = normalise_doi(data_doi), normalise_doi(software_doi)
    for f in include or []:
        if not os.path.isfile(f):
            raise FileNotFoundError(f)
    if sprat_archive_of([os.path.basename(f) for f in include or []]) and is_placeholder(software_doi):
        software_doi = data_doi   # the software is archived in this record: one DOI for both
    os.makedirs(out, exist_ok=True)
    recs = collect.load(records_dir)
    rows = collect.table(recs)
    counts = records.counts(rows)
    org = origin(recs)
    records.write_table(rows, os.path.join(out, "records.csv"), os.path.join(out, "records.jsonl"))
    with open(os.path.join(out, "counts.json"), "w", encoding="utf-8") as fh:
        json.dump(counts, fh, indent=1)
    groups = {"harminv_records": ["harminv"], "spectra_records": ["spectrum", "reference"], "field_records": ["field"],
              "bands_records": ["bands"], "pwe_records": ["pwe"]}
    inventory = {}
    for name, modes in groups.items():
        files = [r["_file"] for r in recs if r["params"]["run"]["mode"] in modes]
        if name == "field_records":
            files += [os.path.join(os.path.dirname(r["_file"]), r["result"]["field_file"]) for r in recs
                      if r["params"]["run"]["mode"] == "field" and r["result"].get("field_file")]
        if files:
            inventory[name] = _tar(os.path.join(out, name + ".tar.gz"), sorted(files), name)
    raw = legacy_info(raw_legacy)
    if raw_legacy and os.path.isdir(raw_legacy):
        files = [p for p in glob.glob(os.path.join(raw_legacy, "**", "*"), recursive=True) if os.path.isfile(p)]
        with tarfile.open(os.path.join(out, "raw_legacy.tar.gz"), "w:gz") as tf:
            for f in sorted(files):
                tf.add(f, arcname=os.path.join("raw_legacy", os.path.relpath(f, raw_legacy)), filter=_owner_free)
        inventory["raw_legacy"] = len(files)
        raw.update(included=True, files=len(files))
    elif raw_legacy:
        log(f"warning: --raw-legacy {raw_legacy} is not a directory; raw_legacy.tar.gz is not written")
    h14b_included = False
    if raw_h14b and os.path.isdir(raw_h14b):
        inventory["raw_h14b"] = _tar_dir(os.path.join(out, "raw_h14b.tar.gz"), raw_h14b, "raw_h14b")
        h14b_included = True
    elif raw_h14b:
        log(f"warning: --raw-h14b {raw_h14b} is not a directory; raw_h14b.tar.gz is not written")
    # analysis outputs
    if tables and os.path.isdir(tables):
        shutil.copytree(tables, os.path.join(out, "analysis"), dirs_exist_ok=True)
    else:
        pipeline.analyze(records_dir, os.path.join(out, "analysis"), log=log)
    root = _repository_root()
    manuscript_registry = False
    src = os.path.join(root, MANUSCRIPT_REGISTRY) if root else None
    if src and os.path.isfile(src):
        shutil.copy2(src, os.path.join(out, "analysis", "manuscript_registry.json"))
        manuscript_registry = True
    else:
        log(f"warning: {MANUSCRIPT_REGISTRY} not found in the repository checkout; "
            "analysis/manuscript_registry.json is not included in the deposit")
    analysis_files = sorted(os.path.relpath(p, os.path.join(out, "analysis"))
                            for p in glob.glob(os.path.join(out, "analysis", "**", "*"), recursive=True) if os.path.isfile(p))
    # frozen predictions
    predictions_included = False
    if root:
        shutil.copytree(os.path.join(root, "predictions"), os.path.join(out, "predictions"), dirs_exist_ok=True)
        shutil.copy2(predictions.CRITERIA_FILE, os.path.join(out, "predictions", "predictions_criteria.json"))
        predictions_included = True
    else:
        log("warning: predictions/frozen not found beside the package or in the current directory (run the command from "
            "the repository checkout); the frozen prediction files are not included in the deposit")
    # Zenodo stores files without folders: the two folders travel as archives
    folder_counts = {}
    for name in ("analysis", "predictions"):
        d = os.path.join(out, name)
        if os.path.isdir(d):
            folder_counts[name] = _tar_dir(os.path.join(out, name + ".tar.gz"), d, name)
            shutil.rmtree(d)
    with open(os.path.join(out, "parameters_manifest.json"), "w", encoding="utf-8") as fh:
        json.dump(parameters_manifest(recs, version), fh, indent=1)
    with open(os.path.join(out, "LICENSE.txt"), "w", encoding="utf-8") as fh:
        fh.write(licence_text(raw, org, sprat_archive=sprat_archive_of([os.path.basename(f) for f in include or []]),
                              h14b=h14b_included))
    included = []
    for f in include or []:
        shutil.copy2(f, os.path.join(out, os.path.basename(f)))
        included.append(os.path.basename(f))
    supplementary = any(f.lower().startswith("supplementary") and f.lower().endswith(".pdf") for f in included)
    sprat_meep = ", ".join(sorted({str(r.get("software", {}).get("meep") or "") for r in recs
                                   if r.get("task", {}).get("source") == "sprat"} - {""}))
    with open(os.path.join(out, "zenodo_metadata.json"), "w", encoding="utf-8") as fh:
        json.dump(zenodo_metadata(counts, version, raw_included=raw["included"], included=included, data_doi=data_doi,
                                  software_doi=software_doi, org=org, raw=raw, manuscript_registry=manuscript_registry,
                                  h14b=h14b_included, arxiv_id=arxiv_id),
                  fh, indent=1)
    with open(os.path.join(out, "DATA_AVAILABILITY.md"), "w", encoding="utf-8") as fh:
        fh.write(data_availability(counts, data_doi=data_doi, software_doi=software_doi, org=org, raw_included=raw["included"],
                                   predictions_included=predictions_included, supplementary=supplementary))
    with open(os.path.join(out, "README.md"), "w", encoding="utf-8") as fh:
        fh.write(readme(counts, inventory, version, included, data_doi=data_doi, software_doi=software_doi, org=org,
                        raw=raw, analysis_files=analysis_files, manuscript_registry=manuscript_registry,
                        predictions_included=predictions_included, sprat_meep=sprat_meep, folder_counts=folder_counts,
                        arxiv_id=arxiv_id, n_rows=len(rows)))
    _check_build_paths(out, [records_dir, tables, raw_legacy, raw_h14b, out, os.getcwd(), root], log=log)
    # manifest
    lines = []
    for path in sorted(glob.glob(os.path.join(out, "**", "*"), recursive=True)):
        if os.path.isfile(path) and os.path.basename(path) != "MANIFEST_sha256.txt":
            lines.append(f"{_sha(path)}  {os.path.relpath(path, out)}")
    with open(os.path.join(out, "MANIFEST_sha256.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    log(f"deposit written to {out}: {counts}; {len(lines)} files in the manifest")
    return dict(counts=counts, inventory=inventory, files=len(lines), origin=org, manuscript_registry=manuscript_registry)


def parameters_manifest(recs: list[dict], version: str) -> dict[str, Any]:
    ref = collect.best(recs, mode="harminv", radius=0.06, row=4, cladding_rows=12, resolution=24, termination_x="pml", guide_periods=25,
                       n=1.33, dx=0.0, dy=0.0, no_rows=True)
    s = ref["structure"] if ref else {}
    p = ref["params"] if ref else {}
    return dict(
        version=version, built=time.strftime("%Y-%m-%d %H:%M:%S"), software=f"sprat {__version__}",
        model=dict(lattice="square", polarisation="TM (E_z)", units="a = c = 1, f = a/lambda",
                   rod_radius_over_a=s.get("lattice", {}).get("rod_radius"), rod_eps=s.get("lattice", {}).get("rod_eps"),
                   a_nm=s.get("lattice", {}).get("a_nm"), target_wavelength_nm=s.get("lattice", {}).get("target_wavelength_nm"),
                   n_analyte_reference=1.33, gap_edges_by_index=records.GAP_EDGES_MPB),
        reference_geometry=s, reference_params=p,
        numerics=dict(fdtd="Meep 1.30 (conda-forge pymeep, serial)", courant=0.5, subpixel_smoothing=True,
                      harminv="Meep Harminv on E_z at the defect centre + (0.13, 0.07) a",
                      Q_limit_rule="pi f_cen t, a sampling rule (conservative: harmonic inversion resolves decays longer than the signal)",
                      validity="margin Q_lim/Q >= 1 and gap position 0.2-0.8 (gap edges interpolated in the MPB table)"),
        analysis=dict(leak_removal="1/Q_barrier = 1/Q - 1/Q_top(c), Q_top(c) = Q_top(8) exp(kappa_top (c - 8))",
                      kappa="OLS slope of ln Q_w(N) over N = 3..6 after leak removal and division of the PML records by the measured G",
                      G="Q_PML / Q_absorber, both leak-removed, at the same N",
                      eta="Hellmann-Feynman: sum (1 - f) eps_a |E_z|^2 / sum <eps> |E_z|^2 with <eps> of the subpixel smoothing reconstructed "
                          "on the E_z grid (the stored map is Meep's harmonic mean of the eigenvalues of the smoothed tensor); "
                          "Q_abs = lambda_r / (2 k S)",
                      complex_bands="physical roots only: eigenvector weight on the two outermost G_y rings below 1e-3 and the replica "
                                    "with |Re k_y| <= pi/a",
                      admission="fits and statistics use the admitted records only; excluded records are shown, marked, in the figures",
                      selection="records are selected by parameters; the longest harmonic-inversion signal wins; the figures of the systematic runs exclude the verification runs"))


def licence_text(raw: dict | None = None, org: dict | None = None, sprat_archive: str | None = None, h14b: bool = False) -> str:
    raw, org = raw or {}, org or dict(imported=1, sprat=0)
    L = ["Creative Commons Attribution 4.0 International (CC BY 4.0)", "",
         "The data in this record are released under CC BY 4.0: https://creativecommons.org/licenses/by/4.0/",
         "Attribution: Hasan Oguz (ORCID 0000-0001-7484-4415), Istanbul Okan University and Pamukkale University.", ""]
    if raw.get("included"):
        where = f"; its text is {raw['code_licence']} inside the archive" if raw.get("code_licence") else ""
        L += ["The scripts inside raw_legacy.tar.gz (the code of the 2026 systematic runs, verification runs and analysis) are released under the "
              f"MIT licence{where}.", ""]
    if h14b:
        L += ["The scripts inside raw_h14b.tar.gz (the job, the comparison and the amended grading of the absorber spectra) are released "
              "under the MIT licence.", ""]
    if sprat_archive:
        L += [f"SPRAT, the software with which the records were {_how(org)}, is included in this record as {sprat_archive} and is "
              f"released under the MIT licence (its LICENSE file); it is maintained at {REPOSITORY_URL}"]
    else:
        L += [f"SPRAT, the software that {_how(org)} the records, is released separately under the MIT licence: {REPOSITORY_URL}"]
    return "\n".join(L) + "\n"


def zenodo_metadata(counts: dict, version: str, raw_included: bool = False, included: list[str] | None = None,
                    data_doi: str = DATA_DOI_PLACEHOLDER, software_doi: str = SOFTWARE_DOI_PLACEHOLDER,
                    org: dict | None = None, raw: dict | None = None, manuscript_registry: bool = False, h14b: bool = False,
                    arxiv_id: str | None = None) -> dict[str, Any]:
    org = org or dict(imported=1, sprat=0)
    raw = dict(raw or {}, included=raw_included)
    included = included or []
    sprat_archive = sprat_archive_of(included)
    if same_record(data_doi, software_doi) or sprat_archive:
        software = f"SPRAT {__version__} (included in this record" + (f" as {sprat_archive}" if sprat_archive else "") + f"; {REPOSITORY_URL})"
    else:
        software = f"SPRAT {__version__} ({REPOSITORY_URL}" + ("" if is_placeholder(software_doi) else f"; https://doi.org/{software_doi}") + ")"
    if org.get("imported") and not org.get("sprat"):
        lead = f"{org['imported']} records were" if org.get("sprat_pwe") else "The records were"
        origin_sentence = (f"{lead} produced on the UHeM Altay cluster by the 2026 scripts of the computational runs (02_kavite.py, with the v4 "
                           f"patch for the verification runs; Meep {LEGACY_MEEP}) and converted into the SPRAT record schema (English keys, "
                           f"sprat-record-1.0) with sprat import-legacy; {software} converts and analyses them and can regenerate each "
                           "record with Meep.")
    elif org.get("sprat") and not org.get("imported"):
        origin_sentence = f"The records (English keys, schema sprat-record-1.0) were produced and are analysed by {software}."
    else:
        origin_sentence = (f"{org.get('imported', 0)} records were produced on the UHeM Altay cluster by the 2026 scripts of the computational runs and "
                           f"converted into the SPRAT record schema with sprat import-legacy; "
                           f"{_plural(org.get('sprat', 0), 'was', 'were')} produced by {software}, which analyses all of them.")
    if org.get("sprat_pwe"):
        origin_sentence += (f" {_plural(org['sprat_pwe'], 'plane-wave record was', 'plane-wave records were', start=True)} computed by the "
                            f"plane-wave layer of SPRAT {org.get('sprat_pwe_version') or __version__}.")
    runs = [f"{counts['harminv']} harmonic-inversion records ({counts['valid']} admitted by the validity criterion of the paper)"]
    if counts.get("spectra"):
        sa = counts.get("spectra_absorber", 0)
        runs.append(f"{counts['spectra']} transmission spectra with their cavity-less reference runs"
                    + (f" ({counts['spectra'] - sa} with the guide ending in the perfectly matched layer, {sa} with the guide continued "
                       "into an absorber)" if sa else ""))
    if counts.get("fields"):
        runs.append(f"{counts['fields']} field maps")
    if counts.get("bands"):
        runs.append(f"{counts['bands']} band computations")
    if counts.get("pwe"):
        runs.append("the plane-wave layer")
    contents = [_join(runs), "the flat table with the validity flags",
                "the analysis outputs with the numbers registry" + (" and the registry of the manuscript" if manuscript_registry else ""),
                "the criterion sets fixed before the runs" + (", one of them amended while its runs were under way," if h14b else "") + " and their grading"]
    for f in included:
        if f.lower().startswith("supplementary"):
            contents.append(f"the supplementary document of the paper ({f})")
        elif f == sprat_archive:
            contents.append(f"the source archive of the software SPRAT {__version__} ({f})")
        else:
            contents.append(f)
    description = (f"Numerical data set of the paper \"{PAPER_TITLE}\": a point-defect microcavity side-coupled to a W1 waveguide in a "
                   "two-dimensional square lattice of silicon rods immersed in an aqueous analyte (TM polarisation; Meep FDTD, MPB "
                   f"and plane-wave expansion). {origin_sentence} The record holds {'; '.join(contents[:-1])}; and {contents[-1]}.")
    if raw_included:
        description += (f" raw_legacy.tar.gz holds {_legacy_what(raw)}: the original files with Turkish keys"
                        + (", whose data dictionary is its README.md," if raw.get("readme") else "") + " and the scripts that wrote them.")
    if h14b:
        description += (" raw_h14b.tar.gz holds the original files of the absorber-terminated spectra with their job files, the criterion "
                        "fixed before those runs, its amendment and the grading.")
    if arxiv_id:
        description += f" The record supplements the preprint arXiv:{arxiv_id}."
    description += " README.md gives the data dictionary and the commands that reproduce the analysis."
    licences = [dict(id="CC-BY-4.0", applies_to="the data, the analysis outputs and the documents of this record")]
    notes = NOTES_FUNDING + " The data are released under CC BY 4.0"
    mit = (([f"the software SPRAT ({sprat_archive})"] if sprat_archive else []) + (["the scripts inside raw_legacy.tar.gz"] if raw_included else [])
           + (["the scripts inside raw_h14b.tar.gz"] if h14b else []))
    if mit:
        licences.append(dict(id="MIT", applies_to=_join(mit)))
        notes += "; " + _join(mit) + " under the MIT licence"
    notes += "."
    if raw_included and raw.get("version") in UNPUBLISHED_LEGACY_VERSIONS:
        notes += (f" raw_legacy.tar.gz holds the unpublished version {raw['version']} of this deposit byte for byte"
                  + (f", apart from the changes listed in its {raw['changes']}." if raw.get("changes") else "."))
    meta: dict[str, Any] = dict(
        upload_type="dataset", title=DATA_TITLE, version=version, creators=[dict(AUTHOR)], description=description,
        keywords=["photonic crystal", "microcavity", "quality factor", "FDTD", "Meep", "plane-wave expansion", "refractive-index sensing",
                  "evanescent Bloch modes", "harmonic inversion", "silicon rods", "pre-specified criteria", "reproducibility", "adiabatic absorber"],
        license="CC-BY-4.0", licenses=licences, access_right="open", language="eng", notes=notes)
    if not is_placeholder(data_doi):
        meta["reserved_doi"] = data_doi
    related = []
    if not is_placeholder(software_doi) and not same_record(data_doi, software_doi):
        related.append(dict(relation="isCompiledBy", identifier=software_doi, scheme="doi", resource_type="software"))
    if arxiv_id:
        related.append(dict(relation="isSupplementTo", identifier=f"arXiv:{arxiv_id}", scheme="arxiv", resource_type="publication-preprint"))
    if related:
        meta["related_identifiers"] = related
    return meta


def data_availability(counts: dict, data_doi: str = DATA_DOI_PLACEHOLDER, software_doi: str = SOFTWARE_DOI_PLACEHOLDER,
                      org: dict | None = None, raw_included: bool = False, predictions_included: bool = True,
                      supplementary: bool = False) -> str:
    return ("# Data availability statement\n\n"
            "The statement printed in the article, matched to the content of this deposit:\n\n> "
            + statement(counts, data_doi, software_doi, org=org, raw_included=raw_included,
                        predictions_included=predictions_included, supplementary=supplementary) + "\n\n"
            "## Counts in this deposit\n\n"
            "| quantity | value | source in this deposit |\n|---|---|---|\n"
            f"| harmonic-inversion records | {counts['harminv']} | rows of records.csv with mode harminv |\n"
            f"| admitted (margin >= 1, gap position 0.2-0.8) | {counts['valid']} | rows with valid = True |\n"
            f"| excluded, margin below 1 | {counts['excluded_margin']} | exclusion_reason = Q > Q_lim |\n"
            f"| excluded, outside the gap window | {counts['excluded_gap']} | exclusion_reason = gap position outside 0.2-0.8 |\n"
            f"| admitted, holding a mode of the finite waveguide (preliminary sweep: N_sep = 1 or r_d >= 0.13a) | "
            f"{counts.get('waveguide_mode_admitted', 0)} | rows with waveguide_mode = True |\n"
            f"| of these, in the dead zone (r_d >= 0.14a) | {counts.get('dead_zone_admitted', 0)} "
            f"| rows with valid = True and radius >= 0.14 |\n"
            f"| admitted, holding a cavity mode | {counts['valid'] - counts.get('waveguide_mode_admitted', 0)} "
            f"| rows with valid = True and waveguide_mode = False |\n"
            f"| clearing the margin alone | {counts['resolved']} | rows with resolved = True |\n"
            f"| spectra / references | {counts['spectra']} / {counts['references']} | spectra_records.tar.gz |\n"
            + (f"| of these, with the guide continued into the absorber | {counts['spectra_absorber']} / "
               f"{counts.get('references_absorber', 0)} | spectra_records.tar.gz, structure.cell.termination_x = absorber |\n"
               if counts.get("spectra_absorber") else "") +
            f"| field maps | {counts['fields']} | field_records.tar.gz |\n"
            f"| band computations | {counts['bands']} | bands_records.tar.gz |\n")


def readme(counts: dict, inventory: dict, version: str, included: list[str] | None = None,
           data_doi: str = DATA_DOI_PLACEHOLDER, software_doi: str = SOFTWARE_DOI_PLACEHOLDER, org: dict | None = None,
           raw: dict | None = None, analysis_files: list[str] | None = None, manuscript_registry: bool = False,
           predictions_included: bool = True, sprat_meep: str = "", folder_counts: dict | None = None,
           arxiv_id: str | None = None, n_rows: int | None = None) -> str:
    org = org or dict(imported=1, sprat=0)
    raw = dict(raw or {}, included="raw_legacy" in inventory)
    record_files = {"harminv_records": "the harmonic-inversion records (JSON, English keys)",
                    "spectra_records": "the transmission spectra and the cavity-less reference runs",
                    "field_records": "the field-map records and their .npz maps",
                    "bands_records": "the band computations and the lattice-constant calibration",
                    "pwe_records": "the plane-wave layer"}
    L = [f"# {DATA_TITLE}", "",
         f"- **Version:** {version}",
         "- **Author:** Hasan Oguz (ORCID 0000-0001-7484-4415), Istanbul Okan University and Pamukkale University",
         f"- **This record:** https://doi.org/{data_doi}",
         (f"- **Software:** SPRAT {__version__}, in this record" + (f" (`{sprat_archive_of(included)}`)" if sprat_archive_of(included) else "")
          + f", maintained at {REPOSITORY_URL}" if (same_record(data_doi, software_doi) or sprat_archive_of(included))
          else f"- **Software:** SPRAT {__version__}, https://doi.org/{software_doi}, maintained at {REPOSITORY_URL}"),
         "- **Licence:** data CC BY 4.0 (LICENSE.txt)" + ("; the scripts inside raw_legacy.tar.gz MIT" if raw["included"] else "")
         + ("; the scripts inside raw_h14b.tar.gz MIT" if "raw_h14b" in inventory else "")
         + ("; SPRAT MIT" if sprat_archive_of(included) else "; SPRAT has its own record under the MIT licence"),
         f"- **Built:** {time.strftime('%Y-%m-%d')} with sprat {__version__}"]
    if arxiv_id:
        L.append(f"- **Preprint:** arXiv:{arxiv_id}")
    L += ["",
         "## Where the records come from", "", origin_text(org, sprat_meep, raw["included"], h14b="raw_h14b" in inventory), "",
         "## Contents", "", "| file | content | count |", "|---|---|---|",
         "| records.csv, records.jsonl | one row per record with the validity flags | "
         + (f"{n_rows} rows, {counts['harminv']} harmonic inversions |" if n_rows else f"{counts['harminv']} harmonic-inversion rows |")]
    for k, v in inventory.items():
        if k == "raw_legacy":
            L.append(f"| raw_legacy.tar.gz | {_legacy_description(raw, version)} | {v} file{'s' if v != 1 else ''} |")
        elif k == "raw_h14b":
            L.append("| raw_h14b.tar.gz | the seven transmission spectra repeated with the guide continued into the absorber and their reference "
                     "runs (Turkish keys, byte for byte), the job files, the criterion fixed before the runs, its amendment and the grading; "
                     f"its README.md describes them | {v} files |")
        else:
            L.append(f"| {k}.tar.gz | {record_files.get(k, k.replace('_', ' '))} | {v} file{'s' if v != 1 else ''} |")
    fc = folder_counts or {}
    shown = [f for f in (analysis_files or []) if f != "manuscript_registry.json"]
    n_an = fc.get("analysis")
    L.append("| analysis.tar.gz | unpacks to `analysis/`: the outputs of `sprat analyze` and `sprat audit`"
             + (": " + ", ".join(shown) if shown else "")
             + ("; and `manuscript_registry.json`, the numbers registry of the manuscript: every number quoted in the paper "
                "with the record or the derivation behind it, or its literature source" if manuscript_registry else "")
             + (f" | {n_an} files |" if n_an else " | - |"))
    if predictions_included:
        n_pr = fc.get("predictions")
        L.append("| predictions.tar.gz | unpacks to `predictions/`: the criteria fixed before the runs, byte for byte, with the "
                 "English criteria file and a README" + (f" | {n_pr} files |" if n_pr else " | - |"))
    L += ["| parameters_manifest.json | every constant of the model, the numerics and the analysis | - |",
          "| counts.json | the counts of the records, those of the table in DATA_AVAILABILITY.md among them | - |",
          "| DATA_AVAILABILITY.md | the data availability statement of the paper and the counts of this deposit | - |",
          "| LICENSE.txt, zenodo_metadata.json | licence, metadata of this record | - |"]
    for f in included or []:
        what = ("the supplementary document of the paper, also available with the article"
                if f.lower().startswith("supplementary") else
                (f"the source archive of SPRAT {__version__}, the software with which the records were converted and analysed "
                 "(MIT licence); `pip install` it or unpack it and follow its README" if f == sprat_archive_of(included)
                 else "additional file"))
        L.append(f"| {f} | {what} | - |")
    L += ["| MANIFEST_sha256.txt | SHA-256 of every file above | - |", "",
          "## Model", "",
          "Square lattice of silicon rods (r = 0.20a, eps = 11.9025) in an analyte of index n_a; W1 waveguide; point defect of radius r_d at "
          "N_sep rows from the guide; TM polarisation; a = 481.4 nm for the 1550 nm target. FDTD: Meep 1.30, resolution 24 in the reference "
          "geometry, Courant 0.5, subpixel smoothing, PML one period, mirror symmetry in x when the defect is on the axis. Q by harmonic "
          "inversion; the sampling rule sets Q_lim = pi f_cen t and the margin Q_lim/Q, and a record is admitted when the margin is at "
          "least 1 and the mode lies between 20 and 80 per cent of the gap (a conservative rule: harmonic inversion resolves decays "
          "longer than the signal).", "",
          "## Record schema (sprat-record-1.0)", "", "| key | meaning |", "|---|---|"]
    L += [f"| `{k}` | {v} |" for k, v in DATA_DICTIONARY]
    names = [k for k in inventory if k not in ("raw_legacy", "raw_h14b")]
    arch = sprat_archive_of(included)
    where = f"SPRAT (`{arch}` in this record, maintained at {REPOSITORY_URL})" if arch else f"SPRAT ({REPOSITORY_URL})"
    install = f"Install it from `{arch}`" if arch else "Install it"
    L += ["", "## Reproducing", "",
          f"{where} reads the records of this deposit directly. {install} (its README gives the commands), "
          "unpack the record archives into one directory and run the analysis:", "",
          "```bash",
          "mkdir records",
          f"for f in {' '.join(names)}; do tar -xzf $f.tar.gz -C records --strip-components=1; done",
          "tar -xzf analysis.tar.gz; tar -xzf predictions.tar.gz",
          "sprat analyze records -o tables",
          "sprat audit records --expected analysis/numbers_registry.json -o tables",
          "sprat figures records -o figures --tables tables",
          "```", "",
          "Zenodo stores files without folders, so `analysis/` and `predictions/` travel as the two archives above.", "",
          "`sprat analyze` recomputes the analysis outputs of `analysis/` from the records: the derived numbers, the numbers registry "
          "and the grading of the criteria fixed before the runs. `sprat audit` recomputes the registered numbers of the paper (the numbers "
          "registry) from the records and compares them, row by row, with the registry of this deposit; run on regenerated records, "
          "it shows how far a regeneration departs from the deposit. `sprat figures` draws the figures of the paper from the records.",
          "",
          "Each record stores the resolved structure and parameters of its run, so `sprat run-one --task <record.json> --out "
          "<new.json>` runs it again with Meep, and `sprat reproduce campaigns/manuscript` ("
          + (f"the parameter files are in `{arch}`" if arch else "in the SPRAT repository") + ") plans and runs "
          "the whole set on a workstation."]
    if manuscript_registry:
        L += ["", "`analysis/manuscript_registry.json` lists every number of the manuscript with its source. "
                  "`sprat audit records --expected analysis/manuscript_registry.json -o tables` compares its record-derived rows "
                  "with the numbers recomputed from the records and lists the literature values and references apart."]
    if raw["included"] and org.get("imported"):
        L += ["", "The original files in `raw_legacy.tar.gz` are the input of `sprat import-legacy`: `tar -xzf raw_legacy.tar.gz` "
                  "followed by `sprat import-legacy raw_legacy -o records` converts them again"
                  + ("; `tar -xzf raw_h14b.tar.gz` and `sprat import-legacy raw_h14b -o records` add the absorber-terminated spectra."
                     if "raw_h14b" in inventory else ".")]
    L.append("")
    return "\n".join(L)


def build_cli(records_dir: str, out: str, raw_legacy: str | None = None, tables: str | None = None, version: str = "3.1.0",
              include: list[str] | None = None, data_doi: str = DATA_DOI_PLACEHOLDER,
              software_doi: str = SOFTWARE_DOI_PLACEHOLDER, raw_h14b: str | None = None, arxiv_id: str | None = None) -> None:
    build(records_dir, out, raw_legacy=raw_legacy, tables=tables, version=version, include=include, data_doi=data_doi,
          software_doi=software_doi, raw_h14b=raw_h14b, arxiv_id=arxiv_id)


__all__ = ["build", "build_cli", "parameters_manifest", "zenodo_metadata", "data_availability", "readme", "statement",
           "licence_text", "normalise_doi", "is_placeholder", "origin", "legacy_info", "REPOSITORY_URL",
           "DATA_DOI_PLACEHOLDER", "SOFTWARE_DOI_PLACEHOLDER", "DATA_TITLE", "PAPER_TITLE"]
