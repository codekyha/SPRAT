"""The manuscript registry of version 10 audited against the registry SPRAT recomputes from the data record: the imported
2026 deposit (raw_legacy), the absorber-terminated spectra (raw_h14b) and the plane-wave record of SPRAT 1.1.0
(tests/fixtures). Every record-derived row agrees within its printed precision, every superseded statement is
consistent, literature values, references and notes are listed as not record-derived and never counted as OK.

The registry of version 9, audited on raw_legacy alone, deviates only on the two quantities SPRAT 1.2.0 renames; the
registries of manuscript versions 7 and 8 and of the unpublished 2.0.2 deposit deviate exactly on the rows that version
9 corrects (V9_CORRECTED) and on the renamed row they carry, and on nothing else. The outputs carry no absolute path.
Skipped unless the deposits are available: SPRAT_LEGACY_DEPOSIT (as for test_legacy_import) for everything, and
SPRAT_H14B_DEPOSIT (the unpacked raw_h14b.tar.gz) in addition for the registry of version 10."""
import json
import os
import shutil

import pytest

HERE = os.path.dirname(__file__)
ROOT = os.path.dirname(HERE)
EXPECTED = os.path.join(ROOT, "campaigns", "manuscript", "expected")
MANUSCRIPT = os.path.join(EXPECTED, "numbers_registry_manuscript_v7.json")
MANUSCRIPT_V8 = os.path.join(EXPECTED, "numbers_registry_manuscript_v8.json")
MANUSCRIPT_V9 = os.path.join(EXPECTED, "numbers_registry_manuscript_v9.json")
MANUSCRIPT_V10 = os.path.join(EXPECTED, "numbers_registry_manuscript_v10.json")
MANUSCRIPT_V10_2 = os.path.join(EXPECTED, "numbers_registry_manuscript_v10_2.json")
PWE_FIXTURE = os.path.join(HERE, "fixtures", "pwe_na1.3300_r0.2.json")
NOT_RECORD = {"literature_v7", "references_added_v6", "references_added_v7", "references_removed_v7", "dataset_note",
              "references_added_v9", "references_removed_v9", "literature_v9", "corrected_v9_note", "corrected_v10_note",
              "literature_v10_2", "references_added_v10_2", "derived_v10_2_note"}
# rows of the earlier registries that SPRAT 1.1.0 no longer reproduces, because manuscript version 9 corrects them:
# the fits over all points (now admitted records only), the loss to water (Q_abs = lambda_r / (2 k S)), the channel
# of the complex band structure (physical roots only; the slower root was a basis artefact), and the analyte energy
# fraction (the stored map is not the permittivity E_z sees); the removed quantities are MISSING
V9_CORRECTED = {
    ("DEVIATION", "derived", "D, effective reflector distance"),
    ("DEVIATION", "derived", "FOM at N_sep = 5, in water"),
    ("DEVIATION", "derived", "Q_abs in water"),
    ("DEVIATION", "derived", "Q_tot in water"),
    ("DEVIATION", "derived", "Q_w at the coupling optimum"),
    ("DEVIATION", "derived", "T_min in water"),
    ("DEVIATION", "derived", "absorber scale factor c"),
    ("DEVIATION", "derived", "kappa_pred, zone-edge channel"),
    ("DEVIATION", "derived", "leak-corrected per-row step of the n_cl = 8 series, radius by radius"),
    ("DEVIATION", "derived", "mean per-row step of the n_cl = 8 series, raw"),
    ("DEVIATION", "derived", "suppression of the Q(r_d) oscillation by the absorber"),
    ("DEVIATION", "field_map", "per-row barrier decay at k_x = beta, rows 3 to 2"),
    ("MISSING", "derived", "S_th from the eta_a bracket"),
    ("MISSING", "derived", "eta_a, analyte energy fraction"),
    ("MISSING", "derived", "kappa_pred, slowest-pair channel"),
    ("MISSING", "field_map", "S_th from the lower bound"),
    ("MISSING", "field_map", "S_th from the lower bound, n_a = 1.45"),
    ("MISSING", "field_map", "energy fraction in mixed pixels"),
    ("MISSING", "field_map", "eta_a, rigorous lower bound"),
}
# the quantities SPRAT 1.2.0 renames: every earlier registry names the first, version 9 also the second
V10_RENAMED_EARLY = {("MISSING", "dataset", "records from the parameter campaign")}
V10_RENAMED = V10_RENAMED_EARLY | {("MISSING", "corrected_v9", "eta_a, threshold mask of the campaign")}
V9_CORRECTED_V7 = V9_CORRECTED | {("DEVIATION", "derived_v7", "coupling optimum"),
                                  ("DEVIATION", "vertical", "Q_perp leaving Q_tot within 3 per cent"),
                                  ("DEVIATION", "vertical", "Q_tot with no vertical channel")}
V9_CORRECTED_V8 = V9_CORRECTED_V7 | {("MISSING", "corrected_v8", "largest deviation of the quadratic law from the thirteen-point "
                                                                 "fine sweep in the reference geometry")}


def _deposit():
    for cand in (os.environ.get("SPRAT_LEGACY_DEPOSIT", ""), os.path.join(HERE, "data", "zenodo_pack_v2")):
        if cand and os.path.isfile(os.path.join(cand, "records_all.csv")):
            return cand
    return None


def _h14b():
    cand = os.environ.get("SPRAT_H14B_DEPOSIT", "")
    return cand if cand and os.path.isdir(os.path.join(cand, "v4_spectra")) else None


pytestmark = pytest.mark.skipif(_deposit() is None, reason="full deposit not available (set SPRAT_LEGACY_DEPOSIT)")
needs_h14b = pytest.mark.skipif(_h14b() is None, reason="absorber spectra not available (set SPRAT_H14B_DEPOSIT)")


def _analyze(d, with_h14b):
    from sprat import legacy
    from sprat.analysis import pipeline
    legacy.import_deposit(_deposit(), str(d / "records"), log=lambda *a: None)
    if with_h14b:
        legacy.import_deposit(_h14b(), str(d / "records"), log=lambda *a: None)
    shutil.copyfile(PWE_FIXTURE, d / "records" / os.path.basename(PWE_FIXTURE))
    return pipeline.analyze(str(d / "records"), str(d / "tables"), log=lambda *a: None)


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    """The records of manuscript version 9: raw_legacy and the plane-wave record."""
    d = tmp_path_factory.mktemp("audit")
    return _analyze(d, False), d


@pytest.fixture(scope="module")
def run_v10(tmp_path_factory):
    """The records of manuscript version 10: raw_legacy, raw_h14b and the plane-wave record."""
    if _h14b() is None:
        pytest.skip("absorber spectra not available (set SPRAT_H14B_DEPOSIT)")
    d = tmp_path_factory.mktemp("audit_v10")
    return _analyze(d, True), d


def _not_ok(cmp):
    return {(r["verdict"], r["section"], r["quantity"]) for r in cmp["rows"] if r["verdict"] != "OK"}


@needs_h14b
def test_manuscript_registry_every_record_row(run_v10):
    from sprat.analysis import audit
    A, _ = run_v10
    cmp = audit.compare_registries(A["registry"], json.load(open(MANUSCRIPT_V10, encoding="utf-8")))
    bad = [r for r in cmp["rows"] if r["verdict"] != "OK"]
    assert not bad, bad
    assert cmp["counts"] == {"OK": 223}
    assert len(cmp["superseded"]) == 21 and all(r["verdict"] == "CONSISTENT" for r in cmp["superseded"])
    assert cmp["not_record_derived_count"] == 29 and cmp["ok"]
    sections = {r["section"] for r in cmp["rows"]}
    assert {"measured", "derived", "second_radius", "field_map", "dataset", "vertical", "barrier_row_v7", "absorption_v7", "derived_v7",
            "corrected_v8", "corrected_v9", "termination", "corrected_v10"} <= sections
    assert not sections & NOT_RECORD
    rows = {r["quantity"]: r for r in cmp["rows"]}
    for q in ("resonance shift over the barrier-row series 0.18a to 0.22a", "ratio of the transverse-displacement coefficient to the row coefficient",
              "reproduction of the N_sep = 5 point by repeats at higher margin", "reproduction of the reference point by a repeat at higher margin",
              "Q_perp leaving Q_tot within 3 per cent", "reflectionless linewidth count for a 1 nm radius error", "coupling optimum",
              "verdict by the criterion fixed before the runs (D-22)", "verdict by the amended criterion (A1)",
              "share of the harmonic-inversion records admitted and holding a cavity mode", "transmission spectra"):
        assert rows[q]["verdict"] == "OK"
    assert rows["verdict by the amended criterion (A1)"]["recomputed"] == "pass"
    assert rows["transmission spectra"]["recomputed"] == 14


def test_registry_of_version_9_on_its_own_records(run):
    """On raw_legacy alone the registry of version 9 is reproduced apart from the two renamed quantities."""
    from sprat.analysis import audit
    A, _ = run
    cmp = audit.compare_registries(A["registry"], json.load(open(MANUSCRIPT_V9, encoding="utf-8")))
    assert _not_ok(cmp) == V10_RENAMED and cmp["counts"] == {"OK": 190, "MISSING": 2}
    # one reason of version 9 is reworded by SPRAT 1.2.0 (the original plane-wave layer); its values are unchanged
    assert cmp["superseded_counts"] == {"CONSISTENT": 20, "DEVIATION": 1}
    dev = [r for r in cmp["superseded"] if r["verdict"] != "CONSISTENT"][0]
    assert (dev["section"], dev["index"]) == ("superseded_v9", 0) and "original plane-wave layer" in dev["recomputed"]["reason"]
    assert dev["recomputed"]["v3_value"] == dev["v3_value"] and dev["recomputed"]["v4_value"] == dev["v4_value"]
    assert cmp["not_record_derived_count"] == 28
    # without the absorber spectra the registry of version 10 misses exactly the section termination and counts 7 spectra
    cmp = audit.compare_registries(A["registry"], json.load(open(MANUSCRIPT_V10, encoding="utf-8")))
    nk = _not_ok(cmp)
    assert {s for _, s, _ in nk} == {"termination", "dataset"} and len([x for x in nk if x[1] == "termination"]) == 29
    assert {x for x in nk if x[1] == "dataset"} == {("DEVIATION", "dataset", "transmission spectra")}


def test_earlier_registries_deviate_only_where_version_9_corrects(run):
    from sprat.analysis import audit
    A, _ = run
    cmp = audit.compare_registries(A["registry"], json.load(open(os.path.join(EXPECTED, "numbers_registry_expected.json"), encoding="utf-8")))
    assert _not_ok(cmp) == V9_CORRECTED | V10_RENAMED_EARLY and cmp["counts"] == {"OK": 63, "DEVIATION": 12, "MISSING": 8}
    assert cmp["superseded_counts"] == {"CONSISTENT": 7, "DEVIATION": 2} and cmp["not_record_derived_count"] == 1
    cmp = audit.compare_registries(A["registry"], json.load(open(MANUSCRIPT, encoding="utf-8")))
    assert _not_ok(cmp) == V9_CORRECTED_V7 | V10_RENAMED_EARLY and cmp["counts"]["OK"] == 81
    cmp = audit.compare_registries(A["registry"], json.load(open(MANUSCRIPT_V8, encoding="utf-8")))
    assert _not_ok(cmp) == V9_CORRECTED_V8 | V10_RENAMED_EARLY and cmp["counts"]["OK"] == 83
    # the two superseded statements that version 9 rewrites (the pure-pixel value is no bound; the suppression factor)
    dev = [r for r in cmp["superseded"] if r["verdict"] != "CONSISTENT"]
    assert len(dev) == 2 and {r["index"] for r in dev} == {4, 7}


@needs_h14b
def test_audit_command_on_a_pack_layout(run_v10, tmp_path, monkeypatch, capsys):
    """The command of the supplementary, run from an unpacked pack: records/ and analysis/manuscript_registry.json."""
    from sprat.analysis import audit
    _, d = run_v10
    shutil.copytree(d / "records", tmp_path / "records")
    (tmp_path / "analysis").mkdir()
    shutil.copyfile(MANUSCRIPT_V10, tmp_path / "analysis" / "manuscript_registry.json")
    monkeypatch.chdir(tmp_path)
    audit.audit_cli("records", "tables", expected="analysis/manuscript_registry.json")
    out = capsys.readouterr().out
    assert "Verdicts: OK 223" in out and "DEVIATION" not in out and "MISSING" not in out
    assert "Superseded statements: CONSISTENT 21" in out and "Not record-derived: 29" in out


def test_outputs_carry_no_absolute_path(run):
    _, d = run
    for name in ("analysis.json", "predictions.json", "ANALYSIS_REPORT.md", "numbers_registry.json", "numbers_registry.md", "records.csv"):
        text = open(os.path.join(d, "tables", name), encoding="utf-8").read()
        assert str(d) not in text and "/tmp/" not in text and "/home/" not in text, name
    A = json.load(open(os.path.join(d, "tables", "analysis.json"), encoding="utf-8"))
    assert A["records_dir"] == "records"
    assert A["predictions"]["criteria_file"] == "sprat/data/predictions_criteria.json"


def test_manuscript_registry_structure():
    """Version 8 adds the section corrected_v8 to version 7; version 9 keeps every section of version 8, carries the
    rows it corrects with their version-8 values (value_v8), and adds corrected_v9 and superseded_v9; version 10 keeps
    every section of version 9, adds termination and corrected_v10, carries its changed rows with their version-9 value
    or derivation, and no longer says "campaign"."""
    v7 = json.load(open(MANUSCRIPT, encoding="utf-8"))
    v8 = json.load(open(MANUSCRIPT_V8, encoding="utf-8"))
    v9 = json.load(open(MANUSCRIPT_V9, encoding="utf-8"))
    v10 = json.load(open(MANUSCRIPT_V10, encoding="utf-8"))
    assert {k: v for k, v in v8.items() if not k.startswith("corrected_v8")} == v7
    assert set(v8) <= set(v9) and {"corrected_v9", "superseded_v9", "corrected_v9_note"} <= set(v9)
    assert len(v9["superseded_v9"]) == 12 and len(v9["corrected_v9"]) >= 80
    changed = [r for rows in v9.values() if isinstance(rows, list) for r in rows if isinstance(r, dict) and "value_v8" in r]
    assert len(changed) == 24
    assert set(v9) <= set(v10) and set(v10) - set(v9) == {"termination", "corrected_v10", "corrected_v10_note"}
    assert len(v10["termination"]) == 29 and len(v10["corrected_v10"]) == 2
    # the statements version 9 corrects are the same; one reason is reworded
    assert [(s["v3_value"], s["v4_value"]) for s in v10["superseded_v9"]] == [(s["v3_value"], s["v4_value"]) for s in v9["superseded_v9"]]
    assert sum(1 for a, b in zip(v9["superseded_v9"], v10["superseded_v9"]) if a["reason"] != b["reason"]) == 1
    rows10 = [r for rows in v10.values() if isinstance(rows, list) for r in rows if isinstance(r, dict)]
    assert sum(1 for r in rows10 if "value_v9" in r) == 2 and sum(1 for r in rows10 if "derivation_v9" in r) == 7
    assert "campaign" not in json.dumps(v10).lower()
    n9 = sum(1 for r in (x for rows in v9.values() if isinstance(rows, list) for x in rows) if isinstance(r, dict) and "quantity" in r)
    n10 = sum(1 for r in rows10 if "quantity" in r)
    assert n10 == n9 + 31


DERIVED_V10_2_ROWS = 14


@needs_h14b
def test_manuscript_registry_v10_2(run_v10):
    """Version 10.2 (SPRAT 1.2.1): the registry of version 10 with the section derived_v10_2, every record-derived row OK."""
    from sprat.analysis import audit
    A, _ = run_v10
    reg = json.load(open(MANUSCRIPT_V10_2, encoding="utf-8"))
    v10 = json.load(open(MANUSCRIPT_V10, encoding="utf-8"))
    assert set(reg) - set(v10) == {"derived_v10_2", "literature_v10_2", "references_added_v10_2", "derived_v10_2_note"}
    assert len(reg["derived_v10_2"]) == DERIVED_V10_2_ROWS
    assert [r["key"] for r in reg["references_added_v10_2"]] == ["fan2002", "delasson2018", "poblet2025"]
    for sec in v10:
        if sec not in ("created", "rule", "dataset"):
            assert reg[sec] == v10[sec], sec
    cmp = audit.compare_registries(A["registry"], reg)
    bad = [r for r in cmp["rows"] if r["verdict"] != "OK"]
    assert not bad, bad
    assert cmp["counts"] == {"OK": 223 + DERIVED_V10_2_ROWS}
    assert cmp["not_record_derived_count"] == 34 and len(cmp["superseded"]) == 21
