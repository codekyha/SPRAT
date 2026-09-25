"""The English criteria file must match the frozen originals value by value, and the frozen files
must carry their recorded checksums."""
import hashlib
import json
import os

from sprat.analysis.predictions import load_criteria

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def test_frozen_files_checksums():
    c = load_criteria()
    for name, s in c["sets"].items():
        path = os.path.join(ROOT, s["frozen_file"])
        assert sha(path) == s["sha256"], name
        if "sealed_copy" in s:
            assert sha(os.path.join(ROOT, s["sealed_copy"])) == s["sealed_sha256"]


def test_first_set_values_match_frozen():
    c = load_criteria()
    frozen = json.load(open(os.path.join(ROOT, c["sets"]["first"]["frozen_file"]), encoding="utf-8"))
    by_id = {o["id"]: o for o in frozen["ongoruler"]}
    for crit in c["criteria"]:
        if crit["set"] != "first":
            continue
        o = by_id[crit["id"]]
        expected_frozen = o["beklenen"]
        if isinstance(expected_frozen, dict) and "tolerans_bagil" in expected_frozen:
            # P8: the tolerance sits inside the expected block in both files
            assert crit["expected"]["eta_a"] == expected_frozen["eta_a"]
            assert crit["expected"]["S_th_nm_per_RIU"] == expected_frozen["S_th_nm_per_RIU"]
            assert crit["expected"]["tolerance_rel"] == expected_frozen["tolerans_bagil"]
        else:
            assert crit["expected"] == expected_frozen, crit["id"]
        if "tolerans_bagil" in o:
            assert crit["tolerance_rel"] == o["tolerans_bagil"], crit["id"]
        if "tolerans" in o:
            assert crit["tolerance_abs"] == o["tolerans"], crit["id"]
        if "tolerans_ln_bagil" in o:
            assert crit["tolerance_ln_rel"] == o["tolerans_ln_bagil"]
    A = c["sets"]["first"]["archive_constants"]
    F = frozen["arsiv_sabitleri"]
    assert (A["Q3_cl12"], A["Q4_cl12"], A["Q5_cl12_margin1.27"], A["kappa_paper"], A["kappa_error"]) == \
           (F["Q3_nk12"], F["Q4_nk12"], F["Q5_nk12_pay1.27"], F["kappa_makale"], F["kappa_hata"])
    assert (A["Q_peak_rd0.065"], A["Q_trough_rd0.085"], A["f_ref"], A["f_cen"]) == (F["Q_pk_rd0.065"], F["Q_tr_rd0.085"], F["f_ref"], F["f_cen"])


def test_second_set_values_match_frozen():
    c = load_criteria()
    frozen = json.load(open(os.path.join(ROOT, c["sets"]["second"]["frozen_file"]), encoding="utf-8"))
    by_id = {o["id"]: o for o in frozen["ongoruler"]}
    for crit in c["criteria"]:
        if crit["set"] != "second":
            continue
        text = by_id[crit["id"]]["beklenen"]
        e = crit["expected"]
        if crit["id"] == "PA1":
            assert text.startswith("[1.64, 1.78]") and e == [1.64, 1.78]
        elif crit["id"] == "PA2":
            assert text == "< 0.02" and e == {"max": 0.02}
        elif crit["id"] == "PA3":
            assert text == "1.723 +- 0.03" and e == 1.723 and crit["tolerance_abs"] == 0.03
        elif crit["id"] == "PA4":
            assert text == "5041.1 +- 1 %" and e == 5041.1 and crit["tolerance_rel"] == 0.01
        elif crit["id"] == "PA5":
            assert text == "2925.9 +- 1 %" and e == 2925.9 and crit["tolerance_rel"] == 0.01
        elif crit["id"] == "PA6":
            assert text == "kappa(0.100a) < kappa(0.060a)" and e == {"less_than": 1.9595}
        else:
            assert crit.get("closed") is True


def test_grader_on_sample_reports_missing(tmp_path):
    """On the small fixture most criteria are MISSING; the grader must not raise."""
    from sprat import legacy
    from sprat.analysis import collect, predictions
    SAMPLE = os.path.join(os.path.dirname(__file__), "fixtures", "legacy_sample")
    legacy.import_deposit(SAMPLE, str(tmp_path / "r"), log=lambda *a: None)
    recs = collect.load(str(tmp_path / "r"))
    g = predictions.grade(recs, None)
    ids = {r["id"] for r in g["rows"]}
    assert {"P1", "P4b", "P10", "PA1", "PB1"} <= ids
    v = {r["id"]: r["verdict"] for r in g["rows"]}
    assert v["P4b"] == "PASS"              # the H1a record is in the sample
    assert v["P1"] == "MISSING" and v["PB1"] == "CLOSED"
    assert "ERROR" not in v.values()
