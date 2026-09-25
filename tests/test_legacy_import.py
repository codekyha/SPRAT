"""Legacy import: the converted records must flatten to the numbers of the deposited records_all.csv.

The small fixture under ``fixtures/legacy_sample`` always runs.  When the full deposit is
available (environment variable SPRAT_LEGACY_DEPOSIT pointing at the unpacked zenodo_pack
directory, or ``tests/data/zenodo_pack_v2``), all 858 rows are checked.
"""
import csv
import json
import math
import os

import pytest

from sprat import legacy, records
from sprat.analysis import collect

HERE = os.path.dirname(__file__)
SAMPLE = os.path.join(HERE, "fixtures", "legacy_sample")


def _deposit():
    for cand in (os.environ.get("SPRAT_LEGACY_DEPOSIT", ""), os.path.join(HERE, "data", "zenodo_pack_v2")):
        if cand and os.path.isfile(os.path.join(cand, "records_all.csv")) and \
                (os.path.exists(os.path.join(cand, "harminv_json.tar.gz")) or os.path.isdir(os.path.join(cand, "harminv_json"))):
            return cand
    return None


def _compare(deposit, out):
    legacy.import_deposit(deposit, out, log=lambda *a: None)
    recs = collect.load(out)
    rows = {r["legacy_file"]: r for r in collect.table(recs)}
    expected = list(csv.DictReader(open(os.path.join(deposit, "records_all.csv"))))
    assert expected
    for e in expected:
        r = rows[e["dosya"]]
        assert math.isclose(float(e["Q"]), r["Q"], rel_tol=1e-12)
        assert math.isclose(float(e["f_r"]), r["f_r"], rel_tol=1e-12)
        assert math.isclose(float(e["pay"]), r["margin"], rel_tol=1e-9)
        assert math.isclose(float(e["gap_konumu"]), r["gap_position"], rel_tol=1e-9)
        assert (e["gecerli"] == "True") == r["valid"]
        assert (e["Q_le_Qlim"] == "True") == r["resolved"]
        assert int(e["nsep"]) == r["row"] and int(e["nkaplama"]) == r["cladding_rows"] and int(e["res"]) == r["resolution"]
        assert math.isclose(float(e["rd"]), r["radius"]) and math.isclose(float(e["na"]), r["n"])
        assert ({"pml": "pml", "emici": "absorber"}[e["sonlandirma"]]) == r["termination_x"]
        assert (e["bariyer_satir_r"].replace(":", ":").replace(",", ";") == r["rows"]) or (not e["bariyer_satir_r"] and not r["rows"])
        assert e["kod_sha256"] == r["code_sha256"]
        assert math.isclose(float(e["duvar_s"]), r["wall_s"], rel_tol=1e-12)
        assert e["dislama_nedeni"] in ("", "Q > Q_lim", "gap konumu 0.2-0.8 disi", "mod yok")
        mapped = {"": "", "Q > Q_lim": "Q > Q_lim", "gap konumu 0.2-0.8 disi": "gap position outside 0.2-0.8", "mod yok": "no mode"}
        assert mapped[e["dislama_nedeni"]] == r["exclusion_reason"]
    return recs, rows


def test_sample_round_trip(tmp_path):
    recs, rows = _compare(SAMPLE, str(tmp_path / "records"))
    assert sum(1 for r in rows.values() if r["mode"] == "harminv") == 6
    v4 = [r for r in recs if r["task"]["source"] == "legacy-verification"]
    assert len(v4) == 3 and all(r["provenance"]["legacy_sha256"] for r in v4)
    h5 = [r for r in v4 if r["task"]["legacy_wp"] == "H5"][0]
    assert h5["structure"]["rows"] == {"2": 0.18}
    abs_ = [r for r in v4 if r["task"]["legacy_wp"] == "H4b"][0]
    assert abs_["structure"]["cell"]["termination_x"] == "absorber" and abs_["structure"]["cell"]["absorber_periods"] == 8.0
    ref = [r for r in recs if r["params"]["run"]["mode"] == "reference"]
    spec = [r for r in recs if r["params"]["run"]["mode"] == "spectrum"]
    assert len(ref) == 1 and len(spec) == 1 and len(spec[0]["result"]["f"]) == len(spec[0]["result"]["flux_out"])
    bands = {r["result"]["task"]: r for r in recs if r["params"]["run"]["mode"] == "bands"}
    assert set(bands) == {"bulk", "sweep", "w1", "calibration"}
    assert bands["bulk"]["result"]["gap"]["exists"] and 0.27 < bands["bulk"]["result"]["gap"]["lower"] < 0.28
    assert bands["calibration"]["result"]["a_nm"] > 400
    assert len(bands["sweep"]["result"]["sweep"]) == 11
    # the record schema round-trips through save/load
    path = tmp_path / "records" / (recs[0]["task"]["label"] + ".json")
    assert records.load(str(path))["schema"] == "sprat-record-1.0"
    manifest = json.load(open(tmp_path / "records" / "_IMPORT_MANIFEST.json"))
    assert manifest["counts"]["harminv"] == 3 and manifest["counts"]["verification"] == 3


def test_gap_table_from_bands_matches_builtin(tmp_path):
    legacy.import_deposit(SAMPLE, str(tmp_path / "r"), log=lambda *a: None)
    recs = collect.load(str(tmp_path / "r"))
    table = records.gap_table_from_records(recs)
    assert table is not None and len(table) == 11
    for (n, lo, hi), (n2, lo2, hi2) in zip(table, records.GAP_EDGES_MPB):
        assert math.isclose(n, n2) and abs(lo - lo2) < 1e-5 and abs(hi - hi2) < 1e-5


@pytest.mark.skipif(_deposit() is None, reason="full deposit not available (set SPRAT_LEGACY_DEPOSIT)")
def test_full_deposit_round_trip(tmp_path):
    recs, rows = _compare(_deposit(), str(tmp_path / "records"))
    flat = collect.table(recs)
    c = records.counts(flat)
    assert c["harminv"] == 858 and c["valid"] == 741 and c["resolved"] == 841
    assert c["excluded_margin"] == 17 and c["excluded_gap"] == 100
    assert c["spectra"] == 7 and c["references"] == 7 and c["fields"] == 2 and c["bands"] == 4 and c["pwe"] == 1
