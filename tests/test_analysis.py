"""Analysis regression: from the imported deposit, every row of the deposited numbers registry is
reproduced (needs the full deposit; see test_legacy_import for how it is found).  The plane-wave
layer and the Fano/coupling functions are tested without the deposit."""
import math
import os

import numpy as np
import pytest

from sprat import pwe
from sprat.analysis import fano

HERE = os.path.dirname(__file__)
ROOT = os.path.dirname(HERE)


def _deposit():
    for cand in (os.environ.get("SPRAT_LEGACY_DEPOSIT", ""), os.path.join(HERE, "data", "zenodo_pack_v2")):
        if cand and os.path.isfile(os.path.join(cand, "records_all.csv")):
            return cand
    return None


def test_fano_fit_recovers_parameters():
    lam = np.linspace(1578, 1583, 601)
    truth = dict(base=0.02, T0=1.0, q=0.4, lam_r=1580.36, fwhm=0.23)
    T = fano.fano(lam, **truth)
    fit = fano.fit_fano(lam, T, 1580.3, 0.25)
    assert math.isclose(fit["lambda_r"], truth["lam_r"], abs_tol=1e-6) and math.isclose(fit["fwhm"], truth["fwhm"], rel_tol=1e-6)
    assert math.isclose(fit["q"], truth["q"], abs_tol=1e-6) and math.isclose(fit["Q"], 1580.36 / 0.23, rel_tol=1e-6)


def test_coupling_branches_identity():
    b = fano.coupling_branches(5000.0, 0.04)
    for br in b.values():
        assert math.isclose(1 / br["Q_coupling"] + 1 / br["Q_intrinsic"], 1 / 5000.0, rel_tol=1e-12)
        r = br["ratio_Qi_over_Qc"]
        assert math.isclose(((1 - r) / (1 + r)) ** 2, 0.04, rel_tol=1e-9)


def test_printed_precision_of_registry_numbers():
    from sprat.analysis import audit
    assert audit.number("9.1e3") == (9100.0, 100.0) and audit.number("0.003 %") == (0.003, 0.001)
    assert audit.number(97500.0) == (97500.0, 100.0) and audit.number(9000) == (9000.0, 1000.0)
    assert audit.number(1.9595)[1] == pytest.approx(1e-4) and audit.number(-8.901e-05)[1] == pytest.approx(1e-8)
    ok, _, last = audit.compare(2.1466, 2.146)
    assert ok and last
    assert not audit.compare(2.148, 2.146)[0] and not audit.compare(857, 858, exact=True)[0]
    assert audit.compare(33.63, "34 %")[0] and audit.compare(9131.4, "9.1e3")[0] and not audit.compare(9300.0, "9.1e3")[0]
    assert audit.compare_text("kappa = 1.89009 +- 0.0192", "kappa = 1.8901 +- 0.0192")[0]
    assert not audit.compare_text("kappa = 1.8901 at r_d", "kappa = 1.8901 +- 0.0192")[0]


def test_audit_lists_non_record_rows_apart():
    from sprat.analysis import audit
    new = dict(derived=[dict(quantity="x", value=1.23)])
    exp = dict(spec="s", derived=[dict(quantity="x", value=1.2)], literature_v7=[dict(key="a", values={"S": "63 +- 9"})],
               references_added_v7=[dict(key="b", doi="10.1/x")], references_removed_v7=["c"], dataset_note="text",
               superseded=[dict(v3_value="k = 1.0", v4_value="k = 2.0", reason="r")])
    cmp = audit.compare_registries(new, exp)
    assert cmp["counts"] == {"OK": 1} and cmp["not_record_derived_count"] == 4
    assert cmp["superseded_counts"] == {"MISSING": 1} and not cmp["ok"]


def test_vertical_channel_algebra():
    from sprat.analysis import barrier
    v = barrier.vertical_channel(4501.77, 9131.42)
    q0 = v["Q_tot_no_vertical"]
    for r in ("3", "10", "50"):
        qp = v["Q_perp_for_reduction"][r]
        assert math.isclose(1 / (1 / q0 + 1 / qp), (1 - float(r) / 100) * q0, rel_tol=1e-12)
    assert math.isclose(v["reduction_pct_at_Q_perp"]["1e4"], 100 * (1 - v["Q_tot_at_Q_perp"]["1e4"] / q0), rel_tol=1e-12)
    assert v["overtakes_absorption_below"] == 9131.42 and v["largest_channel_below"] == 4501.77


@pytest.mark.slow
def test_pwe_quick_reproduces_channels():
    """The plane-wave layer in its quick setting reproduces the physical channels of the full basis (M = 10) to 1e-3:
    the slowest physical root is the zone-edge root; the slower root of the square basis is reported as a truncation
    artefact (the 'least-evanescent pair' of SPRAT 1.0.0 and of the 2026 plane-wave file)."""
    r = pwe.feasibility(quick=True, log=lambda *a: None)
    k = r["kappa"]["reference"]
    assert abs(k["kappa_pred"] - 1.9620) < 2e-3 and k["kappa_pred"] == k["kappa_pred_zone_edge"]
    assert abs(k["kappa_boundary_artefact"] - 1.9015) < 2e-3 and k["kappa_next"] > 6.0
    assert abs(k["kappa_at_kx0"] - 1.3676) < 2e-3 and abs(k["kappa_at_kxpi"] - 2.5039) < 3e-3
    assert abs(k["n_g"] - 2.815) < 0.01
    s = r["kappa"]["second"]
    assert abs(s["kappa_pred"] - 1.7398) < 2e-3 and s["kappa_pred"] == s["kappa_pred_zone_edge"]
    assert abs(s["kappa_boundary_artefact"] - 1.6783) < 2e-3
    assert abs(r["fabry_perot"]["D_over_a"] - 11.984) < 0.02
    assert abs(r["bulk_gap"]["f_lo"] - 0.27506) < 5e-4 and abs(r["bulk_gap"]["f_hi"] - 0.34596) < 5e-4
    assert pwe.kappa_block(r, "reference") is k and pwe.channels_at(r, 0.30461)[1] > 1.9


@pytest.mark.skipif(_deposit() is None, reason="full deposit not available (set SPRAT_LEGACY_DEPOSIT)")
def test_registry_reproduces_deposit(tmp_path):
    import json
    from sprat import legacy
    from sprat.analysis import audit, pipeline
    import shutil
    legacy.import_deposit(_deposit(), str(tmp_path / "records"), log=lambda *a: None)
    shutil.copyfile(os.path.join(HERE, "fixtures", "pwe_na1.3300_r0.2.json"), tmp_path / "records" / "pwe_na1.3300_r0.2.json")
    A = pipeline.analyze(str(tmp_path / "records"), str(tmp_path / "tables"), log=lambda *a: None)
    # the registry of manuscript version 9 on the records it was built from: every row but the two that SPRAT 1.2.0
    # renames (records of the systematic runs, threshold-mask estimate); version 10 is audited in test_audit_manuscript
    expected = json.load(open(os.path.join(ROOT, "campaigns", "manuscript", "expected", "numbers_registry_manuscript_v9.json")))
    cmp = audit.compare_registries(A["registry"], expected)
    bad = {(r["verdict"], r["quantity"]) for r in cmp["rows"] if r["verdict"] != "OK"}
    assert bad == {("MISSING", "records from the parameter campaign"), ("MISSING", "eta_a, threshold mask of the campaign")}, bad
    assert cmp["counts"]["OK"] == 190
    # the unpublished 2.0.2 registry: 63 rows as before, the others corrected by version 9 or renamed (test_audit_manuscript)
    old = json.load(open(os.path.join(ROOT, "campaigns", "manuscript", "expected", "numbers_registry_expected.json")))
    assert audit.compare_registries(A["registry"], old)["counts"]["OK"] == 63
    # headline numbers
    assert abs(A["analyte_sweep"]["at_1p33"]["S"] - 646.52) < 0.01
    assert abs(A["kappa_series"]["kappa"] - 1.8901) < 5e-4
    assert abs(A["barrier"]["series"]["kappa"] - 1.9595) < 5e-5
    assert abs(A["fields"]["reference"]["eta_mask"] - 0.6181) < 5e-4
    assert abs(A["fields"]["reference"]["eta_hf"] - 0.5523) < 5e-4 and A["fields"]["reference"]["reconstruction_max_abs"] < 1e-4
    assert A["counts"]["dead_zone_admitted"] == 436 and A["counts"]["valid"] == 741
    # the modes of the finite waveguide in the preliminary sweep, flagged per record, and the admitted cavity records
    assert A["counts"]["waveguide_mode_admitted"] == 537
    ps = A["preliminary_sweep"]
    assert ps["waveguide_modes"]["admitted"] == 537 and ps["cavity_modes"]["admitted"] == 46
    assert 50 < ps["waveguide_modes"]["Q_min"] and ps["waveguide_modes"]["Q_max"] < 80 < ps["cavity_modes"]["Q_min"]
    assert ps["groups"] == 48 and ps["cavity_groups"]["groups"] == 15 and abs(ps["cavity_groups"]["median_ratio"] - 0.984) < 5e-4
    import csv
    rows = list(csv.DictReader(open(os.path.join(str(tmp_path / "tables"), "records.csv"))))
    assert sum(1 for r in rows if r["waveguide_mode"] == "True") == 537
    assert sum(1 for r in rows if r["valid"] == "True" and r["waveguide_mode"] == "False") == 204
    v = {r["id"]: r["verdict"] for r in A["predictions"]["rows"] if r["id"] in ("P4b", "P12", "P13", "PA1", "PA2", "PA3")}
    assert v == {"P4b": "PASS", "P12": "PASS", "P13": "PASS", "PA1": "PASS", "PA2": "FAIL", "PA3": "PASS"}


@pytest.mark.skipif(_deposit() is None, reason="full deposit not available (set SPRAT_LEGACY_DEPOSIT)")
def test_figures_build(tmp_path):
    pytest.importorskip("matplotlib")
    from sprat import legacy
    from sprat.figures import build
    legacy.import_deposit(_deposit(), str(tmp_path / "records"), log=lambda *a: None)
    done = build.build_all(str(tmp_path / "records"), str(tmp_path / "figures"), log=lambda *a: None)
    assert set(done) == {f"figure{i}" for i in range(1, 8)}
    for files in done.values():
        assert all(os.path.getsize(f) > 1000 for f in files)
