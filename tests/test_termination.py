"""The comparison of the two guide terminations (the absorber-terminated spectra, H14b).

Always run: the fit on a synthetic Fano line, the identity gate, the constants against the frozen criterion and its
amendment, and the index n_a = 1.300 under both terminations from the fixtures (the PML pair of
``fixtures/legacy_sample``, the absorber pair of ``fixtures/h14b_sample``, both byte for byte from the deposit), graded
as the job graded it. With the full deposits (SPRAT_LEGACY_DEPOSIT for raw_legacy and SPRAT_H14B_DEPOSIT for the
unpacked raw_h14b.tar.gz) every index is compared with the comparison the job wrote (``log/fano_S4b_528207.json``) and
with the amended grading (``grading/fano_S4b_528207_A1A2.json``)."""
import json
import math
import os

import numpy as np
import pytest

from sprat import legacy
from sprat.analysis import collect, termination

HERE = os.path.dirname(__file__)
ROOT = os.path.dirname(HERE)
SAMPLE = os.path.join(HERE, "fixtures", "legacy_sample")
H14B_SAMPLE = os.path.join(HERE, "fixtures", "h14b_sample")
FROZEN = os.path.join(ROOT, "predictions", "frozen")

# the row n_a = 1.300 of the comparison the job wrote on the cluster (log/fano_S4b_528207.json of raw_h14b)
JOB_1P300 = dict(pml=dict(q=-0.034182331767607106, lambda_r=1560.7495965528194, FWHM_nm=0.18339286923322495, Q=8510.415934263932,
                          t_end_over_tau=3.4807679048789644, T_max_near=0.9944552726196937, rms=0.018787411889817238),
                 absorber=dict(q=-0.021270133633516666, lambda_r=1560.7498646905601, FWHM_nm=0.30316316593975445,
                               Q=5148.217329940133, t_end_over_tau=4.595532585758078, T_max_near=0.9970114858870986,
                               rms=0.011501774036628803),
                 factor=1.6070576873924238)


def _legacy_deposit():
    cand = os.environ.get("SPRAT_LEGACY_DEPOSIT", "")
    return cand if cand and os.path.isfile(os.path.join(cand, "records_all.csv")) else None


def _h14b_deposit():
    cand = os.environ.get("SPRAT_H14B_DEPOSIT", "")
    return cand if cand and os.path.isfile(os.path.join(cand, "log", "fano_S4b_528207.json")) else None


def test_fit_recovers_a_fano_line():
    lam = np.linspace(1555.0, 1566.0, 601)
    truth = dict(T_bg=-0.01, T0=1.0, q=-0.12, lam_r=1560.75, dlam=0.30)
    fit = termination.fit(lam, termination.fano(lam, **truth))
    assert math.isclose(fit["q"], truth["q"], abs_tol=1e-7) and math.isclose(fit["lambda_r"], truth["lam_r"], abs_tol=1e-7)
    assert math.isclose(fit["FWHM_nm"], truth["dlam"], rel_tol=1e-6) and fit["rms"] < 1e-9


def test_identity_gate():
    a = dict(structure=dict(cell=dict(termination_x="pml", absorber_periods=8.0, guide_periods=25)),
             params=dict(run=dict(label="x", q_est=8327.0), numerics=dict(resolution=24)))
    b = json.loads(json.dumps(a))
    b["structure"]["cell"]["termination_x"] = "absorber"
    b["params"]["run"]["label"] = "x_em8b"
    assert termination.identity(a, b) == []
    b["params"]["numerics"]["resolution"] = 20
    assert termination.identity(a, b) == ["params.numerics.resolution: PML 24, absorber 20"]
    assert termination.identity(a, a) == ["termination of the second run is not the absorber"]


def test_constants_equal_the_frozen_criterion_and_its_amendment():
    import hashlib
    for name, digest in (("predictions_S4b_2026-09-24.json", "e06be4988c516ddbe583d988f7b18a972f8c43905cd795900b57d78ac42bfac4"),
                         ("predictions_S4b_amendment_A1A2_2026-09-24.json", "84d13c85df37dedf2cd93985fb775915e0849300139f2e5f2ba76125b54bb1c3")):
        assert hashlib.sha256(open(os.path.join(FROZEN, name), "rb").read()).hexdigest() == digest, name
    crit = json.load(open(os.path.join(FROZEN, "predictions_S4b_2026-09-24.json"), encoding="utf-8"))
    amend = json.load(open(os.path.join(FROZEN, "predictions_S4b_amendment_A1A2_2026-09-24.json"), encoding="utf-8"))
    assert [float(k) for k in crit["protocol"]["per_index"]] == list(termination.INDICES)
    assert crit["criterion"]["supports"] == "factor >= %g" % termination.SUPPORTS
    assert crit["criterion"]["refutes"].startswith("factor <= %g " % termination.REFUTES)
    assert "within 5 % of 4495.1" in crit["gates"]["S1"] and termination.S1_REFERENCE_Q == 4495.1 and termination.S1_TOLERANCE == 0.05
    assert [float(k) for k in amend["A1"]["graded_indices"]] == list(termination.A1_GRADED)
    assert [float(k) for k in amend["A1"]["reported_not_graded"]] == [1.300, 1.450]
    assert amend["fixed_at"].startswith("2026-09-24T15:46:39+03:00") and "15:46:39" in termination.AMENDMENT
    # the harmonic-inversion Q of the 25a cell that A2 divides by G is the q_est of each archived spectrum run
    assert [v["q_tahmin_cavity"] for v in crit["protocol"]["per_index"].values()] == [8327.0, 7196.0, 5778.0, 4323.0, 3057.0,
                                                                                      2138.0, 1662.0]


def test_index_1p300_under_both_terminations(tmp_path):
    out = str(tmp_path / "records")
    legacy.import_deposit(SAMPLE, out, log=lambda *a: None)
    legacy.import_deposit(H14B_SAMPLE, out, log=lambda *a: None)
    recs = collect.load(out)
    ab = [r for r in recs if r["structure"]["cell"]["termination_x"] == "absorber" and r["params"]["run"]["mode"] in ("spectrum", "reference")]
    assert len(ab) == 2 and all(r["task"]["legacy_wp"] == "H14b" and r["task"]["source"] == "legacy-verification" for r in ab)
    manifest = json.load(open(os.path.join(out, "_IMPORT_MANIFEST.json"), encoding="utf-8"))
    assert len(manifest["deposits"]) == 2 and manifest["counts"]["spectra"] == 2 and manifest["counts"]["references"] == 2
    T = termination.analyse(recs)
    row = T["rows"][0]
    assert row["n"] == 1.300 and row["status"] == "graded" and all(row["gates"].values())
    for side in ("pml", "absorber"):
        for k, v in JOB_1P300[side].items():
            assert math.isclose(row[side][k], v, rel_tol=1e-6, abs_tol=1e-9), (side, k, row[side][k], v)
    assert math.isclose(row["factor"], JOB_1P300["factor"], rel_tol=1e-5) and row["grade"] == "inconclusive"
    assert all(r["status"].startswith("missing") for r in T["rows"][1:])
    assert T["D22"]["verdict"] == "not a pass: inconclusive at n_a = 1.300" and not T["D22"]["complete"]
    assert T["A1"]["verdict"] == "no graded index" and T["S1"] is None and T["A2"] == dict(evaluated=False)
    # the spectrum analysis of section 3.5 keeps to the PML spectra; without an absorber spectrum there is no comparison
    from sprat.analysis import fano
    assert [s["label"] for s in fano.analyse_spectra(recs)] == [row["labels"]["pml"]]
    assert len(fano.analyse_spectra(recs, termination=None)) == 2
    assert termination.analyse([r for r in recs if r["structure"]["cell"]["termination_x"] == "pml"]) is None


@pytest.mark.skipif(_legacy_deposit() is None or _h14b_deposit() is None,
                    reason="full deposits not available (set SPRAT_LEGACY_DEPOSIT and SPRAT_H14B_DEPOSIT)")
def test_every_index_reproduces_the_job(tmp_path):
    out = str(tmp_path / "records")
    legacy.import_deposit(_legacy_deposit(), out, log=lambda *a: None)
    legacy.import_deposit(_h14b_deposit(), out, log=lambda *a: None)
    T = termination.analyse(collect.load(out))
    job = json.load(open(os.path.join(_h14b_deposit(), "log", "fano_S4b_528207.json"), encoding="utf-8"))
    graded = json.load(open(os.path.join(_h14b_deposit(), "grading", "fano_S4b_528207_A1A2.json"), encoding="utf-8"))
    assert len(job["rows"]) == len(T["rows"]) == 7
    for mine, theirs in zip(T["rows"], job["rows"]):
        assert mine["n"] == theirs["na"] and mine["status"] == theirs["status"] == "graded" and mine["grade"] == theirs["grade"]
        for side, key in (("pml", "archive"), ("absorber", "absorber")):
            for k in ("q", "lambda_r", "FWHM_nm", "Q", "t_end_over_tau", "T_max_near", "edge_margin_linewidths"):
                assert math.isclose(mine[side][k], theirs[key][k], rel_tol=1e-6, abs_tol=1e-9), (mine["n"], side, k)
            assert mine[side]["gates"] == theirs[key]["gates"]
        assert math.isclose(mine["factor"], theirs["factor"], rel_tol=1e-5)
    assert T["S1"] is True and graded["S1"] is True and abs(T["S1_value"] - 4530) < 1
    assert T["D22"]["verdict"] == "not a pass: inconclusive at n_a = 1.300" and graded["D22"]["verdict"] == "NOT A PASS: inconclusive at n_a = 1.3000"
    assert T["A1"]["verdict"] == "pass" and graded["A1"]["verdict"] == "PASS" and T["A1"]["complete"]
    for n, g in graded["A2"]["rows"].items():
        mine = T["A2"]["rows"][float(n)]
        assert math.isclose(mine["G"], g["G"], rel_tol=1e-5) and math.isclose(mine["Q_w"], g["Q_w"], rel_tol=1e-5)
    assert math.isclose(T["A2"]["G_1p33"], graded["A2"]["G_1.33_interpolated"], rel_tol=1e-5)
    assert math.isclose(T["A2"]["Q_w_fall"], graded["A2"]["Q_w_fall_1.300_to_1.450"], rel_tol=1e-5)
    assert math.isclose(T["A2"]["Q_cell_fall_harminv"], graded["A2"]["cell_fall_1.300_to_1.450"], rel_tol=1e-12)
    assert abs(T["A2"]["G_1p33_harminv"] - 1.5411) < 5e-5 and abs(T["A2"]["fwhm_cell_opening"] - 5.31) < 5e-3
    s = T["summary"]
    assert s["abs_q_absorber_max"] < 0.024 and 8.9 < s["factor_inner_min"] < 9.0 and 213 < s["factor_inner_max"] < 214
    assert s["T_max_near_absorber_max"] < 0.998 and s["T_max_near_pml_max"] > 1.24
