"""SPRAT 1.2.1 on the full deposits (SPRAT_LEGACY_DEPOSIT for raw_legacy and SPRAT_H14B_DEPOSIT for the unpacked
raw_h14b.tar.gz): the reflectionless quality factor and figure of merit at the rows of table 5, the convergence of the
resonance frequency with resolution, the guided-mode wavevector and the per-row factor at k_x = 0, and the
absorption-limited quality factor with the mixed points of the grid counted either way; every value at the precision
the manuscript prints it, and the registry rows in their sections."""
import os

import pytest

from sprat import legacy
from sprat.analysis import pipeline


def _deposits():
    a, b = os.environ.get("SPRAT_LEGACY_DEPOSIT", ""), os.environ.get("SPRAT_H14B_DEPOSIT", "")
    ok = a and os.path.isfile(os.path.join(a, "records_all.csv")) and b and os.path.isdir(os.path.join(b, "v4_spectra"))
    return (a, b) if ok else None


@pytest.fixture(scope="module")
def analysis(tmp_path_factory):
    dep = _deposits()
    if not dep:
        pytest.skip("full deposits not available (set SPRAT_LEGACY_DEPOSIT and SPRAT_H14B_DEPOSIT)")
    rec = str(tmp_path_factory.mktemp("records"))
    quiet = lambda *a: None                                               # noqa: E731
    legacy.import_deposit(dep[0], rec, log=quiet)
    legacy.import_deposit(dep[1], rec, log=quiet)
    root = os.path.dirname(os.path.dirname(__file__))
    import shutil
    shutil.copy(os.path.join(root, "tests", "fixtures", "pwe_na1.3300_r0.2.json"), os.path.join(rec, "pwe_na1.3300_r0.2.json"))
    return pipeline.analyze(rec, str(tmp_path_factory.mktemp("tables")), log=quiet)


def test_table5_reflectionless(analysis):
    t = analysis["combined"]["physical_values"]["table"]
    assert [r["n"] for r in t] == [1.305, 1.330, 1.350, 1.375, 1.400, 1.425, 1.445]
    assert [round(r["Q_w"]) for r in t] == [4938, 4502, 4114, 3714, 3343, 3009, 2792]
    assert [round(r["FOM_w"]) for r in t] == [2079, 1842, 1644, 1440, 1255, 1090, 979]
    assert [r["source"][:5] for r in t] == ["inter", "harmo", "spect", "spect", "spect", "spect", "inter"]
    assert round(analysis["combined"]["physical_values"]["at_1p33"]["deviation_pct"], 1) == -1.0


def test_frequency_convergence(analysis):
    cv = analysis["convergence"]
    assert round(cv["resolution_f_max_change_pct"], 3) == 0.023
    rc = analysis["barrier"]["resolution_check"]
    assert round(rc["f_change_pct"], 3) == 0.029 and round(abs(rc["dlambda_nm_at_target"]), 2) == 0.46
    assert round(rc["defect_radius_pixels"], 2) == 1.44


def test_channel_counterfactual(analysis):
    ch = analysis["barrier"]["channels"]
    assert round(ch["beta_over_pi"], 4) == 0.5465 and round(ch["per_row_factor_kx0"], 2) == 3.93
    assert round(ch["per_row_factor_zone_edge"], 2) == 7.11


def test_absorption_brackets(analysis):
    ab = analysis["barrier"]["performance"]["absorption"]
    assert float("%.2g" % ab["Q_abs_mixed_as_analyte"]) == 6700 and float("%.2g" % ab["Q_abs_mixed_as_silicon"]) == 9600


def test_registry_rows(analysis):
    reg = analysis["registry"]
    where = {r["quantity"]: (sec, r["value"]) for sec, rows in reg.items() if isinstance(rows, list)
             for r in rows if isinstance(r, dict) and "quantity" in r}
    assert where["per-row factor of the slowest channel at k_x = 0, reference resonance"] == ("derived", 3.93)
    assert where["per-row factor of the slowest channel at k_x = pi/a, reference resonance"] == ("derived", 12.2)
    assert where["guided-mode wavevector beta at the reference resonance"] == ("derived", 0.5465)
    assert where["largest change of f_r from resolution 20 to 32"] == ("systematic_runs", 0.023)
    assert where["change of f_r from resolution 24 to 32 at the reference radius, N_sep = 5"] == ("derived", 0.029)
    assert where["defect radius of the reference geometry in pixels at resolution 24"] == ("derived", 1.44)
    assert where["Q_abs with the mixed points counted wholly as analyte and wholly as silicon"] == ("absorption", [6700.0, 9600.0])
    assert where["spectral Q_w interpolated to n_a = 1.33"] == ("reflectionless", 4455)
    assert [r["quantity"][:30] for r in reg["reflectionless"]] == [
        "reflectionless Q_w at the rows", "reflectionless FOM_w at the ro", "sources of Q_w at the rows of ",
        "spectral Q_w interpolated to n", "deviation of the spectral Q_w "]
