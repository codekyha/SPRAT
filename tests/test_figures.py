"""Figures without the deposit: argument handling, number format, house style, figure 1 from one
synthetic record (byte-identical on a second build) and the skip path when the records lack the data.
The build of all seven figures from the deposit is test_analysis.py::test_figures_build."""
import json
import os

import pytest

pytest.importorskip("matplotlib")

from sprat import records  # noqa: E402
from sprat.figures import build, style  # noqa: E402
from sprat.params import default_params  # noqa: E402
from sprat.structure import default_structure  # noqa: E402


def test_parse_which():
    assert build.parse_which(None) == list(range(1, 8))
    assert build.parse_which("all") == list(range(1, 8))
    assert build.parse_which("7, 1,3") == [1, 3, 7]
    assert build.parse_which("figure4") == [4]
    assert build.parse_which([2, "5", 2]) == [2, 5]
    for bad in ("8", "x", ""):
        with pytest.raises(ValueError):
            build.parse_which(bad)


def test_number_format():
    assert style.num(6927.2, 0) == "6927"
    assert style.num(12345.0, 0) == "12 345"
    assert style.num(-0.1, 2) == "\u22120.10"
    assert style.num(-0.0001, 2) == "0.00"
    assert style.pct(-30.96) == "\u221231%"


def test_style():
    assert style.font_family() in ("Arial", "Liberation Sans", "DejaVu Sans")
    for scheme in ("campaign", "mechanism"):
        rc = style.rc(scheme)
        assert rc["pdf.fonttype"] == 42 and not rc["axes.spines.top"] and rc["xtick.direction"] == "out"
    with pytest.raises(ValueError):
        style.rc("other")


def _one_record(records_dir):
    s, p = default_structure(), default_params()
    assert (s["defect"]["radius"], s["defect"]["row"], s["cell"]["cladding_rows"]) == (0.060, 4, 12)
    rec = records.new_record(s, p, dict(label="reference", source="sprat"))
    rec["result"] = dict(modes=[dict(f=0.30461, Q=6927.0, amplitude=1.0, error=0.0, wavelength_nm=1580.36, fwhm_nm=0.228)],
                         t_used=1000.0, Q_limit=17000.0)
    records.save(rec, os.path.join(records_dir, "reference.json"))


def test_figure1_from_one_record(tmp_path):
    _one_record(str(tmp_path / "records"))
    done = build.build_all(str(tmp_path / "records"), str(tmp_path / "figures"), which="1", log=lambda *a: None)
    assert list(done) == ["figure1"]
    pdf, png = done["figure1"]
    assert os.path.getsize(pdf) > 1000 and os.path.getsize(png) > 1000
    data = open(pdf, "rb").read()
    assert b"Figure 1" in data and b"CreationDate" not in data
    again = build.build_all(str(tmp_path / "records"), str(tmp_path / "again"), which="1", log=lambda *a: None)
    assert open(again["figure1"][0], "rb").read() == data          # reproducible byte for byte
    assert not (tmp_path / "tables").exists()                       # figure 1 alone needs no analysis


def test_missing_data_is_skipped(tmp_path):
    _one_record(str(tmp_path / "records"))
    tables = tmp_path / "tables"
    tables.mkdir()
    (tables / "analysis.json").write_text(json.dumps({"errors": {}}))
    msgs = []
    done = build.build_all(str(tmp_path / "records"), str(tmp_path / "figures"), which="2,5", tables=str(tables),
                           log=msgs.append)
    assert done == {}
    assert any("figure 2" in m and "skipped" in m and "bands" in m for m in msgs)
    assert any("figure 5" in m and "skipped" in m for m in msgs)
