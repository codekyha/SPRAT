import math

import pytest

from sprat._textfile import TextFileError, parse_list, parse_text
from sprat.params import (expand_sweeps, harminv_time, parse_override, parse_params_text, q_limit,
                          resolve, source_time)
from sprat.structure import parse_structure_text, structure_signature, structure_text

REFERENCE = """
[lattice]
type = square
a_nm = 481.4
rod_radius = 0.20
rod_eps = 11.9025

[cell]
guide_periods = 25
cladding_rows = 12
pad_x = 1.5
pad_y = 1.0
termination_x = pml

[waveguide]
type = W1

[defect]
type = rod
radius = 0.060
row = 4
dx = 0
dy = 0
"""


def test_parse_text_sections_keys_comments():
    s = parse_text("[a]\nx = 1 # c\n# only comment\n[B]\ny=2\nbare line\n", "t")
    assert s == {"a": [(2, "x", "1")], "b": [(5, "y", "2"), (6, None, "bare line")]}


def test_value_before_section_is_an_error():
    with pytest.raises(TextFileError):
        parse_text("x = 1\n[a]\n", "t")


def test_structure_defaults_and_validation():
    s = parse_structure_text(REFERENCE)
    assert s["lattice"]["a_nm"] == 481.4 and s["cell"]["absorber_periods"] == 8.0
    assert s["waveguide"]["type"] == "w1" and s["rows"] == {} and s["rods"] == []
    with pytest.raises(TextFileError):
        parse_structure_text(REFERENCE + "\n[cell]\nguide_periods = 3\n")           # duplicate section key
    with pytest.raises(TextFileError):
        parse_structure_text(REFERENCE.replace("row = 4", "row = 13"))              # beyond the cladding
    with pytest.raises(TextFileError):
        parse_structure_text(REFERENCE + "\n[cell]\nfoo = 1\n")                      # unknown key
    with pytest.raises(TextFileError):
        parse_structure_text(REFERENCE + "\n[rows]\n0 = 0.1\n")                      # guide row


def test_structure_rows_and_rods():
    s = parse_structure_text(REFERENCE + "\n[rows]\n2 = 0.18\n-3 = 0.21\n[rods]\nadd 3 -2 0.2\nremove 0 4\n")
    assert s["rows"] == {2: 0.18, -3: 0.21}
    assert s["rods"] == [{"op": "add", "x": 3.0, "y": -2.0, "radius": 0.2}, {"op": "remove", "x": 0.0, "y": 4.0}]
    text = structure_text(s)
    assert parse_structure_text(text) == s


def test_signature_encodes_the_knobs():
    s = parse_structure_text(REFERENCE)
    assert structure_signature(s) == "rd0.0600_dx+0.000_dy+0.000_row4_cl12_nx25_pml"
    s["cell"]["termination_x"] = "absorber"
    s["rows"] = {2: 0.18}
    assert structure_signature(s) == "rd0.0600_dx+0.000_dy+0.000_row4_cl12_nx25_abs8_rows2-0.18"


def test_parse_list_range_and_list():
    assert parse_list("0.050 : 0.070 : 0.005", "float", "t") == [0.05, 0.055, 0.06, 0.065, 0.07]
    assert parse_list("8, 10, 12", "int", "t") == [8, 10, 12]
    with pytest.raises(TextFileError):
        parse_list("0.05:0.072:0.005", "float", "t")


def test_params_defaults_resolution_and_label():
    s = parse_structure_text(REFERENCE)
    p = parse_params_text("[run]\nmode = harminv\nq_est = 6927\n[harminv]\nmargin = 4.5\n")
    r = resolve(s, p)
    assert math.isclose(r["source"]["fcen"], 481.4 / 1550.0)
    assert r["harminv"]["t"] == harminv_time(4.5, 6927, 481.4 / 1550.0) == 32000.0
    assert math.isclose(q_limit(r["source"]["fcen"], r["harminv"]["t"]), math.pi * 481.4 / 1550.0 * 32000.0)
    assert r["run"]["label"] == "harminv_rd0.0600_dx+0.000_dy+0.000_row4_cl12_nx25_pml_na1.3300_res24_pml1_t32000"
    assert source_time(0.06, 5.0) == pytest.approx(166.6667, rel=1e-6)


def test_reference_label_drops_the_defect():
    s = parse_structure_text(REFERENCE)
    p = resolve(s, parse_params_text("[run]\nmode = reference\n"))
    assert p["run"]["label"] == "reference_cl12_nx25_pml_na1.3300_res24_pml1_nf201"


def test_sweep_expansion_product_and_zip():
    s = parse_structure_text(REFERENCE)
    p = parse_params_text("[run]\nmode = harminv\n[sweep]\nstructure.defect.radius = 0.05, 0.06\n"
                          "structure.cell.cladding_rows = 8, 12\nanalyte.n = 1.30, 1.45\n")
    pts = expand_sweeps(s, p)
    assert len(pts) == 8
    assert {pt[0]["defect"]["radius"] for pt in pts} == {0.05, 0.06}
    p2 = parse_params_text("[run]\nmode = harminv\n[sweep]\nstructure.defect.radius = 0.05, 0.06, 0.07\n"
                           "run.q_est = 100, 200, 300\nzip = structure.defect.radius, run.q_est\n")
    pts2 = expand_sweeps(s, p2)
    assert [(pt[0]["defect"]["radius"], pt[1]["run"]["q_est"]) for pt in pts2] == [(0.05, 100.0), (0.06, 200.0), (0.07, 300.0)]
    with pytest.raises(TextFileError):
        parse_params_text("[sweep]\nrun.mode = harminv, spectrum\n")
    with pytest.raises(TextFileError):
        parse_params_text("[sweep]\nstructure.defect.radius = 0.05\nzip = structure.defect.radius, analyte.n\n")


def test_row_sweep_and_override():
    s = parse_structure_text(REFERENCE)
    p = parse_params_text("[run]\nmode = harminv\n[sweep]\nstructure.rows.2 = 0.18, 0.19, 0.21, 0.22\n")
    pts = expand_sweeps(s, p)
    assert [pt[0]["rows"][2] for pt in pts] == [0.18, 0.19, 0.21, 0.22]
    k, v = parse_override("structure.defect.row=5")
    assert (k, v) == ("structure.defect.row", 5)
    with pytest.raises(TextFileError):
        parse_override("structure.defect.rows=5")


def test_field_mode_needs_f_res():
    with pytest.raises(TextFileError):
        parse_params_text("[run]\nmode = field\n")
    parse_params_text("[run]\nmode = field\n[field]\nf_res = 0.3046\n")
