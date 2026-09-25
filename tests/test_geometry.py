"""The new geometry builder must reproduce the rod list of the original FDTD script.

``fixtures/legacy_rods.json`` was generated once from the original ``02_kavite.py``
(``yapi_kur``) of the 2026 campaign for nine geometries; the test compares cell sizes and
rod lists exactly.
"""
import json
import os

import pytest

from sprat import geometry as geo
from sprat.params import parse_params_text, resolve
from sprat.structure import parse_structure_text

HERE = os.path.dirname(__file__)
FIX = json.load(open(os.path.join(HERE, "fixtures", "legacy_rods.json")))


def structure_for(p):
    s = parse_structure_text("[lattice]\n")
    s["cell"].update(guide_periods=p.get("nx", 25), cladding_rows=p.get("nkaplama", 12), pad_x=p.get("pad", 1.5),
                     pad_y=p.get("ypad", 1.0),
                     termination_x="absorber" if p.get("sonlandirma") == "emici" else "pml",
                     absorber_periods=p.get("emici_kalinlik", 8.0))
    s["defect"].update(radius=p.get("rd", 0.060), row=p.get("nsep", 4), dx=p.get("dx", 0.0), dy=p.get("dy", 0.0))
    s["rows"] = {int(k): v for k, v in p.get("bariyer_satirlari", {}).items()}
    return s


@pytest.mark.parametrize("name", sorted(FIX))
def test_rod_list_matches_legacy(name):
    case = FIX[name]
    p = case["params"]
    s = structure_for(p)
    mode = "reference" if p.get("kavite") is False else "harminv"
    par = resolve(s, parse_params_text(f"[run]\nmode = {mode}\n[numerics]\npml = {p.get('dpml', 1.0)}\n"))
    g = geo.build(s, par)
    assert (g.sx, g.sy) == (case["sx"], case["sy"])
    rods = sorted((round(x, 6), round(y, 6), round(r, 6)) for x, y, r in g.rods)
    assert rods == [tuple(r) for r in case["rods"]]


def test_symmetry_rule_and_probe():
    s = parse_structure_text("[lattice]\n")
    par = resolve(s, parse_params_text("[run]\nmode = harminv\n"))
    g = geo.build(s, par)
    assert g.mirror_x and geo.symmetry_sector(s, par, g) == ("even", 1)
    assert geo.probe_point(s, par) == (0.13, 4.07)
    s["defect"]["dx"] = 0.04
    g2 = geo.build(s, par)
    assert not g2.mirror_x and geo.symmetry_sector(s, par, g2) == ("none", 0)
    s["defect"]["dx"] = 0.0
    par["numerics"]["symmetry"] = "odd"
    assert geo.symmetry_sector(s, par, g) == ("odd", -1)
    assert geo.source_and_flux_x(g) == (-15 + 1 + 0.5, -15 + 2.5, 15 - 1.5)


def test_free_rod_edits():
    s = parse_structure_text("[lattice]\n[rods]\nremove 0 4\nadd 0 4 0.1\nremove 7.5 7\n")
    par = resolve(s, parse_params_text("[run]\nmode = harminv\n"))
    g = geo.build(s, par)
    assert any(abs(x) < 1e-9 and abs(y - 4) < 1e-9 and abs(r - 0.1) < 1e-9 for x, y, r in g.rods)
    assert not any(abs(x) < 1e-9 and abs(y - 4) < 1e-9 and abs(r - 0.06) < 1e-9 for x, y, r in g.rods)
    assert g.notes == ["remove 7.5 7.0: no rod there"]
