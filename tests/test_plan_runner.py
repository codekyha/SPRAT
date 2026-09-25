import json
import math
import os

import pytest

from sprat import cost, plan, records, runner
from sprat.params import parse_params_text, resolve
from sprat.structure import default_structure

REF = """[lattice]
a_nm = 481.4
[cell]
guide_periods = 25
cladding_rows = 12
"""


def test_cost_model_reference_run():
    # the reference computation: 30a x 29a, res 24, t = 17 750 plus the 166.7 source time, mirror symmetry
    s = default_structure()
    p = resolve(s, parse_params_text("[run]\nmode = harminv\nq_est = 6927\n[harminv]\nt = 17750\n"))
    est = cost.estimate(s, p)
    assert (est["sx"], est["sy"]) == (30.0, 29.0) and est["symmetric"]
    T = 17750 + 2 * 5 / 0.06
    assert math.isclose(est["pixel_steps"], 0.5 * 2 * 30 * 29 * T * 24 ** 3)
    assert math.isclose(est["rru"], est["pixel_steps"] / 2.16e11, rel_tol=1e-9)
    assert 0.98 < est["rru"] < 1.02                       # one RRU by construction of the reference
    assert math.isclose(cost.seconds(est, 3.19e8), est["pixel_steps"] / 3.19e8)
    assert cost.makespan([10, 10, 10, 10], 2) == 20 and cost.makespan([13.6, 1, 1, 1], 4) == 13.6


def test_memory_and_no_symmetry_for_spectrum():
    s = default_structure()
    p = resolve(s, parse_params_text("[run]\nmode = spectrum\n"))
    est = cost.estimate(s, p)
    assert not est["symmetric"] and est["memory_mb"] > 30


def test_plan_writes_tasks(tmp_path):
    sf, pf = tmp_path / "s.phc", tmp_path / "p.par"
    sf.write_text(REF)
    pf.write_text("[run]\nmode = harminv\nq_est = 500\n[numerics]\nresolution = 8\n[harminv]\nt = 300\n"
                  "[sweep]\nstructure.defect.radius = 0.05, 0.06, 0.07\nstructure.cell.cladding_rows = 6, 8\n")
    tasks = plan.build_tasks(str(sf), str(pf), ["analyte.n=1.40"], records_dir=str(tmp_path / "records"), throughput=3e8)
    assert len(tasks) == 6 and all(t["params"]["analyte"]["n"] == 1.40 for t in tasks)
    assert len({t["label"] for t in tasks}) == 6
    files = plan.write_tasks(tasks, str(tmp_path / "tasks"), 3e8, workers=2)
    assert os.path.exists(files["jsonl"]) and os.path.exists(files["tsv"]) and os.path.exists(files["summary"])
    back = plan.read_tasks(files["jsonl"])
    assert back[0]["label"] == tasks[0]["label"]
    text = open(files["summary"]).read()
    assert "tasks: 6" in text


def test_plan_rejects_duplicate_labels(tmp_path):
    sf, pf = tmp_path / "s.phc", tmp_path / "p.par"
    sf.write_text(REF)
    pf.write_text("[run]\nmode = harminv\nlabel = fixed\n[sweep]\nstructure.defect.radius = 0.05, 0.06\n")
    with pytest.raises(ValueError):
        plan.build_tasks(str(sf), str(pf), records_dir=str(tmp_path / "r"), throughput=3e8)


@pytest.fixture
def fake_env(monkeypatch):
    monkeypatch.setenv("SPRAT_FAKE_FDTD", "1")
    monkeypatch.setenv("SPRAT_FAKE_SLEEP", "0.02")


def test_runner_end_to_end_with_resume(tmp_path, fake_env):
    sf, pf = tmp_path / "s.phc", tmp_path / "p.par"
    sf.write_text(REF)
    pf.write_text("[run]\nmode = harminv\nq_est = 500\ntag = t\n[numerics]\nresolution = 8\n[harminv]\nt = 300\n"
                  "[sweep]\nstructure.defect.radius = 0.05, 0.06, 0.07, 0.08\n")
    rdir = tmp_path / "records"
    tasks = plan.build_tasks(str(sf), str(pf), records_dir=str(rdir), throughput=3e8)
    files = plan.write_tasks(tasks, str(tmp_path / "tasks"), 3e8)
    s = runner.run(files["jsonl"], workers=2, pin_cores=True, log_dir=str(tmp_path / "logs"), poll_s=0.05, quiet=True)
    assert (s.done, s.skipped, s.failed) == (4, 0, 0)
    recs = records.load_dir(str(rdir))
    assert len(recs) == 4 and all(r["schema"] == "sprat-record-1.0" for r in recs)
    r0 = recs[0]
    assert r0["result"]["modes"][0]["Q"] == 500 and r0["provenance"]["throughput"] > 0
    assert r0["task"]["tag"] == "t" and r0["software"]["meep"] == "fake"
    with open(s.joblog) as fh:
        lines = fh.read().strip().splitlines()
    assert len(lines) == 5 and all(l.split("\t")[2] == "0" for l in lines[1:])
    # resume: everything is skipped
    s2 = runner.run(files["jsonl"], workers=2, log_dir=str(tmp_path / "logs"), poll_s=0.05, quiet=True)
    assert (s2.done, s2.skipped) == (0, 4)
    st = runner.status(files["jsonl"], str(tmp_path / "logs"))
    assert st["done"] == 4 and st["pending"] == 0


def test_runner_detects_label_collision(tmp_path, fake_env):
    sf, pf = tmp_path / "s.phc", tmp_path / "p.par"
    sf.write_text(REF)
    pf.write_text("[run]\nmode = harminv\nlabel = same\nq_est = 500\n[numerics]\nresolution = 8\n[harminv]\nt = 300\n")
    rdir = tmp_path / "records"
    tasks = plan.build_tasks(str(sf), str(pf), records_dir=str(rdir), throughput=3e8)
    files = plan.write_tasks(tasks, str(tmp_path / "tasks"), 3e8)
    runner.run(files["jsonl"], workers=1, log_dir=str(tmp_path / "logs"), poll_s=0.05, quiet=True)
    tasks2 = plan.build_tasks(str(sf), str(pf), ["structure.defect.radius=0.07"], records_dir=str(rdir), throughput=3e8)
    files2 = plan.write_tasks(tasks2, str(tmp_path / "tasks2"), 3e8)
    s = runner.run(files2["jsonl"], workers=1, log_dir=str(tmp_path / "logs"), poll_s=0.05, quiet=True)
    assert s.collisions == 1 and s.done == 0


def test_runner_retries_and_reports_failure(tmp_path, fake_env):
    sf, pf = tmp_path / "s.phc", tmp_path / "p.par"
    sf.write_text(REF)
    pf.write_text("[run]\nmode = field\nq_est = 500\n[numerics]\nresolution = 8\n[field]\nf_res = 0.3\n")
    rdir = tmp_path / "records"
    tasks = plan.build_tasks(str(sf), str(pf), records_dir=str(rdir), throughput=3e8)
    tasks[0]["params"]["run"]["mode"] = "broken"          # the worker raises ValueError on an unknown mode
    files = plan.write_tasks(tasks, str(tmp_path / "tasks"), 3e8)
    s = runner.run(files["jsonl"], workers=1, log_dir=str(tmp_path / "logs"), poll_s=0.05, retry=1, quiet=True)
    assert s.failed == 1 and s.done == 0 and s.failed_labels == [tasks[0]["label"]]
    with open(s.joblog) as fh:
        assert len(fh.read().strip().splitlines()) == 3        # header + two attempts


def test_cli_run_one_and_docs(tmp_path, fake_env):
    from sprat.cli import main
    sf, pf = tmp_path / "s.phc", tmp_path / "p.par"
    sf.write_text(REF)
    pf.write_text("[run]\nmode = reference\n[numerics]\nresolution = 8\n[spectrum]\nnfreq = 11\n")
    tasks = plan.build_tasks(str(sf), str(pf), records_dir=str(tmp_path / "records"), throughput=3e8)
    tpath = tmp_path / "t.json"
    tpath.write_text(json.dumps(tasks[0]))
    with pytest.raises(SystemExit) as e:
        main(["run-one", "--task", str(tpath)])
    assert e.value.code == 0
    rec = records.load(tasks[0]["record"])
    assert len(rec["result"]["f"]) == 11 and rec["provenance"]["symmetry"] == "none"
    main(["docs", "-o", str(tmp_path / "docs")])
    assert (tmp_path / "docs" / "PARAMETER_FILE.md").exists()


def test_check_runs():
    from sprat import doctor
    assert doctor.check("records") in (0, 1)


def test_run_one_reruns_a_record_from_itself(tmp_path, monkeypatch):
    """A record carries its structure and parameters; `sprat run-one --task <record>` re-runs it."""
    import json
    from sprat import execute
    from sprat.records import new_record, same_task
    from sprat.structure import default_structure, structure_to_json
    from sprat.params import default_params, resolve
    monkeypatch.setenv("SPRAT_FAKE_FDTD", "1")
    s = default_structure()
    p = default_params(); p.pop("sweep", None); p["run"]["mode"] = "harminv"; p["run"]["label"] = "auto"
    rp = resolve(s, p)
    rec = new_record(s, rp, task=dict(label=rp["run"]["label"], tag="", created="", structure_file="", params_file="",
                                       overrides={}, source="sprat"))
    src = tmp_path / "orig.json"
    src.write_text(json.dumps(dict(rec, structure=structure_to_json(s))), encoding="utf-8")
    out = tmp_path / "rerun.json"
    assert execute.run_task_file(str(src), out=str(out)) == 0
    new = json.load(open(out, encoding="utf-8"))
    assert same_task(rec, new) and new["result"]["modes"][0]["Q"] > 0
