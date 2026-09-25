"""Command-line entry point: ``sprat <command> ...``."""

from __future__ import annotations

import argparse
import json
import os
import sys

from . import __version__


def _p_plan(a):
    from . import plan
    from .cost import read_calibration
    cal = read_calibration(a.records)
    tasks = plan.build_tasks(a.structure, a.params, a.set, records_dir=a.records, throughput=a.throughput or cal["throughput"])
    files = plan.write_tasks(tasks, a.out, a.throughput or cal["throughput"], a.workers)
    print(open(files["summary"], encoding="utf-8").read())
    print(f"wrote {files['jsonl']}, {files['tsv']}, {files['summary']}")


def _p_run(a):
    from . import runner
    s = runner.run(a.tasks, workers=a.workers, pin_cores=not a.no_pin, cores_per_task=a.cores_per_task, order=a.order,
                   overwrite=a.overwrite, retry=a.retry, dry_run=a.dry_run, log_dir=a.logs,
                   memory_guard=not a.no_memory_guard, quiet=a.quiet)
    sys.exit(1 if s.failed else 0)


def _p_run_one(a):
    from .execute import run_task_file
    sys.exit(run_task_file(a.task, overwrite=a.overwrite, out=a.out))


def _p_status(a):
    from . import runner
    print(json.dumps(runner.status(a.tasks, a.logs), indent=1))


def _p_check(a):
    from . import doctor
    sys.exit(doctor.check(a.records))


def _p_calibrate(a):
    from . import doctor
    doctor.calibrate(a.records, resolutions=[int(x) for x in a.resolutions.split(",")], t=a.t)


def _p_collect(a):
    from .analysis import collect
    collect.collect_cli(a.records, a.out)


def _p_analyze(a):
    from .analysis import pipeline
    pipeline.analyze(a.records, a.out, pwe_file=a.pwe, predictions_dir=a.predictions)


def _p_grade(a):
    from .analysis import predictions
    predictions.grade_cli(a.records, a.out, a.predictions)


def _p_audit(a):
    from .analysis import audit
    audit.audit_cli(a.records, a.out, expected=a.expected)


def _p_figures(a):
    from .figures import build
    build.build_all(a.records, a.out, which=a.which, tables=a.tables)


def _p_import_legacy(a):
    from . import legacy
    legacy.import_cli(a.deposit, a.out)


def _p_deposit(a):
    from . import deposit
    deposit.build_cli(a.records, a.out, raw_legacy=a.raw_legacy, tables=a.tables, version=a.version, include=a.include,
                      data_doi=a.data_doi, software_doi=a.software_doi, raw_h14b=a.raw_h14b, arxiv_id=a.arxiv_id)


def _p_reproduce(a):
    from . import reproduce
    reproduce.reproduce_cli(a.campaign, a.records, tier=a.tier, workers=a.workers, dry_run=a.dry_run)


def _p_selftest(a):
    from . import doctor
    sys.exit(doctor.selftest(a.records, reference=a.reference, workers=a.workers))


def _p_docs(a):
    from .params import params_docs
    from .structure import structure_docs
    os.makedirs(a.out, exist_ok=True)
    for name, text in (("STRUCTURE_FILE.md", structure_docs()), ("PARAMETER_FILE.md", params_docs())):
        with open(os.path.join(a.out, name), "w", encoding="utf-8") as fh:
            fh.write(text)
        print("wrote", os.path.join(a.out, name))


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="sprat", description="SPRAT: Side-coupled Photonic-crystal Resonator Analysis Toolkit")
    ap.add_argument("--version", action="version", version=f"sprat {__version__}")
    sub = ap.add_subparsers(dest="command", required=True)

    p = sub.add_parser("plan", help="expand a structure file and a parameter file into a task list with cost estimates")
    p.add_argument("structure"); p.add_argument("params")
    p.add_argument("-o", "--out", default="tasks", help="output prefix (writes <out>.jsonl, <out>.tsv, <out>.summary.md)")
    p.add_argument("--records", default="records", help="records directory (default records)")
    p.add_argument("--set", action="append", default=[], metavar="KEY=VALUE", help="override, e.g. structure.defect.row=5")
    p.add_argument("--throughput", type=float, default=None, help="pixel-steps/s per task for the estimate")
    p.add_argument("--workers", type=int, default=None); p.set_defaults(func=_p_plan)

    p = sub.add_parser("run", help="execute a task list in parallel on this machine")
    p.add_argument("tasks", help="tasks.jsonl written by `sprat plan`")
    p.add_argument("--workers", type=int, default=None, help="concurrent tasks (default: physical cores)")
    p.add_argument("--cores-per-task", type=int, default=1); p.add_argument("--no-pin", action="store_true")
    p.add_argument("--order", choices=["longest", "shortest", "file"], default="longest")
    p.add_argument("--overwrite", action="store_true"); p.add_argument("--retry", type=int, default=1)
    p.add_argument("--dry-run", action="store_true"); p.add_argument("--logs", default="logs")
    p.add_argument("--no-memory-guard", action="store_true"); p.add_argument("--quiet", action="store_true")
    p.set_defaults(func=_p_run)

    p = sub.add_parser("run-one", help="execute one task file, or re-run one record from itself")
    p.add_argument("--task", required=True, help="a task file written by `sprat plan`, or a record (*.json) to re-run from itself")
    p.add_argument("--out", default=None, help="where the new record goes (default <label>_rerun.json beside a record; the task's own path for a task)")
    p.add_argument("--overwrite", action="store_true"); p.set_defaults(func=_p_run_one)

    p = sub.add_parser("status", help="progress of a task list")
    p.add_argument("tasks"); p.add_argument("--logs", default="logs"); p.set_defaults(func=_p_status)

    p = sub.add_parser("check", help="environment check: Python, numpy/scipy, Meep, cores, memory, disk, calibration")
    p.add_argument("--records", default="records"); p.set_defaults(func=_p_check)

    p = sub.add_parser("calibrate", help="measure the throughput of this machine (pixel-steps per second per task)")
    p.add_argument("--records", default="records"); p.add_argument("--resolutions", default="12,24")
    p.add_argument("--t", type=float, default=2000.0, help="harminv time of the calibration runs"); p.set_defaults(func=_p_calibrate)

    p = sub.add_parser("collect", help="flatten the records of a directory into records.csv / records.jsonl")
    p.add_argument("records"); p.add_argument("-o", "--out", default="tables"); p.set_defaults(func=_p_collect)

    p = sub.add_parser("analyze", help="run the whole analysis layer on a records directory")
    p.add_argument("records"); p.add_argument("-o", "--out", default="tables")
    p.add_argument("--pwe", default=None, help="pwe record or legacy pwe JSON to use (default: the newest pwe record)")
    p.add_argument("--predictions", default=None, help="directory with the frozen criterion sets"); p.set_defaults(func=_p_analyze)

    p = sub.add_parser("grade", help="grade the criteria fixed before the runs against the records")
    p.add_argument("records"); p.add_argument("-o", "--out", default="tables"); p.add_argument("--predictions", default=None)
    p.set_defaults(func=_p_grade)

    p = sub.add_parser("audit", help="independent recomputation of the headline numbers, optionally against an expected registry")
    p.add_argument("records"); p.add_argument("-o", "--out", default="tables"); p.add_argument("--expected", default=None)
    p.set_defaults(func=_p_audit)

    p = sub.add_parser("figures", help="build the figures of the paper from the records")
    p.add_argument("records"); p.add_argument("-o", "--out", default="figures")
    p.add_argument("--which", default="all", help="comma-separated figure numbers, or all")
    p.add_argument("--tables", default=None, help="analysis output directory (default <out>/../tables or run analyze)")
    p.set_defaults(func=_p_figures)

    p = sub.add_parser("import-legacy", help="convert the 2026 deposit (Turkish keys) into records")
    p.add_argument("deposit", help="directory of the unpacked Zenodo data pack (records_all.csv, *.tar.gz ...)")
    p.add_argument("-o", "--out", default="records"); p.set_defaults(func=_p_import_legacy)

    p = sub.add_parser("deposit", help="build the Zenodo data pack from a records directory")
    p.add_argument("records"); p.add_argument("-o", "--out", default="deposit")
    p.add_argument("--raw-legacy", default=None, help="directory with the original legacy files to include byte for byte")
    p.add_argument("--raw-h14b", default=None, help="directory with the absorber-terminated spectra, their job files and grading")
    p.add_argument("--arxiv-id", default=None, metavar="ID",
                   help="arXiv identifier of the preprint the record supplements (e.g. 2609.12345); adds the related identifier")
    p.add_argument("--tables", default=None, help="analysis output directory to copy into analysis/ (default: run the analysis)")
    p.add_argument("--version", default="3.1.0", help="version of the data pack (default 3.1.0)")
    p.add_argument("--include", nargs="*", default=None, metavar="FILE",
                   help="extra files to place in the pack before the manifest is written (for example the supplementary PDF)")
    p.add_argument("--data-doi", default="10.5281/zenodo.NNNNNNN", metavar="DOI",
                   help="DOI of the data record (default: the placeholder 10.5281/zenodo.NNNNNNN)")
    p.add_argument("--software-doi", default="10.5281/zenodo.SSSSSSS", metavar="DOI",
                   help="DOI of the SPRAT release (default: the placeholder 10.5281/zenodo.SSSSSSS); "
                        "a different real DOI also becomes the related identifier isCompiledBy of the metadata; "
                        "give the data DOI itself when SPRAT is archived in the same record (--include sprat-<version>.zip); "
                        "with such an archive included and no software DOI given, the data DOI is used")
    p.set_defaults(func=_p_deposit)

    p = sub.add_parser("reproduce", help="plan (and run) a campaign of the campaigns/ directory")
    p.add_argument("campaign", help="campaign directory, e.g. campaigns/manuscript")
    p.add_argument("--records", default="records"); p.add_argument("--tier", default="full", choices=["quick", "core", "full"])
    p.add_argument("--workers", type=int, default=None); p.add_argument("--dry-run", action="store_true"); p.set_defaults(func=_p_reproduce)

    p = sub.add_parser("selftest", help="Meep smoke test at low resolution; --reference adds the regression of the reference record (about 15 minutes)")
    p.add_argument("--records", default="selftest_records"); p.add_argument("--reference", action="store_true")
    p.add_argument("--workers", type=int, default=2); p.set_defaults(func=_p_selftest)

    p = sub.add_parser("docs", help="write the schema reference pages for the two text files")
    p.add_argument("-o", "--out", default="docs"); p.set_defaults(func=_p_docs)
    return ap


def main(argv=None) -> None:
    ap = build_parser()
    a = ap.parse_args(argv)
    try:
        a.func(a)
    except FileNotFoundError as e:                      # a missing input file is a usage error, not a bug
        ap.exit(2, f"sprat {a.command}: file not found: {e.filename or e}\n")
    except ValueError as e:                             # TextFileError and the parameter/structure checks
        ap.exit(2, f"sprat {a.command}: {e}\n")


if __name__ == "__main__":
    main()
