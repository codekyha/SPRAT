"""``sprat analyze``: run the whole analysis layer on a records directory and write the outputs.

Outputs (in the ``tables`` directory): ``records.csv`` / ``records.jsonl`` / ``counts.json``,
``analysis.json`` (every series and derived quantity), ``numbers_registry.json`` and ``.md``,
``predictions_report.md`` / ``predictions.json`` and ``ANALYSIS_REPORT.md``.
"""

from __future__ import annotations

import json
import os
import time
from typing import Any

from .. import records
from . import barrier, census, collect, combined, fano, fields, predictions, registry, sensitivity, termination


def portable_path(path: str) -> str:
    """A path for the outputs: relative to the working directory when inside it, else the last component, so that no
    output carries the absolute path of the machine it was built on."""
    ap = os.path.abspath(path)
    try:
        rel = os.path.relpath(ap, os.getcwd())
    except ValueError:                      # another drive (Windows)
        rel = None
    if rel is not None and rel != os.pardir and not rel.startswith(os.pardir + os.sep) and not os.path.isabs(rel):
        return rel.replace(os.sep, "/")
    return os.path.basename(os.path.normpath(ap))


def _pwe_result(recs: list[dict], pwe_file: str | None) -> dict | None:
    if pwe_file:
        with open(pwe_file, encoding="utf-8") as fh:
            d = json.load(fh)
        return d["result"] if d.get("schema") else d
    cands = [r for r in recs if r["params"]["run"]["mode"] == "pwe"]
    if not cands:
        return None
    cands.sort(key=lambda r: (r["provenance"].get("started") or "", r["task"]["label"]))
    return cands[-1]["result"]


def _pwe_legacy(recs: list[dict]) -> dict | None:
    """The plane-wave layer of the campaign (the earliest pwe record), kept for the superseded statements."""
    cands = [r for r in recs if r["params"]["run"]["mode"] == "pwe"]
    if len(cands) < 2:
        return None
    cands.sort(key=lambda r: (r["provenance"].get("started") or "", r["task"]["label"]))
    return cands[0]["result"]


def _legacy_channels(pwe_l: dict, A: dict) -> dict:
    """The channel values the campaign layer reported (every root kept): at the reference resonance, and on its kappa
    curve at the mean resonance frequency of the second radius."""
    from .. import pwe as pwe_mod
    k = pwe_mod.kappa_block(pwe_l, "reference")
    out = dict(reference=dict(slowest=k["kappa_pred"], zone_edge=k["kappa_pred_zone_edge"], kx0=k.get("kappa_at_kx0"), kxpi=k.get("kappa_at_kxpi")),
               eps_rod=(pwe_l.get("model") or {}).get("eps_rod"))
    fbar = ((A.get("barrier") or {}).get("second_radius") or {}).get("f_mean")
    if fbar:
        cs, ce = pwe_mod.channels_at(pwe_l, fbar)
        out["second"] = dict(f=fbar, slowest=cs, zone_edge=ce)
    return out


def _try(name: str, fn, out: dict, log):
    try:
        out[name] = fn()
    except Exception as e:                 # noqa: BLE001
        out.setdefault("errors", {})[name] = f"{type(e).__name__}: {e}"
        log(f"  {name}: skipped ({type(e).__name__}: {e})")


def analyze(records_dir: str, out: str = "tables", pwe_file: str | None = None, predictions_dir: str | None = None,
            predictions_only: bool = False, criteria_path: str | None = None, log=print) -> dict[str, Any]:
    os.makedirs(out, exist_ok=True)
    recs = collect.load(records_dir)
    rows = collect.table(recs)
    counts = records.counts(rows)
    records.write_table(rows, os.path.join(out, "records.csv"), os.path.join(out, "records.jsonl"))
    with open(os.path.join(out, "counts.json"), "w", encoding="utf-8") as fh:
        json.dump(counts, fh, indent=1)
    log(f"records: {counts}")
    pwe = _pwe_result(recs, pwe_file)
    legacy_pwe = _pwe_legacy(recs)
    A: dict[str, Any] = dict(generated=time.strftime("%Y-%m-%d %H:%M:%S"), records_dir=portable_path(records_dir), counts=counts,
                             pwe_available=pwe is not None)
    _try("census", lambda: census.census(recs), A, log)
    _try("repeats", lambda: census.repeats(recs), A, log)
    _try("band_gap", lambda: census.band_gap(recs, pwe), A, log)
    _try("analyte_sweep", lambda: sensitivity.analyte_sweep(recs), A, log)
    S_meas = A.get("analyte_sweep", {}).get("at_1p33", {}).get("S")
    _try("radius_sweep", lambda: sensitivity.radius_sweep(recs), A, log)
    _try("base_series", lambda: sensitivity.base_series(recs), A, log)
    _try("base_series_all_points", lambda: sensitivity.base_series(recs, include_excluded=True), A, log)
    _try("radius_sweep_all_points", lambda: sensitivity.radius_sweep(recs, include_excluded=True), A, log)
    _try("kappa_series", lambda: sensitivity.kappa_series(recs), A, log)
    _try("resolved_pair", lambda: sensitivity.resolved_pair_kappa(recs), A, log)
    _try("cladding_ladder", lambda: sensitivity.cladding_ladder(recs), A, log)
    _try("displacement", lambda: sensitivity.displacement(recs), A, log)
    _try("preliminary_sweep", lambda: sensitivity.preliminary_sweep(recs), A, log)
    _try("convergence", lambda: sensitivity.convergence(recs), A, log)
    _try("spectra", lambda: fano.analyse_spectra(recs), A, log)
    if "spectra" in A:
        A["spectra_summary"] = fano.summary(A["spectra"])
    _try("termination", lambda: termination.analyse(recs), A, log)
    if A.get("termination") is None:
        A.pop("termination", None)
    else:
        _try("spectra_absorber", lambda: termination.absorber_spectra(recs), A, log)
    _try("fields", lambda: fields.analyse(recs, pwe, S_measured=S_meas), A, log)
    if pwe is not None:
        _try("barrier", lambda: barrier.analyse(recs, pwe, field_eta=(A.get("fields") or {}).get("reference"), S_measured=S_meas), A, log)
    else:
        log("  no pwe record: the barrier analysis and the channel comparisons are skipped (run a mode = pwe task)")
    if legacy_pwe is not None and legacy_pwe is not pwe:
        _try("pwe_legacy_channels", lambda: _legacy_channels(legacy_pwe, A), A, log)
    _try("combined", lambda: combined.combine(A), A, log)
    A["predictions"] = predictions.grade(recs, pwe, A.get("barrier"), A.get("fields"), criteria_path)
    reg = registry.build(A, counts)
    A["registry"] = reg
    if not predictions_only:
        with open(os.path.join(out, "analysis.json"), "w", encoding="utf-8") as fh:
            json.dump({k: v for k, v in A.items() if k not in ("spectra", "spectra_absorber")}, fh, indent=1, default=_default)
        with open(os.path.join(out, "spectra.json"), "w", encoding="utf-8") as fh:
            json.dump(A.get("spectra", []), fh, default=_default)
        if A.get("spectra_absorber"):
            with open(os.path.join(out, "spectra_absorber.json"), "w", encoding="utf-8") as fh:
                json.dump(A["spectra_absorber"], fh, default=_default)
        with open(os.path.join(out, "numbers_registry.json"), "w", encoding="utf-8") as fh:
            json.dump(reg, fh, indent=1, default=_default)
        with open(os.path.join(out, "numbers_registry.md"), "w", encoding="utf-8") as fh:
            fh.write(registry.markdown(reg))
    with open(os.path.join(out, "predictions.json"), "w", encoding="utf-8") as fh:
        json.dump(A["predictions"], fh, indent=1, default=_default)
    with open(os.path.join(out, "predictions_report.md"), "w", encoding="utf-8") as fh:
        fh.write(predictions.report_markdown(A["predictions"]))
    if not predictions_only:
        with open(os.path.join(out, "ANALYSIS_REPORT.md"), "w", encoding="utf-8") as fh:
            fh.write(report(A))
        log(f"wrote {out}/analysis.json, numbers_registry.json/.md, predictions_report.md, ANALYSIS_REPORT.md")
    return A


def _default(o):
    import numpy as np
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return str(o)


def report(A: dict[str, Any]) -> str:
    L = ["# Analysis report", "", f"generated {A['generated']} from {A['records_dir']}", ""]
    c = A["counts"]
    L += ["## Records", "", f"- harmonic-inversion records: {c['harminv']} ({c['valid']} admitted by the criterion margin >= 1 and gap position "
          f"0.2-0.8; {c['excluded_margin']} fail the margin, {c['excluded_gap']} lie outside the gap window)",
          f"- admitted records holding a mode of the finite waveguide (preliminary sweep, N_sep = 1 or r_d >= 0.13a; "
          f"column waveguide_mode): {c.get('waveguide_mode_admitted', 0)}, of which {c.get('dead_zone_admitted', 0)} lie in the dead zone "
          f"r_d >= 0.14a; admitted records holding a cavity mode: "
          f"{c['valid'] - c.get('waveguide_mode_admitted', 0)}",
          f"- spectra {c['spectra']}, cavity-less references {c['references']}, field maps {c['fields']}, band computations {c['bands']}, "
          f"plane-wave results {c['pwe']}", ""]
    sw = A.get("analyte_sweep")
    if sw:
        a = sw["at_1p33"]
        L += ["## Sensitivity (reference geometry, analyte sweep)", "",
              f"- {len(sw['n'])} points, n_a = {sw['n'][0]:.3f} to {sw['n'][-1]:.3f}",
              f"- S(1.33) = {a['S']:.2f} nm/RIU (central difference); single linear slope {sw['chord']['S']:.2f} nm/RIU, R^2 = {sw['chord']['R2']:.6f}",
              f"- Q(1.33) = {a['Q']:.2f}, fwhm = {a['fwhm_nm']:.4f} nm, FOM = {a['FOM']:.1f} 1/RIU, DL (R = fwhm/10) = {a['DL10']:.3e} RIU",
              f"- Q drops by a factor {sw['Q_drop_factor']:.2f} over the sweep", ""]
    ks = A.get("kappa_series")
    if ks:
        L += ["## Raw barrier series of the systematic runs (25a cell, PML)", "",
              "- N_sep " + ", ".join(str(n) for n in ks["N"]) + ": Q = " + ", ".join(f"{q:.1f}" for q in ks["Q"]) +
              " (margins " + ", ".join(f"{m:.2f}" for m in ks["margin"]) + ")",
              f"- kappa = {ks['kappa']:.4f} +- {ks['kappa_sd']:.4f} (per-row factor {ks['per_row_factor']:.2f})", ""]
    ba = A.get("barrier")
    if ba:
        s, l, g = ba["series"], ba["leak"], ba["G"]
        L += ["## Barrier analysis (verification runs)", "",
              f"- three-point solve at clearance 8: kappa = {ba['three_point_solve']['kappa']:.4f}, Q_top(8) = {l['Q_top_8']:.0f}",
              f"- leak law: Q_top(7) = {l['Q_top_7']:.0f}, kappa_top = {l['kappa_top']:.3f} per cladding row",
              f"- cell factor G = Q_PML/Q_abs after leak removal: N = 4 {g['4']:.4f}, N = 5 {g['5']:.4f}",
              "- reflectionless leak-free series Q_w(N): " + ", ".join(f"N{k} {v:.1f}" for k, v in s["Q_w"].items()),
              f"- kappa = {s['kappa']:.4f} +- {s['kappa_sd_monte_carlo']:.4f} (Monte Carlo); per-row factor {s['per_row_factor']:.2f}; "
              f"zone-edge channel {ba['channels']['zone_edge']:.4f} ({s['dev_zone_edge_pct']:+.2f} %)"
              + (f"; next physical channel {ba['channels']['next']:.4f}" if ba['channels'].get('next') else "")
              + (f"; slowest truncation artefact of the square basis {ba['channels']['boundary_artefact']:.4f} (rejected)"
                 if ba['channels'].get('boundary_artefact') else ""), ""]
        mr = ba.get("matched_ratio_fit")
        if mr:
            L += [f"- two-reflector fit of G(r_d): R = {mr['R']:.3f} +- {mr['R_se']:.3f}, D = {mr['D']:.2f} +- {mr['D_se']:.2f} a, rms(ln G) = {mr['rms_lnG']:.3f}; "
                  f"absorber scale c = {mr['scale_c']:.3f} +- {mr['scale_c_se']:.3f}", ""]
        sr = ba.get("second_radius")
        if sr:
            L += [f"- second radius r_d = {sr['radius']}: kappa_absorber = {sr['kappa_absorber']:.4f}, kappa_PML = {sr['kappa_pml']:.4f} "
                  f"(difference {sr['kappa_difference']:.4f}); zone edge at the measured f {sr['channels_at_measured_f']['zone_edge']:.4f} ({sr['dev_zone_edge_pct']:+.2f} %)", ""]
        pf = ba.get("performance")
        if pf:
            L += [f"- in water: Q_abs = {pf['Q_abs']:.0f}, Q_tot(N4) = {pf['Q_tot']['N4']:.0f}, FOM lossless {pf['FOM_lossless_N4']:.0f}, "
                  f"in water {pf['FOM_water_N4']:.0f} 1/RIU; coupling optimum Q_w = {pf['Q_w_opt']:.0f}", ""]
    fe = A.get("fields")
    if fe and fe.get("reference"):
        e = fe["reference"]
        L += ["## Field map at n_a = 1.33", "",
              (f"- Hellmann-Feynman estimate over the permittivity E_z sees: eta_a = {e['eta_hf']:.4f}, S_th = {e['S_hf']:.1f} nm/RIU"
               + (f" ({e['dev_hf_pct']:+.2f} % against the measured {e['S_measured']:.2f})" if "dev_hf_pct" in e else "")
               + f"; pure-analyte points {e['eta_hf_low']:.4f}, with the mixed points {e['eta_hf_high']:.4f}; the stored map is reproduced to "
                 f"{e['reconstruction_max_abs']:.1e}") if e.get("eta_hf") is not None else "- no reconstruction of the smoothed permittivity",
              f"- stored-map estimators of the original analysis: mask {e['eta_mask']:.4f}, linear unmixing {e['eta_lin']:.4f}, "
              f"pure pixels {e['eta_low']:.4f}, geometric {e['eta_geo']:.4f}", ""]
    sp = A.get("spectra_summary")
    if sp:
        L += ["## Spectra", "", f"- {sp['count']} spectra; |q| from {sp['abs_q_min']:.3f} to {sp['abs_q_max']:.3f}; "
              f"fwhm(harminv)/fwhm(fit) from {sp['fwhm_ratio_min']:.4f} to {sp['fwhm_ratio_max']:.4f}", ""]
    pr = A.get("predictions")
    if pr:
        L += ["## Criteria fixed before the runs", "", "- " + ", ".join(f"{k} {v}" for k, v in sorted(pr["verdict_counts"].items())),
              "- see predictions_report.md", ""]
    if A.get("errors"):
        L += ["## Skipped analyses", ""] + [f"- {k}: {v}" for k, v in A["errors"].items()] + [""]
    return "\n".join(L)


__all__ = ["analyze", "report"]
