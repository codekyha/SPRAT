# Record schema `sprat-record-1.0`

One JSON file per run, named after the task label (`records/<label>.json`); a field run also
writes `<label>.npz` beside it. Everything a run needs is inside the record, so any record can be
re-run from itself (`sprat run-one --task <record.json> --out <new.json>`), and two records describe the same task exactly when their `structure` and
`params` blocks agree (labels, tags and the quiet/overwrite flags excluded).

| block | keys | meaning |
|---|---|---|
| `schema` | | `sprat-record-1.0` |
| `software` | `name`, `version`, `meep`, `python`, `code_sha256` | the code that wrote the record; `code_sha256` is the SHA-256 of `sprat/fdtd.py`; imported records carry the checksum of the original script when it was recorded |
| `task` | `label`, `tag`, `created`, `structure_file`, `params_file`, `overrides`, `source` | how the task was defined; `source` is `sprat`, `legacy-campaign` (the systematic runs of 2026; the value is kept as written by earlier versions) or `legacy-verification` (the verification runs); imported records add `legacy_layer` (the deposit directory the file came from) and, for verification runs, `legacy_wp` (the work-package code of 2026, `H1a` to `H13`, and `H14b` for the spectra with the guide continued into the absorber), and `software` gains `legacy_script` (the script that wrote the original file) |
| `structure` | the resolved `.phc` content | `lattice`, `cell`, `waveguide`, `defect`, `rows` (keys are row indices as strings), `rods` |
| `params` | the resolved `.par` content | `run`, `analyte`, `numerics`, `source`, `harminv`, `spectrum`, `field`, `bands`, `pwe` (no `sweep` block: sweeps are expanded before a task exists) |
| `provenance` | `host`, `started`, `wall_s`, `steps`, `sim_time`, `pixels`, `pixel_steps`, `throughput`, `symmetry`, `slot`, `job_id` | where and how long the run took; `pixel_steps` counts half the pixels when a mirror sector was used; imported records add `legacy_file`, `legacy_sha256`, `legacy_package`, `array_id` |
| `result` | mode-specific | see below |

## `result` by mode

**harminv**: `modes` (list of `{f, Q, amplitude, error, wavelength_nm, fwhm_nm}`, strongest amplitude first; `modes[0]` is the reported resonance), `t_used`, `auto_history` (the passes of the two-pass rule), `Q_limit` = $\pi f_{\rm cen} t$.

**spectrum**, **reference**: `f` (frequencies), `flux_in`, `flux_out`, `t_end`, `source_x`, `flux_in_x`, `flux_out_x`. The transmission of a spectrum is `flux_out / flux_out(reference at the same analyte index)`; `flux_out / flux_in` collapses at resonance because the input plane sees the reflected wave.

**field**: `analyte_energy_fraction` (the masked estimator), `wavelength_nm`, `S_first_order_nm_per_RIU` = $\lambda_r \eta / n_a$, `field_file`. The `.npz` holds `Ez` (complex64), `eps` (float32), `sx`, `sy`, `Lx`, `Ly` (the interior without the absorbing layers), `pml`, `edge_x`, `resolution`, `f_res`, `wavelength_nm`, `n`, `rod_eps`; axis 0 is $x$, axis 1 is $y$.

**bands**: `task` (`bulk`, `sweep`, `w1`, `calibration`), `gap` = `{exists, lower, upper, mid, relative_width_pct}`, `a_nm` (bulk and calibration), `tm`/`te` bands (bulk), `sweep` (list of `{n, gap}`), `kx_over_2pi`, `freqs`, `guided_bands` (w1).

**pwe**: `bulk_gap`, `W1` (`kx_over_pi`, `f`), `kappa` (`reference`, `second`: `kappa_pred`, `kappa_pred_zone_edge`, `n_g`, `beta_over_pi`, `kappa_at_kx0`, `kappa_at_kxpi`), `kappa_curve`, `reciprocity`, `fabry_perot`, `barrier_radius`, `loss`, `pred_Q`.

## The flat table (`sprat collect`)

`records.csv` / `records.jsonl`: one row per record with the columns of `sprat.records.FLAT_COLUMNS`; for harmonic-inversion records the validity flags of the paper: `margin` = $Q_{\rm lim}/Q$, `gap_position` = $(f_r - f_{\rm lower})/(f_{\rm upper} - f_{\rm lower})$ with the TM gap edges interpolated at the record's analyte index in the MPB table (`sprat.records.GAP_EDGES_MPB`), `resolved` (margin $\ge 1$), `valid` (resolved and gap position in 0.2 to 0.8), `exclusion_reason` (`Q > Q_lim`, `gap position outside 0.2-0.8`, `no mode`).
