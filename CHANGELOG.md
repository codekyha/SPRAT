# Changelog

## Unreleased (documentation only)

- 29 September 2026: the preprint of the paper, arXiv:2609.31952 (physics.optics), and the published
  record 10.5281/zenodo.22912910 (version 3.1.0) are named in `README.md` (introduction, section 9,
  section 12 with BibTeX entries), `CITATION.cff` (preferred citation: the preprint; the record under
  `references`), `codemeta.json`, `.zenodo.json`, `pyproject.toml` (project URLs) and
  `docs/DATA_AVAILABILITY.md`. No code changed; the version stays 1.2.1 and the archived
  `sprat-1.2.1.zip` of the record is unaffected.

## 1.2.1 (2026-09-25)

The release that manuscript version 10.2 and version 3.1.0 of the data record rest on.
Every number that 1.2.0 computes from the records of the data record is unchanged; 1.2.1 adds the
reflectionless values across the analyte sweep, one analyte index under both guide terminations
in figure 6, the convergence of the resonance frequency with resolution, the channel at the
guided-mode wavevector and the absorption-limited quality factor with the mixed points of the
grid counted either way.

- **Reflectionless values across the analyte sweep.** `combined.physical_values` gives the
  reflectionless quality factor $Q_w$ and figure of merit ${\rm FOM}_w = S Q_w/\lambda_r$ at the
  rows of table 5 and along the sweep: at $n_a$ = 1.330 the harmonic-inversion value with the
  absorber and the cladding leak removed (4502), at the indices of the absorber spectra
  $Q_{\rm cell}/G$ (amendment A2 of the criterion of the absorber spectra), and between them the
  linear interpolation of $Q_w$ in $n_a$, which varies smoothly where $G$ does not. At 1.33 the
  interpolated spectral value (4455) lies 1.0 % below the harmonic-inversion one. Figure 7(c)
  draws $Q_w$, figure 7(d) ${\rm FOM}_w$.
- **Figure 6.** An inset of panel (a) shows $n_a$ = 1.375 under both guide terminations, each on
  its own normalised detuning; `sprat analyze` writes the absorber-terminated spectra to
  `spectra_absorber.json`.
- **Convergence of the resonance frequency.** `convergence.resolution_20_vs_32_frequency`
  (at most 0.023 % between resolutions 20 and 32 at the four radii of table 2) and
  `barrier.resolution_check` (0.029 %, 0.46 nm at 1550 nm, between 24 and 32 at the reference
  radius, $N_{\rm sep}$ = 5; the defect rod is 1.44 pixels in radius at resolution 24).
- **The channel at the guided-mode wavevector.** `barrier.channels` gains $\beta$ at the
  reference resonance ($0.5465\,\pi/a$) and the per-row factor the same frequency gives at
  $k_x = 0$ (3.93) and at $\pi/a$ (12.2), against 7.10 measured.
- **Absorption.** `performance.absorption` gives $Q_{\rm abs}$ with the mixed points of the
  grid counted wholly as analyte ($6.7\times10^3$) or wholly as silicon ($9.6\times10^3$).
- **Registry.** Section `reflectionless` (five rows); the frequency convergence in
  `systematic_runs` and `derived`, the channel in `derived`, the absorption bracket in
  `absorption`. The manuscript registry of version 10.2
  (`campaigns/manuscript/expected/numbers_registry_manuscript_v10_2.json`) is the registry of
  version 10 with the section `derived_v10_2` (14 rows), the three references added and the
  literature value taken from one of them.
- **Texts.** The predictions README and the metadata of the data record say that the third
  criterion set was amended while its runs were under way; the English transcription of the first
  set names the source of P12 and P13 as the frozen file does; the analysis report, the README and
  the docstrings call the analysis of the original scripts "the original analysis"; the entry of
  1.2.0 no longer says that version 3.1.0 of the data record rests on 1.2.0 (it rests on 1.2.1).
- **Tests.** Six more (73 in all): the numbers above on the full deposits, each at the precision
  the manuscript prints it and in its section of the registry, and the audit of the manuscript
  registry of version 10.2 (OK on all 237 record-derived rows).

## 1.2.0 (2026-09-24)

The release that manuscript version 10 rests on. Every record and every number that 1.1.0 computes
from the deposit of 2026 is unchanged; 1.2.0 adds the transmission spectra with the guide continued
into the absorber (with them the counts of spectra and references double) and their comparison with
the spectra of the 25a cell.

- **Absorber-terminated spectra (H14b).** `sprat import-legacy` reads the directory `v4_spectra/` of a
  deposit such as `raw_h14b`: the 14 runs of job 528207, seven spectra and seven cavity-less references
  at $n_a$ = 1.300 to 1.450 with the guide continued into an $8a$ adiabatic absorber (`task.source`
  `legacy-verification`, `legacy_wp` `H14b`). The import manifest lists every deposit imported into one
  directory; `records.counts` gains `spectra_absorber`.
- **`analysis/termination.py`.** The comparison of the two guide terminations: the Fano fit and the
  gates G1 to G5 and S1 of the grading scripts `48_fano_karsilastir_S4b.py` and `49_H14b_A1A2.py`,
  kept line for line; the grades by the criterion fixed before the runs (D-22: not a pass,
  inconclusive at $n_a$ = 1.300, factor 1.61) and by its amendment A1 (pass: factors 8.95 to 213.41 at
  the five indices where $|q|$ in the 25a cell exceeds 0.2), side by side; the cell factor
  $G = Q_{\rm cell}/Q_{\rm absorber}$ (1.653 to 0.607) and the reflectionless $Q_w$ (5037 to 2737, a
  fall of 1.84 against 5.01 in the cell) of amendment A2, with G at n_a = 1.33 from harmonic inversion
  (1.5411) and the opening of the linewidth over the sweep (1.95 without reflections, 5.31 in the cell). On the deposited spectra it gives the numbers
  of the grading scripts bit for bit, and it agrees with the comparison the job wrote on the cluster to
  $3\times10^{-7}$. `analysis.json` gains the key `termination`. The spectrum analysis of section 3.5
  (`spectra_summary`, `spectra.json`) keeps to the spectra of the 25a cell
  (`fano.analyse_spectra(..., termination="pml")`), so every number of 1.1.0 stays.
- **Registry.** The section `termination` (29 rows) and, in `dataset`, the share of the
  harmonic-inversion records admitted and holding a cavity mode (23.78 %) and the number of absorber
  spectra; "transmission spectra" counts 14 when the absorber spectra are present. Two quantities are
  renamed, "records from the systematic runs" and "eta_a, threshold mask of the original analysis",
  and the section of the systematic runs is `systematic_runs`; the derivation of the absorption cap
  says "near 9200", as the abstract does.
- **Manuscript registry of version 10**
  (`campaigns/manuscript/expected/numbers_registry_manuscript_v10.json`): version 9 with the sections
  `termination` and `corrected_v10`, the transmission spectra of the record (14), the absorption cap of
  the abstract (near 9200) and five derivation texts that still carried numbers of version 8 (the two
  vertical-channel rows, the reflectionless cross-term, the two cell-factor estimates), each changed row
  with its version-9 value or derivation; `corrected_v10` holds the share of the admitted cavity records
  and the fall of the quality factor over the analyte sweep (5.01), which section 3.6 printed without a
  row. The audit of the full record gives OK on all 223
  record-derived rows and CONSISTENT on all 21 superseded statements; `sprat deposit` copies it into the
  pack.
- **Criteria fixed before the runs.** `predictions/frozen/` gains the criterion of the
  absorber-terminated spectra and its amendment, byte for byte (`predictions_S4b_2026-09-24.json`,
  `predictions_S4b_amendment_A1A2_2026-09-24.json`); `predictions/README.md` is headed "Criteria fixed
  before the runs" and says that the file names are the names of the day.
- **Deposit.** `sprat deposit --raw-h14b` packs the absorber-terminated spectra with the job files, the
  criterion, the amendment and the grading as `raw_h14b.tar.gz`; `--arxiv-id` adds the preprint the
  record supplements (README; metadata relation *Is supplement to*); a `CHANGES_FROM_*.md` in the
  `--raw-legacy` directory is named in the texts of the pack. The texts of the pack speak of the
  systematic runs and the verification runs.
- **Task list of the data set.** `build_campaign.py --append` added the 14 runs of `raw_h14b` to
  `campaigns/manuscript/tasks.jsonl` (892 tasks; the 878 earlier lines unchanged); `sprat reproduce`
  prints the audit command with the registry of version 10.
- **Tests.** `tests/test_termination.py`: the fit, the identity gate, the constants against the frozen
  criterion and its amendment (with their SHA-256), and the index 1.300 under both terminations from the
  fixtures (`tests/fixtures/h14b_sample`); with `SPRAT_H14B_DEPOSIT`, every index against the comparison
  and the grading of the job. `tests/test_audit_manuscript.py` audits the registry of version 10 on the
  full record, and the registry of version 9 on its own records (the two renamed quantities apart).
  `tests/test_deposit_texts.py` checks the texts of the pack: the statement, the metadata with the
  absorber spectra and the preprint, the licence, and the version of the plane-wave record, which the
  pack now takes from the record (1.1.0) instead of the running SPRAT.
- **Documentation.** README (the absorber comparison, the data record with `raw_h14b.tar.gz`, the
  release order with the arXiv identifier), `docs/REPRODUCING_THE_PAPER.md`, `docs/COST_MODEL.md`,
  `docs/RECORD_SCHEMA.md` and the README of `campaigns/manuscript`; `docs/PARAMETER_FILE.md`
  regenerated from the code (the default measured $\kappa$ values of the plane-wave checks, 1.9595 and
  1.7386, had not reached it in 1.1.0). "Campaign" remains only in identifiers: the directory
  `campaigns/manuscript`, the module `figures/campaign.py`, the tier `campaign` of the task list and the
  record value `legacy-campaign`.
- **AI credit.** README, `CITATION.cff` and `.zenodo.json`: SPRAT was written by the author with coding
  assistance from Claude (Anthropic; Claude Opus 5.5 and Claude Fable 5) and Gemini (Google; Gemini 3.6
  Flash and Gemini 3.1 Pro), used in tandem.

## 1.1.0 (2026-09-24)

Corrections that manuscript version 9 and version 3.1.0 of the data record rest on. No record
changes; every change is in the analysis of the records, and every version-8 value is still
computed alongside for comparison.

- **Complex band structure, physical roots only.** The companion linearisation of the plane-wave
  problem in $k_y$ returns spurious roots in a truncated basis. `pwe.decay_roots` keeps a root when
  less than 1e-3 of its eigenvector weight lies on the two outermost $G_y$ rings and a replica has
  $|{\rm Re}\,k_y| \le \pi/a$; the rejected roots are reported (`boundary_roots`,
  `truncation_check`: three square bases and a circular one). The slowest physical channel at the
  guided-mode wavevector is the zone-edge root (1.9620 at the reference resonance, 1.7288 at the
  second radius); the "least-evanescent pair" of 1.0.0 (1.9015) is an artefact of the square basis.
  The next physical channel (6.23) and the channels at $k_x = 0$ and $\pi/a$ are reported; the
  legacy channels of the 2026 plane-wave file stay available as `pwe_legacy_channels`. A new
  plane-wave record at $n_a = 1.33$, $r = 0.20a$ is part of data record 3.1.0.
- **Admission.** Fits and statistics use the records admitted by the criterion of the paper (margin
  >= 1, gap position 0.2 to 0.8): the quadratic law $f_r(r_d)$ and its transfer to the fine sweep,
  the per-row steps and means of the $n_{\rm cl} = 8$ series and their leak correction, the
  two-reflector fit and the detrended oscillation. The all-points values of 1.0.0 are kept in
  `v8_all_points`, `base_series_all_points` and `radius_sweep_all_points`; the figures draw excluded
  records as open grey symbols.
- **Sampling rule.** $Q_{\rm lim} = \pi f_{\rm cen} t$ is described as a conservative sampling
  rule: harmonic inversion resolves decays longer than the signal, as the repeats at higher margin
  show.
- **Analyte energy fraction.** `fields.meep_smoothing` reconstructs Meep's smoothed permittivity
  tensor from the geometry (exact disc-voxel areas, the surface normal) and reproduces the stored
  map, which is the harmonic mean of the eigenvalues of that tensor, to 1e-5; `fields.hellmann_feynman`
  evaluates the analyte energy fraction with the permittivity $E_z$ sees, weighting each mixed point
  by its analyte fraction (0.5523, $S_{\rm th}$ = 656.3 nm/RIU, 1.5 % above the measurement). The
  stored-map estimators of 1.0.0 are kept and marked as biased upwards; the pure-pixel value is no
  bound.
- **Absorption.** $Q_{\rm abs} = \lambda_r/(2kS)$ from the measured sensitivity (9175), with the
  stored-map and Hellmann-Feynman variants alongside; figures of merit and the coupling metric
  follow.
- **Axial displacement.** `sensitivity.displacement` fits $Q(\delta_x)/Q(0) = 1 - c\,\delta_x^2$ with
  and without a free offset (the decrease of the 25a-cell value is second order, 0.39 % and 0.62 % at
  0.08a); `combined` estimates the change of the cell factor that the resonance shift causes (1.30 %
  and 1.18 %), which exceeds the measured change, and the registry states both.
- **Spectra.** The frequency step, the samples next to the minimum and the far-window transmission
  are reported; the smallest sampled transmission is set by the frequency grid.
- **Census and the modes of the waveguide.** The validity criterion tests the numerics alone. In the
  preliminary sweep ($n_{\rm cl} = 6$, resolution 20) the strongest mode of the 537 admitted records
  at $N_{\rm sep} = 1$ or $r_d \ge 0.13a$ is a low-$Q$ mode of the finite waveguide ($Q$ 50 to 77);
  436 of them lie in the dead zone ($r_d \ge 0.14a$). `records.waveguide_mode` states the rule, the
  new column `waveguide_mode` of `records.csv` flags the records, and the counts, the registry and
  the data availability table give 537 such records and 204 admitted cavity records. The ratio
  $Q(\delta_x = 0.08a)/Q(0)$ of that sweep is reported over its 15 cavity groups (0.984) beside the
  earlier value over all 48 groups.
- **Registry and audit.** New rows for every number of manuscript version 9 and the section
  `superseded_v9` (twelve statements of version 8 with both sides recomputed); the manuscript registry
  of version 9 (`campaigns/manuscript/expected/numbers_registry_manuscript_v9.json`) is audited with
  OK on every record-derived row and CONSISTENT on every superseded statement; `sprat deposit`
  copies it into the pack.
- **Figures.** Figure 3: excluded records as open grey circles, the means of panel (c) filled, the
  label of panel (d) clear of the curve; figure 4: panel (c) draws the physical channel only,
  panel (a) legend and annotations moved; figure 5(b) on an expanded scale with the quadratic fits.
- **Documentation.** README, `docs/REPRODUCING_THE_PAPER.md` and the README of `campaigns/manuscript` describe the
  admission, the physical-root filter, the analyte energy fraction and the plane-wave record needed
  after `sprat import-legacy`; `CITATION.cff` and `.zenodo.json` name SPRAT 1.1.0 and record 3.1.0.

## 1.0.0 (2026-09-23)

First release. English-only rewrite of the original code of 2026 (Turkish identifiers, SLURM and
cluster-specific material removed) as a standalone package:

- structure file (`.phc`) and parameter file (`.par`) with schema validation, sweeps and overrides;
- geometry builder tested against the rod lists of the original script for nine geometries;
- four FDTD modes (harminv, spectrum, reference, field), MPB bands (bulk, sweep, W1, calibration) and the plane-wave layer;
- cost model with per-machine calibration; task planner; workstation parallel runner with pinning, resume, retry, job log and ETA;
- English record schema `sprat-record-1.0`; importer for the deposited 2026 records with a documented key map;
- analysis layer: sensitivity, Fano fits, the series of the systematic runs, leak-corrected barrier analysis with Monte Carlo error budget and Fabry-Perot fits, field-map bounds and row decay, grading of the two frozen prediction sets, numbers registry, audit against the deposited registry and the numbers registry of the manuscript;
- figures 1 to 7 of the paper, in its order, drawn from the records by the package `sprat/figures/` (restored: a pre-release archive lacked it); deposit builder (`--include` for extra files such as the supplementary document); the task list of the data set for regeneration (`campaigns/manuscript`);
- `sprat run-one` re-runs a record from itself (`--task <record.json> --out <new.json>`);
- `.gitattributes` keeps every checkout at LF so that the checksum tests pass on Windows as well;
- `.gitignore` anchors the build outputs (`/records/`, `/tables/`, `/figures/`, `/logs/`, `/deposit*/`, `/build/`, `/dist/`, the task lists) to the repository root, so that the package folder `sprat/figures/` is committed and archived;
- the documentation states where the deposited records come from: the original scripts of 2026 on the UHeM Altay cluster (`02_kavite.py` with the v4 patch, Meep 1.30.0), converted with `sprat import-legacy`; SPRAT converts and analyses them and can regenerate each record with Meep;
- repository URL `https://github.com/codekyha/sprat` in `pyproject.toml`, `CITATION.cff`, `codemeta.json` and the documentation;
- `sprat deposit --data-doi --software-doi`: the DOIs enter the data availability statement, the README and the metadata of the pack; the metadata carries the title of the data record without a version, both licences and only identifiers that are given; the numbers registry of the manuscript is copied into the pack as `analysis/manuscript_registry.json`; the builder warns about absolute paths in the text files of the pack;
- no absolute paths of the build machine in the outputs: the analysis outputs use relative paths, and the import manifest stores the name of the deposit directory only;
- manually started CI workflow `regression`: Meep 1.30 from `environment.yml`, `sprat check` and `sprat selftest --reference`;
- the audit covers every record-derived row of the numbers registry of the manuscript (`numbers_registry_manuscript_v8.json`: 107 rows, 9 superseded statements, literature values and references listed apart); the 2.0.2 registry still gives 83;
- `sprat deposit` packs `analysis/` and `predictions/` as `analysis.tar.gz` and `predictions.tar.gz`, because Zenodo stores files without folders;
- figure 4 draws the error bar of the second defect radius from the computed uncertainty of the analysis;
- one Zenodo record for the software and the data: `sprat deposit --software-doi` equal to `--data-doi` with `--include sprat-<version>.zip` writes the statement, README, licence and metadata for SPRAT archived in the data record (no isCompiledBy relation; with such an archive included and no software DOI given, the data DOI stands for both); `CITATION.cff` carries that DOI as the identifier of SPRAT 1.0.0 and names the record as the preferred citation.

Verified with Meep 1.30.0 built from source: `sprat selftest` and `sprat selftest --reference` pass, and `sprat run-one` reproduces archived records to about 1e-11 in Q.
