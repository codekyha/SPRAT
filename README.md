# SPRAT: Side-coupled Photonic-crystal Resonator Analysis Toolkit

SPRAT is a standalone, fully English software pack for the two-dimensional finite-difference
time-domain (FDTD), plane-wave and band-structure study of a point-defect microcavity
side-coupled to a W1 waveguide in a square lattice of dielectric rods immersed in a liquid
analyte. It re-implements the workflow of the computational runs behind the paper "Discrete
quality-factor control in a side-coupled photonic crystal microcavity: evanescent Bloch
tunnelling and the finite-cell correction" (H. Oguz; preprint arXiv:2609.31952). The records of
that paper's data set were produced on the UHeM Altay cluster by the original scripts of 2026
(`02_kavite.py` with the v4 patch, Meep 1.30.0) and converted into the SPRAT schema with
`sprat import-legacy`; SPRAT converts and analyses them and can regenerate each record with Meep.

* the photonic-crystal **structure is read from a text file** (`.phc`);
* the **run parameters are read from a second text file** (`.par`), including sweeps;
* every run writes **one self-describing JSON record** (English keys, schema `sprat-record-1.0`);
* a **workstation runner** executes hundreds of independent Meep runs in parallel, pinned to
  cores, with resume, retry, a job log and an ETA;
* the **analysis layer** turns records into the numbers of the paper (sensitivity, Fano fits,
  the leak-corrected barrier series and its decay constant, the physical channels of the complex
  band structure, the analyte energy fraction with the permittivity $E_z$ sees, the comparison of
  the transmission spectra under the two guide terminations and the cell factor across the analyte
  sweep, the reflectionless quality factor and figure of merit, the grading of the criteria fixed
  before the runs, a numbers registry, an independent audit);
* the **figures** of the paper are drawn from the records;
* **any record can be regenerated** with Meep from the structure and parameters stored in it,
  and so can the whole data set (`campaigns/manuscript/`); `sprat deposit` packs a records
  directory for Zenodo;
* the deposited 2026 data set (Turkish keys) is **imported** by a documented key map, and the
  audit recomputes the registered numbers of the paper (the numbers registry) from the records
  and compares them with the deposit;
* SPRAT was written by the author with coding assistance from Claude (Anthropic; Claude Opus 5.5
  and Claude Fable 5) and Gemini (Google; Gemini 3.6 Flash and Gemini 3.1 Pro), used in tandem
  (section 12).

Units: the lattice constant $a$ and the vacuum speed of light are 1; frequencies are
$f = a/\lambda$.

---

## 1. Model

Square lattice of silicon rods (radius $r = 0.20a$, $\varepsilon = 11.9025$) in an analyte of
index $n_a$ (1.33 for water). One row of rods is removed (the W1 guide). A point defect of
radius $r_d$ sits $N_{\rm sep}$ rows from the guide, optionally displaced by $(\delta_x, \delta_y)$.
TM polarisation ($E_z$). The lattice constant for the 1550 nm target is $a = f_{\rm mid}\,\lambda_{\rm target} = 481.4$ nm,
with $f_{\rm mid}$ the centre of the TM gap.

FDTD (Meep): resolution 24 points per $a$ in the reference geometry, Courant 0.5, subpixel
smoothing, PML of one period on every side (with the `pml` termination) or an adiabatic
absorber along the guide (`absorber` termination, into which the lattice and the guide
continue). Mirror symmetry in $x$ is used when the defect is on the axis. The quality factor
comes from harmonic inversion of $E_z$ at a probe point next to the defect after the source is
off. With $Q_{\rm lim} = \pi f_{\rm cen}\, t$ (the $Q$ of a mode that decays by a factor $e$ over the
signal length $t$) the margin of a record is $Q_{\rm lim}/Q$. Harmonic inversion resolves decays
longer than the signal, so the rule margin $\ge 1$ is a conservative sampling rule. A record is
admitted when its margin is at least 1 and its resonance lies between 0.2 and 0.8 of the TM gap;
fits and statistics use admitted records only, and the figures show the excluded records of a
sweep as open grey symbols. The criterion tests the numerics alone: in the preliminary sweep
($n_{\rm cl} = 6$, resolution 20) the strongest mode of the 537 admitted records at $N_{\rm sep} = 1$
or $r_d \ge 0.13a$ is a low-$Q$ Fabry-Perot mode of the finite waveguide ($Q$ 50 to 77, median 53);
436 of them lie in the dead zone $r_d \ge 0.14a$, where this holds at every radius. The column
`waveguide_mode` of `records.csv` flags them, they enter no fit, and 204 admitted records hold a
cavity mode (`valid` true and `waveguide_mode` false).

Run modes: `harminv` (complex eigenfrequency: $f_r$, $Q$), `spectrum` (transmission), `reference`
(cavity-less transmission, the normalisation), `field` (DFT map of $E_z$ at $f_{\rm res}$ and the
analyte energy fraction), `bands` (MPB: bulk bands, the gap against the analyte index, the W1
projection, the lattice-constant calibration) and `pwe` (plane-wave expansion: bulk gap, W1
band, complex band structure and the decay channels, Fabry-Perot reading, barrier-row radius,
loss budget).

---

## 2. Installation

### Linux and macOS (Meep from conda-forge)

```bash
conda env create -f environment.yml          # python 3.11, pymeep 1.30 (serial), numpy, scipy, matplotlib, pandas, psutil
conda activate sprat
pip install -e .
sprat check
```

`pymeep=1.30.*` is pinned because the deposited records were produced with Meep 1.30 (by the
original scripts of 2026); a newer Meep runs too (`sprat check` reports the version) and the
regression test tells you whether the numbers moved.

### Windows

Meep's Conda packages do not run on native Windows. Install the **Windows Subsystem for Linux**
(WSL2, Ubuntu), then follow the Linux steps inside it. WSL2 exposes every physical core and the
runner pins tasks to cores there as on Linux. Native Windows can still run everything that does
not need Meep: `pip install -e .[figures,tables]` gives the plane-wave layer, the analysis, the
figures and the deposit builder.

### Analysis only (no Meep)

```bash
pip install -e .[figures,tables]
```

---

## 3. Quick start

```bash
sprat check                                                        # environment, cores, memory, Meep, calibration
sprat calibrate --records records                                  # measures pixel-steps per second on this machine (about 2 minutes)
sprat plan examples/structures/reference_w1_notch.phc examples/params/sweep_defect_radius.par -o tasks --records records
sprat run tasks.jsonl --workers 8                                  # 13 harmonic inversions in parallel
sprat status tasks.jsonl
sprat collect records -o tables                                    # flat table with the validity flags
sprat analyze records -o tables                                    # the derived numbers, the numbers registry, the prediction grades
sprat figures records -o figures --tables tables                   # figures 1 to 7
```

`sprat plan` prints the cost table (tasks, pixel-steps, single-core hours, the wall time for
several worker counts, memory per task) before anything runs; `sprat run --dry-run` prints it
again for the actual task list.

---

## 4. The structure file (`.phc`)

Plain text: `[section]` headers, `key = value`, `#` comments. Every key has a default; unknown
keys are errors. Units are $a$ unless stated. Row indices are counted from the guide row
($j = 0$); positive rows lie on the cavity side. Full reference: `docs/STRUCTURE_FILE.md`.

```ini
[lattice]
type = square
a_nm = 481.4                # or "calibrate": f_mid * target_wavelength_nm from the newest bulk bands record
target_wavelength_nm = 1550
rod_radius = 0.20
rod_eps = 11.9025

[cell]
guide_periods = 25          # n_x
cladding_rows = 12          # n_cl per side
pad_x = 1.5                 # pml termination only
pad_y = 1.0
termination_x = pml         # pml | absorber
absorber_periods = 8        # absorber termination only

[waveguide]
type = W1                   # W1 | none

[defect]
type = rod                  # rod | none
radius = 0.060
row = 4                     # N_sep
dx = 0.0
dy = 0.0

[rows]                      # per-row radius overrides
2 = 0.18

[rods]                      # free-form edits, applied last
add 3.0 -2.0 0.20
remove 0.0 4.0
```

Cell size: $s_x = n_x + 2\,{\rm pad}_x + 2\,d_{\rm PML}$ (pml) or $n_x + 2\,n_{\rm abs}$ (absorber);
$s_y = 2 n_{\rm cl} + 1 + 2\,{\rm pad}_y + 2\,d_{\rm PML}$. The builder was tested against the rod
lists of the original script for nine geometries (`tests/test_geometry.py`).

---

## 5. The parameter file (`.par`)

Same grammar. Full reference: `docs/PARAMETER_FILE.md`.

```ini
[run]
mode = harminv              # harminv | spectrum | reference | field | bands | pwe
label = auto                # record name; auto builds it from the resolved values
tag = sweep_fine            # stored in every record
q_est = 8000                # expected Q: sets the harminv time (with margin) and the stopping ceilings

[analyte]
n = 1.33
k = 0.0

[numerics]
resolution = 24
courant = 0.5
pml = 1.0
subpixel = true
symmetry = auto             # auto | none | even | odd

[source]
fcen = auto                 # a_nm / target_wavelength_nm
fwidth = 0.06
cutoff = 5

[harminv]
t = auto                    # margin * q_est / (pi fcen), rounded up to 100
margin = 4.5
auto_t = false              # the two-pass rule of the original scripts

[spectrum]
nfreq = 201

[field]
f_res = 0.304614            # required for mode = field

[sweep]                     # any key of either file by its dotted path; Cartesian product unless zipped
structure.defect.radius = 0.050 : 0.110 : 0.005
analyte.n = 1.30, 1.33, 1.45
zip = structure.defect.row, structure.cell.cladding_rows
```

The harmonic-inversion rule of the verification runs is

$$t = 100\Big\lceil\frac{m\,Q_{\rm est}}{100\,\pi f_{\rm cen}}\Big\rceil,\qquad Q_{\rm lim} = \pi f_{\rm cen}\,t,\qquad \text{margin} = Q_{\rm lim}/Q .$$

Command-line overrides: `sprat plan s.phc p.par --set structure.defect.row=5 --set harminv.margin=2.5`.
Every record stores the fully resolved structure and parameters, so any record can be re-run
from itself (`sprat run-one --task records/<label>.json --out rerun/<label>.json`) and two
records describe the same task exactly when these blocks agree.

---

## 6. Records

One JSON per run, `records/<label>.json` (a field run also writes `<label>.npz`):

```
schema, software {name, version, meep, python, code_sha256},
task {label, tag, created, structure_file, params_file, overrides},
structure {...}, params {...},
provenance {host, started, wall_s, steps, pixels, pixel_steps, throughput, symmetry, slot, job_id},
result {harminv: modes[{f, Q, amplitude, error, wavelength_nm, fwhm_nm}], t_used, auto_history, Q_limit
        spectrum/reference: f, flux_in, flux_out, t_end
        field: analyte_energy_fraction, wavelength_nm, S_first_order_nm_per_RIU, field_file
        bands: task, gap, freqs, ..., a_nm      pwe: the plane-wave results}
```

`sprat collect` flattens a directory into `records.csv` / `records.jsonl` with the validity
flags of the paper: `resolved` (margin $\ge 1$), `gap_position` (0 to 1 inside the TM gap at the
record's analyte index, from the MPB gap-edge table) and `valid` (resolved and 0.2 to 0.8). The
full schema is in `docs/RECORD_SCHEMA.md`.

---

## 7. Parallel execution on a workstation

`sprat run tasks.jsonl [--workers N] [--cores-per-task 1] [--no-pin] [--order longest|shortest|file] [--overwrite] [--retry 1] [--dry-run]`

* every task runs in its own interpreter, so Meep's memory is released after each task and a
  crash cannot take the batch down;
* `--workers` defaults to the physical cores (not the logical threads: SMT gives Meep nothing);
* on Linux and WSL2 each slot is pinned to its own cores; `OMP_NUM_THREADS=1` for every task;
* tasks start longest-first (the makespan is then the longest task or the total divided by the
  workers, whichever is larger);
* a task whose record exists is skipped; a record that exists for *different* parameters under
  the same label is reported as a collision and not overwritten;
* a failed task is retried once; `logs/<task list>_<time>/joblog.tsv` records every attempt
  (return code, wall time, slot, cores, throughput) and `task_<id>_<attempt>.out` the output;
* the slot count is capped so that the estimated memory of the running tasks fits the available
  memory; the ETA uses the cost model corrected by the observed speed;
* Ctrl-C stops launching new tasks and lets the running ones finish; a second Ctrl-C terminates them.

`sprat status tasks.jsonl` prints done / pending / failed and the remaining single-core hours.

### The cost model

With Courant factor $C$, $\Delta t = C/r$, so a run of $T$ time units on an $s_x \times s_y$
cell at resolution $r$ costs

$$W = (s_x r)(s_y r)\,\frac{T}{\Delta t} = \frac{s_x s_y T r^3}{C}\ \text{pixel-steps},$$

halved with mirror symmetry. The reference computation of the paper ($30a\times29a$, $r = 24$,
$t = 17\,750$) is $2.16\times10^{11}$ pixel-steps and took 677 s on one core of an AMD EPYC
7742 node running 32 such tasks ($3.19\times10^{8}$ pixel-steps/s). `sprat calibrate` runs the
reference geometry at resolutions 12 and 24 with a short signal and writes `calibration.json`;
`plan` and `run` use it. Small cells run at a lower pixel-step rate than the model assumes
(the coarse sweep of the systematic runs ran ten times slower per pixel-step than the production
cells), so treat the estimate as an upper bound on efficiency, not a promise.

---

## 8. Reproducing the paper

Everything below is `sprat analyze records -o tables` and `sprat figures records -o figures`
on a records directory that holds the data set (imported or regenerated). Figures and tables are
numbered as in the paper.

| manuscript item | source (`tables/analysis.json` unless stated) | figure |
|---|---|---|
| geometry of the computational cell | the reference record (its structure block) | 1 |
| TM gap against the analyte index | the band-structure records (`bands`, tasks `bulk` and `sweep`) | 2 |
| table 1, figure 3(a, b) | `base_series` ($n_{\rm cl} = 8$, res 20, $N_{\rm sep} = 2, 3, 4$; quadratic law $f_r(r_d)$), `radius_sweep` (the $0.005a$-step fine sweep) | 3 |
| figure 3(c) | `kappa_series` (raw PML series), `resolved_pair` ($r_d = 0.100a$) | 3 |
| table 2, table 3, figure 3(d) | `convergence` (resolution, table 2), `cladding_ladder` (cladding periods, table 3 and figure 3(d)) | 3 |
| barrier mechanism and finite-cell correction | `barrier` (the three guide terminations, the cell factor $G$ with the two-reflector fit, the channel identification of $\kappa$ against frequency with the plane-wave record, the cladding-leak linearisation; also the three-point solve, the four-point series and $\kappa$, Monte Carlo, second radius, clearance rule, performance in water, tolerance) | 4 |
| table 4, figure 5 | `displacement` | 5 |
| transmission spectrum and linewidth check | `spectra_summary` and `tables/spectra.json` (Fano fits, both normalisations, coupling branches; the spectra of the 25a cell); the inset of figure 6(a) also reads `tables/spectra_absorber.json` (the spectra with the guide continued into the absorber) | 6 |
| the spectra with the guide continued into the absorber: the asymmetry $q$ under both terminations, the grades by the criterion fixed before the runs and by its amendment, the cell factor $G$ and the reflectionless $Q_w$ across the analyte sweep | `termination` | - |
| table 5, sensitivity, FOM, detection limit | `analyte_sweep` | 7 |
| the reflectionless $Q_w$ and ${\rm FOM}_w$ of table 5 and figure 7(c, d) | `combined.physical_values` (the rule, the values at the rows of table 5 and along the sweep, the check at $n_a$ = 1.33) | 7 |
| convergence of the resonance frequency with resolution | `convergence.resolution_20_vs_32_frequency`, `barrier.resolution_check` | - |
| analyte energy fraction (Hellmann-Feynman, with the stored-map estimators of the original analysis) | `fields.eta`, `fields.reference` | - |
| row decay of the stored field | `fields.row_decay` | - |
| criteria fixed before the runs | `predictions` and `tables/predictions_report.md` | - |
| the registered numbers of the paper | `tables/numbers_registry.json` and `.md` | - |

The comparison with published sensors (section S11 of the supplementary material) is not
computed from the records.

**Complex band structure.** The plane-wave layer (`sprat pwe`) linearises the quadratic
eigenproblem in $k_y$ at fixed $f$ and $k_x$ into a companion matrix. In a truncated basis part of
its spectrum is spurious; a root is kept as a Bloch channel when less than $10^{-3}$ of its
eigenvector weight lies on the two outermost rings of $G_y$ and one of its replicas has
$|{\rm Re}\,k_y| \le \pi/a$ (within $5\times10^{-3}$). The slowest physical channel at the
guided-mode wavevector lies at the zone edge; the slower root that versions up to 1.0.0 reported as
the least-evanescent pair is an artefact of the square basis (it moves when the basis becomes
circular) and is reported as `boundary_artefact`. A records directory imported from the 2026
deposit carries only the legacy plane-wave result; run `sprat pwe` (or copy the plane-wave record
of the data record, `pwe_records.tar.gz`) before `sprat analyze` to obtain the physical channels.

**Analyte energy fraction.** Under Meep's subpixel smoothing $E_z$ sees the arithmetic mean
$\langle\varepsilon\rangle$ over a voxel, and the Hellmann-Feynman theorem for the discretised
operator weights each point by its analyte fraction. The permittivity map Meep stores
(`get_array(mp.Dielectric)`) is the harmonic mean of the eigenvalues of the smoothed tensor;
`fields.meep_smoothing` reconstructs the tensor from the geometry, checks it against the stored map
and evaluates the consistent estimate. The estimators that read the stored map as
$\langle\varepsilon\rangle$ are kept for comparison and are biased upwards.

**Selection policy.** The analysis never refers to records by name: it selects them by their
parameters. When a point was run more than once, the record with the longest
harmonic-inversion signal (highest $Q_{\rm lim}$) is used, except for figures 2, 3 and 5 to 7,
which first exclude the verification reruns (`tag = verification`, or the imported v4 set) so
that the systematic runs are reproduced as run. One consequence for figure 3(c): the
$N_{\rm sep} = 5$ point is the 120 000-time-unit rerun of the systematic runs (margin 2.54)
instead of the 60 000 run with margin 1.27 that section 2.4 of the paper quotes; the two agree
to 0.003 % in $Q$. The spectrum analysis of section 3.5 (`spectra_summary`) reads the spectra of
the 25a cell; the spectra with the guide continued into the absorber enter only `termination`.

**Barrier model.** For a defect $N_{\rm sep}$ rows from the guide with clearance $c = n_{\rm cl} - N_{\rm sep}$,

$$\frac{1}{Q} = \frac{e^{-\kappa N_{\rm sep}}}{A\,G_T} + \frac{1}{Q_{\rm top}(c)},\qquad Q_{\rm top}(c) = Q_{\rm top}(8)\,e^{\kappa_{\rm top}(c-8)},$$

with $G_{\rm PML}$ the Fabry-Perot cell factor of the PML-terminated guide and $G_{\rm abs} = 1$.
The three PML points at clearance 8 determine $\kappa$, $Q_{\rm top}(8)$ and $A$ exactly; the
$n_{\rm cl} = 12/13$ pair at $N_{\rm sep} = 5$ gives $\kappa_{\rm top}$; every record is
leak-corrected, the PML records are divided by the measured $G = Q_{\rm PML}/Q_{\rm abs}$, and the
four-point series $Q_w(N_{\rm sep})$ gives $\kappa$ with a Monte Carlo error budget.

**Guide termination.** The seven spectra of the analyte sweep ($n_a$ = 1.300 to 1.450) were
repeated on 24 September 2026 with the guide continued into an $8a$ adiabatic absorber, every
other run parameter unchanged. `termination` fits both sets with the routine of the grading
script (the Fano form $T = T_{\rm bg} + T_0(q+\epsilon)^2/(1+\epsilon^2)$,
$\epsilon = 2(\lambda-\lambda_r)/{\rm FWHM}$, over the whole window, each spectrum normalised by the
cavity-less reference run of its own termination) and applies the gates G1 to G5 and S1. It
grades the factor $|q|_{\rm cell}/|q|_{\rm absorber}$ by the criterion fixed before the runs
(D-22: $\ge 3$ supports the attribution of the residual asymmetry to the guide ends, $\le 1.5$
refutes it, a pass needs support at all seven indices) and by its amendment A1 (fixed before the
first absorber cavity spectrum was written: the ratio is graded at the five indices where
$|q|_{\rm cell} > 0.2$), and reports both verdicts. Amendment A2 gives the cell factor
$G(n_a) = Q_{\rm cell}/Q_{\rm absorber}$ from the two fits and the reflectionless
$Q_w = Q_{\rm cell}^{\rm harminv}/G$. The frozen files are in `predictions/frozen/`; the job files,
the grading scripts and the spectra are in `raw_h14b.tar.gz` of the data record.

### The deposited data set

The data record of the paper (https://doi.org/10.5281/zenodo.22912910, version 3.1.0) holds the
records in the SPRAT schema and the originals (Turkish keys) with the scripts that wrote them:
`raw_legacy.tar.gz`, the deposit of 2026 (version 2.2.0), byte for byte apart from the change
listed in its `CHANGES_FROM_2.2.0.md`, and `raw_h14b.tar.gz`, the 14 runs with the guide continued
into the absorber (job 528207) with the job files, the criterion fixed before the job, its
amendment and the grading scripts. The records were produced on the UHeM Altay cluster by the
original scripts (`02_kavite.py` with the v4 patch, Meep 1.30.0) and converted into the SPRAT
schema with `sprat import-legacy`.

The converted records are read directly from the archives of the data record:

```bash
mkdir records
for f in harminv_records spectra_records field_records bands_records pwe_records; do
    tar -xzf $f.tar.gz -C records --strip-components=1
done
sprat analyze records -o tables
sprat audit records --expected campaigns/manuscript/expected/numbers_registry_manuscript_v10_2.json -o tables
```

The originals can be converted again:

```bash
tar -xzf raw_legacy.tar.gz                              # the 2026 deposit (version 2.2.0), Turkish keys
tar -xzf raw_h14b.tar.gz                                # the absorber-terminated spectra of 24 September 2026
sprat import-legacy raw_legacy -o records               # also takes the packed or unpacked earlier packs
sprat import-legacy raw_h14b -o records
tar -xzf pwe_records.tar.gz -C records --strip-components=1   # the plane-wave record with the physical channels
sprat audit records --expected campaigns/manuscript/expected/numbers_registry_manuscript_v10_2.json -o tables
```

The importer converts the Turkish-keyed files with the map in `docs/LEGACY_KEYS.md`, never
modifies the originals, and stores the original file name and SHA-256 in every converted record.
The audit recomputes the registered numbers of the paper (the numbers registry) from the records
and compares them with the deposit row by row, each number to the precision the registry prints
and counts exactly. `numbers_registry_manuscript_v10_2.json` is the registry of the manuscript,
version 10.2: the registry of version 10 with the section `derived_v10_2` (the reflectionless
quality factor and figure of merit at the rows of table 5, the convergence of the resonance
frequency with resolution, the guided-mode wavevector and the per-row factor at $k_x = 0$, the
absorption-limited quality factor with the mixed points of the grid counted either way), the
three references added in version 10.2 and the literature value taken from one of them. Its
literature values and references do not come from the records and are listed apart; the audit
of the full record gives OK on all 237 record-derived rows and CONSISTENT on all 21 superseded
statements. `numbers_registry_manuscript_v10.json` is the registry of version 10: the registry
of version 9 with the section `termination` (the spectra with the guide continued into the
absorber, both grades and the cell factor across the analyte sweep), the section
`corrected_v10` and derivation texts brought up to date (OK on all 223 record-derived rows).
The registries of versions 9 (`numbers_registry_manuscript_v9.json`), 7
(`numbers_registry_manuscript_v7.json`, the `numbers_registry.json` of the 2.2.0 deposit byte for
byte) and 8 (`numbers_registry_manuscript_v8.json`) and of the unpublished version 2.0.2 of the
data set (`numbers_registry_expected.json`) can be given instead. On the records of the 2.2.0
deposit, the registry of version 9 then differs only in the two quantities that SPRAT 1.2.0
renames ("records from the systematic runs", "eta_a, threshold mask of the original analysis"),
and the earlier registries in exactly the rows that version 9 corrects (all points instead of
admitted ones, the basis artefact instead of the physical channel, the stored map instead of the
permittivity $E_z$ sees) and the first of those names.

### Regenerating the data set with Meep

```bash
sprat calibrate --records records
sprat reproduce campaigns/manuscript --tier quick --records records --dry-run     # cost table only
sprat reproduce campaigns/manuscript --tier full  --records records --workers 16
sprat audit records --expected campaigns/manuscript/expected/numbers_registry_manuscript_v10_2.json -o tables
```

`campaigns/manuscript/README.md` gives the tiers and the cost: 892 tasks; about 87 single-core
hours by the cost model at the calibrated throughput, about 256 core-hours as measured on the
cluster in 2026; the longest task (13.6 h on one core) sets the floor of the wall time.

### Regression

`sprat selftest` runs the four FDTD modes and the bands at resolution 8 in a few minutes;
`sprat selftest --reference` adds the regression of the reference record
($Q = 6926.72$, $f_r = 0.304614$; tolerance $|Q/Q_{\rm ref} - 1| < 10^{-3}$, $|f_r - f_{\rm ref}| < 10^{-5}$).
With Meep 1.30.0 and MPB 1.11.1 built from source on Ubuntu 24.04, `sprat check` reports the
environment complete, `sprat selftest` passes, and `sprat selftest --reference` gives
$Q = 6927.27$ ($|\Delta Q/Q| = 7.9\times10^{-5}$) and $f_r = 0.304614$
($|\Delta f_r| = 8\times10^{-8}$) in 855 s; the small offset in $Q$ comes from the signal length
of the built-in task, $t = 17\,750$, against 17 751.78 in the archived record. `sprat run-one`
repeats archived records: a smoke run with a $6a$ absorber agrees in $Q$ to $6\times10^{-12}$,
the absorber-terminated point at $N_{\rm sep} = 3$ ($Q = 634.03$) in $Q$ to $3.5\times10^{-12}$
and in $f_r$ to $3\times10^{-15}$, and the cavity-less reference run at $n_a = 1.425$ in its flux
spectra to $1.3\times10^{-14}$.
On GitHub the workflow `regression` (started by hand: Actions, regression, Run workflow) installs
Meep 1.30 from `environment.yml` and runs `sprat check` and `sprat selftest --reference`.

---

## 9. Depositing

```bash
sprat deposit records -o deposit_v3 --tables tables --raw-legacy raw_legacy --raw-h14b raw_h14b --version 3.1.0 \
      --include supplementary.pdf sprat-1.2.1.zip \
      --data-doi 10.5281/zenodo.22912910 --software-doi 10.5281/zenodo.22912910 --arxiv-id 2609.31952
```

builds the data pack: flat tables, tarballs of the records by mode, the analysis outputs with
`analysis/manuscript_registry.json` (the numbers registry of the manuscript, copied from
`campaigns/manuscript/expected/`) as `analysis.tar.gz`, the frozen predictions as
`predictions.tar.gz` (Zenodo stores files without folders), `parameters_manifest.json`,
`README.md` with the origin of the records and the data dictionary, `DATA_AVAILABILITY.md` (the
statement for the paper, then the counts), `LICENSE.txt`, `zenodo_metadata.json` and
`MANIFEST_sha256.txt`. `--raw-legacy` packs a directory of original files byte for byte as
`raw_legacy.tar.gz` (for the paper, the unpublished version 2.2.0 of the deposit; a file
`CHANGES_FROM_*.md` in it is named in the texts of the pack), and `--raw-h14b` packs the directory
of the absorber-terminated spectra with their job files, criterion and grading as
`raw_h14b.tar.gz`; `--include` places extra files such as the supplementary document in the pack
before the manifest is written. `--arxiv-id` adds the preprint the record supplements to the
README and, as the related identifier *Is supplement to*, to the metadata.
`--data-doi` and `--software-doi` enter the statement, the README and the metadata; without them
the texts carry the placeholders `10.5281/zenodo.NNNNNNN` and `10.5281/zenodo.SSSSSSS`, and the
metadata carries no identifier. Run the command from the repository checkout: the frozen
predictions and the task list of the data set live in the repository, not in the installed
package. The builder warns about text files of the pack that contain absolute paths of the build
machine.

For the paper, the software and the data share one Zenodo record: `--software-doi` equal to
`--data-doi` (or left out) and `--include sprat-1.2.1.zip` put the source archive of the release
into the pack, and the statement, the README, the licence file and the metadata then describe
SPRAT as part of the record (no *Is compiled by* relation). A separate software record would come from Zenodo's
GitHub integration (`.zenodo.json` holds its metadata); it is not used for the paper, so the
integration stays switched off for this repository.

### Release checklist

1. The record 10.5281/zenodo.22912910 holds the software and the data; its version 3.1.0 is
   published (September 2026).
2. Push the repository to GitHub; leave Zenodo's GitHub integration switched off, otherwise a
   release mints a second DOI. A GitHub release `v1.2.1` can be published as usual.
3. The preprint is listed as arXiv:2609.31952 (physics.optics, September 2026).
4. Build the release archive `sprat-1.2.1.zip` from the files Git commits, and the data pack with
   `--version 3.1.0 --data-doi 10.5281/zenodo.22912910 --software-doi 10.5281/zenodo.22912910
   --raw-legacy raw_legacy --raw-h14b raw_h14b --include supplementary.pdf sprat-1.2.1.zip`, with
   `--arxiv-id 2609.31952`. The files of the published record 3.1.0 were built before the identifier
   existed; the record names the preprint as a related work (*Is supplement to*) in its metadata.
5. Upload the files of the pack into the draft and publish it. A later version of the pack goes
   into a new version of the record: Zenodo then mints a new DOI for it, and a manuscript cites
   that DOI or the concept DOI. The journal DOI (relation *Is supplement to*) is added to the
   record's metadata later.

---

## 10. Repository layout

```
sprat/            the package: _textfile, structure, params, geometry, cost, plan, runner, execute, fdtd, bands, pwe,
                  records, legacy, doctor, reproduce, deposit, cli; analysis/ (collect, fano, sensitivity, barrier,
                  fields, termination, combined, predictions, registry, audit, pipeline); figures/ (style, inputs,
                  schematic, campaign, mechanism, build)
examples/         structure and parameter files
campaigns/        manuscript/: the task list of the data set, the stage files, the expected registries
predictions/      the criteria fixed before the runs (three sets and one amendment), byte for byte
docs/             STRUCTURE_FILE, PARAMETER_FILE, RECORD_SCHEMA, LEGACY_KEYS, COST_MODEL, REPRODUCING_THE_PAPER, DATA_AVAILABILITY
tests/            pytest (parsers, geometry equivalence, cost, runner, legacy import, analysis regression, guide terminations,
                  manuscript audit, figures, predictions)
```

## 11. Tests

```bash
pytest                                   # without Meep: parsers, geometry, planner, runner (synthetic engine), legacy sample, plane-wave layer
SPRAT_LEGACY_DEPOSIT=/path/to/raw_legacy pytest      # adds the full-deposit round trip, the analysis regression and the figures
SPRAT_LEGACY_DEPOSIT=/path/to/raw_legacy SPRAT_H14B_DEPOSIT=/path/to/raw_h14b pytest   # adds the absorber spectra, the reflectionless values, the frequency convergence and the audits of versions 10 and 10.2
```

`SPRAT_LEGACY_DEPOSIT` points at the 2026 deposit: the unpacked `raw_legacy.tar.gz` of the data
record (or `tests/data/zenodo_pack_v2`); `SPRAT_H14B_DEPOSIT` at the unpacked `raw_h14b.tar.gz`.
Without them the comparison of the two terminations is tested at $n_a = 1.300$ on the fixtures. The runner tests use `SPRAT_FAKE_FDTD=1`, a synthetic
engine that writes records with the cost model's steps and pixels; the Meep smoke test runs only
where Meep is installed (`sprat selftest`).

## 12. Citing

SPRAT accompanies the paper "Discrete quality-factor control in a side-coupled photonic crystal
microcavity: evanescent Bloch tunnelling and the finite-cell correction" (H. Oguz, 2026), available
as the preprint arXiv:2609.31952 (https://arxiv.org/abs/2609.31952). SPRAT 1.2.1 is archived
together with the data set of the paper in one Zenodo record, https://doi.org/10.5281/zenodo.22912910
(version 3.1.0; the file `sprat-1.2.1.zip`). Please cite the paper and the record. `CITATION.cff`
carries both, and GitHub's "Cite this repository" button reads it. The journal version, once
published, replaces the preprint in this citation.

```bibtex
@misc{oguz2026discretequalityfactorcontrolsidecoupled,
  title         = {Discrete quality-factor control in a side-coupled photonic crystal microcavity: evanescent Bloch tunnelling and the finite-cell correction},
  author        = {Hasan Oguz},
  year          = {2026},
  eprint        = {2609.31952},
  archivePrefix = {arXiv},
  primaryClass  = {physics.optics},
  url           = {https://arxiv.org/abs/2609.31952}
}

@misc{oguz2026sprat_record,
  title     = {Data, code and supplementary material for: Discrete quality-factor control in a side-coupled photonic crystal microcavity: evanescent Bloch tunnelling and the finite-cell correction},
  author    = {Hasan Oguz},
  year      = {2026},
  publisher = {Zenodo},
  version   = {3.1.0},
  doi       = {10.5281/zenodo.22912910},
  url       = {https://doi.org/10.5281/zenodo.22912910},
  note      = {Zenodo record; holds the data set, the supplementary document and SPRAT 1.2.1 (sprat-1.2.1.zip)}
}
```

Funding of the original computations: Istanbul Okan University BAP OBAP2026010006, Pamukkale
University BAP 2025ALDEP037, UHeM grant 5027772026 (UHeM Altay cluster). Licence: MIT (software);
the data set is CC BY 4.0, and the scripts inside its `raw_legacy.tar.gz` and `raw_h14b.tar.gz`
are MIT.

SPRAT was written by the author with coding assistance from Claude (Anthropic; Claude Opus 5.5
and Claude Fable 5) and Gemini (Google; Gemini 3.6 Flash and Gemini 3.1 Pro), used in tandem under
the author's supervision; the author is responsible for its design, its physics and its
verification.

This software and dataset are provided "as is" without warranty of any kind, express or implied, or commitment to ongoing maintenance or support.
