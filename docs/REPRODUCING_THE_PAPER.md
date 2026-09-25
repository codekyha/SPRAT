# Reproducing the paper

Three routes, in increasing cost. The records of the paper's data set were produced on the UHeM
Altay cluster by the original scripts of 2026 (`02_kavite.py` with the v4 patch, Meep 1.30.0) and
converted into the SPRAT schema with `sprat import-legacy`. The data record
(https://doi.org/10.5281/zenodo.22912910) holds the converted records and the originals with the
scripts that wrote them: `raw_legacy.tar.gz` (the deposit of 2026) and `raw_h14b.tar.gz` (the
spectra with the guide continued into the absorber, with the job files, the criterion fixed
before the job, its amendment and the grading scripts).

## 1. From the deposited records (no Meep, minutes)

Either the converted records of the data record:

```bash
mkdir records
for f in harminv_records spectra_records field_records bands_records pwe_records; do
    tar -xzf $f.tar.gz -C records --strip-components=1
done
```

or the originals, converted again:

```bash
tar -xzf raw_legacy.tar.gz && tar -xzf raw_h14b.tar.gz
sprat import-legacy raw_legacy -o records
sprat import-legacy raw_h14b -o records
tar -xzf pwe_records.tar.gz -C records --strip-components=1   # the plane-wave record of SPRAT 1.1.0
```

then

```bash
sprat analyze records -o tables
sprat figures records -o figures --tables tables
sprat audit records --expected campaigns/manuscript/expected/numbers_registry_manuscript_v10_2.json -o tables
```

`tables/ANALYSIS_REPORT.md` summarises the headline numbers, `tables/numbers_registry.md` lists
the registered numbers with their record or derivation, `tables/predictions_report.md` grades the
frozen predictions, and `tables/AUDIT_REPORT.md` compares the numbers recomputed from the records
with the registry of the manuscript, row by row (its literature values and references are listed
apart, since they do not come from the records).

## 2. Regenerating the records (Meep, a workstation)

```bash
sprat calibrate --records records
sprat reproduce campaigns/manuscript --tier full --records records --workers 16
sprat analyze records -o tables && sprat figures records -o figures --tables tables
sprat audit records --expected campaigns/manuscript/expected/numbers_registry_manuscript_v10_2.json -o tables
```

See `campaigns/manuscript/README.md` for the tiers and the cost. A single record is run again
from itself with `sprat run-one --task <record.json> --out <new.json>`.

## 3. Your own runs

Write a structure file and a parameter file (`examples/`), `sprat plan`, `sprat run`, then the
analysis. The analysis selects records by their parameters, so runs that cover the same points
(the reference geometry, the analyte sweep with the spectra under both guide terminations, the
fine radius sweep, the verification series with both terminations, the field map, the bands and
the plane-wave layer) yield the same tables and figures; runs that cover fewer points yield the
subset the records support and list the skipped blocks in `ANALYSIS_REPORT.md`.

## Which command produces what

Figures and tables are numbered as in the paper.

| item | command | output |
|---|---|---|
| figure 1 (geometry) | `sprat figures` | drawn from the reference record |
| figure 2 (TM gap against the analyte index) | `sprat figures` | drawn from the band-structure records |
| table 1, figures 3(a, b) | `sprat analyze` | `analysis.json: base_series, radius_sweep` |
| figure 3(c) | `sprat analyze` | `kappa_series, resolved_pair` |
| table 2 | `sprat analyze` | `convergence` |
| table 3, figure 3(d) | `sprat analyze` | `cladding_ladder` |
| barrier mechanism and finite-cell correction, figure 4 | `sprat analyze` | `barrier` (panel (c) also reads the plane-wave record) |
| table 4, figure 5 | `sprat analyze` | `displacement` |
| transmission spectrum, figure 6 | `sprat analyze` | `spectra.json`, `analysis.json: spectra_summary`; the inset `spectra_absorber.json` |
| the spectra with the guide continued into the absorber, both grades, the cell factor across the sweep | `sprat analyze` | `analysis.json: termination` |
| table 5, sensitivity, FOM, DL, figure 7 | `sprat analyze` | `analyte_sweep`; the reflectionless columns and curves `analysis.json: combined.physical_values` |
| convergence of the resonance frequency | `sprat analyze` | `convergence`, `barrier.resolution_check` |
| analyte energy fraction, row decay | `sprat analyze` | `fields` |
| criteria fixed before the runs | `sprat grade` or `sprat analyze` | `predictions_report.md`; the third set in `analysis.json: termination` |
| registered numbers | `sprat analyze` | `numbers_registry.json`, `.md` |
| comparison with the deposit | `sprat audit` | `audit.json`, `AUDIT_REPORT.md` |
| figures 1 to 7 | `sprat figures` | `figures/figure1..7.pdf/.png` |
| data pack | `sprat deposit` | `deposit/` |

The comparison with published sensors (section S11 of the supplementary material) is not
computed from the records.
