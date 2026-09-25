# The task list of the data set of the paper

Everything needed to regenerate the data set of the paper, point by point, on a workstation. The
deposited records were produced on the UHeM Altay cluster by the original scripts of 2026
(`02_kavite.py` with the v4 patch, Meep 1.30.0) and converted into the SPRAT schema with
`sprat import-legacy`; the tasks below run the same points with SPRAT.

| file | content |
|---|---|
| `tasks.jsonl` | 892 fully resolved tasks built from the deposited records by `build_campaign.py` (one exact duplicate merged): 858 harmonic inversions, 14 spectra and 14 cavity-less references (seven of each with the guide ending in the PML, seven with it continued into the absorber), 2 field maps, 4 band computations, the plane-wave layer. The first 878 tasks were built on 2026-09-23 from the deposit of 2026; SPRAT 1.2.0 appended the 14 runs of `raw_h14b` with `build_campaign.py --append`, leaving the earlier lines unchanged. Every task carries `legacy_labels`, the names of the deposited files it regenerates. |
| `tasks_summary.md` | tasks and single-core hours per tier |
| `structure_reference.phc`, `structure_absorber.phc` | the reference geometry with the two guide terminations |
| `stages/*.par` | the parameter files of the stages as a reader would write them (the authoritative task set is `tasks.jsonl`, which also carries the odd extra points of the runs of 2026, such as repeated runs with longer signals) |
| `expected/numbers_registry_expected.json` | the numbers registry of the unpublished version 2.0.2 of the data set |
| `expected/numbers_registry_manuscript_v7.json` | the numbers registry of manuscript version 7 (the `numbers_registry.json` of the 2.2.0 deposit, byte for byte) |
| `expected/numbers_registry_manuscript_v8.json` | the numbers registry of manuscript version 8: version 7 plus the section `corrected_v8` (two statements corrected against the records) |
| `expected/numbers_registry_manuscript_v9.json` | the numbers registry of manuscript version 9, recomputed with SPRAT 1.1.0 (sections `corrected_v9` and `superseded_v9`) |
| `expected/numbers_registry_manuscript_v10.json` | the numbers registry of manuscript version 10, recomputed with SPRAT 1.2.0 (sections `termination` and `corrected_v10`) |
| `expected/numbers_registry_manuscript_v10_2.json` | the numbers registry of manuscript version 10.2: version 10 with the section `derived_v10_2` (the reflectionless quality factor and figure of merit at the rows of table 5, the convergence of the resonance frequency, the channel at the guided-mode wavevector, the absorption bracket), the three references added and the literature value quoted; copied into the data pack as `analysis/manuscript_registry.json` |
| `expected/numbers_registry_sprat.json` | the registry recomputed by SPRAT from the imported deposit |
| `expected/records_expected.csv` | the flat table of the 2026 deposit (version 2.2.0, 879 records) in the SPRAT schema |

## Tiers

```
sprat reproduce campaigns/manuscript --tier quick --records records     # bands, plane-wave layer, smoke tests, the reference record
sprat reproduce campaigns/manuscript --tier core  --records records     # every record of the systematic runs behind figures 2, 3 and 5 to 7
sprat reproduce campaigns/manuscript --tier full  --records records     # everything, including the 54 verification runs of the barrier analysis (figure 4) and the 14 absorber-terminated spectrum runs
```

`--dry-run` writes the task list of the tier (`<records>/_campaign/tasks_<tier>.jsonl`, with its `.tsv` and `.summary.md`) and prints the cost table without running anything; this directory itself is never written to. The tier `campaign` of `tasks.jsonl` holds the systematic runs.

## Cost

| tier | tasks | cost-model estimate at the calibrated Altay throughput (single-core hours) | measured on Altay as run (core-hours) |
|---|---|---|---|
| quick | 15 | 0.7 | - |
| core | 824 | 26 | about 133 |
| full | 892 | 87 | about 256 |

The two columns differ because the coarse sweep of 2026 (693 small cells at resolution 20) ran at a much lower pixel-step rate than the production runs (about 3.6e7 against 3.7e8 pixel-steps per second per task), and the cost model uses one calibrated throughput for every cell. The 14 absorber-terminated spectrum runs (job 528207, 24 September 2026) ran as 14 concurrent tasks of four cores each on one node: 14.4 hours of task time, 58 core-hours, 1 h 54 min of wall time; the cost model gives them 6 single-core hours. Run `sprat calibrate` on your machine before trusting the estimate; the longest single task (N_sep = 6, n_cl = 14) is 13.6 h on one Altay core and sets the floor of the wall time whatever the worker count.

## Comparing a regeneration with the deposit

```
sprat audit records --expected campaigns/manuscript/expected/numbers_registry_manuscript_v10_2.json -o tables
```

The audit recomputes the registered numbers of the paper from the regenerated records and lists every deviation from the registry, row by row (`numbers_registry_expected.json` can be given instead). With the same Meep version the measured Q values are expected to agree to about 1e-3 relative (the tolerance of `sprat selftest --reference`) and the derived quantities to their printed precision; another Meep version or CPU can move Q at the 1e-3 level.
