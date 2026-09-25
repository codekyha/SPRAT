# Criteria fixed before the runs

Three sets of criteria were fixed before the runs they judge. The first two were graded afterwards
without being edited; the third was amended once while its runs were under way, after its
reference runs had finished and before its first cavity spectrum was written, and is graded both
as fixed and as amended. The frozen files are kept here byte for byte (Turkish keys in the first
two sets, as written on the day); their names are the names of the day, and "preregistered" in the first one
is such a name. The SHA-256 values of the first two sets are recorded in
`sprat/data/predictions_criteria.json`, the English transcription that the grader reads;
`tests/test_predictions.py` checks that every expected value of the English file equals the value
in the frozen original. The third set and its amendment are read by
`sprat/analysis/termination.py`, whose constants `tests/test_termination.py` checks against them.

| file | frozen on | content |
|---|---|---|
| `frozen/predictions_preregistered_2026-09-17.json` | 2026-09-17 | P1 to P13 (frozen 2026-09-16) and the two convergence criteria H8 and H10 (recorded 2026-09-17), before the first verification run |
| `frozen/predictions_second_set_2026-09-18.json` | 2026-09-18 | PA1 to PA6 for the second-radius runs and PB1 to PB5 for a resolution study, with the grading of 2026-09-20 and the closure note of the PB criteria |
| `frozen/predictions_ek_2026-09-18_sealed.json` | 2026-09-18 | the sealed copy of the second set as frozen, SHA-256 2d1395be...cbd55, before grading |
| `frozen/predictions_S4b_2026-09-24.json` | 2026-09-24 | D-22, the criterion for the transmission spectra with the guide continued into the absorber (H14b), fixed before job 528207 was submitted: at each analyte index the factor \|q\| in the 25a cell over \|q\| with the absorber; 3 or more supports the attribution of the residual asymmetry to the guide ends, 1.5 or less refutes it; a pass needs support at all seven indices, the gates G1 to G5 and S1. SHA-256 e06be498...bfac4 |
| `frozen/predictions_S4b_amendment_A1A2_2026-09-24.json` | 2026-09-24 15:46:39 +03:00 | amendment A1 (the ratio graded at the five indices where \|q\| in the 25a cell exceeds 0.2) and A2 (the cell factor across the analyte sweep, a measurement), fixed after the reference runs had finished (15:00:54) and before the first absorber cavity spectrum was written (15:53:07, +03:00); the time is the one the file records, and its copies carry 15:48:38. SHA-256 84d13c85...bb1c3 |

`sprat grade <records>` (or `sprat analyze`) evaluates every criterion of the first two sets from
records selected by their parameters and writes `predictions_report.md`. `sprat analyze` grades
the third set in `analysis.json`, key `termination`, by D-22 and by A1 side by side: D-22 is not
passed (the factor at $n_a$ = 1.300 is 1.61, inconclusive), A1 is passed (factors of 8.95 or more
at the five graded indices). Every failure is reported first and explained afterwards; no
criterion is changed after its results are seen.
