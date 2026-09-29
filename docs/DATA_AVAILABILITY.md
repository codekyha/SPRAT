# Data availability statement

The statement of the paper, as `sprat deposit ... --raw-legacy raw_legacy --raw-h14b raw_h14b --data-doi 10.5281/zenodo.22912910 --software-doi 10.5281/zenodo.22912910 --include supplementary.pdf sprat-1.2.1.zip` writes it into `DATA_AVAILABILITY.md` of the data pack (version 3.1.0). SPRAT 1.2.1 is archived in the same record as the data. The arXiv version of the paper carries no data clause; the preprint (arXiv:2609.31952) is named in the record's metadata as a related work.

> The data that support the findings of this study are openly available at the following URL/DOI: https://doi.org/10.5281/zenodo.22912910. The deposit holds the harmonic-inversion records with their geometry parameters and validity flags, the transmission spectra with their reference runs, the stored field maps, the band-structure computations, the criteria fixed before the verification runs with their grading, the analysis outputs, the scripts that produced the records and the supplementary document, which is also available with this article. The same record holds SPRAT version 1.2.1, the software with which the records were converted and analysed; SPRAT is maintained at https://github.com/codekyha/sprat.

## Counts in the data pack

| quantity | value | source in the data pack |
|---|---|---|
| harmonic-inversion records | 858 | rows of records.csv with mode harminv |
| admitted (margin >= 1, gap position 0.2-0.8) | 741 | rows with valid = True |
| excluded, margin below 1 | 17 | exclusion_reason = Q > Q_lim |
| excluded, outside the gap window | 100 | exclusion_reason = gap position outside 0.2-0.8 |
| admitted, holding a mode of the finite waveguide (preliminary sweep: N_sep = 1 or r_d >= 0.13a) | 537 | rows with waveguide_mode = True |
| of these, in the dead zone (r_d >= 0.14a) | 436 | rows with valid = True and radius >= 0.14 |
| admitted, holding a cavity mode | 204 | rows with valid = True and waveguide_mode = False |
| clearing the margin alone | 841 | rows with resolved = True |
| spectra / references | 14 / 14 | spectra_records.tar.gz |
| of these, with the guide continued into the absorber | 7 / 7 | spectra_records.tar.gz, structure.cell.termination_x = absorber |
| field maps | 2 | field_records.tar.gz |
| band computations | 4 | bands_records.tar.gz |
