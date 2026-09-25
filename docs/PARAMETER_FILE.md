# Parameter file (.par) reference

## [run]

What to run and how to name it.

| key | type | default | unit | meaning |
|---|---|---|---|---|
| `mode` | harminv / spectrum / reference / field / bands / pwe | 'harminv' |  | harminv: complex eigenfrequency (Q, f_r) by harmonic inversion; spectrum: transmission P_out(f); reference: cavity-less transmission (normalisation); field: DFT field map at f_res and the analyte energy fraction; bands: MPB band structure; pwe: plane-wave-expansion checks |
| `label` | str | 'auto' |  | record name; 'auto' builds it from the resolved values |
| `tag` | str | '' |  | free text stored in every record (the name of a stage or a series) |
| `overwrite` | bool | False |  | rerun even if the record exists |
| `q_est` | float | 2000.0 |  | expected quality factor: sets the harminv time through the margin rule and the run-time ceilings of the spectrum and field modes |
| `quiet` | bool | True |  | silence the Meep log |

## [analyte]

The medium between the rods.

| key | type | default | unit | meaning |
|---|---|---|---|---|
| `n` | float | 1.33 |  | refractive index of the analyte |
| `k` | float | 0.0 |  | extinction coefficient; > 0 adds the corresponding conductivity |

## [numerics]

Discretisation and boundary numerics.

| key | type | default | unit | meaning |
|---|---|---|---|---|
| `resolution` | int | 24 | 1/a | grid points per lattice constant |
| `courant` | float | 0.5 |  | Courant factor S = c dt / dx |
| `pml` | float | 1.0 | a | PML thickness (all sides with pml termination; y sides with the absorber) |
| `subpixel` | bool | True |  | subpixel smoothing of the permittivity |
| `symmetry` | auto / none / even / odd | 'auto' |  | mirror symmetry in x: auto (even when dx = 0 and the structure is mirror-symmetric), none, even, odd |

## [source]

The Gaussian pulse.

| key | type | default | unit | meaning |
|---|---|---|---|---|
| `fcen` | float_or_auto | 'auto' |  | centre frequency a/lambda; auto = a_nm / target_wavelength_nm |
| `fwidth` | float | 0.06 |  | frequency width a/lambda |
| `cutoff` | float | 5.0 |  | Gaussian cutoff in widths (Meep default 5) |

## [harminv]

Harmonic inversion of E_z at a probe point next to the defect.

| key | type | default | unit | meaning |
|---|---|---|---|---|
| `t` | float_or_auto | 'auto' | a/c | signal length after the source is off (Meep time units); auto = margin * q_est / (pi fcen) rounded up to a multiple of 100 |
| `margin` | float | 4.5 |  | target Q_lim / Q_est for the auto rule |
| `t_max` | float | 30000.0 | a/c | ceiling of the second pass when auto_t is on |
| `auto_t` | bool | False |  | two-pass rule: rerun with t = auto_factor * Q / (pi fcen) when the first pass finds a Q that needs more than 1.5 t |
| `auto_factor` | float | 4.0 |  | factor of the two-pass rule |
| `probe_dx` | float | 0.13 | a | probe/source offset from the defect centre along x |
| `probe_dy` | float | 0.07 | a | probe/source offset from the defect centre along y |

## [spectrum]

Transmission spectrum and cavity-less reference.

| key | type | default | unit | meaning |
|---|---|---|---|---|
| `nfreq` | int | 201 |  | number of flux frequencies |
| `dft_tol` | float | 1e-06 |  | relative tolerance of the DFT-decay stopping rule |
| `t_max` | float | 40000.0 | a/c | hard ceiling of the run time after the source |
| `flux_width` | float | 3.0 | a | height of the flux planes across the guide |
| `source_width` | float | 2.0 | a | height of the line source across the guide |

## [field]

DFT field map at one frequency and the analyte energy fraction.

| key | type | default | unit | meaning |
|---|---|---|---|---|
| `f_res` | float | 0.0 |  | frequency of the map (a/lambda); required, normally the resonance of the matching harminv record |
| `t_max` | float | 20000.0 | a/c | hard ceiling of the run time after the source |

## [bands]

MPB band structures (mode = bands).

| key | type | default | unit | meaning |
|---|---|---|---|---|
| `task` | bulk / sweep / w1 / calibration | 'bulk' |  | bulk: TM and TE bands along Gamma-X-M-Gamma and the gap; sweep: TM gap against the analyte index; w1: W1 supercell projection; calibration: the TM gap centre and the lattice constant a = f_mid * target_wavelength_nm only |
| `resolution` | int | 32 |  | MPB grid points per lattice constant |
| `k_points` | int | 24 |  | interpolation points per segment of the k path |
| `num_bands` | int_or_auto | 'auto' |  | bands to compute; auto = 8 (bulk) or 24 (w1) |
| `supercell` | int | 15 |  | rows of the W1 supercell (odd) |
| `sweep_start` | float | 1.0 |  | first analyte index of the sweep |
| `sweep_stop` | float | 1.5 |  | last analyte index of the sweep |
| `sweep_step` | float | 0.05 |  | step of the sweep |

## [pwe]

Plane-wave-expansion checks (mode = pwe): bulk gap, W1 band, complex band structure, Fabry-Perot reading, barrier-row radius, loss budget.

| key | type | default | unit | meaning |
|---|---|---|---|---|
| `quick` | bool | False |  | smaller plane-wave cut-off (M = 7 instead of 10) for a fast check |
| `reference_f` | float | 0.30461 |  | resonance frequency of the reference point (r_d = 0.060a) |
| `second_f` | float | 0.2923 |  | resonance frequency of the second radius (r_d = 0.100a) |
| `kappa_reference` | float | 1.9595 |  | measured kappa at the reference point, for the deviation column |
| `kappa_second` | float | 1.7386 |  | measured kappa at the second radius, for the deviation column |

## [sweep]

Sweeps: `<dotted key> = a:b:step` or `v1, v2, ...`; `zip = key1, key2` pairs two lists element by element instead of taking their product.

Free keys allowed (each value is a str).
