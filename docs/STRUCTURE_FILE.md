# Structure file (.phc) reference

## [lattice]

The periodic lattice and its materials.

| key | type | default | unit | meaning |
|---|---|---|---|---|
| `type` | square | 'square' |  | lattice type (only square is implemented in 1.x) |
| `a_nm` | float_or_calibrate | 481.4 | nm | lattice constant in nanometres, or 'calibrate' to take f_mid * target_wavelength_nm from the most recent bulk band record |
| `target_wavelength_nm` | float | 1550.0 | nm | wavelength that the band-gap centre is calibrated to |
| `rod_radius` | float | 0.2 | a | rod radius r/a |
| `rod_eps` | float | 11.9025 |  | rod permittivity (silicon: 3.45^2) |
| `background` | analyte | 'analyte' |  | medium between the rods: the analyte of the parameter file |

## [cell]

The finite computational cell.

| key | type | default | unit | meaning |
|---|---|---|---|---|
| `guide_periods` | int | 25 |  | number of lattice columns along the guide (n_x) |
| `cladding_rows` | int | 12 |  | rows of rods on each side of the guide (n_cl) |
| `pad_x` | float | 1.5 | a | free analyte between the last column and the PML (PML termination only) |
| `pad_y` | float | 1.0 | a | free analyte between the last row and the PML |
| `termination_x` | pml / absorber | 'pml' |  | guide termination along x: 'pml' (padding + PML) or 'absorber' (lattice and guide continue into an adiabatic absorber) |
| `absorber_periods` | float | 8.0 | a | thickness of the adiabatic absorber (absorber termination only) |

## [waveguide]

The line defect.

| key | type | default | unit | meaning |
|---|---|---|---|---|
| `type` | w1 / none | 'w1' |  | W1: the row j = 0 is removed; none: full lattice |

## [defect]

The point defect (the cavity).

| key | type | default | unit | meaning |
|---|---|---|---|---|
| `type` | rod / none | 'rod' |  | rod: the lattice rod at (0, row) is replaced by a rod of the given radius; none: no cavity |
| `radius` | float | 0.06 | a | defect rod radius r_d/a (0 removes the rod entirely) |
| `row` | int | 4 |  | row of the defect counted from the guide (N_sep) |
| `dx` | float | 0.0 | a | displacement of the defect rod along the guide |
| `dy` | float | 0.0 | a | displacement of the defect rod across the guide |

## [rows]

Per-row radius overrides: `<row index> = <radius/a>`; rows are counted from the guide, positive on the cavity side, negative on the other side.

Free keys allowed (each value is a float).

## [rods]

Free-form edits applied after everything else, one per line: `add x y r` or `remove x y` (units of a).

Bare lines (without `=`) allowed.
