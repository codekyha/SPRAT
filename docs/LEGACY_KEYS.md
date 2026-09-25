# Legacy key map

The original scripts of 2026 wrote their records with Turkish keys.  `sprat import-legacy` converts them with this map; the original files are never modified and every converted record stores `provenance.legacy_file` and `provenance.legacy_sha256`.

| legacy | English | meaning |
|---|---|---|
| `parametreler` | `structure + params` | the parameter block splits into the structure and the parameter file |
| `parametreler.mod` | `params.run.mode` | harminv -> harminv, spektrum -> spectrum, norm -> reference, alan -> field |
| `parametreler.katman` | `params.run.tag` | tarama -> sweep, uretim -> production, elle -> manual |
| `parametreler.etiket` | `task.label` | the record name, kept verbatim |
| `parametreler.a_nm` | `structure.lattice.a_nm` | lattice constant (nm) |
| `parametreler.r` | `structure.lattice.rod_radius` | rod radius r/a |
| `parametreler.rd` | `structure.defect.radius` | defect radius r_d/a |
| `parametreler.dx, dy` | `structure.defect.dx, dy` | defect displacement (a) |
| `parametreler.nsep` | `structure.defect.row` | separating rows N_sep |
| `parametreler.nx` | `structure.cell.guide_periods` | guide periods |
| `parametreler.nkaplama` | `structure.cell.cladding_rows` | cladding rows per side n_cl |
| `parametreler.pad, ypad` | `structure.cell.pad_x, pad_y` | paddings (a) |
| `parametreler.sonlandirma` | `structure.cell.termination_x` | pml -> pml, emici -> absorber |
| `parametreler.emici_kalinlik` | `structure.cell.absorber_periods` | absorber thickness (a) |
| `parametreler.bariyer_satir_r / bariyer_satirlari` | `structure.rows` | 'j:r' -> {j: r} |
| `parametreler.kavite` | `(mode)` | False only for the cavity-less reference run (mode reference) |
| `parametreler.na` | `params.analyte.n` | analyte index |
| `parametreler.kappa` | `params.analyte.k` | analyte extinction coefficient |
| `parametreler.res` | `params.numerics.resolution` | grid points per a |
| `parametreler.dpml` | `params.numerics.pml` | PML thickness (a) |
| `parametreler.simetri` | `params.numerics.symmetry` | oto -> auto, yok -> none, cift -> even, tek -> odd |
| `parametreler.fcen, df` | `params.source.fcen, fwidth` | source centre and width |
| `parametreler.t_harminv` | `params.harminv.t` | requested harminv signal length |
| `parametreler.t_max_harminv` | `params.harminv.t_max` | ceiling of the two-pass rule |
| `parametreler.otomatik_t, otomatik_kat` | `params.harminv.auto_t, auto_factor` | two-pass rule |
| `parametreler.q_tahmin` | `params.run.q_est` | expected Q |
| `parametreler.nfreq, dft_tol, t_max_spektrum` | `params.spectrum.nfreq, dft_tol, t_max` | spectrum settings |
| `parametreler.f_rez, t_alan_max` | `params.field.f_res, t_max` | field-map settings |
| `parametreler.sessiz` | `params.run.quiet` | Meep log silenced |
| `olcum` | `provenance` | the measurement block |
| `olcum.baslangic` | `provenance.started` | start time |
| `olcum.dugum` | `provenance.host` | node name |
| `olcum.is_no, dizi_no` | `provenance.job_id, array_id` | SLURM job and array ids |
| `olcum.adim` | `provenance.steps` | time steps |
| `olcum.zaman` | `provenance.sim_time` | simulated time |
| `olcum.piksel` | `provenance.pixels` | grid points |
| `olcum.piksel_adim` | `provenance.pixel_steps` | work |
| `olcum.piksel_adim_s` | `provenance.throughput` | pixel-steps per second |
| `olcum.duvar_s` | `provenance.wall_s` | wall seconds |
| `olcum.simetri` | `provenance.symmetry` | sector actually used |
| `olcum.kod_sha256` | `software.code_sha256` | SHA-256 of the script that wrote the record (v4 only) |
| `olcum.paket` | `provenance.legacy_package` | 'v4' for the verification runs |
| `sonuc` | `result` | the result block |
| `sonuc.modlar[i].f, Q, genlik, hata, lambda_nm, FWHM_nm` | `result.modes[i].f, Q, amplitude, error, wavelength_nm, fwhm_nm` | harmonic-inversion modes, strongest first |
| `sonuc.t_harminv_kullanilan` | `result.t_used` | signal length actually used |
| `sonuc.otomatik_gecmis` | `result.auto_history` | the passes of the two-pass rule |
| `sonuc.Q_cozunurluk_siniri` | `result.Q_limit` | pi fcen t |
| `sonuc.f, P_in, P_out, t_bitis` | `result.f, flux_in, flux_out, t_end` | spectrum and reference runs |
| `sonuc.f_kesri, lambda_r_nm, S_teori_nm_RIU` | `result.analyte_energy_fraction, wavelength_nm, S_first_order_nm_per_RIU` | field runs |
| `alan_<label>.npz: Ez, eps, sx, sy, dpml, res, f_rez, lambda_nm` | `<label>.npz: Ez, eps, sx, sy, pml, resolution, f_res, wavelength_nm (+ Lx, Ly, edge_x, n, rod_eps)` | field maps; the importer adds the interior size and the medium |
| `bant_bulk_na*.json, bant_tarama.json, bant_w1_na*.json, a_nm.json` | `bands records (task bulk / sweep / w1)` | gap: var -> exists, alt -> lower, ust -> upper, orta -> mid, yuzde -> relative_width_pct |
| `pwe_v4_feasibility_results.json` | `pwe record` | kappa keys 'r_d=0.060a (reference)' -> reference, 'r_d=0.100a' -> second |
