&fortin

  dtemp = 1.2e9
  x_half_max = 2.048d4
  x_half_width = 2048.d0

  dx_model = 10.d0

  dens_base = 3.43e6

  T_star = 1.d8
  T_hi = 2.d8
  T_lo   = 5.d6

  H_star = 2000.d0
  atm_delta  = 50.0

  fuel1_name = "helium-4"
  fuel1_frac = 1.0d0

  ash1_name  = "nickel-56"
  ash1_frac = 1.0d0

  low_density_cutoff = 1.d-4

  cutoff_density = 2.5e4
  max_hse_tagging_level = 3
  max_base_tagging_level = 2

  X_min = 1.e-2

  x_refine_distance = 6.144e4
/

&sponge

  sponge_upper_density = 1.d2
  sponge_lower_density = 1.d0
  sponge_timescale     = 1.d-7

/


&extern

  rtol_spec = 1.d-6
  atol_spec = 1.d-6

  react_boost = 10.0
  use_tables = T
/