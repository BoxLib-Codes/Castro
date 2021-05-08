&fortin

  rho_fuel = 2.d6
  T_fuel = 5.e7

  T_ash = 2.e9
  
  pert_frac = 0.2d0
  pert_delta = 0.05d0

/

&tagging

  denerr = 1.d-7
  dengrad = 0.01
  max_denerr_lev = 5
  max_dengrad_lev = 5

  presserr = 1.d20
  pressgrad = 1.d20
  max_presserr_lev = 5
  max_pressgrad_lev = 5

/

&extern
  rtol_spec = 1.d-8
  atol_spec = 1.d-8

  small_x = 1.e-10

  jacobian = 2
/
