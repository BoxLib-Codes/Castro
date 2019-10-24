# Building and running with the stellarcollapse EOS:

To build with the StellarCollapse EOS in Microphysics, pass
`EOS_DIR=stellarcollapse` as a make flag.

Also, replace the following in the extern namelist in probin.

For weaklib EOS:

```
&extern

  eos_input_is_constant = .true.
  weaklib_eos_table_name = "EquationOfStateTable.h5"
  small_x = 1.d-30
  use_tables = .false.
  use_c12ag_deboer17 = .false.

/
```

For stellarcollapse EOS:

```
&extern
  eos_file = "./Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.h5"
  use_energy_shift = T
/
```

