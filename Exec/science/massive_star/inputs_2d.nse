# ------------------  INPUTS TO MAIN PROGRAM  -------------------

amr.plot_files_output = 1
amr.checkpoint_files_output = 1

max_step = 500000
stop_time = 3600

geometry.is_periodic = 0 0
geometry.coord_sys = 1         # r-z coordinates

geometry.prob_lo   =  0.    0.
geometry.prob_hi   =  1.75e10 1.75e10

amr.n_cell         = 2048 2048

amr.max_level      = 1      # maximum level number allowed

castro.lo_bc       =  3 2
castro.hi_bc       =  2 2

# >>>>>>>>>>>>>  BC FLAGS <<<<<<<<<<<<<<<<
# 0 = Interior           3 = Symmetry
# 1 = Inflow             4 = SlipWall
# 2 = Outflow            5 = NoSlipWall
# >>>>>>>>>>>>>  BC FLAGS <<<<<<<<<<<<<<<<

castro.do_hydro = 1
castro.do_grav  = 1
castro.do_react = 1
castro.do_sponge = 1

castro.ppm_type = 1
castro.ppm_temp_fix = 0

castro.use_flattening = 1

castro.riemann_solver = 1

gravity.gravity_type = MonopoleGrav
gravity.drdxfac = 2
castro.grav_source_type = 4

castro.sponge_upper_density = 1.e3
castro.sponge_lower_density = 1.e2
castro.sponge_timescale     = 1.e-3

castro.cfl            = 0.7     # cfl number for hyperbolic system
castro.init_shrink    = 0.01     # scale back initial timestep by this factor
castro.change_max     = 1.2    # factor by which dt is allowed to change each timestep
castro.sum_interval   = 0       # timesteps between computing and printing volume averages

#castro.dtnuc_e = 0.25
#castro.dtnuc_X = 0.25

amr.ref_ratio       = 2 2 2 2 # refinement ratio
amr.regrid_int      = 10000   # how often to regrid
amr.n_error_buf     = 2 2 2 2 # number of buffer cells in error est
amr.grid_eff        = 0.7     # what constitutes an efficient grid

amr.check_file      = chk     # root name of checkpoint file
amr.check_int       = 100     # number of timesteps between checkpoints
amr.plot_file       = plt     # root name of plot file
castro.plot_per_is_exact = 0
amr.plot_per = 0.1

amr.max_grid_size   = 512       # maximum grid size allowed -- used to control parallelism
amr.blocking_factor = 32       # block factor in grid generation

amr.v               = 1       # control verbosity in Amr.cpp
castro.v            = 1       # control verbosity in Castro.cpp

amr.derive_plot_vars = ALL

castro.small_dens   = 1.e-4
castro.small_temp   = 1.e6

castro.time_integration_method = 0
castro.use_retry = 1
castro.max_subcycles = 16

# problem initialization

problem.model_name =  "15m_500_sec.aprox19.hse.8192"

# refinement

amr.refinement_indicators = denerr

amr.refine.denerr.max_level = 3
amr.refine.denerr.value_greater = 1.e4
amr.refine.denerr.field_name = density


# Microphysics

integrator.rtol_spec = 1.e-6
integrator.atol_spec = 1.e-6

integrator.retry_burn = 0
integrator.abort_on_failure = 0

integrator.jacobian = 2

network.rho_nse = 2.e6
network.T_nse = 3.e9
network.nse_relax_factor = 0.9

