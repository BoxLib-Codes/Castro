# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step = 50000000

stop_time    = 0.07

# PROBLEM SIZE & GEOMETRY
geometry.is_periodic = 0
geometry.coord_sys   = 0        # 0 => cart, 1 => RZ  2=>spherical
geometry.prob_lo     = 0  
geometry.prob_hi     = 8.192e4
amr.n_cell           = 256

# >>>>>>>>>>>>>  BC FLAGS <<<<<<<<<<<<<<<<
# 0 = Interior           3 = Symmetry
# 1 = Inflow             4 = SlipWall
# 2 = Outflow            5 = NoSlipWall
# >>>>>>>>>>>>>  BC FLAGS <<<<<<<<<<<<<<<<
castro.lo_bc       =  3
castro.hi_bc       =  2

# WHICH PHYSICS
castro.do_hydro = 1
castro.do_react = 1

castro.time_integration_method = 2
castro.sdc_order = 4

castro.sdc_solve_for_rhoe = 1
castro.sdc_solver_tol_dens = 1.e-10
castro.sdc_solver_tol_spec = 1.e-10
castro.sdc_solver_tol_ener = 1.e-6
castro.sdc_solver_atol = 1.e-10
castro.sdc_solver = 1
castro.sdc_use_analytic_jac = 1
castro.sdc_solver_relax_factor = 1

castro.use_reconstructed_gamma1 = 1

castro.small_temp = 1.e6
castro.small_dens = 1.e-5

castro.ppm_type = 1

castro.diffuse_temp = 1
castro.diffuse_cutoff_density = 1.e-2


# TIME STEP CONTROL
castro.cfl            = 0.6    # cfl number for hyperbolic system
castro.init_shrink    = 0.1     # scale back initial timestep
castro.change_max     = 1.1     # max time step growth

castro.dtnuc_e        = 0.25


# DIAGNOSTICS & VERBOSITY
castro.sum_interval   = 1       # timesteps between computing mass
amr.data_log          = "toy_flame.log"
castro.v              = 1       # verbosity in Castro.cpp
amr.v                 = 1       # verbosity in Amr.cpp

# REFINEMENT / REGRIDDING
amr.max_level       = 1       # maximum level number allowed
amr.ref_ratio       = 4 4  2 2 2 # refinement ratio
amr.regrid_int      = 2 2 2 2 # how often to regrid
amr.blocking_factor = 64      # block factor in grid generation
amr.max_grid_size   = 64
amr.n_error_buf     = 2 2 2 2 # number of buffer cells in error est

# CHECKPOINT FILES
amr.checkpoint_files_output = 1
amr.check_file      = chk        # root name of checkpoint file
amr.check_int       = 10000       # number of timesteps between checkpoints

# PLOTFILES
amr.plot_file       = plt        # root name of plotfile
amr.plot_per        = 1.e-4
castro.plot_per_is_exact = 1
amr.derive_plot_vars = ALL

#PROBIN FILENAME
amr.probin_file = probin.sdc
