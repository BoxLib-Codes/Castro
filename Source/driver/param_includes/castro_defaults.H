
// This file is automatically created by parse_castro_params.py.  To update
// or add runtime parameters, please edit _cpp_parameters and then run
// mk_params.sh

int         Castro::state_interp_order = 1;
int         Castro::lin_limit_state_interp = 0;
int         Castro::state_nghost = 0;
int         Castro::do_reflux = 1;
int         Castro::update_sources_after_reflux = 1;
int         Castro::use_custom_knapsack_weights = 0;
amrex::Real Castro::difmag = 0.1;
amrex::Real Castro::small_dens = -1.e200;
amrex::Real Castro::small_temp = -1.e200;
amrex::Real Castro::small_pres = -1.e200;
amrex::Real Castro::small_ener = -1.e200;
int         Castro::do_hydro = -1;
int         Castro::time_integration_method = 0;
int         Castro::limit_fourth_order = 1;
int         Castro::initialization_is_cell_average = 0;
int         Castro::use_reconstructed_gamma1 = 0;
int         Castro::add_ext_src = 0;
int         Castro::hybrid_hydro = 0;
int         Castro::ppm_type = 1;
int         Castro::ppm_temp_fix = 0;
int         Castro::ppm_predict_gammae = 0;
int         Castro::ppm_reference_eigenvectors = 0;
int         Castro::plm_iorder = 2;
int         Castro::plm_limiter = 2;
int         Castro::plm_well_balanced = 0;
int         Castro::hybrid_riemann = 0;
int         Castro::riemann_solver = 0;
int         Castro::cg_maxiter = 12;
amrex::Real Castro::cg_tol = 1.0e-5;
int         Castro::cg_blend = 2;
int         Castro::use_eos_in_riemann = 0;
amrex::Real Castro::riemann_speed_limit = 2.99792458e10;
int         Castro::use_flattening = 1;
int         Castro::transverse_use_eos = 0;
int         Castro::transverse_reset_density = 1;
int         Castro::transverse_reset_rhoe = 0;
amrex::Real Castro::dual_energy_eta1 = 1.0e0;
amrex::Real Castro::dual_energy_eta2 = 1.0e-4;
int         Castro::use_pslope = 1;
int         Castro::limit_fluxes_on_small_dens = 0;
int         Castro::density_reset_method = 1;
int         Castro::allow_small_energy = 1;
int         Castro::do_sponge = 0;
int         Castro::sponge_implicit = 1;
int         Castro::source_term_predictor = 0;
int         Castro::first_order_hydro = 0;
std::string Castro::xl_ext_bc_type = "";
std::string Castro::xr_ext_bc_type = "";
std::string Castro::yl_ext_bc_type = "";
std::string Castro::yr_ext_bc_type = "";
std::string Castro::zl_ext_bc_type = "";
std::string Castro::zr_ext_bc_type = "";
int         Castro::hse_zero_vels = 0;
int         Castro::hse_interp_temp = 0;
int         Castro::hse_reflect_vels = 0;
int         Castro::sdc_order = 2;
int         Castro::sdc_quadrature = 0;
int         Castro::sdc_extra = 0;
int         Castro::sdc_solver = 1;
amrex::Real Castro::sdc_solver_tol_dens = 1.e-6;
amrex::Real Castro::sdc_solver_tol_spec = 1.e-6;
amrex::Real Castro::sdc_solver_tol_ener = 1.e-6;
amrex::Real Castro::sdc_solver_atol = 1.e-10;
amrex::Real Castro::sdc_solver_relax_factor = 1.0;
int         Castro::sdc_solve_for_rhoe = 1;
int         Castro::sdc_use_analytic_jac = 1;
amrex::Real Castro::fixed_dt = -1.0;
amrex::Real Castro::initial_dt = -1.0;
amrex::Real Castro::dt_cutoff = 0.0;
amrex::Real Castro::max_dt = 1.e200;
amrex::Real Castro::cfl = 0.8;
amrex::Real Castro::init_shrink = 1.0;
amrex::Real Castro::change_max = 1.1;
int         Castro::plot_per_is_exact = 0;
int         Castro::small_plot_per_is_exact = 0;
int         Castro::use_retry = 0;
amrex::Real Castro::retry_tolerance = 0.02;
amrex::Real Castro::retry_neg_dens_factor = 1.e-1;
amrex::Real Castro::retry_subcycle_factor = 0.5;
int         Castro::use_post_step_regrid = 0;
int         Castro::max_subcycles = 10;
int         Castro::sdc_iters = 2;
amrex::Real Castro::dtnuc_e = 1.e200;
amrex::Real Castro::dtnuc_X = 1.e200;
amrex::Real Castro::dtnuc_X_threshold = 1.e-3;
int         Castro::do_react = -1;
amrex::Real Castro::react_T_min = 0.0;
amrex::Real Castro::react_T_max = 1.e200;
amrex::Real Castro::react_rho_min = 0.0;
amrex::Real Castro::react_rho_max = 1.e200;
int         Castro::disable_shock_burning = 0;
amrex::Real Castro::T_guess = 1.e8;
int         Castro::do_grav = -1;
int         Castro::moving_center = 0;
int         Castro::grav_source_type = 4;
int         Castro::do_rotation = -1;
int         Castro::do_acc = -1;
int         Castro::bndry_func_thread_safe = 1;
int         Castro::grown_factor = 1;
int         Castro::star_at_center = -1;
int         Castro::do_scf_initial_model = 0;
amrex::Real Castro::scf_maximum_density = -1.e6;
amrex::Real Castro::scf_equatorial_radius = -1.e9;
amrex::Real Castro::scf_polar_radius = -1.e9;
amrex::Real Castro::scf_relax_tol = 1.e-3;
int         Castro::scf_max_iterations = 30;
int         Castro::do_special_tagging = 0;
int         Castro::spherical_star = 0;
int         Castro::verbose = 0;
bool        Castro::dump_old = false;
#ifdef AMREX_DEBUG
int         Castro::print_fortran_warnings = 1;
#else
int         Castro::print_fortran_warnings = 0;
#endif
#ifdef AMREX_DEBUG
int         Castro::print_update_diagnostics = 1;
#else
int         Castro::print_update_diagnostics = 0;
#endif
int         Castro::track_grid_losses = 0;
int         Castro::sum_interval = -1;
amrex::Real Castro::sum_per = -1.0e0;
int         Castro::show_center_of_mass = 0;
int         Castro::hard_cfl_limit = 1;
std::string Castro::job_name = "";
int         Castro::output_at_completion = 1;
amrex::Real Castro::reset_checkpoint_time = -1.e200;
int         Castro::reset_checkpoint_step = -1;
#ifdef GRAVITY
int         Castro::use_point_mass = 0;
amrex::Real Castro::point_mass = 0.0;
int         Castro::point_mass_fix_solution = 0;
#endif
#ifdef AMREX_PARTICLES
int         Castro::do_tracer_particles = 0;
#endif
#ifdef ROTATION
amrex::Real Castro::rotational_period = -1.e200;
amrex::Real Castro::rotational_dPdt = 0.0;
int         Castro::rotation_include_centrifugal = 1;
int         Castro::rotation_include_coriolis = 1;
int         Castro::rotation_include_domegadt = 1;
int         Castro::state_in_rotating_frame = 1;
int         Castro::rot_source_type = 4;
int         Castro::implicit_rotation_update = 1;
int         Castro::rot_axis = 3;
#endif
#ifdef DIFFUSION
int         Castro::diffuse_temp = 0;
amrex::Real Castro::diffuse_cutoff_density = -1.e200;
amrex::Real Castro::diffuse_cutoff_density_hi = -1.e200;
amrex::Real Castro::diffuse_cond_scale_fac = 1.0;
#endif
