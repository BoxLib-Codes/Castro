
! This file is automatically created by parse_castro_params.py.  To update
! or add runtime parameters, please edit _cpp_parameters and then run
! mk_params.sh

! This module stores the runtime parameters and integer names for
! indexing arrays.
!
! The Fortran-specific parameters are initialized in set_method_params(),
! and the ones that we are mirroring from C++ and obtaining through the
! ParmParse module are initialized in ca_set_castro_method_params().

module meth_params_module

  use castro_error_module
  use amrex_fort_module, only: rt => amrex_real
  use state_sizes_module, only : nadv, NQAUX, NVAR, NGDNV, NQ, NQSRC

  implicit none

  ! number of ghost cells for the hyperbolic solver
  integer, parameter     :: NHYP    = 4

  ! conservative variables
  integer, allocatable, save :: URHO, UMX, UMY, UMZ, UMR, UML, UMP, UEDEN, UEINT, UTEMP, UFA, UFS, UFX
  integer, allocatable, save :: USHK

  ! primitive variables
  integer, allocatable, save :: QRHO, QU, QV, QW, QPRES, QREINT, QTEMP, QGAME, QGC
  integer, allocatable, save :: QGAMC, QC, QDPDR, QDPDE
#ifdef RADIATION
  integer, allocatable, save :: QGAMCG, QCG, QLAMS
#endif
  integer, allocatable, save :: QFA, QFS, QFX

#ifdef RADIATION
  integer, save :: QRAD, QRADHI, QPTOT, QREITOT
  integer, save :: fspace_type
  logical, save :: do_inelastic_scattering
  logical, save :: comoving

  real(rt)        , save :: flatten_pp_threshold = -1.e0_rt
#endif

  integer, save, allocatable :: npassive
  integer, save, allocatable :: qpass_map(:), upass_map(:)

  ! These are used for the Godunov state
  ! Note that the velocity indices here are picked to be the same value
  ! as in the primitive variable array
  integer, save, allocatable :: GDRHO, GDU, GDV, GDW, GDPRES, GDGAME
#ifdef RADIATION
  integer, save, allocatable :: GDLAMS, GDERADS
#endif

  ! Numerical values corresponding to the gravity types
#ifdef GRAVITY
  integer, save, allocatable :: gravity_type_int
  integer, parameter :: ConstantGrav = 0
  integer, parameter :: MonopoleGrav = 1
  integer, parameter :: PoissonGrav = 2
  integer, parameter :: PrescribedGrav = 3
#endif

  integer         , save :: numpts_1d

  real(rt)        , save, allocatable :: outflow_data_old(:,:)
  real(rt)        , save, allocatable :: outflow_data_new(:,:)
  real(rt)        , save :: outflow_data_old_time
  real(rt)        , save :: outflow_data_new_time
  logical         , save :: outflow_data_allocated
  real(rt)        , save :: max_dist

  ! these flags are for interpreting the EXT_DIR BCs
  integer, parameter :: EXT_UNDEFINED = -1
  integer, parameter :: EXT_HSE = 1
  integer, parameter :: EXT_INTERP = 2

  integer, allocatable, save :: xl_ext, yl_ext, zl_ext, xr_ext, yr_ext, zr_ext

  ! Create versions of these variables on the GPU
  ! the device update is then done in Castro_nd.f90

#ifdef AMREX_USE_CUDA
  attributes(managed) :: URHO, UMX, UMY, UMZ, UMR, UML, UMP, UEDEN, UEINT, UTEMP, UFA, UFS, UFX
  attributes(managed) :: USHK
  attributes(managed) :: QRHO, QU, QV, QW, QPRES, QREINT, QTEMP, QGAME, QGC
  attributes(managed) :: QGAMC, QC, QDPDR, QDPDE
#ifdef RADIATION
  attributes(managed) :: QGAMCG, QCG, QLAMS
#endif
  attributes(managed) :: QFA, QFS, QFX
  attributes(managed) :: npassive
  attributes(managed) :: qpass_map, upass_map
  attributes(managed) :: GDRHO, GDU, GDV, GDW, GDPRES, GDGAME
#ifdef RADIATION
  attributes(managed) :: GDLAMS, GDERADS
#endif
#ifdef GRAVITY
  attributes(managed) :: gravity_type_int
#endif
  attributes(managed) :: xl_ext, yl_ext, zl_ext, xr_ext, yr_ext, zr_ext
#endif

  !$acc declare &
  !$acc create(URHO, UMX, UMY, UMZ, UMR, UML, UMP, UEDEN, UEINT, UTEMP, UFA, UFS,UFX) &
  !$acc create(USHK) &
  !$acc create(QRHO, QU, QV, QW, QPRES, QREINT, QTEMP) &
  !$acc create(QC, QDPDR, QDPDE, QGAMC, QGAME, QGC) &
#ifdef RADIATION
  !$acc create(QGAMCG, QCG, QLAMS) &
  !$acc create(QRAD, QRADHI, QPTOT, QREITOT) &
  !$acc create(fspace_type, do_inelastic_scattering, comoving) &
#endif
  !$acc create(QFA, QFS, QFX) &
  !$acc create(xl_ext, yl_ext, zl_ext, xr_ext, yr_ext, zr_ext)

  ! Begin the declarations of the ParmParse parameters

  real(rt), allocatable, save :: difmag
  real(rt), allocatable, save :: small_dens
  real(rt), allocatable, save :: small_temp
  real(rt), allocatable, save :: small_pres
  real(rt), allocatable, save :: small_ener
  integer,  allocatable, save :: do_hydro
  integer,  allocatable, save :: time_integration_method
  integer,  allocatable, save :: limit_fourth_order
  integer,  allocatable, save :: use_reconstructed_gamma1
  integer,  allocatable, save :: hybrid_hydro
  integer,  allocatable, save :: ppm_type
  integer,  allocatable, save :: ppm_temp_fix
  integer,  allocatable, save :: ppm_predict_gammae
  integer,  allocatable, save :: ppm_reference_eigenvectors
  integer,  allocatable, save :: plm_iorder
  integer,  allocatable, save :: plm_limiter
  integer,  allocatable, save :: plm_well_balanced
  integer,  allocatable, save :: hybrid_riemann
  integer,  allocatable, save :: riemann_solver
  integer,  allocatable, save :: cg_maxiter
  real(rt), allocatable, save :: cg_tol
  integer,  allocatable, save :: cg_blend
  integer,  allocatable, save :: use_eos_in_riemann
  real(rt), allocatable, save :: riemann_speed_limit
  integer,  allocatable, save :: use_flattening
  integer,  allocatable, save :: transverse_use_eos
  integer,  allocatable, save :: transverse_reset_density
  integer,  allocatable, save :: transverse_reset_rhoe
  real(rt), allocatable, save :: dual_energy_eta1
  real(rt), allocatable, save :: dual_energy_eta2
  integer,  allocatable, save :: use_pslope
  integer,  allocatable, save :: limit_fluxes_on_small_dens
  integer,  allocatable, save :: density_reset_method
  integer,  allocatable, save :: allow_small_energy
  integer,  allocatable, save :: do_sponge
  integer,  allocatable, save :: sponge_implicit
  integer,  allocatable, save :: first_order_hydro
  character (len=:), allocatable, save :: xl_ext_bc_type
  character (len=:), allocatable, save :: xr_ext_bc_type
  character (len=:), allocatable, save :: yl_ext_bc_type
  character (len=:), allocatable, save :: yr_ext_bc_type
  character (len=:), allocatable, save :: zl_ext_bc_type
  character (len=:), allocatable, save :: zr_ext_bc_type
  integer,  allocatable, save :: hse_zero_vels
  integer,  allocatable, save :: hse_interp_temp
  integer,  allocatable, save :: hse_reflect_vels
  integer,  allocatable, save :: sdc_order
  integer,  allocatable, save :: sdc_quadrature
  integer,  allocatable, save :: sdc_extra
  integer,  allocatable, save :: sdc_solver
  real(rt), allocatable, save :: sdc_solver_tol_dens
  real(rt), allocatable, save :: sdc_solver_tol_spec
  real(rt), allocatable, save :: sdc_solver_tol_ener
  real(rt), allocatable, save :: sdc_solver_atol
  real(rt), allocatable, save :: sdc_solver_relax_factor
  integer,  allocatable, save :: sdc_solve_for_rhoe
  integer,  allocatable, save :: sdc_use_analytic_jac
  real(rt), allocatable, save :: cfl
  real(rt), allocatable, save :: dtnuc_e
  real(rt), allocatable, save :: dtnuc_X
  real(rt), allocatable, save :: dtnuc_X_threshold
  integer,  allocatable, save :: do_react
  real(rt), allocatable, save :: react_T_min
  real(rt), allocatable, save :: react_T_max
  real(rt), allocatable, save :: react_rho_min
  real(rt), allocatable, save :: react_rho_max
  integer,  allocatable, save :: disable_shock_burning
  real(rt), allocatable, save :: T_guess
  integer,  allocatable, save :: diffuse_temp
  real(rt), allocatable, save :: diffuse_cutoff_density
  real(rt), allocatable, save :: diffuse_cutoff_density_hi
  real(rt), allocatable, save :: diffuse_cond_scale_fac
  integer,  allocatable, save :: do_grav
  integer,  allocatable, save :: grav_source_type
  integer,  allocatable, save :: do_rotation
  real(rt), allocatable, save :: rot_period
  real(rt), allocatable, save :: rot_period_dot
  integer,  allocatable, save :: rotation_include_centrifugal
  integer,  allocatable, save :: rotation_include_coriolis
  integer,  allocatable, save :: rotation_include_domegadt
  integer,  allocatable, save :: state_in_rotating_frame
  integer,  allocatable, save :: rot_source_type
  integer,  allocatable, save :: implicit_rotation_update
  integer,  allocatable, save :: rot_axis
  integer,  allocatable, save :: use_point_mass
  real(rt), allocatable, save :: point_mass
  integer,  allocatable, save :: point_mass_fix_solution
  integer,  allocatable, save :: do_acc
  integer,  allocatable, save :: grown_factor
  integer,  allocatable, save :: track_grid_losses
  character (len=:), allocatable, save :: gravity_type
  real(rt), allocatable, save :: const_grav
  integer,  allocatable, save :: get_g_from_phi

#ifdef AMREX_USE_CUDA
attributes(managed) :: difmag
attributes(managed) :: small_dens
attributes(managed) :: small_temp
attributes(managed) :: small_pres
attributes(managed) :: small_ener
attributes(managed) :: do_hydro
attributes(managed) :: time_integration_method
attributes(managed) :: limit_fourth_order
attributes(managed) :: use_reconstructed_gamma1
attributes(managed) :: hybrid_hydro
attributes(managed) :: ppm_type
attributes(managed) :: ppm_temp_fix
attributes(managed) :: ppm_predict_gammae
attributes(managed) :: ppm_reference_eigenvectors
attributes(managed) :: plm_iorder
attributes(managed) :: plm_limiter
attributes(managed) :: plm_well_balanced
attributes(managed) :: hybrid_riemann
attributes(managed) :: riemann_solver
attributes(managed) :: cg_maxiter
attributes(managed) :: cg_tol
attributes(managed) :: cg_blend
attributes(managed) :: use_eos_in_riemann
attributes(managed) :: riemann_speed_limit
attributes(managed) :: use_flattening
attributes(managed) :: transverse_use_eos
attributes(managed) :: transverse_reset_density
attributes(managed) :: transverse_reset_rhoe
attributes(managed) :: dual_energy_eta1
attributes(managed) :: dual_energy_eta2
attributes(managed) :: use_pslope
attributes(managed) :: limit_fluxes_on_small_dens
attributes(managed) :: density_reset_method
attributes(managed) :: allow_small_energy
attributes(managed) :: do_sponge
attributes(managed) :: sponge_implicit
attributes(managed) :: first_order_hydro






attributes(managed) :: hse_zero_vels
attributes(managed) :: hse_interp_temp
attributes(managed) :: hse_reflect_vels
attributes(managed) :: sdc_order
attributes(managed) :: sdc_quadrature
attributes(managed) :: sdc_extra
attributes(managed) :: sdc_solver
attributes(managed) :: sdc_solver_tol_dens
attributes(managed) :: sdc_solver_tol_spec
attributes(managed) :: sdc_solver_tol_ener
attributes(managed) :: sdc_solver_atol
attributes(managed) :: sdc_solver_relax_factor
attributes(managed) :: sdc_solve_for_rhoe
attributes(managed) :: sdc_use_analytic_jac
attributes(managed) :: cfl
attributes(managed) :: dtnuc_e
attributes(managed) :: dtnuc_X
attributes(managed) :: dtnuc_X_threshold
attributes(managed) :: do_react
attributes(managed) :: react_T_min
attributes(managed) :: react_T_max
attributes(managed) :: react_rho_min
attributes(managed) :: react_rho_max
attributes(managed) :: disable_shock_burning
attributes(managed) :: T_guess
#ifdef DIFFUSION
attributes(managed) :: diffuse_temp
#endif
#ifdef DIFFUSION
attributes(managed) :: diffuse_cutoff_density
#endif
#ifdef DIFFUSION
attributes(managed) :: diffuse_cutoff_density_hi
#endif
#ifdef DIFFUSION
attributes(managed) :: diffuse_cond_scale_fac
#endif
attributes(managed) :: do_grav
attributes(managed) :: grav_source_type
attributes(managed) :: do_rotation
#ifdef ROTATION
attributes(managed) :: rot_period
#endif
#ifdef ROTATION
attributes(managed) :: rot_period_dot
#endif
#ifdef ROTATION
attributes(managed) :: rotation_include_centrifugal
#endif
#ifdef ROTATION
attributes(managed) :: rotation_include_coriolis
#endif
#ifdef ROTATION
attributes(managed) :: rotation_include_domegadt
#endif
#ifdef ROTATION
attributes(managed) :: state_in_rotating_frame
#endif
#ifdef ROTATION
attributes(managed) :: rot_source_type
#endif
#ifdef ROTATION
attributes(managed) :: implicit_rotation_update
#endif
#ifdef ROTATION
attributes(managed) :: rot_axis
#endif
#ifdef GRAVITY
attributes(managed) :: use_point_mass
#endif
#ifdef GRAVITY
attributes(managed) :: point_mass
#endif
#ifdef GRAVITY
attributes(managed) :: point_mass_fix_solution
#endif
attributes(managed) :: do_acc
attributes(managed) :: grown_factor
attributes(managed) :: track_grid_losses

attributes(managed) :: const_grav
attributes(managed) :: get_g_from_phi
#endif

  !$acc declare &
  !$acc create(difmag) &
  !$acc create(small_dens) &
  !$acc create(small_temp) &
  !$acc create(small_pres) &
  !$acc create(small_ener) &
  !$acc create(do_hydro) &
  !$acc create(time_integration_method) &
  !$acc create(limit_fourth_order) &
  !$acc create(use_reconstructed_gamma1) &
  !$acc create(hybrid_hydro) &
  !$acc create(ppm_type) &
  !$acc create(ppm_temp_fix) &
  !$acc create(ppm_predict_gammae) &
  !$acc create(ppm_reference_eigenvectors) &
  !$acc create(plm_iorder) &
  !$acc create(plm_limiter) &
  !$acc create(plm_well_balanced) &
  !$acc create(hybrid_riemann) &
  !$acc create(riemann_solver) &
  !$acc create(cg_maxiter) &
  !$acc create(cg_tol) &
  !$acc create(cg_blend) &
  !$acc create(use_eos_in_riemann) &
  !$acc create(riemann_speed_limit) &
  !$acc create(use_flattening) &
  !$acc create(transverse_use_eos) &
  !$acc create(transverse_reset_density) &
  !$acc create(transverse_reset_rhoe) &
  !$acc create(dual_energy_eta1) &
  !$acc create(dual_energy_eta2) &
  !$acc create(use_pslope) &
  !$acc create(limit_fluxes_on_small_dens) &
  !$acc create(density_reset_method) &
  !$acc create(allow_small_energy) &
  !$acc create(do_sponge) &
  !$acc create(sponge_implicit) &
  !$acc create(first_order_hydro) &
  !$acc create(hse_zero_vels) &
  !$acc create(hse_interp_temp) &
  !$acc create(hse_reflect_vels) &
  !$acc create(sdc_order) &
  !$acc create(sdc_quadrature) &
  !$acc create(sdc_extra) &
  !$acc create(sdc_solver) &
  !$acc create(sdc_solver_tol_dens) &
  !$acc create(sdc_solver_tol_spec) &
  !$acc create(sdc_solver_tol_ener) &
  !$acc create(sdc_solver_atol) &
  !$acc create(sdc_solver_relax_factor) &
  !$acc create(sdc_solve_for_rhoe) &
  !$acc create(sdc_use_analytic_jac) &
  !$acc create(cfl) &
  !$acc create(dtnuc_e) &
  !$acc create(dtnuc_X) &
  !$acc create(dtnuc_X_threshold) &
  !$acc create(do_react) &
  !$acc create(react_T_min) &
  !$acc create(react_T_max) &
  !$acc create(react_rho_min) &
  !$acc create(react_rho_max) &
  !$acc create(disable_shock_burning) &
  !$acc create(T_guess) &
#ifdef DIFFUSION
  !$acc create(diffuse_temp) &
#endif
#ifdef DIFFUSION
  !$acc create(diffuse_cutoff_density) &
#endif
#ifdef DIFFUSION
  !$acc create(diffuse_cutoff_density_hi) &
#endif
#ifdef DIFFUSION
  !$acc create(diffuse_cond_scale_fac) &
#endif
  !$acc create(do_grav) &
  !$acc create(grav_source_type) &
  !$acc create(do_rotation) &
#ifdef ROTATION
  !$acc create(rot_period) &
#endif
#ifdef ROTATION
  !$acc create(rot_period_dot) &
#endif
#ifdef ROTATION
  !$acc create(rotation_include_centrifugal) &
#endif
#ifdef ROTATION
  !$acc create(rotation_include_coriolis) &
#endif
#ifdef ROTATION
  !$acc create(rotation_include_domegadt) &
#endif
#ifdef ROTATION
  !$acc create(state_in_rotating_frame) &
#endif
#ifdef ROTATION
  !$acc create(rot_source_type) &
#endif
#ifdef ROTATION
  !$acc create(implicit_rotation_update) &
#endif
#ifdef ROTATION
  !$acc create(rot_axis) &
#endif
#ifdef GRAVITY
  !$acc create(use_point_mass) &
#endif
#ifdef GRAVITY
  !$acc create(point_mass) &
#endif
#ifdef GRAVITY
  !$acc create(point_mass_fix_solution) &
#endif
  !$acc create(do_acc) &
  !$acc create(grown_factor) &
  !$acc create(track_grid_losses) &
  !$acc create(const_grav) &
  !$acc create(get_g_from_phi)

  ! End the declarations of the ParmParse parameters

  real(rt)        , save :: rot_vec(3)

contains

  subroutine ca_set_castro_method_params() bind(C, name="ca_set_castro_method_params")

    use amrex_parmparse_module, only: amrex_parmparse_build, amrex_parmparse_destroy, amrex_parmparse

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type (amrex_parmparse) :: pp


    allocate(URHO, UMX, UMY, UMZ, UMR, UML, UMP, UEDEN, UEINT, UTEMP, UFA, UFS, UFX)
    allocate(USHK)
    allocate(QRHO, QU, QV, QW, QPRES, QREINT, QTEMP, QGAME, QGC)
    allocate(QGAMC, QC, QDPDR, QDPDE)
#ifdef RADIATION
    allocate(QGAMCG, QCG, QLAMS)
#endif
    allocate(QFA, QFS, QFX)
    allocate(npassive)
    allocate(GDRHO, GDU, GDV, GDW, GDPRES, GDGAME)
#ifdef RADIATION
    allocate(GDLAMS, GDERADS)
#endif
    allocate(xl_ext, yl_ext, zl_ext, xr_ext, yr_ext, zr_ext)

    allocate(difmag)
    difmag = 0.1_rt;
    allocate(small_dens)
    small_dens = -1.e200_rt;
    allocate(small_temp)
    small_temp = -1.e200_rt;
    allocate(small_pres)
    small_pres = -1.e200_rt;
    allocate(small_ener)
    small_ener = -1.e200_rt;
    allocate(do_hydro)
    do_hydro = -1;
    allocate(time_integration_method)
    time_integration_method = 0;
    allocate(limit_fourth_order)
    limit_fourth_order = 1;
    allocate(use_reconstructed_gamma1)
    use_reconstructed_gamma1 = 0;
    allocate(hybrid_hydro)
    hybrid_hydro = 0;
    allocate(ppm_type)
    ppm_type = 1;
    allocate(ppm_temp_fix)
    ppm_temp_fix = 0;
    allocate(ppm_predict_gammae)
    ppm_predict_gammae = 0;
    allocate(ppm_reference_eigenvectors)
    ppm_reference_eigenvectors = 0;
    allocate(plm_iorder)
    plm_iorder = 2;
    allocate(plm_limiter)
    plm_limiter = 2;
    allocate(plm_well_balanced)
    plm_well_balanced = 0;
    allocate(hybrid_riemann)
    hybrid_riemann = 0;
    allocate(riemann_solver)
    riemann_solver = 0;
    allocate(cg_maxiter)
    cg_maxiter = 12;
    allocate(cg_tol)
    cg_tol = 1.0e-5_rt;
    allocate(cg_blend)
    cg_blend = 2;
    allocate(use_eos_in_riemann)
    use_eos_in_riemann = 0;
    allocate(riemann_speed_limit)
    riemann_speed_limit = 2.99792458e10_rt;
    allocate(use_flattening)
    use_flattening = 1;
    allocate(transverse_use_eos)
    transverse_use_eos = 0;
    allocate(transverse_reset_density)
    transverse_reset_density = 1;
    allocate(transverse_reset_rhoe)
    transverse_reset_rhoe = 0;
    allocate(dual_energy_eta1)
    dual_energy_eta1 = 1.0e0_rt;
    allocate(dual_energy_eta2)
    dual_energy_eta2 = 1.0e-4_rt;
    allocate(use_pslope)
    use_pslope = 1;
    allocate(limit_fluxes_on_small_dens)
    limit_fluxes_on_small_dens = 0;
    allocate(density_reset_method)
    density_reset_method = 1;
    allocate(allow_small_energy)
    allow_small_energy = 1;
    allocate(do_sponge)
    do_sponge = 0;
    allocate(sponge_implicit)
    sponge_implicit = 1;
    allocate(first_order_hydro)
    first_order_hydro = 0;
    allocate(character(len=1)::xl_ext_bc_type)
    xl_ext_bc_type = "";
    allocate(character(len=1)::xr_ext_bc_type)
    xr_ext_bc_type = "";
    allocate(character(len=1)::yl_ext_bc_type)
    yl_ext_bc_type = "";
    allocate(character(len=1)::yr_ext_bc_type)
    yr_ext_bc_type = "";
    allocate(character(len=1)::zl_ext_bc_type)
    zl_ext_bc_type = "";
    allocate(character(len=1)::zr_ext_bc_type)
    zr_ext_bc_type = "";
    allocate(hse_zero_vels)
    hse_zero_vels = 0;
    allocate(hse_interp_temp)
    hse_interp_temp = 0;
    allocate(hse_reflect_vels)
    hse_reflect_vels = 0;
    allocate(sdc_order)
    sdc_order = 2;
    allocate(sdc_quadrature)
    sdc_quadrature = 0;
    allocate(sdc_extra)
    sdc_extra = 0;
    allocate(sdc_solver)
    sdc_solver = 1;
    allocate(sdc_solver_tol_dens)
    sdc_solver_tol_dens = 1.e-6_rt;
    allocate(sdc_solver_tol_spec)
    sdc_solver_tol_spec = 1.e-6_rt;
    allocate(sdc_solver_tol_ener)
    sdc_solver_tol_ener = 1.e-6_rt;
    allocate(sdc_solver_atol)
    sdc_solver_atol = 1.e-10_rt;
    allocate(sdc_solver_relax_factor)
    sdc_solver_relax_factor = 1.0_rt;
    allocate(sdc_solve_for_rhoe)
    sdc_solve_for_rhoe = 1;
    allocate(sdc_use_analytic_jac)
    sdc_use_analytic_jac = 1;
    allocate(cfl)
    cfl = 0.8_rt;
    allocate(dtnuc_e)
    dtnuc_e = 1.e200_rt;
    allocate(dtnuc_X)
    dtnuc_X = 1.e200_rt;
    allocate(dtnuc_X_threshold)
    dtnuc_X_threshold = 1.e-3_rt;
    allocate(do_react)
    do_react = -1;
    allocate(react_T_min)
    react_T_min = 0.0_rt;
    allocate(react_T_max)
    react_T_max = 1.e200_rt;
    allocate(react_rho_min)
    react_rho_min = 0.0_rt;
    allocate(react_rho_max)
    react_rho_max = 1.e200_rt;
    allocate(disable_shock_burning)
    disable_shock_burning = 0;
    allocate(T_guess)
    T_guess = 1.e8_rt;
    allocate(do_grav)
    do_grav = -1;
    allocate(grav_source_type)
    grav_source_type = 4;
    allocate(do_rotation)
    do_rotation = -1;
    allocate(do_acc)
    do_acc = -1;
    allocate(grown_factor)
    grown_factor = 1;
    allocate(track_grid_losses)
    track_grid_losses = 0;
#ifdef GRAVITY
    allocate(use_point_mass)
    use_point_mass = 0;
    allocate(point_mass)
    point_mass = 0.0_rt;
    allocate(point_mass_fix_solution)
    point_mass_fix_solution = 0;
#endif
#ifdef ROTATION
    allocate(rot_period)
    rot_period = -1.e200_rt;
    allocate(rot_period_dot)
    rot_period_dot = 0.0_rt;
    allocate(rotation_include_centrifugal)
    rotation_include_centrifugal = 1;
    allocate(rotation_include_coriolis)
    rotation_include_coriolis = 1;
    allocate(rotation_include_domegadt)
    rotation_include_domegadt = 1;
    allocate(state_in_rotating_frame)
    state_in_rotating_frame = 1;
    allocate(rot_source_type)
    rot_source_type = 4;
    allocate(implicit_rotation_update)
    implicit_rotation_update = 1;
    allocate(rot_axis)
    rot_axis = 3;
#endif
#ifdef DIFFUSION
    allocate(diffuse_temp)
    diffuse_temp = 0;
    allocate(diffuse_cutoff_density)
    diffuse_cutoff_density = -1.e200_rt;
    allocate(diffuse_cutoff_density_hi)
    diffuse_cutoff_density_hi = -1.e200_rt;
    allocate(diffuse_cond_scale_fac)
    diffuse_cond_scale_fac = 1.0_rt;
#endif

    call amrex_parmparse_build(pp, "castro")
    call pp%query("difmag", difmag)
    call pp%query("small_dens", small_dens)
    call pp%query("small_temp", small_temp)
    call pp%query("small_pres", small_pres)
    call pp%query("small_ener", small_ener)
    call pp%query("do_hydro", do_hydro)
    call pp%query("time_integration_method", time_integration_method)
    call pp%query("limit_fourth_order", limit_fourth_order)
    call pp%query("use_reconstructed_gamma1", use_reconstructed_gamma1)
    call pp%query("hybrid_hydro", hybrid_hydro)
    call pp%query("ppm_type", ppm_type)
    call pp%query("ppm_temp_fix", ppm_temp_fix)
    call pp%query("ppm_predict_gammae", ppm_predict_gammae)
    call pp%query("ppm_reference_eigenvectors", ppm_reference_eigenvectors)
    call pp%query("plm_iorder", plm_iorder)
    call pp%query("plm_limiter", plm_limiter)
    call pp%query("plm_well_balanced", plm_well_balanced)
    call pp%query("hybrid_riemann", hybrid_riemann)
    call pp%query("riemann_solver", riemann_solver)
    call pp%query("cg_maxiter", cg_maxiter)
    call pp%query("cg_tol", cg_tol)
    call pp%query("cg_blend", cg_blend)
    call pp%query("use_eos_in_riemann", use_eos_in_riemann)
    call pp%query("riemann_speed_limit", riemann_speed_limit)
    call pp%query("use_flattening", use_flattening)
    call pp%query("transverse_use_eos", transverse_use_eos)
    call pp%query("transverse_reset_density", transverse_reset_density)
    call pp%query("transverse_reset_rhoe", transverse_reset_rhoe)
    call pp%query("dual_energy_eta1", dual_energy_eta1)
    call pp%query("dual_energy_eta2", dual_energy_eta2)
    call pp%query("use_pslope", use_pslope)
    call pp%query("limit_fluxes_on_small_dens", limit_fluxes_on_small_dens)
    call pp%query("density_reset_method", density_reset_method)
    call pp%query("allow_small_energy", allow_small_energy)
    call pp%query("do_sponge", do_sponge)
    call pp%query("sponge_implicit", sponge_implicit)
    call pp%query("first_order_hydro", first_order_hydro)
    call pp%query("xl_ext_bc_type", xl_ext_bc_type)
    call pp%query("xr_ext_bc_type", xr_ext_bc_type)
    call pp%query("yl_ext_bc_type", yl_ext_bc_type)
    call pp%query("yr_ext_bc_type", yr_ext_bc_type)
    call pp%query("zl_ext_bc_type", zl_ext_bc_type)
    call pp%query("zr_ext_bc_type", zr_ext_bc_type)
    call pp%query("hse_zero_vels", hse_zero_vels)
    call pp%query("hse_interp_temp", hse_interp_temp)
    call pp%query("hse_reflect_vels", hse_reflect_vels)
    call pp%query("sdc_order", sdc_order)
    call pp%query("sdc_quadrature", sdc_quadrature)
    call pp%query("sdc_extra", sdc_extra)
    call pp%query("sdc_solver", sdc_solver)
    call pp%query("sdc_solver_tol_dens", sdc_solver_tol_dens)
    call pp%query("sdc_solver_tol_spec", sdc_solver_tol_spec)
    call pp%query("sdc_solver_tol_ener", sdc_solver_tol_ener)
    call pp%query("sdc_solver_atol", sdc_solver_atol)
    call pp%query("sdc_solver_relax_factor", sdc_solver_relax_factor)
    call pp%query("sdc_solve_for_rhoe", sdc_solve_for_rhoe)
    call pp%query("sdc_use_analytic_jac", sdc_use_analytic_jac)
    call pp%query("cfl", cfl)
    call pp%query("dtnuc_e", dtnuc_e)
    call pp%query("dtnuc_X", dtnuc_X)
    call pp%query("dtnuc_X_threshold", dtnuc_X_threshold)
    call pp%query("do_react", do_react)
    call pp%query("react_T_min", react_T_min)
    call pp%query("react_T_max", react_T_max)
    call pp%query("react_rho_min", react_rho_min)
    call pp%query("react_rho_max", react_rho_max)
    call pp%query("disable_shock_burning", disable_shock_burning)
    call pp%query("T_guess", T_guess)
    call pp%query("do_grav", do_grav)
    call pp%query("grav_source_type", grav_source_type)
    call pp%query("do_rotation", do_rotation)
    call pp%query("do_acc", do_acc)
    call pp%query("grown_factor", grown_factor)
    call pp%query("track_grid_losses", track_grid_losses)
#ifdef GRAVITY
    call pp%query("use_point_mass", use_point_mass)
    call pp%query("point_mass", point_mass)
    call pp%query("point_mass_fix_solution", point_mass_fix_solution)
#endif
#ifdef ROTATION
    call pp%query("rotational_period", rot_period)
    call pp%query("rotational_dPdt", rot_period_dot)
    call pp%query("rotation_include_centrifugal", rotation_include_centrifugal)
    call pp%query("rotation_include_coriolis", rotation_include_coriolis)
    call pp%query("rotation_include_domegadt", rotation_include_domegadt)
    call pp%query("state_in_rotating_frame", state_in_rotating_frame)
    call pp%query("rot_source_type", rot_source_type)
    call pp%query("implicit_rotation_update", implicit_rotation_update)
    call pp%query("rot_axis", rot_axis)
#endif
#ifdef DIFFUSION
    call pp%query("diffuse_temp", diffuse_temp)
    call pp%query("diffuse_cutoff_density", diffuse_cutoff_density)
    call pp%query("diffuse_cutoff_density_hi", diffuse_cutoff_density_hi)
    call pp%query("diffuse_cond_scale_fac", diffuse_cond_scale_fac)
#endif
    call amrex_parmparse_destroy(pp)


    allocate(character(len=1)::gravity_type)
    gravity_type = "fillme";
    allocate(const_grav)
    const_grav = 0.0_rt;
    allocate(get_g_from_phi)
    get_g_from_phi = 0;

    call amrex_parmparse_build(pp, "gravity")
    call pp%query("gravity_type", gravity_type)
    call pp%query("const_grav", const_grav)
    call pp%query("get_g_from_phi", get_g_from_phi)
    call amrex_parmparse_destroy(pp)



    !$acc update &
    !$acc device(difmag) &
    !$acc device(small_dens) &
    !$acc device(small_temp) &
    !$acc device(small_pres) &
    !$acc device(small_ener) &
    !$acc device(do_hydro) &
    !$acc device(time_integration_method) &
    !$acc device(limit_fourth_order) &
    !$acc device(use_reconstructed_gamma1) &
    !$acc device(hybrid_hydro) &
    !$acc device(ppm_type) &
    !$acc device(ppm_temp_fix) &
    !$acc device(ppm_predict_gammae) &
    !$acc device(ppm_reference_eigenvectors) &
    !$acc device(plm_iorder) &
    !$acc device(plm_limiter) &
    !$acc device(plm_well_balanced) &
    !$acc device(hybrid_riemann) &
    !$acc device(riemann_solver) &
    !$acc device(cg_maxiter) &
    !$acc device(cg_tol) &
    !$acc device(cg_blend) &
    !$acc device(use_eos_in_riemann) &
    !$acc device(riemann_speed_limit) &
    !$acc device(use_flattening) &
    !$acc device(transverse_use_eos) &
    !$acc device(transverse_reset_density) &
    !$acc device(transverse_reset_rhoe) &
    !$acc device(dual_energy_eta1) &
    !$acc device(dual_energy_eta2) &
    !$acc device(use_pslope) &
    !$acc device(limit_fluxes_on_small_dens) &
    !$acc device(density_reset_method) &
    !$acc device(allow_small_energy) &
    !$acc device(do_sponge) &
    !$acc device(sponge_implicit) &
    !$acc device(first_order_hydro) &
    !$acc device(hse_zero_vels) &
    !$acc device(hse_interp_temp) &
    !$acc device(hse_reflect_vels) &
    !$acc device(sdc_order) &
    !$acc device(sdc_quadrature) &
    !$acc device(sdc_extra) &
    !$acc device(sdc_solver) &
    !$acc device(sdc_solver_tol_dens) &
    !$acc device(sdc_solver_tol_spec) &
    !$acc device(sdc_solver_tol_ener) &
    !$acc device(sdc_solver_atol) &
    !$acc device(sdc_solver_relax_factor) &
    !$acc device(sdc_solve_for_rhoe) &
    !$acc device(sdc_use_analytic_jac) &
    !$acc device(cfl) &
    !$acc device(dtnuc_e) &
    !$acc device(dtnuc_X) &
    !$acc device(dtnuc_X_threshold) &
    !$acc device(do_react) &
    !$acc device(react_T_min) &
    !$acc device(react_T_max) &
    !$acc device(react_rho_min) &
    !$acc device(react_rho_max) &
    !$acc device(disable_shock_burning) &
    !$acc device(T_guess) &
#ifdef DIFFUSION
    !$acc device(diffuse_temp) &
#endif
#ifdef DIFFUSION
    !$acc device(diffuse_cutoff_density) &
#endif
#ifdef DIFFUSION
    !$acc device(diffuse_cutoff_density_hi) &
#endif
#ifdef DIFFUSION
    !$acc device(diffuse_cond_scale_fac) &
#endif
    !$acc device(do_grav) &
    !$acc device(grav_source_type) &
    !$acc device(do_rotation) &
#ifdef ROTATION
    !$acc device(rot_period) &
#endif
#ifdef ROTATION
    !$acc device(rot_period_dot) &
#endif
#ifdef ROTATION
    !$acc device(rotation_include_centrifugal) &
#endif
#ifdef ROTATION
    !$acc device(rotation_include_coriolis) &
#endif
#ifdef ROTATION
    !$acc device(rotation_include_domegadt) &
#endif
#ifdef ROTATION
    !$acc device(state_in_rotating_frame) &
#endif
#ifdef ROTATION
    !$acc device(rot_source_type) &
#endif
#ifdef ROTATION
    !$acc device(implicit_rotation_update) &
#endif
#ifdef ROTATION
    !$acc device(rot_axis) &
#endif
#ifdef GRAVITY
    !$acc device(use_point_mass) &
#endif
#ifdef GRAVITY
    !$acc device(point_mass) &
#endif
#ifdef GRAVITY
    !$acc device(point_mass_fix_solution) &
#endif
    !$acc device(do_acc) &
    !$acc device(grown_factor) &
    !$acc device(track_grid_losses) &
    !$acc device(const_grav) &
    !$acc device(get_g_from_phi)


#ifdef GRAVITY
    ! Set the gravity type integer

    allocate(gravity_type_int)

    if (gravity_type == "ConstantGrav") then
       gravity_type_int = ConstantGrav
    else if (gravity_type == "MonopoleGrav") then
       gravity_type_int = MonopoleGrav
    else if (gravity_type == "PoissonGrav") then
       gravity_type_int = PoissonGrav
    else if (gravity_type == "PrescribedGrav") then
       gravity_type_int = PrescribedGrav
    else
       call castro_error("Unknown gravity type")
    end if
#endif

    ! now set the external BC flags
    select case (xl_ext_bc_type)
    case ("hse", "HSE")
       xl_ext = EXT_HSE
    case ("interp", "INTERP")
       xl_ext = EXT_INTERP
    case default
       xl_ext = EXT_UNDEFINED
    end select

    select case (yl_ext_bc_type)
    case ("hse", "HSE")
       yl_ext = EXT_HSE
    case ("interp", "INTERP")
       yl_ext = EXT_INTERP
    case default
       yl_ext = EXT_UNDEFINED
    end select

    select case (zl_ext_bc_type)
    case ("hse", "HSE")
       zl_ext = EXT_HSE
    case ("interp", "INTERP")
       zl_ext = EXT_INTERP
    case default
       zl_ext = EXT_UNDEFINED
    end select

    select case (xr_ext_bc_type)
    case ("hse", "HSE")
       xr_ext = EXT_HSE
    case ("interp", "INTERP")
       xr_ext = EXT_INTERP
    case default
       xr_ext = EXT_UNDEFINED
    end select

    select case (yr_ext_bc_type)
    case ("hse", "HSE")
       yr_ext = EXT_HSE
    case ("interp", "INTERP")
       yr_ext = EXT_INTERP
    case default
       yr_ext = EXT_UNDEFINED
    end select

    select case (zr_ext_bc_type)
    case ("hse", "HSE")
       zr_ext = EXT_HSE
    case ("interp", "INTERP")
       zr_ext = EXT_INTERP
    case default
       zr_ext = EXT_UNDEFINED
    end select

    !$acc update device(xl_ext, yl_ext, zl_ext, xr_ext, yr_ext, zr_ext)


  end subroutine ca_set_castro_method_params


  subroutine ca_finalize_meth_params() bind(C, name="ca_finalize_meth_params")
    implicit none

    deallocate(URHO, UMX, UMY, UMZ, UMR, UML, UMP, UEDEN, UEINT, UTEMP, UFA, UFS, UFX)
    deallocate(USHK)
    deallocate(QRHO, QU, QV, QW, QPRES, QREINT, QTEMP, QGAME, QGC)
    deallocate(QGAMC, QC, QDPDR, QDPDE)
#ifdef RADIATION
    deallocate(QGAMCG, QCG, QLAMS)
#endif
    deallocate(QFA, QFS, QFX)
    deallocate(npassive)
    deallocate(GDRHO, GDU, GDV, GDW, GDPRES, GDGAME)
#ifdef RADIATION
    deallocate(GDLAMS, GDERADS)
    deallocate(xl_ext, yl_ext, zl_ext, xr_ext, yr_ext, zr_ext)
#endif

    if (allocated(difmag)) then
        deallocate(difmag)
    end if
    if (allocated(small_dens)) then
        deallocate(small_dens)
    end if
    if (allocated(small_temp)) then
        deallocate(small_temp)
    end if
    if (allocated(small_pres)) then
        deallocate(small_pres)
    end if
    if (allocated(small_ener)) then
        deallocate(small_ener)
    end if
    if (allocated(do_hydro)) then
        deallocate(do_hydro)
    end if
    if (allocated(time_integration_method)) then
        deallocate(time_integration_method)
    end if
    if (allocated(limit_fourth_order)) then
        deallocate(limit_fourth_order)
    end if
    if (allocated(use_reconstructed_gamma1)) then
        deallocate(use_reconstructed_gamma1)
    end if
    if (allocated(hybrid_hydro)) then
        deallocate(hybrid_hydro)
    end if
    if (allocated(ppm_type)) then
        deallocate(ppm_type)
    end if
    if (allocated(ppm_temp_fix)) then
        deallocate(ppm_temp_fix)
    end if
    if (allocated(ppm_predict_gammae)) then
        deallocate(ppm_predict_gammae)
    end if
    if (allocated(ppm_reference_eigenvectors)) then
        deallocate(ppm_reference_eigenvectors)
    end if
    if (allocated(plm_iorder)) then
        deallocate(plm_iorder)
    end if
    if (allocated(plm_limiter)) then
        deallocate(plm_limiter)
    end if
    if (allocated(plm_well_balanced)) then
        deallocate(plm_well_balanced)
    end if
    if (allocated(hybrid_riemann)) then
        deallocate(hybrid_riemann)
    end if
    if (allocated(riemann_solver)) then
        deallocate(riemann_solver)
    end if
    if (allocated(cg_maxiter)) then
        deallocate(cg_maxiter)
    end if
    if (allocated(cg_tol)) then
        deallocate(cg_tol)
    end if
    if (allocated(cg_blend)) then
        deallocate(cg_blend)
    end if
    if (allocated(use_eos_in_riemann)) then
        deallocate(use_eos_in_riemann)
    end if
    if (allocated(riemann_speed_limit)) then
        deallocate(riemann_speed_limit)
    end if
    if (allocated(use_flattening)) then
        deallocate(use_flattening)
    end if
    if (allocated(transverse_use_eos)) then
        deallocate(transverse_use_eos)
    end if
    if (allocated(transverse_reset_density)) then
        deallocate(transverse_reset_density)
    end if
    if (allocated(transverse_reset_rhoe)) then
        deallocate(transverse_reset_rhoe)
    end if
    if (allocated(dual_energy_eta1)) then
        deallocate(dual_energy_eta1)
    end if
    if (allocated(dual_energy_eta2)) then
        deallocate(dual_energy_eta2)
    end if
    if (allocated(use_pslope)) then
        deallocate(use_pslope)
    end if
    if (allocated(limit_fluxes_on_small_dens)) then
        deallocate(limit_fluxes_on_small_dens)
    end if
    if (allocated(density_reset_method)) then
        deallocate(density_reset_method)
    end if
    if (allocated(allow_small_energy)) then
        deallocate(allow_small_energy)
    end if
    if (allocated(do_sponge)) then
        deallocate(do_sponge)
    end if
    if (allocated(sponge_implicit)) then
        deallocate(sponge_implicit)
    end if
    if (allocated(first_order_hydro)) then
        deallocate(first_order_hydro)
    end if
    if (allocated(xl_ext_bc_type)) then
        deallocate(xl_ext_bc_type)
    end if
    if (allocated(xr_ext_bc_type)) then
        deallocate(xr_ext_bc_type)
    end if
    if (allocated(yl_ext_bc_type)) then
        deallocate(yl_ext_bc_type)
    end if
    if (allocated(yr_ext_bc_type)) then
        deallocate(yr_ext_bc_type)
    end if
    if (allocated(zl_ext_bc_type)) then
        deallocate(zl_ext_bc_type)
    end if
    if (allocated(zr_ext_bc_type)) then
        deallocate(zr_ext_bc_type)
    end if
    if (allocated(hse_zero_vels)) then
        deallocate(hse_zero_vels)
    end if
    if (allocated(hse_interp_temp)) then
        deallocate(hse_interp_temp)
    end if
    if (allocated(hse_reflect_vels)) then
        deallocate(hse_reflect_vels)
    end if
    if (allocated(sdc_order)) then
        deallocate(sdc_order)
    end if
    if (allocated(sdc_quadrature)) then
        deallocate(sdc_quadrature)
    end if
    if (allocated(sdc_extra)) then
        deallocate(sdc_extra)
    end if
    if (allocated(sdc_solver)) then
        deallocate(sdc_solver)
    end if
    if (allocated(sdc_solver_tol_dens)) then
        deallocate(sdc_solver_tol_dens)
    end if
    if (allocated(sdc_solver_tol_spec)) then
        deallocate(sdc_solver_tol_spec)
    end if
    if (allocated(sdc_solver_tol_ener)) then
        deallocate(sdc_solver_tol_ener)
    end if
    if (allocated(sdc_solver_atol)) then
        deallocate(sdc_solver_atol)
    end if
    if (allocated(sdc_solver_relax_factor)) then
        deallocate(sdc_solver_relax_factor)
    end if
    if (allocated(sdc_solve_for_rhoe)) then
        deallocate(sdc_solve_for_rhoe)
    end if
    if (allocated(sdc_use_analytic_jac)) then
        deallocate(sdc_use_analytic_jac)
    end if
    if (allocated(cfl)) then
        deallocate(cfl)
    end if
    if (allocated(dtnuc_e)) then
        deallocate(dtnuc_e)
    end if
    if (allocated(dtnuc_X)) then
        deallocate(dtnuc_X)
    end if
    if (allocated(dtnuc_X_threshold)) then
        deallocate(dtnuc_X_threshold)
    end if
    if (allocated(do_react)) then
        deallocate(do_react)
    end if
    if (allocated(react_T_min)) then
        deallocate(react_T_min)
    end if
    if (allocated(react_T_max)) then
        deallocate(react_T_max)
    end if
    if (allocated(react_rho_min)) then
        deallocate(react_rho_min)
    end if
    if (allocated(react_rho_max)) then
        deallocate(react_rho_max)
    end if
    if (allocated(disable_shock_burning)) then
        deallocate(disable_shock_burning)
    end if
    if (allocated(T_guess)) then
        deallocate(T_guess)
    end if
    if (allocated(diffuse_temp)) then
        deallocate(diffuse_temp)
    end if
    if (allocated(diffuse_cutoff_density)) then
        deallocate(diffuse_cutoff_density)
    end if
    if (allocated(diffuse_cutoff_density_hi)) then
        deallocate(diffuse_cutoff_density_hi)
    end if
    if (allocated(diffuse_cond_scale_fac)) then
        deallocate(diffuse_cond_scale_fac)
    end if
    if (allocated(do_grav)) then
        deallocate(do_grav)
    end if
    if (allocated(grav_source_type)) then
        deallocate(grav_source_type)
    end if
    if (allocated(do_rotation)) then
        deallocate(do_rotation)
    end if
    if (allocated(rot_period)) then
        deallocate(rot_period)
    end if
    if (allocated(rot_period_dot)) then
        deallocate(rot_period_dot)
    end if
    if (allocated(rotation_include_centrifugal)) then
        deallocate(rotation_include_centrifugal)
    end if
    if (allocated(rotation_include_coriolis)) then
        deallocate(rotation_include_coriolis)
    end if
    if (allocated(rotation_include_domegadt)) then
        deallocate(rotation_include_domegadt)
    end if
    if (allocated(state_in_rotating_frame)) then
        deallocate(state_in_rotating_frame)
    end if
    if (allocated(rot_source_type)) then
        deallocate(rot_source_type)
    end if
    if (allocated(implicit_rotation_update)) then
        deallocate(implicit_rotation_update)
    end if
    if (allocated(rot_axis)) then
        deallocate(rot_axis)
    end if
    if (allocated(use_point_mass)) then
        deallocate(use_point_mass)
    end if
    if (allocated(point_mass)) then
        deallocate(point_mass)
    end if
    if (allocated(point_mass_fix_solution)) then
        deallocate(point_mass_fix_solution)
    end if
    if (allocated(do_acc)) then
        deallocate(do_acc)
    end if
    if (allocated(grown_factor)) then
        deallocate(grown_factor)
    end if
    if (allocated(track_grid_losses)) then
        deallocate(track_grid_losses)
    end if
    if (allocated(gravity_type)) then
        deallocate(gravity_type)
    end if
    if (allocated(const_grav)) then
        deallocate(const_grav)
    end if
    if (allocated(get_g_from_phi)) then
        deallocate(get_g_from_phi)
    end if



  end subroutine ca_finalize_meth_params


#ifdef RADIATION
  subroutine ca_init_radhydro_pars(fsp_type_in, do_is_in, com_in,fppt) &
       bind(C, name="ca_init_radhydro_pars")

    use rad_params_module, only : ngroups

    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: fsp_type_in, do_is_in, com_in
    real(rt)        , intent(in) :: fppt


    if (ngroups .eq. 1) then
       fspace_type = 1
    else
       fspace_type = fsp_type_in
    end if

#ifndef AMREX_USE_GPU
    if (fsp_type_in .ne. 1 .and. fsp_type_in .ne. 2) then
       print *, "fspace_type = ", fspace_type
       call castro_error("Unknown fspace_type")
    end if
#endif

    do_inelastic_scattering = (do_is_in .ne. 0)

    if (com_in .eq. 1) then
       comoving = .true.
    else if (com_in .eq. 0) then
       comoving = .false.
    else
#ifndef AMREX_USE_GPU
       call castro_error("Wrong value for comoving")
#endif
    end if

    flatten_pp_threshold = fppt

    !$acc update &
    !$acc device(QRAD, QRADHI, QPTOT, QREITOT) &
    !$acc device(fspace_type) &
    !$acc device(do_inelastic_scattering) &
    !$acc device(comoving)
    !$acc device(flatten_pp_threshold = -1.e0_rt)

  end subroutine ca_init_radhydro_pars
#endif

end module meth_params_module
