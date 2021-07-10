# 21.07

   * The sponge is now applied in a fully implicit manner at the end of
     the CTU advance, rather than using a predictor-corrector approach
     with time centering. This is more consistent with the original form
     of the sponge in Castro. (#1876)

   * Castro can now validate the runtime parameters set in the inputs
     file or on the commandline by setting
     castro.abort_on_invalid_params=1 (#1882)

# 21.06

   * Starting with this release, problem setups written in Fortran are
     no longer supported and will no longer work. Please consult the
     code documentation and example problem setups in Exec/ to understand
     the new problem setup format. If you need help converting a Fortran
     setup to C++, please file an issue. (#1728, #1732)

   * Sponge parameters are now only accepted through the inputs file; the
     &sponge namelist in the probin file is no longer read. (#1731)

   * Ambient parameters are now only accepted through the inputs file; the
     &ambient namelist in the probin file is no longer read. (#1742)

   * The update_sponge_params hook has been removed. (#1716)

   * The Fortran problem-specific source file, ext_src_nd.F90, has been
     removed. Problem-specific sources should be implemented in C++ in
     problem_source.H. (#1856)

   * Support for the legacy tagging scheme based on probin parameters (denerr,
     tempgrad, etc.) has been removed. These can be replaced with equivalent
     tagging criteria constructed in the inputs file; see the docs or examples
     in Exec/ to see how to use `amr.refinement_indicators`. (#1834)

   * The Fortran set_problem_tags hook has been removed. The C++ replacement
     is `problem_tagging()` in `problem_tagging.H`. (#1828)

   * The PrescribedGrav functionality has been removed (not replaced with a C++
     implementation). If you want to obtain the same functionality, you can use
     a problem-defined source term (look for problem_source in the documentation)
     and make the appropriate modification for applying it directly to the state
     (e.g. the momentum source term is rho * g). (#1854)

   * The custom radiation boundary using lo_bcflag and hi_bcflag coupled with
     an implementation of rbndry has been removed. (#1743)

   * We no longer store Reactions_Type in checkpoint files.  This means
     that newer versions of Castro will not restart from old version.
     
# 21.05

   * The parameter use_eos_in_riemann was removed -- we found no
     instances of it being used (#1623)

   * The option castro.apply_sources_consecutively was removed (#1636)

# 21.04

   * For simplified-SDC, we now correctly store only the reactive
     part of the update for rho_enuc, enuc, and rho_omegadot
     plotfile variables (#1602)

# 21.02

   * In axisymmetric geometry, there are additional forces that arise
     due to the changing direction of the unit vectors in the div{rho
     U U} term. The paper by Bernand-Champmartin discusses this. See
     issue #913. This adds those forces.  Note that the Coriolis force
     in 2-d axisymmetry is already consistent with a right-handed
     system despite our internal ordering of the state was r, z,
     theta.  (#923)

   * We can now set any of the Microphysics runtime parameters in the
     inputs file instead of probin.  Each group of parameters has a
     namesapce for the inputs file when set this way
     (e.g. eos.use_coulomb = 1), and the C++ inputs value will take
     precedence over the value set in probin if it is set in both
     places. (#1527)

# 21.01

   * grav_source_type and rot_source_type have been removed. (#580)

   * The minimum C++ standard supported by Castro is now C++17. Most modern compilers
     support C++17; the notable exception is RHEL 7 and its derivatives like CentOS 7,
     where the default compiler is gcc 4.8. In that case a newer compiler must be loaded,
     particularly a version of gcc >= 7.0, for example by installing devtoolset-7 or (if
     running on an HPC cluster that provides modules) using a more recent gcc module. (#1506)

   * There can now be multiple _prob_params files throughout the source
     tree.  We read the problem's file last and that takes precedence over
     any other _prob_params files found. (#1500)

   * The timestep limiter dtnuc_T has been removed. dtnuc_e and dtnuc_X
     are still available for controlling the burning timestep. (#1501)

   * A bug was fixed in the 2nd order true SDC (with reactions) that
     was giving the wrong solution and convergence (#1494).  A second
     bug was fixed in defining the weights for the Radau quadrature
     when using true SDC (#1493)

   * Compiling with the PGI compiler is no longer a requirement for the CUDA build of Castro.
     We recommend using COMP=gnu with a version of gcc that is C++17 compliant (gcc >= 7).

# 20.12

   * An issue with incorrect application of HSE boundary conditions on derived quantities
     is now resolved (#1356). Also, at this point the old Fortran implementations hypfill,
     denfill, ext_hypfill, and ext_denfill have been removed; problem-specific boundary
     conditions should be implemented using the new C++ interface in this release from #1289.

   * The minimum supported Hypre version is now 2.19.0. (#1333)

   * We have switched from a Fortran to a C++ implementation of VODE in Microphysics.
     As a result we have also switched the Strang and simplified SDC burners in Castro
     to use this C++ implementation. Most networks used in Castro have already been
     ported to C++. While networks are not required to have a C++ implementation,
     networks implemented only in Fortran  will not be useable on GPUs, and eventually
     we will use C++ only. (#1313)

   * `problem_checkpoint` and `problem_restart` are moved to C++ from Fortran. See
     Exec/science/wdmerger for an example of the new scheme. `Problem.f90` and `Problem_F.H`
     are now deleted from the code; if you were using these to implement problem-specific
     functionality, you can still manually add these files to the `Make.package` for your
     problem setup. (#1311)

   * For setups using Poisson gravity, tagging is now turned off in locations where
     the fine levels would have been adjacent to a physical boundary. (This previously
     led to an abort.) (#1302)

   * An interface for doing problem tagging in C++ has been added. (#1289)

   * Simplified SDC now only supports the C++ integrators (#1294)

   * MHD problems can now do the magnetic field initialization in C++
     (#1298)

# 20.11

   * The minimum C++ standard supported by Castro is now C++14. Most modern compilers
     support C++14; the notable exception is RHEL 7 and its derivatives like CentOS 7,
     where the default compiler is gcc 4.8. In that case a newer compiler must be loaded,
     particularly a version of gcc >= 5.0, for example by installing devtoolset-7 or (if
     running on an HPC cluster that provides modules) using a more recent gcc module. (#1284)

   * A new option, `castro.retry_small_density_cutoff`, has been added. In some
     cases a small or negative density retry may be triggered on an update that
     moves a zone already close to small_dens just below it. This is not uncommon
     for "ambient"/"fluff" material outside a star. Since these zones are not
     dynamically important anyway, triggering a retry is unnecessary (and possibly
     counterproductive, since it may require a very small timestep to avoid). By
     setting this cutoff value appropriately, the retry will be skipped if the
     density of the zone prior to the update was below the cutoff. (#1273)

# 20.10

   * A new refinement scheme using the inputs file rather than the Fortran
     tagging namelist has been added. (#1243, #1246) As an example, consider:

     ```
     amr.refinement_indicators = dens temp

     amr.refine.dens.max_level = 1
     amr.refine.dens.value_greater = 2.0
     amr.refine.dens.field_name = density

     amr.refine.temp.max_level = 2
     amr.refine.temp.value_less = 1.0
     amr.refine.temp.field_name = Temp
     ```

     `amr.refinement_indicators` is a list of user-defined names for refinement
     schemes. For each defined name, amr.refine.<name> accepts predefined fields
     describing when to tag. In the current implementation, these are `max_level`
     (maximum level to refine to), `start_time` (when to start tagging), `end_time`
     (when to stop tagging), `value_greater` (value above which we refine),
     `value_less` (value below which to refine), `gradient` (absolute value of the
     difference between adjacent cells above which we refine), and `field_name`
     (name of the string defining the field in the code). If a refinement indicator
     is added, either `value_greater`, `value_less`, or `gradient` must be provided.

   * Automatic problem parameter configuration is now available to every
     problem by placing a _prob_params file in your problem directory.
     Examples can be found in most of the problems in Castro/Exec, and you
     can look at the "Setting Up Your Own Problem" section of the documentation
     for more information. This functionality is optional, however note that
     a file containing a Fortran module named "probdata_module" is now
     automatically generated, so if you have a legacy probdata.F90 file
     containing a module with that name it should be renamed. (#1210)

   * The variable "center" (in the `problem` namespace) is now part of this
     automatically generated probdata module; at the present time, the only
     valid way to change the problem center to a value other than zero is in
     amrex_probinit(). (#1222)

   * Initialization of these problem parameters is now done automatically for
     you, so a call to probdata_init() is no longer required in amrex_probinit(). (#1226)

   * Problems may now be initialized in C++ instead of Fortran. Instead of implementing
     amrex_probinit() and ca_initdata(), the problem should implement the analogous
     functions initialize_problem() and initialize_problem_state_data(). If you switch to
     the new C++ initialization, be sure to delete your Prob_nd.F90 file. By default both
     implementations do nothing, so you can pick either one but do not pick both. (#1227)

   * The external heat source term routines have been ported to C++
     (#1191).  Any problem using an external heat source should look
     convert the code over to C++.

   * The interpolate_nd.F90 file has been moved to Util/interpolate and
     is only compiled into Castro if you set USE_INTERPOLATE=TRUE

# 20.09

   * Reactions now work with MHD (#1179)

   * MHD now uses the main slope routine (#1058) The order of the
     slope is now controlled by plm_iorder, just as with hydro.  There
     is an additional option, mhd_limit_characteristic, that
     determines if the limiting is done on the primitive or
     characteristic variables (the default).

# 20.08

   * Rotation_Type has been removed from StateData. (#1128)

   * castro.use_post_step_regrid now unconditionally regrids after
     every timestep on every level. (#898)

   * An issue with gravity.max_solve_level resulting in accesses to invalid data
     (#469, #1118) has been resolved. (#1123)

   * If castro.speed_limit is set to a number greater than zero, this
     will now be strictly enforced on the magnitude of the velocity. (#1115)

   * When using AMR and gravity or rotation, the source terms applied after
     a reflux would have been incorrect if the previous timestep had a retry
     (#1020). This has now been fixed. (#1112)

   * We now have the ability to access the problem-specific runtime
     parameters in C++ (#1093)

# 20.07

   * The master branch has been renamed the main branch. If you have an
     existing clone of Castro, then do the following to update for this
     change. First, do `git checkout master` if you're not already on the
     old master branch. Then do `git pull`. This will gather the updates
     to the repo, but will fail with the message `Your configuration specifies
     to merge with the ref 'refs/heads/master' from the remote, but no such ref
     was fetched.` Then you can simply do `git checkout main` and your local
     repo should automatically switch to that branch and track updates from
     the upstream repo on GitHub. If you like, you can then delete the old
     master branch with `git branch -D master`.

   * The CUDA build no longer has a requirement that amr.blocking_factor
     be a multiple of 8. Though this is recommended for performance reasons,
     it was previously required due to correctness reasons because of the
     use of an AMReX Fortran function, amrex_filccn. As noted in #1048, this
     function is no longer required due to recent changes in Castro (problems
     overriding bc_fill_nd.F90 or bc_ext_fill_nd.F90 do not need to provide an
     initial fill of the ghost zone data before implementing their specific
     boundary conditions; this is now done for you). Calling this function
     may now result in race conditions and correctness issues in the CUDA
     build, so it should be removed from any problem setups. (#1049)

   * The functionality that permitted the rotation rate to change as a
     function of time, castro.rotation_include_domegadt and
     castro.rotational_dPdt, has been removed. (#1045)

   * A CUDA illegal memory access error in Poisson gravity and diffusion
     has been fixed (#1039).

   * The parameter castro.track_grid_losses has been removed. (#1035)

   * The parameter castro.print_fortran_warnings, which no longer had any
     effect, has been removed. (#1036)

   * PPM reconstruction has been added to the MHD solver (#1002)

   * The Reactions_Type StateData has been reworked so that its first
     NumSpec components are rho * omegadot rather than omegadot; then,
     the NumAux auxiliary components are stored, if the network has any
     auxiliary variables; then, rho * enuc is stored (enuc itself is
     removed), and finally the burn weights are stored. The checkpoint
     version has been incremented, so this version of the code cannot
     restart from checkpoints generated with earlier versions of the
     code. (#927)

   * A bug where refluxing between AMR levels resulted in incorrect results
     when a retry occurred in the previous timestep has been fixed. (#1018)

# 20.06

   * The parameter castro.density_reset_method has been removed. A density
     reset now unconditionally sets the density to small_dens, the temperature
     to small_temp, and zeros out the velocities. (#989)

   * A constrained-transport corner transport upwind MHD solver has been
     added.  This can be used by compiling with USE_MPI = TRUE.  Presently
     it only works for a single level (no AMR).  (#307)

   * A burning timestep limiter dtnuc_T has been added which restricts the
     burning from updating the temperature by more than the factor
     dtnuc_T * T / dT/dt. (#972)

   * The reaction weights metric implemented in version 20.05 (#863) has been
     added to the simplified SDC reactions driver. (#930)

   * When using the simplified SDC integration scheme, we now save new-time
     Reactions_Type data to plotfiles. (#929)

# 20.05

   * The parameter use_custom_knapsack_weights and its associated
     functionality have been removed. (#877)

   * We've changed how the runtime parameters are stored.  Previously
     they were static members of their respective class, but this
     prevented their use in lambda-capture functions on GPUs.  Now the
     runtime parameters are grouped into namespaces as extern managed
     data. (#873)

   * We currently have a scheme for storing reactions weightings, which
     are a measure of the number of RHS evaluations during the burn and
     therefore a proxy for the difficulty of the burn. These weights were
     added as separate StateData depending on the runtime option
     use_custom_knapsack_weights. Now, instead we place the weights
     directly in the Reactions_Type StateData as a new component.

     The number of ghost zones in Reactions_Type is increased to 4.

     The checkpoint version has now been incremented; this version of the
     code will not be able to restart from a checkpoint generated by earlier
     versions of the code. (#863)

   * The meaning of dt_cutoff has changed: it is now the fraction of the
     current simulation time which dt may be no smaller than, instead of
     being an absolute measure. We now have set a non-zero default
     (1.e-12) as well. (#865)

   * Backwards compatibility in restarting from a checkpoint is no longer
     supported. Checkpoints from older versions of the code (as determined
     by the checkpoint version in the CastroHeader file in the checkpoint
     directory) cannot be restarted from. (#860)

   * Added an option to do CTU reactions in C++.  A compile flag
     USE_CXX_REACTIONS is added which switches to the C++ integrator
     in Microphysics. Since we will be doing a phased implementation
     of the networks in Microphysics, this is opt-in for now.  (#836)

   * More of the core routines have been ported to C++, including the
     hydro and diffusion timestep estimators (#853) and the sponge
     (#857)

   * AMReX provides CpuBndryFuncFab and GpuBndryFuncFab which are very
     similar to what generic_fill and hypfill did. The AMReX
     implementations are now used. We still have a hypfill and denfill
     function, so that existing problems are not broken, but the main
     one in Source/ no longer calls amrex_filcc (it only has the
     ambient code now). The problems that do override bc_fill_nd.F90
     are thus no longer required to call amrex_filcc. (#837)

   * We now always issue a timestep retry if the density after an
     advance is negative (or less than small_dens). The parameter
     castro.retry_neg_dens_factor is removed. The parameter
     castro.retry_tolerance is also removed as it no longer has
     any effect. (#796)

   * The timestep control parameter castro.change_max now also will
     prevent the timestep by shrinking too much in one timestep
     (previously it would only prevent it from growing too much).
     If change_max is violated in a timestep we will do a retry
     to take more graceful steps. (#844)

   * We now check if the problem setup initialized the density or
     temperature to a value near small_dens or small_temp and abort.
     If this happens, the recourse is to adjust small_dens and
     small_temp to a meaningful value for your problem.  (#822)

   * The src_q multifab was removed and instead we convert the
     conserved state sources to primitive state sources FAB by FAB.
     This saves a lot of memory at the expense of an EOS call. (#829)

   * The plm_well_balanced option was removed.  It was essentially the
     same as use_pslope except it was lower order and only worked with
     constant gravity.  use_pslope now works with both CTU+PLM and
     SDC2+PLM.  A new test problem, hse_convergence, was added to look
     at the behavior of the different reconstruction methods with HSE.

# 20.04

   * A potential undefined flux from the HLL solver when using
     hybrid_riemann has been fixed (#823)

   * The parameter castro.allow_small_energy has been removed. The
     code behavior is now similar to what it would have been with
     allow_small_energy == 0 (the internal energy can never be
     smaller than that allowed by small_temp). (#817)

   * The BC interfaces have been merged and converted to a new FAB
     interface as part of the port to C++. (#819)

   * All boundary fill interfaces other than hypfill and denfill have
     been removed. So, we no longer support overriding the boundary
     conditions for data other than State_Type. Radiation still has
     its own set of custom boundary conditions that can be accessed
     through the inputs file, as described in the docs. (#815)

   * The conversion of the CTU hydrodynamics code to C++ continues.
     The Riemann solvers were converted to C++ (#801) and the
     hybrid momentum routines (#805), the PLM reconstruction (#814),
     the conversion of primitive to conserved variables (#804)

   * We've changed how the backup for retries is done.  Presently if
     use_retry is enabled we make a pre-emptive copy of the StateData
     right at the beginning of the timestep.  Now we only backup when
     we detect that a retry is needed (#812)

# 20.03

   * We now depend on the fundamental constants from Microphysics
     instead of keep our own copy in Castro (#787)

   * We removed the ppm_predict_gammae option for the CTU hydro solver.
     This was not used frequently and did not show much difference with
     the default (rho e) reconstruction. (#780)

   * The Microphysics "extern" parameters are now available in C++

   * We've started converting the CTU hydro solver from Fortran to C++
     (#731).  The PPM reconstruction is now done in C++ (#784).

   * The option ppm_temp_fix = 3 was removed.  This used a
     temperature-based eigensystem for characteristic tracing but was
     never used for production science.

   * If a derived variable has multiple components, all components are now
     added to plotfiles. Previously only the first component was used. (#758)

   * We have updated our workflow when it comes to Castro's dependencies.

     Previously Castro shipped with it a minimal set of microphysics that
     allowed basic problem setups like Sedov to compile, and more advanced
     setups (like ones that include nuclear burning) required downloading
     the starkiller-astro Microphysics repository as an additional step.
     Now, that Microphysics repository is a requirement for using Castro.
     If you are a current user of the Microphysics repository and prefer
     the current workflow where you maintain Microphysics as a separate
     installation from Castro, no change in your workflow is necessary:
     if MICROPHYSICS_HOME is set as an environment variable, Castro will
     use the Microphysics installation in that directory. However we have
     also added Microphysics as a git submodule to Castro, which is now
     the required path if you previously were not using the more advanced
     microphysics (but is also a possibility for those previously using a
     standalone Microphysics installation). To obtain this, you can use
     git submodule update --init --recursive from the top-level directory
     of Castro. The developer team ensures that the version of Microphysics
     that you obtain this way is consistent with the current version of Castro.
     Then, you can keep up to date with the code mostly as normal, except now
     using git pull --recurse-submodules instead of git pull.

     Similarly, AMReX is now maintained as a git submodule rather than as an
     external standalone installation. If you use the same git submodule command
     as above, you'll obtain AMReX. As with Microphysics, you may opt to
     rely on your own installation of AMReX by setting the AMREX_HOME
     environment variable. However you are then responsible for keeping it
     in sync with Castro; if you use the submodule, then you'll get the version
     of AMReX that we have tested to ensure compatibility with the current
     version of Castro. (#651, #760, #762, #765)

   * The names of the conserved state variables in C++ (Density, Xmom, etc.)
     have been changed to match the names in Fortran (URHO, UMX, etc.).
     For user code, this will only affect problem-specific setup code
     like Prob.cpp that references specific state variables. For compatibility,
     we have kept a copy of the old names around that redirect to the
     new names, but the old names are now considered deprecated and will
     be removed in a future release. (#757)

# 20.02

   * Fixed a bug in the nuclear burning timestep estimator when on GPUs
     (#745)

   * rewrote the 4th order SDC hydro driver in C++ to allow code reuse
     with other solvers (#742), and simplified the 2nd order SDC code
     to do dimensional sweeps to reduce memory (#749)

   * The option radiation.integrate_planck has been removed; it was only
     used by one test. By default we always do the full integral of the
     Planck function. (#740)

   * Most of the radiation test problems have been moved over to a new
     opacity directory, rad_power_law, and all of the parameters that
     controlled the behavior of the power law opacity have been moved
     to the extern probin module. We now always expect you to pick a
     specific opacity implementation, so the parameter
     radiation.use_opacity_table_module has been removed. The "null"
     opacity implementation has been previously moved, and the code
     will fail to compile if you attempt to use it; you will need to
     update to rad_power_law. (See the documentation for information
     about how to use this new implementation.)

     Additionally, the code for the multigroup solver was effectively
     previously setting the Rosseland opacity, kappa_r, equal to the
     Planck opacity, kappa_p, if the latter was set but the former was
     not. There was similar unintuitive behavior for the behavior of
     the scattering parameter. Now you will get exactly what you ask
     for in the probin file, given the defaults in the _parameters file
     for the rad_power_law opacity. By default the constant coefficients
     for both are negative, which is invalid, so both must be set to a
     non-negative value for the code to work. Problems that were previously
     setting const_kappa_p but not const_kappa_r should set the latter
     equal to the former to maintain the same code behavior. The analogous
     thing should be done for the exponents (kappa_p_exp_m, kappa_p_exp_n,
     and kappa_p_exp_p). (#725)

   * The parameter radiation.do_real_eos = 0 has been removed, and its
     functionality is now enabled with a new equation of state called
     rad_power_law. This new EOS is only compatible with the pure
     radiation-diffusion tests, not with castro.do_hydro = 1. (#722)

   * We now default to use_retry = 1, instructing Castro to retry a
     step with a smaller dt if there is a CFL violation, burning
     failure, or negative timestep.  For the burning failure, we have
     Castro set the Microphysics parameter abort_on_failure to .false.
     at a high priority (so it overrides the Microphysics default).
     We also check to make sure the combination of parameters makes
     sense at runtime. (#724)

   * The parameter castro.hard_cfl_limit has been removed. (#723)

   * Some unnecessary clean_state calls were removed (#721)

   * Support for neutrino radiation diffusion has been removed.

   * A bug was fixed in the hydro CFL timestep estimator for
     simplified-SDC.  The timestep was more restrictive than it needed
     to be. (#727)

   * A bug was fixed in the simplified-SDC nuclear burning timestep
     estimator (#733)

# 20.01

   * A new option castro.limit_fluxes_on_large_vel has been added. It
     is similar to the existing option limit_fluxes_on_small_dens --
     fluxes are limited to prevent the velocity in any zone from
     getting too high. The largest legal speed is set by
     castro.speed_limit. (#712) This is more general than the previous
     solution proposed by castro.riemann_speed_limit, so that
     parameter has been removed. (#714)

   * The AMR parameter amr.compute_new_dt_on_regrid is now on by
     default. This avoids crashes that result from the CFL number
     being too large after regridding, because we update the
     timestep after seeing that larger velocity. You can still opt
     to set this off if you want to in your inputs file. (#720)

   * We have added calls into Hypre that only exist as of version
     2.15.0, so that is the new minimum requirement for Castro
     radiation. Note that Hypre is now hosted on GitHub at
     https://github.com/hypre-space/hypre.

   * A new option castro.limit_fluxes_on_large_vel has been added. It
     is similar to the existing option limit_fluxes_on_small_dens --
     fluxes are limited to prevent the velocity in any zone from
     getting too high. The largest legal speed is set by
     castro.riemann_speed_limit. (#712)

   * A new option castro.apply_sources_consecutively has been
     added. By default we add all source terms together at once. This
     option, if enabled, adds the sources one at a time, so that each
     source sees the effect of the previously added sources. This can
     matter, as an example, for the sponge source term, which may be
     more effective if it is added after source terms such as gravity
     that update the velocity. (#710)

   * A new option castro.ext_src_implicit has been added. The external
     source terms were previously only implemented as an explicit
     predictor-corrector scheme. The new option, if turned on, changes
     the handling of the external source terms to allow an implicit
     solve. This is done by subtracting the full old-time source and
     adding the full new-time source in the corrector, rather than
     -0.5 and +0.5 of each, respectively. It is still up to the
     individual problem to make sure it is consistent with this scheme
     if the option is turned on. (#709)

   * Add option for using monopole BCs in 3D.  By setting
     gravity.max_multipole_order to a negative number, you can use
     monopole gravity to fill the boundary conditions, rather than the
     multiple BCs. This is useful for debugging purposes.  To make the
     behavior consistent, we now use multipole BCs by default in 2D as
     well. (#716)


# 19.12

   * The use_retry mechanism has been enabled for the simplified
     SDC time integration method. (#695)

   * A case where use_retry could result in a very small last
     subcycle has been avoided. (#701)

   * We no longer allocate memory for sources for the species
     in the conserved state unless PRIM_SPECIES_HAVE_SOURCES is set
     (#699)

   * A subroutine eos_on_host has been added to the EOS module.
     This is a wrapper for the EOS that must be used for CUDA
     builds if the EOS is being called in probinit or other
     places that don't run on the GPU. (#693)

   * We now use VODE90 instead of VODE by default. (#677)

   * A new unit test was added, model_burner, which reads in a 1-d
     initial model and calls the reaction network on it.  This can
     be used to test tolerances, etc.

# 19.11

   * The density flux limiter was simplified and fixes a race condition
     (#646)

   * The SDC algorithm can now use Radau quadrature instead of
     Gauss-Lobatto quadrature. (#666)

   * The option castro.ppm_reference_eigenvectors has been removed.  This
     is now used by default with the CTU PPM solver.

# 19.10

   * The SDC algorithm now implements the burning conditionals
     depending on rho and T (react_rho_min, react_rho_max,
     react_T_min, react_T_max) (#598, #654)

   * The SDC/MOL PLM reconstruction now implements reflecting BCs on
     the interface states (#652, #654)

   * A well-balanced scheme has been added to the piecewise linear SDC
     method, enabled with castro.plm_well_balanaced=1.  At the moment
     it only supports constant gravity. (#294, $654))

   * The weighting of the time-node fluxes stored in the flux registers
     for SDC has been fixed (#654, #658)

   * As before, we can choose the reconstruction with PLM using the
     castro.plm_iorder flag: 1 = piecewise constant, 2 = piecewise
     linear slopes.  Now we added a way to specify the limiter used
     with the linear slopes.  castro.plm_limiter = 1 will use the 2nd
     order MC limiter and castro.plm_limiter = 2 will use the default
     4th order MC limiter (previously there was no way to select the
     2nd order limiter). (#654)

   * The Runge-Kutta based method-of-lines integration method has
     been removed in favor of the SDC integration. (#657)

   * A new way of specifying the problem runtime parameters has been
     introduced.  They can now be specified in a plain text file,
     _prob_params, and at compile time, the probdata_module is
     automatically created.  This automates the creation of the
     probdata variables, the namelist for reading them, setting them
     as managed for CUDA, and adds the ability to output the values to
     a file (like job_info).  This feature is opt-in. You need to set
     USE_PROB_PARAMS in your GNUmakefile and then define the problem
     parameters in a file _prob_params in the problem directory.
     (#234, #619, #673)

   * The time to output is now stored in the job_info file (#365)

   * The SDC time advancement method has been documented

   * The job_info file now reports the number of GPUs being used.

# 19.09

   * You can now type ./Castro.gnu.ex  *describe to see the list of
     modules / compilers the code was built with (#660)

   * The reaction quantities are now computed as proper 4th order
     averages for the plotfile, when using sdc_order = 4 (#647)

   * The velerr tagging now takes the abs() of the velocity component
     to ensure we tag large positive and negative velocities.


# 19.08.1

   * Fix CUDA compilation

   * Remove special treatment of 4th order outflow BCs (see #648)

# 19.08

   * We slightly changed how the characteristic tracing is done for
     the CTU PPM hydro solver  * we now use the limit of the parabola
     as the edge state if the wave is not moving toward the interface
     (#632)

   * The CTU PPM solver now uses a lot less memory by computing the
     integrals under the parabolas as needed instead of precomputing
     and storing them (#624)

   * We created a new error wrapper, castro_error(), to replace the
     AMReX amrex_error().  This will allow us to deal with error when
     on the GPU.

   * The new SDC solver has had substantial improvements:

      * Explicit thermal diffusion is now implemented for both 2nd and
        4th order accurate solvers. (#610)

      * There is a new option (castro.sdc_extra) for taking extra SDC
        iterations.

      * The Newton solver for the SDC update can now subcycle.

      * The sdc_solver_relax_factor option was fixed

      * There is now an absolute tolerance on species used in the
        error check for the SDC solve.

      * Some scalings of terms in the Jacobian were fixed.

      * We now use one-sided stencils for the reconstruction at
        physical boundaries (#633)

   * The flame problem setup now can initialize conservatively from a
     pre-existing flame solution.  The analysis routines have seen
     some improvements for working with 4th order accurate
     simulations.

   * the diffusion_test unit test now works for 4th order problems.

# 19.07

   * There is now a separate set of boundary filling routines for
     source terms, source_fill.F90.  Previously this was handled
     by generic_fill.F90 (#627)

   * The tagging routines have been reformulated so that they run on
     the GPU.  Since the tags_and_untags routine in AMReX is not
     GPU-accelerated, we opt to directly fill in the TagBoxArray in
     the tagging routines. We now pass the TagBox to Fortran as an
     int8_t. This means that the interface to problem_tagging_nd.F90
     has been updated to use integer(1).

     Castro_prob_err_list.H and other related files have been deleted
     as they were not actually used anywhere.

     Castro_error.cpp is now removed and there is no further support
     for writing custom tagging routines. The set of variables that we
     check for tagging is hard-coded in Castro and can be expanded as
     needed. Problem-specific tagging should be done through the
     set_problem_tags functionality. (#611)

   * The dimension-specific code for problem initialization, boundary
     conditions, and external heat terms has been removed, as warned
     in the previous release notice.

# 19.06

   * Deprecation notice: as of the 19.06 release, dimension-specific
     problem setups are deprecated. Presently they are opt-in by
     adding DIMENSION_AGNOSTIC = TRUE to your makefile, and using
     a Prob_nd.F90 file instead of a Prob_{1,2,3}d.F90 file. The
     dimension-agnostic Prob_nd.F90 is used to fill initial data
     for all dimensions. There is always a 3D loop over dimensions,
     and in 1D or 2D the unused dimensions have indices (lo, hi) =
     (0, 0) which is valid in Fortran. The current interface is found
     in Source/problems/Prob_nd.F90. Most of the problems have been
     converted to dimension-agnostic form and any remaining ones will
     be done shortly, so you can use e.g. the Sedov or DustCollapse
     problems to model you own problem on. The dimension agnostic
     problem setup also implies dimension agnostic helper routines
     such as in bc_fill_nd.F90  * any user-facing file ending in a
     1d/2d/3d.f90 file extension is deprecated. In the 19.07 release
     support for these will be removed and problem setups will only
     work if they are dimension agnostic. Please file an issue if you
     need assistance converting your problem.

   * Deprecation notice: as of the 19.06 release, problem-specific
     overrides of Castro_error.cpp, and in general custom tagging
     routines (including Castro_prob_err.cpp and associated files),
     are deprecated. The only supported mechanism for problem-specific
     tagging is through the set_problem_tags function in problem_tagging_nd.F90.
     (There are also dimension-specific versions of this file, but these
     are now deprecated as above.) Please file an issue if you need
     assistance converting your tagging setup to use the problem tagging,
     or if you need more data in that interface to be able to implement
     your tagging scheme. Support will be removed in the 19.07 release.

   * Deprecation notice: as of the 19.06 release, the problem_pre_tagging_hook
     and problem_post_tagging_hook are deprecated. These were not actually
     being used in any problem. These will be removed in the 19.07 release.

# 19.05

   * The dimension agnostic version of the external source term in
     ext_src_nd.F90 has been updated to use the ISO C binding interface,
     and two parameters, time and dt, are now passed by value  * see
     Source/sources/ext_src_nd.F90 for the new interface.

   * problem_derive_nd.f90 has been renamed to problem_derive_nd.F90.

   * The velocity calculated for the interface in the Riemann solve for
     the CGF/CG Riemann solvers can no longer exceed the speed of light.
     The parameter castro.riemann_speed_limit can be set to control the
     speed limit applied in the Riemann solver  * this is useful for
     preventing unphysically large velocities from being created at
     shock fronts or near large density gradients.

   * The algorithm for recalculating source terms after an AMR reflux
     did not set some data needed to correctly calculate positions on
     the fine levels, so position-dependent source terms (rotation,
     hybrid momentum) were being calculated incorrectly. This caused
     some strange effects in AMR simulations with rotation, such as
     binary star orbits getting wider with time and drifting relative
     to the system center of mass. This has now been fixed. (#599)

   * density_reset_method == 3 (which, in the event of a density reset,
     reset to the density at the beginning of the timestep) no longer
     exists. (#538)

   * The behavior of use_retry = 1 with retry_neg_dens_factor changed
     slightly since we now choose the retry timestep based on the difference
     between the (incorrect) negative density and the density it was
     reset to, rather than the old density at the beginning of the step.
     It still does a similar type of timestep limiting, but quantitatively
     the timesteps it chooses will be different. (#538)

   * A sign error was fixed in the hybrid_hydro angular momentum
     algorithm This addresses issue #462.  During commmit 0f09693, a
     change in signs was introduced in add_hybrid_momentum_sources,
     which should be analogous to linear_to_hybrid (#594)

# 19.04

   * The runtime parameter castro.fix_mass_flux has been removed: it is not
     clear what the use case is, and it had no test suite coverage. (#572)

   * Fixed a bug introduced in August 2015 that resulted in incorrect
     real bounds being passed to ca_initdata after a restart for problems
     using a grown domain. This would have resulted in incorrect initialization
     for problems using the grown restart capability if their initialization
     depended on the position on the grid. (#566)

   * Using point-mass gravity no longer requires USE_POINTMASS = TRUE
     in your makefile; USE_GRAV = TRUE is sufficient. However, to
     compensate for this, you must now include castro.use_point_mass = 1
     in your inputs file to enable the point mass. This input parameter
     already existed, but was defaulted to 1 since it only mattered
     if the compile flag was enabled. Now the default is 0.

   * Also, a couple bugs in the point-mass gravity have been fixed.
     The algorithm was not correct in 1D and 2D, and this has been
     resolved. And the point mass value was not being preserved across
     restarts, which is an issue if you're using point_mass_fix_solution
     to update the point mass value as mass accretes to the center of
     the domain. This has been fixed as well.

   * fixed a bug in the source term to (rho e) evolution when using
     MOL or the new SDC integration (#543, #545) and also no longer
     recompute the source terms after reflux for these methods (#549)

   * Dimension agnostic problem setups have had the interface to the
     physical boundary conditions changed (hypfill, denfill, etc.).
     If your problem is dimension agnostic, please consult the new
     interfaces in Source/problems/bc_fill_nd.F90 to understand how
     to convert your problem. The changes are that (1) the "ca_" prefixes
     have been removed from the subroutine names, (2) the "time" parameter
     is passed by value, and (3) the (lo, hi) indices that are the target
     region to update the boundaries on are explicitly passed in as
     the first two arguments. (#546)

   * removed the code to extrapolate diffusion terms to ghost cells
     as it is no longer necessary (#532)

   * we remove enthalpy diffusion (there were no known applications of
     this) and species and velocity diffusion (they were 1-d only).
     None of these routines were regularly tested.  (#534)

   * the problem diagnostics in Castro/Diagnostics have been converted to
     C++ to remain compatible with the AMReX build system.

# 19.03

   * Fixed a minor long-standing bug in the simplified SDC implementation
     involving incorrect indexing. This changes results slightly.

   * a number of tests involving reactions have been moved from
     hydro_tests to reacting_tests (#527)

   * The old spectral deferred corrections method has been renamed
     "simplified" SDC. It is accessed with time_integration_method = 3.
     This still requires building with USE_SDC = TRUE, and when building
     this way, the other time integration methods are unavailable.

   * The 4th order hydro was extended to support general equations
     of state (it is still single level only).  Artificial viscosity
     was also added.

   * A framework for a new spectral deferred corrections method was
     added.  This will allow for higher-order time integration. (#310)

   * Fix a bug where the self_heat parameter was not being initialized
     for the burning timestep limiter (#521).

   * By default, we no longer allocation storage for source terms to
     species in the primitive variable state.  This is set via the
     _variables file, parsed by set_variables.py.  To allow for
     species sources, you need to set PRIM_SPECIES_HAVE_SOURCES.  This
     is done currently for SDC. (#519)

   * renamed QVAR to NQSRC to make it clear that this is the number of
     source terms for the primitive variable state.  We also fixed a
     number of places where QVAR was used instead of NQ. (#517)

   * A new runtime parameter, T_guess, was created.  This is used as
     the initial temperature guess when calling the EOS with inputs
     other than rho, T, for the initial Newton iteration. (#509)

   * The CTU hydrodynamics driver has been rewritten in C++
     (#498)

   * The input parameter castro.do_ctu has been renamed
     castro.time_integration_method. The current legal values
     are 0 (CTU, the default) and 1 (MOL).

   * fixed a bug in the ppm_temp_fix = 1 reconstruction  * we were not
     doing the initial reconstruction of temperature

# 19.02

   * The flux limiter used with the options limit_fluxes_on_low_dens
     was not implemented correctly.  This has been fixed.  (PR #493)

   * the CTU hydro solver no longer does any allocation in any of the
     support routines  * it is all done in the top-level
     driver. (#455)

   * The CTU solver now makes explicit the range of cells looped over
     in the transverse routines (#456)

   * The plotfile quantities divu and magvort were fixed in
     axisymmetric coordinates and diff_term was fixed in all
     geometries.  (#446, 448, 449, 472)

   * abar is a new derived variable (#470)

   * the job_info file now stores domain information (#451), and the
     job_info files is also stored in checkpoints now too (#450)

   * we can now refine on enuc  * the nuclear energy generation rate.
     This is controlled by the parameters (in the &tagging probin
     namespace) enucerr, and max_enucerr_lev.  We also moved the dxnuc
     tagging parameters from inputs to probin, where they are now
     named dxnuc_min (formerly castro.dxnuc), dxnuc_max, and
     max_dxnuc_lev.  (#364, #437, #473)

   * The diffusion cutoff now is a linear ramp instead of a
     discontinous cutoff.  For densities less than
     diffuse_cutoff_density_hi, the transport coefficient is scaled
     linearly until the density reaches diffuse_cutoff_density, where
     it is zero.

# 19.01.4

   * fixed the .zenodo.json

# 19.01

   * The User's Guide is now automatically built from the development
     branch using travis.

   * we now store the job_info file in the checkpoints (#450)

   * we now automatically generate Doxygen docs along with the User's
     Guide and have started adding docstrings throughout the
     code. (#461)

   * The MG solver was optimized a bit (#464)

# 18.12

   * fixed a bug in the CUDA version of the MOL integrator  * it was
     setting the interface states incorrectly.

   * removed ppm_type = 2.  This was not used for science simulations
     because it was never shown to be completely stable.  In the near
     future, the full fourth order method will be merged which will be
     a better replacement.

   * a bug was fixed in the 1-d SDC integration  * we were not
     applying the reactive source terms to the hydrodynamics interface
     prediction.

   * we now apply the source terms for PLM before the transverse
     routines.  This is more consistent with the PPM version.

   * the angular momentum in the plotfile is not computed with respect
     to the center array initialized in the probinit

   * fixed a bug in 2-d with rotation  * we were adding the source
     terms to the out-of-plane velocity twice in the prediction of the
     interface states, resulting in an overall first-order
     approximation there.

   * fixed an issue that resulted in the custom knapsack distribution
     map not generating a useful distribution map for the reactions.
     Also, fixed an edge case where this custom distribution map could
     cause a numerical overflow and code crash. (#382)

   * rewrote the reconstruction routines to eliminate temporary arrays
     to make them GPU friendly

   * the documentation was converted to sphinx

   * changed the interface to the conductivity routines so the conductivity
     is not part of eos_t (#431)


# 18.11

   * we've restructured the CTU solver to no longer use a slab
     approach in 3d.  This is in preparation for offloading the solver
     to GPUs.

   * a few minor bugs were fixed in the transverse routines of the CTU
     solver with radiation

   * simplified conductivity interface to match the eos interface by
     moving the conductivity variable in the arguments into the eos
     type.


# 18.10

   * fixed handling of external BCs (#402)

   * unified some CUDA hydro solver code with the method-of-lines code

   * offloaded gravity source terms, rotation, reactions, and sponge
     to GPUs with CUDA

   * merged the different dimensional versions of the CTU consup
     routine into a single routine (#399)

   * removed the unsafe option "allow_negative_energy"


# 18.09

   * we now only trace under sources that are non-zero, to save
     computational expense.  (#381)

   * we now update T for consistency when we reset a small internal
     energy, e (#384)

   * The parameter dual_energy_update_E_from_e has been removed,
     and the default behavior is now that the total energy (E)
     will not be changed when we reset the internal energy (e).
     This will cause changes in simulation output. (#368)

   * The probin parameter eos_input_is_constant is now true by
     default. This means that when calling the EOS in the mode
     eos_input_re, the energy will not be updated after the EOS
     call (e.g. by the Newton iteration scheme in Helmholtz).
     This will cause changes in simulation output. (#368)

   * the problem-specific runtime parameters (probin) are now
     written to the job_info file (#380)

   * we now skip the initial EOS call prior to the burn  * this
     was redundant because we already did a clean_state (#377)

   * we now support recent versions of hypre (#373)

# 18.08

   * the old use_mlmg_solver parameters were removed  * this
     has been the only multigrid solver in Castro for some time,
     so the parameters had no effect.

   * The parameter dual_energy_eta3 was removed. This had been
     introduced mostly for testing purposes and ended up being
     unnecessary. Also, the EOS call at the beginning of the burn was
     removed; this should provide a modest computational gain. Answers
     may change at the level of the EOS tolerance (#377).

   * The functionality that checks for a regrid at the end of
     the timestep will now apply all tagging criteria, not
     just the t_sound / t_enuc criterion.

   * A bug was fixed that occurred when a retry was triggered in
     the middle of a group of subcycles, rather than at the
     beginning of the advance (#358).

   * A bug with the logic of refluxing when using retries was
     fixed (#357).

   * Tagging can now be done on relative gradients, in addition
     to the existing capability for absolute gradients (#354).
     For example, tempgrad_rel is a relative gradient criterion
     that will tag a zone with temperature T if any adjacent zone
     has a temperature that is different by more than tempgrad_rel * T.
     The tagging is enabled up to a given level with the parameter
     max_tempgrad_rel_lev. The corresponding new tagging criteria
     for other fields are named similarly.

   * Retries can now be done in a dynamic fashion (#179). An
     advance is itself a subcycled advance always, and we keep
     subcycling until we are done with the step. By default we
     still do only a single timestep when use_retry is not
     enabled, but this helps in cases where we have to reject
     the timestep for (say) a CFL violation, and then during
     the retry the CFL criterion is violated again. In the past,
     we would simply have to abort the run if this happened. Now
     we can cut the timestep again and keep going. Additionally,
     if you set abort_on_false to F in your probin file's extern
     parameters, then a burn in Microphysics will not cause an
     abort of the run, and Castro now knows how to deal with that
     by doing a retry and taking a shorter timestep (under the
     logic that most burn failures come from taking too large of
     a hydrodynamic timestep for the burner to be able to keep up).

# 18.07

   * A new GPU (CUDA) hydrodynamics solver (based on the
     method-of-lines solver) has been added, based on the work
     initially done in StarLord.  This is a work in progress, and
     requires the "gpu" branch of AMReX.

   * We removed all dependencies on the AMReX F_BoxLib source, in
     preparation for this source being removed in the future.

   * we now set the number of variables at compile time, by parsing
     the _variables file and interpreting which options are set in the
     preprocessor.  Note that a side effect of this change is that the
     number of radiation groups is now set at compile time instead of
     at runtime.  This change is needed for the GPU port.

     To set the number of radiation groups, set NGROUPS=4, e.g. for
     4 groups, in your problem's GNUmakefile.  Similar options exist
     for neutrinos.

     A related change is that it is now possible to set the number
     of advected quantities (that are not species or EOS auxillary
     fields) via NUMADV in your GNUmakefile.

# 18.06

   * The new multilevel multigrid solvers (MLMG) in the AMReX
     framework are now the default for self-gravity and constructing
     the operator for explicit diffusion.

   * A new test problem, the classic double Mach reflection, was
     added.

   * fix an issue in the retry logic that could sometimes result in
     an overflow of the estimated number of subcycle steps.

   * Improved the behavior of retries when they hit CFL violations
     (#334, #335).

# 18.05

   * Gamma_1 is now a derived variable

   * a new diffusion solver is implemented that uses the new muligrid
     framework in AMReX to compute the diffusive operator.  This can be
     enabled with diffusion.use_mlmg_solver = 1.

# 18.04

   * The job_info file now indicates which runtime parameters were
     changed from their default value (#314)

   * Improvements made to the 4th order hydro solver for single-level

# 18.03

   * The option ppm_trace_sources has been removed  * we now
     always trace on source terms with ppm.  Additionally, all
     sources are now traced, not just momentum sources.

   * The method-of-lines integrator has been rewritten.  It now works
     properly with sources.  Note: it is not intended for multilevel
     yet.  (#288, #287, #286, #164, #291, #137)

# 18.02

   * The approximate state Riemann solvers have been split into two
     parts: the computation of the interface state and the evaluation
     of the fluxes from this interface state.  This gives additional
     flexibililty in using these solvers in other methods.

# 18.01

   * The parameter dtnuc_mode has been removed. This was initially used
     for testing various forms of the burning timestep limiter before a
     final form was settled on.

   * Minor inconsistencies in how the external and diffusion source terms
     were constructed when simultaneously using reactions (#268, #269)
     have been fixed (#271).

   * The deprecated parameter castro.use_colglaz is removed. It was
     deprecated in June 2016 because it was obsoleted by the parameter
     castro.riemann_solver, which can be set to 1 to use the Colella
     and Glaz Riemann solver.

   * The state variable indicies in Fortran are now all defined in a
     single file, Source/driver/_variables.  This makes it much
     clearer and consistent and will allow for autodocumentation and
     clean-ups for GPU acceleration in the future.

# 17.12

   * The sponge can now operate based on pressure. The new parameters
     belong in the sponge namelist, and are named
     sponge_lower_pressure and sponge_upper_pressure. It works on the
     same principle as the density sponge.

   * The sponge can now drive the system to a particular velocity (the
     default is still zero velocity). The new parameters belong in the
     sponge namelist in your probin file, and are named
     sponge_target_{x,y,z}_velocity.

   * The SDC_Source_Type StateData was removed, as its purpose is now
     supplanted by the change to always keep the source terms in
     StateData (see below), and it was thus redundant. This does not
     change code output but does mean that old checkpoints generated
     while using SDC are no longer compatible with the current
     code. Normally we strive to maintain forward compatibility of
     checkpoints, but this change was considered justified because SDC
     is still considered an experimental feature under development and
     to our knowledge no production science runs have yet been done
     using SDC.

   * The parameter gravity.max_solve_level was removed. This was added
     to work around convergence issues in the multigrid solve, but
     those convergence issues have since been fixed, so the parameter
     is no longer used.

   * The Source_Type StateData now holds the actual old- and new-time
     sources (previously it held the source term predictor). This
     StateData is used to fill the sources_for_hydro MultiFab which
     provides the source terms used in the hydrodynamic update. Since
     it is StateData, this operation is done with a
     FillPatch. Consequently the sources_for_hydro data has meaningful
     data in both physical domain ghost zones and ghost zones at a
     coarse-fine interface (previously it only had meaningful data on
     fully interior ghost zones). Checkpoints will now have both old
     and new data. This change does result in a difference in the
     simulation output for simulations with source terms, as the
     answer will be different at the boundaries and at coarse-fine
     interfaces. (#116, #253) A related bug when using SDC was fixed
     too. (#56)

   * The parameter castro.keep_sources_until_end has been removed.

   * Source terms (gravity, rotation, etc.) have now been all
     coalesced into a single MultiFab. This reduces the memory
     footprint of the code. This negates the need for the code
     parameters update_state_between_sources and
     coalesce_update_diagnostics, so they have been removed. This will
     cause a change in results: the previous code behavior was to
     update the state with each successive source term as it was
     applied at the new time. Now every source term will be calculated
     using the same new-time state (the one coming out of the hydro)
     and the source terms are all applied to the state in one shot at
     the end. Note that both this and the old approach are
     second-order accurate in time. For problems that combine multiple
     source terms, this will cause changes that are larger than
     roundoff.  In particular, if rotation is used in conjunction with
     other source terms, changes will be relatively large, because the
     Coriolis force depends on the velocity, so the source term is a
     bit different now that it is seeing a different velocity. (#165,
     #249)

   * As of 17.10, there is a new option castro.plot_per_is_exact. If
     this is set to 1, timesteps will be shortened to exactly hit the
     time interval specified by amr.plot_per. An issue with this
     (#242) was fixed (#243) where an incorrect timestep would be
     taken after a restart if the previous step had been shortened.

   * We can now use the new multigrid solver from AMReX (implemented
     in C++) instead of the older Fortran solver (#241).  This is
     enabled with gravity.use_mlmg_solver=1.  Note, this only works in
     3-d currently.  This has several options:

     gravity.mlmg_max_fmg_iter = 0 : This integer parameter determines
     how many FMG cycles will be performed before switching to
     V-cycle.

     gravity.mlmg_agglomeration = 0 : This boolean flag determines if
     AMR level 0 grids will be agglomerated as the grids are coarsen
     in the multi-grid hierarchy.

     gravity.mlmg_consolidation = 0 : This boolean flag determines if
     grids on an AMR fine level that is in a single-level solve or the
     lowest AMR level of a multi-level composite solve will be
     consolidated as the grids are coarsen in the multi-grid
     hierarchy. Here, consolidation means the grids will be moved to
     fewer MPI processes.

     Numerical experiments have show this scales better than the old
     solver.

   * Apply the sources to the state's ghost zones (#255).  This fixes
     #213  * it ensures the advances give a valid update for the ghost
     zones in the State_Type.


# 17.11.1

   * Minor bug fixes from the 17.11 release. There is a corresponding
     17.11.1 release of AMReX.

# 17.11

   * A bug was fixed in non-Cartesian simulations with AMR
     (1D spherical and 2D cylindrical). The bug was introduced
     around version 17.02 and resulted in incorrect synchronization
     of the pressure term in the momentum equations. The effect
     would have manifested as non-conservation of momentum or
     strange effects at coarse-fine interfaces.

   * The sponge is now always time centered. The option to
     do otherwise was introduced in 17.02, and has now been
     removed. Additionally, the form of the energy source
     term has been corrected for the time centered case,
     and brought into line with how we do the energy source
     term for other sources. (Issue #7, Issue #57)

   * Fixed a bug in the fix for #188.

   * Conductivity_dir has been renamed CONDUCTIVITY_DIR to be consistent
     with EOS_DIR and NETWORK_DIR

   * we no longer get the compositional derivatives as part of the EOS
     call.  If you need this functionality, you need to set the
     preprocessor variable (in Fortran), EXTRA_THERMO

   * you can now use a system's BLAS routines, instead of compiling
     the versions from Microphysics by setting USE_SYSTEM_BLAS=TRUE.
     This then looks at BLAS_LIBRARY for the link line.


# 17.10

   * It is sometimes useful to be able to do some sort of initialization
     phase in your simulation, stop with a checkpoint, and then restart
     (possibly with different options) to do the main phase of the run.
     In this case, you may want to reset the simulation time to zero for
     analysis purposes. The new option castro.reset_checkpoint_time allows
     you to do this: by setting it to the time you want, the checkpoint you
     generate will have this new time. Similarly, castro.reset_checkpoint_step
     allows you to reset the timestep number (for example, to 0). Both options
     only work when you're using amr.checkpoint_on_restart=1, which itself
     requires amr.regrid_on_restart=1. This option is only intended to be used
     for the case where you're generating this checkpoint, so you also need
     to temporarily set max_step and stop_time to the target values you're
     resetting them to, to prevent further steps after the restart. After you
     have the new checkpoint, then you can undo those temporary variables
     and continue your run as usual.

   * A minor error in the gravity source terms was fixed (#109).
     This error should not normally have been observable.

   * fixed a bug in the artifical viscosity in 1-d in
     non-Cartesian geometries (issue #175)

   * the README.md now describes the process to become a
     "core developer" of Castro, and what this means.

   * Network_dir has been renamed NETWORK_DIR and EOS_dir has been
     renamed EOS_DIR.  All of the problem GNUmakefiles have been
     updated.  The old names will continue to work in the near future,
     but users are encouraged to change any of their problems to use
     the new all-caps names (PR #184)

   * the density flux limiting functionality now has a small tolerance
     (#185). It has also been documented (#193).

   * the timestep retry is now conservative  * this was accomplished
     by saving the density fluxes to use in the conservative gravity
     update (#178). Also, a bug in the timestep retry for Poisson
     gravity was fixed (#188).

# 17.09

   * the Source/ directory structure has been reorganized,
     putting the source files into directories by physics and
     eliminating the Src_1d, Src_2d, ... subdirectories

   * the Riemann solvers have been merged into a single
     dimensional-agnostic version in Src_nd.  In 2-d there was an
     issue with how the Godunov state for the CG solver was stored on
     interfaces, which would affect the internal enery evolution.

   * the PLM and PPM reconstruction routines were merged into
     a single dimensional-agnostic version in hydro/

   * the characteristic tracing routines were merged into
     dimensional-agnostic versions in hydro/ and radiation/.  This
     change fixed and outstanding issue  * the PLM reconstruction in
     1-d now uses a reference state. (issue #11)

# 17.08

   * the option castro.limit_fluxes_on_small_dens now only limits
     on density as the name suggests. It originally also limited
     fluxes if the internal energy would go negative, but this
     caused problems in runs with MPI, so it was removed. It was
     not strictly needed anyway, as the normal logic for handling
     negative internal energies is reliable.

   * two errors were fixed in the implementation of the triggered
     regrid at the end of a timestep. The method now correctly conserves
     fluid variables at coarse-fine boundaries.

   * the XGRAPH stuff that output xmgr-compatible 1-d ASCII profiles

   * fixed a bug where the gravity runtime parameters were not
     being properly initialized in the Fortran side of the code.

   * the viscosity routine is now separate from conductivity
     in Microphysics/.  Also, Castro can now use the stellar
     conductivity that is part of StarKiller.

   * the StarKiller-astro Microphysics repo now uses a denser table
     for the Helmholtz EOS (thanks to Frank Timmes).  If you are using
     this EOS, the new table will be soft-linked to your build
     directory automatically.  If you have an old copy laying around,
     it might fail to run, with an I/O error.

# 17.07

   * start of some code cleaning for eventual GPU offload support
     merging from StarLord

   * added a framework for method-of-lines integration of the
     hydrodynamics.  This does not do characteristic tracing, but
     instead does reconstruction through multiple stages of an ODE
     integrator.  At the moment, radiation is not supported.

# 17.06

   * we now require the AMReX library instead of the BoxLib library

   * the Microphysics repository that we rely on for the EOS and
     reaction networks is now part of the starkiller-astro github.
     You can change your clone to refer to this via:

     git remote set-url origin ssh://git@github.com/starkiller-astro/Microphysics

   * a new mechanism for using a stability criterion to trigger
     a regrid at the end of a timestep was added (PR #122)

   * some cleaning of the logic for momentum fluxes and
     limit_hydro_fluxes_on_small_dens (issues #130, #131)

# 17.05

   * some protections added in the retry code


# 17.04

   * rewrote the conservative gravity formulation to work off of the
     potential.  This gives the best conservation with AMR.  This is
     equivalent to the description in Appendix B from Katz et
     al. (2016).


# 17.03

   * the new refluxing method introduced in 16.11 has been removed,
     as it was determined to not provide any benefit in accuracy.

   * new derived plot variables are available when using thermal
     diffusion, including the conductivity, diffusion coefficient, and
     the entire diffusion term to the energy equation.  (issue #104)

   * when all derived variables were stored in the plotfile, we were
     storing the mass fractions twice.  E.g. for he4, we were saving
     "he4" and "X(he4)".  Now only the latter is stored.

   * created a post_simulation() function that is called at the end of
     a simulation.  An example is provided by test_diffusion where we
     output the norm of the error against the analytic solution.
     (issue #107, 108)


# 17.02

   * diagnostic information about how source terms update the state
     has been overhauled and made uniform. All source terms, including
     hydro and resets, print their changes to the state in the same
     format. The parameter print_energy_diagnostics has been renamed
     print_update_diagnostics, and a new parameter
     coalesce_update_diagnostics has been added so that you can
     combine all of the old-time and new-time updates into one print.
     (issue #58)

   * to support both single and double precision, all of the floating
     point declarations use the amrex_real type defined in the
     amrex_fort_module  * this is set to single or double precision at
     compile time.  All constants also now use this type.  For
     brevity, we rename it to 'rt' in the use statement.  (issue #34)

   * the sponge is now time-centered by default (issue #7)

   * the ppm_temp_fix stuff has been documented and made
     consistent across dimensions (issue #25)

   * the job info git information is now the output of git describe,
     this gives more information, including the last tag, how far we
     are from the tag, and an abbreviated hash.  It also indicates if
     your working directory is dirty


# 17.01

   * the radiation-specific version of computeTemp has been removed
     and instead everything goes through the main Castro computeTemp.
     This affects, in particular, how we treated small internal
     energies in radiation. (issue #64)

   * the radiation-specific versions of umeth and consup have been
     merged with the pure hydro routines.  This gives round-off level
     differences.  This also now uses all velocity components in the
     kinetic energy correction for radiation. (issues #66, 70)

   * a minor bug was fixed in the 3-d radiation characteristic tracing,
     regarding which gamma_1 (gamc) is used.


# 16.12a

   * fix a restart bug with radiation that was introduced after 16.10
     (this was cherry-picked from development) (issues #76, 78)

# 16.12

   * BoxLib now requires a C++ 11 compiler by default.  As part of
     this transition, PArrays are replaced by C++ Arrays.
     Additionally, changes to the BoxLib build system mean that we
     only need to supple COMP for the compiler.  FCOMP is now
     ignored.

   * The User's Guide has been updated to reflect the current flow of
     the algorithm.

# 16.11

   * we now distinguish between gravity (which can include a
     constant gravitational acceleration) and self-gravity,
     with the GRAVITY and SELF_GRAVITY preprocessor directives

   * some work on the sync between levels was done  * this will
     be described in a forthcoming paper. The main change by default
     is that after a reflux, we recompute the value of the source terms
     on the affected levels so that the new-time source term knows about
     the updated state due to the flux. For gravity, this resembles
     what the original Castro paper described for a sync source, but this
     is now done in a consistent way for all source terms. This should be fairly
     cheap which is why it is enabled by default, but you can disable it
     (see castro.update_sources_after_reflux). An additional optional
     change is a new strategy for refluxing (see castro.reflux_strategy).
     In the existing standard method, we only reflux after all fine timesteps
     over a coarse timestep have completed. In the new method, we do a
     partial reflux at the end of each fine timestep. This means that
     the coarse state used in time interpolation for the fine level
     is slightly more accurate as we go for later fine timesteps. It
     should also be needed for self-consistently conserving energy for gravity.
     At present it is more expensive than the standard method when there
     are gravity sync solves because there are more of them, but the tradeoff
     is that the simulation is more accurate.

   * the order of computing the temperature and reseting internal
     energy was changed in a few spots. This will change results by default.

   * the radiation-specific source was moved into the Radiation/
     subdirectory

# 16.10

   * the parameter first_order_hydro has been moved from the
     radiation namespace to the castro namespace

   * the problem setups have been moved into sub-directory
     categories to make it easier to read (issue #32)

   * the way we we use tolerances in the multigrid solve for Poisson
     gravity has changed. The old behavior is that you would pass in
     gravity.ml_tol as the relative tolerance on each grid level, and
     absolute tolerances would not be used. This suffered from some
     defects, notably that on fine grids you often had to loosen the
     relative tolerance on each higher level to achieve convergence,
     and in certain cases the scheme would fail completely, for
     example if the fine grids were not covering the mass on the
     grid. We now use an input parameter gravity.abs_tol which
     controls the absolute scale of the tolerance. This can either be
     an array of values, one for each level, or a single scalar
     value. If it is the latter, then the absolute tolerance passed
     into the multigrid scheme is the tolerance multiplied by the
     maximum value of the RHS over the entire domain.  On the coarse
     grid, then, the absolute tolerance is 4*pi*G*rho_max*abs_tol, and
     on fine grids this is multiplied by ref_ratio**2. If you do not
     specify gravity.abs_tol, then a reasonable value is selected for
     the coarse level, and the same scheme is used to give it
     reasonable values on the fine levels as well. The parameter
     gravity.ml_tol has been renamed gravity.rel_tol, and has the same
     meaning as before, but it now defaults to zero. gravity.ml_tol is
     now deprecated, and will be removed in a future release. Note
     that the tolerance used in determining convergence is always the
     less restrictive of the relative and absolute tolerance
     requirements.  gravity.delta_tol has been removed. (issue #43)

   * the radiation hydro solver, that used to live in
     CastroRadiation.git has now been completely integrated into the
     main Castro git repo.  The history was preserved in the
     transition It has also been cleaned up a little (issues #24, #31,
     #33, #48)

     The radiation build variable Network_inputs was renamed
     to NETWORK_INPUTS for consistency.

     The EOSes that used to come with CastroRadiation are available
     in Microphysics.git

   * the gravity and diffusion runtime parameters have been moved
     to the master _cpp_parameters file (issue #42)

   * enthalpy, species, and temperature diffusion are now properly
     time-centered (issue #22), and a bug in the hi boundary
     inflow boundary conditions for diffusion was fixed (issue #41)

   * a flux limiter has been added that limits the size of the hydro
     fluxes if they would cause rho or (rho e) to go negative. This
     can be used with castro.limit_hydro_fluxes_on_small_dens = 1.

   * a bug for single-level problems with Poisson gravity has been
     fixed where the multi-grid tolerance was being set to an
     uninitialized value

   * a flaw in the source terms for the primitive variables in the
     hydro update has been fixed, so that source terms like gravity
     should no longer affect the pressure and (rho e) interface states
     (issue #19)

   * the prediction of source terms to the half-time used in the
     hydrodynamics reconstruction was not being done properly.
     This has been fixed (issue #18)

   * the radiation hydro ppm now implements the ppm_predict_gammae
     option

   * we no longer ship VODE or BLAS with Castro  * these are provided
     by the separate Microphysics git repo

   * the documentation of the architecture of Castro has been
     significantly improved (issues #20, #23, #29, #31)

# 16.09:

   * the PPM tracing routine for radiation was synced up with the pure
     hydro version.  In particular, it now supports ppm_trace_sources,
     implements the reference states and fixes an issue with the
     flattening.

   * The 1-d PPM routine was also updated to support tracing,
     predicting gamma_e instead of (rho e), and an inconsistency in
     the flattening was fixed.

   * the parameters ppm_reference and ppm_reference_edge_limit
     have been removed  * there was no reason to use anything other
     than the defaults

   * the parameter ppm_tau_in_tracing has been removed.  The useful
     part of this is preserved in the ppm_predict_gammae = 1
     functionality, which uses a different set of primitive variables
     (tau, u, p, gamma_e) in the prediction of the interface states.

   * the flux correction for axisymmetric and 1-d spherical coords
     has been fixed.  In particular, there is now a separate flux
     register for the pressure term that enters as a gradient in the
     momentum equation.

   * The sign on the gravitational potential has been flipped to be
     consistent with the usual convention in the physics literature,
     i.e. the potential is negative and we solve del**2 phi = 4 * pi *
     G * rho.

   * Castro_advance.cpp has been significantly cleaned up. Each source
     term (gravity, rotation, diffusion, etc.) has a MultiFab
     associated with it through which it affects the state data. This
     has changed results slightly (typically a relative change no
     larger than the 1e-4 level) due to updates in the order of
     operations and consistency in usage on ghost zones.

   * An iterative solver for coupling between reactions and
     hydrodynamics has been introduced, which you can enable with
     USE_SDC = TRUE in the makefile. The number of iterations done for
     each timestep is controlled with castro.sdc_max_iters.

   * We changed the defaults for the gravity and rotation sources.
     now we do grav_source_type and rot_source_type = 4 by default.
     This is a conservative formulation for the energy equation that
     incorporates the source as a flux in the energy equation.  See
     Katz et al. 2016 for details.

     We also do implicit_rotation_update = 1 by default  * this does a
     slightly better coupling of the Coriolis force in the momentum
     equation by doing an implicit velocity update

     We also make ppm_trace_sources = 1 the default  * this does
     parabolic reconstruction of the momentum sources and traces
     under them when constructing the interface states

   * we now set castro.cg_blend = 2 by default.  This has no effect for
     the default CGF Riemann solver, but for the Colella & Glaz solver
     (castro.riemann_solver = 1), this will augment the secant iteration
     for the pstar find with bisection if we fail to converge.  This
     makes the root find for the star state more robust.

   * a new "fake" setup, riemann_test_zone can be use to send a left /
     right hydro state to the CG Riemann solver for testing  * this acts
     as a unit test for that solver.

   * the default for castro.allow_negative_energy is now 0  * this is
     the safer choice.

   * the default for castro.point_mass_fix_solution was changed to 0
      * this is a more expected behavior for new users.


# 16.08

   * A new parameter gravity.max_multipole_moment_level was added.
     This comes into play when using the multipole solver to compute
     the boundary conditions on the domain for isolated mass
     distributions.  The default behavior in Castro when constructing
     boundary conditions for the gravity solve is to do a multipole
     expansion sum over the density on the coarse grid only.  If you
     increase the value of that new parameter from its default of 0 to
     some higher number, it will use the data from those more refined
     levels in constructing the boundary condition values.

   * The file sponge_nd.f90 in Source/Src_nd/ has been renamed to
     sponge_nd.F90, the file extension change indicating that it can
     now be run through the preprocessor. Please update your local
     name for this file if you're overriding it in your problem setup.

   * The sponge update is usually done in an implicit fashion,
     but you can now instead do an explicit update with
     castro.sponge_implicit == 0.

   * the shock variable is now output if we are running with shock
     detection enabled

   * Microphysics/eos is now Microphysics/EOS

   * a number of changes were done in the Microphysics repo  * see
     Microphysics/CHANGES for a log of those


# 16.07

   * For consistency across the BoxLib suite of astro codes, we've
     renamed the main environment variables.  CASTRO_HOME now replaces
     CASTRO_DIR; MICROPHYSICS_HOME now replaces MICROPHYSICS_DIR.

   * The EOS, network, and conductivity routines have been moved to
     sub-directories or Castro/Microphysics/.  This reflects the way
     the layout in the standalone Microphysics repo as well as that in
     Maestro.

   * Some of the routines in Source/Src_nd/ have been renamed from
     .f90 files to .F90 files so that we can use the preprocessor. If
     you were using any of them (Prob_nd.f90, problem_tagging_nd.f90,
     Rotation_frequency.f90, or ext_src_nd.f90) by having local copies
     in your problem directory that overwrote them, please be sure to
     update the file extension so that Castro will recognize them.

   * If you were using allow_negative_energy == 0, the case where
     (rho*E), the total gas energy of the zone, was negative was
     indirectly covered and it would be reset in this case due to the
     way the logic worked for resetting the internal energy and then
     updating the total energy to be consistent with it. However at
     one point we added an option castro.dual_energy_update_E_from_e
     which disabled that second update and also meant that negative
     (rho*E) was again possible. This possibility has now been
     precluded directly, by resetting (rho*E) the same way if we
     detect that it is negative. This should not change results unless
     you were using castro.dual_energy_update_E_from_e = 1. This is
     also a good time to plug the newer option
     castro.allow_small_energy, which if set to 1 will reset when you
     hit a (rho*e) that is less than the smallest possible energy for
     the (rho, small_temp, X) in that zone. Note that it requires an
     extra EOS call.

   * The default interpolation for coarse zones into fine zones is
     piecewise linear.  There is now an option to use piecewise
     constant instead  * set castro.state_interp_order to 0. Note that
     if you use piecewise linear you can set
     castro.lin_limit_state_interp to 1 if you want to preserve linear
     combinations and therefore guarantee that, say, sum(X) = 1.

   * If you set the new option castro.limit_fluxes_on_small_dens = 1,
     the fluxes will be explicitly limited such that a negative
     density is never created.

   * Along similar lines, there are also new options for how to reset
     a negative density if one should arise. Set
     castro.density_reset_method = 2 to use the average of all
     adjacent zones instead of the default, which is the
     characteristics of the adjacent zone with the highest
     density. Set it to 3 if you want to reset it to the original zone
     state before the hydro update.

   * We have fixed an issue where diffusion did not work correctly if
     add_ext_src = 0.  The diffusion source term is now independent of
     whether you have user-defined source terms.

   * ConvertCheckpoint/ now lives under Util/

   * UsersGuide/ is now Docs/  * this is consistent with the other
     BoxLib codes

   * Burning is no longer done in ghost cells for boundaries with
     neighbors on the same level of refinement.  Instead a ghost cell
     fill is done to fill the like-level neighbor cells.  As a
     consequence of this change, if reset_internal_energy() is invoked
     in a cell, to reset the internal energy to E - K, this reset is
     now reflected in the ghost cells (this is a more consistent
     behavior).  Previously, the energy was never reset in the ghost
     cells.
