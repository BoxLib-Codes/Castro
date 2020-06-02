
#include "Castro.H"
#include "Castro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

#ifdef GRAVITY
#include "Gravity.H"
#endif

#include <climits>
#include <cmath>

using std::string;
using namespace amrex;

Real Castro::advance(Real time, Real dt, int amr_iteration, int amr_ncycle)

// the main driver for a single level.  This will do either the SDC
// algorithm or the Strang-split reactions algorithm.
//
// arguments:
//    time          : the current simulation time
//    dt            : the timestep to advance (e.g., go from time to
//                    time + dt)
//    amr_iteration : where we are in the current AMR subcycle.  Each
//                    level will take a number of steps to reach the
//                    final time of the coarser level below it.  This
//                    counter starts at 1
//    amr_ncycle    : the number of subcycles at this level

{
    BL_PROFILE("Castro::advance()");

    // Save the wall time when we started the step.

    wall_time_start = ParallelDescriptor::second();

    MultiFab::RegionTag amrlevel_tag("AmrLevel_Level_" + std::to_string(level));

    Real dt_new = dt;

    initialize_advance(time, dt, amr_iteration, amr_ncycle);

    // Do the advance.

    if (time_integration_method == CornerTransportUpwind ||
        time_integration_method == SimplifiedSpectralDeferredCorrections) {

        dt_new = std::min(dt_new, subcycle_advance_ctu(time, dt, amr_iteration, amr_ncycle));

#ifndef MHD
#ifndef AMREX_USE_CUDA
#ifdef TRUE_SDC
    } else if (time_integration_method == SpectralDeferredCorrections) {

        for (int iter = 0; iter < sdc_order + sdc_extra; ++iter) {
            sdc_iteration = iter;
            dt_new = do_advance_sdc(time, dt, amr_iteration, amr_ncycle);
        }

#endif  // TRUE_SDC
#endif  // AMREX_USE_CUDA
#endif  //MHD
    }

    // Optionally kill the job at this point, if we've detected a violation.

    if (cfl_violation && !use_retry) {
        amrex::Abort(
            "CFL is too high at this level; go back to a checkpoint and restart with lower CFL "
            "number, or set castro.use_retry = 1");
    }

    // If we didn't kill the job, reset the violation counter.

    cfl_violation = 0;

    if (use_post_step_regrid) {
        check_for_post_regrid(time + dt);
    }

#ifdef AUX_UPDATE
    advance_aux(time, dt);
#endif

#ifdef GRAVITY
#if (BL_SPACEDIM > 1)
    // We do this again here because the solution will have changed
    if ((level == 0) && (spherical_star == 1)) {
        int is_new = 1;
        make_radial_data(is_new);
    }
#endif
#endif

#ifdef GRAVITY
    // Update the point mass.
    if (use_point_mass) {
        pointmass_update(time, dt);
    }
#endif

#ifdef RADIATION
    MultiFab& S_new = get_new_data(State_Type);
    final_radiation_call(S_new, amr_iteration, amr_ncycle);
#endif

#ifdef AMREX_PARTICLES
    advance_particles(amr_iteration, time, dt);
#endif

    finalize_advance(time, dt, amr_iteration, amr_ncycle);

    return dt_new;
}

void Castro::initialize_do_advance(Real time, Real dt, int amr_iteration, int amr_ncycle) {
    BL_PROFILE("Castro::initialize_do_advance()");

    // Reset the CFL violation flag.

    cfl_violation = 0;

#ifdef RADIATION
    // make sure these are filled to avoid check/plot file errors:
    if (do_radiation) {
        get_old_data(Rad_Type).setBndry(0.0);
        get_new_data(Rad_Type).setBndry(0.0);
    } else {
        get_old_data(Rad_Type).setVal(0.0);
        get_new_data(Rad_Type).setVal(0.0);
    }
#endif

    // Reset the grid loss tracking.

    if (track_grid_losses) {
        for (int i = 0; i < n_lost; i++) {
            material_lost_through_boundary_temp[i] = 0.0;
        }
    }

#ifdef GRAVITY
    if (moving_center == 1) {
        define_new_center(get_old_data(State_Type), time);
    }

#if (BL_SPACEDIM > 1)
    if ((level == 0) && (spherical_star == 1)) {
        swap_outflow_data();
        int is_new = 0;
        make_radial_data(is_new);
    }
#endif
#endif

    // For the hydrodynamics update we need to have NUM_GROW ghost
    // zones available, but the state data does not carry ghost
    // zones. So we use a FillPatch using the state data to give us
    // Sborder, which does have ghost zones.

    if (time_integration_method == CornerTransportUpwind ||
        time_integration_method == SimplifiedSpectralDeferredCorrections) {
#ifdef MHD
        MultiFab& Bx_old = get_old_data(Mag_Type_x);
        MultiFab& By_old = get_old_data(Mag_Type_y);
        MultiFab& Bz_old = get_old_data(Mag_Type_z);

        Bx_old_tmp.define(Bx_old.boxArray(), Bx_old.DistributionMap(), 1, NUM_GROW);
        By_old_tmp.define(By_old.boxArray(), By_old.DistributionMap(), 1, NUM_GROW);
        Bz_old_tmp.define(Bz_old.boxArray(), Bz_old.DistributionMap(), 1, NUM_GROW);

        FillPatch(*this, Bx_old_tmp, NUM_GROW, time, Mag_Type_x, 0, 1);
        FillPatch(*this, By_old_tmp, NUM_GROW, time, Mag_Type_y, 0, 1);
        FillPatch(*this, Bz_old_tmp, NUM_GROW, time, Mag_Type_z, 0, 1);
#endif
        // for the CTU unsplit method, we always start with the old
        // state note: a clean_state has already been done on the old
        // state in initialize_advance so we don't need to do another
        // one here
        Sborder.define(grids, dmap, NUM_STATE, NUM_GROW, MFInfo().SetTag("Sborder"));
        const Real prev_time = state[State_Type].prevTime();
        expand_state(Sborder, prev_time, NUM_GROW);

    } else if (time_integration_method == SpectralDeferredCorrections) {

        // we'll handle the filling inside of do_advance_sdc
        Sborder.define(grids, dmap, NUM_STATE, NUM_GROW, MFInfo().SetTag("Sborder"));

    } else {
        amrex::Abort("invalid time_integration_method");
    }

#ifdef SHOCK_VAR
    // Zero out the shock data, and fill it during the advance.
    // For subcycling cases this will always give the shock
    // variable for the latest subcycle, rather than averaging.

    Sborder.setVal(0.0, USHK, 1, Sborder.nGrow());
#endif
}

void Castro::finalize_do_advance(Real time, Real dt, int amr_iteration, int amr_ncycle) {
    BL_PROFILE("Castro::finalize_do_advance()");

#ifdef RADIATION
    if (!do_hydro && Radiation::rad_hydro_combined) {
        MultiFab& Er_old = get_old_data(Rad_Type);
        MultiFab& Er_new = get_new_data(Rad_Type);
        Er_new.copy(Er_old);
    }
#endif

    Sborder.clear();
}

void Castro::initialize_advance(Real time, Real dt, int amr_iteration, int amr_ncycle) {
    BL_PROFILE("Castro::initialize_advance()");

    // Save the current iteration.

    iteration = amr_iteration;

    sub_iteration = 0;
    sub_ncycle = 0;
    dt_subcycle = 1.e200;
    dt_advance = dt;

    keep_prev_state = false;

    // Reset the retry timestep information.

    lastDtRetryLimited = 0;
    lastDtFromRetry = 1.e200;
    in_retry = false;

    if (use_post_step_regrid && level > 0) {

        if (getLevel(level - 1).post_step_regrid && amr_iteration == 1) {

            // If the level below this just triggered a special regrid,
            // the coarse contribution to this level's FluxRegister
            // is no longer valid because the grids have, in general, changed.
            // Zero it out, and add them back using the saved copy of the fluxes.

            getLevel(level - 1).FluxRegCrseInit();
        }
    }

    // Pass some information about the state of the simulation to a Fortran module.

    ca_set_amr_info(level, amr_iteration, amr_ncycle, time, dt);

    // The option of whether to do a multilevel initialization is
    // controlled within the radiation class.  This step belongs
    // before the swap.

#ifdef RADIATION
    if (do_radiation) {
        radiation->pre_timestep(level);
    }

    Erborder.define(grids, dmap, Radiation::nGroups, NUM_GROW);
    lamborder.define(grids, dmap, Radiation::nGroups, NUM_GROW);
#endif

#ifdef GRAVITY
    // If we're on level 0, update the maximum density used in the gravity solver
    // for setting the tolerances. This will be used in all level solves to follow.
    // This must be done before the swap because it relies on the new data.

    if (level == 0 && gravity->get_gravity_type() == "PoissonGrav") {
        gravity->update_max_rhs();
    }
#endif

    // This array holds the sum of all source terms that affect the
    // hydrodynamics.

    sources_for_hydro.define(grids, dmap, NSRC, NUM_GROW);
    sources_for_hydro.setVal(0.0, NUM_GROW);

    // This array holds the source term corrector.

    source_corrector.define(grids, dmap, NSRC, NUM_GROW);
    source_corrector.setVal(0.0, NUM_GROW);

    // Swap the new data from the last timestep into the old state data.

    swap_state_time_levels(dt);

#ifdef GRAVITY
    if (do_grav) {
        gravity->swapTimeLevels(level);
    }
#endif

    // Ensure data is valid before beginning advance. This addresses
    // the fact that we may have new data on this level that was interpolated
    // from a coarser level, and the interpolation in general cannot be
    // trusted to respect the consistency between certain state variables
    // (e.g. UEINT and UEDEN) that we demand in every zone.

    MultiFab& S_old = get_old_data(State_Type);
    clean_state(
#ifdef MHD
        get_old_data(Mag_Type_x), get_old_data(Mag_Type_y), get_old_data(Mag_Type_z),
#endif
        S_old, time, S_old.nGrow());

    // Initialize the previous state data container now, so that we can
    // always ask if it has valid data.

    for (int k = 0; k < num_state_type; ++k) {
        prev_state[k].reset(new StateData());
    }

    // This array holds the hydrodynamics update.
    if (time_integration_method == CornerTransportUpwind ||
        time_integration_method == SimplifiedSpectralDeferredCorrections) {
        hydro_source.define(grids, dmap, NUM_STATE, 0);
    }

    // Allocate space for the primitive variables.

    q.define(grids, dmap, NQ, NUM_GROW);
    q.setVal(0.0);
    qaux.define(grids, dmap, NQAUX, NUM_GROW);

    if (sdc_order == 4) {
        q_bar.define(grids, dmap, NQ, NUM_GROW);
        qaux_bar.define(grids, dmap, NQAUX, NUM_GROW);
#ifdef DIFFUSION
        T_cc.define(grids, dmap, 1, NUM_GROW);
#endif
    }

#ifdef TRUE_SDC
    if (time_integration_method == SpectralDeferredCorrections) {

        MultiFab& S_old = get_old_data(State_Type);
        k_new.resize(SDC_NODES);

        k_new[0].reset(new MultiFab(S_old, amrex::make_alias, 0, NUM_STATE));
        for (int n = 1; n < SDC_NODES; ++n) {
            k_new[n].reset(new MultiFab(grids, dmap, NUM_STATE, 0));
            k_new[n]->setVal(0.0);
        }

        A_old.resize(SDC_NODES);
        for (int n = 0; n < SDC_NODES; ++n) {
            A_old[n].reset(new MultiFab(grids, dmap, NUM_STATE, 0));
            A_old[n]->setVal(0.0);
        }

        A_new.resize(SDC_NODES);
        A_new[0].reset(new MultiFab(*A_old[0], amrex::make_alias, 0, NUM_STATE));
        for (int n = 1; n < SDC_NODES; ++n) {
            A_new[n].reset(new MultiFab(grids, dmap, NUM_STATE, 0));
            A_new[n]->setVal(0.0);
        }

        // We use Sburn a few ways for the SDC integration.  First, we
        // use it to store the initial guess to the nonlinear solve.
        // Second, at the end of the SDC update, we copy the cell-center
        // reaction source into it, including one ghost cell, for later
        // filling of the plotfile.  Finally, we use it as a temporary
        // buffer for when we convert the state to centers while making the
        // source term
        Sburn.define(grids, dmap, NUM_STATE, 2);

#ifdef REACTIONS
        R_old.resize(SDC_NODES);
        for (int n = 0; n < SDC_NODES; ++n) {
            R_old[n].reset(new MultiFab(grids, dmap, NUM_STATE, 0));
            R_old[n]->setVal(0.0);
        }
#endif
    }
#endif

    // Zero out the current fluxes.

    for (int dir = 0; dir < 3; ++dir) {
        fluxes[dir]->setVal(0.0);
        mass_fluxes[dir]->setVal(0.0);
    }

#if (BL_SPACEDIM <= 2)
    if (!Geom().IsCartesian()) {
        P_radial.setVal(0.0);
    }
#endif

#ifdef RADIATION
    if (Radiation::rad_hydro_combined) {
        for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
            rad_fluxes[dir]->setVal(0.0);
        }
    }
#endif
}

void Castro::finalize_advance(Real time, Real dt, int amr_iteration, int amr_ncycle) {
    BL_PROFILE("Castro::finalize_advance()");

    // Add the material lost in this timestep to the cumulative losses.

    if (track_grid_losses) {

        ParallelDescriptor::ReduceRealSum(material_lost_through_boundary_temp, n_lost);

        for (int i = 0; i < n_lost; i++) {
            material_lost_through_boundary_cumulative[i] += material_lost_through_boundary_temp[i];
        }
    }

    if (do_reflux) {
        FluxRegCrseInit();
        FluxRegFineAdd();
    }

    if (time_integration_method == CornerTransportUpwind ||
        time_integration_method == SimplifiedSpectralDeferredCorrections) {
        hydro_source.clear();
    }

    q.clear();
    qaux.clear();

    if (sdc_order == 4) {
        q_bar.clear();
        qaux_bar.clear();
#ifdef DIFFUSION
        T_cc.clear();
#endif
    }

#ifdef RADIATION
    Erborder.clear();
    lamborder.clear();
#endif

    source_corrector.clear();
    sources_for_hydro.clear();

    if (!keep_prev_state) {
        amrex::FillNull(prev_state);
    }

#ifdef TRUE_SDC
    if (time_integration_method == SpectralDeferredCorrections) {
        k_new.clear();
        A_new.clear();
        A_old.clear();
#ifdef REACTIONS
        R_old.clear();
        Sburn.clear();
#endif
    }
#endif

    // Record how many zones we have advanced.

    num_zones_advanced += static_cast<Real>(grids.numPts()) / getLevel(0).grids.numPts();

    Real wall_time = ParallelDescriptor::second() - wall_time_start;
    Real fom_advance = grids.numPts() / wall_time / 1.e6;

    if (verbose >= 1) {
        amrex::Print() << "  Zones advanced per microsecond at this level: " << fom_advance
                       << std::endl
                       << std::endl;
    }
}
