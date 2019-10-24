
#include "Castro.H"
#include "Castro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

#ifdef SELF_GRAVITY
#include "Gravity.H"
#endif

using namespace amrex;

Real
Castro::do_advance_ctu(Real time,
                       Real dt,
                       int  amr_iteration,
                       int  amr_ncycle)
{

    // this routine will advance the old state data (called S_old here)
    // to the new time, for a single level.  The new data is called
    // S_new here.  The update includes reactions (if we are not doing
    // SDC), hydro, and the source terms.

    BL_PROFILE("Castro::do_advance_ctu()");

    const Real prev_time = state[State_Type].prevTime();
    const Real  cur_time = state[State_Type].curTime();

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    // Perform initialization steps.

    initialize_do_advance(time, dt, amr_iteration, amr_ncycle);

#ifndef AMREX_USE_CUDA
    // Check for NaN's.

    check_for_nan(S_old);
#endif

    // Since we are Strang splitting the reactions, do them now


#ifdef REACTIONS
    if (time_integration_method != SimplifiedSpectralDeferredCorrections) {

        // this operates on Sborder (which is initially S_old).  The result
        // of the reactions is added directly back to Sborder.
        strang_react_first_half(prev_time, 0.5 * dt);

    }
#endif

    // Initialize the new-time data. This copy needs to come after the
    // reactions.

    MultiFab::Copy(S_new, Sborder, 0, 0, NUM_STATE, S_new.nGrow());

#ifdef REACTIONS
    if (time_integration_method != SimplifiedSpectralDeferredCorrections) {

        // Do this for the reactions as well, in case we cut the timestep
        // short due to it being rejected.

        MultiFab& R_old = get_old_data(Reactions_Type);
        MultiFab& R_new = get_new_data(Reactions_Type);
        MultiFab::Copy(R_new, R_old, 0, 0, R_new.nComp(), R_new.nGrow());

        // Skip the rest of the advance if the burn was unsuccessful.

        if (burn_success != 1)
            return dt;

    }
#endif

    // Construct the old-time sources from Sborder.  This will already
    // be applied to S_new (with full dt weighting), to be correctly
    // later.  Note -- this does not affect the prediction of the
    // interface state, an explict source will be traced there as
    // needed.

#ifdef SELF_GRAVITY
    construct_old_gravity(amr_iteration, amr_ncycle, prev_time);
#endif

    MultiFab& old_source = get_old_data(Source_Type);

    if (apply_sources()) {

      do_old_sources(old_source, Sborder, prev_time, dt, amr_iteration, amr_ncycle);
      apply_source_to_state(S_new, old_source, dt, 0);
      clean_state(S_new, cur_time, 0);

      // Apply the old sources to the sources for the hydro.
      // Note that we are doing an add here, not a copy,
      // in case we have already started with some source
      // terms (e.g. the source term predictor, or the SDC source).

      if (do_hydro) {
          AmrLevel::FillPatchAdd(*this, sources_for_hydro, NUM_GROW, time, Source_Type, 0, NUM_STATE);
      }

    } else {
      old_source.setVal(0.0, NUM_GROW);

    }


    // Do the hydro update.  We build directly off of Sborder, which
    // is the state that has already seen the burn

    if (do_hydro)
    {
      // Construct the primitive variables.
      cons_to_prim(time);

      // Check for CFL violations.
      check_for_cfl_violation(dt);

      // If we detect one, return immediately.
      if (cfl_violation && hard_cfl_limit)
          return dt;

      construct_ctu_hydro_source(time, dt);
      apply_source_to_state(S_new, hydro_source, dt, 0);
      clean_state(S_new, cur_time, 0);
    }


    // Sync up state after old sources and hydro source.
    frac_change = clean_state(S_new, cur_time, 0);

#ifndef AMREX_USE_CUDA
    // Check for NaN's.

    check_for_nan(S_new);
#endif

    // if we are done with the update do the source correction and
    // then the second half of the reactions

#ifdef SELF_GRAVITY
    // Must define new value of "center" before we call new gravity
    // solve or external source routine
    if (moving_center == 1)
      define_new_center(S_new, time);
#endif

#ifdef SELF_GRAVITY
    // We need to make the new radial data now so that we can use it when we
    // FillPatch in creating the new source.

#if (BL_SPACEDIM > 1)
    if ( (level == 0) && (spherical_star == 1) ) {
      int is_new = 1;
      make_radial_data(is_new);
    }
#endif
#endif

    // Construct and apply new-time source terms.

#ifdef SELF_GRAVITY
    construct_new_gravity(amr_iteration, amr_ncycle, cur_time);
#endif

    MultiFab& new_source = get_new_data(Source_Type);

    if (apply_sources()) {

      do_new_sources(new_source, Sborder, S_new, cur_time, dt, amr_iteration, amr_ncycle);
      apply_source_to_state(S_new, new_source, dt, 0);
      clean_state(S_new, cur_time, 0);

    } else {

      new_source.setVal(0.0, NUM_GROW);

    }

    // If the state has ghost zones, sync them up now
    // since the hydro source only works on the valid zones.

    if (S_new.nGrow() > 0) {
        clean_state(S_new, cur_time, 0);
        expand_state(S_new, cur_time, S_new.nGrow());
    }

#ifdef THORNADO
    create_thornado_source(dt);
#endif

    // Do the second half of the reactions.

#ifdef REACTIONS
    if (time_integration_method != SimplifiedSpectralDeferredCorrections) {

        strang_react_second_half(cur_time - 0.5 * dt, 0.5 * dt);

        // Skip the rest of the advance if the burn was unsuccessful.

        if (burn_success != 1)
            return dt;

    }
#endif

    finalize_do_advance(time, dt, amr_iteration, amr_ncycle);

    return dt;
}

bool
Castro::retry_advance_ctu(Real& time, Real dt, int amr_iteration, int amr_ncycle)
{
    BL_PROFILE("Castro::retry_advance_ctu()");

    Real dt_sub = 1.e200;

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

#ifdef REACTIONS
    MultiFab& R_old = get_old_data(Reactions_Type);
    MultiFab& R_new = get_new_data(Reactions_Type);
#endif

    const Real* dx = geom.CellSize();

    bool do_retry = false;

    // By default, we don't do a retry unless the criteria are violated.

#ifdef _OPENMP
#pragma omp parallel reduction(min:dt_sub)
#endif
    for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

#pragma gpu box(bx)
        ca_check_timestep(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                          BL_TO_FORTRAN_ANYD(S_old[mfi]),
                          BL_TO_FORTRAN_ANYD(S_new[mfi]),
#ifdef REACTIONS
                          BL_TO_FORTRAN_ANYD(R_old[mfi]),
                          BL_TO_FORTRAN_ANYD(R_new[mfi]),
#endif
                          AMREX_REAL_ANYD(dx),
                          dt, AMREX_MFITER_REDUCE_MIN(&dt_sub));

    }

    if (retry_neg_dens_factor > 0.0) {

        // Negative density criterion
        // Reset so that the desired maximum fractional change in density
        // is not larger than retry_neg_dens_factor.

        ParallelDescriptor::ReduceRealMin(frac_change);

        if (frac_change < 0.0)
            dt_sub = std::min(dt_sub, dt * -(retry_neg_dens_factor / frac_change));

    }

    ParallelDescriptor::ReduceRealMin(dt_sub);

    // Do the retry if the suggested timestep is smaller than the actual one.
    // A user-specified tolerance parameter can be used here to prevent
    // retries that are caused by small differences. Note that we are going
    // to intentionally ignore the actual suggested subcycle, and just go with
    // retry_subcycle_factor * the current timestep. The reason is that shrinking
    // the timestep by that factor will substantially change the evolution, and it
    // could be enough to get the simulation to become sane again. If this is the
    // case, we end up saving a lot of timesteps relative to the potentially very
    // small timestep recommended by the above limiters.

    if (dt_sub * (1.0 + retry_tolerance) < std::min(dt, dt_subcycle) || burn_success != 1) {

        do_retry = true;

        dt_subcycle = std::min(dt, dt_subcycle) * retry_subcycle_factor;

        if (verbose && ParallelDescriptor::IOProcessor()) {
            std::cout << std::endl;
            std::cout << "  Timestep " << dt << " rejected at level " << level << "." << std::endl;
            std::cout << "  Performing a retry, with subcycled timesteps of maximum length dt = " << dt_subcycle << std::endl;
            std::cout << std::endl;
        }

        // Restore the original values of the state data.

        for (int k = 0; k < num_state_type; k++) {

            if (prev_state[k]->hasOldData())
                state[k].copyOld(*prev_state[k]);

            if (prev_state[k]->hasNewData())
                state[k].copyNew(*prev_state[k]);

        }

        // Reset the source term predictor.

        if (do_hydro) {
            sources_for_hydro.setVal(0.0, NUM_GROW);
        }

        // Clear the contribution to the fluxes from this step.

        for (int dir = 0; dir < 3; ++dir)
            fluxes[dir]->setVal(0.0);

        for (int dir = 0; dir < 3; ++dir)
            mass_fluxes[dir]->setVal(0.0);

#if (BL_SPACEDIM <= 2)
        if (!Geom().IsCartesian())
            P_radial.setVal(0.0);
#endif

#ifdef RADIATION
        if (Radiation::rad_hydro_combined)
            for (int dir = 0; dir < BL_SPACEDIM; ++dir)
                rad_fluxes[dir]->setVal(0.0);
#endif

        if (time_integration_method == CornerTransportUpwind && source_term_predictor == 1) {

            // Normally the source term predictor is done before the swap,
            // but the prev_state data is saved after the initial swap had
            // been done. So we will temporarily swap the state data back,
            // and reset the time levels.

            // Note that unlike the initial application of the source term
            // predictor before the swap, the old data will have already
            // been allocated when we get to this point. So we want to skip
            // this step if we didn't have old data initially.

            if (prev_state_had_old_data) {

                swap_state_time_levels(0.0);

                const Real dt_old = prev_state_new_time - prev_state_old_time;

                for (int k = 0; k < num_state_type; k++)
                    state[k].setTimeLevel(prev_state_new_time, dt_old, 0.0);

                apply_source_term_predictor();

                swap_state_time_levels(0.0);

                for (int k = 0; k < num_state_type; k++)
                    state[k].setTimeLevel(time + dt_subcycle, dt_subcycle, 0.0);

            }

        }

        if (track_grid_losses)
            for (int i = 0; i < n_lost; i++)
                material_lost_through_boundary_temp[i] = 0.0;

    }

    return do_retry;

}



Real
Castro::subcycle_advance_ctu(const Real time, const Real dt, int amr_iteration, int amr_ncycle)
{
    BL_PROFILE("Castro::subcycle_advance_ctu()");

    // Start the subcycle time off with the main dt,
    // unless we already came in here with an estimate
    // that is different from the initial value we assigned,
    // for example from the post-step regrid algorithm.

    if (dt_subcycle == 1.e200)
        dt_subcycle = dt;

    Real subcycle_time = time;

    Real dt_new = 1.e200;

    sub_iteration = 0;

    // Subcycle until we've reached the target time.
    // Compare against a slightly smaller number to
    // avoid roundoff concerns.

    Real eps = 1.0e-14;

    bool do_swap = false;

    while (subcycle_time < (1.0 - eps) * (time + dt)) {

        sub_iteration += 1;

        if (dt_subcycle < dt_cutoff) {
            if (ParallelDescriptor::IOProcessor()) {
                std::cout << std::endl;
                std::cout << "  The subcycle mechanism requested subcycled timesteps of maximum length dt = " << dt_subcycle << "," << std::endl
                          << "  but this timestep is shorter than the user-defined minimum, " << std::endl
                          << "  castro.dt_cutoff = " << dt_cutoff << ". Aborting." << std::endl;
            }
            amrex::Abort("Error: subcycled timesteps too short.");
        }

        // Shorten the last timestep so that we don't overshoot
        // the ending time.

        if (subcycle_time + dt_subcycle > (time + dt))
            dt_subcycle = (time + dt) - subcycle_time;

        // Check on whether we are going to take too many subcycles.

        int num_subcycles_remaining = int(round(((time + dt) - subcycle_time) / dt_subcycle));

        if (num_subcycles_remaining > max_subcycles) {
            amrex::Print() << std::endl
                           << "  The subcycle mechanism requested " << num_subcycles_remaining << " subcycled timesteps, which is larger than the maximum of " << max_subcycles << "." << std::endl
                           << "  If you would like to override this, increase the parameter castro.max_subcycles." << std::endl;
            amrex::Abort("Error: too many subcycles.");
        }

        // If we get to this point, we survived the sanity checks. Print out the current subcycle iteration.

        if (verbose && ParallelDescriptor::IOProcessor()) {
            std::cout << std::endl;
            std::cout << "  Beginning subcycle " << sub_iteration << " starting at time " << subcycle_time
                      << " with dt = " << dt_subcycle << std::endl;
            std::cout << "  Estimated number of subcycles remaining: " << num_subcycles_remaining << std::endl << std::endl;
        }

        // Swap the time levels. Only do this after the first iteration,
        // and when we are not doing a retry (which handles the swap).

        if (do_swap) {

            // Reset the source term predictor.
            // This must come before the swap.

            if (do_hydro) {
                sources_for_hydro.setVal(0.0, NUM_GROW);
            }

            if (time_integration_method == CornerTransportUpwind && source_term_predictor == 1)
                apply_source_term_predictor();

            swap_state_time_levels(0.0);

#ifdef SELF_GRAVITY
            if (do_grav) {
                gravity->swapTimeLevels(level);
            }
#endif

        }

        // Assume we want to do a swap in the next iteration,
        // unless the retry tells us otherwise.

        do_swap = true;

        // Set the relevant time levels.

        for (int k = 0; k < num_state_type; k++)
            state[k].setTimeLevel(subcycle_time + dt_subcycle, dt_subcycle, 0.0);

        do_advance_ctu(subcycle_time, dt_subcycle, amr_iteration, amr_ncycle);

        if (verbose && ParallelDescriptor::IOProcessor()) {
            std::cout << "  Subcycle completed" << std::endl << std::endl;
        }

        subcycle_time += dt_subcycle;

        // If we have hit a CFL violation during this subcycle, we must abort.

        if (cfl_violation && hard_cfl_limit && !use_retry)
            amrex::Abort("CFL is too high at this level, and we are already inside a retry -- go back to a checkpoint and restart with lower cfl number");

        if (burn_success != 1 && !use_retry)
            amrex::Abort("Burn was unsuccessful");

        // If we're allowing for retries, check for that here.

        if (use_retry) {

            // If we hit a retry, signal that we want to try again.
            // The retry function will handle resetting the state,
            // and updating dt_subcycle.

            if (retry_advance_ctu(subcycle_time, dt_subcycle, amr_iteration, amr_ncycle)) {
                do_swap = false;
                sub_iteration = 0;
                subcycle_time = time;
                lastDtRetryLimited = true;
                lastDtFromRetry = dt_subcycle;
            }

        }

    }

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "  Subcycling complete" << std::endl << std::endl;

    if (sub_iteration > 1) {

        // Finally, copy the original data back to the old state
        // data so that externally it appears like we took only
        // a single timestep. We'll do this as a swap so that
        // we still have the last iteration's old data if we need
        // it later.

        for (int k = 0; k < num_state_type; k++) {

            if (prev_state[k]->hasOldData())
                state[k].replaceOldData(*prev_state[k]);

            state[k].setTimeLevel(time + dt, dt, 0.0);
            prev_state[k]->setTimeLevel(time + dt, dt_subcycle, 0.0);

        }

        // If we took more than one step and are going to do a reflux,
        // keep the data past the end of the step.

        if (do_reflux && update_sources_after_reflux) {

            // Note that since we only want to do this if there's actually a
            // reflux immediately following this, skip this if we're on the
            // finest level and this is not the last iteration.

            if (!(amr_iteration < amr_ncycle && level == parent->finestLevel()))
                keep_prev_state = true;

        }

    }

    // We want to return the subcycled timestep as a suggestion.

    dt_new = std::min(dt_new, dt_subcycle);

    return dt_new;

}
