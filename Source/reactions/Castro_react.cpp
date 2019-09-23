
#include "Castro.H"
#include "Castro_F.H"

#include "AMReX_DistributionMapping.H"

using std::string;
using namespace amrex;

void
Castro::strang_react_first_half(Real time, Real dt)
{
    BL_PROFILE("Castro::strang_react_first_half()");

    // Sanity check: should only be in here if we're doing CTU or MOL.

    if (time_integration_method != CornerTransportUpwind) {
        amrex::Error("Strang reactions are only supported for the CTU and MOL advance.");
    }

    // Get the reactions MultiFab to fill in.

    MultiFab& reactions = get_old_data(Reactions_Type);

    // Ensure we always have valid data, even if we don't do the burn.

    reactions.setVal(0.0);

    if (do_react != 1) return;

    // Get the current state data.

    MultiFab& state = Sborder;

    // Check if we have any zones to burn.

    if (!valid_zones_to_burn(state)) return;

    const int ng = state.nGrow();

    // Reactions are expensive and we would usually rather do a
    // communication step than burn on the ghost zones. So what we
    // will do here is create a mask that indicates that we want to
    // turn on the valid interior zones but NOT on the ghost zones
    // that are interior to the level. However, we DO want to burn on
    // the ghost zones that are on the coarse-fine interfaces, since
    // that is going to be more accurate than interpolating from
    // coarse zones. So we will not mask out those zones, and the
    // subsequent FillBoundary call will not interfere with it.

    iMultiFab& interior_mask = build_interior_boundary_mask(ng);

    // The reaction weightings.

    MultiFab* weights;

    // If we're using the custom knapsack weighting, copy state data
    // to a new MultiFab with a corresponding distribution map.

    MultiFab* state_temp;
    MultiFab* reactions_temp;
    iMultiFab* mask_temp;
    MultiFab* weights_temp;

    // Use managed arrays so that we don't have to worry about
    // deleting things at the end.

    Vector<std::unique_ptr<MultiFab> > temp_data;
    std::unique_ptr<iMultiFab> temp_idata;

    if (use_custom_knapsack_weights) {

	weights = &get_old_data(Knapsack_Weight_Type);

	// Note that we want to use the "old" data here; we've already
	// done a swap on the new data for the old data, so this is
	// really the new-time burn from the last timestep.

	const DistributionMapping& dm = DistributionMapping::makeKnapSack(get_old_data(Knapsack_Weight_Type));
	
	temp_data.push_back(std::unique_ptr<MultiFab>(
				state_temp = new MultiFab(state.boxArray(), dm, state.nComp(), state.nGrow())));
	temp_data.push_back(std::unique_ptr<MultiFab>(
				reactions_temp = new MultiFab(reactions.boxArray(), dm, reactions.nComp(), reactions.nGrow())));
	temp_data.push_back(std::unique_ptr<MultiFab>(
				weights_temp = new MultiFab(weights->boxArray(), dm, weights->nComp(), weights->nGrow())));
	temp_idata.reset(mask_temp = new iMultiFab(interior_mask.boxArray(), dm, interior_mask.nComp(), interior_mask.nGrow()));

	// Copy data from the state. Note that this is a parallel copy
	// from FabArray, and the parallel copy assumes that the data
	// on the ghost zones in state is valid and consistent with
	// the data on the interior zones, since either the ghost or
	// valid zones may end up filling a given destination zone.

	state_temp->copy(state, 0, 0, state.nComp(), state.nGrow(), state.nGrow());

	// Create the mask. We cannot use the interior_mask from above, generated by
	// Castro::build_interior_boundary_mask, because we need it to exist on the
	// current DistributionMap and a parallel copy won't work for the mask.

	int ghost_covered_by_valid = 0;
	int other_cells = 1; // uncovered ghost, valid, and outside domain cells are set to 1

	mask_temp->BuildMask(geom.Domain(), geom.periodicity(),
		             ghost_covered_by_valid, other_cells, other_cells, other_cells);

    }
    else {

	state_temp = &state;
	reactions_temp = &reactions;
	mask_temp = &interior_mask;

	// Create a dummy weight array to pass to Fortran.

	temp_data.push_back(std::unique_ptr<MultiFab>(
				weights_temp = new MultiFab(reactions.boxArray(), reactions.DistributionMap(),
							    reactions.nComp(), reactions.nGrow())));

    }

    if (verbose)
        amrex::Print() << "... Entering burner and doing half-timestep of burning." << std::endl << std::endl;

    react_state(*state_temp, *reactions_temp, *mask_temp, *weights_temp, time, dt, 1, ng);

    if (verbose)
        amrex::Print() << "... Leaving burner after completing half-timestep of burning." << std::endl << std::endl;

    // Note that this FillBoundary *must* occur before we copy any data back
    // to the main state data; it is the only way to ensure that the parallel
    // copy to follow is sensible, because when we're working with ghost zones
    // the valid and ghost zones must be consistent for the parallel copy.

    state_temp->FillBoundary(geom.periodicity());

    // Copy data back to the state data if necessary.

    if (use_custom_knapsack_weights) {

	state.copy(*state_temp, 0, 0, state_temp->nComp(), state_temp->nGrow(), state_temp->nGrow());
	reactions.copy(*reactions_temp, 0, 0, reactions_temp->nComp(), reactions_temp->nGrow(), reactions_temp->nGrow());
	weights->copy(*weights_temp, 0, 0, weights_temp->nComp(), weights_temp->nGrow(), weights_temp->nGrow());

    }

    // Ensure consistency in internal energy and recompute temperature.

    clean_state(state, time, state.nGrow());

}



void
Castro::strang_react_second_half(Real time, Real dt)
{
    BL_PROFILE("Castro::strang_react_second_half()");

    // Sanity check: should only be in here if we're doing CTU or MOL.

    if (time_integration_method != CornerTransportUpwind) {
        amrex::Error("Strang reactions are only supported for the CTU and MOL advance.");
    }

    MultiFab& reactions = get_new_data(Reactions_Type);

    // Ensure we always have valid data, even if we don't do the burn.

    reactions.setVal(0.0);

    if (Knapsack_Weight_Type > 0)
        get_new_data(Knapsack_Weight_Type).setVal(1.0);

    if (do_react != 1) return;

    MultiFab& state = get_new_data(State_Type);

    // Check if we have any zones to burn.

    if (!valid_zones_to_burn(state)) return;

    // To be consistent with other source term types,
    // we are only applying this on the interior zones.

    const int ng = 0;

    iMultiFab& interior_mask = build_interior_boundary_mask(ng);

    MultiFab* weights;

    // Most of the following uses the same logic as strang_react_first_half;
    // look there for explanatory comments.

    MultiFab* state_temp;
    MultiFab* reactions_temp;
    iMultiFab* mask_temp;
    MultiFab* weights_temp;

    Vector<std::unique_ptr<MultiFab> > temp_data;
    std::unique_ptr<iMultiFab> temp_idata;

    if (use_custom_knapsack_weights) {

	weights = &get_new_data(Knapsack_Weight_Type);

	// Here we use the old-time weights filled in during the first-half Strang-split burn.

	const DistributionMapping& dm = DistributionMapping::makeKnapSack(get_old_data(Knapsack_Weight_Type));

	temp_data.push_back(std::unique_ptr<MultiFab>(
				state_temp = new MultiFab(state.boxArray(), dm, state.nComp(), state.nGrow())));
	temp_data.push_back(std::unique_ptr<MultiFab>(
				reactions_temp = new MultiFab(reactions.boxArray(), dm, reactions.nComp(), reactions.nGrow())));
	temp_data.push_back(std::unique_ptr<MultiFab>(
				weights_temp = new MultiFab(weights->boxArray(), dm, weights->nComp(), weights->nGrow())));
	temp_idata.reset(mask_temp = new iMultiFab(interior_mask.boxArray(), dm, interior_mask.nComp(), interior_mask.nGrow()));

	state_temp->copy(state, 0, 0, state.nComp(), state.nGrow(), state.nGrow());

	int ghost_covered_by_valid = 0;
	int other_cells = 1;

	mask_temp->BuildMask(geom.Domain(), geom.periodicity(),
		             ghost_covered_by_valid, other_cells, other_cells, other_cells);

    }
    else {

	state_temp = &state;
	reactions_temp = &reactions;
	mask_temp = &interior_mask;

	temp_data.push_back(std::unique_ptr<MultiFab>(
				weights_temp = new MultiFab(reactions.boxArray(), reactions.DistributionMap(),
							    reactions.nComp(), reactions.nGrow())));

    }

    if (verbose)
        amrex::Print() << "... Entering burner and doing half-timestep of burning." << std::endl << std::endl;

    react_state(*state_temp, *reactions_temp, *mask_temp, *weights_temp, time, dt, 2, ng);

    if (verbose)
        amrex::Print() << "... Leaving burner after completing half-timestep of burning." << std::endl << std::endl;

    state_temp->FillBoundary(geom.periodicity());

    if (use_custom_knapsack_weights) {

	state.copy(*state_temp, 0, 0, state_temp->nComp(), state_temp->nGrow(), state_temp->nGrow());
	reactions.copy(*reactions_temp, 0, 0, reactions_temp->nComp(), reactions_temp->nGrow(), reactions_temp->nGrow());
	weights->copy(*weights_temp, 0, 0, weights_temp->nComp(), weights_temp->nGrow(), weights_temp->nGrow());

    }

    clean_state(state, time + 0.5 * dt, state.nGrow());

}



void
Castro::react_state(MultiFab& s, MultiFab& r, const iMultiFab& mask, MultiFab& w, Real time, Real dt_react, int strang_half, int ngrow)
{

    BL_PROFILE("Castro::react_state()");

    // Sanity check: should only be in here if we're doing CTU or MOL.

    if (time_integration_method != CornerTransportUpwind) {
        amrex::Error("Strang reactions are only supported for the CTU and MOL advance.");
    }

    const Real strt_time = ParallelDescriptor::second();

    // Initialize the weights to the default value (everything is weighted equally).

    w.setVal(1.0);

    // Start off assuming a successful burn.

    burn_success = 1;
    Real burn_failed = 0.0;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(s, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

	const Box& bx = mfi.growntilebox(ngrow);

	// Note that box is *not* necessarily just the valid region!
#pragma gpu box(bx)
	ca_react_state(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
		       BL_TO_FORTRAN_ANYD(s[mfi]),
		       BL_TO_FORTRAN_ANYD(r[mfi]),
		       BL_TO_FORTRAN_ANYD(w[mfi]),
		       BL_TO_FORTRAN_ANYD(mask[mfi]),
		       time, dt_react, strang_half,
	               AMREX_MFITER_REDUCE_SUM(&burn_failed));

    }

    if (burn_failed != 0.0) burn_success = 0;

    ParallelDescriptor::ReduceIntMin(burn_success);

    if (print_update_diagnostics) {

	Real e_added = r.sum(NumSpec + 1);

	if (e_added != 0.0)
            amrex::Print() << "... (rho e) added from burning: " << e_added << std::endl << std::endl;

    }

    if (verbose > 0)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        amrex::Print() << "Castro::react_state() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
	});
#endif
    }

}

// Simplified SDC version

void
Castro::react_state(Real time, Real dt)
{
    BL_PROFILE("Castro::react_state()");

    // Sanity check: should only be in here if we're doing simplified SDC.

    if (time_integration_method != SimplifiedSpectralDeferredCorrections) {
        amrex::Error("This react_state interface is only supported for simplified SDC.");
    }

    const Real strt_time = ParallelDescriptor::second();

    if (verbose)
        amrex::Print() << "... Entering burner and doing full timestep of burning." << std::endl << std::endl;

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    // Build the burning mask, in case the state has ghost zones.

    const int ng = S_new.nGrow();
    const iMultiFab& interior_mask = build_interior_boundary_mask(ng);

    // Create a MultiFab with all of the non-reacting source terms.

    MultiFab A_src(grids, dmap, NUM_STATE, ng);
    sum_of_sources(A_src);

    MultiFab& reactions = get_old_data(Reactions_Type);

    reactions.setVal(0.0);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
    {

	const Box& bx = mfi.growntilebox(ng);

	FArrayBox& uold    = S_old[mfi];
	FArrayBox& unew    = S_new[mfi];
	FArrayBox& a       = A_src[mfi];
	FArrayBox& r       = reactions[mfi];
	const IArrayBox& m = interior_mask[mfi];

	ca_react_state_simplified_sdc(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                                      uold.dataPtr(), ARLIM_3D(uold.loVect()), ARLIM_3D(uold.hiVect()),
                                      unew.dataPtr(), ARLIM_3D(unew.loVect()), ARLIM_3D(unew.hiVect()),
                                      a.dataPtr(), ARLIM_3D(a.loVect()), ARLIM_3D(a.hiVect()),
                                      r.dataPtr(), ARLIM_3D(r.loVect()), ARLIM_3D(r.hiVect()),
                                      m.dataPtr(), ARLIM_3D(m.loVect()), ARLIM_3D(m.hiVect()),
                                      time, dt, sdc_iteration);

    }

    if (ng > 0)
        S_new.FillBoundary(geom.periodicity());

    if (print_update_diagnostics) {

        Real e_added = reactions.sum(NumSpec + 1);

	if (e_added != 0.0)
            amrex::Print() << "... (rho e) added from burning: " << e_added << std::endl << std::endl;

    }

    if (verbose) {

        amrex::Print() << "... Leaving burner after completing full timestep of burning." << std::endl << std::endl;

        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time, IOProc);

        amrex::Print() << "Castro::react_state() time = " << run_time << std::endl << std::endl;
#ifdef BL_LAZY
	});
#endif

    }

}



bool
Castro::valid_zones_to_burn(MultiFab& State)
{

    // The default values of the limiters are 0 and 1.e200, respectively.

    Real small = 1.e-10;
    Real large = 1.e199;

    // Check whether we are limiting on either rho or T.

    bool limit_small_rho = react_rho_min >= small;
    bool limit_large_rho = react_rho_max <= large;

    bool limit_rho = limit_small_rho || limit_large_rho;

    bool limit_small_T = react_T_min >= small;
    bool limit_large_T = react_T_max <= large;

    bool limit_T = limit_small_T || limit_large_T;

    bool limit = limit_rho || limit_T;

    if (!limit) return true;

    // Now, if we're limiting on rho, collect the
    // minimum and/or maximum and compare.

    amrex::Vector<Real> small_limiters;
    amrex::Vector<Real> large_limiters;

    bool local = true;

    Real small_dens = small;
    Real large_dens = large;

    if (limit_small_rho) {
      small_dens = State.min(Density, 0, local);
      small_limiters.push_back(small_dens);
    }

    if (limit_large_rho) {
      large_dens = State.max(Density, 0, local);
      large_limiters.push_back(large_dens);
    }

    Real small_T = small;
    Real large_T = large;

    if (limit_small_T) {
      small_T = State.min(Temp, 0, local);
      small_limiters.push_back(small_T);
    }

    if (limit_large_T) {
      large_T = State.max(Temp, 0, local);
      large_limiters.push_back(large_T);
    }

    // Now do the reductions. We're being careful here
    // to limit the amount of work and communication,
    // because regularly doing this check only makes sense
    // if it is negligible compared to the amount of work
    // needed to just do the burn as normal.

    int small_size = small_limiters.size();

    if (small_size > 0) {
        amrex::ParallelDescriptor::ReduceRealMin(small_limiters.dataPtr(), small_size);

        if (limit_small_rho) {
            small_dens = small_limiters[0];
            if (limit_small_T) {
                small_T = small_limiters[1];
            }
        } else {
            small_T = small_limiters[0];
        }
    }

    int large_size = large_limiters.size();

    if (large_size > 0) {
        amrex::ParallelDescriptor::ReduceRealMax(large_limiters.dataPtr(), large_size);

        if (limit_large_rho) {
            large_dens = large_limiters[0];
            if (limit_large_T) {
                large_T = large_limiters[1];
            }
        } else {
            large_T = large_limiters[1];
        }
    }

    // Finally check on whether min <= rho <= max
    // and min <= T <= max. The defaults are small
    // and large respectively, so if the limiters
    // are not on, these checks will not be triggered.

    if (large_dens >= react_rho_min && small_dens <= react_rho_max &&
        large_T >= react_T_min && small_T <= react_T_max) {
        return true;
    }

    // If we got to this point, we did not survive the limiters,
    // so there are no zones to burn.

    if (verbose > 1)
        amrex::Print() << "  No valid zones to burn, skipping react_state()." << std::endl;

    return false;

}
