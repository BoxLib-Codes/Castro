
#include "Castro.H"
#include "Castro_F.H"

#include "AMReX_DistributionMapping.H"

using std::string;
using namespace amrex;

#ifndef SDC

void
Castro::strang_react_first_half(Real time, Real dt)
{

    // Get the reactions MultiFab to fill in.

    MultiFab& reactions = get_old_data(Reactions_Type);

    reactions.setVal(0.0);

    if (do_react != 1) return;

    // Get the current state data.

    MultiFab& state = Sborder;

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

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "\n" << "... Entering burner and doing half-timestep of burning." << "\n";

    react_state(*state_temp, *reactions_temp, *mask_temp, *weights_temp, time, dt, 1, ng);

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... Leaving burner after completing half-timestep of burning." << "\n";

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

    clean_state(state);

}



void
Castro::strang_react_second_half(Real time, Real dt)
{

    MultiFab& reactions = get_new_data(Reactions_Type);

    reactions.setVal(0.0);

    if (do_react != 1) return;

    MultiFab& state = get_new_data(State_Type);

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

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "\n" << "... Entering burner and doing half-timestep of burning." << "\n";

    react_state(*state_temp, *reactions_temp, *mask_temp, *weights_temp, time, dt, 2, ng);

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... Leaving burner after completing half-timestep of burning." << "\n";

    state_temp->FillBoundary(geom.periodicity());

    if (use_custom_knapsack_weights) {

	state.copy(*state_temp, 0, 0, state_temp->nComp(), state_temp->nGrow(), state_temp->nGrow());
	reactions.copy(*reactions_temp, 0, 0, reactions_temp->nComp(), reactions_temp->nGrow(), reactions_temp->nGrow());
	weights->copy(*weights_temp, 0, 0, weights_temp->nComp(), weights_temp->nGrow(), weights_temp->nGrow());

    }

    int is_new = 1;
    clean_state(is_new);

}



void
Castro::react_state(MultiFab& s, MultiFab& r, const iMultiFab& mask, MultiFab& w, Real time, Real dt_react, int strang_half, int ngrow)
{

    BL_PROFILE("Castro::react_state()");

    const Real strt_time = ParallelDescriptor::second();

    // Initialize the weights to the default value (everything is weighted equally).

    w.setVal(1.0);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(s, true); mfi.isValid(); ++mfi)
    {

	const Box& bx = mfi.growntilebox(ngrow);

	// Note that box is *not* necessarily just the valid region!
	ca_react_state(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		       BL_TO_FORTRAN_3D(s[mfi]),
		       BL_TO_FORTRAN_3D(r[mfi]),
		       BL_TO_FORTRAN_3D(w[mfi]),
		       BL_TO_FORTRAN_3D(mask[mfi]),
		       time, dt_react, strang_half);

    }

    if (verbose) {

	Real e_added = r.sum(NumSpec + 1);

	if (ParallelDescriptor::IOProcessor() && e_added != 0.0)
	    std::cout << "... (rho e) added from burning: " << e_added << std::endl;

    }

    if (verbose > 0)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

	if (ParallelDescriptor::IOProcessor())
	  std::cout << "Castro::react_state() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
	});
#endif
    }

}

#else

// SDC version

void
Castro::react_state(Real time, Real dt)
{
    BL_PROFILE("Castro::react_state()");

    const Real strt_time = ParallelDescriptor::second();

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "\n" << "... Entering burner and doing full timestep of burning." << "\n";

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

	ca_react_state(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		       uold.dataPtr(), ARLIM_3D(uold.loVect()), ARLIM_3D(uold.hiVect()),
		       unew.dataPtr(), ARLIM_3D(unew.loVect()), ARLIM_3D(unew.hiVect()),
		       a.dataPtr(), ARLIM_3D(a.loVect()), ARLIM_3D(a.hiVect()),
		       r.dataPtr(), ARLIM_3D(r.loVect()), ARLIM_3D(r.hiVect()),
		       m.dataPtr(), ARLIM_3D(m.loVect()), ARLIM_3D(m.hiVect()),
		       time, dt, sdc_iteration);

    }

    if (ng > 0)
        S_new.FillBoundary(geom.periodicity());

    if (verbose) {

        Real e_added = reactions.sum(NumSpec + 1);

	if (ParallelDescriptor::IOProcessor() && e_added != 0.0)
	    std::cout << "... (rho e) added from burning: " << e_added << std::endl;

	if (ParallelDescriptor::IOProcessor())
	    std::cout << "... Leaving burner after completing full timestep of burning." << "\n";

        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time, IOProc);

	if (ParallelDescriptor::IOProcessor())
	  std::cout << "Castro::react_state() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
	});
#endif

    }

}

#endif
