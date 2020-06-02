
#include "Castro.H"
#include "Castro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

#ifdef GRAVITY
#include "Gravity.H"
#endif

#ifdef DIFFUSION
#include "Diffusion.H"
#endif

#include <climits>
#include <cmath>

using std::string;
using namespace amrex;

#ifndef MHD
#ifndef AMREX_USE_CUDA
Real Castro::do_advance_sdc(Real time, Real dt, int amr_iteration, int amr_ncycle) {

    // this is the new "formal" SDC integration routine.

    // this does the entire update in time for 1 SDC iteration.

    BL_PROFILE("Castro::do_advance_sdc()");

    const Real prev_time = state[State_Type].prevTime();
    const Real cur_time = state[State_Type].curTime();

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    auto domain_lo = geom.Domain().loVect3d();
    auto domain_hi = geom.Domain().hiVect3d();

    // Perform initialization steps.

    initialize_do_advance(time, dt, amr_iteration, amr_ncycle);

    // Check for NaN's.

    check_for_nan(S_old);

    MultiFab& old_source = get_old_data(Source_Type);
    MultiFab& new_source = get_new_data(Source_Type);

    bool apply_sources_to_state = false;

    // we loop over all nodes, even the last, since we need to compute
    // the advective update source at each node

    for (int m = 0; m < SDC_NODES; m++) {

        current_sdc_node = m;

        // k_new represents carries the solution.  Coming into here, it
        // will be entirely the old state, but we update it on each time
        // node in place.

        Real node_time = time + dt_sdc[m] * dt;

        // fill Sborder with the starting node's info -- we use S_new as
        // our staging area.  Note we need to pass new_time here to the
        // FillPatch so it only pulls from the new MF -- this will not
        // work for multilevel.
        MultiFab::Copy(S_new, *(k_new[m]), 0, 0, S_new.nComp(), 0);
        clean_state(S_new, cur_time, 0);
        expand_state(Sborder, cur_time, NUM_GROW);

        // the next chunk of code constructs the advective term for the
        // current node, m.  First we get the sources, then full hydro
        // source term from the MOL driver.  Note, for m = 0, we only have
        // to do all of this the first iteration, since that state never
        // changes
        if (!(sdc_iteration > 0 && m == 0) &&
            !(sdc_iteration == sdc_order + sdc_extra - 1 && m == SDC_NODES - 1)) {

            // Construct the "old-time" sources from Sborder.  Since we are
            // working from Sborder, this will actually evaluate the sources
            // using the current stage's starting point.

            // TODO: this is not using the density at the current stage
#ifdef GRAVITY
            construct_old_gravity(amr_iteration, amr_ncycle, prev_time);
#endif

            if (apply_sources()) {
#ifndef AMREX_USE_CUDA
                if (sdc_order == 4) {
                    // if we are 4th order, convert to cell-center Sborder -> Sborder_cc
                    // we'll use Sburn for this memory buffer at the moment

                    for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
                        const Box& gbx = mfi.growntilebox(1);

                        make_cell_center(gbx, Sborder.array(mfi), Sburn.array(mfi), domain_lo,
                                         domain_hi);
                    }

                    // we pass in the stage time here
                    do_old_sources(old_source, Sburn, Sburn, node_time, dt, apply_sources_to_state,
                                   amr_iteration, amr_ncycle);

                    // fill the ghost cells for the sources -- note since we have
                    // not defined the new_source yet, we either need to copy this
                    // into new_source for the time-interpolation in the ghost
                    // fill to make sense, or so long as we are not multilevel,
                    // just use the old time (prev_time) in the fill instead of
                    // the node time (time)
                    AmrLevel::FillPatch(*this, old_source, old_source.nGrow(), prev_time,
                                        Source_Type, 0, NSRC);

                    // Now convert to cell averages.  This loop cannot be tiled.
                    FArrayBox tmp;

                    for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
                        const Box& bx = mfi.tilebox();

                        tmp.resize(bx, 1);
                        Elixir elix_tmp = tmp.elixir();
                        auto tmp_arr = tmp.array();

                        make_fourth_in_place(bx, old_source.array(mfi), tmp_arr, domain_lo,
                                             domain_hi);
                    }

                } else {
                    // there is a ghost cell fill hidden in diffusion, so we need
                    // to pass in the time associate with Sborder
                    do_old_sources(old_source, Sborder, Sborder, cur_time, dt,
                                   apply_sources_to_state, amr_iteration, amr_ncycle);
                }

                // note: we don't need a FillPatch on the sources, since they
                // are only used in the valid box in the conservative flux
                // update construction.  The only exception is if we are doing
                // the well-balanced method in the reconstruction of the
                // pressure.
                if (sdc_order == 2 && use_pslope == 1) {
                    AmrLevel::FillPatch(*this, old_source, old_source.nGrow(), prev_time,
                                        Source_Type, 0, NSRC);
                }
#endif

                // store the result in sources_for_hydro -- this is what will
                // be used in the final conservative update
                MultiFab::Copy(sources_for_hydro, old_source, 0, 0, NSRC, old_source.nGrow());

            } else {
                sources_for_hydro.setVal(0.0, 0);
            }

            // Now compute the advective term for the current node -- this
            // will be used to advance us to the next node the new time

            // Construct the primitive variables.
            if (sdc_order == 4) {
                cons_to_prim_fourth(time);
            } else {
                cons_to_prim(time);
            }

            if (do_hydro) {
                // Check for CFL violations.
                check_for_cfl_violation(dt);

                // If we detect one, return immediately.
                if (cfl_violation) return dt;
            }

            // construct the update for the current stage -- this fills
            // A_new[m] with the righthand side for this stage.
            A_new[m]->setVal(0.0);
            construct_mol_hydro_source(time, dt, *A_new[m]);

        }  // end of the m = 0 sdc_iter > 0 check

        // also, if we are the first SDC iteration, we haven't yet stored
        // any old advective terms, so we cannot yet do the quadrature
        // over nodes.  Initialize those now.  Recall, A_new[0] and A_old[0]
        // are aliased.
        if (sdc_iteration == 0 && m == 0) {
            for (int n = 1; n < SDC_NODES; n++) {
                MultiFab::Copy(*(A_old[n]), *(A_new[0]), 0, 0, NUM_STATE, 0);
            }

#ifdef REACTIONS
            // if this is the first node of a new iteration, then we need
            // to compute and store the old reactive source

            // we already have the node state with ghost cells in Sborder,
            // so we can just use that as the starting point
            bool input_is_average = true;
            construct_old_react_source(Sborder, *(R_old[0]), input_is_average);

            // copy to the other nodes -- since the state is the same on all
            // nodes for sdc_iteration == 0
            for (int n = 1; n < SDC_NODES; n++) {
                MultiFab::Copy(*(R_old[n]), *(R_old[0]), 0, 0, R_old[0]->nComp(), 0);
            }
#endif
        }

        // update to the next stage -- this involves computing the
        // integral over the k-1 iteration data.  Note we don't do
        // this if we are on the final node (since there is nothing to
        // update to
        if (m < SDC_NODES - 1) {

            amrex::Print() << "... doing the SDC update, iteration = " << sdc_iteration
                           << " from node " << m << " to " << m + 1 << std::endl;

            do_sdc_update(m, m + 1, dt);  //(dt_sdc[m+1] - dt_sdc[m])*dt);

            // we now have a new value of k_new[m+1], do a clean_state on it
            clean_state(S_new, cur_time, 0);
        }

    }  // node iteration

    // we are done with the integration over all nodes for this iteration

    // note: as of right now, S_new is not updated with the state from
    // the final time node.  This means we can still use S_new as
    // "scratch" until we finally set it.

    if (sdc_iteration != sdc_order + sdc_extra - 1) {
        // store A_old for the next SDC iteration -- don't need to do n=0,
        // since that is unchanged
        for (int n = 1; n < SDC_NODES; n++) {
            MultiFab::Copy(*(A_old[n]), *(A_new[n]), 0, 0, NUM_STATE, 0);
        }
    }

#ifdef REACTIONS
    // we just did the update, so now recompute the "old" reactive
    // source for the next SDC iteration.  We don't need to do this for
    // m = 0, since that state never changes.

    for (int m = 1; m < SDC_NODES; ++m) {
        // TODO: do we need a clean state here?
        MultiFab::Copy(S_new, *(k_new[m]), 0, 0, S_new.nComp(), 0);
        expand_state(Sburn, cur_time, 2);
        bool input_is_average = true;
        construct_old_react_source(Sburn, *(R_old[m]), input_is_average);
    }
#endif

    if (sdc_iteration == sdc_order + sdc_extra - 1) {

        // store the new solution
        MultiFab::Copy(S_new, *(k_new[SDC_NODES - 1]), 0, 0, S_new.nComp(), 0);

        // We need to make source_old and source_new be the source terms at
        // the old and new time.  we never actually evaluated the sources
        // using the new time state (since we just constructed it).  Note:
        // we always use do_old_sources here, since we want the actual
        // source and not a correction.

        // note: we need to have ghost cells here cause some sources (in
        // particular pdivU) need them.  Perhaps it would be easier to just
        // always require State_Type to have 1 ghost cell?

        // TODO: we also need to make these 4th order!
        clean_state(S_old, prev_time, 0);
        expand_state(Sborder, prev_time, Sborder.nGrow());
        do_old_sources(old_source, Sborder, Sborder, prev_time, dt, apply_sources_to_state,
                       amr_iteration, amr_ncycle);
        AmrLevel::FillPatch(*this, old_source, old_source.nGrow(), prev_time, Source_Type, 0, NSRC);

        clean_state(S_new, cur_time, 0);
        expand_state(Sborder, cur_time, Sborder.nGrow());
        do_old_sources(new_source, Sborder, Sborder, cur_time, dt, apply_sources_to_state,
                       amr_iteration, amr_ncycle);
        AmrLevel::FillPatch(*this, new_source, new_source.nGrow(), cur_time, Source_Type, 0, NSRC);
    }

    finalize_do_advance(time, dt, amr_iteration, amr_ncycle);

#ifdef REACTIONS
    // store the reaction information as well.  Note: this will be
    // the instantaneous reactive source from the last burn.  In the
    // future, we might want to do a quadrature over R_old[]

    // At this point, Sburn contains the cell-center reaction source
    // on one ghost-cell.  So we can use this to derive what we need.

    // this is done only for the plotfile
    MultiFab& R_new = get_new_data(Reactions_Type);

    if (sdc_order == 4) {
        // fill ghost cells on S_new -- we'll need these to convert to
        // centers
        Real cur_time = state[State_Type].curTime();
        // we'll use Sborder to expand the state, but we already cleared
        // it at the end of the andance
        Sborder.define(grids, dmap, NUM_STATE, NUM_GROW, MFInfo().SetTag("Sborder"));

        expand_state(Sborder, cur_time, 2);
    }

    FArrayBox U_center;
    FArrayBox R_center;
    FArrayBox tmp;

    // this cannot be tiled
    for (MFIter mfi(R_new); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const Box& obx = mfi.growntilebox(1);

        if (sdc_order == 4) {

            // convert S_new to cell-centers
            U_center.resize(obx, NUM_STATE);
            Elixir elix_u_center = U_center.elixir();
            auto const U_center_arr = U_center.array();

            make_cell_center(obx, Sborder.array(mfi), U_center_arr, domain_lo, domain_hi);

            // pass in the reaction source and state at centers, including one ghost cell
            // and derive everything that is needed including 1 ghost cell
            R_center.resize(obx, R_new.nComp());
            Elixir elix_r_center = R_center.elixir();
            auto const R_center_arr = R_center.array();

            Array4<const Real> const Sburn_arr = Sburn.array(mfi);

            // we don't worry about the difference between centers and averages
            ca_store_reaction_state(obx, Sburn_arr, U_center_arr, R_center_arr);

            // convert R_new from centers to averages in place
            tmp.resize(bx, 1);
            Elixir elix_tmp = tmp.elixir();
            auto const tmp_arr = tmp.array();

            make_fourth_in_place(bx, R_center_arr, tmp_arr, domain_lo, domain_hi);

            // store
            R_new[mfi].copy(R_center, bx, 0, bx, 0, R_new.nComp());

        } else {

            Array4<const Real> const R_old_arr = R_old[SDC_NODES - 1]->array(mfi);
            Array4<const Real> const S_new_arr = S_new.array(mfi);
            Array4<Real> const R_new_arr = R_new.array(mfi);
            // we don't worry about the difference between centers and averages
            ca_store_reaction_state(bx, R_old_arr, S_new_arr, R_new_arr);
        }
    }

    if (sdc_order == 4) {
        Sborder.clear();
    }

#endif  // REACTIONS

    return dt;
}

#endif
#endif
