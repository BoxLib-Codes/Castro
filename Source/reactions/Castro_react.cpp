
#include "Castro.H"
#include "Castro_F.H"

using std::string;
using namespace amrex;

// Strang version

bool
Castro::react_state(MultiFab& s, MultiFab& r, Real time, Real dt)
{
    BL_PROFILE("Castro::react_state()");

    // Sanity check: should only be in here if we're doing CTU.

    if (time_integration_method != CornerTransportUpwind) {
        amrex::Error("Strang reactions are only supported for the CTU advance.");
    }

    const Real strt_time = ParallelDescriptor::second();

    // Start off by assuming a successful burn.

    int burn_success = 1;
#ifndef CXX_REACTIONS
    Real burn_failed = 0.0;
#endif

    if (do_react != 1) {

        // Ensure we always have valid data, even if we don't do the burn.
        r.setVal(0.0, r.nGrow());

        return burn_success;

    }

    // Check if we have any zones to burn.

    if (!valid_zones_to_burn(s)) {

        // Ensure we always have valid data, even if we don't do the burn.
        r.setVal(0.0, r.nGrow());

        return burn_success;

    }

    // If we're not actually doing the burn, interpolate from the level below.

    if (level > castro::reactions_max_solve_level && level > 0) {
        FillCoarsePatch(r, 0, time, Reactions_Type, 0, r.nComp(), r.nGrow());
    }

    const int ng = s.nGrow();

    if (verbose) {
        amrex::Print() << "... Entering burner and doing half-timestep of burning." << std::endl << std::endl;
    }

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
#ifdef CXX_REACTIONS
    using ReduceTuple = typename decltype(reduce_data)::Type;
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(s, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

        const Box& bx = mfi.growntilebox(ng);

        auto U = s.array(mfi);
        auto reactions = r.array(mfi);

        if (level <= castro::reactions_max_solve_level) {

#ifdef CXX_REACTIONS

            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
            {

                burn_t burn_state;

                // Initialize some data for later.

                bool do_burn = true;
                burn_state.success = true;
                Real burn_failed = 0.0_rt;

                // Don't burn on zones inside shock regions, if the relevant option is set.

#ifdef SHOCK_VAR
                if (U(i,j,k,USHK) > 0.0_rt && disable_shock_burning == 1) {
                    do_burn = false;
                }
#endif

                Real rhoInv = 1.0_rt / U(i,j,k,URHO);

                burn_state.rho = U(i,j,k,URHO);
                burn_state.T   = U(i,j,k,UTEMP);
                burn_state.e   = 0.0_rt; // Energy generated by the burn

                for (int n = 0; n < NumSpec; ++n) {
                    burn_state.xn[n] = U(i,j,k,UFS+n) * rhoInv;
                }

#if naux > 0
                for (int n = 0; n < NumAux; ++n) {
                    burn_state.aux[n] = U(i,j,k,UFX+n) * rhoInv;
                }
#endif

                // Ensure we start with no RHS or Jacobian calls registered.

                burn_state.n_rhs = 0;
                burn_state.n_jac = 0;

                // Don't burn if we're outside of the relevant (rho, T) range.

                if (burn_state.T < castro::react_T_min || burn_state.T > castro::react_T_max ||
                    burn_state.rho < castro::react_rho_min || burn_state.rho > castro::react_rho_max) {
                    do_burn = false;
                }

                if (do_burn) {
                    burner(burn_state, dt);
                }

                // If we were unsuccessful, update the failure count.

                if (!burn_state.success) {
                    burn_failed = 1.0_rt;
                }

                if (do_burn) {

                    // Add burning rates to reactions MultiFab, but be
                    // careful because the reactions and state MFs may
                    // not have the same number of ghost cells.

                    if (reactions.contains(i,j,k)) {
                        for (int n = 0; n < NumSpec; ++n) {
                            reactions(i,j,k,n) = U(i,j,k,URHO) * (burn_state.xn[n] - U(i,j,k,UFS+n) * rhoInv) / dt;
                        }
#if naux > 0
                        for (int n = 0; n < NumAux; ++n) {
                            reactions(i,j,k,n+NumSpec) = U(i,j,k,URHO) * (burn_state.aux[n] - U(i,j,k,UFX+n) * rhoInv) / dt;
                        }
#endif
                        reactions(i,j,k,NumSpec+NumAux  ) = U(i,j,k,URHO) * burn_state.e / dt;
                        reactions(i,j,k,NumSpec+NumAux+1) = amrex::max(1.0_rt, static_cast<Real>(burn_state.n_rhs + 2 * burn_state.n_jac));
                    }

                }
                else {

                    if (reactions.contains(i,j,k)) {
                        for (int n = 0; n < NumSpec + NumAux + 1; ++n) {
                            reactions(i,j,k,n) = 0.0_rt;
                        }

                        reactions(i,j,k,NumSpec+NumAux+1) = 1.0_rt;
                    }

                }

                return {burn_failed};

            });

#else

#pragma gpu box(bx)
            ca_react_state(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                           BL_TO_FORTRAN_ANYD(s[mfi]),
                           BL_TO_FORTRAN_ANYD(r[mfi]),
                           time, dt,
                           AMREX_MFITER_REDUCE_SUM(&burn_failed));

#endif

        }

        // Now update the state with the reactions data.

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
        {
            if (U.contains(i,j,k) && reactions.contains(i,j,k)) {
                for (int n = 0; n < NumSpec; ++n) {
                    U(i,j,k,UFS+n) += reactions(i,j,k,n) * dt;
                }
#if naux > 0
                for (int n = 0; n < NumAux; ++n) {
                    U(i,j,k,UFX+n) += reactions(i,j,k,n+NumSpec) * dt;
                }
#endif
                U(i,j,k,UEINT) += reactions(i,j,k,NumSpec+NumAux) * dt;
                U(i,j,k,UEDEN) += reactions(i,j,k,NumSpec+NumAux) * dt;
            }
        });

    }

#ifdef CXX_REACTIONS
    ReduceTuple hv = reduce_data.value();
    Real burn_failed = amrex::get<0>(hv);
#endif

    if (burn_failed != 0.0) {
      burn_success = 0;
    }

    ParallelDescriptor::ReduceIntMin(burn_success);

    if (print_update_diagnostics) {

        Real e_added = r.sum(NumSpec + 1);

        if (e_added != 0.0) {
            amrex::Print() << "... (rho e) added from burning: " << e_added << std::endl << std::endl;
        }

    }

    if (verbose) {
        amrex::Print() << "... Leaving burner after completing half-timestep of burning." << std::endl << std::endl;
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

    return burn_success;

}

// Simplified SDC version

bool
Castro::react_state(Real time, Real dt)
{
    BL_PROFILE("Castro::react_state()");

    // Sanity check: should only be in here if we're doing simplified SDC.

    if (time_integration_method != SimplifiedSpectralDeferredCorrections) {
        amrex::Error("This react_state interface is only supported for simplified SDC.");
    }

    const Real strt_time = ParallelDescriptor::second();

    if (verbose) {
        amrex::Print() << "... Entering burner and doing full timestep of burning." << std::endl << std::endl;
    }

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    // Build the burning mask, in case the state has ghost zones.

    const int ng = S_new.nGrow();
    const iMultiFab& interior_mask = build_interior_boundary_mask(ng);

    // Create a MultiFab with all of the non-reacting source terms.

    MultiFab A_src(grids, dmap, NUM_STATE, ng);
    sum_of_sources(A_src);

    MultiFab& reactions = get_new_data(Reactions_Type);

    reactions.setVal(0.0, reactions.nGrow());

    // Start off assuming a successful burn.

    int burn_success = 1;
    Real burn_failed = 0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:burn_failed)
#endif
    for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

        const Box& bx = mfi.growntilebox(ng);

        FArrayBox& uold    = S_old[mfi];
        FArrayBox& unew    = S_new[mfi];
        FArrayBox& a       = A_src[mfi];
        FArrayBox& r       = reactions[mfi];
        const IArrayBox& m = interior_mask[mfi];

#pragma gpu box(bx)
        ca_react_state_simplified_sdc(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                                      BL_TO_FORTRAN_ANYD(uold),
                                      BL_TO_FORTRAN_ANYD(unew),
                                      BL_TO_FORTRAN_ANYD(a),
                                      BL_TO_FORTRAN_ANYD(r),
                                      BL_TO_FORTRAN_ANYD(m),
                                      time, dt, sdc_iteration,
                                      AMREX_MFITER_REDUCE_SUM(&burn_failed));

    }

    if (burn_failed != 0.0) burn_success = 0;

    ParallelDescriptor::ReduceIntMin(burn_success);

    if (ng > 0) {
        S_new.FillBoundary(geom.periodicity());
    }

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

    if (burn_success) {
        return true;
    } else {
        return false;
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

    if (!limit) {
      return true;
    }

    // Now, if we're limiting on rho, collect the
    // minimum and/or maximum and compare.

    amrex::Vector<Real> small_limiters;
    amrex::Vector<Real> large_limiters;

    bool local = true;

    Real smalldens = small;
    Real largedens = large;

    if (limit_small_rho) {
      smalldens = State.min(URHO, 0, local);
      small_limiters.push_back(smalldens);
    }

    if (limit_large_rho) {
      largedens = State.max(URHO, 0, local);
      large_limiters.push_back(largedens);
    }

    Real small_T = small;
    Real large_T = large;

    if (limit_small_T) {
      small_T = State.min(UTEMP, 0, local);
      small_limiters.push_back(small_T);
    }

    if (limit_large_T) {
      large_T = State.max(UTEMP, 0, local);
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
            smalldens = small_limiters[0];
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
            largedens = large_limiters[0];
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

    if (largedens >= react_rho_min && smalldens <= react_rho_max &&
        large_T >= react_T_min && small_T <= react_T_max) {
        return true;
    }

    // If we got to this point, we did not survive the limiters,
    // so there are no zones to burn.

    if (verbose > 1) {
        amrex::Print() << "  No valid zones to burn, skipping react_state()." << std::endl;
    }

    return false;

}
