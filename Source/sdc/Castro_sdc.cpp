#include <Castro.H>
#include <Castro_F.H>
// #include <Castro_hydro_F.H>
#include <Castro_sdc_util.H>

using namespace amrex;

void
Castro::do_sdc_update(int m_start, int m_end, Real dt)
{

    BL_PROFILE("Castro::do_sdc_update()");

    // NOTE: dt here is the full dt not the dt between time nodes

    // this routine needs to do the update from time node m to m+1
    //
    // We come in with:
    //   A_new[m_start] : this is the advective update at node m_start
    //   A_old[:] : this is the advective source for all nodes at the old iterate
    //
    //   R_old[:] : this is the reaction source for all nodes at the old iterate

    // If we do advection only, then the update is explicit.  If we do
    // reactions, then the update is implicit within a zone.

    // for 4th order reactive SDC, we need to first compute the source, C
    // and do a ghost cell fill on it

    auto domain_lo = geom.Domain().loVect3d();
    auto domain_hi = geom.Domain().hiVect3d();

    // the timestep from m to m+1
    Real dt_m = (dt_sdc[m_end] - dt_sdc[m_start]) * dt;

#ifdef REACTIONS
    // SDC_Source_Type is only defined for 4th order
    MultiFab tmp;
    MultiFab& C_source = (sdc_order == 4) ? get_new_data(SDC_Source_Type) : tmp;

    if (sdc_order == 4)
    {

        // for 4th order reacting flow, we need to create the "source" C
        // as averages and then convert it to cell centers.  The cell-center
        // version needs to have 2 ghost cells
        for (MFIter mfi(*k_new[0]); mfi.isValid(); ++mfi)
        {

            const Box& bx = mfi.tilebox();
            Array4<Real> const& C_source_arr=C_source.array(mfi);

            Array4<const Real> const& A_new_arr=(A_new[m_start])->array(mfi);
            Array4<const Real> const& A_old_0_arr=(A_old[0])->array(mfi);
            Array4<const Real> const& A_old_1_arr=(A_old[1])->array(mfi);
            Array4<const Real> const& A_old_2_arr=(A_old[2])->array(mfi);
            Array4<const Real> const& R_old_0_arr=(R_old[0])->array(mfi);
            Array4<const Real> const& R_old_1_arr=(R_old[1])->array(mfi);
            Array4<const Real> const& R_old_2_arr=(R_old[2])->array(mfi);
            if (sdc_quadrature == 0)
            {

                ca_sdc_compute_C4_lobatto(bx, dt_m, dt, A_new_arr, A_old_0_arr, A_old_1_arr,
                                          A_old_2_arr, R_old_0_arr, R_old_1_arr, R_old_2_arr,
                                          C_source_arr, m_start);

            }
            else
            {

                Array4<const Real> const& A_old_3_arr=(A_old[3])->array(mfi);
                Array4<const Real> const& R_old_3_arr=(R_old[3])->array(mfi);

                ca_sdc_compute_C4_radau(bx, dt_m, dt, A_new_arr, A_old_0_arr, A_old_1_arr,
                                        A_old_2_arr,
                                        A_old_3_arr, R_old_0_arr, R_old_1_arr, R_old_2_arr, R_old_3_arr,
                                        C_source_arr, m_start);

            }
        }

        // need to construct the time for this stage -- but it is not really
        // at a single instance in time.  For single level this does not matter,
        Real time = state[SDC_Source_Type].curTime();
        AmrLevel::FillPatch(*this, C_source, C_source.nGrow(), time,
                            SDC_Source_Type, 0, NUM_STATE);


        // we'll also construct an initial guess for the nonlinear solve,
        // and store this in the Sburn MultiFab.  We'll use S_new as the
        // staging place so we can do a FillPatch
        MultiFab& S_new = get_new_data(State_Type);

        for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {

            const Box& bx = mfi.tilebox();

            Array4<const Real> const& k_new_m_start_arr=
                (k_new[m_start])->array(mfi);
            Array4<const Real> const& k_new_m_end_arr=(k_new[m_end])->array(
                                                                        mfi);
            Array4<const Real> const& A_old_arr=(A_old[m_start])->array(mfi);
            Array4<const Real> const& R_old_arr=(R_old[m_start])->array(mfi);
            Array4<Real> const& S_new_arr=S_new.array(mfi);

            ca_sdc_compute_initial_guess(bx, k_new_m_start_arr, k_new_m_end_arr,
                                         A_old_arr, R_old_arr, S_new_arr,
                                         dt_m, sdc_iteration);


        }

        const Real cur_time = state[State_Type].curTime();
        expand_state(Sburn, cur_time, 2);

    }
#endif

    // main update loop -- we are updating k_new[m_start] to
    // k_new[m_end]

    FArrayBox U_center;
    FArrayBox C_center;
    FArrayBox U_new_center;
    FArrayBox R_new;
    FArrayBox tlap;

    FArrayBox C2;

    for (MFIter mfi(*k_new[0]); mfi.isValid(); ++mfi)
    {

        const Box& bx = mfi.tilebox();
        const Box& bx1 = mfi.growntilebox(1);

#ifdef REACTIONS
        // advection + reactions
        if (sdc_order == 2)
        {

            // second order SDC reaction update -- we don't care about
            // the difference between cell-centers and averages

            // first compute the source term, C -- this differs depending
            // on whether we are Lobatto or Radau
            C2.resize(bx, NUM_STATE);
            Elixir elix_C2 = C2.elixir();
            Array4<Real> const& C2_arr=C2.array();

            Array4<const Real> const& A_new_arr=(A_new[m_start])->array(mfi);
            Array4<const Real> const& A_old_0_arr=(A_old[0])->array(mfi);
            Array4<const Real> const& A_old_1_arr=(A_old[1])->array(mfi);
            Array4<const Real> const& R_old_0_arr=(R_old[0])->array(mfi);
            Array4<const Real> const& R_old_1_arr=(R_old[1])->array(mfi);

            if (sdc_quadrature == 0)
            {

                ca_sdc_compute_C2_lobatto(bx, dt_m, dt, A_new_arr, A_old_0_arr, A_old_1_arr,
                                          R_old_0_arr, R_old_1_arr, C2_arr, m_start);

            }
            else
            {

                Array4<const Real> const& A_old_2_arr=(A_old[2])->array(mfi);
                Array4<const Real> const& R_old_2_arr=(R_old[2])->array(mfi);
                ca_sdc_compute_C2_radau(bx, dt_m, dt, A_new_arr, A_old_0_arr, A_old_1_arr,
                                        A_old_2_arr,
                                        R_old_0_arr, R_old_1_arr, R_old_2_arr, C2_arr, m_start);

            }

            auto k_m = (*k_new[m_start]).array(mfi);
            auto k_n = (*k_new[m_end]).array(mfi);
            auto A_m = (*A_new[m_start]).array(mfi);
            auto A_n = (*A_new[m_end]).array(mfi);
            auto R_m = (*R_old[m_start]).array(mfi);
            auto C_arr = C2.array();

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
            {
                sdc_update_o2(i, j, k, k_m, k_n, A_m, A_n, C_arr, dt_m, sdc_iteration, m_start);

                normalize_species_sdc(i, j, k, k_n);
            });
        }
        else
        {

            // fourth order SDC reaction update -- we need to respect the
            // difference between cell-centers and averages

            Array4<const Real> const& k_new_m_start_arr=
                (k_new[m_start])->array(mfi);
            Array4<Real> const& k_new_m_end_arr=(k_new[m_end])->array(mfi);
            Array4<const Real> const& C_source_arr=C_source.array(mfi);

            // convert the starting U to cell-centered on a fab-by-fab basis
            // -- including one ghost cell
            U_center.resize(bx1, NUM_STATE);
            Elixir elix_u_center = U_center.elixir();
            auto U_center_arr = U_center.array();

            make_cell_center(bx1, Sborder.array(mfi), U_center_arr, domain_lo, domain_hi);

            // sometimes the Laplacian can make the species go negative near discontinuities
            amrex::ParallelFor(bx1,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
            {
                normalize_species_sdc(i, j, k, U_center_arr);
            });

            // convert the C source to cell-centers
            C_center.resize(bx1, NUM_STATE);
            Elixir elix_c_center = C_center.elixir();
            auto C_center_arr = C_center.array();

            make_cell_center(bx1, C_source.array(mfi), C_center_arr, domain_lo, domain_hi);

            // solve for the updated cell-center U using our cell-centered C -- we
            // need to do this with one ghost cell
            U_new_center.resize(bx1, NUM_STATE);
            Elixir elix_u_new_center = U_new_center.elixir();
            auto U_new_center_arr = U_new_center.array();

            // initialize U_new with our guess for the new state, stored as
            // an average in Sburn
            make_cell_center(bx1, Sburn.array(mfi), U_new_center_arr, domain_lo, domain_hi);

            amrex::ParallelFor(bx1,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
            {
                sdc_update_centers_o4(i, j, k, U_center_arr, U_new_center_arr, C_center_arr, dt_m, sdc_iteration);
            });

            // compute R_i and in 1 ghost cell and then convert to <R> in
            // place (only for the interior)
            R_new.resize(bx1, NUM_STATE);
            Elixir elix_R_new = R_new.elixir();
            Array4<Real> const& R_new_arr = R_new.array();

            // ca_instantaneous_react(BL_TO_FORTRAN_BOX(bx1),
            //                        BL_TO_FORTRAN_3D(U_new_center),
            //                        BL_TO_FORTRAN_3D(R_new));
            amrex::ParallelFor(bx1,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
            {
                instantaneous_react(i, j, k, U_new_center_arr, R_new_arr);
            });

            tlap.resize(bx, 1);
            Elixir elix_tlap = tlap.elixir();
            auto const tlap_arr = tlap.array();

            make_fourth_in_place(bx, R_new_arr, tlap_arr, domain_lo, domain_hi);

            // now do the conservative update using this <R> to get <U>
            // We'll also need to pass in <C>
            ca_sdc_conservative_update(bx, dt_m, k_new_m_start_arr, k_new_m_end_arr,
                                       C_source_arr, R_new_arr);

        }
#else
        Array4<const Real> const& k_new_m_start_arr=
            (k_new[m_start])->array(mfi);
        Array4<Real> const& k_new_m_end_arr=(k_new[m_end])->array(mfi);
        Array4<const Real> const& A_new_arr=(A_new[m_start])->array(mfi);
        Array4<const Real> const& A_old_0_arr=(A_old[0])->array(mfi);
        Array4<const Real> const& A_old_1_arr=(A_old[1])->array(mfi);
        // pure advection
        if (sdc_order == 2)
        {

            if (sdc_quadrature == 0)
            {
                ca_sdc_update_advection_o2_lobatto(bx, dt_m, dt, k_new_m_start_arr,
                                                   k_new_m_end_arr,
                                                   A_new_arr, A_old_0_arr, A_old_1_arr,
                                                   m_start);

            }
            else
            {
                Array4<const Real> const& A_old_2_arr=(A_old[2])->array(mfi);
                ca_sdc_update_advection_o2_radau(bx, dt_m, dt, k_new_m_start_arr,
                                                 k_new_m_end_arr,
                                                 A_new_arr, A_old_0_arr, A_old_1_arr, A_old_2_arr,
                                                 m_start);

            }

        }
        else
        {
            Array4<const Real> const& A_old_2_arr=(A_old[2])->array(mfi);
            if (sdc_quadrature == 0)
            {
                ca_sdc_update_advection_o4_lobatto(bx, dt_m, dt, k_new_m_start_arr,
                                                   k_new_m_end_arr,
                                                   A_new_arr, A_old_0_arr, A_old_1_arr, A_old_2_arr,
                                                   m_start);

            }
            else
            {
                Array4<const Real> const& A_old_3_arr=(A_old[3])->array(mfi);
                ca_sdc_update_advection_o4_radau(bx, dt_m, dt, k_new_m_start_arr,
                                                 k_new_m_end_arr,
                                                 A_new_arr, A_old_0_arr, A_old_1_arr, A_old_2_arr,
                                                 A_old_3_arr, m_start);

            }

        }
#endif

    }
}


#ifdef REACTIONS
void
Castro::construct_old_react_source(MultiFab& U_state,
                                   MultiFab& R_source,
                                   const bool input_is_average)
{

    BL_PROFILE("Castro::construct_old_react_source()");

    auto domain_lo = geom.Domain().loVect3d();
    auto domain_hi = geom.Domain().hiVect3d();

    if (sdc_order == 4 && input_is_average)
    {
        // we have cell-averages
        // Note: we cannot tile these operations

        FArrayBox U_center;
        FArrayBox R_center;
        FArrayBox tmp;

        for (MFIter mfi(U_state); mfi.isValid(); ++mfi)
        {

            const Box& bx = mfi.tilebox();
            const Box& obx = mfi.growntilebox(1);

            // Convert to centers
            U_center.resize(obx, NUM_STATE);
            Elixir elix_u_center = U_center.elixir();
            auto const U_center_arr = U_center.array();

            make_cell_center(obx, U_state.array(mfi), U_center_arr, domain_lo, domain_hi);

            // burn, including one ghost cell
            R_center.resize(obx, NUM_STATE);
            Elixir elix_r_center = R_center.elixir();
            auto const R_center_arr = R_center.array();

            // ca_instantaneous_react(BL_TO_FORTRAN_BOX(obx),
            //                        BL_TO_FORTRAN_3D(U_center),
            //                        BL_TO_FORTRAN_3D(R_center));
            amrex::ParallelFor(obx,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
            {
                instantaneous_react(i, j, k, U_center_arr, R_center_arr);
            });

            // at this point, we have the reaction term on centers,
            // including a ghost cell.  Save this into Sburn so we can use
            // it later for the plotfile filling
            Sburn[mfi].copy(R_center, obx, 0, obx, 0, NUM_STATE);

            // convert R to averages (in place)

            tmp.resize(bx, 1);
            Elixir elix_tmp = tmp.elixir();
            auto const tmp_arr = tmp.array();

            make_fourth_in_place(bx, R_center_arr, tmp_arr, domain_lo, domain_hi);

            // copy this to the center
            R_source[mfi].copy(R_center, bx, 0, bx, 0, NUM_STATE);
        }

    }
    else
    {
        // we are cell-centers

        for (MFIter mfi(U_state); mfi.isValid(); ++mfi)
        {

            const Box& bx = mfi.tilebox();

            auto const U_state_arr = U_state.array(mfi);
            auto const R_source_arr = R_source.array(mfi);

            // construct the reactive source term
            // ca_instantaneous_react(BL_TO_FORTRAN_BOX(bx),
            //                        BL_TO_FORTRAN_3D(U_state[mfi]),
            //                        BL_TO_FORTRAN_3D(R_source[mfi]));
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
            {
                instantaneous_react(i, j, k, U_state_arr, R_source_arr);
            }); 
        }
    }
}
#endif
