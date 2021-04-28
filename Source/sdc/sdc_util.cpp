#include <Castro.H>
#include <Castro_F.H>
#include <Castro_sdc_util.H>

using namespace amrex;

void
Castro::ca_sdc_update_advection_o2_lobatto(const Box& bx,
                                           Real dt_m, Real dt,
                                           Array4<const Real> const& k_m,
                                           Array4<Real> const& k_n,
                                           Array4<const Real> const& A_m,
                                           Array4<const Real> const& A_0_old,
                                           Array4<const Real> const& A_1_old,
                                           int m_start)
{
    // update k_m to k_n via advection -- this is a second-order accurate update
    // for the Gauss-Lobatto discretization of the time nodes

    // here, dt_m is the update for this stage, from one time node to the next
    // dt is the update over the whole timestep, n to n+1

    // Gauss-Lobatto / trapezoid

    AMREX_PARALLEL_FOR_4D(bx, k_n.nComp(), i, j, k, n,
    {
        k_n(i,j,k,n) = k_m(i,j,k,n) + 0.5_rt * dt * (A_0_old(i,j,k,n) + A_1_old(i,j,k,n));
    });
}


void
Castro::ca_sdc_update_advection_o2_radau(const Box& bx,
                                         Real dt_m, Real dt,
                                         Array4<const Real> const& k_m,
                                         Array4<Real> const& k_n,
                                         Array4<const Real> const& A_m,
                                         Array4<const Real> const& A_0_old,
                                         Array4<const Real> const& A_1_old,
                                         Array4<const Real> const& A_2_old,
                                         int m_start)
{
    // update k_m to k_n via advection -- this is a second-order accurate update
    // for the Radau discretization of the time nodes

    // here, dt_m is the update for this stage, from one time node to the next
    // dt is the update over the whole timestep, n to n+1

    // Radau

    if (m_start == 0)
    {
        AMREX_PARALLEL_FOR_4D(bx, k_n.nComp(), i, j, k, n,
        {
            k_n(i,j,k,n) = k_m(i,j,k,n) +
                dt_m * (A_m(i,j,k,n) - A_0_old(i,j,k,n)) +
                dt/12.0_rt * (5.0_rt*A_1_old(i,j,k,n) - A_2_old(i,j,k,n));
        });
    }
    else if (m_start == 1)
    {
        AMREX_PARALLEL_FOR_4D(bx, k_n.nComp(), i, j, k, n,
        {
            k_n(i,j,k,n) = k_m(i,j,k,n) +
                dt_m * (A_m(i,j,k,n) - A_1_old(i,j,k,n)) +
                dt/3.0_rt * (A_1_old(i,j,k,n) + A_2_old(i,j,k,n));
        });
    }
}

void
Castro::ca_sdc_update_advection_o4_lobatto(const Box& bx,
                                           Real dt_m, Real dt,
                                           Array4<const Real> const& k_m,
                                           Array4<Real> const& k_n,
                                           Array4<const Real> const& A_m,
                                           Array4<const Real> const& A_0_old,
                                           Array4<const Real> const& A_1_old,
                                           Array4<const Real> const& A_2_old,
                                           int m_start)
{
    // update k_m to k_n via advection -- this is a fourth order accurate update

    // here, dt_m is the update for this stage, from one time node to the next
    // dt is the update over the whole timestep, n to n+1

    // Gauss-Lobatto (Simpsons)

    if (m_start == 0)
    {
        AMREX_PARALLEL_FOR_4D(bx, k_n.nComp(), i, j, k, n,
        {
            k_n(i,j,k,n) = k_m(i,j,k,n) +
                dt_m * (A_m(i,j,k,n) - A_0_old(i,j,k,n)) +
                dt/24.0_rt * (5.0_rt*A_0_old(i,j,k,n) + 8.0_rt*A_1_old(i,j,k,n) - A_2_old(i,j,k,n));
        });
    }
    else if (m_start == 1)
    {
        AMREX_PARALLEL_FOR_4D(bx, k_n.nComp(), i, j, k, n,
        {
            k_n(i,j,k,n) = k_m(i,j,k,n) +
                dt_m * (A_m(i,j,k,n) - A_1_old(i,j,k,n)) +
                dt/24.0_rt * (-A_0_old(i,j,k,n) + 8.0_rt*A_1_old(i,j,k,n) + 5.0_rt*A_2_old(i,j,k,n));
        });
    }
    else
    {
        amrex::Abort("error in ca_sdc_update_advection_o4_lobatto -- should not be here");
    }
}

void
Castro::ca_sdc_update_advection_o4_radau(const Box& bx,
                                         Real dt_m, Real dt,
                                         Array4<const Real> const& k_m,
                                         Array4<Real> const& k_n,
                                         Array4<const Real> const& A_m,
                                         Array4<const Real> const& A_0_old,
                                         Array4<const Real> const& A_1_old,
                                         Array4<const Real> const& A_2_old,
                                         Array4<const Real> const& A_3_old,
                                         int m_start)
{
    // update k_m to k_n via advection -- this is a fourth-order accurate update

    // here, dt_m is the update for this stage, from one time node to the next
    // dt is the update over the whole timestep, n to n+1

    // Gauss-Lobatto (Simpsons)

    if (m_start == 0)
    {
        AMREX_PARALLEL_FOR_4D(bx, k_n.nComp(), i, j, k, n,
        {
            k_n(i,j,k,n) = k_m(i,j,k,n) +
                dt_m * (A_m(i,j,k,n) - A_0_old(i,j,k,n)) +
                dt/1800.0_rt * ((-35.0_rt*std::sqrt(6.0_rt) + 440.0_rt)*A_1_old(i,j,k,n) +
                                (-169.0_rt*std::sqrt(6.0_rt) + 296.0_rt)*A_2_old(i,j,k,n) +
                                (-16.0_rt + 24.0_rt*std::sqrt(6.0_rt))*A_3_old(i,j,k,n));
        });
    }
    else if (m_start == 1)
    {
        AMREX_PARALLEL_FOR_4D(bx, k_n.nComp(), i, j, k, n,
        {
            k_n(i,j,k,n) = k_m(i,j,k,n) +
                dt_m * (A_m(i,j,k,n) - A_1_old(i,j,k,n)) +
                dt/150.0_rt * ((-12.0_rt + 17.0_rt*std::sqrt(6.0_rt))*A_1_old(i,j,k,n) +
                               (12.0_rt + 17.0_rt*std::sqrt(6.0_rt))*A_2_old(i,j,k,n) +
                               (-4.0_rt*std::sqrt(6.0_rt))*A_3_old(i,j,k,n));
        });
    }
    else if (m_start == 2)
    {
        AMREX_PARALLEL_FOR_4D(bx, k_n.nComp(), i, j, k, n,
        {
            k_n(i,j,k,n) = k_m(i,j,k,n) +
                dt_m * (A_m(i,j,k,n) - A_2_old(i,j,k,n)) +
                dt/600.0_rt * ((168.0_rt - 73.0_rt*std::sqrt(6.0_rt))*A_1_old(i,j,k,n) +
                               (120.0_rt + 5.0_rt*std::sqrt(6.0_rt))*A_2_old(i,j,k,n) +
                               (72.0_rt + 8.0_rt*std::sqrt(6.0_rt))*A_3_old(i,j,k,n));
        });
    }
    else
    {
        amrex::Abort("error in ca_sdc_update_advection_o4_radau -- should not be here");
    }
}

#ifdef REACTIONS
void
Castro::ca_sdc_compute_C2_lobatto(const Box& bx,
                                  Real dt_m, Real dt,
                                  Array4<const Real> const& A_m,
                                  Array4<const Real> const& A_0_old,
                                  Array4<const Real> const& A_1_old,
                                  Array4<const Real> const& R_0_old,
                                  Array4<const Real> const& R_1_old,
                                  Array4<Real> const& C,
                                  int m_start)
{
    // compute the source term C for the 2nd order Lobatto update

    // Here, dt_m is the timestep between time-nodes m and m+1

    AMREX_PARALLEL_FOR_4D(bx, C.nComp(), i, j, k, n,
    {
        // construct the source term to the update for 2nd order
        // Lobatto, there is no advective correction, and we have
        // C = - R(U^{m+1,k}) + I_m^{m+1}/dt
        C(i,j,k,n) = -R_1_old(i,j,k,n) +
            0.5_rt * (A_0_old(i,j,k,n) + A_1_old(i,j,k,n)) +
            0.5_rt * (R_0_old(i,j,k,n) + R_1_old(i,j,k,n));
    });

} // end ca_sdc_compute_C2_lobatto

void
Castro::ca_sdc_compute_C2_radau(const Box& bx,
                                Real dt_m, Real dt,
                                Array4<const Real> const& A_m,
                                Array4<const Real> const& A_0_old,
                                Array4<const Real> const& A_1_old,
                                Array4<const Real> const& A_2_old,
                                Array4<const Real> const& R_0_old,
                                Array4<const Real> const& R_1_old,
                                Array4<const Real> const& R_2_old,
                                Array4<Real> const& C,
                                int m_start)
{
    // compute the source term C for the 2nd order Radau update

    // Here, dt_m is the timestep between time-nodes m and m+1

    // construct the source term to the update for 2nd order
    // Radau

    if (m_start == 0)
    {
        AMREX_PARALLEL_FOR_4D(bx, C.nComp(), i, j, k, n,
        {
            C(i,j,k,n) = -R_1_old(i,j,k,n) +
                (A_m(i,j,k,n) - A_0_old(i,j,k,n)) +
                (dt/dt_m) * (1.0_rt/12.0_rt) *
                5.0_rt*(A_1_old(i,j,k,n) + R_1_old(i,j,k,n)) -
                (A_2_old(i,j,k,n) + R_2_old(i,j,k,n));
        });
    }
    else if (m_start == 1)
    {
        AMREX_PARALLEL_FOR_4D(bx, C.nComp(), i, j, k, n,
        {
            C(i,j,k,n) = -R_2_old(i,j,k,n) +
                (A_m(i,j,k,n) - A_1_old(i,j,k,n)) +
                (dt/dt_m) * (1.0_rt/3.0_rt) *
                (A_1_old(i,j,k,n) + R_1_old(i,j,k,n)) +
                (A_2_old(i,j,k,n) + R_2_old(i,j,k,n));
        });
    }
} // end ca_sdc_compute_C2_radau

void
Castro::ca_sdc_compute_C4_lobatto(const Box& bx,
                                  Real dt_m, Real dt,
                                  Array4<const Real> const& A_m,
                                  Array4<const Real> const& A_0_old,
                                  Array4<const Real> const& A_1_old,
                                  Array4<const Real> const& A_2_old,
                                  Array4<const Real> const& R_0_old,
                                  Array4<const Real> const& R_1_old,
                                  Array4<const Real> const& R_2_old,
                                  Array4<Real> const& C,
                                  int m_start)
{
    // compute the 'C' term for the 4th-order solve with reactions
    // note: this 'C' is cell-averages


    // Gauss-Lobatto (Simpsons)

    if (m_start == 0)
    {
        // compute the integral from [t_m, t_{m+1}], normalized by dt_m
        AMREX_PARALLEL_FOR_4D(bx, C.nComp(), i, j, k, n,
        {
            Real integral = 1.0_rt/12.0_rt * (5.0_rt*(A_0_old(i,j,k,n) + R_0_old(i,j,k,n)) +
                                              8.0_rt*(A_1_old(i,j,k,n) + R_1_old(i,j,k,n)) -
                                              (A_2_old(i,j,k,n) + R_2_old(i,j,k,n)));

            C(i,j,k,n) = (A_m(i,j,k,n) - A_0_old(i,j,k,n)) - R_1_old(i,j,k,n) + integral;
        });
    }
    else if (m_start == 1)
    {
        // compute the integral from [t_m, t_{m+1}], normalized by dt_m
        AMREX_PARALLEL_FOR_4D(bx, C.nComp(), i, j, k, n,
        {
            Real integral = 1.0_rt/12.0_rt * (-(A_0_old(i,j,k,n) + R_0_old(i,j,k,n)) +
                                              8.0_rt*(A_1_old(i,j,k,n) + R_1_old(i,j,k,n)) +
                                              5.0_rt*(A_2_old(i,j,k,n) + R_2_old(i,j,k,n)));

            C(i,j,k,n) = (A_m(i,j,k,n) - A_1_old(i,j,k,n)) - R_2_old(i,j,k,n) + integral;
        });
    }
    else
    {
        amrex::Abort("error in ca_sdc_compute_C4 -- should not be here");
    }

} // end subroutine ca_sdc_compute_C4_lobatto


void
Castro::ca_sdc_compute_C4_radau(const Box& bx,
                                Real dt_m, Real dt,
                                Array4<const Real> const& A_m,
                                Array4<const Real> const& A_0_old,
                                Array4<const Real> const& A_1_old,
                                Array4<const Real> const& A_2_old,
                                Array4<const Real> const& A_3_old,
                                Array4<const Real> const& R_0_old,
                                Array4<const Real> const& R_1_old,
                                Array4<const Real> const& R_2_old,
                                Array4<const Real> const& R_3_old,
                                Array4<Real> const& C,
                                int m_start)
{
    // compute the 'C' term for the 4th-order solve with reactions
    // note: this 'C' is cell-averages

    // Gauss-Lobatto (Simpsons)

    if (m_start == 0)
    {
        // compute the integral from [t_m, t_{m+1}], normalized by dt_m
        AMREX_PARALLEL_FOR_4D(bx, C.nComp(), i, j, k, n,
        {
            Real integral = (dt/dt_m) * (1.0_rt/1800.0_rt) *
                ((-35.0_rt*std::sqrt(6.0_rt) + 440.0_rt)*(A_1_old(i,j,k,n) + R_1_old(i,j,k,n)) +
                 (-169.0_rt*std::sqrt(6.0_rt) + 296.0_rt)*(A_2_old(i,j,k,n) + R_2_old(i,j,k,n)) +
                 (-16.0_rt + 24.0_rt*std::sqrt(6.0_rt))*(A_3_old(i,j,k,n) + R_3_old(i,j,k,n)));

            C(i,j,k,n) = (A_m(i,j,k,n) - A_0_old(i,j,k,n)) - R_1_old(i,j,k,n) + integral;
        });
    }
    else if (m_start == 1)
    {
        // compute the integral from [t_m, t_{m+1}], normalized by dt_m
        AMREX_PARALLEL_FOR_4D(bx, C.nComp(), i, j, k, n,
        {
            Real integral = (dt/dt_m) * (1.0_rt/150.0_rt) *
                ((-12.0_rt + 17.0_rt*std::sqrt(6.0_rt))*(A_1_old(i,j,k,n) + R_1_old(i,j,k,n)) +
                 (12.0_rt + 17.0_rt*std::sqrt(6.0_rt))*(A_2_old(i,j,k,n) + R_2_old(i,j,k,n)) +
                 (-4.0_rt*std::sqrt(6.0_rt))*(A_3_old(i,j,k,n) + R_3_old(i,j,k,n)));

            C(i,j,k,n) = (A_m(i,j,k,n) - A_1_old(i,j,k,n)) - R_2_old(i,j,k,n) + integral;
        });
    }
    else if (m_start == 2)
    {
        // compute the integral from [t_m, t_{m+1}], normalized by dt_m
        AMREX_PARALLEL_FOR_4D(bx, C.nComp(), i, j, k, n,
        {
            Real integral = (dt/dt_m) * (1.0_rt/600.0_rt) *
                ((168.0_rt - 73.0_rt*std::sqrt(6.0_rt))*(A_1_old(i,j,k,n) + R_1_old(i,j,k,n)) +
                 (120.0_rt + 5.0_rt*std::sqrt(6.0_rt))*(A_2_old(i,j,k,n) + R_2_old(i,j,k,n)) +
                 (72.0_rt + 8.0_rt*std::sqrt(6.0_rt))*(A_3_old(i,j,k,n) + R_3_old(i,j,k,n)));

            C(i,j,k,n) = (A_m(i,j,k,n) - A_2_old(i,j,k,n)) - R_3_old(i,j,k,n) + integral;
        });
    }
    else
    {
        amrex::Abort("error in ca_sdc_compute_C4 -- should not be here");
    }
} //  end ca_sdc_compute_C4_radau

void
Castro::ca_sdc_conservative_update(const Box& bx, Real const dt_m,
                                   Array4<const Real> const& U_old,
                                   Array4<Real> const& U_new,
                                   Array4<const Real> const& C,
                                   Array4<const Real> const& R_new)
{
    // given <U>_old, <R>_new, and <C>, compute <U>_new

    // now consider the reacting system
    AMREX_PARALLEL_FOR_4D(bx, U_new.nComp(), i, j, k, n,
    {
        U_new(i,j,k,n) = U_old(i,j,k,n) + dt_m * R_new(i,j,k,n) + dt_m * C(i,j,k,n);
    });

    // enforce normalization of the species

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
    {
        normalize_species_sdc(i, j, k, U_new);
    });


} // end subroutine ca_sdc_conservative_update
#endif


void Castro::ca_sdc_compute_initial_guess(const Box& bx,
                                          Array4<const Real> const& U_old,
                                          Array4<const Real> const& U_new,
                                          Array4<const Real> const& A_old,
                                          Array4<const Real> const& R_old,
                                          Array4<Real> const& U_guess,
                                          Real const dt_m, int const sdc_iteration)
{
    // compute the initial guess for the Newton solve
    // Here dt_m is the timestep to update from time node m to m+1

    if (sdc_iteration == 0)
    {
        AMREX_PARALLEL_FOR_4D(bx, U_guess.nComp(), i, j, k, n,
        {
            U_guess(i,j,k,n) = U_old(i,j,k,n) + dt_m * A_old(i,j,k,n) + dt_m * R_old(i,j,k,n);
        });
    }
    else
    {
        AMREX_PARALLEL_FOR_4D(bx, U_guess.nComp(), i, j, k, n,
        {
            U_guess(i,j,k,n) = U_new(i,j,k,n);
        });
    }

}


#ifdef REACTIONS
void Castro::ca_store_reaction_state(const Box& bx,
                                     Array4<const Real> const& R_old,
                                     Array4<const Real> const& state,
                                     Array4<Real> const& R_store)
{
    // copy the data from the last node's reactive source to the state data

    // for R_store we use the indices defined in Castro_setup.cpp for
    // Reactions_Type

    AMREX_PARALLEL_FOR_4D(bx, NumSpec, i, j, k, n,
    {
        R_store(i,j,k,n) = R_old(i,j,k,UFS+n);
    });

#if NAUX_NET > 0
    AMREX_PARALLEL_FOR_4D(bx, NumAux, i, j, k, n,
    {
        R_store(i,j,k,n) = R_old(i,j,k,UFX+n);
    }
#endif

    AMREX_PARALLEL_FOR_3D(bx, i, j, k,
    {
        R_store(i,j,k,NumSpec+NumAux  ) = R_old(i,j,k,UEDEN);
        R_store(i,j,k,NumSpec+NumAux+1) = 0.0; // we're not storing the weights
    });
}

#endif
