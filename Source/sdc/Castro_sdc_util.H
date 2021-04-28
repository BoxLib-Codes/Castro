#ifndef CASTRO_SDC_UTIL_H
#define CASTRO_SDC_UTIL_H

#include <Castro.H>
#ifdef REACTIONS
#include <Castro_react_util.H>
#include <vode_dvode.H>
#include <linpack.H>
#endif

// error codes
constexpr int NEWTON_SUCCESS = 0;
constexpr int SINGULAR_MATRIX = -1;
constexpr int CONVERGENCE_FAILURE = -2;

// solvers
constexpr int NEWTON_SOLVE = 1;
constexpr int VODE_SOLVE = 2;
constexpr int HYBRID_SOLVE = 3;

constexpr int ldjac = NumSpec + 2;


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
normalize_species_sdc(const int i, const int j, const int k,
                      Array4<Real> const& u) {

    Real sum_rhoX = 0.0_rt;

    // we are working on partial densities, so we want to first limit
    // the species to lie between [0, rho]

    sum_rhoX = 0.0_rt;

    for (int n = 0; n < NumSpec; ++n) {
        u(i,j,k,UFS+n) = amrex::max(small_x * u(i,j,k,URHO),
                                    amrex::min(u(i,j,k,URHO), u(i,j,k,UFS+n)));
        sum_rhoX += u(i,j,k,UFS+n);
    }

    for (int n = 0; n < NumSpec; ++n) {
        u(i,j,k,UFS+n) *= u(i,j,k,URHO) / sum_rhoX;
    }
}

#ifdef REACTIONS


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
f_sdc_jac(const Real dt_m,
          GpuArray<Real, NumSpec+2> const& U,
          GpuArray<Real, NumSpec+2>& f,
          RArray2D& Jac,
          GpuArray<Real, NumSpec+2>& f_source,
          GpuArray<Real, 3>& mom_source,
          const Real T_old,
          const Real E_var) {

    // This is used with the Newton solve and returns f and the Jacobian

    GpuArray<Real, NUM_STATE> U_full;
    GpuArray<Real, NUM_STATE> R_full;
    GpuArray<Real, NumSpec+2> R_react;

    Array2D<Real, 0, NumSpec+1, 0, NumSpec+1> dRdw = {0.0_rt};
    Array2D<Real, 0, NumSpec+1, 0, NumSpec+1> dwdU = {0.0_rt};

    Array2D<Real, 0, NumSpec+1, 0, NumSpec+1> dRdU;

    // we are not solving the momentum equations
    // create a full state -- we need this for some interfaces
    U_full[URHO] = U[0];

    for (int n = 0; n < NumSpec; ++n) {
        U_full[UFS+n] = U[1+n];
    }

    if (sdc_solve_for_rhoe == 1) {
        U_full[UEINT] = U[NumSpec+1];
        U_full[UEDEN] = E_var;
    } else {
        U_full[UEDEN] = U[NumSpec+1];
        U_full[UEINT] = E_var;
    }

    for (int n = 0; n < 3; ++n) {
        U_full[UMX+n] = mom_source[n];
    }

    // normalize the species
    auto sum_rhoX = 0.0_rt;
    for (int n = 0; n < NumSpec; ++n) {
        U_full[UFS+n] = amrex::max(small_x, U_full[UFS+n]);
        sum_rhoX += U_full[UFS+n];
    }

    for (int n = 0; n < NumSpec; ++n) {
        U_full[UFS+n] *= U_full[URHO] / sum_rhoX;
    }

    // compute the temperature and species derivatives --
    // maybe this should be done using the burn_state
    // returned by single_zone_react_source, since it is
    // more consistent T from e

    eos_t eos_state;
    eos_state.rho = U_full[URHO];
    eos_state.T = T_old;   // initial guess
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = U_full[UFS+n] / U_full[URHO];
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; ++n) {
        eos_state.aux[n] = U_full[UFX+n] / U_full[URHO];
    }
#endif
    eos_state.e = U_full[UEINT] / U_full[URHO];  // (U_full[UEDEN] - 0.5_rt*sum(U_full(UMX:UMZ))/U_full[URHO])/U_full[URHO]

    eos(eos_input_re, eos_state);

    U_full[UTEMP] = eos_state.T;

    // we'll create a burn_state to pass stuff from the RHS to the Jac function
    burn_t burn_state;

    single_zone_react_source(U_full, R_full, burn_state);

    // store the subset of R used in the Jacobian
    R_react[0] = R_full[URHO];
    for (int n = 0; n < NumSpec; ++n) {
        R_react[n+1] = R_full[UFS+n];
    }
    if (sdc_solve_for_rhoe == 1) {
        R_react[NumSpec+1] = R_full[UEINT];
    } else {
        R_react[NumSpec+1] = R_full[UEDEN];
    }

    for (int n = 0; n < NumSpec+2; ++n) {
        f[n] = U[n] - dt_m * R_react[n] - f_source[n];
    }

    // get dRdw -- this may do a numerical approxiation or use the
    // network's analytic Jac
    single_zone_jac(U_full, burn_state, dRdw);

    // the density row
    dwdU(iwrho, 0) = 1.0_rt;

    // the X_k rows
    for (int m = 1; m < NumSpec+1; ++m) {
        dwdU(iwfs-1+m, 0) = -U[m] / (U[0]*U[0]);
        dwdU(iwfs-1+m, m) = 1.0_rt / U[0];
    }

    auto eos_xderivs = composition_derivatives(eos_state);

    // now the e row -- this depends on whether we are evolving (rho E) or (rho e)
    auto denom = 1.0_rt / eos_state.rho;
    auto xn_sum = 0.0_rt;
    for (int n = 0; n < NumSpec; ++n) {
        xn_sum += eos_state.xn[n] * eos_xderivs.dedX[n];
    }

    if (sdc_solve_for_rhoe == 1) {
        dwdU(iwe, 0) = denom * (xn_sum - eos_state.e);
    } else {
        auto u2_sum = 0.0_rt;
        for (auto n = 0; n < 3; ++n) {
            u2_sum += U_full[UMX+n] * U_full[UMX+n];
        }
        dwdU(iwe, 0) = denom * (xn_sum - eos_state.e -
                               0.5_rt * u2_sum / (eos_state.rho * eos_state.rho));
    }

    for (int m = 0; m < NumSpec; ++m) {
        dwdU(iwe,m+1) = -denom * eos_xderivs.dedX[m];
    }

    dwdU(iwe, NumSpec+1) = denom;

    // construct the Jacobian -- we can get most of the
    // terms from the network itself, but we do not rely on
    // it having derivative wrt density

    // Note: Jac is 1-based, because that is what the linear algebra
    // routines expect but the components, dRdw and dwdU are 0-based!

    for (int n = 1; n <= NumSpec+2; ++n) {
        for (int m = 1; m <= NumSpec+2; ++m) {
            Jac(n, m) = 0.0_rt;
        }
    }

    for (int m = 1; m <= NumSpec+2; ++m) {
        Jac(m, m) = 1.0_rt;
    }

    // auto dRdU = matmul(dRdw, dwdU);
    // dRdU is 0-based

    for (int n = 0; n <= NumSpec+1; ++n) {
        for (int m = 0; m <= NumSpec+1; ++m) {
            dRdU(n,m) = 0.0_rt;
            for (int l = 0; l <= NumSpec+1; ++l) {
                dRdU(n,m) += dRdw(n,l) * dwdU(l,m);
            }
        }
    }

    for (int n = 0; n < NumSpec+2; ++n) {
        for (int m = 0; m < NumSpec+2; ++m) {
            Jac(n+1, m+1) -= dt_m * dRdU(n,m);
        }
    }

}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
sdc_newton_solve(const Real dt_m,
                 GpuArray<Real, NUM_STATE> const& U_old,
                 GpuArray<Real, NUM_STATE> & U_new,
                 GpuArray<Real, NUM_STATE> const& C,
                 const int sdc_iteration,
                 Real& err_out,
                 int& ierr) {
    // the purpose of this function is to solve the system
    // U - dt R(U) = U_old + dt C using a Newton solve.
    //
    // here, U_new should come in as a guess for the new U
    // and will be returned with the value that satisfied the
    // nonlinear function

    RArray2D Jac;

    // we will do the implicit update of only the terms that
    // have reactive sources
    //
    // 0 : rho
    // 1:NumSpec : species
    // NumSpec+1 : (rho E) or (rho e)

    GpuArray<Real, NumSpec+2> U_react;
    GpuArray<Real, NumSpec+2> f_source;
    GpuArray<Real, 3> mom_source;
    GpuArray<Real, NumSpec+2> dU_react;
    GpuArray<Real, NumSpec+2> f;
    RArray1D f_rhs;

    const int MAX_ITER = 100;

    ierr = NEWTON_SUCCESS;

    // the tolerance we are solving to may depend on the
    // iteration
    Real relax_fac = std::pow(sdc_solver_relax_factor, sdc_order - sdc_iteration - 1);
    Real tol_dens = sdc_solver_tol_dens * relax_fac;
    Real tol_spec = sdc_solver_tol_spec * relax_fac;
    Real tol_ener = sdc_solver_tol_ener * relax_fac;

    // update the momenta for this zone -- they don't react
    for (int n = 0; n < 3; ++n) {
        U_new[UMX+n] = U_old[UMX+n] + dt_m * C[UMX+n];
    }

    // now only save the subset that participates in the
    // nonlinear solve -- note: we include the old state in
    // f_source

    // load rpar

    // for the Jacobian solve, we are solving
    //   f(U) = U - dt R(U) - U_old - dt C = 0
    // we define f_source = U_old + dt C so we are solving
    //   f(U) = U - dt R(U) - f_source = 0

    f_source[0] = U_old[URHO] + dt_m * C[URHO];
    for (int n = 0; n < NumSpec; ++n) {
        f_source[1 + n] = U_old[UFS + n] + dt_m * C[UFS + n];
    }
    if (sdc_solve_for_rhoe == 1) {
        f_source[NumSpec+1] = U_old[UEINT] + dt_m * C[UEINT];
    } else {
        f_source[NumSpec+1] = U_old[UEDEN] + dt_m * C[UEDEN];
    }

    // set the momenta to be U_new
    for (int n = 0; n < 3; ++n) {
        mom_source[n] = U_new[UMX+n];
    }

    // temperature will be used as an initial guess in the EOS

    Real T_old = U_old[UTEMP];

    // we should be able to do an update for this somehow?

    Real E_var;

    if (sdc_solve_for_rhoe == 1) {
        E_var = U_new[UEDEN];
    } else {
        E_var = U_new[UEINT];
    }

    // store the subset for the nonlinear solve
    // We use an initial guess if possible
    U_react[0] = U_new[URHO];
    for (int n = 0; n < NumSpec; ++n) {
        U_react[1+n] = U_new[UFS+n];
    }
    if (sdc_solve_for_rhoe == 1) {
        U_react[NumSpec+1] = U_new[UEINT];
    } else {
        U_react[NumSpec+1] = U_new[UEDEN];
    }

#if (INTEGRATOR == 0)

    // do a simple Newton solve

    // iterative loop
    int iter = 0;
    int max_newton_iter = MAX_ITER;

    Real err = 1.e30_rt;
    bool converged = false;

    while (!converged && iter < max_newton_iter) {
        int info = 0;
        f_sdc_jac(dt_m, U_react, f, Jac, f_source, mom_source, T_old, E_var);

        IArray1D ipvt;

        // solve the linear system: Jac dU_react = -f
        dgefa<NumSpec+2>(Jac, ipvt, info);
        if (info != 0) {
            ierr = SINGULAR_MATRIX;
            return;
        }

        for (int n = 1; n <= NumSpec+2; ++n) {
            f_rhs(n) = -f[n-1];
        }

        dgesl<NumSpec+2>(Jac, ipvt, f_rhs);

        for (int n = 0; n < NumSpec+2; ++n) {
            dU_react[n] = f_rhs(n+1);
        }

        // how much of dU_react should we apply?
        Real eta = 1.0_rt;
        for (int n = 0; n < NumSpec+2; ++n) {
            dU_react[n] *= eta;
            U_react[n] += dU_react[n];
        }

        GpuArray<Real, NumSpec+2> eps_tot;

        eps_tot[0] = tol_dens * std::abs(U_react[0]) + sdc_solver_atol;

        // for species, atol is the mass fraction limit, so we
        // multiply by density to get a partial density limit
        for (int n = 0; n < NumSpec; ++n) {
            eps_tot[1 + n] = tol_spec * std::abs(U_react[1 + n]) + sdc_solver_atol * std::abs(U_react[0]);
        }
        eps_tot[NumSpec+1] = tol_ener * std::abs(U_react[NumSpec+1]) + sdc_solver_atol;

        // compute the norm of the weighted error, where the
        // weights are 1/eps_tot
        auto err_sum = 0.0_rt;
        for (int n = 0; n < NumSpec+2; ++n) {
            err_sum += dU_react[n] * dU_react[n] / (eps_tot[n]* eps_tot[n]);
        }
        err = std::sqrt(err_sum / (NumSpec+2));

        if (err < 1.0_rt) {
            converged = true;
        }
        iter++;
    }

    err_out = err;

    if (!converged) {
        ierr = CONVERGENCE_FAILURE;
        return;
    }

#endif

    // update the full U_new
    // if we updated total energy, then correct internal,
    // or vice versa
    U_new[URHO] = U_react[0];
    for (int n = 0; n < NumSpec; ++n) {
        U_new[UFS+n] = U_react[1+n];
    }
    auto v2 = 0.0_rt;
    for (int m = 0; m < 3; ++m) {
        v2 += U_new[UMX+m] * U_new[UMX+m];
    }

    if (sdc_solve_for_rhoe == 1) {
        U_new[UEINT] = U_react[NumSpec+1];
        U_new[UEDEN] = U_new[UEINT] + 0.5_rt * v2 / U_new[URHO];
    } else {
        U_new[UEDEN] = U_react[NumSpec+1];
        U_new[UEINT] = U_new[UEDEN] - 0.5_rt * v2 / U_new[URHO];
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
sdc_newton_subdivide(const Real dt_m,
                     GpuArray<Real, NUM_STATE> const& U_old,
                     GpuArray<Real, NUM_STATE>& U_new,
                     GpuArray<Real, NUM_STATE> const& C,
                     const int sdc_iteration,
                     Real& err_out,
                     int& ierr) {
    // This is the driver for solving the nonlinear update for
    // the reating/advecting system using Newton's method. It
    // attempts to do the solution for the full dt_m requested,
    // but if it fails, will subdivide the domain until it
    // converges or reaches our limit on the number of
    // subintervals.

    const int MAX_NSUB = 64;
    GpuArray<Real, NUM_STATE> U_begin;

    // subdivide the timestep and do multiple Newtons. We come
    // in here with an initial guess for the new solution
    // stored in U_new. That only really makes sense for the
    // case where we have 1 substep. Otherwise, we should just
    // use the old time solution.

    int nsub = 1;
    ierr = CONVERGENCE_FAILURE;

    for (int n = 0; n < NUM_STATE; ++n) {
        U_begin[n] = U_old[n];
    }

    while (nsub < MAX_NSUB && ierr != NEWTON_SUCCESS) {
        if (nsub > 1) {
            for (int n = 0; n < NUM_STATE; ++n) {
                U_new[n] = U_old[n];
            }
        }
        Real dt_sub = dt_m / nsub;

        for (int isub = 0; isub < nsub; ++isub) {
            // normalize species
            Real sum_rhoX = 0.0_rt;
            for (int n = 0; n < NumSpec; ++n) {
                U_begin[UFS + n] = amrex::max(small_x, U_begin[UFS + n]);
                sum_rhoX += U_begin[UFS + n];
            }
            for (int n = 0; n < NumSpec; ++n) {
                U_begin[UFS + n] *= U_begin[URHO] / sum_rhoX;
            }

            sdc_newton_solve(dt_sub, U_begin, U_new, C, sdc_iteration, err_out, ierr);

            for (int n = 0; n < NUM_STATE; ++n) {
                U_begin[n] = U_new[n];
            }
        }
        nsub *= 2;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
sdc_vode_solve(const Real dt_m,
               GpuArray<Real, NUM_STATE> const& U_old,
               GpuArray<Real, NUM_STATE>& U_new,
               GpuArray<Real, NUM_STATE> const& C,
               const int sdc_iteration) {

    // The purpose of this function is to solve the system the
    // approximate system dU/dt = R + C using the VODE ODE
    // integrator.
    // The solution we get here will then be used as the
    // initial guess to the Newton solve on the real system.

    // We will do the implicit update of only the terms that
    // have reactive sources
    //
    // 0 : rho
    // 1:NumSpec : species
    // NumSpec+1 : (rho E) or (rho e)

#if (INTEGRATOR == 0)

    // The tolerance we are solving to may depend on the
    // iteration
    auto relax_fac = std::pow(sdc_solver_relax_factor, sdc_order - sdc_iteration - 1);
    auto tol_dens = sdc_solver_tol_dens * relax_fac;
    auto tol_spec = sdc_solver_tol_spec * relax_fac;
    auto tol_ener = sdc_solver_tol_ener * relax_fac;

    // Update the momenta for this zone -- they don't react
    for (int n = 0; n < 3; ++n) {
        U_new[UMX+n] = U_old[UMX+n] + dt_m * C[UMX+n];
    }

    // Now only save the subset that participates in the
    // nonlinear solve -- note: we include the old state
    // in f_source

    // load rpar

    // If we are solving the system as an ODE, then we are
    // solving
    //   dU/dt = R(U) + C
    // so we simply pass in C
    GpuArray<Real, NumSpec+2> C_react;

    C_react[0] = C[URHO];
    for (int n = 0; n < NumSpec; ++n) {
        C_react[1+n] = C[UFS+n];
    }
    C_react[NumSpec+1] = C[UEINT];

    dvode_t dvode_state;

    burn_t burn_state;

    burn_state.success = true;

    for (int n = 0; n < NumSpec+2; n++) {
        burn_state.f_source[n] = C_react[n];
    }

    for (int n = 0; n < 3; n++) {
        burn_state.mom[n] = U_new[UMX+n];
    }

    // we are always solving for rhoe with the VODE predict
    burn_state.e = U_new[UEDEN];

    burn_state.E_var = U_old[UEDEN];

    // temperature will be used as an initial guess in the EOS
    burn_state.T = U_old[UTEMP];


    // store the subset for the nonlinear solve. We only
    // consider (rho e), not (rho E). This is because at
    // present we do not have a method of updating the
    // velocities during the multistep integration.

    // Also note that the dvode_state is 1-based, but we'll
    // access it as 0-based in our implementation of the
    // RHS routine

    dvode_state.y(1) = U_old[URHO];
    for (int n = 0; n < NumSpec; ++n) {
        dvode_state.y(2+n) = U_old[UFS+n];
    }
    dvode_state.y(NumSpec+2) = U_old[UEINT];

    // // set the maximum number of steps allowed
    // dvode_state.MXSTEP = 25000;

    dvode_state.t = 0.0_rt;
    dvode_state.tout = dt_m;
    dvode_state.HMXI = 1.0_rt / ode_max_dt;

    GpuArray<Real, NumSpec+2> rtol;
    GpuArray<Real, NumSpec+2> atol;

    // relative tolerances
    rtol[0] = tol_dens;
    for (int n = 0; n < NumSpec; ++n) {
        rtol[1+n] = tol_spec;
    }
    rtol[NumSpec+1] = tol_ener;

    // absolute tolerances
    atol[0] = sdc_solver_atol * U_old[URHO];
    for (int n = 0; n < NumSpec; ++n) {
        atol[1+n] = sdc_solver_atol * U_old[URHO]; // this way, atol is the minimum x
    }

    atol[NumSpec+1] = sdc_solver_atol * U_old[UEINT];

    for (int n = 1; n <= NumSpec+2; ++n) {
        dvode_state.atol(n) = atol[n-1];
        dvode_state.rtol(n) = rtol[n-1];
    }

    int istate = dvode(burn_state, dvode_state);

    if (istate < 0) {
        Abort("Vode terminated poorly");
    }

    // update the full U_new

    U_new[URHO] = dvode_state.y(1);
    Real sum_rhoX = 0.0_rt;
    for (int n = 0; n < NumSpec; ++n) {
        U_new[UFS+n] = amrex::max(0.0_rt, dvode_state.y(2+n));
        sum_rhoX += U_new[UFS+n];
    }

    // normalize the abundances to account for integration error

    for (int n = 0; n < NumSpec; ++n) {
        U_new[UFS+n] *= U_new[URHO] / sum_rhoX;
    }

    U_new[UEINT] = dvode_state.y(NumSpec+2);
    auto v2 = 0.0_rt;
    for (int n = 0; n < 3; ++n) {
        v2 += U_new[UMX+n] * U_new[UMX+n];
    }
    U_new[UEDEN] = U_new[UEINT] + 0.5_rt * v2 / U_new[URHO];

    // keep our temperature guess

    U_new[UTEMP] = U_old[UTEMP];

#endif
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
sdc_solve(const Real dt_m,
          GpuArray<Real, NUM_STATE> const& U_old,
          GpuArray<Real, NUM_STATE>& U_new,
          GpuArray<Real, NUM_STATE> const& C,
          const int sdc_iteration) {

    int ierr;
    Real err_out;

    // for debugging
    GpuArray<Real, NUM_STATE> U_orig;

    for (int n = 0; n < NUM_STATE; ++n) {
        U_orig[n] = U_old[n];
    }

    if (sdc_solver == NEWTON_SOLVE) {
        // We are going to assume we already have a good guess
        // for the solve in U_new and just pass the solve onto
        // the main Newton solve
        sdc_newton_subdivide(dt_m, U_old, U_new, C, sdc_iteration, err_out, ierr);

        // failing?
        if (ierr != NEWTON_SUCCESS) {
            Abort("Newton subcycling failed in sdc_solve");
        }
    } else if (sdc_solver == VODE_SOLVE) {
        // Use VODE to do the solution
        sdc_vode_solve(dt_m, U_old, U_new, C, sdc_iteration);

    } else if (sdc_solver == HYBRID_SOLVE) {
        // If it is the first iteration, we will use VODE to
        // predict the solution. Otherwise, we will use Newton.
        if (sdc_iteration == 0) {
            sdc_vode_solve(dt_m, U_old, U_new, C, sdc_iteration);
        }

        // Now U_new is the update that VODE predicts, so we
        // will use that as the initial guess to the Newton solve
        sdc_newton_subdivide(dt_m, U_old, U_new, C, sdc_iteration, err_out, ierr);

        // Failing?
        if (ierr != NEWTON_SUCCESS) {
            Abort("Newton failure in sdc_solve");
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
sdc_solve(const int i, const int j, const int k,
          const Real dt_m,
          Array4<const Real> const& U_old,
          Array4<Real> const& U_new,
          Array4<const Real> const& C,
          const int sdc_iteration) {
    // wrapper for the zone-by-zone version

    GpuArray<Real, NUM_STATE> U_old_zone;
    GpuArray<Real, NUM_STATE> U_new_zone;
    GpuArray<Real, NUM_STATE> C_zone;

    for (int n = 0; n < NUM_STATE; ++n) {
        U_old_zone[n] = U_old(i,j,k,n);
        U_new_zone[n] = U_new(i,j,k,n);
        C_zone[n] = C(i,j,k,n);
    }

    sdc_solve(dt_m, U_old_zone, U_new_zone, C_zone, sdc_iteration);

    for (int n = 0; n < NUM_STATE; ++n) {
        U_new(i,j,k,n) = U_new_zone[n];
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
sdc_update_o2(const int i, const int j, const int k,
              Array4<const Real> const& k_m,
              Array4<Real> const& k_n,
              Array4<const Real> const& A_m,
              Array4<const Real> const& R_m_old,
              Array4<const Real> const& C,
              const Real dt_m,
              const int sdc_iteration, const int m_start) {
    // update k_m to k_n via advection -- this is a second-order accurate update

    // Here, dt_m is the timestep between time-nodes m and m+1

    burn_t burn_state;

    GpuArray<Real, NUM_STATE> U_old;
    GpuArray<Real, NUM_STATE> U_new;
    GpuArray<Real, NUM_STATE> R_full;
    GpuArray<Real, NUM_STATE> C_zone;

    for (int n = 0; n < NUM_STATE; ++n) {
        U_old[n] = k_m(i,j,k,n);
        C_zone[n] = C(i,j,k,n);
    }

    // Only burn if we are within the temperature and density
    // limits for burning
    if (!okay_to_burn(U_old)) {
        for (int n = 0; n < NUM_STATE; ++n) {
            R_full[n] = 0.0_rt;
        }
    } else {

        // This is the full state -- this will be updated as we
        // solve the nonlinear system.  We want to start with a
        // good initial guess.  For later iterations, we should
        // begin with the result from the previous iteration.  For
        // the first iteration, let's try to extrapolate forward
        // in time.
        if (sdc_iteration == 0) {
            for (int n = 0; n < NUM_STATE; ++n) {
                U_new[n] = U_old[n] + dt_m * A_m(i,j,k,n) + dt_m * R_m_old(i,j,k,n);
            }
        } else {
            for (int n = 0; n < NUM_STATE; ++n) {
                U_new[n] = k_n(i,j,k,n);
            }
        }

        sdc_solve(dt_m, U_old, U_new, C_zone, sdc_iteration);

        // we solved our system to some tolerance, but let's be sure we are conservative by
        // reevaluating the reactions and { doing the full step update
        single_zone_react_source(U_new, R_full, burn_state);
    }

    for (int n = 0; n < NUM_STATE; ++n) {
        U_new[n] = U_old[n] + dt_m * R_full[n] + dt_m * C_zone[n];

        // copy back to k_n
        k_n(i,j,k,n) = U_new[n];
    }
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
sdc_update_centers_o4(const int i, const int j, const int k,
                      Array4<const Real> const& U_old,
                      Array4<Real> const& U_new,
                      Array4<const Real> const& C,
                      const Real dt_m,
                      const int sdc_iteration) {
    // Update U_old to U_new on cell-centers.  This is an implicit
    // solve because of reactions.  Here U_old corresponds to time node
    // m and U_new is node m+1.  dt_m is the timestep between m and
    // m+1

    // We come in with U_new being a guess for the updated solution
    if (okay_to_burn(i, j, k, U_old)) {
        sdc_solve(i, j, k, dt_m, U_old, U_new, C, sdc_iteration);
    } else {
        // no reactions, so it is a straightforward update
        for (int n = 0; n < NUM_STATE; ++n) {
            U_new(i,j,k,n) = U_old(i,j,k,n) + dt_m * C(i,j,k,n);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
instantaneous_react(const int i, const int j, const int k,
                    Array4<const Real> const& state,
                    Array4<Real> const& R_source) {

    // convert from cons to prim -- show this be here or in C++-land?
    // or should I do things like we do in burn_state and convert it manually?
    // (in that case, I am not sure if I can assume UTEMP is defined)

    if (okay_to_burn(i, j, k, state)) {
        burn_t burn_state;
        single_zone_react_source(i, j, k, state, R_source, burn_state);
    } else {
        for (int n = 0; n < NUM_STATE; ++n) {
            R_source(i, j, k, n) = 0.0_rt;
        }
    }
}
#endif

#endif
