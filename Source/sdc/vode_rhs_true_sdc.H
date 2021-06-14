#ifndef VODE_RHS_TRUE_SDC_H
#define VODE_RHS_TRUE_SDC_H

#include <network.H>
#include <burn_type.H>
#include <extern_parameters.H>

#ifdef NETWORK_HAS_CXX_IMPLEMENTATION
#include <actual_network.H>
#include <actual_rhs.H>
#else
#include <fortran_to_cxx_actual_rhs.H>
#endif
#include <Castro_react_util.H>

// The f_rhs routine provides the right-hand-side for the DVODE solver.
// This is a generic interface that calls the specific RHS routine in the
// network you're actually using.

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void rhs(const Real time, burn_t& burn_state, dvode_t& vode_state, RArray1D& dUdt)
{

    // note: dUdt is 1-based

    GpuArray<Real, NUM_STATE> U_full;
    GpuArray<Real, NUM_STATE> R_full;

    GpuArray<Real, NumSpec+2> R_react;

    // evaluate R

    // we are not solving the momentum equations
    // create a full state -- we need this for some interfaces
    // note: vode_state % y(:) is 1-based

    U_full[URHO] = vode_state.y(1);
    for (int n = 0; n < NumSpec; n++) {
        U_full[UFS+n] = vode_state.y(2+n);
    }
    U_full[UEINT] = vode_state.y(NumSpec+2);
    U_full[UEDEN] = burn_state.E_var;

    U_full[UMX] = burn_state.mom[0];
    U_full[UMY] = burn_state.mom[1];
    U_full[UMZ] = burn_state.mom[2];

    // initialize the temperature -- a better value will be found when we do the EOS
    // call in single_zone_react_source

    U_full[UTEMP] = burn_state.T;

    // create a temporary burn_t for this call.  It is just going to
    // get loaded up with U_full

    burn_t burn_state_pass;
    single_zone_react_source(U_full, R_full, burn_state_pass);

    // update our temperature for next time

    burn_state.T = burn_state_pass.T;

    R_react[0] = R_full[URHO];
    for (int n = 0; n < NumSpec; n++) {
        R_react[1+n] = R_full[UFS+n];
    }
    R_react[NumSpec+1] = R_full[UEINT];

    // create the RHS -- this is 1-based

    for (int n = 1; n <= NumSpec+2; n++) {
        dUdt(n) = R_react[n-1] + burn_state.f_source[n-1];
    }
}


template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void jac (burn_t& burn_state, dvode_t& vode_state, MatrixType& pd)
{
    // this is the analytic Jacobian.  This routine is called by
    // the VODE when integrator.jacobian != 2.

    // NOTE: the time at which to evaluate the Jacobian is not
    // explicitly passed. VODE always evaluates the analytic
    // Jacobian at vode_state.tn, and we pass vode_state.

    // Call the specific network routine to get the Jacobian.

    // note the Jacobian is 1-based at the moment

    Array2D<Real, 0, NumSpec+1, 0, NumSpec+1> dRdw = {0.0_rt};
    Array2D<Real, 0, NumSpec+1, 0, NumSpec+1> dwdU = {0.0_rt};

    GpuArray<Real, NUM_STATE> U_full;
    GpuArray<Real, NUM_STATE> R_full;
    GpuArray<Real, NumSpec+2> R_react;

    // we are not solving the momentum equations
    // create a full state -- we need this for some interfaces

    U_full[URHO] = vode_state.y(1);
    for (int n = 0; n < NumSpec; n++) {
        U_full[UFS+n] = vode_state.y(2+n);
    }
    U_full[UEINT] = vode_state.y(NumSpec+2);
    U_full[UEDEN] = burn_state.E_var;

    U_full[UMX] = burn_state.mom[0];
    U_full[UMY] = burn_state.mom[1];
    U_full[UMZ] = burn_state.mom[2];

    // compute the temperature and species derivatives --
    // maybe this should be done using the burn_state
    // returned by single_zone_react_source, since it is
    // more consistent T from e

    eos_t eos_state;
    eos_state.rho = U_full[URHO];
    eos_state.T = burn_state.T; // initial guess
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = U_full[UFS+n] / U_full[URHO];
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = U_full[UFX+n] / U_full[URHO];
    }
#endif
    eos_state.e = U_full[UEINT] / U_full[URHO];  //(U_full(UEDEN) - HALF*sum(U_full(UMX:UMZ))/U_full(URHO))/U_full(URHO)


    eos(eos_input_re, eos_state);

    U_full[UTEMP] = eos_state.T;

    burn_t burn_state_pass;
    single_zone_react_source(U_full, R_full, burn_state_pass);

    single_zone_jac(U_full, burn_state_pass, dRdw);

    // construct dwdU

    // the density row

    dwdU(iwrho, 0) = 1.0_rt;

    // the X_k rows

    for (int m = 1; m < NumSpec+1; m++) {
       dwdU(iwfs-1+m, 0) = -vode_state.y(m+1) / (vode_state.y(1) * vode_state.y(1));
       dwdU(iwfs-1+m, m) = 1.0_rt / vode_state.y(1);
    }

    auto eos_xderivs = composition_derivatives(eos_state);

    // now the e row -- this depends on whether we are evolving (rho E) or (rho e)

    Real denom = 1.0_rt / eos_state.rho;
    Real xn_sum = 0.0_rt;
    for (int n = 0; n < NumSpec; ++n) {
        xn_sum += eos_state.xn[n] * eos_xderivs.dedX[n];
    }

    dwdU(iwe, 0) = denom * (xn_sum - eos_state.e);

    for (int m = 0; m < NumSpec; m++) {
        dwdU(iwe, m+1) = -denom * eos_xderivs.dedX[m];
    }

    dwdU(iwe, iwe) = denom;

    // construct the Jacobian -- note: the Jacobian is 1-based

    for (int n = 0; n <= NumSpec+1; ++n) {
        for (int m = 0; m <= NumSpec+1; ++m) {
            pd(n+1, m+1) = 0.0_rt;
            for (int l = 0; l <= NumSpec+1; ++l) {
                pd(n+1, m+1) = dRdw(n,l) * dwdU(l,m);
            }
        }
    }

}
#endif

