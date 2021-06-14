#ifndef CASTRO_REACT_UTIL_H
#define CASTRO_REACT_UTIL_H

#include <Castro.H>
#include <burn_type.H>
#include <numerical_jacobian.H>


constexpr int iwrho = 0;
constexpr int iwfs = 1;
constexpr int iwe = NumSpec+1;

// #include <actual_rhs.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool
okay_to_burn(const int i, const int j, const int k,
             Array4<const Real> const& state) {

    if (state(i,j,k,UTEMP) < react_T_min || state(i,j,k,UTEMP) > react_T_max ||
        state(i,j,k,URHO) < react_rho_min || state(i,j,k,URHO) > react_rho_max) {
            return false;
    }

    return true;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool
okay_to_burn(GpuArray<Real, NUM_STATE> const& state) {

    if (state[UTEMP] < react_T_min || state[UTEMP] > react_T_max ||
        state[URHO] < react_rho_min || state[URHO] > react_rho_max) {
        return false;
    }

    return true;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool
okay_to_burn_type(burn_t const& state) {

    if (state.T < react_T_min || state.T > react_T_max ||
        state.rho < react_rho_min || state.rho > react_rho_max) {
            return false;
    }

    return true;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
single_zone_react_source(GpuArray<Real, NUM_STATE> const& state,
                         GpuArray<Real, NUM_STATE>& R,
                         burn_t& burn_state) {

    // note: we want this burn_state to come out of here so we can
    // reuse it elsewhere

    auto rhoInv = 1.0_rt / state[URHO];

    burn_state.rho = state[URHO];
    burn_state.T = state[UTEMP];
    burn_state.e = state[UEINT] * rhoInv;

    for (int n = 0; n < NumSpec; ++n) {
       burn_state.xn[n] = amrex::max(amrex::min(state[UFS+n] * rhoInv, 1.0_rt), small_x);
    }

#if NAUX_NET > 0
    for (int n = 0; n < NumAux; ++n) {
       burn_state.aux[n] = state[UFX+n] * rhoInv;
    }
#endif

    eos_t eos_state;

    // Ensure that the temperature going in is consistent with the internal energy.
    burn_to_eos(burn_state, eos_state);
    eos(eos_input_re, eos_state);
    eos_to_burn(eos_state, burn_state);

    // eos_get_small_temp(&small_temp);
    burn_state.T = amrex::min(MAX_TEMP, amrex::max(burn_state.T, small_temp));

#ifdef STRANG
    burn_state.self_heat = false;
#endif

    Array1D<Real, 1, neqs> ydot;

    actual_rhs(burn_state, ydot);

    // store the instantaneous R
    for (int n = 0; n < NUM_STATE; ++n) {
        R[n] = 0.0_rt;
    }

    // species rates come back in terms of molar fractions
    for (int n = 0; n < NumSpec; ++n) {
        R[UFS+n] = state[URHO] * aion[n] * ydot(1+n);
    }

    R[UEDEN] = state[URHO] * ydot(net_ienuc);
    R[UEINT] = state[URHO] * ydot(net_ienuc);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
single_zone_react_source(const int i, const int j, const int k,
                         Array4<const Real> const& state,
                         Array4<Real> const& R,
                         burn_t& burn_state) {
    GpuArray<Real, NUM_STATE> state_arr;
    GpuArray<Real, NUM_STATE> R_arr;

    for (int n = 0; n < NUM_STATE; ++n) {
        state_arr[n] = state(i,j,k,n);
        R_arr[n] = R(i,j,k,n);
    }

    single_zone_react_source(state_arr, R_arr, burn_state);

    for (int n = 0; n < NUM_STATE; ++n) {
        R(i,j,k,n) = R_arr[n];
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
single_zone_jac(GpuArray<Real, NUM_STATE> const& state,
                burn_t burn_state,
                Array2D<Real, 0, NumSpec+1, 0, NumSpec+1>& dRdw) {

    // we enter this routine if we are doing the Newton solver or VODE
    // with an analytic Jacobian.
    //
    // For Newton, we select the analytic or numerical Jacobian
    // evalution via sdc_newton_use_analytic_jac.
    //
    // For VODE, we enter this for integrator.jacobian = 1 or 3.  In
    // both cases, the Jacobian is expressed as dR/dU = dR/dw dw/dU,
    // and we always compute dw/dU analytically.  For
    // integrator.jacobian = 1, we will compute dR/dw analytically,
    // using the network.  Otherwise, we will compute it numerically,
    // via differencing.

    // Regardless of the above, the derivatives wrt density are always
    // computed numerically.

#ifdef SIMPLIFIED_SDC
#ifndef AMREX_USE_GPU
    amrex::Error("we shouldn't be here with the simplified SDC method (USE_SIMPLIFIED_SDC=TRUE)");
#endif
#else

    Array1D<Real, 1, neqs> ydot;
    Array1D<Real, 1, neqs> ydot_pert;
    JacNetArray2D Jac;

    const auto eps = 1.e-8_rt;

    actual_rhs(burn_state, ydot);

    // Jac has the derivatives with respect to the native
    // network variables, X, e.  It does not have derivatives with
    // respect to density, so we'll have to compute those ourselves.

    if (sdc_newton_use_analytic_jac == 0 || jacobian == 3) {
        // note the numerical Jacobian will be returned in terms of X
        numerical_jac(burn_state, Jac);
    } else {
        actual_jac(burn_state, Jac);

        // The Jacobian from the nets is in terms of dYdot/dY, but we want
        // it was dXdot/dX, so convert here.
        for (int n = 1; n <= NumSpec; n++) {
            for (int m = 1; m <= neqs; m++) {
                // d(X_n)/dq_m
                Jac(n,m) = Jac(n,m) * aion[n-1];
            }
        }

        for (int m = 1; m <= neqs; m++) {
            for (int n = 1; n <= NumSpec; n++) {
                // d RHS_m / dX_n
                Jac(m,n) = Jac(m,n) * aion_inv[n-1];
            }
        }
    }
#endif

    // at this point, our Jacobian should be entirely in terms of X,
    // not Y.  Let's now fix the rhs terms themselves to be in terms of
    // dX/dt and not dY/dt.
    for (int n = 0; n < NumSpec; ++n) {
        ydot(n+1) *= aion[n];
    }

    // Our jacobian, dR/dw has the form:
    //
    //  /      0                  0                  0                       0          \
    //  | d(rho X1dot)/drho  d(rho X1dot)/dX1   d(rho X1dit)/dX2   ...  d(rho X1dot)/dT |
    //  | d(rho X2dot)/drho  d(rho X2dot)/dX1   d(rho X2dot)/dX2   ...  d(rho X2dot)/dT |
    //  |   ...                                                                         |
    //  \ d(rho Edot)/drho   d(rho Edot)/dX1    d(rho Edot)/dX2    ...  d(rho Edot)/dT  /

    for (int n = 0; n < NumSpec+2; ++n) {
        for (int m = 0; m < NumSpec+2; ++m) {
            dRdw(n,m) = 0.0_rt;
        }
    }

    // now perturb density and call the RHS to compute the derivative wrt rho
    // species rates come back in terms of molar fractions
    burn_t burn_state_pert = burn_state;

    burn_state_pert.rho = burn_state.rho * (1.0_rt + eps);
    burn_state_pert.T = burn_state.T;
    burn_state_pert.e = burn_state.e;
    for (int n = 0; n < NumSpec; ++n) {
        burn_state_pert.xn[n] = burn_state.xn[n];
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; ++n) {
        burn_state_pert.aux[n] = burn_state.aux[n];
    }
#endif

    actual_rhs(burn_state_pert, ydot_pert);

    // make the rates dX/dt and not dY/dt
    for (int n = 0; n < NumSpec; ++n) {
        ydot_pert(n+1) *= aion[n];
    }

    // fill the column of dRdw corresponding to the derivative
    // with respect to rho
    for (int m = 1; m <= NumSpec; ++m) {
       // d( d(rho X_m)/dt)/drho
       dRdw(m, iwrho) = ydot(m) +
            state[URHO] * (ydot_pert(m) - ydot(m))/(eps * burn_state.rho);
    }

    // d( d(rho E)/dt)/drho
    dRdw(NumSpec+1, iwrho) = ydot(net_ienuc) +
         state[URHO] * (ydot_pert(net_ienuc) - ydot(net_ienuc))/(eps * burn_state.rho);

    // fill the columns of dRdw corresponding to each derivative
    // with respect to species mass fraction
    for (int n = 1; n <= NumSpec; ++n) {
       dRdw(0, iwfs-1+n) = 0.0_rt;  // density source

       for (int m = 1; m <= NumSpec; ++m) {
          // d( d(rho X_m)/dt)/dX_n
          dRdw(m, iwfs-1+n) = state[URHO] * Jac(m, n);
       }

       // d( d(rho E)/dt)/dX_n
       dRdw(NumSpec+1, iwfs-1+n) = state[URHO] * Jac(net_ienuc, n);

    }

    // now fill the column corresponding to derivatives with respect to
    // energy -- this column is iwe
    dRdw(0, iwe) = 0.0_rt;

    // d( d(rho X_m)/dt)/de
    for (int m = 1; m <= NumSpec; ++m) {
       dRdw(m, iwe) = state[URHO] * Jac(m, net_ienuc);
    }

    // d( d(rho E)/dt)/de
    dRdw(NumSpec+1, iwe) = state[URHO] * Jac(net_ienuc, net_ienuc);
}

#endif
