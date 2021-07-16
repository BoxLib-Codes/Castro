#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>
#include <network.H>
#ifdef NSE_THERMO
#include <nse.H>
#endif

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{

    int coord_type = geomdata.Coord();

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();
    const Real* probhi = geomdata.ProbHi();

    Real width = problem::w_T * (probhi[0] - problo[0]);

    Real height = probhi[1] - problo[1];

    Real c_T = problo[0] + problem::center_T * (probhi[0] - problo[0]) * (1.0_rt + problem::pert_amplitude * std::sin(2.0_rt * M_PI * 2.0_rt * dx[1] * (static_cast<Real>(j) + 0.5_rt) / height));

    Real xcen = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);

    Real sigma = 1.0_rt / (1.0_rt + std::exp(-(c_T - xcen)/ width));

    state(i,j,k,URHO) = problem::dens_ash + (problem::dens_fuel - problem::dens_ash) * (1.0_rt - sigma);

    state(i,j,k,UTEMP) = problem::T_ash + (problem::T_fuel - problem::T_ash) * (1.0_rt - sigma);

    Real xn[NumSpec] ;

    for (int n = 0; n < NumSpec; n++) {
        xn[n] = problem::xn_ash[n] + (problem::xn_fuel[n] - problem::xn_ash[n]) * (1.0_rt - sigma);
}

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = state(i,j,k,URHO) * xn[n];
        problem::xn[n] = problem::xn_ash[n] + (problem::xn_fuel[n] - problem::xn_ash[n]) * (1.0_rt - sigma);
        state(i,j,k,UFS+n) = state(i,j,k,URHO) * problem::xn[n];
    }

    eos_t eos_state;
    eos_state.rho = state(i,j,k,URHO);
    eos_state.T = state(i,j,k,UTEMP);
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = xn[n];
    }

#ifdef NSE_THERMO
    // set the aux quantities -- we need to do this if we are using the NSE network
    set_nse_aux_from_X(eos_state);

    for (int n = 0; n < NumAux; n++) {
        state(i,j,k,UFX+n) = state(i,j,k,URHO) * eos_state.aux[n];
    }
#endif

    eos(eos_input_rt, eos_state);

    state(i,j,k,UMX) = state(i,j,k,URHO) * (problem::vel - 2 * problem::vel * (1.0_rt - sigma));
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;
    state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state.e;
    state(i,j,k,UEDEN) = state(i,j,k,UEINT) + 0.5_rt * state(i,j,k,UMX) * state(i,j,k,UMX) / state(i,j,k,URHO);

}

#endif
