#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>
#include <network.H>
#ifdef NSE_THERMO
#include <nse.H>
#endif

AMREX_INLINE
void problem_initialize ()
{

    int ifuel1 = network_spec_index(problem::fuel1_name);
    int iash1 = network_spec_index(problem::ash1_name);

    int ifuel2 = network_spec_index(problem::fuel2_name);
    int iash2 = network_spec_index(problem::ash2_name);

    int ifuel3 = network_spec_index(problem::fuel3_name);
    int iash3 = network_spec_index(problem::ash3_name);

    int ifuel4 = network_spec_index(problem::fuel4_name);
    int iash4 = network_spec_index(problem::ash4_name);

    int ifuel5 = network_spec_index(problem::fuel5_name);
    int iash5 = network_spec_index(problem::ash5_name);

    int ifuel6 = network_spec_index(problem::fuel6_name);
    int iash6 = network_spec_index(problem::ash6_name);

    if (iash1 < 0 && iash2 < 0 && iash3 < 0 && iash4 < 0 && iash5 < 0 && iash6 < 0) {
        amrex::Error("no valid ash state defined");
    }

    // fuel state

    for (int n = 0; n < NumSpec; n++) {
        problem::xn_fuel[n] = problem::smallx;
        problem::xn_ash[n] = problem::smallx;
    }

    if (ifuel1 >= 0) {
        problem::xn_fuel[ifuel1] = problem::X_fuel1;
    }

    if (ifuel2 >= 0) {
        problem::xn_fuel[ifuel2] = problem::X_fuel2;
    }

    if (ifuel3 >= 0) {
        problem::xn_fuel[ifuel3] = problem::X_fuel3;
    }

    if (ifuel4 >= 0) {
        problem::xn_fuel[ifuel4] = problem::X_fuel4;
    }

    if (ifuel5 >= 0) {
        problem::xn_fuel[ifuel5] = problem::X_fuel5;
    }

    if (ifuel6 >= 0) {
        problem::xn_fuel[ifuel6] = problem::X_fuel6;
    }

    if (iash1 >= 0) {
        problem::xn_ash[iash1] = problem::X_ash1;
    }

    if (iash2 >= 0) {
        problem::xn_ash[iash2] = problem::X_ash2;
    }

    if (iash3 >= 0) {
        problem::xn_ash[iash3] = problem::X_ash3;
    }

    if (iash4 >= 0) {
        problem::xn_ash[iash4] = problem::X_ash4;
    }

    if (iash5 >= 0) {
        problem::xn_ash[iash5] = problem::X_ash5;
    }

    if (iash6 >= 0) {
        problem::xn_ash[iash6] = problem::X_ash6;
    }


    // normalize

    Real sum_Xf = 0.0_rt;
    Real sum_Xa = 0.0_rt;
    for (int n = 0; n < NumSpec; n++) {
        sum_Xf += problem::xn_fuel[n];
        sum_Xa += problem::xn_ash[n];
    }
    for (int n = 0; n < NumSpec; n++) {
        problem::xn_fuel[n] /= sum_Xf;
        problem::xn_ash[n] /= sum_Xa;
    }

}

#endif
