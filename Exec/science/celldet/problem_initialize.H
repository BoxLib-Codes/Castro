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

    // get the species indices

    problem::ihe4 = network_spec_index("helium-4");
    problem::ic12 = network_spec_index("carbon-12");
    problem::io16 = network_spec_index("oxygen-16");

    if (problem::ihe4 < 0 || problem::ic12 < 0 || problem::io16 < 0) {
        amrex::Error("ERROR: species indices not found");
    }

    // make sure that the carbon fraction falls between 0 and 1

    if (problem::cfrac > 1.0_rt || problem::cfrac < 0.0_rt) {
        amrex::Error("ERROR: cfrac must fall between 0 and 1");
    }

    // make sure that the oxygen fraction falls between 0 and 1

    if (problem::ofrac > 1.0_rt || problem::ofrac < 0.0_rt) {
        amrex::Error("ERROR: ofrac must fall between 0 and 1");
    }

    // make sure that the C/O fraction sums to no more than 1

    if (problem::cfrac + problem::ofrac > 1.0_rt) {
        amrex::Error("ERROR: cfrac + ofrac cannot exceed 1.");
    }

    // set the default mass fractions

    for (int n = 0; n < NumSpec; n++) {
        problem::xn[n] = problem::smallx;
    }
    problem::xn[problem::ic12] = amrex::max(problem::cfrac, problem::smallx);
    problem::xn[problem::io16] = amrex::max(problem::ofrac, problem::smallx);
    problem::xn[problem::ihe4] = 1.0_rt - problem::cfrac - problem::ofrac - (NumSpec - 2) * problem::smallx;

    // Set the ambient material

    problem::ambient_dens = problem::dens_l;
    for (int n = 0; n < NumSpec; n++) {
        problem::ambient_comp[n] = problem::xn[n];
    }

    eos_t eos_state;


    eos_state.rho = problem::ambient_dens;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n]  = problem::ambient_comp[n];
    }

    eos_state.T = problem::T_l;

    eos(eos_input_rt, eos_state);

    problem::ambient_e_l = eos_state.e;

    eos_state.T = problem::T_r;

    eos(eos_input_rt, eos_state);

    problem::ambient_e_r = eos_state.e;

}

#endif
