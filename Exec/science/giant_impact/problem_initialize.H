#ifndef problem_initialize_H
#define problem_initialize_H

#include <model_parser.H>
#include <initial_model.H>
#include <fundamental_constants.H>
#include <prob_parameters.H>

const Real M_earth = 5.9722e27_rt;

// Given the total mass of a binary system, obtain the velocity at a
// given separation assuming the point masses fell in from infinity
// starting from zero velocity. Note this is in the center of mass frame.

AMREX_INLINE
Real freefall_velocity (Real mass, Real distance)
{
    return std::sqrt(2.0_rt * C::Gconst * mass / distance);
}

AMREX_INLINE
void problem_initialize ()
{
    // Set some default values for these quantities;
    // we'll update them soon.

    for (int n = 0; n < 3; ++n) {
        problem::target_initial_center[n] = problem::center[n];
        problem::impactor_initial_center[n] = problem::center[n];

        problem::target_initial_velocity[n] = 0.0_rt;
        problem::impactor_initial_velocity[n] = 0.0_rt;
    }

    // Integrate to reach the desired mass for both target and impactor.

    if (problem::target_mass <= 0.0_rt) {
        amrex::Error("Target mass must be positive");
    }
    problem::target_mass *= M_earth;

    Real target_comp[NumSpec];

    // Arbitrarily put all the composition in the first component

    for (int n = 0; n < NumSpec; ++n) {
        target_comp[n] = 0.0_rt;
    }
    target_comp[0] = 1.0_rt;

    Real target_central_density;

    generate_initial_model(problem::target_mass, problem::target_temperature,
                           target_central_density, problem::target_radius,
                           target_comp, problem::initial_model_dx,
                           0);

    amrex::Print() << std::endl;

    amrex::Print() << "Generated initial model for target of mass " << std::setprecision(3) << problem::target_mass / M_earth
                   << " earth masses, central density " << std::setprecision(3) << std::scientific << target_central_density
                   << " g cm**-3, and radius " << std::setprecision(3) << std::scientific << problem::target_radius << " cm."
                   << std::endl << std::endl;

    Real impactor_comp[NumSpec];

    // Arbitrarily put all the composition in the first component

    for (int n = 0; n < NumSpec; ++n) {
        impactor_comp[n] = 0.0_rt;
    }
    impactor_comp[0] = 1.0_rt;

    if (problem::impactor_mass > 0.0_rt) {
        problem::impactor_mass *= M_earth;

        Real impactor_central_density;

        generate_initial_model(problem::impactor_mass, problem::impactor_temperature,
                               impactor_central_density, problem::impactor_radius,
                               impactor_comp, problem::initial_model_dx,
                               1);

        amrex::Print() << "Generated initial model for impactor of mass " << std::setprecision(3) << problem::impactor_mass / M_earth
                       << " earth masses, central density " << std::setprecision(3) << std::scientific << impactor_central_density
                       << " g cm**-3, and radius " << std::setprecision(3) << std::scientific << problem::impactor_radius << " cm."
                       << std::endl << std::endl;
    }
    else {
        amrex::Print() << "No impactor provided; simulating target only." << std::endl << std::endl;
    }

    // Set ambient properties
    
    for (int n = 0; n < NumSpec; ++n) {
        ambient::ambient_state[UFS+n] = ambient::ambient_state[URHO] * (target_comp[n] + impactor_comp[n]) / 2;
    }

    // We have completed (rho, T, xn) for the ambient state, so we can call the EOS.

    eos_t eos_state;
    eos_state.rho = ambient::ambient_state[URHO];
    eos_state.T   = ambient::ambient_state[UTEMP];
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = ambient::ambient_state[UFS+n] / ambient::ambient_state[URHO];
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; ++n) {
        eos_state.aux[n] = ambient::ambient_state[UFX+n] / ambient::ambient_state[URHO];
    }
#endif

    eos(eos_input_rt, eos_state);

    ambient::ambient_state[UEINT] = ambient::ambient_state[URHO] * eos_state.e;
    ambient::ambient_state[UEDEN] = ambient::ambient_state[UEINT];

    // Set up the initial distances and velocities.

    if (problem::impactor_mass > 0.0) {

        // Scale the initial distances by the radius of the target. This assumes that
        // the target's radius is larger than the impactor's.

        AMREX_ALWAYS_ASSERT(problem::target_radius >= problem::impactor_radius);

        Real initial_distance_x = problem::initial_distance * problem::target_radius;
        Real initial_distance_y = problem::impact_parameter * problem::target_radius;

        // Compute the initial velocity of the impactor with respect to the target.
        // This has two components: the freefall velocity due to the gravitational
        // force (assumed to have started acting at a distance of infinity) and also
        // the initial velocity at infinity.

        Real total_mass = problem::target_mass + problem::impactor_mass;
        Real v_ff = freefall_velocity(total_mass, initial_distance_x);

        // The initial velocity is provided as a fraction of the escape velocity of the
        // system, as defined by Meier, Reinhardt and Stadel (2021). This happens to be
        // the same formula as the freefall velocity evaluated at a distance of the
        // critical radius for impact (the sum of the target and impactor radii).

        Real critical_radius = problem::target_radius + problem::impactor_radius;
        Real v_esc = freefall_velocity(total_mass, critical_radius);

        Real impact_velocity = v_ff + problem::initial_velocity * v_esc;

        // Set the velocity in the center of mass frame. Note that the fraction of the
        // velocity for each object is inversely proportional to the fraction of the system
        // mass it has.

        problem::target_initial_velocity[0] = (problem::impactor_mass / total_mass) * impact_velocity;
        problem::impactor_initial_velocity[0] = -(problem::target_mass / total_mass) * impact_velocity;

        // Set the initial locations in the center of mass frame. As above, the less massive
        // object will be further away from the center of mass.

        problem::target_initial_center[0] = -(problem::impactor_mass / total_mass) * initial_distance_x;
        problem::impactor_initial_center[0] = (problem::target_mass / total_mass) * initial_distance_x;

        problem::target_initial_center[1] = -(problem::impactor_mass / total_mass) * initial_distance_y;
        problem::impactor_initial_center[1] = (problem::target_mass / total_mass) * initial_distance_y;
    }

    // Make sure the domain is big enough to hold both the target and impactor.

    const Real* problo = DefaultGeometry().ProbLo();
    const Real* probhi = DefaultGeometry().ProbHi();

    if (problem::target_initial_center[0] - problem::target_radius <= problo[0] ||
        problem::target_initial_center[1] - problem::target_radius <= problo[1] ||
        problem::target_initial_center[2] - problem::target_radius <= problo[2]) {
        amrex::Error("ERROR: The domain width is too small to include the target.");
    }

    if (problem::impactor_mass > 0.0) {
        if (problem::impactor_initial_center[0] + problem::impactor_radius >= probhi[0] ||
            problem::impactor_initial_center[1] + problem::impactor_radius >= probhi[1] ||
            problem::impactor_initial_center[2] + problem::impactor_radius >= probhi[2]) {
            amrex::Error("ERROR: The domain width is too small to include the impactor.");
        }
    }

    amrex::Print() << "Target location: (" << std::setprecision(3) << std::scientific
                   << problem::target_initial_center[0] << ", "
                   << problem::target_initial_center[1] << ", "
                   << problem::target_initial_center[2] << ")"
                   << std::endl;

    amrex::Print() << "Target velocity: (" << std::setprecision(3) << std::scientific
                   << problem::target_initial_velocity[0] << ", "
                   << problem::target_initial_velocity[1] << ", "
                   << problem::target_initial_velocity[2] << ")"
                   << std::endl;

    if (problem::impactor_mass > 0.0) {
        amrex::Print() << "Impactor location: (" << std::setprecision(3) << std::scientific
                       << problem::impactor_initial_center[0] << ", "
                       << problem::impactor_initial_center[1] << ", "
                       << problem::impactor_initial_center[2] << ")"
                       << std::endl;

        amrex::Print() << "Impactor velocity: (" << std::setprecision(3) << std::scientific
                       << problem::impactor_initial_velocity[0] << ", "
                       << problem::impactor_initial_velocity[1] << ", "
                       << problem::impactor_initial_velocity[2] << ")"
                       << std::endl;

        amrex::Print() << std::endl;
    }
}

#endif
