#ifndef INITIAL_MODEL_H
#define INITIAL_MODEL_H

AMREX_INLINE
void generate_initial_model (Real mass_want, Real surface_temperature_want,
                             Real& central_density, Real& radius,
                             const Real comp[NumSpec], Real dx,
                             int model_index)
{
    // Integrate HSE to obtain an isentropic profile with the desired mass and surface temperature.

    // Check to make sure the target mass makes sense.

    if (mass_want <= 0.0_rt) {
        amrex::Error("Target mass must be positive");
    }

    // Check to make sure the initial temperature makes sense.

    if (surface_temperature_want <= 0.0_rt) {
        amrex::Error("Target surface temperature must be positive");
    }

    if (surface_temperature_want < castro::small_temp) {
        amrex::Error("Target surface temperature is less than small_temp");
    }

    // Our unknowns are the initial central density and temperature.
    // Since we're integrating outward, we need to do a secant iteration
    // over the central density and temperature. rho_c_old is the 'old'
    // guess for the central density and rho_c is the current guess.
    // After two loops, we can start estimating the central density
    // and temperature required to yield our desired characteristics.

    const int max_hse_iter = 250;
    int max_mass_iter = 250;

    Real temp_c, temp_c_old, dtemp_c;
    Real rho_c, rho_c_old, drho_c;
    Real mass, mass_old;
    Real surface_temperature, surface_temperature_old;

    rho_c_old = -1.0_rt;
    rho_c     = 2.0_rt * eos_rho_0;

    temp_c_old = -1.0_rt;
    temp_c     = 10000.0_rt;

    const Real mass_tol = 1.0e-4;
    const Real temp_tol = 1.0e-4;

    bool mass_converged = false;

    model::initial_model_t& model = model::profile(model_index);

    for (int mass_iter = 1; mass_iter <= max_mass_iter; ++mass_iter) {

        bool fluff = false;

        // We start at the center and integrate outward. Initialize
        // the central conditions.

        model.state(0, model::itemp) = temp_c;
        model.state(0, model::idens) = rho_c;
        for (int n = 0; n < NumSpec; ++n) {
            model.state(0, model::ispec + n) = comp[n];
        }

        eos_t eos_state;
        eos_state.rho  = model.state(0, model::idens);
        eos_state.T    = model.state(0, model::itemp);
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = model.state(0, model::ispec + n);
        }

        eos(eos_input_rt, eos_state);

        model.state(0, model::ipres) = eos_state.p;

        model.r(0) = 0.5 * dx;

        int icutoff = NPTS_MODEL;

        // Make the initial guess be completely uniform.

        for (int i = 1; i < NPTS_MODEL; ++i) {
            model.state(i, model::idens) = model.state(0, model::idens);
            model.state(i, model::itemp) = model.state(0, model::itemp);
            model.state(i, model::ipres) = model.state(0, model::ipres);
            for (int n = 0; n < NumSpec; ++n) {
                model.state(i, model::ispec + n) = model.state(0, model::ispec + n);
            }
            model.r(i) = model.r(i-1) + dx;
        }

        // Keep track of the mass enclosed below the current zone.

        Real rl = 0.0_rt;
        Real rr = rl + dx;

        Real M_enclosed = (4.0_rt / 3.0_rt) * M_PI * (std::pow(rr, 3) - std::pow(rl, 3)) * model.state(0, model::idens);
        mass = M_enclosed;

        //-------------------------------------------------------------------------
        // HSE solve
        //-------------------------------------------------------------------------
        for (int i = 1; i < NPTS_MODEL; ++i) {

            rl += dx;
            rr += dx;

            // As the initial guess for the density and temperature, use the underlying zone.

            model.state(i, model::idens) = model.state(i-1, model::idens);
            model.state(i, model::itemp) = model.state(i-1, model::itemp);

            for (int n = 0; n < NumSpec; ++n) {
                model.state(i, model::ispec + n) = comp[n];
                eos_state.xn[n] = model.state(i, model::ispec + n);
            }

            // Call the EOS to get the internal energy and pressure in the underlying zone.

            Real rho_0 = model.state(i-1, model::idens);
            Real T_0 = model.state(i-1, model::itemp);

            eos(eos_input_rt, eos_state);

            Real e_0 = eos_state.e;
            Real p_0 = eos_state.p;

            Real grav = -C::Gconst * M_enclosed / (std::pow(rl, 2));

            //----------------------------------------------------------------------
            // Iteration loop
            //----------------------------------------------------------------------

            // Start off the Newton loop by assuming that the zone has not converged.

            bool converged_hse = false;

            for (int hse_iter = 1; hse_iter <= max_hse_iter; ++hse_iter) {

                if (fluff) {
                    model.state(i, model::idens) = ambient::ambient_state[URHO];
                    eos_state.rho = ambient::ambient_state[URHO];
                    break;
                }

                // The isentropic constraint means we need to simultaneously satisfy
                // the standard HSE constraint on dP/dr and also de/drho = P / rho**2.
                // Using a straightforward discretization of both, the equations to satisfy are:
                //
                // 0 = f(rho, e) = p(rho, e) - rho**2 * (2 * (e - e_0) / (rho - rho_0) - p_0 / rho_0**2)
                // 0 = g(rho, e) = p(rho, e) - p_0 - (rho + rho_0) * g * dr / 2
                //
                // where rho, e are the unknowns and rho_0, e_0 are the values in the zone below.
                // This is a set of two nonlinear equations in two unknowns. We can solve it using
                //
                // J(x_{n}) * (x_{n+1} - x_{n}) = -y(x_{n})
                //
                // where x is the vector [rho, e], y is [f, g], J is the Jacobian
                //
                // [ df/drho  df/de ]
                // [ dg/drho  dg/de ]
                //
                // and n refers to the current iteration of the nonlinear solve.
                //
                // Using Cramer's rule, the solution is:
                //
                // drho = -(f * dg/de   - g * df/de  ) / (df/drho * dg/de - df/de * dg/drho)
                // de   = -(g * df/drho - f * dg/drho) / (df/drho * dg/de - df/de * dg/drho)
                //
                // Note that we are storing temperature, not energy, in the model, so we need to call
                // the EOS to get energy given (rho, T). Similarly, once we decide the new energy,
                // we call the EOS in the inverse mode to get the T to store.

                Real rho = model.state(i, model::idens);
                Real T = model.state(i, model::itemp);

                eos_state.rho = rho;
                eos_state.T   = T;

                eos(eos_input_rt, eos_state);

                Real e = eos_state.e;
                Real p = eos_state.p;

                Real dpde = eos_state.dpde;
                Real dpdr = eos_state.dpdr_e;

                Real denom = rho - rho_0;
                if (denom == 0.0_rt) {
                    denom = 1.e-30_rt;
                }

                Real f = p - rho * rho * (2.0_rt * (e - e_0) / denom - p_0 / (rho_0 * rho_0));
                Real g = p - p_0 - 0.5_rt * dx * (rho + rho_0) * grav;

                Real dfde = dpde - 2.0_rt * rho * rho / denom;
                Real dfdr = dpdr - 2.0_rt * rho * (2.0_rt * (e - e_0) / denom - p_0 / (rho_0 * rho_0)) +
                            2.0_rt * rho * rho * (e - e_0) / (denom * denom);

                Real dgde = dpde;
                Real dgdr = dpdr - 0.5_rt * dx * grav;

                denom = dfdr * dgde - dfde * dgdr;
                if (denom == 0.0_rt) {
                    denom = 1.e-30_rt;
                }

                Real drho = -(f * dgde - g * dfde) / denom;
                Real de   = -(g * dfdr - f * dgdr) / denom;

                // Constrain drho and de so they don't change the state by too much.

                drho = amrex::max(-0.1_rt * rho, amrex::min(drho, 0.1_rt * rho));
                de = amrex::max(-0.1_rt * e, amrex::min(de, 0.1_rt * e));

                // Get the updated temperature from the new (rho, e) pair.

                eos_state.rho = rho + drho;
                eos_state.e = e + de;

                eos(eos_input_re, eos_state);

                Real dT = eos_state.T - T;

                model.state(i, model::idens) += drho;
                model.state(i, model::itemp) += dT;

                // The model stops when rho < rho_0.

                if (model.state(i, model::idens) < eos_rho_0) {
                    icutoff = i;
                    fluff = true;
                }

                if (std::abs(drho) < model::hse_tol * model.state(i, model::idens)) {
                    converged_hse = true;
                    break;
                }

            }

            if (!converged_hse && (!fluff)) {

                std::cout << "Error: zone " <<  i << " did not converge in init_hse()" << std::endl;
                std::cout << model.state(i, model::idens) << " " << model.state(i, model::itemp) << std::endl;
                std::cout << " p " << eos_state.p << " ";
                std::cout << model::hse_tol * model.state(i, model::idens);
                amrex::Error("Error: HSE non-convergence.");

            }

            // Call the EOS to establish the final properties of this zone.

            eos(eos_input_rt, eos_state);

            model.state(i, model::ipres) = eos_state.p;

            // Discretize the mass enclosed as (4 pi / 3) * rho * dr * (rl**2 + rl * rr + rr**2).

            Real dM = (4.0_rt / 3.0_rt) * M_PI * model.state(i, model::idens) * dx *
                      (rr * rr + rl * rr + rl * rl);
            M_enclosed += dM;

            if (i <= icutoff) {
                // Also update the final mass if we're not in the ambient material,
                // as well as the current estimate of the surface temperature.
                mass += dM;
                surface_temperature = model.state(i, model::itemp);
            }

        } // End loop over zones

        radius = model.r(icutoff);

        if (rho_c_old < 0.0_rt) {

            // Not enough iterations yet -- use an arbitrary guess for the next iteration.

            rho_c_old = rho_c;
            rho_c = 0.9_rt * rho_c_old;

            temp_c_old = temp_c;
            temp_c = 0.9_rt * temp_c_old;

        }
        else {

            // Check if we have converged.

            if (std::abs(mass - mass_want) / mass_want < mass_tol &&
                std::abs(surface_temperature_want - surface_temperature) / surface_temperature_want < temp_tol) {
                mass_converged = true;
                break;
            }

            // Do a secant iteration:
            // M_tot = M(rho_c) + dM/drho |_rho_c x drho + ...

            drho_c = (mass_want - mass) / ((mass - mass_old) / (rho_c - rho_c_old));

            rho_c_old = rho_c;
            rho_c = amrex::min(1.1e0_rt * rho_c_old, amrex::max((rho_c + drho_c), 0.9e0_rt * rho_c_old));

            dtemp_c = (surface_temperature_want - surface_temperature) /
                      ((surface_temperature - surface_temperature_old) / (temp_c - temp_c_old));

            temp_c_old = temp_c;
            temp_c = amrex::min(1.1e0_rt * temp_c_old, amrex::max((temp_c + dtemp_c), 0.9e0_rt * temp_c_old));

        }

        mass_old = mass;
        surface_temperature_old = surface_temperature;

    } // End mass constraint loop

    if (!mass_converged && max_mass_iter > 1) {
        amrex::Error("ERROR: Mass did not converge.");
    }

    central_density = model.state(0, model::idens);

    model::initialized = true;
    model::npts = NPTS_MODEL;
}

#endif
