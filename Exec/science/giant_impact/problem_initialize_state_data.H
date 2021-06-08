#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <model_parser.H>
#include <ambient.H>
#include <Castro_util.H>
#include <prob_parameters.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{
    GpuArray<Real, 3> loc;
    position(i, j, k, geomdata, loc);

    Real dist_target = std::sqrt(std::pow(loc[0] - problem::target_initial_center[0], 2) +
                                 std::pow(loc[1] - problem::target_initial_center[1], 2) +
                                 std::pow(loc[2] - problem::target_initial_center[2], 2));

    Real dist_impactor = std::sqrt(std::pow(loc[0] - problem::impactor_initial_center[0], 2) +
                                   std::pow(loc[1] - problem::impactor_initial_center[1], 2) +
                                   std::pow(loc[2] - problem::impactor_initial_center[2], 2));

    const Real* dx = geomdata.CellSize();

    eos_t zone_state;

    if (dist_target < problem::target_radius)
    {
        Real pos[3] = {loc[0] - problem::target_initial_center[0],
                       loc[1] - problem::target_initial_center[1],
                       loc[2] - problem::target_initial_center[2]};

        zone_state.rho = interpolate_3d(pos, dx, model::idens, 1, 0);
        zone_state.T   = interpolate_3d(pos, dx, model::itemp, 1, 0);
        for (int n = 0; n < NumSpec; ++n) {
            zone_state.xn[n] = interpolate_3d(pos, dx, model::ispec + n, 1, 0);
        }

        eos(eos_input_rt, zone_state);

        state(i,j,k,UMX) = problem::target_initial_velocity[0] * zone_state.rho;
        state(i,j,k,UMY) = problem::target_initial_velocity[1] * zone_state.rho;
        state(i,j,k,UMZ) = problem::target_initial_velocity[2] * zone_state.rho;
    }
    else if (dist_impactor < problem::impactor_radius && problem::impactor_mass > 0.0)
    {
        Real pos[3] = {loc[0] - problem::impactor_initial_center[0],
                       loc[1] - problem::impactor_initial_center[1],
                       loc[2] - problem::impactor_initial_center[2]};

        zone_state.rho = interpolate_3d(pos, dx, model::idens, 1, 1);
        zone_state.T   = interpolate_3d(pos, dx, model::itemp, 1, 1);
        for (int n = 0; n < NumSpec; ++n) {
            zone_state.xn[n] = interpolate_3d(pos, dx, model::ispec + n, 1, 1);
        }

        eos(eos_input_rt, zone_state);

        state(i,j,k,UMX) = problem::impactor_initial_velocity[0] * zone_state.rho;
        state(i,j,k,UMY) = problem::impactor_initial_velocity[1] * zone_state.rho;
        state(i,j,k,UMZ) = problem::impactor_initial_velocity[2] * zone_state.rho;
    }
    else
    {
        zone_state.rho = ambient::ambient_state[URHO];
        zone_state.T   = ambient::ambient_state[UTEMP];
        zone_state.e   = ambient::ambient_state[UEINT] / ambient::ambient_state[URHO];
        for (int n = 0; n < NumSpec; ++n) {
            zone_state.xn[n] = ambient::ambient_state[UFS+n] / ambient::ambient_state[URHO];
        }

        state(i,j,k,UMX) = 0.0_rt;
        state(i,j,k,UMY) = 0.0_rt;
        state(i,j,k,UMZ) = 0.0_rt;
    }

    state(i,j,k,URHO)  = zone_state.rho;
    state(i,j,k,UTEMP) = zone_state.T;
    state(i,j,k,UEINT) = zone_state.e * zone_state.rho;
    state(i,j,k,UEDEN) = zone_state.e * zone_state.rho;
    for (int n = 0; n < NumSpec; ++n) {
        state(i,j,k,UFS+n) = zone_state.rho * zone_state.xn[n];
    }

    state(i,j,k,UEDEN) += 0.5_rt * (state(i,j,k,UMX) * state(i,j,k,UMX) +
                                    state(i,j,k,UMY) * state(i,j,k,UMY) +
                                    state(i,j,k,UMZ) * state(i,j,k,UMZ)) / state(i,j,k,URHO);
}

#endif
