#include "Castro.H"
#include "Castro_F.H"

#ifdef DIFFUSION
#include "conductivity.H"
#endif

#ifdef MHD
#include "mhd_util.H"
#endif

using namespace amrex;

Real Castro::estdt_cfl(const Real time) {

    // Courant-condition limited timestep

    GpuArray<Real, 3> center;
    ca_get_center(center.begin());

#ifdef ROTATION
    GpuArray<Real, 3> omega;
    get_omega(time, omega.begin());
#endif

    const auto dx = geom.CellSizeArray();

    ReduceOps<ReduceOpMin> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    const MultiFab& stateMF = get_new_data(State_Type);

    const int ltime_integration_method = time_integration_method;
#ifdef ROTATION
    const int ldo_rotation = do_rotation;
    const int lstate_in_rotating_frame = state_in_rotating_frame;
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(stateMF, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& box = mfi.tilebox();

        auto u = stateMF.array(mfi);

        reduce_op.eval(box, reduce_data,
                       [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept->ReduceTuple {
                           Real rhoInv = 1.0_rt / u(i, j, k, URHO);

                           eos_t eos_state;
                           eos_state.rho = u(i, j, k, URHO);
                           eos_state.T = u(i, j, k, UTEMP);
                           eos_state.e = u(i, j, k, UEINT) * rhoInv;
                           for (int n = 0; n < NumSpec; n++) {
                               eos_state.xn[n] = u(i, j, k, UFS + n) * rhoInv;
                           }
                           for (int n = 0; n < NumAux; n++) {
                               eos_state.aux[n] = u(i, j, k, UFX + n) * rhoInv;
                           }

                           eos(eos_input_re, eos_state);

                           // Compute velocity and then calculate CFL timestep.

                           Real ux = u(i, j, k, UMX) * rhoInv;
                           Real uy = u(i, j, k, UMY) * rhoInv;
                           Real uz = u(i, j, k, UMZ) * rhoInv;

#ifdef ROTATION
                           if (ldo_rotation == 1 && lstate_in_rotating_frame != 1) {
                               Real vel[3];
                               vel[0] = ux;
                               vel[1] = uy;
                               vel[2] = uz;

                               GeometryData geomdata = geom.data();

                               inertial_to_rotational_velocity_c(i, j, k, geomdata, center.begin(),
                                                                 omega.begin(), time, vel);

                               ux = vel[0];
                               uy = vel[1];
                               uz = vel[2];
                           }
#endif

                           Real c = eos_state.cs;

                           Real dt1 = dx[0] / (c + std::abs(ux));

                           Real dt2;
#if AMREX_SPACEDIM >= 2
                           dt2 = dx[1] / (c + std::abs(uy));
#else
                           dt2 = dt1;
#endif

                           Real dt3;
#if AMREX_SPACEDIM == 3
                           dt3 = dx[2] / (c + std::abs(uz));
#else
                           dt3 = dt1;
#endif

                           // The CTU method has a less restrictive timestep than MOL-based
                           // schemes (including the true SDC).  Since the simplified SDC
                           // solver is based on CTU, we can use its timestep.
                           if (ltime_integration_method == 0 || ltime_integration_method == 3) {
                               return {amrex::min(dt1, dt2, dt3)};

                           } else {
                               // method of lines-style constraint is tougher
                               Real dt_tmp = 1.0_rt / dt1;
#if AMREX_SPACEIM >= 2
                               dt_tmp += 1.0_rt / dt2;
#endif
#if AMREX_SPACEDIM == 3
                               dt_tmp += 1.0_rt / dt3;
#endif

                               return 1.0_rt / dt_tmp;
                           }
                       });
    }

    ReduceTuple hv = reduce_data.value();
    Real estdt_hydro = amrex::get<0>(hv);

    return estdt_hydro;
}

#ifdef MHD
Real Castro::estdt_mhd() {

    // MHD timestep limiter
    const auto dx = geom.CellSizeArray();

    ReduceOps<ReduceOpMin> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    const MultiFab& state = get_new_data(State_Type);

    const MultiFab& bx = get_new_data(Mag_Type_x);
    const MultiFab& by = get_new_data(Mag_Type_y);
    const MultiFab& bz = get_new_data(Mag_Type_z);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& box = mfi.tilebox();

        auto u_arr = state.array(mfi);

        auto bx_arr = bx.array(mfi);
        auto by_arr = by.array(mfi);
        auto bz_arr = bz.array(mfi);

        reduce_op.eval(box, reduce_data,
                       [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept->ReduceTuple {
                           Real rhoInv = 1.0_rt / u_arr(i, j, k, URHO);
                           Real bcx = 0.5_rt * (bx_arr(i + 1, j, k) + bx_arr(i, j, k));
                           Real bcy = 0.5_rt * (by_arr(i, j + 1, k) + by_arr(i, j, k));
                           Real bcz = 0.5_rt * (bz_arr(i, j, k + 1) + bz_arr(i, j, k));

                           Real ux = u_arr(i, j, k, UMX) * rhoInv;
                           Real uy = u_arr(i, j, k, UMY) * rhoInv;
                           Real uz = u_arr(i, j, k, UMZ) * rhoInv;

                           eos_t eos_state;
                           eos_state.rho = u_arr(i, j, k, URHO);
                           eos_state.e = u_arr(i, j, k, UEINT) * rhoInv;
                           eos_state.T = u_arr(i, j, k, UTEMP);
                           for (int n = 0; n < NumSpec; n++) {
                               eos_state.xn[n] = u_arr(i, j, k, UFS + n) * rhoInv;
                           }
                           for (int n = 0; n < NumAux; n++) {
                               eos_state.aux[n] = u_arr(i, j, k, UFX + n) * rhoInv;
                           }

                           eos(eos_input_re, eos_state);

                           Real e = u_arr(i, j, k, UEINT) * rhoInv;

                           Real as = eos_state.gam1 * eos_state.p * rhoInv;
                           Real ca = (bcx * bcx + bcy * bcy + bcz * bcz) * rhoInv;

                           Real cx;
                           Real cy;
                           Real cz;

                           if (e > 0_rt) {
                               Real cad = bcx * bcx * rhoInv;
                               eos_soundspeed_mhd(cx, as, ca, cad);

                               cad = bcy * bcy * rhoInv;
                               eos_soundspeed_mhd(cy, as, ca, cad);

                               cad = bcz * bcz * rhoInv;
                               eos_soundspeed_mhd(cz, as, ca, cad);

                           } else {
                               cx = 0.0_rt;
                               cy = 0.0_rt;
                               cz = 0.0_rt;
                           }

                           Real dt1 = dx[0] / (cx + std::abs(ux));

                           Real dt2;
#if AMREX_SPACEDIM >= 2
                           dt2 = dx[1] / (cy + std::abs(uy));
#else
                           dt2 = dt1;
#endif

                           Real dt3;
#if AMREX_SPACEDIM == 3
                           dt3 = dx[2] / (cz + std::abs(uz));
#else
                           dt3 = dt1;
#endif

                           return {amrex::min(dt1, dt2, dt3)};
                       });
    }

    ReduceTuple hv = reduce_data.value();
    Real estdt_mhd = amrex::get<0>(hv);

    return estdt_mhd;
}
#endif

#ifdef DIFFUSION

Real Castro::estdt_temp_diffusion(void) {

    // Diffusion-limited timestep
    //
    // dt < 0.5 dx**2 / D
    // where D = k/(rho c_v), and k is the conductivity

    const auto dx = geom.CellSizeArray();

    ReduceOps<ReduceOpMin> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    const MultiFab& stateMF = get_new_data(State_Type);

    const Real ldiffuse_cutoff_density = diffuse_cutoff_density;
    const Real lmax_dt = max_dt;
    const Real lcfl = cfl;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(stateMF, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& box = mfi.tilebox();

        auto ustate = stateMF.array(mfi);

        reduce_op.eval(box, reduce_data,
                       [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept->ReduceTuple {
                           if (ustate(i, j, k, URHO) > ldiffuse_cutoff_density) {

                               Real rho_inv = 1.0_rt / ustate(i, j, k, URHO);

                               // we need cv
                               eos_t eos_state;
                               eos_state.rho = ustate(i, j, k, URHO);
                               eos_state.T = ustate(i, j, k, UTEMP);
                               eos_state.e = ustate(i, j, k, UEINT) * rho_inv;
                               for (int n = 0; n < NumSpec; n++) {
                                   eos_state.xn[n] = ustate(i, j, k, UFS + n) * rho_inv;
                               }
                               for (int n = 0; n < NumAux; n++) {
                                   eos_state.aux[n] = ustate(i, j, k, UFX + n) * rho_inv;
                               }

                               eos(eos_input_re, eos_state);

                               // we also need the conductivity
                               conductivity(eos_state);

                               // maybe we should check (and take action) on negative cv here?
                               Real D = eos_state.conductivity * rho_inv / eos_state.cv;

                               Real dt1 = 0.5_rt * dx[0] * dx[0] / D;

                               Real dt2;
#if AMREX_SPACEDIM >= 2
                               dt2 = 0.5_rt * dx[1] * dx[1] / D;
#else
                               dt2 = dt1;
#endif

                               Real dt3;
#if AMREX_SPACEDIM >= 3
                               dt3 = 0.5_rt * dx[2] * dx[2] / D;
#else
                               dt3 = dt1;
#endif

                               return {amrex::min(dt1, dt2, dt3)};

                           } else {
                               return lmax_dt / lcfl;
                           }
                       });
    }

    ReduceTuple hv = reduce_data.value();
    Real estdt_diff = amrex::get<0>(hv);

    return estdt_diff;
}
#endif
