#include <Castro.H>
#include <Castro_F.H>

#ifdef DIFFUSION
#include <conductivity.H>
#endif

#ifdef MHD
#include <mhd_util.H>
#endif

#ifdef ROTATION
#include <Rotation.H>
#endif

#ifdef REACTIONS
#ifdef NETWORK_HAS_CXX_IMPLEMENTATION
#include <actual_rhs.H>
#else
#include <fortran_to_cxx_actual_rhs.H>
#endif
#endif

#ifdef RADIATION
#include <Radiation.H>
#endif

using namespace amrex;

Real
Castro::estdt_cfl(const Real time)
{

  // Courant-condition limited timestep

#ifdef ROTATION
  GeometryData geomdata = geom.data();
#endif

  const auto dx = geom.CellSizeArray();

  ReduceOps<ReduceOpMin> reduce_op;
  ReduceData<Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  const MultiFab& stateMF = get_new_data(State_Type);

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(stateMF, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const Box& box = mfi.tilebox();

    auto u = stateMF.array(mfi);

    reduce_op.eval(box, reduce_data,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
    {

      Real rhoInv = 1.0_rt / u(i,j,k,URHO);

      eos_rep_t eos_state;
      eos_state.rho = u(i,j,k,URHO);
      eos_state.T = u(i,j,k,UTEMP);
      eos_state.e = u(i,j,k,UEINT) * rhoInv;
      for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = u(i,j,k,UFS+n) * rhoInv;
      }
#if NAUX_NET > 0
      for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = u(i,j,k,UFX+n) * rhoInv;
      }
#endif

      eos(eos_input_re, eos_state);

      // Compute velocity and then calculate CFL timestep.

      Real ux = u(i,j,k,UMX) * rhoInv;
      Real uy = u(i,j,k,UMY) * rhoInv;
      Real uz = u(i,j,k,UMZ) * rhoInv;

#ifdef ROTATION
      if (castro::do_rotation == 1 && castro::state_in_rotating_frame != 1) {
        GpuArray<Real, 3> vel;
        vel[0] = ux;
        vel[1] = uy;
        vel[2] = uz;

        inertial_to_rotational_velocity(i, j, k, geomdata, time, vel);

        ux = vel[0];
        uy = vel[1];
        uz = vel[2];
      }
#endif

      Real c = eos_state.cs;

      Real dt1 = dx[0]/(c + std::abs(ux));

      Real dt2;
#if AMREX_SPACEDIM >= 2
      dt2 = dx[1]/(c + std::abs(uy));
#else
      dt2 = dt1;
#endif

      Real dt3;
#if AMREX_SPACEDIM == 3
      dt3 = dx[2]/(c + std::abs(uz));
#else
      dt3 = dt1;
#endif

      // The CTU method has a less restrictive timestep than MOL-based
      // schemes (including the true SDC).  Since the simplified SDC
      // solver is based on CTU, we can use its timestep.
      if (castro::time_integration_method == 0 || castro::time_integration_method == 3) {
        return {amrex::min(dt1, dt2, dt3)};

      } else {
        // method of lines-style constraint is tougher
        Real dt_tmp = 1.0_rt/dt1;
#if AMREX_SPACEIM >= 2
        dt_tmp += 1.0_rt/dt2;
#endif
#if AMREX_SPACEDIM == 3
        dt_tmp += 1.0_rt/dt3;
#endif

        return 1.0_rt/dt_tmp;
      }

    });

  }

  ReduceTuple hv = reduce_data.value();
  Real estdt_hydro = amrex::get<0>(hv);

  return estdt_hydro;

}

#ifdef MHD
Real
Castro::estdt_mhd()
{

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
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
    {

      Real rhoInv = 1.0_rt / u_arr(i,j,k,URHO);
      Real bcx = 0.5_rt * (bx_arr(i+1,j,k) + bx_arr(i,j,k));
      Real bcy = 0.5_rt * (by_arr(i,j+1,k) + by_arr(i,j,k));
      Real bcz = 0.5_rt * (bz_arr(i,j,k+1) + bz_arr(i,j,k));

      Real ux = u_arr(i,j,k,UMX) * rhoInv;
      Real uy = u_arr(i,j,k,UMY) * rhoInv;
      Real uz = u_arr(i,j,k,UMZ) * rhoInv;

      eos_rep_t eos_state;
      eos_state.rho = u_arr(i,j,k,URHO);
      eos_state.e = u_arr(i,j,k,UEINT) * rhoInv;
      eos_state.T = u_arr(i,j,k,UTEMP);
      for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = u_arr(i,j,k,UFS+n) * rhoInv;
      }
#if NAUX_NET > 0
      for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = u_arr(i,j,k,UFX+n) * rhoInv;
      }
#endif

      eos(eos_input_re, eos_state);

      Real e  = u_arr(i,j,k,UEINT) * rhoInv;

      Real as = eos_state.gam1 * eos_state.p * rhoInv;
      Real ca = (bcx*bcx + bcy*bcy + bcz*bcz) * rhoInv;

      Real cx;
      Real cy;
      Real cz;

      if (e > 0_rt) {
        Real cad = bcx*bcx * rhoInv;
        eos_soundspeed_mhd(cx, as, ca, cad);

        cad = bcy*bcy * rhoInv;
        eos_soundspeed_mhd(cy, as, ca, cad);

        cad = bcz*bcz * rhoInv;
        eos_soundspeed_mhd(cz, as, ca, cad);

      } else {
        cx = 0.0_rt;
        cy = 0.0_rt;
        cz = 0.0_rt;
      }

      Real dt1 = dx[0]/(cx + std::abs(ux));

      Real dt2;
#if AMREX_SPACEDIM >= 2
      dt2 = dx[1]/(cy + std::abs(uy));
#else
      dt2 = dt1;
#endif

      Real dt3;
#if AMREX_SPACEDIM == 3
      dt3 = dx[2]/(cz + std::abs(uz));
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

Real
Castro::estdt_temp_diffusion(void)
{

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
                   [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
                   {

                     if (ustate(i,j,k,URHO) > ldiffuse_cutoff_density) {

                       Real rho_inv = 1.0_rt/ustate(i,j,k,URHO);

                       // we need cv
                       eos_t eos_state;

                       eos_state.rho = ustate(i,j,k,URHO);
                       eos_state.T = ustate(i,j,k,UTEMP);
                       eos_state.e = ustate(i,j,k,UEINT) * rho_inv;
                       for (int n = 0; n < NumSpec; n++) {
                         eos_state.xn[n]  = ustate(i,j,k,UFS+n) * rho_inv;
                       }
#if NAUX_NET > 0
                       for (int n = 0; n < NumAux; n++) {
                         eos_state.aux[n] = ustate(i,j,k,UFX+n) * rho_inv;
                       }
#endif

                       eos(eos_input_re, eos_state);

                       // we also need the conductivity
                       conductivity(eos_state);

                       // maybe we should check (and take action) on negative cv here?
                       Real D = eos_state.conductivity * rho_inv / eos_state.cv;

                       Real dt1 = 0.5_rt * dx[0]*dx[0] / D;

                       Real dt2;
#if AMREX_SPACEDIM >= 2
                       dt2 = 0.5_rt * dx[1]*dx[1] / D;
#else
                       dt2 = dt1;
#endif

                       Real dt3;
#if AMREX_SPACEDIM >= 3
                       dt3 = 0.5_rt * dx[2]*dx[2] / D;
#else
                       dt3 = dt1;
#endif

                       return {amrex::min(dt1, dt2, dt3)};

                     } else {
                       return lmax_dt/lcfl;
      }
    });
  }

  ReduceTuple hv = reduce_data.value();
  Real estdt_diff = amrex::get<0>(hv);

  return estdt_diff;
}
#endif

#ifdef REACTIONS
Real
Castro::estdt_burning()
{

    if (castro::dtnuc_e > 1.e199_rt && castro::dtnuc_X > 1.e199_rt) return 1.e200_rt;

    ReduceOps<ReduceOpMin> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    MultiFab& S_new = get_new_data(State_Type);

    const MultiFab& R_new = get_new_data(Reactions_Type);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();

        const auto S = S_new[mfi].array();
        const auto R = R_new[mfi].array();

        const auto dx = geom.CellSizeArray();

        // Set a floor on the minimum size of a derivative. This floor
        // is small enough such that it will result in no timestep limiting.

        const Real derivative_floor = 1.e-50_rt;

        // We want to limit the timestep so that it is not larger than
        // dtnuc_e * (e / (de/dt)).  If the timestep factor dtnuc is
        // equal to 1, this says that we don't want the
        // internal energy to change by any more than its current
        // magnitude in the next timestep.
        //
        // If dtnuc is less than one, it controls the fraction we will
        // allow the internal energy to change in this timestep due to
        //  nuclear burning, provided that our instantaneous estimate
        // of the energy release is representative of the full timestep.
        //
        // We also do the same thing for the species, using a timestep
        // limiter dtnuc_X * (X_k / (dX_k/dt)). To prevent changes
        // due to trace isotopes that we probably are not interested in,
        // only apply the limiter to species with an abundance greater
        // than a user-specified threshold.
        //
        // To estimate de/dt and dX/dt, we are going to call the RHS of the
        // burner given the current state data. We need to do an EOS
        // call before we do the RHS call so that we have accurate
        // values for the thermodynamic data like abar, zbar, etc.
        // But we will call in (rho, T) mode, which is inexpensive.

        reduce_op.eval(box, reduce_data,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            Real rhoInv = 1.0_rt / S(i,j,k,URHO);

            burn_t state;

            state.rho = S(i,j,k,URHO);
            state.T   = S(i,j,k,UTEMP);
            state.e   = S(i,j,k,UEINT) * rhoInv;
            for (int n = 0; n < NumSpec; ++n) {
                state.xn[n] = S(i,j,k,UFS+n) * rhoInv;
            }
#if NAUX_NET > 0
            for (int n = 0; n < NumAux; ++n) {
                state.aux[n] = S(i,j,k,UFX+n) * rhoInv;
            }
#endif

            if (state.T < castro::react_T_min || state.T > castro::react_T_max ||
                state.rho < castro::react_rho_min || state.rho > castro::react_rho_max) {
                return {1.e200_rt};
            }

            Real e    = state.e;
            Real T    = amrex::max(state.T, castro::small_temp);
            Real X[NumSpec];
            for (int n = 0; n < NumSpec; ++n) {
                X[n] = amrex::max(state.xn[n], small_x);
            }

            eos(eos_input_rt, state);

#ifdef STRANG
            state.self_heat = true;
#endif
            Array1D<Real, 1, neqs> ydot;
            actual_rhs(state, ydot);

            Real dedt = ydot(net_ienuc);
            Real dXdt[NumSpec];
            for (int n = 0; n < NumSpec; ++n) {
                dXdt[n] = ydot(n+1) * aion[n];
            }

            // Apply a floor to the derivatives. This ensures that we don't
            // divide by zero; it also gives us a quick method to disable
            // the timestep limiting, because the floor is small enough
            // that the implied timestep will be very large, and thus
            // ignored compared to other limiters.

            //dedt = amrex::max(std::abs(dedt), derivative_floor);
            dedt = amrex::max(std::abs(R(i,j,k,0) / S(i,j,k,URHO)), derivative_floor);

            for (int n = 0; n < NumSpec; ++n) {
                if (X[n] >= castro::dtnuc_X_threshold) {
                    dXdt[n] = amrex::max(std::abs(dXdt[n]), derivative_floor);
                } else {
                    dXdt[n] = derivative_floor;
                }
            }

            Real dt_tmp = 1.e200_rt;

#ifdef NSE
            // we need to use the eos_state interface here because for
            // SDC, if we come in with a burn_t, it expects to
            // evaluate the NSE criterion based on the conserved state.

            eos_t eos_state;
            burn_to_eos(state, eos_state);

            if (!in_nse(eos_state)) {
#endif
                dt_tmp = dtnuc_e * e / dedt;
#ifdef NSE
            }
#endif
            for (int n = 0; n < NumSpec; ++n) {
                dt_tmp = amrex::min(dt_tmp, dtnuc_X * (X[n] / dXdt[n]));
            }

            return {dt_tmp};
        });

    }

    ReduceTuple hv = reduce_data.value();
    Real estdt = amrex::get<0>(hv);

    return estdt;
}
#endif

#ifdef RADIATION
Real
Castro::estdt_rad ()
{
    auto dx = geom.CellSizeArray();

    const MultiFab& stateMF = get_new_data(State_Type);
    const MultiFab& radMF = get_new_data(Rad_Type);

    // Compute radiation + hydro limited timestep.

    ReduceOps<ReduceOpMin> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(stateMF, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& tbox = mfi.tilebox();
        const Box& vbox = mfi.validbox();

        FArrayBox gPr;
        gPr.resize(tbox);
        radiation->estimate_gamrPr(stateMF[mfi], radMF[mfi], gPr, dx.data(), vbox);

        auto u = stateMF[mfi].array();
        auto gPr_arr = gPr.array();

        // Note that we synchronize in the call to estimate_gamrPr. Ideally
        // we would merge these into one loop later (probably by making that
        // call occur on a zone-by-zone basis) so that we can simplify.

        reduce_op.eval(tbox, reduce_data,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            Real rhoInv = 1.0_rt / u(i,j,k,URHO);

            eos_t eos_state;
            eos_state.rho = u(i,j,k,URHO);
            eos_state.T   = u(i,j,k,UTEMP);
            eos_state.e   = u(i,j,k,UEINT) * rhoInv;
            for (int n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = u(i,j,k,UFS+n) * rhoInv;
            }
#if NAUX_NET > 0
            for (int n = 0; n < NumAux; ++n) {
                eos_state.aux[n] = u(i,j,k,UFX+n) * rhoInv;
            }
#endif

            eos(eos_input_re, eos_state);

            Real c = eos_state.cs;
            c = std::sqrt(c * c + gPr_arr(i,j,k) * rhoInv);

            Real ux = u(i,j,k,UMX) * rhoInv;
            Real uy = u(i,j,k,UMY) * rhoInv;
            Real uz = u(i,j,k,UMZ) * rhoInv;

            Real dt1 = dx[0] / (c + std::abs(ux));
#if AMREX_SPACEDIM >= 2
            Real dt2 = dx[1] / (c + std::abs(uy));
#else
            Real dt2 = std::numeric_limits<Real>::max();
#endif
#if AMREX_SPACEDIM == 3
            Real dt3 = dx[2] / (c + std::abs(uz));
#else
            Real dt3 = std::numeric_limits<Real>::max();
#endif

            Real dt_min = amrex::min(dt1, dt2, dt3);

            return {dt_min};
        });

        Gpu::synchronize();
    }

    ReduceTuple hv = reduce_data.value();
    Real estdt = amrex::get<0>(hv);
    estdt = std::min(estdt, max_dt / cfl);

    return estdt;
}
#endif
