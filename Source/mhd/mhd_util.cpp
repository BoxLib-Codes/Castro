#include "Castro.H"
#include "Castro_F.H"
#include "Castro_util.H"
#include "Castro_hydro_F.H"

#include "mhd_util.H"

using namespace amrex;

void
Castro::check_for_mhd_cfl_violation(const Box& bx,
                                    const Real dt,
                                    Array4<Real const> const& q_arr,
                                    Array4<Real const> const& qaux_arr) {

  auto dx = geom.CellSizeArray();

  Real dtdx = dt / dx[0];

#if AMREX_SPACEDIM >= 2
  Real dtdy = dt / dx[1];
#else
  Real dtdy = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
  Real dtdz = dt / dx[2];
#else
  Real dtdz = 0.0_rt;
#endif

  ReduceOps<ReduceOpMax> reduce_op;
  ReduceData<Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  reduce_op.eval(bx, reduce_data,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
  {

    //sound speed for ideal mhd
    Real as2 = qaux_arr(i,j,k,QC) * qaux_arr(i,j,k,QC);
    Real rho_inv = 1.0_rt / q_arr(i,j,k,QRHO);

    Real ca = (q_arr(i,j,k,QMAGX) * q_arr(i,j,k,QMAGX) +
               q_arr(i,j,k,QMAGY) * q_arr(i,j,k,QMAGY) +
               q_arr(i,j,k,QMAGZ) * q_arr(i,j,k,QMAGZ)) * rho_inv;

    Real cx = 0.0;
    Real cad = q_arr(i,j,k,QMAGX) * q_arr(i,j,k,QMAGX) * rho_inv;
    eos_soundspeed_mhd(cx, as2, ca, cad);

    Real cy = 0.0;
    cad = q_arr(i,j,k,QMAGY) * q_arr(i,j,k,QMAGY) * rho_inv;
    eos_soundspeed_mhd(cy, as2, ca, cad);

    Real cz = 0.0;
    cad = q_arr(i,j,k,QMAGZ) * q_arr(i,j,k,QMAGZ) * rho_inv;
    eos_soundspeed_mhd(cz, as2, ca, cad);

    Real courx = (cx + std::abs(q_arr(i,j,k,QU))) * dtdx;
    Real coury = (cy + std::abs(q_arr(i,j,k,QV))) * dtdy;
    Real courz = (cz + std::abs(q_arr(i,j,k,QW))) * dtdz;

#ifndef AMREX_USE_CUDA
    if (verbose == 1) {

      if (courx > 1.0_rt) {
        std::cout << std::endl;
        std::cout << "Warning:: CFL violation in check_for_mhd_cfl_violation" << std::endl;
        std::cout << ">>> ... (u+c) * dt / dx > 1 " << courx << std::endl;
        std::cout << ">>> ... at cell i = " << i << " j = " << j << " k = " << k << std::endl;
        std::cout << ">>> ... u = " << q_arr(i,j,k,QU) << " c = " << cx << std::endl;
        std::cout << ">>> ... B = " << q_arr(i,j,k,QMAGX) << " " << q_arr(i,j,k,QMAGY) << " " << q_arr(i,j,k,QMAGZ) << std::endl;
        std::cout << ">>> ... density = " << q_arr(i,j,k,QRHO) << std::endl;
        std::cout << ">>> ... internal e = " << q_arr(i,j,k,QREINT) << std::endl;
        std::cout << ">>> ... pressure = " << q_arr(i,j,k,QPRES) << std::endl;
      }

      if (coury > 1.0_rt) {
        std::cout << std::endl;
        std::cout << "Warning:: CFL violation in check_for_mhd_cfl_violation" << std::endl;
        std::cout << ">>> ... (v+c) * dt / dx > 1 " << coury << std::endl;
        std::cout << ">>> ... at cell i = " << i << " j = " << j << " k = " << k << std::endl;
        std::cout << ">>> ... v = " << q_arr(i,j,k,QV) << " c = " << cy << std::endl;
        std::cout << ">>> ... B = " << q_arr(i,j,k,QMAGX) << " " << q_arr(i,j,k,QMAGY) << " " << q_arr(i,j,k,QMAGZ) << std::endl;
        std::cout << ">>> ... density = " << q_arr(i,j,k,QRHO) << std::endl;
        std::cout << ">>> ... internal e = " << q_arr(i,j,k,QREINT) << std::endl;
        std::cout << ">>> ... pressure = " << q_arr(i,j,k,QPRES) << std::endl;
      }

      if (courz > 1.0_rt) {
        std::cout << std::endl;
        std::cout << "Warning:: CFL violation in check_for_mhd_cfl_violation" << std::endl;
        std::cout << ">>> ... (w+c) * dt / dx > 1 " << courz << std::endl;
        std::cout << ">>> ... at cell i = " << i << " j = " << j << " k = " << k << std::endl;
        std::cout << ">>> ... w = " << q_arr(i,j,k,QW) << " c = " << cz << std::endl;
        std::cout << ">>> ... B = " << q_arr(i,j,k,QMAGX) << " " << q_arr(i,j,k,QMAGY) << " " << q_arr(i,j,k,QMAGZ) << std::endl;
        std::cout << ">>> ... density = " << q_arr(i,j,k,QRHO) << std::endl;
        std::cout << ">>> ... internal e = " << q_arr(i,j,k,QREINT) << std::endl;
        std::cout << ">>> ... pressure = " << q_arr(i,j,k,QPRES) << std::endl;

      }

    }
#endif

    return {amrex::max(courx, coury, courz)};

  });

  ReduceTuple hv = reduce_data.value();
  Real courno = amrex::get<0>(hv);

  if (courno > 1.0) {
    amrex::Print() << "WARNING -- EFFECTIVE CFL AT LEVEL " << level << " IS " << courno << std::endl << std::endl;

    cfl_violation = 1;
  }

}


void
Castro::consup_mhd(const Box& bx,
                   Array4<Real> const& update,
                   Array4<Real const> const& flux0,
                   Array4<Real const> const& flux1,
                   Array4<Real const> const& flux2) {

  // do the conservative update and store - div{F} in update.  Note,
  // in contrast to the CTU hydro case, we don't add the pdivu term to
  // (rho e) here, because we included that in the hydro source terms.

  const auto dx = geom.CellSizeArray();

  Real dxinv = 1.0_rt/dx[0];
#if AMREX_SPACEDIM >= 2
  Real dyinv = 1.0_rt/dx[1];
#endif
#if AMREX_SPACEDIM == 3
  Real dzinv = 1.0_rt/dx[2];
#endif

  amrex::ParallelFor(bx, NUM_STATE,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k, int n) noexcept
  {

    if (n == UTEMP) {
      update(i,j,k,n) = 0.0_rt;
#ifdef SHOCK_VAR
    } else if (n == USHK) {
      update(i,j,k,n) = 0.0_rt;
#endif
    } else {
      update(i,j,k,n) = (flux0(i,j,k,n) - flux0(i+1,j,k,n)) * dxinv;
#if AMREX_SPACEDIM >= 2
      update(i,j,k,n) += (flux1(i,j,k,n) - flux1(i,j+1,k,n)) * dyinv;
#endif
#if AMREX_SPACEDIM == 3
      update(i,j,k,n) += (flux2(i,j,k,n) - flux2(i,j,k+1,n)) * dzinv;
#endif
    }

  });

}


void
Castro::PrimToCons(const Box& bx,
                   Array4<Real const> const& q_arr,
                   Array4<Real> const& u_arr) {

  // calculate the conserved variables from the primitive

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
  {

    u_arr(i,j,k,URHO) = q_arr(i,j,k,QRHO);
    u_arr(i,j,k,UMX) = q_arr(i,j,k,QRHO)*q_arr(i,j,k,QU);
    u_arr(i,j,k,UMY) = q_arr(i,j,k,QRHO)*q_arr(i,j,k,QV);
    u_arr(i,j,k,UMZ) = q_arr(i,j,k,QRHO)*q_arr(i,j,k,QW);

    eos_t eos_state;
    eos_state.rho = q_arr(i,j,k,QRHO);
    eos_state.p = q_arr(i,j,k,QPRES);
    eos_state.T = 100.0_rt; // dummy initial T.
    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n] = q_arr(i,j,k,QFS+n);
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; n++) {
      eos_state.aux[n] = q_arr(i,j,k,QFX+n);
    }
#endif

    eos(eos_input_rp, eos_state);

    u_arr(i,j,k,UEDEN) = eos_state.rho * eos_state.e +
      + 0.5_rt * q_arr(i,j,k,QRHO) * (q_arr(i,j,k,QU) * q_arr(i,j,k,QU) +
                                      q_arr(i,j,k,QV) * q_arr(i,j,k,QV) +
                                      q_arr(i,j,k,QW) * q_arr(i,j,k,QW)) +
      + 0.5_rt * (q_arr(i,j,k,QMAGX) * q_arr(i,j,k,QMAGX) +
                  q_arr(i,j,k,QMAGY) * q_arr(i,j,k,QMAGY) +
                  q_arr(i,j,k,QMAGZ) * q_arr(i,j,k,QMAGZ));

    u_arr(i,j,k,UEINT) = eos_state.rho * eos_state.e;
    u_arr(i,j,k,UTEMP) = eos_state.T;

    u_arr(i,j,k,UMAGX) = q_arr(i,j,k,QMAGX);
    u_arr(i,j,k,UMAGY) = q_arr(i,j,k,QMAGY);
    u_arr(i,j,k,UMAGZ) = q_arr(i,j,k,QMAGZ);

    // species
    for (int n = 0; n < NumSpec; n++) {
      u_arr(i,j,k,UFS+n) = q_arr(i,j,k,QRHO) * q_arr(i,j,k,QFS+n);
    }

  });
}


void
Castro::prim_half(const Box& bx,
                  Array4<Real> const& q2D,
                  Array4<Real const> const& q_arr,
                  Array4<Real const> const& flxx,
                  Array4<Real const> const& flxy,
                  Array4<Real const> const& flxz,
                  const Real dt) {

  // Find the 2D corrected primitive variables and n+1/2

  auto dx = geom.CellSizeArray();

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
  {

    Array1D<Real, 0, NQ-1> q_zone;
    Array1D<Real, 0, NUM_STATE+2> U_zone;

    // we come in with q, so we need to compute this to U
    for (int n = 0; n < NQ; n++) {
      q_zone(n) = q_arr(i,j,k,n);
    }

    Real gam1;
    PToC(q_zone, U_zone, gam1);

    // now do the conservative update of U

    for (int n = 0; n < NUM_STATE+3; n++) {
      U_zone(n) -= 0.5 * dt * (flxx(i+1,j,k,n) - flxx(i,j,k,n)) / dx[0];
#if AMREX_SPACEDIM >= 2
      U_zone(n) -= 0.5 * dt * (flxy(i,j+1,k,n) - flxy(i,j,k,n)) / dx[1];
#endif
#if AMREX_SPACEDIM == 3
      U_zone(n) -= 0.5 * dt * (flxz(i,j,k+1,n) - flxz(i,j,k,n)) / dx[2];
#endif
    }

    // now convert this to q
    ConsToPrim(q_zone.arr, U_zone.arr);

    // Right below eq. 48
    for (int n = 0; n < NQ; n++) {
      q2D(i,j,k,n) = q_zone(n);
    }
  });
}
