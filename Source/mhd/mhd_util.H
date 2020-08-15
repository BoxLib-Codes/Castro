#ifndef CASTRO_MHD_SOUND_SPEED_H
#define CASTRO_MHD_SOUND_SPEED_H

#include "Castro_util.H"
#include "mhd_sizes.H"

using namespace amrex;

// sound speed calc for ideal MHD
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
eos_soundspeed_mhd(Real& c, Real as, Real ca, Real bd) {

  // inputs:
  // as = P_g * gam1/rho
  // ca = B^2/rho
  // bd = B_direction^2/rho

  // Fast Magneto-Sonic Wave
  c = 0.5_rt * ((as + ca) + std::sqrt((as + ca)*(as + ca) - 4.0_rt*as*bd));
  c = std::sqrt(c);
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
electric(Real* q_zone, Real& E_zone, const int comp) {

  // this takes the cell-center primitive state, q_zone, and computes the cell-center
  // electric field, E_zone, using Faraday's law:
  // E = -v X B

  if (comp == 0) {
    E_zone = -q_zone[QV] * q_zone[QMAGZ] + q_zone[QW] * q_zone[QMAGY];

  } else if (comp == 1) {
    E_zone = -q_zone[QW] * q_zone[QMAGX] + q_zone[QU] * q_zone[QMAGZ];

  } else {
    E_zone = -q_zone[QU] * q_zone[QMAGY] + q_zone[QV] * q_zone[QMAGX];
  }

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
ConsToPrim(Real* q_zone, Real* U_zone) {

  // calculate the primitive variables from the conserved
  // U has NUM_STATE+3 components
  // q has NQ components


  q_zone[QRHO] = U_zone[URHO];
  q_zone[QU] = U_zone[UMX] / q_zone[QRHO];
  q_zone[QV] = U_zone[UMY] / q_zone[QRHO];
  q_zone[QW] = U_zone[UMZ] / q_zone[QRHO];
  q_zone[QREINT] = U_zone[UEDEN] -
    0.5_rt * q_zone[QRHO] * (q_zone[QU] * q_zone[QU] +
                             q_zone[QV] * q_zone[QV] +
                             q_zone[QW] * q_zone[QW]) -
    0.5_rt * (U_zone[UMAGX] * U_zone[UMAGX] +
              U_zone[UMAGY] * U_zone[UMAGY] +
              U_zone[UMAGZ] * U_zone[UMAGZ]);

  // species

  for (int n = 0; n < NumSpec; n++) {
    q_zone[QFS+n] = U_zone[UFS+n] / U_zone[URHO];
  }

  eos_t eos_state;

  eos_state.rho = q_zone[QRHO];
  eos_state.e = q_zone[QREINT] / eos_state.rho;
  for (int n = 0; n < NumSpec; n++) {
    eos_state.xn[n] = q_zone[QFS+n];
  }
#if NAUX_NET > 0
  for (int n = 0; n < NumAux; n++) {
    eos_state.aux[n] = q_zone[QFX+n];
  }
#endif
  eos_state.T = 100.0_rt; // initial guess

  eos(eos_input_re, eos_state);

  q_zone[QTEMP] = eos_state.T;
  q_zone[QPRES] = eos_state.p;

  q_zone[QMAGX] = U_zone[UMAGX];
  q_zone[QMAGY] = U_zone[UMAGY];
  q_zone[QMAGZ] = U_zone[UMAGZ];

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
PToC(Array1D<Real, 0, NQ-1>& q,
     Array1D<Real, 0, NUM_STATE+2>& U,
     Real& gam1) {

  // first the conserved state
  U(URHO) = q(QRHO);
  U(UMX) = q(QRHO) * q(QU);
  U(UMY) = q(QRHO) * q(QV);
  U(UMZ) = q(QRHO) * q(QW);

  eos_t eos_state;

  eos_state.rho = q(QRHO);
  eos_state.p = q(QPRES);
  for (int n = 0; n < NumSpec; n++) {
    eos_state.xn[n] = q(QFS+n);
  }
#if NAUX_NET > 0
  for (int n = 0; n < NumAux; n++) {
    eos_state.aux[n] = q(QFX+n);
  }
#endif
  eos_state.T = 100.0;  // dummy initial guess

  eos(eos_input_rp, eos_state);

  Real B2 = q(QMAGX) * q(QMAGX) + q(QMAGY) * q(QMAGY) + q(QMAGZ) * q(QMAGZ);

  U(UEINT) = eos_state.rho * eos_state.e;
  U(UEDEN) = U(UEINT) +
    0.5_rt * q(QRHO) * (q(QU) * q(QU) + q(QV) * q(QV) + q(QW) * q(QW)) + 0.5_rt * B2;

  U(UMAGX) = q(QMAGX);
  U(UMAGY) = q(QMAGY);
  U(UMAGZ) = q(QMAGZ);

  for (int n = 0; n < NumSpec; n++) {
    U(UFS+n) = U(URHO) * q(QFS+n);
  }

  U(UTEMP) = eos_state.T;

  gam1 = eos_state.gam1;

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int
epsilon_ijk(const int i, const int j, const int k) {

  int val;

  if (i == j || j == k || i == k) {
    val = 0;
  } else if ((i == 0 && j == 1 && k == 2) ||
             (i == 1 && j == 2 && k == 0) ||
             (i == 2 && j == 0 && k == 1)) {
    val = 1;
  } else {
    val = -1;
  }

  return val;
}

#endif

