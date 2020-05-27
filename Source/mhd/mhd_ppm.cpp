#include "Castro.H"
#include "Castro_F.H"

#include "ppm.H"

using namespace amrex;

#include "mhd_eigen.H"
#include "mhd_slope.H"

void
Castro::ppm_mhd(const Box& bx,
                const int idir,
                Array4<Real const> const& q_arr,
                Array4<Real const> const& qaux,
                Array4<Real const> const& flatn,
                Array4<Real const> const& Bx,
                Array4<Real const> const& By,
                Array4<Real const> const& Bz,
                Array4<Real> const& qleft,
                Array4<Real> const& qright,
                Array4<Real const> const& srcQ,
                const Real dt) {

  // these loops are over cell-centers and for each cell-center, we find the left and
  // right interface states

  const auto dx = geom.CellSizeArray();

  Real dtdx = dt/dx[idir];

  // these are the characteristic variables for this direction
  int cvars[NEIGN];

  if (idir == 0) {

    // component (Bx) is omitted

    cvars[IEIGN_RHO] = QRHO;
    cvars[IEIGN_U] = QU;
    cvars[IEIGN_V] = QV;
    cvars[IEIGN_W] = QW;
    cvars[IEIGN_P] = QPRES;
    cvars[IEIGN_BT] = QMAGY;
    cvars[IEIGN_BTT] = QMAGZ;

  } else if (idir == 1) {

    // component (By) is omitted

    cvars[IEIGN_RHO] = QRHO;
    cvars[IEIGN_U] = QU;
    cvars[IEIGN_V] = QV;
    cvars[IEIGN_W] = QW;
    cvars[IEIGN_P] = QPRES;
    cvars[IEIGN_BT] = QMAGX;
    cvars[IEIGN_BTT] = QMAGZ;

  } else {

    // component (Bz) is omitted

    cvars[IEIGN_RHO] = QRHO;
    cvars[IEIGN_U] = QU;
    cvars[IEIGN_V] = QV;
    cvars[IEIGN_W] = QW;
    cvars[IEIGN_P] = QPRES;
    cvars[IEIGN_BT] = QMAGX;
    cvars[IEIGN_BTT] = QMAGY;

  }


  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
  {

    // compute the eigenvectors and eigenvalues for this coordinate direction

    Array1D<Real, 0, NQ-1> q_zone;
    for (int n = 0; n < NQ; n++) {
      q_zone(n) = q_arr(i,j,k,n);
    }

    Real as = qaux(i,j,k,QC);

    Array1D<Real, 0, NEIGN-1> lam;

    evals(lam, as, q_zone, idir);

    Array2D<Real, 0, NEIGN-1, 0, NEIGN-1> leig;
    Array2D<Real, 0, NEIGN-1, 0, NEIGN-1> reig;

    if (idir == 0) {
      evecx(leig, reig, as, q_zone);

    } else if (idir == 1) {
      evecy(leig, reig, as, q_zone);

    } else {
      evecz(leig, reig, as, q_zone);
    }


    // do the parabolic reconstruction and compute the integrals under
    // the characteristic waves
    Real s[5];
    Real flat = flatn(i,j,k);
    Real sm;
    Real sp;

    Real Ip[NEIGN][NEIGN];  // first component is var, second is wave
    Real Im[NEIGN][NEIGN];


    for (int n = 0; n < NEIGN; n++) {

      int v = cvars[n];

      if (idir == 0) {
        s[im2] = q_arr(i-2,j,k,v);
        s[im1] = q_arr(i-1,j,k,v);
        s[i0]  = q_arr(i,j,k,v);
        s[ip1] = q_arr(i+1,j,k,v);
        s[ip2] = q_arr(i+2,j,k,v);

      } else if (idir == 1) {
        s[im2] = q_arr(i,j-2,k,v);
        s[im1] = q_arr(i,j-1,k,v);
        s[i0]  = q_arr(i,j,k,v);
        s[ip1] = q_arr(i,j+1,k,v);
        s[ip2] = q_arr(i,j+2,k,v);

      } else {
        s[im2] = q_arr(i,j,k-2,v);
        s[im1] = q_arr(i,j,k-1,v);
        s[i0]  = q_arr(i,j,k,v);
        s[ip1] = q_arr(i,j,k+1,v);
        s[ip2] = q_arr(i,j,k+2,v);

      }

      ppm_reconstruct(s, flat, sm, sp);

      for (int ii = 0; ii < NEIGN; ii++) {
        Real Ipt;
        Real Imt;
        ppm_int_profile_single(sm, sp, s[i0], lam(ii), dtdx, Ipt, Imt);

        Ip[n][ii] = Ipt;
        Im[n][ii] = Imt;
      }
    }

    // MHD Source Terms -- from the Miniati paper, Eq. 32 and 33
    Real smhd[NEIGN];

    smhd[IEIGN_RHO] = 0.0_rt;
    smhd[IEIGN_U] = q_zone(QMAGX) / q_zone(QRHO);
    smhd[IEIGN_V] = q_zone(QMAGY) / q_zone(QRHO);
    smhd[IEIGN_W] = q_zone(QMAGZ) / q_zone(QRHO);
    smhd[IEIGN_P] = q_zone(QMAGX) * q_zone(QU) +
                    q_zone(QMAGY) * q_zone(QV) +
                    q_zone(QMAGZ) * q_zone(QW);

    if (idir == 0) {
      smhd[IEIGN_BT] = q_zone(QV);
      smhd[IEIGN_BTT] = q_zone(QW);

      // cross-talk of normal magnetic field direction
      for (int n = 0; n < NEIGN; n++) {
        smhd[n] = smhd[n] * (Bx(i+1,j,k) - Bx(i,j,k)) / dx[idir];
      }

    } else if (idir == 1) {
      smhd[IEIGN_BT] = q_zone(QU);
      smhd[IEIGN_BTT] = q_zone(QW);

      // cross-talk of normal magnetic field direction
      for (int n = 0; n < NEIGN; n++) {
        smhd[n] = smhd[n] * (By(i,j+1,k) - By(i,j,k)) / dx[idir];
      }

    } else {
      smhd[IEIGN_BT] = q_zone(QU);
      smhd[IEIGN_BTT] = q_zone(QV);

      // cross-talk of normal magnetic field direction
      for (int n = 0; n < NEIGN; n++) {
        smhd[n] = smhd[n] * (Bz(i,j,k+1) - Bz(i,j,k)) / dx[idir];
      }
    }

    // Perform the characteristic projection.  Since we are using
    // Using HLLD, we sum over all eigenvalues -- see the discussion
    // after Eq. 31

    // right state at i-1/2

    // Im is the integral from the left edge, so we take as the
    // reference state the fastest wave moving to the left

    Real summ_m[NEIGN] = {0.0_rt};

    // loop over the waves
    for (int ii = 0; ii < NEIGN; ii++) {

      Real LdQ = 0.0_rt;

      // loop over variables in Im[n][ii]
      for (int n = 0; n < NEIGN; n++) {
        LdQ += leig(ii,n) * (Im[n][0] - Im[n][ii]);
      }

      // add the contribution of this wave to each variable
      for (int n = 0; n < NEIGN; n++) {
        summ_m[n] += LdQ * reig(n,ii);
      }
    }

    qright(i,j,k,QRHO) = amrex::max(small_dens,
                                    Im[IEIGN_RHO][0] - summ_m[IEIGN_RHO] + 0.5_rt*dt*smhd[IEIGN_RHO]);
    qright(i,j,k,QU) = Im[IEIGN_U][0] - summ_m[IEIGN_U] + 0.5_rt*dt*smhd[IEIGN_U];
    qright(i,j,k,QV) = Im[IEIGN_V][0] - summ_m[IEIGN_V] + 0.5_rt*dt*smhd[IEIGN_V];
    qright(i,j,k,QW) = Im[IEIGN_W][0] - summ_m[IEIGN_W] + 0.5_rt*dt*smhd[IEIGN_W];
    qright(i,j,k,QPRES) = amrex::max(small_pres,
                                     Im[IEIGN_P][0] - summ_m[IEIGN_P] + 0.5_rt*dt*smhd[IEIGN_P]);

    if (idir == 0) {
      qright(i,j,k,QMAGX) = Bx(i,j,k); // Bx stuff
      qright(i,j,k,QMAGY) = Im[IEIGN_BT][0] -  summ_m[IEIGN_BT] + 0.5_rt*dt*smhd[IEIGN_BT];
      qright(i,j,k,QMAGZ) = Im[IEIGN_BTT][0] - summ_m[IEIGN_BTT] + 0.5_rt*dt*smhd[IEIGN_BTT];

    } else if (idir == 1) {
      qright(i,j,k,QMAGX) = Im[IEIGN_BT][0] - summ_m[IEIGN_BT] + 0.5_rt*dt*smhd[IEIGN_BT];
      qright(i,j,k,QMAGY) = By(i,j,k); // By stuff
      qright(i,j,k,QMAGZ) = Im[IEIGN_BTT][0] - summ_m[IEIGN_BTT] + 0.5_rt*dt*smhd[IEIGN_BTT];

    } else {
      qright(i,j,k,QMAGX) = Im[IEIGN_BT][0] - summ_m[IEIGN_BT] + 0.5_rt*dt*smhd[IEIGN_BT];
      qright(i,j,k,QMAGY) = Im[IEIGN_BTT][0] - summ_m[IEIGN_BTT] + 0.5_rt*dt*smhd[IEIGN_BTT];
      qright(i,j,k,QMAGZ) = Bz(i,j,k); // Bz stuff
    }


    // left state at i+1/2

    // Ip is the integral from the right edge, so we take as the
    // reference state the fastest wave moving to the right

    Real summ_p[NEIGN] = {0.0_rt};

    // loop over the waves
    for (int ii = 0; ii < NEIGN; ii++) {

      Real LdQ = 0.0_rt;

      // loop over variables in Im[n][ii]
      for (int n = 0; n < NEIGN; n++) {
        LdQ += leig(ii,n) * (Ip[n][NEIGN-1] - Ip[n][ii]);
      }

      // add the contribution of this wave to each variable
      for (int n = 0; n < NEIGN; n++) {
        summ_p[n] += LdQ * reig(n,ii);
      }
    }

    if (idir == 0) {
      qleft(i+1,j,k,QRHO) = amrex::max(small_dens,
                                       Ip[IEIGN_RHO][NEIGN-1] - summ_p[IEIGN_RHO] + 0.5_rt*dt*smhd[IEIGN_RHO]);
      qleft(i+1,j,k,QU) = Ip[IEIGN_U][NEIGN-1] - summ_p[IEIGN_U] + 0.5_rt*dt*smhd[IEIGN_U];
      qleft(i+1,j,k,QV) = Ip[IEIGN_V][NEIGN-1] - summ_p[IEIGN_V] + 0.5_rt*dt*smhd[IEIGN_V];
      qleft(i+1,j,k,QW) = Ip[IEIGN_W][NEIGN-1] - summ_p[IEIGN_W] + 0.5_rt*dt*smhd[IEIGN_W];
      qleft(i+1,j,k,QPRES) = amrex::max(small_pres,
                                        Ip[IEIGN_P][NEIGN-1] - 0.5_rt*summ_p[IEIGN_P] + 0.5_rt*dt*smhd[IEIGN_P]);

      qleft(i+1,j,k,QMAGX) = Bx(i+1,j,k); // Bx stuff
      qleft(i+1,j,k,QMAGY) = Ip[IEIGN_BT][NEIGN-1] - summ_p[IEIGN_BT] + 0.5_rt*dt*smhd[IEIGN_BT];
      qleft(i+1,j,k,QMAGZ) = Ip[IEIGN_BTT][NEIGN-1] - summ_p[IEIGN_BTT] + 0.5_rt*dt*smhd[IEIGN_BTT];

    } else if (idir == 1) {
      qleft(i,j+1,k,QRHO) = amrex::max(small_dens,
                                       Ip[IEIGN_RHO][NEIGN-1] - summ_p[IEIGN_RHO] + 0.5_rt*dt*smhd[IEIGN_RHO]);
      qleft(i,j+1,k,QU) = Ip[IEIGN_U][NEIGN-1] - summ_p[IEIGN_U] + 0.5_rt*dt*smhd[IEIGN_U];
      qleft(i,j+1,k,QV) = Ip[IEIGN_V][NEIGN-1] - summ_p[IEIGN_V] + 0.5_rt*dt*smhd[IEIGN_V];
      qleft(i,j+1,k,QW) = Ip[IEIGN_W][NEIGN-1] - summ_p[IEIGN_W] + 0.5_rt*dt*smhd[IEIGN_W];
      qleft(i,j+1,k,QPRES) = amrex::max(small_pres,
                                        Ip[IEIGN_P][NEIGN-1] - summ_p[IEIGN_P] + 0.5_rt*dt*smhd[IEIGN_P]);

      qleft(i,j+1,k,QMAGX) = Ip[IEIGN_BT][NEIGN-1] - summ_p[IEIGN_BT] + 0.5_rt*dt*smhd[IEIGN_BT];
      qleft(i,j+1,k,QMAGY) = By(i,j+1,k); // By stuff
      qleft(i,j+1,k,QMAGZ) = Ip[IEIGN_BTT][NEIGN-1] - summ_p[IEIGN_BTT] + 0.5_rt*dt*smhd[IEIGN_BTT];

    } else {
      qleft(i,j,k+1,QRHO) = amrex::max(small_dens,
                                       Ip[IEIGN_RHO][NEIGN-1] - summ_p[IEIGN_RHO] + 0.5_rt*dt*smhd[IEIGN_RHO]);
      qleft(i,j,k+1,QU) = Ip[IEIGN_U][NEIGN-1] - summ_p[IEIGN_U] + 0.5_rt*dt*smhd[IEIGN_U];
      qleft(i,j,k+1,QV) = Ip[IEIGN_V][NEIGN-1] - summ_p[IEIGN_V] + 0.5_rt*dt*smhd[IEIGN_V];
      qleft(i,j,k+1,QW) = Ip[IEIGN_W][NEIGN-1] - summ_p[IEIGN_W] + 0.5_rt*dt*smhd[IEIGN_W];
      qleft(i,j,k+1,QPRES) = amrex::max(small_pres,
                                        Ip[IEIGN_P][NEIGN-1] - summ_p[IEIGN_P] + 0.5_rt*dt*smhd[IEIGN_P]);

      qleft(i,j,k+1,QMAGX) = Ip[IEIGN_BT][NEIGN-1] - summ_p[IEIGN_BT] + 0.5_rt*dt*smhd[IEIGN_BT];
      qleft(i,j,k+1,QMAGY) = Ip[IEIGN_BTT][NEIGN-1] - summ_p[IEIGN_BTT] + 0.5_rt*dt*smhd[IEIGN_BTT];
      qleft(i,j,k+1,QMAGZ) = Bz(i,j,k+1); // Bz stuff
    }


    // species
    for (int n = 0; n < NumSpec; n++) {

      int v = QFS+n;

      if (idir == 0) {
        s[im2] = q_arr(i-2,j,k,v);
        s[im1] = q_arr(i-1,j,k,v);
        s[i0]  = q_arr(i,j,k,v);
        s[ip1] = q_arr(i+1,j,k,v);
        s[ip2] = q_arr(i+2,j,k,v);

      } else if (idir == 1) {
        s[im2] = q_arr(i,j-2,k,v);
        s[im1] = q_arr(i,j-1,k,v);
        s[i0]  = q_arr(i,j,k,v);
        s[ip1] = q_arr(i,j+1,k,v);
        s[ip2] = q_arr(i,j+2,k,v);

      } else {
        s[im2] = q_arr(i,j,k-2,v);
        s[im1] = q_arr(i,j,k-1,v);
        s[i0]  = q_arr(i,j,k,v);
        s[ip1] = q_arr(i,j,k+1,v);
        s[ip2] = q_arr(i,j,k+2,v);

      }

      Real Ips;
      Real Ims;

      Real un;

      if (idir == 0) {
        un = q_arr(i,j,k,QU);

      } else if (idir == 1) {
        un = q_arr(i,j,k,QV);

      } else {
        un = q_arr(i,j,k,QW);
      }

      ppm_reconstruct(s, flat, sm, sp);
      ppm_int_profile_single(sm, sp, s[i0], un, dtdx, Ips, Ims);

      if (idir == 0) {
        qleft(i+1,j,k,QFS+n) = Ips;
      } else if (idir == 1) {
        qleft(i,j+1,k,QFS+n) = Ips;
      } else {
        qleft(i,j,k+1,QFS+n) = Ips;
      }
      qright(i,j,k,QFS+n) = Ims;
    }

    // rho e
    eos_t eos_state;

    if (idir == 0) {
      eos_state.rho = qleft(i+1,j,k,QRHO);
      eos_state.p = qleft(i+1,j,k,QPRES);
      eos_state.T = q_arr(i,j,k,QTEMP); // some initial guess?
      for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = qleft(i+1,j,k,QFS+n);
      }
      eos(eos_input_rp, eos_state);
      qleft(i+1,j,k,QREINT) = eos_state.e * eos_state.rho;

    } else if (idir == 1) {
      eos_state.rho = qleft(i,j+1,k,QRHO);
      eos_state.p = qleft(i,j+1,k,QPRES);
      eos_state.T = q_arr(i,j,k,QTEMP); // some initial guess?
      for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = qleft(i,j+1,k,QFS+n);
      }
      eos(eos_input_rp, eos_state);
      qleft(i,j+1,k,QREINT) = eos_state.e * eos_state.rho;

    } else {
      eos_state.rho = qleft(i,j,k+1,QRHO);
      eos_state.p = qleft(i,j,k+1,QPRES);
      eos_state.T = q_arr(i,j,k,QTEMP); // some initial guess?
      for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = qleft(i,j,k+1,QFS+n);
      }
      eos(eos_input_rp, eos_state);
      qleft(i,j,k+1,QREINT) = eos_state. e * eos_state.rho;
    }

    eos_state.rho = qright(i,j,k,QRHO);
    eos_state.p = qright(i,j,k,QPRES);
    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n] = qright(i,j,k,QFS+n);
    }

    eos(eos_input_rp, eos_state);
    qright(i,j,k,QREINT) = eos_state.e * eos_state.rho;

    // add source terms
    if (idir == 0) {
      qleft(i+1,j,k,QRHO) = amrex::max(small_dens,
                                       qleft(i+1,j,k,QRHO) + 0.5_rt*dt*srcQ(i,j,k,QRHO));
      qleft(i+1,j,k,QU) = qleft(i+1,j,k,QU) + 0.5_rt*dt*srcQ(i,j,k,QU);
      qleft(i+1,j,k,QV) = qleft(i+1,j,k,QV) + 0.5_rt*dt*srcQ(i,j,k,QV);
      qleft(i+1,j,k,QW) = qleft(i+1,j,k,QW) + 0.5_rt*dt*srcQ(i,j,k,QW);
      qleft(i+1,j,k,QPRES) = qleft(i+1,j,k,QPRES) + 0.5_rt*dt*srcQ(i,j,k,QPRES);
      qleft(i+1,j,k,QREINT) = qleft(i+1,j,k,QREINT) + 0.5_rt*dt*srcQ(i,j,k,QREINT);

    } else if (idir == 1) {
      qleft(i,j+1,k,QRHO) = amrex::max(small_dens,
                                       qleft(i,j+1,k,QRHO) + 0.5_rt*dt*srcQ(i,j,k,QRHO));
      qleft(i,j+1,k,QU) = qleft(i,j+1,k,QU) + 0.5_rt*dt*srcQ(i,j,k,QU);
      qleft(i,j+1,k,QV) = qleft(i,j+1,k,QV) + 0.5_rt*dt*srcQ(i,j,k,QV);
      qleft(i,j+1,k,QW) = qleft(i,j+1,k,QW) + 0.5_rt*dt*srcQ(i,j,k,QW);
      qleft(i,j+1,k,QPRES) = qleft(i,j+1,k,QPRES) + 0.5_rt*dt*srcQ(i,j,k,QPRES);
      qleft(i,j+1,k,QREINT) = qleft(i,j+1,k,QREINT) + 0.5_rt*dt*srcQ(i,j,k,QREINT);

    } else {
      qleft(i,j,k+1,QRHO) = max(small_dens, qleft(i,j,k+1,QRHO) + 0.5_rt*dt*srcQ(i,j,k,QRHO));
      qleft(i,j,k+1,QU) = qleft(i,j,k+1,QU) + 0.5_rt*dt*srcQ(i,j,k,QU);
      qleft(i,j,k+1,QV) = qleft(i,j,k+1,QV) + 0.5_rt*dt*srcQ(i,j,k,QV);
      qleft(i,j,k+1,QW) = qleft(i,j,k+1,QW) + 0.5_rt*dt*srcQ(i,j,k,QW);
      qleft(i,j,k+1,QPRES) = qleft(i,j,k+1,QPRES) + 0.5_rt*dt*srcQ(i,j,k,QPRES);
      qleft(i,j,k+1,QREINT) = qleft(i,j,k+1,QREINT) + 0.5_rt*dt*srcQ(i,j,k,QREINT);
    }

    qright(i,j,k,QRHO) = max(small_dens, qright(i,j,k,QRHO) + 0.5_rt*dt*srcQ(i,j,k,QRHO));
    qright(i,j,k,QU) = qright(i,j,k,QU) + 0.5_rt*dt*srcQ(i,j,k,QU);
    qright(i,j,k,QV) = qright(i,j,k,QV) + 0.5_rt*dt*srcQ(i,j,k,QV);
    qright(i,j,k,QW) = qright(i,j,k,QW) + 0.5_rt*dt*srcQ(i,j,k,QW);
    qright(i,j,k,QPRES) = qright(i,j,k,QPRES) + 0.5_rt*dt*srcQ(i,j,k,QPRES);
    qright(i,j,k,QREINT) = qright(i,j,k,QREINT) + 0.5_rt*dt*srcQ(i,j,k,QREINT);

  });
}
