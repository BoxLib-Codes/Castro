#include "Castro.H"
#include "Castro_F.H"
#include "Rotation.H"
#include "Castro_util.H"
#include "math.H"

void
Castro::rsrc(const Box& bx,
             Array4<Real const> const& phi,
             Array4<Real const> const& rot,
             Array4<Real const> const& uold,
             Array4<Real> const& source, 
             Array4<Real const> const& vol,
             const Real dt, const Real time) {

  GpuArray<Real, 3> center;
  ca_get_center(center.begin());

  GeometryData geomdata = geom.data();

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
  {

    Real Sr[3] = {};
    Real src[NSRC] = {};

    // Temporary array for seeing what the new state would be if the update were applied here.

    Real snew[NUM_STATE] = {};

    GpuArray<Real, 3> loc;
    position(i, j, k, geomdata, loc);

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      loc[dir] -= center[dir];
    }

    Real rho = uold(i,j,k,URHO);
    Real rhoInv = 1.0_rt / rho;

    for (int n = 0; n < NUM_STATE; n++) {
      snew[n] = uold(i,j,k,n);
    }

    Real old_ke = 0.5_rt * (snew[UMX] * snew[UMX] + snew[UMY] * snew[UMY] + snew[UMZ] * snew[UMZ]) * rhoInv;

    for (int n = 0; n < 3; n++) {
      Sr[n] = rho * rot(i,j,k,n);
    }

    src[UMX] = Sr[UMX];
    src[UMY] = Sr[UMY];
    src[UMZ] = Sr[UMZ];

    snew[UMX] += dt * src[UMX];
    snew[UMY] += dt * src[UMY];
    snew[UMZ] += dt * src[UMZ];

#ifdef HYBRID_MOMENTUM
    if (state_in_rotating_frame == 1) {
      set_hybrid_momentum_source(loc, src(UMR:UMP), Sr);

      snew[UMR] += dt * src[UMR];
      snew[UML] += dt * src[UML];
      snew[UMP] += dt * src[UMP];
    }
#endif

    // Kinetic energy source: this is v . the momentum source.
    // We don't apply in the case of the conservative energy
    // formulation.

    Real SrE;

    if (rot_source_type == 1 || rot_source_type == 2) {

      SrE = uold(i,j,k,UMX) * rhoInv * Sr[UMX] +
            uold(i,j,k,UMY) * rhoInv * Sr[UMY] +
            uold(i,j,k,UMZ) * rhoInv * Sr[UMZ];

    } else if (rot_source_type == 3) {

      Real new_ke = 0.5_rt * (snew[UMX] * snew[UMX] + snew[UMY] * snew[UMY] + snew[UMZ] * snew[UMZ]) * rhoInv;
      SrE = new_ke - old_ke;

    } else if (rot_source_type == 4) {

      // The conservative energy formulation does not strictly require
      // any energy source-term here, because it depends only on the
      // fluid motions from the hydrodynamical fluxes which we will only
      // have when we get to the 'corrector' step. Nevertheless we add a
      // predictor energy source term in the way that the other methods
      // do, for consistency. We will fully subtract this predictor value
      // during the corrector step, so that the final result is correct.
      // Here we use the same approach as rot_source_type == 2.

      SrE = uold(i,j,k,UMX) * rhoInv * Sr[UMX] +
            uold(i,j,k,UMY) * rhoInv * Sr[UMY] +
            uold(i,j,k,UMZ) * rhoInv * Sr[UMZ];

    } else {
#ifndef AMREX_USE_GPU
      amrex::Error("Error:: rotation_sources_nd.F90 :: invalid rot_source_type");
#endif
    }

    src[UEDEN] += SrE;

    snew[UEDEN] += dt * src[UEDEN];

    // Add to the outgoing source array.

    for (int n = 0; n < NSRC; n++) {
      source(i,j,k,n) += src[n];
    }

  });

}


void
Castro::corrrsrc(const Box& bx,
                 Array4<Real const> const& phi_old,
                 Array4<Real const> const& phi_new,
                 Array4<Real const> const& rold,
                 Array4<Real const> const& rnew,
                 Array4<Real const> const& uold,
                 Array4<Real const> const& unew,
                 Array4<Real> const& source,
                 Array4<Real const> const& flux1,
                 Array4<Real const> const& flux2,
                 Array4<Real const> const& flux3,
                 const Real dt, const Real time,
                 Array4<Real const> const& vol) {

  // Corrector step for the rotation source terms. This is applied
  // after the hydrodynamics update to fix the time-level n
  // prediction and add the time-level n+1 data.  This subroutine
  // exists outside of the Fortran module above because it needs to
  // be called directly from C++.

  // phi_old and phi_new are used to compute the time centered rotational potential

  // rold and rnew are the old and new time rotational acceleration

  // uold and unew are the old and new time state data

  // source is the source term to send back

  // flux1, flux2, and flux3 are the hydrodynamical mass fluxes


  // Rotation source options for how to add the work to (rho E):
  // rot_source_type =
  // 1: Standard version ("does work")
  // 2: Modification of type 1 that updates the momentum before constructing the energy corrector
  // 3: Puts all rotational work into KE, not (rho e)
  // 4: Conservative energy formulation

  // Note that the time passed to this function
  // is the new time at time-level n+1.

  GpuArray<Real, 3> center;
  ca_get_center(center.begin());

  GpuArray<Real, 3> omega_old;
  get_omega(time-dt, omega_old);

  GpuArray<Real, 3> omega_new;
  get_omega(time, omega_new);

  GpuArray<Real, 3> domegadt_old;
  get_domegadt(time-dt, domegadt_old);

  GpuArray<Real, 3> domegadt_new;
  get_domegadt(time, domegadt_new);

  GeometryData geomdata = geom.data();

  Real dt_omega[3];

  Array2D<Real, 0, 2, 0, 2> dt_omega_matrix = {};

  if (implicit_rotation_update == 1) {

    // Don't do anything here if we've got the Coriolis force disabled.

    if (rotation_include_coriolis == 1) {

      // If the state variables are in the inertial frame, then we are doing
      // an implicit solve using (dt / 2) multiplied by the standard Coriolis term.
      // If not, then the rotation source term to the linear momenta (Equations 16
      // and 17 in Byerly et al., 2014) still retains a Coriolis-like form, with
      // the only difference being that the magnitude is half as large. Consequently
      // we can still do an implicit solve in that case.

      if (state_in_rotating_frame == 1) {

        for (int idir = 0; idir < 3; idir++) {
          dt_omega[idir] = dt * omega_new[idir];
        }

      } else {

        for (int idir = 0; idir < 3; idir++) {
          dt_omega[idir] = 0.5_rt * dt * omega_new[idir];
        }

      }

    } else {

      for (int idir = 0; idir < 3; idir++) {
        dt_omega[idir] = 0.0_rt;
      }

    }


    dt_omega_matrix(0, 0) = 1.0_rt + dt_omega[0] * dt_omega[0];
    dt_omega_matrix(0, 1) = dt_omega[0] * dt_omega[1] + dt_omega[2];
    dt_omega_matrix(0, 2) = dt_omega[0] * dt_omega[2] - dt_omega[1];

    dt_omega_matrix(1, 0) = dt_omega[1] * dt_omega[0] - dt_omega[2];
    dt_omega_matrix(1, 1) = 1.0_rt + dt_omega[1] * dt_omega[1];
    dt_omega_matrix(1, 2) = dt_omega[1] * dt_omega[2] + dt_omega[0];

    dt_omega_matrix(2, 0) = dt_omega[2] * dt_omega[0] + dt_omega[1];
    dt_omega_matrix(2, 1) = dt_omega[2] * dt_omega[1] - dt_omega[0];
    dt_omega_matrix(2, 2) = 1.0_rt + dt_omega[2] * dt_omega[2];

    for (int l = 0; l < 3; l++) {
      for (int m = 0; m < 3; m++) {
        dt_omega_matrix(l, m) /= (1.0_rt + dt_omega[0] * dt_omega[0] +
                                           dt_omega[1] * dt_omega[1] +
                                           dt_omega[2] * dt_omega[2]);
      }
    }

  }

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
  {

    Real Sr_old[3] = {};
    Real Sr_new[3] = {};
    Real Srcorr[3] = {};
    Real src[NSRC] = {};

    // Temporary array for seeing what the new state would be if the update were applied here.

    Real snew[NUM_STATE] = {};

    GpuArray<Real, 3> loc;
    position(i, j, k, geomdata, loc);

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      loc[dir] -= center[dir];
    }

    Real rhoo = uold(i,j,k,URHO);
    Real rhooinv = 1.0_rt / uold(i,j,k,URHO);

    Real rhon = unew(i,j,k,URHO);
    Real rhoninv = 1.0_rt / unew(i,j,k,URHO);

    for (int n = 0; n < NUM_STATE; n++) {
      snew[n] = unew(i,j,k,n);
    }

    Real old_ke = 0.5_rt * (snew[UMX] * snew[UMX] + snew[UMY] * snew[UMY] + snew[UMZ] * snew[UMZ]) * rhoninv;

    // Define old source terms

    Real vold[3];

    vold[0] = uold(i,j,k,UMX) * rhooinv;
    vold[1] = uold(i,j,k,UMY) * rhooinv;
    vold[2] = uold(i,j,k,UMZ) * rhooinv;

    for (int n = 0; n < 3; n++) {
      Sr_old[n] = rhoo * rold(i,j,k,n);
    }

    Real SrE_old = vold[0] * Sr_old[0] + vold[1] * Sr_old[1] + vold[2] * Sr_old[2];

    // Define new source terms

    GpuArray<Real, 3> vnew;

    vnew[0] = unew(i,j,k,UMX) * rhoninv;
    vnew[1] = unew(i,j,k,UMY) * rhoninv;
    vnew[2] = unew(i,j,k,UMZ) * rhoninv;

    for (int n = 0; n < 3; n++) {
      Sr_new[n] = rhon * rnew(i,j,k,n);
    }

    Real SrE_new = vnew[0] * Sr_new[0] + vnew[1] * Sr_new[1] + vnew[2] * Sr_new[2];

    // Define correction terms

    for (int n = 0; n < 3; n++) {
      Srcorr[n] = 0.5_rt * (Sr_new[n] - Sr_old[n]);
    }

    if (implicit_rotation_update == 1) {

      // Coupled/implicit momentum update (wdmerger paper I; Section 2.4)
      // http://adsabs.harvard.edu/abs/2016ApJ...819...94K

      // Do the full corrector step with the old contribution (subtract 1/2 times the old term) and do
      // the non-Coriolis parts of the new contribution (add 1/2 of the new term).

      Real acc[3];
      bool coriolis = false;
      rotational_acceleration(loc, vnew, omega_new, domegadt_new, coriolis, acc);

      Real new_mom_tmp[3];
      for (int n = 0; n < 3; n++) {
        new_mom_tmp[n] = unew(i,j,k,UMX+n) - 0.5_rt * Sr_old[n] * dt + 0.5_rt * rhon * acc[n] * dt;
      }

      // The following is the general solution to the 3D coupled system,
      // assuming that the rotation vector has components along all three
      // axes, obtained using Cramer's rule (the coefficient matrix is
      // defined above). In practice the user will probably only be using
      // one axis for rotation; if it's the z-axis, then this reduces to
      // Equations 25 and 26 in the wdmerger paper. Note that this will
      // have the correct form regardless of whether the state variables are
      // measured in the rotating frame or not; we handled that in the construction
      // of the dt_omega_matrix. It also has the correct form if we have disabled
      // the Coriolis force entirely; at that point it reduces to the identity matrix.

      Real new_mom[3];

      // new_mom = matmul(dt_omega_matrix, new_mom)

      for (int l = 0; l < 3; l++) {
        for (int m = 0; m < 3; m++) {
          new_mom[l] = dt_omega_matrix(l,m) * new_mom_tmp[m];
        }
      }

      // Obtain the effective source term; remember that we're ultimately going
      // to multiply the source term by dt to get the update to the state.

      for (int n = 0; n < 3; n++) {
        Srcorr[n] = (new_mom[n] - unew(i,j,k,UMX+n)) / dt;
      }

    }

    // Correct momenta

    src[UMX] = Srcorr[0];
    src[UMY] = Srcorr[1];
    src[UMZ] = Srcorr[2];

    snew[UMX] += dt * src[UMX];
    snew[UMY] += dt * src[UMY];
    snew[UMZ] += dt * src[UMZ];

#ifdef HYBRID_MOMENTUM
    // The source terms vanish if the state variables are measured in the
    // inertial frame; see wdmerger paper III.

    if (state_in_rotating_frame == 1) {
      set_hybrid_momentum_source(loc, src(UMR:UMP), Srcorr);

      snew[UMR] += dt * src[UMR];
      snew[UML] += dt * src[UML];
      snew[UMP] += dt * src[UMP];
    }
#endif

    // Correct energy

    Real SrEcorr;

    if (rot_source_type == 1) {

      // If rot_source_type == 1, then we calculated SrEcorr before updating the velocities.

      SrEcorr = 0.5_rt * (SrE_new - SrE_old);

    } else if (rot_source_type == 2) {

      // For this source type, we first update the momenta
      // before we calculate the energy source term.

      GpuArray<Real, 3> vnew;

      vnew[0] = snew[UMX] * rhoninv;
      vnew[1] = snew[UMY] * rhoninv;
      vnew[2] = snew[UMZ] * rhoninv;

      Real acc[3];
      bool coriolis = true;
      rotational_acceleration(loc, vnew, omega_new, domegadt_new, coriolis, acc);

      Sr_new[0] = rhon * acc[0];
      Sr_new[1] = rhon * acc[1];
      Sr_new[2] = rhon * acc[2];

      Real SrE_new = vnew[0] * Sr_new[0] + vnew[1] * Sr_new[1] + vnew[2] * Sr_new[2];

      SrEcorr = 0.5_rt * (SrE_new - SrE_old);

    } else if (rot_source_type == 3) {

      // Instead of calculating the energy source term explicitly,
      // we simply update the kinetic energy.

      Real new_ke = 0.5_rt * (snew[UMX] * snew[UMX] + snew[UMY] * snew[UMY] + snew[UMZ] * snew[UMZ]) * rhoninv;
      SrEcorr = new_ke - old_ke;

    } else if (rot_source_type == 4) {

      // Conservative energy update

      // First, subtract the predictor step we applied earlier.

      SrEcorr = - SrE_old;

      // The change in the gas energy is equal in magnitude to, and opposite in sign to,
      // the change in the rotational potential energy, rho * phi.
      // This must be true for the total energy, rho * E_gas + rho * phi, to be conserved.
      // Consider as an example the zone interface i+1/2 in between zones i and i + 1.
      // There is an amount of mass drho_{i+1/2} leaving the zone. From this zone's perspective
      // it starts with a potential phi_i and leaves the zone with potential phi_{i+1/2} =
      // (1/2) * (phi_{i-1}+phi_{i}). Therefore the new rotational energy is equal to the mass
      // change multiplied by the difference between these two potentials.
      // This is a generalization of the cell-centered approach implemented in
      // the other source options, which effectively are equal to
      // SrEcorr = - drho(i,j,k) * phi(i,j,k),
      // where drho(i,j,k) = HALF * (unew(i,j,k,URHO) - uold(i,j,k,URHO)).

      // Note that in the hydrodynamics step, the fluxes used here were already
      // multiplied by dA and dt, so dividing by the cell volume is enough to
      // get the density change (flux * dt * dA / dV). We then divide by dt
      // so that we get the source term and not the actual update, which will
      // be applied later by multiplying by dt.

      Real phi = 0.5_rt * (phi_new(i,j,k) + phi_old(i,j,k));

      Real phixl = 0.5_rt * (phi_new(i-1,j,k) + phi_old(i-1,j,k));
      Real phixr = 0.5_rt * (phi_new(i+1,j,k) + phi_old(i+1,j,k));
      Real phiyl = 0.5_rt * (phi_new(i,j-dg1,k) + phi_old(i,j-dg1,k));
      Real phiyr = 0.5_rt * (phi_new(i,j+dg1,k) + phi_old(i,j+dg1,k));
      Real phizl = 0.5_rt * (phi_new(i,j,k-dg2) + phi_old(i,j,k-dg2));
      Real phizr = 0.5_rt * (phi_new(i,j,k+dg2) + phi_old(i,j,k+dg2));

      SrEcorr = SrEcorr - (0.5_rt / dt) * ( flux1(i    ,j,k) * (phi - phixl) -
                                            flux1(i+1  ,j,k) * (phi - phixr) +
                                            flux2(i,    j,k) * (phi - phiyl) -
                                            flux2(i,j+dg1,k) * (phi - phiyr) +
                                            flux3(i,j,k    ) * (phi - phizl) -
                                            flux3(i,j,k+dg2) * (phi - phizr) ) / vol(i,j,k);

      // Correct for the time rate of change of the potential, which acts
      // purely as a source term. This is only necessary for this source type;
      // it is captured automatically for the others since the time rate of change
      // of omega also appears in the velocity source term.

      GpuArray<Real, 3> cp_tmp;

      cross_product(domegadt_old, loc, cp_tmp);
      for (int idir = 0; idir < 3; idir++) {
        Sr_old[idir] = - rhoo * cp_tmp[idir];
      }

      cross_product(domegadt_new, loc, cp_tmp);
      for (int idir = 0; idir < 3; idir++) {
        Sr_new[idir] = - rhon * cp_tmp[idir];
      }

      Real vnew[3];

      vnew[0] = snew[UMX] * rhoninv;
      vnew[1] = snew[UMY] * rhoninv;
      vnew[2] = snew[UMZ] * rhoninv;

      SrEcorr = SrEcorr + 0.5_rt * (vold[0] * Sr_old[0] + vold[1] * Sr_old[1] + vold[2] * Sr_old[2]) +
                                   (vnew[0] * Sr_new[0] + vnew[1] * Sr_new[1] + vnew[2] * Sr_new[2]);

    } else {
#ifndef AMREX_USE_GPU
      amrex::Error("Error:: rotation_sources_nd.F90 :: invalid rot_source_type");
#endif
    }

    src[UEDEN] = SrEcorr;

    // Add to the outgoing source array.

    for (int n = 0; n < NSRC; n++) {
      source(i,j,k,n) += src[n];
    }

  });

}

