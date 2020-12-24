#include <Castro.H>
#include <Castro_F.H>
#include <Castro_util.H>
#include <Rotation.H>
#ifdef HYBRID_MOMENTUM
#include <hybrid.H>
#endif

void
Castro::rsrc(const Box& bx,
             Array4<Real const> const& uold,
             Array4<Real> const& source, 
             const Real dt) {

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
      loc[dir] -= problem::center[dir];
    }

    Real rho = uold(i,j,k,URHO);
    Real rhoInv = 1.0_rt / rho;

    for (int n = 0; n < NUM_STATE; n++) {
      snew[n] = uold(i,j,k,n);
    }

    Real old_ke = 0.5_rt * (snew[UMX] * snew[UMX] + snew[UMY] * snew[UMY] + snew[UMZ] * snew[UMZ]) * rhoInv;

    GpuArray<Real, 3> v;

    v[0] = uold(i,j,k,UMX) * rhoInv;
    v[1] = uold(i,j,k,UMY) * rhoInv;
    v[2] = uold(i,j,k,UMZ) * rhoInv;

    bool coriolis = true;
    rotational_acceleration(loc, v, coriolis, Sr);

    for (int n = 0; n < 3; n++) {
        Sr[n] = rho * Sr[n];
    }

    src[UMX] = Sr[0];
    src[UMY] = Sr[1];
    src[UMZ] = Sr[2];

    snew[UMX] += dt * src[UMX];
    snew[UMY] += dt * src[UMY];
    snew[UMZ] += dt * src[UMZ];

#ifdef HYBRID_MOMENTUM
    if (state_in_rotating_frame == 1) {

      GpuArray<Real, 3> linear_momentum;
      linear_momentum[0] = src[UMX];
      linear_momentum[1] = src[UMY];
      linear_momentum[2] = src[UMZ];

      GpuArray<Real, 3> hybrid_source;
      set_hybrid_momentum_source(loc, linear_momentum, hybrid_source);

      snew[UMR] += dt * hybrid_source[0];
      snew[UML] += dt * hybrid_source[1];
      snew[UMP] += dt * hybrid_source[2];

      src[UMR] = hybrid_source[0];
      src[UML] = hybrid_source[1];
      src[UMP] = hybrid_source[2];

    }
#endif

    // Kinetic energy source: this is v . the momentum source.
    // We don't apply in the case of the conservative energy
    // formulation.

    Real SrE;

    if (rot_source_type == 1 || rot_source_type == 2) {

      SrE = uold(i,j,k,UMX) * rhoInv * Sr[0] +
            uold(i,j,k,UMY) * rhoInv * Sr[1] +
            uold(i,j,k,UMZ) * rhoInv * Sr[2];

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

      SrE = uold(i,j,k,UMX) * rhoInv * Sr[0] +
            uold(i,j,k,UMY) * rhoInv * Sr[1] +
            uold(i,j,k,UMZ) * rhoInv * Sr[2];

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
                 Array4<Real const> const& uold,
                 Array4<Real const> const& unew,
                 Array4<Real> const& source,
                 Array4<Real const> const& flux0,
                 Array4<Real const> const& flux1,
                 Array4<Real const> const& flux2,
                 const Real dt,
                 Array4<Real const> const& vol) {

  // Corrector step for the rotation source terms. This is applied
  // after the hydrodynamics update to fix the time-level n
  // prediction and add the time-level n+1 data.  This subroutine
  // exists outside of the Fortran module above because it needs to
  // be called directly from C++.

  // phi_old and phi_new are used to compute the time centered rotational potential

  // uold and unew are the old and new time state data

  // source is the source term to send back

  // flux0, flux1, and flux2 are the hydrodynamical mass fluxes


  // Rotation source options for how to add the work to (rho E):
  // rot_source_type =
  // 1: Standard version ("does work")
  // 2: Modification of type 1 that updates the momentum before constructing the energy corrector
  // 3: Puts all rotational work into KE, not (rho e)
  // 4: Conservative energy formulation

  // Note that the time passed to this function
  // is the new time at time-level n+1.

  GeometryData geomdata = geom.data();

  Real hdtInv = 0.5_rt / dt;

  GpuArray<Real, 3> dx;
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      dx[i] = geom.CellSizeArray()[i];
  }
  for (int i = AMREX_SPACEDIM; i < 3; ++i) {
      dx[i] = 0.0_rt;
  }

  auto omega = get_omega();

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
          dt_omega[idir] = dt * omega[idir];
        }

      } else {

        for (int idir = 0; idir < 3; idir++) {
          dt_omega[idir] = 0.5_rt * dt * omega[idir];
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
      loc[dir] -= problem::center[dir];
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

    GpuArray<Real, 3> vold;

    vold[0] = uold(i,j,k,UMX) * rhooinv;
    vold[1] = uold(i,j,k,UMY) * rhooinv;
    vold[2] = uold(i,j,k,UMZ) * rhooinv;

    bool coriolis = true;
    rotational_acceleration(loc, vold, coriolis, Sr_old);

    for (int n = 0; n < 3; n++) {
        Sr_old[n] = rhoo * Sr_old[n];
    }

    Real SrE_old = vold[0] * Sr_old[0] + vold[1] * Sr_old[1] + vold[2] * Sr_old[2];


    // Define new source terms

    GpuArray<Real, 3> vnew;

    vnew[0] = unew(i,j,k,UMX) * rhoninv;
    vnew[1] = unew(i,j,k,UMY) * rhoninv;
    vnew[2] = unew(i,j,k,UMZ) * rhoninv;

    rotational_acceleration(loc, vnew, coriolis, Sr_new);

    for (int n = 0; n < 3; n++) {
        Sr_new[n] = rhon * Sr_new[n];
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
      coriolis = false;
      rotational_acceleration(loc, vnew, coriolis, acc);

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

      Real new_mom[3] = {}; 

      // new_mom = matmul(dt_omega_matrix, new_mom)

      for (int l = 0; l < 3; l++) {
        for (int m = 0; m < 3; m++) {
          new_mom[l] += dt_omega_matrix(l,m) * new_mom_tmp[m];
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
      GpuArray<Real, 3> hybrid_source;

      GpuArray<Real, 3> linear_momentum;
      linear_momentum[0] = src[UMX];
      linear_momentum[1] = src[UMY];
      linear_momentum[2] = src[UMZ];

      set_hybrid_momentum_source(loc, linear_momentum, hybrid_source);

      snew[UMR] += dt * hybrid_source[0];
      snew[UML] += dt * hybrid_source[1];
      snew[UMP] += dt * hybrid_source[2];

      src[UMR] = hybrid_source[0];
      src[UML] = hybrid_source[1];
      src[UMP] = hybrid_source[2];
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
      coriolis = true;
      rotational_acceleration(loc, vnew, coriolis, acc);

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

      // See corrgsrc for an explanation of this algorithm. We can implement this identically
      // to the conservative gravity, swapping out the gravitational acceleration for the
      // rotational acceleration. We only need to calculate the rotational acceleration on
      // edges, which is straightforward to do. The Coriolis term does not play a role in
      // the energy update (it does no work), so we turn it off and pass in a temporary
      // velocity array which will be ignored.

      Real edge_Sr[3][2];

      // Loop over directions and edges

      for (int dir = 0; dir < 3; ++dir) {
          for (int edge = 0; edge <= 1; ++edge) {

              bool ccx = dir == 0 ? true : false;
              bool ccy = dir == 1 ? true : false;
              bool ccz = dir == 2 ? true : false;

              int ie = dir == 0 ? i + edge : i;
              int je = dir == 1 ? j + edge : j;
              int ke = dir == 2 ? k + edge : k;

              position(ie, je, ke, geomdata, loc, ccx, ccy, ccz);

              for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                  loc[n] -= problem::center[n];
              }

              GpuArray<Real, 3> temp_vel{};
              Real temp_Sr[3];

              coriolis = false;
              rotational_acceleration(loc, temp_vel, coriolis, temp_Sr);

              edge_Sr[dir][edge] = temp_Sr[dir];

          }
      }

      SrEcorr += hdtInv * (flux0(i      ,j,k) * edge_Sr[0][0] * dx[0] +
                           flux0(i+1*dg0,j,k) * edge_Sr[0][1] * dx[0] +
                           flux1(i,j      ,k) * edge_Sr[1][0] * dx[1] +
                           flux1(i,j+1*dg1,k) * edge_Sr[1][1] * dx[1] +
                           flux2(i,j,k      ) * edge_Sr[2][0] * dx[2] +
                           flux2(i,j,k+1*dg2) * edge_Sr[2][1] * dx[2]) / vol(i,j,k);

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

