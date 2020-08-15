#include "Castro.H"
#include "Castro_F.H"
#include "Castro_util.H"
#include "Castro_hydro_F.H"

#include "mhd_util.H"

using namespace amrex;

void
Castro::get_inplane_Bs_transverse_flux(const int face, const int comp, const int tdir,
                                       const int i, const int j, const int k,
                                       Array4<Real const> const& Ex,
                                       Array4<Real const> const& Ey,
                                       Array4<Real const> const& Ez,
                                       Real& Fl, Real& Fr) {

  // this is for the 2 B components that are on a face (e.g., y and z on the x faces)
  // face is the coordinate direction of the face
  // comp is the component of the B field we are dealing with
  // tdir is the transverse direction (orthogonal to comp)

  // only the fields tangental to the face direction (in-plane)
  BL_ASSERT(face != comp);

  // for a field pointing in direction comp, there is no transverse update in the comp direction
  BL_ASSERT(comp != tdir);


  // we are working on a face, and considering B field component in
  // the plane of the face.  there are 2 transverse flux differences
  // that we need to account for.

  if (face == 0) {
    // on the x-face, the y and z components are in plane

    if (comp == 1) {
      // By on the x-face should get corrected by an x and z flux difference

      if (tdir == 0) {
        // Fx_{i+1/2,j,k}(By) = - Ez_{i+1/2,j,k}
        //                    = -1/2 [ Ez_{i+1/2,j-1/2,k} + Ez_{i+1/2,j+1/2,k} ]

        Fr = -0.5_rt * (Ez(i+1,j,k) + Ez(i+1,j+1,k));
        Fl = -0.5_rt * (Ez(i,j,k) + Ez(i,j+1,k));

      } else if (tdir == 2) {
        // Fz_{i,j,k+1/2}(By) = Ex_{i,j,k+1/2}
        //                    = 1/2 [ Ex_{i,j-1/2,k+1/2} + Ex_{i,j+1/2,k+1/2} ]

        Fr = 0.5_rt * (Ex(i,j,k+1) + Ex(i,j+1,k+1));
        Fl = 0.5_rt * (Ex(i,j,k) + Ex(i,j+1,k));

      }

    } else if (comp == 2) {
      // Bz on the x-face should get corrected by an x and y flux difference

      if (tdir == 0) {
        // Fx_{i+1/2,j,k}(Bz) = Ey_{i+1/2,j,k}
        //                    = 1/2 [ Ey_{i+1/2,j,k-1/2} + Ey_{i+1/2,j,k+1/2} ]

        Fr = 0.5_rt * (Ey(i+1,j,k) + Ey(i+1,j,k+1));
        Fl = 0.5_rt * (Ey(i,j,k) + Ey(i,j,k+1));

      } else if (tdir == 1) {
        // Fy_{i,j+1/2,k}(Bz) = -Ex_{i,j+1/2,k}
        //                    = -1/2 [ Ex_{i,j+1/2,k-1/2} + Ex_{i,j+1/2,k+1/2} ]

        Fr = -0.5_rt * (Ex(i,j+1,k) + Ex(i,j+1,k+1));
        Fl = -0.5_rt * (Ex(i,j,k) + Ex(i,j,k+1));

      }

    }

  } else if (face == 1) {

    if (comp == 0) {
      // Bx on the y-face should get corrected by an y and z flux difference

      if (tdir == 1) {
        // Fy_{i,j+1/2,k}(Bx) = Ez_{i,j+1/2,k}
        //                    = 1/2 [ Ez_{i-1/2,j+1/2,k} + Ez_{i+1/2,j+1/2,k} ]

        Fr = 0.5_rt * (Ez(i,j+1,k) + Ez(i+1,j+1,k));
        Fl = 0.5_rt * (Ez(i,j,k) + Ez(i+1,j,k));

      } else if (tdir == 2) {
        // Fz_{i,j,k+1/2}(Bx) = -Ey_{i,j,k+1/2}
        //                    = -1/2 (Ey_{i-1/2,j,k+1/2} + Ey_{i+1/2,j,k+1/2})

        Fr = -0.5_rt * (Ey(i,j,k+1) + Ey(i+1,j,k+1));
        Fl = -0.5_rt * (Ey(i,j,k) + Ey(i+1,j,k));
      }

    } else if (comp == 2) {
      // Bz on the y-face should get corrected by an x and y flux difference

      if (tdir == 0) {
        // Fx_{i+1/2,j,k}(Bz) = Ey_{i+1/2,j,k}
        //                    = 1/2 [ Ey_{i+1/2,j,k-1/2} + Ey_{i+1/2,j,k+1/2} ]

        Fr = 0.5_rt * (Ey(i+1,j,k) + Ey(i+1,j,k+1));
        Fl = 0.5_rt * (Ey(i,j,k) + Ey(i,j,k+1));

      } else if (tdir == 1) {
        // Fy_{i,j+1/2,k}(Bz) = -Ex_{i,j+1/2,k}
        //                    = -1/2 [ Ex_{i,j+1/2,k-1/2} + Ex_{i,j+1/2,k+1/2} ]

        Fr = -0.5_rt * (Ex(i,j+1,k) + Ex(i,j+1,k+1));
        Fl = -0.5_rt * (Ex(i,j,k) + Ex(i,j,k+1));
      }

    }

  } else {  // face == 2

    if (comp == 0) {
      //  Bx on the z-face should get corrected by a y and z flux difference

      if (tdir == 1) {
        // Fy_{i,j+1/2,k}(Bx) = Ez_{i,j+1/2,k}
        //                    = 1/2 [ Ez_{i-1/2,j+1/2,k} + Ez_{i+1/2,j+1/2,k}

        Fr = 0.5_rt * (Ez(i,j+1,k) + Ez(i+1,j+1,k));
        Fl = 0.5_rt * (Ez(i,j,k) + Ez(i+1,j,k));

      } else if (tdir == 2) {
        // Fz_{i,j,k+1/2}(Bx) = -Ey_{i,j,k+1/2}
        //                    = -1/2 [ Ey_{i-1/2,j,k+1/2} + Ey_{i+1/2,j,k+1/2} ]

        Fr = -0.5_rt * (Ey(i,j,k+1) + Ey(i+1,j,k+1));
        Fl = -0.5_rt * (Ey(i,j,k) + Ey(i+1,j,k));
      }

    } else if (comp == 1) {
      // By on the z-face should get corrected by a x and z flux difference

      if (tdir == 0) {
        // Fx_{i+1/2,j,k}(By) = -Ez_{i+1/2,j,k}
        //                    = -1/2 [ Ez_{i+1/2,j-1/2,k} + Ez_{i+1/2,j+1/2,k} ]

        Fr = -0.5_rt * (Ez(i+1,j,k) + Ez(i+1,j+1,k));
        Fl = -0.5_rt * (Ez(i,j,k) + Ez(i,j+1,k));

      } else if (tdir == 2) {
        // Fz_{i,j,k+1/2}(By) = Ex_{i,j,k+1/2}
        //                    = 1/2 [ Ex_{i,j-1/2,k+1/2} + Ex_{i,j+1/2,k+1/2} ]

        Fr = 0.5_rt * (Ex(i,j,k+1) + Ex(i,j+1,k+1));
        Fl = 0.5_rt * (Ex(i,j,k) + Ex(i,j+1,k));
      }

    }

  }

}

void
Castro::corner_couple(const Box& bx,
                      Array4<Real> const& qr_out,
                      Array4<Real> const& ql_out,
                      Array4<Real const> const& ur,
                      Array4<Real const> const& ul,
                      Array4<Real const> const& flxd2,
                      Array4<Real const> const& Ed1,
                      Array4<Real const> const& Ed3,
                      Array4<Real const> const& Ex,
                      Array4<Real const> const& Ey,
                      Array4<Real const> const& Ez,
                      const int d1, const int d2, const int d3,
                      const Real dt) {

  // take conservative interface states ul and ur and update them
  // with with the transverse flux difference (corner coupling) to
  // produce ul_out and ur_out
  //
  // This implements MM step 3 of the CTU algorithm.
  // the normal direction (for the interface states) is d1
  // the transverse direction (for the flux difference) is d2

  GpuArray<Real, 3> dx;
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      dx[i] = geom.CellSizeArray()[i];
  }
  for (int i = AMREX_SPACEDIM; i < 3; ++i) {
      dx[i] = 0.0_rt;
  }

  // cl and cr are the offsets to the indices for the conserved state fluxes
  // they will be offset in d2 to capture the flux difference

  int cl[3] = {};
  int cr[3] = {};

  // b is for indexing into the magnetic field that is in the d1
  // (normal) direction

  int b[3] = {};


  // update the state on interface direction d1 with the input flux in direction d2

  // for the flux difference, F_r - F_l, we need to shift the indices in the first flux (F_r)
  // in d2 to get a difference across the interface.  We also need to shift by a zone in d1
  // for the left interface.  cr(:) and cl(:) will hold these shifts.

  // the first term of the flxd2 substraction is shifted by 1 on the direction d2
  cr[d2] = 1;

  // for the normal B component
  b[d2] = 1;


  // ell is the lower-left E in both directions
  // (e.g. Ez_{i-1/2,j-1/2,k})

  int sgn = epsilon_ijk(d1, d2, d3);
  Real cdtdx = dt/(3.0_rt * dx[d1]);

  // UMAGD1 corresponds to d1, and UMAGD2 to d2, UMAGD3 to d3
  int UMAGD1 = UMAGX + d1;
  int UMAGD2 = UMAGX + d2;
  int UMAGD3 = UMAGX + d3;

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
  {

    // first the conserved state

    // right interface (e.g. U_{i-1/2,j,k,R} or the "-" state in MM notation)
    // MM Eq. 37

    Real utmp[NUM_STATE+3];
    Real qtmp[NQ];

    for (int n = 0; n < NUM_STATE; n++) {
      if (n != UTEMP) {
        utmp[n] = ur(i,j,k,n) -
          cdtdx*(flxd2(i+cr[0],j+cr[1],k+cr[2],n) - flxd2(i,j,k,n));

      } else {
        utmp[n] = 0.0_rt;
      }
    }

    // now magnetic fields

    // right state on the interface (e.g. B_{i-1/2,j,k,R} or `-` in MM notation)

    // d1 -- this is perpendicular to the face, MM Eq. 38
    // (note MM Eq. 38 has a sign error) e.g., for d1 = x and
    // d2 = y, this gets updated as
    //
    // Bx|y_{i-1/2,j,k,R} = Bx_{i-1/2,j,k) -
    //     dt/3dx (Ez_{i-1/2,j+1/2,k) - Ez_{i-1/2,j-1/2,k})
    //
    // we use b[] to captured the j+1/2 indexing into Ez

    utmp[UMAGD1] = ur(i,j,k,UMAGD1) - sgn * cdtdx *
      (Ed3(i+b[0],j+b[1],k+b[2]) - Ed3(i,j,k));

    // d3 -- this is in the plane of the face, MM Eq. 39
    // e.g.,g for d1 = x, and d2 = y, this gets updated as
    //
    // Bz|y_{i-1/2,j,k,R} = Bz_{i-1/2,j,k} + 1/2 dt/3dx
    //     (Ex_{i,j+1/2,k+1/2} - Ex_{i,j-1/2,k+1/2} +
    //      Ex_{i,j+1/2,k-1/2} - Ex_{i,j-1/2,k-1/2})
    //
    // we use err(:) for the first E term, elr(:) for the
    // second, and erl(:) for the third

    Real Fl;
    Real Fr;

    get_inplane_Bs_transverse_flux(d1, d3, d2, i, j, k,
                                   Ex, Ey, Ez,
                                   Fl, Fr);

    utmp[UMAGD3] = ur(i,j,k,UMAGD3) - cdtdx * (Fr - Fl);

    // the component pointing in the transverse update direction, d2, is unchanged
    utmp[UMAGD2] = ur(i,j,k,UMAGD2);

    ConsToPrim(qtmp, utmp);

    for (int n = 0; n < NQ; n++) {
      qr_out(i,j,k,n) = qtmp[n];
    }
  });


  // left interface (e.g., U_{i-1/2,j,k,L} or the "+" state in MM notation)
  // note: this uses information one zone to the left in d1

  cl[d1] -= 1;
  cr[d1] -= 1;

  // The in-plane B component at B_{i-1/2,j,k,L} uses the information one zone to the left
  // in direction d1


  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
  {

    // conservative state

    Real utmp[NUM_STATE+3];
    Real qtmp[NQ];

    for (int n = 0; n < NUM_STATE; n++) {
      if (n != UTEMP) {
        utmp[n] = ul(i,j,k,n) -
          cdtdx*(flxd2(i+cr[0],j+cr[1],k+cr[2],n) -
                 flxd2(i+cl[0],j+cl[1],k+cl[2],n));

      } else {
        utmp[n] = 0.0_rt;
      }
    }

    // left state on the interface (e.g. B_{i-1/2,j,k,L} or `+` in MM notation)

    int ii = i;
    int jj = j;
    int kk = k;

    if (d1 == 0) {
      ii -= 1;
    } else if (d1 == 1) {
      jj -= 1;
    } else {
      kk -= 1;
    }

    utmp[UMAGD1] = ul(i,j,k,UMAGD1) - sgn * cdtdx *
      (Ed3(i+b[0],j+b[1],k+b[2]) - Ed3(i,j,k));

    Real Fl;
    Real Fr;

    get_inplane_Bs_transverse_flux(d1, d3, d2, ii, jj, kk,
                                   Ex, Ey, Ez,
                                   Fl, Fr);

    utmp[UMAGD3] = ul(i,j,k,UMAGD3) - cdtdx * (Fr - Fl);

    utmp[UMAGD2] = ul(i,j,k,UMAGD2);

    ConsToPrim(qtmp, utmp);

    for (int n = 0; n < NQ; n++) {
      ql_out(i,j,k,n) = qtmp[n];
    }
  });
}


void
Castro::half_step(const Box& bx,
                  Array4<Real> const& qr_out,
                  Array4<Real> const& ql_out,
                  Array4<Real const> const& ur,
                  Array4<Real const> const& ul,
                  Array4<Real const> const& flxd1,
                  Array4<Real const> const& flxd2,
                  Array4<Real const> const& Ed,
                  Array4<Real const> const& Ed1,
                  Array4<Real const> const& Ed2,
                  Array4<Real const> const& Ex,
                  Array4<Real const> const& Ey,
                  Array4<Real const> const& Ez,
                  const int d, const int d1, const int d2, const Real dt) {

  // Final transverse flux corrections to the conservative state

  GpuArray<Real, 3> dx;
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      dx[i] = geom.CellSizeArray()[i];
  }
  for (int i = AMREX_SPACEDIM; i < 3; ++i) {
      dx[i] = 0.0_rt;
  }

  // c1l, c1r are for indexing flxd1 offsets, c2l, c2r are for flxd2

  int c1l[3] = {};
  int c1r[3] = {};
  int c2l[3] = {};
  int c2r[3] = {};

  int a1[3] = {};
  int a2[3] = {};

  Real hdtdx = 0.5_rt * dt / dx[d];
  int sgn = -1 * epsilon_ijk(d, d1, d2);

  // UMAGD1 corresponds to d1, and UMAGD2 to d2, UMAGD3 to d3
  int UMAGD = UMAGX + d;
  int UMAGD1 = UMAGX + d1;
  int UMAGD2 = UMAGX + d2;


  c1r[d1] = 1;  // add +1 to the d1 direction in the first flxd1 term of the subtraction
  c2r[d2] = 1;  // add +1 to the d2 direction in the first flxd2 term of the subtraction

  // for the normal component of B
  a1[d2] = 1;  // shift on first term of Ed1 substraction, in d2 direction
  a2[d1] = 1;  // shift on first term of Ed2 substraction, in d1 direction


  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
  {

    // first the conservative state

    // right interface (e.g. U_{i-1/2,j,k,R} or the "-" state in MM notation)
    // MM Eq. 44

    Real utmp[NUM_STATE+3];
    Real qtmp[NQ];

    for (int n = 0; n < NUM_STATE; n++) {
      if (n != UTEMP) {
        utmp[n] = ur(i,j,k,n) -
          hdtdx * (flxd1(i+c1r[0],j+c1r[1],k+c1r[2],n) - flxd1(i,j,k,n)) -
          hdtdx * (flxd2(i+c2r[0],j+c2r[1],k+c2r[2],n) - flxd2(i,j,k,n));

      } else {
        utmp[n] = 0.0_rt;
      }
    }

    // now the magnetic field

    // right state on the interface (e.g. B_{i-1/2,j,k,R} or `-` in MM notation)

    // Bd -- this is perpendicular to the face. Note MM eq.45
    // in Miniati has a sign error in the epsilon term

    utmp[UMAGD] = ur(i,j,k,UMAGD) - sgn * hdtdx *
      ((Ed1(i+a1[0],j+a1[1],k+a1[2]) - Ed1(i,j,k)) -
       (Ed2(i+a2[0],j+a2[1],k+a2[2]) - Ed2(i,j,k)));

    // Bd1 -- this is one of the components of B in the plane of the face d
    // Eq.46 in Miniati

    Real F1l;
    Real F1r;
    Real F2l;
    Real F2r;

    get_inplane_Bs_transverse_flux(d, d1, d, i, j, k,
                                   Ex, Ey, Ez,
                                   F1l, F1r);

    get_inplane_Bs_transverse_flux(d, d1, d2, i, j, k,
                                   Ex, Ey, Ez,
                                   F2l, F2r);

    utmp[UMAGD1] = ur(i,j,k,UMAGD1) - hdtdx * ((F1r - F1l) + (F2r - F2l));

    get_inplane_Bs_transverse_flux(d, d2, d, i, j, k,
                                   Ex, Ey, Ez,
                                   F1l, F1r);

    get_inplane_Bs_transverse_flux(d, d2, d1, i, j, k,
                                   Ex, Ey, Ez,
                                   F2l, F2r);

    utmp[UMAGD2] = ur(i,j,k,UMAGD2) - hdtdx * ((F1r - F1l) + (F2r - F2l));

    // convert to primitive

    ConsToPrim(qtmp, utmp);

    for (int n = 0; n < NQ; n++) {
      qr_out(i,j,k,n) = qtmp[n];
    }
  });


  // for the left state U components on the face dir, the flux
  // difference is in the zone to the left

    c1r[d] -= 1;
    c1l[d] -= 1;
    c2r[d] -= 1;
    c2l[d] -= 1;

    // left state on the interface (e.g., B_{i-1/2,j,k,L} or `+` in MM notation)

    // The in-plane B component at B_{i-1/2,j,k,L} uses the information one zone to the left
    // in direction d

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
  {

    // left interface (e.g., U_{i+1/2,j,k,L} or the "+" state in MM notation)

    Real utmp[NUM_STATE+3];
    Real qtmp[NQ];

    // first the conserved state

    for (int n = 0; n < NUM_STATE; n++) {
      if (n != UTEMP) {
        utmp[n] = ul(i,j,k,n) -
          hdtdx * (flxd1(i+c1r[0],j+c1r[1],k+c1r[2],n) -
                   flxd1(i+c1l[0],j+c1l[1],k+c1l[2],n)) -
          hdtdx * (flxd2(i+c2r[0],j+c2r[1],k+c2r[2],n) -
                   flxd2(i+c2l[0],j+c2l[1],k+c2l[2],n));

      } else {
        utmp[n] = 0.0_rt;
      }
    }

    // now the B fields

    // Bd -- this is perpendicular to the face (MM Eq. 45 with sign fix)

    // this is the same face as the right state, so the update the identical
    utmp[UMAGD] = ul(i,j,k,UMAGD) - sgn * hdtdx *
      ((Ed1(i+a1[0],j+a1[1],k+a1[2]) - Ed1(i,j,k)) -
       (Ed2(i+a2[0],j+a2[1],k+a2[2]) - Ed2(i,j,k)));


    int ii = i;
    int jj = j;
    int kk = k;

    if (d == 0) {
      ii -= 1;
    } else if (d == 1) {
      jj -= 1;
    } else {
      kk -= 1;
    }

    Real F1l;
    Real F1r;
    Real F2l;
    Real F2r;

    // Bd1 -- first component on face d, eq. 46 in Miniati

    get_inplane_Bs_transverse_flux(d, d1, d, ii, jj, kk,
                                   Ex, Ey, Ez,
                                   F1l, F1r);

    get_inplane_Bs_transverse_flux(d, d1, d2, ii, jj, kk,
                                   Ex, Ey, Ez,
                                   F2l, F2r);

    utmp[UMAGD1] = ul(i,j,k,UMAGD1) - hdtdx * ((F1r - F1l) + (F2r - F2l));

    // Bd2 -- second component on face d, eq. 46 in Miniati

    get_inplane_Bs_transverse_flux(d, d2, d, ii, jj, kk,
                                   Ex, Ey, Ez,
                                   F1l, F1r);

    get_inplane_Bs_transverse_flux(d, d2, d1, ii, jj, kk,
                                   Ex, Ey, Ez,
                                   F2l, F2r);

    utmp[UMAGD2] = ul(i,j,k,UMAGD2) - hdtdx * ((F1r - F1l) + (F2r - F2l));

    // convert to primitive

    ConsToPrim(qtmp, utmp);

    for (int n = 0; n < NQ; n++) {
      ql_out(i,j,k,n) = qtmp[n];
    }
  });

}
