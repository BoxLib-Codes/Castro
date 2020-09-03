
// do the advection in frequency space for the radiation energy.
// see Paper III, section 2.4

constexpr int rk_order = 3;
constexpr bool use_WENO = false;

constexpr Real cfl = 0.5_rt;

constexpr Real onethird = 1.0_rt/3.0_rt;
constexpr Real twothirds = 2.0_rt/3.0_rt;
constexpr Real onesixth = 1.0_rt/6.0_rt;

// RK5
constexpr Real B1 = 0.5_rt;
constexpr Real B2 = 1.0_rt/16.0_rt;
constexpr Real B3 = 0.5_rt;
constexpr Real B4 = 9.0_rt/16.0_rt;
constexpr Real B5 = 8.0_rt/7.0_rt;
constexpr Real B6 = 7.0_rt/90.0_rt;

constexpr Real C20 = 5.0_rt/8.0_rt;
constexpr Real C21 = 3.0_rt/8.0_rt;

constexpr Real C40 = 17.0_rt/8.0_rt;
constexpr Real C41 = 9.0_rt/8.0_rt;
constexpr Real C42 = -3.0_rt;
constexpr Real C43 = 0.75_rt;

constexpr Real C50 = -5.0_rt/21.0_rt;
constexpr Real C51 = 2.0_rt/7.0_rt;
constexpr Real C52 = 0.0_rt;
constexpr Real C53 = 4.0_rt;
constexpr Real C54 = -64.0_rt/21.0_rt;

constexpr Real C60 = -8.0_rt/27.0_rt;
constexpr Real C61 = -1.0_rt/5.0_rt;
constexpr Real C62 = 32.0_rt/45.0_rt;
constexpr Real C63 = -32.0_rt/45.0_rt;
constexpr Real C64 = 32.0_rt/27.0_rt;
constexpr Real C65 = 14.0_rt/45.0_rt;


void advect_in_fspace(Real* ustar, Real* af, const Real dt, Real& nstep_fsp) {

  call update_one_species(nGroups, ustar, af, dlognu, dt, nstep_fsp)

}

void update_one_species(const int n, Real* u, Real* a, Real* dx, const Real tend, int& nstepmax) {

  Real dt = 1.e50_rt;
  Real acfl;

  int nstep;

  for (int i = 0; i < n; i++) {
    acfl = 1.e-50_rt + std::abs(a[i]);
    dt = amrex::min(dt, dx[i]/acfl*cfl);
  }

  if (dt >= tend) {
    nstep = 1;
    dt = tend;
  } else {
    nstep = std::ceil(tend/dt);
    dt = tend / static_cast<Real>(nstep);
  }

  Real u1[n];
  Real u2[n];
  Real u3[n];
  Real u4[n];
  Real u5[n];
  Real dudt_tmp[n];

  for (int istep = 1; istep < nstep; istep++) {
    if (rk_order == 5) {
      // RK5
      dudt(u, a, dx, n, dudt_tmp);
      for (int g = 0; g < n; g++) {
        u1[g] = u[g] + B1 * dt * dudt_tmp[g];
      }
      dudt(u1, a, dx, n, dudt_tmp);
      for (int g = 0; g < n; g++) {
        u2[g] = (C20*u[g] + C21*u1[g]) + B2 * dt * dudt_tmp[g];
      }
      dudt(u2, a, dx, n, dudt_tmp);
      for (int g = 0; g < n; g++) {
        u3[g] = u[g] + B3 * dt * dudt_tmp[g];
      }
      dudt(u3, a, dx, n, dudt_tmp);
      for (int g = 0; g < n; g++) {
        u4[g] = (C40*u[g] + C41*u1[g] + C42*u2[g] + C43*u3[g]) + B4 * dt * dudt_tmp[g];
      }
      dudt(u4, a, dx, n, dudt_tmp);
      for (int g = 0; g < n; g++) {
        u5[g] = (C50*u[g] + C51*u1[g] + C52*u2[g] + C53*u3[g] + C54*u4[g]) + B5 * dt * dudt_tmp[g];
      }
      dudt(u5, a, dx, n, dudt_tmp);
      for (int g = 0; g < n; g++) {
        u[g] = (C60*u[g] + C61*u1[g] + C62*u2[g] + C63*u3[g] + C64*u4[g] + C65*u5[g]) + B6 * dt * dudt_tmp[g];
      }

    } else if (rk_order == 4) {
      // RK4
      dudt(u, a, dx, n, dudt_tmp);
      for (int g = 0; g < n; g++) {
        u1[g] = u[g] + 0.5_rt * dt * dudt_tmp[g];
      }
      dudt(u1, a, dx, n, dudt_tmp);
      for (int g = 0; g < n; g++) {
        u2[g] = u[g] + 0.5_rt * dt * dudt_tmp[g];
      }
      dudt(u2, a, dx, n, dudt_tmp);
      for (int g = 0; g < n; g++) {
        u3[g] = u[g] + dt * dudt_tmp[g];
      }
      dudt(u3, a, dx, n, dudt_tmp);
      for (int g = 0; g < n; g++) {
        u[g] = onethird * (u1 + 2.0_rt * u2 + u3 - u[g]) + onesixth * dt * dudt_tmp[g];
      }

    } else if (rk_order == 3) {
      // RK3
      dudt(u, a, dx, n, dudt_tmp);
      for (int g = 0; g < n; g++) {
        u1[g] = u[g] + dt * dudt_tmp[g];
      }
      dudt(u1, a, dx, n, dudt_tmp);
      for (int g = 0; g < n; g++) {
        u1[g] = 0.75e0_rt * u[g] + 0.25_rt * (u1[g] + dt * dudt_tmp[g]);
      }
      dudt(u1, a, dx, n, dudt_tmp);
      for (int g = 0; g < n; g++) {
        u[g] = onethird * u[g] + twothirds * (u1[g] + dt * dudt_tmp[g]);
      }

    } else {
      // first-order
      dudt(u, a, dx, n, dudt_tmp);
      for (int g = 0; g < n; g++) {
        u[g] += dt * dudt_tmp[g];
      }
    }
  }

  nstepmax = amrex::max(nstepmax, nstep);
}


void dudt(Real* u, Real* a, Real* dx, const int n, Real* dudt_tmp) {

  // compute the RHS for the g'th energy group out of n note, since
  // this is advection in frequency space, we will need to assume some
  // boundary conditions on the stencil used to reconstruct over nu

  if (use_WENO) {

       ag(-2) = -a(1)
       ag(-1) = -a(0)
       ag(0:n-1) = a(0:n-1)
       ag(n) = -ag(n-1)
       ag(n+1) = -ag(n-2)

       ug(-2) = u(1)
       ug(-1) = u(0)
       ug(0:n-1) = u(0:n-1)
       ug(n) = ug(n-1)
       ug(n+1) = ug(n-2)

       fg = ag*ug
       ag = abs(ag)

       f(0) = 0.e0_rt
       do i=1,n-1
          alpha = maxval(ag(i-3:i+2))
          fp = 0.5e0_rt * (fg(i-3:i+1) + alpha*ug(i-3:i+1))
          fm = 0.5e0_rt * (fg(i-2:i+2) - alpha*ug(i-2:i+2))
          call weno5(fp(1),fp(2),fp(3),fp(4),fp(5),fpw)
          call weno5(fm(5),fm(4),fm(3),fm(2),fm(1),fmw)
          f(i) = fpw + fmw
       end do
       f(n) = 0.e0_rt

    else

       ag(-1) = -a(0)
       ag(0:n-1) = a(0:n-1)
       ag(n) = -a(n-1)

       ug(-1) = u(0)
       ug(0:n-1) = u(0:n-1)
       ug(n) = u(n-1)

       f(0) = 0.e0_rt
       do i=1,n-1
          r = (ug(i-1)-ug(i-2)) / (ug(i)-ug(i-1) + 1.e-50_rt)
          ul = ug(i-1) + 0.5e0_rt * (ug(i)-ug(i-1)) * MC(r)

          r = (ag(i-1)-ag(i-2)) / (ag(i)-ag(i-1) + 1.e-50_rt)
          al = ag(i-1) + 0.5e0_rt * (ag(i)-ag(i-1)) * MC(r)

          fl = al*ul

          r = (ug(i) - ug(i-1)) / (ug(i+1) - ug(i) + 1.e-50_rt)
          ur = ug(i) - 0.5e0_rt * (ug(i+1) - ug(i)) * MC(r)

          r = (ag(i) - ag(i-1)) / (ag(i+1) - ag(i) + 1.e-50_rt)
          ar = ag(i) - 0.5e0_rt * (ag(i+1) - ag(i)) * MC(r)

          fr = ar*ur

          a_plus = max(0.e0_rt, al, ar)
          a_minus = max(0.e0_rt, -al, -ar)
          f(i) = (a_plus*fl + a_minus*fr - a_plus*a_minus*(ur-ul)) &
                 / (a_plus + a_minus + 1.e-50_rt)
       end do
       f(n) = 0.e0_rt

    end if

    do i = 0, n-1
       dudt(i) = (f(i) - f(i+1)) / dx(i)
    end do

  end function dudt


  function MC(r) result(MCr)

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: r
    real(rt) :: MCr

    !$gpu

    MCr = max(0.e0_rt, min(2.e0_rt*r, 0.5e0_rt*(1.e0_rt+r), 2.e0_rt))

  end function MC


  subroutine weno5(vm2, vm1, v, vp1, vp2, v_weno5)

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: vm2, vm1, v, vp1, vp2
    real(rt), intent(out) :: v_weno5

    real(rt), parameter :: epsw=1.0e-6_rt, b1=13.e0_rt/12.e0_rt, b2=1.e0_rt/6.e0_rt

    real(rt) :: djm1, ejm1, dj, ej, djp1, ejp1, dis0, dis1, dis2, &
                q30, q31, q32, d01, d02, a1ba0, a2ba0, w0, w1, w2

    !$gpu

    djm1 = vm2 - 2.e0_rt*vm1 + v
    ejm1 = vm2 - 4.e0_rt*vm1 + 3.e0_rt*v
    dj   = vm1 - 2.e0_rt*v + vp1
    ej   = vm1 - vp1
    djp1 = v - 2.e0_rt*vp1 + vp2
    ejp1 = 3.e0_rt*v - 4.e0_rt*vp1 + vp2

    dis0 = b1*djm1*djm1 + 0.25e0_rt*ejm1*ejm1 + epsw
    dis1 = b1*dj*dj     + 0.25e0_rt*ej*ej     + epsw
    dis2 = b1*djp1*djp1 + 0.25e0_rt*ejp1*ejp1 + epsw

    q30 = 2.e0_rt*vm2 - 7.e0_rt*vm1 + 11.e0_rt*v
    q31 = -vm1 + 5.e0_rt*v + 2.e0_rt*vp1
    q32 = 2.e0_rt*v + 5.e0_rt*vp1 - vp2

    d01 = dis0 / dis1
    d02 = dis0 / dis2
    a1ba0 = 6.e0_rt * d01 * d01
    a2ba0 = 3.e0_rt * d02 * d02
    w0 = 1.e0_rt / (1.e0_rt + a1ba0 + a2ba0)
    w1 = a1ba0 * w0
    w2 = 1.e0_rt - w0 - w1

    if (w0.lt.1.0e-10_rt) w0 = 0.e0_rt
    if (w1.lt.1.0e-10_rt) w1 = 0.e0_rt
    if (w2.lt.1.0e-10_rt) w2 = 0.e0_rt

    v_weno5 = b2*(w0*q30 + w1*q31 + w2*q32)

  end subroutine weno5
