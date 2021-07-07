#include <Castro.H>
#include <Castro_util.H>

#ifdef RADIATION
#include <Radiation.H>
#include <fluxlimiter.H>
#endif

using namespace amrex;

// add the transverse flux difference in direction idir_t to the
// interface states in direction idir_n

void
Castro::trans_single(const Box& bx,
                     int idir_t, int idir_n,
                     Array4<Real const> const& qm,
                     Array4<Real> const& qmo,
                     Array4<Real const> const& qp,
                     Array4<Real> const& qpo,
                     Array4<Real const> const& qaux_arr,
                     Array4<Real const> const& flux_t,
#ifdef RADIATION
                     Array4<Real const> const& rflux_t,
#endif
                     Array4<Real const> const& q_t,
#if AMREX_SPACEDIM == 2
                     Array4<Real const> const& area_t,
                     Array4<Real const> const& vol,
#endif
                     Real hdt, Real cdtdx)
{
    // Evaluate the transverse terms for both
    // the minus and plus states.

    actual_trans_single(bx, idir_t, idir_n, -1,
                        qm, qmo,
                        qaux_arr,
                        flux_t,
#ifdef RADIATION
                        rflux_t,
#endif
                        q_t,
#if AMREX_SPACEDIM == 2
                        area_t,
                        vol,
#endif
                        hdt, cdtdx);

    actual_trans_single(bx, idir_t, idir_n, 0,
                        qp, qpo,
                        qaux_arr,
                        flux_t,
#ifdef RADIATION
                        rflux_t,
#endif
                        q_t,
#if AMREX_SPACEDIM == 2
                        area_t,
                        vol,
#endif
                        hdt, cdtdx);
}


void
Castro::actual_trans_single(const Box& bx,
                            int idir_t, int idir_n, int d,
                            Array4<Real const> const& q_arr,
                            Array4<Real> const& qo_arr,
                            Array4<Real const> const& qaux_arr,
                            Array4<Real const> const& flux_t,
#ifdef RADIATION
                            Array4<Real const> const& rflux_t,
#endif
                            Array4<Real const> const& q_t,
#if AMREX_SPACEDIM == 2
                            Array4<Real const> const& area_t,
                            Array4<Real const> const& vol,
#endif
                            Real hdt, Real cdtdx)
{

    //       qm|qp
    //         |
    // --------+--------
    //   i-1       i
    //        i-1/2
    //
    // the qp state will see the transverse flux in zone i
    // the qm state will see the transverse flux in zone i-1

    // we account for this with the 'd' variable:
    // d = 0 will do qp and d = -1 will do qm

    // idir_t is the transverse direction and we set il,jl,kl
    // and ir,jr,kr to be the face-centered indices needed for
    // the transverse flux difference

    int coord = geom.Coord();

    bool reset_density = transverse_reset_density;
    bool reset_rhoe = transverse_reset_rhoe;
    Real small_p = small_pres;

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {

        // We are handling the states at the interface of
        // (i, i+1) in the x-direction, and similarly for
        // the y- and z- directions.

        int il = i;
        int jl = j;
        int kl = k;

        int ir = i;
        int jr = j;
        int kr = k;

        // set the face indices in the transverse direction

        if (idir_t == 0) {
          ir = i+1;
          jr = j;
          kr = k;

        } else if (idir_t == 1) {
          ir = i;
          jr = j+1;
          kr = k;

        } else {
          ir = i;
          jr = j;
          kr = k+1;
        }

        // We're handling both the plus and minus states;
        // for the minus state we're shifting one zone to
        // the left in our chosen direction.

        if (idir_n == 0) {
          il += d;
          ir += d;

        } else if (idir_n == 1) {
          jl += d;
          jr += d;

        } else {
          kl += d;
          kr += d;
        }

        // construct the conserved variables

        Real U_int[NUM_STATE];

        U_int[URHO] = q_arr(i,j,k,QRHO);
        U_int[UMX] = U_int[URHO] * q_arr(i,j,k,QU);
        U_int[UMY] = U_int[URHO] * q_arr(i,j,k,QV);
        U_int[UMZ] = U_int[URHO] * q_arr(i,j,k,QW);
        U_int[UEDEN] = q_arr(i,j,k,QREINT) +
            0.5_rt * (U_int[UMX] * U_int[UMX] +
                      U_int[UMY] * U_int[UMY] +
                      U_int[UMZ] * U_int[UMZ]) / U_int[URHO];
        U_int[UEINT] = q_arr(i,j,k,QREINT);

        for (int ipassive = 0; ipassive < npassive; ipassive++) {
            int n = upassmap(ipassive);
            int nqp = qpassmap(ipassive);

            U_int[n] = U_int[URHO] * q_arr(i,j,k,nqp);
        }

        // TODO: how does the hybrid stuff fit in here?

        // now do the transverse flux update

#if AMREX_SPACEDIM == 2
        const Real volinv = 1.0_rt / vol(il,jl,kl);
#endif

        for (int n = 0; n < NUM_STATE; n++) {

            if (n == UTEMP) {
                continue;
            }

#ifdef SHOCK_VAR
            if (n == USHK) {
                continue;
            }
#endif

#if AMREX_SPACEDIM == 2
            U_int[n] += - hdt * (area_t(ir,jr,kr) * flux_t(ir,jr,kr,n) -
                                 area_t(il,jl,kl) * flux_t(il,jl,kl,n)) * volinv;
#else
            U_int[n] += - cdtdx * (flux_t(ir,jr,kr,n) - flux_t(il,jl,kl,n));
#endif

        }


        Real pgp  = q_t(ir,jr,kr,GDPRES);
        Real pgm  = q_t(il,jl,kl,GDPRES);
        Real ugp  = q_t(ir,jr,kr,GDU+idir_t);
        Real ugm  = q_t(il,jl,kl,GDU+idir_t);

#if AMREX_SPACEDIM == 2
        Real du = area_t(ir,jr,kr) * ugp - area_t(il,jl,kl) * ugm;
#else
        Real du = ugp - ugm;
#endif
        Real pav = 0.5_rt * (pgp + pgm);

#if AMREX_SPACEDIM == 2

        // For the x-momentum, we need to consider the grad p term.
        // Our y-interface equation for (rho u) is:
        //
        //  d(rho u)/dt + d(rho u v)/dy = - 1/r d(r rho u u)/dr - dp/dr
        //
        // in cylindrical coords -- note that the p term is not
        // in a divergence for UMX in the x-direction, so there
        // are no area factors.  For this geometry, we do not
        // include p in our definition of the flux in the
        // x-direction, for we need to fix this now.

        if (idir_t == 0 && !mom_flux_has_p(0, idir_t, coord)) {
            U_int[UMX] += - cdtdx * (pgp - pgm);
        }
#endif

        // For the dual energy approach, we need to correct the
        // internal energy equation with the p dV term

#if AMREX_SPACEDIM == 2
        U_int[UEINT] += - hdt * pav * du * volinv;
#else
        U_int[UEINT] += - cdtdx * pav * du;
#endif

        // Reset to original value if adding transverse terms made
        // density negative

        bool reset_state = false;
        if (reset_density == 1 && U_int[URHO] < 0.0_rt) {
            reset_state = true;
        }

        if (! reset_state) {

            // Convert back to primitive form

            qo_arr(i,j,k,QRHO) = U_int[URHO];
            Real rhoinv = 1.0_rt / U_int[URHO];

            qo_arr(i,j,k,QU) = U_int[UMX] * rhoinv;
            qo_arr(i,j,k,QV) = U_int[UMY] * rhoinv;
            qo_arr(i,j,k,QW) = U_int[UMZ] * rhoinv;

            Real kineng = 0.5_rt * qo_arr(i,j,k,QRHO) *
                (qo_arr(i,j,k,QU) * qo_arr(i,j,k,QU) +
                 qo_arr(i,j,k,QV) * qo_arr(i,j,k,QV) +
                 qo_arr(i,j,k,QW) * qo_arr(i,j,k,QW));

            if ((U_int[UEDEN] - kineng) > castro::dual_energy_eta1 * U_int[UEDEN]) {
                qo_arr(i,j,k,QREINT) = U_int[UEDEN] - kineng;
            } else {
                qo_arr(i,j,k,QREINT) = U_int[UEINT];
            }

            qo_arr(i,j,k,QREINT) = amrex::max(qo_arr(i,j,k,QREINT), small_dens * small_ener);

            // Load passively advected quatities into q
            for (int ipassive = 0; ipassive < npassive; ipassive++) {
                int n  = upassmap(ipassive);
                int iq = qpassmap(ipassive);
                qo_arr(i,j,k,iq) = U_int[n] * rhoinv;
            }

            eos_rep_t eos_state;

            eos_state.T = T_guess; // the input T may not be valid
            eos_state.rho = qo_arr(i,j,k,QRHO);
            eos_state.e = qo_arr(i,j,k,QREINT) * rhoinv;
            for (int n = 0; n < NumSpec; n++) {
                eos_state.xn[n]  = qo_arr(i,j,k,QFS+n);
            }
#if NAUX_NET > 0
            for (int n = 0; n < NumAux; n++) {
                eos_state.aux[n] = qo_arr(i,j,k,QFX+n);
            }
#endif

            eos(eos_input_re, eos_state);

            qo_arr(i,j,k,QTEMP) = eos_state.T;
            qo_arr(i,j,k,QPRES) = amrex::max(eos_state.p, small_pres);

        } else {

            for (int n = 0; n < NQ; n++) {
                qo_arr(i,j,k,n) = q_arr(i,j,k,n);
            }

        }
    });

}



void
Castro::trans_final(const Box& bx,
                    int idir_n, int idir_t1, int idir_t2,
                    Array4<Real const> const& qm,
                    Array4<Real> const& qmo,
                    Array4<Real const> const& qp,
                    Array4<Real> const& qpo,
                    Array4<Real const> const& qaux_arr,
                    Array4<Real const> const& flux_t1,
#ifdef RADIATION
                    Array4<Real const> const& rflux_t1,
#endif
                    Array4<Real const> const& flux_t2,
#ifdef RADIATION
                    Array4<Real const> const& rflux_t2,
#endif
                    Array4<Real const> const& q_t1,
                    Array4<Real const> const& q_t2,
                    Real cdtdx_t1, Real cdtdx_t2)
{
    // Evaluate the transverse terms for both
    // the minus and plus states.

    actual_trans_final(bx, idir_n, idir_t1, idir_t2, -1,
                       qm, qmo,
                       qaux_arr,
                       flux_t1,
#ifdef RADIATION
                       rflux_t1,
#endif
                       flux_t2,
#ifdef RADIATION
                       rflux_t2,
#endif
                       q_t1,
                       q_t2,
                       cdtdx_t1, cdtdx_t2);

    actual_trans_final(bx, idir_n, idir_t1, idir_t2, 0,
                       qp, qpo,
                       qaux_arr,
                       flux_t1,
#ifdef RADIATION
                       rflux_t1,
#endif
                       flux_t2,
#ifdef RADIATION
                       rflux_t2,
#endif
                       q_t1,
                       q_t2,
                       cdtdx_t1, cdtdx_t2);

}



void
Castro::actual_trans_final(const Box& bx,
                           int idir_n, int idir_t1, int idir_t2, int d,
                           Array4<Real const> const& q_arr,
                           Array4<Real> const& qo_arr,
                           Array4<Real const> const& qaux_arr,
                           Array4<Real const> const& flux_t1,
#ifdef RADIATION
                           Array4<Real const> const& rflux_t1,
#endif
                           Array4<Real const> const& flux_t2,
#ifdef RADIATION
                           Array4<Real const> const& rflux_t2,
#endif
                           Array4<Real const> const& q_t1,
                           Array4<Real const> const& q_t2,
                           Real cdtdx_t1, Real cdtdx_t2)
{

    bool reset_density = transverse_reset_density;
    bool reset_rhoe = transverse_reset_rhoe;
    Real small_p = small_pres;

#ifdef RADIATION
    int fspace_t = Radiation::fspace_advection_type;
    int comov = Radiation::comoving;
    int limiter = Radiation::limiter;
    int closure = Radiation::closure;
#endif

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {

        // the normal state
        int iln = i;
        int jln = j;
        int kln = k;

        // the first transverse state
        int il_t1 = i;
        int jl_t1 = j;
        int kl_t1 = k;

        int ir_t1 = i;
        int jr_t1 = j;
        int kr_t1 = k;

        // the second transverse state
        int il_t2 = i;
        int jl_t2 = j;
        int kl_t2 = k;

        int ir_t2 = i;
        int jr_t2 = j;
        int kr_t2 = k;

        if (idir_n == 0) {

            // x is the normal direction

            // y is the first transverse state
            ir_t1 += d;
            jr_t1 += 1;

            // z is the second transverse state
            ir_t2 += d;
            kr_t2 += 1;

            // offset for the plus/minus state
            iln += d;
            il_t1 += d;
            il_t2 += d;

        }
        else if (idir_n == 1) {

            // y is the normal direction

            // x is the first transverse state
            ir_t1 += 1;
            jr_t1 += d;

            // z is the second transverse state
            jr_t2 += d;
            kr_t2 += 1;

            // offset for the plus/minus state
            jln += d;
            jl_t1 += d;
            jl_t2 += d;

        }
        else {

            // z is the normal direction

            // x is the first transverse state
            ir_t1 += 1;
            kr_t1 += d;

            // y is the second transverse state
            jr_t2 += 1;
            kr_t2 += d;

            // offset for the plus/minus state
            kln += d;
            kl_t1 += d;
            kl_t2 += d;

        }

        // construct the conserved variables

        Real U_int[NUM_STATE];

        U_int[URHO] = q_arr(i,j,k,QRHO);
        U_int[UMX] = U_int[URHO] * q_arr(i,j,k,QU);
        U_int[UMY] = U_int[URHO] * q_arr(i,j,k,QV);
        U_int[UMZ] = U_int[URHO] * q_arr(i,j,k,QW);
        U_int[UEDEN] = q_arr(i,j,k,QREINT) +
            0.5_rt * (U_int[UMX] * U_int[UMX] +
                      U_int[UMY] * U_int[UMY] +
                      U_int[UMZ] * U_int[UMZ]) / U_int[URHO];
        U_int[UEINT] = q_arr(i,j,k,QREINT);

        for (int ipassive = 0; ipassive < npassive; ipassive++) {
            int n = upassmap(ipassive);
            int nqp = qpassmap(ipassive);

            U_int[n] = U_int[URHO] * q_arr(i,j,k,nqp);
        }

        // TODO: how does the hybrid stuff fit in here?

        // now do the transverse flux update

        for (int n = 0; n < NUM_STATE; n++) {

            if (n == UTEMP) {
                continue;
            }

#ifdef SHOCK_VAR
            if (n == USHK) {
                continue;
            }
#endif

            U_int[n] += - cdtdx_t1 * (flux_t1(ir_t1,jr_t1,kr_t1,n) -
                                      flux_t1(il_t1,jl_t1,kl_t1,n))
                        - cdtdx_t2 * (flux_t2(ir_t2,jr_t2,kr_t2,n) -
                                      flux_t2(il_t2,jl_t2,kl_t2,n));

        }

        // Add the transverse differences to the normal states for the
        // fluid variables.

        Real pgt1p  = q_t1(ir_t1,jr_t1,kr_t1,GDPRES);
        Real pgt1m  = q_t1(il_t1,jl_t1,kl_t1,GDPRES);
        Real ugt1p  = q_t1(ir_t1,jr_t1,kr_t1,GDU+idir_t1);
        Real ugt1m  = q_t1(il_t1,jl_t1,kl_t1,GDU+idir_t1);

        Real pgt2p  = q_t2(ir_t2,jr_t2,kr_t2,GDPRES);
        Real pgt2m  = q_t2(il_t2,jl_t2,kl_t2,GDPRES);
        Real ugt2p  = q_t2(ir_t2,jr_t2,kr_t2,GDU+idir_t2);
        Real ugt2m  = q_t2(il_t2,jl_t2,kl_t2,GDU+idir_t2);

        Real pt1av = 0.5_rt * (pgt1p + pgt1m);
        Real dut1 = ugt1p - ugt1m;

        Real pt2av = 0.5_rt * (pgt2p + pgt2m);
        Real dut2 = ugt2p - ugt2m;

        // For the dual energy approach, we need to correct the
        // internal energy equation with the p dV term

        U_int[UEINT] += - cdtdx_t1 * pt1av * dut1 - cdtdx_t2 * pt2av * dut2;

        // Reset to original value if adding transverse terms made
        // density negative

        bool reset_state = false;
        if (reset_density == 1 && U_int[URHO] < 0.0_rt) {
            reset_state = true;
        }


        if (! reset_state) {

            // Convert back to primitive form

            qo_arr(i,j,k,QRHO) = U_int[URHO];
            Real rhoinv = 1.0_rt / U_int[URHO];

            qo_arr(i,j,k,QU) = U_int[UMX] * rhoinv;
            qo_arr(i,j,k,QV) = U_int[UMY] * rhoinv;
            qo_arr(i,j,k,QW) = U_int[UMZ] * rhoinv;

            Real kineng = 0.5_rt * qo_arr(i,j,k,QRHO) *
                (qo_arr(i,j,k,QU) * qo_arr(i,j,k,QU) +
                 qo_arr(i,j,k,QV) * qo_arr(i,j,k,QV) +
                 qo_arr(i,j,k,QW) * qo_arr(i,j,k,QW));

            if ((U_int[UEDEN] - kineng) > castro::dual_energy_eta1 * U_int[UEDEN]) {
                qo_arr(i,j,k,QREINT) = U_int[UEDEN] - kineng;
            } else {
                qo_arr(i,j,k,QREINT) = U_int[UEINT];
            }

            // Load passively advected quatities into q
            for (int ipassive = 0; ipassive < npassive; ipassive++) {
                int n  = upassmap(ipassive);
                int iq = qpassmap(ipassive);
                qo_arr(i,j,k,iq) = U_int[n] * rhoinv;
            }

            eos_rep_t eos_state;

            eos_state.T = T_guess; // the input T may not be valid
            eos_state.rho = qo_arr(i,j,k,QRHO);
            eos_state.e = qo_arr(i,j,k,QREINT) * rhoinv;
            for (int n = 0; n < NumSpec; n++) {
                eos_state.xn[n]  = qo_arr(i,j,k,QFS+n);
            }
#if NAUX_NET > 0
            for (int n = 0; n < NumAux; n++) {
                eos_state.aux[n] = qo_arr(i,j,k,QFX+n);
            }
#endif

            eos(eos_input_re, eos_state);

            qo_arr(i,j,k,QTEMP) = eos_state.T;
            qo_arr(i,j,k,QPRES) = amrex::max(eos_state.p, small_pres);

        } else {

            for (int n = 0; n < NQ; n++) {
                qo_arr(i,j,k,n) = q_arr(i,j,k,n);
            }

        }

    });

}
