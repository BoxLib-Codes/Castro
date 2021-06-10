#include <Castro.H>
#include <Castro_F.H>

#include <riemann_solvers.H>

#ifdef RADIATION
#include <Radiation.H>
#endif

#include <cmath>

#include <eos.H>
using namespace amrex;

void
Castro::cmpflx_plus_godunov(const Box& bx,
                            Array4<Real> const& qm,
                            Array4<Real> const& qp,
                            Array4<Real> const& flx,
#ifdef RADIATION
                            Array4<Real> const& rflx,
#endif
                            Array4<Real> const& qgdnv,
                            Array4<Real const> const& qaux_arr,
                            Array4<Real const> const& shk,
                            const int idir, const bool store_full_state) {

    // note: bx is not necessarily the limits of the valid (no ghost
    // cells) domain, but could be hi+1 in some dimensions.  We rely on
    // the caller to specify the interfaces over which to solve the
    // Riemann problems

    // Solve Riemann problem to get the fluxes

    // store_full_state determines what is put into qgdnv.  if
    // store_full_state is True, we put all NQ variables into the
    // qgdnv.  if store_full_state is False, we only store the NGDNV
    // needed elsewhere in the algorithm.


    // Note: because the NQ variables do not include lambda and the
    // NGDNV do, we only support store_full_state = false for
    // Radiation

#ifdef RADIATION
    if (store_full_state == true) {
        amrex::Error("cannot store full interface state with radiation");
    }
#endif

    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();

    // do we want to force the flux to zero at the boundary?
    const bool special_bnd_lo = (lo_bc[idir] == Symmetry ||
                                 lo_bc[idir] == SlipWall ||
                                 lo_bc[idir] == NoSlipWall);
    const bool special_bnd_hi = (hi_bc[idir] == Symmetry ||
                                 hi_bc[idir] == SlipWall ||
                                 hi_bc[idir] == NoSlipWall);

    auto coord = geom.Coord();

    GeometryData geomdata = geom.data();

    const auto domlo = geom.Domain().loVect3d();
    const auto domhi = geom.Domain().hiVect3d();

#ifdef RADIATION
    int fspace_t = Radiation::fspace_advection_type;
    int comov = Radiation::comoving;
    int limiter = Radiation::limiter;
    int closure = Radiation::closure;
#endif

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {


        if (riemann_solver == 0 || riemann_solver == 1) {
            // approximate state Riemann solvers

            // first find the interface state on the current interface

            RiemannState qint;

            riemann_state(i, j, k, idir,
                          qm, qp, qaux_arr,
                          qint,
                          geomdata,
                          special_bnd_lo, special_bnd_hi,
                          domlo, domhi);

            // now use the interface state to compute and store the flux

            compute_flux_q(i, j, k, idir,
                           geomdata,
                           qint, flx,
#ifdef RADIATION
                           rflx,
                           fspace_t, comov,
                           limiter, closure,
#endif
                           qgdnv, store_full_state);

            // now do the passives -- we didn't include them in qint, qgdnv, or flux above

            // the passives are always just upwinded, so we do that here
            // regardless of the solver

            Real sgnm = std::copysign(1.0_rt, qint.un);
            if (qint.un == 0.0_rt) {
                sgnm = 0.0_rt;
            }

            Real fp = 0.5_rt*(1.0_rt + sgnm);
            Real fm = 0.5_rt*(1.0_rt - sgnm);

            // we don't want to store all the values of X, so we are going
            // to reconstruct in 2 steps to implement the CMA method

            Real sum_X = 0.0_rt;

            for (int n = 0; n < NumSpec; n++) {
                sum_X += amrex::max(0.0_rt,
                                    amrex::min(1.0_rt,
                                               fp * qm(i,j,k,QFS+n) + fm * qp(i,j,k,QFS+n)));
            }

            for (int ipassive = 0; ipassive < npassive; ipassive++) {
                int nqp = qpassmap(ipassive);
                int n  = upassmap(ipassive);

                Real X_int = fp * qm(i,j,k,nqp) + fm * qp(i,j,k,nqp);

                if (nqp >= QFS && nqp < QFS+NumSpec) {
                    // normalize
                    X_int = amrex::max(0.0_rt, amrex::min(1.0_rt, X_int));
                    X_int = X_int / sum_X;
                }

                flx(i,j,k,n) = flx(i,j,k,URHO) * X_int;

                if (store_full_state) {
                    qgdnv(i,j,k,nqp) = X_int;
                }
            }

        } else if (riemann_solver == 2) {
            // HLLC
            HLLC(i, j, k, idir,
                 qm, qp,
                 qaux_arr,
                 flx,
                 qgdnv, store_full_state,
                 geomdata,
                 special_bnd_lo, special_bnd_hi,
                 domlo, domhi);

#ifndef AMREX_USE_GPU
        } else {
            amrex::Error("ERROR: invalid value of riemann_solver");
#endif
        }

        if (hybrid_riemann == 1) {
            // correct the fluxes using an HLL scheme if we are in a shock
            // and doing the hybrid approach

            int is_shock = 0;

            if (idir == 0) {
                is_shock = static_cast<int>(shk(i-1,j,k) + shk(i,j,k));
            } else if (idir == 1) { 
                is_shock = static_cast<int>(shk(i,j-1,k) + shk(i,j,k));
            } else {
                is_shock = static_cast<int>(shk(i,j,k-1) + shk(i,j,k));
            }

            if (is_shock >= 1) {

                Real cl;
                Real cr;
                if (idir == 0) {
                    cl = qaux_arr(i-1,j,k,QC);
                    cr = qaux_arr(i,j,k,QC);
                } else if (idir == 1) { 
                    cl = qaux_arr(i,j-1,k,QC);
                    cr = qaux_arr(i,j,k,QC);
                } else {
                    cl = qaux_arr(i,j,k-1,QC);
                    cr = qaux_arr(i,j,k,QC);
                }

                Real ql_zone[NQ];
                Real qr_zone[NQ];
                Real flx_zone[NUM_STATE];

                for (int n = 0; n < NQ; n++) {
                    ql_zone[n] = qm(i,j,k,n);
                    qr_zone[n] = qp(i,j,k,n);
                }

                // pass in the current flux -- the
                // HLL solver will overwrite this
                // if necessary
                for (int n = 0; n < NUM_STATE; n++) {
                    flx_zone[n] = flx(i,j,k,n);
                }

                HLL(ql_zone, qr_zone, cl, cr,
                    idir, coord,
                    flx_zone);

                for (int n = 0; n < NUM_STATE; n++) {
                    flx(i,j,k,n) = flx_zone[n];
                }
            }
        }
    });

}


