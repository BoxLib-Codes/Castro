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

    auto coord = geom.Coord();

    GeometryData geomdata = geom.data();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {

        GpuArray<Real, NQ> qint_local;
#ifdef RADIATION
        GpuArrau<Real, NGROUPS> lambda_int_local;
#endif

        if (riemann_solver == 0 || riemann_solver == 1) {
            // approximate state Riemann solvers

            // first find the interface state on the current interface

            riemann_state(i, j, k, idir,
                          qm, qp, qaux_arr,
                          qint_local,
#ifdef RADIATION
                          lambda_int_local,
#endif
                          geomdata,
                          lo_bc, hi_bc);

            // now use the interface state to compute and store the flux

            compute_flux_q(i, j, k, idir,
                           geomdata,
                           qint_local, flx,
#ifdef RADIATION
                           lambda_int_local, rflx,
#endif
                           qgdnv, store_full_state);


        } else if (riemann_solver == 2) {
            // HLLC
            HLLC(i, j, k, idir,
                 qm, qp,
                 qaux_arr,
                 flx,
                 qgdnv, store_full_state,
                 geomdata,
                 lo_bc, hi_bc);

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


