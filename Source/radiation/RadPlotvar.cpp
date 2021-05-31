
#include <Radiation.H>
#include <fluxlimiter.H>

#include <RAD_F.H>

using namespace amrex;

void Radiation::save_lambda_in_plotvar(int level, const Array<MultiFab,AMREX_SPACEDIM>& lambda)
{
    int nlambda = lambda[0].nComp();
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*plotvar[level],true); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

        auto lamx = lambda[0][mfi].array();
#if AMREX_SPACEDIM >= 2
        auto lamy = lambda[1][mfi].array();
#endif
#if AMREX_SPACEDIM == 3
        auto lamz = lambda[2][mfi].array();
#endif

        auto lamc = (*plotvar[level])[mfi].array();

        int dcomp = icomp_lambda;

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {
            for (int n = 0; n < nlambda; ++n) {
#if AMREX_SPACEDIM == 1
                lamc(i,j,k,dcomp+n) = (lamx(i,j,k,n) + lamx(i+1,j,k,n)) * 0.5_rt;
#elif AMREX_SPACEDIM == 2
                lamc(i,j,k,dcomp+n) = (lamx(i,j,k,n) + lamx(i+1,j,k,n) +
                                       lamy(i,j,k,n) + lamy(i,j+1,k,n)) * 0.25_rt;
#else
                lamc(i,j,k,dcomp+n) = (lamx(i,j,k,n) + lamx(i+1,j,k,n) +
                                       lamy(i,j,k,n) + lamy(i,j+1,k,n) +
                                       lamz(i,j,k,n) + lamz(i,j,k+1,n)) * (1.0_rt / 6.0_rt);
#endif
            }
        });
    }
}

void Radiation::save_lab_Er_in_plotvar(int level, const MultiFab& Snew, 
                                       const MultiFab& Ecom, const MultiFab& F, int iflx)
{
    int nflx = F.nComp();

    int ifx = iflx;

#if AMREX_SPACEDIM >= 2
    int ify = iflx + NGROUPS;
#else
    int ify = 0;
#endif

#if AMREX_SPACEDIM == 3
    int ifz = iflx + NGROUPS * NGROUPS;
#else
    int ifz = 0;
#endif

    int ier = icomp_lab_Er;

    Real c2 = 1.e0_rt / (Radiation::clight * Radiation::clight);

    GpuArray<Real, NGROUPS> dlognu = {0.0};

    if (NGROUPS > 1) {
        ca_get_dlognu(dlognu.begin());
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*plotvar[level],true); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

        auto Snew_arr = Snew[mfi].array();
        auto Ecom_arr = Ecom[mfi].array();
        auto F_arr = F[mfi].array();
        auto Elab_arr = (*plotvar[level])[mfi].array();

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {
            Real rhoInv = 1.0_rt / Snew_arr(i,j,k,URHO);

            Real vxc2 = Snew_arr(i,j,k,UMX) * rhoInv * c2;

#if AMREX_SPACEDIM >= 2
            Real vyc2 = Snew_arr(i,j,k,UMY) * rhoInv * c2;
#else
            Real vyc2 = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
            Real vzc2 = Snew_arr(i,j,k,UMZ) * rhoInv * c2;
#else
            Real vzc2 = 0.0_rt;
#endif

            for (int g = 0; g < NGROUPS; ++g) {
                Elab_arr(i,j,k,g+ier) = Ecom_arr(i,j,k,g) + 2.0_rt * (vxc2 * F_arr(i,j,k,ifx+g) +
                                                                      vyc2 * F_arr(i,j,k,ify+g) +
                                                                      vzc2 * F_arr(i,j,k,ifz+g));
            }

            if (NGROUPS > 1) {
                Array1D<Real, -1, NGROUPS> nufnux, nufnuy, nufnuz;

                for (int g = 0; g < NGROUPS; ++g) {
                    nufnux(g) = F_arr(i,j,k,ifx+g) / dlognu[g];
                    nufnuy(g) = F_arr(i,j,k,ify+g) / dlognu[g];
                    nufnuz(g) = F_arr(i,j,k,ifz+g) / dlognu[g];
                }

                nufnux(-1) = -nufnux(0);
                nufnuy(-1) = -nufnuy(0);
                nufnuz(-1) = -nufnuz(0);
                nufnux(NGROUPS) = -nufnux(NGROUPS-1);
                nufnuy(NGROUPS) = -nufnuy(NGROUPS-1);
                nufnuz(NGROUPS) = -nufnuz(NGROUPS-1);

                for (int g = 0; g < NGROUPS; ++g) {
                    Elab_arr(i,j,k,g+ier) -= vxc2 * 0.5_rt * (nufnux(g+1) - nufnux(g-1)) +
                                             vyc2 * 0.5_rt * (nufnuy(g+1) - nufnuy(g-1)) +
                                             vzc2 * 0.5_rt * (nufnuz(g+1) - nufnuz(g-1));
                }
            }
        });
    }
}

void Radiation::save_flux_in_plotvar(int level, const MultiFab& Snew, 
                                     const Array<MultiFab,AMREX_SPACEDIM>& lambda,
                                     const MultiFab& Er, const MultiFab& Fr, int iflx,
                                     const Real lab_factor)
{
    int nflx = Fr.nComp();
    int nlambda = lambda[0].nComp();

    GpuArray<Real, NGROUPS> dlognu = {0.0};

    if (NGROUPS > 1) {
        ca_get_dlognu(dlognu.begin());
    }

    int limiter = Radiation::limiter;
    int closure = Radiation::closure;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*plotvar[level],true); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

        auto Snew_arr = Snew[mfi].array();
        auto Er_arr = Er[mfi].array();
        auto Fr_arr = Fr[mfi].array();
        auto Fo_arr = (*plotvar[level])[mfi].array();

        auto lamx = lambda[0][mfi].array();
#if AMREX_SPACEDIM >= 2
        auto lamy = lambda[1][mfi].array();
#endif
#if AMREX_SPACEDIM == 3
        auto lamz = lambda[2][mfi].array();
#endif

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {
            Array1D<Real, 0, NGROUPS-1> Eddf = {0.0};

            for (int g = 0; g < NGROUPS; ++g) {
                int ilam = amrex::min(g, nlambda - 1);

#if AMREX_SPACEDIM == 1
                Real lamcc = 0.5_rt * (lamx(i,j,k,ilam) + lamx(i+1,j,k,ilam));
#elif AMREX_SPACEDIM >= 2
                Real lamcc = 0.25_rt * (lamx(i,j,k,ilam) + lamx(i+1,j,k,ilam) +
                                        lamy(i,j,k,ilam) + lamy(i,j+1,k,ilam));
#else
                Real lamcc = (1.0_rt / 6.0_rt) * (lamx(i,j,k,ilam) + lamx(i+1,j,k,ilam) +
                                                  lamy(i,j,k,ilam) + lamy(i,j+1,k,ilam) +
                                                  lamz(i,j,k,ilam) + lamz(i,j,k+1,ilam));
#endif

                Eddf(g) = Edd_factor(lamcc, limiter, closure);
            }

            int ifix = iflx;
            int ifox = nflx;

#if AMREX_SPACEDIM >= 2
            int ifiy = iflx + NGROUPS;
            int ifoy = nflx + NGROUPS;
#endif

#if AMREX_SPACEDIM == 3
            int ifiz = iflx + NGROUPS * NGROUPS;
            int ifoz = nflx + NGROUPS * NGROUPS;
#endif

            Real rhoInv = 1.0_rt / Snew_arr(i,j,k,URHO);

            // lab_factor should be +1 for lab frame, -1 for comoving frame
            Real vx = Snew_arr(i,j,k,UMX) * rhoInv * lab_factor;
#if AMREX_SPACEDIM >= 2
            Real vy = Snew_arr(i,j,k,UMY) * rhoInv * lab_factor;
#endif
#if AMREX_SPACEDIM == 3
            Real vz = Snew_arr(i,j,k,UMZ) * rhoInv * lab_factor;
#endif

            Array1D<Real, 0, NGROUPS-1> vdotpx;
#if AMREX_SPACEDIM >= 2
            Array1D<Real, 0, NGROUPS-1> vdotpy;
#endif
#if AMREX_SPACEDIM == 3
            Array1D<Real, 0, NGROUPS-1> vdotpz;
#endif

            for (int g = 0; g < NGROUPS; ++g) {
                Real f1 = (1.0_rt - Eddf(g));
                Real f2 = (3.0_rt * Eddf(g) - 1.0_rt);
                Real foo = 1.0_rt / std::sqrt(Fr_arr(i,j,k,ifix+g) * Fr_arr(i,j,k,ifix+g) +
#if AMREX_SPACEDIM >= 2
                                              Fr_arr(i,j,k,ifiy+g) * Fr_arr(i,j,k,ifiy+g) +
#endif
#if AMREX_SPACEDIM == 3
                                              Fr_arr(i,j,k,ifiz+g) * Fr_arr(i,j,k,ifiz+g) +
#endif
                                              1.e-50_rt);

                Real nx = Fr_arr(i,j,k,ifix+g) * foo;
#if AMREX_SPACEDIM >= 2
                Real ny = Fr_arr(i,j,k,ifiy+g) * foo;
#endif
#if AMREX_SPACEDIM == 3
                Real nz = Fr_arr(i,j,k,ifiz+g) * foo;
#endif

                Real vdotn = vx * nx;
#if AMREX_SPACEDIM >= 2
                vdotn += vy * ny;
#endif
#if AMREX_SPACEDIM == 3
                vdotn += vz * vz;
#endif

                vdotpx(g) = 0.5_rt * Er_arr(i,j,k,g) * (f1 * vx + f2 * vdotn * nx);
#if AMREX_SPACEDIM >= 2
                vdotpy(g) = 0.5_rt * Er_arr(i,j,k,g) * (f1 * vy + f2 * vdotn * ny);
#endif
#if AMREX_SPACEDIM == 3
                vdotpz(g) = 0.5_rt * Er_arr(i,j,k,g) * (f1 * vz + f2 * vdotn * nz);
#endif

                Fo_arr(i,j,k,ifox+g) = Fr_arr(i,j,k,ifix+g) + vx * Er_arr(i,j,k,g) + vdotpx(g);
#if AMREX_SPACEDIM >= 2
                Fo_arr(i,j,k,ifoy+g) = Fr_arr(i,j,k,ifiy+g) + vy * Er_arr(i,j,k,g) + vdotpy(g);
#endif
#if AMREX_SPACEDIM == 3
                Fo_arr(i,j,k,ifoz+g) = Fr_arr(i,j,k,ifiz+g) + vz * Er_arr(i,j,k,g) + vdotpz(g);
#endif
            }

            if (NGROUPS > 1) {
                Array1D<Real, -1, NGROUPS> nuvpnux;
#if AMREX_SPACEDIM >= 2
                Array1D<Real, -1, NGROUPS> nuvpnuy;
#endif
#if AMREX_SPACEDIM == 3
                Array1D<Real, -1, NGROUPS> nuvpnuz;
#endif

                for (int g = 0; g < NGROUPS; ++g) {
                    nuvpnux(g) = vdotpx(g) / dlognu[g];
#if AMREX_SPACEDIM >= 2
                    nuvpnuy(g) = vdotpy(g) / dlognu[g];
#endif
#if AMREX_SPACEDIM == 3
                    nuvpnuz(g) = vdotpz(g) / dlognu[g];
#endif
                }

                nuvpnux(-1) = -nuvpnux(0);
                nuvpnux(NGROUPS) = -nuvpnux(NGROUPS-1);

#if AMREX_SPACEDIM >= 2
                nuvpnuy(-1) = -nuvpnuy(0);
                nuvpnuy(NGROUPS) = -nuvpnuy(NGROUPS-1);
#endif

#if AMREX_SPACEDIM == 3
                nuvpnuz(-1) = -nuvpnuz(0);
                nuvpnuz(NGROUPS) = -nuvpnuz(NGROUPS-1);
#endif

                for (int g = 0; g < NGROUPS; ++g) {
                    Fo_arr(i,j,k,ifox+g) -= 0.5_rt * (nuvpnux(g+1)-nuvpnux(g-1));
#if AMREX_SPACEDIM >= 2
                    Fo_arr(i,j,k,ifoy+g) -= 0.5_rt * (nuvpnuy(g+1)-nuvpnuy(g-1));
#endif
#if AMREX_SPACEDIM == 3
                    Fo_arr(i,j,k,ifoz+g) -= 0.5_rt * (nuvpnuz(g+1)-nuvpnuz(g-1));
#endif
                }
            }
        });
    }
}
