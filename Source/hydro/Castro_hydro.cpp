#include "Castro.H"
#include "Castro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

using namespace amrex;

#ifndef AMREX_USE_CUDA
void
Castro::construct_hydro_source(Real time, Real dt)
{

  BL_PROFILE("Castro::construct_hydro_source()");

  const Real strt_time = ParallelDescriptor::second();

  // this constructs the hydrodynamic source (essentially the flux
  // divergence) using the CTU framework for unsplit hydrodynamics

  if (verbose && ParallelDescriptor::IOProcessor())
    std::cout << "... Entering hydro advance" << std::endl << std::endl;

  hydro_source.setVal(0.0);

  int finest_level = parent->finestLevel();

  const Real *dx = geom.CellSize();

  const int* domain_lo = geom.Domain().loVect();
  const int* domain_hi = geom.Domain().hiVect();

  MultiFab& S_new = get_new_data(State_Type);

#ifdef RADIATION
  MultiFab& Er_new = get_new_data(Rad_Type);

  if (!Radiation::rad_hydro_combined) {
    amrex::Abort("Castro::construct_hydro_source -- we don't implement a mode where we have radiation, but it is not coupled to hydro");
  }

  int nstep_fsp = -1;
#endif

  // note: the radiation consup currently does not fill these
  Real mass_lost       = 0.;
  Real xmom_lost       = 0.;
  Real ymom_lost       = 0.;
  Real zmom_lost       = 0.;
  Real eden_lost       = 0.;
  Real xang_lost       = 0.;
  Real yang_lost       = 0.;
  Real zang_lost       = 0.;

  BL_PROFILE_VAR("Castro::advance_hydro_ca_umdrv()", CA_UMDRV);

#if AMREX_SPACEDIM >= 2
  MultiFab ftmp1;
  ftmp1.define(grids, dmap, NUM_STATE, 1);

  MultiFab ftmp2;
  ftmp2.define(grids, dmap, NUM_STATE, 1);

#ifdef RADIATION
  MultiFab rftmp1;
  rftmp1.define(grids, dmap, Radiation::nGroups, 1);

  MultiFab rftmp2;
  rftmp2.define(grids, dmap, Radiation::nGroups, 1);
#endif

  MultiFab qgdnvtmp1;
  qgdnvtmp1.define(grids, dmap, NGDNV, 2);

  MultiFab qgdnvtmp2;
  qgdnvtmp2.define(grids, dmap, NGDNV, 2);

  MultiFab ql;
  ql.define(grids, dmap, NQ, 2);

  MultiFab qr;
  qr.define(grids, dmap, NQ, 2);

#endif

#if AMREX_SPACEDIM == 3
  MultiFab qmxy;
  qmxy.define(grids, dmap, NQ, 2);

  MultiFab qpxy;
  qpxy.define(grids, dmap, NQ, 2);

  MultiFab qmxz;
  qmxz.define(grids, dmap, NQ, 2);

  MultiFab qpxz;
  qpxz.define(grids, dmap, NQ, 2);

  MultiFab qmyx;
  qmyx.define(grids, dmap, NQ, 2);

  MultiFab qpyx;
  qpyx.define(grids, dmap, NQ, 2);

  MultiFab qmyz;
  qmyz.define(grids, dmap, NQ, 2);

  MultiFab qpyz;
  qpyz.define(grids, dmap, NQ, 2);

  MultiFab qmzx;
  qmzx.define(grids, dmap, NQ, 2);

  MultiFab qpzx;
  qpzx.define(grids, dmap, NQ, 2);

  MultiFab qmzy;
  qmzy.define(grids, dmap, NQ, 2);

  MultiFab qpzy;
  qpzy.define(grids, dmap, NQ, 2);
#endif

#ifdef _OPENMP
#ifdef RADIATION
#pragma omp parallel reduction(max:nstep_fsp) \
                     reduction(+:mass_lost,xmom_lost,ymom_lost,zmom_lost) \
                     reduction(+:eden_lost,xang_lost,yang_lost,zang_lost)
#endif
#else
#pragma omp parallel
#ifdef _OPENMP
#pragma omp parallel reduction(+:mass_lost,xmom_lost,ymom_lost,zmom_lost) \
                     reduction(+:eden_lost,xang_lost,yang_lost,zang_lost)
#endif
#endif
  {

#ifdef RADIATION
    int priv_nstep_fsp = -1;
#endif

    for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

      // the valid region box
      const Box& bx = mfi.tilebox();

      const Box& obx = amrex::grow(bx, 1);

      AsyncFab flatn(obx, 1);

#ifdef RADIATION
      AsyncFab flatg(obx, 1);
#endif
      
      // compute the flattening coefficient

      if (first_order_hydro == 1) {
        flatn.hostFab().setVal(0.0, obx);
      } else if (use_flattening == 1) {
#ifdef RADIATION
        ca_rad_flatten(ARLIM_3D(obx.loVect()), ARLIM_3D(obx.hiVect()),
                       BL_TO_FORTRAN_ANYD(q[mfi]),
                       BL_TO_FORTRAN_ANYD(flatn.hostFab()),
                       BL_TO_FORTRAN_ANYD(flatg.hostFab()));
#else
        ca_uflatten(ARLIM_3D(obx.loVect()), ARLIM_3D(obx.hiVect()),
                    BL_TO_FORTRAN_ANYD(q[mfi]),
                    BL_TO_FORTRAN_ANYD(flatn.hostFab()), QPRES+1);
#endif
      } else {
        flatn.hostFab().setVal(1.0, obx);
      }

#ifdef RADIATION
      flatg.clear();
#endif

      const Box& xbx = amrex::surroundingNodes(bx, 0);
      const Box& ybx = amrex::surroundingNodes(bx, 1);
      const Box& zbx = amrex::surroundingNodes(bx, 2);

      AsyncFab dq(obx, AMREX_SPACEDIM*NQ);
      AsyncFab Ip(obx, AMREX_SPACEDIM*3*NQ);
      AsyncFab Im(obx, AMREX_SPACEDIM*3*NQ);
      AsyncFab Ip_src(obx, AMREX_SPACEDIM*3*QVAR);
      AsyncFab Im_src(obx, AMREX_SPACEDIM*3*QVAR);
      AsyncFab Ip_gc(obx, AMREX_SPACEDIM*3);
      AsyncFab Im_gc(obx,AMREX_SPACEDIM*3);
      AsyncFab sm(obx, AMREX_SPACEDIM);
      AsyncFab sp(obx, AMREX_SPACEDIM);
      AsyncFab shk(obx, 1);
      AsyncFab qxm(xbx, NQ);
      AsyncFab qxp(xbx, NQ);
#if AMREX_SPACEDIM >= 2
      AsyncFab qym(ybx, NQ);
      AsyncFab qyp(ybx, NQ);
#endif
#if AMREX_SPACEDIM == 3
      AsyncFab qzm(zbx, NQ);
      AsyncFab qzp(zbx, NQ);
#endif

      ctu_normal_states(ARLIM_3D(obx.loVect()), ARLIM_3D(obx.hiVect()),
                        ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                        BL_TO_FORTRAN_ANYD(q[mfi]),
                        BL_TO_FORTRAN_ANYD(flatn.hostFab()),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(src_q[mfi]),
                        BL_TO_FORTRAN_ANYD(shk.hostFab()),
                        BL_TO_FORTRAN_ANYD(Ip.hostFab()),
                        BL_TO_FORTRAN_ANYD(Im.hostFab()),
                        BL_TO_FORTRAN_ANYD(Ip_src.hostFab()),
                        BL_TO_FORTRAN_ANYD(Im_src.hostFab()),
                        BL_TO_FORTRAN_ANYD(Ip_gc.hostFab()),
                        BL_TO_FORTRAN_ANYD(Im_gc.hostFab()),
                        BL_TO_FORTRAN_ANYD(dq.hostFab()),
                        BL_TO_FORTRAN_ANYD(sm.hostFab()),
                        BL_TO_FORTRAN_ANYD(sp.hostFab()),
                        BL_TO_FORTRAN_ANYD(qxm.hostFab()),
                        BL_TO_FORTRAN_ANYD(qxp.hostFab()),
#if AMREX_SPACEDIM >= 2
                        BL_TO_FORTRAN_ANYD(qym.hostFab()),
                        BL_TO_FORTRAN_ANYD(qyp.hostFab()),
#endif
#if AMREX_SPACEDIM == 3
                        BL_TO_FORTRAN_ANYD(qzm.hostFab()),
                        BL_TO_FORTRAN_ANYD(qzp.hostFab()),
#endif
                        ZFILL(dx), dt,
#if (AMREX_SPACEDIM < 3)
                        BL_TO_FORTRAN_ANYD(dLogArea[0][mfi]),
#endif
                        ARLIM_3D(domain_lo), ARLIM_3D(domain_hi));

      // Clear local Fabs
      dq.clear();
      Ip.clear();
      Im.clear();
      Ip_src.clear();
      Im_src.clear();
      Ip_gc.clear();
      Im_gc.clear();
      sm.clear();
      sp.clear();
      flatn.clear();

      AsyncFab div(obx, 1);
      
      // compute divu -- we'll use this later when doing the artifical viscosity
      divu(ARLIM_3D(obx.loVect()), ARLIM_3D(obx.hiVect()),
           BL_TO_FORTRAN_ANYD(q[mfi]),
           ZFILL(dx),
           BL_TO_FORTRAN_ANYD(div.hostFab()));

      const Box& tbx = amrex::grow(bx, 2);
      

      AsyncFab q_int(tbx, NQ);
#ifdef RADIATION
      AsyncFab lambda_int(tbx, Radiation::nGroups);
#endif

      Array<AsyncFab, AMREX_SPACEDIM> flux{AMREX_D_DECL(AsyncFab(xbx, NUM_STATE),
                                                        AsyncFab(ybx, NUM_STATE),
                                                        AsyncFab(zbx, NUM_STATE))};
      Array<AsyncFab, AMREX_SPACEDIM> qe{AMREX_D_DECL(AsyncFab(xbx, NGDNV),
                                                      AsyncFab(ybx, NGDNV),
                                                      AsyncFab(zbx, NGDNV))};
#ifdef RADIATION
      Array<AsyncFab, AMREX_SPACEDIM> rad_flux{AMREX_D_DECL(AsyncFab(xbx, Radiation::nGroups),
                                                            AsyncFab(ybx, Radiation::nGroups),
                                                            AsyncFab(zbx, Radiation::nGroups))};
#endif

#if AMREX_SPACEDIM <= 2
      AsyncFab pradial(xbx, 1);
#endif

#if AMREX_SPACEDIM == 1
      const Box& nxbx = mfi.nodaltilebox(0);

      cmpflx_plus_godunov(ARLIM_3D(nxbx.loVect()), ARLIM_3D(nxbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qxm.hostFab()),
                          BL_TO_FORTRAN_ANYD(qxp.hostFab()), 1, 1,
                          BL_TO_FORTRAN_ANYD(flux[0].hostFab()),
                          BL_TO_FORTRAN_ANYD(q_int.hostFab()),
  #ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rad_flux[0].hostFab()),
                          BL_TO_FORTRAN_ANYD(lambda_int.hostFab()),
  #endif
                          BL_TO_FORTRAN_ANYD(qe[0].hostFab()),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk.hostFab()),
                          1, ARLIM_3D(domain_lo), ARLIM_3D(domain_hi));

#endif // 1-d



#if AMREX_SPACEDIM == 2

      const amrex::Real hdt = 0.5*dt;
      const amrex::Real hdtdx = 0.5*dt/dx[0];
      const amrex::Real hdtdy = 0.5*dt/dx[1];

      // compute F^x
      // [lo(1), lo(2)-1, 0], [hi(1)+1, hi(2)+1, 0]
      const Box& cxbx = mfi.grownnodaltilebox(0, IntVect(AMREX_D_DECL(0,1,0)));

      // ftmp1 = fx
      // rftmp1 = rfx
      // qgdnvtmp1 = qgdnxv
      cmpflx_plus_godunov(ARLIM_3D(cxbx.loVect()), ARLIM_3D(cxbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qxm.hostFab()),
                          BL_TO_FORTRAN_ANYD(qxp.hostFab()), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp1[mfi]),
                          BL_TO_FORTRAN_ANYD(q_int.hostFab()),
  #ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp1[mfi]),
                          BL_TO_FORTRAN_ANYD(lambda_int.hostFab()),
  #endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp1[mfi]),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk.hostFab()),
                          1, ARLIM_3D(domain_lo), ARLIM_3D(domain_hi));

      // compute F^y
      // [lo(1)-1, lo(2), 0], [hi(1)+1, hi(2)+1, 0]
      const Box& cybx = mfi.grownnodaltilebox(1,IntVect(AMREX_D_DECL(1,0,0)));

      // ftmp2 = fy
      // rftmp2 = rfy
      cmpflx_plus_godunov(ARLIM_3D(cybx.loVect()), ARLIM_3D(cybx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qym.hostFab()),
                          BL_TO_FORTRAN_ANYD(qyp.hostFab()), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp2[mfi]),
                          BL_TO_FORTRAN_ANYD(q_int.hostFab()),
  #ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp2[mfi]),
                          BL_TO_FORTRAN_ANYD(lambda_int.hostFab()),
  #endif
                          BL_TO_FORTRAN_ANYD(qe[1].hostFab()),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk.hostFab()),
                          2, ARLIM_3D(domain_lo), ARLIM_3D(domain_hi));

      // add the transverse flux difference in y to the x states

      // [lo(1), lo(2), 0], [hi(1)+1, hi(2), 0]
      const Box& nxbx = mfi.nodaltilebox(0);

      // ftmp2 = fy
      // rftmp2 = rfy
      transy_on_xstates(ARLIM_3D(nxbx.loVect()), ARLIM_3D(nxbx.hiVect()),
                        BL_TO_FORTRAN_ANYD(qxm.hostFab()),
                        BL_TO_FORTRAN_ANYD(ql[mfi]),
                        BL_TO_FORTRAN_ANYD(qxp.hostFab()),
                        BL_TO_FORTRAN_ANYD(qr[mfi]),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp2[mfi]),
  #ifdef RADIATION
                        BL_TO_FORTRAN_ANYD(rftmp2[mfi]),
  #endif
                        BL_TO_FORTRAN_ANYD(qe[1].hostFab()),
                        hdtdy);

      // solve the final Riemann problem axross the x-interfaces

      cmpflx_plus_godunov(ARLIM_3D(nxbx.loVect()), ARLIM_3D(nxbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(ql[mfi]),
                          BL_TO_FORTRAN_ANYD(qr[mfi]), 1, 1,
                          BL_TO_FORTRAN_ANYD(flux[0].hostFab()),
                          BL_TO_FORTRAN_ANYD(q_int.hostFab()),
  #ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rad_flux[0].hostFab()),
                          BL_TO_FORTRAN_ANYD(lambda_int.hostFab()),
  #endif
                          BL_TO_FORTRAN_ANYD(qe[0].hostFab()),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk.hostFab()),
                          1, ARLIM_3D(domain_lo), ARLIM_3D(domain_hi));

      // add the transverse flux difference in x to the y states
      // [lo(1), lo(2)-1, 0], [hi(1), hi(2)+1, 0]
      const Box& nybx = mfi.nodaltilebox(1);

      // ftmp1 = fx
      // rftmp1 = rfx
      // qgdnvtmp1 = qgdnvx

      transx_on_ystates(ARLIM_3D(nybx.loVect()), ARLIM_3D(nybx.hiVect()),
                        BL_TO_FORTRAN_ANYD(qym.hostFab()),
                        BL_TO_FORTRAN_ANYD(ql[mfi]),
                        BL_TO_FORTRAN_ANYD(qyp.hostFab()),
                        BL_TO_FORTRAN_ANYD(qr[mfi]),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp1[mfi]),
  #ifdef RADIATION
                        BL_TO_FORTRAN_ANYD(rftmp1[mfi]),
  #endif
                        BL_TO_FORTRAN_ANYD(qgdnvtmp1[mfi]),
                        BL_TO_FORTRAN_ANYD(area[0][mfi]),
                        BL_TO_FORTRAN_ANYD(volume[mfi]),
                        hdt, hdtdx);

      // solve the final Riemann problem axross the y-interfaces

      cmpflx_plus_godunov(ARLIM_3D(nybx.loVect()), ARLIM_3D(nybx.hiVect()),
                          BL_TO_FORTRAN_ANYD(ql[mfi]),
                          BL_TO_FORTRAN_ANYD(qr[mfi]), 1, 1,
                          BL_TO_FORTRAN_ANYD(flux[1].hostFab()),
                          BL_TO_FORTRAN_ANYD(q_int.hostFab()),
  #ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rad_flux[1].hostFab()),
                          BL_TO_FORTRAN_ANYD(lambda_int.hostFab()),
  #endif
                          BL_TO_FORTRAN_ANYD(qe[1].hostFab()),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk.hostFab()),
                          2, ARLIM_3D(domain_lo), ARLIM_3D(domain_hi));
  #endif // 2-d



  #if AMREX_SPACEDIM == 3

      const amrex::Real hdt = 0.5*dt;

      const amrex::Real hdtdx = 0.5*dt/dx[0];
      const amrex::Real hdtdy = 0.5*dt/dx[1];
      const amrex::Real hdtdz = 0.5*dt/dx[2];

      const amrex::Real cdtdx = dt/dx[0]/3.0;
      const amrex::Real cdtdy = dt/dx[1]/3.0;
      const amrex::Real cdtdz = dt/dx[2]/3.0;

      // compute F^x
      // [lo(1), lo(2)-1, lo(3)-1], [hi(1)+1, hi(2)+1, hi(3)+1]
      const Box& cxbx = mfi.grownnodaltilebox(0, IntVect(AMREX_D_DECL(0,1,1)));

      // ftmp1 = fx
      // rftmp1 = rfx
      // qgdnvtmp1 = qgdnxv
      cmpflx_plus_godunov(ARLIM_3D(cxbx.loVect()), ARLIM_3D(cxbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qxm.hostFab()),
                          BL_TO_FORTRAN_ANYD(qxp.hostFab()), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp1[mfi]),
                          BL_TO_FORTRAN_ANYD(q_int.hostFab()),
  #ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp1[mfi]),
                          BL_TO_FORTRAN_ANYD(lambda_int.hostFab()),
  #endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp1[mfi]),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk.hostFab()),
                          1, ARLIM_3D(domain_lo), ARLIM_3D(domain_hi));

      // [lo(1), lo(2), lo(3)-1], [hi(1), hi(2)+1, hi(3)+1]
      const Box& txybx = mfi.grownnodaltilebox(1, IntVect(AMREX_D_DECL(0,0,1)));

      // ftmp1 = fx
      // rftmp1 = rfx
      // qgdnvtmp1 = qgdnvx
      transx_on_ystates(ARLIM_3D(txybx.loVect()), ARLIM_3D(txybx.hiVect()),
                        BL_TO_FORTRAN_ANYD(qym.hostFab()),
                        BL_TO_FORTRAN_ANYD(qmyx[mfi]),
                        BL_TO_FORTRAN_ANYD(qyp.hostFab()),
                        BL_TO_FORTRAN_ANYD(qpyx[mfi]),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp1[mfi]),
  #ifdef RADIATION
                        BL_TO_FORTRAN_ANYD(rftmp1[mfi]),
  #endif
                        BL_TO_FORTRAN_ANYD(qgdnvtmp1[mfi]),
                        hdt, cdtdx);

      // [lo(1), lo(2)-1, lo(3)], [hi(1), hi(2)+1, hi(3)+1]
      const Box& txzbx = mfi.grownnodaltilebox(2, IntVect(AMREX_D_DECL(0,1,0)));

      transx_on_zstates(ARLIM_3D(txzbx.loVect()), ARLIM_3D(txzbx.hiVect()),
                        BL_TO_FORTRAN_ANYD(qzm.hostFab()),
                        BL_TO_FORTRAN_ANYD(qmzx[mfi]),
                        BL_TO_FORTRAN_ANYD(qzp.hostFab()),
                        BL_TO_FORTRAN_ANYD(qpzx[mfi]),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp1[mfi]),
  #ifdef RADIATION
                        BL_TO_FORTRAN_ANYD(rftmp1[mfi]),
  #endif
                        BL_TO_FORTRAN_ANYD(qgdnvtmp1[mfi]),
                        hdt, cdtdx);

      // compute F^y
      // [lo(1)-1, lo(2), lo(3)-1], [hi(1)+1, hi(2)+1, hi(3)+1]
      const Box& cybx = mfi.grownnodaltilebox(1,IntVect(AMREX_D_DECL(1,0,1)));

      // ftmp1 = fy
      // rftmp1 = rfy
      // qgdnvtmp1 = qgdnvy
      cmpflx_plus_godunov(ARLIM_3D(cybx.loVect()), ARLIM_3D(cybx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qym.hostFab()),
                          BL_TO_FORTRAN_ANYD(qyp.hostFab()), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp1[mfi]),
                          BL_TO_FORTRAN_ANYD(q_int.hostFab()),
  #ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp1[mfi]),
                          BL_TO_FORTRAN_ANYD(lambda_int.hostFab()),
  #endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp1[mfi]),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk.hostFab()),
                          2, ARLIM_3D(domain_lo), ARLIM_3D(domain_hi));

      // [lo(1), lo(2), lo(3)-1], [hi(1)+1, hi(2), lo(3)+1]
      const Box& tyxbx = mfi.grownnodaltilebox(0, IntVect(AMREX_D_DECL(0,0,1)));

      // ftmp1 = fy
      // rftmp1 = rfy
      // qgdnvtmp1 = qgdnvy
      transy_on_xstates(ARLIM_3D(tyxbx.loVect()), ARLIM_3D(tyxbx.hiVect()),
                        BL_TO_FORTRAN_ANYD(qxm.hostFab()),
                        BL_TO_FORTRAN_ANYD(qmxy[mfi]),
                        BL_TO_FORTRAN_ANYD(qxp.hostFab()),
                        BL_TO_FORTRAN_ANYD(qpxy[mfi]),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp1[mfi]),
  #ifdef RADIATION
                        BL_TO_FORTRAN_ANYD(rftmp1[mfi]),
  #endif
                        BL_TO_FORTRAN_ANYD(qgdnvtmp1[mfi]),
                        cdtdy);

      // [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2), lo(3)+1]
      const Box& tyzbx = mfi.grownnodaltilebox(2, IntVect(AMREX_D_DECL(1,0,0)));

      // ftmp1 = fy
      // rftmp1 = rfy
      // qgdnvtmp1 = qgdnvy
      transy_on_zstates(ARLIM_3D(tyzbx.loVect()), ARLIM_3D(tyzbx.hiVect()),
                        BL_TO_FORTRAN_ANYD(qzm.hostFab()),
                        BL_TO_FORTRAN_ANYD(qmzy[mfi]),
                        BL_TO_FORTRAN_ANYD(qzp.hostFab()),
                        BL_TO_FORTRAN_ANYD(qpzy[mfi]),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp1[mfi]),
  #ifdef RADIATION
                        BL_TO_FORTRAN_ANYD(rftmp1[mfi]),
  #endif
                        BL_TO_FORTRAN_ANYD(qgdnvtmp1[mfi]),
                        cdtdy);

      // compute F^z
      // [lo(1)-1, lo(2)-1, lo(3)], [hi(1)+1, hi(2)+1, hi(3)+1]
      const Box& czbx = mfi.grownnodaltilebox(2, IntVect(AMREX_D_DECL(1,1,0)));

      // ftmp1 = fz
      // rftmp1 = rfz
      // qgdnvtmp1 = qgdnvz
      cmpflx_plus_godunov(ARLIM_3D(czbx.loVect()), ARLIM_3D(czbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qzm.hostFab()),
                          BL_TO_FORTRAN_ANYD(qzp.hostFab()), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp1[mfi]),
                          BL_TO_FORTRAN_ANYD(q_int.hostFab()),
  #ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp1[mfi]),
                          BL_TO_FORTRAN_ANYD(lambda_int.hostFab()),
  #endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp1[mfi]),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk.hostFab()),
                          3, ARLIM_3D(domain_lo), ARLIM_3D(domain_hi));

      // [lo(1)-1, lo(2)-1, lo(3)], [hi(1)+1, hi(2)+1, lo(3)]
      const Box& tzxbx = mfi.grownnodaltilebox(0, IntVect(AMREX_D_DECL(0,1,0)));

      // ftmp1 = fz
      // rftmp1 = rfz
      // qgdnvtmp1 = qgdnvz
      transz_on_xstates(ARLIM_3D(tzxbx.loVect()), ARLIM_3D(tzxbx.hiVect()),
                        BL_TO_FORTRAN_ANYD(qxm.hostFab()),
                        BL_TO_FORTRAN_ANYD(qmxz[mfi]),
                        BL_TO_FORTRAN_ANYD(qxp.hostFab()),
                        BL_TO_FORTRAN_ANYD(qpxz[mfi]),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp1[mfi]),
  #ifdef RADIATION
                        BL_TO_FORTRAN_ANYD(rftmp1[mfi]),
  #endif
                        BL_TO_FORTRAN_ANYD(qgdnvtmp1[mfi]),
                        cdtdz);

      // [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2)+1, lo(3)]
      const Box& tzybx = mfi.grownnodaltilebox(1, IntVect(AMREX_D_DECL(1,0,0)));

      // ftmp1 = fz
      // rftmp1 = rfz
      // qgdnvtmp1 = qgdnvz
      transz_on_ystates(ARLIM_3D(tzybx.loVect()), ARLIM_3D(tzybx.hiVect()),
                        BL_TO_FORTRAN_ANYD(qym.hostFab()),
                        BL_TO_FORTRAN_ANYD(qmyz[mfi]),
                        BL_TO_FORTRAN_ANYD(qyp.hostFab()),
                        BL_TO_FORTRAN_ANYD(qpyz[mfi]),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp1[mfi]),
  #ifdef RADIATION
                        BL_TO_FORTRAN_ANYD(rftmp1[mfi]),
  #endif
                        BL_TO_FORTRAN_ANYD(qgdnvtmp1[mfi]),
                        cdtdz);

      // we now have q?zx, q?yx, q?zy, q?xy, q?yz, q?xz

      //
      // Use qx?, q?yz, q?zy to compute final x-flux
      //

      // compute F^{y|z}
      // [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2)+1, hi(3)]
      const Box& cyzbx = mfi.grownnodaltilebox(1, IntVect(AMREX_D_DECL(1,0,0)));

      // ftmp1 = fyz
      // rftmp1 = rfyz
      // qgdnvtmp1 = qgdnvyz
      cmpflx_plus_godunov(ARLIM_3D(cyzbx.loVect()), ARLIM_3D(cyzbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qmyz[mfi]),
                          BL_TO_FORTRAN_ANYD(qpyz[mfi]), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp1[mfi]),
                          BL_TO_FORTRAN_ANYD(q_int.hostFab()),
  #ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp1[mfi]),
                          BL_TO_FORTRAN_ANYD(lambda_int.hostFab()),
  #endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp1[mfi]),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk.hostFab()),
                          2, ARLIM_3D(domain_lo), ARLIM_3D(domain_hi));

      // compute F^{z|y}
      // [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)+1]
      const Box& czybx = mfi.grownnodaltilebox(2, IntVect(AMREX_D_DECL(1,0,0)));

      // ftmp2 = fzy
      // rftmp2 = rfzy
      // qgdnvtmp2 = qgdnvzy
      cmpflx_plus_godunov(ARLIM_3D(czybx.loVect()), ARLIM_3D(czybx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qmzy[mfi]),
                          BL_TO_FORTRAN_ANYD(qpzy[mfi]), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp2[mfi]),
                          BL_TO_FORTRAN_ANYD(q_int.hostFab()),
  #ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp2[mfi]),
                          BL_TO_FORTRAN_ANYD(lambda_int.hostFab()),
  #endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp2[mfi]),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk.hostFab()),
                          3, ARLIM_3D(domain_lo), ARLIM_3D(domain_hi));

      // compute the corrected x interface states and fluxes
      // [lo(1), lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)]
      const Box& fcxbx = mfi.nodaltilebox(0);

      transyz(ARLIM_3D(fcxbx.loVect()), ARLIM_3D(fcxbx.hiVect()),
              BL_TO_FORTRAN_ANYD(qxm.hostFab()),
              BL_TO_FORTRAN_ANYD(ql[mfi]),
              BL_TO_FORTRAN_ANYD(qxp.hostFab()),
              BL_TO_FORTRAN_ANYD(qr[mfi]),
              BL_TO_FORTRAN_ANYD(qaux[mfi]),
              BL_TO_FORTRAN_ANYD(ftmp1[mfi]),
  #ifdef RADIATION
              BL_TO_FORTRAN_ANYD(rftmp1[mfi]),
  #endif
              BL_TO_FORTRAN_ANYD(ftmp2[mfi]),
  #ifdef RADIATION
              BL_TO_FORTRAN_ANYD(rftmp2[mfi]),
  #endif
              BL_TO_FORTRAN_ANYD(qgdnvtmp1[mfi]),
              BL_TO_FORTRAN_ANYD(qgdnvtmp2[mfi]),
              hdt, hdtdy, hdtdz);

      cmpflx_plus_godunov(ARLIM_3D(cxbx.loVect()), ARLIM_3D(cxbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(ql[mfi]),
                          BL_TO_FORTRAN_ANYD(qr[mfi]), 1, 1,
                          BL_TO_FORTRAN_ANYD(flux[0].hostFab()),
                          BL_TO_FORTRAN_ANYD(q_int.hostFab()),
  #ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rad_flux[0].hostFab()),
                          BL_TO_FORTRAN_ANYD(lambda_int.hostFab()),
  #endif
                          BL_TO_FORTRAN_ANYD(qe[0].hostFab()),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk.hostFab()),
                          1, ARLIM_3D(domain_lo), ARLIM_3D(domain_hi));

      //
      // Use qy?, q?zx, q?xz to compute final y-flux
      //

      // compute F^{z|x}
      // [lo(1), lo(2)-1, lo(3)], [hi(1), hi(2)+1, hi(3)+1]
      const Box& czxbx = mfi.grownnodaltilebox(2, IntVect(AMREX_D_DECL(0,1,0)));

      // ftmp1 = fzx
      // rftmp1 = rfzx
      // qgdnvtmp1 = qgdnvzx
      cmpflx_plus_godunov(ARLIM_3D(czxbx.loVect()), ARLIM_3D(czxbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qmzx[mfi]),
                          BL_TO_FORTRAN_ANYD(qpzx[mfi]), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp1[mfi]),
                          BL_TO_FORTRAN_ANYD(q_int.hostFab()),
  #ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp1[mfi]),
                          BL_TO_FORTRAN_ANYD(lambda_int.hostFab()),
  #endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp1[mfi]),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk.hostFab()),
                          3, ARLIM_3D(domain_lo), ARLIM_3D(domain_hi));

      // compute F^{x|z}
      // [lo(1), lo(2)-1, lo(3)], [hi(1)+1, hi(2)+1, hi(3)]
      const Box& cxzbx = mfi.grownnodaltilebox(0, IntVect(AMREX_D_DECL(0,1,0)));

      // ftmp2 = fxz
      // rftmp2 = rfxz
      // qgdnvtmp2 = qgdnvxz
      cmpflx_plus_godunov(ARLIM_3D(cxzbx.loVect()), ARLIM_3D(cxzbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qmxz[mfi]),
                          BL_TO_FORTRAN_ANYD(qpxz[mfi]), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp2[mfi]),
                          BL_TO_FORTRAN_ANYD(q_int.hostFab()),
  #ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp2[mfi]),
                          BL_TO_FORTRAN_ANYD(lambda_int.hostFab()),
  #endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp2[mfi]),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk.hostFab()),
                          1, ARLIM_3D(domain_lo), ARLIM_3D(domain_hi));

      // Compute the corrected y interface states and fluxes
      // [lo(1), lo(2), lo(3)], [hi(1), hi(2)+1, hi(3)]
      const Box& fcybx = mfi.nodaltilebox(1);

      transxz(ARLIM_3D(fcybx.loVect()), ARLIM_3D(fcybx.hiVect()),
              BL_TO_FORTRAN_ANYD(qym.hostFab()),
              BL_TO_FORTRAN_ANYD(ql[mfi]),
              BL_TO_FORTRAN_ANYD(qyp.hostFab()),
              BL_TO_FORTRAN_ANYD(qr[mfi]),
              BL_TO_FORTRAN_ANYD(qaux[mfi]),
              BL_TO_FORTRAN_ANYD(ftmp2[mfi]),
  #ifdef RADIATION
              BL_TO_FORTRAN_ANYD(rftmp2[mfi]),
  #endif
              BL_TO_FORTRAN_ANYD(ftmp1[mfi]),
  #ifdef RADIATION
              BL_TO_FORTRAN_ANYD(rftmp1[mfi]),
  #endif
              BL_TO_FORTRAN_ANYD(qgdnvtmp2[mfi]),
              BL_TO_FORTRAN_ANYD(qgdnvtmp1[mfi]),
              hdt, hdtdx, hdtdz);

      // Compute the final F^y
      // [lo(1), lo(2), lo(3)], [hi(1), hi(2)+1, hi(3)]
      cmpflx_plus_godunov(ARLIM_3D(cybx.loVect()), ARLIM_3D(cybx.hiVect()),
                          BL_TO_FORTRAN_ANYD(ql[mfi]),
                          BL_TO_FORTRAN_ANYD(qr[mfi]), 1, 1,
                          BL_TO_FORTRAN_ANYD(flux[1].hostFab()),
                          BL_TO_FORTRAN_ANYD(q_int.hostFab()),
  #ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rad_flux[1].hostFab()),
                          BL_TO_FORTRAN_ANYD(lambda_int.hostFab()),
  #endif
                          BL_TO_FORTRAN_ANYD(qe[1].hostFab()),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk.hostFab()),
                          2, ARLIM_3D(domain_lo), ARLIM_3D(domain_hi));

      //
      // Use qz?, q?xy, q?yx to compute final z-flux
      //

      // compute F^{x|y}
      // [lo(1), lo(2), lo(3)-1], [hi(1)+1, hi(2), hi(3)+1]
      const Box& cxybx = mfi.grownnodaltilebox(0, IntVect(AMREX_D_DECL(0,0,1)));

      // ftmp1 = fxy
      // rftmp1 = rfxy
      // qgdnvtmp1 = qgdnvxy
      cmpflx_plus_godunov(ARLIM_3D(cxybx.loVect()), ARLIM_3D(cxybx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qmxy[mfi]),
                          BL_TO_FORTRAN_ANYD(qpxy[mfi]), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp1[mfi]),
                          BL_TO_FORTRAN_ANYD(q_int.hostFab()),
  #ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp1[mfi]),
                          BL_TO_FORTRAN_ANYD(lambda_int.hostFab()),
  #endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp1[mfi]),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk.hostFab()),
                          1, ARLIM_3D(domain_lo), ARLIM_3D(domain_hi));

      // compute F^{y|x}
      // [lo(1), lo(2), lo(3)-1], [hi(1), hi(2)+dg(2), hi(3)+1]
      const Box& cyxbx = mfi.grownnodaltilebox(1, IntVect(AMREX_D_DECL(0,0,1)));

      // ftmp2 = fyx
      // rftmp2 = rfyx
      // qgdnvtmp2 = qgdnvyx
      cmpflx_plus_godunov(ARLIM_3D(cyxbx.loVect()), ARLIM_3D(cyxbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qmyx[mfi]),
                          BL_TO_FORTRAN_ANYD(qpyx[mfi]), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp2[mfi]),
                          BL_TO_FORTRAN_ANYD(q_int.hostFab()),
  #ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp2[mfi]),
                          BL_TO_FORTRAN_ANYD(lambda_int.hostFab()),
  #endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp2[mfi]),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk.hostFab()),
                          2, ARLIM_3D(domain_lo), ARLIM_3D(domain_hi));

      // compute the corrected z interface states and fluxes
      // [lo(1), lo(2), lo(3)], [hi(1), hi(2), hi(3)+1]
      const Box& fczbx =  mfi.nodaltilebox(2);

      transxy(ARLIM_3D(fczbx.loVect()), ARLIM_3D(fczbx.hiVect()),
              BL_TO_FORTRAN_ANYD(qzm.hostFab()),
              BL_TO_FORTRAN_ANYD(ql[mfi]),
              BL_TO_FORTRAN_ANYD(qzp.hostFab()),
              BL_TO_FORTRAN_ANYD(qr[mfi]),
              BL_TO_FORTRAN_ANYD(qaux[mfi]),
              BL_TO_FORTRAN_ANYD(ftmp1[mfi]),
  #ifdef RADIATION
              BL_TO_FORTRAN_ANYD(rftmp1[mfi]),
  #endif
              BL_TO_FORTRAN_ANYD(ftmp2[mfi]),
  #ifdef RADIATION
              BL_TO_FORTRAN_ANYD(rftmp2[mfi]),
  #endif
              BL_TO_FORTRAN_ANYD(qgdnvtmp1[mfi]),
              BL_TO_FORTRAN_ANYD(qgdnvtmp2[mfi]),
              hdt, hdtdx, hdtdy);

      // compute the final z fluxes F^z
      // [lo(1), lo(2), lo(3)], [hi(1), hi(2), hi(3)+1]

      cmpflx_plus_godunov(ARLIM_3D(czbx.loVect()), ARLIM_3D(czbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(ql[mfi]),
                          BL_TO_FORTRAN_ANYD(qr[mfi]), 1, 1,
                          BL_TO_FORTRAN_ANYD(flux[2].hostFab()),
                          BL_TO_FORTRAN_ANYD(q_int.hostFab()),
  #ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rad_flux[2].hostFab()),
                          BL_TO_FORTRAN_ANYD(lambda_int.hostFab()),
  #endif
                          BL_TO_FORTRAN_ANYD(qe[2].hostFab()),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk.hostFab()),
                          3, ARLIM_3D(domain_lo), ARLIM_3D(domain_hi));

  #endif // 3-d



      // clean the fluxes

      for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

          const Box& nbx = mfi.nodaltilebox(idir);

          int idir_f = idir + 1;

          ctu_clean_fluxes(ARLIM_3D(nbx.loVect()), ARLIM_3D(nbx.hiVect()),
                           idir_f,
                           BL_TO_FORTRAN_ANYD(Sborder[mfi]),
                           BL_TO_FORTRAN_ANYD(q[mfi]),
                           BL_TO_FORTRAN_ANYD(flux[idir].hostFab()),
      #ifdef RADIATION
                           BL_TO_FORTRAN_ANYD(Erborder[mfi]),
                           BL_TO_FORTRAN_ANYD(rad_flux[idir].hostFab()),
      #endif
                           BL_TO_FORTRAN_ANYD(area[idir][mfi]),
                           BL_TO_FORTRAN_ANYD(volume[mfi]),
                           BL_TO_FORTRAN_ANYD(div.hostFab()),
                           ZFILL(dx), dt);

      }



      AsyncFab pdivu(bx, 1);

      // conservative update

      ctu_consup(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                 BL_TO_FORTRAN_ANYD(Sborder[mfi]),
                 BL_TO_FORTRAN_ANYD(q[mfi]),
                 BL_TO_FORTRAN_ANYD(shk.hostFab()),
                 BL_TO_FORTRAN_ANYD(S_new[mfi]),
                 BL_TO_FORTRAN_ANYD(hydro_source[mfi]),
                 D_DECL(BL_TO_FORTRAN_ANYD(flux[0].hostFab()),
                        BL_TO_FORTRAN_ANYD(flux[1].hostFab()),
                        BL_TO_FORTRAN_ANYD(flux[2].hostFab())),
  #ifdef RADIATION
                 BL_TO_FORTRAN_ANYD(Erborder[mfi]),
                 BL_TO_FORTRAN_ANYD(Er_new[mfi]),
                 D_DECL(BL_TO_FORTRAN_ANYD(rad_flux[0].hostFab()),
                        BL_TO_FORTRAN_ANYD(rad_flux[1].hostFab()),
                        BL_TO_FORTRAN_ANYD(rad_flux[2].hostFab())),
                 &priv_nstep_fsp,
  #endif
                 D_DECL(BL_TO_FORTRAN_ANYD(qe[0].hostFab()),
                        BL_TO_FORTRAN_ANYD(qe[1].hostFab()),
                        BL_TO_FORTRAN_ANYD(qe[2].hostFab())),
                 D_DECL(BL_TO_FORTRAN_ANYD(area[0][mfi]),
                        BL_TO_FORTRAN_ANYD(area[1][mfi]),
                        BL_TO_FORTRAN_ANYD(area[2][mfi])),
                 BL_TO_FORTRAN_ANYD(volume[mfi]),
                 BL_TO_FORTRAN_ANYD(pdivu.hostFab()),
                 ZFILL(dx), dt);

      pdivu.clear();

  #ifdef RADIATION
      nstep_fsp = std::max(nstep_fsp, priv_nstep_fsp);
  #endif

      for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

        const Box& nbx = mfi.nodaltilebox(idir);

        scale_flux(ARLIM_3D(nbx.loVect()), ARLIM_3D(nbx.hiVect()),
  #if AMREX_SPACEDIM == 1
                   BL_TO_FORTRAN_ANYD(qe[idir].hostFab()),
  #endif
                   BL_TO_FORTRAN_ANYD(flux[idir].hostFab()),
                   BL_TO_FORTRAN_ANYD(area[idir][mfi]), dt);

  #ifdef RADIATION
        scale_rad_flux(ARLIM_3D(nbx.loVect()), ARLIM_3D(nbx.hiVect()),
                       BL_TO_FORTRAN_ANYD(rad_flux[idir].hostFab()),
                       BL_TO_FORTRAN_ANYD(area[idir][mfi]), dt);
  #endif

  #if AMREX_SPACEDIM <= 2
        if (idir == 0) {
          // get the scaled radial pressure -- we need to treat this specially
          // TODO: we should be able to do this entirely in C++, but we need to
          // know the value of mom_flux_has_p
          store_pradial(ARLIM_3D(nbx.loVect()), ARLIM_3D(nbx.hiVect()),
                        BL_TO_FORTRAN_ANYD(qe[idir].hostFab()),
                        BL_TO_FORTRAN_ANYD(pradial.hostFab()), dt);
        }
  #endif
        // Store the fluxes from this advance.

        // For normal integration we want to add the fluxes from this advance
        // since we may be subcycling the timestep. But for SDC integration
        // we want to copy the fluxes since we expect that there will not be
        // subcycling and we only want the last iteration's fluxes.
  #ifndef SDC
        (*fluxes[idir])[mfi].plus(flux[idir].hostFab(), nbx, 0, 0, NUM_STATE);
  #ifdef RADIATION
        (*rad_fluxes[idir])[mfi].plus(rad_flux[idir].hostFab(), nbx, 0, 0, Radiation::nGroups);
  #endif
  #else
        (*fluxes[idir])[mfi].copy(flux[idir].hostFab(), nbx, 0, nbx, 0, NUM_STATE);
  #ifdef RADIATION
        (*rad_fluxes[idir])[mfi].copy(rad_flux[idir].hostFab() nbx,0, nbx, 0, Radiation::nGroups);
  #endif
  #endif
        (*mass_fluxes[idir])[mfi].copy(flux[idir].hostFab(), nbx, Density, nbx, 0, 1);
      } // idir loop

  #if (AMREX_SPACEDIM <= 2)
      if (!Geometry::IsCartesian()) {
  #ifndef SDC
        P_radial[mfi].plus(pradial.hostFab(), mfi.nodaltilebox(0), 0, 0, 1);
  #else
        P_radial[mfi].copy(pradial.hostFab(), mfi.nodaltilebox(0), 0, mfi.nodaltilebox(0), 0, 1);
  #endif
      }
  #endif

      if (track_grid_losses == 1) {

          const Box& bx = mfi.tilebox();

          ca_track_grid_losses(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                               D_DECL(BL_TO_FORTRAN_ANYD(flux[0].hostFab()),
                                      BL_TO_FORTRAN_ANYD(flux[1].hostFab()),
                                      BL_TO_FORTRAN_ANYD(flux[2].hostFab())),
                               mass_lost, xmom_lost, ymom_lost, zmom_lost,
                               eden_lost, xang_lost, yang_lost, zang_lost);
      }

    } // MFIter loop

  } // OMP loop


  BL_PROFILE_VAR_STOP(CA_UMDRV);

#ifdef RADIATION
  if (radiation->verbose>=1) {
#ifdef BL_LAZY
    Lazy::QueueReduction( [=] () mutable {
#endif
       ParallelDescriptor::ReduceIntMax(nstep_fsp, ParallelDescriptor::IOProcessorNumber());
       if (ParallelDescriptor::IOProcessor() && nstep_fsp > 0) {
         std::cout << "Radiation f-space advection on level " << level
                   << " takes as many as " << nstep_fsp;
         if (nstep_fsp == 1) {
           std::cout<< " substep.\n";
         }
         else {
           std::cout<< " substeps.\n";
         }
       }
#ifdef BL_LAZY
     });
#endif
  }
#else
  // Flush Fortran output

  if (verbose)
    flush_output();

  if (track_grid_losses)
    {
      material_lost_through_boundary_temp[0] += mass_lost;
      material_lost_through_boundary_temp[1] += xmom_lost;
      material_lost_through_boundary_temp[2] += ymom_lost;
      material_lost_through_boundary_temp[3] += zmom_lost;
      material_lost_through_boundary_temp[4] += eden_lost;
      material_lost_through_boundary_temp[5] += xang_lost;
      material_lost_through_boundary_temp[6] += yang_lost;
      material_lost_through_boundary_temp[7] += zang_lost;
    }

  if (print_update_diagnostics)
    {

      bool local = true;
      Vector<Real> hydro_update = evaluate_source_change(hydro_source, dt, local);

#ifdef BL_LAZY
      Lazy::QueueReduction( [=] () mutable {
#endif
         ParallelDescriptor::ReduceRealSum(hydro_update.dataPtr(), hydro_update.size(), ParallelDescriptor::IOProcessorNumber());

         if (ParallelDescriptor::IOProcessor())
           std::cout << std::endl << "  Contributions to the state from the hydro source:" << std::endl;

         print_source_change(hydro_update);

#ifdef BL_LAZY
      });
#endif
    }
#endif

  if (verbose && ParallelDescriptor::IOProcessor())
    std::cout << "... Leaving hydro advance" << std::endl << std::endl;

  if (verbose > 0)
    {
      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
      Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
      Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

	if (ParallelDescriptor::IOProcessor())
	  std::cout << "Castro::construct_hydro_source() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
	});
#endif
    }

}
#endif



void
Castro::construct_mol_hydro_source(Real time, Real dt)
{

  BL_PROFILE("Castro::construct_mol_hydro_source()");

  // this constructs the hydrodynamic source (essentially the flux
  // divergence) using method of lines integration.  The output, as a
  // update to the state, is stored in the k_mol array of multifabs.

  const Real strt_time = ParallelDescriptor::second();

  if (verbose && ParallelDescriptor::IOProcessor())
    std::cout << "... hydro MOL stage " << mol_iteration << std::endl;


  // we'll add each stage's contribution to -div{F(U)} as we compute them
  if (mol_iteration == 0) {
    hydro_source.setVal(0.0);
  }


  const Real *dx = geom.CellSize();

  MultiFab& S_new = get_new_data(State_Type);

  MultiFab& k_stage = *k_mol[mol_iteration];

#ifdef RADIATION
  MultiFab& Er_new = get_new_data(Rad_Type);

  if (!Radiation::rad_hydro_combined) {
    amrex::Abort("Castro::construct_mol_hydro_source -- we don't implement a mode where we have radiation, but it is not coupled to hydro");
  }

  int nstep_fsp = -1;
#endif

  BL_PROFILE_VAR("Castro::advance_hydro_ca_umdrv()", CA_UMDRV);

  const int* domain_lo = geom.Domain().loVect();
  const int* domain_hi = geom.Domain().hiVect();

#ifndef AMREX_USE_CUDA

#ifdef _OPENMP
#ifdef RADIATION
#pragma omp parallel reduction(max:nstep_fsp)
#endif
#endif
  {

    FArrayBox flux[AMREX_SPACEDIM];
#if (AMREX_SPACEDIM <= 2)
    FArrayBox pradial(Box::TheUnitBox(),1);
#endif
#ifdef RADIATION
    FArrayBox rad_flux[AMREX_SPACEDIM];
#endif

#ifdef RADIATION
    int priv_nstep_fsp = -1;
#endif
    // The fourth order stuff cannot do tiling because of the Laplacian corrections
    for (MFIter mfi(S_new, (fourth_order) ? no_tile_size : hydro_tile_size); mfi.isValid(); ++mfi)
      {
	const Box& bx  = mfi.tilebox();

	const int* lo = bx.loVect();
	const int* hi = bx.hiVect();

	FArrayBox &statein  = Sborder[mfi];
	FArrayBox &stateout = S_new[mfi];

	FArrayBox &source_in  = sources_for_hydro[mfi];

	// the output of this will be stored in the correct stage MF
	FArrayBox &source_out = k_stage[mfi];
	FArrayBox &source_hydro_only = hydro_source[mfi];

#ifdef RADIATION
	FArrayBox &Er = Erborder[mfi];
	FArrayBox &lam = lamborder[mfi];
	FArrayBox &Erout = Er_new[mfi];
#endif

	// All cate fabs for fluxes
	for (int i = 0; i < AMREX_SPACEDIM ; i++)  {
	  const Box& bxtmp = amrex::surroundingNodes(bx,i);
	  flux[i].resize(bxtmp,NUM_STATE);
#ifdef RADIATION
	  rad_flux[i].resize(bxtmp,Radiation::nGroups);
#endif
	}

#if (AMREX_SPACEDIM <= 2)
	if (!Geometry::IsCartesian()) {
	  pradial.resize(amrex::surroundingNodes(bx,0),1);
	}
#endif
        if (fourth_order) {
          ca_fourth_single_stage
            (ARLIM_3D(lo), ARLIM_3D(hi), &time, ARLIM_3D(domain_lo), ARLIM_3D(domain_hi),
             &(b_mol[mol_iteration]),
             BL_TO_FORTRAN_ANYD(statein), 
             BL_TO_FORTRAN_ANYD(stateout),
             BL_TO_FORTRAN_ANYD(q[mfi]),
             BL_TO_FORTRAN_ANYD(q_bar[mfi]),
             BL_TO_FORTRAN_ANYD(qaux[mfi]),
             BL_TO_FORTRAN_ANYD(source_in),
             BL_TO_FORTRAN_ANYD(source_out),
             BL_TO_FORTRAN_ANYD(source_hydro_only),
             ZFILL(dx), &dt,
             D_DECL(BL_TO_FORTRAN_ANYD(flux[0]),
                    BL_TO_FORTRAN_ANYD(flux[1]),
                    BL_TO_FORTRAN_ANYD(flux[2])),
             D_DECL(BL_TO_FORTRAN_ANYD(area[0][mfi]),
                    BL_TO_FORTRAN_ANYD(area[1][mfi]),
                    BL_TO_FORTRAN_ANYD(area[2][mfi])),
#if (AMREX_SPACEDIM < 3)
             BL_TO_FORTRAN_ANYD(pradial),
             BL_TO_FORTRAN_ANYD(dLogArea[0][mfi]),
#endif
             BL_TO_FORTRAN_ANYD(volume[mfi]),
             verbose);

        } else {
          ca_mol_single_stage
            (ARLIM_3D(lo), ARLIM_3D(hi), &time, ARLIM_3D(domain_lo), ARLIM_3D(domain_hi),
             &(b_mol[mol_iteration]),
             BL_TO_FORTRAN_ANYD(statein), 
             BL_TO_FORTRAN_ANYD(stateout),
             BL_TO_FORTRAN_ANYD(q[mfi]),
             BL_TO_FORTRAN_ANYD(qaux[mfi]),
             BL_TO_FORTRAN_ANYD(source_in),
             BL_TO_FORTRAN_ANYD(source_out),
             BL_TO_FORTRAN_ANYD(source_hydro_only),
             ZFILL(dx), &dt,
             D_DECL(BL_TO_FORTRAN_ANYD(flux[0]),
                    BL_TO_FORTRAN_ANYD(flux[1]),
                    BL_TO_FORTRAN_ANYD(flux[2])),
             D_DECL(BL_TO_FORTRAN_ANYD(area[0][mfi]),
                    BL_TO_FORTRAN_ANYD(area[1][mfi]),
                    BL_TO_FORTRAN_ANYD(area[2][mfi])),
#if (AMREX_SPACEDIM < 3)
             BL_TO_FORTRAN_ANYD(pradial),
             BL_TO_FORTRAN_ANYD(dLogArea[0][mfi]),
#endif
             BL_TO_FORTRAN_ANYD(volume[mfi]),
             verbose);
        }

	// Store the fluxes from this advance -- we weight them by the
	// integrator weight for this stage
	for (int i = 0; i < AMREX_SPACEDIM ; i++) {
	  (*fluxes    [i])[mfi].saxpy(b_mol[mol_iteration], flux[i], 
				      mfi.nodaltilebox(i), mfi.nodaltilebox(i), 0, 0, NUM_STATE);
#ifdef RADIATION
	  (*rad_fluxes[i])[mfi].saxpy(b_mol[mol_iteration], rad_flux[i], 
				      mfi.nodaltilebox(i), mfi.nodaltilebox(i), 0, 0, Radiation::nGroups);
#endif
	}

#if (AMREX_SPACEDIM <= 2)
	if (!Geometry::IsCartesian()) {
	  P_radial[mfi].saxpy(b_mol[mol_iteration], pradial,
                              mfi.nodaltilebox(0), mfi.nodaltilebox(0), 0, 0, 1);
	}
#endif
      } // MFIter loop

#ifdef RADIATION
    nstep_fsp = std::max(nstep_fsp, priv_nstep_fsp);
#endif
  }  // end of omp parallel region

#else
  // CUDA version

#ifndef RADIATION

  MultiFab flatn;
  flatn.define(grids, dmap, 1, 1);

  MultiFab div;
  div.define(grids, dmap, 1, 1);

  MultiFab qm;
  qm.define(grids, dmap, AMREX_SPACEDIM*NQ, 2);

  MultiFab qp;
  qp.define(grids, dmap, AMREX_SPACEDIM*NQ, 2);

  MultiFab shk;
  shk.define(grids, dmap, 1, 1);


  MultiFab flux[AMREX_SPACEDIM];
  MultiFab qe[AMREX_SPACEDIM];
  MultiFab qi[AMREX_SPACEDIM];

  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      flux[i].define(getEdgeBoxArray(i), dmap, NUM_STATE, 0);
      qe[i].define(getEdgeBoxArray(i), dmap, NGDNV, 0);
      qi[i].define(getEdgeBoxArray(i), dmap, NQ, 0);
  }


#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

      const Box& obx = mfi.growntilebox(1);

      // Compute divergence of velocity field.
#pragma gpu
      divu(AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
           BL_TO_FORTRAN_ANYD(q[mfi]),
           AMREX_REAL_ANYD(dx),
           BL_TO_FORTRAN_ANYD(div[mfi]));

      // Compute flattening coefficient for slope calculations.
#pragma gpu
      ca_uflatten
          (AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
           BL_TO_FORTRAN_ANYD(q[mfi]),
           BL_TO_FORTRAN_ANYD(flatn[mfi]), QPRES+1);

      // Do PPM reconstruction to the zone edges.
      int put_on_edges = 1;

#pragma gpu
      ca_ppm_reconstruct
          (AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()), put_on_edges,
           BL_TO_FORTRAN_ANYD(q[mfi]), NQ, 1, NQ,
           BL_TO_FORTRAN_ANYD(flatn[mfi]),
           BL_TO_FORTRAN_ANYD(qm[mfi]),
           BL_TO_FORTRAN_ANYD(qp[mfi]), NQ, 1, NQ);

      // Compute the shk variable
#pragma gpu
      ca_shock
        (AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
         BL_TO_FORTRAN_ANYD(q[mfi]),
         BL_TO_FORTRAN_ANYD(shk[mfi]),
         AMREX_REAL_ANYD(dx));

  } // MFIter loop



#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

      for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

          const Box& ebx = mfi.nodaltilebox(idir);

          int idir_f = idir + 1;

#pragma gpu
          ca_construct_flux_cuda
              (AMREX_INT_ANYD(ebx.loVect()), AMREX_INT_ANYD(ebx.hiVect()),
               AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
               AMREX_REAL_ANYD(dx), dt,
               idir_f,
               BL_TO_FORTRAN_ANYD(Sborder[mfi]),
               BL_TO_FORTRAN_ANYD(div[mfi]),
               BL_TO_FORTRAN_ANYD(qaux[mfi]),
               BL_TO_FORTRAN_ANYD(shk[mfi]),
               BL_TO_FORTRAN_ANYD(qm[mfi]),
               BL_TO_FORTRAN_ANYD(qp[mfi]),
               BL_TO_FORTRAN_ANYD(qi[idir][mfi]),
               BL_TO_FORTRAN_ANYD(flux[idir][mfi]),
               BL_TO_FORTRAN_ANYD(area[idir][mfi]));

#pragma gpu
          ca_store_godunov_state
            (AMREX_INT_ANYD(ebx.loVect()), AMREX_INT_ANYD(ebx.hiVect()),
             BL_TO_FORTRAN_ANYD(qi[idir][mfi]),
             BL_TO_FORTRAN_ANYD(qe[idir][mfi]));

          // Store the fluxes from this advance -- we weight them by the
          // integrator weight for this stage
#ifdef AMREX_USE_CUDA
          Cuda::Device::synchronize();  // because saxpy below is run on cpu
#endif
          (*fluxes[idir])[mfi].saxpy(b_mol[mol_iteration], flux[idir][mfi], ebx, ebx, 0, 0, NUM_STATE);

      }

  } // MFIter loop


#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

      const Box& bx = mfi.tilebox();

#pragma gpu
      ca_construct_hydro_update_cuda
          (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
           AMREX_REAL_ANYD(dx), dt,
           BL_TO_FORTRAN_ANYD(qe[0][mfi]),
           BL_TO_FORTRAN_ANYD(qe[1][mfi]),
           BL_TO_FORTRAN_ANYD(qe[2][mfi]),
           BL_TO_FORTRAN_ANYD(flux[0][mfi]),
           BL_TO_FORTRAN_ANYD(flux[1][mfi]),
           BL_TO_FORTRAN_ANYD(flux[2][mfi]),
           BL_TO_FORTRAN_ANYD(area[0][mfi]),
           BL_TO_FORTRAN_ANYD(area[1][mfi]),
           BL_TO_FORTRAN_ANYD(area[2][mfi]),
           BL_TO_FORTRAN_ANYD(volume[mfi]),
           BL_TO_FORTRAN_ANYD(sources_for_hydro[mfi]),
           BL_TO_FORTRAN_ANYD(k_stage[mfi]));

  } // MFIter loop

#endif // RADIATION

#endif // CUDA check

  BL_PROFILE_VAR_STOP(CA_UMDRV);

  // Flush Fortran output

  if (verbose)
    flush_output();


  if (print_update_diagnostics)
    {

      bool local = true;
      Vector<Real> hydro_update = evaluate_source_change(k_stage, dt, local);

#ifdef BL_LAZY
      Lazy::QueueReduction( [=] () mutable {
#endif
	  ParallelDescriptor::ReduceRealSum(hydro_update.dataPtr(), hydro_update.size(), ParallelDescriptor::IOProcessorNumber());

	  if (ParallelDescriptor::IOProcessor())
	    std::cout << std::endl << "  Contributions to the state from the hydro source:" << std::endl;

	  print_source_change(hydro_update);

#ifdef BL_LAZY
	});
#endif
    }

    if (verbose > 0)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

	if (ParallelDescriptor::IOProcessor())
	  std::cout << "Castro::construct_mol_hydro_source() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
	});
#endif
    }


}



void
Castro::cons_to_prim(const Real time)
{

#ifdef RADIATION
    AmrLevel::FillPatch(*this, Erborder, NUM_GROW, time, Rad_Type, 0, Radiation::nGroups);

    MultiFab lamborder(grids, dmap, Radiation::nGroups, NUM_GROW);
    if (radiation->pure_hydro) {
      lamborder.setVal(0.0, NUM_GROW);
    }
    else {
      radiation->compute_limiter(level, grids, Sborder, Erborder, lamborder);
    }
#endif

    MultiFab& S_new = get_new_data(State_Type);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

        const Box& qbx = mfi.growntilebox(NUM_GROW);

        // Convert the conservative state to the primitive variable state.
        // This fills both q and qaux.

#pragma gpu
        ca_ctoprim(AMREX_INT_ANYD(qbx.loVect()), AMREX_INT_ANYD(qbx.hiVect()),
                   BL_TO_FORTRAN_ANYD(Sborder[mfi]),
#ifdef RADIATION
                   BL_TO_FORTRAN_ANYD(Erborder[mfi]),
                   BL_TO_FORTRAN_ANYD(lamborder[mfi]),
#endif
                   BL_TO_FORTRAN_ANYD(q[mfi]),
                   BL_TO_FORTRAN_ANYD(qaux[mfi]));

        // Convert the source terms expressed as sources to the conserved state to those
        // expressed as sources for the primitive state.
#ifndef AMREX_USE_CUDA
        if (do_ctu) {
          ca_srctoprim(BL_TO_FORTRAN_BOX(qbx),
                       BL_TO_FORTRAN_ANYD(q[mfi]),
                       BL_TO_FORTRAN_ANYD(qaux[mfi]),
                       BL_TO_FORTRAN_ANYD(sources_for_hydro[mfi]),
                       BL_TO_FORTRAN_ANYD(src_q[mfi]));
        }
#endif

#ifndef RADIATION

        // Add in the reactions source term; only done in SDC.

#ifdef SDC
#ifdef REACTIONS
        MultiFab& SDC_react_source = get_new_data(SDC_React_Type);

        if (do_react)
	    src_q[mfi].plus(SDC_react_source[mfi],qbx,qbx,0,0,QVAR);
#endif
#endif
#endif
      
    }

}


#ifndef AMREX_USE_CUDA
void
Castro::cons_to_prim_fourth(const Real time)
{
  // convert the conservative state cell averages to primitive cell
  // averages with 4th order accuracy

    MultiFab& S_new = get_new_data(State_Type);

    // we don't support radiation here
#ifdef RADIATION
    amrex::Abort("radiation not supported to fourth order");
#else
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

      const Box& qbx = mfi.growntilebox(NUM_GROW);
      const Box& qbxm1 = mfi.growntilebox(NUM_GROW-1);

      // note: these conversions are using a growntilebox, so it
      // will include ghost cells

      // convert U_avg to U_cc -- this will use a Laplacian
      // operation and will result in U_cc defined only on
      // NUM_GROW-1 ghost cells at the end.
      FArrayBox U_cc;
      U_cc.resize(qbx, NUM_STATE);

      ca_make_cell_center(BL_TO_FORTRAN_BOX(qbxm1),
                          BL_TO_FORTRAN_FAB(Sborder[mfi]),
                          BL_TO_FORTRAN_FAB(U_cc));

      // convert U_avg to q_bar -- this will be done on all NUM_GROW
      // ghost cells.
      FArrayBox qaux_bar;
      qaux_bar.resize(qbx, NQAUX);

      ca_ctoprim(BL_TO_FORTRAN_BOX(qbx),
                 BL_TO_FORTRAN_ANYD(Sborder[mfi]),
                 BL_TO_FORTRAN_ANYD(q_bar[mfi]),
                 BL_TO_FORTRAN_ANYD(qaux_bar));

      // this is what we should construct the flattening coefficient
      // from

      // convert U_cc to q_cc (we'll store this temporarily in q,
      // qaux).  This will remain valid only on the NUM_GROW-1 ghost
      // cells.
      ca_ctoprim(BL_TO_FORTRAN_BOX(qbxm1),
                 BL_TO_FORTRAN_ANYD(U_cc),
                 BL_TO_FORTRAN_ANYD(q[mfi]),
                 BL_TO_FORTRAN_ANYD(qaux[mfi]));
    }


    // check for NaNs
    check_for_nan(q);
    check_for_nan(q_bar);


#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

      const Box& qbxm1 = mfi.growntilebox(NUM_GROW-1);

      // now convert q, qaux into 4th order accurate averages
      // this will create q, qaux in NUM_GROW-1 ghost cells, but that's
      // we need here

      ca_make_fourth_average(BL_TO_FORTRAN_BOX(qbxm1),
                             BL_TO_FORTRAN_FAB(q[mfi]),
                             BL_TO_FORTRAN_FAB(q_bar[mfi]));

      // not sure if we need to convert qaux this way, or if we can
      // just evaluate it (we may not need qaux at all actually)

    }

    check_for_nan(q_bar);
#endif // RADIATION
}
#endif



void
Castro::check_for_cfl_violation(const Real dt)
{

    Real courno = -1.0e+200;

    const Real *dx = geom.CellSize();

    MultiFab& S_new = get_new_data(State_Type);

#ifdef _OPENMP
#pragma omp parallel reduction(max:courno)
#endif
    for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

#pragma gpu
        ca_compute_cfl(BL_TO_FORTRAN_BOX(bx),
                       BL_TO_FORTRAN_ANYD(q[mfi]),
                       BL_TO_FORTRAN_ANYD(qaux[mfi]),
                       dt, AMREX_REAL_ANYD(dx), AMREX_MFITER_REDUCE_MAX(&courno), print_fortran_warnings);

    }

    ParallelDescriptor::ReduceRealMax(courno);

    if (courno > 1.0) {
        amrex::Print() << "WARNING -- EFFECTIVE CFL AT LEVEL " << level << " IS " << courno << std::endl << std::endl;

        cfl_violation = 1;
    }

}
