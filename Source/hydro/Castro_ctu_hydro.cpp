#include <Castro.H>
#include <Castro_util.H>
#include <Castro_F.H>
#include <Castro_hydro.H>

#ifdef RADIATION
#include <Radiation.H>
#endif

#ifdef HYBRID_MOMENTUM
#include <hybrid.H>
#endif

using namespace amrex;

void
Castro::construct_ctu_hydro_source(Real time, Real dt)
{

#ifndef TRUE_SDC

  BL_PROFILE("Castro::construct_ctu_hydro_source()");

  const Real strt_time = ParallelDescriptor::second();

  // this constructs the hydrodynamic source (essentially the flux
  // divergence) using the CTU framework for unsplit hydrodynamics

  if (verbose && ParallelDescriptor::IOProcessor())
    std::cout << "... Entering construct_ctu_hydro_source()" << std::endl << std::endl;

#ifdef HYBRID_MOMENTUM
  GeometryData geomdata = geom.data();
#endif

#if AMREX_SPACEDIM == 2
  int coord = geom.Coord();
#endif

  const Real *dx = geom.CellSize();

  MultiFab& S_new = get_new_data(State_Type);

  // we will treat the hydro source as any other source term

#ifdef RADIATION
  MultiFab& Er_new = get_new_data(Rad_Type);

  if (!Radiation::rad_hydro_combined) {
    amrex::Abort("Castro::construct_ctu_hydro_source -- we don't implement a mode where we have radiation, but it is not coupled to hydro");
  }

  int nstep_fsp = -1;

  AmrLevel::FillPatch(*this, Erborder, NUM_GROW, time, Rad_Type, 0, Radiation::nGroups);

  MultiFab lamborder(grids, dmap, Radiation::nGroups, NUM_GROW);
  if (radiation->pure_hydro) {
      lamborder.setVal(0.0, NUM_GROW);
  }
  else {
      radiation->compute_limiter(level, grids, Sborder, Erborder, lamborder);
  }
#endif

#ifdef _OPENMP
#ifdef RADIATION
#pragma omp parallel reduction(max:nstep_fsp)
#else
#pragma omp parallel
#endif
#endif
  {

#ifdef RADIATION
    int priv_nstep_fsp = -1;
#endif

    // Declare local storage now. This should be done outside the MFIter loop,
    // and then we will resize the Fabs in each MFIter loop iteration. Then,
    // we apply an Elixir to ensure that their memory is saved until it is no
    // longer needed (only relevant for the asynchronous case, usually on GPUs).

    FArrayBox flatn;
#ifdef RADIATION
    FArrayBox flatg;
#endif
    FArrayBox shk;
    FArrayBox q, qaux;
    FArrayBox src_q;
    FArrayBox qxm, qxp;
#if AMREX_SPACEDIM >= 2
    FArrayBox qym, qyp;
#endif
#if AMREX_SPACEDIM == 3
    FArrayBox qzm, qzp;
#endif
    FArrayBox div;
#if AMREX_SPACEDIM >= 2
    FArrayBox ftmp1, ftmp2;
#ifdef RADIATION
    FArrayBox rftmp1, rftmp2;
#endif
    FArrayBox qgdnvtmp1, qgdnvtmp2;
    FArrayBox ql, qr;
#endif
    FArrayBox flux[AMREX_SPACEDIM], qe[AMREX_SPACEDIM];
#ifdef RADIATION
    FArrayBox rad_flux[AMREX_SPACEDIM];
#endif
#if AMREX_SPACEDIM <= 2
    FArrayBox pradial;
#endif
#if AMREX_SPACEDIM == 3
    FArrayBox qmyx, qpyx;
    FArrayBox qmzx, qpzx;
    FArrayBox qmxy, qpxy;
    FArrayBox qmzy, qpzy;
    FArrayBox qmxz, qpxz;
    FArrayBox qmyz, qpyz;
#endif

#ifdef AMREX_USE_GPU
    size_t starting_size = MultiFab::queryMemUsage("AmrLevel_Level_" + std::to_string(level));
    size_t current_size = starting_size;
#endif

    MultiFab& old_source = get_old_data(Source_Type);

    for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

      size_t fab_size = 0;

      // the valid region box
      const Box& bx = mfi.tilebox();

      const Box& obx = amrex::grow(bx, 1);

      flatn.resize(obx, 1);
      Elixir elix_flatn = flatn.elixir();
      fab_size += flatn.nBytes();

#ifdef RADIATION
      flatg.resize(obx, 1);
      Elixir elix_flatg = flatg.elixir();
      fab_size += flatg.nBytes();
#endif

      // If we are oversubscribing the GPU, performance of the hydro will be constrained
      // due to its heavy memory requirements. We can help the situation by prefetching in
      // all the data we will need, and then prefetching it out at the end. This at least
      // improves performance by mitigating the number of unified memory page faults.

      // An empirical threshold on NVIDIA GPUs is that we're probably oversubscribing if
      // there are less than 10 MB left.

      bool oversubscribed = false;

#ifdef AMREX_USE_GPU
      if (Gpu::Device::freeMemAvailable() < 10000000) {
          oversubscribed = true;
      }
#endif

      if (oversubscribed) {
          volume[mfi].prefetchToDevice();
          Sborder[mfi].prefetchToDevice();
          S_new[mfi].prefetchToDevice();
          for (int i = 0; i < AMREX_SPACEDIM; ++i) {
              area[i][mfi].prefetchToDevice();
              (*fluxes[i])[mfi].prefetchToDevice();
          }
#if AMREX_SPACEDIM < 3
          dLogArea[0][mfi].prefetchToDevice();
          P_radial[mfi].prefetchToDevice();
#endif
#ifdef RADIATION
          Erborder[mfi].prefetchToDevice();
          Er_new[mfi].prefetchToDevice();
#endif
      }

      // Compute the primitive variables (both q and qaux) from
      // the conserved variables.

      const Box& qbx = amrex::grow(bx, NUM_GROW);
      const Box& qbx3 = amrex::grow(bx, 3);

      q.resize(qbx, NQ);
      Elixir elix_q = q.elixir();
      fab_size += q.nBytes();
      Array4<Real> const q_arr = q.array();

      qaux.resize(qbx, NQAUX);
      Elixir elix_qaux = qaux.elixir();
      fab_size += qaux.nBytes();
      Array4<Real> const qaux_arr = qaux.array();

      ctoprim(qbx, time, Sborder.array(mfi),
#ifdef RADIATION
              Erborder.array(mfi), lamborder.array(mfi),
#endif
              q_arr, qaux_arr);



      Array4<Real const> const areax_arr = area[0].array(mfi);
#if AMREX_SPACEDIM >= 2
      Array4<Real const> const areay_arr = area[1].array(mfi);
#endif
#if AMREX_SPACEDIM == 3
      Array4<Real const> const areaz_arr = area[2].array(mfi);
#endif

      Array4<Real> const vol_arr = volume.array(mfi);

#if AMREX_SPACEDIM < 3
      Array4<Real const> const dLogArea_arr = (dLogArea[0]).array(mfi);
#endif

      // compute the flattening coefficient

      Array4<Real> const flatn_arr = flatn.array();
#ifdef RADIATION
      Array4<Real> const flatg_arr = flatg.array();
#endif

      if (first_order_hydro == 1) {
        amrex::ParallelFor(obx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {
          flatn_arr(i,j,k) = 0.0;
        });
      } else if (use_flattening == 1) {

        uflatten(obx, q_arr, flatn_arr, QPRES);

#ifdef RADIATION
        uflatten(obx, q_arr, flatg_arr, QPTOT);

        Real flatten_pp_thresh = radiation::flatten_pp_threshold;

        amrex::ParallelFor(obx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {
          flatn_arr(i,j,k) = flatn_arr(i,j,k) * flatg_arr(i,j,k);

          if (flatten_pp_thresh > 0.0) {
            if ( q_arr(i-1,j,k,QU) + q_arr(i,j-1,k,QV) + q_arr(i,j,k-1,QW) >
                 q_arr(i+1,j,k,QU) + q_arr(i,j+1,k,QV) + q_arr(i,j,k+1,QW) ) {

              if (q_arr(i,j,k,QPRES) < flatten_pp_thresh * q_arr(i,j,k,QPTOT)) {
                flatn_arr(i,j,k) = 0.0;
              }
            }
          }
        });
#endif

      } else {
        amrex::ParallelFor(obx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {
          flatn_arr(i,j,k) = 1.0;
        });
      }

      const Box& xbx = amrex::surroundingNodes(bx, 0);
      const Box& gxbx = amrex::grow(xbx, 1);
#if AMREX_SPACEDIM >= 2
      const Box& ybx = amrex::surroundingNodes(bx, 1);
      const Box& gybx = amrex::grow(ybx, 1);
#endif
#if AMREX_SPACEDIM == 3
      const Box& zbx = amrex::surroundingNodes(bx, 2);
      const Box& gzbx = amrex::grow(zbx, 1);
#endif

      shk.resize(obx, 1);
      Elixir elix_shk = shk.elixir();
      fab_size += shk.nBytes();

      Array4<Real> const shk_arr = shk.array();

      // Multidimensional shock detection
      // Used for the hybrid Riemann solver

#ifdef SHOCK_VAR
      bool compute_shock = true;
#else
      bool compute_shock = false;
#endif

      if (hybrid_riemann == 1 || compute_shock) {
        shock(obx, q_arr, shk_arr);
      }
      else {
        amrex::ParallelFor(obx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {
          shk_arr(i,j,k) = 0.0;
        });
      }

      // get the primitive variable hydro sources

      src_q.resize(qbx3, NQSRC);
      Elixir elix_src_q = src_q.elixir();
      fab_size += src_q.nBytes();
      Array4<Real> const src_q_arr = src_q.array();

      Array4<Real> const old_src_arr = old_source.array(mfi);
      Array4<Real> const src_corr_arr = source_corrector.array(mfi);

      src_to_prim(qbx3, dt, q_arr, old_src_arr, src_corr_arr, src_q_arr);

#ifndef RADIATION
#ifdef SIMPLIFIED_SDC
#ifdef REACTIONS
        // Add in the reactions source term; only done in simplified SDC.

        if (time_integration_method == SimplifiedSpectralDeferredCorrections) {

            MultiFab& SDC_react_source = get_new_data(Simplified_SDC_React_Type);

            if (do_react)
              src_q.plus<RunOn::Device>(SDC_react_source[mfi], qbx3, qbx3, 0, 0, NQSRC);

        }
#endif
#endif
#endif


      // work on the interface states

      qxm.resize(obx, NQ);
      Elixir elix_qxm = qxm.elixir();
      fab_size += shk.nBytes();

      qxp.resize(obx, NQ);
      Elixir elix_qxp = qxp.elixir();
      fab_size += qxp.nBytes();

      Array4<Real> const qxm_arr = qxm.array();
      Array4<Real> const qxp_arr = qxp.array();

#if AMREX_SPACEDIM >= 2
      qym.resize(obx, NQ);
      Elixir elix_qym = qym.elixir();
      fab_size += qym.nBytes();

      qyp.resize(obx, NQ);
      Elixir elix_qyp = qyp.elixir();
      fab_size += qyp.nBytes();

      Array4<Real> const qym_arr = qym.array();
      Array4<Real> const qyp_arr = qyp.array();

#endif

#if AMREX_SPACEDIM == 3
      qzm.resize(obx, NQ);
      Elixir elix_qzm = qzm.elixir();
      fab_size += qzm.nBytes();

      qzp.resize(obx, NQ);
      Elixir elix_qzp = qzp.elixir();
      fab_size += qzp.nBytes();

      Array4<Real> const qzm_arr = qzm.array();
      Array4<Real> const qzp_arr = qzp.array();

#endif

      if (ppm_type == 0) {

        ctu_plm_states(obx, bx,
                       q_arr,
                       flatn_arr,
                       qaux_arr,
                       src_q_arr,
                       qxm_arr, qxp_arr,
#if AMREX_SPACEDIM >= 2
                       qym_arr, qyp_arr,
#endif
#if AMREX_SPACEDIM == 3
                       qzm_arr, qzp_arr,
#endif
#if (AMREX_SPACEDIM < 3)
                       dLogArea_arr,
#endif
                       dt);

      } else {

#ifdef RADIATION
        ctu_ppm_rad_states(obx, bx,
                           q_arr, flatn_arr, qaux_arr, src_q_arr,
                           qxm_arr, qxp_arr,
#if AMREX_SPACEDIM >= 2
                           qym_arr, qyp_arr,
#endif
#if AMREX_SPACEDIM == 3
                           qzm_arr, qzp_arr,
#endif
#if AMREX_SPACEDIM < 3
                           dLogArea_arr,
#endif
                           dt);
#else

        ctu_ppm_states(obx, bx,
                       q_arr, flatn_arr, qaux_arr, src_q_arr,
                       qxm_arr, qxp_arr,
#if AMREX_SPACEDIM >= 2
                       qym_arr, qyp_arr,
#endif
#if AMREX_SPACEDIM == 3
                       qzm_arr, qzp_arr,
#endif
#if AMREX_SPACEDIM < 3
                       dLogArea_arr,
#endif
                       dt);
#endif

      }

      div.resize(obx, 1);
      Elixir elix_div = div.elixir();
      fab_size += div.nBytes();
      auto div_arr = div.array();

      // compute divu -- we'll use this later when doing the artifical viscosity
      divu(obx, q_arr, div_arr);

      flux[0].resize(gxbx, NUM_STATE);
      Elixir elix_flux_x = flux[0].elixir();
      fab_size += flux[0].nBytes();
      Array4<Real> const flux0_arr = (flux[0]).array();

      qe[0].resize(gxbx, NGDNV);
      Elixir elix_qe_x = qe[0].elixir();
      auto qex_arr = qe[0].array();
      fab_size += qe[0].nBytes();

#ifdef RADIATION
      rad_flux[0].resize(gxbx, Radiation::nGroups);
      Elixir elix_rad_flux_x = rad_flux[0].elixir();
      fab_size += rad_flux[0].nBytes();
      auto rad_flux0_arr = (rad_flux[0]).array();
#endif

#if AMREX_SPACEDIM >= 2
      flux[1].resize(gybx, NUM_STATE);
      Elixir elix_flux_y = flux[1].elixir();
      fab_size += flux[1].nBytes();
      Array4<Real> const flux1_arr = (flux[1]).array();

      qe[1].resize(gybx, NGDNV);
      Elixir elix_qe_y = qe[1].elixir();
      auto qey_arr = qe[1].array();
      fab_size += qe[1].nBytes();

#ifdef RADIATION
      rad_flux[1].resize(gybx, Radiation::nGroups);
      Elixir elix_rad_flux_y = rad_flux[1].elixir();
      fab_size += rad_flux[1].nBytes();
      auto const rad_flux1_arr = (rad_flux[1]).array();
#endif
#endif

#if AMREX_SPACEDIM == 3
      flux[2].resize(gzbx, NUM_STATE);
      Elixir elix_flux_z = flux[2].elixir();
      fab_size += flux[2].nBytes();
      Array4<Real> const flux2_arr = (flux[2]).array();

      qe[2].resize(gzbx, NGDNV);
      Elixir elix_qe_z = qe[2].elixir();
      auto qez_arr = qe[2].array();
      fab_size += qe[2].nBytes();

#ifdef RADIATION
      rad_flux[2].resize(gzbx, Radiation::nGroups);
      Elixir elix_rad_flux_z = rad_flux[2].elixir();
      fab_size += rad_flux[2].nBytes();
      auto const rad_flux2_arr = (rad_flux[2]).array();
#endif
#endif

#if AMREX_SPACEDIM <= 2
      if (!Geom().IsCartesian()) {
          pradial.resize(xbx, 1);
      }
      Elixir elix_pradial = pradial.elixir();
      fab_size += pradial.nBytes();
#endif

#if AMREX_SPACEDIM == 1

#ifdef PRIM_SPECIES_HAVE_SOURCES
      // if we are doing species sources, add them here

      add_species_source_to_states(xbx, 0, dt,
                                   qxm_arr, qxp_arr, src_q_arr);
#endif

      // compute the fluxes through the x-interface

      cmpflx_plus_godunov(xbx,
                          qxm_arr, qxp_arr,
                          flux0_arr,
#ifdef RADIATION
                          rad_flux0_arr,
#endif
                          qex_arr,
                          qaux_arr,
                          shk_arr,
                          0, false);

#endif // 1-d



#if AMREX_SPACEDIM >= 2
      ftmp1.resize(obx, NUM_STATE);
      Elixir elix_ftmp1 = ftmp1.elixir();
      auto ftmp1_arr = ftmp1.array();
      fab_size += ftmp1.nBytes();

      ftmp2.resize(obx, NUM_STATE);
      Elixir elix_ftmp2 = ftmp2.elixir();
      auto ftmp2_arr = ftmp2.array();
      fab_size += ftmp2.nBytes();

#ifdef RADIATION
      rftmp1.resize(obx, Radiation::nGroups);
      Elixir elix_rftmp1 = rftmp1.elixir();
      auto rftmp1_arr = rftmp1.array();
      fab_size += rftmp1.nBytes();

      rftmp2.resize(obx, Radiation::nGroups);
      Elixir elix_rftmp2 = rftmp2.elixir();
      auto rftmp2_arr = rftmp2.array();
      fab_size += rftmp2.nBytes();
#endif

      qgdnvtmp1.resize(obx, NGDNV);
      Elixir elix_qgdnvtmp1 = qgdnvtmp1.elixir();
      auto qgdnvtmp1_arr = qgdnvtmp1.array();
      fab_size += qgdnvtmp1.nBytes();

#if AMREX_SPACEDIM == 3
      qgdnvtmp2.resize(obx, NGDNV);
      Elixir elix_qgdnvtmp2 = qgdnvtmp2.elixir();
      auto qgdnvtmp2_arr = qgdnvtmp2.array();
      fab_size += qgdnvtmp2.nBytes();
#endif

      ql.resize(obx, NQ);
      Elixir elix_ql = ql.elixir();
      auto ql_arr = ql.array();
      fab_size += ql.nBytes();

      qr.resize(obx, NQ);
      Elixir elix_qr = qr.elixir();
      auto qr_arr = qr.array();
      fab_size += qr.nBytes();
#endif



#if AMREX_SPACEDIM == 2

      const amrex::Real hdt = 0.5*dt;
      const amrex::Real hdtdx = 0.5*dt/dx[0];
      const amrex::Real hdtdy = 0.5*dt/dx[1];

      // compute F^x
      // [lo(1), lo(2)-1, 0], [hi(1)+1, hi(2)+1, 0]
      const Box& cxbx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0,1,0)));

      // ftmp1 = fx
      // rftmp1 = rfx
      // qgdnvtmp1 = qgdnxv
      cmpflx_plus_godunov(cxbx,
                          qxm_arr, qxp_arr,
                          ftmp1_arr,
#ifdef RADIATION
                          rftmp1_arr,
#endif
                          qgdnvtmp1_arr,
                          qaux_arr, shk_arr,
                          0, false);

      // compute F^y
      // [lo(1)-1, lo(2), 0], [hi(1)+1, hi(2)+1, 0]
      const Box& cybx = amrex::grow(ybx, IntVect(AMREX_D_DECL(1,0,0)));

      // ftmp2 = fy
      // rftmp2 = rfy
      cmpflx_plus_godunov(cybx,
                          qym_arr, qyp_arr,
                          ftmp2_arr,
#ifdef RADIATION
                          rftmp2_arr,
#endif
                          qey_arr,
                          qaux_arr, shk_arr,
                          1, false);

      // add the transverse flux difference in y to the x states
      // [lo(1), lo(2), 0], [hi(1)+1, hi(2), 0]

      // ftmp2 = fy
      // rftmp2 = rfy
      trans_single(xbx, 1, 0,
                   qxm_arr, ql_arr,
                   qxp_arr, qr_arr,
                   qaux_arr,
                   ftmp2_arr,
#ifdef RADIATION
                   rftmp2_arr,
#endif
                   qey_arr,
                   areay_arr,
                   vol_arr,
                   hdt, hdtdy);

      //reset_edge_state_thermo(xbx, ql.array());

      //reset_edge_state_thermo(xbx, qr.array());

      // solve the final Riemann problem axross the x-interfaces

#ifdef PRIM_SPECIES_HAVE_SOURCES
      add_species_source_to_states(xbx, 0, dt,
                                   ql_arr, qr_arr, src_q_arr);

#endif

      cmpflx_plus_godunov(xbx,
                          ql_arr, qr_arr,
                          flux0_arr,
#ifdef RADIATION
                          rad_flux0_arr,
#endif
                          qex_arr,
                          qaux_arr, shk_arr,
                          0, false);

      // add the transverse flux difference in x to the y states
      // [lo(1), lo(2), 0], [hi(1), hi(2)+1, 0]

      // ftmp1 = fx
      // rftmp1 = rfx
      // qgdnvtmp1 = qgdnvx

      trans_single(ybx, 0, 1,
                   qym_arr, ql_arr,
                   qyp_arr, qr_arr,
                   qaux_arr,
                   ftmp1_arr,
#ifdef RADIATION
                   rftmp1_arr,
#endif
                   qgdnvtmp1_arr,
                   areax_arr,
                   vol_arr,
                   hdt, hdtdx);

      //reset_edge_state_thermo(ybx, ql.array());

      //reset_edge_state_thermo(ybx, qr.array());


      // solve the final Riemann problem axross the y-interfaces

#ifdef PRIM_SPECIES_HAVE_SOURCES
      add_species_source_to_states(ybx, 1, dt,
                                   ql_arr, qr_arr, src_q_arr);

#endif

      cmpflx_plus_godunov(ybx,
                          ql_arr, qr_arr,
                          flux1_arr,
#ifdef RADIATION
                          rad_flux1_arr,
#endif
                          qey_arr,
                          qaux_arr, shk_arr,
                          1, false);
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
      const Box& cxbx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0,1,1)));

      // ftmp1 = fx
      // rftmp1 = rfx
      // qgdnvtmp1 = qgdnxv
      cmpflx_plus_godunov(cxbx,
                          qxm_arr, qxp_arr,
                          ftmp1_arr,
#ifdef RADIATION
                          rftmp1_arr,
#endif
                          qgdnvtmp1_arr,
                          qaux_arr, shk_arr,
                          0, false);


      // [lo(1), lo(2), lo(3)-1], [hi(1), hi(2)+1, hi(3)+1]
      const Box& tyxbx = amrex::grow(ybx, IntVect(AMREX_D_DECL(0,0,1)));

      qmyx.resize(tyxbx, NQ);
      Elixir elix_qmyx = qmyx.elixir();
      auto qmyx_arr = qmyx.array();
      fab_size += qmyx.nBytes();

      qpyx.resize(tyxbx, NQ);
      Elixir elix_qpyx = qpyx.elixir();
      auto qpyx_arr = qpyx.array();
      fab_size += qpyx.nBytes();

      // ftmp1 = fx
      // rftmp1 = rfx
      // qgdnvtmp1 = qgdnvx
      trans_single(tyxbx, 0, 1,
                   qym_arr, qmyx_arr,
                   qyp_arr, qpyx_arr,
                   qaux_arr,
                   ftmp1_arr,
#ifdef RADIATION
                   rftmp1_arr,
#endif
                   qgdnvtmp1_arr,
                   hdt, cdtdx);

      //reset_edge_state_thermo(tyxbx, qmyx.array());

      //reset_edge_state_thermo(tyxbx, qpyx.array());

      // [lo(1), lo(2)-1, lo(3)], [hi(1), hi(2)+1, hi(3)+1]
      const Box& tzxbx = amrex::grow(zbx, IntVect(AMREX_D_DECL(0,1,0)));

      qmzx.resize(tzxbx, NQ);
      Elixir elix_qmzx = qmzx.elixir();
      auto qmzx_arr = qmzx.array();
      fab_size += qmzx.nBytes();

      qpzx.resize(tzxbx, NQ);
      Elixir elix_qpzx = qpzx.elixir();
      auto qpzx_arr = qpzx.array();
      fab_size += qpzx.nBytes();

      trans_single(tzxbx, 0, 2,
                   qzm_arr, qmzx_arr,
                   qzp_arr, qpzx_arr,
                   qaux_arr,
                   ftmp1_arr,
#ifdef RADIATION
                   rftmp1_arr,
#endif
                   qgdnvtmp1_arr,
                   hdt, cdtdx);

      //reset_edge_state_thermo(tzxbx, qmzx.array());

      //reset_edge_state_thermo(tzxbx, qpzx.array());

      // compute F^y
      // [lo(1)-1, lo(2), lo(3)-1], [hi(1)+1, hi(2)+1, hi(3)+1]
      const Box& cybx = amrex::grow(ybx, IntVect(AMREX_D_DECL(1,0,1)));

      // ftmp1 = fy
      // rftmp1 = rfy
      // qgdnvtmp1 = qgdnvy
      cmpflx_plus_godunov(cybx,
                          qym_arr, qyp_arr,
                          ftmp1_arr,
#ifdef RADIATION
                          rftmp1_arr,
#endif
                          qgdnvtmp1_arr,
                          qaux_arr, shk_arr,
                          1, false);

      // [lo(1), lo(2), lo(3)-1], [hi(1)+1, hi(2), lo(3)+1]
      const Box& txybx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0,0,1)));

      qmxy.resize(txybx, NQ);
      Elixir elix_qmxy = qmxy.elixir();
      auto qmxy_arr = qmxy.array();
      fab_size += qmxy.nBytes();

      qpxy.resize(txybx, NQ);
      Elixir elix_qpxy = qpxy.elixir();
      auto qpxy_arr = qpxy.array();
      fab_size += qpxy.nBytes();

      // ftmp1 = fy
      // rftmp1 = rfy
      // qgdnvtmp1 = qgdnvy
      trans_single(txybx, 1, 0,
                   qxm_arr, qmxy_arr,
                   qxp_arr, qpxy_arr,
                   qaux_arr,
                   ftmp1_arr,
#ifdef RADIATION
                   rftmp1_arr,
#endif
                   qgdnvtmp1_arr,
                   hdt, cdtdy);

      //reset_edge_state_thermo(txybx, qmxy.array());

      //reset_edge_state_thermo(txybx, qpxy.array());

      // [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2), lo(3)+1]
      const Box& tzybx = amrex::grow(zbx, IntVect(AMREX_D_DECL(1,0,0)));

      qmzy.resize(tzybx, NQ);
      Elixir elix_qmzy = qmzy.elixir();
      auto qmzy_arr = qmzy.array();
      fab_size += qmzy.nBytes();

      qpzy.resize(tzybx, NQ);
      Elixir elix_qpzy = qpzy.elixir();
      auto qpzy_arr = qpzy.array();
      fab_size += qpzy.nBytes();

      // ftmp1 = fy
      // rftmp1 = rfy
      // qgdnvtmp1 = qgdnvy
      trans_single(tzybx, 1, 2,
                   qzm_arr, qmzy_arr,
                   qzp_arr, qpzy_arr,
                   qaux_arr,
                   ftmp1_arr,
#ifdef RADIATION
                   rftmp1_arr,
#endif
                   qgdnvtmp1_arr,
                   hdt, cdtdy);

      //reset_edge_state_thermo(tzybx, qmzy.array());

      //reset_edge_state_thermo(tzybx, qpzy.array());

      // compute F^z
      // [lo(1)-1, lo(2)-1, lo(3)], [hi(1)+1, hi(2)+1, hi(3)+1]
      const Box& czbx = amrex::grow(zbx, IntVect(AMREX_D_DECL(1,1,0)));

      // ftmp1 = fz
      // rftmp1 = rfz
      // qgdnvtmp1 = qgdnvz
      cmpflx_plus_godunov(czbx,
                          qzm_arr, qzp_arr,
                          ftmp1_arr,
#ifdef RADIATION
                          rftmp1_arr,
#endif
                          qgdnvtmp1_arr,
                          qaux_arr, shk_arr,
                          2, false);

      // [lo(1)-1, lo(2)-1, lo(3)], [hi(1)+1, hi(2)+1, lo(3)]
      const Box& txzbx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0,1,0)));

      qmxz.resize(txzbx, NQ);
      Elixir elix_qmxz = qmxz.elixir();
      auto qmxz_arr = qmxz.array();
      fab_size += qmxz.nBytes();

      qpxz.resize(txzbx, NQ);
      Elixir elix_qpxz = qpxz.elixir();
      auto qpxz_arr = qpxz.array();
      fab_size += qpxz.nBytes();

      // ftmp1 = fz
      // rftmp1 = rfz
      // qgdnvtmp1 = qgdnvz
      trans_single(txzbx, 2, 0,
                   qxm_arr, qmxz_arr,
                   qxp_arr, qpxz_arr,
                   qaux_arr,
                   ftmp1_arr,
#ifdef RADIATION
                   rftmp1_arr,
#endif
                   qgdnvtmp1_arr,
                   hdt, cdtdz);

      //reset_edge_state_thermo(txzbx, qmxz.array());

      //reset_edge_state_thermo(txzbx, qpxz.array());

      // [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2)+1, lo(3)]
      const Box& tyzbx = amrex::grow(ybx, IntVect(AMREX_D_DECL(1,0,0)));

      qmyz.resize(tyzbx, NQ);
      Elixir elix_qmyz = qmyz.elixir();
      auto qmyz_arr = qmyz.array();
      fab_size += qmyz.nBytes();

      qpyz.resize(tyzbx, NQ);
      Elixir elix_qpyz = qpyz.elixir();
      auto qpyz_arr = qpyz.array();
      fab_size += qpyz.nBytes();

      // ftmp1 = fz
      // rftmp1 = rfz
      // qgdnvtmp1 = qgdnvz
      trans_single(tyzbx, 2, 1,
                   qym_arr, qmyz_arr,
                   qyp_arr, qpyz_arr,
                   qaux_arr,
                   ftmp1_arr,
#ifdef RADIATION
                   rftmp1_arr,
#endif
                   qgdnvtmp1_arr,
                   hdt, cdtdz);

      //reset_edge_state_thermo(tyzbx, qmyz.array());

      //reset_edge_state_thermo(tyzbx, qpyz.array());

      // we now have q?zx, q?yx, q?zy, q?xy, q?yz, q?xz

      //
      // Use qx?, q?yz, q?zy to compute final x-flux
      //

      // compute F^{y|z}
      // [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2)+1, hi(3)]
      const Box& cyzbx = amrex::grow(ybx, IntVect(AMREX_D_DECL(1,0,0)));

      // ftmp1 = fyz
      // rftmp1 = rfyz
      // qgdnvtmp1 = qgdnvyz
      cmpflx_plus_godunov(cyzbx,
                          qmyz_arr, qpyz_arr,
                          ftmp1_arr,
#ifdef RADIATION
                          rftmp1_arr,
#endif
                          qgdnvtmp1_arr,
                          qaux_arr, shk_arr,
                          1, false);

      // compute F^{z|y}
      // [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)+1]
      const Box& czybx = amrex::grow(zbx, IntVect(AMREX_D_DECL(1,0,0)));

      // ftmp2 = fzy
      // rftmp2 = rfzy
      // qgdnvtmp2 = qgdnvzy
      cmpflx_plus_godunov(czybx,
                          qmzy_arr, qpzy_arr,
                          ftmp2_arr,
#ifdef RADIATION
                          rftmp2_arr,
#endif
                          qgdnvtmp2_arr,
                          qaux_arr, shk_arr,
                          2, false);

      // compute the corrected x interface states and fluxes
      // [lo(1), lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)]

      trans_final(xbx, 0, 1, 2,
                  qxm_arr, ql_arr,
                  qxp_arr, qr_arr,
                  qaux_arr,
                  ftmp1_arr,
#ifdef RADIATION
                  rftmp1_arr,
#endif
                  ftmp2_arr,
#ifdef RADIATION
                  rftmp2_arr,
#endif
                  qgdnvtmp1_arr,
                  qgdnvtmp2_arr,
                  hdtdy, hdtdz);

      //reset_edge_state_thermo(xbx, ql.array());

      //reset_edge_state_thermo(xbx, qr.array());

#ifdef PRIM_SPECIES_HAVE_SOURCES
      add_species_source_to_states(xbx, 0, dt,
                                   ql_arr, qr_arr, src_q_arr);

#endif


      cmpflx_plus_godunov(xbx,
                          ql_arr, qr_arr,
                          flux0_arr,
#ifdef RADIATION
                          rad_flux0_arr,
#endif
                          qex_arr,
                          qaux_arr, shk_arr,
                          0, false);

      //
      // Use qy?, q?zx, q?xz to compute final y-flux
      //

      // compute F^{z|x}
      // [lo(1), lo(2)-1, lo(3)], [hi(1), hi(2)+1, hi(3)+1]
      const Box& czxbx = amrex::grow(zbx, IntVect(AMREX_D_DECL(0,1,0)));

      // ftmp1 = fzx
      // rftmp1 = rfzx
      // qgdnvtmp1 = qgdnvzx
      cmpflx_plus_godunov(czxbx,
                          qmzx_arr, qpzx_arr,
                          ftmp1_arr,
#ifdef RADIATION
                          rftmp1_arr,
#endif
                          qgdnvtmp1_arr,
                          qaux_arr, shk_arr,
                          2, false);

      // compute F^{x|z}
      // [lo(1), lo(2)-1, lo(3)], [hi(1)+1, hi(2)+1, hi(3)]
      const Box& cxzbx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0,1,0)));

      // ftmp2 = fxz
      // rftmp2 = rfxz
      // qgdnvtmp2 = qgdnvxz
      cmpflx_plus_godunov(cxzbx,
                          qmxz_arr, qpxz_arr,
                          ftmp2_arr,
#ifdef RADIATION
                          rftmp2_arr,
#endif
                          qgdnvtmp2_arr,
                          qaux_arr, shk_arr,
                          0, false);

      // Compute the corrected y interface states and fluxes
      // [lo(1), lo(2), lo(3)], [hi(1), hi(2)+1, hi(3)]

      trans_final(ybx, 1, 0, 2,
                  qym_arr, ql_arr,
                  qyp_arr, qr_arr,
                  qaux_arr,
                  ftmp2_arr,
#ifdef RADIATION
                  rftmp2_arr,
#endif
                  ftmp1_arr,
#ifdef RADIATION
                  rftmp1_arr,
#endif
                  qgdnvtmp2_arr,
                  qgdnvtmp1_arr,
                  hdtdx, hdtdz);

      //reset_edge_state_thermo(ybx, ql.array());

      //reset_edge_state_thermo(ybx, qr.array());

#ifdef PRIM_SPECIES_HAVE_SOURCES
      add_species_source_to_states(ybx, 1, dt,
                                   ql_arr, qr_arr, src_q_arr);

#endif


      // Compute the final F^y
      // [lo(1), lo(2), lo(3)], [hi(1), hi(2)+1, hi(3)]
      cmpflx_plus_godunov(ybx,
                          ql_arr, qr_arr,
                          flux1_arr,
#ifdef RADIATION
                          rad_flux1_arr,
#endif
                          qey_arr,
                          qaux_arr, shk_arr,
                          1, false);

      //
      // Use qz?, q?xy, q?yx to compute final z-flux
      //

      // compute F^{x|y}
      // [lo(1), lo(2), lo(3)-1], [hi(1)+1, hi(2), hi(3)+1]
      const Box& cxybx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0,0,1)));

      // ftmp1 = fxy
      // rftmp1 = rfxy
      // qgdnvtmp1 = qgdnvxy
      cmpflx_plus_godunov(cxybx,
                          qmxy_arr, qpxy_arr,
                          ftmp1_arr,
#ifdef RADIATION
                          rftmp1_arr,
#endif
                          qgdnvtmp1_arr,
                          qaux_arr, shk_arr,
                          0, false);

      // compute F^{y|x}
      // [lo(1), lo(2), lo(3)-1], [hi(1), hi(2)+dg(2), hi(3)+1]
      const Box& cyxbx = amrex::grow(ybx, IntVect(AMREX_D_DECL(0,0,1)));

      // ftmp2 = fyx
      // rftmp2 = rfyx
      // qgdnvtmp2 = qgdnvyx
      cmpflx_plus_godunov(cyxbx,
                          qmyx_arr, qpyx_arr,
                          ftmp2_arr,
#ifdef RADIATION
                          rftmp2_arr,
#endif
                          qgdnvtmp2_arr,
                          qaux_arr, shk_arr,
                          1, false);

      // compute the corrected z interface states and fluxes
      // [lo(1), lo(2), lo(3)], [hi(1), hi(2), hi(3)+1]

      trans_final(zbx, 2, 0, 1,
                  qzm_arr, ql_arr,
                  qzp_arr, qr_arr,
                  qaux_arr,
                  ftmp1_arr,
#ifdef RADIATION
                  rftmp1_arr,
#endif
                  ftmp2_arr,
#ifdef RADIATION
                  rftmp2_arr,
#endif
                  qgdnvtmp1_arr,
                  qgdnvtmp2_arr,
                  hdtdx, hdtdy);

      //reset_edge_state_thermo(zbx, ql.array());

      //reset_edge_state_thermo(zbx, qr.array());

#ifdef PRIM_SPECIES_HAVE_SOURCES
      add_species_source_to_states(zbx, 2, dt,
                                   ql_arr, qr_arr, src_q_arr);

#endif

      // compute the final z fluxes F^z
      // [lo(1), lo(2), lo(3)], [hi(1), hi(2), hi(3)+1]

      cmpflx_plus_godunov(zbx,
                          ql_arr, qr_arr,
                          flux2_arr,
#ifdef RADIATION
                          rad_flux2_arr,
#endif
                          qez_arr,
                          qaux_arr, shk_arr,
                          2, false);

#endif // 3-d



      // clean the fluxes

      for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

          const Box& nbx = amrex::surroundingNodes(bx, idir);

          Array4<Real> const flux_arr = (flux[idir]).array();
          Array4<Real const> const uin_arr = Sborder.array(mfi);

          // Zero out shock and temp fluxes -- these are physically meaningless here
          amrex::ParallelFor(nbx,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
          {
              flux_arr(i,j,k,UTEMP) = 0.e0;
#ifdef SHOCK_VAR
              flux_arr(i,j,k,USHK) = 0.e0;
#endif
          });

          apply_av(nbx, idir, div_arr, uin_arr, flux_arr);

#ifdef RADIATION
          Array4<Real> const rad_flux_arr = (rad_flux[idir]).array();
          Array4<Real const> const Erin_arr = Erborder.array(mfi);

          apply_av_rad(nbx, idir, div_arr, Erin_arr, rad_flux_arr);
#endif

          if (limit_fluxes_on_small_dens == 1) {
              limit_hydro_fluxes_on_small_dens
                  (nbx, idir,
                   Sborder.array(mfi),
                   q.array(),
                   volume.array(mfi),
                   flux[idir].array(),
                   area[idir].array(mfi),
                   dt);
          }

          if (limit_fluxes_on_large_vel == 1) {
              limit_hydro_fluxes_on_large_vel
                  (nbx, idir,
                   Sborder.array(mfi),
                   q.array(),
                   volume.array(mfi),
                   flux[idir].array(),
                   area[idir].array(mfi),
                   dt);
          }

          normalize_species_fluxes(nbx, flux_arr);

      }



      // conservative update
      Array4<Real> const update_arr = S_new.array(mfi);

      Array4<Real> const flx_arr = (flux[0]).array();
      Array4<Real> const qx_arr = (qe[0]).array();

#if AMREX_SPACEDIM >= 2
      Array4<Real> const fly_arr = (flux[1]).array();
      Array4<Real> const qy_arr = (qe[1]).array();
#endif

#if AMREX_SPACEDIM == 3
      Array4<Real> const flz_arr = (flux[2]).array();
      Array4<Real> const qz_arr = (qe[2]).array();
#endif

      consup_hydro(bx,
#ifdef SHOCK_VAR
                   shk_arr,
#endif
                   update_arr,
                   flx_arr, qx_arr,
#if AMREX_SPACEDIM >= 2
                   fly_arr, qy_arr,
#endif
#if AMREX_SPACEDIM == 3
                   flz_arr, qz_arr,
#endif
                   dt);


#ifdef HYBRID_MOMENTUM
      auto dx_arr = geom.CellSizeArray();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
      {

          GpuArray<Real, 3> loc;

          position(i, j, k, geomdata, loc);

          for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
              loc[dir] -= problem::center[dir];

          Real R = amrex::max(std::sqrt(loc[0] * loc[0] + loc[1] * loc[1]), R_min);
          Real RInv = 1.0_rt / R;

          update_arr(i,j,k,UMR) -= dt * (loc[0] * RInv) * (qx_arr(i+1,j,k,GDPRES) - qx_arr(i,j,k,GDPRES)) / dx_arr[0];
#if AMREX_SPACEDIM >= 2
          update_arr(i,j,k,UMR) -= dt * (loc[1] * RInv) * (qy_arr(i,j+1,k,GDPRES) - qy_arr(i,j,k,GDPRES)) / dx_arr[1];
#endif
      });
#endif

#ifdef RADIATION
      ctu_rad_consup(bx,
                     update_arr,
                     Erborder.array(mfi),
                     Er_new.array(mfi),
                     (rad_flux[0]).array(),
                     (qe[0]).array(),
                     (area[0]).array(mfi),
#if AMREX_SPACEDIM >= 2
                     (rad_flux[1]).array(),
                     (qe[1]).array(),
                     (area[1]).array(mfi),
#endif
#if AMREX_SPACEDIM == 3
                     (rad_flux[2]).array(),
                     (qe[2]).array(),
                     (area[2]).array(mfi),
#endif
                     priv_nstep_fsp,
                     volume.array(mfi),
                     dt);

      nstep_fsp = std::max(nstep_fsp, priv_nstep_fsp);
#endif


      for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

        const Box& nbx = amrex::surroundingNodes(bx, idir);

        Array4<Real> const flux_arr = (flux[idir]).array();
        Array4<Real const> const area_arr = (area[idir]).array(mfi);

        scale_flux(nbx,
#if AMREX_SPACEDIM == 1
                   qex_arr,
#endif
                   flux_arr, area_arr, dt);

#ifdef RADIATION
        Array4<Real> const rad_flux_arr = (rad_flux[idir]).array();
        scale_rad_flux(nbx, rad_flux_arr, area_arr, dt);
#endif

        if (idir == 0) {
#if AMREX_SPACEDIM <= 2
            Array4<Real> pradial_fab = pradial.array();
#endif

            // get the scaled radial pressure -- we need to treat this specially
#if AMREX_SPACEDIM <= 2

#if AMREX_SPACEDIM == 1
            if (!Geom().IsCartesian()) {
#elif AMREX_SPACEDIM == 2
            if (!mom_flux_has_p(0, 0, coord)) {
#endif
                amrex::ParallelFor(nbx,
                [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
                {
                    pradial_fab(i,j,k) = qex_arr(i,j,k,GDPRES) * dt;
                });
            }

#endif
        }

        // Store the fluxes from this advance. For simplified SDC integration we
        // only need to do this on the last iteration.

        bool add_fluxes = true;

        if (time_integration_method == SimplifiedSpectralDeferredCorrections &&
            sdc_iteration != sdc_iters - 1) {
            add_fluxes = false;
        }

        if (add_fluxes) {

            Array4<Real> const flux_fab = (flux[idir]).array();
            Array4<Real> fluxes_fab = (*fluxes[idir]).array(mfi);

            amrex::ParallelFor(mfi.nodaltilebox(idir), NUM_STATE,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k, int n)
            {
                fluxes_fab(i,j,k,n) += flux_fab(i,j,k,n);
            });

#ifdef RADIATION
            Array4<Real> const rad_flux_fab = (rad_flux[idir]).array();
            Array4<Real> rad_fluxes_fab = (*rad_fluxes[idir]).array(mfi);

            amrex::ParallelFor(mfi.nodaltilebox(idir), Radiation::nGroups,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k, int n)
            {
                rad_fluxes_fab(i,j,k,n) += rad_flux_fab(i,j,k,n);
            });
#endif

#if AMREX_SPACEDIM <= 2

#if AMREX_SPACEDIM == 1
            if (idir == 0 && !Geom().IsCartesian()) {
#elif AMREX_SPACEDIM == 2
            if (idir == 0 && !mom_flux_has_p(0, 0, coord)) {
#endif
                Array4<Real> pradial_fab = pradial.array();
                Array4<Real> P_radial_fab = P_radial.array(mfi);

                amrex::ParallelFor(mfi.nodaltilebox(0),
                [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
                {
                    P_radial_fab(i,j,k,0) += pradial_fab(i,j,k,0);
                });
            }

#endif

        } // add_fluxes

        Array4<Real> const flux_fab = (flux[idir]).array();
        Array4<Real> mass_fluxes_fab = (*mass_fluxes[idir]).array(mfi);

        amrex::ParallelFor(mfi.nodaltilebox(idir),
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {
            // This is a copy, not an add, since we need mass_fluxes to be
            // only this subcycle's data when we evaluate the gravitational
            // forces.

            mass_fluxes_fab(i,j,k,0) = flux_fab(i,j,k,URHO);
        });

      } // idir loop

#ifdef AMREX_USE_GPU
      // Check if we're going to run out of memory in the next MFIter iteration.
      // If so, do a synchronize here so that we don't oversubscribe GPU memory.
      // Note that this will capture the case where we started with more memory
      // than what the GPU has, on the logic that even in that case, it makes
      // sense to not further pile on the oversubscription demands.

      // This could (and should) be generalized in the future to operate with
      // more granularity than the MFIter loop boundary. We would have potential
      // synchronization points prior to each of the above kernel launches, and
      // we would check whether the sum of all previously allocated fabs would
      // result in oversubscription, including any contributions from a partial
      // MFIter loop. A further optimization would be to not apply a device
      // synchronize, but rather to use CUDA events to poll on a check about
      // whether enough memory has freed up to begin the next iteration, and then
      // immediately proceed to the next kernel when there's enough space for it.

      current_size += fab_size;
      if (current_size + fab_size >= Gpu::Device::totalGlobalMem()) {
          Gpu::Device::synchronize();
          current_size = starting_size;
      }
#endif

      if (oversubscribed) {
          volume[mfi].prefetchToHost();
          Sborder[mfi].prefetchToHost();
          S_new[mfi].prefetchToHost();
          for (int i = 0; i < AMREX_SPACEDIM; ++i) {
              area[i][mfi].prefetchToHost();
              (*fluxes[i])[mfi].prefetchToHost();
          }
#if AMREX_SPACEDIM < 3
          dLogArea[0][mfi].prefetchToHost();
          P_radial[mfi].prefetchToHost();
#endif
#ifdef RADIATION
          Erborder[mfi].prefetchToHost();
          Er_new[mfi].prefetchToHost();
#endif
      }


    } // MFIter loop

  } // OMP loop

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
#endif

  if (verbose && ParallelDescriptor::IOProcessor())
    std::cout << "... Leaving construct_ctu_hydro_source()" << std::endl << std::endl;

  if (verbose > 0)
    {
      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
      Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
      Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        if (ParallelDescriptor::IOProcessor())
          std::cout << "Castro::construct_ctu_hydro_source() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
        });
#endif
    }

#endif

}
