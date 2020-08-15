#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void
Castro::just_the_mhd(Real time, Real dt)
{
      if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... mhd ...!!! " << std::endl << std::endl;

      hydro_source.setVal(0.0);

      const auto dx = geom.CellSizeArray();

      MultiFab& S_new = get_new_data(State_Type);
      MultiFab& Bx_new= get_new_data(Mag_Type_x);
      MultiFab& By_new= get_new_data(Mag_Type_y);
      MultiFab& Bz_new= get_new_data(Mag_Type_z);


      //MultiFab electric[BL_SPACEDIM];
      //for (int j = 0; j < BL_SPACEDIM; j++)
      //{
      //  electric[j].define(getEdgeBoxArray(j), dmap, 1, 0);
      //  electric[j].setVal(0.0);
      //}


      BL_ASSERT(NUM_GROW == 6);


#ifdef _OPENMP
#pragma omp parallel
#endif
    {

      FArrayBox flux[AMREX_SPACEDIM], E[AMREX_SPACEDIM];

      FArrayBox q;
      FArrayBox qaux;
      FArrayBox srcQ;

      FArrayBox flatn;
      FArrayBox flatg;

      FArrayBox qleft[AMREX_SPACEDIM];
      FArrayBox qright[AMREX_SPACEDIM];

      FArrayBox flxx1D;
      FArrayBox flxy1D;
      FArrayBox flxz1D;

      FArrayBox ux_left;
      FArrayBox ux_right;
      FArrayBox uy_left;
      FArrayBox uy_right;
      FArrayBox uz_left;
      FArrayBox uz_right;

      FArrayBox qtmp_left;
      FArrayBox qtmp_right;

      FArrayBox flx_xy;
      FArrayBox flx_xz;

      FArrayBox flx_yx;
      FArrayBox flx_yz;

      FArrayBox flx_zx;
      FArrayBox flx_zy;

      FArrayBox q2D;

      FArrayBox div;

      for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {

          const Box& bx = mfi.tilebox();
          const Box& obx = amrex::grow(bx, 1);
          const Box& gbx = amrex::grow(bx, 2);

          // box with NUM_GROW ghost cells for PPM stuff
          const Box& bx_gc = amrex::grow(bx, NUM_GROW);

          FArrayBox &statein  = Sborder[mfi];
          auto u_arr = statein.array();

          FArrayBox &source_in  = sources_for_hydro[mfi];
          auto src_arr = source_in.array();

          FArrayBox &hydro_update = hydro_source[mfi];
          auto update_arr = hydro_update.array();

          FArrayBox& Bx  = Bx_old_tmp[mfi];
          auto Bx_arr = Bx.array();
          auto elix_Bx = Bx.elixir();

          FArrayBox& By  = By_old_tmp[mfi];
          auto By_arr = By.array();
          auto elix_By = By.elixir();

          FArrayBox& Bz  = Bz_old_tmp[mfi];
          auto Bz_arr = Bz.array();
          auto elix_Bz = Bz.elixir();

          FArrayBox& Bxout = Bx_new[mfi];
          auto Bxo_arr = Bxout.array();

          FArrayBox& Byout = By_new[mfi];
          auto Byo_arr = Byout.array();

          FArrayBox& Bzout = Bz_new[mfi];
          auto Bzo_arr = Bzout.array();


          // allocate fabs for fluxes and electric field

          const Box& nbx = amrex::surroundingNodes(bx, 0);
          const Box& nby = amrex::surroundingNodes(bx, 1);
          const Box& nbz = amrex::surroundingNodes(bx, 2);

          const Box& nbxf = amrex::grow(nbx, IntVect(0, 1, 1));
          const Box& nbyf = amrex::grow(nby, IntVect(1, 0, 1));
          const Box& nbzf = amrex::grow(nbz, IntVect(1, 1, 0));

          // need to revisit these box sizes
          const Box& nbxe = amrex::grow(nbx, IntVect(2, 3, 3));
          const Box& nbye = amrex::grow(nby, IntVect(3, 2, 3));
          const Box& nbze = amrex::grow(nbz, IntVect(3, 3, 2));

          flux[0].resize(nbxf, NUM_STATE+3);
          auto flxx_arr = flux[0].array();
          auto elix_flxx = flux[0].elixir();

          E[0].resize(nbxe);
          auto Ex_arr = E[0].array();
          auto elix_Ex = E[0].elixir();

          flux[1].resize(nbyf, NUM_STATE+3);
          auto flxy_arr = flux[1].array();
          auto elix_flxy = flux[1].elixir();

          E[1].resize(nbye);
          auto Ey_arr = E[1].array();
          auto elix_Ey = E[1].elixir();

          flux[2].resize(nbzf, NUM_STATE+3);
          auto flxz_arr = flux[2].array();
          auto elix_flxz = flux[2].elixir();

          E[2].resize(nbze);
          auto Ez_arr = E[2].array();
          auto elix_Ez = E[2].elixir();


          // Calculate primitives based on conservatives
          q.resize(bx_gc, NQ);
          auto q_arr = q.array();
          auto elix_q = q.elixir();

          qaux.resize(bx_gc, NQAUX);
          auto qaux_arr = qaux.array();
          auto elix_qaux = qaux.elixir();

          srcQ.resize(bx_gc, NQSRC);
          auto src_q_arr = srcQ.array();
          auto elix_src_q = srcQ.elixir();

          ctoprim(bx_gc, time,
                  u_arr,
                  Bx_arr, By_arr, Bz_arr,
                  q_arr, qaux_arr);

          src_to_prim(bx_gc, q_arr, src_arr, src_q_arr);

          check_for_mhd_cfl_violation(bx, dt, q_arr, qaux_arr);

          // we need to compute the flattening coefficient for every zone
          // center where we do reconstruction

          const Box& bxi = amrex::grow(bx, IntVect(3, 3, 3));

          flatn.resize(bxi, 1);
          auto flatn_arr = flatn.array();
          auto elix_flatn = flatn.elixir();

          flatg.resize(bxi, 1);
          auto flatg_arr = flatg.array();
          auto elix_flatg = flatg.elixir();

          if (use_flattening == 0) {
            amrex::ParallelFor(bxi,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
            {
              flatn_arr(i,j,k) = 1.0;
            });

          } else {

            uflatten(bxi, q_arr, flatn_arr, QPRES);
            uflatten(bxi, q_arr, flatg_arr, QPTOT);

            amrex::ParallelFor(bxi,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
            {
              flatn_arr(i,j,k) = flatn_arr(i,j,k) * flatg_arr(i,j,k);
            });

          }

          // Interpolate Cell centered values to faces
          qleft[0].resize(bx_gc, NQ);
          auto qx_left_arr = qleft[0].array();
          auto elix_qx_left = qleft[0].elixir();

          qright[0].resize(bx_gc, NQ);
          auto qx_right_arr = qright[0].array();
          auto elix_qx_right = qright[0].elixir();

          qleft[1].resize(bx_gc, NQ);
          auto qy_left_arr = qleft[1].array();
          auto elix_qy_left = qleft[1].elixir();

          qright[1].resize(bx_gc, NQ);
          auto qy_right_arr = qright[1].array();
          auto elix_qy_right = qright[1].elixir();

          qleft[2].resize(bx_gc, NQ);
          auto qz_left_arr = qleft[2].array();
          auto elix_qz_left = qleft[2].elixir();

          qright[2].resize(bx_gc, NQ);
          auto qz_right_arr = qright[2].array();
          auto elix_qz_right = qright[2].elixir();


          for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {

            if (ppm_type == 0) {
              plm(bxi, idir,
                  q_arr, qaux_arr, flatn_arr,
                  Bx_arr, By_arr, Bz_arr,
                  qleft[idir].array(), qright[idir].array(),
                  src_q_arr, dt);

            } else {
              ppm_mhd(bxi, idir,
                      q_arr, qaux_arr, flatn_arr,
                      Bx_arr, By_arr, Bz_arr,
                      qleft[idir].array(), qright[idir].array(),
                      src_q_arr, dt);
            }
          }

          // Corner Couple and find the correct fluxes + electric fields

          // Do the corner coupling and the CT updates

          // MM CTU Step 1
          // Calculate Flux 1D, eq.35

          // x-dir
          // [lo(1)-2, lo(2)-3, lo(3)-3] [hi(1)+3, hi(2)+3, hi(3)+3]
          const Box& bfx = amrex::grow(nbx, IntVect(2, 3, 3));

          flxx1D.resize(bfx, NUM_STATE+3);
          auto flxx1D_arr = flxx1D.array();
          auto elix_flxx1D = flxx1D.elixir();

          hlld(bfx, qleft[0].array(), qright[0].array(), flxx1D_arr, 0);

          // y-dir
          // [lo(1)-3, lo(2)-2, lo(3)-3] [hi(1)+3, hi(2)+3, hi(3)+3]
          const Box& bfy = amrex::grow(nby, IntVect(3, 2, 3));

          flxy1D.resize(bfy, NUM_STATE+3);
          auto flxy1D_arr = flxy1D.array();
          auto elix_flxy1D = flxy1D.elixir();

          hlld(bfy, qleft[1].array(), qright[1].array(), flxy1D_arr, 1);

          // z-dir
          // [lo(1)-3, lo(2)-3, lo(3)-2] [hi(1)+3, hi(2)+3, hi(3)+3]
          const Box& bfz = amrex::grow(nbz, IntVect(3, 3, 2));

          flxz1D.resize(bfz, NUM_STATE+3);
          auto flxz1D_arr = flxz1D.array();
          auto elix_flxz1D = flxz1D.elixir();

          hlld(bfz, qleft[2].array(), qright[2].array(), flxz1D_arr, 2);


          // Prim to Cons

          ux_left.resize(gbx, NUM_STATE+3);
          auto ux_left_arr = ux_left.array();
          auto elix_ux_left = ux_left.elixir();

          ux_right.resize(gbx, NUM_STATE+3);
          auto ux_right_arr = ux_right.array();
          auto elix_ux_right = ux_right.elixir();

          PrimToCons(gbx, qx_left_arr, ux_left_arr);
          PrimToCons(gbx, qx_right_arr, ux_right_arr);

          uy_left.resize(gbx, NUM_STATE+3);
          auto uy_left_arr = uy_left.array();
          auto elix_uy_left = uy_left.elixir();

          uy_right.resize(gbx, NUM_STATE+3);
          auto uy_right_arr = uy_right.array();
          auto elix_uy_right = uy_right.elixir();

          PrimToCons(gbx, qy_left_arr, uy_left_arr);
          PrimToCons(gbx, qy_right_arr, uy_right_arr);

          uz_left.resize(gbx, NUM_STATE+3);
          auto uz_left_arr = uz_left.array();
          auto elix_uz_left = uz_left.elixir();

          uz_right.resize(gbx, NUM_STATE+3);
          auto uz_right_arr = uz_right.array();
          auto elix_uz_right = uz_right.elixir();

          PrimToCons(gbx, qz_left_arr, uz_left_arr);
          PrimToCons(gbx, qz_right_arr, uz_right_arr);

          // MM CTU Step 2
          // Use "1D" fluxes To interpolate Temporary Edge Centered Electric Fields, eq.36

          // [lo(1)-2, lo(2)-2, lo(3)-2][hi(1)+2, hi(2)+3, hi(3)+3]
          Box eebx = amrex::grow(bx, 2);
          eebx.growHi(1);
          eebx.growHi(2);

          electric_edge_x(eebx, q_arr, Ex_arr, flxy1D_arr, flxz1D_arr);


          // [lo(1)-2, lo(2)-2, lo(3)-2][hi(1)+3, hi(2)+2, hi(3)+3]
          Box eeby = amrex::grow(bx, 2);
          eeby.growHi(0);
          eeby.growHi(2);

          electric_edge_y(eeby, q_arr, Ey_arr, flxx1D_arr, flxz1D_arr);


          // [lo(1)-2, lo(2)-2, lo(3)-2][hi(1)+3, hi(2)+3, hi(3)+2]
          Box eebz = amrex::grow(bx, 2);
          eebz.growHi(0);
          eebz.growHi(1);

          electric_edge_z(eebz, q_arr, Ez_arr, flxx1D_arr, flxy1D_arr);


          // MM CTU Steps 3, 4, and 5
          // Corner Couple, eq. 37, 38 and 39 Correct Conservative vars using Transverse Fluxes

          // X direction

          // affected by Y Flux
          // [lo(1)-1, lo(2)-2, lo(3)-2] [hi(1)+2, hi(2)+2, hi(2)+2]
          const Box& ccbx = amrex::grow(nbx, IntVect(1, 2, 2));

          qtmp_left.resize(gbx, NQ);
          auto qtmp_left_arr = qtmp_left.array();
          auto elix_qtmp_left = qtmp_left.elixir();

          qtmp_right.resize(gbx, NQ);
          auto qtmp_right_arr = qtmp_right.array();
          auto elix_qtmp_right = qtmp_right.elixir();

          corner_couple(ccbx,
                        qtmp_right_arr, qtmp_left_arr,
                        ux_right_arr, ux_left_arr,
                        flxy1D_arr, Ex_arr, Ez_arr,
                        Ex_arr, Ey_arr, Ez_arr,
                        0, 1, 2, dt);

          // Calculate Flux 2D eq. 40
          // F^{x|y}
          flx_xy.resize(ccbx, NUM_STATE+3);
          auto flx_xy_arr = flx_xy.array();
          auto elix_flx_xy = flx_xy.elixir();

          hlld(ccbx, qtmp_left_arr, qtmp_right_arr, flx_xy_arr, 0);

          // affected by Z Flux
          corner_couple(ccbx,
                        qtmp_right_arr, qtmp_left_arr,
                        ux_right_arr, ux_left_arr,
                        flxz1D_arr, Ex_arr, Ey_arr,
                        Ex_arr, Ey_arr, Ez_arr,
                        0, 2, 1, dt);

          // F^{x|z}
          flx_xz.resize(ccbx, NUM_STATE+3);
          auto flx_xz_arr = flx_xz.array();
          auto elix_flx_xz = flx_xz.elixir();

          hlld(ccbx, qtmp_left_arr, qtmp_right_arr, flx_xz_arr, 0);


          // Y direction

          // affected by X Flux
          // [lo(1)-2, lo(2)-1, lo(3)-2] [hi(1)+2, hi(2)+2, hi(3)+2]
          const Box& ccby = amrex::grow(nby, IntVect(2, 1, 2));

          corner_couple(ccby,
                        qtmp_right_arr, qtmp_left_arr,
                        uy_right_arr, uy_left_arr,
                        flxx1D_arr, Ey_arr, Ez_arr,
                        Ex_arr, Ey_arr, Ez_arr,
                        1, 0, 2, dt);

          // F^{y|x}
          flx_yx.resize(ccby, NUM_STATE+3);
          auto flx_yx_arr = flx_yx.array();
          auto elix_flx_yx = flx_yx.elixir();

          hlld(ccby, qtmp_left_arr, qtmp_right_arr, flx_yx_arr, 1);

          // affected by Z Flux

          corner_couple(ccby,
                        qtmp_right_arr, qtmp_left_arr,
                        uy_right_arr, uy_left_arr,
                        flxz1D_arr, Ey_arr, Ex_arr,
                        Ex_arr, Ey_arr, Ez_arr,
                        1, 2, 0, dt);

          // F^{y|z}
          flx_yz.resize(ccby, NUM_STATE+3);
          auto flx_yz_arr = flx_yz.array();
          auto elix_flx_yz = flx_yz.elixir();

          hlld(ccby, qtmp_left_arr, qtmp_right_arr, flx_yz_arr, 1);

          // Z direction

          // affected by X Flux
          // [lo(1)-2, lo(2)-2, lo(3)-1] [hi(1)+2, hi(2)+2, hi(3)+2]
          const Box& ccbz = amrex::grow(nbz, IntVect(2, 2, 1));

          corner_couple(ccbz,
                        qtmp_right_arr, qtmp_left_arr,
                        uz_right_arr, uz_left_arr,
                        flxx1D_arr, Ez_arr, Ey_arr,
                        Ex_arr, Ey_arr, Ez_arr,
                        2, 0, 1, dt);

          // F^{z|x}
          flx_zx.resize(ccbz, NUM_STATE+3);
          auto flx_zx_arr = flx_zx.array();
          auto elix_flx_zx = flx_zx.elixir();

          hlld(ccbz, qtmp_left_arr, qtmp_right_arr, flx_zx_arr, 2);

          // affected by Y Flux

          corner_couple(ccbz,
                        qtmp_right_arr, qtmp_left_arr,
                        uz_right_arr, uz_left_arr,
                        flxy1D_arr, Ez_arr, Ex_arr,
                        Ex_arr, Ey_arr, Ez_arr,
                        2, 1, 0, dt);

          // F^{z|y}
          flx_zy.resize(ccbz, NUM_STATE+3);
          auto flx_zy_arr = flx_zy.array();
          auto elix_flx_zy = flx_zy.elixir();

          hlld(ccbz, qtmp_left_arr, qtmp_right_arr, flx_zy_arr, 2);


          // MM CTU Step 6
          // Use Averaged 2D fluxes to interpolate temporary Edge Centered Electric Fields, reuse "flx1D"
          // eq. 42 and 43

          amrex::ParallelFor(ccbx, NUM_STATE+3,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k, int n) noexcept
          {
            flxx1D_arr(i,j,k,n) = 0.5_rt * (flx_xy_arr(i,j,k,n) + flx_xz_arr(i,j,k,n));
          });

          amrex::ParallelFor(ccby, NUM_STATE+3,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k, int n) noexcept
          {
            flxy1D_arr(i,j,k,n) = 0.5_rt * (flx_yx_arr(i,j,k,n) + flx_yz_arr(i,j,k,n));
          });

          amrex::ParallelFor(ccbz, NUM_STATE+3,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k, int n) noexcept
          {
            flxz1D_arr(i,j,k,n) = 0.5_rt * (flx_zx_arr(i,j,k,n) + flx_zy_arr(i,j,k,n));
          });


          // eq. 41
          // [lo(1)-1, lo(2)-1, lo(3)-1][hi(1)+1, hi(2)+2, hi(3)+2]
          Box eebx2 = amrex::grow(bx, 1);
          eebx2.growHi(1);
          eebx2.growHi(2);

          electric_edge_x(eebx2, q_arr, Ex_arr, flxy1D_arr, flxz1D_arr);


          // [lo(1)-1, lo(2)-1, lo(3)-1][hi(1)+2, hi(2)+1, hi(3)+2]
          Box eeby2 = amrex::grow(bx, 1);
          eeby2.growHi(0);
          eeby2.growHi(2);

          electric_edge_y(eeby2, q_arr, Ey_arr, flxx1D_arr, flxz1D_arr);


          // [lo(1)-1, lo(2)-1, lo(3)-1][hi(1)+2, hi(2)+2, hi(3)+1]
          Box eebz2 = amrex::grow(bx, 1);
          eebz2.growHi(0);
          eebz2.growHi(1);

          electric_edge_z(eebz2, q_arr, Ez_arr, flxx1D_arr, flxy1D_arr);


          // MM CTU Step 7, 8, and 9
          // Half Step conservative vars eq.44, eq.45, eq.46
          // Here we reuse qtmp_left/right to denote the half-time conservative state

          // for x direction
          // [lo(1), lo(2)-1, lo(3)-1][hi(1)+1, hi(2)+1, hi(3)+1]
          const Box& nbx1 = amrex::grow(nbx, IntVect(0, 1, 1));

          half_step(nbx1,
                    qtmp_right_arr, qtmp_left_arr,
                    ux_right_arr, ux_left_arr,
                    flx_yz_arr, flx_zy_arr,
                    Ex_arr, Ey_arr, Ez_arr,
                    Ex_arr, Ey_arr, Ez_arr,
                    0, 1, 2, dt);

          // Final Fluxes eq.47

          // We need to compute these on a box 1 larger in the transverse directions
          // than we'd need for hydro alone due to the electric update

          hlld(nbx1, qtmp_left_arr, qtmp_right_arr, flux[0].array(), 0);

          // for y direction
          const Box& nby1 = amrex::grow(nby, IntVect(1, 0, 1));

          half_step(nby1,
                    qtmp_right_arr, qtmp_left_arr,
                    uy_right_arr, uy_left_arr,
                    flx_xz_arr, flx_zx_arr,
                    Ey_arr, Ex_arr, Ez_arr,
                    Ex_arr, Ey_arr, Ez_arr,
                    1, 0, 2, dt);

          hlld(nby1, qtmp_left_arr, qtmp_right_arr, flux[1].array(), 1);

          // for z direction
          const Box& nbz1 = amrex::grow(nbz, IntVect(1, 1, 0));

          half_step(nbz1,
                    qtmp_right_arr, qtmp_left_arr,
                    uz_right_arr, uz_left_arr,
                    flx_xy_arr, flx_yx_arr,
                    Ez_arr, Ex_arr, Ey_arr,
                    Ex_arr, Ey_arr, Ez_arr,
                    2, 0, 1, dt);

          hlld(nbz1, qtmp_left_arr, qtmp_right_arr, flux[2].array(), 2);


          // MM CTU Step 10
          // Primitive update eq. 48
          q2D.resize(obx, NQ);
          auto q2D_arr = q2D.array();
          auto elix_q2D = q2D.elixir();

          prim_half(obx, q2D_arr, q_arr,
                    flxx1D_arr, flxy1D_arr, flxz1D_arr, dt);

          // Final Electric Field Update eq.48

          // [lo(1), lo(2), lo(3)][hi(1), hi(2)+1, hi(3)+1]
          Box eebxf = mfi.tilebox();
          eebxf.growHi(1, 1);
          eebxf.growHi(2, 1);

          electric_edge_x(eebxf, q2D_arr, Ex_arr, flxy_arr, flxz_arr);

          // [lo(1), lo(2), lo(3)][hi(1)+1, hi(2), hi(3)+1]
          Box eebyf = mfi.tilebox();
          eebyf.growHi(0, 1);
          eebyf.growHi(2, 1);

          electric_edge_y(eebyf, q2D_arr, Ey_arr, flxx_arr, flxz_arr);

          // [lo(1), lo(2), lo(3)][hi(1)+1, hi(2)+1 ,hi(3)]
          Box eebzf = mfi.tilebox();
          eebzf.growHi(0, 1);
          eebzf.growHi(1, 1);

          electric_edge_z(eebzf, q2D_arr, Ez_arr, flxx_arr, flxy_arr);

          // clean the final fluxes

          div.resize(obx, 1);
          Elixir elix_div = div.elixir();
          auto div_arr = div.array();

          // compute divu -- we'll use this later when doing the artifical viscosity
          divu(obx, q_arr, div_arr);

          for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

            const Box& nbx = amrex::surroundingNodes(bx, idir);

            Array4<Real> const flux_arr = (flux[idir]).array();

            // Zero out shock and temp fluxes -- these are physically meaningless here
            amrex::ParallelFor(nbx,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
            {
              flux_arr(i,j,k,UTEMP) = 0.e0;
#ifdef SHOCK_VAR
              flux_arr(i,j,k,USHK) = 0.e0;
#endif
            });

            apply_av(nbx, idir, div_arr, u_arr, flux_arr);

            normalize_species_fluxes(nbx, flux_arr);

          }


          // Conservative update

          consup_mhd(bx, update_arr, flxx_arr, flxy_arr, flxz_arr);

          // magnetic update

          Real dtdx = dt / dx[0];

          amrex::ParallelFor(nbx,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
          {
            Bxo_arr(i,j,k) = Bx_arr(i,j,k) + dtdx *
              ((Ey_arr(i,j,k+1) - Ey_arr(i,j,k)) - (Ez_arr(i,j+1,k) - Ez_arr(i,j,k)));
          });

#if AMREX_SPACEDIM >= 2
          dtdx = dt / dx[1];
#else
          dtdx = 0.0_rt;
#endif

          amrex::ParallelFor(nby,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
          {
            Byo_arr(i,j,k) = By_arr(i,j,k) + dtdx *
              ((Ez_arr(i+1,j,k) - Ez_arr(i,j,k)) - (Ex_arr(i,j,k+1) - Ex_arr(i,j,k)));
          });

#if AMREX_SPACEDIM == 3
          dtdx = dt / dx[2];
#else
          dtdx = 0.0_rt;
#endif

          amrex::ParallelFor(nbz,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
          {
            Bzo_arr(i,j,k) = Bz_arr(i,j,k) + dtdx *
              ((Ex_arr(i,j+1,k) - Ex_arr(i,j,k)) - (Ey_arr(i+1,j,k) - Ey_arr(i,j,k)));
          });


          // not sure if this is needed

          //Ex(ex_lo(1):ex_hi(1),ex_lo(2):ex_hi(2), ex_lo(3):ex_hi(3)) = Extemp(ex_lo(1):ex_hi(1),ex_lo(2):ex_hi(2),ex_lo(3):ex_hi(3))
          //Ey(ey_lo(1):ey_hi(1),ey_lo(2):ey_hi(2), ey_lo(3):ey_hi(3)) = Eytemp(ey_lo(1):ey_hi(1),ey_lo(2):ey_hi(2),ey_lo(3):ey_hi(3))
          //Ez(ez_lo(1):ez_hi(1),ez_lo(2):ez_hi(2), ez_lo(3):ez_hi(3)) = Eztemp(ez_lo(1):ez_hi(1),ez_lo(2):ez_hi(2),ez_lo(3):ez_hi(3))

          // Store the fluxes from this advance.

          // For normal integration we want to add the fluxes from this advance
          // since we may be subcycling the timestep. But for simplified SDC integration
          // we want to copy the fluxes since we expect that there will not be
          // subcycling and we only want the last iteration's fluxes.

          for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {

            Array4<Real> const flux_fab = (flux[idir]).array();
            Array4<Real> fluxes_fab = (*fluxes[idir]).array(mfi);
            const int numcomp = NUM_STATE;

            if (time_integration_method == SimplifiedSpectralDeferredCorrections) {

              AMREX_HOST_DEVICE_FOR_4D(mfi.nodaltilebox(idir), numcomp, i, j, k, n,
              {
                fluxes_fab(i,j,k,n) = flux_fab(i,j,k,n);
              });

            } else {

              AMREX_HOST_DEVICE_FOR_4D(mfi.nodaltilebox(idir), numcomp, i, j, k, n,
              {
                fluxes_fab(i,j,k,n) += flux_fab(i,j,k,n);
              });

            }


            Array4<Real> mass_fluxes_fab = (*mass_fluxes[idir]).array(mfi);

            AMREX_HOST_DEVICE_FOR_4D(mfi.nodaltilebox(idir), 1, i, j, k, n,
            {
              mass_fluxes_fab(i,j,k,0) = flux_fab(i,j,k,URHO);
            });

          } // idir loop

        }

    }

}

