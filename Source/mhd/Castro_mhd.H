
    void just_the_mhd(amrex::Real time, amrex::Real dt);

    void
    check_for_mhd_cfl_violation(const amrex::Box& bx, const amrex::Real dt,
                                amrex::Array4<amrex::Real const> const& q_arr,
                                amrex::Array4<amrex::Real const> const& qaux_arr);

    void
    plm(const amrex::Box& bx,
        const int idir,
        amrex::Array4<amrex::Real const> const& s,
        amrex::Array4<amrex::Real const> const& qaux,
        amrex::Array4<amrex::Real const> const& flatn,
        amrex::Array4<amrex::Real const> const& Bx,
        amrex::Array4<amrex::Real const> const& By,
        amrex::Array4<amrex::Real const> const& Bz,
        amrex::Array4<amrex::Real> const& qleft,
        amrex::Array4<amrex::Real> const& qright,
        amrex::Array4<amrex::Real const> const& srcQ,
        const amrex::Real dt);

    void
    ppm_mhd(const amrex::Box& bx,
            const int idir,
            amrex::Array4<amrex::Real const> const& q_arr,
            amrex::Array4<amrex::Real const> const& qaux,
            amrex::Array4<amrex::Real const> const& flatn,
            amrex::Array4<amrex::Real const> const& Bx,
            amrex::Array4<amrex::Real const> const& By,
            amrex::Array4<amrex::Real const> const& Bz,
            amrex::Array4<amrex::Real> const& qleft,
            amrex::Array4<amrex::Real> const& qright,
            amrex::Array4<amrex::Real const> const& srcQ,
            const amrex::Real dt);

    void
    consup_mhd(const amrex::Box& bx,
               amrex::Array4<amrex::Real> const& update,
               amrex::Array4<amrex::Real const> const& flux0,
               amrex::Array4<amrex::Real const> const& flux1,
               amrex::Array4<amrex::Real const> const& flux2);

    void
    PrimToCons(const amrex::Box& bx,
               amrex::Array4<amrex::Real const> const& q_arr,
               amrex::Array4<amrex::Real> const& u_arr);

    void
    prim_half(const amrex::Box& bx,
              amrex::Array4<amrex::Real> const& q2D,
              amrex::Array4<amrex::Real const> const& q_arr,
              amrex::Array4<amrex::Real const> const& flxx,
              amrex::Array4<amrex::Real const> const& flxy,
              amrex::Array4<amrex::Real const> const& flxz,
              const amrex::Real dt);


    void
    electric_edge_x(const amrex::Box& bx,
                    amrex::Array4<amrex::Real const> const& q_arr,
                    amrex::Array4<amrex::Real> const& E,
                    amrex::Array4<amrex::Real const> const& flxy,
                    amrex::Array4<amrex::Real const> const& flxz);

    void
    electric_edge_y(const amrex::Box& bx,
                    amrex::Array4<amrex::Real const> const& q_arr,
                    amrex::Array4<amrex::Real> const& E,
                    amrex::Array4<amrex::Real const> const& flxx,
                    amrex::Array4<amrex::Real const> const& flxz);

    void
    electric_edge_z(const amrex::Box& bx,
                    amrex::Array4<amrex::Real const> const& q_arr,
                    amrex::Array4<amrex::Real> const& E,
                    amrex::Array4<amrex::Real const> const& flxx,
                    amrex::Array4<amrex::Real const> const& flxy);

    void
    corner_couple(const amrex::Box& bx,
                  amrex::Array4<amrex::Real> const& qr_out,
                  amrex::Array4<amrex::Real> const& ql_out,
                  amrex::Array4<amrex::Real const> const& ur,
                  amrex::Array4<amrex::Real const> const& ul,
                  amrex::Array4<amrex::Real const> const& flxd2,
                  amrex::Array4<amrex::Real const> const& Ed1,
                  amrex::Array4<amrex::Real const> const& Ed3,
                  amrex::Array4<amrex::Real const> const& Ex,
                  amrex::Array4<amrex::Real const> const& Ey,
                  amrex::Array4<amrex::Real const> const& Ez,
                  const int d1, const int d2, const int d3,
                  const amrex::Real dt);

    void
    half_step(const amrex::Box& bx,
              amrex::Array4<amrex::Real> const& qr_out,
              amrex::Array4<amrex::Real> const& ql_out,
              amrex::Array4<amrex::Real const> const& ur,
              amrex::Array4<amrex::Real const> const& ul,
              amrex::Array4<amrex::Real const> const& flxd1,
              amrex::Array4<amrex::Real const> const& flxd2,
              amrex::Array4<amrex::Real const> const& Ed,
              amrex::Array4<amrex::Real const> const& Ed1,
              amrex::Array4<amrex::Real const> const& Ed2,
              amrex::Array4<amrex::Real const> const& Ex,
              amrex::Array4<amrex::Real const> const& Ey,
              amrex::Array4<amrex::Real const> const& Ez,
              const int d, const int d1, const int d2, const amrex::Real dt);

    void
    get_inplane_Bs_transverse_flux(const int face, const int comp, const int tdir,
                                   const int i, const int j, const int k,
                                   amrex::Array4<amrex::Real const> const& Ex,
                                   amrex::Array4<amrex::Real const> const& Ey,
                                   amrex::Array4<amrex::Real const> const& Ez,
                                   amrex::Real& Fl, amrex::Real& Fr);

    void
    hlld(const amrex::Box& bx,
         amrex::Array4<amrex::Real const> const& qleft,
         amrex::Array4<amrex::Real const> const& qright,
         amrex::Array4<amrex::Real> const& flx,
         const int dir);


