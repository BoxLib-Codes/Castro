#ifndef CASTRO_RADIATION_H
#define CASTRO_RADIATION_H

#include <AMReX_AmrLevel.H>

#include <Castro.H>
#include <RadTypes.H>
#include <RadBndry.H>
#include <MGRadBndry.H>
#include <RadSolve.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_Array.H>

#include <radiation_params.H>

///
/// @class Radiation
/// @brief
///
class Radiation {

public:

  enum Solver_Type { InvalidSolver = -1,
                     SingleGroupSolver = 0,
                     // ProtoMultiGroupSolver = 1,  // deleted
                     //  NeutrinoSolver = 2,  // deleted
                     //  SGFLDHybrid = 3,  // deleted
                     //  SGFLDImplicit = 4, // deleted
                     SGFLDSolver = 5,
                     MGFLDSolver = 6};

  static Solver_Type SolverType;

  static amrex::Real radtoE, Etorad, radfluxtoF;

  static int do_multigroup;
  static int nGroups;

  static int rad_hydro_combined;
  static int comoving;
  static int Er_Lorentz_term;
  static int fspace_advection_type;
  static int do_inelastic_scattering;    ///< Compton scattering for photon only

  static int plot_lambda, plot_kappa_p, plot_kappa_r;
  static int plot_lab_Er, plot_lab_flux, plot_com_flux;
  static int icomp_lambda, icomp_kp, icomp_kr;
  static int icomp_lab_Er, icomp_lab_Fr, icomp_com_Fr;
  static int nplotvar;   // extra plot variables
  static amrex::Vector<std::string> plotvar_names;

  static int filter_lambda_T, filter_lambda_S;
  static int filter_prim_int, filter_prim_T, filter_prim_S;

  static int accelerate;        ///< controls multigroup convergence acceleration

  int verbose;

  int maxiter;           ///< iteration limit for implicit update loop
  int miniter;
  int convergence_check_type; ///< 0 (default): both 2 and 3,
                              ///< 1: rhoe,
                              ///< 2: residue of Eq. rhoe,
                              ///< 3: T
  amrex::Real relInTol, absInTol; ///< tolerance for inner iternation of J equation
  int maxInIter;           ///< iteration limit for inner iteration of J equation
  int minInIter;
  int skipAccelAllowed;   ///< Skip acceleration if it doesn't help
  int matter_update_type; ///< 0: conservative  1: non-conservative  2: C and NC interwoven
                          ///< The last outer iteration is always conservative.
  int n_bisect;  ///< Bisection after n_bisect iterations
  amrex::Real dedT_fac; ///< Make dedT larger for safety in Newton iteration
  int inner_convergence_check;
  amrex::Real delta_e_rat_dt_tol; ///< energy change tolerance for adjusting timestep
  amrex::Real delta_T_rat_dt_tol; ///< T change tolerance for adjusting timestep
  static int limiter;          ///< 0: no limiter, 2: Lev-Pom limiter
                               ///< 12: Bruenn, 22: square root, 32: Minerbo
  static int closure;      ///< 0: \f$f = \lambda\f$, 1: \f$f = 1/3\f$,
                           ///< 2: \f$f = 1-2\lambda\f$, 3: \f$f = \lambda+(\lambda R)^2\f$
                           ///< 4: \f$f = 1/3 + 2/3 (\lambda R)^2\f$
  int update_planck;     ///< after this number of iterations, lag planck
  int update_rosseland;  ///< after this number of iterations, lag rosseland
  int update_opacity;
  int update_limiter;    ///< after this number of iterations, lag limiter
  int inner_update_limiter; ///< This is for MGFLD solver.
                            ///< Stop updating limiter after ? inner iterations
                            ///< 0 means lagging by one outer iteration
  amrex::Real dT;               ///< temperature step for derivative estimate
  amrex::Real underfac;         ///< factor controlling progressive underrelaxation
  int do_sync;               ///< perform sync (if false zero out sync source)
  int do_kappa_stm_emission;

  static void read_static_params();

///
/// Physical constants for use in multigroup, made static and public
/// so that other modules such as RadSolve can get to them easily.
/// Values are set in the Radiation constructor.
///
  static amrex::Real convert_MeV_erg, clight, hPlanck, kBoltz, Avogadro;

///
/// Values used internally for speed of light, Stefan-Boltzmann constant,
/// radiation constant.  (For contrived test problems, these may differ
/// from the proper physical values, so be careful if you use them!)
///
  static amrex::Real c, sigma, aRad;

  static int current_group_number;
  static std::string current_group_name;

  static int pure_hydro;

  amrex::Vector<std::unique_ptr<amrex::MultiFab> > plotvar;


///
/// @param Parent
/// @param castro
/// @param restart
///
  Radiation(amrex::Amr* Parent, class Castro* castro, int restart = 0);
  ~Radiation() { }


///
/// @param level
/// @param grids
/// @param dmap
///
  void regrid(int level, const amrex::BoxArray& grids,
              const amrex::DistributionMapping& dmap);


///
/// @param level
///
  void close(int level);


///
/// @param level
/// @param grids
/// @param dmap
/// @param dir
/// @param is
///
  void restart(int level, const amrex::BoxArray& grids,
               const amrex::DistributionMapping& dmap,
               const std::string& dir, std::istream& is);


///
/// @param level
/// @param dir
/// @param os
/// @param how
///
  void checkPoint(int level,
                  const std::string& dir,
                  std::ostream&  os,
                  amrex::VisMF::How     how);


///
/// access to group information

///
/// @param igroup
///
  void set_current_group(int igroup);

/// @param i
///
  amrex::Real group_center(int i) {
    return nugroup[i];
  }


///
/// @param i
///
  amrex::Real group_weight(int i) {
    return dnugroup[i];
  }



///
/// diagnostics for adjusting timestep
/// @param lev
///
  amrex::Real deltaEnergyRatMax(int lev) {
    return delta_e_rat_level[lev];
  }


///
/// @param lev
///
  amrex::Real deltaTRatMax(int lev) {
    return delta_T_rat_level[lev];
  }

  amrex::Real deltaEnergyTol() {
    return delta_e_rat_dt_tol;
  }

  amrex::Real deltaTTol() {
    return delta_T_rat_dt_tol;
  }

///
/// access to conservation flux register for energy diagnostics
///
/// @param level
///
  amrex::FluxRegister* consRegister(int level) {
      return flux_cons_old[level].get();
  }



///
/// @param level
///
  void post_init(int level);


///
/// @param level
///
  void pre_timestep(int level);


///
/// @param level
/// @param ncycle
///
  void init_flux(int level, int ncycle);


///
/// @param level
///
  void deferred_sync_setup(int level);


///
/// @param level
/// @param rhs
/// @param indx
///
  void deferred_sync(int level, amrex::MultiFab& rhs, int indx);


///
/// @param level
/// @param iteration
/// @param ncycle
///
  void single_group_update(int level, int iteration, int ncycle);


///
/// @param level
/// @param iteration
/// @param ncycle
///
  void MGFLD_implicit_update(int level, int iteration, int ncycle);


///
/// @param level
///
  void inelastic_scattering(int level);


///
/// @param os
///
  void write_groups(ostream& os);


///
/// @param level
///
  void reflux(int level);


///
/// @param level
/// @param grids
/// @param &Sborder
/// @param &Erborder
/// @param &lamborder
///
  void compute_limiter(int level, const amrex::BoxArray& grids,
                       const amrex::MultiFab &Sborder,
                       const amrex::MultiFab &Erborder,
                       amrex::MultiFab &lamborder);


///
/// @param state
/// @param Er
/// @param gPr
/// @param dx
/// @param box
///
  void estimate_gamrPr(const amrex::FArrayBox& state, const amrex::FArrayBox& Er,
                       amrex::FArrayBox& gPr, const amrex::Real* dx, const amrex::Box& box);


///
/// @param estdt
/// @param level
///
  void EstTimeStep(amrex::Real& estdt, int level);


///
/// @param level
/// @param State
///
  void filter_prim(int level, amrex::MultiFab& State);

  // internal utility functions:


///
/// @param exch
/// @param Er
/// @param fkp
/// @param igroup
///
  void compute_exchange(amrex::MultiFab& exch, amrex::MultiFab& Er, amrex::MultiFab& fkp,
                        int igroup = -1);


///
/// @param eta
/// @param etainv
/// @param state
/// @param temp
/// @param fkp
/// @param Er
/// @param delta_t
/// @param c
/// @param underrel
/// @param lag_planck
/// @param igroup
///
  void compute_eta(amrex::MultiFab& eta, amrex::MultiFab& etainv, amrex::MultiFab& state,
                   amrex::MultiFab& temp, amrex::MultiFab& fkp, amrex::MultiFab& Er,
                   amrex::Real delta_t, amrex::Real c, amrex::Real underrel, int lag_planck,
                   int igroup = 0);


///
/// @param relative
/// @param absolute
/// @param frhoes
/// @param frhoem
/// @param eta
/// @param etainv
/// @param dflux_old
/// @param dflux_new
/// @param exchange
/// @param delta_t
///
  void internal_energy_update(amrex::Real& relative, amrex::Real& absolute,
                              amrex::MultiFab& frhoes, amrex::MultiFab& frhoem,
                              amrex::MultiFab& eta, amrex::MultiFab& etainv,
                              amrex::MultiFab& dflux_old, amrex::MultiFab& dflux_new,
                              amrex::MultiFab& exchange, amrex::Real delta_t);


///
/// @param relative
/// @param absolute
/// @param frhoes
/// @param frhoem
/// @param eta
/// @param etainv
/// @param dflux_old
/// @param dflux_new
/// @param exchange
/// @param Dterm
/// @param delta_t
///
  void internal_energy_update(amrex::Real& relative, amrex::Real& absolute,
                              amrex::MultiFab& frhoes, amrex::MultiFab& frhoem,
                              amrex::MultiFab& eta, amrex::MultiFab& etainv,
                              amrex::MultiFab& dflux_old, amrex::MultiFab& dflux_new,
                              amrex::MultiFab& exchange, amrex::MultiFab& Dterm,
                              amrex::Real delta_t);


///
/// @param relative
/// @param absolute
/// @param frhoes
/// @param frhoem
/// @param eta
/// @param etainv
/// @param Er_new
/// @param dflux_old
/// @param dflux_new
/// @param temp
/// @param fkp
/// @param state
/// @param delta_t
///
  void nonconservative_energy_update(amrex::Real& relative, amrex::Real& absolute,
                                     amrex::MultiFab& frhoes, amrex::MultiFab& frhoem,
                                     amrex::MultiFab& eta, amrex::MultiFab& etainv,
                                     amrex::MultiFab& Er_new,
                                     amrex::MultiFab& dflux_old, amrex::MultiFab& dflux_new,
                                     amrex::MultiFab& temp, amrex::MultiFab& fkp,
                                     amrex::MultiFab& state, amrex::Real delta_t);


///
/// @param state
/// @param frhoes
///
  void state_update(amrex::MultiFab& state,  amrex::MultiFab& frhoes);


///
/// @param bd
/// @param Er
/// @param time
/// @param level
///
  void getBndryData(RadBndry& bd, amrex::MultiFab& Er, amrex::Real time, int level);

///
/// multigroup version
///
/// @param mgbd
/// @param Er
/// @param time
/// @param level
///
  void getBndryDataMG(MGRadBndry& mgbd, amrex::MultiFab& Er, amrex::Real time, int level);

///
/// @param mgbd
/// @param Er
/// @param level
///
  void getBndryDataMG_ga(MGRadBndry& mgbd, amrex::MultiFab& Er, int level);


///
/// @param bdry
/// @param level
/// @param time
///
  void filBndry(amrex::BndryRegister& bdry, int level, amrex::Real time);

  // Flux limiter functions, potentially for use by all update modules

///
/// Computes the scaled gradient for use in flux limiters & Eddington factors
///
/// @param level
/// @param amrex::Array<amrex::MultiFab
/// @param R
/// @param kappa_r
/// @param kcomp
/// @param Er
/// @param igroup
/// @param limiter
/// @param nGrow_Er
/// @param Rcomp
///
  void scaledGradient(int level,
                      amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& R,
                      amrex::MultiFab& kappa_r, int kcomp,
                      amrex::MultiFab& Er, int igroup,
                      int limiter, int nGrow_Er=0, int Rcomp=0);

///
/// On input, lambda should contain scaled gradient.
/// On output this will be overwritten with the flux limiter.
///
/// @param level
/// @param amrex::Array<amrex::MultiFab
/// @param lambda
/// @param limiter
/// @param lamcomp
///
  void fluxLimiter(int level,
                   amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& lambda,
                   int limiter, int lamcomp=0);

///
/// Fab versions of conversion functions.
///

///
/// @param c_v
/// @param temp
/// @param state
/// @param reg
///
  void get_c_v(amrex::FArrayBox& c_v, amrex::FArrayBox& temp,
               amrex::FArrayBox& state, const amrex::Box& reg);


///
/// MultiFab versions of conversion functions call Fab versions for eos data
///
/// @param frhoe
/// @param state
///
  void get_frhoe(amrex::MultiFab& frhoe, amrex::MultiFab& state);

///
/// Planck mean is only needed in the grid interiors, but must be
/// computed separately for each fluid (fkp.nVar() == N fluids).
/// temp contains frhoe on input
///
/// @param fkp
/// @param temp
/// @param state
/// @param igroup
/// @param delta_t
///
  void get_planck_and_temp(amrex::MultiFab& fkp, amrex::MultiFab& temp, amrex::MultiFab& state,
                           int igroup = 0, amrex::Real delta_t = 0.0);

///
/// Rosseland mean is needed on an extended region (nGrow == 1).
/// Fills patches with state information on extended region
///
/// @param kappa_r
/// @param hclaw
/// @param igroup
///
  void get_rosseland(amrex::MultiFab& kappa_r, amrex::AmrLevel* hclaw, int igroup = 0);

///
/// For the interior opacities, temperature is assumed to be known (input)
///
/// Updates Rosseland mean in interior of level, leaves ghost values
/// at physical or coarse-fine interfaces untouched
///
/// @param kappa_r
/// @param temp
/// @param state
/// @param geom
/// @param igroup
///
  void update_rosseland_from_temp(amrex::MultiFab& kappa_r,
                                  amrex::MultiFab& temp, amrex::MultiFab& state,
                                  const amrex::Geometry& geom,
                                  int igroup = 0);

///
/// RadHydro
///
/// @param kappa_r
/// @param v
/// @param dcf
/// @param delta_t
/// @param c
/// @param hclaw
/// @param igroup
///
  void get_rosseland_v_dcf(amrex::MultiFab& kappa_r, amrex::MultiFab& v, amrex::MultiFab& dcf,
                           amrex::Real delta_t, amrex::Real c,
                           amrex::AmrLevel* hclaw,
                           int igroup = 0);

///
/// multigroup methods
///
/// @param verbose
///
  void get_groups(int verbose);


///
/// @param f
/// @param indx
///
  void extrapolateBorders(amrex::MultiFab& f, int indx);


///
/// @param kappa_r
/// @param state
///
  void SGFLD_compute_rosseland(amrex::MultiFab& kappa_r, const amrex::MultiFab& state);

///
/// @param kappa_r
/// @param state
///
  void SGFLD_compute_rosseland(amrex::FArrayBox& kappa_r, const amrex::FArrayBox& state);

///
/// @param dcf
/// @param etainv
/// @param kp
/// @param kr
/// @param geom
///
  void update_dcf(amrex::MultiFab& dcf, amrex::MultiFab& etainv, amrex::MultiFab& kp, amrex::MultiFab& kr,
                  const amrex::Geometry& geom);

///
/// <MGFLD>
///
  int use_dkdT;


///
/// <Shestakov-Bolstad>
///
  int use_WiensLaw;
  amrex::Real Tf_Wien;

///
/// </ Shestakov-Bolstad>
///

///
/// @param relative_in
/// @param absolute_in
/// @param error_er
/// @param Er_new
/// @param Er_pi
/// @param kappa_p
/// @param etaTz
/// @param temp_new
/// @param delta_t
///
  void check_convergence_er(amrex::Real& relative_in, amrex::Real& absolute_in, amrex::Real& error_er,
                            const amrex::MultiFab& Er_new, const amrex::MultiFab& Er_pi,
                            const amrex::MultiFab& kappa_p,
                            const amrex::MultiFab& etaTz,
                            const amrex::MultiFab& temp_new,
                            amrex::Real delta_t);

///
/// @param rhoe_new
/// @param rhoe_star
/// @param rhoe_step
/// @param Er_new
/// @param temp_new
/// @param temp_star
/// @param rho
/// @param kappa_p
/// @param jg
/// @param dedT
/// @param rel_rhoe
/// @param abs_rhoe
/// @param rel_FT
/// @param abs_FT
/// @param rel_T
/// @param abs_T
/// @param delta_t
///
  void check_convergence_matt(const amrex::MultiFab& rhoe_new, const amrex::MultiFab& rhoe_star,
                              const amrex::MultiFab& rhoe_step, const amrex::MultiFab& Er_new,
                              const amrex::MultiFab& temp_new, const amrex::MultiFab& temp_star,
                              const amrex::MultiFab& rho,
                              const amrex::MultiFab& kappa_p, const amrex::MultiFab& jg,
                              const amrex::MultiFab& dedT,
                              amrex::Real& rel_rhoe, amrex::Real& abs_rhoe,
                              amrex::Real& rel_FT,   amrex::Real& abs_FT,
                              amrex::Real& rel_T,    amrex::Real& abs_T,
                              amrex::Real delta_t);

///
/// @param coupT
/// @param kappa_p
/// @param Er_pi
/// @param jg
///
  void compute_coupling(amrex::MultiFab& coupT,
                        const amrex::MultiFab& kappa_p, const amrex::MultiFab& Er_pi,
                        const amrex::MultiFab& jg);

///
/// @param etaT
/// @param etaTz
/// @param eta1
/// @param djdT
/// @param dkdT
/// @param dedT
/// @param Er_star
/// @param rho
/// @param delta_t
/// @param ptc_tau
///
  void compute_etat(amrex::MultiFab& etaT, amrex::MultiFab& etaTz,
                    amrex::MultiFab& eta1, amrex::MultiFab& djdT,
                    const amrex::MultiFab& dkdT, const amrex::MultiFab& dedT,
                    const amrex::MultiFab& Er_star, const amrex::MultiFab& rho,
                    amrex::Real delta_t, amrex::Real ptc_tau);

///
/// @param S_new
/// @param temp_new
/// @param temp_star
/// @param kappa_p
/// @param kappa_r
/// @param jg
/// @param djdT
/// @param dkdT
/// @param dedT
/// @param level
/// @param it
/// @param ngrow
///
  void eos_opacity_emissivity(const amrex::MultiFab& S_new,
                              const amrex::MultiFab& temp_new,
                              const amrex::MultiFab& temp_star,
                              amrex::MultiFab& kappa_p, amrex::MultiFab& kappa_r, amrex::MultiFab& jg,
                              amrex::MultiFab& djdT, amrex::MultiFab& dkdT, amrex::MultiFab& dedT,
                              int level, int it, int ngrow);

///
/// @param Er_new
/// @param Er_pi
/// @param kappa_p
/// @param kappa_r
/// @param etaT
/// @param eta1
/// @param mugT
/// @param lambda
/// @param solver
/// @param mgbd
/// @param grids
/// @param level
/// @param time
/// @param delta_t
/// @param ptc_tau
///
  void gray_accel(amrex::MultiFab& Er_new, amrex::MultiFab& Er_pi,
                  amrex::MultiFab& kappa_p, amrex::MultiFab& kappa_r,
                  amrex::MultiFab& etaT, amrex::MultiFab& eta1,
                  amrex::MultiFab& mugT,
                  amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& lambda,
                  RadSolve* solver, MGRadBndry& mgbd,
                  const amrex::BoxArray& grids, int level, amrex::Real time, amrex::Real delta_t, amrex::Real ptc_tau);

///
/// @param Er_new
/// @param Er_pi
/// @param kappa_p
/// @param etaT
/// @param mugT
/// @param delta_t
/// @param ptc_tau
///
  void local_accel(amrex::MultiFab& Er_new, const amrex::MultiFab& Er_pi,
                   const amrex::MultiFab& kappa_p,
                   const amrex::MultiFab& etaT,
                   const amrex::MultiFab& mugT,
                   amrex::Real delta_t, amrex::Real ptc_tau);

///
/// @param state
/// @param rhoe
/// @param temp
/// @param grids
/// @param derat
/// @param dT
/// @param level
///
  void state_energy_update(amrex::MultiFab& state, const amrex::MultiFab& rhoe,
                           const amrex::MultiFab& temp, const amrex::BoxArray& grids,
                           amrex::Real& derat, amrex::Real& dT, int level);

///
/// @param rhoe_new
/// @param temp_new
/// @param Er_new
/// @param Er_pi
/// @param rhoe_star
/// @param rhoe_step
/// @param etaT
/// @param etaTz
/// @param eta1
/// @param coupT
/// @param kappa_p
/// @param jg
/// @param mugT
/// @param S_new
/// @param level
/// @param delta_t
/// @param ptc_tau
/// @param it
/// @param conservative_update
///
  void update_matter(amrex::MultiFab& rhoe_new, amrex::MultiFab& temp_new,
                     const amrex::MultiFab& Er_new, const amrex::MultiFab& Er_pi,
                     const amrex::MultiFab& rhoe_star,
                     const amrex::MultiFab& rhoe_step,
                     const amrex::MultiFab& etaT, const amrex::MultiFab& etaTz,
                     const amrex::MultiFab& eta1,
                     const amrex::MultiFab& coupT,
                     const amrex::MultiFab& kappa_p, const amrex::MultiFab& jg,
                     const amrex::MultiFab& mugT,
                     const amrex::MultiFab& S_new,
                     int level, amrex::Real delta_t,
                     amrex::Real ptc_tau, int it, bool conservative_update);

///
/// @param rhoe_new
/// @param temp_new
/// @param rhoe_star
/// @param temp_star
/// @param S_new
/// @param grids
/// @param level
///
  void bisect_matter(amrex::MultiFab& rhoe_new, amrex::MultiFab& temp_new,
                     const amrex::MultiFab& rhoe_star, const amrex::MultiFab& temp_star,
                     const amrex::MultiFab& S_new, const amrex::BoxArray& grids, int level);

///
/// for the hyperbolic solver
///
/// @param kappa_r
/// @param state
///
  void MGFLD_compute_rosseland(amrex::MultiFab& kappa_r, const amrex::MultiFab& state);

///
/// @param kappa_r
/// @param state
///
  void MGFLD_compute_rosseland(amrex::FArrayBox& kappa_r, const amrex::FArrayBox& state);


///
/// @param kappa_s
/// @param state
///
  void MGFLD_compute_scattering(amrex::FArrayBox& kappa_s, const amrex::FArrayBox& state);


///
/// @param rhs
/// @param dt
/// @param level
///
  void rhstoEr(amrex::MultiFab& rhs, amrex::Real dt, int level);

///
/// </ MGFLD>
///
/// @param level
/// @param amrex::Array<amrex::MultiFab
/// @param lambda
///
  void save_lambda_in_plotvar(int level, const amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>& lambda);

///
/// @param level
/// @param Snew
/// @param Ecom
/// @param F
/// @param iflx
///
  void save_lab_Er_in_plotvar(int level, const amrex::MultiFab& Snew,
                              const amrex::MultiFab& Ecom, const amrex::MultiFab& F, int iflx);

///
/// @param level
/// @param Snew
/// @param amrex::Array<amrex::MultiFab
/// @param lambda
/// @param Er
/// @param F
/// @param iflx
/// @param lab_factor
///
  void save_flux_in_plotvar(int level, const amrex::MultiFab& Snew,
                            const amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>& lambda,
                            const amrex::MultiFab& Er, const amrex::MultiFab& F,
                            int iflx, const amrex::Real lab_factor = 1.0);

///
/// @param level
/// @param Snew
/// @param amrex::Array<amrex::MultiFab
/// @param lambda
/// @param Er
/// @param F
/// @param iflx
///
  void save_com_flux_in_plotvar(int level, const amrex::MultiFab& Snew,
                                const amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>& lambda,
                                const amrex::MultiFab& Er, const amrex::MultiFab& F, int iflx);

  amrex::Vector<amrex::Real> xnu, nugroup, dnugroup;

protected:

  amrex::Amr* parent;
  amrex::BCRec rad_bc;          ///< types defined in LO_BCTYPES.H, not BCTYPES.H
  amrex::Real reltol, abstol;   ///< tolerances for implicit update loop
  amrex::Vector<amrex::Real> delta_e_rat_level, delta_T_rat_level;


///
/// arrays of flux registers
///   flux_cons[level] is a amrex::FluxRegister, which contains fluxes for all groups
///
  amrex::Vector <std::unique_ptr<amrex::FluxRegister> > flux_cons;
  amrex::Vector <std::unique_ptr<amrex::FluxRegister> > flux_cons_old;
  amrex::Vector <std::unique_ptr<amrex::FluxRegister> > flux_trial;


///
/// for deferred sync
///
  amrex::Vector<amrex::Real> delta_t_old;


///
/// divergence of flux
///
  amrex::Vector<std::unique_ptr<amrex::MultiFab> > dflux;

  std::string group_units;
  amrex::Real group_print_factor;

  amrex::IntVect spot;

};

#endif
