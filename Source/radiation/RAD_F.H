#ifndef CASTRO_RAD_F_H
#define CASTRO_RAD_F_H

#include <AMReX_ArrayLim.H>

#ifdef BL_FORT_USE_UNDERSCORE
#  define   FORT_RADBNDRY2 rbndr2_
#  define   FORT_INIT_OPACITY_TABLE init_opacity_table_

#elif defined(BL_FORT_USE_UPPERCASE)

#  define   FORT_RADBNDRY2 RBNDR2
#  define   FORT_INIT_OPACITY_TABLE INIT_OPACITY_TABLE

#endif

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
void ca_init_radhydro_pars
    (const int & fspace_advection_type, const int & do_inelastic_scattering,
     const int & comoving);

void ca_initradconstants
  (const amrex::Real& pi, amrex::Real& clight, amrex::Real& hplanck,
   amrex::Real& kboltz, amrex::Real& stefbol, amrex::Real& avogadro,
   amrex::Real& mev2erg);

void ca_initsinglegroup
  (const int& ngroups);

void ca_get_dlognu
  (amrex::Real* dlognu_out);

void ca_get_nugroup
  (amrex::Real* nugroup);

void ca_get_dnugroup
  (amrex::Real* dnugroup);

#ifdef __cplusplus
}
#endif


BL_FORT_PROC_DECL(CA_INITGROUPS,ca_initgroups)
  (const amrex::Real* nugroup, const amrex::Real* dnugroup,
   const int& ngroups, const int& ng0, const int& ng1);

BL_FORT_PROC_DECL(CA_INITGROUPS2,ca_initgroups2)
  (const amrex::Real* nugroup, const amrex::Real* dnugroup,
   const amrex::Real* xnu, const int& ngroups);

BL_FORT_PROC_DECL(CA_INITGROUPS3,ca_initgroups3)
  (const amrex::Real* nugroup, const amrex::Real* dnugroup, const amrex::Real* dlognugroup,
   const amrex::Real* xnu, const int& ngroups, const int& ng0, const int& ng1);

BL_FORT_PROC_DECL(CA_COMPUTE_KAPKAP, ca_compute_kapkap)
     (BL_FORT_FAB_ARG(kapkap), const BL_FORT_FAB_ARG(kap_r)); 

BL_FORT_PROC_DECL(CA_COMPUTE_LAMBORDER, ca_compute_lamborder)
   (const BL_FORT_FAB_ARG(Er), const BL_FORT_FAB_ARG(kap), 
    BL_FORT_FAB_ARG(lam),  const amrex::Real* dx, const int* ngrow, const int* limiter,
    const int* filter_lambda_T, const int* filter_lambda_S);
 
#ifdef __cplusplus
extern "C" {
#endif
  void ca_correct_dterm
    (D_DECL(BL_FORT_FAB_ARG(dfx),
            BL_FORT_FAB_ARG(dfy),
            BL_FORT_FAB_ARG(dfz)),
     const amrex::Real* re, const amrex::Real* rc);
#ifdef __cplusplus
}
#endif

// <MGFLD>
#ifdef __cplusplus
extern "C" {
#endif
void ca_initfluxlimiter
   (const int* limiter, const int* closure);
#ifdef __cplusplus
}
#endif

extern "C"
{
  void ca_rad_source
    (const int* lo, const int* hi, 
     BL_FORT_FAB_ARG_3D(rhs),
     const amrex::Real* dx, amrex::Real dt,
     amrex::Real time, int igroup);

}

#ifdef __cplusplus
extern "C" {
#endif
  void ca_filt_prim
    (const int lo[], const int hi[],
     BL_FORT_FAB_ARG(Stmp),
     BL_FORT_FAB_ARG(Snew),
     BL_FORT_FAB_ARG(mask),
     const int* filter_prim_T, const int* filter_prim_S,
     const int* domlo, const int* domhi,
     const amrex::Real* delta, const amrex::Real* xlo, const amrex::Real* problo,
     const amrex::Real* time, const int* level);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C"
{
#endif

  void ca_flux_face2center
   (const int lo[], const int hi[],
    BL_FORT_FAB_ARG(test), 
    const BL_FORT_FAB_ARG(flx),
    const amrex::Real* x, const int* xlo, const int* xhi, 
    const int* ntest, const int* idim, const int* iflx);

#ifdef __cplusplus
}
#endif
// </ MGFLD>

#ifdef __cplusplus
extern "C" {
#endif
  void sphc(amrex::Real* r, amrex::Real* s, ARLIM_P(reglo), ARLIM_P(reghi), const amrex::Real* dx);

  void sphe(amrex::Real* r, amrex::Real* s, const int& idim, 
            ARLIM_P(edgeboxlo), ARLIM_P(edgeboxhi), const amrex::Real* dx);

  void ceupdterm(ARLIM_P(reglo), ARLIM_P(reghi), 
                 amrex::Real& relres, amrex::Real& absres,
                 BL_FORT_FAB_ARG(frhoes), 
                 amrex::Real* frhoem, amrex::Real* eta, amrex::Real* etainv,
                 amrex::Real* dflux_old, amrex::Real* dflux_new,
                 amrex::Real* exch, amrex::Real* dterm, 
                 const amrex::Real& delta_t, const amrex::Real& theta);

  void nceup(ARLIM_P(reglo), ARLIM_P(reghi), 
             amrex::Real& relres, amrex::Real& absres,
             BL_FORT_FAB_ARG(frhoes), 
             amrex::Real* frhoem, amrex::Real* eta, amrex::Real* etainv,
             BL_FORT_FAB_ARG(er), 
             amrex::Real* dflux_old, amrex::Real* dflux_new,
             amrex::Real* temp, amrex::Real* fkp, amrex::Real* cv, 
             BL_FORT_FAB_ARG(state), 
             const amrex::Real& sigma, const amrex::Real& c,
             const amrex::Real& delta_t, const amrex::Real& theta);

  void fkpn(const int* lo, const int* hi,
            BL_FORT_FAB_ARG_3D(fkp),       
            amrex::Real con, amrex::Real em, amrex::Real en,
            amrex::Real ep, amrex::Real xnu, amrex::Real tf,
            BL_FORT_FAB_ARG_3D(temp), 
            BL_FORT_FAB_ARG_3D(state),
            amrex::Real temp_offset);

  void rosse1(const int* lo, const int* hi,
              amrex::Real con, amrex::Real em, amrex::Real en,
              amrex::Real ep, amrex::Real nu,
              amrex::Real tf, amrex::Real kfloor,
              BL_FORT_FAB_ARG_3D(state),
              BL_FORT_FAB_ARG_3D(kappar));

  void rosse1s(const int* lo, const int* hi,
               amrex::Real con, amrex::Real em, amrex::Real en,
               amrex::Real ep,  amrex::Real scon,
               amrex::Real sem, amrex::Real sen,
               amrex::Real sep, amrex::Real nu,
               amrex::Real tf, amrex::Real kfloor,
               BL_FORT_FAB_ARG_3D(state),
               BL_FORT_FAB_ARG_3D(kappar));
  
  void rfface(BL_FORT_FAB_ARG(fine), 
              BL_FORT_FAB_ARG(crse), 
              const int& idim, const int* irat);

  void FORT_RADBNDRY2(amrex::Real* bf, ARLIM_P(blo), ARLIM_P(bhi), 
                      int* tfab, ARLIM_P(dlo), ARLIM_P(dhi),
                      const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real& time);

  void bextrp(BL_FORT_FAB_ARG(f), 
              ARLIM_P(reglo), ARLIM_P(reghi));

  void FORT_INIT_OPACITY_TABLE(const int& iverb);

#ifdef __cplusplus
}
#endif

#endif
