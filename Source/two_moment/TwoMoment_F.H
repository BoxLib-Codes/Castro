#ifndef _TwoMoment_F_H_
#define _TwoMoment_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
  void InitThornado(const int* nDimsX, const int* nDimsE, const int* swE,
                    const amrex::Real* eL, const amrex::Real* eR,
                    const amrex::Real* zoomE, const int* nSpecies);

  void ca_get_rad_ncomp(int* rad_ncomp);

  void ca_get_thornado_node_averages(int* n_node_avgs);

  void ca_init_thornado_data(
     const int& level, const amrex::Real& time,
     const int* lo, const int* hi,
     const int& n_rad_comp,
     BL_FORT_FAB_ARG_3D(thornado_data),
     BL_FORT_FAB_ARG_3D(fluid_data),
     const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* xhi);

  void InitThornado_Patch(const int* len, int* swX,
                          const amrex::Real* grid_lo, const amrex::Real* grid_hi);

  void FreeThornado_Patch();

  void compute_thornado_timestep(const amrex::Real* dx, 
                                       amrex::Real& dt); 

  void call_to_thornado(const int lo[], const int hi[], amrex::Real* dt_sub,
                        BL_FORT_FAB_ARG_3D( U_F  ), const int* ns,
                        BL_FORT_FAB_ARG_3D(dU_F  ), const int* nds,
                        BL_FORT_FAB_ARG_3D( U_R_o), const int* n_uro,
                        BL_FORT_FAB_ARG_3D( U_R_n), const int* n_urn,
                        int* n_moments, int* my_ngrow);

  void ca_der_avgs_per_E
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(J)        ,const int* nvar,
     const BL_FORT_FAB_ARG_3D(U_R),const int* ncomp);

  void ca_der_J
    (BL_FORT_FAB_ARG_3D(J),const int* nvar,
     const BL_FORT_FAB_ARG_3D(U_R),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

  void ca_der_Hx
    (BL_FORT_FAB_ARG_3D(H_x),const int* nvar,
     const BL_FORT_FAB_ARG_3D(U_R),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

  void ca_der_Hy
    (BL_FORT_FAB_ARG_3D(H_y),const int* nvar,
     const BL_FORT_FAB_ARG_3D(U_R),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

  void ca_der_Hz
    (BL_FORT_FAB_ARG_3D(H_z),const int* nvar,
     const BL_FORT_FAB_ARG_3D(U_R),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

  void ca_dg_refine(const int* fblo, const int* fbhi,
                    BL_FORT_FAB_ARG_ANYD(fine),
                    const BL_FORT_FAB_ARG_ANYD(crse),
                    const int* nFineX,
                    const int* nvar);

  void ca_dg_coarsen(const int* cblo, const int* cbhi,
                     const BL_FORT_FAB_ARG_ANYD(fine),
                     BL_FORT_FAB_ARG_ANYD(crse),
                     const int* nFineX,
                     const int* nvar);

#ifdef __cplusplus
}
#endif

#endif
