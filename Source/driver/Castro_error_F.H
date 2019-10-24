#ifndef _Castro_error_F_H_
#define _Castro_error_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void ca_denerror
    (const int* lo, const int* hi,
     int8_t* tag, const int* tag_lo, const int* tag_hi,
     BL_FORT_FAB_ARG_3D(den), const int ncomp,
     const amrex::Real* dx, const amrex::Real* problo,
     const int8_t tagval, const int8_t clearval,
     const amrex::Real time, const int level);

  void ca_ye_error
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
      BL_FORT_FAB_ARG_3D(ye),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* problo,
     const amrex::Real* time, const int* level);
  
  void ca_enterror
    (const int* lo, const int* hi,
     int8_t* tag, const int* tag_lo, const int* tag_hi,
     BL_FORT_FAB_ARG_3D(ent), const int ncomp,
     const amrex::Real* dx, const amrex::Real* problo,
     const int8_t tagval, const int8_t clearval,
     const amrex::Real time, const int level);

  void ca_velerror
    (const int* lo, const int* hi,
     int8_t* tag, const int* tag_lo, const int* tag_hi,
     BL_FORT_FAB_ARG_3D(vel), const int ncomp,
     const amrex::Real* dx, const amrex::Real* problo,
     const int8_t tagval, const int8_t clearval,
     const amrex::Real time, const int level);

  void ca_temperror
    (const int* lo, const int* hi,
     int8_t* tag, const int* tag_lo, const int* tag_hi,
     BL_FORT_FAB_ARG_3D(temp), const int ncomp,
     const amrex::Real* dx, const amrex::Real* problo,
     const int8_t tagval, const int8_t clearval,
     const amrex::Real time, const int level);

  void ca_presserror
    (const int* lo, const int* hi,
     int8_t* tag, const int* tag_lo, const int* tag_hi,
     BL_FORT_FAB_ARG_3D(press), const int ncomp,
     const amrex::Real* dx, const amrex::Real* problo,
     const int8_t tagval, const int8_t clearval,
     const amrex::Real time, const int level);

  void ca_nucerror
    (const int* lo, const int* hi,
     int8_t* tag, const int* tag_lo, const int* tag_hi,
     BL_FORT_FAB_ARG_3D(state), const int ncomp,
     const amrex::Real* dx, const amrex::Real* problo,
     const int8_t tagval, const int8_t clearval,
     const amrex::Real time, const int level);

  void ca_enucerror
    (const int* lo, const int* hi,
     int8_t* tag, const int* tag_lo, const int* tag_hi,
     BL_FORT_FAB_ARG_3D(state), const int ncomp,
     const amrex::Real* dx, const amrex::Real* problo,
     const int8_t tagval, const int8_t clearval,
     const amrex::Real time, const int level);

#ifdef RADIATION
  void ca_raderror
    (const int* lo, const int* hi,
     int8_t* tag, const int* tag_lo, const int* tag_hi,
     BL_FORT_FAB_ARG_3D(rad), const int ncomp,
     const amrex::Real* dx, const amrex::Real* problo,
     const int8_t tagval, const int8_t clearval,
     const amrex::Real time, const int level);
#endif

  void set_problem_tags
    (const int* lo, const int* hi,
     int8_t* tag, const int* tag_lo, const int* tag_hi,
     const BL_FORT_FAB_ARG_3D(state),
     const amrex::Real* dx, const amrex::Real* problo,
     const int8_t tagval, const int8_t clearval,
     const amrex::Real time, const int level);

#ifdef __cplusplus
}
#endif

#endif
