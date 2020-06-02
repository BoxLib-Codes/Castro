#ifndef CASTRO_DERIVE_H
#define CASTRO_DERIVE_H

#include <AMReX_BLFort.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

#ifdef __cplusplus
extern "C" {
#endif

void ca_derpres(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
                amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_dereint1(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                 const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
                 amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_dereint2(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                 const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
                 amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derlogden(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                  const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
                  amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_deruplusc(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                  const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
                  amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_deruminusc(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                   const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
                   amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_dersoundspeed(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                      const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
                      amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_dergamma1(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                  const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
                  amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_dermachnumber(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                      const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
                      amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derentropy(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                   const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
                   amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

#ifdef DIFFUSION
void ca_dercond(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
                amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derdiffcoeff(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
                     amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derdiffterm(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                    const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
                    amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);
#endif

void ca_derenuctimescale(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                         const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
                         amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_dervel(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
               const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
               amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_dermagvel(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                  const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
                  amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_dermaggrav(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                   const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
                   amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derradialvel(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
                     amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_dercircvel(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                   const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
                   amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_dermagmom(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                  const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
                  amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derangmomx(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                   const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                   amrex::Real time, const int* bcrec, int level);

void ca_derangmomy(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                   const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                   amrex::Real time, const int* bcrec, int level);

void ca_derangmomz(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                   const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                   amrex::Real time, const int* bcrec, int level);

void ca_derkineng(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                  const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                  amrex::Real time, const int* bcrec, int level);

void ca_dernull(BL_FORT_FAB_ARG_3D(der), const int* nvar, const BL_FORT_FAB_ARG_3D(data),
                const int* ncomp, const int* lo, const int* hi, const int* domain_lo,
                const int* domain_hi, const amrex::Real* delta, const amrex::Real* xlo,
                const amrex::Real* time, const amrex::Real* dt, const int* bcrec, const int* level,
                const int* grid_no);

void ca_derspec(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
                amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derabar(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
                amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_dermagvort(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                   const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
                   amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derdivu(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
                amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derstate(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                 const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
                 amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

#ifdef MHD
void ca_dermagcenx(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                   const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
                   amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_dermagceny(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                   const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
                   amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_dermagcenz(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                   const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
                   amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derex(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
              const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata, amrex::Real /*time*/,
              const int* /*bcrec*/, int /*level*/);

void ca_derey(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
              const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata, amrex::Real /*time*/,
              const int* /*bcrec*/, int /*level*/);

void ca_derez(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
              const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata, amrex::Real /*time*/,
              const int* /*bcrec*/, int /*level*/);

void ca_derdivb(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
                amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

#endif

#ifdef __cplusplus
}
#endif

/* problem-specific includes */
#include <Problem_Derive.H>

#endif
