#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

#ifndef CASTRO_MGUTILS_H
#define CASTRO_MGUTILS_H

using namespace amrex;

void apply_metric(const Box& bx, Array4<Real> const rhs, const Box& rbx, Array4<Real> const ecx,
                  const Box& xbx,

#if AMREX_SPACEDIM >= 2
                  Array4<Real> const ecy, const Box& ybx,
#endif
                  const Real* dx, const int coord_type);

void do_weight_cc(const Box& bx, Array4<Real> const cc, const Real* dx, const int coord_type);

void do_unweight_cc(const Box& bx, Array4<Real> const cc, const Real* dx, const int coord_type);

void do_unweight_edges(const Box& bx, Array4<Real> const ec, const int idir, const Real* dx,
                       const int coord_type);

#endif
