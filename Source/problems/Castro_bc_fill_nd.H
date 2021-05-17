#ifndef CASTRO_BC_FILL_ND_H
#define CASTRO_BC_FILL_ND_H

#include <AMReX_BLFort.H>

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

void ca_statefill(amrex::Box const& bx, amrex::FArrayBox& data,
                  const int dcomp, const int numcomp,
                  amrex::Geometry const& geom, const amrex::Real time,
                  const amrex::Vector<amrex::BCRec>& bcr, const int bcomp,
                  const int scomp);

void
hse_fill(const amrex::Box& bx, amrex::Array4<amrex::Real> const& adv,
         amrex::Geometry const& geom, const amrex::Vector<amrex::BCRec>& bcr,
         const amrex::Real time);

///
/// Fill the boundaries with the ambient state
///
/// @param bx    box
/// @param state state
/// @param geom  geometry
/// @param bcr   boundary conditions
///
void
ambient_fill(const amrex::Box& bx, amrex::Array4<amrex::Real> const& state,
             amrex::Geometry const& geom, const amrex::Vector<amrex::BCRec>& bcr);

///
/// Fill the boundaries with the ambient state (density only)
///
/// @param bx    box
/// @param state state
/// @param geom  geometry
/// @param bcr   boundary conditions
///
void
ambient_denfill(const amrex::Box& bx, amrex::Array4<amrex::Real> const& state,
                amrex::Geometry const& geom, const amrex::Vector<amrex::BCRec>& bcr);

#endif
