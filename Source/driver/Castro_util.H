#ifndef CASTRO_UTIL_H
#define CASTRO_UTIL_H

#include "AMReX_Geometry.H"

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int upassmap(int ipassive) {
    if (ipassive < NumAdv) {
        return UFA + ipassive;
    } else if (ipassive < NumAdv + NumSpec) {
        return UFS + ipassive - NumAdv;
    } else {
        return UFX + ipassive - NumAdv - NumSpec;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int qpassmap(int ipassive) {
    if (ipassive < NumAdv) {
        return QFA + ipassive;
    } else if (ipassive < NumAdv + NumSpec) {
        return QFS + ipassive - NumAdv;
    } else {
        return QFX + ipassive - NumAdv - NumSpec;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool mom_flux_has_p(const int mom_dir, const int flux_dir,
                                                             const int coord) {

    // Keep track of which components of the momentum flux have pressure

#if AMREX_SPACEDIM <= 2

    if (mom_dir == 0 && flux_dir == 0) {

        if (AMREX_SPACEDIM == 1 || (AMREX_SPACEDIM == 2 && coord == CoordSys::RZ)) {
            return false;
        } else {
            return true;
        }

    } else {

        return (mom_dir == flux_dir);
    }

#else  // AMREX_SPACEDIM == 3; Cartesian

    return (mom_dir == flux_dir);

#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void position(int i, int j, int k,
                                                       GeometryData const& geomdata,
                                                       GpuArray<Real, 3>& loc, bool ccx = true,
                                                       bool ccy = true, bool ccz = true) {
    // Given 3D indices (i,j,k), return the cell-centered spatial position.
    // Optionally we can also be edge-centered in any of the directions.

    int idx[3] = {i, j, k};

    bool cc[3] = {ccx, ccy, ccz};

    Real offset[AMREX_SPACEDIM];

    const Real* problo = geomdata.ProbLo();
    const Real* probhi = geomdata.ProbHi();
    const Real* dx = geomdata.CellSize();
    const int* domlo = geomdata.domain.loVect();
    const int* domhi = geomdata.domain.hiVect();

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (cc[dir]) {
            // If we're cell-centered, we want to be in the middle of the zone.
            offset[dir] = problo[dir] + 0.5_rt * dx[dir];
        } else {
            // We are edge-centered.
            offset[dir] = problo[dir];
        }
    }

    // Be careful when using periodic boundary conditions. In that case,
    // we need to loop around to the other side of the domain.

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {

        if (geomdata.isPeriodic(dir)) {
            if (idx[dir] < domlo[dir]) {
                offset[dir] += probhi[dir] - problo[dir];
            }

            // For the hi end we need to deal correctly with
            // edge-centered indexing.
            if ((cc[dir] && idx[dir] > domhi[dir]) || (!cc[dir] && idx[dir] > domhi[dir] + 1)) {
                offset[dir] += problo[dir] - probhi[dir];
            }
        }
    }

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        loc[dir] = offset[dir] + idx[dir] * dx[dir];
    }

    for (int dir = AMREX_SPACEDIM; dir < 3; ++dir) {
        loc[dir] = 0.0_rt;
    }
}

#endif
