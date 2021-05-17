
#include <AMReX_BLFort.H>
#include <Castro.H>
#include <Castro_bc_fill_nd.H>
#include <Castro_generic_fill.H>
#include <problem_bc_fill.H>

using namespace amrex;

void ca_statefill(Box const& bx, FArrayBox& data,
                  const int dcomp, const int numcomp,
                  Geometry const& geom, const Real time,
                  const Vector<BCRec>& bcr, const int bcomp,
                  const int scomp)
{
    // Here dcomp is the component in the destination array that we
    // are filling and bcr is a vector of length ncomp which are the
    // BC values corresponding to components dcomp to dcomp + ncomp -
    // 1

    // First, fill all the BC data using the default routines.
    // We replace inflow with outflow in the generic fill to ensure that
    // valid data is always present.

    Vector<BCRec> bcr_noinflow{bcr};
    for (int i = 0; i < bcr_noinflow.size(); ++i) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            if (bcr_noinflow[i].lo(dir) == EXT_DIR) {
                bcr_noinflow[i].setLo(dir, FOEXTRAP);
            }
            if (bcr_noinflow[i].hi(dir) == EXT_DIR) {
                bcr_noinflow[i].setHi(dir, FOEXTRAP);
            }
        }
    }

    GpuBndryFuncFab<CastroGenericFill> gpu_bndry_func(CastroGenericFill{});
    gpu_bndry_func(bx, data, dcomp, numcomp, geom, time, bcr_noinflow, bcomp, scomp);

    // At this point, if we filling anything other than the full state (numcomp == NUM_STATE),
    // we assume that we are doing a derive or some other routine which is not actually setting
    // boundary conditions on the state data. So we immediately return. This means that the derive
    // data may not be exactly what the user wants at the physical boundary, but the impact of that
    // is usually negligible.

    if (numcomp != NUM_STATE) {
        return;
    }

    // Fill ambient BCs.

    ambient_fill(bx, data.array(dcomp), geom, bcr);

    // Now we consider external BCs (HSE).  Note, if we are at a
    // corner where two (or three) faces want to do HSE, we may run
    // into a situation that the data is not valid in the corner where
    // we start the integration.  We'll abort, for now, if we run into
    // this case.
    //
    // The future fix is to first call ext_fill on the ghost cells
    // that are not corners and then call it a second time on just the
    // corners.

#if AMREX_SPACEDIM == 2
    if ((bcr[URHO].lo(0) == EXT_DIR && bcr[URHO].lo(1) == EXT_DIR) ||
        (bcr[URHO].lo(0) == EXT_DIR && bcr[URHO].hi(1) == EXT_DIR) ||
        (bcr[URHO].hi(0) == EXT_DIR && bcr[URHO].lo(1) == EXT_DIR) ||
        (bcr[URHO].hi(0) == EXT_DIR && bcr[URHO].hi(1) == EXT_DIR)) {
        amrex::Error("Error: external boundaries meeting at a corner not supported");
    }
#endif

#if AMREX_SPACEDIM == 3
    if ((bcr[URHO].lo(0) == EXT_DIR &&           // xl, yl, zl corner
         (bcr[URHO].lo(1) == EXT_DIR || bcr[URHO].lo(2) == EXT_DIR)) ||
        (bcr[URHO].lo(1) == EXT_DIR && bcr[URHO].lo(2) == EXT_DIR) ||
        (bcr[URHO].lo(0) == EXT_DIR &&           // xl, yr, zl corner
         (bcr[URHO].hi(1) == EXT_DIR || bcr[URHO].lo(2) == EXT_DIR)) ||
        (bcr[URHO].hi(1) == EXT_DIR && bcr[URHO].lo(2) == EXT_DIR) ||
        (bcr[URHO].lo(0) == EXT_DIR &&           // xl, yl, zr corner
         (bcr[URHO].lo(1) == EXT_DIR || bcr[URHO].hi(2) == EXT_DIR)) ||
        (bcr[URHO].lo(1) == EXT_DIR && bcr[URHO].hi(2) == EXT_DIR) ||
        (bcr[URHO].lo(0) == EXT_DIR &&           // xl, yr, zr corner
         (bcr[URHO].hi(1) == EXT_DIR || bcr[URHO].hi(2) == EXT_DIR)) ||
        (bcr[URHO].hi(1) == EXT_DIR && bcr[URHO].hi(2) == EXT_DIR) ||
        (bcr[URHO].hi(0) == EXT_DIR &&           // xr, yl, zl corner
         (bcr[URHO].lo(1) == EXT_DIR || bcr[URHO].lo(2) == EXT_DIR)) ||
        (bcr[URHO].lo(1) == EXT_DIR && bcr[URHO].lo(2) == EXT_DIR) ||
        (bcr[URHO].hi(0) == EXT_DIR &&           // xr, yr, zl corner
         (bcr[URHO].hi(1) == EXT_DIR || bcr[URHO].lo(2) == EXT_DIR)) ||
        (bcr[URHO].hi(1) == EXT_DIR && bcr[URHO].lo(2) == EXT_DIR) ||
        (bcr[URHO].hi(0) == EXT_DIR &&           // xr, yl, zr corner
         (bcr[URHO].lo(1) == EXT_DIR || bcr[URHO].hi(2) == EXT_DIR)) ||
        (bcr[URHO].lo(1) == EXT_DIR && bcr[URHO].hi(2) == EXT_DIR) ||
        (bcr[URHO].hi(0) == EXT_DIR &&           // xr, yr, zr corner
         (bcr[URHO].hi(1) == EXT_DIR || bcr[URHO].hi(2) == EXT_DIR)) ||
        (bcr[URHO].hi(1) == EXT_DIR && bcr[URHO].hi(2) == EXT_DIR)) {
        amrex::Error("Error: external boundaries meeting at a corner not supported");
    }
#endif

#ifdef GRAVITY
    hse_fill(bx, data.array(), geom, bcr, time);
#endif

    // Finally, override with problem-specific boundary conditions.

    const auto state = data.array();

    // Copy BCs to an Array1D so they can be passed by value to the ParallelFor.

    Array1D<BCRec, 0, NUM_STATE - 1> bcs;
    for (int n = 0; n < numcomp; ++n) {
        bcs(n) = bcr[n];
    }

    const auto geomdata = geom.data();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
        problem_bc_fill(i, j, k, state, time, bcs, geomdata);
    });
}

