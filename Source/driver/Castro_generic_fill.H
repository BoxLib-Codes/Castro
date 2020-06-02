#ifndef CASTRO_GENERIC_FILL_H
#define CASTRO_GENERIC_FILL_H

#include <AMReX_BLFort.H>
#include <Castro.H>

#ifdef AMREX_USE_CUDA
#include <AMReX_Arena.H>
#include <cuda_runtime_api.h>
#endif

#ifdef AMREX_USE_CUDA
static void set_bc_launch_config() {
    amrex::Gpu::Device::setNumThreadsMin(Castro::minBCThreads(0), Castro::minBCThreads(1),
                                         Castro::minBCThreads(2));
}

static void clean_bc_launch_config() { amrex::Gpu::Device::setNumThreadsMin(1, 1, 1); }

// Return a pointer to bc valid for use in Fortran. For the CPU this is a no-op.

static int* prepare_bc(const int* bc, const int nvar) {
    int* bc_f = (int*)amrex::The_Arena()->alloc(AMREX_SPACEDIM * 2 * nvar * sizeof(int));
    AMREX_GPU_SAFE_CALL(cudaMemcpyAsync(bc_f, bc, AMREX_SPACEDIM * 2 * nvar * sizeof(int),
                                        cudaMemcpyHostToDevice, amrex::Gpu::Device::cudaStream()));
    return bc_f;
}

static void clean_bc(int* bc_f) {
    amrex::Gpu::Device::streamSynchronize();
    amrex::The_Arena()->free(bc_f);
}
#endif

struct CastroGenericFill {
    AMREX_GPU_DEVICE
    void operator()(const IntVect& iv, Array4<Real> const& data, const int dcomp, const int numcomp,
                    GeometryData const& geom, const Real time, const BCRec* bcr, const int bcomp,
                    const int orig_comp) const {
        // Dummy routine that does nothing for inflow boundaries.
        // We assume that there are no inflow boundaries for a
        // generic fill, since we overwrote them with first-order
        // extrapolation in variableSetUp().
    }
};

namespace {
    static CastroGenericFill castro_generic_fill_func;
}

void ca_generic_fill(amrex::Box const& bx, amrex::FArrayBox& data, const int dcomp,
                     const int numcomp, amrex::Geometry const& geom, const amrex::Real time,
                     const amrex::Vector<amrex::BCRec>& bcr, const int bcomp, const int scomp);

#endif
