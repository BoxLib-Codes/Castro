
#include <climits>

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_INTERP_F.H>
#include <AMReX_Interpolater.H>

#include "DG_Interpolater.H"

namespace amrex {

DG_Interpolater dg_interp;

DG_Interpolater::DG_Interpolater ()
{}

DG_Interpolater::~DG_Interpolater ()
{}

Box
DG_Interpolater::CoarseBox (const Box&     fine,
                                   const IntVect& ratio)
{
    Box crse = amrex::coarsen(fine,ratio);
    // we don't need slopes
    //    crse.grow(1);
    return crse;
}

Box
DG_Interpolater::CoarseBox (const Box& fine,
                                   int        ratio)
{
    Box crse(amrex::coarsen(fine,ratio));
    // we don't need slopes
    //     crse.grow(1);
    return crse;
}

void
DG_Interpolater::interp (const FArrayBox&     crse,
                         int                  crse_comp,
                         FArrayBox&           fine,
                         int                  fine_comp,
                         int                  ncomp,
                         const Box&           fine_region,
                         const IntVect&       ratio,
                         const Geometry&      crse_geom,
                         const Geometry&      fine_geom,
                         Vector<BCRec> const& bcr,
                         int                  actual_comp,
                         int                  actual_state,
                         RunOn                gpu_or_cpu)
{
    BL_PROFILE("DG_Interpolater::interp()");
    BL_ASSERT(bcr.size() >= ncomp);

    //
    // Make box which is intersection of fine_region and domain of fine.
    //
    Box target_fine_region = fine_region & fine.box();

    const int* fblo   = target_fine_region.loVect();
    const int* fbhi   = target_fine_region.hiVect();

    const int* ratioV = ratio.getVect();

    ca_dg_refine(AMREX_ARLIM_ANYD(fblo), AMREX_ARLIM_ANYD(fbhi),
                 BL_TO_FORTRAN_ANYD(fine),
                 BL_TO_FORTRAN_ANYD(crse),
                 AMREX_ARLIM_ANYD(ratioV), &ncomp);
}

}
