#ifndef DG_INTERPOLATER_H_
#define DG_INTERPOLATER_H_

#include <AMReX_Box.H>
#include <AMReX_BCRec.H>
#include <AMReX_REAL.H>
#include <AMReX_Interpolater.H>
#include "Castro.H"
#include "TwoMoment_F.H"

namespace amrex {

class Geometry;
class FArrayBox;
class InterpolaterBoxCoarsener;

class DG_Interpolater
    :
    public Interpolater
{
public:
    //
    // The constructor.
    //
    DG_Interpolater ();
    //
    // The destructor.
    //
    virtual ~DG_Interpolater () override;
    //
    // Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box& fine,
                           int        ratio) override;
    //
    // Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box&     fine,
                           const IntVect& ratio) override;
    //
    // Coarse to fine interpolation in space.
    //
    virtual void interp (const FArrayBox&     crse,
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
                         RunOn                gpu_or_cpu) override;

protected:
};

extern DG_Interpolater dg_interp;

}
#endif /*DG_INTERPOLATER_H_*/

