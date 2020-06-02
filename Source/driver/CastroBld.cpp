
#include "AMReX_LevelBld.H"
#include "Castro.H"

using namespace amrex;

class CastroBld : public LevelBld {
    virtual void variableSetUp() override;
    virtual void variableCleanUp() override;
    virtual AmrLevel* operator()() override;
    virtual AmrLevel* operator()(Amr& papa, int lev, const Geometry& level_geom, const BoxArray& ba,
                                 const DistributionMapping& dm, Real time) override;
};

CastroBld Castro_bld;

LevelBld* getLevelBld() { return &Castro_bld; }

void CastroBld::variableSetUp() { Castro::variableSetUp(); }

void CastroBld::variableCleanUp() { Castro::variableCleanUp(); }

AmrLevel* CastroBld::operator()() { return new Castro; }

AmrLevel* CastroBld::operator()(Amr& papa, int lev, const Geometry& level_geom, const BoxArray& ba,
                                const DistributionMapping& dm, Real time) {
    return new Castro(papa, lev, level_geom, ba, dm, time);
}
