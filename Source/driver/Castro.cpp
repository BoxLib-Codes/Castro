#ifndef WIN32
#include <unistd.h>
#endif

#include <iomanip>

#include <algorithm>
#include <cstdio>
#include <ctime>
#include <iostream>
#include <string>
#include <vector>

#include <AMReX_CONSTANTS.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_TagBox.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>
#include <Castro.H>
#include <Castro_F.H>
#include <Castro_error_F.H>

#ifdef RADIATION
#include "RAD_F.H"
#include "Radiation.H"
#endif

#ifdef AMREX_PARTICLES
#include <AMReX_Particles.H>
#endif

#ifdef GRAVITY
#include "Gravity.H"
#endif

#ifdef DIFFUSION
#include "Diffusion.H"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef AMREX_USE_CUDA
#include <cuda_profiler_api.h>
#endif

#include <extern_parameters.H>

using namespace amrex;

bool Castro::signalStopJob = false;

std::vector<std::string> Castro::err_list_names;
std::vector<int> Castro::err_list_ng;
int Castro::num_err_list_default = 0;
int Castro::radius_grow = 1;
BCRec Castro::phys_bc;
int Castro::NUM_GROW = -1;

int Castro::lastDtPlotLimited = 0;
Real Castro::lastDtBeforePlotLimiting = 0.0;

Real Castro::num_zones_advanced = 0.0;

Vector<std::string> Castro::source_names;

Vector<int> Castro::upass_map;
Vector<int> Castro::qpass_map;

#ifdef TRUE_SDC
int Castro::SDC_NODES;
Vector<Real> Castro::dt_sdc;
Vector<Real> Castro::node_weights;
#endif

#ifdef AMREX_USE_CUDA
int Castro::numBCThreadsMin[3] = {1, 1, 1};
#endif

// the sponge parameters are controlled by Fortran, so
// this just initializes them before we grab their values
// from Fortran
#include <sponge_defaults.H>

#ifdef GRAVITY
// the gravity object
Gravity* Castro::gravity = 0;
#endif

#ifdef DIFFUSION
// the diffusion object
Diffusion* Castro::diffusion = 0;
#endif

#ifdef RADIATION
int Castro::do_radiation = -1;

// the radiation object
Radiation* Castro::radiation = 0;
#endif

std::string Castro::probin_file = "probin";

#if BL_SPACEDIM == 1
#ifndef AMREX_USE_CUDA
IntVect Castro::hydro_tile_size(1024);
#else
IntVect Castro::hydro_tile_size(1048576);
#endif
IntVect Castro::no_tile_size(1024);
#elif BL_SPACEDIM == 2
#ifndef AMREX_USE_CUDA
IntVect Castro::hydro_tile_size(1024, 16);
#else
IntVect Castro::hydro_tile_size(1048576, 1048576);
#endif
IntVect Castro::no_tile_size(1024, 1024);
#else
#ifndef AMREX_USE_CUDA
IntVect Castro::hydro_tile_size(1024, 16, 16);
#else
IntVect Castro::hydro_tile_size(1048576, 1048576, 1048576);
#endif
IntVect Castro::no_tile_size(1024, 1024, 1024);
#endif

// this will be reset upon restart
Real Castro::previousCPUTimeUsed = 0.0;

Real Castro::startCPUTime = 0.0;

int Castro::SDC_Source_Type = -1;
int Castro::num_state_type = 0;

namespace amrex {
    extern int compute_new_dt_on_regrid;
}

// Castro::variableSetUp is in Castro_setup.cpp
// variableCleanUp is called once at the end of a simulation
void Castro::variableCleanUp() {
#ifdef GRAVITY
    if (gravity != 0) {
        if (verbose > 1 && ParallelDescriptor::IOProcessor()) {
            std::cout << "Deleting gravity in variableCleanUp..." << '\n';
        }
        delete gravity;
        gravity = 0;
    }
#endif

#ifdef DIFFUSION
    if (diffusion != 0) {
        if (verbose > 1 && ParallelDescriptor::IOProcessor()) {
            std::cout << "Deleting diffusion in variableCleanUp..." << '\n';
        }
        delete diffusion;
        diffusion = 0;
    }
#endif

#ifdef RADIATION
    if (radiation != 0) {
        int report = (verbose || radiation->verbose);
        if (report && ParallelDescriptor::IOProcessor()) {
            std::cout << "Deleting radiation in variableCleanUp..." << '\n';
        }
        delete radiation;
        radiation = 0;
        if (report && ParallelDescriptor::IOProcessor()) {
            std::cout << "                                        done" << std::endl;
        }
    }
#endif

#ifdef AMREX_PARTICLES
    delete TracerPC;
    TracerPC = 0;
#endif

    desc_lst.clear();

    ca_finalize_meth_params();

    ca_network_finalize();

    eos_finalize();
    ca_eos_finalize();

#ifdef SPONGE
    sponge_finalize();
#endif
    amrinfo_finalize();
}

void Castro::read_params() {
    static bool done = false;

    if (done) return;

    done = true;

    ParmParse pp("castro");

    using namespace castro;

#include <castro_queries.H>

    // Get boundary conditions
    Vector<int> lo_bc(BL_SPACEDIM), hi_bc(BL_SPACEDIM);
    pp.getarr("lo_bc", lo_bc, 0, BL_SPACEDIM);
    pp.getarr("hi_bc", hi_bc, 0, BL_SPACEDIM);
    for (int i = 0; i < BL_SPACEDIM; i++) {
        phys_bc.setLo(i, lo_bc[i]);
        phys_bc.setHi(i, hi_bc[i]);
    }

    const Geometry& dgeom = DefaultGeometry();

    //
    // Check phys_bc against possible periodic geometry
    // if periodic, must have internal BC marked.
    //
    if (dgeom.isAnyPeriodic()) {
        //
        // Do idiot check.  Periodic means interior in those directions.
        //
        for (int dir = 0; dir < BL_SPACEDIM; dir++) {
            if (dgeom.isPeriodic(dir)) {
                if (lo_bc[dir] != Interior) {
                    std::cerr << "Castro::read_params:periodic in direction " << dir
                              << " but low BC is not Interior\n";
                    amrex::Error();
                }
                if (hi_bc[dir] != Interior) {
                    std::cerr << "Castro::read_params:periodic in direction " << dir
                              << " but high BC is not Interior\n";
                    amrex::Error();
                }
            }
        }
    } else {
        //
        // Do idiot check.  If not periodic, should be no interior.
        //
        for (int dir = 0; dir < BL_SPACEDIM; dir++) {
            if (lo_bc[dir] == Interior) {
                std::cerr << "Castro::read_params:interior bc in direction " << dir
                          << " but not periodic\n";
                amrex::Error();
            }
            if (hi_bc[dir] == Interior) {
                std::cerr << "Castro::read_params:interior bc in direction " << dir
                          << " but not periodic\n";
                amrex::Error();
            }
        }
    }

    if (dgeom.IsRZ() && (lo_bc[0] != Symmetry)) {
        std::cerr
            << "ERROR:Castro::read_params: must set r=0 boundary condition to Symmetry for r-z\n";
        amrex::Error();
    }

#if (BL_SPACEDIM == 1)
    if (dgeom.IsSPHERICAL()) {
        if ((lo_bc[0] != Symmetry) && (dgeom.ProbLo(0) == 0.0)) {
            std::cerr << "ERROR:Castro::read_params: must set r=0 boundary condition to Symmetry "
                         "for spherical\n";
            amrex::Error();
        }
    }
#elif (BL_SPACEDIM == 2)
    if (dgeom.IsSPHERICAL()) {
        amrex::Abort("We don't support spherical coordinate systems in 2D");
    }
#elif (BL_SPACEDIM == 3)
    if (dgeom.IsRZ()) {
        amrex::Abort("We don't support cylindrical coordinate systems in 3D");
    } else if (dgeom.IsSPHERICAL()) {
        amrex::Abort("We don't support spherical coordinate systems in 3D");
    }
#endif

#ifdef HYBRID_MOMENTUM
    // We do not support hybrid advection when using the HLLC solver.

    if (riemann_solver == 2) {
        amrex::Abort("HLLC Riemann solver unsupported when using hybrid momentum.");
    }
#endif

#ifdef REACTIONS
#ifdef SIMPLIFIED_SDC
    if (jacobian == 1) {
        amrex::Abort("Simplified SDC requires the numerical Jacobian now (jacobian = 2)");
    }
#endif
#endif
    // sanity checks

    if (grown_factor < 1) {
        amrex::Error("grown_factor must be integer >= 1");
    }

    if (cfl <= 0.0 || cfl > 1.0) {
        amrex::Error("Invalid CFL factor; must be between zero and one.");
    }

    // SDC does not support CUDA yet
#ifdef AMREX_USE_CUDA
    if (time_integration_method == SpectralDeferredCorrections) {
        amrex::Error("CUDA SDC is currently disabled.");
    }
#endif

    // Simplified SDC currently requires USE_SIMPLIFIED_SDC to be defined.
    // Also, if we have USE_SIMPLIFIED_SDC defined, we can't use the other
    // time integration_methods, because only the SDC burner
    // interface is available in Microphysics in this case.
#ifndef SIMPLIFIED_SDC
    if (time_integration_method == SimplifiedSpectralDeferredCorrections) {
        amrex::Error("Simplified SDC currently requires USE_SIMPLIFIED_SDC=TRUE when compiling.");
    }
#else
    if (time_integration_method != SimplifiedSpectralDeferredCorrections) {
        amrex::Error(
            "When building with USE_SIMPLIFIED_SDC=TRUE, only simplified SDC can be used.");
    }
#endif

#ifndef TRUE_SDC
    if (time_integration_method == SpectralDeferredCorrections) {
        amrex::Error("True SDC currently requires USE_TRUE_SDC=TRUE when compiling.");
    }
#else
    if (time_integration_method != SpectralDeferredCorrections) {
        amrex::Error("When building with USE_TRUE_SDC=TRUE, only true SDC can be used.");
    }
#endif

    if (hybrid_riemann == 1 && BL_SPACEDIM == 1) {
        std::cerr << "hybrid_riemann only implemented in 2- and 3-d\n";
        amrex::Error();
    }

    if (hybrid_riemann == 1 && (dgeom.IsSPHERICAL() || dgeom.IsRZ())) {
        std::cerr << "hybrid_riemann should only be used for Cartesian coordinates\n";
        amrex::Error();
    }

    // Make sure not to call refluxing if we're not actually doing any hydro.
    if (do_hydro == 0) {
        do_reflux = 0;
    }

    if (max_dt < fixed_dt) {
        std::cerr << "cannot have max_dt < fixed_dt\n";
        amrex::Error();
    }

#ifdef AMREX_PARTICLES
    read_particle_params();
#endif

#ifdef RADIATION
    pp.get("do_radiation", do_radiation);

    // Some radiation parameters are initialized here because they
    // may be used in variableSetUp, well before the call to the
    // Radiation constructor,

    if (do_radiation) {
        Radiation::read_static_params();
    }

    // radiation is only supported with CTU
    if (do_radiation && time_integration_method != CornerTransportUpwind) {
        amrex::Error("Radiation is currently only supported for CTU time advancement.");
    }
#endif

#ifdef ROTATION
    if (do_rotation) {
        if (rotational_period <= 0.0) {
            std::cerr << "Error:Castro::Rotation enabled but rotation period less than zero\n";
            amrex::Error();
        }
    }
    if (dgeom.IsRZ()) rot_axis = 2;
#if (BL_SPACEDIM == 1)
    if (do_rotation) {
        std::cerr << "ERROR:Castro::Rotation not implemented in 1d\n";
        amrex::Error();
    }
#endif
#endif

    // SCF initial model construction can only be done if both
    // rotation and gravity have been compiled in.

#if (!defined(GRAVITY) || !defined(ROTATION))
    if (do_scf_initial_model) {
        amrex::Error(
            "SCF initial model construction is only permitted if USE_GRAV=TRUE and "
            "USE_ROTATION=TRUE at compile time.");
    }
#endif

#ifdef AMREX_USE_CUDA
    if (do_scf_initial_model) {
        amrex::Error(
            "SCF initial model construction is currently not permitted if USE_CUDA=TRUE at compile "
            "time.");
    }
#endif

    StateDescriptor::setBndryFuncThreadSafety(bndry_func_thread_safe);

    ParmParse ppa("amr");
    ppa.query("probin_file", probin_file);

    Vector<int> tilesize(BL_SPACEDIM);
    if (pp.queryarr("hydro_tile_size", tilesize, 0, BL_SPACEDIM)) {
        for (int i = 0; i < BL_SPACEDIM; i++) {
            hydro_tile_size[i] = tilesize[i];
        }
    }

    // Override Amr defaults. Note: this function is called after Amr::Initialize()
    // in Amr::InitAmr(), right before the ParmParse checks, so if the user opts to
    // override our overriding, they can do so.

    compute_new_dt_on_regrid = 1;
}

Castro::Castro() : prev_state(num_state_type) {}

Castro::Castro(Amr& papa, int lev, const Geometry& level_geom, const BoxArray& bl,
               const DistributionMapping& dm, Real time)
    : AmrLevel(papa, lev, level_geom, bl, dm, time), prev_state(num_state_type) {
    MultiFab::RegionTag amrlevel_tag("AmrLevel_Level_" + std::to_string(lev));

    buildMetrics();

    initMFs();

    for (int i = 0; i < n_lost; i++) {
        material_lost_through_boundary_cumulative[i] = 0.0;
        material_lost_through_boundary_temp[i] = 0.0;
    }

    // Coterminous AMR boundaries are not supported in Castro if we're doing refluxing.

    if (do_hydro && do_reflux) {
        for (int ilev = 0; ilev <= parent->maxLevel(); ++ilev) {
            if (parent->nErrorBuf(ilev) == 0) {
                amrex::Error("n_error_buf = 0 is unsupported when using hydro.");
            }
        }
    }

#ifdef AMREX_USE_CUDA
    // Enforce our requirement on the blocking factor for CUDA. See Castro::variableSetUp() for details.
    for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
        for (int ilev = 0; ilev <= parent->maxLevel(); ++ilev) {
            if (parent->blockingFactor(ilev)[dim] % numBCThreadsMin[dim] != 0) {
                amrex::Error("Using CUDA requires a blocking factor that is a multiple of 8.");
            }
        }
    }
#endif

    // initialize all the new time level data to zero
    for (int k = 0; k < num_state_type; k++) {
        MultiFab& data = get_new_data(k);
        data.setVal(0.0, data.nGrow());
    }

#ifdef GRAVITY

    if (do_grav) {
        // gravity is a static object, only alloc if not already there
        if (gravity == 0) {
            gravity = new Gravity(parent, parent->finestLevel(), &phys_bc, URHO);
        }

        // Passing numpts_1d at level 0
        if (!level_geom.isAllPeriodic() && gravity != 0) {
            int numpts_1d = get_numpts();

            // For 1D, we need to add ghost cells to the numpts
            // given to us by Castro.

#if (BL_SPACEDIM == 1)
            numpts_1d += 2 * NUM_GROW;
#endif

            gravity->set_numpts_in_gravity(numpts_1d);
        }

        gravity->install_level(level, this, volume, area.data());

        if (verbose && level == 0 && ParallelDescriptor::IOProcessor()) {
            std::cout << "Setting the gravity type to " << gravity->get_gravity_type() << std::endl;
        }

#ifdef GRAVITY
        if (gravity->get_gravity_type() == "PoissonGrav" && gravity->NoComposite() != 0 &&
            gravity->NoSync() == 0) {
            std::cerr << "Error: not meaningful to have gravity.no_sync == 0 without having "
                         "gravity.no_composite == 0.";
            amrex::Error();
        }
#endif
    }

#endif

#ifdef DIFFUSION
    // diffusion is a static object, only alloc if not already there
    if (diffusion == 0) {
        diffusion = new Diffusion(parent, &phys_bc);
    }

    diffusion->install_level(level, this, volume, area.data());
#endif

#ifdef RADIATION
    if (do_radiation) {
        if (radiation == 0) {
            // radiation is a static object, only alloc if not already there
            radiation = new Radiation(parent, this);
        }
        radiation->regrid(level, grids, dmap);

        rad_solver.reset(new RadSolve(parent, level, grids, dmap));
    }
#endif
}

Castro::~Castro() {
#ifdef RADIATION
    if (radiation != 0) {
        //radiation->cleanup(level);
        radiation->close(level);
    }
#endif
}

void Castro::buildMetrics() {
    const int ngrd = grids.size();

    radius.resize(ngrd);

    const Real* dx = geom.CellSize();

    for (int i = 0; i < ngrd; i++) {
        const Box& b = grids[i];
        int ilo = b.smallEnd(0) - radius_grow;
        int ihi = b.bigEnd(0) + radius_grow;
        int len = ihi - ilo + 1;

        radius[i].resize(len);

        Real* rad = radius[i].dataPtr();

        if (Geom().IsCartesian()) {
            for (int j = 0; j < len; j++) {
                rad[j] = 1.0;
            }
        } else {
            RealBox gridloc = RealBox(grids[i], geom.CellSize(), geom.ProbLo());

            const Real xlo = gridloc.lo(0) + (0.5 - radius_grow) * dx[0];

            for (int j = 0; j < len; j++) {
                rad[j] = xlo + j * dx[0];
            }
        }
    }

    volume.clear();
    volume.define(grids, dmap, 1, NUM_GROW);
    geom.GetVolume(volume);

    for (int dir = 0; dir < BL_SPACEDIM; dir++) {
        area[dir].clear();
        area[dir].define(getEdgeBoxArray(dir), dmap, 1, NUM_GROW);
        geom.GetFaceArea(area[dir], dir);
    }
    for (int dir = BL_SPACEDIM; dir < 3; dir++) {
        area[dir].clear();
        area[dir].define(grids, dmap, 1, 0);
        area[dir].setVal(0.0);
    }

    dLogArea[0].clear();
#if (BL_SPACEDIM <= 2)
    geom.GetDLogA(dLogArea[0], grids, dmap, 0, NUM_GROW);
#endif

    if (level == 0) setGridInfo();

    wall_time_start = 0.0;
}

// Initialize the MultiFabs and flux registers that live as class members.

void Castro::initMFs() {
    fluxes.resize(3);

    for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
        fluxes[dir].reset(new MultiFab(getEdgeBoxArray(dir), dmap, NUM_STATE, 0));
    }

    for (int dir = BL_SPACEDIM; dir < 3; ++dir) {
        fluxes[dir].reset(new MultiFab(get_new_data(State_Type).boxArray(), dmap, NUM_STATE, 0));
    }

    mass_fluxes.resize(3);

    for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
        mass_fluxes[dir].reset(new MultiFab(getEdgeBoxArray(dir), dmap, 1, 0));
    }

    for (int dir = BL_SPACEDIM; dir < 3; ++dir) {
        mass_fluxes[dir].reset(new MultiFab(get_new_data(State_Type).boxArray(), dmap, 1, 0));
    }

#if (BL_SPACEDIM <= 2)
    if (!Geom().IsCartesian()) {
        P_radial.define(getEdgeBoxArray(0), dmap, 1, 0);
    }
#endif

#ifdef RADIATION
    if (Radiation::rad_hydro_combined) {
        rad_fluxes.resize(BL_SPACEDIM);
        for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
            rad_fluxes[dir].reset(new MultiFab(getEdgeBoxArray(dir), dmap, Radiation::nGroups, 0));
        }
    }
#endif

    if (do_reflux && level > 0) {

        flux_reg.define(grids, dmap, crse_ratio, level, NUM_STATE);
        flux_reg.setVal(0.0);

#if (BL_SPACEDIM < 3)
        if (!Geom().IsCartesian()) {
            pres_reg.define(grids, dmap, crse_ratio, level, 1);
            pres_reg.setVal(0.0);
        }
#endif

#ifdef RADIATION
        if (Radiation::rad_hydro_combined) {
            rad_flux_reg.define(grids, dmap, crse_ratio, level, Radiation::nGroups);
            rad_flux_reg.setVal(0.0);
        }
#endif

#ifdef GRAVITY
        if (do_grav && gravity->get_gravity_type() == "PoissonGrav" && gravity->NoSync() == 0) {
            phi_reg.define(grids, dmap, crse_ratio, level, 1);
            phi_reg.setVal(0.0);
        }
#endif
    }

    // Set the flux register scalings.

    if (do_reflux) {

        flux_crse_scale = -1.0;
        flux_fine_scale = 1.0;

        // The fine pressure scaling depends on dimensionality,
        // as the dimensionality determines the number of
        // adjacent zones. In 1D the face is a point so
        // there's only one fine neighbor for a given coarse
        // face; in 2D there's crse_ratio[1] faces adjacent
        // to a face perpendicular to the radial dimension;
        // and in 3D there would be crse_ratio**2, though
        // we do not separate the pressure out in 3D. Note
        // that the scaling by dt has already been handled
        // in the construction of the P_radial array.

        // The coarse pressure scaling is the same as for the
        // fluxes, we want the total refluxing contribution
        // over the full set of fine timesteps to equal P_radial.

#if (BL_SPACEDIM == 1)
        pres_crse_scale = -1.0;
        pres_fine_scale = 1.0;
#elif (BL_SPACEDIM == 2)
        pres_crse_scale = -1.0;
        pres_fine_scale = 1.0 / crse_ratio[1];
#endif
    }

    post_step_regrid = 0;

    lastDtRetryLimited = false;
    lastDtFromRetry = 1.e200;

    lastDt = 1.e200;
}

void Castro::setTimeLevel(Real time, Real dt_old, Real dt_new) {
    AmrLevel::setTimeLevel(time, dt_old, dt_new);
}

void Castro::setGridInfo() {

    // Send refinement data to Fortran. We do it here
    // because now the grids have been initialized and
    // we need this data for setting up the problem.
    // Note that this routine will always get called
    // on level 0, even if we are doing a restart,
    // so it is safe to put this here.

    if (level == 0) {

        const int max_level = parent->maxLevel();
        const int nlevs = max_level + 1;

        Vector<Real> dx_level(3 * nlevs);
        Vector<int> domlo_level(3 * nlevs);
        Vector<int> domhi_level(3 * nlevs);
        Vector<int> ref_ratio_to_f(3 * nlevs);
        Vector<int> n_error_buf_to_f(nlevs);
        Vector<int> blocking_factor_to_f(nlevs);

        const Real* dx_coarse = geom.CellSize();

        const int* domlo_coarse = geom.Domain().loVect();
        const int* domhi_coarse = geom.Domain().hiVect();

        for (int dir = 0; dir < 3; dir++) {
            if (dir < BL_SPACEDIM) {
                dx_level[dir] = dx_coarse[dir];

                domlo_level[dir] = domlo_coarse[dir];
                domhi_level[dir] = domhi_coarse[dir];
            } else {
                dx_level[dir] = 0.0;

                domlo_level[dir] = 0;
                domhi_level[dir] = 0;
            }

            // Refinement ratio and error buffer on finest level are meaningless,
            // and we want them to be zero on the finest level because some
            // of the algorithms depend on this feature.

            ref_ratio_to_f[dir + 3 * (nlevs - 1)] = 0;
            n_error_buf_to_f[nlevs - 1] = 0;
        }

        for (int lev = 0; lev <= max_level; lev++) {
            blocking_factor_to_f[lev] = parent->blockingFactor(lev)[0];
        }

        for (int lev = 1; lev <= max_level; lev++) {
            IntVect ref_ratio = parent->refRatio(lev - 1);

            // Note that we are explicitly calculating here what the
            // data would be on refined levels rather than getting the
            // data directly from those levels, because some potential
            // refined levels may not exist at the beginning of the simulation.

            for (int dir = 0; dir < 3; dir++) {
                if (dir < BL_SPACEDIM) {
                    dx_level[3 * lev + dir] = dx_level[3 * (lev - 1) + dir] / ref_ratio[dir];
                    int ncell =
                        (domhi_level[3 * (lev - 1) + dir] - domlo_level[3 * (lev - 1) + dir] + 1) *
                        ref_ratio[dir];
                    domlo_level[3 * lev + dir] = domlo_level[dir];
                    domhi_level[3 * lev + dir] = domlo_level[3 * lev + dir] + ncell - 1;
                    ref_ratio_to_f[3 * (lev - 1) + dir] = ref_ratio[dir];
                } else {
                    dx_level[3 * lev + dir] = 0.0;
                    domlo_level[3 * lev + dir] = 0;
                    domhi_level[3 * lev + dir] = 0;
                    ref_ratio_to_f[3 * (lev - 1) + dir] = 0;
                }
            }

            n_error_buf_to_f[lev - 1] = parent->nErrorBuf(lev - 1);
        }

        ca_set_grid_info(max_level, dx_level.data(), domlo_level.data(), domhi_level.data(),
                         ref_ratio_to_f.data(), n_error_buf_to_f.data(),
                         blocking_factor_to_f.data());
    }
}

void Castro::initData() {
    BL_PROFILE("Castro::initData()");

    //
    // Loop over grids, call FORTRAN function to init with data.
    //
    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    MultiFab& S_new = get_new_data(State_Type);
    Real cur_time = state[State_Type].curTime();

    S_new.setVal(0.);

    // make sure dx = dy = dz -- that's all we guarantee to support
#if (BL_SPACEDIM == 2)
    const Real SMALL = 1.e-13;
    if (fabs(dx[0] - dx[1]) > SMALL * dx[0]) {
        amrex::Abort("We don't support dx != dy");
    }
#elif (BL_SPACEDIM == 3)
    const Real SMALL = 1.e-13;
    if ((fabs(dx[0] - dx[1]) > SMALL * dx[0]) || (fabs(dx[0] - dx[2]) > SMALL * dx[0])) {
        amrex::Abort("We don't support dx != dy != dz");
    }
#endif

    ca_set_amr_info(level, -1, -1, -1.0, -1.0);

    if (verbose && ParallelDescriptor::IOProcessor()) {
        std::cout << "Initializing the data at level " << level << std::endl;
    }

#ifdef MHD
    MultiFab& Bx_new = get_new_data(Mag_Type_x);
    Bx_new.setVal(0.0);

    MultiFab& By_new = get_new_data(Mag_Type_y);
    By_new.setVal(0.0);

    MultiFab& Bz_new = get_new_data(Mag_Type_z);
    Bz_new.setVal(0.0);

#endif

    // Don't profile for this code, since there will be a lot of host
    // activity and GPU page faults that we're uninterested in.
#ifdef AMREX_USE_CUDA
    AMREX_GPU_SAFE_CALL(cudaProfilerStop());
#endif

#ifdef RADIATION
    // rad quantities are in the state even if (do_radiation == 0)
    MultiFab& Rad_new = get_new_data(Rad_Type);
    Rad_new.setVal(0.);
#endif

#ifdef REACTIONS
    MultiFab& React_new = get_new_data(Reactions_Type);
    React_new.setVal(0.);
#endif

#ifdef SIMPLIFIED_SDC
#ifdef REACTIONS
    if (time_integration_method == SimplifiedSpectralDeferredCorrections) {
        MultiFab& react_src_new = get_new_data(Simplified_SDC_React_Type);
        react_src_new.setVal(0.0, NUM_GROW);
    }
#endif
#endif

#ifdef MAESTRO_INIT
    MAESTRO_init();
#else
    {

#ifdef MHD
        int nbx = Bx_new.nComp();
        int nby = By_new.nComp();
        int nbz = Bz_new.nComp();

        Bx_new.setVal(0.0);
        By_new.setVal(0.0);
        Bz_new.setVal(0.0);

        for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
            RealBox gridloc(grids[mfi.index()], geom.CellSize(), geom.ProbLo());
            const Box& box = mfi.validbox();
            const int* lo = box.loVect();
            const int* hi = box.hiVect();

            BL_FORT_PROC_CALL(CA_INITMAG, ca_initmag)
            (level, cur_time, lo, hi, nbx, BL_TO_FORTRAN_3D(Bx_new[mfi]), nby,
             BL_TO_FORTRAN_3D(By_new[mfi]), nbz, BL_TO_FORTRAN_3D(Bz_new[mfi]), dx, gridloc.lo(),
             gridloc.hi());
        }

#endif  //MHD

#ifdef AMREX_USE_CUDA
        for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
#ifdef GPU_COMPATIBLE_PROBLEM
            // Prefetch data to the device to avoid page faults while we're initializing.
            S_new.prefetchToDevice(mfi);
#else
            // Prefetch data to the host (and then back to the device at the end)
            // to avoid expensive page faults while the initialization is done.
            S_new.prefetchToHost(mfi);
#endif
        }
#endif

        for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
            const Box& box = mfi.validbox();
            const int* lo = box.loVect();
            const int* hi = box.hiVect();

#ifdef GPU_COMPATIBLE_PROBLEM

#pragma gpu box(box)
            ca_initdata(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi), BL_TO_FORTRAN_ANYD(S_new[mfi]),
                        AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(prob_lo));

#else
            RealBox gridloc = RealBox(grids[mfi.index()], geom.CellSize(), geom.ProbLo());

            BL_FORT_PROC_CALL(CA_INITDATA, ca_initdata)
            (level, cur_time, ARLIM_3D(lo), ARLIM_3D(hi), NUM_STATE, BL_TO_FORTRAN_ANYD(S_new[mfi]),
             ZFILL(dx), ZFILL(gridloc.lo()), ZFILL(gridloc.hi()));

#endif
        }

        // it is not a requirement that the problem setup defines the
        // temperature, so we do that here _and_ ensure that we are
        // within any small limits
        computeTemp(
#ifdef MHD
            Bx_new, By_new, Bz_new,
#endif
            S_new, cur_time, 0);

        ReduceOps<ReduceOpSum, ReduceOpSum> reduce_op;
        ReduceData<int, int> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();

            auto S_arr = S_new.array(mfi);

            Real lsmall_temp = small_temp;
            Real lsmall_dens = small_dens;

            reduce_op.eval(bx, reduce_data,
                           [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept->ReduceTuple {
                               // if the problem tried to initialize a thermodynamic
                               // state that is at or below small_temp, then we abort.
                               // This is dangerous and we should recommend a smaller
                               // small_temp
                               int T_failed = 0;
                               if (S_arr(i, j, k, UTEMP) < lsmall_temp * 1.001) {
                                   T_failed = 1;
                               }

                               int rho_failed = 0;
                               if (S_arr(i, j, k, URHO) < lsmall_dens * 1.001) {
                                   rho_failed = 1;
                               }

                               return {T_failed, rho_failed};
                           });
        }

        ReduceTuple hv = reduce_data.value();
        int init_failed_T = amrex::get<0>(hv);
        int init_failed_rho = amrex::get<1>(hv);

        if (init_failed_rho != 0) {
            amrex::Error("Error: initial data has rho <~ small_dens");
        }

        if (init_failed_T != 0) {
            amrex::Error("Error: initial data has T <~ small_temp");
        }

#ifdef AMREX_USE_CUDA
#ifndef GPU_COMPATIBLE_PROBLEM
        for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
            S_new.prefetchToDevice(mfi);
        }
#endif
#endif

#ifdef HYBRID_MOMENTUM
        // Generate the initial hybrid momenta based on this user data.

        linear_to_hybrid_momentum(S_new, 0);
#endif

        // Verify that the sum of (rho X)_i = rho at every cell

        for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();

            auto S_arr = S_new.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept {
                Real spec_sum = 0.0_rt;
                for (int n = 0; n < NumSpec; n++) {
                    spec_sum += S_arr(i, j, k, UFS + n);
                }
                if (std::abs(S_arr(i, j, k, URHO) - spec_sum) > 1.e-8_rt * S_arr(i, j, k, URHO)) {
#ifndef AMREX_USE_CUDA
                    std::cout << "Sum of (rho X)_i vs rho at (i,j,k): " << i << " " << j << " " << k
                              << " " << spec_sum << " " << S_arr(i, j, k, URHO) << std::endl;
#endif
                    amrex::Error("Error: failed check of initial species summing to 1");
                }
            });
        }

#ifdef TRUE_SDC
        if (initialization_is_cell_average == 0) {
            // we are assuming that the initialization was done to cell-centers

            // Enforce that the total and internal energies are consistent.
            enforce_consistent_e(
#ifdef MHD
                Bx_new, By_new, Bz_new,
#endif
                S_new);

            // For fourth-order, we need to convert to cell-averages now.
            // (to second-order, these are cell-averages, so we're done in that case).

#ifndef AMREX_USE_CUDA
            if (sdc_order == 4) {
                Sborder.define(grids, dmap, NUM_STATE, NUM_GROW);
                AmrLevel::FillPatch(*this, Sborder, NUM_GROW, cur_time, State_Type, 0, NUM_STATE);

                // note: this cannot be tiled
                auto domain_lo = geom.Domain().loVect3d();
                auto domain_hi = geom.Domain().hiVect3d();

                FArrayBox tmp;

                for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
                    const Box& box = mfi.validbox();

                    tmp.resize(box, 1);
                    Elixir elix_tmp = tmp.elixir();
                    auto tmp_arr = tmp.array();

                    make_fourth_in_place(box, Sborder.array(mfi), tmp_arr, domain_lo, domain_hi);
                }

                // now copy back the averages
                MultiFab::Copy(S_new, Sborder, 0, 0, NUM_STATE, 0);
                Sborder.clear();
            }
#endif
        } else {

            Sborder.define(grids, dmap, NUM_STATE, NUM_GROW);
            AmrLevel::FillPatch(*this, Sborder, NUM_GROW, cur_time, State_Type, 0, NUM_STATE);

            // convert to centers -- not tile safe
            auto domain_lo = geom.Domain().loVect3d();
            auto domain_hi = geom.Domain().hiVect3d();

            FArrayBox tmp;

            for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
                const Box& box = mfi.growntilebox(2);

                tmp.resize(box, 1);
                Elixir elix_tmp = tmp.elixir();
                auto tmp_arr = tmp.array();

                make_cell_center_in_place(box, Sborder.array(mfi), tmp_arr, domain_lo, domain_hi);
            }

            // reset the energy -- do this in one ghost cell so we can average in place below
            for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
                const Box& box = mfi.growntilebox(1);

                auto S_arr = Sborder.array(mfi);

                amrex::ParallelFor(box, [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept {
                    Real rhoInv = 1.0_rt / S_arr(i, j, k, URHO);
                    Real u = S_arr(i, j, k, UMX) * rhoInv;
                    Real v = S_arr(i, j, k, UMY) * rhoInv;
                    Real w = S_arr(i, j, k, UMZ) * rhoInv;

                    eos_t eos_state;
                    eos_state.rho = S_arr(i, j, k, URHO);
                    eos_state.T = S_arr(i, j, k, UTEMP);
                    eos_state.e = S_arr(i, j, k, UEINT) * rhoInv - 0.5_rt * (u * u + v * v + w * w);
                    for (int n = 0; n < NumSpec; n++) {
                        eos_state.xn[n] = S_arr(i, j, k, UFS + n) * rhoInv;
                    }
                    for (int n = 0; n < NumAux; n++) {
                        eos_state.aux[n] = S_arr(i, j, k, UFX + n) * rhoInv;
                    }

                    eos(eos_input_re, eos_state);

                    S_arr(i, j, k, UTEMP) = eos_state.T;

                    S_arr(i, j, k, UEINT) = eos_state.rho * eos_state.e;
                });
            }

            // convert back to averages -- not tile safe
            for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
                const Box& box = mfi.validbox();

                tmp.resize(box, 1);
                Elixir elix_tmp = tmp.elixir();
                auto tmp_arr = tmp.array();

                make_fourth_in_place(box, Sborder.array(mfi), tmp_arr, domain_lo, domain_hi);
            }

            // now copy back the averages for UEINT and UTEMP only
            MultiFab::Copy(S_new, Sborder, UEINT, UEINT, 1, 0);
            MultiFab::Copy(S_new, Sborder, UTEMP, UTEMP, 1, 0);
            Sborder.clear();
        }
#else
        // Enforce that the total and internal energies are consistent.
        enforce_consistent_e(
#ifdef MHD
            Bx_new, By_new, Bz_new,
#endif
            S_new);
#endif

        // Do a FillPatch so that we can get the ghost zones filled.

        int ng = S_new.nGrow();

        if (ng > 0) {
            AmrLevel::FillPatch(*this, S_new, ng, cur_time, State_Type, 0, S_new.nComp());
        }
    }

    clean_state(
#ifdef MHD
        Bx_new, By_new, Bz_new,
#endif
        S_new, cur_time, S_new.nGrow());

#ifdef RADIATION
    if (do_radiation) {
        for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
            int i = mfi.index();

            if (radiation->verbose > 2) {
                std::cout << "Calling RADINIT at level " << level << ", grid " << i << std::endl;
            }

            const Box& box = mfi.validbox();
            const int* lo = box.loVect();
            const int* hi = box.hiVect();

#ifdef GPU_COMPATIBLE_PROBLEM

#pragma gpu box(box)
            ca_initrad(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi), BL_TO_FORTRAN_ANYD(Rad_new[mfi]),
                       AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(prob_lo));

#else
            RealBox gridloc(grids[mfi.index()], geom.CellSize(), geom.ProbLo());

            BL_FORT_PROC_CALL(CA_INITRAD, ca_initrad)
            (level, cur_time, ARLIM_3D(lo), ARLIM_3D(hi), Radiation::nGroups,
             BL_TO_FORTRAN_ANYD(Rad_new[mfi]), ZFILL(dx), ZFILL(gridloc.lo()), ZFILL(gridloc.hi()));

#endif
        }
    }
#endif  // RADIATION

#endif  // MAESTRO_INIT

#ifdef GRAVITY
#if (BL_SPACEDIM > 1)
    if ((level == 0) && (spherical_star == 1)) {
        const int nc = S_new.nComp();
        const int n1d = get_numpts();
        allocate_outflow_data(&n1d, &nc);
        int is_new = 1;
        make_radial_data(is_new);
    }
#endif

    MultiFab& G_new = get_new_data(Gravity_Type);
    G_new.setVal(0.);

    MultiFab& phi_new = get_new_data(PhiGrav_Type);
    phi_new.setVal(0.);
#endif

    MultiFab& source_new = get_new_data(Source_Type);
    source_new.setVal(0., source_new.nGrow());

#ifdef ROTATION
    MultiFab& rot_new = get_new_data(Rotation_Type);
    rot_new.setVal(0.);

    MultiFab& phirot_new = get_new_data(PhiRot_Type);
    phirot_new.setVal(0.);
#endif

#ifdef AMREX_PARTICLES
    if (level == 0) {
        init_particles();
    }
#endif

#ifdef AMREX_USE_CUDA
    AMREX_GPU_SAFE_CALL(cudaProfilerStart());
#endif

    if (verbose && ParallelDescriptor::IOProcessor()) {
        std::cout << "Done initializing the level " << level << " data " << std::endl;
    }
}

void Castro::init(AmrLevel& old) {
    BL_PROFILE("Castro::init(old)");

    Castro* oldlev = (Castro*)&old;

    //
    // Create new grid data by fillpatching from old.
    //
    Real dt_new = parent->dtLevel(level);
    Real cur_time = oldlev->state[State_Type].curTime();
    Real prev_time = oldlev->state[State_Type].prevTime();
    Real dt_old = cur_time - prev_time;
    setTimeLevel(cur_time, dt_old, dt_new);

    for (int s = 0; s < num_state_type; ++s) {
        MultiFab& state_MF = get_new_data(s);
        FillPatch(old, state_MF, state_MF.nGrow(), cur_time, s, 0, state_MF.nComp());
    }
}

//
// This version inits the data on a new level that did not
// exist before regridding.
//
void Castro::init() {
    BL_PROFILE("Castro::init()");

    Real dt = parent->dtLevel(level);
    Real cur_time = getLevel(level - 1).state[State_Type].curTime();
    Real prev_time = getLevel(level - 1).state[State_Type].prevTime();

    Real dt_old = (cur_time - prev_time) / (Real)parent->MaxRefRatio(level - 1);

    Real time = cur_time;

    // If we just triggered a regrid, we need to account for the fact that
    // the data on the coarse level has already been advanced.

    if (getLevel(level - 1).post_step_regrid) {
        time = prev_time;
    }

    setTimeLevel(time, dt_old, dt);

    for (int s = 0; s < num_state_type; ++s) {
        MultiFab& state_MF = get_new_data(s);
        FillCoarsePatch(state_MF, 0, time, s, 0, state_MF.nComp());
    }
}

Real Castro::initialTimeStep() {
    Real dummy_dt = 0.0;
    Real init_dt = 0.0;

    if (initial_dt > 0.0) {
        init_dt = initial_dt;
    } else {
        init_dt = init_shrink * estTimeStep(dummy_dt);
    }

    return init_dt;
}

Real Castro::estTimeStep(Real dt_old) {
    BL_PROFILE("Castro::estTimeStep()");

    if (fixed_dt > 0.0) {
        return fixed_dt;
    }

    ca_set_amr_info(level, -1, -1, -1.0, -1.0);

    Real estdt = max_dt;

    Real time = state[State_Type].curTime();

    const Real* dx = geom.CellSize();

    std::string limiter = "castro.max_dt";

    // Start the hydro with the max_dt value, but divide by CFL
    // to account for the fact that we multiply by it at the end.
    // This ensures that if max_dt is more restrictive than the hydro
    // criterion, we will get exactly max_dt for a timestep.

    Real estdt_hydro = max_dt / cfl;

    if (do_hydro) {

#ifdef RADIATION
        if (Radiation::rad_hydro_combined) {

            const MultiFab& stateMF = get_new_data(State_Type);

            // Compute radiation + hydro limited timestep.

#ifdef _OPENMP
#pragma omp parallel reduction(min : estdt_hydro)
#endif
            {
                Real dt = max_dt / cfl;

                const MultiFab& radMF = get_new_data(Rad_Type);
                FArrayBox gPr;

                for (MFIter mfi(stateMF, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const Box& tbox = mfi.tilebox();
                    const Box& vbox = mfi.validbox();

                    gPr.resize(tbox);
                    radiation->estimate_gamrPr(stateMF[mfi], radMF[mfi], gPr, dx, vbox);

                    ca_estdt_rad(tbox.loVect(), tbox.hiVect(), BL_TO_FORTRAN(stateMF[mfi]),
                                 BL_TO_FORTRAN(gPr), dx, &dt);
                }
                estdt_hydro = std::min(estdt_hydro, dt);
            }

        } else {
#endif

#ifdef MHD
            estdt_hydro = estdt_mhd();
#else
        estdt_hydro = estdt_cfl(time);
#endif

#ifdef RADIATION
        }
#endif

        ParallelDescriptor::ReduceRealMin(estdt_hydro);
        estdt_hydro *= cfl;
        if (verbose) {
            amrex::Print() << "...estimated hydro-limited timestep at level " << level << ": "
                           << estdt_hydro << std::endl;
        }

        // Determine if this is more restrictive than the maximum timestep limiting

        if (estdt_hydro < estdt) {
            limiter = "hydro";
            estdt = estdt_hydro;
        }
    }

#ifdef DIFFUSION
    // Diffusion-limited timestep
    // Note that the diffusion uses the same CFL safety factor
    // as the main hydrodynamics timestep limiter.

    Real estdt_diffusion = max_dt / cfl;

    if (diffuse_temp) {
        estdt_diffusion = estdt_temp_diffusion();
    }

    ParallelDescriptor::ReduceRealMin(estdt_diffusion);
    estdt_diffusion *= cfl;
    if (verbose) {
        amrex::Print() << "...estimated diffusion-limited timestep at level " << level << ": "
                       << estdt_diffusion << std::endl;
    }

    // Determine if this is more restrictive than the hydro limiting

    if (estdt_diffusion < estdt) {
        limiter = "diffusion";
        estdt = estdt_diffusion;
    }
#endif  // diffusion

#ifdef REACTIONS
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& R_new = get_new_data(Reactions_Type);

    // Dummy value to start with
    Real estdt_burn = max_dt;

    if (do_react) {

        // Compute burning-limited timestep.

#ifdef _OPENMP
#pragma omp parallel reduction(min : estdt_burn)
#endif
        {
            Real dt = max_dt;

            for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
                const Box& box = mfi.validbox();

#pragma gpu box(box)
                ca_estdt_burning(AMREX_INT_ANYD(box.loVect()), AMREX_INT_ANYD(box.hiVect()),
                                 BL_TO_FORTRAN_ANYD(S_new[mfi]), BL_TO_FORTRAN_ANYD(R_new[mfi]),
                                 AMREX_REAL_ANYD(dx), AMREX_MFITER_REDUCE_MIN(&dt));
            }

            estdt_burn = std::min(estdt_burn, dt);
        }

        ParallelDescriptor::ReduceRealMin(estdt_burn);

        if (verbose && estdt_burn < max_dt) {
            amrex::Print() << "...estimated burning-limited timestep at level " << level << ": "
                           << estdt_burn << std::endl;
        }

        // Determine if this is more restrictive than the hydro limiting

        if (estdt_burn < estdt) {
            limiter = "burning";
            estdt = estdt_burn;
        }
    }
#endif

#ifdef RADIATION
    if (do_radiation) radiation->EstTimeStep(estdt, level);
#endif

    if (verbose) {
        amrex::Print() << "Castro::estTimeStep (" << limiter << "-limited) at level " << level
                       << ":  estdt = " << estdt << '\n'
                       << std::endl;
    }

    return estdt;
}

void Castro::computeNewDt(int finest_level, int sub_cycle, Vector<int>& n_cycle,
                          const Vector<IntVect>& ref_ratio, Vector<Real>& dt_min,
                          Vector<Real>& dt_level, Real stop_time, int post_regrid_flag) {
    BL_PROFILE("Castro::computeNewDt()");

    //
    // We are at the start of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0) return;

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++) {
        Castro& adv_level = getLevel(i);
        dt_min[i] = adv_level.estTimeStep(dt_level[i]);
    }

    if (fixed_dt <= 0.0) {
        if (post_regrid_flag == 1) {
            //
            // Limit dt's by pre-regrid dt
            //
            for (int i = 0; i <= finest_level; i++) {
                dt_min[i] = std::min(dt_min[i], dt_level[i]);
            }
        } else {
            //
            // Limit dt's by change_max * old dt,
            // if we didn't limit the last timestep
            // to hit a plotfile interval.
            //
            if (lastDtPlotLimited) {

                dt_min[0] = std::min(dt_min[0], lastDtBeforePlotLimiting);

                lastDtPlotLimited = 0;
                lastDtBeforePlotLimiting = 0.0;

            } else {

                for (int i = 0; i <= finest_level; i++) {
                    if (verbose && ParallelDescriptor::IOProcessor()) {
                        if (dt_min[i] > change_max * dt_level[i]) {
                            std::cout << "Castro::compute_new_dt : limiting dt at level " << i
                                      << '\n';
                            std::cout << " ... new dt computed: " << dt_min[i] << '\n';
                            std::cout << " ... but limiting to: " << change_max * dt_level[i]
                                      << " = " << change_max << " * " << dt_level[i] << '\n';
                        }
                    }
                    dt_min[i] = std::min(dt_min[i], change_max * dt_level[i]);
                }
            }
        }
    }

    //
    // If we limited the last step by a retry,
    // apply that here if the retry-recommended
    // timestep is smaller than what we calculated.
    //
    for (int i = 0; i <= finest_level; ++i) {
        if (getLevel(i).lastDtRetryLimited == 1) {
            if (getLevel(i).lastDtFromRetry < dt_min[i]) {
                if (verbose && ParallelDescriptor::IOProcessor()) {
                    std::cout << " ... limiting dt at level " << i
                              << " to: " << getLevel(i).lastDtFromRetry
                              << " = retry-limited timestep\n";
                }
                dt_min[i] = getLevel(i).lastDtFromRetry;
            }
        }
    }

    //
    // Find the minimum over all levels
    //
    for (int i = 0; i <= finest_level; i++) {
        n_factor *= n_cycle[i];
        dt_0 = std::min(dt_0, n_factor * dt_min[i]);
    }

    //
    // Optionally, limit dt's by the value of
    // plot_per or small_plot_per.
    //
    if (plot_per_is_exact) {

        const Real plot_per = parent->plotPer();

        if (plot_per > 0.0) {

            const Real cur_time = state[State_Type].curTime();

            // Calculate the new dt by comparing to the dt needed to get
            // to the next multiple of plot_per.

            const Real dtMod = std::fmod(cur_time, plot_per);

            Real newPlotDt;

            // Note that if we are just about exactly on a multiple of plot_per,
            // then we need to be careful to avoid floating point issues.

            if (std::abs(dtMod - plot_per) <= std::numeric_limits<Real>::epsilon() * cur_time) {
                newPlotDt = plot_per + (plot_per - dtMod);
            } else {
                newPlotDt = plot_per - dtMod;
            }

            if (newPlotDt < dt_0) {
                lastDtPlotLimited = 1;
                lastDtBeforePlotLimiting = dt_0;
                dt_0 = newPlotDt;

                // Avoid taking timesteps that are so small that
                // they may cause problems in the hydrodynamics.

                const Real epsDt = 1.e-4 * lastDtBeforePlotLimiting;
                dt_0 = std::max(dt_0, epsDt);

                if (verbose) {
                    amrex::Print()
                        << " ... limiting dt to " << dt_0 << " to hit the next plot interval.\n";
                }
            }
        }
    }

    if (small_plot_per_is_exact) {

        const Real small_plot_per = parent->smallplotPer();

        if (small_plot_per > 0.0) {

            const Real cur_time = state[State_Type].curTime();

            // Same logic as for plot_per_is_exact.

            const Real dtMod = std::fmod(cur_time, small_plot_per);

            Real newSmallPlotDt;

            if (std::abs(dtMod - small_plot_per) <=
                std::numeric_limits<Real>::epsilon() * cur_time) {
                newSmallPlotDt = small_plot_per + (small_plot_per - dtMod);
            } else {
                newSmallPlotDt = small_plot_per - dtMod;
            }

            if (newSmallPlotDt < dt_0) {
                lastDtPlotLimited = 1;
                lastDtBeforePlotLimiting = dt_0;
                dt_0 = newSmallPlotDt;

                const Real epsDt = 1.e-4 * lastDtBeforePlotLimiting;
                dt_0 = std::max(dt_0, epsDt);

                if (verbose)
                    amrex::Print() << " ... limiting dt to " << dt_0
                                   << " to hit the next smallplot interval.\n";
            }
        }
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = std::numeric_limits<Real>::epsilon();
    Real cur_time = state[State_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) >= (stop_time - eps)) {
            dt_0 = stop_time - cur_time;
            if (verbose)
                amrex::Print() << " ... limiting dt to " << dt_0 << " to hit the stop_time.\n";
        }
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++) {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0 / n_factor;
    }
}

void Castro::computeInitialDt(int finest_level, int sub_cycle, Vector<int>& n_cycle,
                              const Vector<IntVect>& ref_ratio, Vector<Real>& dt_level,
                              Real stop_time) {
    BL_PROFILE("Castro::computeInitialDt()");

    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0) {
        return;
    }

    int i;

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    // TODO/DEBUG: This will need to change for optimal subcycling.
    for (i = 0; i <= finest_level; i++) {
        dt_level[i] = getLevel(i).initialTimeStep();
        n_factor *= n_cycle[i];
        dt_0 = std::min(dt_0, n_factor * dt_level[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001 * dt_0;
    Real cur_time = state[State_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps)) {
            dt_0 = stop_time - cur_time;
        }
    }

    n_factor = 1;
    for (i = 0; i <= finest_level; i++) {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0 / n_factor;
    }
}

void Castro::post_timestep(int iteration) {
    BL_PROFILE("Castro::post_timestep()");

    // Pass some information about the state of the simulation to a Fortran module.

    ca_set_amr_info(level, iteration, -1, -1.0, -1.0);

    //
    // Integration cycle on fine level grids is complete
    // do post_timestep stuff here.
    //
    int finest_level = parent->finestLevel();

#ifdef RADIATION
    if (do_radiation && (level < finest_level)) {
        // computeTemp is not needed before or after this call because
        // setup for deferred sync does not touch state, only flux registers.
        radiation->deferred_sync_setup(level);

        if (do_reflux) {
            radiation->reflux(level);
            // Since radiation->reflux does not touch the fluid state,
            // we do need to recompute Temp here.
        }
    }
#endif

    // Now do the refluxing. If we're using gravity it
    // will also do the sync solve associated with the reflux.

    if (do_reflux && level < parent->finestLevel()) {
        reflux(level, level + 1);
    }

    // Ensure consistency with finer grids.

    if (level < finest_level) avgDown();

#ifdef MHD
    MultiFab& Bx_new = get_new_data(Mag_Type_x);
    MultiFab& By_new = get_new_data(Mag_Type_y);
    MultiFab& Bz_new = get_new_data(Mag_Type_z);
#endif

    // Clean up any aberrant state data generated by the reflux and average-down,
    // and then update quantities like temperature to be consistent.
    MultiFab& S_new = get_new_data(State_Type);
    clean_state(
#ifdef MHD
        Bx_new, By_new, Bz_new,
#endif
        S_new, state[State_Type].curTime(), S_new.nGrow());

    // Flush Fortran output

    if (verbose) flush_output();

#ifdef DO_PROBLEM_POST_TIMESTEP

    // Provide a hook for the user to do things after all of
    // the normal updates have been applied. The user is
    // responsible for any actions after this point, like
    // doing a computeTemp call if they change the state data.

    problem_post_timestep();

#endif

    if (level == 0) {
        int nstep = parent->levelSteps(0);
        Real dtlev = parent->dtLevel(0);
        Real cumtime = parent->cumTime() + dtlev;

        bool sum_int_test = false;

        if (sum_interval > 0) {

            if (nstep % sum_interval == 0) {
                sum_int_test = true;
            }
        }

        bool sum_per_test = false;

        if (sum_per > 0.0) {

            const int num_per_old = floor((cumtime - dtlev) / sum_per);
            const int num_per_new = floor((cumtime) / sum_per);

            if (num_per_old != num_per_new) {
                sum_per_test = true;
            }
        }

        if (sum_int_test || sum_per_test) {
            sum_integrated_quantities();
        }

#ifdef GRAVITY
        if (moving_center) {
            write_center();
        }
#endif
    }

#ifdef RADIATION
    // diagnostic stuff

    if (level == 0) {
        do_energy_diagnostics();
    }
#endif

#ifdef AMREX_PARTICLES
    if (TracerPC) {
        const int ncycle = parent->nCycle(level);
        //
        // Don't redistribute/timestamp on the final subiteration except on the coarsest grid.
        //
        if (iteration < ncycle || level == 0) {
            int ngrow = (level == 0) ? 0 : iteration;

            TracerPC->Redistribute(level, parent->finestLevel(), ngrow);

            TimestampParticles(ngrow + 1);
        }
    }
#endif
}

void Castro::post_restart() {
    BL_PROFILE("Castro::post_restart()");

    Real cur_time = state[State_Type].curTime();

#ifdef AMREX_PARTICLES
    ParticlePostRestart(parent->theRestartFile());
#endif

#ifdef GRAVITY
    if (do_grav) {
        if (level == 0) {
            // Passing numpts_1d at level 0
            int numpts_1d = get_numpts();

#if (BL_SPACEDIM == 1)
            numpts_1d += 2 * NUM_GROW;
#endif

            gravity->set_numpts_in_gravity(numpts_1d);

            for (int lev = 0; lev <= parent->finestLevel(); lev++) {
                AmrLevel& this_level = getLevel(lev);
                Castro& cs_level = getLevel(lev);
                gravity->install_level(lev, &this_level, cs_level.Volume(), cs_level.Area());
            }

            if (moving_center == 1) {
                MultiFab& S_new = get_new_data(State_Type);
                define_new_center(S_new, cur_time);
            }

            gravity->set_mass_offset(cur_time);

            if (gravity->get_gravity_type() == "PoissonGrav") {
                if (gravity->NoComposite() != 1) {
                    int use_previous_phi = 1;

                    // Update the maximum density, used in setting the solver tolerance.

                    gravity->update_max_rhs();

                    gravity->multilevel_solve_for_new_phi(0, parent->finestLevel(),
                                                          use_previous_phi);
                    if (gravity->test_results_of_solves() == 1) gravity->test_composite_phi(level);
                }
            }

            if (grown_factor > 1) post_grown_restart();
        }
    }
#endif

#ifdef ROTATION
    MultiFab& phirot_new = get_new_data(PhiRot_Type);
    MultiFab& rot_new = get_new_data(Rotation_Type);
    MultiFab& S_new = get_new_data(State_Type);
    if (do_rotation) {
        fill_rotation_field(phirot_new, rot_new, S_new, cur_time);
    } else {
        phirot_new.setVal(0.0);
        rot_new.setVal(0.0);
    }
#endif

#ifdef DIFFUSION
    // diffusion is a static object, only alloc if not already there
    if (diffusion == 0) diffusion = new Diffusion(parent, &phys_bc);

    if (level == 0)
        for (int lev = 0; lev <= parent->finestLevel(); lev++) {
            AmrLevel& this_level = getLevel(lev);
            Castro& cs_level = getLevel(lev);
            diffusion->install_level(lev, &this_level, cs_level.Volume(), cs_level.Area());
        }
#endif

#ifdef DO_PROBLEM_POST_RESTART
    problem_post_restart();
#endif
}

void Castro::postCoarseTimeStep(Real cumtime) {
    // postCoarseTimeStep() is only called by level 0.
    BL_ASSERT(level == 0);
    AmrLevel::postCoarseTimeStep(cumtime);
#ifdef GRAVITY
    if (do_grav) gravity->set_mass_offset(cumtime, 0);
#endif
}

void Castro::check_for_post_regrid(Real time) {

    BL_PROFILE("Castro::check_for_post_regrid()");

    // Check whether we have any zones at this time signifying that they
    // need to be tagged that do not have corresponding zones on the
    // fine level.

    if (level < parent->maxLevel()) {

        TagBoxArray tags(grids, dmap);

        for (int i = 0; i < err_list_names.size(); ++i) {
            apply_tagging_func(tags, time, i);
        }

        apply_problem_tags(tags, time);

        // Globally collate the tags.

        Vector<IntVect> tvec;

        tags.collate(tvec);

        // If we requested any tags at all, we have a potential trigger for a regrid.

        int num_tags = tvec.size();

        bool missing_on_fine_level = false;

        if (num_tags > 0) {

            if (level == parent->finestLevel()) {

                // If there is no level above us at all, we know a regrid is needed.

                missing_on_fine_level = true;

            } else {

                // If there is a level above us, we need to check whether
                // every tagged zone has a corresponding entry in the fine
                // grid. If not, we need to trigger a regrid so we can get
                // those other zones refined too.

                const BoxArray& fgrids = getLevel(level + 1).boxArray();

                for (int i = 0; i < tvec.size(); ++i) {

                    Box c_bx(tvec[i], tvec[i]);
                    Box f_bx = c_bx.refine(parent->refRatio(level));

                    if (!fgrids.contains(f_bx)) {
                        missing_on_fine_level = true;
                        break;
                    }
                }
            }
        }

        if (missing_on_fine_level) {
            post_step_regrid = 1;

            if (amrex::ParallelDescriptor::IOProcessor()) {
                std::cout << "\n"
                          << "Current refinement pattern insufficient at level " << level << ".\n"
                          << "Performing a regrid to obtain more refinement.\n";
            }
        }
    }
}

void Castro::post_regrid(int lbase, int new_finest) {

    BL_PROFILE("Castro::post_regrid()");

    fine_mask.clear();

#ifdef AMREX_PARTICLES
    if (TracerPC && level == lbase) {
        TracerPC->Redistribute(lbase);
    }
#endif

#ifdef GRAVITY
    if (do_grav) {

        if (use_post_step_regrid && getLevel(lbase).post_step_regrid &&
            gravity->get_gravity_type() == "PoissonGrav") {

            if (level > lbase) {

                // In the case where we're coming here during a regrid that occurs
                // after a timestep, we only want to interpolate the gravitational
                // field from the old time. The state data will already have been
                // filled, so all we need to do is interpolate the grad_phi data.

                // Instantiate a bare physical BC function for grad_phi. It doesn't do anything
                // since the fine levels for Poisson gravity do not touch the physical boundary.

                GradPhiPhysBCFunct gp_phys_bc;

                // We need to use a nodal interpolater.

                Interpolater* gp_interp = &node_bilinear_interp;

                Vector<MultiFab*> grad_phi_coarse =
                    amrex::GetVecOfPtrs(gravity->get_grad_phi_prev(level - 1));
                Vector<MultiFab*> grad_phi_fine =
                    amrex::GetVecOfPtrs(gravity->get_grad_phi_curr(level));

                Real time = getLevel(lbase).get_state_data(Gravity_Type).prevTime();

                // For the BCs, we will use the Gravity_Type BCs for convenience, but these will
                // not do anything because we do not fill on physical boundaries.

                const Vector<BCRec>& gp_bcs = getLevel(level).get_desc_lst()[Gravity_Type].getBCs();

                for (int n = 0; n < BL_SPACEDIM; ++n) {
                    amrex::InterpFromCoarseLevel(*grad_phi_fine[n], time, *grad_phi_coarse[n], 0, 0,
                                                 1, parent->Geom(level - 1), parent->Geom(level),
                                                 gp_phys_bc, 0, gp_phys_bc, 0,
                                                 parent->refRatio(level - 1), gp_interp, gp_bcs, 0);
                }
            }

        } else {

            const Real cur_time = state[State_Type].curTime();
            if ((level == lbase) && cur_time > 0.) {
                if (gravity->get_gravity_type() == "PoissonGrav" && (gravity->NoComposite() != 1)) {
                    int use_previous_phi = 1;

                    // Update the maximum density, used in setting the solver tolerance.

                    if (level == 0) {
                        gravity->update_max_rhs();
                    }

                    gravity->multilevel_solve_for_new_phi(level, new_finest, use_previous_phi);
                }
            }
        }
    }
#endif
}

void Castro::post_init(Real stop_time) {
    BL_PROFILE("Castro::post_init()");

    if (level > 0) {
        return;
    }

    //
    // Average data down from finer levels
    // so that conserved data is consistent between levels.
    //
    int finest_level = parent->finestLevel();
    for (int k = finest_level - 1; k >= 0; k--) {
        getLevel(k).avgDown();
    }

#ifdef GRAVITY

    if (do_grav) {

        Real cur_time = state[State_Type].curTime();

        if (gravity->get_gravity_type() == "PoissonGrav") {

            // Update the maximum density, used in setting the solver tolerance.

            gravity->update_max_rhs();

            // Calculate offset before first multilevel solve.
            gravity->set_mass_offset(cur_time);

            if (gravity->NoComposite() != 1) {
                gravity->multilevel_solve_for_new_phi(level, finest_level);
                if (gravity->test_results_of_solves() == 1) {
                    gravity->test_composite_phi(level);
                }
            }
        }

        // Make this call just to fill the initial state data.
        for (int k = 0; k <= parent->finestLevel(); k++) {
            BoxArray ba = getLevel(k).boxArray();
            MultiFab& grav_new = getLevel(k).get_new_data(Gravity_Type);
            gravity->get_new_grav_vector(k, grav_new, cur_time);
        }
    }
#endif

#ifdef ROTATION
    MultiFab& phirot_new = get_new_data(PhiRot_Type);
    MultiFab& rot_new = get_new_data(Rotation_Type);
    MultiFab& S_new = get_new_data(State_Type);
    if (do_rotation) {
        Real cur_time = state[State_Type].curTime();
        fill_rotation_field(phirot_new, rot_new, S_new, cur_time);
    } else {
        phirot_new.setVal(0.0);
        rot_new.setVal(0.0);
    }
#endif

#ifdef RADIATION
    if (do_radiation) {
        // The option of whether to do a multilevel initialization is
        // controlled within the radiation class.

        radiation->post_init(level);

        for (int k = finest_level - 1; k >= 0; k--) {
            getLevel(k).avgDown(Rad_Type);
        }

        do_energy_diagnostics();

        // re-estimate the initial timestep using the initialized
        // fine-level data and radiation field.
        //if (level == 0)
        //post_init_estDT();
    }
#endif

    // Allow the user to define their own post_init functions.

#ifdef DO_PROBLEM_POST_INIT

    problem_post_init();

#endif

    // If we're doing SCF initialization, do it here.

#ifdef GRAVITY
#ifdef ROTATION
#ifndef AMREX_USE_CUDA
    if (do_scf_initial_model) {
        scf_relaxation();
    }
#endif
#endif
#endif

    int nstep = parent->levelSteps(0);
    Real dtlev = parent->dtLevel(0);
    Real cumtime = parent->cumTime();
    if (cumtime != 0.0) cumtime += dtlev;

    bool sum_int_test = false;

    if (sum_interval > 0) {

        if (nstep % sum_interval == 0) {
            sum_int_test = true;
        }
    }

    bool sum_per_test = false;

    if (sum_per > 0.0) {

        const int num_per_old = floor((cumtime - dtlev) / sum_per);
        const int num_per_new = floor((cumtime) / sum_per);

        if (num_per_old != num_per_new) {
            sum_per_test = true;
        }
    }

    if (sum_int_test || sum_per_test) {
        sum_integrated_quantities();
    }

#ifdef GRAVITY
    if (level == 0 && moving_center == 1) {
        write_center();
    }
#endif
}

void Castro::post_grown_restart() {

    BL_PROFILE("Castro::post_grown_restart()");

    if (level > 0) return;

#ifdef GRAVITY
    if (do_grav) {
        int finest_level = parent->finestLevel();
        Real cur_time = state[State_Type].curTime();

        if (gravity->get_gravity_type() == "PoissonGrav") {

            // Update the maximum density, used in setting the solver tolerance.

            gravity->update_max_rhs();

            // Calculate offset before first multilevel solve.
            gravity->set_mass_offset(cur_time);

            if (gravity->NoComposite() != 1) {
                gravity->multilevel_solve_for_new_phi(level, finest_level);
                if (gravity->test_results_of_solves() == 1) {
                    gravity->test_composite_phi(level);
                }
            }
        }

        // Make this call just to fill the initial state data.
        for (int k = 0; k <= parent->finestLevel(); k++) {
            MultiFab& grav_new = getLevel(k).get_new_data(Gravity_Type);
            gravity->get_new_grav_vector(k, grav_new, cur_time);
        }
    }
#endif

#ifdef ROTATION
    MultiFab& phirot_new = get_new_data(PhiRot_Type);
    MultiFab& rot_new = get_new_data(Rotation_Type);
    MultiFab& S_new = get_new_data(State_Type);
    if (do_rotation) {
        Real cur_time = state[State_Type].curTime();
        fill_rotation_field(phirot_new, rot_new, S_new, cur_time);
    } else {
        phirot_new.setVal(0.0);
        rot_new.setVal(0.0);
    }
#endif

#ifdef RADIATION
    if (do_radiation) {
        // The option of whether to do a multilevel initialization is
        // controlled within the radiation class.

        radiation->post_init(level);

        int finest_level = parent->finestLevel();

        for (int k = finest_level - 1; k >= 0; k--) {
            getLevel(k).avgDown(Rad_Type);
        }

        do_energy_diagnostics();

        // re-estimate the initial timestep using the initialized
        // fine-level data and radiation field.
        //if (level == 0)
        //post_init_estDT();
    }
#endif
}

int Castro::okToContinue() {
    if (level > 0) {
        return 1;
    }

    int test = 1;

    if (signalStopJob) {
        test = 0;
        if (ParallelDescriptor::IOProcessor()) {
            std::cout << " Signalling a stop of the run due to signalStopJob = true." << std::endl;
        }
    } else if (parent->dtLevel(level) < dt_cutoff * parent->cumTime()) {
        test = 0;
        if (ParallelDescriptor::IOProcessor()) {
            std::cout << " Signalling a stop of the run because dt < dt_cutoff * time."
                      << std::endl;
        }
    }

    return test;
}

#ifdef AUX_UPDATE
void Castro::advance_aux(Real time, Real dt) {
    BL_PROFILE("Castro::advance_aux()");

    if (verbose && ParallelDescriptor::IOProcessor()) {
        std::cout << "... special update for auxiliary variables \n";
    }

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_old, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& box = mfi.tilebox();
        FArrayBox& old_fab = S_old[mfi];
        FArrayBox& new_fab = S_new[mfi];
        void ca_auxupdate(BL_TO_FORTRAN(old_fab), BL_TO_FORTRAN(new_fab), box.loVect(),
                          box.hiVect(), &dt);
    }
}
#endif

void Castro::FluxRegCrseInit() {

    if (level == parent->finestLevel()) {
        return;
    }

    Castro& fine_level = getLevel(level + 1);

    for (int i = 0; i < BL_SPACEDIM; ++i) {
        fine_level.flux_reg.CrseInit(*fluxes[i], i, 0, 0, NUM_STATE, flux_crse_scale);
    }

#if (BL_SPACEDIM <= 2)
    if (!Geom().IsCartesian()) {
        fine_level.pres_reg.CrseInit(P_radial, 0, 0, 0, 1, pres_crse_scale);
    }
#endif

#ifdef RADIATION
    if (Radiation::rad_hydro_combined) {
        for (int i = 0; i < BL_SPACEDIM; ++i) {
            fine_level.rad_flux_reg.CrseInit(*rad_fluxes[i], i, 0, 0, Radiation::nGroups,
                                             flux_crse_scale);
        }
    }
#endif
}

void Castro::FluxRegFineAdd() {

    BL_PROFILE("Castro::FluxRegFineAdd()");

    if (level == 0) {
        return;
    }

    for (int i = 0; i < BL_SPACEDIM; ++i) {
        flux_reg.FineAdd(*fluxes[i], i, 0, 0, NUM_STATE, flux_fine_scale);
    }

#if (BL_SPACEDIM <= 2)
    if (!Geom().IsCartesian()) {
        getLevel(level).pres_reg.FineAdd(P_radial, 0, 0, 0, 1, pres_fine_scale);
    }
#endif

#ifdef RADIATION
    if (Radiation::rad_hydro_combined) {
        for (int i = 0; i < BL_SPACEDIM; ++i) {
            getLevel(level).rad_flux_reg.FineAdd(*rad_fluxes[i], i, 0, 0, Radiation::nGroups,
                                                 flux_fine_scale);
        }
    }
#endif
}

void Castro::reflux(int crse_level, int fine_level) {
    BL_PROFILE("Castro::reflux()");

    BL_ASSERT(fine_level > crse_level);

    const Real strt = ParallelDescriptor::second();

#ifdef GRAVITY
    int nlevs = fine_level - crse_level + 1;

    Vector<std::unique_ptr<MultiFab> > drho(nlevs);
    Vector<std::unique_ptr<MultiFab> > dphi(nlevs);

    if (do_grav && gravity->get_gravity_type() == "PoissonGrav" && gravity->NoSync() == 0) {

        for (int lev = crse_level; lev <= fine_level; ++lev) {

            const auto& amrlevel = getLevel(lev);
            const auto& ba = amrlevel.boxArray();
            const auto& dm = amrlevel.DistributionMap();

            drho[lev - crse_level].reset(new MultiFab(ba, dm, 1, 0));
            dphi[lev - crse_level].reset(new MultiFab(ba, dm, 1, 0));

            drho[lev - crse_level]->setVal(0.0);
            dphi[lev - crse_level]->setVal(0.0);
        }
    }
#endif

    FluxRegister* reg;

    for (int lev = fine_level; lev > crse_level; --lev) {

        reg = &getLevel(lev).flux_reg;

        Castro& crse_lev = getLevel(lev - 1);
        Castro& fine_lev = getLevel(lev);

        MultiFab& crse_state = crse_lev.get_new_data(State_Type);

        // Clear out the data that's not on coarse-fine boundaries so that this register only
        // modifies the fluxes on coarse-fine interfaces.

        reg->ClearInternalBorders(crse_lev.geom);

        // Trigger the actual reflux on the coarse level now.

        reg->Reflux(crse_state, crse_lev.volume, 1.0, 0, 0, NUM_STATE, crse_lev.geom);

        // Store the density change, for the gravity sync.

#ifdef GRAVITY
        int ilev = lev - crse_level - 1;

        if (do_grav && gravity->get_gravity_type() == "PoissonGrav" && gravity->NoSync() == 0) {
            reg->Reflux(*drho[ilev], crse_lev.volume, 1.0, 0, URHO, 1, crse_lev.geom);
            amrex::average_down(*drho[ilev + 1], *drho[ilev], 0, 1, getLevel(lev).crse_ratio);
        }
#endif

        // Also update the coarse fluxes MultiFabs using the reflux data. This should only make
        // a difference if we re-evaluate the source terms later.

        Vector<std::unique_ptr<MultiFab> > temp_fluxes(3);

        if (update_sources_after_reflux) {

            for (int i = 0; i < BL_SPACEDIM; ++i) {
                temp_fluxes[i].reset(new MultiFab(
                    crse_lev.fluxes[i]->boxArray(), crse_lev.fluxes[i]->DistributionMap(),
                    crse_lev.fluxes[i]->nComp(), crse_lev.fluxes[i]->nGrow()));
                temp_fluxes[i]->setVal(0.0);
            }
            for (OrientationIter fi; fi; ++fi) {
                const FabSet& fs = (*reg)[fi()];
                int idir = fi().coordDir();
                fs.copyTo(*temp_fluxes[idir], 0, 0, 0, temp_fluxes[idir]->nComp());
            }
            for (int i = 0; i < BL_SPACEDIM; ++i) {
                MultiFab::Add(*crse_lev.fluxes[i], *temp_fluxes[i], 0, 0,
                              crse_lev.fluxes[i]->nComp(), 0);
                MultiFab::Add(*crse_lev.mass_fluxes[i], *temp_fluxes[i], URHO, 0, 1, 0);
                temp_fluxes[i].reset();
            }
        }

        // We no longer need the flux register data, so clear it out.

        reg->setVal(0.0);

#if (BL_SPACEDIM <= 2)
        if (!Geom().IsCartesian()) {

            reg = &getLevel(lev).pres_reg;

            MultiFab dr(crse_lev.grids, crse_lev.dmap, 1, 0);
            dr.setVal(crse_lev.geom.CellSize(0));

            reg->ClearInternalBorders(crse_lev.geom);

            reg->Reflux(crse_state, dr, 1.0, 0, UMX, 1, crse_lev.geom);

            if (update_sources_after_reflux) {

                temp_fluxes[0].reset(
                    new MultiFab(crse_lev.P_radial.boxArray(), crse_lev.P_radial.DistributionMap(),
                                 crse_lev.P_radial.nComp(), crse_lev.P_radial.nGrow()));
                temp_fluxes[0]->setVal(0.0);

                for (OrientationIter fi; fi; ++fi) {
                    const FabSet& fs = (*reg)[fi()];
                    int idir = fi().coordDir();
                    if (idir == 0) {
                        fs.copyTo(*temp_fluxes[idir], 0, 0, 0, temp_fluxes[idir]->nComp());
                    }
                }

                MultiFab::Add(crse_lev.P_radial, *temp_fluxes[0], 0, 0, crse_lev.P_radial.nComp(),
                              0);
                temp_fluxes[0].reset();
            }

            reg->setVal(0.0);
        }
#endif

#ifdef RADIATION

        // This follows the same logic as the pure hydro fluxes; see above for details.

        if (Radiation::rad_hydro_combined) {

            reg = &getLevel(lev).rad_flux_reg;

            reg->ClearInternalBorders(crse_lev.geom);

            reg->Reflux(crse_lev.get_new_data(Rad_Type), crse_lev.volume, 1.0, 0, 0,
                        Radiation::nGroups, crse_lev.geom);

            if (update_sources_after_reflux) {

                for (int i = 0; i < BL_SPACEDIM; ++i) {
                    temp_fluxes[i].reset(new MultiFab(crse_lev.rad_fluxes[i]->boxArray(),
                                                      crse_lev.rad_fluxes[i]->DistributionMap(),
                                                      crse_lev.rad_fluxes[i]->nComp(),
                                                      crse_lev.rad_fluxes[i]->nGrow()));
                    temp_fluxes[i]->setVal(0.0);
                }
                for (OrientationIter fi; fi; ++fi) {
                    const FabSet& fs = (*reg)[fi()];
                    int idir = fi().coordDir();
                    fs.copyTo(*temp_fluxes[idir], 0, 0, 0, temp_fluxes[idir]->nComp());
                }
                for (int i = 0; i < BL_SPACEDIM; ++i) {
                    MultiFab::Add(*crse_lev.rad_fluxes[i], *temp_fluxes[i], 0, 0,
                                  crse_lev.rad_fluxes[i]->nComp(), 0);
                    temp_fluxes[i].reset();
                }
            }

            reg->setVal(0.0);
        }

#endif

#ifdef GRAVITY
        if (do_grav && gravity->get_gravity_type() == "PoissonGrav" && gravity->NoSync() == 0) {

            reg = &getLevel(lev).phi_reg;

            // Note that the scaling by the area here is corrected for by dividing by the
            // cell volume in the reflux. In this way we get a discrete divergence that
            // is analogous to the divergence of the flux in the hydrodynamics. See Equation
            // 37 in the Castro I paper. The dimensions of dphi are therefore actually
            // phi / cm**2, which makes it correct for the RHS of the Poisson equation.

            for (int i = 0; i < BL_SPACEDIM; ++i) {
                reg->CrseInit(*(gravity->get_grad_phi_curr(lev - 1)[i]), crse_lev.area[i], i, 0, 0,
                              1, -1.0);
                reg->FineAdd(*(gravity->get_grad_phi_curr(lev)[i]), fine_lev.area[i], i, 0, 0, 1,
                             1.0);
            }

            reg->Reflux(*dphi[ilev], crse_lev.volume, 1.0, 0, 0, 1, crse_lev.geom);

            amrex::average_down(*dphi[ilev + 1], *dphi[ilev], 0, 1, getLevel(lev).crse_ratio);

            reg->setVal(0.0);
        }
#endif
    }

    // Do the sync solve across all levels.

#ifdef GRAVITY
    if (do_grav && gravity->get_gravity_type() == "PoissonGrav" && gravity->NoSync() == 0) {
        gravity->gravity_sync(crse_level, fine_level, amrex::GetVecOfPtrs(drho),
                              amrex::GetVecOfPtrs(dphi));
    }
#endif

    // Now subtract the new-time updates to the state data,
    // recompute it, and add it back. This corrects for the fact
    // that the new-time data was calculated the first time around
    // using a state that hadn't yet been refluxed. Note that this
    // needs to come after the gravity sync solve because the gravity
    // source depends on an up-to-date value of phi. We'll do this
    // on the fine level in addition to the coarser levels, because
    // global sources like gravity or source terms that rely on
    // ghost zone fills like diffusion depend on the data in the
    // coarser levels.

    if (update_sources_after_reflux &&
        (time_integration_method == CornerTransportUpwind ||
         time_integration_method == SimplifiedSpectralDeferredCorrections)) {

        for (int lev = fine_level; lev >= crse_level; --lev) {

            MultiFab& S_old = getLevel(lev).get_old_data(State_Type);

            MultiFab& S_new = getLevel(lev).get_new_data(State_Type);
#ifdef MHD
            MultiFab& Bx_new = getLevel(lev).get_new_data(Mag_Type_x);
            MultiFab& By_new = getLevel(lev).get_new_data(Mag_Type_y);
            MultiFab& Bz_new = getLevel(lev).get_new_data(Mag_Type_z);
#endif
            MultiFab& source = getLevel(lev).get_new_data(Source_Type);
            Real time = getLevel(lev).state[State_Type].curTime();
            Real dt_advance_local =
                getLevel(lev)
                    .dt_advance;  // Note that this may be shorter than the full timestep due to subcycling.
            Real dt_amr = parent->dtLevel(lev);  // The full timestep expected by the Amr class.

            ca_set_amr_info(lev, -1, -1, time, dt_advance_local);

            if (getLevel(lev).apply_sources()) {

                getLevel(lev).apply_source_to_state(S_new, source, -dt_advance_local, 0);
                getLevel(lev).clean_state(
#ifdef MHD
                    Bx_new, By_new, Bz_new,
#endif
                    S_new, time, 0);
            }

            // Temporarily restore the last iteration's old data for the purposes of recalculating the corrector.
            // This is only necessary if we've done subcycles on that level.

            if (use_retry && dt_advance_local < dt_amr && getLevel(lev).keep_prev_state) {

                for (int k = 0; k < num_state_type; k++) {

                    if (getLevel(lev).prev_state[k]->hasOldData()) {

                        // Use the new-time data as a temporary buffer. Ideally this would be done
                        // as a pointer swap, but we cannot assume that the distribution mapping
                        // is the same between the current state and the original state, due to
                        // possible regrids that have occurred in between.

                        MultiFab& old = getLevel(lev).get_old_data(k);
                        MultiFab::Copy(getLevel(lev).prev_state[k]->newData(), old, 0, 0,
                                       old.nComp(), old.nGrow());
                        MultiFab::Copy(old, getLevel(lev).prev_state[k]->oldData(), 0, 0,
                                       old.nComp(), old.nGrow());

                        getLevel(lev).state[k].setTimeLevel(time, dt_advance_local, 0.0);
                        getLevel(lev).prev_state[k]->setTimeLevel(time, dt_amr, 0.0);
                    }
                }
            }

            if (getLevel(lev).apply_sources()) {
                bool apply_sources_to_state = true;
                getLevel(lev).do_new_sources(
#ifdef MHD
                    Bx_new, By_new, Bz_new,
#endif
                    source, S_old, S_new, time, dt_advance_local, apply_sources_to_state);
            }

            if (use_retry && dt_advance_local < dt_amr && getLevel(lev).keep_prev_state) {

                for (int k = 0; k < num_state_type; k++) {

                    if (getLevel(lev).prev_state[k]->hasOldData()) {

                        // Now retrieve the original old time data.

                        MultiFab& old = getLevel(lev).get_old_data(k);
                        MultiFab::Copy(old, getLevel(lev).prev_state[k]->newData(), 0, 0,
                                       old.nComp(), old.nGrow());

                        getLevel(lev).state[k].setTimeLevel(time, dt_amr, 0.0);
                        getLevel(lev).prev_state[k]->setTimeLevel(time, dt_advance_local, 0.0);
                    }
                }

                // Now deallocate the old data, it is no longer needed.

                if (lev == 0 || lev > level) {
                    amrex::FillNull(getLevel(lev).prev_state);
                    getLevel(lev).keep_prev_state = false;
                }
            }
        }
    }

    if (verbose) {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real end = ParallelDescriptor::second() - strt;

#ifdef BL_LAZY
        Lazy::QueueReduction([=]() mutable {
#endif
            ParallelDescriptor::ReduceRealMax(end, IOProc);
            if (ParallelDescriptor::IOProcessor()) {
                std::cout << "Castro::reflux() at level " << level << " : time = " << end
                          << std::endl;
            }
#ifdef BL_LAZY
        });
#endif
    }
}

void Castro::avgDown() {
    BL_PROFILE("Castro::avgDown()");

    if (level == parent->finestLevel()) return;

    for (int k = 0; k < num_state_type; k++) {
        avgDown(k);
    }
}

void Castro::normalize_species(MultiFab& S_new, int ng) {
    BL_PROFILE("Castro::normalize_species()");

    Real lsmall_x = small_x;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.growntilebox(ng);

        auto u = S_new.array(mfi);

        // Ensure the species mass fractions are between small_x and 1,
        // then normalize them so that they sum to 1.

        amrex::ParallelFor(bx, [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept {
            Real rhoX_sum = 0.0_rt;

            for (int n = 0; n < NumSpec; ++n) {
                u(i, j, k, UFS + n) = amrex::max(lsmall_x * u(i, j, k, URHO),
                                                 amrex::min(u(i, j, k, URHO), u(i, j, k, UFS + n)));
                rhoX_sum += u(i, j, k, UFS + n);
            }

            Real fac = u(i, j, k, URHO) / rhoX_sum;

            for (int n = 0; n < NumSpec; ++n) {
                u(i, j, k, UFS + n) *= fac;
            }
        });
    }
}

void Castro::enforce_consistent_e(
#ifdef MHD
    MultiFab& Bx, MultiFab& By, MultiFab& Bz,
#endif
    MultiFab& S) {

    BL_PROFILE("Castro::enforce_consistent_e()");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& box = mfi.tilebox();

        auto S_arr = S.array(mfi);

#ifdef MHD
        auto Bx_arr = Bx.array(mfi);
        auto By_arr = By.array(mfi);
        auto Bz_arr = Bz.array(mfi);
#endif

        ParallelFor(box, [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept {
            Real rhoInv = 1.0_rt / S_arr(i, j, k, URHO);
            Real u = S_arr(i, j, k, UMX) * rhoInv;
            Real v = S_arr(i, j, k, UMY) * rhoInv;
            Real w = S_arr(i, j, k, UMZ) * rhoInv;

            S_arr(i, j, k, UEDEN) =
                S_arr(i, j, k, UEINT) + 0.5_rt * S_arr(i, j, k, URHO) * (u * u + v * v + w * w);

#ifdef MHD
            Real bx_cell_c = 0.5_rt * (Bx_arr(i, j, k) + Bx_arr(i + 1, j, k));
            Real by_cell_c = 0.5_rt * (By_arr(i, j, k) + By_arr(i, j + 1, k));
            Real bz_cell_c = 0.5_rt * (Bz_arr(i, j, k) + Bz_arr(i, j, k + 1));

            S_arr(i, j, k, UEDEN) +=
                0.5_rt * (bx_cell_c * bx_cell_c + by_cell_c * by_cell_c + bz_cell_c * bz_cell_c);
#endif
        });
    }
}

void Castro::enforce_min_density(MultiFab& state_in, int ng) {

    BL_PROFILE("Castro::enforce_min_density()");

    // This routine sets the density in state_in to be larger than the
    // density floor.  Note that it will operate everywhere on state_in,
    // including ghost zones.

    MultiFab reset_source;

    if (print_update_diagnostics) {

        // Before we do anything, make a copy of the state.

        reset_source.define(state_in.boxArray(), state_in.DistributionMap(), state_in.nComp(), 0);

        MultiFab::Copy(reset_source, state_in, 0, 0, state_in.nComp(), 0);
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state_in, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.growntilebox(ng);

        do_enforce_minimum_density(bx, state_in.array(mfi), verbose);
    }

    if (print_update_diagnostics) {

        // Evaluate what the effective reset source was.

        MultiFab::Subtract(reset_source, state_in, 0, 0, state_in.nComp(), 0);

        bool local = true;
        Vector<Real> reset_update = evaluate_source_change(reset_source, 1.0, local);

#ifdef BL_LAZY
        Lazy::QueueReduction([=]() mutable {
#endif
            ParallelDescriptor::ReduceRealSum(reset_update.dataPtr(), reset_update.size(),
                                              ParallelDescriptor::IOProcessorNumber());

            if (ParallelDescriptor::IOProcessor()) {
                if (std::abs(reset_update[0]) != 0.0) {
                    std::cout << std::endl
                              << "  Contributions to the state from negative density resets:"
                              << std::endl;

                    print_source_change(reset_update);
                }
            }

#ifdef BL_LAZY
        });
#endif
    }
}

void Castro::avgDown(int state_indx) {
    BL_PROFILE("Castro::avgDown(state_indx)");

    if (level == parent->finestLevel()) return;

    Castro& fine_lev = getLevel(level + 1);

    const Geometry& fgeom = fine_lev.geom;
    const Geometry& cgeom = geom;

    MultiFab& S_crse = get_new_data(state_indx);
    MultiFab& S_fine = fine_lev.get_new_data(state_indx);

    amrex::average_down(S_fine, S_crse, fgeom, cgeom, 0, S_fine.nComp(), fine_ratio);
}

void Castro::allocOldData() {
    MultiFab::RegionTag amrlevel_tag("AmrLevel_Level_" + std::to_string(level));
    MultiFab::RegionTag statedata_tag("StateData_Level_" + std::to_string(level));
    for (int k = 0; k < num_state_type; k++) state[k].allocOldData();
}

void Castro::removeOldData() { AmrLevel::removeOldData(); }

void Castro::errorEst(TagBoxArray& tags, int clearval, int tagval, Real time, int n_error_buf,
                      int ngrow) {
    BL_PROFILE("Castro::errorEst()");

    ca_set_amr_info(level, -1, -1, -1.0, -1.0);

    Real ltime = time;

    // If we are forcing a post-timestep regrid,
    // note that we need to use the new time here,
    // not the old time.

    if (post_step_regrid) {
        ltime = get_state_data(State_Type).curTime();
    }

    // Apply each of the specified tagging functions.

    for (int j = 0; j < num_err_list_default; j++) {
        apply_tagging_func(tags, ltime, j);
    }

    // Now we'll tag any user-specified zones using the full state array.

    apply_problem_tags(tags, ltime);
}

void Castro::apply_problem_tags(TagBoxArray& tags, Real time) {

    BL_PROFILE("Castro::apply_problem_tags()");

    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();

    MultiFab& S_new = get_new_data(State_Type);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        for (MFIter mfi(tags); mfi.isValid(); ++mfi) {
            // tile box
            const Box& bx = mfi.validbox();

            TagBox& tagfab = tags[mfi];

            const int8_t tagval = (int8_t)TagBox::SET;
            const int8_t clearval = (int8_t)TagBox::CLEAR;

#ifdef GPU_COMPATIBLE_PROBLEM
#pragma gpu box(bx)
            set_problem_tags(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                             (int8_t*)BL_TO_FORTRAN_ANYD(tagfab), BL_TO_FORTRAN_ANYD(S_new[mfi]),
                             AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(prob_lo), tagval, clearval, time,
                             level);
#else
            set_problem_tags(AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
                             (int8_t*)BL_TO_FORTRAN_ANYD(tagfab), BL_TO_FORTRAN_ANYD(S_new[mfi]),
                             AMREX_ZFILL(dx), AMREX_ZFILL(prob_lo), tagval, clearval, time, level);
#endif
        }
    }
}

void Castro::apply_tagging_func(TagBoxArray& tags, Real time, int j) {

    BL_PROFILE("Castro::apply_tagging_func()");

    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();

    auto mf = derive(err_list_names[j], time, err_list_ng[j]);

    BL_ASSERT(mf);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        for (MFIter mfi(tags); mfi.isValid(); ++mfi) {
            // FABs
            FArrayBox& datfab = (*mf)[mfi];
            TagBox& tagfab = tags[mfi];

            // tile box
            const Box& bx = mfi.validbox();

            const int* lo = bx.loVect();
            const int* hi = bx.hiVect();
            //
            const int ncomp = datfab.nComp();

            const int8_t tagval = (int8_t)TagBox::SET;
            const int8_t clearval = (int8_t)TagBox::CLEAR;

            if (err_list_names[j] == "density") {
#pragma gpu box(bx)
                ca_denerror(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                            (int8_t*)BL_TO_FORTRAN_ANYD(tagfab), BL_TO_FORTRAN_ANYD(datfab), ncomp,
                            AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(prob_lo), tagval, clearval, time,
                            level);
            } else if (err_list_names[j] == "Temp") {
#pragma gpu box(bx)
                ca_temperror(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                             (int8_t*)BL_TO_FORTRAN_ANYD(tagfab), BL_TO_FORTRAN_ANYD(datfab), ncomp,
                             AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(prob_lo), tagval, clearval, time,
                             level);
            } else if (err_list_names[j] == "pressure") {
#pragma gpu box(bx)
                ca_presserror(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                              (int8_t*)BL_TO_FORTRAN_ANYD(tagfab), BL_TO_FORTRAN_ANYD(datfab),
                              ncomp, AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(prob_lo), tagval,
                              clearval, time, level);
            } else if (err_list_names[j] == "x_velocity" || err_list_names[j] == "y_velocity" ||
                       err_list_names[j] == "z_velocity") {
#pragma gpu box(bx)
                ca_velerror(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                            (int8_t*)BL_TO_FORTRAN_ANYD(tagfab), BL_TO_FORTRAN_ANYD(datfab), ncomp,
                            AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(prob_lo), tagval, clearval, time,
                            level);
            }
#ifdef REACTIONS
            else if (err_list_names[j] == "t_sound_t_enuc") {
#pragma gpu box(bx)
                ca_nucerror(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                            (int8_t*)BL_TO_FORTRAN_ANYD(tagfab), BL_TO_FORTRAN_ANYD(datfab), ncomp,
                            AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(prob_lo), tagval, clearval, time,
                            level);
            } else if (err_list_names[j] == "enuc") {
#pragma gpu box(bx)
                ca_enucerror(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                             (int8_t*)BL_TO_FORTRAN_ANYD(tagfab), BL_TO_FORTRAN_ANYD(datfab), ncomp,
                             AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(prob_lo), tagval, clearval, time,
                             level);
            }
#endif
#ifdef RADIATION
            else if (err_list_names[j] == "rad") {
#pragma gpu box(bx)
                ca_raderror(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                            (int8_t*)BL_TO_FORTRAN_ANYD(tagfab), BL_TO_FORTRAN_ANYD(datfab), ncomp,
                            AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(prob_lo), tagval, clearval, time,
                            level);
            }
#endif
        }
    }
}

std::unique_ptr<MultiFab> Castro::derive(const std::string& name, Real time, int ngrow) {

    BL_PROFILE("Castro::derive()");

#ifdef AMREX_PARTICLES
    return ParticleDerive(name, time, ngrow);
#else
    return AmrLevel::derive(name, time, ngrow);
#endif
}

void Castro::derive(const std::string& name, Real time, MultiFab& mf, int dcomp) {

    BL_PROFILE("Castro::derive()");

    AmrLevel::derive(name, time, mf, dcomp);
}

void Castro::amrinfo_init() { ca_amrinfo_init(); }

void Castro::amrinfo_finalize() { ca_amrinfo_finalize(); }

void Castro::extern_init() {
    // initialize the external runtime parameters -- these will
    // live in the probin

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "reading extern runtime parameters ..." << std::endl;
    }

    const int probin_file_length = probin_file.length();
    Vector<int> probin_file_name(probin_file_length);

    for (int i = 0; i < probin_file_length; i++) {
        probin_file_name[i] = probin_file[i];
    }

    // read them in in Fortran
    ca_extern_init(probin_file_name.dataPtr(), &probin_file_length);

    // grab them from Fortran to C++
    init_extern_parameters();
}

void Castro::reset_internal_energy(const Box& bx,
#ifdef MHD
                                   Array4<Real> const Bx, Array4<Real> const By,
                                   Array4<Real> const Bz,
#endif
                                   Array4<Real> const u) {
    Real lsmall_temp = small_temp;
    Real ldual_energy_eta2 = dual_energy_eta2;

    amrex::ParallelFor(bx, [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept {
        Real rhoInv = 1.0_rt / u(i, j, k, URHO);
        Real Up = u(i, j, k, UMX) * rhoInv;
        Real Vp = u(i, j, k, UMY) * rhoInv;
        Real Wp = u(i, j, k, UMZ) * rhoInv;
        Real ke = 0.5_rt * (Up * Up + Vp * Vp + Wp * Wp);

        eos_t eos_state;

        eos_state.rho = u(i, j, k, URHO);
        eos_state.T = lsmall_temp;
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = u(i, j, k, UFS + n) * rhoInv;
        }
        for (int n = 0; n < NumAux; ++n) {
            eos_state.aux[n] = u(i, j, k, UFX + n) * rhoInv;
        }

        eos(eos_input_rt, eos_state);

        Real small_e = eos_state.e;

#ifdef MHD
        Real bx_cell_c = 0.5_rt * (Bx(i, j, k) + Bx(i + 1, j, k));
        Real by_cell_c = 0.5_rt * (By(i, j, k) + By(i, j + 1, k));
        Real bz_cell_c = 0.5_rt * (Bz(i, j, k) + Bz(i, j, k + 1));

        Real B_ener =
            0.5_rt * (bx_cell_c * bx_cell_c + by_cell_c * by_cell_c + bz_cell_c * bz_cell_c);
#else
        Real B_ener = 0.0_rt;
#endif

        // Ensure the internal energy is at least as large as this minimum
        // from the EOS; the same holds true for the total energy.

        u(i, j, k, UEINT) = amrex::max(u(i, j, k, UEINT), u(i, j, k, URHO) * small_e);
        u(i, j, k, UEDEN) =
            amrex::max(u(i, j, k, UEDEN), u(i, j, k, URHO) * (small_e + ke) + B_ener);

        // Apply the dual energy criterion: get e from E if (E - K) > eta * E.

        Real rho_eint = u(i, j, k, UEDEN) - u(i, j, k, URHO) * ke - B_ener;

        if (rho_eint > ldual_energy_eta2 * u(i, j, k, UEDEN)) {
            u(i, j, k, UEINT) = rho_eint;
        }
    });
}

void Castro::reset_internal_energy(
#ifdef MHD
    MultiFab& Bx, MultiFab& By, MultiFab& Bz,
#endif
    MultiFab& S_new, int ng)

{

    BL_PROFILE("Castro::reset_internal_energy()");

    MultiFab old_state;

    // Make a copy of the state so we can evaluate how much changed.

    if (print_update_diagnostics) {
        old_state.define(S_new.boxArray(), S_new.DistributionMap(), S_new.nComp(), 0);
        MultiFab::Copy(old_state, S_new, 0, 0, S_new.nComp(), 0);
    }

    // Ensure (rho e) isn't too small or negative
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.growntilebox(ng);

        reset_internal_energy(bx,
#ifdef MHD
                              Bx.array(mfi), By.array(mfi), Bz.array(mfi),
#endif
                              S_new.array(mfi));
    }

    if (print_update_diagnostics) {
        // Evaluate what the effective reset source was.

        MultiFab reset_source(S_new.boxArray(), S_new.DistributionMap(), S_new.nComp(), 0);

        MultiFab::Copy(reset_source, S_new, 0, 0, S_new.nComp(), 0);

        MultiFab::Subtract(reset_source, old_state, 0, 0, old_state.nComp(), 0);

        bool local = true;
        Vector<Real> reset_update = evaluate_source_change(reset_source, 1.0, local);

#ifdef BL_LAZY
        Lazy::QueueReduction([=]() mutable {
#endif
            ParallelDescriptor::ReduceRealSum(reset_update.dataPtr(), reset_update.size(),
                                              ParallelDescriptor::IOProcessorNumber());

            if (ParallelDescriptor::IOProcessor()) {
                if (std::abs(reset_update[UEINT]) != 0.0) {
                    std::cout << std::endl
                              << "  Contributions to the state from negative energy resets:"
                              << std::endl;

                    print_source_change(reset_update);
                }
            }

#ifdef BL_LAZY
        });
#endif
    }
}

void Castro::computeTemp(
#ifdef MHD
    MultiFab& Bx, MultiFab& By, MultiFab& Bz,
#endif

    MultiFab& State, Real time, int ng)

{

    BL_PROFILE("Castro::computeTemp()");

    MultiFab Stemp;

    // for 4th order, the only variables that may change here are Temp
    // and Eint.  To ensure that we don't modify values unless there is
    // a reset for Eint, we want to store the Laplacian that we use for
    // the cell-average -> cell-center conversion so we can use the same
    // Laplacian to convert back, resulting only in roundoff changes if
    // Eint is not modified.  We have to store it here, since we
    // overwrite the grown state as we work.
    MultiFab Eint_lap;

#ifdef TRUE_SDC
    if (sdc_order == 4) {

        // we need to make the data live at cell-centers first

        // fill Stemp with S_new.
        // we only need 2 ghost cells here, then the make_cell_center
        // makes 1 ghost cell a valid center, we compute its temp, and
        // then the final average results only in interior temps valid
        Stemp.define(State.boxArray(), State.DistributionMap(), NUM_STATE, 2);
        expand_state(Stemp, time, Stemp.nGrow());

        // store the Laplacian term for the internal energy
        Eint_lap.define(State.boxArray(), State.DistributionMap(), 1, 0);

        // convert to cell centers -- this will result in Stemp being
        // cell centered only on 1 ghost cells
        auto domain_lo = geom.Domain().loVect3d();
        auto domain_hi = geom.Domain().hiVect3d();

        FArrayBox tmp;

        for (MFIter mfi(Stemp); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.growntilebox(1);
            const Box& bx0 = mfi.tilebox();
            const int idx = mfi.tileIndex();

            compute_lap_term(bx0, Stemp.array(mfi), Eint_lap.array(mfi), UEINT, domain_lo,
                             domain_hi);

            tmp.resize(bx, 1);
            Elixir elix_tmp = tmp.elixir();
            auto tmp_arr = tmp.array();

            make_cell_center_in_place(bx, Stemp.array(mfi), tmp_arr, domain_lo, domain_hi);
        }
    }
#endif

#ifdef TRUE_SDC
    if (sdc_order == 4) {
        // we need to enforce minimum density here, since the conversion
        // from cell-average to centers could have made rho < 0 near steep
        // gradients
        enforce_min_density(Stemp, Stemp.nGrow());
        reset_internal_energy(Stemp, Stemp.nGrow());
    } else {
#endif
        reset_internal_energy(
#ifdef MHD
            Bx, By, Bz,
#endif
            State, ng);
#ifdef TRUE_SDC
    }
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(State, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        int num_ghost = ng;
#ifdef TRUE_SDC
        if (sdc_order == 4) {
            // only one ghost cell is at cell-centers
            num_ghost = 1;
        }
#endif

        const Box& bx = mfi.growntilebox(num_ghost);

#ifdef TRUE_SDC
        FArrayBox& u_fab = (sdc_order == 4) ? Stemp[mfi] : State[mfi];
#else
        FArrayBox& u_fab = State[mfi];
#endif

        Array4<Real> const u = u_fab.array();

        amrex::ParallelFor(bx, [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept {
            Real rhoInv = 1.0_rt / u(i, j, k, URHO);

            eos_t eos_state;

            eos_state.rho = u(i, j, k, URHO);
            eos_state.T = u(i, j, k, UTEMP);  // Initial guess for the EOS
            eos_state.e = u(i, j, k, UEINT) * rhoInv;
            for (int n = 0; n < NumSpec; ++n) eos_state.xn[n] = u(i, j, k, UFS + n) * rhoInv;
            for (int n = 0; n < NumAux; ++n) eos_state.aux[n] = u(i, j, k, UFX + n) * rhoInv;

            eos(eos_input_re, eos_state);

            u(i, j, k, UTEMP) = eos_state.T;
        });

        if (clamp_ambient_temp == 1) {
#pragma gpu box(bx)
            ca_clamp_temp(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                          BL_TO_FORTRAN_ANYD(u_fab));
        }
    }

#ifdef TRUE_SDC
    if (sdc_order == 4) {

        // we need to copy back from Stemp into S_new, making it
        // cell-average in the process.  For temperature, we will
        // construct the Laplacian from the new state and use for the
        // correction, since the temperature was just updated.  For the
        // internal energy, most zones will not have been reset, so we
        // don't want to modify them from what they were before, therefore
        // we use the stored Laplacian computed above to convert back to
        // cell-averages -- this is 4th-order and will be a no-op for
        // those zones where e wasn't changed.

        auto domain_lo = geom.Domain().loVect3d();
        auto domain_hi = geom.Domain().hiVect3d();

        FArrayBox tmp;

        for (MFIter mfi(Stemp); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.tilebox();
            const int idx = mfi.tileIndex();

            tmp.resize(bx, 1);
            Elixir elix_tmp = tmp.elixir();
            auto tmp_arr = tmp.array();

            // only temperature
            make_fourth_in_place_n(bx, Stemp.array(mfi), UTEMP, tmp_arr, domain_lo, domain_hi);
        }

        // correct UEINT
        MultiFab::Add(Stemp, Eint_lap, 0, UEINT, 1, 0);

        // copy back UTEMP and UEINT -- those are the only things that
        // should have changed.
        MultiFab::Copy(State, Stemp, UTEMP, UTEMP, 1, 0);
        MultiFab::Copy(State, Stemp, UEINT, UEINT, 1, 0);

        // now that we redid these, redo the ghost fill -- technically,
        // only need this for UTEMP and UEINT, and only if ng > 0
        if (ng > 0) {
            AmrLevel::FillPatch(*this, State, State.nGrow(), time, State_Type, 0, NUM_STATE);
        }

        Stemp.clear();
    }
#endif
}

void Castro::create_source_corrector() {

    BL_PROFILE("Castro::create_source_corrector()");

    if (time_integration_method == CornerTransportUpwind && source_term_predictor == 1) {

        // Optionally predict the source terms to t + dt/2,
        // which is the time-level n+1/2 value, To do this we use a
        // lagged predictor estimate: dS/dt_n = (S_n - S_{n-1}) / dt, so
        // S_{n+1/2} = S_n + (dt / 2) * dS/dt_n. We'll add the S_n
        // terms later; now we add the second term. We defer the
        // multiplication by dt / 2 until the actual advance, since
        // we may be subcycling and thus not know yet what the
        // advance timestep is.

        // Note that since for dS/dt we want (S^{n+1} - S^{n}) / dt,
        // we only need to take twice the new-time source term from the
        // last timestep, since in the predictor-corrector approach,
        // the new-time source term is 1/2 * S^{n+1} - 1/2 * S^{n}.
        // This is untrue in general for the non-momentum sources,
        // so for safety we'll only apply it to the momentum sources.

        // Even though we're calculating this predictor from the last
        // timestep, we've already done the swap, so the "new" data
        // from the last timestep is currently residing in the "old"
        // StateData. (As a corollary, this operation must be done
        // prior to updating any of the source StateData.) Since the
        // dt comes from the last timestep, which is no longer equal
        // to the difference between prevTime and curTime, we rely
        // on our recording of the last dt from the previous advance.

        const Real time = get_state_data(Source_Type).prevTime();

        AmrLevel::FillPatch(*this, source_corrector, NUM_GROW, time, Source_Type, UMX, 3, UMX);

        source_corrector.mult(2.0 / lastDt, NUM_GROW);

    } else if (time_integration_method == SimplifiedSpectralDeferredCorrections) {

        // If we're doing simplified SDC, time-center the source term (using the
        // current iteration's old sources and the last iteration's new
        // sources). Since the "new-time" sources are just the corrector step
        // of the predictor-corrector formalism, we want to add the full
        // value of the "new-time" sources to the old-time sources to get a
        // time-centered value. Note that, as above, the "new" data from the
        // last step is currently residing in the "old" StateData since we
        // have already done the swap.

        const Real time = get_state_data(Source_Type).prevTime();

        AmrLevel::FillPatch(*this, source_corrector, NUM_GROW, time, Source_Type, 0, NSRC);
    }
}

void Castro::swap_state_time_levels(const Real dt) {

    BL_PROFILE("Castro::swap_state_time_levels()");

    MultiFab::RegionTag statedata_tag("StateData_Level_" + std::to_string(level));
    MultiFab::RegionTag amrlevel_tag("AmrLevel_Level_" + std::to_string(level));

    for (int k = 0; k < num_state_type; k++) {

        // The following is a hack to make sure that we only
        // ever have new data for certain state types that only
        // ever need new time data; by doing a swap now, we'll
        // guarantee that allocOldData() does nothing. We do
        // this because we never need the old data, so we
        // don't want to allocate memory for it.

#ifdef SIMPLIFIED_SDC
#ifdef REACTIONS
        if (time_integration_method == SimplifiedSpectralDeferredCorrections) {
            if (k == Simplified_SDC_React_Type) {
                state[k].swapTimeLevels(0.0);
            }
        }
#endif
#endif

#ifdef TRUE_SDC
#ifdef REACTIONS
        if (time_integration_method == SpectralDeferredCorrections && sdc_order == 4 &&
            k == SDC_Source_Type) {
            state[k].swapTimeLevels(0.0);
        }
#endif
#endif
        state[k].allocOldData();

        state[k].swapTimeLevels(dt);
    }
}

#ifdef GRAVITY
int Castro::get_numpts() {
    int numpts_1d;

    Box bx(geom.Domain());
    long nx = bx.size()[0];

#if (BL_SPACEDIM == 1)
    numpts_1d = nx;
#elif (BL_SPACEDIM == 2)
    long ny = bx.size()[1];
    Real ndiagsq = Real(nx * nx + ny * ny);
    numpts_1d = int(sqrt(ndiagsq)) + 2 * NUM_GROW;
#elif (BL_SPACEDIM == 3)
    long ny = bx.size()[1];
    long nz = bx.size()[2];
    Real ndiagsq = Real(nx * nx + ny * ny + nz * nz);
    numpts_1d = int(sqrt(ndiagsq)) + 2 * NUM_GROW;
#endif

    if (verbose && ParallelDescriptor::IOProcessor()) {
        std::cout << "Castro::numpts_1d at level  " << level << " is " << numpts_1d << std::endl;
    }

    return numpts_1d;
}

void Castro::make_radial_data(int is_new) {
#if (BL_SPACEDIM > 1)

    BL_PROFILE("Castro::make_radial_data()");

    // We only call this for level = 0
    BL_ASSERT(level == 0);

    int numpts_1d = get_numpts();

    Vector<Real> radial_vol(numpts_1d, 0);

    const Real* dx = geom.CellSize();
    Real dr = dx[0];

    if (is_new == 1) {
        MultiFab& S = get_new_data(State_Type);
        const int nc = S.nComp();
        Vector<Real> radial_state(numpts_1d * nc, 0);
        for (MFIter mfi(S); mfi.isValid(); ++mfi) {
            Box bx(mfi.validbox());
            ca_compute_avgstate(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()), ZFILL(dx), &dr, &nc,
                                BL_TO_FORTRAN_ANYD(S[mfi]), radial_state.dataPtr(),
                                BL_TO_FORTRAN_ANYD(volume[mfi]), radial_vol.dataPtr(),
                                ZFILL(geom.ProbLo()), &numpts_1d);
        }

        ParallelDescriptor::ReduceRealSum(radial_vol.dataPtr(), numpts_1d);
        ParallelDescriptor::ReduceRealSum(radial_state.dataPtr(), numpts_1d * nc);

        int first = 0;
        int np_max = 0;
        for (int i = 0; i < numpts_1d; i++) {
            if (radial_vol[i] > 0.) {
                for (int j = 0; j < nc; j++) {
                    radial_state[nc * i + j] /= radial_vol[i];
                }
            } else if (first == 0) {
                np_max = i;
                first = 1;
            }
        }

        Vector<Real> radial_state_short(np_max * nc, 0);

        for (int i = 0; i < np_max; i++) {
            for (int j = 0; j < nc; j++) {
                radial_state_short[nc * i + j] = radial_state[nc * i + j];
            }
        }

        const Real new_time = state[State_Type].curTime();
        set_new_outflow_data(radial_state_short.dataPtr(), &new_time, &np_max, &nc);
    } else {
        MultiFab& S = get_old_data(State_Type);
        const int nc = S.nComp();
        Vector<Real> radial_state(numpts_1d * nc, 0);
        for (MFIter mfi(S); mfi.isValid(); ++mfi) {
            Box bx(mfi.validbox());
            ca_compute_avgstate(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()), ZFILL(dx), &dr, &nc,
                                BL_TO_FORTRAN_ANYD(S[mfi]), radial_state.dataPtr(),
                                BL_TO_FORTRAN_ANYD(volume[mfi]), radial_vol.dataPtr(),
                                ZFILL(geom.ProbLo()), &numpts_1d);
        }

        ParallelDescriptor::ReduceRealSum(radial_vol.dataPtr(), numpts_1d);
        ParallelDescriptor::ReduceRealSum(radial_state.dataPtr(), numpts_1d * nc);

        int first = 0;
        int np_max = 0;
        for (int i = 0; i < numpts_1d; i++) {
            if (radial_vol[i] > 0.) {
                for (int j = 0; j < nc; j++) {
                    radial_state[nc * i + j] /= radial_vol[i];
                }
            } else if (first == 0) {
                np_max = i;
                first = 1;
            }
        }

        Vector<Real> radial_state_short(np_max * nc, 0);

        for (int i = 0; i < np_max; i++) {
            for (int j = 0; j < nc; j++) {
                radial_state_short[nc * i + j] = radial_state[nc * i + j];
            }
        }

        const Real old_time = state[State_Type].prevTime();
        set_old_outflow_data(radial_state_short.dataPtr(), &old_time, &np_max, &nc);
    }

#endif
}

void Castro::define_new_center(MultiFab& S, Real time) {
    BL_PROFILE("Castro::define_new_center()");

    Real center[3];
    const Real* dx = geom.CellSize();

    IntVect max_index = S.maxIndex(URHO, 0);
    Box bx(max_index, max_index);
    bx.grow(1);
    BoxArray ba(bx);
    int owner = ParallelDescriptor::IOProcessorNumber();
    DistributionMapping dm{Vector<int>(1, owner)};
    MultiFab mf(ba, dm, 1, 0);

    // Define a cube 3-on-a-side around the point with the maximum density
    FillPatch(*this, mf, 0, time, State_Type, URHO, 1);

    int mi[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) {
        mi[i] = max_index[i];
    }

    // Find the position of the "center" by interpolating from data at cell centers
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        ca_find_center(mf[mfi].dataPtr(), &center[0], ARLIM_3D(mi), ZFILL(dx),
                       ZFILL(geom.ProbLo()));
    }
    // Now broadcast to everyone else.
    ParallelDescriptor::Bcast(&center[0], BL_SPACEDIM, owner);

    // Make sure if R-Z that center stays exactly on axis
    if (Geom().IsRZ()) {
        center[0] = 0;
    }

    ca_set_center(ZFILL(center));
}

void Castro::write_center() {
    int ndatalogs = parent->NumDataLogs();

    if ((moving_center == 1) && (ndatalogs > 0) && ParallelDescriptor::IOProcessor()) {
        std::ostream& data_logc = parent->DataLog(0);

        int nstep = parent->levelSteps(0);
        Real time = state[State_Type].curTime();

        Real center[3];
        ca_get_center(center);

        if (time == 0.0) {
            data_logc << std::setw(8) << "   nstep";
            data_logc << std::setw(14) << "         time  ";
            data_logc << std::setw(14) << "         center" << std::endl;
            ;
        }

        data_logc << std::setw(8) << nstep;
        data_logc << std::setw(14) << std::setprecision(6) << time;
        data_logc << std::setw(14) << std::setprecision(6) << center[0];
#if (BL_SPACEDIM >= 2)
        data_logc << std::setw(14) << std::setprecision(6) << center[1];
#endif
#if (BL_SPACEDIM == 3)
        data_logc << std::setw(14) << std::setprecision(6) << center[2];
#endif
        data_logc << std::endl;
    }
}
#endif

Real Castro::getCPUTime() {

    int numCores = ParallelDescriptor::NProcs();
#ifdef _OPENMP
    numCores = numCores * omp_get_max_threads();
#endif

    Real T = numCores * (ParallelDescriptor::second() - startCPUTime) + previousCPUTimeUsed;

    return T;
}

MultiFab& Castro::build_fine_mask() {
    BL_PROFILE("Castro::build_fine_mask()");

    BL_ASSERT(level > 0);  // because we are building a mask for the coarser level

    if (fine_mask.empty()) {
        fine_mask = makeFineMask(parent->boxArray(level - 1), parent->DistributionMap(level - 1),
                                 parent->boxArray(level), crse_ratio,
                                 1.0,   // coarse
                                 0.0);  // fine
    }

    return fine_mask;
}

iMultiFab& Castro::build_interior_boundary_mask(int ng) {
    BL_PROFILE("Castro::build_interior_boundary_mask()");

    for (int i = 0; i < ib_mask.size(); ++i) {
        if (ib_mask[i]->nGrow() == ng) {
            return *ib_mask[i];
        }
    }

    //  If we got here, we need to build a new one
    ib_mask.push_back(std::unique_ptr<iMultiFab>(new iMultiFab(grids, dmap, 1, ng)));

    iMultiFab& imf = *ib_mask.back();

    int ghost_covered_by_valid = 0;
    int other_cells = 1;  // uncovered ghost, valid, and outside domain cells are set to 1

    imf.BuildMask(geom.Domain(), geom.periodicity(), ghost_covered_by_valid, other_cells,
                  other_cells, other_cells);

    return imf;
}

// Fill a version of the state with ng ghost zones from the state data.

void Castro::expand_state(MultiFab& S, Real time, int ng) {
    BL_PROFILE("Castro::expand_state()");

    BL_ASSERT(S.nGrow() >= ng);

    AmrLevel::FillPatch(*this, S, ng, time, State_Type, 0, NUM_STATE);
}

void Castro::check_for_nan(MultiFab& state_in, int check_ghost) {
    BL_PROFILE("Castro::check_for_nan()");

    int ng = 0;
    if (check_ghost == 1) {
        ng = state_in.nGrow();
    }

    if (state_in.contains_nan(URHO, state_in.nComp(), ng, true)) {
        for (int i = 0; i < state_in.nComp(); i++) {
            if (state_in.contains_nan(URHO + i, 1, ng, true)) {
                std::string abort_string = std::string("State has NaNs in the ") +
                                           desc_lst[State_Type].name(i) +
                                           std::string(" component::check_for_nan()");
                amrex::Abort(abort_string.c_str());
            }
        }
    }
}

// Given State_Type state data, perform a number of cleaning steps to make
// sure the data is sensible.

void Castro::clean_state(
#ifdef MHD
    MultiFab& bx, MultiFab& by, MultiFab& bz,
#endif
    MultiFab& state_in, Real time, int ng) {

    BL_PROFILE("Castro::clean_state()");

    // Enforce a minimum density.

    enforce_min_density(state_in, ng);

    // Ensure all species are normalized.

    normalize_species(state_in, ng);

    // Sync the linear and hybrid momenta.

#ifdef HYBRID_MOMENTUM
    if (hybrid_hydro) {
        hybrid_to_linear_momentum(state_in, ng);
    }
#endif

    // Compute the temperature (note that this will also reset
    // the internal energy for consistency with the total energy).

    computeTemp(
#ifdef MHD
        bx, by, bz,
#endif
        state_in, time, ng);
}
