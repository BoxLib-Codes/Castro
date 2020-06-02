
#ifndef WIN32
#include <unistd.h>
#endif

#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>

#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include "Castro.H"
#include "Castro_F.H"
#include "Castro_io.H"

#ifdef RADIATION
#include "Radiation.H"
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

#include "AMReX_buildInfo.H"

using std::string;
using namespace amrex;

// Castro maintains an internal checkpoint version numbering system.
// This allows us to ensure that we don't attempt to restart from a
// checkpoint that is incompatible with the current code. The version
// number is stored in the CastroHeader file inside a checkpoint.
// The following were the changes made in updating version numbers:
// 0: all checkpoints before we began the numbering system (July 27, 2015)
// 1: PhiGrav_Type was added to the checkpoint
// 2: Source_Type was added to the checkpoint
// 3: A ReactHeader file was generated and the maximum de/dt was stored there
// 4: Reactions_Type added to checkpoint; ReactHeader functionality deprecated
// 5: Simplified_SDC_Source_Type and Simplified_SDC_React_Type added to checkpoint
// 6: Simplified_SDC_Source_Type removed from Castro
// 7: A weights field was added to Reactions_Type; number of ghost zones increased to NUM_GROW

namespace {
    int input_version = -1;
    int current_version = 7;
}  // namespace

// I/O routines for Castro

void Castro::restart(Amr& papa, istream& is, bool bReadSpecial) {
    // Let's check Castro checkpoint version first;
    // trying to read from checkpoint; if nonexisting, set it to 0.
    if (input_version == -1) {
        if (ParallelDescriptor::IOProcessor()) {
            std::ifstream CastroHeaderFile;
            std::string FullPathCastroHeaderFile = papa.theRestartFile();
            FullPathCastroHeaderFile += "/CastroHeader";
            CastroHeaderFile.open(FullPathCastroHeaderFile.c_str(), std::ios::in);
            if (CastroHeaderFile.good()) {
                char foo[256];
                // first line: Checkpoint version: ?
                CastroHeaderFile.getline(foo, 256, ':');
                CastroHeaderFile >> input_version;
                CastroHeaderFile.close();
            } else {
                input_version = 0;
            }
        }
        ParallelDescriptor::Bcast(&input_version, 1, ParallelDescriptor::IOProcessorNumber());
    }

    // Check that the version number matches.

    if (input_version != current_version) {
        amrex::Error("Checkpoint format incompatible with current code");
    }

    // Check if there's a file in the header indicating that the
    // previous timestep was limited to hit a plot interval. If so,
    // read in the value, so that the timestep after this restart is
    // limited appropriately.

    if (ParallelDescriptor::IOProcessor()) {

        std::ifstream dtHeaderFile;
        std::string FullPathdtHeaderFile = papa.theRestartFile();
        FullPathdtHeaderFile += "/dtHeader";
        dtHeaderFile.open(FullPathdtHeaderFile.c_str(), std::ios::in);

        if (dtHeaderFile.good()) {

            lastDtPlotLimited = 1;
            dtHeaderFile >> lastDtBeforePlotLimiting;
            dtHeaderFile.close();
        }
    }

    // Check if we have the same state variables

    if (ParallelDescriptor::IOProcessor()) {

        std::ifstream StateListFile;
        std::string FullPathStateList = papa.theRestartFile();
        FullPathStateList += "/state_names.txt";
        StateListFile.open(FullPathStateList.c_str(), std::ios::in);

        if (StateListFile.good()) {
            std::string var;
            for (int n = 0; n < NUM_STATE; n++) {
                StateListFile >> var;
                if (var != desc_lst[State_Type].name(n)) {
                    amrex::Error("state variables do not agree");
                }
            }
            StateListFile.close();
        }
    }

    ParallelDescriptor::Bcast(&lastDtPlotLimited, 1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&lastDtBeforePlotLimiting, 1,
                              ParallelDescriptor::IOProcessorNumber());

    // also need to mod checkPoint function to store the new version in a text file

    AmrLevel::restart(papa, is, bReadSpecial);

    buildMetrics();

    initMFs();

    // get the elapsed CPU time to now;
    if (level == 0 && ParallelDescriptor::IOProcessor()) {
        // get ellapsed CPU time
        std::ifstream CPUFile;
        std::string FullPathCPUFile = parent->theRestartFile();
        FullPathCPUFile += "/CPUtime";
        CPUFile.open(FullPathCPUFile.c_str(), std::ios::in);

        if (CPUFile.good()) {

            CPUFile >> previousCPUTimeUsed;
            CPUFile.close();
        }

        std::cout << "read CPU time: " << previousCPUTimeUsed << "\n";
    }

    if (track_grid_losses && level == 0) {

        // get the current value of the diagnostic quantities
        std::ifstream DiagFile;
        std::string FullPathDiagFile = parent->theRestartFile();
        FullPathDiagFile += "/Diagnostics";
        DiagFile.open(FullPathDiagFile.c_str(), std::ios::in);

        if (DiagFile.good()) {

            for (int i = 0; i < n_lost; i++) {
                DiagFile >> material_lost_through_boundary_cumulative[i];
                material_lost_through_boundary_temp[i] = 0.0;
            }

            DiagFile.close();
        }
    }

#ifdef GRAVITY
    if (use_point_mass && level == 0) {

        // get the current value of the point mass
        std::ifstream PMFile;
        std::string FullPathPMFile = parent->theRestartFile();
        FullPathPMFile += "/point_mass";
        PMFile.open(FullPathPMFile.c_str(), std::ios::in);

        if (PMFile.good()) {
            PMFile >> point_mass;
            set_pointmass(&point_mass);
            PMFile.close();
        }
    }
#endif

    if (level == 0) {
        // get problem-specific stuff -- note all processors do this,
        // eliminating the need for a broadcast
        std::string dir = parent->theRestartFile();

        char* dir_for_pass = new char[dir.size() + 1];
        std::copy(dir.begin(), dir.end(), dir_for_pass);
        dir_for_pass[dir.size()] = '\0';

        int len = dir.size();

        Vector<int> int_dir_name(len);
        for (int j = 0; j < len; j++) {
            int_dir_name[j] = (int)dir_for_pass[j];
        }

        problem_restart(int_dir_name.dataPtr(), &len);

        delete[] dir_for_pass;
    }

    const Real* dx = geom.CellSize();

    if ((grown_factor > 1) && (parent->maxLevel() < 1)) {
        std::cout << "grown_factor is " << grown_factor << std::endl;
        std::cout << "max_level is " << parent->maxLevel() << std::endl;
        amrex::Error("Must have max_level > 0 if doing special restart with grown_factor");
    }

    if (grown_factor > 1 && level == 0) {
        if (verbose && ParallelDescriptor::IOProcessor()) {
            std::cout << "Doing special restart with grown_factor " << grown_factor << std::endl;
        }
        MultiFab& S_new = get_new_data(State_Type);

        Box orig_domain;
        if (star_at_center == 0) {
            orig_domain = amrex::coarsen(geom.Domain(), grown_factor);
        } else if (star_at_center == 1) {

            Box domain(geom.Domain());
            int lo = 0, hi = 0;
            if (geom.IsRZ()) {
                if (grown_factor != 2) {
                    amrex::Abort("Must have grown_factor = 2");
                }

                int d = 0;
                int dlen = domain.size()[d];
                lo = 0;
                hi = dlen / 2;
                orig_domain.setSmall(d, lo);
                orig_domain.setBig(d, hi);

                d = 1;
                dlen = domain.size()[d];
                lo = dlen / 4;
                hi = 3 * dlen / 4 - 1;
                orig_domain.setSmall(d, lo);
                orig_domain.setBig(d, hi);

            } else {
                for (int d = 0; d < BL_SPACEDIM; d++) {
                    int dlen = domain.size()[d];
                    if (grown_factor == 2) {
                        lo = dlen / 4;
                        hi = 3 * dlen / 4 - 1;
                    } else if (grown_factor == 3) {
                        lo = (dlen) / 3;
                        hi = 2 * (dlen) / 3 - 1;
                    } else {
                        amrex::Abort("Must have grown_factor = 2 or 3");
                    }
                    orig_domain.setSmall(d, lo);
                    orig_domain.setBig(d, hi);
                }
            }
        } else {
            if (ParallelDescriptor::IOProcessor()) {
                std::cout << "... invalid value of star_at_center: " << star_at_center << std::endl;
            }
            amrex::Abort();
        }

        for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {

            const Real* prob_lo = geom.ProbLo();
            const Box& bx = mfi.validbox();
            const int* lo = bx.loVect();
            const int* hi = bx.hiVect();

            if (!orig_domain.contains(bx)) {

#ifdef GPU_COMPATIBLE_PROBLEM

#pragma gpu box(bx)
                ca_initdata(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi), BL_TO_FORTRAN_ANYD(S_new[mfi]),
                            AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(prob_lo));

#else

                RealBox gridloc = RealBox(grids[mfi.index()], geom.CellSize(), geom.ProbLo());

                Real cur_time = state[State_Type].curTime();

                BL_FORT_PROC_CALL(CA_INITDATA, ca_initdata)
                (level, cur_time, ARLIM_3D(lo), ARLIM_3D(hi), NUM_STATE,
                 BL_TO_FORTRAN_ANYD(S_new[mfi]), ZFILL(dx), ZFILL(gridloc.lo()),
                 ZFILL(gridloc.hi()));

#endif
            }
        }
    }

    if (grown_factor > 1 && level == 1) {
        getLevel(0).avgDown();
    }

#ifdef GRAVITY
#if (BL_SPACEDIM > 1)
    if ((level == 0) && (spherical_star == 1)) {
        MultiFab& S_new = get_new_data(State_Type);
        const int nc = S_new.nComp();
        const int n1d = get_numpts();
        allocate_outflow_data(&n1d, &nc);
        int is_new = 1;
        make_radial_data(is_new);
    }
#endif

    if (do_grav && level == 0) {
        BL_ASSERT(gravity == 0);
        gravity = new Gravity(parent, parent->finestLevel(), &phys_bc, URHO);
    }
#endif

#ifdef DIFFUSION
    if (level == 0) {
        BL_ASSERT(diffusion == 0);
        diffusion = new Diffusion(parent, &phys_bc);
    }
#endif

#ifdef RADIATION
    if (do_radiation) {
        if (radiation == 0) {
            // radiation is a static object, only alloc if not already there
            int rad_restart = 1;  // disables quasi-steady initialization
            radiation = new Radiation(parent, this, rad_restart);
        }
        radiation->regrid(level, grids, dmap);
        radiation->restart(level, grids, dmap, parent->theRestartFile(), is);

        rad_solver.reset(new RadSolve(parent, level, grids, dmap));
    }
#endif

    // If we want, we can restart the checkpoint at a new time.

    if (reset_checkpoint_time > -1.e199) {

        if (!parent->RegridOnRestart()) {
            amrex::Error("reset_checkpoint_time only makes sense when amr.regrid_on_restart=1");
        }

        const Real cur_time = get_state_data(State_Type).curTime();
        const Real prev_time = get_state_data(State_Type).prevTime();
        const Real dt = cur_time - prev_time;

        parent->setStartTime(reset_checkpoint_time);
        parent->setCumTime(reset_checkpoint_time);

        for (int n = 0; n < num_state_type; ++n) {
            StateData& cur_state = get_state_data(n);
            cur_state.setOldTimeLevel(reset_checkpoint_time - dt);
            cur_state.setNewTimeLevel(reset_checkpoint_time);
        }
    }

    if (reset_checkpoint_step > -1) {

        if (!parent->RegridOnRestart()) {
            amrex::Error("reset_checkpoint_step only makes sense when amr.regrid_on_restart=1");
        }

        parent->setLevelSteps(level, reset_checkpoint_step);
        parent->setLevelCount(level, reset_checkpoint_step);
    }
}

void Castro::set_state_in_checkpoint(Vector<int>& state_in_checkpoint) {
    for (int i = 0; i < num_state_type; ++i) {
        state_in_checkpoint[i] = 1;
    }
}

void Castro::checkPoint(const std::string& dir, std::ostream& os, VisMF::How how,
                        bool dump_old_default) {

    const Real io_start_time = ParallelDescriptor::second();

    AmrLevel::checkPoint(dir, os, how, dump_old);

    const Real io_time = ParallelDescriptor::second() - io_start_time;

#ifdef RADIATION
    if (do_radiation) {
        radiation->checkPoint(level, dir, os, how);
    }
#endif

#ifdef AMREX_PARTICLES
    ParticleCheckPoint(dir);
#endif

    if (level == 0 && ParallelDescriptor::IOProcessor()) {
        {
            std::ofstream CastroHeaderFile;
            std::string FullPathCastroHeaderFile = dir;
            FullPathCastroHeaderFile += "/CastroHeader";
            CastroHeaderFile.open(FullPathCastroHeaderFile.c_str(), std::ios::out);

            CastroHeaderFile << "Checkpoint version: " << current_version << std::endl;
            CastroHeaderFile.close();

            writeJobInfo(dir, io_time);

            // output the list of state variables, so we can do a sanity check on restart
            std::ofstream StateListFile;
            std::string FullPathStateList = dir;
            FullPathStateList += "/state_names.txt";
            StateListFile.open(FullPathStateList.c_str(), std::ios::out);

            for (int n = 0; n < NUM_STATE; n++) {
                StateListFile << desc_lst[State_Type].name(n) << "\n";
            }
            StateListFile.close();
        }

        // If we have limited this last timestep to hit a plot interval,
        // store the timestep we took before limiting. After the restart,
        // this dt will be read in and used to limit the next timestep
        // appropriately, rather than the shortened timestep.

        if (lastDtPlotLimited == 1) {

            std::ofstream dtHeaderFile;
            std::string FullPathdtHeaderFile = dir;
            FullPathdtHeaderFile += "/dtHeader";
            dtHeaderFile.open(FullPathdtHeaderFile.c_str(), std::ios::out);

            dtHeaderFile << lastDtBeforePlotLimiting << std::endl;
            dtHeaderFile.close();
        }

        {
            // store elapsed CPU time
            std::ofstream CPUFile;
            std::string FullPathCPUFile = dir;
            FullPathCPUFile += "/CPUtime";
            CPUFile.open(FullPathCPUFile.c_str(), std::ios::out);

            CPUFile << std::setprecision(17) << getCPUTime();
            CPUFile.close();
        }

        if (track_grid_losses) {

            // store diagnostic quantities
            std::ofstream DiagFile;
            std::string FullPathDiagFile = dir;
            FullPathDiagFile += "/Diagnostics";
            DiagFile.open(FullPathDiagFile.c_str(), std::ios::out);

            for (int i = 0; i < n_lost; i++) {
                DiagFile << std::setprecision(17) << material_lost_through_boundary_cumulative[i]
                         << std::endl;
            }

            DiagFile.close();
        }

#ifdef GRAVITY
        if (use_point_mass) {

            // store current value of the point mass
            std::ofstream PMFile;
            std::string FullPathPMFile = dir;
            FullPathPMFile += "/point_mass";
            PMFile.open(FullPathPMFile.c_str(), std::ios::out);

            PMFile << std::setprecision(17) << point_mass << std::endl;

            PMFile.close();
        }
#endif

        {
            // store any problem-specific stuff
            char* dir_for_pass = new char[dir.size() + 1];
            std::copy(dir.begin(), dir.end(), dir_for_pass);
            dir_for_pass[dir.size()] = '\0';

            int len = dir.size();

            Vector<int> int_dir_name(len);
            for (int j = 0; j < len; j++) {
                int_dir_name[j] = (int)dir_for_pass[j];
            }

            problem_checkpoint(int_dir_name.dataPtr(), &len);

            delete[] dir_for_pass;
        }
    }
}

std::string Castro::thePlotFileType() const {
    //
    // Increment this whenever the writePlotFile() format changes.
    //
    static const std::string the_plot_file_type("HyperCLaw-V1.1");

    return the_plot_file_type;
}

void Castro::setPlotVariables() {
    AmrLevel::setPlotVariables();

    // Don't add the Source_Type data to the plotfile, we only
    // want to store it in the checkpoints.

    for (int i = 0; i < desc_lst[Source_Type].nComp(); i++) {
        parent->deleteStatePlotVar(desc_lst[Source_Type].name(i));
    }

#ifdef SIMPLIFIED_SDC
#ifdef REACTIONS
    if (time_integration_method == SimplifiedSpectralDeferredCorrections) {
        for (int i = 0; i < desc_lst[Simplified_SDC_React_Type].nComp(); i++) {
            parent->deleteStatePlotVar(desc_lst[Simplified_SDC_React_Type].name(i));
        }
    }
#endif
#endif

    ParmParse pp("castro");

    bool plot_X;

    if (pp.query("plot_X", plot_X)) {
        if (plot_X) {
            // Get the species names from the network model.
            //
            for (int i = 0; i < NumSpec; i++) {
                string spec_string = "X(" + short_spec_names_cxx[i] + ")";
                parent->addDerivePlotVar(spec_string);
            }
        }
    }
}

void Castro::writeJobInfo(const std::string& dir, const Real io_time) {

    // job_info file with details about the run
    std::ofstream jobInfoFile;
    std::string FullPathJobInfoFile = dir;
    FullPathJobInfoFile += "/job_info";
    jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

    std::string PrettyLine = std::string(78, '=') + "\n";
    std::string OtherLine = std::string(78, '-') + "\n";
    std::string SkipSpace = std::string(8, ' ');

    // job information
    jobInfoFile << PrettyLine;
    jobInfoFile << " Castro Job Information\n";
    jobInfoFile << PrettyLine;

    jobInfoFile << "job name: " << job_name << "\n\n";
    jobInfoFile << "inputs file: " << inputs_name << "\n\n";

    jobInfoFile << "number of MPI processes: " << ParallelDescriptor::NProcs() << "\n";
#ifdef _OPENMP
    jobInfoFile << "number of threads:       " << omp_get_max_threads() << "\n";
#endif
    jobInfoFile << "\n";
    jobInfoFile << "hydro tile size:         " << hydro_tile_size << "\n";

    jobInfoFile << "\n";
    jobInfoFile << "CPU time used since start of simulation (CPU-hours): " << getCPUTime() / 3600.0;

    jobInfoFile << "\n\n";

    // plotfile information
    jobInfoFile << PrettyLine;
    jobInfoFile << " Plotfile Information\n";
    jobInfoFile << PrettyLine;

    time_t now = time(0);

    // Convert now to tm struct for local timezone
    tm* localtm = localtime(&now);
    jobInfoFile << "output date / time: " << asctime(localtm);

    char currentDir[FILENAME_MAX];
    if (getcwd(currentDir, FILENAME_MAX)) {
        jobInfoFile << "output dir:         " << currentDir << "\n";
    }

    jobInfoFile << "I/O time (s):       " << io_time << "\n";

    jobInfoFile << "\n\n";

#ifdef AMREX_USE_GPU
    // This output assumes for simplicity that every rank uses the
    // same type of GPU.

    jobInfoFile << PrettyLine;
    jobInfoFile << "GPU Information:       "
                << "\n";
    jobInfoFile << PrettyLine;

    jobInfoFile << "GPU model name: " << Gpu::Device::deviceName() << "\n";
    jobInfoFile << "Number of GPUs used: " << Gpu::Device::numDevicesUsed() << "\n";

    jobInfoFile << "\n\n";
#endif

    // build information
    jobInfoFile << PrettyLine;
    jobInfoFile << " Build Information\n";
    jobInfoFile << PrettyLine;

    jobInfoFile << "build date:    " << buildInfoGetBuildDate() << "\n";
    jobInfoFile << "build machine: " << buildInfoGetBuildMachine() << "\n";
    jobInfoFile << "build dir:     " << buildInfoGetBuildDir() << "\n";
    jobInfoFile << "AMReX dir:     " << buildInfoGetAMReXDir() << "\n";

    jobInfoFile << "\n";

    jobInfoFile << "COMP:          " << buildInfoGetComp() << "\n";
    jobInfoFile << "COMP version:  " << buildInfoGetCompVersion() << "\n";

    jobInfoFile << "\n";

    jobInfoFile << "C++ compiler:  " << buildInfoGetCXXName() << "\n";
    jobInfoFile << "C++ flags:     " << buildInfoGetCXXFlags() << "\n";

    jobInfoFile << "\n";

    jobInfoFile << "Fortran comp:  " << buildInfoGetFName() << "\n";
    jobInfoFile << "Fortran flags: " << buildInfoGetFFlags() << "\n";

    jobInfoFile << "\n";

    jobInfoFile << "Link flags:    " << buildInfoGetLinkFlags() << "\n";
    jobInfoFile << "Libraries:     " << buildInfoGetLibraries() << "\n";

    jobInfoFile << "\n";

    for (int n = 1; n <= buildInfoGetNumModules(); n++) {
        jobInfoFile << buildInfoGetModuleName(n) << ": " << buildInfoGetModuleVal(n) << "\n";
    }

    jobInfoFile << "\n";

    const char* githash1 = buildInfoGetGitHash(1);
    const char* githash2 = buildInfoGetGitHash(2);
    const char* githash3 = buildInfoGetGitHash(3);
    if (strlen(githash1) > 0) {
        jobInfoFile << "Castro       git describe: " << githash1 << "\n";
    }
    if (strlen(githash2) > 0) {
        jobInfoFile << "AMReX        git describe: " << githash2 << "\n";
    }
    if (strlen(githash3) > 0) {
        jobInfoFile << "Microphysics git describe: " << githash3 << "\n";
    }

    const char* buildgithash = buildInfoGetBuildGitHash();
    const char* buildgitname = buildInfoGetBuildGitName();
    if (strlen(buildgithash) > 0) {
        jobInfoFile << buildgitname << " git describe: " << buildgithash << "\n";
    }

    jobInfoFile << "\n\n";

    // grid information
    jobInfoFile << PrettyLine;
    jobInfoFile << " Grid Information\n";
    jobInfoFile << PrettyLine;

    int f_lev = parent->finestLevel();

    for (int i = 0; i <= f_lev; i++) {
        jobInfoFile << " level: " << i << "\n";
        jobInfoFile << "   number of boxes = " << parent->numGrids(i) << "\n";
        jobInfoFile << "   maximum zones   = ";
        for (int n = 0; n < BL_SPACEDIM; n++) {
            jobInfoFile << parent->Geom(i).Domain().length(n) << " ";
            //jobInfoFile << parent->Geom(i).ProbHi(n) << " ";
        }
        jobInfoFile << "\n\n";
    }

    jobInfoFile << " Boundary conditions\n";
    Vector<int> lo_bc_out(BL_SPACEDIM), hi_bc_out(BL_SPACEDIM);
    ParmParse pp("castro");
    pp.getarr("lo_bc", lo_bc_out, 0, BL_SPACEDIM);
    pp.getarr("hi_bc", hi_bc_out, 0, BL_SPACEDIM);

    // these names correspond to the integer flags setup in the
    // Castro_setup.cpp
    const char* names_bc[] = {"interior", "inflow",   "outflow",
                              "symmetry", "slipwall", "noslipwall"};

    jobInfoFile << "   -x: " << names_bc[lo_bc_out[0]] << "\n";
    jobInfoFile << "   +x: " << names_bc[hi_bc_out[0]] << "\n";
    if (BL_SPACEDIM >= 2) {
        jobInfoFile << "   -y: " << names_bc[lo_bc_out[1]] << "\n";
        jobInfoFile << "   +y: " << names_bc[hi_bc_out[1]] << "\n";
    }
    if (BL_SPACEDIM == 3) {
        jobInfoFile << "   -z: " << names_bc[lo_bc_out[2]] << "\n";
        jobInfoFile << "   +z: " << names_bc[hi_bc_out[2]] << "\n";
    }

    jobInfoFile << "\n\n";

    jobInfoFile << " Domain geometry info\n";

    Real center[3];
    ca_get_center(center);

    jobInfoFile << "     center: " << center[0] << " , " << center[1] << " , " << center[2] << "\n";
    jobInfoFile << "\n";

    jobInfoFile << "     geometry.is_periodic: ";
    for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
        jobInfoFile << geom.isPeriodic(idir) << " ";
    }
    jobInfoFile << "\n";

    jobInfoFile << "     geometry.coord_sys:   " << geom.Coord() << "\n";

    jobInfoFile << "     geometry.prob_lo:     ";
    for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
        jobInfoFile << geom.ProbLo(idir) << " ";
    }
    jobInfoFile << "\n";

    jobInfoFile << "     geometry.prob_hi:     ";
    for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
        jobInfoFile << geom.ProbHi(idir) << " ";
    }
    jobInfoFile << "\n";

    jobInfoFile << "     amr.n_cell:           ";
    const int* domain_lo = geom.Domain().loVect();
    const int* domain_hi = geom.Domain().hiVect();
    for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
        jobInfoFile << domain_hi[idir] - domain_lo[idir] + 1 << " ";
    }
    jobInfoFile << "\n";

    int max_level = parent->maxLevel();
    jobInfoFile << "     amr.max_level:        " << max_level << "\n";

    jobInfoFile << "     amr.ref_ratio:        ";
    for (int lev = 1; lev <= max_level; lev++) {
        IntVect ref_ratio = parent->refRatio(lev - 1);
        jobInfoFile << ref_ratio[0] << " ";
    }
    jobInfoFile << "\n";

    jobInfoFile << "\n\n";

    // species info
    int mlen = 20;

    jobInfoFile << PrettyLine;
    jobInfoFile << " Species Information\n";
    jobInfoFile << PrettyLine;

    jobInfoFile << std::setw(6) << "index" << SkipSpace << std::setw(mlen + 1) << "name"
                << SkipSpace << std::setw(7) << "A" << SkipSpace << std::setw(7) << "Z"
                << "\n";
    jobInfoFile << OtherLine;

    for (int i = 0; i < NumSpec; i++) {
        jobInfoFile << std::setw(6) << i << SkipSpace << std::setw(mlen + 1) << std::setfill(' ')
                    << short_spec_names_cxx[i] << SkipSpace << std::setw(7) << aion[i] << SkipSpace
                    << std::setw(7) << zion[i] << "\n";
    }
    jobInfoFile << "\n\n";

    // runtime parameters
    jobInfoFile << PrettyLine;
    jobInfoFile << " Inputs File Parameters\n";
    jobInfoFile << PrettyLine;

#include "castro_job_info_tests.H"
#ifdef AMREX_PARTICLES
#include "particles_job_info_tests.H"
#endif

#ifdef GRAVITY
    gravity->output_job_info_params(jobInfoFile);
#endif
#ifdef DIFFUSION
    diffusion->output_job_info_params(jobInfoFile);
#endif

    jobInfoFile.close();

    // now the external parameters
    const int jobinfo_file_length = FullPathJobInfoFile.length();
    Vector<int> jobinfo_file_name(jobinfo_file_length);

    for (int i = 0; i < jobinfo_file_length; i++) {
        jobinfo_file_name[i] = FullPathJobInfoFile[i];
    }

    runtime_pretty_print(jobinfo_file_name.dataPtr(), &jobinfo_file_length);

#ifdef PROB_PARAMS
    prob_params_pretty_print(jobinfo_file_name.dataPtr(), &jobinfo_file_length);
#endif
}

void Castro::writeBuildInfo() {
    std::string PrettyLine = std::string(78, '=') + "\n";
    std::string OtherLine = std::string(78, '-') + "\n";
    std::string SkipSpace = std::string(8, ' ');

    // build information
    std::cout << PrettyLine;
    std::cout << " Castro Build Information\n";
    std::cout << PrettyLine;

    std::cout << "build date:    " << buildInfoGetBuildDate() << "\n";
    std::cout << "build machine: " << buildInfoGetBuildMachine() << "\n";
    std::cout << "build dir:     " << buildInfoGetBuildDir() << "\n";
    std::cout << "AMReX dir:     " << buildInfoGetAMReXDir() << "\n";

    std::cout << "\n";

    std::cout << "COMP:          " << buildInfoGetComp() << "\n";
    std::cout << "COMP version:  " << buildInfoGetCompVersion() << "\n";

    std::cout << "\n";

    std::cout << "C++ compiler:  " << buildInfoGetCXXName() << "\n";
    std::cout << "C++ flags:     " << buildInfoGetCXXFlags() << "\n";

    std::cout << "\n";

    std::cout << "Fortran comp:  " << buildInfoGetFName() << "\n";
    std::cout << "Fortran flags: " << buildInfoGetFFlags() << "\n";

    std::cout << "\n";

    std::cout << "Link flags:    " << buildInfoGetLinkFlags() << "\n";
    std::cout << "Libraries:     " << buildInfoGetLibraries() << "\n";

    std::cout << "\n";

    for (int n = 1; n <= buildInfoGetNumModules(); n++) {
        std::cout << buildInfoGetModuleName(n) << ": " << buildInfoGetModuleVal(n) << "\n";
    }

    std::cout << "\n";

    const char* githash1 = buildInfoGetGitHash(1);
    const char* githash2 = buildInfoGetGitHash(2);
    const char* githash3 = buildInfoGetGitHash(3);
    if (strlen(githash1) > 0) {
        std::cout << "Castro       git describe: " << githash1 << "\n";
    }
    if (strlen(githash2) > 0) {
        std::cout << "AMReX        git describe: " << githash2 << "\n";
    }
    if (strlen(githash3) > 0) {
        std::cout << "Microphysics git describe: " << githash3 << "\n";
    }

    const char* buildgithash = buildInfoGetBuildGitHash();
    const char* buildgitname = buildInfoGetBuildGitName();
    if (strlen(buildgithash) > 0) {
        std::cout << buildgitname << " git describe: " << buildgithash << "\n";
    }

    std::cout << "\n\n";
}

void Castro::writePlotFile(const std::string& dir, ostream& os, VisMF::How how) {
    plotFileOutput(dir, os, how, 0);
}

void Castro::writeSmallPlotFile(const std::string& dir, ostream& os, VisMF::How how) {
    plotFileOutput(dir, os, how, 1);
}

void Castro::plotFileOutput(const std::string& dir, ostream& os, VisMF::How how,
                            const int is_small) {
#ifdef AMREX_PARTICLES
    ParticlePlotFile(dir);
#endif

    //
    // The list of indices of State to write to plotfile.
    // first component of pair is state_type,
    // second component of pair is component # within the state_type
    //
    std::vector<std::pair<int, int> > plot_var_map;
    for (int typ = 0; typ < desc_lst.size(); typ++) {
        for (int comp = 0; comp < desc_lst[typ].nComp(); comp++) {
            if (((parent->isStatePlotVar(desc_lst[typ].name(comp)) && is_small == 0) ||
                 (parent->isStateSmallPlotVar(desc_lst[typ].name(comp)) && is_small == 1)) &&
                desc_lst[typ].getType() == IndexType::TheCellType()) {
                plot_var_map.push_back(std::pair<int, int>(typ, comp));
            }
        }
    }

    int num_derive = 0;
    std::list<std::string> derive_names;
    const std::list<DeriveRec>& dlist = derive_lst.dlist();

    for (auto it = dlist.begin(); it != dlist.end(); ++it) {
        if ((parent->isDerivePlotVar(it->name()) && is_small == 0) ||
            (parent->isDeriveSmallPlotVar(it->name()) && is_small == 1)) {
#ifdef AMREX_PARTICLES
            if (it->name() == "particle_count" || it->name() == "total_particle_count") {
                if (Castro::theTracerPC()) {
                    derive_names.push_back(it->name());
                    num_derive = num_derive + it->numDerive();
                }
            } else
#endif
            {
                derive_names.push_back(it->name());
                num_derive = num_derive + it->numDerive();
            }
        }
    }

    int n_data_items = plot_var_map.size() + num_derive;

#ifdef RADIATION
    if (Radiation::nplotvar > 0) n_data_items += Radiation::nplotvar;
#endif

    Real cur_time = state[State_Type].curTime();

    if (level == 0 && ParallelDescriptor::IOProcessor()) {
        //
        // The first thing we write out is the plotfile type.
        //
        os << thePlotFileType() << '\n';

        if (n_data_items == 0) {
            amrex::Error("Must specify at least one valid data item to plot");
        }

        os << n_data_items << '\n';

        //
        // Names of variables -- first state, then derived
        //
        for (int i = 0; i < plot_var_map.size(); i++) {
            int typ = plot_var_map[i].first;
            int comp = plot_var_map[i].second;
            os << desc_lst[typ].name(comp) << '\n';
        }

        for (auto it = derive_names.begin(); it != derive_names.end(); ++it) {
            const DeriveRec* rec = derive_lst.get(*it);
            if (rec->numDerive() > 1) {
                for (int i = 0; i < rec->numDerive(); ++i) {
                    os << rec->variableName(0) + '_' + std::to_string(i) + '\n';
                }
            } else {
                os << rec->variableName(0) << '\n';
            }
        }

#ifdef RADIATION
        for (int i = 0; i < Radiation::nplotvar; ++i) {
            os << Radiation::plotvar_names[i] << '\n';
        }
#endif

        os << BL_SPACEDIM << '\n';
        os << parent->cumTime() << '\n';
        int f_lev = parent->finestLevel();
        os << f_lev << '\n';
        for (int i = 0; i < BL_SPACEDIM; i++) {
            os << geom.ProbLo(i) << ' ';
        }
        os << '\n';
        for (int i = 0; i < BL_SPACEDIM; i++) {
            os << geom.ProbHi(i) << ' ';
        }
        os << '\n';
        for (int i = 0; i < f_lev; i++) {
            os << parent->refRatio(i)[0] << ' ';
        }
        os << '\n';
        for (int i = 0; i <= f_lev; i++) {
            os << parent->Geom(i).Domain() << ' ';
        }
        os << '\n';
        for (int i = 0; i <= f_lev; i++) {
            os << parent->levelSteps(i) << ' ';
        }
        os << '\n';
        for (int i = 0; i <= f_lev; i++) {
            for (int k = 0; k < BL_SPACEDIM; k++) {
                os << parent->Geom(i).CellSize()[k] << ' ';
            }
            os << '\n';
        }
        os << (int)geom.Coord() << '\n';
        os << "0\n";  // Write bndry data.

#ifdef RADIATION
        if (do_radiation && Radiation::do_multigroup) {
            std::ofstream groupfile;
            std::string FullPathGroupFile = dir;
            FullPathGroupFile += "/RadiationGroups";
            groupfile.open(FullPathGroupFile.c_str(), std::ios::out);

            radiation->write_groups(groupfile);

            groupfile.close();
        }
#endif
    }
    // Build the directory to hold the MultiFab at this level.
    // The name is relative to the directory containing the Header file.
    //
    static const std::string BaseName = "/Cell";
    char buf[64];
    sprintf(buf, "Level_%d", level);
    std::string Level = buf;
    //
    // Now for the full pathname of that directory.
    //
    std::string FullPath = dir;
    if (!FullPath.empty() && FullPath[FullPath.size() - 1] != '/') {
        FullPath += '/';
    }
    FullPath += Level;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor()) {
        if (!amrex::UtilCreateDirectory(FullPath, 0755)) {
            amrex::CreateDirectoryFailed(FullPath);
        }
    }
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor()) {
        os << level << ' ' << grids.size() << ' ' << cur_time << '\n';
        os << parent->levelSteps(level) << '\n';

        for (int i = 0; i < grids.size(); ++i) {
            RealBox gridloc = RealBox(grids[i], geom.CellSize(), geom.ProbLo());
            for (int n = 0; n < BL_SPACEDIM; n++) {
                os << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
            }
        }
        //
        // The full relative pathname of the MultiFabs at this level.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
        if (n_data_items > 0) {
            std::string PathNameInHeader = Level;
            PathNameInHeader += BaseName;
            os << PathNameInHeader << '\n';
        }
    }
    //
    // We combine all of the multifabs -- state, derived, etc -- into one
    // multifab -- plotMF.
    // NOTE: we are assuming that each state variable has one component,
    // but a derived variable is allowed to have multiple components.
    int cnt = 0;
    const int nGrow = 0;
    MultiFab plotMF(grids, dmap, n_data_items, nGrow);
    MultiFab* this_dat = 0;
    //
    // Cull data from state variables -- use no ghost cells.
    //
    for (int i = 0; i < plot_var_map.size(); i++) {
        int typ = plot_var_map[i].first;
        int comp = plot_var_map[i].second;
        this_dat = &state[typ].newData();
        MultiFab::Copy(plotMF, *this_dat, comp, cnt, 1, nGrow);
        cnt++;
    }
    //
    // Cull data from derived variables.
    //
    if (dlist.size() > 0) {
        for (auto it = dlist.begin(); it != dlist.end(); ++it) {
            if ((parent->isDerivePlotVar(it->name()) && is_small == 0) ||
                (parent->isDeriveSmallPlotVar(it->name()) && is_small == 1)) {

                auto derive_dat = derive(it->variableName(0), cur_time, nGrow);
                MultiFab::Copy(plotMF, *derive_dat, 0, cnt, it->numDerive(), nGrow);
                cnt = cnt + it->numDerive();
            }
        }
    }

#ifdef RADIATION
    if (Radiation::nplotvar > 0) {
        MultiFab::Copy(plotMF, *(radiation->plotvar[level]), 0, cnt, Radiation::nplotvar, 0);
        cnt += Radiation::nplotvar;
    }
#endif

    //
    // Use the Full pathname when naming the MultiFab.
    //
    std::string TheFullPath = FullPath;
    TheFullPath += BaseName;

    const Real io_start_time = ParallelDescriptor::second();

    if (amrex::AsyncOut::UseAsyncOut()) {
        VisMF::AsyncWrite(std::move(plotMF), TheFullPath);
    } else {
        VisMF::Write(plotMF, TheFullPath, how, true);
    }

    const Real io_time = ParallelDescriptor::second() - io_start_time;

    if (level == 0 && ParallelDescriptor::IOProcessor()) {
        writeJobInfo(dir, io_time);
    }

    if (track_grid_losses && level == 0) {

        // store diagnostic quantities
        std::ofstream DiagFile;
        std::string FullPathDiagFile = dir;
        FullPathDiagFile += "/Diagnostics";
        DiagFile.open(FullPathDiagFile.c_str(), std::ios::out);

        for (int i = 0; i < n_lost; i++) {
            DiagFile << std::setprecision(17) << material_lost_through_boundary_cumulative[i]
                     << std::endl;
        }

        DiagFile.close();
    }

#ifdef GRAVITY
    if (use_point_mass && level == 0) {

        // store current value of the point mass
        std::ofstream PMFile;
        std::string FullPathPMFile = dir;
        FullPathPMFile += "/point_mass";
        PMFile.open(FullPathPMFile.c_str(), std::ios::out);

        PMFile << std::setprecision(17) << point_mass << std::endl;

        PMFile.close();
    }
#endif
}
