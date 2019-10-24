
#ifndef _Castro_H_
#define _Castro_H_

#include <AMReX_BC_TYPES.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_ErrorList.H>
#include <AMReX_FluxRegister.H>

#ifdef BL_LAZY
#include <AMReX_Lazy.H>
#endif

#ifdef AMREX_PARTICLES
#include "AMReX_AmrParticles.H"
#endif

#include <memory>
#include <iostream>

using std::istream;
using std::ostream;


enum StateType { State_Type = 0,
#ifdef RADIATION
                 Rad_Type,
#endif
#ifdef SELF_GRAVITY
                 PhiGrav_Type,
                 Gravity_Type,
#endif
#ifdef ROTATION
		 PhiRot_Type,
		 Rotation_Type,
#endif
		 Source_Type,
#ifdef REACTIONS
                 Reactions_Type,
#endif
#ifdef REACTIONS
		 Simplified_SDC_React_Type
#endif
#ifdef THORNADO
		 Thornado_Type,
		 Thornado_Fluid_Source_Type,
		 Thornado_Rad_Source_Type
#endif
};

// Create storage for all source terms.

enum sources { ext_src = 0,
               thermo_src,
#ifdef SPONGE
               sponge_src,
#endif
#ifdef DIFFUSION
               diff_src,
#endif
#ifdef HYBRID_MOMENTUM
               hybrid_src,
#endif
#ifdef GRAVITY
               grav_src,
#endif
#ifdef ROTATION
               rot_src,
#endif
	       num_src };


// time integration method

enum int_method { CornerTransportUpwind = 0,
                  UnusedTimeIntegration,
                  SpectralDeferredCorrections,
                  SimplifiedSpectralDeferredCorrections
                };

///
/// @class Castro
///
/// @brief ``AmrLevel``-derived class for hyperbolic conservation equations for stellar media
///
class Castro
    :
    public amrex::AmrLevel
{
public:

///
/// Default constructor.  Builds invalid object.
///
    Castro ();

///
/// The basic constructor.
///
/// @param papa     parent ``amrex::Amr`` object
/// @param lev      which level are we on
/// @param level_geom level geometry
/// @param bl       ``BoxArray`` object
/// @param dm       the mapping
/// @param time     current time
///
    Castro (amrex::Amr&            papa,
            int             lev,
            const amrex::Geometry& level_geom,
            const amrex::BoxArray& bl,
            const amrex::DistributionMapping& dm,
            amrex::Real            time);

///
/// The destructor.
///
    virtual ~Castro () override;

///
/// Restart from a checkpoint file.
///
/// @param papa             parent ``amrex::Amr`` object
/// @param is
/// @param bReadSpecial
///
    virtual void restart (amrex::Amr&     papa,
                          istream& is,
			  bool bReadSpecial = false) override;

///
/// This is called only when we restart from an old checkpoint.
/// Sets the Vector `state_in_checkpoint` depending on the input
/// version and type.
///
/// @param state_in_checkpoint      Vector of integer flags
///                                 describing the state
///
    virtual void set_state_in_checkpoint (amrex::Vector<int>& state_in_checkpoint) override;

///
/// Call ``amrex::AmrLevel::checkPoint`` and then add radiation info
///
/// @param dir          Directory to store checkpoint in
/// @param os           ``std::ostream`` object
/// @param how          `V`isMF::How`` object
/// @param dump_old     boolean
///
    virtual void checkPoint(const std::string& dir,
                            std::ostream&      os,
                            amrex::VisMF::How         how,
                            bool               dump_old) override;

///
/// A string written as the first item in writePlotFile() at
/// level zero. It is so we can distinguish between different
/// types of plot files. For Castro it has the form: Castro-Vnnn
///
    virtual std::string thePlotFileType () const override;

    virtual void setPlotVariables () override;

///
/// Write a plotfile to specified directory.
///
/// @param dir  directory to write to
/// @param os   output stream
/// @param how  ``VisMF::How`` object
///
    virtual void writePlotFile (const std::string& dir,
                                ostream&       os,
                                amrex::VisMF::How     how) override;

///
/// @param dir  directory to write to
/// @param os   output stream
/// @param how  ``VisMF::How`` object
///
    virtual void writeSmallPlotFile (const std::string& dir,
				     ostream&       os,
				     amrex::VisMF::How     how) override;

///
/// @param dir      directory to write to
/// @param os       output stream
/// @param how      ``VisMF::How`` object
/// @param is_small integer
///
    void plotFileOutput(const std::string& dir,
                        ostream& os,
                        amrex::VisMF::How how,
                        const int is_small);


///
/// Write job info to file
///
/// @param dir       directory to write to
/// @param io_time   how long it took to generate the plotfile/checkpoint
///
    void writeJobInfo (const std::string& dir, const amrex::Real io_time);


///
/// Dump build info 
///
    static void writeBuildInfo ();


///
/// Define data descriptors.
///
/// ::variableSetUp is called in the constructor of Amr.cpp, so
/// it should get called every time we start or restart a job
///
    static void variableSetUp ();

///
/// Define tagging functions.
///
    static void ErrorSetUp ();

///
/// Cleanup data descriptors at end of run.
///
/// ::variableSetUp is in ``Castro_setup.cpp``
/// ``variableCleanUp`` is called once at the end of a simulation
///
    static void variableCleanUp ();

///
/// Initialize grid data at problem start-up.
///
    virtual void initData () override;

#ifdef AMREX_PARTICLES

///
/// Read particle-related inputs
///
    static void read_particle_params ();

///
/// Initialize particle locations and velocities (and strengths if relevant)
///
    void init_particles ();

///
/// Write particles in checkpoint directories
///
/// @param dir  directory containing checkpoint
///
    void ParticleCheckPoint (const std::string& dir);

///
/// Write particles in plotfile directories
///
/// @param dir  directory containing plotfile
///
    void ParticlePlotFile (const std::string& dir);

///
/// How to initialize at restart
///
/// @param restart_file     path to restart file
///
    void ParticlePostRestart (const std::string& restart_file);

///
/// Derived quantities associated with particles
///
/// @param name         name of particle variable to derive
/// @param time         current time
/// @param ngrow        Number of ghost cells
///
    std::unique_ptr<amrex::MultiFab> ParticleDerive (const std::string& name,
						     amrex::Real           time,
						     int            ngrow);

///
/// Timestamp particles
///
/// @param ngrow
///
    void TimestampParticles (int ngrow);

///
/// Advance the particles by dt
///
/// @param iteration    where we are in the current AMR subcycle
/// @param time         current time
/// @param dt           timestep
///
    void advance_particles (int iteration, amrex::Real time, amrex::Real dt);

#endif

#ifdef MAESTRO_INIT
    void MAESTRO_init ();
#endif

#ifdef THORNADO
    static int init_thornado();
    void init_thornado_data ();
    void init_thornado_on_patch();
    void create_thornado_source(amrex::Real dt);
#endif

///
///
///
    amrex::MultiFab* Area ();

///
/// Area of the multifab with surface normal along ``dir``
///
/// @param dir  direction of normal
///
    amrex::MultiFab& Area (int dir);

///
/// The volume of the multifab.
///
    amrex::MultiFab& Volume ();

///
/// Return the n'th fluxes MultiFab.
///
/// @param dir  direction in which to find the fluxes along
///
    amrex::MultiFab& Fluxes (int dir);

///
/// Set time levels of state data.
///
/// @param time         current time
/// @param dt_old       old timestep
/// @param dt_new       new timestep
///
    virtual void setTimeLevel (amrex::Real time,
                               amrex::Real dt_old,
                               amrex::Real dt_new) override;

///
/// Initialize data on this level from another Castro (during regrid).
///
/// @param old  amrex::AmrLevel object to use to initialize data
///
    virtual void init (amrex::AmrLevel& old) override;

///
/// Initialize data on this level after regridding if old level did not
/// previously exist
///
    virtual void init () override;

///
/// Proceed with next timestep?
///
    virtual int okToContinue () override;

///
/// Advance grids at this level in time.
///
/// The main driver for a single level.  This will do either the SDC
/// algorithm or the Strang-split reactions algorithm.
///
/// @param time     the current simulation time
/// @param dt       the timestep to advance (e.g., go from time to
///                    time + dt)
/// @param iteration    where we are in the current AMR subcycle.  Each
///                    level will take a number of steps to reach the
///                    final time of the coarser level below it.  This
///                    counter starts at 1
/// @param ncycle   the number of subcycles at this level
///
    virtual amrex::Real advance (amrex::Real time,
                                 amrex::Real dt,
                                 int  iteration,
                                 int  ncycle) override;

///
/// This routine will advance the old state data (called ``S_old`` here)
/// to the new time, for a single level.  The new data is called
/// ``S_new`` here.  The update includes reactions (if we are not doing
/// SDC), hydro, and the source terms. The advance is doing using the CTU
/// (corner transport upwind) approach. We also use this interface for the
/// simplified SDC update.
///
/// @param time     the current simulation time
/// @param dt       the timestep to advance (e.g., go from time to
///                    time + dt)
/// @param amr_iteration    where we are in the current AMR subcycle.  Each
///                    level will take a number of steps to reach the
///                    final time of the coarser level below it.  This
///                    counter starts at 1
/// @param amr_ncycle   the number of subcycles at this level
///
    amrex::Real do_advance_ctu(amrex::Real time,
                               amrex::Real dt,
                               int  amr_iteration,
                               int  amr_ncycle);


#ifndef AMREX_USE_CUDA
    amrex::Real do_advance_sdc (amrex::Real time,
                                amrex::Real dt,
                                int  amr_iteration,
                                int  amr_ncycle);
#endif

///
/// Should we retry advancing the simulation? By default, we
/// don't do a retry unless the criteria are violated. This is
/// used only for the CTU advance.
///
/// @param time     the current simulation time
/// @param dt       the timestep to advance (e.g., go from time to
///                    time + dt)
/// @param amr_iteration    where we are in the current AMR subcycle.  Each
///                    level will take a number of steps to reach the
///                    final time of the coarser level below it.  This
///                    counter starts at 1
/// @param amr_ncycle   the number of subcycles at this level
///
    bool retry_advance_ctu(amrex::Real& time, amrex::Real dt, int amr_iteration, int amr_ncycle);

///
/// Subcyles until we've reached the target time, ``time`` + ``dt``.
/// The last timestep will be shortened if needed so that
/// it does not overshoot the ending time. This is used only for
/// the CTU advance.
///
/// @param time     the current simulation time
/// @param dt       the timestep to advance (e.g., go from time to
///                    time + dt)
/// @param amr_iteration    where we are in the current AMR subcycle.  Each
///                    level will take a number of steps to reach the
///                    final time of the coarser level below it.  This
///                    counter starts at 1
/// @param amr_ncycle   the number of subcycles at this level
///
    amrex::Real subcycle_advance_ctu(amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle);

///
/// This is run at the start of an ::advance on a single level.
/// It sets up all the necessary arrays and data containers,
/// ensures the data is valid and makes copies of the MultiFabs
/// in the old and new state in case a retry is performed.
///
/// @param time     the current simulation time
/// @param dt       the timestep to advance (e.g., go from time to
///                    time + dt)
/// @param amr_iteration    where we are in the current AMR subcycle.  Each
///                    level will take a number of steps to reach the
///                    final time of the coarser level below it.  This
///                    counter starts at 1
/// @param amr_ncycle   the number of subcycles at this level
///
    void initialize_advance(amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle);

///
/// If required, adds material lost in timestep to cumulative
/// losses, does refluxing, and clears arrays/data containers
/// initialized at start of the level's advance.
///
/// @param time     the current simulation time
/// @param dt       the timestep to advance (e.g., go from time to
///                    time + dt)
/// @param amr_iteration    where we are in the current AMR subcycle.  Each
///                    level will take a number of steps to reach the
///                    final time of the coarser level below it.  This
///                    counter starts at 1
/// @param amr_ncycle   the number of subcycles at this level
///
    void finalize_advance(amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle);

///
/// Performed at the start of ::do_advance, resets flags and
/// ensures required ghost zones are filled.
///
/// @param time     the current simulation time
/// @param dt       the timestep to advance (e.g., go from time to
///                    time + dt)
/// @param amr_iteration    where we are in the current AMR subcycle.  Each
///                    level will take a number of steps to reach the
///                    final time of the coarser level below it.  This
///                    counter starts at 1
/// @param amr_ncycle   the number of subcycles at this level
///
    void initialize_do_advance(amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle);

///
/// Performed at the end of ::do_advance, cleans up temporary data
/// created by ::initialize_do_advance
///
/// @param time     the current simulation time
/// @param dt       the timestep to advance (e.g., go from time to
///                    time + dt)
/// @param amr_iteration    where we are in the current AMR subcycle.  Each
///                    level will take a number of steps to reach the
///                    final time of the coarser level below it.  This
///                    counter starts at 1
/// @param amr_ncycle   the number of subcycles at this level
///
    void finalize_do_advance(amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle);

///
/// Check for NaNs in the given MultiFab
///
/// @param state        MultiFab to check
/// @param check_ghost  do we check the ghost cells?
///
    void check_for_nan(amrex::MultiFab& state, int check_ghost=0);

///
/// Returns true if flag corresponding to source type ``src`` is set.
///
/// @param src      integer, index corresponding to source type
///
    bool source_flag(int src);

///
/// Returns whether any sources are actually applied.
///
    bool apply_sources();

    static int get_output_at_completion();


///
/// Construct source terms at old time
///
/// @param source   MultiFab to save sources to
/// @param state    Current state
/// @param time     the current simulation time
/// @param dt       the timestep to advance (e.g., go from time to
///                    time + dt)
/// @param amr_iteration    where we are in the current AMR subcycle.  Each
///                    level will take a number of steps to reach the
///                    final time of the coarser level below it.  This
///                    counter starts at 1
/// @param amr_ncycle   the number of subcycles at this level
///
    void do_old_sources(amrex::MultiFab& source, amrex::MultiFab& state, amrex::Real time, amrex::Real dt, int amr_iteration = -1, int amr_ncycle = -1);

///
/// Construct the old-time source
///
/// @param src      integer corresponding to source type to construct
/// @param source   MultiFab to save sources to
/// @param state    State data
/// @param time     the current simulation time
/// @param dt       the timestep to advance (e.g., go from time to
///                    time + dt)
/// @param amr_iteration    where we are in the current AMR subcycle.  Each
///                    level will take a number of steps to reach the
///                    final time of the coarser level below it.  This
///                    counter starts at 1
/// @param amr_ncycle   the number of subcycles at this level
///
    void construct_old_source(int src, amrex::MultiFab& source, amrex::MultiFab& state, amrex::Real time, amrex::Real dt, int amr_iteration = -1, int amr_ncycle = -1);


///
/// Construct new time sources
///
/// @param source   MultiFab to save sources to
/// @param state_old    Old state
/// @param state_new    New state
/// @param time     the current simulation time
/// @param dt       the timestep to advance (e.g., go from time to
///                    time + dt)
/// @param amr_iteration    where we are in the current AMR subcycle.  Each
///                    level will take a number of steps to reach the
///                    final time of the coarser level below it.  This
///                    counter starts at 1
/// @param amr_ncycle   the number of subcycles at this level
///
    void do_new_sources(amrex::MultiFab& source, amrex::MultiFab& state_old, amrex::MultiFab& state_new, amrex::Real time, amrex::Real dt, int amr_iteration = -1, int amr_ncycle = -1);

///
/// Construct the new-time sources.
///
/// @param src          Integer corresponding to type of source to construct
/// @param source       MultiFab to save source to
/// @param state_old    Old state
/// @param state_new    New state
/// @param time         the current simulation time
/// @param dt       t   he timestep to advance (e.g., go from time to
///                    time + dt)
/// @param amr_iteration    where we are in the current AMR subcycle.  Each
///                    level will take a number of steps to reach the
///                    final time of the coarser level below it.  This
///                    counter starts at 1
/// @param amr_ncycle   the number of subcycles at this level
///
    void construct_new_source(int src, amrex::MultiFab& source, amrex::MultiFab& state_old, amrex::MultiFab& state_new, amrex::Real time, amrex::Real dt, int amr_iteration = -1, int amr_ncycle = -1);

///
/// Evaluate diagnostics quantities describing the effect of an
/// update on the state. The optional parameter local determines
/// whether we want to do this calculation globally over all processes
/// or locally just on this processor. The latter is useful if you
/// are evaluating the contribution from multiple source changes at once
/// and want to amortize the cost of the parallel reductions.
/// Note that the resultant output is volume-weighted.
///
/// @param update   MultiFab containing source
/// @param dt       timestep
/// @param local    boolean, do we want to do this calculation
///                 globally over all processes
///                 or locally just on this processor?
///
    amrex::Vector<amrex::Real> evaluate_source_change(amrex::MultiFab& update, amrex::Real dt, bool local = false);

///
/// Print the change due to a given source term update.
/// We assume here that the input array is lined up with
/// the ``NUM_STATE`` components of ``State_Type`` because we are
/// interested in printing changes to energy, mass, etc.
///
/// @param update       Real Vector of changes
///
    void print_source_change(amrex::Vector<amrex::Real> update);

///
/// For the old-time or new-time sources update, evaluate the change in the state
/// for all source terms, then print the results.
///
/// @param dt       timestep
/// @param is_new   boolean, do we get the new or old state data?
///
    void print_all_source_changes(amrex::Real dt, bool is_new);


///
/// Obtain the sum of all source terms.
///
/// this computes advective_source + 1/2 (old source + new source)
///
/// @note the advective source is defined as \f$-\nabla\cdot F\f$
///
/// the time-centering is accomplished since new source is defined
/// to be 1/2 (new source - old source) generally.
///
/// @param source       MultiFab of source terms to be summed
///
    void sum_of_sources(amrex::MultiFab& source);


///
/// Center source terms in time, storing in ``S_new``
///
/// @param S_new    MultiFab to store centered sources in
/// @param src_old  Old source term
/// @param src_new  New source term
/// @param dt       timestep
///
    void time_center_source_terms (amrex::MultiFab& S_new,
				   amrex::MultiFab& src_old,
				   amrex::MultiFab& src_new,
				   amrex::Real dt = 1.0);

#ifdef SPONGE

///
/// Construct sponge source terms at old timestep
///
/// @param source   MultiFab to save source terms to
/// @param state    Old state
/// @param time     current time
/// @param dt       timestep
///
    void construct_old_sponge_source(amrex::MultiFab& source, amrex::MultiFab& state, amrex::Real time, amrex::Real dt);


///
/// Construct sponge source terms at new timestep
///
/// @param source       MultiFab to save source terms to
/// @param state_old    Old state
/// @param state_new    New state
/// @param time         current time
/// @param dt           timestep
///
    void construct_new_sponge_source(amrex::MultiFab& source, amrex::MultiFab& state_old, amrex::MultiFab& state_new, amrex::Real time, amrex::Real dt);

///
/// Allocate sponge parameters
///
    static void sponge_init();

///
/// Deallocate sponge parameters
///
    static void sponge_finalize();
#endif


///
/// Create external source terms at old timestep
///
/// @param source   MultiFab to save source terms to
/// @param state    Old state
/// @param time     current time
/// @param dt       timestep
///
    void construct_old_ext_source(amrex::MultiFab& source, amrex::MultiFab& state, amrex::Real time, amrex::Real dt);


///
/// Create external source terms at new timestep
///
/// @param source       MultiFab to save source terms to
/// @param state_old    Old state
/// @param state_new    New state
/// @param time         current time
/// @param dt           timestep
///
    void construct_new_ext_source(amrex::MultiFab& source, amrex::MultiFab& state_old, amrex::MultiFab& state_new, amrex::Real time, amrex::Real dt);


///
/// Fill ``ext_src`` with external sources
///
/// @param time     current time
/// @param dt       timestep
/// @param S_old    Old state
/// @param S_new    New state
/// @param ext_src MultiFab to save sources to
///
    void fill_ext_source(amrex::Real time, amrex::Real dt, amrex::MultiFab& S_old, amrex::MultiFab& S_new, amrex::MultiFab& ext_src);


///
/// Construct thermal sources at old timestep
///
/// @param source   MultiFab to save sources to
/// @param state    Old state
/// @param time     current time
/// @param dt       timestep
///
    void construct_old_thermo_source(amrex::MultiFab& source, amrex::MultiFab& state,
                                     amrex::Real time, amrex::Real dt);


///
/// Construct thermal sources at new timestep
///
/// @param source       MultiFab to save sources to
/// @param state_old    Old state
/// @param state_new    New state
/// @param time         current time
/// @param dt           timestep
///
    void construct_new_thermo_source(amrex::MultiFab& source, amrex::MultiFab& state_old,
                                     amrex::MultiFab& state_new, amrex::Real time, amrex::Real dt);


///
/// Fill ``ext_src`` with thermal sources
///
/// @param time     current time
/// @param dt       timestep
/// @param S_old    Old state
/// @param S_new    New state
/// @param ext_src  MultiFab to fill with sources
///
    void fill_thermo_source(amrex::Real time, amrex::Real dt,
                            amrex::MultiFab& S_old, amrex::MultiFab& S_new,
                            amrex::MultiFab& ext_src);

#ifdef GRAVITY
#ifdef SELF_GRAVITY

///
/// Construct gravitational field at old timestep
///
/// @param amr_iteration    where we are in the current AMR subcycle
/// @param amr_ncycle       the number of subcycles at this level
/// @param time             current time
///
    void construct_old_gravity(int amr_iteration, int amr_ncycle,
			       amrex::Real time);


///
/// Construct gravitational field at new timestep
///
/// @param amr_iteration    where we are in the current AMR subcycle
/// @param amr_ncycle       the number of subcycles at this level
/// @param time             current time
///
    void construct_new_gravity(int amr_iteration, int amr_ncycle,
			       amrex::Real time);
#endif

///
/// Construct gravitational source terms at old timestep
///
/// @param source   MultiFab to save sources to
/// @param state    Old state
/// @param time     current time
/// @param dt       timestep
///
    void construct_old_gravity_source(amrex::MultiFab& source, amrex::MultiFab& state, amrex::Real time, amrex::Real dt);


///
/// Construct gravitational source terms at new timestep
///
/// @param source       MultiFab to save sources to
/// @param state_old    Old state
/// @param state_new    New state
/// @param time         current time
/// @param dt           timestep
///
    void construct_new_gravity_source(amrex::MultiFab& source, amrex::MultiFab& state_old, amrex::MultiFab& state_new, amrex::Real time, amrex::Real dt);
#endif


///
/// Estimate time step.
///
/// @param dt_old   previous timestep
///
    amrex::Real estTimeStep (amrex::Real dt_old);

///
/// Compute initial time step.
///
    amrex::Real initialTimeStep ();

///
/// Compute initial ``dt``.
///
/// @param finest_level     Index of finest level
/// @param sub_cycle
/// @param n_cycle          the number of subcycles at this level
/// @param ref_ratio        Vector of refinement ratios between levels
/// @param dt_level         Real Vector, will contain initial timestep at each level.
/// @param stop_time        End time of simulation
///
    virtual void computeInitialDt (int                   finest_level,
                                   int                   sub_cycle,
                                   amrex::Vector<int>&           n_cycle,
                                   const amrex::Vector<amrex::IntVect>& ref_ratio,
                                   amrex::Vector<amrex::Real>&          dt_level,
                                   amrex::Real                  stop_time) override;

///
/// Compute new ``dt``.
///
/// @param finest_level     Index of finest level
/// @param sub_cycle
/// @param n_cycle          the number of subcycles at this level
/// @param ref_ratio        Vector of refinement ratios between levels
/// @param dt_min
/// @param dt_level
/// @param stop_time        End time of simulation
/// @param post_regrid_flag Have we done regridding yet?
///
    virtual void computeNewDt (int                   finest_level,
                               int                   sub_cycle,
                               amrex::Vector<int>&           n_cycle,
                               const amrex::Vector<amrex::IntVect>& ref_ratio,
                               amrex::Vector<amrex::Real>&          dt_min,
                               amrex::Vector<amrex::Real>&          dt_level,
                               amrex::Real                  stop_time,
                               int                   post_regrid_flag) override;

///
/// Allocate data at old time.
///
    virtual void allocOldData () override;

///
/// Remove data at old time.
///
    virtual void removeOldData () override;

///
/// Passes some data about the grid to a Fortran module.
///
/// Send refinement data to Fortran. We do it here
/// because now the grids have been initialized and
/// we need this data for setting up the problem.
///
/// @note this routine will always get called
/// on level 0, even if we are doing a restart,
/// so it is safe to put this here.
///
    void setGridInfo();

///
/// Print information about energy budget.
///
    void do_energy_diagnostics();

///
/// Do work after timestep().
///
/// @param iteration    where we are in the current AMR subcycle
///
    virtual void post_timestep (int iteration) override;

///
/// Contains operations to be done only after a full coarse timestep.
///
/// @param cumtime
///
    virtual void postCoarseTimeStep (amrex::Real cumtime) override;

///
/// Check whether we want to do a regrid.
///
/// @param time     current time
///
    void check_for_post_regrid (amrex::Real time);

///
/// Do work after regrid().
///
/// @param lbase
/// @param new_finest
///
    virtual void post_regrid (int lbase,
                              int new_finest) override;

///
/// Do work after a restart().
///
    virtual void post_restart () override;

///
/// Do work after init().
///
/// @param stop_time    End time of simulation
///
    virtual void post_init (amrex::Real stop_time) override;

#ifdef GRAVITY
#ifdef ROTATION
///
/// Initialize a model using the self-consistent field method.
///
    void scf_relaxation();

///
/// Perform the Hachisu SCF method.
///
    void do_hscf_solve();
#endif
#endif

///
/// Do work after restart with ``grown_factor`` > 1
///
    void post_grown_restart ();

#ifdef DO_PROBLEM_POST_SIMULATION

///
/// Do work at the end of the simulation - before the last outputs
///
/// @param amr_level
///
    static void problem_post_simulation (amrex::Vector<std::unique_ptr<AmrLevel> >& amr_level);
#endif


///
/// Error estimation for regridding.
///
/// @param tb           TagBoxArray of tags
/// @param clearval     integer
/// @param tagval       integer
/// @param time         current time
/// @param n_error_buf  integer
/// @param ngrow        Number of ghost cells
///
    virtual void errorEst (amrex::TagBoxArray& tb,
                           int          clearval,
                           int          tagval,
                           amrex::Real  time,
			   int          n_error_buf = 0,
			   int          ngrow = 0) override;


///
/// Apply the ::set_problem_tags routine.
///
/// @param tags         TagBoxArray of tags
/// @param time         current time
///
    void apply_problem_tags (amrex::TagBoxArray& tags, amrex::Real time);


///
/// Apply a given tagging function.
///
/// @param tags         TagBoxArray of tags
/// @param time         current time
/// @param j            integer
///
    void apply_tagging_func (amrex::TagBoxArray& tags, amrex::Real time, int j);


///
/// Returns a MultiFab containing the derived data for this level.
/// If ngrow>0 the MultiFab is built on the appropriately grown BoxArray.
///
/// @param name     Name of derived data
/// @param time     Current time
/// @param ngrow    Number of ghost cells
///
    virtual std::unique_ptr<amrex::MultiFab> derive (const std::string& name,
						     amrex::Real        time,
						     int                ngrow) override;

///
/// This version of derive() fills the dcomp'th component of mf with the derived quantity.
///
/// @param name         Name of quantity to derive
/// @param time         current time
/// @param mf           MultiFab to store derived quantity in
/// @param dcomp        index of component of `mf` to fill with derived quantity
///
    virtual void derive (const std::string& name,
			 amrex::Real               time,
			 amrex::MultiFab&          mf,
			 int                dcomp) override;

    static int numGrow();

    static int numState();

#ifdef REACTIONS
///
/// React state through a single reaction timestep, saving reaction source terms to ``reactions`` MultiFab. This is the CTU version.
///
/// @param state        Current state
/// @param reactions    MultiFab to save reaction sources to
/// @param mask
/// @param weights
/// @param time         current time
/// @param dt_react     reaction timestep
/// @param strang_half  Is this a Strang-split half timestep?
/// @param ngrow        Number of ghost cells
///
    void react_state(amrex::MultiFab& state,
		     amrex::MultiFab& reactions,
		     const amrex::iMultiFab& mask,
		     amrex::MultiFab& weights,
		     amrex::Real time,
		     amrex::Real dt_react,
		     int strang_half, int  ngrow = 0);

///
/// Simplified SDC version of react_state. Reacts the current state through a single timestep.
///
/// @param time     current time
/// @param dt       timestep
///
    void react_state(amrex::Real time, amrex::Real dt);

///
/// React through first half timestep
///
/// @param time     current time
/// @param dt       timestep
///
    void strang_react_first_half(amrex::Real time, amrex::Real dt);


///
/// React through second half stimestep
///
/// @param time     current time
/// @param dt       timestep
///
    void strang_react_second_half(amrex::Real time, amrex::Real dt);

///
/// Obtain the effective source term due to reactions on the primitive variables
///
/// @param source   MultiFab to save source term to
/// @param dt       timestep
///
    void get_react_source_prim(amrex::MultiFab& source, amrex::Real time, amrex::Real dt);

///
/// @param state
/// @param react_mf
/// @param time     current time
/// @param dt       timestep
///
    void strang_chem (amrex::MultiFab&  state,
                      amrex::MultiFab&  react_mf,
                      amrex::Real       time,
                      amrex::Real       dt);

///
/// Are there any zones in ``State`` that can burn?
///
/// @param State    State MultiFab
///
    bool valid_zones_to_burn(amrex::MultiFab& State);
#endif

#ifdef ROTATION

///
/// Construct rotation source terms at old time
///
/// @param source   MultiFab to save source terms to
/// @param state    Old state
/// @param time     current time
/// @param dt       timestep
///
    void construct_old_rotation_source(amrex::MultiFab& source, amrex::MultiFab& state, amrex::Real time, amrex::Real dt);


///
/// Construct rotation source terms at new time.
///
/// @param source       MultiFab to save source terms to
/// @param state_old    Old State
/// @param state_new    New State
/// @param time         current time
/// @param dt           timestep
///
    void construct_new_rotation_source(amrex::MultiFab& source, amrex::MultiFab& state_old, amrex::MultiFab& state_new, amrex::Real time, amrex::Real dt);


///
/// Fill rotational potential and acceleration
///
/// @param phi      Rotational potential
/// @param rot      Rotational acceleration
/// @param state    Current state
/// @param time     current time
///
    void fill_rotation_field(amrex::MultiFab& phi, amrex::MultiFab& rot, amrex::MultiFab& state, amrex::Real time);
#endif

#ifdef RADIATION

///
/// @param S_new
/// @param iteration    where we are in the current AMR subcycle
/// @param ncycle       the number of subcycles at this level
///
    void final_radiation_call (amrex::MultiFab& S_new,
                               int iteration,
                               int ncycle);
#endif


///
/// Reset the internal energy
///
/// @param State    Current state
/// @param ng       number of ghost cells
///
    void reset_internal_energy (amrex::MultiFab& State, int ng);



///
/// Compute the current temperature
///
/// @param state    the state to operate on
/// @param time     current time
/// @param ng       number of ghost cells
///
    void computeTemp(amrex::MultiFab& state, amrex::Real time, int ng);


///
/// Optionally predict the source terms to ``t + dt/2``,
/// which is the time-level ``n+1/2`` value, To do this we use a
/// lagged predictor estimate: \f$ \frac{dS}{dt_n} = (S_n - S_{n-1}) / dt \f$, so
/// \f$ S_{n+1/2} = S_n + (dt / 2) \frac{dS}{dt_n} \f$. We'll add the ``S_n``
/// terms later; now we add the second term. We defer the
/// multiplication by ``dt / 2`` to ::initialize_do_advance since
/// for a retry we may not yet know at this point what the
/// advance timestep is.
///
/// Note that since for ``dS/dt`` we want \f$ (S^{n+1} - S^{n}) / dt \f$,
/// we only need to take twice the new-time source term, since in
/// the predictor-corrector approach, the new-time source term is
/// \f$ \frac{1}{2} S^{n+1} - \frac{1}{2} S^{n} \f$. This is
/// untrue in general for the
/// non-momentum sources, so for safety we'll only apply it to the
/// momentum sources.
///
/// Note that if the old data doesn't exist yet (e.g. it is
/// the first step of the simulation) ``FillPatch`` will just
/// return the new data, so this is a valid operation and
/// the result will be zero, so there is no source term
/// prediction in the first step.
///
    void apply_source_term_predictor();


///
/// This is a hack to make sure that we only
/// ever have new data for certain state types that only
/// ever need new time data; by doing a swap now, we'll
/// guarantee that allocOldData() does nothing. We do
/// this because we never need the old data, so we
/// don't want to allocate memory for it.
///
/// @param dt   timestep
///
    void swap_state_time_levels (const amrex::Real dt);

#ifdef DIFFUSION

///
/// Construct diffusion source at old time
///
/// @param source   MultiFab to save source terms to
/// @param state    Old state
/// @param time     current time
/// @param dt       timestep
///
    void construct_old_diff_source(amrex::MultiFab& source, amrex::MultiFab& state, amrex::Real time, amrex::Real dt);

///
/// Construct diffusion source at new time
///
/// @param source       MultiFab to save source terms to
/// @param state_old    Old state
/// @param state_new    New state
/// @param time         current time
/// @param dt           timestep
///
    void construct_new_diff_source(amrex::MultiFab& source, amrex::MultiFab& state_old, amrex::MultiFab& state_new, amrex::Real time, amrex::Real dt);


///
/// Get thermal conductivity diffusion term at given time
///
/// @param time     current time
/// @param state    Current state
/// @param DiffTerm MultiFab to save term to
///
    void getTempDiffusionTerm (amrex::Real time, amrex::MultiFab& state, amrex::MultiFab& DiffTerm);

///
/// Get enthalpy conductivity diffusion term at given time
///
/// @param time     current time
/// @param state    Current state
/// @param DiffTerm MultiFab to save term to
///
    void getEnthDiffusionTerm (amrex::Real time, amrex::MultiFab& state, amrex::MultiFab& DiffTerm);
#if (BL_SPACEDIM == 1)

///
/// Get species' diffusion terms at given time
///
/// @param time     current time
/// @param state    Current state
/// @param DiffTerm MultiFab to save term to
///
    void getSpecDiffusionTerm (amrex::Real time, amrex::MultiFab& state, amrex::MultiFab& DiffTerm);
#endif
#if (BL_SPACEDIM == 1)

///
/// Get viscous term at given time
///
/// @param time                     current time
/// @param state                    Current state
/// @param ViscousTermforMomentum   Viscous term for momentum equation
/// @param ViscousTermforEnergy     Viscous term for energy equation
///
    void getViscousTerm (amrex::Real time, amrex::MultiFab& state, amrex::MultiFab& ViscousTermforMomentum, amrex::MultiFab& ViscousTermforEnergy);

///
///
///
/// @param time         current time
/// @param state        Current state
/// @param ViscousTerm  MultiFab to save term to
///
    void getFirstViscousTerm (amrex::Real time, amrex::MultiFab& state, amrex::MultiFab& ViscousTerm);

///
///
///
/// @param time         current time
/// @param state        Current state
/// @param ViscousTerm  MultiFab to save term to
///
    void getSecndViscousTerm (amrex::Real time, amrex::MultiFab& state, amrex::MultiFab& ViscousTerm);

///
/// @param time                     current time
/// @param state                    Current state
/// @param ViscousTermforEnergy     MultiFab to save term to
///
    void getViscousTermForEnergy (amrex::Real time, amrex::MultiFab& state, amrex::MultiFab& ViscousTermforEnergy);
#endif

///
/// Calculate temperature or enthalpty diffusion terms and add to ``ext_src`` (multiplied by ``mult_factor``).
///
/// @param ext_src      Source terms to add diffusion sources to
/// @param state        Current state
/// @param DiffTerm     MultiFab to save diffusion sources to
/// @param t            Current time
/// @param mult_factor  Real, factor to multiply diffusion
///                     sources by before adding to other sources.
///
    void add_temp_diffusion_to_source (amrex::MultiFab& ext_src, amrex::MultiFab& state, amrex::MultiFab& DiffTerm, amrex::Real t, amrex::Real mult_factor = 1.0);
#if (BL_SPACEDIM == 1)

///
/// Calculate species' diffusion source terms and add to ``ext_src``
///
/// @param ext_src      MultiFab of source terms to add diffusion terms to
/// @param state        Current state
/// @param DiffTerm     MultiFab to save diffusion terms to
/// @param t            Current time
/// @param mult_factor  Real, factor to multiply diffusion term
///                     by before addding to other sources
///
    void add_spec_diffusion_to_source (amrex::MultiFab& ext_src, amrex::MultiFab& state, amrex::MultiFab& DiffTerm, amrex::Real t, amrex::Real mult_factor = 1.0);
#endif
#if (BL_SPACEDIM == 1)

///
/// Calculate viscous source terms and add to ``ext_src``
///
/// @param ext_src                  MultiFab to add viscous sources to
/// @param state                    Current state
/// @param ViscousTermforMomentum   MultiFab to save viscous term in momentum equation to
/// @param ViscousTermforEnergy     MultiFab to save viscous term in energy equation to
/// @param t                        Current time
/// @param mult_factor              Real factor to multiply
///                                 diffusion terms by before
///                                 adding to other sources.
///
    void add_viscous_term_to_source   (amrex::MultiFab& ext_src, amrex::MultiFab& state, amrex::MultiFab& ViscousTermforMomentum,
                                       amrex::MultiFab& ViscousTermforEnergy, amrex::Real t, amrex::Real mult_factor = 1.0);
#endif
#endif

#ifdef ROTATION

///
/// Calculate rotation source term and add to ``ext_src``
///
/// @param ext_src          MultiFab to add rotation term to
/// @param RotationTerm     MultiFab to save rotation term to
/// @param old_time         Old time
///
  void add_rotation_to_source(amrex::MultiFab& ext_src, amrex::MultiFab& RotationTerm, amrex::Real old_time);

///
/// Center rotation source term in time
///
/// @param S_new
/// @param OldRotationTerm
/// @param cur_time         Current time
/// @param dt               timestep
///
  void time_center_rotation(amrex::MultiFab& S_new, amrex::MultiFab& OldRotationTerm, amrex::Real cur_time, amrex::Real dt);
#endif


///
/// Volume weighted sum of given quantity
///
/// @param name         Name of quantity
/// @param time         current time
/// @param local        boolean, is sum local (over each patch) or over entire MultiFab?
/// @param finemask     boolean, should we build a mask to exclude finer levels?
///
    amrex::Real volWgtSum (const std::string& name, amrex::Real time, bool local=false, bool finemask=true);


///
/// Volume weight sum of (given quantity) squared
///
/// @param name     Name of quantity
/// @param time     current time
/// @param local    boolean, is sum local (over each patch) or over entire MultiFab?
///
    amrex::Real volWgtSquaredSum (const std::string& name, amrex::Real time, bool local=false);


///
/// Sum weighted by volume multiplied by distance from center in given direction
///
/// @param name     Name of quantity
/// @param time     current time
/// @param idir     Axis along which to compute distance from center
/// @param local    boolean, is sum local (over each patch) or over entire MultiFab?
///
    amrex::Real locWgtSum (const std::string& name, amrex::Real time, int idir, bool local=false);


///
/// Sum weighted by volume multiplied by product of distances from center in two given directions
///
/// @param name     Name of quantity
/// @param time     current time
/// @param idir1    First direction
/// @param idir2    Second direction
/// @param local    boolean, is sum local (over each patch) or over entire MultiFab?
///
    amrex::Real locWgtSum2D (const std::string& name, amrex::Real time, int idir1, int idir2, bool local=false);


///
/// Volume weighted sum of given component of multifab
///
/// @param mf       MultiFab to sum
/// @param comp     Component to sum
/// @param local    boolean, is sum local (over each patch) or over entire MultiFab?
///
    amrex::Real volWgtSumMF (const amrex::MultiFab& mf, int comp, bool local=false);


///
/// This function is a clone of ::locWgtSum except that it only sums over one half of the domain.
/// The lower half corresponds to side == 0, and the upper half corresponds to side == 1.
/// The argument idir (x == 0, y == 1, z == 2) gives the direction to location weight by,
/// and the argument bdir gives the direction along which to bisect.
///
/// @param name     Name of quantity
/// @param time     current time
/// @param idir     direction to weight by
/// @param side     lower half (0) or upper half (1)
/// @param bdir     Direction along which to bisect
/// @param local    boolean, is sum local (over each patch) or over entire MultiFab?
///
    amrex::Real locWgtSumOneSide (const std::string& name, amrex::Real time, int idir, int side, int bdir, bool local=false);


///
/// This function is a clone of ::volWgtSum except it computes the result only on half of the domain.
/// The lower half corresponds to side == 0 and the upper half corresponds to side == 1.
/// The argument bdir gives the direction along which to bisect.
///
/// @note Only designed to work in Cartesian geometries.
///
/// @param name     Name of quantity
/// @param time     current time
/// @param side     lower half (0) or upper half (1)
/// @param bdir     Direction along which to bisect
/// @param local    boolean, is sum local (over each patch) or over entire MultiFab?
///
    amrex::Real volWgtSumOneSide (const std::string& name, amrex::Real time, int side, int bdir, bool local=false);


///
/// Volume weighted sum of the product of two quantities
///
/// @param name1    Name of first quantity
/// @param name2    Name of second quantity
/// @param time     current time
/// @param local    boolean, is sum local (over each patch) or over entire MultiFab?
///
    amrex::Real volProductSum (const std::string& name1, const std::string& name2, amrex::Real time, bool local=false);


///
/// Location weighted sum of (quantity) squared
///
/// @param name     Name of quantity
/// @param time     current time
/// @param idir     Direction along which to weight
/// @param local    boolean, is sum local (over each patch) or over entire MultiFab?
///
    amrex::Real locSquaredSum (const std::string& name, amrex::Real time, int idir, bool local=false);

#ifdef GRAVITY
///
/// Are we using point mass gravity?
///
    int using_point_mass ();

///
/// Return the value of the point mass.
///
    amrex::Real get_point_mass ();
#endif

#ifdef HYBRID_MOMENTUM

///
/// Construct hybrid source terms at old time
///
/// @param source   MultiFab to save source to
/// @param state    Old state
/// @param time     current time
/// @param dt       timestep
///
    void construct_old_hybrid_source(amrex::MultiFab& source, amrex::MultiFab& state, amrex::Real time, amrex::Real dt);


///
/// Construct hybrid source terms at new time
///
/// @param source       MultiFab to save source to
/// @param state_old    Old state
/// @param state_new    New state
/// @param time         current time
/// @param dt           timestep
///
    void construct_new_hybrid_source(amrex::MultiFab& source, amrex::MultiFab& state_old, amrex::MultiFab& state_new, amrex::Real time, amrex::Real dt);


///
/// Fill ``source`` with hybrid source terms
///
/// @param state        Current state
/// @param source       MultiFab to save source to
/// @param mult_factor
///
    void fill_hybrid_hydro_source(amrex::MultiFab& state, amrex::MultiFab& source, const amrex::Real mult_factor);


///
/// Synchronize linear momentum with hybrid momentum
///
/// @param state    Current state
/// @param ng       Number of ghost cells
///
    void hybrid_to_linear_momentum(amrex::MultiFab& state, int ng = 0);


///
/// Synchronize hybrid momentum with linear momentum
///
/// @param state    Current state
/// @param ng       Number of ghost cells
///
    void linear_to_hybrid_momentum(amrex::MultiFab& state, int ng = 0);
#endif


///
/// @param force
/// @param sources
/// @param state
///
    void add_force_to_sources (amrex::MultiFab& force, amrex::MultiFab& sources, amrex::MultiFab& state);


///
/// Add source terms to ``state``
///
/// @param state    State
/// @param source   Source terms
/// @param dt       timestep
/// @param ng       number of ghost cells
///
    void apply_source_to_state (amrex::MultiFab& state, amrex::MultiFab& source, amrex::Real dt, int ng);


///
/// Fill a version of the state with ``ng`` ghost zones from the state data.
/// After we do the fill, optionally clean the data using the iclean parameter.
///
/// @param S        MultiFab to be filled
/// @param time     current time
/// @param ng       number of ghost cells
///
    void expand_state(amrex::MultiFab& S, amrex::Real time, int ng);

#ifdef SELF_GRAVITY

///
/// Calculate state in radial direction by averaging
///
/// @param is_new   do we use new (0) or old (1) data?
///
    void make_radial_data (int is_new);
#endif

#ifdef AUX_UPDATE

///
/// Calculate auxiliary variables at new timestep
///
/// @param time     current time
/// @param dt       timestep
///
    void advance_aux(amrex::Real time, amrex::Real dt);
#endif

#ifndef AMREX_USE_CUDA
    void do_sdc_update(int m1, int m2, amrex::Real dt);
#endif

#ifdef REACTIONS
/// Take an input conserved state U_state, evaluate the instantaneous
/// reaction rates, and store the result in R_source/
///
/// Special consideration is needed for fourth_order accuracy:
///
///  * If we come in with U_state representing an average
///    (input_is_average = true), then, for fourth_order, we will
///    convert U_state to cell-centers (consuming one ghost cell),
///    burn on the cell-centers (including one ghost cell) and then
///    average back to cell-centers (consuming another ghosst cell).
///    This requires >= 2 ghost cells for U_state, which we assume
///    are valid.
///
///  * If we come in with U_state representing cell centers
///    (input_is_average = false), then we just return R_source
///    at the cell-centers.
///
    void construct_old_react_source(amrex::MultiFab& U_state,
                                    amrex::MultiFab& R_source,
                                    const bool input_is_average);
#endif

// Hydrodynamics
#include <Castro_hydro.H>

#ifdef AMREX_PARTICLES
    static amrex::AmrTracerParticleContainer* theTracerPC () { return TracerPC; }
#endif

    static int       NUM_STATE;

    /// conservative variable index keys
    static int       Density, Xmom, Ymom, Zmom, Eden, Eint, Temp;
#ifdef HYBRID_MOMENTUM
    static int       Rmom, Lmom, Pmom;
#endif
#ifdef SHOCK_VAR
    static int       Shock;
#endif

    static int       FirstAdv,  NumAdv;
    static int       FirstSpec, NumSpec;
    static int       FirstAux,  NumAux;

    /// primitive variable index keys
    static int       QRHO, QU, QV, QW;
    static int       QGAME, QGC, QPRES, QREINT, QTEMP;
    static int       QFA, QFS, QFX;
#ifdef MHD
    static int       QMAGX, QMAGY, QMAGZ;
#endif
#ifdef RADIATION
    static int       QPTOT, QREITOT, QRAD;
#endif

    // Godunov state index keys
    static int GDRHO, GDU, GDV, GDW;
    static int GDPRES, GDGAME;
#ifdef RADIATION
    static int GDLAMS, GDERADS;
#endif

///
/// variable counts
///
    static int       NQSRC;
    static int       NQAUX;
    static int       NQ;
    static int       NGDNV;


    static amrex::Vector<std::string> source_names;
#ifdef THORNADO
    static amrex::Vector<std::string> thornado_plotvar_names;
#endif


///
/// This MultiFab is on the coarser level.  This is useful for the coarser level
///     to mask out the finer level.  We only build this when it is needed.
///     This coarse MultiFab has to live on the fine level because it must be updated
///     even when only the fine level changes.
///
    amrex::MultiFab fine_mask;
    amrex::MultiFab& build_fine_mask();


///
/// A record of how many cells we have advanced throughout the simulation.
///     This is saved as a real because we will be storing the number of zones
///     advanced as a ratio with the number of zones on the coarse grid (to
///     help prevent the number from getting too large), and that may be a
///     non-integer number.
///
    static amrex::Real num_zones_advanced;

#ifdef THORNADO
    static int THORNADO_NMOMENTS;
    static int THORNADO_RAD_NDOF;

    static int         thornado_nSpecies;

    // These are the lower and upper limits of the energy groups
    static int         thornado_swE;
    static amrex::Real thornado_eL;
    static amrex::Real thornado_eR;
    static int         thornado_ndimse;
    static amrex::Real thornado_zoome;
    static int         thornado_nplotvar;
    static int         thornado_plot_node_averages;
#endif

#ifdef AMREX_USE_CUDA
///
/// Minimum CUDA threadblock size for BC fills.
///
    static int numBCThreadsMin[3];

    static int minBCThreads(int dir) { return numBCThreadsMin[dir]; }
#endif

protected:


///
/// Build a mask that ghost cells overlapping with interior cells in the same multifab
///     are set to 0, whereas others are set to 1.
///
    amrex::Vector<std::unique_ptr<amrex::iMultiFab> > ib_mask;

///
/// @param ng number of ghost cells
///
    amrex::iMultiFab& build_interior_boundary_mask (int ng);


#ifdef THORNADO
    //
    // Read thornado-related inputs
    //
    static void read_thornado_params();

    std::unique_ptr<amrex::MultiFab> get_thornado_plotMF ();
#endif

#ifdef SELF_GRAVITY
    int get_numpts();


///
/// Difference between composite and level solve for phi.
///
    amrex::MultiFab comp_minus_level_phi;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > comp_minus_level_grad_phi;
#endif


///
/// A state array with ghost zones.
///
    amrex::MultiFab Sborder;

#ifdef RADIATION
    amrex::MultiFab Erborder;
    amrex::MultiFab lamborder;
#endif


///
/// A state array for the some reaction stuff with SDC
///
    amrex::MultiFab Sburn;


///
/// The primitive variable state array.
///
    amrex::MultiFab q;


///
/// we need a second version when doing fourth order
///
    amrex::MultiFab q_bar;


///
/// for diffusion in fourth order we need the cell-center temperature
//
#ifdef DIFFUSION
    amrex::MultiFab T_cc;
#endif

///
/// The auxiliary primitive variable state array.
///
    amrex::MultiFab qaux;
    amrex::MultiFab qaux_bar;


///
/// The source terms in primitive form.
///
    amrex::MultiFab src_q;


///
/// Source terms to the hydrodynamics solve.
///
    amrex::MultiFab sources_for_hydro;


///
/// Source term representing hydrodynamics update - this should just be
///     \f$-\nabla\cdot \{F(U)\}\f$ - no source terms
///
    amrex::MultiFab hydro_source;


///
/// Hydrodynamic (and radiation) fluxes.
///
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > fluxes;
#if (BL_SPACEDIM <= 2)
    amrex::MultiFab         P_radial;
#endif
#ifdef RADIATION
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > rad_fluxes;
#endif

    amrex::Vector<std::unique_ptr<amrex::MultiFab> > mass_fluxes;

    amrex::FluxRegister flux_reg;
#if (BL_SPACEDIM <= 2)
    amrex::FluxRegister pres_reg;
#endif
#ifdef RADIATION
    amrex::FluxRegister rad_flux_reg;
#endif
#ifdef GRAVITY
    amrex::FluxRegister phi_reg;
#endif

    // Does the momentum flux include the pressure term?
    // First index is the coordinate direction, second
    // index is the momentum direction.
    amrex::Array<amrex::Array<bool, 3>, 3> mom_flux_has_p;

///
/// Scalings for the flux registers.
///
    amrex::Real flux_crse_scale;
    amrex::Real flux_fine_scale;

#if (BL_SPACEDIM <= 2)
    amrex::Real pres_crse_scale;
    amrex::Real pres_fine_scale;
#endif


///
/// Did we trigger a CFL violation?
///
    int cfl_violation;


///
/// Did we have a successful burn?
///
    int burn_success;


///
/// State data to hold if we want to do a retry.
///
    amrex::Vector<std::unique_ptr<amrex::StateData> > prev_state;


///
/// Save the simulation time of the prev_state.
///
    amrex::Real prev_state_old_time;
    amrex::Real prev_state_new_time;
    bool prev_state_had_old_data;


///
/// Flag for indicating that we want to save prev_state until the reflux.
///
    bool keep_prev_state;



    //
    // Storage for the SDC time integration

    // this is the new iterations solution at the time nodes
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > k_new;

    // this is the old value of the advective update at the
    // nodes of the time integration
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > A_old;

    // this is the new value of the advective update at the
    // nodes of the time integration
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > A_new;

    // this is the old value of the reaction source at the
    // nodes of the time integration
#ifdef REACTIONS
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > R_old;
#endif

    static int SDC_NODES;
    static amrex::Vector<amrex::Real> dt_sdc;
    static amrex::Vector<amrex::Real> node_weights;


///
/// Wall time that we started the timestep
///
    amrex::Real wall_time_start;


///
/// Call extern/networks/.../network.f90::network_init()
///
    static void extern_init ();

    static void network_init ();

    static void network_finalize ();

    static void eos_finalize ();

    static void amrinfo_init();

    static void amrinfo_finalize();

    static void read_params ();


///
/// @param lev
///
    Castro& getLevel (int lev);

    void FluxRegCrseInit();
    void FluxRegFineAdd();


///
/// Do refluxing
///
/// @param crse_level
/// @param fine_level
///
    void reflux (int crse_level, int fine_level);


///
/// Normalize species fractions so they sum to 1
///
/// @param S_new    State
/// @param ng       number of ghost cells
///
    void normalize_species (amrex::MultiFab& S_new, int ng);


///
/// Enforces
/// \f[
/// \rho E = \rho e + \frac{1}{2} \rho (u^2 + v^2 + w^2)
/// \f]
///
/// @param S  State
///
    void enforce_consistent_e (amrex::MultiFab& S);

///
/// This routine sets the density in ``state`` to be larger than the density floor.
/// Note that it will operate everywhere on ``state``, including ghost zones.
///
/// The return value is the the negative fractional change in the state that has the
/// largest magnitude. If there is no reference state, this is meaningless.
///
/// @param state        state
/// @param ng           number of ghost cells
///
    amrex::Real enforce_min_density (amrex::MultiFab& state, int ng);

///
/// Given ``State_Type`` state data, perform a number of cleaning steps to make
/// sure the data is sensible. The return value is the same as the return
/// value of ::enforce_min_density.
///
/// @param state    State data
/// @param time     current time
/// @param ng       number of ghost cells
///
    amrex::Real clean_state (amrex::MultiFab& state, amrex::Real time, int ng);

///
/// Average new state from ``level+1`` down to ``level``
///
    void avgDown ();


///
/// Average given component of new state from ``level+1`` down to ``level``
///
/// @param state_indx   Index of state to average down
///
    void avgDown (int state_indx);

#ifdef THORNADO
    void average_down_thornado_data(const amrex::MultiFab& S_fine, amrex::MultiFab& S_crse, 
                                    int ncomp, const amrex::IntVect& ratio);
#endif

    void buildMetrics ();


///
/// Initialize the MultiFabs and flux registers that live as class members.
///
    void initMFs ();

///
/// integrate derived quantities over domain
///
/// @param name     Name of quantity to integrate
/// @param time     current time
/// @param local
///
    amrex::Real sumDerive (const std::string& name, amrex::Real time, bool local=false);

///
/// Calculates volume weight sums of variables and prints to screen
///
    void sum_integrated_quantities ();

    void write_info ();

#ifdef SELF_GRAVITY

///
/// Find the position of the "center" by interpolating from data at cell centers
///
/// @param S        Current state
/// @param time     current time
///
    void define_new_center (amrex::MultiFab& S, amrex::Real time);
    void write_center ();
#endif

#ifdef GRAVITY
///
/// @param time     current time
/// @param dt       timestep
///
    void pointmass_update(amrex::Real time, amrex::Real dt);
#endif

    void stopJob ();

///
/// Returns true if we're on the GPU and the total
/// memory on this level oversubscribes GPU memory.
///
    bool oversubscribing() {
#ifdef AMREX_USE_GPU
        return (amrex::MultiFab::queryMemUsage("AmrLevel_Level_" + std::to_string(level)) >= amrex::Gpu::Device::totalGlobalMem());
#else
        return false;
#endif
    }

///
/// The data.
///
    amrex::MultiFab             volume;
    amrex::MultiFab             area[3];
    amrex::MultiFab             dLogArea[1];
    amrex::Vector< amrex::Vector<amrex::Real> > radius;


///
/// Keep track of which AMR iteration we're on.
///
    int iteration;

    //
    // Static data members.
    //
#include <castro_params.H>

#ifdef AMREX_PARTICLES
#include "particles_params.H"
#endif

    static bool      signalStopJob;
    static int       radius_grow;
    static std::vector<std::string> err_list_names;
    static std::vector<int> err_list_ng;
    static int              num_err_list_default;
    static amrex::BCRec     phys_bc;
    static int       NUM_GROW;

    static int         lastDtPlotLimited;
    static amrex::Real lastDtBeforePlotLimiting;

    int lastDtRetryLimited;
    amrex::Real lastDtFromRetry;


///
/// for keeping track of mass changes from negative density resets
///
    static amrex::Real      frac_change;

///
/// For keeping track of fluid quantities lost at physical grid boundaries.
/// This should persist through restarts, but right now only on level 0.
///
    static const int n_lost = 8;
    amrex::Real             material_lost_through_boundary_cumulative[n_lost];
    amrex::Real             material_lost_through_boundary_temp[n_lost];


///
/// for keeping track of the amount of CPU time used -- this will persist
///     after restarts
///
    static amrex::Real      previousCPUTimeUsed;
    static amrex::Real      startCPUTime;

///
/// Get total CPU time
///
    static amrex::Real getCPUTime();

    bool             FillPatchedOldState_ok;


///
/// There can be only one Gravity object, it covers all levels:
///
    static class Gravity *gravity;


///
/// There can be only one Diffusion object, it covers all levels:
///
    static class Diffusion *diffusion;

#ifdef RADIATION

///
/// permits radiation to be turned on and off without recompiling:
///
    static int do_radiation;


///
/// There can be only one Radiation object, it covers all levels:
///
    static class Radiation *radiation;

    friend class Radiation;
#endif

#ifdef AMREX_PARTICLES
    static amrex::AmrTracerParticleContainer* TracerPC;
#endif

///
/// Name of the probin file and its length.
///
    static std::string probin_file;

    static amrex::IntVect hydro_tile_size;
    static amrex::IntVect no_tile_size;

    static int Knapsack_Weight_Type;
    static int SDC_Source_Type;
    static int num_state_type;


    // counters for various retries in Castro

///
/// subcycle capability
///
    int sub_iteration;
    int sub_ncycle;
    bool do_subcycle;
    amrex::Real dt_subcycle;
    amrex::Real dt_advance;


///
/// sdc
///
    int sdc_iteration;
    int current_sdc_node;



/* problem-specific includes */
#include <Problem.H>

};

//
// Inlines.
//

inline
int
Castro::numGrow()
{
    return NUM_GROW;
}

inline
int
Castro::numState()
{
    return NUM_STATE;
}

inline
amrex::MultiFab*
Castro::Area ()
{
    return area;
}


///
/// @param dir
///
inline
amrex::MultiFab&
Castro::Area (int dir)
{
    return area[dir];
}

inline
amrex::MultiFab&
Castro::Volume ()
{
    return volume;
}


///
/// @param lev
///
inline
Castro&
Castro::getLevel (int lev)
{
    return *(Castro *) &parent->getLevel(lev);
}

inline
void
Castro::stopJob()
{
  signalStopJob = true;


}

inline
int
Castro::get_output_at_completion()
{
    return output_at_completion;
}

#ifdef GRAVITY
inline
amrex::Real
Castro::get_point_mass ()
{
    return point_mass;
}

inline
int
Castro::using_point_mass ()
{
    return use_point_mass;
}
#endif

#endif /*_Castro_H_*/
