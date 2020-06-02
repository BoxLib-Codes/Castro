#include <cstdio>

#include <AMReX_ParmParse.H>
#include "AMReX_LevelBld.H"
#include "Castro.H"
#include "Castro_F.H"
#include "Castro_bc_fill_nd.H"
#include "Castro_bc_fill_nd_F.H"
#include "Castro_generic_fill.H"
#include "Derive.H"
#include "eos.H"
#ifdef RADIATION
#include "RAD_F.H"
#include "Radiation.H"
#endif
#include <Problem_Derive_F.H>

#include "AMReX_buildInfo.H"

using std::string;
using namespace amrex;

static Box the_same_box(const Box& b) { return b; }
static Box grow_box_by_one(const Box& b) { return amrex::grow(b, 1); }

typedef StateDescriptor::BndryFunc BndryFunc;

//
// Components are:
//  Interior, Inflow, Outflow,  Symmetry,     SlipWall,     NoSlipWall
//
static int scalar_bc[] = {INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN};

static int norm_vel_bc[] = {INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_ODD, REFLECT_ODD, REFLECT_ODD};

static int tang_vel_bc[] = {INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN};

#ifdef MHD
static int mag_field_bc[] = {INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, FOEXTRAP, HOEXTRAP};
#endif

static void set_scalar_bc(BCRec& bc, const BCRec& phys_bc) {
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < BL_SPACEDIM; i++) {
        bc.setLo(i, scalar_bc[lo_bc[i]]);
        bc.setHi(i, scalar_bc[hi_bc[i]]);
    }
}

static void set_x_vel_bc(BCRec& bc, const BCRec& phys_bc) {
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0, norm_vel_bc[lo_bc[0]]);
    bc.setHi(0, norm_vel_bc[hi_bc[0]]);
#if (BL_SPACEDIM >= 2)
    bc.setLo(1, tang_vel_bc[lo_bc[1]]);
    bc.setHi(1, tang_vel_bc[hi_bc[1]]);
#endif
#if (BL_SPACEDIM == 3)
    bc.setLo(2, tang_vel_bc[lo_bc[2]]);
    bc.setHi(2, tang_vel_bc[hi_bc[2]]);
#endif
}

static void set_y_vel_bc(BCRec& bc, const BCRec& phys_bc) {
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0, tang_vel_bc[lo_bc[0]]);
    bc.setHi(0, tang_vel_bc[hi_bc[0]]);
#if (BL_SPACEDIM >= 2)
    bc.setLo(1, norm_vel_bc[lo_bc[1]]);
    bc.setHi(1, norm_vel_bc[hi_bc[1]]);
#endif
#if (BL_SPACEDIM == 3)
    bc.setLo(2, tang_vel_bc[lo_bc[2]]);
    bc.setHi(2, tang_vel_bc[hi_bc[2]]);
#endif
}

static void set_z_vel_bc(BCRec& bc, const BCRec& phys_bc) {
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0, tang_vel_bc[lo_bc[0]]);
    bc.setHi(0, tang_vel_bc[hi_bc[0]]);
#if (BL_SPACEDIM >= 2)
    bc.setLo(1, tang_vel_bc[lo_bc[1]]);
    bc.setHi(1, tang_vel_bc[hi_bc[1]]);
#endif
#if (BL_SPACEDIM == 3)
    bc.setLo(2, norm_vel_bc[lo_bc[2]]);
    bc.setHi(2, norm_vel_bc[hi_bc[2]]);
#endif
}

#ifdef MHD
static void set_mag_field_bc(BCRec& bc, const BCRec& phys_bc) {
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < BL_SPACEDIM; i++) {
        bc.setLo(i, mag_field_bc[lo_bc[i]]);
        bc.setHi(i, mag_field_bc[hi_bc[i]]);
    }
}
#endif

// In some cases we want to replace inflow boundaries with
// first-order extrapolation boundaries. This is intended to
// be used for state data that the user is not going to
// provide inflow boundary conditions for, like gravity
// and reactions, and it works in conjunction with the
// generic_fill boundary routine.

static void replace_inflow_bc(BCRec& bc) {
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (bc.lo(dir) == EXT_DIR) {
            bc.setLo(dir, FOEXTRAP);
        }
        if (bc.hi(dir) == EXT_DIR) {
            bc.setHi(dir, FOEXTRAP);
        }
    }
}

void Castro::variableSetUp() {

    // Castro::variableSetUp is called in the constructor of Amr.cpp, so
    // it should get called every time we start or restart a job

    // initialize the start time for our CPU-time tracker
    startCPUTime = ParallelDescriptor::second();

    // Output the git commit hashes used to build the executable.

    if (ParallelDescriptor::IOProcessor()) {

        const char* castro_hash = buildInfoGetGitHash(1);
        const char* amrex_hash = buildInfoGetGitHash(2);
        const char* microphysics_hash = buildInfoGetGitHash(3);
        const char* buildgithash = buildInfoGetBuildGitHash();
        const char* buildgitname = buildInfoGetBuildGitName();

        if (strlen(castro_hash) > 0) {
            std::cout << "\n"
                      << "Castro git describe: " << castro_hash << "\n";
        }
        if (strlen(amrex_hash) > 0) {
            std::cout << "AMReX git describe: " << amrex_hash << "\n";
        }
        if (strlen(microphysics_hash) > 0) {
            std::cout << "Microphysics git describe: " << microphysics_hash << "\n";
        }
        if (strlen(buildgithash) > 0) {
            std::cout << buildgitname << " git describe: " << buildgithash << "\n";
        }

        std::cout << "\n";
    }

    BL_ASSERT(desc_lst.size() == 0);

    // read the C++ parameters that are set in inputs and do other
    // initializations (e.g., set phys_bc)
    read_params();

    // Read in the input values to Fortran.
    ca_set_castro_method_params();

    // Initialize the runtime parameters for any of the external
    // microphysics (these are the parameters that are in the &extern
    // block of the probin file)
    extern_init();

    // set small positive values of the "small" quantities if they are
    // negative this mirrors the logic in ca_set_method_params for
    // Fortran
    if (small_dens < 0.0_rt) {
        small_dens = 1.e-200_rt;
    }

    if (small_temp < 0.0_rt) {
        small_temp = 1.e-200_rt;
    }

    if (small_pres < 0.0_rt) {
        small_pres = 1.e-200_rt;
    }

    if (small_ener < 0.0_rt) {
        small_ener = 1.e-200_rt;
    }

    // Initialize the network
    ca_network_init();
#ifdef CXX_REACTIONS
    network_init();
#endif

    // Initialize the EOS
    ca_eos_init();
    eos_init();

    // Ensure that Castro's small variables are consistent
    // with the minimum permitted by the EOS, and vice versa.

    Real new_min_T = std::max(small_temp, EOSData::mintemp);
    small_temp = new_min_T;
    EOSData::mintemp = new_min_T;

    Real new_min_rho = std::max(small_dens, EOSData::mindens);
    small_dens = new_min_rho;
    EOSData::mindens = new_min_rho;

    // some consistency checks on the parameters
#ifdef REACTIONS
    int abort_on_failure;
    ca_get_abort_on_failure(&abort_on_failure);

#ifdef TRUE_SDC
    // for TRUE_SDC, we don't support retry, so we need to ensure that abort_on_failure = T
    if (use_retry) {
        amrex::Warning("use_retry = 1 is not supported with true SDC.  Disabling");
        use_retry = 0;
    }
    if (!abort_on_failure) {
        amrex::Warning("abort_on_failure = F not supported with true SDC.  Resetting");
        abort_on_failure = 1;
        ca_set_abort_on_failure(&abort_on_failure);
    }
#else
    if (!use_retry && !abort_on_failure) {
        amrex::Error("use_retry = 0 and abort_on_failure = F is dangerous and not supported");
    }
#endif
#endif

#ifdef REACTIONS
    // Initialize the burner
    burner_init();
#endif

    // Initialize the amr info
    amrinfo_init();

    const int dm = BL_SPACEDIM;

    // NUM_GROW is the number of ghost cells needed for the hyperbolic portions
#ifdef MHD
    NUM_GROW = 6;
#else
    NUM_GROW = 4;
#endif

    const Real run_strt = ParallelDescriptor::second();

    // set the conserved, primitive, aux, and godunov indices in Fortran
    ca_set_method_params(dm);

    // setup the passive maps -- this follows the same logic as the
    // Fortran versions in ca_set_method_params
    int ipassive = 0;

    upass_map.resize(npassive);
    qpass_map.resize(npassive);

    for (int iadv = 0; iadv < NumAdv; ++iadv) {
        upass_map[ipassive] = UFA + iadv;
        qpass_map[ipassive] = QFA + iadv;
        ++ipassive;
    }

    for (int ispec = 0; ispec < NumSpec; ++ispec) {
        upass_map[ipassive] = UFS + ispec;
        qpass_map[ipassive] = QFS + ispec;
        ++ipassive;
    }

    for (int iaux = 0; iaux < NumAux; ++iaux) {
        upass_map[ipassive] = UFX + iaux;
        qpass_map[ipassive] = QFX + iaux;
        ++ipassive;
    }

    Real run_stop = ParallelDescriptor::second() - run_strt;

    ParallelDescriptor::ReduceRealMax(run_stop, ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "\nTime in ca_set_method_params: " << run_stop << '\n';
    }

    const Geometry& dgeom = DefaultGeometry();

    const int coord_type = dgeom.Coord();

    // Get the center variable from the inputs and pass it directly to Fortran.
    Vector<Real> center(BL_SPACEDIM, 0.0);
    ParmParse ppc("castro");
    ppc.queryarr("center", center, 0, BL_SPACEDIM);

    ca_set_problem_params(dm, phys_bc.lo(), phys_bc.hi(), Interior, Inflow, Outflow, Symmetry,
                          SlipWall, NoSlipWall, coord_type, dgeom.ProbLo(), dgeom.ProbHi(),
                          center.dataPtr());

    // Read in the parameters for the tagging criteria
    // and store them in the Fortran module.

    const int probin_file_length = probin_file.length();
    Vector<int> probin_file_name(probin_file_length);

    for (int i = 0; i < probin_file_length; i++) {
        probin_file_name[i] = probin_file[i];
    }

    ca_get_tagging_params(probin_file_name.dataPtr(), &probin_file_length);

#ifdef SPONGE
    // Initialize the sponge

    sponge_init();

    // Read in the parameters for the sponge
    // and store them in the Fortran module.

    ca_read_sponge_params(probin_file_name.dataPtr(), &probin_file_length);

    // bring the sponge parameters into C++
    ca_get_sponge_params(sponge_lower_factor, sponge_upper_factor, sponge_lower_radius,
                         sponge_upper_radius, sponge_lower_density, sponge_upper_density,
                         sponge_lower_pressure, sponge_upper_pressure, sponge_target_velocity,
                         sponge_timescale);

#endif

    // Read in the ambient state parameters.

    ca_get_ambient_params(probin_file_name.dataPtr(), &probin_file_length);

    Interpolater* interp;

    if (state_interp_order == 0) {
        interp = &pc_interp;
    } else {
        if (lin_limit_state_interp == 1) {
            interp = &lincc_interp;
        } else {
            interp = &cell_cons_interp;
        }
    }

    // Note that the default is state_data_extrap = false,
    // store_in_checkpoint = true.  We only need to put these in
    // explicitly if we want to do something different,
    // like not store the state data in a checkpoint directory
    bool state_data_extrap = false;
    bool store_in_checkpoint;

#if defined(RADIATION)
    // Radiation should always have at least one ghost zone.
    int ngrow_state = std::max(1, state_nghost);
#else
    int ngrow_state = state_nghost;
#endif

    BL_ASSERT(ngrow_state >= 0);

    store_in_checkpoint = true;
    desc_lst.addDescriptor(State_Type, IndexType::TheCellType(), StateDescriptor::Point,
                           ngrow_state, NUM_STATE, interp, state_data_extrap, store_in_checkpoint);

#ifdef MHD
    store_in_checkpoint = true;
    IndexType xface(IntVect{AMREX_D_DECL(1, 0, 0)});
    desc_lst.addDescriptor(Mag_Type_x, xface, StateDescriptor::Point, 0, 1, interp,
                           state_data_extrap, store_in_checkpoint);
    IndexType yface(IntVect{AMREX_D_DECL(0, 1, 0)});
    desc_lst.addDescriptor(Mag_Type_y, yface, StateDescriptor::Point, 0, 1, interp,
                           state_data_extrap, store_in_checkpoint);
    IndexType zface(IntVect{AMREX_D_DECL(0, 0, 1)});
    desc_lst.addDescriptor(Mag_Type_z, zface, StateDescriptor::Point, 0, 1, interp,
                           state_data_extrap, store_in_checkpoint);
#endif

#ifdef GRAVITY
    store_in_checkpoint = true;
    desc_lst.addDescriptor(PhiGrav_Type, IndexType::TheCellType(), StateDescriptor::Point, 1, 1,
                           &cell_cons_interp, state_data_extrap, store_in_checkpoint);

    store_in_checkpoint = false;
    desc_lst.addDescriptor(Gravity_Type, IndexType::TheCellType(), StateDescriptor::Point, NUM_GROW,
                           3, &cell_cons_interp, state_data_extrap, store_in_checkpoint);
#endif

    // Source terms -- for the CTU method, because we do characteristic
    // tracing on the source terms, we need NUM_GROW ghost cells to do
    // the reconstruction.  For SDC, on the other hand, we only
    // need 1 (for the fourth-order stuff). Simplified SDC uses the CTU
    // advance, so it behaves the same way as CTU here.

    store_in_checkpoint = true;
    int source_ng;
    if (time_integration_method == CornerTransportUpwind ||
        time_integration_method == SimplifiedSpectralDeferredCorrections) {
        source_ng = NUM_GROW;
    } else if (time_integration_method == SpectralDeferredCorrections) {
        if (sdc_order == 2 && use_pslope) {
            source_ng = NUM_GROW;
        } else {
            source_ng = 1;
        }
    } else {
        amrex::Error("Unknown time_integration_method");
    }
    desc_lst.addDescriptor(Source_Type, IndexType::TheCellType(), StateDescriptor::Point, source_ng,
                           NSRC, &cell_cons_interp, state_data_extrap, store_in_checkpoint);

#ifdef ROTATION
    store_in_checkpoint = false;
    desc_lst.addDescriptor(PhiRot_Type, IndexType::TheCellType(), StateDescriptor::Point, 1, 1,
                           &cell_cons_interp, state_data_extrap, store_in_checkpoint);

    store_in_checkpoint = false;
    desc_lst.addDescriptor(Rotation_Type, IndexType::TheCellType(), StateDescriptor::Point,
                           NUM_GROW, 3, &cell_cons_interp, state_data_extrap, store_in_checkpoint);
#endif

#ifdef REACTIONS
    // Components 0:Numspec-1         are      omegadot_i
    // Component    NumSpec            is      enuc =      (eout-ein)
    // Component    NumSpec+1          is  rho_enuc= rho * (eout-ein)
    // Component    NumSpec+2          is      weights ~ number of RHS calls
    store_in_checkpoint = true;
    desc_lst.addDescriptor(Reactions_Type, IndexType::TheCellType(), StateDescriptor::Point,
                           NUM_GROW, NumSpec + 3, &cell_cons_interp, state_data_extrap,
                           store_in_checkpoint);
#endif

#ifdef SIMPLIFIED_SDC
#ifdef REACTIONS
    // For simplified SDC, we want to store the reactions source.

    if (time_integration_method == SimplifiedSpectralDeferredCorrections) {

        store_in_checkpoint = true;
        desc_lst.addDescriptor(Simplified_SDC_React_Type, IndexType::TheCellType(),
                               StateDescriptor::Point, NUM_GROW, NQSRC, &cell_cons_interp,
                               state_data_extrap, store_in_checkpoint);
    }
#endif
#endif

    Vector<BCRec> bcs(NUM_STATE);
    Vector<std::string> name(NUM_STATE);

    BCRec bc;
    set_scalar_bc(bc, phys_bc);
    bcs[URHO] = bc;
    name[URHO] = "density";

    set_x_vel_bc(bc, phys_bc);
    bcs[UMX] = bc;
    name[UMX] = "xmom";

    set_y_vel_bc(bc, phys_bc);
    bcs[UMY] = bc;
    name[UMY] = "ymom";

    set_z_vel_bc(bc, phys_bc);
    bcs[UMZ] = bc;
    name[UMZ] = "zmom";

#ifdef HYBRID_MOMENTUM
    set_scalar_bc(bc, phys_bc);
    bcs[UMR] = bc;
    name[UMR] = "rmom";

    set_scalar_bc(bc, phys_bc);
    bcs[UML] = bc;
    name[UML] = "lmom";

    set_scalar_bc(bc, phys_bc);
    bcs[UMP] = bc;
    name[UMP] = "pmom";

#endif
    set_scalar_bc(bc, phys_bc);
    bcs[UEDEN] = bc;
    name[UEDEN] = "rho_E";

    set_scalar_bc(bc, phys_bc);
    bcs[UEINT] = bc;
    name[UEINT] = "rho_e";

    set_scalar_bc(bc, phys_bc);
    bcs[UTEMP] = bc;
    name[UTEMP] = "Temp";

    for (int i = 0; i < NumAdv; ++i) {
        char buf[64];
        sprintf(buf, "adv_%d", i);
        set_scalar_bc(bc, phys_bc);
        bcs[UFA + i] = bc;
        name[UFA + i] = string(buf);
    }

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << NumSpec << " Species: " << std::endl;
        for (int i = 0; i < NumSpec; i++) {
            std::cout << short_spec_names_cxx[i] << ' ' << ' ';
        }
        std::cout << std::endl;
    }

    for (int i = 0; i < NumSpec; ++i) {
        set_scalar_bc(bc, phys_bc);
        bcs[UFS + i] = bc;
        name[UFS + i] = "rho_" + short_spec_names_cxx[i];
    }

    // Get the auxiliary names from the network model.
    std::vector<std::string> aux_names;
    for (int i = 0; i < NumAux; i++) {
        aux_names.push_back(short_aux_names_cxx[i]);
    }

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << NumAux << " Auxiliary Variables: " << std::endl;
        for (int i = 0; i < NumAux; i++) {
            std::cout << aux_names[i] << ' ' << ' ';
        }
        std::cout << std::endl;
    }

    for (int i = 0; i < NumAux; ++i) {
        set_scalar_bc(bc, phys_bc);
        bcs[UFX + i] = bc;
        name[UFX + i] = "rho_" + aux_names[i];
    }

#ifdef SHOCK_VAR
    set_scalar_bc(bc, phys_bc);
    bcs[USHK] = bc;
    name[USHK] = "Shock";
#endif

    BndryFunc stateBndryFunc(ca_statefill);
    stateBndryFunc.setRunOnGPU(true);

    desc_lst.setComponent(State_Type, URHO, name, bcs, stateBndryFunc);

    BndryFunc genericBndryFunc(ca_generic_fill);
    genericBndryFunc.setRunOnGPU(true);

#ifdef MHD
    set_mag_field_bc(bc, phys_bc);
    desc_lst.setComponent(Mag_Type_x, 0, "b_x", bc, BndryFunc(ca_face_fillx));
    desc_lst.setComponent(Mag_Type_y, 0, "b_y", bc, BndryFunc(ca_face_filly));
    desc_lst.setComponent(Mag_Type_z, 0, "b_z", bc, BndryFunc(ca_face_fillz));
#endif

#ifdef GRAVITY
    set_scalar_bc(bc, phys_bc);
    replace_inflow_bc(bc);
    desc_lst.setComponent(PhiGrav_Type, 0, "phiGrav", bc, genericBndryFunc);
    set_x_vel_bc(bc, phys_bc);
    replace_inflow_bc(bc);
    desc_lst.setComponent(Gravity_Type, 0, "grav_x", bc, genericBndryFunc);
    set_y_vel_bc(bc, phys_bc);
    replace_inflow_bc(bc);
    desc_lst.setComponent(Gravity_Type, 1, "grav_y", bc, genericBndryFunc);
    set_z_vel_bc(bc, phys_bc);
    replace_inflow_bc(bc);
    desc_lst.setComponent(Gravity_Type, 2, "grav_z", bc, genericBndryFunc);
#endif

#ifdef ROTATION
    set_scalar_bc(bc, phys_bc);
    replace_inflow_bc(bc);
    desc_lst.setComponent(PhiRot_Type, 0, "phiRot", bc, genericBndryFunc);
    set_x_vel_bc(bc, phys_bc);
    replace_inflow_bc(bc);
    desc_lst.setComponent(Rotation_Type, 0, "rot_x", bc, genericBndryFunc);
    set_y_vel_bc(bc, phys_bc);
    replace_inflow_bc(bc);
    desc_lst.setComponent(Rotation_Type, 1, "rot_y", bc, genericBndryFunc);
    set_z_vel_bc(bc, phys_bc);
    replace_inflow_bc(bc);
    desc_lst.setComponent(Rotation_Type, 2, "rot_z", bc, genericBndryFunc);
#endif

    // Source term array will use source fill

    Vector<BCRec> source_bcs(NSRC);
    Vector<std::string> state_type_source_names(NSRC);

    for (int i = 0; i < NSRC; ++i) {
        state_type_source_names[i] = name[i] + "_source";
        source_bcs[i] = bcs[i];

        // Replace any instances of inflow boundaries with first-order extrapolation.
        replace_inflow_bc(source_bcs[i]);
    }

    desc_lst.setComponent(Source_Type, URHO, state_type_source_names, source_bcs, genericBndryFunc);

#ifdef REACTIONS
    std::string name_react;
    for (int i = 0; i < NumSpec; ++i) {
        set_scalar_bc(bc, phys_bc);
        replace_inflow_bc(bc);
        name_react = "omegadot_" + short_spec_names_cxx[i];
        desc_lst.setComponent(Reactions_Type, i, name_react, bc, genericBndryFunc);
    }
    desc_lst.setComponent(Reactions_Type, NumSpec, "enuc", bc, genericBndryFunc);
    desc_lst.setComponent(Reactions_Type, NumSpec + 1, "rho_enuc", bc, genericBndryFunc);
    desc_lst.setComponent(Reactions_Type, NumSpec + 2, "weights", bc, genericBndryFunc);
#endif

#ifdef SIMPLIFIED_SDC
#ifdef REACTIONS
    if (time_integration_method == SimplifiedSpectralDeferredCorrections) {
        for (int i = 0; i < NQSRC; ++i) {
            char buf[64];
            sprintf(buf, "sdc_react_source_%d", i);
            set_scalar_bc(bc, phys_bc);
            replace_inflow_bc(bc);

            desc_lst.setComponent(Simplified_SDC_React_Type, i, std::string(buf), bc,
                                  genericBndryFunc);
        }
    }
#endif
#endif

#ifdef RADIATION
    int ngrow = 1;
    int ncomp = Radiation::nGroups;
    desc_lst.addDescriptor(Rad_Type, IndexType::TheCellType(), StateDescriptor::Point, ngrow, ncomp,
                           interp);
    set_scalar_bc(bc, phys_bc);
    replace_inflow_bc(bc);

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "Radiation::nGroups = " << Radiation::nGroups << std::endl;
    }

    char rad_name[10];
    if (!Radiation::do_multigroup) {
        desc_lst.setComponent(Rad_Type, Rad, "rad", bc, genericBndryFunc);
    } else {
        for (int i = 0; i < Radiation::nGroups; i++) {
            sprintf(rad_name, "rad%d", i);
            desc_lst.setComponent(Rad_Type, i, rad_name, bc, genericBndryFunc);
        }
    }
#endif

    // some optional State_Type's -- since these depend on the value of
    // runtime parameters, we don't add these to the enum, but instead
    // add them to the count of State_Type's if we will use them

#ifdef REACTIONS
    if (time_integration_method == SpectralDeferredCorrections && sdc_order == 4) {

        // we are doing 4th order reactive SDC.  We need 2 ghost cells here
        SDC_Source_Type = desc_lst.size();

        store_in_checkpoint = false;
        desc_lst.addDescriptor(SDC_Source_Type, IndexType::TheCellType(), StateDescriptor::Point, 2,
                               NUM_STATE, interp, state_data_extrap, store_in_checkpoint);

        // this is the same thing we do for the sources, but now we use the generic fill
        Vector<BCRec> sdc_source_bcs(NUM_STATE);
        Vector<std::string> sdc_source_names(NUM_STATE);

        for (int i = 0; i < NUM_STATE; ++i) {
            sdc_source_names[i] = name[i] + "_sdc_source";
            sdc_source_bcs[i] = bcs[i];

            // Replace any instances of inflow boundaries with first-order extrapolation.
            replace_inflow_bc(sdc_source_bcs[i]);
        }

        desc_lst.setComponent(SDC_Source_Type, URHO, sdc_source_names, sdc_source_bcs,
                              genericBndryFunc);
    }
#endif

    num_state_type = desc_lst.size();

    //
    // DEFINE DERIVED QUANTITIES
    //
    // Pressure
    //
    derive_lst.add("pressure", IndexType::TheCellType(), 1, ca_derpres, the_same_box);
    derive_lst.addComponent("pressure", desc_lst, State_Type, URHO, NUM_STATE);

    //
    // Kinetic energy
    //
    derive_lst.add("kineng", IndexType::TheCellType(), 1, ca_derkineng, the_same_box);
    derive_lst.addComponent("kineng", desc_lst, State_Type, URHO, 1);
    derive_lst.addComponent("kineng", desc_lst, State_Type, UMX, 3);

    //
    // Sound speed (c)
    //
    derive_lst.add("soundspeed", IndexType::TheCellType(), 1, ca_dersoundspeed, the_same_box);
    derive_lst.addComponent("soundspeed", desc_lst, State_Type, URHO, NUM_STATE);

    //
    // Gamma_1
    //
    derive_lst.add("Gamma_1", IndexType::TheCellType(), 1, ca_dergamma1, the_same_box);
    derive_lst.addComponent("Gamma_1", desc_lst, State_Type, URHO, NUM_STATE);

    //
    // Mach number(M)
    //
    derive_lst.add("MachNumber", IndexType::TheCellType(), 1, ca_dermachnumber, the_same_box);
    derive_lst.addComponent("MachNumber", desc_lst, State_Type, URHO, NUM_STATE);

#if (BL_SPACEDIM == 1)
    //
    // Wave speed u+c
    //
    derive_lst.add("uplusc", IndexType::TheCellType(), 1, ca_deruplusc, the_same_box);
    derive_lst.addComponent("uplusc", desc_lst, State_Type, URHO, NUM_STATE);

    //
    // Wave speed u-c
    //
    derive_lst.add("uminusc", IndexType::TheCellType(), 1, ca_deruminusc, the_same_box);
    derive_lst.addComponent("uminusc", desc_lst, State_Type, URHO, NUM_STATE);
#endif

    //
    // Gravitational forcing
    //
#ifdef GRAVITY
    //    derive_lst.add("rhog",IndexType::TheCellType(),1,
    //                   BL_FORT_PROC_CALL(CA_RHOG,ca_rhog),the_same_box);
    //    derive_lst.addComponent("rhog",desc_lst,State_Type, URHO, 1);
    //    derive_lst.addComponent("rhog",desc_lst,Gravity_Type,0,BL_SPACEDIM);
#endif

    //
    // Entropy (S)
    //
    derive_lst.add("entropy", IndexType::TheCellType(), 1, ca_derentropy, the_same_box);
    derive_lst.addComponent("entropy", desc_lst, State_Type, URHO, NUM_STATE);

#ifdef DIFFUSION
    if (diffuse_temp) {
        //
        // thermal conductivity (k_th)
        //
        derive_lst.add("thermal_cond", IndexType::TheCellType(), 1, ca_dercond, the_same_box);
        derive_lst.addComponent("thermal_cond", desc_lst, State_Type, URHO, NUM_STATE);

        //
        // thermal diffusivity (k_th/(rho c_v))
        //
        derive_lst.add("diff_coeff", IndexType::TheCellType(), 1, ca_derdiffcoeff, the_same_box);
        derive_lst.addComponent("diff_coeff", desc_lst, State_Type, URHO, NUM_STATE);

        //
        // diffusion term (the divergence of thermal flux)
        //
        derive_lst.add("diff_term", IndexType::TheCellType(), 1, ca_derdiffterm, grow_box_by_one);
        derive_lst.addComponent("diff_term", desc_lst, State_Type, URHO, NUM_STATE);
    }
#endif

    //
    // Vorticity
    //
    derive_lst.add("magvort", IndexType::TheCellType(), 1, ca_dermagvort, grow_box_by_one);
    // Here we exploit the fact that UMX = URHO + 1
    //   in order to use the correct interpolation.
    if (UMX != URHO + 1) {
        amrex::Error("We are assuming UMX = URHO + 1 in Castro_setup.cpp");
    }
    derive_lst.addComponent("magvort", desc_lst, State_Type, URHO, 4);

    //
    // Div(u)
    //
    derive_lst.add("divu", IndexType::TheCellType(), 1, ca_derdivu, grow_box_by_one);
    derive_lst.addComponent("divu", desc_lst, State_Type, URHO, 1);
    derive_lst.addComponent("divu", desc_lst, State_Type, UMX, 3);

    //
    // Internal energy as derived from rho*E, part of the state
    //
    derive_lst.add("eint_E", IndexType::TheCellType(), 1, ca_dereint1, the_same_box);
    derive_lst.addComponent("eint_E", desc_lst, State_Type, URHO, NUM_STATE);

    //
    // Internal energy as derived from rho*e, part of the state
    //
    derive_lst.add("eint_e", IndexType::TheCellType(), 1, ca_dereint2, the_same_box);
    derive_lst.addComponent("eint_e", desc_lst, State_Type, URHO, NUM_STATE);

    //
    // Log(density)
    //
    derive_lst.add("logden", IndexType::TheCellType(), 1, ca_derlogden, the_same_box);
    derive_lst.addComponent("logden", desc_lst, State_Type, URHO, NUM_STATE);

    derive_lst.add("StateErr", IndexType::TheCellType(), 3, ca_derstate, grow_box_by_one);
    derive_lst.addComponent("StateErr", desc_lst, State_Type, URHO, 1);
    derive_lst.addComponent("StateErr", desc_lst, State_Type, UTEMP, 1);
    derive_lst.addComponent("StateErr", desc_lst, State_Type, UFS, 1);

    //
    // X from rhoX
    //
    for (int i = 0; i < NumSpec; i++) {
        std::string spec_string = "X(" + short_spec_names_cxx[i] + ")";
        derive_lst.add(spec_string, IndexType::TheCellType(), 1, ca_derspec, the_same_box);
        derive_lst.addComponent(spec_string, desc_lst, State_Type, URHO, 1);
        derive_lst.addComponent(spec_string, desc_lst, State_Type, UFS + i, 1);
    }

    //
    // Abar
    //
    derive_lst.add("abar", IndexType::TheCellType(), 1, ca_derabar, the_same_box);
    derive_lst.addComponent("abar", desc_lst, State_Type, URHO, 1);
    derive_lst.addComponent("abar", desc_lst, State_Type, UFS, NumSpec);

    //
    // Velocities
    //
    derive_lst.add("x_velocity", IndexType::TheCellType(), 1, ca_dervel, the_same_box);
    derive_lst.addComponent("x_velocity", desc_lst, State_Type, URHO, 1);
    derive_lst.addComponent("x_velocity", desc_lst, State_Type, UMX, 1);

    derive_lst.add("y_velocity", IndexType::TheCellType(), 1, ca_dervel, the_same_box);
    derive_lst.addComponent("y_velocity", desc_lst, State_Type, URHO, 1);
    derive_lst.addComponent("y_velocity", desc_lst, State_Type, UMY, 1);

    derive_lst.add("z_velocity", IndexType::TheCellType(), 1, ca_dervel, the_same_box);
    derive_lst.addComponent("z_velocity", desc_lst, State_Type, URHO, 1);
    derive_lst.addComponent("z_velocity", desc_lst, State_Type, UMZ, 1);

#ifdef REACTIONS
    //
    // Nuclear energy generation timescale t_e == e / edot
    // Sound-crossing time t_s == dx / c_s
    // Ratio of these is t_s_t_e == t_s / t_e
    //
    derive_lst.add("t_sound_t_enuc", IndexType::TheCellType(), 1, ca_derenuctimescale,
                   the_same_box);
    derive_lst.addComponent("t_sound_t_enuc", desc_lst, State_Type, URHO, NUM_STATE);
    derive_lst.addComponent("t_sound_t_enuc", desc_lst, Reactions_Type, NumSpec, 1);
#endif

    derive_lst.add("magvel", IndexType::TheCellType(), 1, ca_dermagvel, the_same_box);
    derive_lst.addComponent("magvel", desc_lst, State_Type, URHO, 1);
    derive_lst.addComponent("magvel", desc_lst, State_Type, UMX, 3);

    derive_lst.add("radvel", IndexType::TheCellType(), 1, ca_derradialvel, the_same_box);
    derive_lst.addComponent("radvel", desc_lst, State_Type, URHO, 1);
    derive_lst.addComponent("radvel", desc_lst, State_Type, UMX, 3);

    derive_lst.add("circvel", IndexType::TheCellType(), 1, ca_dercircvel, the_same_box);
    derive_lst.addComponent("circvel", desc_lst, State_Type, URHO, 1);
    derive_lst.addComponent("circvel", desc_lst, State_Type, UMX, 3);

    derive_lst.add("magmom", IndexType::TheCellType(), 1, ca_dermagmom, the_same_box);
    derive_lst.addComponent("magmom", desc_lst, State_Type, UMX, 3);

    derive_lst.add("angular_momentum_x", IndexType::TheCellType(), 1, ca_derangmomx, the_same_box);
    derive_lst.addComponent("angular_momentum_x", desc_lst, State_Type, URHO, 1);
    derive_lst.addComponent("angular_momentum_x", desc_lst, State_Type, UMX, 3);

    derive_lst.add("angular_momentum_y", IndexType::TheCellType(), 1, ca_derangmomy, the_same_box);
    derive_lst.addComponent("angular_momentum_y", desc_lst, State_Type, URHO, 1);
    derive_lst.addComponent("angular_momentum_y", desc_lst, State_Type, UMX, 3);

    derive_lst.add("angular_momentum_z", IndexType::TheCellType(), 1, ca_derangmomz, the_same_box);
    derive_lst.addComponent("angular_momentum_z", desc_lst, State_Type, URHO, 1);
    derive_lst.addComponent("angular_momentum_z", desc_lst, State_Type, UMX, 3);

#ifdef GRAVITY
    derive_lst.add("maggrav", IndexType::TheCellType(), 1, ca_dermaggrav, the_same_box);
    derive_lst.addComponent("maggrav", desc_lst, Gravity_Type, 0, 3);
#endif

#ifdef AMREX_PARTICLES
    //
    // We want a derived type that corresponds to the number of particles
    // in each cell.  We only intend to use it in plotfiles for debugging
    // purposes.  We'll just use the DERNULL since don't do anything in
    // fortran for now.  We'll actually set the values in writePlotFile().
    //
    derive_lst.add("particle_count", IndexType::TheCellType(), 1, ca_dernull, the_same_box);
    derive_lst.addComponent("particle_count", desc_lst, State_Type, URHO, 1);

    derive_lst.add("total_particle_count", IndexType::TheCellType(), 1, ca_dernull, the_same_box);
    derive_lst.addComponent("total_particle_count", desc_lst, State_Type, URHO, 1);
#endif

#ifdef RADIATION
    if (Radiation::do_multigroup) {
        derive_lst.add("Ertot", IndexType::TheCellType(), 1, ca_derertot, the_same_box);
        derive_lst.addComponent("Ertot", desc_lst, Rad_Type, 0, Radiation::nGroups);
    }
#endif

#ifdef MHD
    //Electric Field at the face
    //
    //Magentic Field Cell Centered
    //x component
    derive_lst.add("B_x", IndexType::TheCellType(), 1, ca_dermagcenx, the_same_box);
    derive_lst.addComponent("B_x", desc_lst, Mag_Type_x, 0, 1);
    //y component
    derive_lst.add("B_y", IndexType::TheCellType(), 1, ca_dermagceny, the_same_box);
    derive_lst.addComponent("B_y", desc_lst, Mag_Type_y, 0, 1);
    //z component
    derive_lst.add("B_z", IndexType::TheCellType(), 1, ca_dermagcenz, the_same_box);
    derive_lst.addComponent("B_z", desc_lst, Mag_Type_z, 0, 1);

    //Electric Field
    //x component
    derive_lst.add("E_x", IndexType::TheCellType(), 1, ca_derex, the_same_box);
    derive_lst.addComponent("E_x", desc_lst, Mag_Type_y, 0, 1);
    derive_lst.addComponent("E_x", desc_lst, Mag_Type_z, 0, 1);
    derive_lst.addComponent("E_x", desc_lst, State_Type, Density, 1);  //For velocities
    derive_lst.addComponent("E_x", desc_lst, State_Type, Ymom, 1);
    derive_lst.addComponent("E_x", desc_lst, State_Type, Zmom, 1);

    //y component
    derive_lst.add("E_y", IndexType::TheCellType(), 1, ca_derey, the_same_box);
    derive_lst.addComponent("E_y", desc_lst, Mag_Type_x, 0, 1);
    derive_lst.addComponent("E_y", desc_lst, Mag_Type_z, 0, 1);
    derive_lst.addComponent("E_y", desc_lst, State_Type, Density, 1);  //For velocities
    derive_lst.addComponent("E_y", desc_lst, State_Type, Xmom, 1);
    derive_lst.addComponent("E_y", desc_lst, State_Type, Zmom, 1);

    //z component
    derive_lst.add("E_z", IndexType::TheCellType(), 1, ca_derez, the_same_box);
    derive_lst.addComponent("E_z", desc_lst, Mag_Type_x, 0, 1);
    derive_lst.addComponent("E_z", desc_lst, Mag_Type_y, 0, 1);
    derive_lst.addComponent("E_z", desc_lst, State_Type, Density, 1);  //For velocities
    derive_lst.addComponent("E_z", desc_lst, State_Type, Xmom, 1);
    derive_lst.addComponent("E_z", desc_lst, State_Type, Ymom, 1);

    //Divergence of B
    derive_lst.add("Div_B", IndexType::TheCellType(), 1, ca_derdivb, the_same_box);
    derive_lst.addComponent("Div_B", desc_lst, Mag_Type_x, 0, 1);
    derive_lst.addComponent("Div_B", desc_lst, Mag_Type_y, 0, 1);
    derive_lst.addComponent("Div_B", desc_lst, Mag_Type_z, 0, 1);

#endif

    for (int i = 0; i < NumAux; i++) {
        derive_lst.add(aux_names[i], IndexType::TheCellType(), 1, ca_derspec, the_same_box);
        derive_lst.addComponent(aux_names[i], desc_lst, State_Type, URHO, 1);
        derive_lst.addComponent(aux_names[i], desc_lst, State_Type, UFX + i, 1);
    }

    //
    // Problem-specific adds
#include <Problem_Derives.H>

    //
    // DEFINE ERROR ESTIMATION QUANTITIES
    //
    err_list_names.push_back("density");
    err_list_ng.push_back(1);

    err_list_names.push_back("Temp");
    err_list_ng.push_back(1);

    err_list_names.push_back("pressure");
    err_list_ng.push_back(1);

    err_list_names.push_back("x_velocity");
    err_list_ng.push_back(1);

#if (BL_SPACEDIM >= 2)
    err_list_names.push_back("y_velocity");
    err_list_ng.push_back(1);
#endif

#if (BL_SPACEDIM == 3)
    err_list_names.push_back("z_velocity");
    err_list_ng.push_back(1);
#endif

#ifdef REACTIONS
    err_list_names.push_back("t_sound_t_enuc");
    err_list_ng.push_back(0);

    err_list_names.push_back("enuc");
    err_list_ng.push_back(0);
#endif

#ifdef RADIATION
    if (do_radiation && !Radiation::do_multigroup) {
        err_list_names.push_back("rad");
        err_list_ng.push_back(1);
    }
#endif

    // Save the number of built-in functions; this will help us
    // distinguish between those, and the ones the user is about to add.

    num_err_list_default = err_list_names.size();

    //
    // Construct an array holding the names of the source terms.
    //

    source_names.resize(num_src);

    // Fill with an empty string to initialize.

    for (int n = 0; n < num_src; ++n) {
        source_names[n] = "";
    }

    source_names[ext_src] = "user-defined external";
    source_names[thermo_src] = "pdivU source";
#ifdef SPONGE
    source_names[sponge_src] = "sponge";
#endif

#ifdef DIFFUSION
    source_names[diff_src] = "diffusion";
#endif

#ifdef HYBRID_MOMENTUM
    source_names[hybrid_src] = "hybrid";
#endif

#ifdef GRAVITY
    source_names[grav_src] = "gravity";
#endif

#ifdef ROTATION
    source_names[rot_src] = "rotation";
#endif

#ifdef AMREX_USE_CUDA
    // Set the minimum number of threads needed per
    // threadblock to do BC fills with CUDA. We will
    // force this to be 8. The reason is that it is
    // not otherwise guaranteed for our thread blocks
    // to be aligned with the grid in such a way that
    // the synchronization logic in amrex_filccn works
    // out. We need at least NUM_GROW + 1 threads in a
    // block for CTU. If we used this minimum of 5, we
    // would hit cases where this doesn't work since
    // our blocking_factor is usually a power of 2, and
    // the thread blocks would not be aligned to guarantee
    // that the threadblocks containing the ghost zones
    // contained all of the ghost zones, as well as the
    // required interior zone. And for reflecting BCs,
    // we need NUM_GROW * 2 == 8 threads anyway. This logic
    // then requires that blocking_factor be a multiple
    // of 8. It is a little wasteful for SDC and for
    // problems that only have outflow BCs, but the BC
    // fill is not the expensive part of the algorithm
    // for our production science problems anyway, so
    // we ignore this extra cost in favor of safety.

    for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
        numBCThreadsMin[dim] = 8;
    }
#endif

#ifdef TRUE_SDC
    if (sdc_quadrature == 0) {
        // Gauss-Lobatto

        if (sdc_order == 2) {
            // trapezoid
            SDC_NODES = 2;

            dt_sdc.resize(SDC_NODES);
            dt_sdc = {0.0, 1.0};

            node_weights.resize(SDC_NODES);
            node_weights = {0.5, 0.5};

        } else if (sdc_order == 4) {
            // Simpsons
            SDC_NODES = 3;

            dt_sdc.resize(SDC_NODES);
            dt_sdc = {0.0, 0.5, 1.0};

            node_weights.resize(SDC_NODES);
            node_weights = {1.0 / 6.0, 4.0 / 6.0, 1.0 / 6.0};

        } else {
            amrex::Error("invalid value of sdc_order");
        }

    } else if (sdc_quadrature == 1) {
        // Radau

        if (sdc_order == 2) {
            SDC_NODES = 3;

            dt_sdc.resize(SDC_NODES);
            dt_sdc = {0.0, 1.0 / 3.0, 1.0};

            node_weights.resize(SDC_NODES);
            node_weights = {0.0, 3.0 / 4.0, 1.0 / 4.0};

        } else if (sdc_order == 4) {
            SDC_NODES = 4;

            dt_sdc.resize(SDC_NODES);
            dt_sdc = {0.0, (4.0 - std::sqrt(6.0)) / 10.0, (4.0 + std::sqrt(6.0)) / 10.0, 1.0};

            node_weights.resize(SDC_NODES);
            node_weights = {0.0, (16.0 - std::sqrt(6.0)) / 36.0, (16.0 + std::sqrt(6.0)) / 36.0,
                            1.0 / 9.0};

        } else {
            amrex::Error("invalid value of sdc_order");
        }

        node_weights.resize(SDC_NODES);
        node_weights = {1.0 / 6.0, 4.0 / 6.0, 1.0 / 6.0};

    } else {
        amrex::Error("invalid value of sdc_quadrature");
    }
#endif
}
