// functions specific for MGFLDSolver (SolverType 6)

#include <AMReX_LO_BCTYPES.H>

#include <Radiation.H>
#include <RadSolve.H>

#include <Castro_F.H>

#include <RAD_F.H>

#include <iostream>
#include <iomanip>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace amrex;

void Radiation::MGFLD_implicit_update(int level, int iteration, int ncycle)
{ 
  BL_PROFILE("Radiation::MGFLD_implicit_update");
  if (verbose) {
      amrex::Print() << "Radiation MGFLD implicit update, level " << level << "..." << std::endl;
  }

  BL_ASSERT(Radiation::nGroups > 0);

  int fine_level =  parent->finestLevel();

  // allocation and intialization
  Castro *castro = dynamic_cast<Castro*>(&parent->getLevel(level));
  const BoxArray& grids = castro->boxArray();
  const DistributionMapping& dmap = castro->DistributionMap();
  Real delta_t = parent->dtLevel(level);

  int ngrow = 1; 

  Real time = castro->get_state_data(Rad_Type).curTime();
  Real oldtime = castro->get_state_data(Rad_Type).prevTime();

  MultiFab& S_new = castro->get_new_data(State_Type);
  AmrLevel::FillPatch(*castro,S_new,ngrow,time,State_Type,0,S_new.nComp(),0); 

  Array<MultiFab, AMREX_SPACEDIM> lambda;
  if (limiter > 0) {
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
        lambda[idim].define(castro->getEdgeBoxArray(idim), dmap, nGroups, 0);
    }

    if (inner_update_limiter == -1) {
      MultiFab& Er_lag = castro->get_old_data(Rad_Type);
      Er_lag.setBndry(-1.0);
      Er_lag.FillBoundary(parent->Geom(level).periodicity());

      MultiFab& S_lag = castro->get_old_data(State_Type);
      for (FillPatchIterator fpi(*castro,S_lag,ngrow,oldtime,State_Type,
                                 0,S_lag.nComp()); fpi.isValid(); ++fpi) {
          S_lag[fpi].copy<RunOn::Device>(fpi());
      }

      MultiFab kpr_lag(grids,dmap,nGroups,1);
      MGFLD_compute_rosseland(kpr_lag, S_lag); 

      for (int igroup=0; igroup<nGroups; ++igroup) {
        scaledGradient(level, lambda, kpr_lag, igroup, Er_lag, igroup, limiter, 1, igroup);
        // lambda now contains scaled gradient
        fluxLimiter(level, lambda, limiter, igroup);
        // lambda now contains flux limiter
      }
    }
  }
  else {
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
        lambda[idim].define(castro->getEdgeBoxArray(idim), dmap, 1, 0);
      lambda[idim].setVal(1./3.);
    }    
  }

  // Er_new: work copy
  // Er_old: the input state of the implicit update
  // Er_pi: previous inner iteration
  // Er_star: previous outer iteration
  //
  MultiFab& Er_new = castro->get_new_data(Rad_Type);
  {
    MultiFab rhs(grids,dmap,1,0);
    for (int igroup=0; igroup<nGroups; igroup++) {
      rhs.setVal(0.0);
      deferred_sync(level, rhs, igroup);
      rhstoEr(rhs, delta_t, level);
      MultiFab::Add(Er_new, rhs, 0, igroup, 1, 0);
    }
  }
  MultiFab Er_old(grids, dmap , Er_new.nComp(), 0);
  MultiFab::Copy(Er_old, Er_new, 0, 0, Er_new.nComp(), 0);
  MultiFab Er_pi(grids,dmap,nGroups,1);
  MultiFab Er_star(grids, dmap, nGroups, 1);
  Er_pi.setBndry(-1.0); // later we may use it to compute limiter
  Er_star.setBndry(-1.0); // later we may use it to compute limiter

  MultiFab rhoe_new(grids,dmap,1,0);
  MultiFab rhoe_old(grids,dmap,1,0);
  MultiFab rhoe_star(grids,dmap,1,0);

  MultiFab rho(grids,dmap,1,1);
  MultiFab temp_new(grids,dmap,1,1); // ghost cell for kappa_r
  MultiFab temp_star(grids,dmap,1,0);

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi) {

      const Box &gbx = mfi.growntilebox(1);
      const Box &bx  = mfi.tilebox();

      rho[mfi].copy<RunOn::Device>(S_new[mfi],gbx, URHO, gbx,0,1);

      rhoe_new[mfi].copy<RunOn::Device>(S_new[mfi], bx, UEINT, bx,0,1);
      rhoe_old[mfi].copy<RunOn::Device>(rhoe_new[mfi], bx);

      temp_new[mfi].copy<RunOn::Device>(S_new[mfi],gbx, UTEMP, gbx,0,1);
    
  }

  // Planck mean and Rosseland 
  MultiFab kappa_p(grids,dmap,nGroups,1);
  MultiFab kappa_r(grids,dmap,nGroups,1); 

  // emissivity, j_g = \int j_nu dnu
  // j_nu = 4 pi /c * \eta_0^{th} = \kappa_0 * B_\nu (assuming LTE),
  // where B_\nu is the usual Planck function \times 4 pi / c
  MultiFab jg(grids,dmap,nGroups,1);    
  MultiFab djdT(grids,dmap,nGroups,1);  
  MultiFab dkdT(grids,dmap,nGroups,1);  
  MultiFab etaT(grids,dmap,1,0);
  MultiFab etaTz(grids,dmap,1,0);
  MultiFab eta1(grids,dmap,1,0); // eta1 = 1 - etaT

  MultiFab& mugT = djdT;

  MultiFab dedT(grids,dmap,1,0);

  MultiFab coupT(grids,dmap,1,0); // \sum{\kappa E - j}

  // multigroup boundary object
  MGRadBndry mgbd(grids,dmap, nGroups, castro->Geom());
  getBndryDataMG(mgbd, Er_new, time, level);

  bool have_Sanchez_Pomraning = false;
  int lo_bc[3]={0}, hi_bc[3]={0};
  for (int idim=0; idim<AMREX_SPACEDIM; idim++) {
    lo_bc[idim] = rad_bc.lo(idim);
    hi_bc[idim] = rad_bc.hi(idim);
    if (lo_bc[idim] == LO_SANCHEZ_POMRANING || 
        hi_bc[idim] == LO_SANCHEZ_POMRANING) {
      have_Sanchez_Pomraning = true;
    }
  }

  RadSolve* const solver = castro->rad_solver.get();

  Real relative_in, absolute_in, error_er;
  Real rel_rhoe, abs_rhoe;
  Real rel_T, abs_T;
  Real rel_FT, abs_FT;

  // point to flux_trial (or NULL) for appropriate levels
  FluxRegister* flux_in = (level < fine_level) ? flux_trial[level+1].get() : nullptr;
  FluxRegister* flux_out = (level > 0) ? flux_trial[level].get() : nullptr;

  Array<MultiFab, AMREX_SPACEDIM> Flux;
  for (int n = 0; n < AMREX_SPACEDIM; n++) {
      Flux[n].define(castro->getEdgeBoxArray(n), dmap, 1, 0);
  }

  std::unique_ptr<MultiFab> flxsave;
  MultiFab* flxcc;
  int icomp_flux = -1;
  if (plot_com_flux) {
      flxcc = plotvar[level].get();
      icomp_flux = icomp_com_Fr;
  } else if (plot_lab_Er || plot_lab_flux) {
      flxsave.reset(new MultiFab(grids, dmap, nGroups*AMREX_SPACEDIM, 0));
      flxcc = flxsave.get();
      icomp_flux = 0;
  } 

  // Er_step: starting state of the inner iteration (e.g., ^(2))
  // There used to be an extra velocity term update
  MultiFab& Er_step = Er_old;
  MultiFab& rhoe_step = rhoe_old;

  Real reltol_in = relInTol;
  Real ptc_tau = 0.0;  // not being used 

  // nonlinear loop for all groups
  int it = 0;
  bool conservative_update = false;
  bool outer_ready = false;
  bool converged = false;
  bool inner_converged = false;
  do {
    it++;

    if (it == 1) {
      eos_opacity_emissivity(S_new, temp_new,
                             temp_star, // input
                             kappa_p, kappa_r, jg, 
                             djdT, dkdT, dedT, // output
                             level, it, 1); 
      // It's OK that temp_star does not have a valid value for it==1
    }

    MultiFab::Copy(rhoe_star, rhoe_new, 0, 0, 1, 0);
    MultiFab::Copy(temp_star, temp_new, 0, 0, 1, 0);
    MultiFab::Copy(Er_star, Er_new, 0, 0, nGroups, 0);

    if (limiter>0 && inner_update_limiter==0) {
      Er_star.FillBoundary(parent->Geom(level).periodicity());

      for (int igroup=0; igroup<nGroups; ++igroup) {
        scaledGradient(level, lambda, kappa_r, igroup, Er_star, igroup, limiter, 1, igroup);
        // lambda now contains scaled gradient
        fluxLimiter(level, lambda, limiter, igroup);
        // lambda now contains flux limiter
      }
    }
    
    // djdT is both input and output
    compute_etat(etaT, etaTz,
                 eta1, djdT, 
                 dkdT, dedT,
                 Er_star, rho, 
                 delta_t, ptc_tau);

    // After this, djdT contains mugT.

    // The inner loops does not update rhoe and T
    int innerIteration = 0;
    inner_converged = false;
    Real relative_in_prev = 1.e200, absolute_in_prev = 1.e200;
    bool accel_allowed = true;
    do {
      innerIteration++;

      MultiFab::Copy(Er_pi, Er_new, 0, 0, nGroups, 0);

      if (limiter>0 && inner_update_limiter>0) { 
        if (innerIteration <= inner_update_limiter) {
          Er_pi.FillBoundary(parent->Geom(level).periodicity());
          
          for (int igroup=0; igroup<nGroups; ++igroup) {
            scaledGradient(level, lambda, kappa_r, igroup, Er_pi, igroup, limiter, 1, igroup);
            // lambda now contains scaled gradient
            fluxLimiter(level, lambda, limiter, igroup);
            // lambda now contains flux limiter
          }
        }
      }

      compute_coupling(coupT, kappa_p, Er_pi, jg);

      for (int igroup=0; igroup<nGroups; ++igroup) {

        set_current_group(igroup);

        // setup and solve linear system

        // set boundary condition
        solver->levelBndry(mgbd, igroup);
        
        solver->levelACoeffs(level, kappa_p, delta_t, c, igroup, ptc_tau);

        int lamcomp = (limiter==0) ? 0 : igroup;
        solver->levelBCoeffs(level, lambda, kappa_r, igroup, c, lamcomp);

        if (have_Sanchez_Pomraning) {
          solver->levelSPas(level, lambda, igroup, lo_bc, hi_bc);
        }

        { // src and rhd block
                  
          MultiFab rhs(grids,dmap,1,0);

          solver->levelRhs(level, rhs, jg, mugT,
                           coupT, etaT,
                           Er_step, rhoe_step, Er_star, rhoe_star,
                           delta_t, igroup, it, ptc_tau);

          // solve Er equation and put solution in Er_new(igroup)
          solver->levelSolve(level, Er_new, igroup, rhs, 0.01);
        } // end src and rhs block

        solver->levelFlux(level, Flux, Er_new, igroup);
        solver->levelFluxReg(level, flux_in, flux_out, Flux, igroup);
          
        if (icomp_flux >= 0) 
            solver->levelFluxFaceToCenter(level, Flux, *flxcc, icomp_flux+igroup);

      } // end loop over groups
      
      // Check for convergence *before* acceleration step:
      check_convergence_er(relative_in, absolute_in, error_er, Er_new, Er_pi,
                           kappa_p, etaTz, temp_new, delta_t);

      if (verbose >= 2) {
        int oldprec = std::cout.precision(3);
        amrex::Print() << "Outer = " << it << ", Inner = " << innerIteration
                       << ", inner err =  " << std::setw(8) << relative_in << " (rel),  " 
                       << std::setw(8) << absolute_in << " (abs)" << std::endl;
        std::cout.precision(oldprec);
      }

      if (relative_in < 1.e-15) {
        inner_converged = true;
      }
      else if (innerIteration < minInIter) {
        inner_converged = false;
      }
      else if ( (relative_in <= reltol_in || absolute_in <= absInTol) 
                && error_er <= reltol ) {
        inner_converged = true;
      }

      if (!inner_converged) {
        Real accel_fac=1.+1.e-6;
        if (skipAccelAllowed &&
            relative_in>accel_fac*relative_in_prev && 
            absolute_in>accel_fac*absolute_in_prev) {
          accel_allowed = false;
          if (relative_in>10.*relative_in_prev && 
              absolute_in>10.*absolute_in_prev) {
            MultiFab::Copy(Er_new, Er_star, 0, 0, nGroups, 0);
          }
        }
        relative_in_prev = relative_in;
        absolute_in_prev = absolute_in;

        if (accel_allowed) {
          if (accelerate == 1) {
            local_accel(Er_new, Er_pi, kappa_p, etaT,
                        mugT, delta_t, ptc_tau);
          } 
          else if (accelerate == 2) {
            gray_accel(Er_new, Er_pi, kappa_p, kappa_r, 
                       etaT, eta1, mugT,
                       lambda, solver, mgbd, grids, level, time, delta_t, ptc_tau);
          } 
        }
      }

    } while(!inner_converged && innerIteration < maxInIter); 

    if (verbose == 1) {
      int oldprec = std::cout.precision(3);
      amrex::Print() << "Outer = " << it << ", Inner = " << innerIteration
                     << ", inner tol =  " << std::setw(8) << relative_in << "  " 
                     << std::setw(8) << absolute_in << std::endl;
      std::cout.precision(oldprec);
    }

    // update rhoe and T
    if (matter_update_type == 0) {
      conservative_update = true;
    }
    else if (matter_update_type == 1) {
      if (outer_ready || it >= maxiter) {
        conservative_update = true;
      }
      else {
        conservative_update = false;
      }
    }
    else if (matter_update_type == 2) {
      if (it%2 == 0) {
        conservative_update = true;
      }
      else {
        conservative_update = false;
      }
    }
    else if (matter_update_type == 3) {
      if (outer_ready || it >= maxiter) {
        conservative_update = true;
      }
      else {
        if (it%2 == 0) {
          conservative_update = true;
        }
        else {
          conservative_update = false;
        }
      }
    }
    else {
      conservative_update = true;
    }

    update_matter(rhoe_new, temp_new, Er_new, Er_pi,
                  rhoe_star, rhoe_step,
                  etaT, etaTz, eta1,
                  coupT,
                  kappa_p, jg, mugT,
                  S_new, level, delta_t, ptc_tau, it, conservative_update);

    eos_opacity_emissivity(S_new, temp_new,
                           temp_star, // input
                           kappa_p, kappa_r, jg, 
                           djdT, dkdT, dedT, // output
                           level, it+1, 0);

    check_convergence_matt(rhoe_new, rhoe_star, rhoe_step, Er_new,
                           temp_new, temp_star, 
                           rho, kappa_p, jg, dedT,
                           rel_rhoe, abs_rhoe, rel_FT, abs_FT, rel_T, abs_T,
                           delta_t);

    Real relative_out, absolute_out;

    switch (convergence_check_type) {
    case 1:
      relative_out = rel_rhoe;
      absolute_out = abs_rhoe;
      break;
    case 2:
      relative_out = rel_FT;
      absolute_out = abs_FT;
      break;
    case 3:
      relative_out = rel_T;
      absolute_out = abs_T;
      break;
    default:
      relative_out = rel_T;
      if (conservative_update) {
        //      relative_out = (relative_out > rel_rhoe) ? relative_out : rel_rhoe;
        relative_out = (relative_out > rel_FT  ) ? relative_out : rel_FT;
      }
      absolute_out = abs_T;
    }

    if (verbose >= 2) {
      int oldprec = std::cout.precision(4);
      amrex::Print() << "Update Errors for      rhoe,        FT,         T" 
                     << std::endl;
      amrex::Print() << "       Relative = " << std::setw(9) << rel_rhoe << ", " 
                     << std::setw(9) << rel_FT << ", " << std::setw(9) << rel_T 
                     << std::endl;
      amrex::Print() << "       Absolute = " << std::setw(9) << abs_rhoe << ", " 
                     << std::setw(9) << abs_FT << ", " << std::setw(9) << abs_T 
                     << std::endl;
      std::cout.precision(oldprec);
    }

    if (it < miniter) {
      converged = false;
    }
    else if (relative_out <= reltol || absolute_out <= abstol) { 
      //      || rel_rhoe < 1.e-15) {
      converged = true;
    }
    else {
      converged = false;
    }

    //    if (relative_out <= 10.*reltol) {
    if (relative_out <= reltol) {
      outer_ready = true;
    }

    if (!converged && it > n_bisect) {
      bisect_matter(rhoe_new, temp_new,
                    rhoe_star, temp_star,
                    S_new, grids, level);

      eos_opacity_emissivity(S_new, temp_new,
                             temp_star, // input
                             kappa_p, kappa_r, jg, 
                             djdT, dkdT, dedT, // output
                             level, it+1, 0);
    }
   
  } while ( ((!converged || !inner_converged) && it<maxiter)
            || !conservative_update);

  if (verbose == 1) {
    int oldprec = std::cout.precision(4);
    amrex::Print() << "Update Errors for      rhoe,        FT,         T" 
                   << std::endl;
    amrex::Print() << "       Relative = " << std::setw(9) << rel_rhoe << ", " 
                   << std::setw(9) << rel_FT << ", " << std::setw(9) << rel_T 
                   << std::endl;
    amrex::Print() << "       Absolute = " << std::setw(9) << abs_rhoe << ", " 
                   << std::setw(9) << abs_FT << ", " << std::setw(9) << abs_T 
                   << std::endl;
    std::cout.precision(oldprec);
  }

  if (!converged) {
      amrex::Abort("Implicit Update Failed to Converge");
  }

  // update flux registers

  flux_in = (level < fine_level) ? flux_trial[level+1].get() : nullptr;
  if (flux_in) {
      for (OrientationIter face; face; ++face) {
          Orientation ori = face();
          (*flux_cons[level+1])[ori].linComb(1.0, -1.0,
                                             (*flux_in)[ori], 0, 0, nGroups);
      }
  }

  flux_out = (level > 0) ? flux_trial[level].get() : nullptr;
  if (flux_out) {
      for (OrientationIter face; face; ++face) {
          Orientation ori = face();
          (*flux_cons[level])[ori].linComb(1.0, 1.0 / ncycle,
                                        (*flux_out)[ori], 0, 0, nGroups);
    }
  }

  if (level == fine_level && iteration == 1) {

    // We've now advanced the finest level, so we've just passed the
    // last time we might want to reflux from any existing flux_cons_old
    // on any level until the next sync.  Setting the stored delta_t
    // to 0.0 is a hack equivalent to marking each corresponding
    // register as "not to be used".

    for (int flev = level; flev > 0; flev--) {
      delta_t_old[flev-1] = 0.0;
    }
  }

  Real derat = 0.0;
  Real dTrat = 0.0;

  // update state with new fluid energy and temperature
  state_energy_update(S_new, rhoe_new, temp_new, grids, derat, dTrat, level);

  delta_e_rat_level[level] = derat;
  delta_T_rat_level[level] = dTrat;

  if (verbose >= 2) {
      amrex::Print() << "Delta T      Ratio = " << dTrat << std::endl;
  }

  if (plot_lambda) {
      save_lambda_in_plotvar(level, lambda);
  }

  if (plot_kappa_p) {
      MultiFab::Copy(*plotvar[level], kappa_p, 0, icomp_kp, nGroups, 0);
  }

  if (plot_kappa_r) {
      MultiFab::Copy(*plotvar[level], kappa_r, 0, icomp_kr, nGroups, 0);
  }

  if (plot_lab_Er) {
      save_lab_Er_in_plotvar(level, S_new, Er_new, *flxcc, icomp_flux);
  }

  // if (plot_com_flux) {
  //     already done when calling solver->levelFluxFaceToCenter()
  // }

  if (plot_lab_flux) {
      save_flux_in_plotvar(level, S_new, lambda, Er_new, *flxcc, icomp_flux);
  }

  if (verbose) {
      amrex::Print() << "                                     done" << std::endl;
  }
}
