#include <iomanip>

#include <Castro.H>
#include <Castro_F.H>

#ifdef GRAVITY
#include <Gravity.H>
#include <Gravity_F.H>
#endif

Real
Castro::volWgtSum (const std::string& name,
                   Real               time)
{
    Real        sum     = 0.0;
    const Real* dx      = geom.CellSize();
    MultiFab*   mf      = derive(name,time,0);

    BL_ASSERT(mf != 0);

    BoxArray baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        if (level < parent->finestLevel())
        {
            std::vector< std::pair<int,Box> > isects = baf.intersections(grids[mfi.index()]);

            for (int ii = 0; ii < isects.size(); ii++)
            {
                fab.setVal(0,isects[ii].second,0,fab.nComp());
            }
        }
        Real s = 0.0;
        const Box& box  = mfi.validbox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();
#if(BL_SPACEDIM < 3) 
        const Real* rad = radius[mfi.index()].dataPtr();
        int irlo        = lo[0]-radius_grow;
        int irhi        = hi[0]+radius_grow;
#endif

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //
#if(BL_SPACEDIM == 1) 
	BL_FORT_PROC_CALL(CA_SUMMASS,ca_summass)
            (BL_TO_FORTRAN(fab),lo,hi,dx,&s,rad,irlo,irhi);
#elif(BL_SPACEDIM == 2)
	BL_FORT_PROC_CALL(CA_SUMMASS,ca_summass)
            (BL_TO_FORTRAN(fab),lo,hi,dx,&s,rad,irlo,irhi);
#elif(BL_SPACEDIM == 3)
	BL_FORT_PROC_CALL(CA_SUMMASS,ca_summass)
            (BL_TO_FORTRAN(fab),lo,hi,dx,&s);
#endif
        sum += s;
    }

    delete mf;

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
Castro::volWgtSquaredSum (const std::string& name,
                          Real               time)
{
    Real        sum     = 0.0;
    const Real* dx      = geom.CellSize();
    MultiFab*   mf      = derive(name,time,0);

    BL_ASSERT(mf != 0);

    BoxArray baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        if (level < parent->finestLevel())
        {
            std::vector< std::pair<int,Box> > isects = baf.intersections(grids[mfi.index()]);

            for (int ii = 0; ii < isects.size(); ii++)
            {
                fab.setVal(0,isects[ii].second,0,fab.nComp());
            }
        }
        Real s = 0.0;
        const Box& box  = mfi.validbox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();
#if(BL_SPACEDIM < 3) 
        const Real* rad = radius[mfi.index()].dataPtr();
        int irlo        = lo[0]-radius_grow;
        int irhi        = hi[0]+radius_grow;
#endif

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //
#if(BL_SPACEDIM == 1) 
	BL_FORT_PROC_CALL(CA_SUMSQUARED,ca_summass)
            (BL_TO_FORTRAN(fab),lo,hi,dx,&s,rad,irlo,irhi);
#elif(BL_SPACEDIM == 2)
	BL_FORT_PROC_CALL(CA_SUMSQUARED,ca_summass)
            (BL_TO_FORTRAN(fab),lo,hi,dx,&s,rad,irlo,irhi);
#elif(BL_SPACEDIM == 3)
	BL_FORT_PROC_CALL(CA_SUMSQUARED,ca_summass)
            (BL_TO_FORTRAN(fab),lo,hi,dx,&s);
#endif
        sum += s;
    }

    delete mf;

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
Castro::locWgtSum (const std::string& name,
                   Real               time,
                   int                idir)
{
    Real        sum     = 0.0;
    const Real* dx      = geom.CellSize();
    MultiFab*   mf      = derive(name,time,0);

    BL_ASSERT(mf != 0);

    BoxArray baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        if (level < parent->finestLevel())
        {
            std::vector< std::pair<int,Box> > isects = baf.intersections(grids[mfi.index()]);

            for (int ii = 0; ii < isects.size(); ii++)
            {
                fab.setVal(0,isects[ii].second,0,fab.nComp());
            }
        }
        Real s = 0.0;
        const Box& box  = mfi.validbox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();
#if (BL_SPACEDIM < 3)
        const Real* rad = radius[mfi.index()].dataPtr();
        int irlo        = lo[0]-radius_grow;
        int irhi        = hi[0]+radius_grow;
#endif

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //
#if (BL_SPACEDIM == 1) 
	BL_FORT_PROC_CALL(CA_SUMLOCMASS,ca_sumlocmass)
            (BL_TO_FORTRAN(fab),lo,hi,geom.ProbLo(),dx,&s,rad,irlo,irhi,idir);
#elif (BL_SPACEDIM == 2)
        int geom_flag = Geometry::IsRZ() ? 1 : 0;
        if (idir == 0 && geom_flag == 1) {
            s = 0.0;
        } else {
	   BL_FORT_PROC_CALL(CA_SUMLOCMASS,ca_sumlocmass)
               (BL_TO_FORTRAN(fab),lo,hi,geom.ProbLo(),dx,&s,rad,irlo,irhi,idir);
        }
#else
	BL_FORT_PROC_CALL(CA_SUMLOCMASS,ca_sumlocmass)
            (BL_TO_FORTRAN(fab),lo,hi,geom.ProbLo(),dx,&s,idir);
#endif
        sum += s;
    }

    delete mf;

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
Castro::volWgtSumMF (MultiFab* mf, int comp) 
{
    Real        sum     = 0.0;
    const Real* dx      = geom.CellSize();

    BL_ASSERT(mf != 0);

    BoxArray baf;

    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        Real s = 0.0;
        const Box& box  = mfi.validbox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();
#if(BL_SPACEDIM < 3) 
        const Real* rad = radius[mfi.index()].dataPtr();
        int irlo        = lo[0]-radius_grow;
        int irhi        = hi[0]+radius_grow;
#endif

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //
#if(BL_SPACEDIM == 1) 
	BL_FORT_PROC_CALL(CA_SUMMASS,ca_summass)
            (BL_TO_FORTRAN(fab),lo,hi,dx,&s,rad,irlo,irhi);
#elif(BL_SPACEDIM == 2)
	BL_FORT_PROC_CALL(CA_SUMMASS,ca_summass)
            (BL_TO_FORTRAN(fab),lo,hi,dx,&s,rad,irlo,irhi);
#elif(BL_SPACEDIM == 3)
	BL_FORT_PROC_CALL(CA_SUMMASS,ca_summass)
            (BL_TO_FORTRAN(fab),lo,hi,dx,&s);
#endif
        sum += s;
    }

    delete mf;

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
Castro::volWgtSumOneSide (const std::string& name,
                          Real               time, 
                          int                side,
                          int                bdir)
{
    // This function is a clone of volWgtSum except it computes the result only on half of the domain.
    // The lower half corresponds to side == 0 and the upper half corresponds to side == 1.
    // The argument bdir gives the direction along which to bisect.
    // ONLY WORKS IN THREE DIMENSIONS.

    Real        sum     = 0.0;
    const Real* dx      = geom.CellSize();
    MultiFab*   mf      = derive(name,time,0);
    const int* domlo    = geom.Domain().loVect(); 
    const int* domhi    = geom.Domain().hiVect();

    BL_ASSERT(mf != 0);

    BoxArray baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        if (level < parent->finestLevel())
        {
            std::vector< std::pair<int,Box> > isects = baf.intersections(grids[mfi.index()]);

            for (int ii = 0; ii < isects.size(); ii++)
            {
                fab.setVal(0,isects[ii].second,0,fab.nComp());
            }
        }
        Real s = 0.0;
        const Box& box  = mfi.validbox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();

        int hiLeft[3]       = { *hi, *(hi+1), *(hi+2) };
        hiLeft[bdir]        = *(domhi+bdir) / 2;
        int loRight[3]      = { *lo, *(lo+1), *(lo+2) };
        loRight[bdir]       = *(domhi+bdir) / 2 + 1;

        const int* hiLeftPtr  = hiLeft;
        const int* loRightPtr = loRight;

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //
        
        if ( side == 0 && *(lo + bdir) <= *(hiLeftPtr + bdir) ) {
          if ( *(hi + bdir) <= *(hiLeftPtr + bdir) )
            BL_FORT_PROC_CALL(CA_SUMMASS,ca_summass)
              (BL_TO_FORTRAN(fab),lo,hi,dx,&s);
          else
            BL_FORT_PROC_CALL(CA_SUMMASS,ca_summass)
              (BL_TO_FORTRAN(fab),lo,hiLeftPtr,dx,&s);
	}  
        else if ( side == 1 && *(hi + bdir) >= *(loRightPtr + bdir) ) {
          if ( *(lo + bdir) >= *(loRightPtr + bdir) )
            BL_FORT_PROC_CALL(CA_SUMMASS,ca_summass)
              (BL_TO_FORTRAN(fab),lo,hi,dx,&s);
          else
            BL_FORT_PROC_CALL(CA_SUMMASS,ca_summass)
              (BL_TO_FORTRAN(fab),loRightPtr,hi,dx,&s);
	}
        
        sum += s;
		
    }

    delete mf;

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
Castro::locWgtSumOneSide (const std::string& name,
                          Real               time,
                          int                idir, 
                          int                side,
                          int                bdir)
{
  // This function is a clone of locWgtSum except that it only sums over one half of the domain.
  // The lower half corresponds to side == 0, and the upper half corresponds to side == 1.
  // The argument idir (x == 0, y == 1, z == 2) gives the direction to location weight by,
  // and the argument bdir gives the direction along which to bisect.
  // ONLY WORKS IN THREE DIMENSIONS.

    Real sum            = 0.0;
    const Real* dx      = geom.CellSize();
    MultiFab*   mf      = derive(name,time,0); 
    const int* domlo = geom.Domain().loVect(); 
    const int* domhi = geom.Domain().hiVect(); 

    BL_ASSERT(mf != 0);

    BoxArray baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        if (level < parent->finestLevel())
        {
            std::vector< std::pair<int,Box> > isects = baf.intersections(grids[mfi.index()]);

            for (int ii = 0; ii < isects.size(); ii++)
            {
                fab.setVal(0,isects[ii].second,0,fab.nComp());
            }
        }
        Real s = 0.0;
        const Box& box  = mfi.validbox();
        const int* lo   = box.loVect();
         const int* hi   = box.hiVect();

        int hiLeft[3]       = { *hi, *(hi+1), *(hi+2) };
        hiLeft[bdir]        = *(domhi+bdir) / 2;
        int loRight[3]      = { *lo, *(lo+1), *(lo+2) };
        loRight[bdir]       = (*(domhi+bdir) / 2) + 1;

        const int* hiLeftPtr  = hiLeft;
        const int* loRightPtr = loRight;
                
        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        // 
        
        if ( (side == 0) && (*(lo + bdir) <= *(hiLeftPtr + bdir)) ) {
          if ( *(hi + bdir) <= *(hiLeftPtr + bdir) )
            BL_FORT_PROC_CALL(CA_SUMLOCMASS,ca_sumlocmass)
              (BL_TO_FORTRAN(fab),lo,hi,geom.ProbLo(),dx,&s,idir);
          else
            BL_FORT_PROC_CALL(CA_SUMLOCMASS,ca_sumlocmass)
              (BL_TO_FORTRAN(fab),lo,hiLeftPtr,geom.ProbLo(),dx,&s,idir);
	}  
        else if ( (side == 1) && (*(hi + bdir) >= *(loRightPtr + bdir)) ) {
          if ( *(lo + bdir) >= *(loRightPtr + bdir) )
            BL_FORT_PROC_CALL(CA_SUMLOCMASS,ca_sumlocmass)
              (BL_TO_FORTRAN(fab),lo,hi,geom.ProbLo(),dx,&s,idir);
          else
            BL_FORT_PROC_CALL(CA_SUMLOCMASS,ca_sumlocmass)
              (BL_TO_FORTRAN(fab),loRightPtr,hi,geom.ProbLo(),dx,&s,idir);
	}
  
        sum += s;
        
    }

    delete mf;

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;

}

#ifdef GRAVITY
Real
Castro::volProductSum (const std::string& name1, 
                       const std::string& name2,
                       Real time)
{
    Real        sum     = 0.0;
    const Real* dx      = geom.CellSize();
    MultiFab*   mf1;
    MultiFab*   mf2;

    if ( name1 == "phi" )
      mf1 = gravity->get_phi_curr(level);
    else
      mf1 = derive(name1,time,0);
    
    if ( name2 == "phi" )
      mf2 = gravity->get_phi_curr(level);
    else
      mf2 = derive(name2,time,0);

    BL_ASSERT(mf1 != 0);
    BL_ASSERT(mf2 != 0);

    BoxArray baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    for (MFIter mfi(*mf1); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab1 = (*mf1)[mfi];
        FArrayBox& fab2 = (*mf2)[mfi];

        if (level < parent->finestLevel())
        {
            std::vector< std::pair<int,Box> > isects = baf.intersections(grids[mfi.index()]);

            for (int ii = 0; ii < isects.size(); ii++)
            {
                fab1.setVal(0,isects[ii].second,0,fab1.nComp());
                fab2.setVal(0,isects[ii].second,0,fab2.nComp());
            }
        }
        Real s = 0.0;
        const Box& box  = mfi.validbox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();

	BL_FORT_PROC_CALL(CA_SUMPRODUCT,ca_sumproduct)
	  (BL_TO_FORTRAN(fab1),BL_TO_FORTRAN(fab2),lo,hi,dx,&s);
        
        sum += s;
    }

    if ( name1 != "phi" )
      delete mf1;
    if ( name2 != "phi" )
      delete mf2;

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}
#endif

Real
Castro::locSquaredSum (const std::string& name,
                       Real               time,
                       int                idir)
{
    Real        sum     = 0.0;
    const Real* dx      = geom.CellSize();
    MultiFab*   mf      = derive(name,time,0);

    BL_ASSERT(mf != 0);

    BoxArray baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        if (level < parent->finestLevel())
        {
            std::vector< std::pair<int,Box> > isects = baf.intersections(grids[mfi.index()]);

            for (int ii = 0; ii < isects.size(); ii++)
            {
                fab.setVal(0,isects[ii].second,0,fab.nComp());
            }
        }
        Real s = 0.0;
        const Box& box  = mfi.validbox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();

	BL_FORT_PROC_CALL(CA_SUMLOCSQUAREDMASS,ca_sumlocsquaredmass)
            (BL_TO_FORTRAN(fab),lo,hi,geom.ProbLo(),dx,&s,idir);

        sum += s;
    }

    delete mf;

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}
