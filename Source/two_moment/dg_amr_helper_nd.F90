module dg_amr_helper_module

  implicit none

contains

  subroutine dg_refine(fblo, fbhi, &
                       fine, fine_lo, fine_hi, &
                       crse, crse_lo, crse_hi, &
                       nFineX, nvar) &
                       bind(c,name='ca_dg_refine')

    use amrex_error_module, only: amrex_error
    use amrex_fort_module, only: rt => amrex_real
    use ProgramHeaderModule           , only: nNodesX, nNodesE, nE, nDOF
    use RadiationFieldsModule         , only: nSpecies, nCR
    use TwoMoment_MeshRefinementModule, only: Refine_TwoMoment

    implicit none

    ! arguments
    integer :: crse_lo(3), crse_hi(3)
    integer :: fine_lo(3), fine_hi(3)
    integer :: fblo(3), fbhi(3)
    integer :: nFineX(3)
    integer :: nvar

    real(rt) :: fine(fine_lo(1):fine_hi(1), fine_lo(2):fine_hi(2), fine_lo(3):fine_hi(3), 0:nvar-1)
    real(rt) :: crse(crse_lo(1):crse_hi(1), crse_lo(2):crse_hi(2), crse_lo(3):crse_hi(3), 0:nvar-1)

    ! local variables
    integer :: ic , jc , kc
    integer :: icc, jcc, kcc
    integer :: ifc, jfc, kfc
    integer :: ispec, imom, id, ie, ii, ioff
    integer :: nFineX_local(3)
    integer :: cblo(3), cbhi(3)

#if (AMREX_SPACEDIM == 2)
    real(rt) :: U_FineCoarse(1:nDOF, 1:nE, 1:nFineX(1), 1:nFineX(2), 1:1)
#elif (AMREX_SPACEDIM == 3)
    real(rt) :: U_FineCoarse(1:nDOF, 1:nE, 1:nFineX(1), 1:nFineX(2), 1:nFineX(3))
#endif

    real(rt) :: U_CrseCoarse(1:nDOF, 1:nE)

    ! The C++ sends in zero in 2D but Refine_TwoMoment expects 1
    nFineX_local(:) = nFineX(:)
#if (AMREX_SPACEDIM == 2)
    nFineX_local(3) = 1
#endif

    if ( nvar .ne. nNodesE * nE * nCR * nSpecies ) &
       call amrex_error("Mismatch of nvar with nNodesE*nE*nCR*nSpecies in dg_refine")
#if (AMREX_SPACEDIM == 2)
    if ( nDOF .ne. nNodesE * nNodesX(1) * nNodesX(2) ) &
       call amrex_error("Mismatch of nDOF with nNodesE*nNodesX(1)*nNodesX(2) in dg_refine")
#elif (AMREX_SPACEDIM == 3)
    if ( nDOF .ne. nNodesE * nNodesX(1) * nNodesX(2) * nNodesX(3) ) &
       call amrex_error("Mismatch of nDOF with nNodesE*nNodesX(1)*nNodesX(2)*nNodesX(3) in dg_refine")
#endif

    ! print *,'FINE DIMS ', fine_lo(:), fine_hi(:)
    ! print *,'CRSE DIMS ', crse_lo(:), crse_hi(:)
    ! print *,'FB LO/HI  ',    fblo(:),    fbhi(:)
 
    cblo(1) = fblo(1) / nFineX(1)
    cbhi(1) = fbhi(1) / nFineX(1)
    cblo(2) = fblo(2) / nFineX(2)
    cbhi(2) = fbhi(2) / nFineX(2)
#if (AMREX_SPACEDIM == 3)
    cblo(3) = fblo(3) / nFineX(3)
    cbhi(3) = fbhi(3) / nFineX(3)
#endif

    ! loop over coarse grid under fine region...
#if (AMREX_SPACEDIM == 2)
    do kc = crse_lo(3), crse_hi(3)
#elif (AMREX_SPACEDIM == 3)
    do kc = cblo(3), cbhi(3), 2
#endif
       do jc = cblo(2), cbhi(2), 2
       do ic = cblo(1), cbhi(1), 2

             do ispec = 1, nSpecies
                do imom = 1, nCR
                   do ie = 1, nE

                          ioff = (ispec-1)*(nCR*nE*nNodesE) + (imom-1)*(nE*nNodesE) + (ie-1)*nNodesE
   
                          do id = 1, nNodesE
                             ii   = ioff + (id-1)
                             U_CrseCoarse(id,ie) = crse(ic,jc,kc,ii)
                          end do
   
                          do id = nNodesE+1, 2*nNodesE
                             ii   = ioff + (id-nNodesE-1)
                             U_CrseCoarse(id,ie) = crse(ic+1,jc,kc,ii)
                          end do
   
                          do id = 2*nNodesE+1, 3*nNodesE
                             ii   = ioff + (id-2*nNodesE-1)
                             U_CrseCoarse(id,ie) = crse(ic,jc+1,kc,ii)
                          end do
   
                          do id = 3*nNodesE+1, 4*nNodesE
                          ii   = ioff + (id-3*nNodesE-1)
                             U_CrseCoarse(id,ie) = crse(ic+1,jc+1,kc,ii)
                          end do

#if (AMREX_SPACEDIM == 3)
                          do id = 4*nNodesE+1, 5*nNodesE
                          ii   = ioff + (id-4*nNodesE-1)
                             U_CrseCoarse(id,ie) = crse(ic,jc,kc+1,ii)
                          end do

                          do id = 5*nNodesE+1, 6*nNodesE
                          ii   = ioff + (id-5*nNodesE-1)
                             U_CrseCoarse(id,ie) = crse(ic+1,jc,kc+1,ii)
                          end do

                          do id = 6*nNodesE+1, 7*nNodesE
                          ii   = ioff + (id-6*nNodesE-1)
                             U_CrseCoarse(id,ie) = crse(ic,jc+1,kc+1,ii)
                          end do

                          do id = 7*nNodesE+1, 8*nNodesE
                          ii   = ioff + (id-7*nNodesE-1)
                             U_CrseCoarse(id,ie) = crse(ic+1,jc+1,kc+1,ii)
                          end do
#endif

                   enddo ! ie

                   ! Coarsen U_FineCoarse --> U_CrseCoarse
                   call Refine_TwoMoment(nE, nFineX, U_CrseCoarse, U_FineCoarse)

                   do ie = 1, nE

                      ioff = (ispec-1)*(nCR*nE*nNodesE) + (imom-1)*(nE*nNodesE) + (ie-1)*nNodesE

                      do kcc = 1,nFineX_local(3)
                      do jcc = 1,nFineX_local(2)
                      do icc = 1,nFineX_local(1)

                         ifc = 2*(ic+icc-1)
                         jfc = 2*(jc+jcc-1)
                         kfc = 2*(kc+kcc-1)

                         do id = 1, nNodesE
                            ii   = ioff + (id-1)
                            fine(ifc,jfc,kfc,ii) = U_FineCoarse(id,ie,icc,jcc,kcc)
                         end do
         
                         do id = nNodesE+1, 2*nNodesE
                            ii   = ioff + (id-nNodesE-1)
                            fine(ifc+1,jfc,kfc,ii) = U_FineCoarse(id,ie,icc,jcc,kcc)
                         end do
      
                         do id = 2*nNodesE+1, 3*nNodesE
                            ii   = ioff + (id-2*nNodesE-1)
                            fine(ifc,jfc+1,kfc,ii) = U_FineCoarse(id,ie,icc,jcc,kcc)
                         end do
         
                         do id = 3*nNodesE+1, 4*nNodesE
                            ii   = ioff + (id-3*nNodesE-1)
                            fine(ifc+1,jfc+1,kfc,ii) = U_FineCoarse(id,ie,icc,jcc,kcc)
                         end do

#if (AMREX_SPACEDIM == 3)
                         do id = 4*nNodesE+1, 5*nNodesE
                            ii   = ioff + (id-4*nNodesE-1)
                            fine(ifc,jfc,kfc+1,ii) = U_FineCoarse(id,ie,icc,jcc,kcc)
                         end do

                         do id = 5*nNodesE+1, 6*nNodesE
                            ii   = ioff + (id-5*nNodesE-1)
                            fine(ifc+1,jfc,kfc+1,ii) = U_FineCoarse(id,ie,icc,jcc,kcc)
                         end do
      
                         do id = 6*nNodesE+1, 7*nNodesE
                            ii   = ioff + (id-6*nNodesE-1)
                            fine(ifc,jfc+1,kfc+1,ii) = U_FineCoarse(id,ie,icc,jcc,kcc)
                         end do
         
                         do id = 7*nNodesE+1, 8*nNodesE
                            ii   = ioff + (id-7*nNodesE-1)
                            fine(ifc+1,jfc+1,kfc+1,ii) = U_FineCoarse(id,ie,icc,jcc,kcc)
                         end do
#endif
                      end do
                      end do
                      end do

                   enddo ! ie

             enddo ! imom
          enddo ! ispec

       enddo ! ic
       enddo ! jc
    enddo ! kc

  end subroutine dg_refine

  subroutine dg_coarsen(cblo, cbhi, &
                        fine, fine_lo, fine_hi, &
                        crse, crse_lo, crse_hi, &
                        nFineX, nvar) &
                        bind(c,name='ca_dg_coarsen')

    use amrex_error_module, only: amrex_error
    use amrex_fort_module, only: rt => amrex_real
    use ProgramHeaderModule           , only: nNodesX, nNodesE, nE, nDOF
    use RadiationFieldsModule         , only: nSpecies, nCR
    use TwoMoment_MeshRefinementModule, only: Coarsen_TwoMoment    
    
    implicit none

    ! arguments
    integer :: crse_lo(3), crse_hi(3)
    integer :: fine_lo(3), fine_hi(3)
    integer :: cblo(3), cbhi(3)
    integer :: nFineX(3)
    integer :: nvar

    real(rt) :: fine(fine_lo(1):fine_hi(1), fine_lo(2):fine_hi(2), fine_lo(3):fine_hi(3), 0:nvar-1)
    real(rt) :: crse(crse_lo(1):crse_hi(1), crse_lo(2):crse_hi(2), crse_lo(3):crse_hi(3), 0:nvar-1)

    ! local variables
    integer :: ic , jc , kc
    integer :: icc, jcc, kcc
    integer :: ifc, jfc, kfc
    integer :: ispec, imom, id, ie, ii, ioff
    integer :: nFineX_local(3)

#if (AMREX_SPACEDIM == 2)
    real(rt) :: U_FineCoarse(1:nDOF, 1:nE, 1:nFineX(1), 1:nFineX(2), 1:1)
#elif (AMREX_SPACEDIM == 3)
    real(rt) :: U_FineCoarse(1:nDOF, 1:nE, 1:nFineX(1), 1:nFineX(2), 1:nFineX(3))
#endif
    real(rt) :: U_CrseCoarse(1:nDOF, 1:nE)

    ! The C++ sends in zero but Coarsen_TwoMoment expects 1
    nFineX_local(:) = nFineX(:)
#if (AMREX_SPACEDIM == 2)
    nFineX_local(3) = 1
#endif

    ! print *,'CRSE LO ', crse_lo(:)
    ! print *,'CRSE HI ', crse_hi(:)

    ! print *,'FINE LO ', fine_lo(:)
    ! print *,'FINE HI ', fine_hi(:)

    if ( nvar .ne. nNodesE * nE * nCR * nSpecies ) &
       call amrex_error("Mismatch of nvar with nNodesE*nE*nCR*nSpecies in dg_coarsen")
#if (AMREX_SPACEDIM == 2)
    if ( nDOF .ne. nNodesE * nNodesX(1) * nNodesX(2) ) &
       call amrex_error("Mismatch of nDOF with nNodesE*nNodesX(1)*nNodesX(2) in dg_refine")
#elif (AMREX_SPACEDIM == 3)
    if ( nDOF .ne. nNodesE * nNodesX(1) * nNodesX(2) * nNodesX(3) ) &
       call amrex_error("Mismatch of nDOF with nNodesE*nNodesX(1)*nNodesX(2)*nNodesX(3) in dg_refine")
#endif

    ! loop over coarse grid -- note that we must skip here because we need to construct an intermediate
    !      array at half of this resolution
#if (AMREX_SPACEDIM == 2)
    do kc = crse_lo(3), crse_hi(3)
#elif (AMREX_SPACEDIM == 3)
    do kc = cblo(3), cbhi(3), 2
#endif
    do jc = cblo(2), cbhi(2), 2
    do ic = cblo(1), cbhi(1), 2

          ! loop over species and moments ...
          do ispec = 1, nSpecies
             do imom = 1, nCR
                do ie = 1, nE

                   ioff = (ispec-1)*(nCR*nE*nNodesE) + (imom-1)*(nE*nNodesE) + (ie-1)*nNodesE

                   do kcc = 1,nFineX_local(3)
                   do jcc = 1,nFineX_local(2)
                   do icc = 1,nFineX_local(1)

                       ifc = 2*(ic+icc-1)
                       jfc = 2*(jc+jcc-1)
                       kfc = 2*(kc+kcc-1)
#if (AMREX_SPACEDIM == 2)
                       kfc = 0
#endif

                       do id = 1, nNodesE
                          ii = ioff + (id-1)
                          U_FineCoarse(id,ie,icc,jcc,kcc) = fine(ifc,jfc,kfc,ii)
                       end do

                       do id = nNodesE+1, 2*nNodesE
                          ii = ioff+ (id-nNodesE-1)
                          U_FineCoarse(id,ie,icc,jcc,kcc) = fine(ifc+1,jfc,kfc,ii)
                       end do

                       do id = 2*nNodesE+1, 3*nNodesE
                          ii = ioff + (id-2*nNodesE-1)
                          U_FineCoarse(id,ie,icc,jcc,kcc) = fine(ifc,jfc+1,kfc,ii)
                       end do

                       do id = 3*nNodesE+1, 4*nNodesE
                          ii = ioff + (id-3*nNodesE-1)
                          U_FineCoarse(id,ie,icc,jcc,kcc) = fine(ifc+1,jfc+1,kfc,ii)
                       end do

#if (AMREX_SPACEDIM == 3)
                       do id = 4*nNodesE+1, 5*nNodesE
                          ii = ioff + (id-4*nNodesE-1)
                          U_FineCoarse(id,ie,icc,jcc,kcc) = fine(ifc+1,jfc+1,kfc,ii)
                          U_FineCoarse(id,ie,icc,jcc,kcc) = fine(ifc,jfc,kfc+1,ii)
                       end do

                       do id = 5*nNodesE+1, 6*nNodesE
                          ii = ioff + (id-5*nNodesE-1)
                          U_FineCoarse(id,ie,icc,jcc,kcc) = fine(ifc+1,jfc,kfc+1,ii)
                       end do

                       do id = 6*nNodesE+1, 7*nNodesE
                          ii = ioff + (id-6*nNodesE-1)
                          U_FineCoarse(id,ie,icc,jcc,kcc) = fine(ifc,jfc+1,kfc+1,ii)
                       end do

                       do id = 7*nNodesE+1, 8*nNodesE
                          ii = ioff + (id-7*nNodesE-1)
                          U_FineCoarse(id,ie,icc,jcc,kcc) = fine(ifc+1,jfc+1,kfc+1,ii)
                       end do
#endif

                   end do ! icc
                   end do ! jcc
                   end do ! kcc

                enddo ! ie
                   
                ! Coarsen U_FineCoarse --> U_CrseCoarse
                call Coarsen_TwoMoment(nE, nFineX_local, U_CrseCoarse, U_FineCoarse)

                do ie = 1, nE

                   ioff = (ispec-1)*(nCR*nE*nNodesE) + (imom-1)*(nE*nNodesE) + (ie-1)*nNodesE

                   do id = 1, nNodesE
                      ii = ioff + (id-1)
                      crse(ic,jc,kc,ii) = U_CrseCoarse(id,ie)
                   end do
   
                   do id = nNodesE+1, 2*nNodesE
                      ii = ioff + (id-nNodesE-1)
                      crse(ic+1,jc,kc,ii) = U_CrseCoarse(id,ie)
                   end do
   
                   do id = 2*nNodesE+1, 3*nNodesE
                      ii = ioff + (id-2*nNodesE-1)
                      crse(ic,jc+1,kc,ii) = U_CrseCoarse(id,ie)
                   end do
   
                   do id = 3*nNodesE+1, 4*nNodesE
                      ii = ioff + (id-3*nNodesE-1)
                      crse(ic+1,jc+1,kc,ii) = U_CrseCoarse(id,ie)
                   end do
#if (AMREX_SPACEDIM == 3)
                   do id = 4*nNodesE+1, 5*nNodesE
                      ii = ioff + (id-4*nNodesE-1)
                      crse(ic,jc,kc+1,ii) = U_CrseCoarse(id,ie)
                   end do
   
                   do id = 5*nNodesE+1, 6*nNodesE
                      ii = ioff + (id-5*nNodesE-1)
                      crse(ic+1,jc,kc+1,ii) = U_CrseCoarse(id,ie)
                   end do
   
                   do id = 6*nNodesE+1, 7*nNodesE
                      ii = ioff + (id-6*nNodesE-1)
                      crse(ic,jc+1,kc+1,ii) = U_CrseCoarse(id,ie)
                   end do
   
                   do id = 7*nNodesE+1, 8*nNodesE
                      ii = ioff + (id-7*nNodesE-1)
                      crse(ic+1,jc+1,kc+1,ii) = U_CrseCoarse(id,ie)
                   end do
#endif

                enddo ! ie

             enddo ! imom
          enddo ! ispec

    enddo ! ic
    enddo ! jc
    enddo ! kc

  end subroutine dg_coarsen
  
end module dg_amr_helper_module
