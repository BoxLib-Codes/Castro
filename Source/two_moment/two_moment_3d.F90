  subroutine compute_thornado_timestep(dx, dt) &
                              bind(C, name="compute_thornado_timestep")

    use amrex_fort_module, only : rt => amrex_real
    use TimeSteppingModule_Castro, only : ComputeTimeStep_TwoMoment

    real(rt), intent(in ) :: dx(3)
    real(rt), intent(out) :: dt

    real(rt) :: dX_CGS(3)

    dX_CGS(:) = dx(:)

    call ComputeTimeStep_TwoMoment( dX_CGS, dt )

  end subroutine compute_thornado_timestep

  subroutine call_to_thornado(lo, hi, dt, &
                              S , s_lo, s_hi, ns , &
                              dS, d_lo, d_hi, nds, &
                              U_R_o, U_R_o_lo, U_R_o_hi, n_uro, &
                              dR   ,    dr_lo,    dr_hi, n_urn, &
                              n_moments, ng) &
                              bind(C, name="call_to_thornado")

    use amrex_constants_module, only : fourth, half, zero, one, two
    use amrex_fort_module, only : rt => amrex_real
    use amrex_error_module, only : amrex_abort
    use meth_params_module, only : URHO,UMX,UMY,UMZ,UEINT,UEDEN,UFX,UFS
    use ProgramHeaderModule, only : nE, nNodesE, swE
    use FluidFieldsModule, only : uCF, nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne
    use FluidFieldsModule, only : CreateFluidFields, DestroyFluidFields
    use GeometryFieldsModule, only : uGF
    use GeometryFieldsModuleE, only : uGE
    use RadiationFieldsModule, only : CreateRadiationFields,DestroyRadiationFields,nSpecies, uCR
    use TimeSteppingModule_Castro, only : Update_IMEX_PDARS
    use UnitsModule, only : Gram, Centimeter, Second,  AtomicMassUnit, Erg

    use SubcellReconstructionModule, only : ProjectionMatrix

    implicit none
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) ::  s_lo(3),  s_hi(3)
    integer, intent(in) ::  d_lo(3),  d_hi(3)
    integer, intent(in) ::  U_R_o_lo(3),  U_R_o_hi(3)
    integer, intent(in) ::     dr_lo(3),     dr_hi(3)
    integer, intent(in) ::  ns, nds, n_uro, n_urn
    integer, intent(in) ::  n_moments
    integer, intent(in) :: ng
    real(rt), intent(in) :: dt

    ! Here we expect  n_rad_comp = 20 x 16 x 6 x 4 (energy x dof x species x moments)

    ! Conserved fluid state (rho, rho u, rho v, rho w, rho E, rho e...)
    real(rt), intent(inout) ::  S(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),ns) 
    real(rt), intent(inout) :: dS(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nds)

    ! Old and new radiation state
    real(rt), intent(inout) ::  U_R_o(U_R_o_lo(1): U_R_o_hi(1),  U_R_o_lo(2): U_R_o_hi(2),   U_R_o_lo(3): U_R_o_hi(3), 0:n_uro-1) 
    real(rt), intent(inout) ::     dR(   dr_lo(1):    dr_hi(1),     dr_lo(2):    dr_hi(2),      dr_lo(3):    dr_hi(3), 0:n_urn-1)

    ! Temporary variables
    integer  :: i,j,k,ioff,nu
    integer  :: ic,jc,kc
    integer  :: ind
    integer  :: ii,id,ie,im,is
    integer  :: u_lo(3), u_hi(3)
    real(rt) :: conv_dens, conv_mom, conv_enr, conv_ne

    real(rt), allocatable :: u0(:,:,:,:,:)

    integer  :: nX(3)
    integer  :: swX(3)

    if (ng.ne. 2) &
      call amrex_abort("Need two ghost cells in call_to_thornado!")

    conv_dens = Gram / Centimeter**3
    conv_mom  = Gram / Centimeter**2 / Second
    conv_enr  = Gram / Centimeter / Second**2
    conv_enr  = Erg / Centimeter**3
    conv_ne   = 1.d0 / Centimeter**3

    nX(:)  = (hi(:) - lo(:) + 1) / 2
    swX(:) = ng

    u_lo(1) = 1    - swX(1)
    u_lo(2) = 1    - swX(2)
    u_lo(3) = 1    - swX(3)
    u_hi(1) = nX(1)+ swX(1)
    u_hi(2) = nX(2)+ swX(2)
    u_hi(3) = nX(3)+ swX(3)
    nu      = 8

    call CreateFluidFields ( nX, swX, Verbose_Option = .FALSE. )

    call CreateRadiationFields ( nX, swX, nE, swE, nSpecies_Option = nSpecies, Verbose_Option = .FALSE. )

    ! ************************************************************************************
    ! Interpolate from the Castro "S" arrays into Thornado "uCF" arrays
    ! ************************************************************************************

    allocate( u0(nu,u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nCF) )

    call interpolate_fluid (lo, hi, &
                            S , s_lo, s_hi, ns , &
                            u0, u_lo, u_hi, nu, &
                            ng)

    ! ************************************************************************************
    ! Copy from the Castro U_R arrays into Thornado arrays from InitThornado_Patch
    ! ************************************************************************************
    do kc = u_lo(3),u_hi(3)
    do jc = u_lo(2),u_hi(2)
    do ic = u_lo(1),u_hi(1)

         ! The uCR array was allocated in CreateRadiationdFields_Conserved with 
         ! ALLOCATE &
         !   ( uCR(1:nDOF, &
         !         1-swE:nE+swE, &
         !         1-swX(1):nX(1)+swX(1), &
         !         1-swX(2):nX(2)+swX(2), &
         !         1-swX(3):nX(3)+swX(3), &
         !         1:nCR, 1:nSpecies) )

         ! U_R_o spatial indices start at lo - (number of ghost zones)
         !   uCR spatial indices start at 1 - (number of ghost zones)
         i = lo(1) + 2*(ic-1)
         j = lo(2) + 2*(jc-1)
         k = lo(3) + 2*(kc-1)

         do is = 1, nSpecies
         do im = 1, n_moments
         do ie = 1, nE

         ioff = (is-1)*(n_moments*nE*nNodesE) + (im-1)*(nE*nNodesE) + (ie-1)*nNodesE

         do id = 1, nNodesE
            ii   = ioff + (id-1)
            uCR(id,ie,ic,jc,kc,im,is) = U_R_o(i,j,k,ii)
         end do

         do id = nNodesE+1, 2*nNodesE
            ii   = ioff + (id-nNodesE-1)
            uCR(id,ie,ic,jc,kc,im,is) = U_R_o(i+1,j,k,ii)
         end do

         do id = 2*nNodesE+1, 3*nNodesE
            ii   = ioff + (id-2*nNodesE-1)
            uCR(id,ie,ic,jc,kc,im,is) = U_R_o(i,j+1,k,ii)
         end do

         do id = 3*nNodesE+1, 4*nNodesE
            ii   = ioff + (id-3*nNodesE-1)
            uCR(id,ie,ic,jc,kc,im,is) = U_R_o(i+1,j+1,k,ii)
         end do

         do id = 4*nNodesE+1,5*nNodesE
            ii   = ioff + (id-4*nNodesE-1)
            uCR(id,ie,ic,jc,kc,im,is) = U_R_o(i,j,k+1,ii)
         end do

         do id = 5*nNodesE+1, 6*nNodesE
            ii   = ioff + (id-5*nNodesE-1)
            uCR(id,ie,ic,jc,kc,im,is) = U_R_o(i+1,j,k+1,ii)
         end do

         do id = 6*nNodesE+1, 7*nNodesE
            ii   = ioff + (id-6*nNodesE-1)
            uCR(id,ie,ic,jc,kc,im,is) = U_R_o(i,j+1,k+1,ii)
         end do

         do id = 7*nNodesE+1, 8*nNodesE
            ii   = ioff + (id-7*nNodesE-1)
            uCR(id,ie,ic,jc,kc,im,is) = U_R_o(i+1,j+1,k+1,ii)
         end do

         end do
         end do
         end do

    end do
    end do
    end do

    ! ************************************************************************************
    ! Call the time stepper that lives in the thornado repo
    ! ************************************************************************************

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: uCF, uCR, uGE, uGF )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( uCF, uCR, uGE, uGF )
#endif

    call Update_IMEX_PDARS(dt*Second, uCF, uCR)

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: uCF, uCR ) &
    !$OMP MAP( release: uGE, uGF )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT( uCF, uCR ) &
    !$ACC DELETE( uGE, uGF )
#endif

    ! Zero out dS so we can accumulate weighted average in it
     dS(:,:,:,:) = 0.d0

    ! ************************************************************************************
    ! Copy back from the thornado arrays into Castro arrays
    ! ************************************************************************************
    do kc = 1, nX(3)
    do jc = 1, nX(2)
    do ic = 1, nX(1)

         ! We store dS as a source term which we can add to S outside of this routine
         ! We now use the weighting from thornado to convert the four node values back
         !  into a single averaged value
         !
         ! Update_IMEX_PC2 doesn't currently change the fluid density or momentum
         !
         !   S spatial indices start at lo - (number of ghost zones)
         ! uCF spatial indices start at 1  - (number of ghost zones)
         i = lo(1) + 2*(ic-1)
         j = lo(2) + 2*(jc-1)
         k = lo(3) + 2*(kc-1)

!!$         dS(i  ,j  ,k  ,URHO ) = dS(i  ,j  ,k  ,URHO ) + (uCF(1,ic,jc,kc,iCF_D ) - u0(1,ic,jc,kc,iCF_D ) )
!!$         dS(i  ,j  ,k  ,UMX  ) = dS(i  ,j  ,k  ,UMX  ) + (uCF(1,ic,jc,kc,iCF_S1) - u0(1,ic,jc,kc,iCF_S1) )
!!$         dS(i  ,j  ,k  ,UMY  ) = dS(i  ,j  ,k  ,UMY  ) + (uCF(1,ic,jc,kc,iCF_S2) - u0(1,ic,jc,kc,iCF_S2) )
!!$         dS(i  ,j  ,k  ,UMZ  ) = dS(i  ,j  ,k  ,UMZ  ) + (uCF(1,ic,jc,kc,iCF_S3) - u0(1,ic,jc,kc,iCF_S3) )
!!$         dS(i  ,j  ,k  ,UEDEN) = dS(i  ,j  ,k  ,UEDEN) + (uCF(1,ic,jc,kc,iCF_E ) - u0(1,ic,jc,kc,iCF_E ) )
!!$         dS(i  ,j  ,k  ,UFX  ) = dS(i  ,j  ,k  ,UFX  ) + (uCF(1,ic,jc,kc,iCF_Ne) - u0(1,ic,jc,kc,iCF_Ne) )
!!$
!!$         dS(i+1,j  ,k  ,URHO ) = dS(i+1,j  ,k  ,URHO ) + (uCF(2,ic,jc,kc,iCF_D ) - u0(2,ic,jc,kc,iCF_D ) )
!!$         dS(i+1,j  ,k  ,UMX  ) = dS(i+1,j  ,k  ,UMX  ) + (uCF(2,ic,jc,kc,iCF_S1) - u0(2,ic,jc,kc,iCF_S1) )
!!$         dS(i+1,j  ,k  ,UMY  ) = dS(i+1,j  ,k  ,UMY  ) + (uCF(2,ic,jc,kc,iCF_S2) - u0(2,ic,jc,kc,iCF_S2) )
!!$         dS(i+1,j  ,k  ,UMZ  ) = dS(i+1,j  ,k  ,UMZ  ) + (uCF(2,ic,jc,kc,iCF_S3) - u0(2,ic,jc,kc,iCF_S3) )
!!$         dS(i+1,j  ,k  ,UEDEN) = dS(i+1,j  ,k  ,UEDEN) + (uCF(2,ic,jc,kc,iCF_E ) - u0(2,ic,jc,kc,iCF_E ) )
!!$         dS(i+1,j  ,k  ,UFX  ) = dS(i+1,j  ,k  ,UFX  ) + (uCF(2,ic,jc,kc,iCF_Ne) - u0(2,ic,jc,kc,iCF_Ne) )
!!$
!!$         dS(i  ,j+1,k  ,URHO ) = dS(i  ,j+1,k  ,URHO ) + (uCF(3,ic,jc,kc,iCF_D ) - u0(3,ic,jc,kc,iCF_D ) )
!!$         dS(i  ,j+1,k  ,UMX  ) = dS(i  ,j+1,k  ,UMX  ) + (uCF(3,ic,jc,kc,iCF_S1) - u0(3,ic,jc,kc,iCF_S1) )
!!$         dS(i  ,j+1,k  ,UMY  ) = dS(i  ,j+1,k  ,UMY  ) + (uCF(3,ic,jc,kc,iCF_S2) - u0(3,ic,jc,kc,iCF_S2) )
!!$         dS(i  ,j+1,k  ,UMZ  ) = dS(i  ,j+1,k  ,UMZ  ) + (uCF(3,ic,jc,kc,iCF_S3) - u0(3,ic,jc,kc,iCF_S3) )
!!$         dS(i  ,j+1,k  ,UEDEN) = dS(i  ,j+1,k  ,UEDEN) + (uCF(3,ic,jc,kc,iCF_E ) - u0(3,ic,jc,kc,iCF_E ) )
!!$         dS(i  ,j+1,k  ,UFX  ) = dS(i  ,j+1,k  ,UFX  ) + (uCF(3,ic,jc,kc,iCF_Ne) - u0(3,ic,jc,kc,iCF_Ne) )
!!$
!!$         dS(i+1,j+1,k  ,URHO ) = dS(i+1,j+1,k  ,URHO ) + (uCF(4,ic,jc,kc,iCF_D ) - u0(4,ic,jc,kc,iCF_D ) )
!!$         dS(i+1,j+1,k  ,UMX  ) = dS(i+1,j+1,k  ,UMX  ) + (uCF(4,ic,jc,kc,iCF_S1) - u0(4,ic,jc,kc,iCF_S1) )
!!$         dS(i+1,j+1,k  ,UMY  ) = dS(i+1,j+1,k  ,UMY  ) + (uCF(4,ic,jc,kc,iCF_S2) - u0(4,ic,jc,kc,iCF_S2) )
!!$         dS(i+1,j+1,k  ,UMZ  ) = dS(i+1,j+1,k  ,UMZ  ) + (uCF(4,ic,jc,kc,iCF_S3) - u0(4,ic,jc,kc,iCF_S3) )
!!$         dS(i+1,j+1,k  ,UEDEN) = dS(i+1,j+1,k  ,UEDEN) + (uCF(4,ic,jc,kc,iCF_E ) - u0(4,ic,jc,kc,iCF_E ) )
!!$         dS(i+1,j+1,k  ,UFX  ) = dS(i+1,j+1,k  ,UFX  ) + (uCF(4,ic,jc,kc,iCF_Ne) - u0(4,ic,jc,kc,iCF_Ne) )
!!$
!!$         dS(i  ,j  ,k+1,URHO ) = dS(i  ,j  ,k+1,URHO ) + (uCF(5,ic,jc,kc,iCF_D ) - u0(5,ic,jc,kc,iCF_D ) )
!!$         dS(i  ,j  ,k+1,UMX  ) = dS(i  ,j  ,k+1,UMX  ) + (uCF(5,ic,jc,kc,iCF_S1) - u0(5,ic,jc,kc,iCF_S1) )
!!$         dS(i  ,j  ,k+1,UMY  ) = dS(i  ,j  ,k+1,UMY  ) + (uCF(5,ic,jc,kc,iCF_S2) - u0(5,ic,jc,kc,iCF_S2) )
!!$         dS(i  ,j  ,k+1,UMZ  ) = dS(i  ,j  ,k+1,UMZ  ) + (uCF(5,ic,jc,kc,iCF_S3) - u0(5,ic,jc,kc,iCF_S3) )
!!$         dS(i  ,j  ,k+1,UEDEN) = dS(i  ,j  ,k+1,UEDEN) + (uCF(5,ic,jc,kc,iCF_E ) - u0(5,ic,jc,kc,iCF_E ) )
!!$         dS(i  ,j  ,k+1,UFX  ) = dS(i  ,j  ,k+1,UFX  ) + (uCF(5,ic,jc,kc,iCF_Ne) - u0(5,ic,jc,kc,iCF_Ne) )
!!$
!!$         dS(i+1,j  ,k+1,URHO ) = dS(i+1,j  ,k+1,URHO ) + (uCF(6,ic,jc,kc,iCF_D ) - u0(6,ic,jc,kc,iCF_D ) )
!!$         dS(i+1,j  ,k+1,UMX  ) = dS(i+1,j  ,k+1,UMX  ) + (uCF(6,ic,jc,kc,iCF_S1) - u0(6,ic,jc,kc,iCF_S1) )
!!$         dS(i+1,j  ,k+1,UMY  ) = dS(i+1,j  ,k+1,UMY  ) + (uCF(6,ic,jc,kc,iCF_S2) - u0(6,ic,jc,kc,iCF_S2) )
!!$         dS(i+1,j  ,k+1,UMZ  ) = dS(i+1,j  ,k+1,UMZ  ) + (uCF(6,ic,jc,kc,iCF_S3) - u0(6,ic,jc,kc,iCF_S3) )
!!$         dS(i+1,j  ,k+1,UEDEN) = dS(i+1,j  ,k+1,UEDEN) + (uCF(6,ic,jc,kc,iCF_E ) - u0(6,ic,jc,kc,iCF_E ) )
!!$         dS(i+1,j  ,k+1,UFX  ) = dS(i+1,j  ,k+1,UFX  ) + (uCF(6,ic,jc,kc,iCF_Ne) - u0(6,ic,jc,kc,iCF_Ne) )
!!$
!!$         dS(i  ,j+1,k+1,URHO ) = dS(i  ,j+1,k+1,URHO ) + (uCF(7,ic,jc,kc,iCF_D ) - u0(7,ic,jc,kc,iCF_D ) )
!!$         dS(i  ,j+1,k+1,UMX  ) = dS(i  ,j+1,k+1,UMX  ) + (uCF(7,ic,jc,kc,iCF_S1) - u0(7,ic,jc,kc,iCF_S1) )
!!$         dS(i  ,j+1,k+1,UMY  ) = dS(i  ,j+1,k+1,UMY  ) + (uCF(7,ic,jc,kc,iCF_S2) - u0(7,ic,jc,kc,iCF_S2) )
!!$         dS(i  ,j+1,k+1,UMZ  ) = dS(i  ,j+1,k+1,UMZ  ) + (uCF(7,ic,jc,kc,iCF_S3) - u0(7,ic,jc,kc,iCF_S3) )
!!$         dS(i  ,j+1,k+1,UEDEN) = dS(i  ,j+1,k+1,UEDEN) + (uCF(7,ic,jc,kc,iCF_E ) - u0(7,ic,jc,kc,iCF_E ) )
!!$         dS(i  ,j+1,k+1,UFX  ) = dS(i  ,j+1,k+1,UFX  ) + (uCF(7,ic,jc,kc,iCF_Ne) - u0(7,ic,jc,kc,iCF_Ne) )
!!$
!!$         dS(i+1,j+1,k+1,URHO ) = dS(i+1,j+1,k+1,URHO ) + (uCF(8,ic,jc,kc,iCF_D ) - u0(8,ic,jc,kc,iCF_D ) )
!!$         dS(i+1,j+1,k+1,UMX  ) = dS(i+1,j+1,k+1,UMX  ) + (uCF(8,ic,jc,kc,iCF_S1) - u0(8,ic,jc,kc,iCF_S1) )
!!$         dS(i+1,j+1,k+1,UMY  ) = dS(i+1,j+1,k+1,UMY  ) + (uCF(8,ic,jc,kc,iCF_S2) - u0(8,ic,jc,kc,iCF_S2) )
!!$         dS(i+1,j+1,k+1,UMZ  ) = dS(i+1,j+1,k+1,UMZ  ) + (uCF(8,ic,jc,kc,iCF_S3) - u0(8,ic,jc,kc,iCF_S3) )
!!$         dS(i+1,j+1,k+1,UEDEN) = dS(i+1,j+1,k+1,UEDEN) + (uCF(8,ic,jc,kc,iCF_E ) - u0(8,ic,jc,kc,iCF_E ) )
!!$         dS(i+1,j+1,k+1,UFX  ) = dS(i+1,j+1,k+1,UFX  ) + (uCF(8,ic,jc,kc,iCF_Ne) - u0(8,ic,jc,kc,iCF_Ne) )

         do ind = 1, 8

            dS(i  ,j  , k  ,URHO ) = dS(i  ,j  ,k  ,URHO ) &
              + ProjectionMatrix(1,ind) * ( uCF(ind,ic,jc,kc,iCF_D ) - u0(ind,ic,jc,kc,iCF_D ) )
            dS(i  ,j  , k  ,UMX  ) = dS(i  ,j  ,k  ,UMX  ) &
              + ProjectionMatrix(1,ind) * ( uCF(ind,ic,jc,kc,iCF_S1) - u0(ind,ic,jc,kc,iCF_S1) )
            dS(i  ,j  , k  ,UMY  ) = dS(i  ,j  ,k  ,UMY  ) &
              + ProjectionMatrix(1,ind) * ( uCF(ind,ic,jc,kc,iCF_S2) - u0(ind,ic,jc,kc,iCF_S2) )
            dS(i  ,j  , k  ,UMZ  ) = dS(i  ,j  ,k  ,UMZ  ) &
              + ProjectionMatrix(1,ind) * ( uCF(ind,ic,jc,kc,iCF_S3) - u0(ind,ic,jc,kc,iCF_S3) )
            dS(i  ,j  , k  ,UEDEN) = dS(i  ,j  ,k  ,UEDEN) &
              + ProjectionMatrix(1,ind) * ( uCF(ind,ic,jc,kc,iCF_E ) - u0(ind,ic,jc,kc,iCF_E ) )
            dS(i  ,j  , k  ,UFX  ) = dS(i  ,j  ,k  ,UFX  ) &
              + ProjectionMatrix(1,ind) * ( uCF(ind,ic,jc,kc,iCF_Ne) - u0(ind,ic,jc,kc,iCF_Ne) )

            dS(i+1,j  , k  ,URHO ) = dS(i+1,j  ,k  ,URHO ) &
              + ProjectionMatrix(2,ind) * ( uCF(ind,ic,jc,kc,iCF_D ) - u0(ind,ic,jc,kc,iCF_D ) )
            dS(i+1,j  , k  ,UMX  ) = dS(i+1,j  ,k  ,UMX  ) &
              + ProjectionMatrix(2,ind) * ( uCF(ind,ic,jc,kc,iCF_S1) - u0(ind,ic,jc,kc,iCF_S1) )
            dS(i+1,j  , k  ,UMY  ) = dS(i+1,j  ,k  ,UMY  ) &
              + ProjectionMatrix(2,ind) * ( uCF(ind,ic,jc,kc,iCF_S2) - u0(ind,ic,jc,kc,iCF_S2) )
            dS(i+1,j  , k  ,UMZ  ) = dS(i+1,j  ,k  ,UMZ  ) &
              + ProjectionMatrix(2,ind) * ( uCF(ind,ic,jc,kc,iCF_S3) - u0(ind,ic,jc,kc,iCF_S3) )
            dS(i+1,j  , k  ,UEDEN) = dS(i+1,j  ,k  ,UEDEN) &
              + ProjectionMatrix(2,ind) * ( uCF(ind,ic,jc,kc,iCF_E ) - u0(ind,ic,jc,kc,iCF_E ) )
            dS(i+1,j  , k  ,UFX  ) = dS(i+1,j  ,k  ,UFX  ) &
              + ProjectionMatrix(2,ind) * ( uCF(ind,ic,jc,kc,iCF_Ne) - u0(ind,ic,jc,kc,iCF_Ne) )

            dS(i  ,j+1, k  ,URHO ) = dS(i  ,j+1,k  ,URHO ) &
              + ProjectionMatrix(3,ind) * ( uCF(ind,ic,jc,kc,iCF_D ) - u0(ind,ic,jc,kc,iCF_D ) )
            dS(i  ,j+1, k  ,UMX  ) = dS(i  ,j+1,k  ,UMX  ) &
              + ProjectionMatrix(3,ind) * ( uCF(ind,ic,jc,kc,iCF_S1) - u0(ind,ic,jc,kc,iCF_S1) )
            dS(i  ,j+1, k  ,UMY  ) = dS(i  ,j+1,k  ,UMY  ) &
              + ProjectionMatrix(3,ind) * ( uCF(ind,ic,jc,kc,iCF_S2) - u0(ind,ic,jc,kc,iCF_S2) )
            dS(i  ,j+1, k  ,UMZ  ) = dS(i  ,j+1,k  ,UMZ  ) &
              + ProjectionMatrix(3,ind) * ( uCF(ind,ic,jc,kc,iCF_S3) - u0(ind,ic,jc,kc,iCF_S3) )
            dS(i  ,j+1, k  ,UEDEN) = dS(i  ,j+1,k  ,UEDEN) &
              + ProjectionMatrix(3,ind) * ( uCF(ind,ic,jc,kc,iCF_E ) - u0(ind,ic,jc,kc,iCF_E ) )
            dS(i  ,j+1, k  ,UFX  ) = dS(i  ,j+1,k  ,UFX  ) &
              + ProjectionMatrix(3,ind) * ( uCF(ind,ic,jc,kc,iCF_Ne) - u0(ind,ic,jc,kc,iCF_Ne) )

            dS(i+1,j+1, k  ,URHO ) = dS(i+1,j+1,k  ,URHO ) &
              + ProjectionMatrix(4,ind) * ( uCF(ind,ic,jc,kc,iCF_D ) - u0(ind,ic,jc,kc,iCF_D ) )
            dS(i+1,j+1, k  ,UMX  ) = dS(i+1,j+1,k  ,UMX  ) &
              + ProjectionMatrix(4,ind) * ( uCF(ind,ic,jc,kc,iCF_S1) - u0(ind,ic,jc,kc,iCF_S1) )
            dS(i+1,j+1, k  ,UMY  ) = dS(i+1,j+1,k  ,UMY  ) &
              + ProjectionMatrix(4,ind) * ( uCF(ind,ic,jc,kc,iCF_S2) - u0(ind,ic,jc,kc,iCF_S2) )
            dS(i+1,j+1, k  ,UMZ  ) = dS(i+1,j+1,k  ,UMZ  ) &
              + ProjectionMatrix(4,ind) * ( uCF(ind,ic,jc,kc,iCF_S3) - u0(ind,ic,jc,kc,iCF_S3) )
            dS(i+1,j+1, k  ,UEDEN) = dS(i+1,j+1,k  ,UEDEN) &
              + ProjectionMatrix(4,ind) * ( uCF(ind,ic,jc,kc,iCF_E ) - u0(ind,ic,jc,kc,iCF_E ) )
            dS(i+1,j+1, k  ,UFX  ) = dS(i+1,j+1,k  ,UFX  ) &
              + ProjectionMatrix(4,ind) * ( uCF(ind,ic,jc,kc,iCF_Ne) - u0(ind,ic,jc,kc,iCF_Ne) )

            dS(i  ,j  , k+1,URHO ) = dS(i  ,j  ,k+1,URHO ) &
              + ProjectionMatrix(5,ind) * ( uCF(ind,ic,jc,kc,iCF_D ) - u0(ind,ic,jc,kc,iCF_D ) )
            dS(i  ,j  , k+1,UMX  ) = dS(i  ,j  ,k+1,UMX  ) &
              + ProjectionMatrix(5,ind) * ( uCF(ind,ic,jc,kc,iCF_S1) - u0(ind,ic,jc,kc,iCF_S1) )
            dS(i  ,j  , k+1,UMY  ) = dS(i  ,j  ,k+1,UMY  ) &
              + ProjectionMatrix(5,ind) * ( uCF(ind,ic,jc,kc,iCF_S2) - u0(ind,ic,jc,kc,iCF_S2) )
            dS(i  ,j  , k+1,UMZ  ) = dS(i  ,j  ,k+1,UMZ  ) &
              + ProjectionMatrix(5,ind) * ( uCF(ind,ic,jc,kc,iCF_S3) - u0(ind,ic,jc,kc,iCF_S3) )
            dS(i  ,j  , k+1,UEDEN) = dS(i  ,j  ,k+1,UEDEN) &
              + ProjectionMatrix(5,ind) * ( uCF(ind,ic,jc,kc,iCF_E ) - u0(ind,ic,jc,kc,iCF_E ) )
            dS(i  ,j  , k+1,UFX  ) = dS(i  ,j  ,k+1,UFX  ) &
              + ProjectionMatrix(5,ind) * ( uCF(ind,ic,jc,kc,iCF_Ne) - u0(ind,ic,jc,kc,iCF_Ne) )

            dS(i+1,j  , k+1,URHO ) = dS(i+1,j  ,k+1,URHO ) &
              + ProjectionMatrix(6,ind) * ( uCF(ind,ic,jc,kc,iCF_D ) - u0(ind,ic,jc,kc,iCF_D ) )
            dS(i+1,j  , k+1,UMX  ) = dS(i+1,j  ,k+1,UMX  ) &
              + ProjectionMatrix(6,ind) * ( uCF(ind,ic,jc,kc,iCF_S1) - u0(ind,ic,jc,kc,iCF_S1) )
            dS(i+1,j  , k+1,UMY  ) = dS(i+1,j  ,k+1,UMY  ) &
              + ProjectionMatrix(6,ind) * ( uCF(ind,ic,jc,kc,iCF_S2) - u0(ind,ic,jc,kc,iCF_S2) )
            dS(i+1,j  , k+1,UMZ  ) = dS(i+1,j  ,k+1,UMZ  ) &
              + ProjectionMatrix(6,ind) * ( uCF(ind,ic,jc,kc,iCF_S3) - u0(ind,ic,jc,kc,iCF_S3) )
            dS(i+1,j  , k+1,UEDEN) = dS(i+1,j  ,k+1,UEDEN) &
              + ProjectionMatrix(6,ind) * ( uCF(ind,ic,jc,kc,iCF_E ) - u0(ind,ic,jc,kc,iCF_E ) )
            dS(i+1,j  , k+1,UFX  ) = dS(i+1,j  ,k+1,UFX  ) &
              + ProjectionMatrix(6,ind) * ( uCF(ind,ic,jc,kc,iCF_Ne) - u0(ind,ic,jc,kc,iCF_Ne) )

            dS(i  ,j+1, k+1,URHO ) = dS(i  ,j+1,k+1,URHO ) &
              + ProjectionMatrix(7,ind) * ( uCF(ind,ic,jc,kc,iCF_D ) - u0(ind,ic,jc,kc,iCF_D ) )
            dS(i  ,j+1, k+1,UMX  ) = dS(i  ,j+1,k+1,UMX  ) &
              + ProjectionMatrix(7,ind) * ( uCF(ind,ic,jc,kc,iCF_S1) - u0(ind,ic,jc,kc,iCF_S1) )
            dS(i  ,j+1, k+1,UMY  ) = dS(i  ,j+1,k+1,UMY  ) &
              + ProjectionMatrix(7,ind) * ( uCF(ind,ic,jc,kc,iCF_S2) - u0(ind,ic,jc,kc,iCF_S2) )
            dS(i  ,j+1, k+1,UMZ ) = dS(i  ,j+1,k+1,UMZ   ) &
              + ProjectionMatrix(7,ind) * ( uCF(ind,ic,jc,kc,iCF_S3) - u0(ind,ic,jc,kc,iCF_S3) )
            dS(i  ,j+1, k+1,UEDEN) = dS(i  ,j+1,k+1,UEDEN) &
              + ProjectionMatrix(7,ind) * ( uCF(ind,ic,jc,kc,iCF_E ) - u0(ind,ic,jc,kc,iCF_E ) )
            dS(i  ,j+1, k+1,UFX  ) = dS(i  ,j+1,k+1,UFX  ) &
              + ProjectionMatrix(7,ind) * ( uCF(ind,ic,jc,kc,iCF_Ne) - u0(ind,ic,jc,kc,iCF_Ne) )

            dS(i+1,j+1, k+1,URHO ) = dS(i+1,j+1,k+1,URHO ) &
              + ProjectionMatrix(8,ind) * ( uCF(ind,ic,jc,kc,iCF_D ) - u0(ind,ic,jc,kc,iCF_D ) )
            dS(i+1,j+1, k+1,UMX  ) = dS(i+1,j+1,k+1,UMX  ) &
              + ProjectionMatrix(8,ind) * ( uCF(ind,ic,jc,kc,iCF_S1) - u0(ind,ic,jc,kc,iCF_S1) )
            dS(i+1,j+1, k+1,UMY  ) = dS(i+1,j+1,k+1,UMY  ) &
              + ProjectionMatrix(8,ind) * ( uCF(ind,ic,jc,kc,iCF_S2) - u0(ind,ic,jc,kc,iCF_S2) )
            dS(i+1,j+1, k+1,UMZ  ) = dS(i+1,j+1,k+1,UMZ  ) &
              + ProjectionMatrix(8,ind) * ( uCF(ind,ic,jc,kc,iCF_S3) - u0(ind,ic,jc,kc,iCF_S3) )
            dS(i+1,j+1, k+1,UEDEN) = dS(i+1,j+1,k+1,UEDEN) &
              + ProjectionMatrix(8,ind) * ( uCF(ind,ic,jc,kc,iCF_E ) - u0(ind,ic,jc,kc,iCF_E ) )
            dS(i+1,j+1, k+1,UFX  ) = dS(i+1,j+1,k+1,UFX  ) &
              + ProjectionMatrix(8,ind) * ( uCF(ind,ic,jc,kc,iCF_Ne) - u0(ind,ic,jc,kc,iCF_Ne) )

         end do

    end do
    end do
    end do

    do k = lo(3),hi(3)
    do j = lo(2),hi(2)
    do i = lo(1),hi(1)

!        dS(i,j,k,URHO ) = dS(i,j,k,URHO ) / conv_dens
!        dS(i,j,k,UMX  ) = dS(i,j,k,UMX  ) / conv_mom
!        dS(i,j,k,UMY  ) = dS(i,j,k,UMY  ) / conv_mom
!        dS(i,j,k,UMZ  ) = dS(i,j,k,UMZ  ) / conv_mom
         dS(i,j,k,UEDEN) = dS(i,j,k,UEDEN) / conv_enr
         dS(i,j,k,UFX  ) = dS(i,j,k,UFX  ) / conv_ne

         dS(i,j,k,UEINT) = dS(i,j,k,UEDEN)     ! TRUE IFF NO MOMENTUM SOURCE TERMS

         ! Store electron molar fraction * density in the species
         dS(i,j,k,UFS)   = dS(i,j,k,UFX) * AtomicMassUnit / Gram

    end do
    end do
    end do

    do kc = 1, nX(3)
    do jc = 1, nX(2)
    do ic = 1, nX(1)

         ! U_R_o spatial indices start at lo - (number of ghost zones)
         !   uCR spatial indices start at 1 - (number of ghost zones)
         i = lo(1) + 2*(ic-1)
         j = lo(2) + 2*(jc-1)
         k = lo(3) + 2*(kc-1)

         do is = 1, nSpecies
         do im = 1, n_moments
         do ie = 1, nE

         ioff = (is-1)*(n_moments*nE*nNodesE) + (im-1)*(nE*nNodesE) + (ie-1)*nNodesE

         do id = 1, nNodesE
            ii   = ioff + (id-1)
            dR(i,j,k,ii) = uCR(id,ie,ic,jc,kc,im,is) - U_R_o(i,j,k,ii)
         end do

         do id = nNodesE+1, 2*nNodesE
            ii   = ioff + (id-nNodesE-1)
            dR(i+1,j,k,ii) = uCR(id,ie,ic,jc,kc,im,is) - U_R_o(i+1,j,k,ii)
         end do

         do id = 2*nNodesE+1, 3*nNodesE
            ii   = ioff + (id-2*nNodesE-1)
            dR(i,j+1,k,ii) = uCR(id,ie,ic,jc,kc,im,is) - U_R_o(i,j+1,k,ii)
         end do

         do id = 3*nNodesE+1, 4*nNodesE
            ii   = ioff + (id-3*nNodesE-1)
            dR(i+1,j+1,k,ii) = uCR(id,ie,ic,jc,kc,im,is) - U_R_o(i+1,j+1,k,ii)
         end do

         do id = 4*nNodesE+1,5*nNodesE
            ii   = ioff + (id-4*nNodesE-1)
            dR(i,j,k+1,ii) = uCR(id,ie,ic,jc,kc,im,is) - U_R_o(i,j,k+1,ii)
         end do

         do id = 5*nNodesE+1, 6*nNodesE
            ii   = ioff + (id-5*nNodesE-1)
            dR(i+1,j,k+1,ii) = uCR(id,ie,ic,jc,kc,im,is) - U_R_o(i+1,j,k+1,ii)
         end do

         do id = 6*nNodesE+1, 7*nNodesE
            ii   = ioff + (id-6*nNodesE-1)
            dR(i,j+1,k+1,ii) = uCR(id,ie,ic,jc,kc,im,is) - U_R_o(i,j+1,k+1,ii)
         end do

         do id = 7*nNodesE+1, 8*nNodesE
            ii   = ioff + (id-7*nNodesE-1)
            dR(i+1,j+1,k+1,ii) = uCR(id,ie,ic,jc,kc,im,is) - U_R_o(i+1,j+1,k+1,ii)
         end do

         end do
         end do
         end do

    end do
    end do
    end do

    deallocate(u0)

    call DestroyFluidFields
    call DestroyRadiationFields

  end subroutine call_to_thornado
