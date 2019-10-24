
  subroutine interpolate_fluid (lo, hi, &
                                S , s_lo, s_hi, ns , &
                                u0, u_lo, u_hi, nu, &
                                ng) 

    use amrex_constants_module, only : fourth, half, zero, one, two
    use amrex_fort_module, only : rt => amrex_real
    use amrex_error_module, only : amrex_abort
    use meth_params_module, only : URHO,UMX,UMY,UMZ,UEDEN,UFX
    use SubcellReconstructionModule, only: ReconstructionMatrix
    use FluidFieldsModule    , only : uCF, nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne
    use EquationOfStateModule_TABLE, only : ComputeTemperatureFromSpecificInternalEnergy_TABLE, &
                                            ComputeThermodynamicStates_Primitive_TABLE
    use UnitsModule, only : Gram, Centimeter, Second, AtomicMassUnit, Erg, Kelvin

    use ReferenceElementModuleX, only: NodesX_q

    implicit none
    integer, intent(in) :: lo(2), hi(2)
    integer, intent(in) :: s_lo(2),  s_hi(2)
    integer, intent(in) :: u_lo(3),  u_hi(3)
    integer, intent(in) :: ns, nu
    integer, intent(in) :: ng

    ! Conserved fluid state (rho, rho u, rho v, rho w, rho E, rho e...)
    real(rt), intent(inout) ::  S(s_lo(1):s_hi(1),s_lo(2):s_hi(2),ns) 

    real(rt), intent(inout) :: u0(nu,u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nCF)

    ! Temporary variables
    integer  :: i,j
    integer  :: ic,jc,kc
    integer  :: ind
    real(rt) :: conv_dens, conv_mom, conv_enr, conv_ne

    real(rt) :: x(4)
    real(rt) :: y(4)
    real(rt) :: z(4)

    integer  :: interp_type

    ! No interpolation:
    !interp_type = 0
    ! Use reconstruction matrix:
    interp_type = 1

    if (s_lo(1) .gt. (lo(1)-ng-1) .or. s_hi(1) .lt. (hi(1)+ng+1) .or. &
        s_lo(2) .gt. (lo(2)-ng-1) .or. s_hi(2) .lt. (hi(2)+ng+1)) then
       print *,'WRONG NUMBER OF GHOST CELLS!'
       stop
    end if

    ! Conversion to thornado units
    conv_dens = Gram / Centimeter**3
    conv_mom  = Gram / Centimeter**2 / Second
    conv_enr  = Erg / Centimeter**3
    conv_ne   = 1.d0 / Centimeter**3

    ! These are the locations of the DG nodes in the space [-.5:.5]
    ! Note that in 2D with 2 nodes in each direction, they are ordered 
    ! ind = 1:  (Lo,Lo)
    ! ind = 2:  (Hi,Lo)
    ! ind = 3:  (Lo,Hi)
    ! ind = 4:  (Hi,Hi)
     do ind = 1, 4
        x(ind) = NodesX_q(1,ind)
        y(ind) = NodesX_q(2,ind)
        z(ind) = NodesX_q(3,ind)
     end do

     ! The uCF array was allocated in CreateFluidFieldsConserved with 
     !     ALLOCATE( uCF &
     !      (1:nDOFX, &
     !       1-swX(1):nX(1)+swX(1), &
     !       1-swX(2):nX(2)+swX(2), &
     !       1-swX(3):nX(3)+swX(3), &
     !       1:nCF) )

    ! ************************************************************************************
    ! No interpolation
    ! ************************************************************************************
    if (interp_type .eq. 0) then
       do jc = u_lo(2),u_hi(2)
       do ic = u_lo(1),u_hi(1)

          !   S spatial indices start at lo - (number of ghost zones)
          ! uCF spatial indices start at 1 - (number of ghost zones)
          i = lo(1) + 2*(ic-1)
          j = lo(2) + 2*(jc-1)

          ! In 2-d, kc = 1
          kc = 1

          ! Thornado uses units where c = G = k = 1, Meter = 1

          uCF(1,ic,jc,kc,iCF_D)  = S(i,j,URHO)  * conv_dens
          uCF(1,ic,jc,kc,iCF_S1) = S(i,j,UMX)   * conv_mom
          uCF(1,ic,jc,kc,iCF_S2) = S(i,j,UMY)   * conv_mom
          uCF(1,ic,jc,kc,iCF_S3) = S(i,j,UMZ)   * conv_mom
          uCF(1,ic,jc,kc,iCF_E)  = S(i,j,UEDEN) * conv_enr
          uCF(1,ic,jc,kc,iCF_Ne) = S(i,j,UFX)   * conv_ne

          uCF(2,ic,jc,kc,iCF_D)  = S(i+1,j,URHO)  * conv_dens
          uCF(2,ic,jc,kc,iCF_S1) = S(i+1,j,UMX)   * conv_mom
          uCF(2,ic,jc,kc,iCF_S2) = S(i+1,j,UMY)   * conv_mom
          uCF(2,ic,jc,kc,iCF_S3) = S(i+1,j,UMZ)   * conv_mom
          uCF(2,ic,jc,kc,iCF_E)  = S(i+1,j,UEDEN) * conv_enr
          uCF(2,ic,jc,kc,iCF_Ne) = S(i+1,j,UFX)   * conv_ne

          uCF(3,ic,jc,kc,iCF_D)  = S(i,j+1,URHO)  * conv_dens
          uCF(3,ic,jc,kc,iCF_S1) = S(i,j+1,UMX)   * conv_mom
          uCF(3,ic,jc,kc,iCF_S2) = S(i,j+1,UMY)   * conv_mom
          uCF(3,ic,jc,kc,iCF_S3) = S(i,j+1,UMZ)   * conv_mom
          uCF(3,ic,jc,kc,iCF_E)  = S(i,j+1,UEDEN) * conv_enr
          uCF(3,ic,jc,kc,iCF_Ne) = S(i,j+1,UFX)   * conv_ne

          uCF(4,ic,jc,kc,iCF_D)  = S(i+1,j+1,URHO)  * conv_dens
          uCF(4,ic,jc,kc,iCF_S1) = S(i+1,j+1,UMX)   * conv_mom
          uCF(4,ic,jc,kc,iCF_S2) = S(i+1,j+1,UMY)   * conv_mom
          uCF(4,ic,jc,kc,iCF_S3) = S(i+1,j+1,UMZ)   * conv_mom
          uCF(4,ic,jc,kc,iCF_E)  = S(i+1,j+1,UEDEN) * conv_enr
          uCF(4,ic,jc,kc,iCF_Ne) = S(i+1,j+1,UFX)   * conv_ne

          ! Make this copy so we can create dS on nodes, not after S has been
          !      averaged back to cell centers
          u0(1:4,ic,jc,kc,iCF_D ) = uCF(1:4,ic,jc,kc,iCF_D )
          u0(1:4,ic,jc,kc,iCF_S1) = uCF(1:4,ic,jc,kc,iCF_S1)
          u0(1:4,ic,jc,kc,iCF_S2) = uCF(1:4,ic,jc,kc,iCF_S2)
          u0(1:4,ic,jc,kc,iCF_S3) = uCF(1:4,ic,jc,kc,iCF_S3)
          u0(1:4,ic,jc,kc,iCF_E ) = uCF(1:4,ic,jc,kc,iCF_E )
          u0(1:4,ic,jc,kc,iCF_Ne) = uCF(1:4,ic,jc,kc,iCF_Ne)

       end do
       end do
    ! ************************************************************************************
    ! Use reconstruction matrix
    ! ************************************************************************************
    elseif (interp_type .eq. 1) then
       do jc = u_lo(2),u_hi(2)
       do ic = u_lo(1),u_hi(1)

          !   S spatial indices start at lo - (number of ghost zones)
          ! uCF spatial indices start at 1 - (number of ghost zones)
          i = lo(1) + 2*(ic-1)
          j = lo(2) + 2*(jc-1)

          ! In 2-d, kc = 1
          kc = 1

          ! Thornado uses units where c = G = k = 1, Meter = 1

          do ind = 1, 4

             uCF(ind,ic,jc,kc,iCF_D) &
               = ( ReconstructionMatrix  (ind,1) * S(i,  j,  URHO) &
                   + ReconstructionMatrix(ind,2) * S(i+1,j,  URHO) &
                   + ReconstructionMatrix(ind,3) * S(i,  j+1,URHO) &
                   + ReconstructionMatrix(ind,4) * S(i+1,j+1,URHO) ) * conv_dens

             uCF(ind,ic,jc,kc,iCF_S1) &
               = ( ReconstructionMatrix  (ind,1) * S(i,  j,  UMX) &
                   + ReconstructionMatrix(ind,2) * S(i+1,j,  UMX) &
                   + ReconstructionMatrix(ind,3) * S(i,  j+1,UMX) &
                   + ReconstructionMatrix(ind,4) * S(i+1,j+1,UMX) ) * conv_mom

             uCF(ind,ic,jc,kc,iCF_S2) &
               = ( ReconstructionMatrix  (ind,1) * S(i,  j,  UMY) &
                   + ReconstructionMatrix(ind,2) * S(i+1,j,  UMY) &
                   + ReconstructionMatrix(ind,3) * S(i,  j+1,UMY) &
                   + ReconstructionMatrix(ind,4) * S(i+1,j+1,UMY) ) * conv_mom

             uCF(ind,ic,jc,kc,iCF_S3) &
               = ( ReconstructionMatrix  (ind,1) * S(i,  j,  UMZ) &
                   + ReconstructionMatrix(ind,2) * S(i+1,j,  UMZ) &
                   + ReconstructionMatrix(ind,3) * S(i,  j+1,UMZ) &
                   + ReconstructionMatrix(ind,4) * S(i+1,j+1,UMZ) ) * conv_mom

             uCF(ind,ic,jc,kc,iCF_E) &
               = ( ReconstructionMatrix  (ind,1) * S(i,  j,  UEDEN) &
                   + ReconstructionMatrix(ind,2) * S(i+1,j,  UEDEN) &
                   + ReconstructionMatrix(ind,3) * S(i,  j+1,UEDEN) &
                   + ReconstructionMatrix(ind,4) * S(i+1,j+1,UEDEN) ) * conv_enr

             uCF(ind,ic,jc,kc,iCF_Ne) &
               = ( ReconstructionMatrix  (ind,1) * S(i,  j,  UFX) &
                   + ReconstructionMatrix(ind,2) * S(i+1,j,  UFX) &
                   + ReconstructionMatrix(ind,3) * S(i,  j+1,UFX) &
                   + ReconstructionMatrix(ind,4) * S(i+1,j+1,UFX) ) * conv_ne

          end do

          ! Make this copy so we can create dS on nodes, not after S has been
          !      averaged back to cell centers
          u0(1:4,ic,jc,kc,iCF_D ) = uCF(1:4,ic,jc,kc,iCF_D )
          u0(1:4,ic,jc,kc,iCF_S1) = uCF(1:4,ic,jc,kc,iCF_S1)
          u0(1:4,ic,jc,kc,iCF_S2) = uCF(1:4,ic,jc,kc,iCF_S2)
          u0(1:4,ic,jc,kc,iCF_S3) = uCF(1:4,ic,jc,kc,iCF_S3)
          u0(1:4,ic,jc,kc,iCF_E ) = uCF(1:4,ic,jc,kc,iCF_E )
          u0(1:4,ic,jc,kc,iCF_Ne) = uCF(1:4,ic,jc,kc,iCF_Ne)

       end do
       end do
    end if

  end subroutine interpolate_fluid
