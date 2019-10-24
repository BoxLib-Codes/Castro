
  subroutine interpolate_fluid (lo, hi, &
                                S , s_lo, s_hi, ns , &
                                u0, u_lo, u_hi, nu, &
                                ng) 

    use amrex_constants_module, only : fourth, half, zero, one, two
    use amrex_fort_module, only : rt => amrex_real
    use amrex_error_module, only : amrex_abort
    use meth_params_module, only : URHO,UMX,UMY,UMZ,UEDEN,UFX,UTEMP
    use SubcellReconstructionModule, only : ReconstructionMatrix
    use FluidFieldsModule, only : uCF, nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne
    use EquationOfStateModule_TABLE, only : ComputeTemperatureFromSpecificInternalEnergy_TABLE, &
                                            ComputeThermodynamicStates_Primitive_TABLE
    use UnitsModule, only : Gram, Centimeter, Second, AtomicMassUnit, Erg, Kelvin

    use ReferenceElementModuleX, only: NodesX_q

    implicit none
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: s_lo(3),  s_hi(3)
    integer, intent(in) :: u_lo(3),  u_hi(3)
    integer, intent(in) :: ns, nu
    integer, intent(in) :: ng

    ! Conserved fluid state (rho, rho u, rho v, rho w, rho E, rho e...)
    real(rt), intent(inout) ::  S(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),ns) 

    real(rt), intent(inout) :: u0(nu,u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nCF) 

    ! Temporary variables
    integer  :: i,j,k,n
    integer  :: ic,jc,kc
    integer  :: ind
    real(rt) :: conv_dens, conv_mom, conv_enr, conv_ne

    ! For interpolation
    real(rt) :: xslope(ns), yslope(ns), zslope(ns)
    real(rt) :: xyslope(ns), xzslope(ns), yzslope(ns)
    real(rt) :: Sval(ns)
    real(rt) :: dlft(ns), drgt(ns), dcen(ns), dlim, dsgn

    real(rt) :: rho_in(1)
    real(rt) :: Ye_in(1)
    real(rt) :: E_in(1)
    real(rt) :: T_in(1)

    real(rt) :: T_out(1)
    real(rt) :: Ne_out(1)
    real(rt) :: Epervol_out(1)
    real(rt) :: Epermass_out(1)

    real(rt) :: fac(ns),fac_all
    real(rt) :: delta_S_slope(ns)
    real(rt) :: delta_S_val(ns)
    real(rt) :: S_lll(ns),S_llh(ns),S_lhl(ns),S_lhh(ns)
    real(rt) :: S_hll(ns),S_hlh(ns),S_hhl(ns),S_hhh(ns)

    real(rt) :: x(8)
    real(rt) :: y(8)
    real(rt) :: z(8)

    integer  :: interp_type

    ! No interpolation
    !interp_type = 0
    ! Use reconstruction matrix:
    interp_type = 1

    ! Conversion to thornado units
    conv_dens = Gram / Centimeter**3
    conv_mom  = Gram / Centimeter**2 / Second
    conv_enr  = Erg / Centimeter**3
    conv_ne   = 1.d0 / Centimeter**3

    ! These are the locations of the DG nodes in the space [-.5:.5]
    ! Note that in 3D with 2 nodes in each direction, they are ordered 
    ! ind = 1:  (Lo,Lo,Lo)
    ! ind = 2:  (Hi,Lo,Lo)
    ! ind = 3:  (Lo,Hi,Lo)
    ! ind = 4:  (Hi,Hi,Lo)
    ! ind = 5:  (Lo,Lo,Hi)
    ! ind = 6:  (Hi,Lo,Hi)
    ! ind = 7:  (Lo,Hi,Hi)
    ! ind = 8:  (Hi,Hi,Hi)
     do ind = 1, nu
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

       do kc = u_lo(3),u_hi(3)
       do jc = u_lo(2),u_hi(2)
       do ic = u_lo(1),u_hi(1)

          !   S spatial indices start at lo - (number of ghost zones)
          ! uCF spatial indices start at 1 - (number of ghost zones)
          i = lo(1) + 2*(ic-1)
          j = lo(2) + 2*(jc-1)
          k = lo(3) + 2*(kc-1)

          ! Thornado uses units where c = G = k = 1, Meter = 1
          uCF(1,ic,jc,kc,iCF_D)  = S(i,j,k,URHO)  * conv_dens
          uCF(1,ic,jc,kc,iCF_S1) = S(i,j,k,UMX)   * conv_mom
          uCF(1,ic,jc,kc,iCF_S2) = S(i,j,k,UMY)   * conv_mom
          uCF(1,ic,jc,kc,iCF_S3) = S(i,j,k,UMZ)   * conv_mom
          uCF(1,ic,jc,kc,iCF_E)  = S(i,j,k,UEDEN) * conv_enr
          uCF(1,ic,jc,kc,iCF_Ne) = S(i,j,k,UFX)   * conv_ne

          uCF(2,ic,jc,kc,iCF_D)  = S(i+1,j,k,URHO)  * conv_dens
          uCF(2,ic,jc,kc,iCF_S1) = S(i+1,j,k,UMX)   * conv_mom
          uCF(2,ic,jc,kc,iCF_S2) = S(i+1,j,k,UMY)   * conv_mom
          uCF(2,ic,jc,kc,iCF_S3) = S(i+1,j,k,UMZ)   * conv_mom
          uCF(2,ic,jc,kc,iCF_E)  = S(i+1,j,k,UEDEN) * conv_enr
          uCF(2,ic,jc,kc,iCF_Ne) = S(i+1,j,k,UFX)   * conv_ne

          uCF(3,ic,jc,kc,iCF_D)  = S(i,j+1,k,URHO)  * conv_dens
          uCF(3,ic,jc,kc,iCF_S1) = S(i,j+1,k,UMX)   * conv_mom
          uCF(3,ic,jc,kc,iCF_S2) = S(i,j+1,k,UMY)   * conv_mom
          uCF(3,ic,jc,kc,iCF_S3) = S(i,j+1,k,UMZ)   * conv_mom
          uCF(3,ic,jc,kc,iCF_E)  = S(i,j+1,k,UEDEN) * conv_enr
          uCF(3,ic,jc,kc,iCF_Ne) = S(i,j+1,k,UFX)   * conv_ne

          uCF(4,ic,jc,kc,iCF_D)  = S(i+1,j+1,k,URHO)  * conv_dens
          uCF(4,ic,jc,kc,iCF_S1) = S(i+1,j+1,k,UMX)   * conv_mom
          uCF(4,ic,jc,kc,iCF_S2) = S(i+1,j+1,k,UMY)   * conv_mom
          uCF(4,ic,jc,kc,iCF_S3) = S(i+1,j+1,k,UMZ)   * conv_mom
          uCF(4,ic,jc,kc,iCF_E)  = S(i+1,j+1,k,UEDEN) * conv_enr
          uCF(4,ic,jc,kc,iCF_Ne) = S(i+1,j+1,k,UFX)   * conv_ne

          uCF(5,ic,jc,kc,iCF_D)  = S(i,j,k+1,URHO)  * conv_dens
          uCF(5,ic,jc,kc,iCF_S1) = S(i,j,k+1,UMX)   * conv_mom
          uCF(5,ic,jc,kc,iCF_S2) = S(i,j,k+1,UMY)   * conv_mom
          uCF(5,ic,jc,kc,iCF_S3) = S(i,j,k+1,UMZ)   * conv_mom
          uCF(5,ic,jc,kc,iCF_E)  = S(i,j,k+1,UEDEN) * conv_enr
          uCF(5,ic,jc,kc,iCF_Ne) = S(i,j,k+1,UFX)   * conv_ne

          uCF(6,ic,jc,kc,iCF_D)  = S(i+1,j,k+1,URHO)  * conv_dens
          uCF(6,ic,jc,kc,iCF_S1) = S(i+1,j,k+1,UMX)   * conv_mom
          uCF(6,ic,jc,kc,iCF_S2) = S(i+1,j,k+1,UMY)   * conv_mom
          uCF(6,ic,jc,kc,iCF_S3) = S(i+1,j,k+1,UMZ)   * conv_mom
          uCF(6,ic,jc,kc,iCF_E)  = S(i+1,j,k+1,UEDEN) * conv_enr
          uCF(6,ic,jc,kc,iCF_Ne) = S(i+1,j,k+1,UFX)   * conv_ne

          uCF(7,ic,jc,kc,iCF_D)  = S(i,j+1,k+1,URHO)  * conv_dens
          uCF(7,ic,jc,kc,iCF_S1) = S(i,j+1,k+1,UMX)   * conv_mom
          uCF(7,ic,jc,kc,iCF_S2) = S(i,j+1,k+1,UMY)   * conv_mom
          uCF(7,ic,jc,kc,iCF_S3) = S(i,j+1,k+1,UMZ)   * conv_mom
          uCF(7,ic,jc,kc,iCF_E)  = S(i,j+1,k+1,UEDEN) * conv_enr
          uCF(7,ic,jc,kc,iCF_Ne) = S(i,j+1,k+1,UFX)   * conv_ne

          uCF(8,ic,jc,kc,iCF_D)  = S(i+1,j+1,k+1,URHO)  * conv_dens
          uCF(8,ic,jc,kc,iCF_S1) = S(i+1,j+1,k+1,UMX)   * conv_mom
          uCF(8,ic,jc,kc,iCF_S2) = S(i+1,j+1,k+1,UMY)   * conv_mom
          uCF(8,ic,jc,kc,iCF_S3) = S(i+1,j+1,k+1,UMZ)   * conv_mom
          uCF(8,ic,jc,kc,iCF_E)  = S(i+1,j+1,k+1,UEDEN) * conv_enr
          uCF(8,ic,jc,kc,iCF_Ne) = S(i+1,j+1,k+1,UFX)   * conv_ne

          ! Make this copy so we can create dS on nodes, not after S has been
          !      averaged back to cell centers
          u0(1:8,ic,jc,kc,iCF_D ) = uCF(1:8,ic,jc,kc,iCF_D )
          u0(1:8,ic,jc,kc,iCF_S1) = uCF(1:8,ic,jc,kc,iCF_S1)
          u0(1:8,ic,jc,kc,iCF_S2) = uCF(1:8,ic,jc,kc,iCF_S2)
          u0(1:8,ic,jc,kc,iCF_S3) = uCF(1:8,ic,jc,kc,iCF_S3)
          u0(1:8,ic,jc,kc,iCF_E ) = uCF(1:8,ic,jc,kc,iCF_E )
          u0(1:8,ic,jc,kc,iCF_Ne) = uCF(1:8,ic,jc,kc,iCF_Ne)

       end do
       end do
       end do
    ! ************************************************************************************
    ! Use reconstruction matrix
    ! ************************************************************************************
    elseif (interp_type .eq. 1) then

       do kc = u_lo(3),u_hi(3)
       do jc = u_lo(2),u_hi(2)
       do ic = u_lo(1),u_hi(1)

          !   S spatial indices start at lo - (number of ghost zones)
          ! uCF spatial indices start at 1  - (number of ghost zones)
          i = lo(1) + 2*(ic-1)
          j = lo(2) + 2*(jc-1)
          k = lo(3) + 2*(kc-1)

          ! Thornado uses units where c = G = k = 1, Meter = 1

          do ind = 1, 8

            uCF(ind,ic,jc,kc,iCF_D) &
              = ( ReconstructionMatrix  (ind,1) * S(i,  j,  k,  URHO) &
                  + ReconstructionMatrix(ind,2) * S(i+1,j,  k,  URHO) &
                  + ReconstructionMatrix(ind,3) * S(i,  j+1,k,  URHO) &
                  + ReconstructionMatrix(ind,4) * S(i+1,j+1,k,  URHO) &
                  + ReconstructionMatrix(ind,5) * S(i,  j,  k+1,URHO) &
                  + ReconstructionMatrix(ind,6) * S(i+1,j,  k+1,URHO) &
                  + ReconstructionMatrix(ind,7) * S(i,  j+1,k+1,URHO) &
                  + ReconstructionMatrix(ind,8) * S(i+1,j+1,k+1,URHO) ) * conv_dens

            uCF(ind,ic,jc,kc,iCF_S1) &
              = ( ReconstructionMatrix  (ind,1) * S(i,  j,  k,  UMX) &
                  + ReconstructionMatrix(ind,2) * S(i+1,j,  k,  UMX) &
                  + ReconstructionMatrix(ind,3) * S(i,  j+1,k,  UMX) &
                  + ReconstructionMatrix(ind,4) * S(i+1,j+1,k,  UMX) &
                  + ReconstructionMatrix(ind,5) * S(i,  j,  k+1,UMX) &
                  + ReconstructionMatrix(ind,6) * S(i+1,j,  k+1,UMX) &
                  + ReconstructionMatrix(ind,7) * S(i,  j+1,k+1,UMX) &
                  + ReconstructionMatrix(ind,8) * S(i+1,j+1,k+1,UMX) ) * conv_mom

            uCF(ind,ic,jc,kc,iCF_S2) &
              = ( ReconstructionMatrix  (ind,1) * S(i,  j,  k,  UMY) &
                  + ReconstructionMatrix(ind,2) * S(i+1,j,  k,  UMY) &
                  + ReconstructionMatrix(ind,3) * S(i,  j+1,k,  UMY) &
                  + ReconstructionMatrix(ind,4) * S(i+1,j+1,k,  UMY) &
                  + ReconstructionMatrix(ind,5) * S(i,  j,  k+1,UMY) &
                  + ReconstructionMatrix(ind,6) * S(i+1,j,  k+1,UMY) &
                  + ReconstructionMatrix(ind,7) * S(i,  j+1,k+1,UMY) &
                  + ReconstructionMatrix(ind,8) * S(i+1,j+1,k+1,UMY) ) * conv_mom

            uCF(ind,ic,jc,kc,iCF_S3) &
              = ( ReconstructionMatrix  (ind,1) * S(i,  j,  k,  UMZ) &
                  + ReconstructionMatrix(ind,2) * S(i+1,j,  k,  UMZ) &
                  + ReconstructionMatrix(ind,3) * S(i,  j+1,k,  UMZ) &
                  + ReconstructionMatrix(ind,4) * S(i+1,j+1,k,  UMZ) &
                  + ReconstructionMatrix(ind,5) * S(i,  j,  k+1,UMZ) &
                  + ReconstructionMatrix(ind,6) * S(i+1,j,  k+1,UMZ) &
                  + ReconstructionMatrix(ind,7) * S(i,  j+1,k+1,UMZ) &
                  + ReconstructionMatrix(ind,8) * S(i+1,j+1,k+1,UMZ) ) * conv_mom

            uCF(ind,ic,jc,kc,iCF_E) &
              = ( ReconstructionMatrix  (ind,1) * S(i,  j,  k,  UEDEN) &
                  + ReconstructionMatrix(ind,2) * S(i+1,j,  k,  UEDEN) &
                  + ReconstructionMatrix(ind,3) * S(i,  j+1,k,  UEDEN) &
                  + ReconstructionMatrix(ind,4) * S(i+1,j+1,k,  UEDEN) &
                  + ReconstructionMatrix(ind,5) * S(i,  j,  k+1,UEDEN) &
                  + ReconstructionMatrix(ind,6) * S(i+1,j,  k+1,UEDEN) &
                  + ReconstructionMatrix(ind,7) * S(i,  j+1,k+1,UEDEN) &
                  + ReconstructionMatrix(ind,8) * S(i+1,j+1,k+1,UEDEN) ) * conv_enr

            uCF(ind,ic,jc,kc,iCF_Ne) &
              = ( ReconstructionMatrix  (ind,1) * S(i,  j,  k,  UFX) &
                  + ReconstructionMatrix(ind,2) * S(i+1,j,  k,  UFX) &
                  + ReconstructionMatrix(ind,3) * S(i,  j+1,k,  UFX) &
                  + ReconstructionMatrix(ind,4) * S(i+1,j+1,k,  UFX) &
                  + ReconstructionMatrix(ind,5) * S(i,  j,  k+1,UFX) &
                  + ReconstructionMatrix(ind,6) * S(i+1,j,  k+1,UFX) &
                  + ReconstructionMatrix(ind,7) * S(i,  j+1,k+1,UFX) &
                  + ReconstructionMatrix(ind,8) * S(i+1,j+1,k+1,UFX) ) * conv_ne

          end do

          ! Make this copy so we can create dS on nodes, not after S has been
          !      averaged back to cell centers
          u0(1:8,ic,jc,kc,iCF_D ) = uCF(1:8,ic,jc,kc,iCF_D )
          u0(1:8,ic,jc,kc,iCF_S1) = uCF(1:8,ic,jc,kc,iCF_S1)
          u0(1:8,ic,jc,kc,iCF_S2) = uCF(1:8,ic,jc,kc,iCF_S2)
          u0(1:8,ic,jc,kc,iCF_S3) = uCF(1:8,ic,jc,kc,iCF_S3)
          u0(1:8,ic,jc,kc,iCF_E ) = uCF(1:8,ic,jc,kc,iCF_E )
          u0(1:8,ic,jc,kc,iCF_Ne) = uCF(1:8,ic,jc,kc,iCF_Ne)

       end do
       end do
       end do

    end if

  end subroutine interpolate_fluid
