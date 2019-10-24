module derive_thornado_module

  use amrex_fort_module, only : rt => amrex_real
  use RadiationFieldsModule, only : nSpecies
  use ProgramHeaderModule, only : nE, nNodesE

  implicit none

  public

contains
  
! All subroutines in this file must be threadsafe because they are called
! inside OpenMP parallel regions.

  subroutine ca_der_avgs_per_E(lo,hi,&
                               U_R_avg,j_lo,j_hi,nv, &
                               U_R    ,d_lo,d_hi,nc) &
    bind(C, name="ca_der_avgs_per_E")

    use ReferenceElementModule, only: weights_q, NodeNumberTable
    use GeometryFieldsModuleE , only: iGE_Ep2, uGE

    implicit none 

    integer, intent(in)     :: lo(3), hi(3)
    integer, intent(in)     :: j_lo(3), j_hi(3), nv
    integer, intent(in)     :: d_lo(3), d_hi(3), nc
    real(rt), intent(inout) :: U_R_avg(j_lo(1):j_hi(1),j_lo(2):j_hi(2),j_lo(3):j_hi(3),0:nv-1)
    real(rt), intent(in)    :: U_R    (d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),0:nc-1)

    integer          :: i, j, k, is, ie, id, im, ii, icomp, n_moments = 4
    integer          :: iNodeE
    integer          :: n_each
    real(rt)         :: weight_esquared, sum_esquared

    U_R_avg(:,:,:,:) = 0.0e0_rt  ! zero it out
 
    ! Divide by the number of moments
    n_each = nv / n_moments

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
          do is = 1, nSpecies

            do ie = 1, nE

              sum_esquared = 0.d0

              do id = 1, nNodesE ! radiation degrees of freedom

                iNodeE = NodeNumberTable(1,id)
  
                weight_esquared = uGE(iNodeE,ie,iGE_Ep2)
  
                sum_esquared = sum_esquared + weights_q(id) * weight_esquared
  
                im = 1  ! J is first moment
                ii = (is-1)*(n_moments*nE*nNodesE) + &
                     (im-1)*(nE*nNodesE) + &
                     (ie-1)*nNodesE + (id-1)
                icomp = (is-1)*nSpecies + ie-1
                U_R_avg(i,j,k,icomp) = U_R_avg(i,j,k,icomp) + U_R(i,j,k,ii) * weights_q(id) * weight_esquared
  
                im = 2  ! H_x is second moment
                ii = (is-1)*(n_moments*nE*nNodesE) + &
                     (im-1)*(nE*nNodesE) + &
                     (ie-1)*nNodesE + (id-1)
                icomp = (is-1)*nSpecies + n_each+ie-1
                U_R_avg(i,j,k,icomp) = U_R_avg(i,j,k,icomp) + U_R(i,j,k,ii) * weights_q(id) * weight_esquared
  
                im = 3  ! H_x is second moment
                ii = (is-1)*(n_moments*nE*nNodesE) + &
                     (im-1)*(nE*nNodesE) + &
                     (ie-1)*nNodesE + (id-1)
                icomp = (is-1)*nSpecies + 2*n_each+ie-1
                U_R_avg(i,j,k,icomp) = U_R_avg(i,j,k,icomp) + U_R(i,j,k,ii) * weights_q(id) * weight_esquared
  
                im = 4  ! H_z is second moment
                ii = (is-1)*(n_moments*nE*nNodesE) + &
                     (im-1)*(nE*nNodesE) + &
                     (ie-1)*nNodesE + (id-1)
                icomp = (is-1)*nSpecies + 3*n_each+ie-1
                U_R_avg(i,j,k,icomp) = U_R_avg(i,j,k,icomp) + U_R(i,j,k,ii) * weights_q(id) * weight_esquared
  
              end do

              icomp = (is-1)*nSpecies + ie-1
              U_R_avg(i,j,k,icomp) = U_R_avg(i,j,k,icomp) / sum_esquared

              icomp = (is-1)*nSpecies + n_each+ie-1
              U_R_avg(i,j,k,icomp) = U_R_avg(i,j,k,icomp) / sum_esquared

              icomp = (is-1)*nSpecies + 2*n_each+ie-1
              U_R_avg(i,j,k,icomp) = U_R_avg(i,j,k,icomp) / sum_esquared

              icomp = (is-1)*nSpecies + 3*n_each+ie-1
              U_R_avg(i,j,k,icomp) = U_R_avg(i,j,k,icomp) / sum_esquared

            end do

          end do
          end do
       end do
    end do

  end subroutine ca_der_avgs_per_E

  subroutine ca_der_J(J_avg,j_lo,j_hi,nv, &
                      U_R,d_lo,d_hi,nc,lo,hi,domlo, &
                      domhi,delta,xlo,time,dt,bc,level,grid_no) &
                      bind(C, name="ca_der_J")

    use ProgramHeaderModule   , only: nNodesE
    use ReferenceElementModule, only: weights_q, NodeNumberTable
    use GeometryFieldsModuleE , only: iGE_Ep2, uGE

    implicit none 

    integer, intent(in)     :: lo(3), hi(3)
    integer, intent(in)     :: j_lo(3), j_hi(3), nv
    integer, intent(in)     :: d_lo(3), d_hi(3), nc
    integer, intent(in)     :: domlo(3), domhi(3)
    integer, intent(in)     :: bc(3,2,nc)
    real(rt), intent(in)    :: delta(3), xlo(3), time, dt
    real(rt), intent(inout) :: J_avg(j_lo(1):j_hi(1),j_lo(2):j_hi(2),j_lo(3):j_hi(3),nv)
    real(rt), intent(in)    :: U_R(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),0:nc)
    integer, intent(in)     :: level, grid_no

    integer          :: i, j, k, is, ie, id, im, ii, icomp, n_moments = 4
    integer          :: iNodeE
    real(rt)         :: weight_esquared, sum_esquared

    J_avg(:,:,:,:) = 0.0e0_rt  ! zero it out

    im = 1  ! J is first moment

    if (nv .ne. 1) then
       print *,'NV NOT ONE IN CA_DER_J ', nv
       stop
    end if

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

            sum_esquared = 0.d0

            do is = 1, nSpecies
            do ie = 1, nE

              ii = (is-1)*(n_moments*nE*nNodesE) + &
                   (im-1)*(nE*nNodesE) + (ie-1)*nNodesE

              do id = 1, nNodesE

                 iNodeE = NodeNumberTable(1,id)

                 weight_esquared = uGE(iNodeE,ie,iGE_Ep2)

                 sum_esquared = sum_esquared + weights_q(id) * weight_esquared

                 icomp = ii + (id-1)

                 J_avg(i,j,k,1) = J_avg(i,j,k,1) + U_R(i,j,k,icomp) * weights_q(id) * weight_esquared

            end do
            end do
            end do

            J_avg(i,j,k,1) = J_avg(i,j,k,1) / sum_esquared

          end do
       end do
    end do

  end subroutine ca_der_J

  subroutine ca_der_Hx(Hx_avg,hx_lo,hx_hi,nv, &
                       U_R,d_lo,d_hi,nc,lo,hi,domlo, &
                       domhi,delta,xlo,time,dt,bc,level,grid_no) &
                       bind(C, name="ca_der_Hx")

    use ProgramHeaderModule   , only: nNodesE
    use ReferenceElementModule, only: weights_q, NodeNumberTable
    use GeometryFieldsModuleE , only: iGE_Ep2, uGE

    implicit none 

    integer, intent(in)     :: lo(3), hi(3)
    integer, intent(in)     :: hx_lo(3), hx_hi(3), nv
    integer, intent(in)     :: d_lo(3), d_hi(3), nc
    integer, intent(in)     :: domlo(3), domhi(3)
    integer, intent(in)     :: bc(3,2,nc)
    real(rt), intent(in)    :: delta(3), xlo(3), time, dt
    real(rt), intent(inout) :: Hx_avg(hx_lo(1):hx_hi(1),hx_lo(2):hx_hi(2),hx_lo(3):hx_hi(3),nv)
    real(rt), intent(in)    :: U_R(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),0:nc)
    integer, intent(in)     :: level, grid_no

    integer          :: i, j, k, is, ie, id, im, ii, icomp, n_moments = 4
    integer          :: iNodeE
    real(rt)         :: weight_esquared, sum_esquared

    Hx_avg(:,:,:,:) = 0.0e0_rt  ! zero it out

    im = 2  ! Hx is second moment

    if (nv .ne. 1) then
       print *,'NV NOT ONE IN CA_DER_HX ', nv
       stop
    end if

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

            sum_esquared = 0.d0

            do is = 1, nSpecies
            do ie = 1, nE

              ii = (is-1)*(n_moments*nE*nNodesE) + &
                   (im-1)*(nE*nNodesE) + (ie-1)*nNodesE

              do id = 1, nNodesE

                 iNodeE = NodeNumberTable(1,id)

                 weight_esquared = uGE(iNodeE,ie,iGE_Ep2)

                 sum_esquared = sum_esquared + weights_q(id) * weight_esquared

                 icomp = ii + (id-1)

                 Hx_avg(i,j,k,1) = Hx_avg(i,j,k,1) + U_R(i,j,k,icomp) * weights_q(id) * weight_esquared

            end do
            end do
            end do

            Hx_avg(i,j,k,1) = Hx_avg(i,j,k,1) / sum_esquared

          end do
       end do
    end do

  end subroutine ca_der_Hx


  subroutine ca_der_Hy(Hy_avg,hy_lo,hy_hi,nv, &
                      U_R,d_lo,d_hi,nc,lo,hi,domlo, &
                      domhi,delta,xlo,time,dt,bc,level,grid_no) &
                      bind(C, name="ca_der_Hy")

    use ProgramHeaderModule   , only: nNodesE
    use ReferenceElementModule, only: weights_q, NodeNumberTable
    use GeometryFieldsModuleE , only: iGE_Ep2, uGE

    implicit none 

    integer, intent(in)     :: lo(3), hi(3)
    integer, intent(in)     :: hy_lo(3), hy_hi(3), nv
    integer, intent(in)     :: d_lo(3), d_hi(3), nc
    integer, intent(in)     :: domlo(3), domhi(3)
    integer, intent(in)     :: bc(3,2,nc)
    real(rt), intent(in)    :: delta(3), xlo(3), time, dt
    real(rt), intent(inout) :: Hy_avg(hy_lo(1):hy_hi(1),hy_lo(2):hy_hi(2),hy_lo(3):hy_hi(3),nv)
    real(rt), intent(in)    :: U_R(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),0:nc)
    integer, intent(in)     :: level, grid_no

    integer          :: i, j, k, is, ie, id, im, ii, icomp, n_moments = 4
    integer          :: iNodeE
    real(rt)         :: weight_esquared, sum_esquared

    Hy_avg(:,:,:,:) = 0.0e0_rt  ! zero it out

    im = 3  ! Hy is third moment

    if (nv .ne. 1) then
       print *,'NV NOT ONE IN CA_DER_HY ', nv
       stop
    end if

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

            sum_esquared = 0.d0

            do is = 1, nSpecies
            do ie = 1, nE

              ii = (is-1)*(n_moments*nE*nNodesE) + &
                   (im-1)*(nE*nNodesE) + (ie-1)*nNodesE

              do id = 1, nNodesE

                 iNodeE = NodeNumberTable(1,id)

                 weight_esquared = uGE(iNodeE,ie,iGE_Ep2)

                 sum_esquared = sum_esquared + weights_q(id) * weight_esquared

                 icomp = ii + (id-1)

                 Hy_avg(i,j,k,1) = Hy_avg(i,j,k,1) + U_R(i,j,k,icomp) * weights_q(id) * weight_esquared

            end do
            end do
            end do

            Hy_avg(i,j,k,1) = Hy_avg(i,j,k,1) / sum_esquared

          end do
       end do
    end do

  end subroutine ca_der_Hy



  subroutine ca_der_Hz(Hz_avg,hz_lo,hz_hi,nv, &
                      U_R,d_lo,d_hi,nc,lo,hi,domlo, &
                      domhi,delta,xlo,time,dt,bc,level,grid_no) &
                      bind(C, name="ca_der_Hz")

    use ProgramHeaderModule   , only: nNodesE
    use ReferenceElementModule, only: weights_q, NodeNumberTable
    use GeometryFieldsModuleE , only: iGE_Ep2, uGE

    implicit none 

    integer, intent(in)     :: lo(3), hi(3)
    integer, intent(in)     :: hz_lo(3), hz_hi(3), nv
    integer, intent(in)     :: d_lo(3), d_hi(3), nc
    integer, intent(in)     :: domlo(3), domhi(3)
    integer, intent(in)     :: bc(3,2,nc)
    real(rt), intent(in)    :: delta(3), xlo(3), time, dt
    real(rt), intent(inout) :: Hz_avg(hz_lo(1):hz_hi(1),hz_lo(2):hz_hi(2),hz_lo(3):hz_hi(3),nv)
    real(rt), intent(in)    :: U_R(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),0:nc)
    integer, intent(in)     :: level, grid_no

    integer          :: i, j, k, is, ie, id, im, ii, icomp, n_moments = 4
    integer          :: iNodeE
    real(rt)         :: weight_esquared, sum_esquared

    Hz_avg(:,:,:,:) = 0.0e0_rt  ! zero it out

    im = 4  ! Hz is fourth moment

    if (nv .ne. 1) then
       print *,'NV NOT ONE IN CA_DER_HZ ', nv
       stop
    end if

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

            sum_esquared = 0.d0

            do is = 1, nSpecies
            do ie = 1, nE

              ii = (is-1)*(n_moments*nE*nNodesE) + &
                   (im-1)*(nE*nNodesE) + (ie-1)*nNodesE

              do id = 1, nNodesE

                 iNodeE = NodeNumberTable(1,id)

                 weight_esquared = uGE(iNodeE,ie,iGE_Ep2)

                 sum_esquared = sum_esquared + weights_q(id) * weight_esquared

                 icomp = ii + (id-1)

                 Hz_avg(i,j,k,1) = Hz_avg(i,j,k,1) + U_R(i,j,k,icomp) * weights_q(id) * weight_esquared

            end do
            end do
            end do

            Hz_avg(i,j,k,1) = Hz_avg(i,j,k,1) / sum_esquared

          end do
       end do
    end do

  end subroutine ca_der_Hz

end module derive_thornado_module
