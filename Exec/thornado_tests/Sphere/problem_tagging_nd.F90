module problem_tagging_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  public

contains

  subroutine set_problem_tags(lo, hi, &
                              tag, tag_lo, tag_hi, &
                              state, state_lo, state_hi, &
                              dx, problo, &
                              set, clear, time, level) &
                              bind(C, name="set_problem_tags")

    use amrex_constants_module, only: HALF
    use meth_params_module, only: URHO, NVAR, UFX
    use probdata_module, only: ye_err, ye_grad, ye_grad_rel, max_ye_err_lev, max_ye_grad_lev, max_ye_grad_rel_lev
    use prob_params_module, only: dg

    implicit none

    integer,    intent(in   ) :: lo(3), hi(3)
    integer,    intent(in   ) :: tag_lo(3), tag_hi(3)
    integer,    intent(in   ) :: state_lo(3), state_hi(3)
    integer(1), intent(inout) :: tag(tag_lo(1):tag_hi(1),tag_lo(2):tag_hi(2),tag_lo(3):tag_hi(3))
    real(rt),   intent(in   ) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)
    real(rt),   intent(in   ) :: problo(3), dx(3)
    integer(1), intent(in   ), value :: set, clear
    integer,    intent(in   ), value :: level
    real(rt),   intent(in   ), value :: time

    integer  :: i, j, k
    real(rt) :: ax, ay, az

    !$gpu

    !     Tag on regions of high Ye
    if (level .lt. max_ye_err_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (state(i,j,k,UFX) .ge. ye_err) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

    !     Tag on regions of high Ye gradient
    if (level .lt. max_ye_grad_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ax = ABS(state(i+1*dg(1),j,k,UFX) - state(i,j,k,UFX))
                ay = ABS(state(i,j+1*dg(2),k,UFX) - state(i,j,k,UFX))
                az = ABS(state(i,j,k+1*dg(3),UFX) - state(i,j,k,UFX))
                ax = MAX(ax,ABS(state(i,j,k,UFX) - state(i-1*dg(1),j,k,UFX)))
                ay = MAX(ay,ABS(state(i,j,k,UFX) - state(i,j-1*dg(2),k,UFX)))
                az = MAX(az,ABS(state(i,j,k,UFX) - state(i,j,k-1*dg(3),UFX)))

                if (MAX(ax,ay,az) .ge. ye_grad) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

    !     Tag on regions of high Ye gradient
    if (level .lt. max_ye_grad_rel_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ax = ABS(state(i+1*dg(1),j,k,UFX) - state(i,j,k,UFX))
                ay = ABS(state(i,j+1*dg(2),k,UFX) - state(i,j,k,UFX))
                az = ABS(state(i,j,k+1*dg(3),UFX) - state(i,j,k,UFX))
                ax = MAX(ax,ABS(state(i,j,k,UFX) - state(i-1*dg(1),j,k,UFX)))
                ay = MAX(ay,ABS(state(i,j,k,UFX) - state(i,j-1*dg(2),k,UFX)))
                az = MAX(az,ABS(state(i,j,k,UFX) - state(i,j,k-1*dg(3),UFX)))

                if (MAX(ax,ay,az) .ge. ABS(ye_grad_rel * state(i,j,k,UFX))) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine set_problem_tags

end module problem_tagging_module
