module castro_util_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains



  function position(i, j, k, ccx, ccy, ccz)
    ! Given 3D indices (i,j,k), return the cell-centered spatial position.
    ! Optionally we can also be edge-centered in any of the directions.
    !

    use amrinfo_module, only: amr_level
    use prob_params_module, only: problo, probhi, physbc_lo, physbc_hi, dx_level, &
         domlo_level, domhi_level, Interior
    use amrex_constants_module, only: ZERO, HALF
    use amrex_fort_module, only: rt => amrex_real

    ! Input arguments

    integer :: i, j, k
    logical, optional :: ccx, ccy, ccz

    ! Local variables
    real(rt) :: position(3), dx(3), offset(3)
    integer  :: idx(3)
    logical  :: cc(3)
    integer  :: domlo(3), domhi(3)
    integer  :: dir

    !$gpu

    idx = [ i, j, k ]

    dx(:) = dx_level(:,amr_level)
    domlo = domlo_level(:,amr_level)
    domhi = domhi_level(:,amr_level)

    offset(:) = problo(:)

    cc(:) = .true.

    if (present(ccx)) then
       cc(1) = ccx
    endif

    if (present(ccy)) then
       cc(2) = ccy
    endif

    if (present(ccz)) then
       cc(3) = ccz
    endif

    do dir = 1, 3
       if (cc(dir)) then
          ! If we're cell-centered, we want to be in the middle of the zone.

          offset(dir) = offset(dir) + HALF * dx(dir)
       else
          ! Take care of the fact that for edge-centered indexing,
          ! we actually range from (domlo, domhi+1).

          domhi(dir) = domhi(dir) + 1
       endif
    enddo

    ! Be careful when using periodic boundary conditions. In that case,
    ! we need to loop around to the other side of the domain.

    do dir = 1, 3
       if      (physbc_lo(dir) .eq. Interior .and. idx(dir) .lt. domlo(dir)) then
          offset(dir) = offset(dir) + (probhi(dir) - problo(dir))
       else if (physbc_hi(dir) .eq. Interior .and. idx(dir) .gt. domhi(dir)) then
          offset(dir) = offset(dir) + (problo(dir) - probhi(dir))
       endif
    enddo

    position(:) = offset(:) + dble(idx(:)) * dx(:)

  end function position





  subroutine ca_enforce_consistent_e(lo,hi,state,s_lo,s_hi) bind(c,name='ca_enforce_consistent_e')
    ! Enforces (rho E) = (rho e) + 1/2 rho (u^2 + v^2 + w^2)

    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT
    use amrex_constants_module, only: HALF, ONE
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    ! Local variables
    integer  :: i,j,k
    real(rt) :: u, v, w, rhoInv

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rhoInv = ONE / state(i,j,k,URHO)
             u = state(i,j,k,UMX) * rhoInv
             v = state(i,j,k,UMY) * rhoInv
             w = state(i,j,k,UMZ) * rhoInv

             state(i,j,k,UEDEN) = state(i,j,k,UEINT) + &
                  HALF * state(i,j,k,URHO) * (u*u + v*v + w*w)

          end do
       end do
    end do

  end subroutine ca_enforce_consistent_e


  subroutine ca_recompute_energetics(lo, hi, state, s_lo, s_hi) bind(c,name='ca_recompute_energetics')
    ! Recomputes T and (rho e) from (rho E)

    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UTEMP, UFS, UFX, UEINT
    use amrex_constants_module, only: HALF, ONE
    use amrex_fort_module, only: rt => amrex_real
    use eos_type_module, only : eos_t, eos_input_re
    use eos_module, only : eos
    use network, only : nspec, naux
    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    ! Local variables
    integer  :: i,j,k
    real(rt) :: u, v, w, rhoInv

    type(eos_t) :: eos_state

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rhoInv = ONE / state(i,j,k,URHO)
             u = state(i,j,k,UMX) * rhoInv
             v = state(i,j,k,UMY) * rhoInv
             w = state(i,j,k,UMZ) * rhoInv

             eos_state % rho = state(i,j,k,URHO)
             eos_state % T = state(i,j,k,UTEMP)
             eos_state % e = state(i,j,k,UEINT) * rhoInv - HALF * (u*u + v*v + w*w)
             eos_state % xn(:) = state(i,j,k,UFS:UFS-1+nspec) * rhoInv
             eos_state % aux(:) = state(i,j,k,UFX:UFX+naux-1) * rhoInv

             call eos(eos_input_re, eos_state)

             state(i,j,k,UTEMP) = eos_state % T

             state(i,j,k,UEINT) = eos_state % rho * eos_state % e

          end do
       end do
    end do

  end subroutine ca_recompute_energetics


  subroutine ca_reset_internal_e(lo,hi,u,u_lo,u_hi,verbose) bind(c,name='ca_reset_internal_e')

    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_rt
    use network, only: nspec, naux
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UFX, &
         UTEMP, small_temp, allow_small_energy, &
         dual_energy_eta2
    use amrex_constants_module, only: ZERO, HALF, ONE
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in), value :: verbose
    integer, intent(in) :: u_lo(3), u_hi(3)
    real(rt), intent(inout) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)

    ! Local variables
    integer  :: i,j,k
    real(rt) :: Up, Vp, Wp, ke, rho_eint, eden, small_e, eint_new, rhoInv

    type (eos_t) :: eos_state

    !$gpu

    ! Reset internal energy

    ! First, check if the internal energy variable is
    ! smaller than the internal energy computed via
    ! a call to the EOS using the small temperature.
    ! If so, reset it using the current temperature,
    ! assuming it is at least as large as small_temp.

    if (allow_small_energy .eq. 0) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                rhoInv = ONE / u(i,j,k,URHO)
                Up = u(i,j,k,UMX) * rhoInv
                Vp = u(i,j,k,UMY) * rhoInv
                Wp = u(i,j,k,UMZ) * rhoInv
                ke = HALF * (Up**2 + Vp**2 + Wp**2)
                eden = u(i,j,k,UEDEN) * rhoInv

                eos_state % rho = u(i,j,k,URHO)
                eos_state % T   = small_temp
                eos_state % xn  = u(i,j,k,UFS:UFS+nspec-1) * rhoInv
                eos_state % aux = u(i,j,k,UFX:UFX+naux-1) * rhoInv

                call eos(eos_input_rt, eos_state)

                small_e = eos_state % e

                ! If E < small_e, reset it so that it's equal to internal + kinetic.

                if (eden < small_e) then

                   if (u(i,j,k,UEINT) * rhoInv < small_e) then

                      eos_state % T = max(u(i,j,k,UTEMP), small_temp)

                      call eos(eos_input_rt, eos_state)

                      u(i,j,k,UEINT) = u(i,j,k,URHO) * eos_state % e

                   endif

                   u(i,j,k,UEDEN) = u(i,j,k,UEINT) + u(i,j,k,URHO) * ke

                else

                   rho_eint = u(i,j,k,UEDEN) - u(i,j,k,URHO) * ke

                   ! Reset (e from e) if it's greater than eta * E.

                   if (rho_eint .gt. ZERO .and. rho_eint / u(i,j,k,UEDEN) .gt. dual_energy_eta2) then

                      u(i,j,k,UEINT) = rho_eint

                   endif

                   if (u(i,j,k,UEINT) * rhoInv < small_e) then

                      eos_state % T = max(u(i,j,k,UTEMP), small_temp)

                      call eos(eos_input_rt, eos_state)

                      u(i,j,k,UEINT) = u(i,j,k,URHO) * eos_state % e
                      u(i,j,k,UTEMP) = eos_state % T

                   endif

                endif

             enddo
          enddo
       enddo

    else

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                rhoInv = ONE/u(i,j,k,URHO)
                Up = u(i,j,k,UMX) * rhoInv
                Vp = u(i,j,k,UMY) * rhoInv
                Wp = u(i,j,k,UMZ) * rhoInv
                ke = HALF * (Up**2 + Vp**2 + Wp**2)

                if (u(i,j,k,UEDEN) < ZERO) then

                   if (u(i,j,k,UEINT) < ZERO) then

                      eos_state % rho   = u(i,j,k,URHO)
                      eos_state % T     = small_temp
                      eos_state % xn(:) = u(i,j,k,UFS:UFS+nspec-1) * rhoInv
                      eos_state % aux(1:naux) = u(i,j,k,UFX:UFX+naux-1) * rhoInv

                      call eos(eos_input_rt, eos_state)

                      u(i,j,k,UEINT) = u(i,j,k,URHO) * eos_state % e

                   endif

                   u(i,j,k,UEDEN) = u(i,j,k,UEINT) + u(i,j,k,URHO) * ke

                else

                   rho_eint = u(i,j,k,UEDEN) - u(i,j,k,URHO) * ke

                   ! Reset (e from e) if it's greater than eta * E.
                   if (rho_eint .gt. ZERO .and. rho_eint / u(i,j,k,UEDEN) .gt. dual_energy_eta2) then

                      u(i,j,k,UEINT) = rho_eint

                      ! If not resetting and little e is negative ...
                   else if (u(i,j,k,UEINT) .le. ZERO) then

                      eos_state % rho   = u(i,j,k,URHO)
                      eos_state % T     = small_temp
                      eos_state % xn(:) = u(i,j,k,UFS:UFS+nspec-1) * rhoInv
                      eos_state % aux(1:naux) = u(i,j,k,UFX:UFX+naux-1) * rhoInv

                      call eos(eos_input_rt, eos_state)

                      eint_new = eos_state % e

#ifndef AMREX_USE_CUDA
                      if (verbose .gt. 0) then
                         print *,'   '
                         print *,'>>> Warning: Castro_util.F90::reset_internal_energy  ',i,j,k
                         print *,'>>> ... resetting neg. e from EOS using small_temp'
                         print *,'>>> ... from ',u(i,j,k,UEINT)/u(i,j,k,URHO),' to ', eint_new
                         print *,'    '
                      end if
#endif

                      u(i,j,k,UEINT) = u(i,j,k,URHO) * eint_new
                      u(i,j,k,UTEMP) = eos_state % T

                   endif

                end if
             enddo
          enddo
       enddo

    endif

  end subroutine ca_reset_internal_e


  subroutine ca_compute_temp(lo,hi,state,s_lo,s_hi) bind(c,name='ca_compute_temp')

    use network, only: nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_input_re, eos_t
    use meth_params_module, only: NVAR, URHO, UEINT, UTEMP, &
         UFS, UFX
    use amrex_constants_module, only: ZERO, ONE
    use castro_error_module
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer , intent(in   ) :: lo(3),hi(3)
    integer , intent(in   ) :: s_lo(3),s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer  :: i,j,k
    real(rt) :: rhoInv

    type (eos_t) :: eos_state

    !$gpu

    ! First check the inputs for validity.

#ifndef AMREX_USE_CUDA
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             if (state(i,j,k,URHO) <= ZERO) then
                print *,'   '
                print *,'>>> Error: Castro_util.F90::ca_compute_temp ',i,j,k
                print *,'>>> ... negative density ',state(i,j,k,URHO)
                print *,'    '
                call castro_error("Error:: compute_temp_nd.f90")
             end if

             if (state(i,j,k,UEINT) <= ZERO) then
                print *,'   '
                print *,'>>> Warning: Castro_util.F90::ca_compute_temp ',i,j,k
                print *,'>>> ... negative (rho e) ',state(i,j,k,UEINT)
                print *,'   '
                call castro_error("Error:: compute_temp_nd.f90")
             end if

          enddo
       enddo
    enddo
#endif

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rhoInv = ONE / state(i,j,k,URHO)

             eos_state % rho = state(i,j,k,URHO)
             eos_state % T   = state(i,j,k,UTEMP) ! Initial guess for the EOS
             eos_state % e   = state(i,j,k,UEINT) * rhoInv
             eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) * rhoInv
             eos_state % aux = state(i,j,k,UFX:UFX+naux-1) * rhoInv

             call eos(eos_input_re, eos_state)

             state(i,j,k,UTEMP) = eos_state % T

          enddo
       enddo
    enddo

  end subroutine ca_compute_temp



  subroutine ca_check_initial_species(lo, hi, state, state_lo, state_hi) bind(c,name='ca_check_initial_species')

    use network           , only: nspec
    use meth_params_module, only: NVAR, URHO, UFS
#ifndef AMREX_USE_CUDA
    use castro_error_module, only: castro_error
#endif
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: state_lo(3), state_hi(3)
    real(rt), intent(in) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)

    ! Local variables
    integer  :: i, j, k
    real(rt) :: spec_sum

#ifndef NUCLEAR_EOS

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             spec_sum = sum(state(i,j,k,UFS:UFS+nspec-1))

#ifndef AMREX_USE_CUDA
             if (abs(state(i,j,k,URHO)-spec_sum) .gt. 1.e-8_rt * state(i,j,k,URHO)) then

                print *,'Sum of (rho X)_i vs rho at (i,j,k): ',i,j,k,spec_sum,state(i,j,k,URHO)
                call castro_error("Error:: Failed check of initial species summing to 1")

             end if
#endif

          enddo
       enddo
    enddo

#endif

  end subroutine ca_check_initial_species



  subroutine ca_normalize_species(lo, hi, u, u_lo, u_hi) bind(c,name='ca_normalize_species')

    use network, only: nspec
    use meth_params_module, only: NVAR, URHO, UFS
    use amrex_constants_module, only: ONE
    use extern_probin_module, only: small_x
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    real(rt), intent(inout) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)

    ! Local variables
    integer  :: i, j, k
    real(rt) :: xn(nspec)

    !$gpu

#ifndef NUCLEAR_EOS
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             xn = u(i,j,k,UFS:UFS+nspec-1)

             xn = max(small_x * u(i,j,k,URHO), min(u(i,j,k,URHO), xn))

             xn = u(i,j,k,URHO) * (xn / sum(xn))

             u(i,j,k,UFS:UFS+nspec-1) = xn

          enddo
       enddo
    enddo
#endif

  end subroutine ca_normalize_species





  function position_to_index(loc) result(index)
    ! Given 3D spatial coordinates, return the cell-centered zone indices closest to it.
    ! Optionally we can also be edge-centered in any of the directions.

    use amrinfo_module, only: amr_level
    use prob_params_module, only: dx_level, dim
    use amrex_fort_module, only: rt => amrex_real

    real(rt), intent(in) :: loc(3)

    integer :: index(3)

    index(1:dim)   = NINT(loc(1:dim) / dx_level(1:dim,amr_level))
    index(dim+1:3) = 0

  end function position_to_index





  function area(i, j, k, dir)
    ! Given 3D indices (i,j,k) and a direction dir, return the
    ! area of the face perpendicular to direction d. We assume
    ! the coordinates perpendicular to the dir axies are edge-centered.
    ! Note that Castro has no support for angular coordinates, so
    ! this function only provides Cartesian in 1D/2D/3D, Cylindrical (R-Z)
    ! in 2D, and Spherical in 1D.

    use amrinfo_module, only: amr_level
    use amrex_constants_module, only: ZERO, ONE, TWO, M_PI, FOUR
    use prob_params_module, only: dim, coord_type, dx_level
    use castro_error_module
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer, intent(in) :: i, j, k, dir

    real(rt) :: area

    logical :: cc(3) = .true.
    real(rt) :: dx(3), loc(3)

    ! Force edge-centering along the direction of interest

    cc(dir) = .false.

    dx = dx_level(:,amr_level)

    if (coord_type .eq. 0) then

       ! Cartesian (1D/2D/3D)

       if (dim .eq. 1) then

          select case (dir)

          case (1)
             area = ONE
          case default
             area = ZERO

          end select

       else if (dim .eq. 2) then

          select case (dir)

          case (1)
             area = dx(2)
          case (2)
             area = dx(1)
          case default
             area = ZERO

          end select

       else if (dim .eq. 3) then

          select case (dir)

          case (1)
             area = dx(2) * dx(3)
          case (2)
             area = dx(1) * dx(3)
          case (3)
             area = dx(1) * dx(2)
          case default
             area = ZERO

          end select

       endif

    else if (coord_type .eq. 1) then

       ! Cylindrical (2D only)

       ! Get edge-centered position

       loc = position(i,j,k,cc(1),cc(2),cc(3))

       if (dim .eq. 2) then

          select case (dir)

          case (1)
             area = TWO * M_PI * loc(1) * dx(2)
          case (2)
             area = TWO * M_PI * loc(1) * dx(1)
          case default
             area = ZERO

          end select

#ifndef AMREX_USE_CUDA
       else

          call castro_error("Cylindrical coordinates only supported in 2D.")
#endif

       endif

    else if (coord_type .eq. 2) then

       ! Spherical (1D only)

       ! Get edge-centered position

       loc = position(i,j,k,cc(1),cc(2),cc(3))

       if (dim .eq. 1) then

          select case (dir)

          case (1)
             area = FOUR * M_PI * loc(1)**2
          case default
             area = ZERO

          end select

#ifndef AMREX_USE_CUDA
       else

          call castro_error("Spherical coordinates only supported in 1D.")
#endif

       endif

    endif

  end function area





  function volume(i, j, k)
    ! Given 3D cell-centered indices (i,j,k), return the volume of the zone.
    ! Note that Castro has no support for angular coordinates, so
    ! this function only provides Cartesian in 1D/2D/3D, Cylindrical (R-Z)
    ! in 2D, and Spherical in 1D.

    use amrinfo_module, only: amr_level
    use castro_error_module
    use amrex_constants_module, only: ZERO, HALF, FOUR3RD, TWO, M_PI
    use prob_params_module, only: dim, coord_type, dx_level
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer, intent(in) :: i, j, k

    real(rt) :: volume

    real(rt) :: dx(3), loc_l(3), loc_r(3)

    dx = dx_level(:,amr_level)

    if (coord_type .eq. 0) then

       ! Cartesian (1D/2D/3D)

       if (dim .eq. 1) then

          volume = dx(1)

       else if (dim .eq. 2) then

          volume = dx(1) * dx(2)

       else if (dim .eq. 3) then

          volume = dx(1) * dx(2) * dx(3)

       endif

    else if (coord_type .eq. 1) then

       ! Cylindrical (2D only)

       ! Get inner and outer radii

       loc_l = position(i  ,j,k,ccx=.true.)
       loc_r = position(i+1,j,k,ccx=.true.)

       if (dim .eq. 2) then

          volume = TWO * M_PI * (HALF * (loc_l(1) + loc_r(1))) * dx(1) * dx(2)

#ifndef AMREX_USE_CUDA
       else

          call castro_error("Cylindrical coordinates only supported in 2D.")
#endif

       endif

    else if (coord_type .eq. 2) then

       ! Spherical (1D only)

       ! Get inner and outer radii

       loc_l = position(i  ,j,k,ccx=.true.)
       loc_r = position(i+1,j,k,ccx=.true.)

       if (dim .eq. 1) then

          volume = FOUR3RD * M_PI * (loc_r(1)**3 - loc_l(1)**3)

#ifndef AMREX_USE_CUDA
       else

          call castro_error("Spherical coordinates only supported in 1D.")
#endif

       endif

    endif

  end function volume




  logical function is_domain_corner(idx) result(is_corner)
    ! Given an index, determine whether it is on a domain corner or not.

    use prob_params_module, only: domlo_level, domhi_level
    use amrinfo_module, only: amr_level

    implicit none

    integer, intent(in) :: idx(3)

    integer :: domlo(3), domhi(3)

    is_corner = .false.

    domlo(:) = domlo_level(:, amr_level)
    domhi(:) = domhi_level(:, amr_level)

    if (idx(1) < domlo(1) .or. idx(1) > domhi(1)) then
       if (idx(2) < domlo(2) .or. idx(2) > domhi(2)) then
          if (idx(3) < domlo(3) .or. idx(3) > domhi(3)) then
             is_corner = .true.
          end if
       end if
    end if

  end function is_domain_corner





  subroutine ca_get_center(center_out) bind(C, name="ca_get_center")
    ! Get the current center of the problem.  This may not be the
    ! center of the domain, due to any problem symmetries.

    use prob_params_module, only: center
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(inout) :: center_out(3)

    center_out = center

  end subroutine ca_get_center





  subroutine ca_set_center(center_in) bind(C, name="ca_set_center")
    !
    ! .. note::
    !    Binds to C function ``ca_set_center``

    use prob_params_module, only: center
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: center_in(3)

    center = center_in

  end subroutine ca_set_center





  subroutine ca_find_center(data,new_center,icen,dx,problo) &
       bind(C, name="ca_find_center")
    !
    ! .. note::
    !    Binds to C function ``ca_find_center``

    use amrex_constants_module, only: ZERO, HALF, TWO
    use prob_params_module, only: dg, dim
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(inout) :: data(-1:1,-1*dg(2):1*dg(2),-1*dg(3):1*dg(3))
    real(rt), intent(  out) :: new_center(3)
    real(rt), intent(in   ) :: dx(3),problo(3)

    real(rt) :: a,b,x,y,z,cen
    integer  :: icen(3)
    integer  :: i,j,k

    if (dim .eq. 1) then

       ! In 1-D it only make sense to have the center at the origin
       new_center = ZERO

    else if (dim .ge. 2) then

       ! We do this to take care of precision issues
       cen = data(0,0,0)
       do k = -1*dg(3),1*dg(3)
          do j = -1*dg(2),1*dg(2)
             do i = -1*dg(1),1*dg(1)
                data(i,j,k) = data(i,j,k) - cen
             end do
          end do
       end do

       ! This puts the "center" at the cell center
       new_center(1:dim) = problo(1:dim) +  (icen(1:dim)+HALF) * dx(1:dim)

       ! Fit parabola y = a x^2  + b x + c through three points
       ! a = 1/2 ( y_1 + y_-1)
       ! b = 1/2 ( y_1 - y_-1)
       ! x_vertex = -b / 2a

       ! ... in x-direction
       a = HALF * (data(1,0,0) + data(-1,0,0)) - data(0,0,0)
       b = HALF * (data(1,0,0) - data(-1,0,0)) - data(0,0,0)
       x = -b / (TWO*a)
       new_center(1) = new_center(1) +  x*dx(1)

       ! ... in y-direction
       a = HALF * (data(0,1,0) + data(0,-1,0)) - data(0,0,0)
       b = HALF * (data(0,1,0) - data(0,-1,0)) - data(0,0,0)
       y = -b / (TWO*a)
       new_center(2) = new_center(2) +  y*dx(2)

       if (dim .eq. 3) then

          ! ... in z-direction
          a = HALF * (data(0,0,1) + data(0,0,-1)) - data(0,0,0)
          b = HALF * (data(0,0,1) - data(0,0,-1)) - data(0,0,0)
          z = -b / (TWO*a)
          new_center(3) = new_center(3) +  z*dx(3)

       endif

    endif

  end subroutine ca_find_center





  subroutine ca_compute_avgstate(lo,hi,dx,dr,nc,&
       state,s_lo,s_hi,radial_state, &
       vol,v_lo,v_hi,radial_vol, &
       problo,numpts_1d) &
       bind(C, name="ca_compute_avgstate")
    !
    ! .. note::
    !    Binds to C function ``ca_compute_avgstate``

    use meth_params_module, only: URHO, UMX, UMY, UMZ
    use prob_params_module, only: center, dim
    use amrex_constants_module, only: HALF
    use castro_error_module
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3),hi(3),nc
    real(rt), intent(in   ) :: dx(3),dr,problo(3)

    integer,  intent(in   ) :: numpts_1d
    real(rt), intent(inout) :: radial_state(nc,0:numpts_1d-1)
    real(rt), intent(inout) :: radial_vol(0:numpts_1d-1)

    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc)

    integer,  intent(in   ) :: v_lo(3), v_hi(3)
    real(rt), intent(in   ) :: vol(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))

    integer  :: i,j,k,n,index
    real(rt) :: x,y,z,r
    real(rt) :: x_mom,y_mom,z_mom,radial_mom

#ifndef AMREX_USE_CUDA
    if (dim .eq. 1) call castro_error("Error: cannot do ca_compute_avgstate in 1D.")
#endif

    !
    ! Do not OMP this.
    !
    do k = lo(3), hi(3)
       z = problo(3) + (dble(k)+HALF) * dx(3) - center(3)
       do j = lo(2), hi(2)
          y = problo(2) + (dble(j)+HALF) * dx(2) - center(2)
          do i = lo(1), hi(1)
             x = problo(1) + (dble(i)+HALF) * dx(1) - center(1)
             r = sqrt(x**2 + y**2 + z**2)
             index = int(r/dr)
#ifndef AMREX_USE_CUDA
             if (index .gt. numpts_1d-1) then
                print *,'COMPUTE_AVGSTATE: INDEX TOO BIG ',index,' > ',numpts_1d-1
                print *,'AT (i,j,k) ',i,j,k
                print *,'R / DR ',r,dr
                call castro_error("Error:: Castro_util.F90 :: ca_compute_avgstate")
             end if
#endif
             radial_state(URHO,index) = radial_state(URHO,index) &
                  + vol(i,j,k)*state(i,j,k,URHO)
             !
             ! Store the radial component of the momentum in the
             ! UMX, UMY and UMZ components for now.
             !
             x_mom = state(i,j,k,UMX)
             y_mom = state(i,j,k,UMY)
             z_mom = state(i,j,k,UMZ)
             radial_mom = x_mom * (x/r) + y_mom * (y/r) + z_mom * (z/r)
             radial_state(UMX,index) = radial_state(UMX,index) + vol(i,j,k)*radial_mom
             radial_state(UMY,index) = radial_state(UMY,index) + vol(i,j,k)*radial_mom
             radial_state(UMZ,index) = radial_state(UMZ,index) + vol(i,j,k)*radial_mom

             do n = UMZ+1,nc
                radial_state(n,index) = radial_state(n,index) + vol(i,j,k)*state(i,j,k,n)
             end do
             radial_vol(index) = radial_vol(index) + vol(i,j,k)
          enddo
       enddo
    enddo

  end subroutine ca_compute_avgstate




  function linear_to_angular_momentum(loc, mom) result(ang_mom)

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: loc(3), mom(3)

    real(rt) :: ang_mom(3)

    !$gpu

    ang_mom(1) = loc(2) * mom(3) - loc(3) * mom(2)
    ang_mom(2) = loc(3) * mom(1) - loc(1) * mom(3)
    ang_mom(3) = loc(1) * mom(2) - loc(2) * mom(1)

  end function linear_to_angular_momentum

end module castro_util_module
