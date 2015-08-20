subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  use eos_module, only : gamma_const
  use bl_error_module
  implicit none

  integer :: init, namlen
  integer :: name(namlen)
  double precision :: problo(3), probhi(3)

  integer untin,i

  namelist /fortin/ frac, &
       rho_1, rho_2, p0_base

  ! Build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen = 256
  character probin*(maxlen)

  if (namlen .gt. maxlen) then
     call bl_error('probin file name too long')
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do
         
  ! set namelist defaults here
  frac = 0.5d0
  rho_1 = 1.0d0
  rho_2 = 2.0d0
  p0_base = 5.0d0

  ! Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)


  ! set local variable defaults
  split(1) = frac*(problo(1)+probhi(1))
  split(2) = frac*(problo(2)+probhi(2))
  split(3) = frac*(problo(3)+probhi(3))
  
  L_x = probhi(1) - problo(1)

end subroutine PROBINIT


! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize data on each grid.  
! ::: 
! ::: NOTE:  all arrays have one cell of ghost zones surrounding
! :::        the grid interior.  Values in these cells need not
! :::        be set here.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: level     => amr level of grid
! ::: time      => time at which to init data             
! ::: lo,hi     => index limits of grid interior (cell centered)
! ::: nstate    => number of state components.  You should know
! :::		   this already!
! ::: state     <=  Scalar array
! ::: delta     => cell size
! ::: xlo,xhi   => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------
subroutine ca_initdata(level,time,lo,hi,nscal, &
                       state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
                       delta,xlo,xhi)

  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, &
       UEDEN, UEINT, UFS, UTEMP
  use bl_constants_module, only: ZERO, HALF, M_PI
  use eos_module, only : gamma_const
  
  implicit none
        
  integer :: level, nscal
  integer :: lo(3), hi(3)
  integer :: state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
  double precision :: xlo(3), xhi(3), time, delta(3)
  double precision :: state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)
  
  integer :: i,j,k
  double precision :: x,y,z,r2d,pres,presmid,pertheight
  
  presmid  = p0_base - rho_1*split(3)
        
  state(:,:,:,UMX)   = ZERO
  state(:,:,:,UMY)   = ZERO
  state(:,:,:,UMZ)   = ZERO
  state(:,:,:,UTEMP) = ZERO

  do k = lo(3), hi(3)
     z = (k+HALF)*delta(3)

     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           
           if (z .lt. split(3)) then
              pres = p0_base - rho_1*z
              state(i,j,k,UEDEN) = pres / (gamma_const - 1.0d0)
              state(i,j,k,UEINT) = pres / (gamma_const - 1.0d0)
           else
              pres = presmid - rho_2*(z-split(3))
              state(i,j,k,UEDEN) = pres / (gamma_const - 1.0d0)
              state(i,j,k,UEINT) = pres / (gamma_const - 1.0d0)
           end if
           
        enddo
     enddo
  enddo

  do k = lo(3), hi(3)
     z = (k+HALF)*delta(3)
        
     do j = lo(2), hi(2)
        y = (j+HALF)*delta(2)
        
        do i = lo(1), hi(1)
           x = (i+HALF)*delta(1)
     
           r2d = min(sqrt((x-split(1))**2+(y-split(2))**2), 0.5d0*L_x)
           pertheight = 0.5d0 - 0.01d0*cos(2.0d0*M_PI*r2d/L_x)
           state(i,j,k,URHO) = rho_1 + ((rho_2-rho_1)/2.0d0)* &
                (1+tanh((z-pertheight)/0.005d0))
           state(i,j,k,UFS) = state(i,j,k,URHO)
           
        enddo
     enddo
  enddo

end subroutine ca_initdata

