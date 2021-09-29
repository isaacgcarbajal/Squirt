module Globals
  use Parameters
  implicit none
  
  !--------------------------------------------------------------------------
  ! Number of equations in the system
  !--------------------------------------------------------------------------
  integer, parameter :: neq = 4
  
  !--------------------------------------------------------------------------
  ! Spatial global constants
  !--------------------------------------------------------------------------
	real*8, parameter :: dx = (xMax - xMin)/nx
	real*8, parameter :: dy = (yMax - yMin)/ny

  !--------------------------------------------------------------------------
  ! Variables and Fluxes
  !--------------------------------------------------------------------------
  real*8 :: prim(neq, nx+2, ny+2)     ! Primitive variables
  real*8 ::    u(neq, nx+2, ny+2)     ! Current conserved variables
  real*8 ::   up(neq, nx+2, ny+2)     ! Advanced conserved variables
  real*8 ::    f(neq, nx+2, ny+2)     ! Physical Fluxes in X
  real*8 ::    g(neq, nx+2, ny+2)     ! Physical Fluxes in Y

  !--------------------------------------------------------------------------
  ! Time and iteration variables
  !--------------------------------------------------------------------------
  integer :: currentIteration     ! Current number of iterations
  integer :: nextFileNumber       ! Number of the next output
  real*8 :: currentTime           ! Current time in the simulation
  real*8 :: timeToPrint           ! Time to next output
  real*8 :: dt                    ! Current valid time step

  !--------------------------------------------------------------------------
  ! Variables to meassure the excecution time
  !--------------------------------------------------------------------------
  integer :: clockStart, clockCount, clockRate, clockMax

end module Globals
