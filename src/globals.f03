module Globals
  use Parameters
  implicit none
  
  !--------------------------------------------------------------------------
  ! Presicion of the real variables
  !--------------------------------------------------------------------------
  integer, parameter ::                             &
      sp = kind(1.0),                               &
      dp = selected_real_kind(2*precision(1.0_sp)), &
      qp = selected_real_kind(2*precision(1.0_dp))
  
  ! real precision (rp)
  integer, parameter :: rp = dp
  
  !--------------------------------------------------------------------------
  ! Number of equations in the system
  !--------------------------------------------------------------------------
  integer, parameter :: neq = 4
  
  !--------------------------------------------------------------------------
  ! Spatial global constants
  !--------------------------------------------------------------------------
  real(rp), parameter :: dx = (xMax - xMin)/nx
  real(rp), parameter :: dy = (yMax - yMin)/ny

  !--------------------------------------------------------------------------
  ! Variables and Fluxes
  !--------------------------------------------------------------------------
  real(rp) :: prim(neq, -1:nx+2, -1:ny+2)     ! Primitive variables
  real(rp) ::    u(neq, -1:nx+2, -1:ny+2)     ! Current conserved variables
  real(rp) ::   up(neq, -1:nx+2, -1:ny+2)     ! Advanced conserved variables
  real(rp) ::    f(neq, -1:nx+2, -1:ny+2)     ! Physical Fluxes in X
  real(rp) ::    g(neq, -1:nx+2, -1:ny+2)     ! Physical Fluxes in Y

  !--------------------------------------------------------------------------
  ! Time and iteration variables
  !--------------------------------------------------------------------------
  integer :: currentIteration     ! Current number of iterations
  integer :: nextFileNumber       ! Number of the next output
  real(rp) :: currentTime         ! Current time in the simulation
  real(rp) :: timeToPrint         ! Time to next output
  real(rp) :: dt                  ! Current valid time step

  !--------------------------------------------------------------------------
  ! Variables to meassure the excecution time
  !--------------------------------------------------------------------------
  integer :: clockStart, clockCount, clockRate, clockMax


end module Globals
