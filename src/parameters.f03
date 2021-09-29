module Parameters
  implicit none
  
  !--------------------------------------------------------------------------
  ! Space parameters
  !--------------------------------------------------------------------------
  ! Number of nodes in X and Y
  integer, parameter :: nx = 400
  integer, parameter :: ny = 400
  ! Domain delimitation in X
  real*8, parameter :: xMin = 0.0
  real*8, parameter :: xMax = 1.0
  ! Domain delimitation in Y
  real*8, parameter :: yMin = 0.0
  real*8, parameter :: yMax = 1.0

  !--------------------------------------------------------------------------
  ! Time parameters
  !--------------------------------------------------------------------------
  real*8, parameter :: endTime = 10.0
  real*8, parameter :: dtToPrint = 0.1

  !--------------------------------------------------------------------------
  ! Physical parameters
  !--------------------------------------------------------------------------
  ! Corresponding Gas Gamma parameter
  real*8, parameter :: gasGamma = 1.4

  !--------------------------------------------------------------------------
  ! Numerical parameters
  !--------------------------------------------------------------------------
  ! Courant–Friedrichs–Lewy paramter
  real*8, parameter :: CFL = 0.99 
  ! Artifitial viscosity
  real*8, parameter :: eta = 0.001
  
  !--------------------------------------------------------------------------
  ! Problem specific parameters
  !--------------------------------------------------------------------------
  real*8, parameter :: A=0.01
  
  !--------------------------------------------------------------------------
  ! Boundary conditions:
  ! 1 .- Dirichlet
  ! 2 .- Newmann BC
  ! 3 .- Periodic
  !--------------------------------------------------------------------------
  !integer :: topBC    = 2
  !integer :: bottomBC = 2
  !integer :: leftBC   = 3
  !integer :: rightBC  = 3

end module Parameters