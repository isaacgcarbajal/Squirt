module Parameters
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
  ! Space parameters
  !--------------------------------------------------------------------------
  ! Number of nodes in X and Y
  integer, parameter :: nx = 256
  integer, parameter :: ny = 256
  ! Domain delimitation in X
  real(rp), parameter :: xMin = 0.0
  real(rp), parameter :: xMax = 1.0
  ! Domain delimitation in Y
  real(rp), parameter :: yMin = 0.0
  real(rp), parameter :: yMax = 1.0

  !--------------------------------------------------------------------------
  ! Time parameters
  !--------------------------------------------------------------------------
  real(rp), parameter :: endTime = 10.0
  real(rp), parameter :: dtToPrint = 0.1

  !--------------------------------------------------------------------------
  ! Physical parameters
  !--------------------------------------------------------------------------
  ! Corresponding Gas Gamma parameter
  real(rp), parameter :: gasGamma = 1.4

  !--------------------------------------------------------------------------
  ! Numerical parameters
  !--------------------------------------------------------------------------
  ! Courant–Friedrichs–Lewy paramter
  real(rp), parameter :: CFL = 0.99 
  ! Artifitial viscosity
  real(rp), parameter :: eta = 0.001
  
  !--------------------------------------------------------------------------
  ! Problem specific parameters
  !--------------------------------------------------------------------------
  real(rp), parameter :: A=0.01
  

end module Parameters
