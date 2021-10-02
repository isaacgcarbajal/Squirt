module Solver
  
  contains
  
  !============================================================================
  ! stepAuxiliarConserved:
  !   Subroutine to update the value of the array of auxiliar conserved
  !   variables for a given "dt" assuming that the fluxes "f" and "g" have
  !   already been calculated.
  ! Input:
  !   dt -> The time step of the update.
  !============================================================================
  subroutine stepAuxiliarConserved(dt)
    use Globals, only: neq, nx, ny, dx, dy, u, up, f, g
    implicit none
    
    real*8, intent(in) :: dt
    
    integer :: i, j
    do j=1,ny
      do i=1,nx
      
        up(:,i,j) = u(:,i,j) - dt*(f(:,i,j) - f(:,i-1,j))/dx  &
                             - dt*(g(:,i,j) - g(:,i,j-1))/dy
      end do
    end do
    
  end subroutine stepAuxiliarConserved

  !============================================================================
  ! solutionStep:
  !   Subroutine to make a full step of the simulation. Given an initial state
  !   it calculates the next value for the conserved variables and updates the
  !   current time and current iteration of the simulation.
  !============================================================================
  subroutine solutionStep()
    use Globals, only: dt, u, up, currentIteration, currentTime
    use HydroCore, only: updatePrimitivesWith, fullPrimitive2conserved, addArtifitialViscosity
    use Boundaries, only: boundaryConditionsI, boundaryConditionsII
    use HLLC, only: fullHLLCFlux
    implicit none
    
    integer :: halfStep
    
    halfStep = 1
    call fullHLLCFlux(halfStep)
    call stepAuxiliarConserved(dt/2)
    call boundaryConditionsII(up)
    call updatePrimitivesWith(up)
    
    halfStep = 2
    call fullHLLCFlux(halfStep)
    call stepAuxiliarConserved(dt)
    call addArtifitialViscosity()
    call boundaryConditionsI(u)
    call updatePrimitivesWith(u)

    currentTime = currentTime + dt
    currentIteration = currentIteration + 1

  end subroutine solutionStep

end module Solver
