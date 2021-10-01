module Solver
  
  contains
  
  subroutine stepConserved(dt)
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
    
  end subroutine stepConserved

  subroutine solutionStep()
    use Globals, only: dt, u, up, currentIteration, currentTime
    use HydroCore, only: updatePrimitivesWith, fullPrimitive2conserved, addArtifitialViscosity
    use Boundaries, only: boundaryConditionsI, boundaryConditionsII
    use HLLC, only: fullHLLCFlux
    implicit none
    
    call fullHLLCFlux(1)
    call stepConserved(dt/2)
    call boundaryConditionsII(up)
    call updatePrimitivesWith(up)
    
    call fullHLLCFlux(2)
    call stepConserved(dt)
    call addArtifitialViscosity()
    call boundaryConditionsI(u)
    call updatePrimitivesWith(u)

    currentTime = currentTime + dt
    currentIteration = currentIteration + 1

  end subroutine solutionStep

end module Solver
