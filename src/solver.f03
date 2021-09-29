module Solver
  
  contains
  
  subroutine stepConserved(dt)
    use Globals, only: neq, nx, ny, dx, dy, u, up, f, g
    implicit none
    
    real*8, intent(in) :: dt
    
    integer :: i, j
    do j=2,ny+1
      do i=2,nx+1
      
        up(:,i,j) = u(:,i,j) - dt*(f(:,i,j) - f(:,i-1,j))/dx  &
                             - dt*(g(:,i,j) - g(:,i,j-1))/dy
      end do
    end do
    
  end subroutine stepConserved

  subroutine solutionStep()
    use Globals, only: dt, u, up, currentIteration, currentTime
    use HydroCore, only: updatePrimitivesWith, fullPrimitive2conserved, addArtifitialViscosity
    use Boundaries, only: boundaryConditions
    use HLLC, only: fullHLLCFlux
    implicit none
    
    call fullHLLCFlux()
    call stepConserved(dt/2)
    call boundaryConditions(up)
    call updatePrimitivesWith(up)
    
    call fullHLLCFlux()
    call stepConserved(dt)
    call addArtifitialViscosity()
    call boundaryConditions(u)
    call updatePrimitivesWith(u)

    currentTime = currentTime + dt
    currentIteration = currentIteration + 1

  end subroutine solutionStep

end module Solver
