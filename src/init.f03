module Init

  contains
  
  subroutine initialConditions()
  
    use Globals, only: neq, nx, ny, xMin, yMin, dx, dy,               &
                       currentTime, currentIteration, nextFileNumber, &
                       timeToPrint, A, gasGamma, prim, u
    use HydroCore, only: fullPrimitive2conserved
    use Utilities, only: writeToDisk
    implicit none
  
    call initializePrimitives()
    call fullPrimitive2conserved()
    call initializeTimeVariables()
    
    call writeToDisk()
    
  end subroutine initialConditions
  
  subroutine initializePrimitives()
  
    use Globals, only: nx, ny, yMin, dy, A, prim
    implicit none
    
    integer :: i,j
    do j=-1,ny+2
      do i=-1,nx+2
        if(yMin+(j-1)*dy <= 0.5) then
          prim(1,i,j) = 2.0
          prim(2,i,j) = 0.5
          prim(3,i,j) = 0.0
          prim(4,i,j) = 2.5
        else
          prim(1,i,j) = 1.0
          prim(2,i,j) =-0.5
          prim(3,i,j) = 0.0
          prim(4,i,j) = 2.5
        end if
        ! Applying the perturbation
        prim(2,i,j) = prim(2,i,j) + A*(rand() - 0.5)
        prim(3,i,j) = prim(3,i,j) + A*(rand() - 0.5)
      end do
    end do
    
  end subroutine initializePrimitives
  
  subroutine initializeTimeVariables()
    use Globals, only: currentTime, currentIteration, nextFileNumber, timeToPrint
    implicit none
    
    currentTime = 0.0
    currentIteration = 0
    nextFileNumber = 0
    timeToPrint = 0.0
  
  end subroutine initializeTimeVariables


end module Init
