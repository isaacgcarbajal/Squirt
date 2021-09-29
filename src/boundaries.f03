module Boundaries

  contains
  
  subroutine boundaryConditions(u)
  
    use Globals, only: neq, nx, ny
    real*8, intent(inout) :: u(neq, nx+2, ny+2)
    
    integer ieq, i, j
    
    ! Top and bottom Neumann
    do i=2, nx+1
      do ieq=1, neq
      
        !call applyNeumann(u(ieq,i,1), u(ieq,i,2))
        !call applyNeumann(u(ieq,i,ny+2), u(ieq,i,ny+1))
        
        u(ieq,i,1) = u(ieq,i,2)
        u(ieq,i,ny+2) = u(ieq,i,ny+1)
        
      end do
    end do
    
    ! Left and right Periodic
    do j=2, ny+1
      do ieq=1, neq
      
        !call applyPeriodic(u(ieq,1,j), u(ieq,2,j), u(ieq,nx+1,j), u(ieq,nx+2,j))
        
        u(ieq,1,j) = u(ieq,nx+1,j)
        u(ieq,nx+2,j) = u(ieq,2,j)
      
      end do
    end do
    
  end subroutine boundaryConditions
  
  subroutine applyNeumann(ghost, first)
    
    implicit none
    real*8, intent(out) :: ghost
    real*8, intent(in)  :: first
    
    ghost = first
  
  end subroutine applyNeumann

  subroutine applyPeriodic(ghostL, first, last, ghostR)
    
    implicit none
    real*8, intent(out) :: ghostL, ghostR
    real*8, intent(in) :: first, last
    
    ghostL = last
    ghostR = first
    
  end subroutine applyPeriodic

end module Boundaries
