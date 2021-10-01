module Boundaries

  contains
  
  subroutine boundaryConditionsI(u)
  
    use Globals, only: neq, nx, ny
    real*8, intent(inout) :: u(neq, -1:nx+2, -1:ny+2)
    
    integer ieq, i, j
    
    ! Top and bottom Neumann
    do i=1, nx
      do ieq=1, neq
      
        call applyNeumann(u(ieq,i,0), u(ieq,i,1))
        call applyNeumann(u(ieq,i,ny+1), u(ieq,i,ny))
        
      end do
    end do
    
    ! Left and right Periodic
    do j=1, ny
      do ieq=1, neq
      
        call applyPeriodic(u(ieq,0,j), u(ieq,1,j), u(ieq,nx,j), u(ieq,nx+1,j))
      
      end do
    end do
    
  end subroutine boundaryConditionsI
  
  subroutine boundaryConditionsII(u)
  
    use Globals, only: neq, nx, ny
    real*8, intent(inout) :: u(neq, -1:nx+2, -1:ny+2)
    
    integer ieq, i, j
    
    ! Top and bottom Neumann
    do i=1, nx
      do ieq=1, neq
      
        call applyNeumann(u(ieq,i,0), u(ieq,i,1))
        call applyNeumann(u(ieq,i,-1), u(ieq,i,0))
        
        call applyNeumann(u(ieq,i,ny+1), u(ieq,i,ny))
        call applyNeumann(u(ieq,i,ny+2), u(ieq,i,ny+1))
        
      end do
    end do
    
    ! Left and right Periodic
    do j=1, ny
      do ieq=1, neq
      
        call applyPeriodic(u(ieq,-1,j), u(ieq,1,j), u(ieq,nx-2,j), u(ieq,nx+1,j))
        call applyPeriodic(u(ieq, 0,j), u(ieq,2,j), u(ieq,nx  ,j), u(ieq,nx+2,j))
      
      end do
    end do
    
  end subroutine boundaryConditionsII
  
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
