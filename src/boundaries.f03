module Boundaries

  contains
  
  !============================================================================
  ! boundaryConditionsI:
  !   Subroutine to set a first order boundary conditions on a one ghost cell
  !   per side of the domain.
  ! Input:
  !   u(neq, -1:nx+2, -1:ny+2) -> Array in which the boundary condicions will
  !                               be imposed.
  !============================================================================
  subroutine boundaryConditionsI(u)
  
    use Globals, only: rp, neq, nx, ny
    real(rp), intent(inout) :: u(neq, -1:nx+2, -1:ny+2)
    
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
  
  !============================================================================
  ! boundaryConditionsII:
  !   Subroutine to set a second order boundary conditions on two ghost cells
  !   per side of the domain.
  ! Input:
  !   u(neq, -1:nx+2, -1:ny+2) -> Array in which the boundary condicions will
  !                               be imposed.
  !============================================================================
  subroutine boundaryConditionsII(u)
  
    use Globals, only: rp, neq, nx, ny
    real(rp), intent(inout) :: u(neq, -1:nx+2, -1:ny+2)
    
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
  
  !============================================================================
  ! applyNeumann:
  !   Subroutine to impose a Neumann boundary condition in one cell.
  ! Input:
  !   physical -> The physical cell to be copied
  ! Output:
  !   ghost    -> The ghost cell filled
  !============================================================================
  subroutine applyNeumann(ghost, physical)
    
    use Globals, only: rp
    implicit none
    real(rp), intent(out) :: ghost
    real(rp), intent(in)  :: physical
    
    ghost = physical
  
  end subroutine applyNeumann

  !============================================================================
  ! applyPeriodic:
  !   Subroutine to impose a periodic boundary condition over two cells.
  ! Input:
  !   physicalL -> The physical cell in the "left" to be copied
  !   physicalR -> The physical cell in the "right" to be copied
  ! Output:
  !   ghostL    -> The ghost "left" cell filled
  !   ghostR    -> The ghost "right" cell filled
  !============================================================================
  subroutine applyPeriodic(ghostL, physicalL, physicalR, ghostR)
    
    use Globals, only: rp
    implicit none
    real(rp), intent(out) :: ghostL, ghostR
    real(rp), intent(in) :: physicalL, physicalR
    
    ghostL = physicalR
    ghostR = physicalL
    
  end subroutine applyPeriodic

end module Boundaries
