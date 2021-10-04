module HydroCore

  contains
  
  !============================================================================
  ! conserved2primitive:
  !   Subroutine to calculate the primitive variables given the conservated in
  !   one cell of the domain.
  ! Input:
  !   uu(neq) -> Conserved variables in the cell.
  ! Output:
  !   pp(neq) -> Primitive variables in the cell.
  !============================================================================
  subroutine conserved2primitive(uu,pp)
  
    use Globals, only: rp, neq, nx, ny, gasGamma
    implicit none
    real(rp), intent(in)  :: uu(neq)
    real(rp), intent(out) :: pp(neq)
    
    pp(1) = uu(1)
    pp(2) = uu(2)/uu(1)
    pp(3) = uu(3)/uu(1)
    pp(4) = (gasGamma - 1.0)*(uu(4) - 0.5*(uu(2)*pp(2) + uu(3)*pp(3)))
  
  end subroutine conserved2primitive
  
  !============================================================================
  ! fullConserved2primitive:
  !   Subroutine to calculate the primitive variables on the entire domain
  !   given the conservated variables. It makes use of the conserved2primitive
  !   subroutine to calculate the variables in each cell.
  !============================================================================
  subroutine fullConserved2primitive()
  
    use Globals, only: nx, ny, u, prim
    implicit none
    
    integer :: i, j
    do j=-1,ny+2
      do i=-1,nx+2
      
        call conserved2primitive(u(:,i,j), prim(:,i,j))
      
      end do
    end do
    
  end subroutine fullConserved2primitive
  
  !============================================================================
  ! primitive2conserved:
  !   Subroutine to calculate the conservated variables given the primitives in
  !   one cell of the domain.
  ! Input:
  !   pp(neq) -> Primitive variables in the cell.
  ! Output:
  !   uu(neq) -> Conserved variables in the cell.
  !============================================================================
  subroutine primitive2conserved(pp,uu)
  
    use Globals, only: rp, neq, nx, ny, gasGamma
    implicit none
    real(rp), intent(in)  :: pp(neq)
    real(rp), intent(out) :: uu(neq)
    
    uu(1) = pp(1)
    uu(2) = pp(1)*pp(2)
    uu(3) = pp(1)*pp(3)
    uu(4) = 0.5*pp(1)*(pp(2)**2 + pp(3)**2) + pp(4)/(gasGamma - 1.0)
  
  end subroutine primitive2conserved
  
  !============================================================================
  ! fullPrimitive2conserved:
  !   Subroutine to calculate the conserved variables on the entire domain
  !   given the primitive variables. It makes use of the primitive2conserved
  !   subroutine to calculate the variables in each cell.
  !============================================================================
  subroutine fullPrimitive2conserved()
  
    use Globals, only: nx, ny, u, prim
    implicit none
    
    integer :: i, j
    do j=-1,ny+2
      do i=-1,nx+2
      
        call primitive2conserved(prim(:,i,j), u(:,i,j))
      
      end do
    end do
    
  end subroutine fullPrimitive2conserved

  !============================================================================
  ! updatePrimitivesWith:
  !   Subroutine to update the value of the global array of primitive variables
  !   using an array of conservated variables as an input.
  ! Input:
  !   u(neq,-1:nx+2,-1:ny+2) -> Array of conservated variables in the entire
  !                             domain.
  !============================================================================
  subroutine updatePrimitivesWith(u)
    
    use Globals, only: rp, neq, nx, ny, prim, gasGamma
    implicit none
    
    real(rp), intent(in) :: u(neq,-1:nx+2,-1:ny+2)

    integer :: i, j
    do j=-1,ny+2
      do i=-1,nx+2
      
        prim(1,i,j) = u(1,i,j)
        prim(2,i,j) = u(2,i,j)/u(1,i,j)
        prim(3,i,j) = u(3,i,j)/u(1,i,j)
        prim(4,i,j) = (gasGamma - 1.0)*(u(4,i,j) - 0.5*(u(2,i,j)*prim(2,i,j) + u(3,i,j)*prim(3,i,j)))
      
      end do
    end do    
    
  end subroutine updatePrimitivesWith

  !============================================================================
  ! soundSpeed:
  !   Subroutine to calculate the sound speed for an ideal gas given the
  !   density and the pressure.
  ! Input:
  !   rho -> Local density
  !   P   -> Local pressure
  ! Output:
  !   cs  -> Sound speed
  !============================================================================
  subroutine soundSpeed(rho, P, cs)
    
    use Globals, only: rp, gasGamma
    implicit none
    real(rp), intent(in) :: rho, P
    real(rp), intent(out) :: cs
    
    cs = sqrt(gasGamma*P/rho)
    
  end subroutine soundSpeed

  !============================================================================
  ! validTimeStep:
  !   Subroutine to calculate the current valid time step for the simulation.
  !   Is related to how fast the information travels across the domain.
  !   This subroutine updates the value of the global variable "dt"
  !============================================================================
  subroutine validTimeStep()
    
    use Globals, only: rp, nx, ny, prim, dt, CFL, dx, dy
    implicit none
    
    real(rp) :: cs
    real(rp) :: maxSpeedX = 0.0
    real(rp) :: maxSpeedY = 0.0
    
    integer :: i, j
    do j=-1,ny+2
      do i=-1,nx+2
      
        call soundSpeed(prim(1,i,j), prim(4,i,j), cs)
        
        maxSpeedX = max(maxSpeedX, cs+abs(prim(2,i,j)))
        maxSpeedY = max(maxSpeedY, cs+abs(prim(3,i,j)))
        
      end do
    end do    
    
    if (maxSpeedX > maxSpeedY) then
      dt = CFL*dx/(sqrt(2.0)*maxSpeedX)
    else
      dt = CFL*dy/(sqrt(2.0)*maxSpeedY)
    end if
    
  end subroutine validTimeStep
  
  !============================================================================
  ! primitive2EulerFluxes:
  !   Subroutine to calculate the fluxes of the Euler equations via the
  !   primitive variables in one cell.
  ! Input:
  !   pp(neq) -> Primitive variables on the cell.
  ! Output:
  !   ff(neq) -> Euler flux on the cell.
  !============================================================================
  subroutine primitive2EulerFluxes(pp,ff)
  
    use Globals, only: rp, neq, gasGamma
    implicit none
    
    real(rp), intent(in)  :: pp(neq)
    real(rp), intent(out) :: ff(neq)
    
    ff(1) = pp(1)*pp(2)
    ff(2) = pp(1)*pp(2)**2 + pp(4)
    ff(3) = pp(1)*pp(2)*pp(3)
    ff(4) = pp(2)*(0.5*pp(1)*(pp(2)**2+pp(3)**2)+pp(4)*gasGamma/(gasGamma-1.0))
    
  end subroutine primitive2EulerFluxes

  !============================================================================
  ! addArtifitialViscosity:
  !   Subroutine to copy the auxiliar array of conserved variables "up" into
  !   the array of conserved variables "u" adding artifitial viscosity moduled
  !   by the parameter "eta".
  !============================================================================
  subroutine addArtifitialViscosity()
  
    use Globals, only: nx, ny, u, up, eta
    implicit none
    
    integer :: i, j
    do j=1,ny
      do i=1,nx
      
        u(:,i,j) = up(:,i,j) + eta*(  up(:,i-1,j) + up(:,i+1,j)  &
                                    + up(:,i,j-1) + up(:,i,j+1)  &
                                    - 4*up(:,i,j))
      end do
    end do
  
  end subroutine addArtifitialViscosity

  !============================================================================
  ! slopeLimiter:
  !   Subroutine to restrict the values of the primitives on the left and right
  !   sides of a cell. It is necessary two ghost cells on each frontier of the
  !   domain. Makes use of the "minmod" limiter.
  ! Input:
  !   ppll(neq) -> Primitive variables on the left of the left side of the cell
  !   ppl(neq)  -> Primitive variables on the left side of the cell
  !   ppr(neq)  -> Primitive variables on the right side of the cell
  !   pprr(neq) -> Primitive variables on the right of the right side of the
  !                cell
  !============================================================================
  subroutine slopeLimiter(ppll, ppl, ppr, pprr)
  
    use Globals, only: rp, neq
    implicit none
    real(rp), intent(in)   :: ppll(neq), pprr(neq)
    real(rp), intent(inout):: ppl(neq), ppr(neq)
    
    real(rp) :: deltaL, deltaM, deltaR
    real(rp) :: avrgL, avrgR
    
    integer::ieq
      
    do ieq=1,neq
      
      deltaL = ppl(ieq)  - ppll(ieq)
      deltaM = ppr(ieq)  - ppl(ieq)
      deltaR = pprr(ieq) - ppr(ieq)
      
      avrgL = average(deltaL, deltaM)
      avrgR = average(deltaM, deltaR)
      
      ppl(ieq) = ppl(ieq) + 0.5*avrgL
      ppR(ieq) = ppR(ieq) - 0.5*avrgR
    
    end do
    
    contains
    
    real(rp) function average(a,b)
      implicit none
      real(rp) :: s

      real(rp), intent(in) :: a, b
    
      ! minmod limiter
      s = sign(1.0D0,a)
      average = s*max(0.0,min(abs(a),s*b))
    
    end function average
  
  end subroutine slopeLimiter

end module HydroCore
