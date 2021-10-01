module HydroCore

  contains
  
  subroutine conserved2primitive(uu,pp)
  
    use Globals, only: neq, nx, ny, gasGamma
    implicit none
    real*8, intent(in)  :: uu(neq)
    real*8, intent(out) :: pp(neq)
    
    pp(1) = uu(1)
    pp(2) = uu(2)/uu(1)
    pp(3) = uu(3)/uu(1)
    pp(4) = (gasGamma - 1.0)*(uu(4) - 0.5*(uu(2)*pp(2) + uu(3)*pp(3)))
  
  end subroutine conserved2primitive
  
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

  subroutine primitive2conserved(pp,uu)
  
    use Globals, only: neq, nx, ny, gasGamma
    implicit none
    real*8, intent(in)  :: pp(neq)
    real*8, intent(out) :: uu(neq)
    
    uu(1) = pp(1)
    uu(2) = pp(1)*pp(2)
    uu(3) = pp(1)*pp(3)
    uu(4) = 0.5*pp(1)*(pp(2)**2 + pp(3)**2) + pp(4)/(gasGamma - 1.0)
  
  end subroutine primitive2conserved
  
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

  subroutine updatePrimitivesWith(u)
    
    use Globals, only: neq, nx, ny, prim, gasGamma
    implicit none
    
    real*8, intent(in) :: u(neq,-1:nx+2,-1:ny+2)

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

  subroutine soundSpeed(rho, P, cs)
    
    use Globals, only: gasGamma
    implicit none
    real*8, intent(in) :: rho, P
    real*8, intent(out) :: cs
    
    cs = sqrt(gasGamma*P/rho)
    
  end subroutine soundSpeed

  subroutine validTimeStep()
    
    use Globals, only: nx, ny, prim, dt, CFL, dx, dy
    implicit none
    
    real*8 :: cs
    real*8 :: maxSpeedX = 0.0
    real*8 :: maxSpeedY = 0.0
    
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
  
  subroutine primitive2EulerFluxes(pp,ff)
  
    use Globals, only: neq, gasGamma
    implicit none
    
    real*8, intent(in)  :: pp(neq)
    real*8, intent(out) :: ff(neq)
    
    ff(1) = pp(1)*pp(2)
    ff(2) = pp(1)*pp(2)**2 + pp(4)
    ff(3) = pp(1)*pp(2)*pp(3)
    ff(4) = pp(2)*(0.5*pp(1)*(pp(2)**2+pp(3)**2)+pp(4)*gasGamma/(gasGamma-1.0))
    
  end subroutine primitive2EulerFluxes

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

  subroutine slopeLimiter(ppll, ppl, ppr, pprr)
  
    use Globals, only: neq
    implicit none
    real*8, intent(in)   :: ppll(neq), pprr(neq)
    real*8, intent(inout):: ppl(neq), ppr(neq)
    
    real*8 :: deltaL, deltaM, deltaR
    real*8 :: avrgL, avrgR
    
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
    
    real function average(a,b)
      implicit none
      real*8 :: s

      real*8, intent(in) :: a, b
    
      ! minmod limiter
      s = sign(1.0D0,a)
      average = s*max(0.0,min(abs(a),s*b))
    
    end function average
  
  end subroutine slopeLimiter

end module HydroCore
