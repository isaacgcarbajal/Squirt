module HLLC

  contains
  
  subroutine primitive2HLLCFlux(ppL, ppR, ff)
  
    use Globals, only: neq, gasGamma
    use HydroCore, only: soundSpeed, primitive2EulerFluxes, primitive2conserved
    implicit none
    real*8, intent(in) :: ppL(neq), ppR(neq)
    real*8, intent(out):: ff(neq)
    
    real*8 :: csL, csR, sL, sR, sst
    real*8 :: rhost, Ek
    
    real*8 :: uu(neq), uuk(neq)
    
    call soundSpeed(ppL(1), ppL(4), csL)
    call soundSpeed(ppR(1), ppR(4), csR)
    
    sL = min(ppL(2)-csL, ppR(2)-csR)
    sR = max(ppL(2)+csL, ppR(2)+csR)
    
    if(sL > 0.0) then
      call primitive2EulerFluxes(ppL,ff)
      return
    end if
    
    if(sR < 0.0) then
      call primitive2EulerFluxes(ppR,ff)
      return
    end if
    
    sst = (ppR(4) - ppL(4) + ppL(1)*ppL(2)*(sL-ppL(2)) - ppR(1)*ppR(2)*(sR-ppR(2)))   &
          / (ppL(1)*(sL-ppL(2)) - ppR(1)*(sR-ppR(2)))
          
    if(sst >= 0.0) then
    
      rhost = ppL(1)*(sL - ppL(2))/(sL - sst)
      Ek = 0.5*ppL(1)*(ppL(2)**2 + ppL(3)**2) + ppL(4)/(gasGamma-1)
    
      uuk(1) = rhost
      uuk(2) = rhost*sst
      uuk(3) = rhost*ppL(3)
      uuk(4) = rhost*( Ek/ppL(1)+(sst-ppL(2))*(sst+ppL(4)/(ppL(1)*(sL-ppL(2)))) )
      
      call primitive2EulerFluxes(ppL, ff)
      call primitive2conserved(ppL, uu)
      
      ff(:) = ff(:) + sL*(uuk(:) - uu(:))
      
      return
    end if
    
    if(sst <= 0.0) then
    
      rhost = ppR(1)*(sR - ppR(2))/(sR - sst)
      Ek = 0.5*ppR(1)*(ppR(2)**2 + ppR(3)**2) + ppR(4)/(gasGamma-1)
    
      uuk(1) = rhost
      uuk(2) = rhost*sst
      uuk(3) = rhost*ppR(3)
      uuk(4) = rhost*( Ek/ppR(1)+(sst-ppR(2))*(sst+ppR(4)/(ppR(1)*(sR-ppR(2)))) )
      
      call primitive2EulerFluxes(ppR, ff)
      call primitive2conserved(ppR, uu)
      
      ff(:) = ff(:) + sR*(uuk(:) - uu(:))
      
      return
    end if
  
  end subroutine primitive2HLLCFlux

  subroutine fullHLLCFlux()
    
    use Globals, only: neq, nx, ny, prim, F, G
    use Utilities, only: swapXY
    implicit none
    
    real*8 :: ppL(neq), ppR(neq)
    real*8 :: tmpFlux(neq)
    
    integer :: i, j
    do j=1,ny+1
      do i=1,nx+1
        
        ! Fluxes in X
        ppL(:) = prim(:, i, j)
        ppR(:) = prim(:, i+1, j)
        
        call primitive2HLLCFlux(ppL, ppR, tmpFlux)
        
        F(:,i,j) = tmpFlux(:)
        
        ! Fluxes in Y
        ppL(:) = prim(:, i, j)
        ppR(:) = prim(:, i, j+1)
        call swapXY(ppL)
        call swapXY(ppR)
        
        call primitive2HLLCFlux(ppL, ppR, tmpFlux)
        call swapXY(tmpFlux)
        
        G(:,i,j) = tmpFlux(:)

        
      end do
    end do
    
  end subroutine fullHLLCFlux

end module HLLC
