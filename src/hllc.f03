module HLLC

  contains
  
  !============================================================================
  ! primitive2HLLCFlux:
  !   Subroutine to calculate the HLLC flux from the primitive vatiables in one
  !   cell.
  ! Input:
  !   ppL(neq) -> Primitive variables on the left side of the cell
  !   ppR(neq) -> Primitive variables on the right side of the cell
  ! Output:
  !   ff(neq)  -> The HLLC fluxes in the cell
  !============================================================================
  subroutine primitive2HLLCFlux(ppL, ppR, ff)
  
    use Globals, only: rp, neq, gasGamma
    use HydroCore, only: soundSpeed, primitive2EulerFluxes, primitive2conserved
    implicit none
    real(rp), intent(in) :: ppL(neq), ppR(neq)
    real(rp), intent(out):: ff(neq)
    
    real(rp) :: csL, csR, sL, sR, sst
    real(rp) :: rhost, Ek
    
    real(rp) :: uu(neq), uuk(neq)
    
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

  !============================================================================
  ! fullHLLCFlux:
  !   Subroutine to calculate the HLLC flux over the entire domain. In the case
  !   of the first half time step there is no need to apply the slope limiters,
  !   but in the second they are needed in order to preserve the second order
  !   accuracy.
  ! Input:
  !   halfStep -> the number of "half step" of the solution. Is used to select
  !               a case inside the subroutine.
  !============================================================================
  subroutine fullHLLCFlux(halfStep)
    
    use Globals, only: rp, neq, nx, ny, prim, F, G
    use Utilities, only: swapXY
    use HydroCore, only: slopeLimiter
    implicit none
    
    integer, intent(in) :: halfStep
    
    real(rp) :: ppL(neq), ppR(neq)
    real(rp) :: ppLL(neq), ppRR(neq)
    real(rp) :: tmpFlux(neq)
    
    integer :: i, j
    
    select case(halfStep)
    
    case(1)
      do j=0,ny
        do i=0,nx
          
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
      
    case(2)
      do j=0,ny
        do i=0,nx
          
          ! Fluxes in X
          ppLL(:) = prim(:, i-1, j)
          ppL(:)  = prim(:, i  , j)
          ppR(:)  = prim(:, i+1, j)
          ppRR(:) = prim(:, i+2, j)
          call slopeLimiter(ppLL, ppL, ppR, ppRR)
          
          call primitive2HLLCFlux(ppL, ppR, tmpFlux)
          
          F(:,i,j) = tmpFlux(:)
          
          ! Fluxes in Y
          ppLL(:) = prim(:, i, j-1)
          ppL(:)  = prim(:, i, j  )
          ppR(:)  = prim(:, i, j+1)
          ppRR(:) = prim(:, i, j+2)
          call swapXY(ppLL)
          call swapXY(ppL)
          call swapXY(ppR)
          call swapXY(ppRR)
          call slopeLimiter(ppLL, ppL, ppR, ppRR)
          
          call primitive2HLLCFlux(ppL, ppR, tmpFlux)
          call swapXY(tmpFlux)
          
          G(:,i,j) = tmpFlux(:)
    
          
        end do
      end do
      
    end select
      
  end subroutine fullHLLCFlux

end module HLLC
