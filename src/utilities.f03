module Utilities

  contains
  
  !============================================================================
  ! writeToDisk:
  !   Subroutine to write to disk a binary of the current array of primitives.
  !   Also it updates the value of the next file and the time to  print the
  !   next file.
  !   It prints to screen the file name of the binary file.
  !============================================================================
  subroutine writeToDisk()
    use Globals, only: neq, nx, ny, nextFileNumber, timeToPrint, dtToPrint, prim
    implicit none

    character (len=80) :: fileName
    integer :: i,j, ieq

    ! Generate the name of the output file
    write(fileName, "(a,i3.3,a)") "HLLC_", nextFileNumber, ".dat"

    ! Open the file
    open(unit=10, file=fileName, access="stream")

    ! Write the values of the primitives to the file
    do ieq=1,neq
      do i=1,nx
        do j=1,ny
        
          write(10) prim(ieq, i, j)
          
        end do
      end do
    end do

    ! Close the file
    close(10)

    write(*,'(a,a)') "Written ", trim(fileName)

    ! Update iteration variables
    nextFileNumber = nextFileNumber + 1
    timeToPrint = timeToPrint + dtToPrint
  
  end subroutine writeToDisk

  !============================================================================
  ! initializeClockVariables:
  !   Subroutine to initialize the time variables via the "system_clock"
  !   function in order to meassure the physical time taken by the program to
  !   finish.
  !============================================================================
  subroutine initializeClockVariables()
    
    use Globals, only: clockStart, clockRate, clockMax
    implicit none
    
    call system_clock(clockStart, clockRate, clockMax)
  
  end subroutine initializeClockVariables
  
  !============================================================================
  ! finishClockVariables:
  !   Subroutine to calculate the elapsed time from the beginning of the
  !   simulation to the end. Also prints in screen the number of iterations
  !   and the total physical elapsed time.
  !============================================================================
  subroutine finishClockVariables()
    use Globals, only: clockStart, clockRate, clockMax, clockCount, currentIteration
  
    call system_clock(clockCount, clockRate, clockMax)
    write(*, '(a, i5,a,f8.3,a)') "Where calculated ", currentIteration, " iterations in ",   &
                                  (clockCount-clockStart)/real(clockRate), " s"
  end subroutine finishClockVariables

  !============================================================================
  ! swapXY:
  !   Subroutine to swap the value of the variables. It is used in order to
  !   calculate the fluxes with only one "x" function, swapping the X and Y
  !   components of the simulation and then swapping again once the fluxes have
  !   been calculated.
  !============================================================================
  subroutine swapXY(var)
    use Globals, only: rp, neq
    real(rp), intent(inout) :: var(neq)
    real(rp) :: aux
    
    aux = var(2)
    var(2) = var(3)
    var(3) = aux    
    
  end subroutine swapXY

end module Utilities
