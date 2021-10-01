module Utilities

  contains
  
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

  subroutine initializeClockVariables()
    
    use Globals, only: clockStart, clockRate, clockMax
    implicit none
    
    call system_clock(clockStart, clockRate, clockMax)
  
  end subroutine initializeClockVariables
  
  subroutine finishClockVariables()
    use Globals, only: clockStart, clockRate, clockMax, clockCount, currentIteration
  
    call system_clock(clockCount, clockRate, clockMax)
    write(*, '(a, i5,a,f8.3,a)') "Where calculated ", currentIteration, " iterations in ",   &
                                  (clockCount-clockStart)/real(clockRate), " s"
  end subroutine finishClockVariables

  subroutine swapXY(var)
    use Globals, only: neq
    real*8, intent(inout) :: var(neq)
    real*8 :: aux
    
    aux = var(2)
    var(2) = var(3)
    var(3) = aux    
    
  end subroutine swapXY

end module Utilities
