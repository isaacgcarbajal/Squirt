program Main

  use Globals, only: currentTime, endTime, timeToPrint
  use Init, only: initialConditions
  use Utilities, only: writeToDisk, initializeClockVariables, finishClockVariables
  use Solver, only: solutionStep
  use HydroCore, only: validTimeStep
  implicit none
  
  call initialConditions()
  
  call initializeClockVariables()
  do while(currentTime < endTime)

    call validTimeStep()
    
    call solutionStep()

    if(currentTime >= timeToPrint) then
      call writeToDisk()
    end if
    
  end do
  call finishClockVariables()

end program Main
