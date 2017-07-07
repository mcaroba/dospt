module good_bye

  use read_input

  implicit none

  contains

  subroutine print_good_bye()
!******************************************************
! Print message at the end depending on execution
! success
  if(execution_error == 0)then
    write(*,*)'                                       |'
    write(*,*)'Program successfully executed          |'
    write(*,*)'                                       |'
    write(*,*)'End of execution:                      |'
    call timestring(time)
    write(*,'(1X,A39,A)')time,'|'
  else if(execution_error > 0)then
    write(*,*)'                                       |'
    write(*,*)'       !!!!!!  ERROR  !!!!!!           |'
    write(*,*)'                                       |'
    write(*,*)' !!! Program could not be executed !!! |'
    write(*,*)'                                       |'
    write(*,*)'The following error was detected:      |'
    write(*,*)'                                       |'
    if(execution_error == 98)then
      write(*,*)'*) Bad input file                      |'
    else if(execution_error == 97)then
      write(*,*)'*) Bad trajectory file                 |'
    end if
    write(*,*)'                                       |'
    write(*,*)'Look above for "ERROR" messages to find|'
    write(*,*)'out what went wrong.                   |'
  end if
  write(*,*)'_______________________________________/'
!******************************************************
  end subroutine

end module
