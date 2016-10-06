module good_bye

  use read_input

  implicit none

  contains

  subroutine print_good_bye()
!******************************************************
! Print message at the end depending on execution
! success
   if(execution_error == 97)then
97   write(*,*)'                                       |'
     write(*,*)'Program successfully executed          |'
     write(*,*)'                                       |'
     write(*,*)'End of execution:                      |'
     call timestring(time)
     write(*,'(1X,A39,A)')time,'|'
   else if(execution_error == 98)then
98   write(*,*)'                                       |'
     write(*,*)'Program could not be executed          |'
     write(*,*)'Check for errors                       |'
   end if
99 write(*,*)'_______________________________________/'
!******************************************************
  end subroutine

end module
