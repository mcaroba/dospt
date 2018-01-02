!=================================================================================================
!=================================================================================================
subroutine sort_supergroups(ngroups, nsupergroups, ngroups_in_supergroup, group_in_supergroup)

  implicit none

  integer :: ngroups, nsupergroups, iostatus, temp_int, i, j
  integer :: group_start, group_end
  integer, allocatable :: ngroups_in_supergroup(:), group_in_supergroup(:,:)
  character*32 :: sg_string
  character*1 :: sg_char

  open(unit=10, file="supergroups", status="old")
  nsupergroups = 0
  do
    read(10,*, iostat=iostatus)
    if(iostatus /= 0)exit
    nsupergroups = nsupergroups + 1
  end do
  rewind(10)

  allocate( ngroups_in_supergroup(1:nsupergroups) )
  ngroups_in_supergroup = 0
  allocate( group_in_supergroup(1:nsupergroups, 1:ngroups) )

  do i=1,nsupergroups
    sg_string = ""
    group_start = 0
    do
      read(10,'(A)', advance="no", iostat=iostatus) sg_char
      if(sg_char == "-")then
        read(sg_string,*) group_start
        sg_string = ""
      else if( sg_char == " " .and. sg_string == "" )then
        continue
      else if( sg_char == "," .or. sg_char == " " .or. iostatus < 0 )then
        read(sg_string,*) group_end
        if( group_start == 0 )then
          group_start = group_end
        end if
        do j=group_start,group_end
          ngroups_in_supergroup(i) = ngroups_in_supergroup(i) + 1
          group_in_supergroup(i, ngroups_in_supergroup(i)) = j
        end do
        sg_string = ""
        group_start = 0
      else
        write(sg_string,*) trim(sg_string) // sg_char
      end if
      if( iostatus /= 0 )then
        exit
      end if
    end do
  end do
  rewind(10)
  close(10)

end subroutine
!=================================================================================================
!=================================================================================================
