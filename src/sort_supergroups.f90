!=================================================================================================
!=================================================================================================
subroutine sort_supergroups(ngroups, nsupergroups, ngroups_in_supergroup, group_in_supergroup)

  implicit none

  integer :: ngroups, nsupergroups, iostatus, temp_int, i, j, k
  integer :: group_start, group_end, nlines
  integer, allocatable :: ngroups_in_supergroup(:), group_in_supergroup(:,:), supergroup_in_this_line(:)
  character*1024 :: sg_string
  character*1 :: sg_char

! Initial parsing of supergroups file
  open(unit=10, file="supergroups", status="old")
  nsupergroups = 0
  nlines = 0
  do
    read(10, '(A)', iostat=iostatus) sg_string
    sg_char = adjustl(sg_string)
    if(iostatus /= 0)then
      exit
    else if( sg_char /= "!" .and. sg_char /= "#" .and. sg_char /= " " )then
      nsupergroups = nsupergroups + 1
    end if
    nlines = nlines + 1
  end do
  rewind(10)
  allocate( supergroup_in_this_line(1:nlines) )
  supergroup_in_this_line = 0
  nsupergroups = 0
  nlines = 0
  do
    read(10, '(A)', iostat=iostatus) sg_string
    sg_char = adjustl(sg_string)
    if(iostatus /= 0)then
      exit
    else if( sg_char /= "!" .and. sg_char /= "#" .and. sg_char /= " " )then
      nlines = nlines + 1
      nsupergroups = nsupergroups + 1
      supergroup_in_this_line(nlines) = nsupergroups
    else
      nlines = nlines + 1
    end if
  end do
  rewind(10)


  allocate( ngroups_in_supergroup(1:nsupergroups) )
  ngroups_in_supergroup = 0
  allocate( group_in_supergroup(1:nsupergroups, 1:ngroups) )

  do k = 1, nlines
    if( supergroup_in_this_line(k) == 0 )then
!     This is a comment or blank line
      read(10, *)
      cycle
    else
!     This is a line with supergroup info in it
      i = supergroup_in_this_line(k)
    end if
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
