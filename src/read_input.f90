module read_input

  use constants

  implicit none

! Variables for input parameters (read from "input" file):
! L(1:3) = lattice vectors (only orthrhombic cells currently supported)
! tau = trajectory length in ps
! n = number of snapshots in the trajectory
! T = temperature
! smooth = whether to filter DoS for individual groups
! sigma_nu = characteristic filtering length in ps^-1
  character*64 :: keyword
  character*1 :: keyword_first
  character*40 :: time
  character*3 :: mode = "gro"
  character :: crap
  real*8 :: L(1:3) = 0., V = 0., tau = 0., T = 1.d-10, twobykT, sigma_nu = 0.d0
  logical :: smooth = .false., estimate_vel = .false., renormalize = .false., vibrational_gas = .false.
  logical :: f_opt = .true., f_rot_opt = .true., write_dos = .true.
  character*3 :: smixture = "mol", hs_formalism
  integer :: natoms, n = 0, iostatus, execution_error

! Variables for grouping interface
  integer :: ngroups, nmasses
  integer, allocatable :: natoms_in_group(:), atoms_in_group(:,:)
  real*8, allocatable :: mass_types_values(:), mass_group(:), Sgroup(:,:,:)
  character*16, allocatable :: mass_types(:)
  real*8, allocatable :: symmetry_number_group(:)

! Supergroup interface
  integer :: nsupergroups = 0
  integer, allocatable :: ngroups_in_supergroup(:), group_in_supergroup(:,:)
  logical :: calc_supergroups

! Looping variables
  integer :: i, i2, j, j2, ios


  contains

  subroutine print_welcome_and_read_input()

  interface
    subroutine sort_supergroups(ngroups, nsupergroups, ngroups_in_supergroup, group_in_supergroup)
      integer, intent(in) :: ngroups
      integer, intent(out) :: nsupergroups
      integer, allocatable, intent(out) :: ngroups_in_supergroup(:), group_in_supergroup(:,:)
    end subroutine
  end interface

!******************************************************
! Print welcome message, starting time of execution and
! read in input options and parameters
  write(*,*)' '
  write(*,*)'____________________________________________________________'
  write(*,*)'                                                            |'
  write(*,*)'  :MMMMMMMMO          .    MMMMMMMM .. MMMMMMMM$.MMMMMMMMMMM|'
  write(*,*)' .MMM ...MMMM..   +MN=   .MMM..........MMM...MMM ....MMM....|'
  write(*,*)'..MMM...  MMM...MMMMMMMM  MMMMN.......MMM .. MMM.....MMM....|'
  write(*,*)'..MMM...  MMM..MMM.  MMM. .MMMMMMM ...MMMMMMMMM.....MMM ....|'
  write(*,*)'.MMM ... MMM,.MMM.  .MMM . .  .ZMMM.. MMMMMO,.......MMM.....|'
  write(*,*)'.MMM . MMMM=..MMM ..MMM7 ..  ..8MMM..MMM ...........MMM.....|'
  write(*,*)'MMMMMMMMMM. .. MMMMMMM   MMMMMMMMD...MMM ..........MMM=.... |'
  write(*,*)'                                                            |'
  write(*,*)'                                                            |'
  write(*,*)'                    You are using the                       |'
  write(*,*)'     Density of States (two-) Phase Thermodynamics, v0.1    |'
  write(*,*)'                     DoSPT v0.1-alpha                       |'
  write(*,*)'                          ...                               |'
  write(*,*)'                Written by Miguel A. Caro                   |'
  write(*,*)'                          ...                               |'
  write(*,*)'                    mcaroba@gmail.com                       |'
  write(*,*)'                          ...                               |'
  write(*,*)'                 Last updated Sep 2016                      |'
  write(*,*)'                                                            |'
  write(*,*)'....................................... ____________________/'

! Prints start of execution time:
  write(*,*)'                                       |'
  write(*,*)'Start of execution:                    |'
  call timestring(time)
  write(*,'(1X,A39,A)')time,'|'
  write(*,*)'                                       |'
  write(*,*)'.......................................|'

! Read input file:
  open(unit=10,file='input',status='old',iostat=iostatus)
  write(*,*)'                                       |'
  write(*,*)'Checking input file...                 |'
  if(iostatus/=0)then
    close(10)
    write(*,*)'ERROR: input file could not be found   |'
    write(*,*)'                                       |'
    write(*,*)'.......................................|'
    execution_error = 98
    return
  else
    do while(iostatus==0)
      read(10,*,iostat=iostatus)keyword
      keyword = trim(keyword)
      if(iostatus/=0)then
        exit
      end if
      keyword_first=keyword
      if(keyword_first=='#')then
        continue
      else if(keyword=='points')then
        backspace(10)
        read(10,*,iostat=iostatus) crap, crap, n
      else if(keyword=='cell')then
        backspace(10)
        read(10,*,iostat=iostatus) crap, crap, L(1), L(2), L(3)
      else if(keyword=='tau')then
        backspace(10)
        read(10,*,iostat=iostatus) crap, crap, tau
      else if(keyword=='temperature')then
        backspace(10)
        read(10,*,iostat=iostatus) crap, crap, T
      else if(keyword=='smooth')then
        backspace(10)
        read(10,*,iostat=iostatus) crap, crap, smooth
      else if(keyword=='sigma_nu')then
        backspace(10)
        read(10,*,iostat=iostatus) crap, crap, sigma_nu
      else if(keyword=='format')then
        backspace(10)
        read(10,*,iostat=iostatus) crap, crap, mode
      else if(keyword=='estimate_velocities')then
        backspace(10)
        read(10,*,iostat=iostatus) crap, crap, estimate_vel
      else if(keyword=='renormalize_dos')then
        backspace(10)
        read(10,*,iostat=iostatus) crap, crap, renormalize
      else if(keyword=='vibrational_gas')then
        backspace(10)
        read(10,*,iostat=iostatus) crap, crap, vibrational_gas
      else if(keyword=='entropy_mixture')then
        backspace(10)
        read(10,*,iostat=iostatus) crap, crap, smixture
      else if(keyword=='f_opt')then
        backspace(10)
        read(10,*,iostat=iostatus) crap, crap, f_opt
      else if(keyword=='f_rot_opt')then
        backspace(10)
        read(10,*,iostat=iostatus) crap, crap, f_rot_opt
      else if(keyword=='write_dos')then
        backspace(10)
        read(10,*,iostat=iostatus) crap, crap, write_dos
      else if(keyword=='hs_formalism')then
!       Define the HS formalism, overriding f_opt and f_rot_opt flags
        backspace(10)
        read(10,*,iostat=iostatus) crap, crap, hs_formalism
        if(hs_formalism == "ome" .or. hs_formalism == "lin")then
          f_rot_opt = .false.
          f_opt = .false.
        else if(hs_formalism == "scs")then
          f_rot_opt = .true.
          f_opt = .true.
        else
          iostatus=1
          write(*,*)'ERROR - invalid value for keyword: ', keyword
        end if
      else
        iostatus=1
        write(*,*)'ERROR - invalid keyword: ', keyword
      end if
    if(iostatus/=0)then
      write(*,*)'ERROR: bad input file                  |'
      write(*,*)'                                       |'
      write(*,*)'.......................................|'
      execution_error = 98
      return
    end if
    end do
  end if
  close(10)
  write(*,*)'                                       |'
  write(*,'(1X,A)')'You specified the following options:   |'
  write(*,*)'                                       |'
  write(*,'(1X,A,I15,A)')'No. of snapshots = ', n, '     |'
  write(*,'(1X,A,F14.3,A)')'Trajectory period = ', tau, ' ps  |'
  write(*,'(1X,A,F7.3,A,F7.3,A,F7.3,A)')'Cell = [', L(1), ', ',  L(2), ', ', L(3), '] nm  |'
  write(*,'(1X,A,F20.3,A)')'Temperature = ', T, ' K   |'
  if(smooth)then
    if(sigma_nu == 0.)then
      sigma_nu = 5.d0/tau
     end if
    write(*,'(1X,A,F15.3,A)')'Filtering length = ', sigma_nu, ' 1/ps|'
  end if
  write(*,*)'                                       |'
  write(*,*)'.......................................|'
  twobykT = 2.d0 / kB / T
  V = L(1) * L(2) * L(3)
!******************************************************






!******************************************************
! Read in group information
  open(unit=10, file="groups", status="old", iostat=iostatus)
    write(*,*)'                                       |'
    write(*,*)'Checking groups file...                |'
    if(iostatus/=0)then
      close(10)
      write(*,*)'ERROR: groups file could not be found  |'
      write(*,*)'                                       |'
      write(*,*)'.......................................|'
      execution_error = 98
      return
    end if
! Number of groups in the system
    read(10,*) natoms, ngroups
    allocate( natoms_in_group(1:ngroups) )
    allocate( mass_group(1:ngroups) )
    allocate( symmetry_number_group(1:ngroups) )
    allocate( Sgroup(1:ngroups,1:(n+1)/2,1:3) )
    do i=1,ngroups
      read(10,*) natoms_in_group(i), symmetry_number_group(i)
      read(10,*)
    end do
! Identify largest group for memory allocation purposes
    j=0
    do i=1,ngroups
      if( natoms_in_group(i) > j)then
        j = natoms_in_group(i)
      end if
    end do
    allocate( atoms_in_group(1:ngroups, 1:j) )
! Read in labels for atoms in each group
    rewind(10)
    read(10,*)
    do i=1,ngroups
      read(10,*)
      read(10,*) ( atoms_in_group(i,j), j = 1, natoms_in_group(i) )
    end do
! Check if same atom shows in more than one group, or more than
! once in the same group
    do i=1,ngroups
      do j=1,natoms_in_group(i)
        do i2=1,ngroups
          do j2=1,natoms_in_group(i2)
!           Same atom more than once in a group
            if(i == i2 .and. j /= j2 .and. atoms_in_group(i,j) == atoms_in_group(i2,j2))then
              write(*,*)       '                                       |'
              write(*,'(1X,A)')'ERROR:                                 |'
              write(*,*)       '                                       |'
              write(*,'(1X,A,I8,A)')'Atom ', atoms_in_group(i,j), ' is present more than once|'
              write(*,'(1X,A,I8,A)')'in group ', i, '. You cannot do that! |'
              write(*,*)'                                       |'
              write(*,*)'.......................................|'
              execution_error = 98
              return
            end if
!           Same atom in different groups
            if(i /= i2 .and. atoms_in_group(i,j) == atoms_in_group(i2,j2))then
              write(*,*)       '                                       |'
              write(*,'(1X,A)')'ERROR:                                 |'
              write(*,*)       '                                       |'
              write(*,'(1X,A,I8,A)')'Atom ', atoms_in_group(i,j), ' is present in more than  |'
              write(*,'(1X,A)')'one group. Not currently implemented.  |'
              write(*,*)'                                       |'
              write(*,*)'.......................................|'
              execution_error = 98
              return
            end if
          end do
        end do
      end do
    end do
  close(10)
  write(*,*)       '                                       |'
  write(*,'(1X,A)')'Group information:                     |'
  write(*,*)       '                                       |'
  write(*,'(1X,A,I19,A)')'No. of atoms = ', natoms, '     |'
  write(*,'(1X,A,I18,A)')'No. of groups = ', ngroups, '     |'
  write(*,*)'                                       |'
  write(*,*)'.......................................|'
!******************************************************






!******************************************************
! Read in supergroup information
  open(unit=10, file="supergroups", status="old", iostat=iostatus)
  if(iostatus == 0)then
    close(10)
!   If file exits calculate supergroup stuff
    calc_supergroups = .true.
    write(*,*)'                                       |'
    write(*,*)'Checking supergroups file...           |'
    call sort_supergroups(ngroups, nsupergroups, ngroups_in_supergroup, group_in_supergroup)
    close(10)
    write(*,*)       '                                       |'
    write(*,'(1X,A)')'Supergroup information:                |'
    write(*,*)       '                                       |'
    write(*,'(1X,A,I13,A)')'No. of supergroups = ', nsupergroups, '     |'
    write(*,*)'                                       |'
    write(*,*)'.......................................|'
  else
    close(10)
  end if
!******************************************************






!******************************************************
! Read in masses from user's library
  open(unit=10, file="masses", status="old", iostat=iostatus)
    write(*,*)'                                       |'
    write(*,*)'Checking masses file...                |'
    if(iostatus/=0)then
      close(10)
      write(*,*)'ERROR: groups file could not be found  |'
      write(*,*)'                                       |'
      write(*,*)'.......................................|'
      execution_error = 98
      return
    end if
    nmasses=0
    do
      read(10,*,iostat=ios)
      if (ios/=0) exit     
      nmasses=nmasses+1
    end do
    allocate( mass_types(1:nmasses) )
    allocate( mass_types_values(1:nmasses) )
    write(*,*)'                                       |'
    write(*,*)'Masses were found for:                 |'
    write(*,*)'                                       |'
    write(*,*)'    Element (tag);      Mass (amu)     |'
    rewind(10)
    do i=1,nmasses
      read(10,*) mass_types(i), mass_types_values(i)
      write(*,'(1X,A17,F17.2,A)') trim(mass_types(i)), mass_types_values(i),'     |'
    end do
  close(10)
  write(*,*)'                                       |'
  write(*,*)'.......................................|'
!******************************************************
  end subroutine

  end module
