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
  integer, allocatable :: natoms_in_group(:), atoms_in_group(:,:), atom_belongs_to_group(:)
  real*8, allocatable :: mass_types_values(:)
  character*16, allocatable :: mass_types(:)
  real*8, allocatable :: symmetry_number_group(:)
  logical :: check_topology = .true.

! Supergroup interface
  integer :: nsupergroups = 0
  integer, allocatable :: ngroups_in_supergroup(:), group_in_supergroup(:,:), group_belongs_to_supergroup(:)
  logical :: calc_supergroups

! Looping variables
  integer :: i, i2, j, j2, ios, k, k2

! Topology interface
  integer :: nbond_types
  character*16, allocatable :: bond_type(:,:)
  real*8, allocatable :: bond_cutoffs(:,:)
  integer, allocatable :: nspecies_in_topology(:), topology_in_supergroup(:), neach_species_in_topology(:,:)
  integer :: ntopology_types
  real*8, allocatable :: symmetry_number_of_topology(:)
  character*16, allocatable :: species_in_topology(:,:)

  contains

  subroutine print_welcome_and_read_input()
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
  write(*,*)'   Density of States (two-) Phase Thermodynamics, v0.1.1    |'
  write(*,*)'                   DoSPT v0.1.1-alpha                       |'
  write(*,*)'                    http://dospt.org                        |'
  write(*,*)'                          ...                               |'
  write(*,*)'                Written by Miguel A. Caro                   |'
  write(*,*)'                          ...                               |'
  write(*,*)'                    mcaroba@gmail.com                       |'
  write(*,*)'                          ...                               |'
  write(*,*)'                 Last updated Feb 2017                      |'
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

! The different input files (other than input) must be read in this order

! Read the masses file
call read_masses()

! Read the groups file
call read_groups()

! Read the supergroups file
call read_supergroups()

! Read the topology file
call read_topology()
end subroutine










subroutine read_groups()
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
    allocate( symmetry_number_group(1:ngroups) )
    allocate( atom_belongs_to_group(1:natoms) )
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
! This variable takes up the atom number and returns the group number (for the group that atom belongs to)
    do i=1, ngroups
      do j=1, natoms_in_group(i)
        atom_belongs_to_group(atoms_in_group(i,j)) = i
      end do
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
end subroutine






subroutine read_supergroups()
  interface
    subroutine sort_supergroups(ngroups, nsupergroups, ngroups_in_supergroup, group_in_supergroup)
      integer, intent(in) :: ngroups
      integer, intent(out) :: nsupergroups
      integer, allocatable, intent(out) :: ngroups_in_supergroup(:), group_in_supergroup(:,:)
    end subroutine
  end interface
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
  if( calc_supergroups )then
! This variable takes up the group number and returns the supergroup number
    allocate( group_belongs_to_supergroup(1:ngroups) )
    do i=1, nsupergroups
      do j=1, ngroups_in_supergroup(i)
        group_belongs_to_supergroup(group_in_supergroup(i,j)) = i
      end do
    end do
  end if
!******************************************************
end subroutine





subroutine read_topology()
!******************************************************
! Read in topology from file
  open(unit=10, file="topology", status="old", iostat=iostatus)
    write(*,*)'                                       |'
    write(*,*)'Checking topology file...              |'
    if(iostatus/=0)then
      close(10)
      write(*,*)'WARNING: topology file could not be    |'
      write(*,*)'found. You can safely disregard this   |'
      write(*,*)'message unless you expect bond breaking|'
      write(*,*)'in your simulation.                    |'
      write(*,*)'                                       |'
      write(*,*)'.......................................|'
      check_topology = .false.
    else
      nbond_types = 0
      do
        read(10,*,iostat=ios) keyword
        if (ios/=0) exit
        if( keyword == "bond" )then
          nbond_types = nbond_types + 1
        end if
      end do
      allocate( bond_type(1:nbond_types, 1:2) )
      allocate( bond_cutoffs(1:nbond_types, 1:2) )
      nbond_types = 0
      rewind(10)
      do
        read(10,*,iostat=ios) keyword
        if (ios/=0) exit
        if( keyword == "bond" )then
          nbond_types = nbond_types + 1
          backspace(10)
! WARNING: YOU SHOULD FIX THIS
! The code should assert that bond_cutoffs(:,2) >= bond_cutoffs(:,1), that is, the bond creation
! should never happen at a larger cutoff than the bond breaking, otherwise the
! code will go crazy
          read(10,*) keyword, bond_type(nbond_types,1:2), bond_cutoffs(nbond_types,1:2)
        end if
      end do
!     Check new group topologies and asignment to supergroups
      rewind(10)
      ntopology_types = 0
      do
        read(10,*,iostat=ios) keyword
        if (ios/=0) exit
        if( keyword == "group" )then
          ntopology_types = ntopology_types + 1
        end if
      end do
      allocate( nspecies_in_topology(1:ntopology_types) )
      allocate( symmetry_number_of_topology(1:ntopology_types) )
      allocate( topology_in_supergroup(1:ntopology_types) )
      rewind(10)
!     k is the largest number of different species in a group; for memory allocation purposes
      k = 0
      j2 = 0
      do
        read(10,*,iostat=ios) keyword
        if (ios/=0) exit
        if( keyword == "group" )then
          j2 = j2 + 1
          backspace(10)
          j = 0
          k2 = 0
          do while( keyword /= "sym" )
            j = j + 1
            read(10, *) ( keyword, i = 1, j + 1)
            backspace(10)
            if( keyword /= "sym" )then
              do i2 = 1, nmasses
                if( keyword == mass_types(i2) )then
                  k2 = k2 + 1
                end if
              end do
            end if
          end do
          nspecies_in_topology(j2) = k2
          if( k2 > k )then
            k = k2
          end if
          read(10,*)
        end if
      end do
      allocate( species_in_topology(1:ntopology_types, 1:k) )
      allocate( neach_species_in_topology(1:ntopology_types, 1:k) )
      rewind(10)
      i = 0
      do
        read(10,*,iostat=ios) keyword
        if (ios/=0) exit
        if( keyword == "group" )then
          i = i + 1
          backspace(10)
          do j = 1, nspecies_in_topology(i)
            read(10,*) keyword, (keyword, k = 1, 2*(j-1)), species_in_topology(i, j), neach_species_in_topology(i, j)
            backspace(10)
          end do
          read(10,*) (keyword, k = 1, 2*nspecies_in_topology(i) + 1), keyword, symmetry_number_of_topology(i), &
                     keyword, topology_in_supergroup(i)
        end if
      end do
!     Increase number of supergroups
      do i = 1, ntopology_types
        if( topology_in_supergroup(i) > nsupergroups )then
          nsupergroups = topology_in_supergroup(i)
        end if
      end do
    end if
  close(10)
  do i = 1, nbond_types
    if( i == 1 )then
      write(*,*)'                                       |'
      write(*,*)'The following bond types were found in |'
      write(*,*)'the topology file:                     |'
      write(*,*)'                                       |'
    end if
      write(*,'(1X,A4,A,A4,A,F5.3,A,F5.3,A)') trim(bond_type(i,1)), ' <-> ', adjustl(bond_type(i,2)), &
           ' Cuts.: (', bond_cutoffs(i,1), ', ', bond_cutoffs(i,2), ') nm |'
  end do
  write(*,*)'                                       |'
  write(*,*)'.......................................|'
!******************************************************
  end subroutine


  


  subroutine read_masses()
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
