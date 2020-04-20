module topology

  use read_input
  use trajectory

  implicit none

  integer, allocatable :: birth_time(:), death_time(:), lifetime(:)
  logical :: error = .false., topology_has_changed = .false., interpolate_dos = .false.
  real*8, allocatable :: weight_group(:), weight_supergroup(:)


  contains

!=================================================================================================
!=================================================================================================
subroutine sort_out_topology()
!******************************************************
! Rebuild topology in case bonds break during the dynamics
  if( check_topology )then
    call rebuild_topology(n, nrebuild_top, natoms, ngroups, nsupergroups, L, positions, species, &
                          atoms_in_group, natoms_in_group, symmetry_number_group, atom_belongs_to_group, &
                          group_in_supergroup, ngroups_in_supergroup, group_belongs_to_supergroup, &
                          nspecies_in_topology, topology_in_supergroup, neach_species_in_topology, &
                          ntopology_types, symmetry_number_of_topology, species_in_topology, &
                          nbond_types, bond_type, bond_cutoffs, birth_time, death_time, topology_has_changed, &
                          nconstraints_in_group, nconstraints_of_topology )
  else
    allocate( birth_time(1:ngroups) )
    allocate( death_time(1:ngroups) )
    birth_time = 1
    death_time = n + 1
  end if
! Make sure that we don't add unnecessary supergroups unless they are required
  nsupergroups = size(ngroups_in_supergroup,1)
  allocate( lifetime(1:ngroups) )
  allocate( weight_group(1:ngroups) )
  allocate( weight_supergroup(1:nsupergroups) )
  weight_supergroup = 0.d0
  lifetime = death_time - birth_time
  do i = 1, ngroups
    weight_group(i) = dfloat(lifetime(i)) / dfloat(n)
  end do
  do i = 1, nsupergroups
    do j = 1, ngroups_in_supergroup(i)
      weight_supergroup(i) = weight_supergroup(i) + weight_group(group_in_supergroup(i,j))
    end do
  end do
! If the topology has not been rebuilt, then we return the nsupergroups variable to its original
! state
  if( .not. topology_has_changed )then
    nsupergroups = size(ngroups_in_supergroup, 1)
  end if
!******************************************************
end subroutine
!=================================================================================================
!=================================================================================================






!=================================================================================================
!=================================================================================================
subroutine rebuild_topology(nsteps, dstep, natoms, ngroups, nsupergroups, L_cell, positions, species, &
           atoms_in_group, natoms_in_group, symmetry_number_group, atom_belongs_to_group, &
           group_in_supergroup, ngroups_in_supergroup, group_belongs_to_supergroup, &
           nspecies_in_topology, topology_in_supergroup, neach_species_in_topology, &
           ntopology_types, symmetry_number_of_topology, species_in_topology, &
           nbond_types, bond_type, bond_cutoffs, birth_time, death_time, topology_has_changed, &
           nconstraints_in_group, nconstraints_of_topology )

  implicit none

! In and out variables
  integer :: nsteps, natoms, ngroups, nbond_types, nsupergroups, dstep
  real*8 :: L_cell(1:3), bond_cutoffs(:,:), update_bar
  real*8, allocatable :: symmetry_number_group(:)
  real*8, allocatable :: nconstraints_in_group(:)
  real*4 :: positions(:,:,:)
  character*16 :: species(:), bond_type(:,:)
  integer, allocatable :: atoms_in_group(:,:), natoms_in_group(:)
  integer, allocatable :: atom_belongs_to_group(:), death_time(:), birth_time(:)
  integer, allocatable :: ngroups_in_supergroup(:), group_in_supergroup(:,:)
  integer :: nspecies_in_topology(:), topology_in_supergroup(:), neach_species_in_topology(:,:)
  integer :: ntopology_types
  real*8 :: symmetry_number_of_topology(:)
  real*8 :: nconstraints_of_topology(:)
  character*16 :: species_in_topology(:,:)
  logical :: topology_has_changed

! Internal variables
  integer :: i, j, k, iostatus, i2, j2, l, i3, step, ios, k2
  integer :: new_groups
  logical, allocatable :: group_still_exists(:), group_still_exists_temp(:)
  logical :: rebuild_groups = .false., broken_message = .false., repeat_broken_message = .true.
  integer, allocatable :: atoms_in_group_temp(:,:), natoms_in_group_temp(:)
  integer, allocatable :: atom_belongs_to_new_group(:), natoms_in_new_group(:)
  real*8 :: d
  real*8, allocatable :: symmetry_number_group_temp(:)
  real*8, allocatable :: nconstraints_in_group_temp(:)
  character*32 :: temp_char
  integer, allocatable :: neach_species_in_topology_temp(:,:), group_belongs_to_supergroup(:), ngroups_in_supergroup_temp(:), &
                          group_in_supergroup_temp(:,:), group_belongs_to_supergroup_temp(:), birth_time_temp(:), &
                          death_time_temp(:)
  logical :: topology_match, broken_bond, atom_reallocated, print_bar = .false., match
  logical, allocatable :: atom_visited(:), atoms_bonded(:,:), atoms_bonded_prev(:,:)


! There are no new groups to begin with
  new_groups = 0
! List of existing groups
  allocate( group_still_exists(1:ngroups) )
  group_still_exists = .true.
  allocate( birth_time(1:ngroups) )
  birth_time = 1
  allocate( death_time(1:ngroups) )
  death_time = nsteps + 1


! Initialize these arrays
  allocate( atoms_in_group_temp(1:1,1:1) )
  allocate( natoms_in_group_temp(1:1) )
  allocate( symmetry_number_group_temp(1:1) )
  allocate( nconstraints_in_group_temp(1:1) )
  nconstraints_in_group_temp = 0.d0
  allocate( group_still_exists_temp(1:1) )
  allocate( birth_time_temp(1:1) )
  allocate( death_time_temp(1:1) )


! Bonding search arrays initialization
  allocate( atom_visited(1:natoms) )
  atom_visited = .false.
  allocate( atoms_bonded(1:natoms, 1:natoms) )
  atoms_bonded = .false.
  allocate( atoms_bonded_prev(1:natoms, 1:natoms) )
  atoms_bonded_prev = .false.





  call system("rm -rf topol; mkdir -p topol")







! Check initial consistency between groups and topology files
  step = 1
  call build_bond_list(natoms, atom_belongs_to_group, species, bond_type, bond_cutoffs, positions, &
                       L_cell, step, atoms_bonded_prev, atoms_bonded, .true.)
  allocate( atom_belongs_to_new_group(1:natoms) )
  do i = 1, natoms
!   Find all the atoms which can be connected to i
    atom_belongs_to_new_group = 0
    atom_visited = .false.
    atom_belongs_to_new_group(i) = 1
    atom_visited(i) = .true.
    do j = 1, natoms
      call find_neighbors(i, j, natoms, atoms_bonded, atom_visited, atom_belongs_to_new_group, 1)
    end do
!   Check that the new group built from the topology has the same number of atoms as the original one from the groups file
    match = .true.
    k = 0
    do j = 1, natoms
      if( atom_belongs_to_new_group(j) == 1 )then
        k = k + 1
      end if
    end do
    if( k /= natoms_in_group(atom_belongs_to_group(i)) )then
      match = .false.
    else
!     If the number of atoms coincides, make sure that the atoms in i's group as given in the groups file coincide with those
!     obtain from connecting bonded atoms
      loop8: do k = 1, natoms_in_group(atom_belongs_to_group(i))
        i2 = atoms_in_group(atom_belongs_to_group(i), k)
        match = .false.
        loop7: do j2 = 1, natoms
          if( atom_belongs_to_new_group(j2) == 1 )then
            if( i2 == j2 )then
              match = .true.
              exit loop7
            end if
          end if
        end do loop7
!       If we couldn't find a match get out
        if( .not. match )then
          exit loop8
        end if
      end do loop8
    end if
    if( .not. match )then
      k = atom_belongs_to_group(i)
      topology_has_changed = .true.
      group_still_exists(k) = .false.
      rebuild_groups = .true.
      broken_message = .true.
      death_time(k) = step
    end if
  end do
  deallocate( atom_belongs_to_new_group )

  if( broken_message .and. step == 1 )then
    write(*,*)'                                       |'
    write(*,*)'WARNING: your initial groups file is   |'
    write(*,*)'not consistent with the specified      |'
    write(*,*)'topology. I am rebuilding your groups. |'
    broken_message = .false.
  end if





  write(*,*)'                                       |'
  write(*,*)'Checking for bond breaking according to|'
  write(*,*)'specified topology...                  |'
  write(*,*)'                                       |'
  write(*,*)'Progress:                              |'
  write(*,*)'                                       |'
  write(*,'(1X,A)',advance='no')'['
  update_bar = dfloat(nsteps)/36.d0
  k2 = 0



  do step = 1, nsteps, dstep

!   Update progress bar every nsteps/36 iterations
    if(dfloat(step) > dfloat(k2)*update_bar)then
      write(*,'(A)', advance='no')'='
      k2 = k2 + 1
    end if

!   Compare list of bonded atoms now and before
    atoms_bonded_prev = atoms_bonded
    if( step > 1 )then
      call build_bond_list(natoms, atom_belongs_to_group, species, bond_type, bond_cutoffs, positions, &
                           L_cell, step, atoms_bonded_prev, atoms_bonded, .false.)
      do i = 1, natoms
        do j = i+1, natoms
!         If bonding info has changed, kill groups involved
          if( atoms_bonded(i,j) .neqv. atoms_bonded_prev(i,j) )then
            topology_has_changed = .true.
            rebuild_groups = .true.
            broken_message = .true.
            k = atom_belongs_to_group(i)
            group_still_exists(k) = .false.
            death_time(k) = step
            k = atom_belongs_to_group(j)
            group_still_exists(k) = .false.
            death_time(k) = step
          end if
        end do
      end do
    end if

!   Now we create the new groups according to the new bonding info
    if( rebuild_groups )then
      new_groups = 0
      allocate( atom_belongs_to_new_group(1:natoms) )
      atom_belongs_to_new_group = 0
      atom_visited = .false.
      do i = 1, natoms
        if( .not. atom_visited(i) .and. .not. group_still_exists(atom_belongs_to_group(i)) .and. &
            death_time(atom_belongs_to_group(i)) == step )then
          new_groups = new_groups + 1
          atom_belongs_to_new_group(i) = new_groups
          atom_visited(i) = .true.
          do j = 1, natoms
            call find_neighbors(i, j, natoms, atoms_bonded, atom_visited, atom_belongs_to_new_group, new_groups)
          end do
        end if
      end do

!     Check how many atoms are in each of the new groups
      allocate( natoms_in_new_group(1:new_groups) )
      do i = 1, new_groups
        k = 0
        do j = 1, natoms
          if( atom_belongs_to_new_group(j) == i )then
            k = k + 1
!           Also make sure all of these groups are properly terminated
            if( group_still_exists(atom_belongs_to_group(j)) )then
              group_still_exists(atom_belongs_to_group(j)) = .false.
              death_time(atom_belongs_to_group(j)) = step
            end if
          end if
        end do
        natoms_in_new_group(i) = k
      end do


!     Now we need to allocate some arrays for the new groups. The arrays that
!     depend on the number of groups need to be extended by + new_groups

!     Identify largest group for memory allocation purposes
      j2=0
      do i=1,ngroups
        if( natoms_in_group(i) > j2)then
          j2 = natoms_in_group(i)
        end if
      end do
      do i=1,new_groups
        if( natoms_in_new_group(i) > j2)then
          j2 = natoms_in_new_group(i)
        end if
      end do
      deallocate( atoms_in_group_temp )
      deallocate( natoms_in_group_temp )
      deallocate( symmetry_number_group_temp )
      deallocate( nconstraints_in_group_temp )
      deallocate( group_still_exists_temp )
      deallocate( birth_time_temp )
      deallocate( death_time_temp )

      allocate( atoms_in_group_temp(1:ngroups+new_groups,1:j2) )
!     Allocate the other arrays
      allocate( natoms_in_group_temp(1:ngroups+new_groups) )
      allocate( symmetry_number_group_temp(1:ngroups+new_groups) )
!     If the symmetry number is not provided we assume sigma = 1
      symmetry_number_group_temp(1:ngroups+new_groups) = 1.d0
      allocate( nconstraints_in_group_temp(1:ngroups+new_groups) )
!     If the number of constraints is not provided we set it to zero
      nconstraints_in_group_temp(1:ngroups+new_groups) = 0.d0
      allocate( group_still_exists_temp(1:ngroups+new_groups) )
      allocate( birth_time_temp(1:ngroups+new_groups) )
      allocate( death_time_temp(1:ngroups+new_groups) )

!     With all the temp arrays allocated, now we copy the current arrays before adding the new info
      i = size(atoms_in_group, 2)
      atoms_in_group_temp(1:ngroups,1:i) = atoms_in_group(1:ngroups,1:i)
      natoms_in_group_temp(1:ngroups) = natoms_in_group(1:ngroups)
      symmetry_number_group_temp(1:ngroups) = symmetry_number_group(1:ngroups)
      nconstraints_in_group_temp(1:ngroups) = nconstraints_in_group(1:ngroups)
      group_still_exists_temp(1:ngroups) = group_still_exists(1:ngroups)
      death_time_temp(1:ngroups) = death_time(1:ngroups)
      birth_time_temp(1:ngroups) = birth_time(1:ngroups)


!     Now we add the new info
      do i = 1, natoms
        if( atom_belongs_to_new_group(i) /= 0 )then
          atom_belongs_to_group(i) = ngroups + atom_belongs_to_new_group(i)
        end if
      end do
      group_still_exists_temp(ngroups+1 : ngroups+new_groups) = .true.
      do k = 1, new_groups
        natoms_in_group_temp(ngroups + k) = natoms_in_new_group(k)
        death_time_temp(ngroups + k) = nsteps + 1
        birth_time_temp(ngroups + k) = step
        j = 0
        do i = 1, natoms
          if( atom_belongs_to_new_group(i) == k )then
            j = j + 1
            atoms_in_group_temp(ngroups + k, j) = i
          end if
        end do
      end do

!     Reallocate arrays with new number of groups
      deallocate( atoms_in_group )
      deallocate( natoms_in_group )
      deallocate( symmetry_number_group )
      deallocate( nconstraints_in_group )
      deallocate( group_still_exists )
      deallocate( death_time )
      deallocate( birth_time )
      ngroups = ngroups + new_groups
      allocate( atoms_in_group(1:ngroups,1:j2) )
      allocate( natoms_in_group(1:ngroups) )
      allocate( symmetry_number_group(1:ngroups) )
      allocate( nconstraints_in_group(1:ngroups) )
      allocate( group_still_exists(1:ngroups) )
      allocate( death_time(1:ngroups) )
      allocate( birth_time(1:ngroups) )
      atoms_in_group = atoms_in_group_temp
      natoms_in_group = natoms_in_group_temp
      symmetry_number_group = symmetry_number_group_temp
      nconstraints_in_group = nconstraints_in_group_temp
      group_still_exists = group_still_exists_temp
      death_time = death_time_temp
      birth_time = birth_time_temp







!     Reallocate supergroup variables
      allocate( group_belongs_to_supergroup_temp(1:ngroups) )
      group_belongs_to_supergroup_temp(1:ngroups-new_groups) = group_belongs_to_supergroup(1:ngroups-new_groups)
      deallocate( group_belongs_to_supergroup )
      allocate( group_belongs_to_supergroup(1:ngroups) )
      group_belongs_to_supergroup = group_belongs_to_supergroup_temp
      deallocate( group_belongs_to_supergroup_temp )
!
      i = size(ngroups_in_supergroup, 1)
      allocate( ngroups_in_supergroup_temp(1:nsupergroups) )
      ngroups_in_supergroup_temp = 0
      ngroups_in_supergroup_temp(1:i) = ngroups_in_supergroup(1:i)
      deallocate( ngroups_in_supergroup )
      allocate( ngroups_in_supergroup(1:nsupergroups) )
      ngroups_in_supergroup = ngroups_in_supergroup_temp
      deallocate( ngroups_in_supergroup_temp )
!
      i = size(group_in_supergroup, 1)
      j = size(group_in_supergroup, 2)
      allocate( group_in_supergroup_temp(1:nsupergroups, 1:ngroups) )
      group_in_supergroup_temp(1:i, 1:j) = group_in_supergroup(1:i, 1:j)
      deallocate( group_in_supergroup )
      allocate( group_in_supergroup(1:nsupergroups, 1:ngroups) )
      group_in_supergroup = group_in_supergroup_temp
      deallocate( group_in_supergroup_temp )
!
!     We need to identify the new groups with any of the given topologies
      allocate( neach_species_in_topology_temp(1:ntopology_types, 1:l) )
      do i = ngroups-new_groups+1, ngroups
        neach_species_in_topology_temp = neach_species_in_topology
        do k = 1, ntopology_types
          topology_match = .true.
!         Find this atom in this topology:
          loop3: do j = 1, natoms_in_group(i)
            j2 = atoms_in_group(i,j)
            do i2 = 1, nspecies_in_topology(k)
!             Decrease count in temp variable every time we find one match
              if( species(j2) == species_in_topology(k,i2) )then
                neach_species_in_topology_temp(k, i2) = neach_species_in_topology_temp(k, i2) - 1
                cycle loop3
              end if
            end do
          end do loop3
!         To match a new group with a topology, we must ensure that all counters are exactly zero
          do i2 = 1, nspecies_in_topology(k)
            if( neach_species_in_topology_temp(k, i2) /= 0 )then
              topology_match = .false.
            end if
          end do
          if( topology_match )then
            group_belongs_to_supergroup(i) = topology_in_supergroup(k)
            symmetry_number_group(i) = symmetry_number_of_topology(k)
            nconstraints_in_group(i) = nconstraints_of_topology(k)
            ngroups_in_supergroup(group_belongs_to_supergroup(i)) = ngroups_in_supergroup(group_belongs_to_supergroup(i)) + 1
            group_in_supergroup(group_belongs_to_supergroup(i), ngroups_in_supergroup(group_belongs_to_supergroup(i))) = i
            exit
          end if
        end do
        if( .not. topology_match )then
! FIX THIS: ADD WARNING MESSAGE IF GROUP CAN'T BE ASSIGNED TO A TOPOLOGY
          continue
        end if
      end do
      deallocate( neach_species_in_topology_temp )


!     This writes the list of groups at each point in the trajectory when there is rebuilding of topology
      write(temp_char,*) step
      open(unit=10, file="topol/groups." // trim(adjustl(temp_char)), status="unknown")
      write(10,*) natoms, ngroups
      do i=1,ngroups
        write(10,*) natoms_in_group(i), symmetry_number_group(i), group_still_exists(i), nconstraints_in_group(i)
        write(10,*) ( atoms_in_group(i,j), j = 1, natoms_in_group(i) )
      end do
      close(10)

!     Reinitialize/deallocate some variables
      rebuild_groups = .false.
      deallocate( atom_belongs_to_new_group )
      deallocate( natoms_in_new_group )
    end if

  end do

  do i = k2, 34
    write(*,'(A1)',advance='no')'='
  end do
  write(*,'(A3)')'] |'

  if( broken_message )then
    write(*,*)'                                       |'
    write(*,*)'WARNING: I have detected bonds breaking|'
    write(*,*)'during your simulation. Check out the  |'
    write(*,*)'topol directory to find out at which   |'
    write(*,*)'time step(s) it is happening.          |'
    write(*,*)'                                       |'
    write(*,*)'If you have a complicated system (e.g.,|'
    write(*,*)'highly mobile protons) this can get    |'
    write(*,*)'messy.                                 |'
  end if

  write(*,*)'                                       |'
  write(*,*)'.......................................|'

! COMMENT THIS PRINTING OUT IN PUBLIC VERSION OF CODE
!do i = 1,ngroups
!write(*,*) i, group_belongs_to_supergroup(i), ngroups_in_supergroup(group_belongs_to_supergroup(i)), &
!symmetry_number_group(i), birth_time(i), death_time(i)
!end do
!do i = 1, nsupergroups
!do j = 1, ngroups_in_supergroup(i)
!write(*,*) i, j, group_in_supergroup(i,j), ngroups_in_supergroup(i)
!end do
!end do

end subroutine












subroutine build_bond_list(natoms, atom_belongs_to_group, species, bond_type, bond_cutoffs, positions, &
                           L_cell, step, atoms_bonded_prev, atoms_bonded, initialize)

  implicit none

  logical :: atoms_bonded(:,:), atoms_bonded_prev(:,:), initialize
  integer :: natoms, i2, j2, l, step, nbond_types, atom_belongs_to_group(:)
  character*16 :: species(:), bond_type(:,:)
  real*4 :: positions(:,:,:)
  real*8 :: L_cell(:), d, bond_cutoffs(:,:)

  nbond_types = size(bond_type, 1)

  atoms_bonded = .false.
  do i2 = 1, natoms
    do j2 = 1, natoms
!     Make sure i2 and j2 are different atoms
      if( i2 /= j2 )then
        do l = 1, nbond_types
          if( (species(i2) == bond_type(l, 1) .and. species(j2) == bond_type(l, 2)) .or. &
              (species(i2) == bond_type(l, 2) .and. species(j2) == bond_type(l, 1)) )then
            call get_distance( (/ positions(i2,1:3,step) /), (/ positions(j2,1:3,step) /), L_cell, d)
!           Are i2 and j2 within bond-forming distance?
            if( d < bond_cutoffs(l, 1) )then
              atoms_bonded(i2,j2) = .true.
              atoms_bonded(j2,i2) = .true.
!           Are i2 and j2 below bond-breaking distance and were they bonded before?
            else if( d < bond_cutoffs(l, 2) .and. initialize .and. atom_belongs_to_group(i2) == atom_belongs_to_group(j2) )then
              atoms_bonded(i2,j2) = .true.
              atoms_bonded(j2,i2) = .true.
!           Are i2 and j2 below bond-breaking distance and were they bonded before?
            else if( d < bond_cutoffs(l, 2) .and. atoms_bonded_prev(i2,j2) )then
              atoms_bonded(i2,j2) = .true.
              atoms_bonded(j2,i2) = .true.
            end if
          end if
        end do
      end if
    end do
  end do
end subroutine














subroutine reallocate_arrays_after_topology_rebuild(ngroups, group_still_exists, birth_time, death_time)

  implicit none

! Input variables
  integer :: ngroups, birth_time(:), death_time(:)
  logical :: group_still_exists(:)

! Internal variables
  integer :: k, noriginal_groups

! Check for goups which were not alive for the whole duration of the dynamics and/or which were created
! after reading the trajectory file. Info for these groups pertaining several arrays (partial volumes,
! etc.) need to be reallocated to accommodate all the groups.
  do k = 1, ngroups
!   We need to operate on all the groups that do not exist anymore or which were created when the
!   topology was rebuilt
    if( .not. group_still_exists(k) )then
! COMMENT THIS PRINTING OUT IN PUBLIC VERSION OF CODE
      write(*,*) k
    end if
  end do

end subroutine













subroutine get_distance(pos1, pos2, l, d)

! Returns distance d between atoms i and j according
! to minimum image convention

  implicit none

  real*8 :: d, l(1:3)
  real*4 :: x, y, z, lx, ly, lz, pos1(1:3), pos2(1:3)

  lx = real(l(1))
  ly = real(l(2))
  lz = real(l(3))  

  x = pos2(1) - pos1(1)
  y = pos2(2) - pos1(2)
  z = pos2(3) - pos1(3)

  if( abs(x) > lx/2.)then
    x = x - sign(lx, x)
  end if
  if( abs(y) > ly/2.)then
    y = y - sign(ly, y)
  end if
  if( abs(z) > lz/2.)then
    z = z - sign(lz, z)
  end if

  d = sqrt( x**2 + y**2 + z**2)

end subroutine








recursive subroutine find_neighbors(i, j, natoms, atoms_bonded, atom_visited, atom_belongs_to_new_group, new_groups)

  implicit none

  integer :: i, j, k, natoms, new_groups, atom_belongs_to_new_group(:)
  logical :: atom_visited(:), atoms_bonded(:,:)

  if( .not. atom_visited(j) .and. i /= j .and. atoms_bonded(i,j) )then
    atom_visited(j) = .true.
    atom_belongs_to_new_group(j) = new_groups
    do k = 1, natoms
      call find_neighbors(j, k, natoms, atoms_bonded, atom_visited, atom_belongs_to_new_group, new_groups)
    end do
  end if
end subroutine



end module
