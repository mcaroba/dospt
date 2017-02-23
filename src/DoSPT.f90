program DOSPT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                                            !!!!!
!!!                                                                !!!
!      :MMMMMMMMO          .    MMMMMMMM .. MMMMMMMM$.MMMMMMMMMMM    !
!     .MMM ...MMMM..   +MN=   .MMM..........MMM...MMM ....MMM....    !
!    ..MMM...  MMM...MMMMMMMM  MMMMN.......MMM .. MMM.....MMM....    !
!    ..MMM...  MMM..MMM.  MMM. .MMMMMMM ...MMMMMMMMM.....MMM ....    !
!    .MMM ... MMM,.MMM.  .MMM . .  .ZMMM.. MMMMMO,.......MMM.....    !
!    .MMM . MMMM=..MMM ..MMM7 ..  ..8MMM..MMM ...........MMM.....    !
!    MMMMMMMMMM. .. MMMMMMM   MMMMMMMMD...MMM ..........MMM=....     !
!                                                                    !
!                          DoSPT v0.1.1                              !
!                                                                    !
!  The following distribution of Fortran routines for thermodynamic  !
!      properties calculation, collectively known as DoSPT, has      !
!                          been written by                           !
!                                                                    !
!                           Miguel A. Caro                           !
!           Dept. of Electrical Engineering and Automation           !
!       COMP Centre of Excellence in Computational Nanoscience       !
!                 Aalto University, Espoo, Finland                   !
!                         mcaroba@gmail.com                          !
!                                                                    !
! They are provided for free, in the hope that they will be useful,  !
!    but with no warranty whatsoever, under the Creative Commons     !
!         Attribution-NonCommercial-ShareAlike license               !
!      http://creativecommons.org/licenses/by-nc-sa/3.0/             !
!                                                                    !
!   When publishing work that makes use of the present distribution  !
!                   please have a look and cite                      !
!                                                                    !
!        Miguel A. Caro, Tomi Laurila and Olga Lopez-Acevedo         !
!   "Accurate schemes for calculation of thermodynamic properties    !
!      of liquid mixtures from molecular dynamics simulations"       !
!              J. Chem. Phys. 145, 244504 (2016)                     !
!                                                                    !
! For an in-depth account of the theory underlying the 2PT method    !
!  method, please read (and cite, as appropriate) the original work: !
!                                                                    !
!           S.-T. Lin, M. Blanco, and W. A. Goddard III              !
! "The two-phase model for calculating thermodynamic properties of   !
!  liquids from molecular dynamics: Validation for the phase diagram !
!                     of Lennard-Jones fluids"                       !
!                 J. Chem. Phys. 119, 11792 (2003)                   !
!                                                                    !
!!!          Distribution last updated on 21 Feb. 2016             !!!
!!!!!                                                            !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!******************************************************
! Modules used in this program
  use constants
  use read_input
  use good_bye
!******************************************************






!******************************************************
! This block declares all the variables needed during
! program execution
  implicit none

  integer :: k3, k4, n_volumes, Ng, niter, di_volumes, n_prev, imax
  real*8, allocatable :: wsave(:), r(:), rw(:), r2(:), Nsg
  real*4, allocatable :: velocities(:,:,:), positions(:,:,:)
  real*8 :: temp(1:6), update_bar, res, D_rot_real, D_ideal, D_rot_ideal, sigma_eff, ngroups_eff = 0.d0
  real*8 ::  nu_i(3), gas_integrand(3), solid_integrand(3), dos_subinterval(3), solid_weight(3), quadrature(3)

  real*8, allocatable :: Stotal(:,:), m(:)


  integer :: n_particles
  real*8, allocatable :: mp(:), rp(:,:), vp(:,:)
  real*4, allocatable :: v_rot(:,:,:), v_cm(:,:,:), v_vib(:,:,:)
  real*4, allocatable :: w(:,:,:), eig_group_inst(:,:,:), rxv(:,:,:), volumes(:,:), dist(:,:), e_kin_tot(:), e_kin(:,:)
  real*8, allocatable :: v_rot_group(:,:), v_cm_group(:), v_vib_group(:,:)
  real*8, allocatable :: w_group(:), rxv_group(:,:), volume_group(:), dist_group(:)
  real*8, allocatable :: x(:), y(:), z(:), volumes_temp(:), Sgroup(:,:,:)

  integer :: nparticles
  real*8, allocatable :: eig_group(:,:), symmetry_number_supergroup(:)
  character*16 :: mass_label
  character*16, allocatable :: species(:)
  real*8, allocatable :: delta_group(:,:), degf_group(:,:), degf_supergroup(:,:), weight_group(:)
  real*8, allocatable :: f_group(:,:), f_supergroup(:,:), y_group(:,:), z_group(:,:), mass_group(:)
  real*8, allocatable :: SHSbykB_group(:,:), SHSbykB_supergroup(:,:), SRbykB_group(:), SRbykB_supergroup(:)
  real*8, allocatable :: entropy_group(:,:), entropy_supergroup(:,:), entropy_supergroup_gas(:,:)
  real*8, allocatable :: entropy_supergroup_solid(:,:), degf_supergroup_gas(:,:), degf_supergroup_solid(:,:)
  real*8, allocatable :: entropy1PT_group(:,:), entropy1PT_supergroup(:,:), eig_supergroup(:,:)
  real*8 :: omega

  integer, allocatable :: birth_time(:), death_time(:), lifetime(:)
  logical :: error = .false., topology_has_changed = .false., interpolate_dos = .false.
  real*8, allocatable :: Ssupergroup(:,:,:), mass_supergroup(:), volume_supergroup(:)
  real*8, allocatable :: delta_supergroup(:,:), y_supergroup(:,:), z_supergroup(:,:)
  real*8, allocatable :: natoms_in_supergroup(:), ngroups_in_supergroup_eff(:)


! Filtering stuff must be single precision
  real*4, allocatable :: smooth_x(:), smooth_y(:), smooth_ys(:), smooth_rw(:), smooth_res(:)
  real*4 :: smooth_f

  real*8, allocatable :: sigma_group(:,:), sigma_supergroup(:,:)

! Binary entropy of mixing variables
  real*8 :: v_dash, alpha, beta, m_frac
!******************************************************





!******************************************************
! This block provides the subroutine interface, which is
! required for several subroutines that require variable
! allocation for some of the input parameters
interface

subroutine get_mass(mass_label, mass, nmasses, mass_types, mass_types_values)
  integer, intent(in) :: nmasses
  character*16, intent(in) :: mass_label
  real*8, intent(out) :: mass
  real*8, intent(in)  :: mass_types_values(:)
  character*16, intent(in)  :: mass_types(:)
end subroutine

subroutine get_decomposed_velocities(m, r, v, natoms_in_group, v_rot, v_cm, v_vib, mass, w_group, &
                                     rxv_group, dist_group, eig)
  integer, intent(in) :: natoms_in_group
  real*8, intent(in) :: m(:), r(:,:), v(:,:)
  real*8, intent(out) :: v_rot(:,:), v_cm(:), v_vib(:,:), mass, w_group(:), rxv_group(:,:), dist_group(:), eig(:)
end subroutine

subroutine restore_replica(L, r, natoms_in_group)
  real*8, intent(inout) :: r(:,:)
  real*8, intent(in) :: L(1:3)
  integer, intent(in) :: natoms_in_group
end subroutine

subroutine read_trajectory(n, natoms, tau, L, mode, positions, velocities, estimate_vel, error, &
                           m, nmasses, mass_types, mass_types_values, volumes, volumes_temp, species, di_volumes)
  integer, intent(in) :: n, natoms, nmasses, di_volumes
  real*8, intent(inout) :: m(:)
  real*8, intent(in)  :: mass_types_values(:)
  character*16, intent(in)  :: mass_types(:)
  character*16, intent(out)  :: species(:)
  real*8, intent(in) :: L(1:3), volumes_temp(:), tau
  real*4, intent(inout) :: positions(:,:,:), velocities(:,:,:), volumes(:,:)
  logical, intent(in) :: estimate_vel
  logical, intent(inout) :: error
  character*3, intent(in) :: mode
end subroutine

subroutine get_omega(this_group, ngroups, n_in_group, volume_group, sigma_group, mass_group, mode, omega)
  real*8, intent(out) :: omega
  real*8, intent(in) :: volume_group(:), sigma_group(:), mass_group(:), n_in_group(:)
  integer, intent(in) :: ngroups, this_group
  character*1, intent(in) :: mode
end subroutine

subroutine fluidicity_calculator(nsupergroups, ngroups_in_supergroup, N_DoF, s0, M, T, V, V_voro, niter, res, sigma, f)
  real*8, intent(in) :: T, V, s0(:), N_DoF(:), M(:), V_voro(:), ngroups_in_supergroup(:)
  integer, intent(in) :: nsupergroups
  real*8, intent(out) :: res, f(:), sigma(:)
  integer, intent(out) :: niter
end subroutine

subroutine get_compressibility(nsupergroups, ngroups_in_supergroup, V, sigma, f, xi, z)
  integer, intent(in) :: nsupergroups
  real*8, intent(in) :: V, sigma(:), f(:), ngroups_in_supergroup(:)
  real*8, intent(out) :: z, xi
end subroutine

subroutine get_real_rotational_diffusivity(tau, n, N_DoF, ngroups, group_in_supergroup, birth_time, death_time, j2, w, D)
  integer, intent(in) :: n, ngroups, group_in_supergroup(:,:), j2, birth_time(:), death_time(:)
  real*8, intent(in) :: tau, N_DoF
  real*4, intent(in) :: w(:,:,:)
  real*8, intent(out) :: D
end subroutine

subroutine rebuild_topology(n, natoms, ngroups, nsupergroups, L, positions, species, &
           atoms_in_group, natoms_in_group, symmetry_number_group, atom_belongs_to_group, &
           group_in_supergroup, ngroups_in_supergroup, group_belongs_to_supergroup, &
           nspecies_in_topology, topology_in_supergroup, neach_species_in_topology, &
           ntopology_types, symmetry_number_of_topology, species_in_topology, &
           nbond_types, bond_type, bond_cutoffs, birth_time, death_time, topology_has_changed)
  integer, intent(in) :: n, natoms, nbond_types, ntopology_types
  integer, intent(inout) :: ngroups, nsupergroups
  real*8, intent(in) :: L(1:3), bond_cutoffs(:,:), symmetry_number_of_topology(:)
  real*8, intent(inout) :: symmetry_number_group(:)
  real*4, intent(in) :: positions(:,:,:)
  character*16, intent(in) :: species(:), bond_type(:,:), species_in_topology(:,:)
  integer, intent(inout) :: atoms_in_group(:,:), natoms_in_group(:), atom_belongs_to_group(:)
  integer, intent(out) :: death_time(:), birth_time(:)
  integer, intent(inout) :: ngroups_in_supergroup(:), group_in_supergroup(:,:), group_belongs_to_supergroup(:)
  integer, intent(in) :: nspecies_in_topology(:), topology_in_supergroup(:), neach_species_in_topology(:,:)
  logical, intent(out) :: topology_has_changed
end subroutine

subroutine interpolate_1D(r2, r1, n2, n1, tau2, tau1, f)
  real*8, intent(inout) :: r1(:)
  real*8, intent(in) :: tau1, tau2, f, r2(:)
  integer, intent(in) :: n2, n1
end subroutine

end interface
!******************************************************






!******************************************************
! Constants and units initialization (constants.f90)
  call initialize_constants()

! Print welcome and read in all input files except for trajectory (read_input.f90)
  call print_welcome_and_read_input()

!
!  call print_good_bye()
!******************************************************





!******************************************************
! Allocate arrays for trajectory input and FFT routine
  allocate( wsave(1:2*n+15) )
  allocate( r(1:n) )
  allocate( r2(1:n) )
!  allocate( rw(1:n) )
  allocate( velocities(1:natoms, 1:3, 1:n) )
  allocate( e_kin_tot(1:n) )
  allocate( positions(1:natoms, 1:3, 1:n) )
  allocate( Stotal(1:(n+1)/2, 1:3) )
  allocate( m(1:natoms) )
  allocate( species(1:natoms) )
! Initialize wsave for FFT routine
  call dffti(n, wsave)
! FIX THIS: di_volumes SHOULD BE AVAILABLE TO THE USER AS AN INPUT VARIABLE
! Only calculate Voronoi cell volumes every 100 time steps
  di_volumes = 100
  n_volumes = 1 + (n - 1)/di_volumes
  allocate( volumes(1:natoms, 1:n_volumes) )
  allocate( volumes_temp(1:natoms) ) 
!******************************************************





!******************************************************
! Read trajectory file
  call read_trajectory(n, natoms, tau, L, mode, positions, velocities, estimate_vel, error, &
                       m, nmasses, mass_types, mass_types_values, volumes, volumes_temp, species, di_volumes)
  if(error)then
    execution_error = 98
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! STOP EXECUTION IF THERE IS AN ERROR !!!!!
  end if
!******************************************************





!******************************************************
! Rebuild topology in case bonds break during the dynamics
  if( check_topology )then
    call rebuild_topology(n, natoms, ngroups, nsupergroups, L, positions, species, &
                          atoms_in_group, natoms_in_group, symmetry_number_group, atom_belongs_to_group, &
                          group_in_supergroup, ngroups_in_supergroup, group_belongs_to_supergroup, &
                          nspecies_in_topology, topology_in_supergroup, neach_species_in_topology, &
                          ntopology_types, symmetry_number_of_topology, species_in_topology, &
                          nbond_types, bond_type, bond_cutoffs, birth_time, death_time, topology_has_changed)
  else
    allocate( birth_time(1:ngroups) )
    allocate( death_time(1:ngroups) )
    birth_time = 1
    death_time = n + 1
  end if
  allocate( lifetime(1:ngroups) )
  allocate( weight_group(1:ngroups) )
  lifetime = death_time - birth_time
  do i = 1, ngroups
    weight_group(i) = dfloat(lifetime(i)) / dfloat(n)
  end do
! If the topology has not been rebuilt, then we return the nsupergroups variable to its original
! state
  if( .not. topology_has_changed )then
    nsupergroups = size(ngroups_in_supergroup, 1)
  end if
!******************************************************






!******************************************************
! Allocate property arrays
! The indexing is given by the group index in order of
! appearance in the "groups" file, followed by the atom
! index in order of appearance in the list of atoms for
! that group
! 
! These quantities are given for each atom
  allocate( v_rot(1:natoms, 1:3, 1:n) )
  allocate( v_cm(1:natoms, 1:3, 1:n) )
  allocate( v_vib(1:natoms, 1:3, 1:n) )
  allocate( rxv(1:natoms, 1:3, 1:n) )
  allocate( dist(1:natoms, 1:n) )
  allocate( x(1:natoms) )
  allocate( y(1:natoms) )
  allocate( z(1:natoms) )
! Principal moments of inertia
  allocate( eig_group(1:ngroups, 1:3) )
  eig_group = 0.d0
! These for groups
  allocate( w(1:ngroups, 1:3, 1:n) )
!  allocate( eig_group_inst(1:ngroups, 1:3, 1:n) )
  allocate( e_kin(1:3, 1:n) )
  allocate( volume_group(1:ngroups) )
  allocate( sigma_group(1:ngroups, 1:3) )
  allocate( delta_group(1:ngroups,1:3) )
  allocate( f_group(1:ngroups,1:3) )
  allocate( y_group(1:ngroups,1:3) )
  allocate( z_group(1:ngroups,1:3) )
  allocate( SHSbykB_group(1:ngroups,1:3) )
  allocate( SRbykB_group(1:ngroups) )
  allocate( entropy_group(1:ngroups, 1:3) )
  allocate( entropy1PT_group(1:ngroups, 1:3) )
  allocate( mass_group(1:ngroups) )
  allocate( Sgroup(1:ngroups,1:(n+1)/2,1:3) )
  mass_group = 0.d0
  volume_group = 0.d0
  sigma_group = 0.d0
  entropy_group = 0.d0
  entropy1PT_group = 0.d0
!******************************************************





!******************************************************
! FIX FILTERING TO WORK WITH TOPOLOGY RECONSTRUCTION. In the meantime print warning message.
  if( smooth .and. topology_has_changed )then
    write(*,*)'                                       |'
    write(*,*)'WARNING: usage of DoS filtering togeth-|'
    write(*,*)'er with topology reconstruction is not |'
    write(*,*)'currently implemented. Expect nonsense |'
    write(*,*)'from this simulation!                  |'
    write(*,*)'                                       |'
    write(*,*)'.......................................|'
  end if
! Take care of the filtering variables
  if(smooth)then
    allocate( smooth_x(1:(n+1)/2) )
    allocate( smooth_y(1:(n+1)/2) )
    allocate( smooth_ys(1:(n+1)/2) )
    allocate( smooth_rw(1:(n+1)/2) )
    allocate( smooth_res(1:(n+1)/2) )
    smooth_f = 4.*sigma_nu*tau/float(n)
!   The FFT array elements must be regularly spaced, thus the
!   contents of smooth_x are irrelevant as long as they also are
!   regularly spaced
    do i=1,(n+1)/2
      smooth_x(i) = float(i)
    end do
  end if
!******************************************************





!******************************************************
! Calculate the trajectory-averaged Voronoi cell volumes
! for each group
  do j= 1, ngroups
    k = 0
    do i = 1, n, di_volumes
!     Transform trajectory index into volume index
      i2 = 1 + i/di_volumes
!     Add the volumes only during the time when the group is alive
      if( birth_time(j) <= i .and. death_time(j) > i )then
        k = k + 1
        do j2=1,natoms_in_group(j)
          volume_group(j) = volume_group(j) + volumes(atoms_in_group(j,j2),i2)
        end do
      end if
    end do
!   If a group is short-lived or the Voronoi sampling is low, it could be that
!   there are no Voronoi volumes to use (k = 0). In that case, we compute one
!   Voronoi volume for the group at a time exactly in between birth and death
    if( k > 0 )then
      volume_group(j) = volume_group(j) / dfloat(k)
    else
      i2 = 1 + (birth_time(j) + death_time(j)) / (2*di_volumes)
      do j2=1,natoms_in_group(j)
        volume_group(j) = volume_group(j) + volumes(atoms_in_group(j,j2),i2)
      end do
    end if
  end do
!******************************************************





!******************************************************
! Calculate translational, rotational and vibrational
! velocities for the atoms in each group
  write(*,*)'                                       |'
  write(*,*)'Carrying out velocity decomposition... |'
  write(*,*)'                                       |'
  write(*,*)'Progress:                              |'
  write(*,*)'                                       |'
  write(*,'(1X,A)',advance='no')'['
  update_bar = dfloat(ngroups*n)/36.d0
  k2 = 1
  do j=1,ngroups
!   Allocate variables
    allocate( mp(1:natoms_in_group(j)) )
    allocate( rp(1:natoms_in_group(j), 1:3) )
    allocate( vp(1:natoms_in_group(j), 1:3) )
    allocate( v_rot_group(1:natoms_in_group(j), 1:3) )
    allocate( v_cm_group(1:3) )
    allocate( w_group(1:3) )
    allocate( rxv_group(1:natoms_in_group(j), 1:3) )
    allocate( dist_group(1:natoms_in_group(j)) )
    allocate( v_vib_group(1:natoms_in_group(j), 1:3) )
!   Assign a mass
    do j2=1,natoms_in_group(j)
      mp(j2) = m( atoms_in_group(j,j2) )
    end do
    do i=1,n
!     Update progress bar every ngroups*n/36 iterations
      if(dfloat((j-1)*n+i) > dfloat(k2)*update_bar)then
        write(*,'(A)', advance='no')'='
        k2 = k2 + 1
      end if
!     Perform operations only during the time the group is alive
      if( birth_time(j) <= i .and. death_time(j) > i )then
!       Read in molecule info
        do j2=1,natoms_in_group(j)
          do k=1,3
            rp(j2,k) = positions(atoms_in_group(j,j2), k, i)
            vp(j2,k) = velocities(atoms_in_group(j,j2), k, i)
          end do
        end do
!       Get decomposed values for present molecule
        call restore_replica(L, rp, natoms_in_group(j))
        call get_decomposed_velocities(mp, rp, vp, natoms_in_group(j), v_rot_group, v_cm_group, v_vib_group, mass_group(j), &
                                       w_group, rxv_group, dist_group, temp(1:3))
        w(j, 1:3, i) = w_group(1:3)
!       Calculate average principal moments of inertia (amu * nm^2)
        eig_group(j,1:3) = eig_group(j,1:3) + temp(1:3)/dfloat(lifetime(j))
!       and instantaneous ones
!        eig_group_inst(j,1:3,i) = temp(1:3)
!       Pass values to atoms
        do j2=1,natoms_in_group(j)
          do k=1,3
            v_rot(atoms_in_group(j,j2), k, i) = v_rot_group(j2, k)
            v_cm(atoms_in_group(j,j2), k, i) = v_cm_group(k)
            v_vib(atoms_in_group(j,j2), k, i) = v_vib_group(j2, k)
            rxv(atoms_in_group(j,j2), k, i) = rxv_group(j2, k)
          end do
          dist(atoms_in_group(j,j2), i) = dist_group(j2)
        end do
      end if
    end do
!   Deallocate some variables
    deallocate(mp, rp, vp)
    deallocate(v_rot_group, v_cm_group, v_vib_group, w_group, rxv_group, dist_group)
  end do
  deallocate(positions)
  write(*,'(A4)')'=] |'
  write(*,*)'                                       |'
  write(*,*)'.......................................|'
!******************************************************





!******************************************************
! Calculate kinetic energy at each point in the trajectory from the velocities
! units are amu * nm^2 / ps^2, and we convert them to eV
  e_kin = 0.
  e_kin_tot = 0.
  do j=1,natoms
    do i=1,n
      e_kin_tot(i) = e_kin_tot(i) + 0.5 * m(j) * dot_product(velocities(j,1:3,i), velocities(j,1:3,i))
      e_kin(1,i) = e_kin(1,i) + 0.5 * m(j) * dot_product(v_cm(j,1:3,i), v_cm(j,1:3,i))
      e_kin(2,i) = e_kin(2,i) + 0.5 * m(j) * dot_product(v_rot(j,1:3,i), v_rot(j,1:3,i))
      e_kin(3,i) = e_kin(3,i) + 0.5 * m(j) * dot_product(v_vib(j,1:3,i), v_vib(j,1:3,i))
    end do
  end do
  e_kin_tot = e_kin_tot * amu * nm**2 / ps**2 / eV
  e_kin = e_kin * amu * nm**2 / ps**2 / eV
  open(unit=10, file="e_kin", status="unknown")
  do i=1,n
    write(10,*) dfloat(i-1)/dfloat(n-1)*tau, e_kin_tot(i), e_kin(1:3,i)
  end do
  close(10)
  deallocate(e_kin)
  deallocate(e_kin_tot)
  deallocate(velocities)
!******************************************************





!******************************************************
! Calculate density of states
! ATTENTION, this is a parallel block
  write(*,*)'                                       |'
  write(*,*)'Calculating density of states...       |'
  write(*,*)'                                       |'
! Check prime number factorization for n and print warning if high factors are found
! FIX THIS: WHEN GROUPS BREAK AND THE TOPOLOGY IS REBUILT SOMETIMES THE lifetime(j) YIELDS HIGH FACTOR
! DECOMPOSITION (THIS IS INDEPENDENT OF n). FIGURE OUT HOW TO DEAL WITH THAT
  call factorize(n)
! Print progress bar only if we have more than 36 degrees of freedom
  k3 = 0
  k4 = 0
  do j=1,ngroups
    do j2=1,natoms_in_group(j)
      k3 = k3 + 1
    end do
  end do
  if(k3*9 >= 36)then
    write(*,*)'Progress:                              |'
    write(*,*)'                                       |'
    write(*,'(1X,A)',advance='no')'['
  end if
  update_bar = dfloat(9*k3)/36.d0
  k3 = 1
!
  n_prev = n
  Sgroup = 0.d0
!$omp parallel do firstprivate(wsave) private(r,rw,r2,i,i2,j,j2,k,k2,n_prev,smooth_y,smooth_ys,smooth_rw,smooth_res) &
!$omp& shared(Sgroup,k3,k4)
  do j=1,ngroups 
    do j2=1,natoms_in_group(j)
      do k=1,3
!     k2 now is the 3 components (cm, rot and vib)
      do k2=1,3
!       Print progress bar
!$omp critical
        k4 = k4 + 1
        if(float(k4) > float(k3)*update_bar .and. update_bar > 0.9999d0)then
          write(*,'(A)', advance='no')'='
          k3 = k3 + 1
        end if
!$omp end critical
!       End printing progress bar
!       We need to Fourier-transform the arrays according to the life cycle of the different groups
        r = 0.d0
        do i = birth_time(j), death_time(j)-1
          i2 = i + 1 - birth_time(j)
!         Translational
          if(k2 == 1)then
            r(i2) = v_cm(atoms_in_group(j,j2), k, i)
!         Rotational
          else if(k2 == 2)then
            r(i2) = rxv(atoms_in_group(j,j2), k, i) / dist(atoms_in_group(j,j2), i)
!         Vibrational
          else if(k2 == 3)then
            r(i2) = v_vib(atoms_in_group(j,j2), k, i)
          end if
        end do
!       Ask for the lifetime to be at least 10 time steps before calculating DoS
!       to prevent transitioning groups from messing things up. DoS = 0 otherwise
        if( lifetime(j) > 10 )then
!         Call FFT routine
!         For groups which do not live for the whole dynamics, we use a different array
          if( lifetime(j) /= n )then
            deallocate(r2)
            allocate( r2(1:lifetime(j)) )
            r2(1:lifetime(j)) = r(1:lifetime(j))
            if( n_prev /= lifetime(j) )then
              call dffti(lifetime(j), wsave)
              n_prev = lifetime(j)
            end if
            call dfftf(lifetime(j), r2, wsave)
            interpolate_dos = .true.
!         For groups which live for the whole dynamics we're happy
          else
            if( n_prev /= n )then
              call dffti(n, wsave)
              n_prev = n
            end if
            call dfftf(n, r, wsave)
            interpolate_dos = .false.
          end if
!         Save DoS for this group
          if( interpolate_dos )then
!           We perform the interpolation on the density of states
            r2(1) = r2(1)**2
            do i = 2, (lifetime(j)+1)/2
              r2(i) = r2(2*i-2)**2+r2(2*i-1)**2
            end do
!           This is the density of states on the r2 mesh
            r2 = r2 * tau*dfloat(lifetime(j))/dfloat(n) / dfloat(lifetime(j))**2
!           We interpolate from the r2 mesh to the r mesh. Remember the period for the interpolation routine
!           is actually in frequency domain: period_nu = (n-1)/2/tau
            call interpolate_1D(r2, r, (lifetime(j)+1)/2, (n+1)/2, &
                                dfloat(lifetime(j)-1)/2.d0/(tau*dfloat(lifetime(j))/dfloat(n)), &
                                dfloat(n-1)/2.d0/tau, 1.d0)
            Sgroup(j,1,k2) = Sgroup(j,1,k2) + m(atoms_in_group(j,j2)) * r(1)
            do i=2,(n+1)/2
              Sgroup(j,i,k2) = Sgroup(j,i,k2) + m(atoms_in_group(j,j2)) * r(i)
            end do
          else
            Sgroup(j,1,k2) = Sgroup(j,1,k2) + m(atoms_in_group(j,j2)) * r(1)**2*tau/dfloat(n)**2
            do i=2,(n+1)/2
              Sgroup(j,i,k2) = Sgroup(j,i,k2) + m(atoms_in_group(j,j2)) * (r(2*i-2)**2+r(2*i-1)**2)*tau/dfloat(n)**2
            end do
          end if
        end if
      end do
      end do
    end do
! If filtering is enabled, perform the LOWESS smoothing now
    if(smooth)then
      do k2=1,3
!       Keep in mind that the lowess routine works in single precision (real*4)
        smooth_y(1:(n+1)/2) = Sgroup(j,1:(n+1)/2,k2)
        call lowess(smooth_x, smooth_y, (n+1)/2, smooth_f, 2, 0., smooth_ys, smooth_rw, smooth_res)
        Sgroup(j,1:(n+1)/2,k2) = smooth_ys(1:(n+1)/2)
      end do
    end if
  end do
! Calculate total DoS
  Stotal = 0.d0
  do j=1,ngroups
    do k2=1,3
      Stotal(1,k2) = Stotal(1,k2) + weight_group(j) * Sgroup(j,1,k2)
      do i=2,(n+1)/2
        Stotal(i,k2) = Stotal(i,k2) + weight_group(j) * Sgroup(j,i,k2)
      end do
    end do
  end do
! Calculate supergroup DoS and masses/volumes
  if(calc_supergroups)then
    allocate( Ssupergroup(1:nsupergroups,1:(n+1)/2,1:3) )
    allocate( mass_supergroup(1:nsupergroups) )
    allocate( symmetry_number_supergroup(1:nsupergroups) )
    allocate( volume_supergroup(1:nsupergroups) )
    allocate( sigma_supergroup(1:nsupergroups, 1:3) )
    allocate( natoms_in_supergroup(1:nsupergroups) )
    allocate( ngroups_in_supergroup_eff(1:nsupergroups) )
    allocate( eig_supergroup(1:nsupergroups,1:3) )
    sigma_supergroup = 0.d0
    Ssupergroup = 0.d0
    mass_supergroup = 0.d0
    symmetry_number_supergroup = 0.d0
    volume_supergroup = 0.d0
    natoms_in_supergroup = 0.d0
    eig_supergroup = 0.d0
    ngroups_in_supergroup_eff = 0.d0
    ngroups_eff = 0.d0
    do j=1,nsupergroups
      do j2=1,ngroups_in_supergroup(j)
        ngroups_in_supergroup_eff(j) = ngroups_in_supergroup_eff(j) + weight_group(group_in_supergroup(j,j2))
      end do
      ngroups_eff = ngroups_eff + ngroups_in_supergroup_eff(j)
    end do
    do j=1,nsupergroups
      do j2=1,ngroups_in_supergroup(j)
        eig_supergroup(j,1:3) = eig_supergroup(j,1:3) + weight_group(group_in_supergroup(j,j2))* &
                                eig_group(group_in_supergroup(j,j2),1:3)/ngroups_in_supergroup_eff(j)
        symmetry_number_supergroup(j) = symmetry_number_supergroup(j) + weight_group(group_in_supergroup(j,j2))* &
                                        symmetry_number_group(group_in_supergroup(j,j2))/ngroups_in_supergroup_eff(j)
        mass_supergroup(j) = mass_supergroup(j) + weight_group(group_in_supergroup(j,j2))* &
                             mass_group(group_in_supergroup(j,j2))
        volume_supergroup(j) = volume_supergroup(j) + weight_group(group_in_supergroup(j,j2))* &
                               volume_group(group_in_supergroup(j,j2))
        natoms_in_supergroup(j) = natoms_in_supergroup(j) + weight_group(group_in_supergroup(j,j2))* &
                                  dfloat(natoms_in_group(group_in_supergroup(j,j2)))
        do k2=1,3
          Ssupergroup(j,1,k2) = Ssupergroup(j,1,k2) + weight_group(group_in_supergroup(j,j2))* &
                                Sgroup(group_in_supergroup(j,j2),1,k2)
          do i=2,(n+1)/2
            Ssupergroup(j,i,k2) = Ssupergroup(j,i,k2) + weight_group(group_in_supergroup(j,j2))* &
                                  Sgroup(group_in_supergroup(j,j2),i,k2)
          end do
        end do
      end do
!     Mass of the supergroup is average mass of its groups. The groups should be equivalent
!     (e.g., all water molecules) and so this should not be necessary, but it's done for
!     consistency in case the groups are not equivalent (for some reason)
      mass_supergroup(j) = mass_supergroup(j) / ngroups_in_supergroup_eff(j)
    end do
  end if
  if(update_bar > 0.9999d0)then
    write(*,'(A4)')'=] |'
    write(*,*)'                                       |'
    write(*,*)'.......................................|'
  else
    write(*,*)'Done.                                  |'
    write(*,*)'                                       |'
    write(*,*)'.......................................|'
  end if
!******************************************************





!******************************************************
! Write DoS to file
  if(write_dos)then
  open(unit=10, file="dos", status="unknown")
!
  write(*,*)'                                       |'
  write(*,*)'Writing DoS to file "dos"...           |'
  write(*,*)'                                       |'
  write(*,*)'Progress:                              |'
  write(*,*)'                                       |'
  write(*,'(1X,A)',advance='no')'['
  update_bar = dfloat(ngroups*imax)/36.d0
  k2 = 1
!
  write(10,*) "# Total density of states for all the groups involved"
! FIX THIS. MAKE IT POSSIBLE FOR THE USER TO CHANGE MAX CUTOFF FOR DoS OUTPUT
! Choose upper cutoff for DoS printing
  if( .true. )then
    imax = int(150.d0 * tau) - 1
  else
    imax = (n+1)/2
  end if
  do i=1,imax
    write(10, *) dfloat(i-1)/tau, conv1 * twobykT * Stotal(i,1), conv1 * twobykT * Stotal(i,2), conv1 * twobykT * Stotal(i,3)
  end do
  do j=1,ngroups
    write(10,*) 
    write(10,*) "# Density of states (in ps) for group", j
    do i=1,imax
!     Update progress bar every ngroups*n/36 iterations
      if(dfloat((j-1)*imax+i) > dfloat(k2)*update_bar)then
        write(*,'(A)', advance='no')'='
        k2 = k2 + 1
      end if
      write(10, *) dfloat(i-1)/tau, conv1 * twobykT * Sgroup(j,i,1), conv1 * twobykT * Sgroup(j,i,2), &
                                    conv1 * twobykT * Sgroup(j,i,3)
    end do
  end do
  write(*,'(A4)')'=] |'
  write(*,*)'                                       |'
  write(*,*)'.......................................|'
  close(10)
! Same for supergroups
  if(calc_supergroups)then
    open(unit=10, file="dos_sg", status="unknown")
!
    write(*,*)'                                       |'
    write(*,*)'Writing DoS to file "dos_sg"...        |'
    write(*,*)'                                       |'
    write(*,*)'Progress:                              |'
    write(*,*)'                                       |'
    write(*,'(1X,A)',advance='no')'['
    update_bar = dfloat(nsupergroups*imax)/36.d0
    k2 = 1
!
    do j=1,nsupergroups
      write(10,*) "# Density of states (in ps) for supergroup", j
      do i=1,imax
!     Update progress bar every ngroups*n/36 iterations
        if(dfloat((j-1)*imax+i) > dfloat(k2)*update_bar)then
          write(*,'(A)', advance='no')'='
          k2 = k2 + 1
        end if
        write(10, *) dfloat(i-1)/tau, conv1 * twobykT * Ssupergroup(j,i,1), conv1 * twobykT * Ssupergroup(j,i,2), &
                                      conv1 * twobykT * Ssupergroup(j,i,3)
      end do
      write(10,*)
    end do
    write(*,'(A4)')'=] |'
    write(*,*)'                                       |'
    write(*,*)'.......................................|'
    close(10)
  end if
  end if
!******************************************************





!******************************************************
! Get diffusivity and fluidicity values for each group
  write(*,*)'                                       |'
  write(*,*)'Calculating partitioning properties    |'
  write(*,*)'(diffusivity, fluidicity, etc.) for    |'
  write(*,*)'all groups and supergroups and writing |'
  write(*,*)'to file "fluidicity"...                |'
  write(*,*)'                                       |'
  allocate( degf_group(1:ngroups, 1:3) )
  degf_group = 0.d0
! Delta is obtained differently for the different components. The (2*S0)/(9M)
! term (J. Chem. Theory and Comp. 7, 1893 (2011), Eq. (14)) is substituted
! by (2*S0)/(3*degrees_of_freedom). We do not make any assumptions about how
! the degrees of freedom are distributed between translational, rotational
! and vibrational motion. Instead, we calculate the effective degrees of freedom
! as the integral of the density of states.
  do j=1,ngroups
    do k=1,3
      do i=1,(n+1)/2
!       Calculate the integral of the DoS:
        degf_group(j,k) = degf_group(j,k) + conv1 * twobykT * Sgroup(j,i,k) / tau
      end do
! FIX THIS: IN GENERAL, IT DOES NOT MAKE A LOT OF SENSE TO COMPUTE THESE PROPERTIES FOR SINGLE GROUPS
! A REFACTORING OF THE CODE SHOULD BE CARRIED OUT TO GIVE THE INDIVIDUAL CONTRIBUTION OF EACH GROUP
! TO ITS SUPERGROUP'S ENTROPY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(1 == 1)then
!       This menthod of calculating the fluidicity is essentially equivalent to Lin's method,
!       except for the inclusion of the "omega term" heuristically in the definition of the
!       normalized diffusivity "delta". This method is now deprecated and the more sophisticated
!       method based on numerical minimization of a penalty function should be employed.
!       For a monocomponent system, the results should be identical.
!       Available for testing and debugging purposes.
!       Calculate Delta:
        call get_omega(j, ngroups, dfloat(1+0*natoms_in_group), volume_group, sigma_group(1:ngroups,k), mass_group, "v", omega)
        call get_delta(conv1 * twobykT * Sgroup(j,1,k), T, conv2, 1.d0, &
                       degf_group(j,k), mass_group(j), V, delta_group(j,k), omega)
!       And f:
        call find_f(delta_group(j,k), f_group(j,k))
!       Packing fraction and compressibility:
        y_group(j,k) = f_group(j,k)**(5.d0/2.d0) / delta_group(j,k)**(3.d0/2.d0)
        z_group(j,k) = (1.d0 + y_group(j,k) + y_group(j,k)**2d0 - y_group(j,k)**3.d0) / (1.d0 - y_group(j,k))**3.d0
        sigma_group(j,k) = (6.d0/pi * y_group(j,k) / f_group(j,k) * V / dfloat(1) )**(1.d0/3.d0)
      end if
    end do
  end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(1 == 0)then
!   Get fluidicity from calculator. This is a bad idea for single groups, because the simultaneous optimization
!   of the hard-sphere diameters of many groups will be very inefficient. For single groups fluidicity, etc.
!   is not well defined anyway, and not used for entropy calculations. We can therefore use the heuristic approach above.
!   For single groups, use this code only for testing and debugging.
    do k = 1, 3
      call fluidicity_calculator(ngroups, dfloat(1+0*natoms_in_group), degf_group(1:ngroups,k), &
                                 conv1 * twobykT * Sgroup(1:ngroups,1,k), &
                                 mass_group, T, V, volume_group, niter, res, sigma_group(1:ngroups,k), f_group(1:ngroups,k))
      do j = 1, ngroups
!       Get partial compressibility and partial packing fraction
        call get_compressibility(1, dfloat(1+0*natoms_in_group(j:j)), volume_group(j), sigma_group(j:j, k), &
                                 f_group(j:j, k), y_group(j,k), z_group(j,k))
      end do
    end do
  end if
! Write to file
  open(unit=10, file="fluidicity_g", status="unknown")
  write(10,*) "# Group no.; effective deg. of freedom [CM; rot; vib]; Voronoi volume (nm^3); &
               &fluidicity f [CM; rot; vib]"
  do j=1,ngroups
    write(10,'(I8,2X,F10.4,2X,F10.4,2X,F10.4,2X,F10.4,2X,F6.4,2X,F6.4,2X,F6.4)') &
          j, degf_group(j,1:3), volume_group(j), f_group(j,1:3)
  end do
  close(10)
! Repeat the same as above for supergroups this time
  if(calc_supergroups)then
    allocate( degf_supergroup(1:nsupergroups, 1:3) )
    allocate( delta_supergroup(1:nsupergroups, 1:3) )
    allocate( f_supergroup(1:nsupergroups, 1:3) )
    allocate( y_supergroup(1:nsupergroups, 1:3) )
    allocate( z_supergroup(1:nsupergroups, 1:3) )
    degf_supergroup = 0.d0
    do j=1,nsupergroups
      do k=1,3
        do i=1,(n+1)/2
!         Calculate the integral of the Dos:
          degf_supergroup(j,k) = degf_supergroup(j,k) + conv1 * twobykT * Ssupergroup(j,i,k) / tau
        end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if((k == 1 .and. .not. f_opt) .or. (k == 2 .and. .not. f_rot_opt) .or. (k == 3))then
!         This menthod of calculating the fluidicity is essentially equivalent to Lin's method,
!         except for the inclusion of the "omega term" heuristically in the definition of the
!         normalized diffusivity "delta". This method is now deprecated and the more sophisticated
!         method based on numerical minimization of a penalty function should be employed.
!         For a monocomponent system, the results should be identical.
!         Available for testing and debugging purposes.
!         Vibrational fluidicity is not a well-defined concept and is calculated only for curious
!         purposes using this method. Vibrational fluidicity is not used to calculate vibrational entropy.
!         Calculate Delta:
          if(hs_formalism == "lin")then
            omega = 1
            temp(1) = volume_supergroup(j)
          else
            call get_omega(j, nsupergroups, ngroups_in_supergroup_eff, volume_supergroup, &
                           sigma_supergroup(1:nsupergroups,k), mass_supergroup, "v", omega)
            temp(1) = V
          end if
          call get_delta(conv1 * twobykT * Ssupergroup(j,1,k), T, conv2, ngroups_in_supergroup_eff(j), &
                         degf_supergroup(j,k), mass_supergroup(j), temp(1), delta_supergroup(j,k), omega)
!         And f:
          call find_f(delta_supergroup(j,k), f_supergroup(j,k))
!         Packing fraction and compressibility:
          y_supergroup(j,k) = f_supergroup(j,k)**(5.d0/2.d0) / delta_supergroup(j,k)**(3.d0/2.d0)
          z_supergroup(j,k) = (1.d0 + y_supergroup(j,k) + y_supergroup(j,k)**2d0 - y_supergroup(j,k)**3.d0) / &
                              (1.d0 - y_supergroup(j,k))**3.d0
          sigma_supergroup(j,k) = (6.d0/pi * y_supergroup(j,k) / f_supergroup(j,k) * volume_supergroup(j) / &
                                  ngroups_in_supergroup_eff(j) )**(1.d0/3.d0)
        end if
      end do
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Get fluidicity from calculator.
    if(f_opt)then
!     Translational part
      k = 1
      call fluidicity_calculator(nsupergroups, ngroups_in_supergroup_eff, degf_supergroup(1:nsupergroups,k), &
                                 conv1 * twobykT * Ssupergroup(1:nsupergroups,1,k), mass_supergroup, T, V, &
                                 volume_supergroup, niter, res, sigma_supergroup(1:nsupergroups,k), &
                                 f_supergroup(1:nsupergroups,k))
      do j = 1, nsupergroups
!       Get partial compressibility and partial packing fraction
        call get_compressibility(1, ngroups_in_supergroup_eff(j:j), volume_supergroup(j), sigma_supergroup(j:j, k), &
                                 f_supergroup(j:j, k), y_supergroup(j,k), z_supergroup(j,k))
      end do
    end if
    if(f_rot_opt)then
!     Rotational fluidicity optimized requires input from the translational part (HS diameters are inherited)
      k = 2
      sigma_supergroup(1:nsupergroups, k) = sigma_supergroup(1:nsupergroups, 1)
!     Obtain rotational fluidicity
      do j = 1, nsupergroups
!       Real rotational diffusion coefficient -> D_rot_real
        call get_real_rotational_diffusivity(tau, n, degf_supergroup(j,k), ngroups_in_supergroup(j), &
                                             group_in_supergroup, birth_time, death_time, j, w, D_rot_real)
!        call get_diffusivity(Ssupergroup(j,1,k), T, mass_supergroup(j), degf_supergroup(j,k), D_rot_real)
!       Ideal translational diffusion coefficient -> D_ideal
        call get_omega(j, nsupergroups, ngroups_in_supergroup_eff, volume_supergroup, &
                       sigma_supergroup(1:nsupergroups, k), mass_supergroup, "v", omega)
        call get_zero_pressure_diffusivity(sigma_supergroup(j,k), omega, T, mass_supergroup(j), &
                                           ngroups_in_supergroup_eff(j), V, D_ideal)
!       Ideal rotational diffusion coefficient -> D_rot_ideal
        call get_ideal_rotational_diffusivity(sigma_supergroup(j,k), T, mass_supergroup(j), &
                                              ngroups_in_supergroup_eff(j), volume_supergroup(j), D_ideal, D_rot_ideal)
        D_rot_ideal = D_rot_ideal / 9.d0 * ( eig_supergroup(j,1) + eig_supergroup(j,2) + eig_supergroup(j,3) ) * &
                    ( 1.d0/eig_supergroup(j,1) + 1.d0/eig_supergroup(j,2) + 1.d0/eig_supergroup(j,3) )
!       Fluidicity
        if( degf_supergroup(j,k) < 1.d-10 )then
          f_supergroup(j, k) = 0.d0
        else
          f_supergroup(j, k) = D_rot_real / D_rot_ideal
        end if
      end do
    end if
!   Write to file
    open(unit=10, file="fluidicity", status="unknown")
    write(10,*) "# Supergroup no.; effective deg. of freedom [CM; rot; vib]; Voronoi volume (nm^3); &
                 &fluidicity f [CM; rot; vib]; HS diameters (nm) [CM, rot, vib]"
    do j=1,nsupergroups
      write(10,'(I8,2X,F10.4,2X,F10.4,2X,F10.4,2X,F10.4,2X,F6.4,2X,F6.4,2X,F6.4,2X,F10.4,2X,F10.4,2X,F10.4)') &
            j, degf_supergroup(j,1:3), volume_supergroup(j), f_supergroup(j,1:3), sigma_supergroup(j,1:3)
    end do
  end if
  close(10)
  write(*,*)'Done.                                  |'
  write(*,*)'                                       |'
  write(*,*)'.......................................|'
!******************************************************





!******************************************************
! Calculate thermodynamic properties
!
  write(*,*)'                                       |'
  write(*,*)'Calculating entropy and writing to file|'
  write(*,*)'"entropy"...                           |'
  write(*,*)'                                       |'
  write(*,*)'Progress:                              |'
  write(*,*)'                                       |'
  write(*,'(1X,A)',advance='no')'['
  update_bar = dfloat((ngroups+nsupergroups)*3*(n+1)/2)/35.d0
  k2 = 0
  k3 = 0
!
! Calculate gas weighting functions
! For all the groups
  do j=1,ngroups
    do k = 1,3
!     Calculate SHS divided by kB
      Ng = 1
! FIX THIS. AGAIN (SEE ABOVE) IT MAKES LITTLE SENSE TO CALCULATE THIS FOR INDIVIDUAL GROUPS. RETHINK THIS PART OF THE CODE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!temp(1) = V
temp(1) = volume_group(j)
      SHSbykB_group(j,k) = 5.d0/2.d0 &
         + dlog(conv3 * (2.d0 * pi * mass_group(j) * kB * T / h**2.d0)**(3.d0/2.d0) * &
           temp(1) * z_group(j,k) / dfloat(Ng) / f_group(j,k)) &
         + y_group(j,k) * (3.d0 * y_group(j,k) - 4.d0) / (1.d0 - y_group(j,k))**2.d0
!     If fluidicity is zero truncate SHS to zero
      if(f_group(j,k) < 1.d-10)then
        SHSbykB_group(j,k) = 0.d0
      end if
!     If k = 2, calculate also rotational entropy SR divided by kB
      if(k == 2)then
!       Store T divided by rotational temperature in temp variable
        temp(1:3) = conv4 * T / h**2.d0 * 8.d0 * pi**2.d0 * kB * eig_group(j,1:3)
        SRbykB_group(j) = - dlog(pi * symmetry_number_group(j))
        do i=1,3
          if(temp(i) > 1.d0)then
            SRbykB_group(j) = SRbykB_group(j) + 0.5 + 0.5 * dlog(pi * temp(i))
          end if
        end do
!       If fluidicity is zero truncate SR to zero
        if(f_group(j,k) < 1.d-10)then
          SRbykB_group(j) = 0.d0
        end if
      end if
    end do
  end do
! For all the supergroups
  if(calc_supergroups)then
    allocate( SHSbykB_supergroup(1:nsupergroups, 1:3) )
    allocate( SRbykB_supergroup(1:nsupergroups) )
    do j=1,nsupergroups
      do k = 1,3
        Nsg = ngroups_in_supergroup_eff(j)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!temp(1) = V
temp(1) = volume_supergroup(j)
        SHSbykB_supergroup(j,k) = 5.d0/2.d0 &
           + dlog(conv3 * (2.d0 * pi * mass_supergroup(j) * kB * T / h**2.d0)**(3.d0/2.d0) * &
             temp(1) * z_supergroup(j,k) / Nsg / f_supergroup(j,k)) &
           + y_supergroup(j,k) * (3.d0 * y_supergroup(j,k) - 4.d0) / (1.d0 - y_supergroup(j,k))**2.d0
!       If fluidicity is zero truncate SHS to zero
        if(f_supergroup(j,k) < 1.d-10)then
          SHSbykB_supergroup(j,k) = 0.d0
        end if
!       If k = 2, calculate also rotational entropy SR divided by kB
        if(k == 2)then
!         Store T divided by rotational temperature in temp variable
          temp(1:3) = conv4 * T / h**2.d0 * 8.d0 * pi**2.d0 * kB * eig_supergroup(j,1:3)
          SRbykB_supergroup(j) = - dlog(pi * symmetry_number_supergroup(j))
          do i=1,3
            if(temp(i) > 1.d0)then
              SRbykB_supergroup(j) = SRbykB_supergroup(j) + 0.5 + 0.5 * dlog(pi * temp(i))
            end if
          end do
!         If fluidicity is zero truncate SR to zero
          if(f_supergroup(j,k) < 1.d-10)then
            SRbykB_supergroup(j) = 0.d0
          end if
        end if
      end do
    end do
  end if
! Calculate entropy in eV / K
! For groups
  open(unit=10, file="entropy", status="unknown")
  open(unit=20, file="entropy1PT", status="unknown")
  write(10,*)"# Group; Entropies in eV/K [trans, rot, vib]; Total entropy in eV/K; Total entropy in J/K"
  write(20,*)"# Group; Entropies in eV/K [trans, rot, vib]; Total entropy in eV/K; Total entropy in J/K"
  do j=1,ngroups
    do k=1,3
!     The integral is performed using a quadrature rule for each subinterval (no. of subintervals is
!     no. of points minus one)
      do i=1,(n+1)/2-1
!       Update progress bar
        k2 = k2 + 1
        if(float(k2) > float(k3)*update_bar)then
          write(*,'(A)', advance='no')'='
          k3 = k3 + 1
        end if
!
!       DoS and frequency in the subinterval. Quadrature coefficients
        dos_subinterval = (/ Sgroup(j,i,k), (Sgroup(j,i,k) + Sgroup(j,i+1,k))/2.d0, Sgroup(j,i+1,k) /)
        nu_i = (/ dfloat(i-1) / tau, (dfloat(i-1) + 0.5d0) / tau, dfloat(i) / tau /)
        quadrature = (/ 1.d0/6.d0, 4.d0/6.d0, 1.d0/6.d0 /)
!
!       We skip the gas part if the no. of degrees of freedom is zero or the fluidicity is zero (i.e., numerically very small)
!       SHS must also be positive
        if(degf_group(j,k) < 1.d-10 .or. f_group(j,k) < 1.d-10 .or. SHSbykB_group(j,k) < 0.d0)then
          gas_integrand = 0.d0
        else if(k == 3 .and. .not. vibrational_gas)then
          gas_integrand = 0.d0
        else
          temp(1) = degf_group(j,k)
          gas_integrand = pi * conv1 * twobykT * Sgroup(j,1,k) * nu_i / 2.d0 / f_group(j,k) / temp(1)
          gas_integrand = conv1 * twobykT * Sgroup(j,1,k) / (1.d0 + gas_integrand**2.d0)
!         The gas DoS cannot be larger than the total DoS
          do j2=1,3
            if(gas_integrand(j2) > conv1 * twobykT * dos_subinterval(j2))then
              gas_integrand(j2) = conv1 * twobykT * dos_subinterval(j2)
            end if
          end do
        end if
!       Calculate gas entropy
        if(k == 1 .or. k == 3)then
          entropy_group(j,k) = entropy_group(j,k) + dot_product(gas_integrand,quadrature) * SHSbykB_group(j,k) / 3.d0
        else if(k == 2)then
          entropy_group(j,k) = entropy_group(j,k) + dot_product(gas_integrand,quadrature) * SRbykB_group(j) / 3.d0
        end if
!
!       Solid part
        solid_weight = h * nu_i / kB / T
        solid_weight = solid_weight / (dexp(solid_weight) -1.d0) - dlog(1.d0 - dexp(-solid_weight))
!       At zero frequency the solid DoS must be zero, otherwise integral may diverge
        if( i == 1 )then
          solid_weight(1) = 0.d0
        end if
!       Add solid part if the remaining DoS is positive
        solid_integrand = conv1 * twobykT * dos_subinterval - gas_integrand
        do j2=1,3
          if(solid_integrand(j2) > 0.d0)then
            entropy_group(j,k) = entropy_group(j,k) + solid_integrand(j2) * solid_weight(j2) * quadrature(j2)
          end if
        end do
!       Calculate 1PT part
        entropy1PT_group(j,k) = entropy1PT_group(j,k) + &
                                dot_product(conv1*twobykT*dos_subinterval*solid_weight,quadrature)
      end do
      entropy_group(j,k) = entropy_group(j,k) * kB / tau
      entropy1PT_group(j,k) = entropy1PT_group(j,k) * kB / tau
    end do
!   Entropy of mixture
    if( smixture == "vol")then
      temp(1) = kB * dfloat(1) * dlog(V / volume_group(j))
    else if( smixture == "mol")then
      temp(1) = kB * dfloat(1) * dlog(dfloat(ngroups) / 1.d0)
    end if
    entropy_group(j,1) = entropy_group(j,1) + temp(1)
    entropy1PT_group(j,1) = entropy1PT_group(j,1) + temp(1)
    temp(5) = entropy_group(j,1) + entropy_group(j,2) + entropy_group(j,3)
    write(10,'(I8,2X,ES15.8,2X,ES15.8,2X,ES15.8,2X,ES15.8,2X,ES15.8)') j, entropy_group(j,1:3), temp(5), temp(5) * eVtoJ
    temp(5) = entropy1PT_group(j,1) + entropy1PT_group(j,2) + entropy1PT_group(j,3)
    write(20,'(I8,2X,ES15.8,2X,ES15.8,2X,ES15.8,2X,ES15.8,2X,ES15.8)') j, entropy1PT_group(j,1:3), temp(5), temp(5) * eVtoJ
  end do
! For supergroups
  if(calc_supergroups)then
    allocate( entropy_supergroup(1:nsupergroups, 1:3) )
    allocate( entropy_supergroup_gas(1:nsupergroups, 1:3) )
    allocate( entropy_supergroup_solid(1:nsupergroups, 1:3) )
    allocate( degf_supergroup_gas(1:nsupergroups, 1:3) )
    allocate( degf_supergroup_solid(1:nsupergroups, 1:3) )
    allocate( entropy1PT_supergroup(1:nsupergroups, 1:3) )
    entropy_supergroup = 0.d0
    entropy_supergroup_gas = 0.d0
    entropy_supergroup_solid = 0.d0
    degf_supergroup_gas = 0.d0
    degf_supergroup_solid = 0.d0
    entropy1PT_supergroup = 0.d0
    open(unit=30, file="entropy_gas", status="unknown")
    open(unit=40, file="entropy_solid", status="unknown")
    write(10,*)
    write(20,*)
    write(30,*)
    write(40,*)
    write(10,*)"# Supergroup; Entropies in eV/K [trans, rot, vib]; Total entropy in eV/K; Total entropy in J/K"
    write(20,*)"# Supergroup; Entropies in eV/K [trans, rot, vib]; Total entropy in eV/K; Total entropy in J/K"
    write(30,*)"# Supergroup; Entropies in eV/K [trans, rot, vib]; number of DoF [trans, rot, vib]"
    write(40,*)"# Supergroup; Entropies in eV/K [trans, rot, vib]; number of DoF [trans, rot, vib]"
    do j=1,nsupergroups
      do k=1,3
!       The integral is performed using a quadrature rule for each subinterval (no. of subintervals is
!       no. of points minus one)
        do i=1,(n+1)/2-1
!         Update progress bar
          k2 = k2 + 1
          if(float(k2) > float(k3)*update_bar)then
            write(*,'(A)', advance='no')'='
            k3 = k3 + 1
          end if
!
!         DoS and frequency in the subinterval. Quadrature coefficients
          dos_subinterval = (/ Ssupergroup(j,i,k), (Ssupergroup(j,i,k) + Ssupergroup(j,i+1,k))/2.d0, Ssupergroup(j,i+1,k) /)
          nu_i = (/ dfloat(i-1) / tau, (dfloat(i-1) + 0.5d0) / tau, dfloat(i) / tau /)
          quadrature = (/ 1.d0/6.d0, 4.d0/6.d0, 1.d0/6.d0 /)
!
!         We skip the gas part if the no. of degrees of freedom is zero or the fluidicity is zero (i.e., numerically very small)
!         SHS must also be positive
          if(degf_supergroup(j,k) < 1.d-10 .or. f_supergroup(j,k) < 1.d-10 .or. &
             (k /= 2 .and. SHSbykB_supergroup(j,k) < 0.d0) .or. (k == 2 .and. SRbykB_supergroup(j) < 0.d0) )then
            gas_integrand = 0.d0
          else if(k == 3 .and. .not. vibrational_gas)then
            gas_integrand = 0.d0
          else
            temp(1) = degf_supergroup(j,k)
            gas_integrand = pi * conv1 * twobykT * Ssupergroup(j,1,k) * nu_i / 2.d0 / f_supergroup(j,k) / temp(1)
            gas_integrand = conv1 * twobykT * Ssupergroup(j,1,k) / (1.d0 + gas_integrand**2.d0)
!           The gas DoS cannot be larger than the total DoS
            do j2=1,3
              if(gas_integrand(j2) > conv1 * twobykT * dos_subinterval(j2))then
                gas_integrand(j2) = conv1 * twobykT * dos_subinterval(j2)
              end if
            end do
          end if
!         Calculate gas entropy
          if(k == 1 .or. k == 3)then
            entropy_supergroup_gas(j,k) = entropy_supergroup_gas(j,k) + dot_product(gas_integrand,quadrature) * &
                                          SHSbykB_supergroup(j,k) / 3.d0
            entropy_supergroup(j,k) = entropy_supergroup(j,k) + dot_product(gas_integrand,quadrature) * &
                                      SHSbykB_supergroup(j,k) / 3.d0
            degf_supergroup_gas(j,k) = degf_supergroup_gas(j,k) + dot_product(gas_integrand,quadrature)
          else if(k == 2)then
            entropy_supergroup_gas(j,k) = entropy_supergroup_gas(j,k) + dot_product(gas_integrand,quadrature) * &
                                          SRbykB_supergroup(j) / 3.d0
            entropy_supergroup(j,k) = entropy_supergroup(j,k) + dot_product(gas_integrand,quadrature) * &
                                      SRbykB_supergroup(j) / 3.d0
            degf_supergroup_gas(j,k) = degf_supergroup_gas(j,k) + dot_product(gas_integrand,quadrature)
          end if
!
!         Solid part
          solid_weight = h * nu_i / kB / T
          solid_weight = solid_weight / (dexp(solid_weight) -1.d0) - dlog(1.d0 - dexp(-solid_weight))
!         At zero frequency the solid DoS must be zero, otherwise integral may diverge
          if( i == 1 )then
            solid_weight(1) = 0.d0
          end if
!         Add solid part if the remaining DoS is positive
          solid_integrand = conv1 * twobykT * dos_subinterval - gas_integrand
          do j2=1,3
            if(solid_integrand(j2) > 0.d0)then
              entropy_supergroup_solid(j,k) = entropy_supergroup_solid(j,k) + &
                                              solid_integrand(j2) * solid_weight(j2) * quadrature(j2)
              entropy_supergroup(j,k) = entropy_supergroup(j,k) + solid_integrand(j2) * solid_weight(j2) * quadrature(j2)
              degf_supergroup_solid(j,k) = degf_supergroup_solid(j,k) + solid_integrand(j2) * quadrature(j2)
            end if
          end do
!         Calculate 1PT part
          entropy1PT_supergroup(j,k) = entropy1PT_supergroup(j,k) + &
                                       dot_product(conv1*twobykT*dos_subinterval*solid_weight,quadrature)
        end do
        entropy_supergroup(j,k) = entropy_supergroup(j,k) * kB / tau
        entropy_supergroup_gas(j,k) = entropy_supergroup_gas(j,k) * kB / tau
        entropy_supergroup_solid(j,k) = entropy_supergroup_solid(j,k) * kB / tau
        entropy1PT_supergroup(j,k) = entropy1PT_supergroup(j,k) * kB / tau
        degf_supergroup_gas(j,k) = degf_supergroup_gas(j,k) / tau
        degf_supergroup_solid(j,k) = degf_supergroup_solid(j,k) / tau
      end do
      ! Entropy of mixture
      if( smixture == "vol")then
        temp(1) = kB * ngroups_in_supergroup_eff(j) * dlog(V / volume_supergroup(j))
      else if( smixture == "mol")then
        temp(1) = kB * ngroups_in_supergroup_eff(j) * dlog(ngroups_eff / ngroups_in_supergroup_eff(j))
!     This entropy of mixing should only be used for binary mixtures !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! UNTESTED!!!!!
      else if( nsupergroups == 2 .and. smixture == "bin")then
        v_dash = volume_supergroup(2) / volume_supergroup(1) * ngroups_in_supergroup_eff(1) / ngroups_in_supergroup_eff(2)
        m_frac = ngroups_in_supergroup_eff(1) / ngroups_eff
        alpha = (dabs(v_dash - 1.d0) * (1.d0 - 2.d0 * m_frac) + v_dash + 1.d0) / 2.d0
        v_dash = 1.d0 / v_dash
        beta = (dabs(v_dash - 1.d0) * (1.d0 - 2.d0 * m_frac) + v_dash + 1.d0) / 2.d0
        if(j == 1)then
          temp(1) = kB * ngroups_in_supergroup_eff(1) * dlog( (ngroups_in_supergroup_eff(1) + &
                    ngroups_in_supergroup_eff(2)*alpha) / ngroups_in_supergroup_eff(1) )
        else if(j == 2)then
          temp(1) = kB * ngroups_in_supergroup_eff(2) * dlog( (ngroups_in_supergroup_eff(2) + &
                    ngroups_in_supergroup_eff(1)*beta) / ngroups_in_supergroup_eff(2) )
        end if
      end if
      entropy_supergroup(j,1) = entropy_supergroup(j,1) + temp(1)
      entropy_supergroup_gas(j,1) = entropy_supergroup_gas(j,1) + temp(1)
      entropy1PT_supergroup(j,1) = entropy1PT_supergroup(j,1) + temp(1)
      temp(5) = entropy_supergroup(j,1) + entropy_supergroup(j,2) + entropy_supergroup(j,3)
      write(10,'(I8,2X,ES15.8,2X,ES15.8,2X,ES15.8,2X,ES15.8,2X,ES15.8)') j, entropy_supergroup(j,1:3), temp(5), temp(5) * eVtoJ
      temp(5) = entropy1PT_supergroup(j,1) + entropy1PT_supergroup(j,2) + entropy1PT_supergroup(j,3)
      write(20,'(I8,2X,ES15.8,2X,ES15.8,2X,ES15.8,2X,ES15.8,2X,ES15.8)') j, entropy1PT_supergroup(j,1:3), temp(5), temp(5) * eVtoJ
!     Write out gas entropy and gas DoF
      write(30,'(I8,2X,ES15.8,2X,ES15.8,2X,ES15.8,2X,ES15.8,2X,ES15.8,2X,ES15.8)') j, &
            entropy_supergroup_gas(j,1:3), degf_supergroup_gas(j,1:3)
!     Write out solid entrop and solid DoF
      write(40,'(I8,2X,ES15.8,2X,ES15.8,2X,ES15.8,2X,ES15.8,2X,ES15.8,2X,ES15.8)') j, &
            entropy_supergroup_solid(j,1:3), degf_supergroup_solid(j,1:3)
    end do
    close(30)
    close(40)
!   If degrees of freedom adjustment is enabled, output the info
!   Write to file
    if(renormalize)then
      write(10,*)
      write(20,*)
      write(10,*)"# Entropy results with DoS normalized to the expected degrees of freedom"
      write(20,*)"# Entropy results with DoS normalized to the expected degrees of freedom"
      write(10,*)"# Supergroup; Entropies in eV/K [trans, rot, vib]; Total entropy in eV/K; Total entropy in J/K"
      write(20,*)"# Supergroup; Entropies in eV/K [trans, rot, vib]; Total entropy in eV/K; Total entropy in J/K"
!     Calculate expected number of DoF
      do j = 1, nsupergroups
        temp(1:3) = 0.d0
        do j2 = 1, ngroups_in_supergroup(j)
!         3 translational DoF per group in this supergroup
          temp(1) = temp(1) + weight_group(group_in_supergroup(j,j2))*3.d0
!         Total number of DoF for this group (takes constrains into account)
          temp(4) = weight_group(group_in_supergroup(j,j2))* &
                    dot_product(degf_group(group_in_supergroup(j,j2),1:3), (/1.d0, 1.d0, 1.d0/))
!         For diatomic molecules
          if(natoms_in_group(group_in_supergroup(j,j2)) == 2)then
!           We add 2 rotational DoF per linear molecule
            temp(2) = temp(2) + weight_group(group_in_supergroup(j,j2))*2.d0
!           We add the rest as vibrational DoF (if the bond is constrained this is correctly zero)
            temp(3) = temp(3) + temp(4) - weight_group(group_in_supergroup(j,j2))*5.d0
!         For molecules with 3 or more atoms
          else if(natoms_in_group(group_in_supergroup(j,j2)) > 2)then
            temp(2) = temp(2) + weight_group(group_in_supergroup(j,j2))*3.d0
            temp(3) = temp(3) + temp(4) - weight_group(group_in_supergroup(j,j2))*6.d0
          end if
        end do
!       Make sure we don't divide by zero if there are e.g. no vibrational and/or rotational DoF
!       Also if numerically the no. of DoF is negative, this will set entropy to zero
        do j2 = 1, 3
          if(degf_supergroup(j,j2) > 1.d-10)then
            temp(j2) = temp(j2) / degf_supergroup(j,j2)
          else
            temp(j2) = 0.d0
          end if
        end do
        temp(5) = temp(1)*entropy_supergroup(j,1) + temp(2)*entropy_supergroup(j,2) + temp(3)*entropy_supergroup(j,3)
        write(10,'(I8,2X,ES15.8,2X,ES15.8,2X,ES15.8,2X,ES15.8,2X,ES15.8)') j, temp(1:3)*entropy_supergroup(j,1:3), &
                                                                           temp(5), temp(5) * eVtoJ
        temp(5) = temp(1)*entropy1PT_supergroup(j,1) + temp(2)*entropy1PT_supergroup(j,2) + temp(3)*entropy1PT_supergroup(j,3)
        write(20,'(I8,2X,ES15.8,2X,ES15.8,2X,ES15.8,2X,ES15.8,2X,ES15.8)') j, temp(1:3)*entropy1PT_supergroup(j,1:3), &
                                                                           temp(5), temp(5) * eVtoJ
      end do
    end if
  end if
  close(10)
  close(20)
  write(*,'(A4)')'=] |'
  write(*,*)'                                       |'
  write(*,*)'.......................................|'
!******************************************************











end program















subroutine interpolate_1D(r2, r1, n2, n1, tau2, tau1, f)

! This subroutine linearly interpolates 1-D array r2 on the mesh given by array 1-D array r1
! assuming that the first element of both arrays is given at zero
! tau1 and tau2 are the size of the domains, that is tau1 = dt1*n1 and tau2 = dt2*n2, where
! dt is the spacing between data points and n is the number of data points.
! That is, data point n is given at tau = dt*(n-1)
! f is a normalizing factor, where data points in r2 are 1/f times too large. The interpolated
! data will be returned as f*r2 on the r1 mesh.
! Array r1 will be overwritten
  implicit none

  real*8 :: r1(:), r2(:), tau1, tau2, dt1, dt2, f
  integer :: n1, n2, t1, t2

  dt1 = tau1 / dfloat(n1)
  dt2 = tau2 / dfloat(n2)
  r1 = 0.d0

  do t1 = 0, n1-1
!   For each t1, we need to identify the closest t2 on the left and the right
!   This one is always on the left or at the same point
    t2 = int( dfloat(t1) * dt1 / dt2 )
!   t2 cannot be larger than n2-2
    if( t2 <= n2-2 )then
!     If points coincide exactly (including at zero), do not interpolate. We allow for rounding tolerance
       if( dabs(t2*dt2-t1*dt1) < 1.d-6 )then
        r1(t1+1) = r2(t2+1)
!     Otherwise, inter-/extrapolate from r2 to r1
      else
        r1(t1+1) = r2(t2+1) + (r2(t2+2) - r2(t2+1))/dt2 * (dfloat(t1)*dt1 - dfloat(t2)*dt2)
      end if
    end if
  end do

  r1 = r1*f

end subroutine













subroutine get_mass(mass_label, mass, nmasses, mass_types, mass_types_values)

  implicit none

  integer :: nmasses, i
  character*16 :: mass_label
  real*8 :: mass
  real*8 :: mass_types_values(:)
  character*16 :: mass_types(:)

  do i=1,nmasses
    if(trim(mass_label) == trim(mass_types(i)))then
      mass = mass_types_values(i)
      return
    end if
  end do

  write(*,*) 'ERROR: Mass for atom ', trim(mass_label), ' not found in the "masses" file!'

end subroutine





subroutine get_decomposed_velocities(m, r, v, natoms_in_group, v_rot, v_cm, v_vib, mass, w_group, &
                                     rxv_group, dist_group, eig)

  implicit none

  integer :: natoms_in_group, i, j, k
  real*8 :: m(:), r(:,:), v(:,:), mass
  real*8 :: v_rot(:,:), v_cm(:), v_vib(:,:), temp, L(1:3) = 0.d0
  real*8 :: Ine(1:3,1:3) = 0.d0, I_inv(1:3,1:3), CM(1:3) = 0.d0
  real*8 :: delta(1:3,1:3) = 0.d0, w(1:3) = 0.d0, eig(:)
  real*8 :: temp1, temp2
  real*8 :: w_group(:), rxv_group(:,:), dist_group(:)

  v_cm = 0.d0
  v_rot = 0.d0
  v_vib = 0.d0

! Kronecker delta
  delta(1,1) = 1.d0
  delta(2,2) = 1.d0
  delta(3,3) = 1.d0

! Get center of mass
  mass = 0.d0
  CM = 0.d0
  do k = 1, natoms_in_group
    mass = mass + m(k)
    do i = 1, 3
      CM(i) = CM(i) + r(k, i) * m(k)
    end do
  end do
  CM = CM / mass

! Get center of mass velocity (same for every atom)
  do k = 1, natoms_in_group
    do i = 1, 3
      v_cm(i) = v_cm(i) + v(k, i) * m(k)
    end do
  end do
  v_cm = v_cm / mass

! If we only have one atom in the system everything is translational.
! Set the values and back to main program (otherwise inverse of inertia will yield
! division by 0 error)
  if(natoms_in_group == 1)then
    v_rot = 0.d0
    v_vib = 0.d0
    w_group = 0.d0
    rxv_group = 0.d0
    dist_group = 1.d0
    return
  end if

! Refer everything to inertial frame where the center of mass is at rest
  do k = 1, natoms_in_group
    do i = 1, 3
      v(k,i) = v(k,i) - v_cm(i)
      r(k,i) = r(k,i) - CM(i)
    end do
    dist_group(k) = dsqrt(r(k,1)**2 + r(k,2)**2 + r(k,3)**2)
!   Readjust dist if it's too close to zero to prevent division by zero
    if(dist_group(k) < 1.d-8)then
      dist_group(k) = 1.d-8
    end if
  end do

! If we only have two atoms in the system, there is one vibrational degree of
! freedom and two rotational degrees of freedom, besides the three translational
! degrees of freedom. The moment of inertia matrix does not have an inverse in
! such a case and we need to do something else. Here we think of a basis where
! u1 is a unit vector along the direction connecting both atoms, and u2, u3
! lie in a plane perpendicular to u1. Vibrational velocities are the u1 components
! and rotational velocities are given as the instantaenous angular velocities
! around u2 and u3
  if(natoms_in_group == 2)then
!   w is the same if we use r1, v1 or r2, v2 because we're in the CM frame
    call cross_product(r(1,1:3), v(1,1:3), w(1:3))
    w(1:3) = w(1:3) / dot_product(r(1,1:3),r(1,1:3))
!   Get v_rot for each atom as v_rot = w x r
    call cross_product(w(1:3), r(1,1:3), v_rot(1,1:3))
    call cross_product(w(1:3), r(2,1:3), v_rot(2,1:3))
!   Get r x v
    call cross_product(r(1,1:3), v_rot(1,1:3), rxv_group(1,1:3))
    call cross_product(r(2,1:3), v_rot(2,1:3), rxv_group(2,1:3))
!   Vibrational velocity
    v_vib(1,1:3) = v(1,1:3) - v_rot(1,1:3)
    v_vib(2,1:3) = v(2,1:3) - v_rot(2,1:3)
    return
  end if

! Calculate tensor of inertia
  Ine=0.d0
  do k = 1, natoms_in_group
    do i = 1, 3
      do j = 1, 3
        Ine(i,j) = Ine(i,j) + m(k) * ( delta(i,j)*( r(k,1)**2 + r(k,2)**2 + r(k,3)**2 ) &
                   - r(k,i)*r(k,j) )
      end do
    end do
  end do
! And inverse of the tensor of inertia which I precomputed analytically to save time
  temp = - Ine(1,3)**2*Ine(2,2) + 2.d0*Ine(1,2)*Ine(1,3)*Ine(2,3) - Ine(1,1)*Ine(2,3)**2 - Ine(1,2)**2*Ine(3,3) &
         + Ine(1,1)*Ine(2,2)*Ine(3,3)
  I_inv(1,1) = (- Ine(2,3)**2 + Ine(2,2)*Ine(3,3) ) / temp
  I_inv(1,2) = (- Ine(1,2)*Ine(3,3) + Ine(1,3)*Ine(2,3) ) / temp
  I_inv(1,3) = (- Ine(1,3)*Ine(2,2) + Ine(1,2)*Ine(2,3) ) / temp
  I_inv(2,1) = I_inv(1,2)
  I_inv(2,2) = (- Ine(1,3)**2 + Ine(1,1)*Ine(3,3) ) / temp
  I_inv(2,3) = (- Ine(2,3)*Ine(1,1) + Ine(1,2)*Ine(1,3) ) / temp
  I_inv(3,1) = I_inv(1,3)
  I_inv(3,2) = I_inv(2,3)
  I_inv(3,3) = (- Ine(1,2)**2 + Ine(1,1)*Ine(2,2) ) / temp

! Calculate principal moments of inertia
  eig = 0.d0
  call get_eigenvalues_3x3(Ine,eig)

! Angular momentum (velocities and position wrt CM)
  L=0.d0
  rxv_group = 0.d0
  do k = 1, natoms_in_group
    rxv_group(k,1) = r(k,2)*v(k,3) - r(k,3)*v(k,2)
    rxv_group(k,2) = - r(k,1)*v(k,3) + r(k,3)*v(k,1)
    rxv_group(k,3) = r(k,1)*v(k,2) - r(k,2)*v(k,1)
    L(1) = L(1) + m(k) * rxv_group(k,1)
    L(2) = L(2) + m(k) * rxv_group(k,2)
    L(3) = L(3) + m(k) * rxv_group(k,3)
  end do

! Angular velocity
  w = 0.d0
  do i = 1, 3
    do j = 1, 3
      w(i) = w(i) + I_inv(i,j) * L(j)
    end do
  end do
  w_group = w

! Rotational velocity of each atom
  do k = 1, natoms_in_group
    v_rot(k,1) = w(2)*r(k,3) - w(3)*r(k,2)
    v_rot(k,2) = - w(1)*r(k,3) + w(3)*r(k,1)
    v_rot(k,3) = w(1)*r(k,2) - w(2)*r(k,1)
  end do

! Rewrite rxv using only rotational velocity now
  rxv_group = 0.d0
  do k = 1, natoms_in_group
    rxv_group(k,1) = r(k,2)*v_rot(k,3) - r(k,3)*v_rot(k,2)
    rxv_group(k,2) = - r(k,1)*v_rot(k,3) + r(k,3)*v_rot(k,1)
    rxv_group(k,3) = r(k,1)*v_rot(k,2) - r(k,2)*v_rot(k,1)
  end do

! Vibrational velocity of each atom
  do k = 1, natoms_in_group
    do i = 1, 3
      v_vib(k,i) = v(k,i) - v_rot(k,i)
    end do
  end do

end subroutine





subroutine restore_replica(L, r, natoms_in_group)

  implicit none

  real*8 :: r(:,:)
  real*8, allocatable :: new_r(:,:)
  real*8 :: L(1:3), d
  integer :: i, j, j2, natoms_in_group

! Brings r(j,i) j>1 to the same set of periodic replicas
! that r(1,i) belongs to, according to a minimum image convention 

  allocate( new_r(1:natoms_in_group,1:3) )

  do j=2,natoms_in_group
  do i=1,3
    d=1.d10
    do j2=-1,1
      if(dabs(r(j,i) + L(i)*dfloat(j2) - r(1,i)) < d)then
        d = dabs(r(j,i) + L(i)*dfloat(j2) - r(1,i))
        new_r(j,i) = r(j,i) + L(i)*dfloat(j2)
      end if
    end do
    r(j,i) = new_r(j,i)
  end do
  end do

end subroutine





subroutine get_delta(s0, T, conv2, nparticles, degf_group, mass_group, V, delta, &
                     omega)

  implicit none

  real*8 :: delta, s0, conv2, mass_group, T, V, degf_group
  real*8 :: nparticles
  real*8 :: pi, kB, omega, n1, n2

  pi = dacos(-1.d0)
! Boltzmann's constant in eV/K
  kB = 8.6173324d-5

  n1 = degf_group
  n2 = nparticles

! Return delta
  if( degf_group < 1.d-5 )then
    delta = 0.d0
  else
    delta = conv2 * 2.d0 * s0 / 3.d0 / n1 * dsqrt(pi * kB * T / mass_group) &
                    * (n2 / V)**(1.d0/3.d0) * (6.d0 / pi)**(2.d0/3.d0) * omega**(1.d0/3.d0)
  end if
end subroutine





!subroutine get_omega(this_group, ngroups, natoms_in_group, volume_group, mass_group, omega)
!
!  implicit none
!
!  real*8 :: omega
!  real*8, allocatable :: volume_group(:), mass_group(:)
!  integer :: ngroups, j, this_group
!  integer, allocatable :: natoms_in_group(:)
!
!  omega = 0.d0
!
!  do j=1, ngroups
!    omega = omega + 1.d0/4.d0/dsqrt(2.d0) * natoms_in_group(j) / natoms_in_group(this_group) * &
!            (1.d0 + (volume_group(j)/volume_group(this_group))**(1.d0/3.d0))**2.d0 * &
!            dsqrt(1.d0 + mass_group(this_group)/mass_group(j))
!  end do
!
!end subroutine











subroutine find_f(delta_in, f_out)

  implicit none

  integer :: i, limit=1000
  real*8 :: tol, p, f, f1, f2, delta, root, mid, delta_in, f_out

! This function must be defined before doing anything else
  p(f) =  2.d0 * delta**(-9./2.) * f**(15./2.) &
         -6.d0 * delta**(-3.)    * f**(5.) &
         -1.d0 * delta**(-3./2.) * f**(7./2.) &
         +6.d0 * delta**(-3./2.) * f**(5./2.) &
         +2.d0 * f -2.d0

! If delta is negative (or zero) set fluidicity to zero and get out
  if(delta_in <= 0.)then
    f_out = 0.d0
    return
  end if

  delta = delta_in

! f must belong to the [0,1] interval so the safest guesses
! are 0 and 1
  f1 = 0.d0 ; f2 = 1.d0 ; tol = 1.d-15

! p(f1) must be negative and p(f2) must be positive
  if( p(f1) > 0.d0 )then
    mid = f1
    f1 = f2
    f2 = mid
  end if

  i = 1
  do while (i < limit)
    mid = (f2 + f1)/2.d0
    if( dabs(p(mid)) < tol )then
      exit
    end if
    if( p(mid) < 0.d0 )then
      f1 = mid
    else
      f2 = mid
    end if
    i = i+1
  end do

  f_out = mid

end subroutine




subroutine cross_product(u,v,cross)

  implicit none

  real*8 :: u(1:3), v(1:3), cross(1:3)

  cross(1:3) = (/ u(2)*v(3)-u(3)*v(2), &
                 -u(1)*v(3)+u(3)*v(1), &
                  u(1)*v(2)-u(2)*v(1) /)

end subroutine
