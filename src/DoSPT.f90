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
!                            DoSPT v0.1                              !
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
!                      Miguel A. Caro et al.                         !
!              [Implementation paper citation info]                  !
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
!!!          Distribution last updated on 09 Sep. 2016             !!!
!!!!!                                                            !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!******************************************************
! This block declares all the variables needed during
! program execution
  implicit none

  integer :: i, i2, natoms, j, k, j2, ios, k2, k3, k4, n_volumes, Ng, Nsg, niter
  real*8, allocatable :: wsave(:), r(:), rw(:)
  real*4, allocatable :: velocities(:,:,:), positions(:,:,:)
  real*8 :: temp(1:6), update_bar, res, D_rot_real, D_ideal, D_rot_ideal, sigma_eff
  real*8 ::  nu_i(3), gas_integrand(3), solid_integrand(3), dos_subinterval(3), solid_weight(3), quadrature(3)
  character :: crap
  real*8, allocatable :: Stotal(:,:), m(:)
  character*40 :: time

  integer :: n_particles
  real*8, allocatable :: mp(:), rp(:,:), vp(:,:)
  real*4, allocatable :: v_rot(:,:,:), v_cm(:,:,:), v_vib(:,:,:)
  real*4, allocatable :: w(:,:,:), rxv(:,:,:), volumes(:,:), dist(:,:)
  real*8, allocatable :: v_rot_group(:,:), v_cm_group(:), v_vib_group(:,:)
  real*8, allocatable :: w_group(:), rxv_group(:,:), volume_group(:), dist_group(:)
  real*8, allocatable :: x(:), y(:), z(:), volumes_temp(:)

! Variables for grouping interface
  integer :: ngroups, nparticles, nmasses
  integer, allocatable :: natoms_in_group(:), atoms_in_group(:,:)
  real*8, allocatable :: mass_types_values(:), mass_group(:), Sgroup(:,:,:), eig_group(:,:)
  character*16, allocatable :: mass_types(:)
  character*16 :: mass_label
  real*8, allocatable :: delta_group(:,:), degf_group(:,:), degf_supergroup(:,:)
  real*8, allocatable :: f_group(:,:), f_supergroup(:,:), y_group(:,:), z_group(:,:)
  real*8, allocatable :: SHSbykB_group(:,:), SHSbykB_supergroup(:,:), SRbykB_group(:), SRbykB_supergroup(:)
  real*8, allocatable :: entropy_group(:,:), entropy_supergroup(:,:), entropy_supergroup_gas(:,:)
  real*8, allocatable :: entropy_supergroup_solid(:,:), degf_supergroup_gas(:,:), degf_supergroup_solid(:,:)
  real*8, allocatable :: entropy1PT_group(:,:), entropy1PT_supergroup(:,:), eig_supergroup(:,:)
  real*8, allocatable :: symmetry_number_group(:), symmetry_number_supergroup(:)
  real*8 :: omega, conv4

! Supergroup interface
  integer :: nsupergroups = 0
  integer, allocatable :: ngroups_in_supergroup(:), group_in_supergroup(:,:), natoms_in_supergroup(:)
  logical :: calc_supergroups, error = .false.
  real*8, allocatable :: Ssupergroup(:,:,:), mass_supergroup(:), volume_supergroup(:)
  real*8, allocatable :: delta_supergroup(:,:), y_supergroup(:,:), z_supergroup(:,:)

! Filtering stuff must be single precision
  real*4, allocatable :: smooth_x(:), smooth_y(:), smooth_ys(:), smooth_rw(:), smooth_res(:)
  real*4 :: smooth_f

  real*8, allocatable :: sigma_group(:,:), sigma_supergroup(:,:)

! Binary entropy of mixing variables
  real*8 :: v_dash, alpha, beta, m_frac

! Variables for input parameters (read from "input" file):
! L(1:3) = lattice vectors (only orthrhombic cells currently supported)
! tau = trajectory length in ps
! n = number of snapshots in the trajectory
! T = temperature
! smooth = whether to filter DoS for individual groups
! sigma_nu = characteristic filtering length in ps^-1
  real*8 :: L(1:3) = 0., V = 0., tau = 0., T = 1.d-10, twobykT, sigma_nu = 0.d0
  integer :: n = 0, iostatus
  character*3 :: mode = "gro"
  character*64 :: keyword
  character*1 :: keyword_first
  logical :: smooth = .false., estimate_vel = .false., renormalize = .false., vibrational_gas = .false.
  logical :: f_opt = .true., f_rot_opt = .true., write_dos = .true.
  character*3 :: smixture = "mol", hs_formalism

! Constants and units
  real*8 :: pi, kB, h, ps, nm, eV, amu, conv1, conv2, conv3, eVtoJ
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

subroutine sort_supergroups(ngroups, nsupergroups, ngroups_in_supergroup, group_in_supergroup)
  integer, intent(in) :: ngroups
  integer, intent(out) :: nsupergroups
  integer, allocatable, intent(out) :: ngroups_in_supergroup(:), group_in_supergroup(:,:)
end subroutine

subroutine read_trajectory(n, natoms, tau, L, mode, positions, velocities, estimate_vel, error, &
                           m, nmasses, mass_types, mass_types_values, volumes, volumes_temp)
  integer, intent(in) :: n, natoms, nmasses
  real*8, intent(inout) :: m(:)
  real*8, intent(in)  :: mass_types_values(:)
  character*16, intent(in)  :: mass_types(:)
  real*8, intent(in) :: L(1:3), volumes_temp(:), tau
  real*4, intent(inout) :: positions(:,:,:), velocities(:,:,:), volumes(:,:)
  logical, intent(in) :: estimate_vel
  logical, intent(inout) :: error
  character*3, intent(in) :: mode
end subroutine

subroutine get_omega(this_group, ngroups, n_in_group, volume_group, sigma_group, mass_group, mode, omega)
  real*8, intent(out) :: omega
  real*8, intent(in) :: volume_group(:), sigma_group(:), mass_group(:)
  integer, intent(in) :: ngroups, this_group
  integer, intent(in) :: n_in_group(:)
  character*1, intent(in) :: mode
end subroutine

subroutine fluidicity_calculator(nsupergroups, ngroups_in_supergroup, N_DoF, s0, M, T, V, V_voro, niter, res, sigma, f)
  real*8, intent(in) :: T, V, s0(:), N_DoF(:), M(:), V_voro(:)
  integer, intent(in) :: nsupergroups, ngroups_in_supergroup(:)
  real*8, intent(out) :: res, f(:), sigma(:)
  integer, intent(out) :: niter
end subroutine

subroutine get_compressibility(nsupergroups, ngroups_in_supergroup, V, sigma, f, xi, z)
  integer, intent(in) :: nsupergroups, ngroups_in_supergroup(:)
  real*8, intent(in) :: V, sigma(:), f(:)
  real*8, intent(out) :: z, xi
end subroutine

subroutine get_real_rotational_diffusivity(tau, n, N_DoF, ngroups, group_in_supergroup, j2, w, D)
  integer, intent(in) :: n, ngroups, group_in_supergroup(:,:), j2
  real*8, intent(in) :: tau, N_DoF
  real*4, intent(in) :: w(:,:,:)
  real*8, intent(out) :: D
end subroutine

end interface
!******************************************************





!******************************************************
! In this block several constants and conversion factors
! and fundamental constants are defined
  pi = dacos(-1.d0)
! Boltzmann's constant in eV/K
  kB = 8.6173324d-5
! Planck's constant in eV*ps
  h = 4.135667662d-3
! Different constants
  ps = 1.d-12
  nm = 1.d-9
  eV = 1.602176565d-19
  amu = 1.660538921d-27
! To transform from eV to J
  eVtoJ = 96486.9d0
! Conversion factors
! Sum(s_j^k)(0) for each atom group comes in nm^2 / ps, we then
! multiply by amu and divide by eV to obtain s0 for each atom group.
! Therefore, the units of s0 are (nm^2 * amu) / (eV * ps), but they
! should be returned in ps. We define the conversion factor:
  conv1 = nm**2 * amu / eV / ps**2
! Delta should be dimensionless, but it's calculated as
! (ps * sqrt(eV / amu) / nm), hence we define the conversion factor:
  conv2 = ps * dsqrt(eV / amu) / nm
! The quantity inside the SHS logarithm will be calculated in
! (amu / eV)^3/2 * nm^3 / ps^3 but it should be dimensionless. We
! define this conversion factor:
  conv3 = (amu / eV)**(3.d0/2.d0) * nm**3.d0 / ps**3.d0
! T/SR will be calculated in (amu nm^2) / (eV ps^2) but should be
! dimensionless:
  conv4 = (amu * nm**2.d0) / (eV * ps**2.d0)
!******************************************************





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
    goto 98
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
      goto 98
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
      goto 98
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
              goto 98
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
              goto 98
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
      goto 98
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
! Only calculate Voronoi cell volumes every 100 time steps
  n_volumes = 1 + n/100
  allocate( volumes(1:natoms, 1:n_volumes) ) 
  allocate( volumes_temp(1:natoms) )
  allocate( x(1:natoms) )
  allocate( y(1:natoms) )
  allocate( z(1:natoms) )
! Principal moments of inertia
  allocate( eig_group(1:ngroups, 1:3) )
  eig_group = 0.d0
! These for groups
  allocate( w(1:ngroups, 1:3, 1:n) )
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
  volume_group = 0.d0
  sigma_group = 0.d0
  entropy_group = 0.d0
  entropy1PT_group = 0.d0
!******************************************************





!******************************************************
! Allocate arrays for trajectory input and FFT routine
  allocate( wsave(1:2*n+15) )
  allocate( r(1:n) )
!  allocate( rw(1:n) )
  allocate( velocities(1:natoms, 1:3, 1:n) )
  allocate( positions(1:natoms, 1:3, 1:n) )
  allocate( Stotal(1:(n+1)/2, 1:3) )
  allocate( m(1:natoms) )
! Initialize wsave for FFT routine
  call dffti(n, wsave)
!******************************************************





!******************************************************
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
! Read trajectory file
  call read_trajectory(n, natoms, tau, L, mode, positions, velocities, estimate_vel, error, &
                       m, nmasses, mass_types, mass_types_values, volumes, volumes_temp)
  if(error)then
    goto 98
  end if
!******************************************************





!******************************************************
! Calculate the trajectory-averaged Voronoi cell volumes
! for each group
  do j=1,ngroups
    do j2=1,natoms_in_group(j)
      do i=1,n_volumes
        volume_group(j) = volume_group(j) + volumes(atoms_in_group(j,j2),i)
      end do
    end do
    volume_group(j) = volume_group(j) / dfloat(n_volumes)
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
!     Read in molecule info
      do j2=1,natoms_in_group(j)
        do k=1,3
          rp(j2,k) = positions(atoms_in_group(j,j2), k, i)
          vp(j2,k) = velocities(atoms_in_group(j,j2), k, i)
        end do
      end do
!     Get decomposed values for present molecule
      call restore_replica(L, rp, natoms_in_group(j))
      call get_decomposed_velocities(mp, rp, vp, natoms_in_group(j), v_rot_group, v_cm_group, v_vib_group, mass_group(j), &
                                     w_group, rxv_group, dist_group, temp(1:3))
      w(j, 1:3, i) = w_group(1:3)
!     Calculate average principal moments of inertia (amu * nm^2)
      eig_group(j,1:3) = eig_group(j,1:3) + temp(1:3)/dfloat(n)
!     Pass values to atoms
      do j2=1,natoms_in_group(j)
        do k=1,3
          v_rot(atoms_in_group(j,j2), k, i) = v_rot_group(j2, k)
          v_cm(atoms_in_group(j,j2), k, i) = v_cm_group(k)
          v_vib(atoms_in_group(j,j2), k, i) = v_vib_group(j2, k)
          rxv(atoms_in_group(j,j2), k, i) = rxv_group(j2, k)
        end do
        dist(atoms_in_group(j,j2), i) = dist_group(j2)
      end do
    end do
!   Deallocate some variables
    deallocate(mp, rp, vp)
    deallocate(v_rot_group, v_cm_group, v_vib_group, w_group, rxv_group, dist_group)
  end do
  deallocate(positions)
  deallocate(velocities)
  write(*,'(A4)')'=] |'
  write(*,*)'                                       |'
  write(*,*)'.......................................|'
!******************************************************





!******************************************************
! Calculate density of states
! ATTENTION, this is a parallel block
  write(*,*)'                                       |'
  write(*,*)'Calculating density of states...       |'
  write(*,*)'                                       |'
! Check prime number factorization for n and print warning if high factors are found
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
  Sgroup = 0.d0
!$omp parallel do firstprivate(wsave) private(r,rw,i,j,j2,k,k2,smooth_y,smooth_ys,smooth_rw,smooth_res) shared(Sgroup,k3,k4)
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
        do i=1,n
!         Translational
          if(k2 == 1)then
            r(i) = v_cm(atoms_in_group(j,j2), k, i)
!         Rotational
          else if(k2 == 2)then
            r(i) = rxv(atoms_in_group(j,j2), k, i) / dist(atoms_in_group(j,j2), i)
!         Vibrational
          else if(k2 == 3)then
            r(i) = v_vib(atoms_in_group(j,j2), k, i)
          end if
        end do
!       Call FFT routine
        call dfftf(n, r, wsave)
!       Save DoS for this group
        Sgroup(j,1,k2) = Sgroup(j,1,k2) + m(atoms_in_group(j,j2)) * r(1)**2*tau/dfloat(n)**2
        do i=2,(n+1)/2
          Sgroup(j,i,k2) = Sgroup(j,i,k2) + m(atoms_in_group(j,j2)) * (r(2*i-2)**2+r(2*i-1)**2)*tau/dfloat(n)**2
        end do
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
      Stotal(1,k2) = Stotal(1,k2) + Sgroup(j,1,k2)
      do i=2,(n+1)/2
        Stotal(i,k2) = Stotal(i,k2) + Sgroup(j,i,k2)
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
    allocate( eig_supergroup(1:nsupergroups,1:3) )
    sigma_supergroup = 0.d0
    Ssupergroup = 0.d0
    mass_supergroup = 0.d0
    symmetry_number_supergroup = 0.d0
    volume_supergroup = 0.d0
    natoms_in_supergroup = 0
    eig_supergroup = 0.d0
    do j=1,nsupergroups
      do j2=1,ngroups_in_supergroup(j)
        eig_supergroup(j,1:3) = eig_supergroup(j,1:3) + &
                                 eig_group(group_in_supergroup(j,j2),1:3)/dfloat(ngroups_in_supergroup(j))
        symmetry_number_supergroup(j) = symmetry_number_supergroup(j) + &
                                        symmetry_number_group(group_in_supergroup(j,j2))/dfloat(ngroups_in_supergroup(j))
        mass_supergroup(j) = mass_supergroup(j) + mass_group(group_in_supergroup(j,j2))
        volume_supergroup(j) = volume_supergroup(j) + volume_group(group_in_supergroup(j,j2))
        natoms_in_supergroup(j) = natoms_in_supergroup(j) + natoms_in_group(group_in_supergroup(j,j2))
        do k2=1,3
          Ssupergroup(j,1,k2) = Ssupergroup(j,1,k2) + Sgroup(group_in_supergroup(j,j2),1,k2)
          do i=2,(n+1)/2
            Ssupergroup(j,i,k2) = Ssupergroup(j,i,k2) + Sgroup(group_in_supergroup(j,j2),i,k2)
          end do
        end do
      end do
!     Mass of the supergroup is average mass of its groups. The groups should be equivalent
!     (e.g., all water molecules) and so this should not be necessary, but it's done for
!     consistency in case the groups are not equivalent (for some reason)
      mass_supergroup(j) = mass_supergroup(j) / dfloat(ngroups_in_supergroup(j))
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
  update_bar = dfloat(ngroups*(n+1)/2)/36.d0
  k2 = 1
!
  write(10,*) "# Total density of states for all the groups involved"
  do i=1,(n+1)/2
    write(10, *) dfloat(i-1)/tau, conv1 * twobykT * Stotal(i,1), conv1 * twobykT * Stotal(i,2), conv1 * twobykT * Stotal(i,3)
  end do
  do j=1,ngroups
    write(10,*) 
    write(10,*) "# Density of states (in ps) for group", j
    do i=1,(n+1)/2
!     Update progress bar every ngroups*n/36 iterations
      if(dfloat((j-1)*(n+1)/2+i) > dfloat(k2)*update_bar)then
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
    update_bar = dfloat(nsupergroups*(n+1)/2)/36.d0
    k2 = 1
!
    do j=1,nsupergroups
      write(10,*) "# Density of states (in ps) for supergroup", j
      do i=1,(n+1)/2
!     Update progress bar every ngroups*n/36 iterations
        if(dfloat((j-1)*(n+1)/2+i) > dfloat(k2)*update_bar)then
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(1 == 1)then
!       This menthod of calculating the fluidicity is essentially equivalent to Lin's method,
!       except for the inclusion of the "omega term" heuristically in the definition of the
!       normalized diffusivity "delta". This method is now deprecated and the more sophisticated
!       method based on numerical minimization of a penalty function should be employed.
!       For a monocomponent system, the results should be identical.
!       Available for testing and debugging purposes.
!       Calculate Delta:
        call get_omega(j, ngroups, 1+0*natoms_in_group, volume_group, sigma_group(1:ngroups,k), mass_group, "v", omega)
        call get_delta(conv1 * twobykT * Sgroup(j,1,k), T, conv2, natoms_in_group(j), 1, &
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
      call fluidicity_calculator(ngroups, 1+0*natoms_in_group, degf_group(1:ngroups,k), conv1 * twobykT * Sgroup(1:ngroups,1,k), &
                                 mass_group, T, V, volume_group, niter, res, sigma_group(1:ngroups,k), f_group(1:ngroups,k))
      do j = 1, ngroups
!       Get partial compressibility and partial packing fraction
        call get_compressibility(1, 1+0*natoms_in_group(j:j), volume_group(j), sigma_group(j:j, k), &
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
            call get_omega(j, nsupergroups, ngroups_in_supergroup, volume_supergroup, &
                           sigma_supergroup(1:nsupergroups,k), mass_supergroup, "v", omega)
            temp(1) = V
          end if
          call get_delta(conv1 * twobykT * Ssupergroup(j,1,k), T, conv2, natoms_in_supergroup(j), ngroups_in_supergroup(j), &
                         degf_supergroup(j,k), mass_supergroup(j), temp(1), delta_supergroup(j,k), omega)
!         And f:
          call find_f(delta_supergroup(j,k), f_supergroup(j,k))
!         Packing fraction and compressibility:
          y_supergroup(j,k) = f_supergroup(j,k)**(5.d0/2.d0) / delta_supergroup(j,k)**(3.d0/2.d0)
          z_supergroup(j,k) = (1.d0 + y_supergroup(j,k) + y_supergroup(j,k)**2d0 - y_supergroup(j,k)**3.d0) / &
                              (1.d0 - y_supergroup(j,k))**3.d0
          sigma_supergroup(j,k) = (6.d0/pi * y_supergroup(j,k) / f_supergroup(j,k) * volume_supergroup(j) / &
                                  dfloat(ngroups_in_supergroup(j)) )**(1.d0/3.d0)
        end if
      end do
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Get fluidicity from calculator.
    if(f_opt)then
!     Translational part
      k = 1
      call fluidicity_calculator(nsupergroups, ngroups_in_supergroup, degf_supergroup(1:nsupergroups,k), &
                                 conv1 * twobykT * Ssupergroup(1:nsupergroups,1,k), mass_supergroup, T, V, &
                                 volume_supergroup, niter, res, sigma_supergroup(1:nsupergroups,k), &
                                 f_supergroup(1:nsupergroups,k))
      do j = 1, nsupergroups
!       Get partial compressibility and partial packing fraction
        call get_compressibility(1, ngroups_in_supergroup(j:j), volume_supergroup(j), sigma_supergroup(j:j, k), &
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
                                             group_in_supergroup, j, w, D_rot_real)
!        call get_diffusivity(Ssupergroup(j,1,k), T, mass_supergroup(j), degf_supergroup(j,k), D_rot_real)
!       Ideal translational diffusion coefficient -> D_ideal
        call get_omega(j, nsupergroups, ngroups_in_supergroup, volume_supergroup, &
                       sigma_supergroup(1:nsupergroups, k), mass_supergroup, "v", omega)
        call get_zero_pressure_diffusivity(sigma_supergroup(j,k), omega, T, mass_supergroup(j), &
                                           dfloat(ngroups_in_supergroup(j)), V, D_ideal)
!       Ideal rotational diffusion coefficient -> D_rot_ideal
        call get_ideal_rotational_diffusivity(sigma_supergroup(j,k), T, mass_supergroup(j), &
                                              dfloat(ngroups_in_supergroup(j)), volume_supergroup(j), D_ideal, D_rot_ideal)
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
        Nsg = ngroups_in_supergroup(j)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!temp(1) = V
temp(1) = volume_supergroup(j)
        SHSbykB_supergroup(j,k) = 5.d0/2.d0 &
           + dlog(conv3 * (2.d0 * pi * mass_supergroup(j) * kB * T / h**2.d0)**(3.d0/2.d0) * &
             temp(1) * z_supergroup(j,k) / dfloat(Nsg) / f_supergroup(j,k)) &
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
        temp(1) = kB * dfloat(ngroups_in_supergroup(j)) * dlog(V / volume_supergroup(j))
      else if( smixture == "mol")then
        temp(1) = kB * dfloat(ngroups_in_supergroup(j)) * dlog(dfloat(ngroups) / dfloat(ngroups_in_supergroup(j)))
!     This entropy of mixing should only be used for binary mixtures !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! UNTESTED!!!!!
      else if( nsupergroups == 2 .and. smixture == "bin")then
        v_dash = volume_supergroup(2) / volume_supergroup(1) * dfloat(ngroups_in_supergroup(1)) / dfloat(ngroups_in_supergroup(2))
        m_frac = dfloat(ngroups_in_supergroup(1)) / dfloat(ngroups)
        alpha = (dabs(v_dash - 1.d0) * (1.d0 - 2.d0 * m_frac) + v_dash + 1.d0) / 2.d0
        v_dash = 1.d0 / v_dash
        beta = (dabs(v_dash - 1.d0) * (1.d0 - 2.d0 * m_frac) + v_dash + 1.d0) / 2.d0
        if(j == 1)then
          temp(1) = kB * dfloat(ngroups_in_supergroup(1)) * dlog( (dfloat(ngroups_in_supergroup(1)) + &
                    dfloat(ngroups_in_supergroup(2))*alpha) / dfloat(ngroups_in_supergroup(1)) )
        else if(j == 2)then
          temp(1) = kB * dfloat(ngroups_in_supergroup(2)) * dlog( (dfloat(ngroups_in_supergroup(2)) + &
                    dfloat(ngroups_in_supergroup(1))*beta) / dfloat(ngroups_in_supergroup(2)) )
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
          temp(1) = temp(1) + 3.d0
!         Total number of DoF for this group (takes constrains into account)
          temp(4) = dot_product(degf_group(group_in_supergroup(j,j2),1:3), (/1.d0, 1.d0, 1.d0/))
!         For diatomic molecules
          if(natoms_in_group(group_in_supergroup(j,j2)) == 2)then
!           We add 2 rotational DoF per linear molecule
            temp(2) = temp(2) + 2.d0
!           We add the rest as vibrational DoF (if the bond is constrained this is correctly zero)
            temp(3) = temp(3) + temp(4) - 5.d0
!         For molecules with 3 or more atoms
          else if(natoms_in_group(group_in_supergroup(j,j2)) > 2)then
            temp(2) = temp(2) + 3.d0
            temp(3) = temp(3) + temp(4) - 6.d0
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





!******************************************************
! Print message at the end depending on execution
! success
97 write(*,*)'                                       |'
   write(*,*)'Program successfully executed          |'
   write(*,*)'                                       |'
   write(*,*)'End of execution:                      |'
   call timestring(time)
   write(*,'(1X,A39,A)')time,'|'
   goto 99
98 write(*,*)'                                       |'
   write(*,*)'Program could not be executed          |'
   write(*,*)'Check for errors                       |'
99 write(*,*)'_______________________________________/'
!******************************************************





end program












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





subroutine get_delta(s0, T, conv2, natoms, nparticles, degf_group, mass_group, V, delta, &
                     omega)

  implicit none

  real*8 :: delta, s0, conv2, mass_group, T, V, degf_group
  integer :: natoms, nparticles
  real*8 :: pi, kB, omega, n1, n2

  pi = dacos(-1.d0)
! Boltzmann's constant in eV/K
  kB = 8.6173324d-5

  n1 = degf_group
  n2 = dfloat(nparticles)

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
    do
      read(10, '(I8)', advance="no", iostat=iostatus) temp_int
      if(iostatus == 0)then
        ngroups_in_supergroup(i) = ngroups_in_supergroup(i) + 1
        group_in_supergroup(i, ngroups_in_supergroup(i)) = temp_int
      else if(iostatus == -2)then
        ngroups_in_supergroup(i) = ngroups_in_supergroup(i) + 1
        group_in_supergroup(i, ngroups_in_supergroup(i)) = temp_int
        exit
      else
        backspace(10)
        sg_string = ""
        do
          read(10,'(A)', advance="no", iostat=iostatus) sg_char
          if(sg_char == "-")then
            read(sg_string,*) group_start
            sg_string = ""
          else if(sg_string == "," .or. iostatus /= 0)then
            read(sg_string,*) group_end
            exit
          else
            write(sg_string,*) trim(sg_string) // sg_char
          end if
        end do
        do j=group_start,group_end
          ngroups_in_supergroup(i) = ngroups_in_supergroup(i) + 1
          group_in_supergroup(i, ngroups_in_supergroup(i)) = j
        end do
        exit
      end if
    end do
  end do
  rewind(10)
  close(10)

end subroutine





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
