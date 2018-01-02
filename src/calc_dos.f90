module calc_dos

! This module does the following things:
!
! * Allocate arrays, including smoothing
! * Calculate trajectory-averaged volumes
! * Carry out velocity decomposition
! * Calculate decomposed kinetic energy
! * Calculate density of states
! * Write density of states to file

  use constants
  use read_input
  use trajectory
  use topology
  use misc2

  implicit none

  integer :: k3, k4, n_prev, imax

  real*8, allocatable :: wsave(:), r(:), rw(:), r2(:)
  real*8 :: temp(1:6), update_bar, ngroups_eff = 0.d0
  real*8, allocatable :: Stotal(:,:)
  real*8, allocatable :: mp(:), rp(:,:), vp(:,:)
  real*8, allocatable :: v_rot_group(:,:), v_cm_group(:), v_vib_group(:,:)
  real*8, allocatable :: w_group(:), rxv_group(:,:), volume_group(:), dist_group(:)
  real*8, allocatable :: x(:), y(:), z(:), Sgroup(:,:,:)
  real*8, allocatable :: eig_group(:,:), symmetry_number_supergroup(:), mass_group(:), eig_supergroup(:,:)
  real*8, allocatable :: Ssupergroup(:,:,:), mass_supergroup(:), volume_supergroup(:)
  real*8, allocatable :: natoms_in_supergroup(:), ngroups_in_supergroup_eff(:)

  real*4, allocatable :: v_rot(:,:,:), v_cm(:,:,:), v_vib(:,:,:)
  real*4, allocatable :: w(:,:,:), eig_group_inst(:,:,:), rxv(:,:,:), dist(:,:), e_kin_tot(:), e_kin(:,:)
! Filtering stuff must be single precision
  real*4, allocatable :: smooth_x(:), smooth_y(:), smooth_ys(:), smooth_rw(:), smooth_res(:)
  real*4 :: smooth_f

! Volume exclusion
  real*8, allocatable :: V_apparent(:)
  real*8, allocatable :: N_apparent(:)


  contains


!=================================================================================================
!=================================================================================================
subroutine get_dos()
!******************************************************
! Allocate arrays for trajectory input and FFT routine
  allocate( wsave(1:2*n+15) )
  allocate( r(1:n) )
  allocate( r2(1:n) )
!  allocate( rw(1:n) )
  allocate( e_kin_tot(1:n) )
  allocate( Stotal(1:(n+1)/2, 1:3) )
! Initialize wsave for FFT routine
  call dffti(n, wsave)
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
  allocate( mass_group(1:ngroups) )
  allocate( Sgroup(1:ngroups,1:(n+1)/2,1:3) )
  mass_group = 0.d0
  volume_group = 0.d0
!******************************************************





!******************************************************
! FIX FILTERING TO WORK WITH TOPOLOGY RECONSTRUCTION. In the meantime print warning message.
  if( smooth .and. topology_has_changed )then
    write(*,*)'                                       |'
    write(*,*)'WARNING: usage of DoS filtering togeth-|  <-- WARNING'
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
    allocate( natoms_in_supergroup(1:nsupergroups) )
    allocate( ngroups_in_supergroup_eff(1:nsupergroups) )
    allocate( eig_supergroup(1:nsupergroups,1:3) )
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
!   If volume exclusion is enabled, we handle it here
    allocate( V_apparent(1:nsupergroups) )
    allocate( N_apparent(1:nsupergroups) )
    V_apparent = V
    N_apparent = ngroups_eff
    allocate( exclude_volume(1:nsupergroups) )
    exclude_volume = .false.
    if( exclude_volume_logical )then
      call exclude_volume_assign()
      do j = 1, nsupergroups
        if( exclude_volume(j) )then
          V_apparent = V_apparent - volume_supergroup(j)
          N_apparent = N_apparent - ngroups_in_supergroup_eff(j)
        end if
      end do
      do j = 1, nsupergroups
        if( exclude_volume(j) )then
          V_apparent(j) = volume_supergroup(j)
          N_apparent(j) = ngroups_in_supergroup_eff(j)
        end if
      end do
    end if
!   End of volume exclusion block
!
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
  if( print_mi )then
    write(*,*)'                                       |'
    write(*,*)'The average principal moments of       |'
    write(*,*)'inertia (in amu*nm^2) are:             |'
    write(*,*)'                                       |'
    write(*,*)'Supgrp        I_1        I_2        I_3|'
    write(*,*)'------ ---------- ---------- ----------|'
    do j = 1, nsupergroups
      write(*,'(1X,I6,1X,F10.4,1X,F10.4,1X,F10.4,A1)') j, eig_supergroup(j,1:3), '|'
    end do
    write(*,*)'                                       |'
    write(*,*)'.......................................|'
  end if
!******************************************************





!******************************************************
! Write DoS to file
! FIX THIS. MAKE IT POSSIBLE FOR THE USER TO CHANGE MAX CUTOFF FOR DoS OUTPUT
! Choose upper cutoff for DoS printing
  if( .true. )then
!   This corresponds to a 150 ps^-1 cutoff for output frequency
    imax = int(150.d0 * tau) - 1
!   imax can't be larger than the array size
    if( imax > (n+1)/2 )then
      imax = (n+1)/2
    end if
  else
    imax = (n+1)/2
  end if
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
end subroutine
!=================================================================================================
!=================================================================================================










!=================================================================================================
!=================================================================================================
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
!=================================================================================================
!=================================================================================================











!=================================================================================================
!=================================================================================================
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
!=================================================================================================
!=================================================================================================










!=================================================================================================
!=================================================================================================
subroutine exclude_volume_assign()

  integer :: group_start, group_end
  character*1 :: sg_char
  character*32 :: sg_string

  sg_string = ""
  group_start = 0
  do i = 1, len(exclude_volume_string)
    read(exclude_volume_string(i:i), '(A)', iostat=iostatus) sg_char
    if(sg_char == "-")then
      read(sg_string,*) group_start
      sg_string = ""
    else if( sg_char == " " .and. sg_string == "" )then
      continue
    else if(sg_char == "," .or. sg_char == " " .or. iostatus < 0)then
      read(sg_string,*) group_end
      if( group_start == 0 )then
        group_start = group_end
      end if
      do j=group_start,group_end
        exclude_volume(j) = .true.
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
end subroutine
!=================================================================================================
!=================================================================================================







end module
