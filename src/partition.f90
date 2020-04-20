module partition

  use constants
  use read_input
  use trajectory
! use topology
  use calc_dos
  use fluidicity

  implicit none


  real*8 :: omega
  real*8, allocatable :: delta_group(:,:), degf_group(:,:), degf_supergroup(:,:)
  real*8, allocatable :: f_group(:,:), f_supergroup(:,:), y_group(:,:), z_group(:,:)
  real*8, allocatable :: delta_supergroup(:,:), y_supergroup(:,:), z_supergroup(:,:)
  real*8, allocatable :: sigma_group(:,:), sigma_supergroup(:,:)
  real*8 :: res, D_rot_real, D_ideal, D_rot_ideal
  integer :: niter



  contains


subroutine dof_partition()
! Allocate arrays
  allocate( delta_group(1:ngroups,1:3) )
  allocate( f_group(1:ngroups,1:3) )
  allocate( y_group(1:ngroups,1:3) )
  allocate( z_group(1:ngroups,1:3) )
  allocate( sigma_group(1:ngroups, 1:3) )
  sigma_group = 0.d0

  if(calc_supergroups)then
    allocate( sigma_supergroup(1:nsupergroups, 1:3) )
    sigma_supergroup = 0.d0
  end if

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
        call get_omega(j, ngroups, dfloat(1+0*natoms_in_group), volume_group, sigma_group(1:ngroups,k), mass_group, "v", omega, &
                       (/ .false. /))
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
                                 mass_group, T, (/ V /), volume_group, niter, res, &
                                 sigma_group(1:ngroups,k), f_group(1:ngroups,k), (/ .false. /))
      do j = 1, ngroups
!       Get partial compressibility and partial packing fraction
        call get_compressibility(1, dfloat(1+0*natoms_in_group(j:j)), volume_group(j:J), sigma_group(j:j, k), &
                                 f_group(j:j, k), y_group(j,k), z_group(j,k), (/ .false. /))
      end do
    end do
  end if
! NOTE: I'M NOT OUTPUTTING FLUIDICITIES FOR SINGLE GROUPS ANYMORE, SINCE IT DOESN'T MAKE ANY SENSE
! Write to file
  open(unit=10, file="fluidicity_g", status="unknown")
  write(10,*) "# Group no.; effective deg. of freedom [CM; rot; vib]; Voronoi volume (nm^3); &
               &weight; supergroup"
  do j=1,ngroups
    write(10,'(I8,2X,F10.4,2X,F10.4,2X,F10.4,2X,F10.4,2X,F6.4,2X,I8)') &
          j, degf_group(j,1:3), volume_group(j), weight_group(j), group_belongs_to_supergroup(j)
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
        if( (k == 1 .and. .not. f_opt) .or. &
            (k == 1 .and. exclude_volume(j)) .or. &
            (k == 2 .and. .not. f_rot_opt) .or. &
            (k == 3))then
!         This menthod of calculating the fluidicity is essentially equivalent to Lin's method,
!         except for the inclusion of the "omega term" heuristically in the definition of the
!         normalized diffusivity "delta". This method is now deprecated and the more sophisticated
!         method based on numerical minimization of a penalty function should be employed.
!         For a monocomponent system, the results should be identical.
!         Available for testing and debugging purposes.
!         Vibrational fluidicity is not a well-defined concept and is calculated only for curious
!         purposes using this method. Vibrational fluidicity is not used to calculate vibrational entropy.
!         Calculate Delta:
          if(hs_formalism == "lin" .or. exclude_volume(j))then
            omega = 1.d0
            temp(1) = volume_supergroup(j)
          else
            call get_omega(j, nsupergroups, ngroups_in_supergroup_eff, volume_supergroup, &
                           sigma_supergroup(1:nsupergroups,k), mass_supergroup, "v", omega, exclude_volume)
!            temp(1) = V
            temp(1) = V_apparent(j)
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
                                 conv1 * twobykT * Ssupergroup(1:nsupergroups,1,k), mass_supergroup, T, V_apparent, &
                                 volume_supergroup, niter, res, sigma_supergroup(1:nsupergroups,k), &
                                 f_supergroup(1:nsupergroups,k), exclude_volume)
      do j = 1, nsupergroups
!       Get partial compressibility and partial packing fraction
!       This should not be used for excluded volumes -> use Lin's formalism instead
        if( .not. exclude_volume(j) )then
          call get_compressibility(1, ngroups_in_supergroup_eff(j:j), volume_supergroup(j:j), sigma_supergroup(j:j, k), &
                                   f_supergroup(j:j, k), y_supergroup(j,k), z_supergroup(j,k), exclude_volume(j:j))
        end if
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
                       sigma_supergroup(1:nsupergroups, k), mass_supergroup, "v", omega, exclude_volume)
!        temp(1) = V
        temp(1) = V_apparent(j)
        call get_zero_pressure_diffusivity(sigma_supergroup(j,k), omega, T, mass_supergroup(j), &
                                           ngroups_in_supergroup_eff(j), temp(1), D_ideal)
!       Ideal rotational diffusion coefficient -> D_rot_ideal
        call get_ideal_rotational_diffusivity(sigma_supergroup(j,k), T, mass_supergroup(j), &
                                              ngroups_in_supergroup_eff(j), volume_supergroup(j), D_ideal, D_rot_ideal)
!       The ratios between temperature and rotational temperatures are used to check whether the
!       molecules are linear or not
        temp(1:3) = conv4 * T / h**2.d0 * 8.d0 * pi**2.d0 * kB * eig_supergroup(j,1:3)
        if( temp(1) < 1.d0 )then
!         We have a monatomic particle
          f_supergroup(j, k) = 0.d0
          cycle
        else if( temp(3) < 1.d0 )then
!         We have a linear molecule 
          D_rot_ideal = D_rot_ideal / 4.d0 * ( eig_supergroup(j,1) + eig_supergroup(j,2) ) * &
                      ( 1.d0/eig_supergroup(j,1) + 1.d0/eig_supergroup(j,2) )
        else
          D_rot_ideal = D_rot_ideal / 9.d0 * ( eig_supergroup(j,1) + eig_supergroup(j,2) + eig_supergroup(j,3) ) * &
                      ( 1.d0/eig_supergroup(j,1) + 1.d0/eig_supergroup(j,2) + 1.d0/eig_supergroup(j,3) )
        end if
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
                 &fluidicity f [CM; rot; vib]; HS diameters (nm) [CM, rot, vib]; weight"
    do j=1,nsupergroups
      write(10,'(I8,2X,F10.4,2X,F10.4,2X,F10.4,2X,F10.4,2X,F7.4,2X,F7.4,2X,F7.4,2X,F10.4,2X,F10.4,2X,F10.4,2X,F12.5)') &
            j, degf_supergroup(j,1:3), volume_supergroup(j), f_supergroup(j,1:3), sigma_supergroup(j,1:3), weight_supergroup(j)
    end do
  end if
  close(10)
  write(*,*)'Done.                                  |'
  write(*,*)'                                       |'
  write(*,*)'.......................................|'
!******************************************************
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


end module
