module fluidicity

contains

subroutine fluidicity_calculator(nsupergroups, ngroups_in_supergroup, N_DoF, s0, M, T, V_apparent, &
                                 V_voro, niter, res, sigma, f, exclude_volume)

  implicit none

  real*8 :: T, V_apparent(:), pi, res_min, res, dsigma
  real*8 :: z, xi, ngroups_in_supergroup(:)
  integer :: nsupergroups, i, j, N_sigma, N_dim, k, niter, level, N_levels
  real*8 :: f_tol, step, dV, res_prev, start_scan, end_scan
  real*8 :: s0(:), N_DoF(:), M(:), f(:), V_voro(:), sigma(:), alpha
  real*8, allocatable :: omega(:), D(:), D0(:), sigma0(:), y(:), sigma_min(:), &
                         f_prev(:), grad(:), res_plus(:), res_minus(:), sigma_temp(:), sigma_prev(:)
  logical :: converged, step_too_large, exclude_volume(:)
!
! Variables for estimation of partial volumes using an approach based on vdW radii
! real*8, allocatable :: vdw_radii(:), pos(:,:)


! IMPORTANT: this subroutine should ignore excluded volumes. For those, sigma and f inputs should be the same as
! the outputs


  pi = dacos(-1.d0)

  alpha = 1.d0

!  allocate( s0(1:nsupergroups) )
!  allocate( M(1:nsupergroups) )
!  allocate( ngroups_in_supergroup(1:nsupergroups) )
!  allocate( N_DoF(1:nsupergroups) )
!  allocate( V_voro(1:nsupergroups) )
!  allocate( f(1:nsupergroups) )

! Allocate variables needed for the minimization
  allocate( omega(1:nsupergroups) )
!  allocate( sigma(1:nsupergroups) )
  allocate( sigma0(1:nsupergroups) )
  allocate( sigma_min(1:nsupergroups) )
  allocate( sigma_temp(1:nsupergroups) )
  allocate( sigma_prev(1:nsupergroups) )
  allocate( f_prev(1:nsupergroups) )
  allocate( grad(1:nsupergroups) )
  allocate( res_plus(1:nsupergroups) )
  allocate( res_minus(1:nsupergroups) )


! Initial estimate for sigma based on Voronoi volumes
  sigma0(1:nsupergroups) = (6.d0 * V_voro(1:nsupergroups) / ngroups_in_supergroup(1:nsupergroups) / pi)**(1./3.)

! Look in the interval 0.1*sigma0 to 2*sigma0 for a good staring guess
! for the minimization. We look at N_dim^N_sg points in this interval with exponential
! increments. In each (sigma) dimension, the number of points used is 16
!
! Points per dimension
  N_dim = 16
! Levels of refinement
  N_levels = 16
! Make sure they're not too many if a lot of supergroups are present in the system
  do while (N_dim**nsupergroups > 10000000)
    N_dim = N_dim / 2
  end do
  start_scan = 0.1d0
  end_scan = 2.d0
  do level = 1, N_levels
!   Set initial values
    do j = 1, nsupergroups
      if( .not. exclude_volume(j) )then
        sigma(j) = start_scan * sigma0(j)
      end if
    end do
    sigma0(1:nsupergroups) = sigma(1:nsupergroups)
!   Total number of points
    N_sigma = N_dim**nsupergroups
!   Increment of the exponent
    dsigma = (dlog10(end_scan) - dlog10(start_scan))/ log10(dfloat(N_dim)) / dfloat(N_dim - 1)
!   Initial residual values at the minimum
    sigma_min(1:nsupergroups) = sigma(1:nsupergroups)
    res_min = 1.d20
!   Loop through points
    do i = 1, N_sigma
!     Set up all the sigma values along each dimension
      do j = 1, nsupergroups
        if( .not. exclude_volume(j) )then
          k = mod( int((i-1)/N_dim**(j-1)), N_dim )
          sigma(j) = sigma0(j) * N_dim**(dfloat(k)*dsigma)
        end if
      end do
!     Calculate the residual for given sigma array
      call HS_res(nsupergroups, ngroups_in_supergroup, s0, T, M, N_DoF, V_apparent, V_voro, sigma, f, alpha, res, exclude_volume)
      if(res < res_min)then
        res_min = res
        sigma_min = sigma
      end if
    end do
    sigma0 = sigma_min
    start_scan = 0.5d0 + start_scan/2.d0
    end_scan = 0.5 + end_scan/2.d0
  end do

! Write out convergence towards minimum
  open(unit=10,file="log_opt",status="unknown")
  write(10,*) 0, sigma_min

! Approach minimum using a gradient descent. Written by Miguel Caro from scratch,
! so likely poorly optimized. Not a problem since this step will be very short compared
! with the overall execution time of the code.
! Convergence parameters
  step = 1.d0
  f_tol = 1.d-5
  do j = 1, nsupergroups
    if( .not. exclude_volume(j) )then
      f(j) = 1.d0
    end if
  end do
  f_prev = 0.d0
  i = 0
  res_prev = res_min
  sigma_prev = sigma_min
  converged = .false.
  do while (.not. converged)
    i = i + 1
    f_prev = f
    step_too_large = .true.
!   Calculate gradient of residual at sigma_min
    do while (step_too_large)
      do j = 1, nsupergroups
!       Excluded volumes should not contribute
        if( .not. exclude_volume(j) )then
          sigma_temp = sigma_min
!         Calculate partial derivatives dres/dsigma in the vicinity of sigma_min
          sigma_temp(j) = sigma_min(j)*(1.d0+step*0.01d0)
          call HS_res(nsupergroups, ngroups_in_supergroup, s0, T, M, N_DoF, V_apparent, V_voro, sigma_temp, f, alpha, &
                      res_plus(j), exclude_volume)
          sigma_temp(j) = sigma_min(j)*(1.d0-step*0.01d0)
          call HS_res(nsupergroups, ngroups_in_supergroup, s0, T, M, N_DoF, V_apparent, V_voro, sigma_temp, f, alpha, &
                      res_minus(j), exclude_volume)
        end if
      end do
!     Estimate gradient or res in the vicinity of sigma_min
      do j = 1, nsupergroups
!       Excluded volumes should not contribute
        if( .not. exclude_volume(j) )then
          grad(j) = (res_plus(j) - res_minus(j)) / sigma_min(j) / 0.02d0 / step
        else
          grad(j) = 0.d0
        end if
      end do
!     Get current residual
      call HS_res(nsupergroups, ngroups_in_supergroup, s0, T, M, N_DoF, V_apparent, V_voro, sigma_min, f, alpha, res_prev, &
                  exclude_volume)
!     Estimate new sigma from gradient descent and get new residual
      sigma_min(1:nsupergroups) = sigma_min(1:nsupergroups) - grad(1:nsupergroups) * step
      call HS_res(nsupergroups, ngroups_in_supergroup, s0, T, M, N_DoF, V_apparent, V_voro, sigma_min, f, alpha, res, &
                  exclude_volume)
!     Decrease step if divergence is observed and recalculate quantities
!     Also make sure there are no negative sigma values (this happens sometimes, I don't know why)
      if( res > res_prev .or. any(sigma_min < 0.d0) )then
        sigma_min = sigma_prev
        step = step/2.d0
      else
        sigma_prev = sigma_min
        res_prev = res
        step_too_large = .false.
      end if
    end do
    write(10,*) i, sigma_min, f
!   Check if we've achieved convergence
    converged = .true.
    do j = 1, nsupergroups
      if( dabs(f(j) - f_prev(j)) > f_tol )then
        converged = .false.
        exit
      end if
    end do
  end do
  niter = i
  do j = 1, nsupergroups
!   Excluded volumes should not contribute
    if( .not. exclude_volume(j) )then
      sigma(j) = sigma_min(j)
    end if
  end do
  close(10)

end subroutine











subroutine get_omega(this_group, ngroups, natoms_in_group, volume_group, sigma_group, mass_group, mode, omega, exclude_volume)

  implicit none

  real*8 :: omega
  real*8 :: volume_group(:), sigma_group(:), mass_group(:)
  integer :: ngroups, j, j2, this_group
  real*8 :: natoms_in_group(:)
  character*1 :: mode
  logical :: exclude_volume(:), exclude_default

  if( ngroups > size(exclude_volume) )then
    exclude_default = .true.
  else
    exclude_default = .false.
  end if

  omega = 0.d0

  if( exclude_default )then
    j2 = 1
  else
    j2 = this_group
  end if
  if( .not. exclude_volume(j2) )then
    if( mode == "v" )then
      do j=1, ngroups
        if( exclude_default )then
          j2 = 1
        else
          j2 = j
        end if
        if( .not. exclude_volume(j2) )then
          omega = omega + 1.d0/4.d0/dsqrt(2.d0) * natoms_in_group(j) / natoms_in_group(this_group) * &
                  (1.d0 + (volume_group(j)/volume_group(this_group))**(1.d0/3.d0))**2.d0 * &
                  dsqrt(1.d0 + mass_group(this_group)/mass_group(j))
        end if
      end do
    else if( mode == "s" )then
      do j=1, ngroups
        if( exclude_default )then
          j2 = 1
        else
          j2 = j
        end if
        if( .not. exclude_volume(j2) )then
          omega = omega + 1.d0/4.d0/dsqrt(2.d0) * natoms_in_group(j) / natoms_in_group(this_group) * &
                  (1.d0 + sigma_group(j)/sigma_group(this_group))**2.d0 * &
                  dsqrt(1.d0 + mass_group(this_group)/mass_group(j))
        end if
      end do
    end if
  else
    omega = 1.d0
  end if

end subroutine














subroutine get_compressibility(nsupergroups, ngroups_in_supergroup, V_apparent, sigma, f, xi, z, exclude_volume)

!***********************************************************************************
! This subroutine calculates the compressibility factor of a mixture of hard
! spheres according to the Mansoori-Carnahan-Starling-Leland expression [JCP 54,
! 1523 (1971)]. If you want to use the original expression simply set the fluidicity
! to one (f=1) when you call this subroutine.
!***********************************************************************************

  implicit none

  integer, intent(in) :: nsupergroups
  real*8, intent(in) :: V_apparent(:), sigma(:), f(:), ngroups_in_supergroup(:)
  real*8, intent(out) :: z, xi
  integer :: i, j, k
  real*8 :: pi, ngroups, sum, y1, y2, y3
  real*8, allocatable :: xi_i(:), delta(:,:)
  logical :: exclude_volume(:)

  pi = dacos(-1.d0)

  allocate( xi_i(1:nsupergroups) )
  allocate( delta(1:nsupergroups, 1:nsupergroups) )

  do i = 1, nsupergroups
!   Excluded volumes should not contribute
    if( .not. exclude_volume(i) )then
      xi_i(i) = pi / 6.d0 * f(i) * ngroups_in_supergroup(i) / V_apparent(i) * sigma(i)**3
    end if
  end do

  ngroups = 0.d0
  xi = 0.d0
  do i = 1, nsupergroups
!   Excluded volumes should not contribute
    if( .not. exclude_volume(i) )then
      xi = xi + xi_i(i)
      ngroups = ngroups + f(i) * ngroups_in_supergroup(i)
    end if
  end do

  delta = 0.d0
  do i = 1, nsupergroups
!   Excluded volumes should not contribute
    if( .not. exclude_volume(i) )then
      do j = i+1, nsupergroups
!       Excluded volumes should not contribute
        if( .not. exclude_volume(j) )then
          delta(i,j) = dsqrt(xi_i(i) * xi_i(j)) / xi * (sigma(i) - sigma(j))**2 / sigma(i) / sigma(j) * &
                       dsqrt(f(i) * ngroups_in_supergroup(i) * f(j) * ngroups_in_supergroup(j)) / &
                       ngroups
          delta(j,i) = delta(i,j)
        end if
      end do
    end if
  end do

  y1 = 0.d0
  do i = 1, nsupergroups
!   Excluded volumes should not contribute
    if( .not. exclude_volume(i) )then
      do j = i+1, nsupergroups
!       Excluded volumes should not contribute
        if( .not. exclude_volume(j) )then
          y1 = y1 + delta(i,j) * (sigma(i) + sigma(j)) / dsqrt(sigma(i) * sigma(j))
        end if
      end do
    end if
  end do

  y2 = 0.d0
  do i = 1, nsupergroups
!   Excluded volumes should not contribute
    if( .not. exclude_volume(i) )then
      do j = i+1, nsupergroups
!       Excluded volumes should not contribute
        if( .not. exclude_volume(j) )then
          sum = 0.d0
          do k = 1, nsupergroups
!           Excluded volumes should not contribute
            if( .not. exclude_volume(k) )then
              sum = sum + xi_i(k) / xi * dsqrt(sigma(i) * sigma(j)) / sigma(k)
            end if
          end do
          y2 = y2 + delta(i,j) * sum
        end if
      end do
    end if
  end do

  y3 = 0.d0
  do i = 1, nsupergroups
!   Excluded volumes should not contribute
    if( .not. exclude_volume(i) )then
      y3 = y3 + (xi_i(i) / xi)**(2.d0/3.d0) * (f(i) * ngroups_in_supergroup(i) / ngroups)**(1.d0/3.d0)
    end if
  end do
  y3 = y3**3

  z = ( (1.d0 + xi + xi**2) - 3.d0*xi * (y1 + y2*xi) - y3*xi**3) / (1.d0 - xi)**3
end subroutine














subroutine HS_res(nsupergroups, ngroups_in_supergroup, s0, T, M, N_DoF, V_apparent, V_partial, sigma, f, alpha, res, &
                  exclude_volume)

!***********************************************************************************
! This subroutine calculates the residual difference between a) the ratio of real
! diffusivity to zero-pressure diffusivity and b) the same ratio predicted by
! Enskog theory from the compressibity of the hard-sphere fluid.
!***********************************************************************************

  implicit none

  real*8 ::  s0(:), N_DoF(:), M(:), f(:), sigma(:), V_partial(:)
  real*8 :: xi, res, pi, V_apparent(:), kB, T, z, temp, alpha, ngroups_in_supergroup(:)
  integer :: nsupergroups
  integer :: i, j
  real*8, allocatable :: D(:), D0(:), z_i(:), xi_i(:), omega(:)
  logical :: exclude_volume(:)

  pi = dacos(-1.d0)

  allocate( D(1:nsupergroups) )
  allocate( D0(1:nsupergroups) )
  allocate( z_i(1:nsupergroups) )
  allocate( xi_i(1:nsupergroups) )
  allocate( omega(1:nsupergroups) )

  do i = 1, nsupergroups
    if( .not. exclude_volume(i) )then
      call get_omega(i, nsupergroups, ngroups_in_supergroup, V_partial, sigma, M, "sigma", omega(i), exclude_volume)
    end if
  end do

  res = 0.d0
  do i = 1, nsupergroups
!   Excluded volumes should not contribute to the penalty function
    if( .not. exclude_volume(i) )then
!     Get fluidicity from D(N) and D0(N)
      call get_diffusivity(s0(i), T, M(i), N_DoF(i), D(i))
      call get_zero_pressure_diffusivity(sigma(i), omega(i), T, M(i), ngroups_in_supergroup(i), V_apparent(i), D0(i))
      f(i) = D(i)/D0(i)
!     Get D(f*N) and D0(f*N)
      call get_diffusivity(s0(i), T, M(i), f(i)*N_DoF(i), D(i))
      call get_zero_pressure_diffusivity(sigma(i), omega(i), T, M(i), f(i)*ngroups_in_supergroup(i), V_apparent(i), D0(i))
!     Get partial compressibility
      call get_compressibility(1, ngroups_in_supergroup(i:i), V_partial(i:i), sigma(i:i), f(i:i), xi_i(i), z_i(i), &
                               exclude_volume(i:i))
!     Update residual with the difference between the actual D(f*N)/D0(f*N) ratio and that
!     predicted by Enskog's theory for hard spheres
      res = res + V_partial(i)/V_apparent(i) * (D(i)/D0(i) - 4.d0*xi_i(i)/(z_i(i)-1.d0))**2
!      res = res + 1.d0/dfloat(nsupergroups) * (D(i)/D0(i) - 4.d0*xi_i(i)/(z_i(i)-1.d0))**2
    end if
  end do

! Obtain compressibility of the hard-sphere mixture with the obtained fluidicities
  call get_compressibility(nsupergroups, ngroups_in_supergroup, V_apparent, sigma, f, xi, z, exclude_volume)

! Update the residual with N_sg times the difference between the MCSL HS compressibility
! and that obtained from summing the partial compressibilities
  temp = 0.d0
  do i = 1, nsupergroups
!   Excluded volumes should not contribute to the penalty function
    if( .not. exclude_volume(i) )then
      temp = temp + V_partial(i)/V_apparent(i) * z_i(i)
    end if
  end do

  res = res + alpha * (z - temp)**2

end subroutine






subroutine get_diffusivity(s0, T, M, N_DoF, D)

!***********************************************************************************
! This subroutine calculates the real diffusivity of the
! system based on the MD information. Note that you have to pass the number
! of degrees of freedom as argument (N_DoF), rather than the number of particles (N)
! If you want to pass 3*N as argument, make sure to pass 3.d0*dfloat(N) instead!
! Units information. The subroutine expects the following units:
!
! [s0] = ps
! [T] = K
! [M] = amu
!
! D is returned in:
!
! [D] = nm^2 / ps
!***********************************************************************************

  implicit none

  real*8 :: D, s0, N_DoF, kB, T, M
  real*8 :: ps, eV, amu

  ps = 1.d-12
  eV = 1.602176565d-19
  amu = 1.660538921d-27
! Boltzmann's constant in eV/K
  kB = 8.6173324d-5

  D = s0 * kB * T / 4.d0 / M / N_DoF
! Get right units. D above is in ps * eV / amu, and we want it back in nm^2 / ps.
! We first transform to SI (m^2 / s) and then to nm^2 / ps
  D = D * ps * eV / amu * 1.d6

end subroutine




subroutine get_zero_pressure_diffusivity(sigma, omega, T, M, N, V, D0)

!***********************************************************************************
! This subroutine calculates the hard-sphere diffusivity of the
! system in the zero pressure limit.
!
! N is the particle number in this case and V is the total volume. Even though
! N is originally an integer, it should be passed as a real*8
!
! Units information. The subroutine expects the following units:
!
! [sigma] = nm
! [V] = nm^3
! [T] = K
! [M] = amu
!
! D0 is returned in:
!
! [D0] = nm^2 / ps
!***********************************************************************************

  implicit none

  real*8 :: D0, V, kB, T, sigma, omega, M, N
  real*8 :: nm, amu, eV, pi

  nm = 1.d-9
  eV = 1.602176565d-19
  amu = 1.660538921d-27
  pi = dacos(-1.d0)
! Boltzmann's constant in eV/K
  kB = 8.6173324d-5

  D0 = 3.d0/8.d0 * V / N / sigma**2 / omega * dsqrt(kB * T / pi / M)
! Get right units. D0 above is in nm * sqrt(eV / amu), and we want it back in nm^2 / ps.
! We first transform to SI (m^2 / s) and then to nm^2 / ps
  D0 = D0 * nm * dsqrt(eV / amu) * 1.d6

end subroutine







subroutine get_ideal_rotational_diffusivity(sigma, T, M, N, V, D0_trn, D0)

!***********************************************************************************
! This subroutine calculates the ideal rotational diffusivity of the
! system for spheres.
!
! N is the particle number in this case and V is the partial volume. Even though
! N is originally an integer, it should be passed as a real*8
!
! Units information. The subroutine expects the following units:
!
! [sigma] = nm
! [V] = nm^3
! [T] = K
! [M] = amu
! [D0_trn] = nm^2 / ps
!
! D0 is returned in:
!
! [D0] = 1 / ps
!***********************************************************************************

  implicit none

  real*8 :: D0, V, kB, T, sigma, M, N, D0_trn
  real*8 :: nm, amu, eV, ps

  nm = 1.d-9
  eV = 1.602176565d-19
  amu = 1.660538921d-27
  ps = 1.d-12
! Boltzmann's constant in eV/K
  kB = 8.6173324d-5

  D0 = 6.d0 / 5.d0 * kB * T / (N / V * M * sigma**3) / D0_trn
! Get right units. D0 above is in eV / amu * ps / nm^2, and we want it back in 1 / ps.
! We first transform to SI (1 / s) and then to 1 / ps
  D0 = D0 * eV / amu * ps / nm**2 * 1.d-12

end subroutine






subroutine get_real_rotational_diffusivity(tau, n, N_DoF, ngroups, group_in_supergroup, birth_time, death_time, j2, w, D)

!***********************************************************************************
! This subroutine calculates the real rotational diffusivity of the
! system from MD information.
!
! Units information. The subroutine expects the following units:
!
! [tau] = ps
! [n] = number
! [N_DoF] = number
! [w] = 1 / ps
! [ngroups] = number
! [group_in_supergroup] = index
! [j2] = index
!
! D is returned in:
!
! [D] = 1 / ps
!***********************************************************************************

  implicit none

  real*8 :: D, N_DoF, tau, S
  real*4 :: w(:,:,:)
  integer :: n, ngroups, j, k, i, j2, group_in_supergroup(:,:), birth_time(:), death_time(:)

  D = 0.d0
  do j = 1, ngroups
    do k = 1, 3
       S = 0.d0
       do i = 1, n
         if( i >= birth_time(group_in_supergroup(j2,j)) .and. i < death_time(group_in_supergroup(j2,j)) )then
           S = S + w(group_in_supergroup(j2,j), k, i)
         end if
       end do
       D = D + S**2
     end do
   end do
 D = tau / (2.d0 * dfloat(n)**2 * N_DoF) * D

end subroutine

end module
