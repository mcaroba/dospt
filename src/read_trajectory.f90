!=================================================================================================
!=================================================================================================
subroutine read_trajectory(n, natoms, tau, L, mode, positions, velocities, estimate_vel, error, &
                           m, nmasses, mass_types, mass_types_values, volumes, volumes_temp)

  implicit none

  real*8 :: L(1:3)
  real*4 :: positions(:,:,:), velocities(:,:,:), volumes(:,:)
  real*8 :: r1(1:3), r2(1:3), r3(1:3), r4(1:3), r5(1:3)
  real*8 :: r0(1:3), v0(1:3), dt, tau
  integer :: n, natoms, i, j, t1, t2, t3, t4, t5, t0, iostatus
  character*64 :: crap, format_gromacs1, format_gromacs2
  character*8 :: filename
  character*3 :: mode
  logical :: estimate_vel, error
  real*8 :: update_bar, volumes_temp(:)
  integer :: k, j2, ndec, k2, j3
  character*1 :: cha

  integer :: nmasses
  character*16 :: mass_label
  real*8 :: m(:)
  real*8 :: mass_types_values(:)
  character*16 :: mass_types(:)

  real*8, allocatable :: x(:), y(:), z(:)

! Quartic approximation can be hard-switched on here - only for testing and debugging
  logical :: quartic = .false.


!******************************************************
! The get_mass subroutine requires an interface
interface
subroutine get_mass(mass_label, mass, nmasses, mass_types, mass_types_values)
  integer, intent(in) :: nmasses
  character*16, intent(in) :: mass_label
  real*8, intent(out) :: mass
  real*8, intent(in)  :: mass_types_values(:)
  character*16, intent(in)  :: mass_types(:)
end subroutine
end interface
!******************************************************


  allocate( x(1:natoms) )
  allocate( y(1:natoms) )
  allocate( z(1:natoms) )

  dt = tau / dfloat(n-1)

! Filename and other options
  if(mode == "gro")then
    filename = "traj.gro"
  else if(mode == "xyz")then
    filename = "traj.xyz"
    estimate_vel = .true.
  end if

! Reads in trajectory file and prints messages
  open(unit=10, file=filename, status="old", iostat=iostatus)
  write(*,*)'                                       |'
  write(*,'(A,1X,A,1X,A)')' Reading', filename, 'file...               |'
  if(iostatus/=0)then
    close(10)
    write(*,'(A,1X,A,1X,A)')' ERROR:', filename, 'file could not be found|'
    write(*,*)'                                       |'
    write(*,*)'.......................................|'
    error = .true.
    return
  end if
  write(*,*)'                                       |'
  write(*,*)'... and calculating Voronoi volumes... |'
  write(*,*)'                                       |'
  write(*,*)'Progress (I/O usually takes longest):  |'
  write(*,*)'                                       |'
  write(*,'(1X,A)',advance='no')'['
  update_bar = dfloat(n)/36.d0
  k = 1
  j2 = 1
!
  do i=1,n
!   Update progress bar every n/36 iterations
    if(dfloat(i) > dfloat(k)*update_bar)then
      write(*,'(A)', advance='no')'='
      k = k + 1
    end if
!
!   xyz format
!   No input velocities, positions given in A, then converted to nm
    if(mode == "xyz")then
      read(10,*)
      read(10,*)
      do j=1,natoms
        read(10,*) mass_label, positions(j,1,i), positions(j,2,i), positions(j,3,i)
        if(i == 1)then
          call get_mass(mass_label, m(j), nmasses, mass_types, mass_types_values)
        end if
        positions(j,1:3,i) = 0.1 * positions(j,1:3,i)
        x(j) = positions(j,1,i)
        y(j) = positions(j,2,i)
        z(j) = positions(j,3,i)
      end do
!
!   gromacs format
!   Input velocities posssible. Positions in nm and velocities in nm/ps
    else if(mode == "gro")then
      read(10,*)
      read(10,*)
!     Prevent missing spaces between numbers from raising a runtime reading error by
!     learning the precision of Gromacs output. The usual Gromacs format is (i5,2a5,i5,3f8.3,3f8.4)
!     However, the user can define more decimal places for positions and velocities, e.g.
!     for DoSPT it should be at least 5 decimals for positions and 5+1 for velocities (velocities
!     are written by Gromacs always with an extra decimal). Here we figure out how many decimal places
!     are being used so that we can generate the correct format for read(,). All of this is necessary
!     bacause sometimes Gromacs' format does not leave a space between e.g. two velocities.
      if(i == 1)then
        read(10,*) crap, crap, crap, crap
        k2 = len_trim(crap)
        crap = trim(crap)
        j3 = 0
        do j = 1, k2
          j3 = j3 + 1
          if(crap(j:j) == ".")then
            ndec = k2 - j3
            exit
          end if
        end do
        backspace(10)
        write(format_gromacs1, '(A,I2,A,I2,A,I2,A,I2,A)') '(I5,2A5,I5,3F', 5+ndec, '.', ndec, ',3F', 5+ndec, '.', 1+ndec, ')'
        write(format_gromacs2, '(A,I2,A,I2,A)') '(I5,2A5,I5,3F', 5+ndec, '.', ndec, ')'
        format_gromacs1 = trim(format_gromacs1)
      end if
      do j=1,natoms
        if(estimate_vel)then
          read(10, format_gromacs2) k2, crap, crap, k2, positions(j,1:3,i)
        else
          read(10, format_gromacs1) k2, crap, crap, k2, positions(j,1:3,i), velocities(j,1:3,i)
        end if
        if(i == 1)then
          read(crap,*) mass_label
          call get_mass(mass_label, m(j), nmasses, mass_types, mass_types_values)
        end if
        x(j) = positions(j,1,i)
        y(j) = positions(j,2,i)
        z(j) = positions(j,3,i)
      end do
      read(10,*)
    end if
!
!   Get Voronoi cell volumes every 100 time steps
    if( i == 1 )then
      call voronoi_volumes(natoms, L, x, y, z, volumes_temp)
      volumes(1:natoms,j2) = volumes_temp(1:natoms)
      j2 = j2 + 1
    else if( mod(i,100) == 0 )then
      call voronoi_volumes(natoms, L, x, y, z, volumes_temp)
      volumes(1:natoms,j2) = volumes_temp(1:natoms)
      j2 = j2 + 1
    end if
  end do
  close(10)
  write(*,'(A4)')'=] |'
  write(*,*)'                                       |'
  write(*,*)'.......................................|'


! Estimate velocities from positions
  if(estimate_vel)then
    write(*,*)'                                       |'
    write(*,*)'Estimating velocities from             |'
    write(*,*)'trajectory...                          |'
    do i=1,n
      if(i == 1)then
        t1 = i
        t2 = i+1
        t3 = i+2
        t0 = 1
      else if(i == n)then
        t1 = i-2
        t2 = i-1
        t3 = i
        t0 = n
      else
        t1 = i-1
        t2 = i
        t3 = i+1
        t0 = i
      end if
!     Quartic approximation - only for testing and debugging!
      if(quartic)then
        if(i == 1)then
          t1 = i
          t2 = i+1
          t3 = i+2
          t4 = i+3
          t5 = i+4
          k = 1
        else if(i == 2)then
          t1 = i-1
          t2 = i
          t3 = i+1
          t4 = i+2
          t5 = i+3
          k = 2
        else if(i == n-1)then
          t1 = i-3
          t2 = i-2
          t3 = i-1
          t4 = i
          t5 = i+1
          k = 4
        else if(i == n)then
          t1 = i-4
          t2 = i-3
          t3 = i-2
          t4 = i-1
          t5 = i
          k = 5
        else
          t1 = i-2
          t2 = i-1
          t3 = i
          t4 = i+1
          t5 = i+2
          k = 3
        end if
      end if
      do j=1,natoms
        r1(1:3) = positions(j,1:3,t1)
        r2(1:3) = positions(j,1:3,t2)
        r3(1:3) = positions(j,1:3,t3)
!       Quartic approximation - only for testing and debugging!
        if(quartic)then
          r4(1:3) = positions(j,1:3,t4)
          r5(1:3) = positions(j,1:3,t5)
        end if
        call restore_replica_traj(L, r1, r2, r3)
!       Quartic approximation - only for testing and debugging!
        if(quartic)then
          call restore_replica_traj(L, r1, r4, r5)
        end if
        call estimate_r_v(dfloat(t1)*dt, dfloat(t2)*dt, dfloat(t3)*dt, r1, r2, r3, dfloat(t0)*dt, r0, v0)
!       Quartic approximation - only for testing and debugging!
        if(quartic)then
          call estimate_v_quartic(r1, r2, r3, r4, r5, dt, k, v0)
        end if
        velocities(j,1:3,i) = v0(1:3)
      end do
    end do
    write(*,*)'                                       |'
    write(*,*)'Done.                                  |'
    write(*,*)'                                       |'
    write(*,*)'.......................................|'
  end if
end subroutine
!=================================================================================================
!=================================================================================================





!=================================================================================================
!=================================================================================================
subroutine restore_replica_traj(L, r1, r2, r3)

  implicit none

  real*8 :: r1(1:3), r2(1:3), r3(1:3)
  real*8 :: L(1:3), d21, d31
  integer :: i, j
  real*8 :: new_r2(1:3), new_r3(1:3)

! Brings r2 and r3 to the same set of periodic replicas
! that r1 belongs to, according to a minimum image convention 

  do i=1,3
    d21=1.d10
    d31=1.d10
    do j=-1,1
      if(dabs(r2(i) + L(i)*dfloat(j) - r1(i)) < d21)then
        d21 = dabs(r2(i) + L(i)*dfloat(j) - r1(i))
        new_r2(i) = r2(i) + L(i)*dfloat(j)
      end if
      if(dabs(r3(i) + L(i)*dfloat(j) - r1(i)) < d31)then
        d31 = dabs(r3(i) + L(i)*dfloat(j) - r1(i))
        new_r3(i) = r3(i) + L(i)*dfloat(j)
      end if
    end do
  end do

  do i=1,3
    r2(i) = new_r2(i)
    r3(i) = new_r3(i)
  end do

end subroutine
!=================================================================================================
!=================================================================================================





!=================================================================================================
!=================================================================================================
subroutine estimate_r_v(t1, t2, t3, r1, r2, r3, t0, r0, v0)

  implicit none

  real*8 :: r1(1:3), r2(1:3), r3(1:3)
  real*8 :: t1, t2, t3, t0, r0(1:3), v0(1:3)
  real*8 :: a(1:3), b(1:3), c(1:3), temp
  integer :: i

! Performs a quadratic approximation given points r1, r2 and r3
! evaluated at t1, t2, t3. It then returns r0 and v0 evaluated
! at t0 using the quadratic expression. Obviously, t0 should
! be in the interval [min(t1,t2,t3):max(t1,t2,t3)] for this
! to be a good approximation

! Calculate the interpolation coefficients for the expression
! r(i) = a(i)*t^2 + b(i)*t + c(i) and evaluate it at t0

  temp = (t1 - t2)*(t1 - t3)*(t2 - t3)
  do i=1,3
    a(i) = t3*(-r1(i) + r2(i)) + t2*(r1(i) - r3(i)) + t1*(-r2(i) + r3(i))
    a(i) = a(i) / temp
    b(i) = t3**2*(r1(i) - r2(i)) + t1**2*(r2(i) - r3(i)) + t2**2*(-r1(i) + r3(i))
    b(i) = b(i) / temp
    c(i) = t1*t3*(-t1 + t3)*r2(i) + t2**2*(t3*r1(i) - t1*r3(i)) + t2*(-t3**2*r1(i) + t1**2*r3(i))
    c(i) = c(i) / temp

    r0(i) = a(i)*t0**2 + b(i)*t0 + c(i)
    v0(i) = 2.d0*a(i)*t0 + b(i)
  end do

end subroutine
!=================================================================================================
!=================================================================================================





!=================================================================================================
!=================================================================================================
subroutine estimate_v_quartic(r1, r2, r3, r4, r5, dt, k, v0)

! This subroutine is provided only for testing and debugging purposes
! The quartic approximation does not perform well if the integration
! time step is too long

  implicit none

  real*8 :: r1(1:3), r2(1:3), r3(1:3), r4(1:3), r5(1:3)
  real*8 :: dt, v0(1:3)
  integer :: k

! Performs a quartic approximation given points r1, r2, r3
! r4 and r5, provided the time step dt is constant.
! It then returns v0 evaluated at the same time as rj
! It is assumed that r1, r2, r3, r4 and r5 are given sequentially
! for increasing time, i.e. if r3 is given at t then r2 is given
! at t-dt etc.

  if(k == 1)then
    v0 = (- 25.d0*r1 + 48.d0*r2 - 36.d0*r3 + 16.d0*r4 - 3.d0*r5) / 12.d0 / dt
  else if(k == 2)then
    v0 = (- 3.d0*r1 - 10.d0*r2 + 18.d0*r3 - 6.d0*r4 + r5) / 12.d0 / dt
  else if(k == 3)then
    v0 = (r1 - 8.d0*r2 + 8.d0*r4 - r5) / 12.d0 / dt
  else if(k == 4)then
    v0 = (- r1 + 6.d0*r2 - 18.d0*r3 + 10.d0*r4 + 3.d0*r5) / 12.d0 / dt
  else if(k == 5)then
    v0 = (+ 3.d0*r1 - 16.d0*r2 + 36.d0*r3 - 48.d0*r4 + 25.d0*r5) / 12.d0 / dt
  end if

end subroutine
!=================================================================================================
!=================================================================================================
