!=================================================================================================
!=================================================================================================
subroutine spheres_volume(pos, r, N, N_k, V)

!***********************************************************************************
! This subroutine calculates the total volume of a set of N spheres that might be
! overlapping. The ovelapping volume is counted only once. This is based on a
! Monte-Carlo-type method. Written by Miguel Caro from scratch, so most likely not
! properly optimized.
!***********************************************************************************

  implicit none

  real*8, intent(in) :: pos(:,:), r(:)
  real*8, intent(out) :: V
  integer, intent(in) :: N, N_k
  integer :: i, j, k, N_in, N_total, seed, N_spheres, N_even, N_odd, N_out
  real*8 :: pi, V_box, rmin(1:3), rmax(1:3), rp(1:3), dV, V12, V3, temp, min
  real*8 :: x1, x2, y1, y2, z1, z2, V_in, V_out
  real*8, allocatable :: faces(:, :)
  real*8, allocatable :: random(:, :), volume_of_box(:), rbox(:,:)
  integer, allocatable :: counter(:), counter_in(:), counter_out(:)
  integer :: loc, N_boxes, i2, j2, k2
  logical :: already_checked, point_in
  logical, allocatable :: overlap(:,:)

  allocate( random(1:N_k, 1:3) )

  pi = dacos(-1.d0)

! Trivial solution for one sphere
  if(N == 1)then
    V = 4.d0 / 3.d0 * pi * r(1)**3
    return
  end if

! If there are only two spheres, or if the spheres overlap only two at a time,
! one can calculate the overlapping volume analytically. This means the value
! is exact and obtained much more cheaply than with the Monte-Carlo procedure
!
! Solution for only two spheres:
  if(N == 2)then
    call spheres_intersection(pos(1,1:3), pos(2,1:3), r(1), r(2), dV)
    V = 4.d0 / 3.d0 * pi * (r(1)**3 + r(2)**3) - dV
    return
  end if

!**********************************************************
! WARNING:
! This method only performs well if the number of points is large compared to the number of
! boxes. To improve the efficiency, I would need to code the box generation in such a way
! that some of the very small boxes are combined into larger ones. One can test for example
! that for non-overlapping spheres the calculated number of boxes can be larger than the
! number of spheres, which is suboptimal. For real systems some boxes might be very small.
! The other method below based on one single box around each sphere is much better than
! this for most purposes, and converges quicker than the brute-force Monte Carlo.
! Use at your own risk.
!
if(1 == 0)then ! change condition to run
! Create system of boxes
! Sort the faces
  allocate( faces(1:2*N, 1:3) )
  do k = 1, 3
  do i = 1, N
    faces(i, k) = pos(i, k) - r(i)
    faces(N + i, k) = pos(i, k) + r(i)
  end do
  end do
  do k = 1, 3
  do i = 1, 2*N-1
    min = faces(i, k)
    loc = i
    do j = i+1, 2*N
      if(faces(j, k) < min)then
        min = faces(j, k)
        loc = j
      end if
    end do
    temp = faces(i, k)
    faces(i, k) = faces(loc, k)
    faces(loc, k) = temp
  end do
  end do
! Store box info for those boxes that are not empty
  allocate( overlap(1:8*(N-1)**3, 1:N) )
  allocate( volume_of_box(1:8*(N-1)**3) )
  allocate( rbox(1:8*(N-1)**3, 1:6) )
  N_boxes = 0
  V_box = 0.d0
  do i2 = 1, 2*N-1
    x1 = faces(i2, 1)
    x2 = faces(i2+1, 1)
  do j2 = 1, 2*N-1
    y1 = faces(j2, 2)
    y2 = faces(j2+1, 2)
  do k2 = 1, 2*N-1
    z1 = faces(k2, 3)
    z2 = faces(k2+1, 3)
    do i=1,N
      rmax(1:3) = pos(i,1:3) + r(i)
      rmin(1:3) = pos(i,1:3) - r(i)
!     Check that there is overlap between this box and this sphere
      if( rmax(1) > x1 .and. rmin(1) < x2 .and. rmax(2) > y1 .and. rmin(2) < y2 .and. rmax(3) > z1 .and. rmin(3) < z2)then
        if( .not. already_checked )then
          N_boxes = N_boxes + 1
          rbox(N_boxes, 1:6) = (/ x1, x2, y1, y2, z1, z2 /)
          volume_of_box(N_boxes) = (x2-x1)*(y2-y1)*(z2-z1)
          V_box = V_box + (x2-x1)*(y2-y1)*(z2-z1)
          already_checked = .true.
        end if
        overlap(N_boxes, i) = .true.
      end if
    end do
    already_checked = .false.
  end do
  end do
  end do

! Start Monte Carlo sampling. We choose the number of points placed in each box from the
! relative volume of that box compared to the volume of all the boxes combined (approx.)
  call init_random_seed()
  call random_number(random)
  N_total = 0
  N_in = 0
  do i = 1, N_boxes
    do k = 1, int(dfloat(N_k)*volume_of_box(i)/V_box)
!   Generate random point
      rp(1) = rbox(i, 2) + random(N_total + k,1)*(rbox(i, 1)-rbox(i, 2))
      rp(2) = rbox(i, 4) + random(N_total + k,2)*(rbox(i, 3)-rbox(i, 4))
      rp(3) = rbox(i, 6) + random(N_total + k,3)*(rbox(i, 5)-rbox(i, 6))
      N_total = N_total + 1
!     WARNING: this loop should be optimized by building up a list of spheres which overlap with
!     this box and only looping through those
      do j=1,N
!       If point is inside increment N_in and go for next point
        if( overlap(i, j) )then
          if( dot_product(rp(1:3)-pos(j,1:3),rp(1:3)-pos(j,1:3)) < r(j)**2 )then
            N_in = N_in + 1
            exit
          end if
        end if
      end do
    end do
  end do
! Get volume
  V = V_box / dfloat(N_total) * dfloat(N_in)
end if ! change condition to run
!**********************************************************


!**********************************************************
! A box is built around each sphere in this method. Then random points are
! placed in each box and the overlapping regions of the boxes is taken into
! account. This method converges faster than the one above and the one
! below. USE THIS ONE.
!
! With typical atomic (vd Waals) volumes, typical overlaps, etc. about 100k sampling
! points per atom will lead to errors below 0.1% in the calculated molecular
! volumes.
!
if(1 == 1)then ! change condition to (not) run
! Generate collection of random numbers
  call init_random_seed()
  call random_number(random)

! In itialize the counters
  allocate( counter_in(1:N) )
  allocate( counter_out(1:N) )
  counter_in = 0
  counter_out = 0

! Calculate total volume of all boxes put together
  allocate( volume_of_box(1:N) )
  V_box = 0.d0
  do i = 1, N
    volume_of_box(i) = 8.d0*r(i)**3
    V_box = V_box + volume_of_box(i)
  end do

  N_total = 0
! Loop through all boxes
  do i = 1, N
!   Loop through a number of random points proportional to the volume of the box
    do k = 1, int(dfloat(N_k)*volume_of_box(i)/V_box)
      N_total = N_total + 1
!     Generate random point within this particular box
      rp(1:3) = pos(i, 1:3)-r(i) + random(N_total,1:3)*2.d0*r(i)

!     Check two things: whether the point is inside/outside the overlapping spheres and
!     how many boxes the point is contained in
      point_in = .false.
      N_boxes = 0
      do j = 1, N
!       Check whether point is inside the spheres
        if( dot_product(rp(1:3)-pos(j,1:3),rp(1:3)-pos(j,1:3)) < r(j)**2 .and. (.not. point_in) )then
          point_in = .true.
        end if
!       If point is inside the boxes increment N_boxes
        if( pos(j, 1) - r(j) < rp(1) .and. rp(1) < pos(j, 1) + r(j) .and. &
            pos(j, 2) - r(j) < rp(2) .and. rp(2) < pos(j, 2) + r(j) .and. &
            pos(j, 3) - r(j) < rp(3) .and. rp(3) < pos(j, 3) + r(j) )then
          N_boxes = N_boxes + 1
        end if
      end do
!     Update the counters
      if( point_in )then
        counter_in(N_boxes) = counter_in(N_boxes) + 1
      else
        counter_out(N_boxes) = counter_out(N_boxes) + 1
      end if
    end do
  end do
! Get volume with proper weighting of the points that fell in several boxes a the same time
  V_in = 0.d0
  V_out = 0.d0
  do i = 1, N
    V_in = V_in + V_box/dfloat(N_total) * dfloat(counter_in(i)) / dfloat(i)
    V_out = V_out + V_box/dfloat(N_total) * dfloat(counter_out(i)) / dfloat(i)
  end do
  V = V_in
end if ! change condition to run
!**********************************************************


!**********************************************************
! WARNING: this code below is a brute-force Monte Carlo which
! does not perform the smart partitioning of space into nice boxes
! Unfortunately more work needs to go into the boxes code before it's safe
! to use for all systems.
!
if(1 == 0)then ! change condition to run
! For more than 3 spheres
! Get V using pure Monte Carlo
! Create a box around the overlapping spheres
  rmin(1:3) = pos(1,1:3) - r(1)
  rmax(1:3) = pos(1,1:3) + r(1)
  do i=2,N
    do j=1,3
      if(pos(i,j) - r(i) < rmin(j))then
        rmin(j) = pos(i,j) - r(i)
      end if
      if(pos(i,j) + r(i) > rmax(j))then
        rmax(j) = pos(i,j) + r(i)
      end if
    end do
  end do

! Randomly place points in the box, then check whether they are inside or outside
! The approximated volume is V_box / N_total * N_in
  V_box = (rmax(1)-rmin(1))*(rmax(2)-rmin(2))*(rmax(3)-rmin(3))

  call init_random_seed()
  call random_number(random)

  N_total = 0
  N_in = 0
  do k = 1, N_k
!   Generate random point
    rp(1:3) = rmin(1:3) + random(k,1:3)*(rmax(1:3)-rmin(1:3))

    N_total = N_total + 1
    do j=1,N
!     If point is inside increment N_in and go for next point
      if( dot_product(rp(1:3)-pos(j,1:3),rp(1:3)-pos(j,1:3)) < r(j)**2 )then
        N_in = N_in + 1
        exit
      end if
    end do
  end do
! Get volume
  V = V_box / dfloat(N_total) * dfloat(N_in)
end if ! change condition to run
!**********************************************************


!**********************************************************
! WARNING!!!
!
! This stuff below is an unsuccessful attempt at solving the problem using
! a hybrid analytical/numerical routine. I didn't manage to get the coefficients
! for the numerical corrections right. It would be interesting to finish the derivations
! since it would improve the accuracy of the results. The code is not useable as is.
! You have been warned!
!
! The region in space that belongs to any two spheres at most (that is,
! it does not belong at the same time to any three or more spheres) can be
! calculated analytically as V_1,2 = S_i (V_i) - S_i,j>i (V_ij).
! The remaining volume, that is V_3 = S_i,j>i,k>i (V_ijk), can be calculated numerically
! and added as a correction: V = V_1,2 + V_3. Note that V_1,2 can be negative
!
! The strategy is create a box around the set of spheres, and use a Monte Carlo
! sampling method to estimate the volume V_3 from the number of random points in space
! that fall inside at least 3 different spheres.

if(1 == 0)then ! change condition to run
!  Get V12
   V12 = 0.d0
   do i = 1, N
     V12 = V12 + 4.d0 / 3.d0 * pi * r(i)**3
     do j = i+1, N
       call spheres_intersection(pos(i,1:3), pos(j,1:3), r(i), r(j), dV)
       V12 = V12 - dV
     end do
   end do

! Get V3 using the Monte Carlo procedure
   call init_random_seed()
   call random_number(random)

  allocate( counter(0:N) )
  counter = 0

  N_in = 0
  do k = 1, N_k
!   Generate random point
    rp(1:3) = rmin(1:3) + random(k,1:3)*(rmax(1:3)-rmin(1:3))

!   Check the number of spheres inside of which this point resides
    N_spheres = 0
    do j=1,N
!     If point is inside increment N_spheres
      if( dot_product(rp(1:3)-pos(j,1:3),rp(1:3)-pos(j,1:3)) < r(j)**2 )then
        N_spheres = N_spheres + 1
      end if
    end do
    counter(N_spheres) = counter(N_spheres) + 1
  end do

! The sum below should work, but it doesn't. I haven't figured out why
  N_in = counter(3) + 3*counter(4) + 6*counter(5) + 10*counter(6) ! + etc.

! Get volume
  V3 = V_box / dfloat(N_k) * dfloat(N_in)

  write(*,*) N_k, V, V12+V3, V12, V3

end if ! change condition to run
!**********************************************************

end subroutine
!=================================================================================================
!=================================================================================================








!=================================================================================================
!=================================================================================================
subroutine spheres_intersection(pos1, pos2, r1, r2, V)

!***********************************************************************************
! This subroutine calculates the volume of the intersection of two spheres. The
! volume is returned in the same units (cubed) as the input radii/positions
!***********************************************************************************

  implicit none

  real*8, intent(in) :: pos1(1:3), pos2(1:3), r1, r2
  real*8, intent(out) :: V
  real*8 :: d, pi, rbig, rsmall, hbig, hsmall, Vbig, Vsmall

  pi = dacos(-1.d0)

  if(r1 > r2)then
    rbig = r1
    rsmall = r2
  else
    rbig = r2
    rsmall = r1
  end if

! Distance between spheres
  d = dsqrt(dot_product(pos1(1:3)-pos2(1:3), pos1(1:3)-pos2(1:3)))
! Cap heights
  hbig = (rsmall - rbig + d)*(rsmall + rbig - d) / 2.d0 / d
  hsmall = (rbig - rsmall + d)*(rsmall + rbig - d) / 2.d0 / d
  if(hsmall > rsmall)then
    hsmall = 2.d0*rsmall - hsmall
  end if
! Cap volumes
  Vbig = pi / 3.d0 * hbig**2 * (3.d0 * rbig - hbig)
  Vsmall = pi / 3.d0 * hsmall**2 * (3.d0 * rsmall - hsmall)

  V = 0.d0

! Cases
! Spheres are not touching:
  if( d >= rbig + rsmall )then
    V = 0.d0
! Spheres are touching and small cap is less than half a sphere
  else if( hsmall < rsmall )then
    V = Vbig + Vsmall
! Spheres are touching and small cap is more than (or equal to) half a sphere 
  else if( hsmall >= rsmall .and. d + rsmall > rbig )then
    V = 4.d0 * pi / 3.d0 * rsmall**3 - Vsmall + Vbig
! Small sphere is completely inside large sphere:
  else if( d + rsmall <= rbig)then
    V = 4.d0 * pi / 3.d0 * rsmall**3
  end if

end subroutine
!=================================================================================================
!=================================================================================================








!=================================================================================================
!=================================================================================================
subroutine init_random_seed()

!***********************************************************************************
! This subroutine was copy-pasted from the gfortran online documentation. It
! initializes the random number generator.
!***********************************************************************************

  use iso_fortran_env, only: int64
  implicit none
  integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt(8), pid
  integer(int64) :: t

   call random_seed(size = n)
   allocate(seed(n))
   ! First try if the OS provides a random number generator
   open(newunit=un, file="/dev/urandom", access="stream", &
     form="unformatted", action="read", status="old", iostat=istat)
   if (istat == 0) then
     read(un) seed
     close(un)
   else
! Fallback to XOR:ing the current time and pid. The PID is
! useful in case one launches multiple instances of the same
! program in parallel.
     call system_clock(t)
     if (t == 0) then
       call date_and_time(values=dt)
       t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
           + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
           + dt(3) * 24_int64 * 60 * 60 * 1000 &
           + dt(5) * 60 * 60 * 1000 &
           + dt(6) * 60 * 1000 + dt(7) * 1000 &
           + dt(8)
     end if
     pid = getpid()
     t = ieor(t, int(pid, kind(t)))
     do i = 1, n
       seed(i) = lcg(t)
     end do
   end if
   call random_seed(put=seed)
   contains
   ! This simple PRNG might not be good enough for real work, but is
   ! sufficient for seeding a better PRNG.
     function lcg(s)
       integer :: lcg
       integer(int64) :: s
       if (s == 0) then
         s = 104729
       else
         s = mod(s, 4294967296_int64)
       end if
         s = mod(s * 279470273_int64, 4294967291_int64)
         lcg = int(mod(s, int(huge(0), int64)), kind(0))
       end function lcg
end subroutine init_random_seed
!=================================================================================================
!=================================================================================================
