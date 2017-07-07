module misc2

  contains

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









subroutine cross_product(u,v,cross)

  implicit none

  real*8 :: u(1:3), v(1:3), cross(1:3)

  cross(1:3) = (/ u(2)*v(3)-u(3)*v(2), &
                 -u(1)*v(3)+u(3)*v(1), &
                  u(1)*v(2)-u(2)*v(1) /)

end subroutine

end module
