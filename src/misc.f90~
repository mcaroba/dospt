! Various utility subroutines used by the main code. They're used for
! info printing purposes mostly. These subroutines do not require an
! explicit interface.
! These routines where adpated but *not* originally written by
! Miguel A. Caro.


subroutine timestring ( string )

!*******************************************************************************
!
!! TIMESTRING writes the current YMDHMS date into a string.
!
!  Example:
!
!    STRING = 'May 31 2001   9:45:54.872 AM'
!
!  Modified:
!
!    15 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) STRING, contains the date information.
!    A character length of 40 should always be sufficient.
!
  implicit none
!
  character ( len = 8 ) ampm
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = * ) string
  character ( len = 10 ) time
  integer values(8)
  integer y
  character ( len = 5 ) zone
!
  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( string, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end






subroutine  factorize(n)

! Prime factorization code based on Ching-Kuang Shene's code:
! http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap04/factor.html
!
! All the text output stuff written by Miguel Caro

  implicit none

! We preset the array sizes to handle 100 different factors which
! should be more than enough for any practical purpose
  integer :: factor(1:100) = 0, factor_exponent(1:100) = 0
  integer :: input, divisor, count, n_diff_factors = 0, i, spaces, n
  character*64 :: temp(4), format

  input = n

  count = 0
! Remove all factors of 2
  do
    if (mod(input,2) /= 0 .or. input == 1) exit
    count = count + 1
    n_diff_factors = 1
    factor(1) = 2
    factor_exponent(1) = factor_exponent(1) + 1
    input = input / 2
  end do

! Remove odd factors now
  divisor = 3
  do
!   If a factor is too large, exit and done
    if (divisor > input) exit
!   Try this factor repeatedly
    do
      if (mod(input, divisor) /= 0 .or. input == 1)  exit
      count = count + 1
      if(divisor /= factor(n_diff_factors) )then
        n_diff_factors = n_diff_factors + 1
        factor(n_diff_factors) = divisor
      end if
      factor_exponent(n_diff_factors) = factor_exponent(n_diff_factors) + 1
      input = input / divisor
    end do
!   Move to next odd number
    divisor = divisor + 2
  end do


! Print error message if we find factors larger than 10
  if(factor(n_diff_factors) > 50)then
    write(*,*)'          *****************            |'
    write(*,*)'WARNING:                               |'
    write(*,*)'The prime factorization of your number |'
    write(*,*)'of trajectory snapshots yiels very high|'
    write(*,*)'factors. This can make the calculation |'
    write(*,*)'of the DoS (which relies on Fast       |'
    write(*,*)'Fourier Transform) grossly inefficient.|'
    write(*,*)'                                       |'
    write(*,*)'The factorization currently yields:    |'
    write(*,*)'                                       |'
    write(*,'(1X,A)',advance='no')'n ='

    write(temp(1),*) factor(1)
    write(temp(2),*) factor_exponent(1)
    write(temp(3),'(A)') trim(adjustl(temp(1))) // '^' // trim(adjustl(temp(2))) 

    count = 0
    do i=2,n_diff_factors
      write(temp(1),*) factor(i)
      write(temp(2),*) factor_exponent(i)
      write(temp(4),'(A)') '* ' // trim(adjustl(temp(1))) // '^' // trim(adjustl(temp(2)))
      if( len(trim(adjustl(temp(3)))) + len(trim(adjustl(temp(4)))) > 33 )then
        spaces = 35 - len(trim(adjustl(temp(3))))
        if(count == 0)then
          write(format,'(A,I2,A)') "(1X,A,", spaces, "X,A)"
          write(*,format) trim(adjustl(temp(3))), '|'
          count = 1
        else
          write(format,'(A,I2,A)') "(5X,A,", spaces, "X,A)"
          write(*,format) trim(adjustl(temp(3))), '|'
        end if
        write(temp(3),'(5X,A)') ' ' // trim(adjustl(temp(4))) 
      else
        write(temp(3),'(A)') trim(adjustl(temp(3))) // ' ' // trim(adjustl(temp(4)))
      end if
    end do

    spaces = 35 - len(trim(adjustl(temp(3))))
    if(count == 0)then
      write(format,'(A,I2,A)') "(1X,A,", spaces, "X,A)"
      write(*,format) trim(adjustl(temp(3))), '|'
      count = 1
    else
      write(format,'(A,I2,A)') "(5X,A,", spaces, "X,A)"
      write(*,format) trim(adjustl(temp(3))), '|'
    end if
    write(*,*)'                                       |'
    write(*,*)'If you notice low performance during   |'
    write(*,*)'the DoS calculation step, try reducing |'
    write(*,*)'or increasing n so that large prime    |'
    write(*,*)'factors are removed.                   |'
    write(*,*)'          *****************            |'
    write(*,*)'                                       |'
  end if

end subroutine






subroutine get_eigenvalues_3x3(A,eig)

! I copied this routine from Wikipedia and modified parts of it
! It assumes matrix A is symmetric

  implicit none

  real*8 :: A(3,3), eig(3), I(3,3) = 0.d0, B(3,3)
  real*8 :: p1, q, p2, p, r, detB, phi, pi

  pi = dacos(-1.d0)

  I(1,1) = 1.d0
  I(2,2) = 1.d0
  I(3,3) = 1.d0

  p1 = A(1,2)**2 + A(1,3)**2 + A(2,3)**2
  if (p1 == 0.) then
    ! A is diagonal.
    eig(1) = A(1,1)
    eig(2) = A(2,2)
    eig(3) = A(3,3)
  else
    q = (A(1,1) + A(2,2) + A(3,3))/3.d0
    p2 = (A(1,1) - q)**2 + (A(2,2) - q)**2 + (A(3,3) - q)**2 + 2.d0 * p1
    p = dsqrt(p2 / 6.d0)
    B = (1.d0 / p) * (A - q * I)
    detB = -B(1,3)**2*B(2,2) + 2.d0*B(1,2)*B(1,3)*B(2,3) - B(1,1)*B(2,3)**2 - B(1,2)**2*B(3,3) + B(1,1)*B(2,2)*B(3,3)
    r = detB / 2.d0
  end if

  ! In exact arithmetic for a symmetric matrix  -1 <= r <= 1
  ! but computation error can leave it slightly outside this range.
  if (r <= -1.) then
    phi = pi / 3.d0
  else if (r >= 1.) then
    phi = 0.d0
  else
    phi = dacos(r) / 3.d0
  end if

  ! the eigenvalues satisfy eig3 <= eig2 <= eig1
  eig(1) = q + 2.d0 * p * dcos(phi)
  eig(3) = q + 2.d0 * p * dcos(phi + (2.d0*pi/3.d0))
  ! since trace(A) = eig1 + eig2 + eig3
  eig(2) = 3.d0 * q - eig(1) - eig(3)

end subroutine
