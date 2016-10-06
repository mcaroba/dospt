module constants

  implicit none

! Constants and units
  real*8 :: pi, kB, h, ps, nm, eV, amu, conv1, conv2, conv3, eVtoJ, conv4

  contains

  subroutine initialize_constants()
!******************************************************
! In this block several constants and conversion factors
! and fundamental constants are defined
  pi = dacos(-1.d0)
! Boltzmann's constant in eV/K
  kB = 8.6173324d-5
! Planck's constant in eV*ps
  h = 4.135667662d-3
! Different units given in their SI equivalent
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
  end subroutine

end module constants
