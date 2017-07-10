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
!                           DoSPT v0.2                               !
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
!!!          Distribution last updated on 10 July 2017             !!!
!!!!!                                                            !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!******************************************************
! Modules used in this program
  use constants
  use read_input
  use trajectory
  use topology
  use calc_dos
  use partition
  use thermodynamic
  use good_bye
!******************************************************



!******************************************************
! This block declares all the variables needed during
! program execution
  implicit none
!******************************************************






!******************************************************
! Error code. Variable "execution_error" will store an
! integer. The code is as follows:
!
! 98 = Program could not read input file
! 97 = Program could not read trajectory file
!******************************************************







!=================================================================================================
!=================================================================================================
! This is the program workflow. All the operations performed by the code are put here
! and appear as a sequential call to subroutines contained in different modules.




!******************************************************
! Constants and units initialization (constants.f90)
  call initialize_constants()

! This module takes care of declaring the following global variables:
!
! pi, kB, h, ps, nm, eV, amu, conv1, conv2, conv3, eVtoJ, conv4
!******************************************************



!******************************************************
! Print welcome and read in all input files except for trajectory (read_input.f90)
  call print_welcome_and_read_input()

! This module takes care of declaring the following global variables:
!
! 
!******************************************************



!******************************************************
! Read trajectory file (read_trajectory.f90)
  if(execution_error == 0)then
    call allocate_trajectory_arrays()
    call read_trajectory(n, natoms, tau, L, mode, positions, velocities, estimate_vel, error, &
                         m, nmasses, mass_types, mass_types_values, volumes, volumes_temp, species, di_volumes)
    if(error)then
      execution_error = 97
    end if
  end if

! This module takes care of declaring the following global variables:
!
! 
!******************************************************



!******************************************************
! Rebuild the topology according to topology file (rebuild_topology.f90)
  if(execution_error == 0)then
    call sort_out_topology()
  end if

! This module takes care of declaring the following global variables:
!
! 
!******************************************************



!******************************************************
! Calculate the DoS and do a bunch of other stuff (calc_dos.f90)
  if(execution_error == 0)then
    call get_dos()
  end if

! This module takes care of declaring the following global variables:
!
! 
!******************************************************



!******************************************************
! Carry out the partition of the DoFs (partition.f90)
  if(execution_error == 0)then
    call dof_partition()
  end if

! This module takes care of declaring the following global variables:
!
! 
!******************************************************



!******************************************************
! Compute thermodynamic properties (thermodynamic.f90)
  if(execution_error == 0)then
    call get_entropy()
  end if

! This module takes care of declaring the following global variables:
!
! 
!******************************************************



!******************************************************
! Print timing and good by / error messages
  call print_good_bye()
!******************************************************
!=================================================================================================
!=================================================================================================


end program
