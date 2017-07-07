module thermodynamic

  use constants
  use read_input
  use trajectory
  use topology
  use calc_dos
  use partition

  implicit none

  integer :: Ng
  real*8 :: Nsg, sigma_eff
  real*8 :: nu_i(3), gas_integrand(3), solid_integrand(3), dos_subinterval(3), solid_weight(3), quadrature(3)
  real*8, allocatable :: SHSbykB_group(:,:), SHSbykB_supergroup(:,:), SRbykB_group(:), SRbykB_supergroup(:)
  real*8, allocatable :: entropy_group(:,:), entropy_supergroup(:,:), entropy_supergroup_gas(:,:)
  real*8, allocatable :: entropy_supergroup_solid(:,:), degf_supergroup_gas(:,:), degf_supergroup_solid(:,:)
  real*8, allocatable :: entropy1PT_group(:,:), entropy1PT_supergroup(:,:)
! Binary entropy of mixing variables
  real*8 :: v_dash, alpha, beta, m_frac

  contains

subroutine get_entropy()
  allocate( SHSbykB_group(1:ngroups,1:3) )
  allocate( SRbykB_group(1:ngroups) )
  allocate( entropy_group(1:ngroups, 1:3) )
  allocate( entropy1PT_group(1:ngroups, 1:3) )

  entropy_group = 0.d0
  entropy1PT_group = 0.d0


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
! FIX THIS. AGAIN (SEE ABOVE) IT MAKES LITTLE SENSE TO CALCULATE THIS FOR INDIVIDUAL GROUPS. RETHINK THIS PART OF THE CODE
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
        Nsg = ngroups_in_supergroup_eff(j)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!temp(1) = V
temp(1) = volume_supergroup(j)
        SHSbykB_supergroup(j,k) = 5.d0/2.d0 &
           + dlog(conv3 * (2.d0 * pi * mass_supergroup(j) * kB * T / h**2.d0)**(3.d0/2.d0) * &
             temp(1) * z_supergroup(j,k) / Nsg / f_supergroup(j,k)) &
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
        temp(1) = kB * ngroups_in_supergroup_eff(j) * dlog(V / volume_supergroup(j))
      else if( smixture == "mol")then
        temp(1) = kB * ngroups_in_supergroup_eff(j) * dlog(ngroups_eff / ngroups_in_supergroup_eff(j))
!     This entropy of mixing should only be used for binary mixtures !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! UNTESTED!!!!!
      else if( nsupergroups == 2 .and. smixture == "bin")then
        v_dash = volume_supergroup(2) / volume_supergroup(1) * ngroups_in_supergroup_eff(1) / ngroups_in_supergroup_eff(2)
        m_frac = ngroups_in_supergroup_eff(1) / ngroups_eff
        alpha = (dabs(v_dash - 1.d0) * (1.d0 - 2.d0 * m_frac) + v_dash + 1.d0) / 2.d0
        v_dash = 1.d0 / v_dash
        beta = (dabs(v_dash - 1.d0) * (1.d0 - 2.d0 * m_frac) + v_dash + 1.d0) / 2.d0
        if(j == 1)then
          temp(1) = kB * ngroups_in_supergroup_eff(1) * dlog( (ngroups_in_supergroup_eff(1) + &
                    ngroups_in_supergroup_eff(2)*alpha) / ngroups_in_supergroup_eff(1) )
        else if(j == 2)then
          temp(1) = kB * ngroups_in_supergroup_eff(2) * dlog( (ngroups_in_supergroup_eff(2) + &
                    ngroups_in_supergroup_eff(1)*beta) / ngroups_in_supergroup_eff(2) )
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
          temp(1) = temp(1) + weight_group(group_in_supergroup(j,j2))*3.d0
!         Total number of DoF for this group (takes constrains into account)
          temp(4) = weight_group(group_in_supergroup(j,j2))* &
                    dot_product(degf_group(group_in_supergroup(j,j2),1:3), (/1.d0, 1.d0, 1.d0/))
!         For diatomic molecules
          if(natoms_in_group(group_in_supergroup(j,j2)) == 2)then
!           We add 2 rotational DoF per linear molecule
            temp(2) = temp(2) + weight_group(group_in_supergroup(j,j2))*2.d0
!           We add the rest as vibrational DoF (if the bond is constrained this is correctly zero)
            temp(3) = temp(3) + temp(4) - weight_group(group_in_supergroup(j,j2))*5.d0
!         For molecules with 3 or more atoms
          else if(natoms_in_group(group_in_supergroup(j,j2)) > 2)then
            temp(2) = temp(2) + weight_group(group_in_supergroup(j,j2))*3.d0
            temp(3) = temp(3) + temp(4) - weight_group(group_in_supergroup(j,j2))*6.d0
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









end subroutine


end module
