!**********************************
! This module calculates energies and forces
! associated with intra-molecular bond, angle,
! and dihedral potentials
!*********************************

module bonded_interactions
  use routines


contains

  !*****************************************
  ! this subroutine returns the intra-molecular
  ! bond, angle, dihedral energy and force of the entire system
  !******************************************
  subroutine intra_molecular_energy_force( system_data, molecule_data, atom_data )
    use global_variables
    use omp_lib
    type(system_data_type) , intent(inout)  :: system_data
    type(molecule_data_type), dimension(:), intent(in) :: molecule_data
    type(atom_data_type) , intent(inout) :: atom_data

    !******** this is a local data structure with pointers that will be set
    ! to subarrays of atom_data arrays for the specific atoms in the molecule
    type(single_molecule_data_type) :: single_molecule_data

    integer :: i_mole, i_mole_type
    real*8  :: E_bond_local, E_angle_local , E_dihedral_local
    real*8  :: E_bond, E_angle , E_dihedral

    ! can't use subobjects in OpenMP reduction clause, so use temporary local variables here
    E_bond=0; E_angle=0; E_dihedral=0

    ! note that not all molecules have dihedrals, and those that do are probably
    ! in continuous indices, so don't split loop into big chunks here

    !    call OMP_SET_NUM_THREADS(n_threads)
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(n_threads, atom_data, molecule_data, system_data ) REDUCTION(+:E_bond, E_angle, E_dihedral)
    !$OMP DO SCHEDULE(DYNAMIC,1)
    do i_mole = 1, system_data%n_mole

       ! this is used to look up bonded force field parameters for this molecule
       i_mole_type = molecule_data(i_mole)%molecule_type_index

       ! set pointers for single_molecule_data to target molecule subarrays of atom_data
       call return_molecule_block( single_molecule_data , molecule_data(i_mole)%n_atom, molecule_data(i_mole)%atom_index, atom_xyz=atom_data%xyz, atom_force=atom_data%force, atom_type_index=atom_data%atom_type_index )

       ! bond energy, force for molecule
       call intra_molecular_bond_energy_force( E_bond_local, single_molecule_data, i_mole_type )

       ! angle energy, force for molecule
       call intra_molecular_angle_energy_force( E_angle_local, single_molecule_data, i_mole_type )

       ! dihedral energy, force for molecule
       call intra_molecular_dihedral_energy_force( E_dihedral_local, single_molecule_data, i_mole_type ) 

       E_bond     =   E_bond + E_bond_local
       E_angle    =   E_angle + E_angle_local
       E_dihedral =   E_dihedral + E_dihedral_local

       call dissociate_single_molecule_data(single_molecule_data)

    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    system_data%E_bond = E_bond
    system_data%E_angle = E_angle
    system_data%E_dihedral = E_dihedral


  end subroutine intra_molecular_energy_force



  !*************************************
  ! this subroutine returns the intra-molecular
  ! bond energy of the entire system
  ! note there is no pbc image check, because
  ! as is convention in this program, molecules
  ! are not split up over periodic boundaries
  !*************************************
  subroutine intra_molecular_bond_energy_force( E_bond, single_molecule_data , i_mole_type ) 
    use global_variables
    real*8, intent(out) :: E_bond
    type(single_molecule_data_type), intent(inout):: single_molecule_data
    integer, intent(in) :: i_mole_type

    integer :: i_atom_type, j_atom_type, i_atom, j_atom, i_bond
    real*8, dimension(3) :: r_ij, force_ij
    real*8 :: r_mag, E_bond_ij

    E_bond = 0d0
    ! loop over bond list for this molecule
    do i_bond=1, size(molecule_type_data(i_mole_type)%bond_list)
       i_atom = molecule_type_data(i_mole_type)%bond_list(i_bond)%i_atom
       j_atom = molecule_type_data(i_mole_type)%bond_list(i_bond)%j_atom

       ! bond between these two atoms
       i_atom_type = single_molecule_data%atom_type_index(i_atom)
       j_atom_type = single_molecule_data%atom_type_index(j_atom)

       r_ij(:) = single_molecule_data%xyz(:,i_atom) - single_molecule_data%xyz(:,j_atom)
       r_mag = dsqrt( dot_product( r_ij , r_ij ) )

       call pairwise_bond_energy_force( E_bond_ij, force_ij, i_atom_type, j_atom_type, r_ij , r_mag )

       E_bond = E_bond + E_bond_ij 

       single_molecule_data%force(:,i_atom) = single_molecule_data%force(:,i_atom) + force_ij(:)
       single_molecule_data%force(:,j_atom) = single_molecule_data%force(:,j_atom) - force_ij(:)

    enddo

  end subroutine intra_molecular_bond_energy_force



  !*************************************
  ! this subroutine calculates the pairwise bond energy
  ! and force between two atomtypes.  Currently, 
  ! either a Harmonic(1), G-96 fourth power(2), or Morse(3) function type is employed
  !*************************************
  subroutine pairwise_bond_energy_force( E_bond_ij, force_ij, i_atom_type, j_atom_type, r_ij , r_mag )
    use global_variables
    real*8, intent(out) :: E_bond_ij
    real*8,dimension(3),intent(out) :: force_ij
    integer,intent(in) :: i_atom_type, j_atom_type
    real*8,dimension(3), intent(in) :: r_ij
    real*8, intent(in) :: r_mag

    integer :: bondtype
    real*8 :: b0 , kb , D, beta

    bondtype = atype_bond_type(i_atom_type,j_atom_type)


    Select Case( bondtype )
    Case(1)
       ! Harmonic
       b0 = atype_bond_parameter(i_atom_type,j_atom_type,1)
       kb = atype_bond_parameter(i_atom_type,j_atom_type,2)
       E_bond_ij = 0.5d0 * kb * ( r_mag - b0 ) ** 2 
       force_ij(:) = - kb * ( r_mag - b0 ) * r_ij(:) / r_mag
    Case(2)
       ! This is the G-96 fourth power potential
       b0 = atype_bond_parameter(i_atom_type,j_atom_type,1)
       kb = atype_bond_parameter(i_atom_type,j_atom_type,2)
       E_bond_ij = 0.25d0 * kb * ( r_mag**2 - b0**2 ) ** 2 
       force_ij(:) = - kb * ( r_mag**2 - b0**2 ) * r_ij(:)
    Case(3)
       ! Morse
       D = atype_bond_parameter(i_atom_type,j_atom_type,1)
       beta = atype_bond_parameter(i_atom_type,j_atom_type,2)
       b0 = atype_bond_parameter(i_atom_type,j_atom_type,3)

       E_bond_ij = D * ( 1d0 - exp(-beta*(r_mag-b0)) ) ** 2
       force_ij(:) = -2d0 * D * beta * exp(-beta*(r_mag-b0)) * ( 1d0 - exp(-beta*(r_mag-b0)) )* r_ij(:) / r_mag

    Case default
       write(*,*) "bond type isn't implemented! Please select either 1 (harmonic)"
       write(*,*) "or 3 (Morse)"
       stop

    End Select

  end subroutine pairwise_bond_energy_force



  !*************************************
  ! this subroutine returns the intra-molecular
  ! bond energy of the entire system
  ! note there is no pbc image check, because
  ! as is convention in this program, molecules
  ! are not split up over periodic boundaries
  !*************************************
  subroutine intra_molecular_angle_energy_force( E_angle, single_molecule_data, i_mole_type ) 
    use global_variables
    real*8, intent(out) :: E_angle
    type(single_molecule_data_type), intent(inout) :: single_molecule_data
    integer, intent(in) :: i_mole_type

    integer :: i_atom_type, j_atom_type, k_atom_type, i_atom, j_atom, k_atom, i_angle
    real*8, dimension(3) :: r_ij, r_kj, f_ij , f_kj
    real*8 :: E_angle_ijk

    E_angle=0d0

    ! loop over angle list for this molecule
    do i_angle=1, size(molecule_type_data(i_mole_type)%angle_list)
       i_atom = molecule_type_data(i_mole_type)%angle_list(i_angle)%i_atom
       j_atom = molecule_type_data(i_mole_type)%angle_list(i_angle)%j_atom
       k_atom = molecule_type_data(i_mole_type)%angle_list(i_angle)%k_atom

       i_atom_type = single_molecule_data%atom_type_index(i_atom)
       j_atom_type = single_molecule_data%atom_type_index(j_atom)
       k_atom_type = single_molecule_data%atom_type_index(k_atom)

       r_ij(:) = single_molecule_data%xyz(:,i_atom) - single_molecule_data%xyz(:,j_atom)
       r_kj(:) = single_molecule_data%xyz(:,k_atom) - single_molecule_data%xyz(:,j_atom)

       call trimer_angle_energy_force( E_angle_ijk, f_ij, f_kj, i_atom_type, j_atom_type, k_atom_type , r_ij, r_kj)

       E_angle = E_angle + E_angle_ijk
       single_molecule_data%force(:,i_atom) = single_molecule_data%force(:,i_atom) + f_ij(:)
       single_molecule_data%force(:,k_atom) = single_molecule_data%force(:,k_atom) + f_kj(:)
       single_molecule_data%force(:,j_atom) = single_molecule_data%force(:,j_atom)  - f_ij(:) - f_kj(:)
    end do

  end subroutine intra_molecular_angle_energy_force



  !**************************************
  ! this computes the angle potential between three atoms
  !**************************************
  subroutine trimer_angle_energy_force( E_angle_ijk, f_ij, f_kj, i_atom_type, j_atom_type, k_atom_type , r_ij, r_kj )
    use global_variables
    integer, intent(in) :: i_atom_type, j_atom_type, k_atom_type
    real*8,dimension(3), intent(in) :: r_ij, r_kj
    real*8,dimension(3), intent(out) :: f_ij, f_kj
    real*8, intent(out)  :: E_angle_ijk

    integer :: angletype
    real*8 :: rij_mag, rkj_mag, theta, th0 , cth, fac, cosine, cosine0
    real*8,parameter :: small=1D-4

    rij_mag = dsqrt( dot_product( r_ij , r_ij ) )
    rkj_mag = dsqrt( dot_product( r_kj , r_kj ) )
    cosine = dot_product( r_ij , r_kj ) / rij_mag / rkj_mag

    angletype = atype_angle_type(i_atom_type,j_atom_type,k_atom_type)

    Select Case( angletype )
    Case(1)
       ! this is harmonic angle potential

       ! be careful about numerical instabilities here! if cosine is 1.0000001, acos will return NaN
       if ( cosine < -0.999999999d0 ) then
          theta = constants%pi
       elseif ( cosine > 0.999999999d0 ) then
          theta = 0d0
       else
          theta = acos( cosine )
       endif

       th0 = atype_angle_parameter(i_atom_type, j_atom_type, k_atom_type,1)
       cth = atype_angle_parameter(i_atom_type, j_atom_type, k_atom_type,2)

       E_angle_ijk = 0.5d0 * cth * ( theta - th0 ) ** 2
       ! force

       ! don't divide by zero here, because of the large force constant,
       ! the only time cosine will be 1 is if theta=theta0=pi
       if ( abs(theta - th0 ) < small ) then
          fac=0d0
       else
          fac = cth * ( theta - th0 ) / dsqrt(1d0 - cosine**2)
       end if

       f_ij(:) = fac * ( r_kj(:) / rij_mag / rkj_mag - cosine * r_ij(:) / rij_mag ** 2 )
       f_kj(:) = fac * ( r_ij(:) / rij_mag / rkj_mag - cosine * r_kj(:) / rkj_mag ** 2 )

    Case(2)
       ! this is cosine based angle potential
       th0 = atype_angle_parameter(i_atom_type, j_atom_type, k_atom_type,1)
       cth = atype_angle_parameter(i_atom_type, j_atom_type, k_atom_type,2)
       ! th0 should be in radians here
       cosine0 = cos( th0 )

       E_angle_ijk = 0.5d0 * cth * ( cosine - cosine0 ) ** 2
       ! force
       fac = -cth * ( cosine - cosine0 ) 
       f_ij(:) = fac * ( r_kj(:) / rij_mag / rkj_mag - cosine * r_ij(:) / rij_mag ** 2 )
       f_kj(:) = fac * ( r_ij(:) / rij_mag / rkj_mag - cosine * r_kj(:) / rkj_mag ** 2 )

    case default
       write(*,*) "requested angle type potential not implemented"
       stop
    End Select

  end subroutine trimer_angle_energy_force




  !*************************************
  ! this subroutine returns the intra-molecular
  ! dihedral energy of the entire system
  ! note there is no pbc image check, because
  ! as is convention in this program, molecules
  ! are not split up over periodic boundaries
  !*************************************
  subroutine intra_molecular_dihedral_energy_force( E_dihedral, single_molecule_data, i_mole_type )
    use global_variables
    real*8, intent(out) :: E_dihedral
    type(single_molecule_data_type), intent(inout) :: single_molecule_data
    integer, intent(in) :: i_mole_type

    integer :: i_atom_type, j_atom_type, k_atom_type, l_atom_type, i_atom, j_atom, k_atom, l_atom, i_dihedral
    real*8,dimension(3) :: r_ji, r_kj, r_lk, f_ji, f_kj, f_lk
    real*8 :: E_dihedral_ijkl


    E_dihedral=0d0
    ! loop over dihedral list for this molecule
    do i_dihedral=1,size(molecule_type_data(i_mole_type)%dihedral_list)
       i_atom = molecule_type_data(i_mole_type)%dihedral_list(i_dihedral)%i_atom
       j_atom = molecule_type_data(i_mole_type)%dihedral_list(i_dihedral)%j_atom
       k_atom = molecule_type_data(i_mole_type)%dihedral_list(i_dihedral)%k_atom
       l_atom = molecule_type_data(i_mole_type)%dihedral_list(i_dihedral)%l_atom


       i_atom_type = single_molecule_data%atom_type_index(i_atom)
       j_atom_type = single_molecule_data%atom_type_index(j_atom)
       k_atom_type = single_molecule_data%atom_type_index(k_atom)
       l_atom_type = single_molecule_data%atom_type_index(l_atom)

       r_ji(:) = single_molecule_data%xyz(:,j_atom) - single_molecule_data%xyz(:,i_atom)
       r_kj(:) = single_molecule_data%xyz(:,k_atom) - single_molecule_data%xyz(:,j_atom)
       r_lk(:) = single_molecule_data%xyz(:,l_atom) - single_molecule_data%xyz(:,k_atom)

       call quartet_dihedral_energy_force( E_dihedral_ijkl, f_ji, f_kj, f_lk, i_atom_type, j_atom_type, k_atom_type , l_atom_type, r_ji, r_kj, r_lk )

       E_dihedral = E_dihedral + E_dihedral_ijkl

       single_molecule_data%force(:,i_atom) = single_molecule_data%force(:,i_atom) - f_ji(:)
       single_molecule_data%force(:,j_atom) = single_molecule_data%force(:,j_atom) + f_ji(:) - f_kj(:)
       single_molecule_data%force(:,k_atom) = single_molecule_data%force(:,k_atom) + f_kj(:) - f_lk(:)
       single_molecule_data%force(:,l_atom) = single_molecule_data%force(:,l_atom) + f_lk(:)
        
    enddo
  end subroutine intra_molecular_dihedral_energy_force



  !**********************************
  ! this subroutine computes the dihedral energy
  ! and force for the quartet of atoms i,j,k,l
  !**********************************
  subroutine quartet_dihedral_energy_force( E_dihedral_ijkl, f_ji, f_kj, f_lk, i_atom_type, j_atom_type, k_atom_type , l_atom_type, r_ji, r_kj, r_lk )
    use global_variables
    real*8, intent(out) :: E_dihedral_ijkl
    real*8,dimension(3), intent(out) :: f_ji, f_kj, f_lk
    real*8,dimension(3), intent(in) :: r_ji, r_kj, r_lk
    integer, intent(in) :: i_atom_type, j_atom_type, k_atom_type, l_atom_type

    real*8 :: rji_mag, rkj_mag, rlk_mag, rji_mag2, rkj_mag2, rlk_mag2, n_mult
    real*8 :: a_dot_b , a_dot_a , b_dot_b, xi, xi0, kxi, cosine, cosine2
    real*8 :: sqrt_a_dot_a, sqrt_b_dot_b
    real*8 :: dot_rkj_rji , dot_rlk_rkj , dot_rlk_rji, fac
    real*8 :: c0, c1, c2, c3, c4, c5
    real*8, dimension(3) :: d_a_dot_b_drji , d_a_dot_b_drkj , d_a_dot_b_drlk
    real*8, dimension(3) :: d_a_dot_a_drji , d_a_dot_a_drkj , d_a_dot_a_drlk
    real*8, dimension(3) :: d_b_dot_b_drji , d_b_dot_b_drkj , d_b_dot_b_drlk
    integer :: dihedraltype, shift
    real*8,parameter :: small=1D-4

    rji_mag2 = dot_product( r_ji , r_ji )
    rkj_mag2 = dot_product( r_kj , r_kj )
    rlk_mag2 = dot_product( r_lk , r_lk )
    rji_mag = dsqrt( rji_mag2 )
    rkj_mag = dsqrt( rkj_mag2 )
    rlk_mag = dsqrt( rlk_mag2 )
    dot_rkj_rji = dot_product( r_kj, r_ji )
    dot_rlk_rkj = dot_product( r_lk, r_kj )
    dot_rlk_rji = dot_product( r_lk, r_ji )

    ! we don't use cross products here, because we use two triple product formulas
    ! to convert the dot product of two cross products into terms that only
    ! depend on dot products

    ! define alpha as r_ji x r_kj and beta as r_kj x r_lk
    ! then we want cos(xi) = alpha .dot. beta / mag(alpha) / mag_beta

    !alpha .dot. beta
    a_dot_b = dot_rkj_rji * dot_rlk_rkj - dot_rlk_rji * rkj_mag2
    ! alpha .dot. alpha
    a_dot_a = rji_mag2 * rkj_mag2 - dot_rkj_rji ** 2
    ! beta .dot. beta
    b_dot_b = rlk_mag2 * rkj_mag2 - dot_rlk_rkj ** 2

    ! need gradients of a_dot_b , a_dot_a, and b_dot_b w.r.t to rji, rkj, rlk
    ! calculate these
    d_a_dot_b_drji(:) = r_kj(:) * dot_rlk_rkj - r_lk(:) * rkj_mag2
    d_a_dot_b_drkj(:) = r_ji(:) * dot_rlk_rkj + dot_rkj_rji * r_lk(:) - dot_rlk_rji * 2d0 * r_kj(:)
    d_a_dot_b_drlk(:) = dot_rkj_rji * r_kj(:) - r_ji(:) * rkj_mag2

    d_a_dot_a_drji(:) = rkj_mag2 * 2d0 * r_ji(:) - 2d0 * dot_rkj_rji * r_kj(:)
    d_a_dot_a_drkj(:) = rji_mag2 * 2d0 * r_kj(:) - 2d0 * dot_rkj_rji * r_ji(:)
    d_a_dot_a_drlk(:) = 0d0

    d_b_dot_b_drji(:) = 0d0
    d_b_dot_b_drkj(:) = rlk_mag2 * 2d0 * r_kj(:) - 2d0 * dot_rlk_rkj * r_lk(:)
    d_b_dot_b_drlk(:) = rkj_mag2 * 2d0 * r_lk(:) - 2d0 * dot_rlk_rkj * r_kj(:)

    ! now get dihedral angle xi
    sqrt_a_dot_a = dsqrt( a_dot_a )
    sqrt_b_dot_b = dsqrt( b_dot_b )
    cosine = a_dot_b / sqrt_a_dot_a  / sqrt_b_dot_b

    ! be careful about numerical instabilities here! if cosine is 1.0000001, acos will return NaN
    if ( cosine < -0.999999999d0 ) then
       xi = constants%pi
    elseif ( cosine > 0.999999999d0 ) then
       xi = 0d0
    else
       xi = acos( cosine )
    endif



    dihedraltype = atype_dihedral_type(i_atom_type,j_atom_type,k_atom_type,l_atom_type)

    Select Case( dihedraltype )
    Case(1)
       ! proper dihedral
       ! note that the angle of the proper dihedral corresponds to zero in the cis configuration
       ! our definition of alpha and beta is consistent with this

       ! xi0 should be in radians here
       xi0 = atype_dihedral_parameter(i_atom_type, j_atom_type, k_atom_type,l_atom_type,1)
       kxi = atype_dihedral_parameter(i_atom_type, j_atom_type, k_atom_type,l_atom_type,2)
       n_mult = atype_dihedral_parameter(i_atom_type, j_atom_type, k_atom_type,l_atom_type,3)

       E_dihedral_ijkl = kxi * ( 1d0 + cos( n_mult * xi - xi0 ) )

       ! now force

       cosine2 = cosine**2
       ! don't divide by zero here
       if ( abs(cosine2 - 1d0 ) < small ) then
          ! make sure shift is zero or 180, so that sine is zero
          if ( (abs( xi0 ) < small ) .or. (abs( xi0 - constants%pi ) < small ) ) then
             fac = 0d0
          else
             ! sine isn't zero, undefined force
             stop "undefined dihedral force"
          endif
       else
          fac = kxi * -sin( n_mult * xi - xi0 ) * n_mult / dsqrt(1d0 - cosine2)
       end if

       f_ji(:) = fac * ( d_a_dot_b_drji(:) / sqrt_a_dot_a  / sqrt_b_dot_b  &
            & - 0.5d0 * a_dot_b / a_dot_a**1.5d0 / sqrt_b_dot_b  * d_a_dot_a_drji(:) &
            & - 0.5d0 * a_dot_b / sqrt_a_dot_a  / b_dot_b**1.5d0 * d_b_dot_b_drji(:) )

       f_kj(:) = fac * ( d_a_dot_b_drkj(:) / sqrt_a_dot_a  / sqrt_b_dot_b  &
            & - 0.5d0 * a_dot_b / a_dot_a**1.5d0 / sqrt_b_dot_b  * d_a_dot_a_drkj(:) &
            & - 0.5d0 * a_dot_b / sqrt_a_dot_a  / b_dot_b**1.5d0 * d_b_dot_b_drkj(:) )

       f_lk(:) = fac * ( d_a_dot_b_drlk(:) / sqrt_a_dot_a  / sqrt_b_dot_b  &
            & - 0.5d0 * a_dot_b / a_dot_a**1.5d0 / sqrt_b_dot_b * d_a_dot_a_drlk(:) &
            & - 0.5d0 * a_dot_b / sqrt_a_dot_a  / b_dot_b**1.5d0 * d_b_dot_b_drlk(:) )


    Case(2)
       ! improper dihedral

       ! acos has a range from 0 to pi, however for the dihedral potential, we want the
       ! smallest angle between the planes, which has a range from 0 to pi / 2
       ! so if xi > pi / 2 , we need to shift it by 180, (and then take absolute values)
       ! this effects the forces as follows: if cosine was decreasing (more negative) in the range
       ! pi / 2 < xi < pi , then xi was increasing. but an increasing xi in this range
       ! leads to a decreasing xi after the shift.  Therefore, if angle was shifted,
       ! multiply forces by -1.
       if ( xi > ( constants%pi / 2d0 ) ) then
          xi = abs( xi - constants%pi )
          shift=1
       else
          shift=0
       end if

       ! xi0 should be in radians here
       xi0 = atype_dihedral_parameter(i_atom_type, j_atom_type, k_atom_type,l_atom_type,1)
       kxi = atype_dihedral_parameter(i_atom_type, j_atom_type, k_atom_type,l_atom_type,2)

       E_dihedral_ijkl = 0.5d0 * kxi * ( xi - xi0 ) ** 2

       ! now force

       ! don't divide by zero here, because of the large force constant for impropers,
       ! the only time cosine will be 1 is if xi=xi0=0
       if ( abs(xi - xi0 ) < small ) then
          fac=0d0
       else
          fac = kxi * ( xi - xi0 ) / dsqrt(1d0 - cosine**2)
       end if


       f_ji(:) = fac * ( d_a_dot_b_drji(:) / sqrt_a_dot_a  / sqrt_b_dot_b  &
            & - 0.5d0 * a_dot_b / a_dot_a**1.5d0 / sqrt_b_dot_b  * d_a_dot_a_drji(:) &
            & - 0.5d0 * a_dot_b / sqrt_a_dot_a  / b_dot_b**1.5d0 * d_b_dot_b_drji(:) )

       f_kj(:) = fac * ( d_a_dot_b_drkj(:) / sqrt_a_dot_a  / sqrt_b_dot_b  &
            & - 0.5d0 * a_dot_b / a_dot_a**1.5d0 / sqrt_b_dot_b  * d_a_dot_a_drkj(:) &
            & - 0.5d0 * a_dot_b / sqrt_a_dot_a  / b_dot_b**1.5d0 * d_b_dot_b_drkj(:) )

       f_lk(:) = fac * ( d_a_dot_b_drlk(:) / sqrt_a_dot_a  / sqrt_b_dot_b  &
            & - 0.5d0 * a_dot_b / a_dot_a**1.5d0 / sqrt_b_dot_b  * d_a_dot_a_drlk(:) &
            & - 0.5d0 * a_dot_b / sqrt_a_dot_a  / b_dot_b**1.5d0 * d_b_dot_b_drlk(:) )

       ! if angle was shifted, change sign of forces
       if ( shift == 1 ) then
          f_ji = -f_ji
          f_kj = -f_kj
          f_lk = -f_lk
       end if
    Case(3)
       !Calculate the Ryckaert-Bellemans dihedral potential.
    
       !6 Parameters for the Ryckaert-Bellemans cosine expansion
       c0 = atype_dihedral_parameter(i_atom_type, j_atom_type,k_atom_type,l_atom_type,1)
       c1 = atype_dihedral_parameter(i_atom_type, j_atom_type,k_atom_type,l_atom_type,2)
       c2 = atype_dihedral_parameter(i_atom_type, j_atom_type,k_atom_type,l_atom_type,3)
       c3 = atype_dihedral_parameter(i_atom_type, j_atom_type,k_atom_type,l_atom_type,4)
       c4 = atype_dihedral_parameter(i_atom_type, j_atom_type,k_atom_type,l_atom_type,5)
       c5 = atype_dihedral_parameter(i_atom_type, j_atom_type,k_atom_type,l_atom_type,6)

       !Dihedral definition for Ryckaert-Bellemans is opposite that of the standard
       !definition (cis dihedral angle is 0) and the angle is +180 compared to the 
       !standard potential. This makes cosine equal to negative cosine. We adjust 
       !this by putting a negative term in front of the odd-powered terms for the 
       !energy function.


       E_dihedral_ijkl = c0 - c1*cosine + c2*cosine**2 - c3*cosine**3 + &
       c4*cosine**4 - c5*cosine**5

       !force

       fac = c1 - 2d0*c2*cosine + 3d0*c3*cosine**2-4d0*c4* &
            & cosine**3 + 5d0*c5*cosine**4

       f_ji(:) = fac * ( d_a_dot_b_drji(:) / sqrt_a_dot_a  / sqrt_b_dot_b  &
            & - 0.5d0 * a_dot_b / a_dot_a**1.5d0 / sqrt_b_dot_b  *d_a_dot_a_drji(:) &
            & - 0.5d0 * a_dot_b / sqrt_a_dot_a  / b_dot_b**1.5d0 *d_b_dot_b_drji(:) )

       f_kj(:) = fac * ( d_a_dot_b_drkj(:) / sqrt_a_dot_a  / sqrt_b_dot_b  &
            & - 0.5d0 * a_dot_b / a_dot_a**1.5d0 / sqrt_b_dot_b  *d_a_dot_a_drkj(:) &
            & - 0.5d0 * a_dot_b / sqrt_a_dot_a  / b_dot_b**1.5d0 *d_b_dot_b_drkj(:) )
       
       f_lk(:) = fac * ( d_a_dot_b_drlk(:) / sqrt_a_dot_a  / sqrt_b_dot_b  &
            & - 0.5d0 * a_dot_b / a_dot_a**1.5d0 / sqrt_b_dot_b  *d_a_dot_a_drlk(:) &
            & - 0.5d0 * a_dot_b / sqrt_a_dot_a  / b_dot_b**1.5d0 *d_b_dot_b_drlk(:) )


    End Select

  end subroutine quartet_dihedral_energy_force

  !*****************************************
  ! This subroutine generates an exclusions list
  ! "pair_exclusions" for every molecule type
  ! for intra-molecular non-bonded interactions
  ! based on the setting of n_exclusions global variable parameter
  ! atom pairs connected by n_exclusions bonds or less will be
  ! excluded
  ! exclusions will be marked with a "1"
  !
  ! in addition to exclusions, we may have a special treatment
  ! of 1-4 interactions, and therefore we need to identify 1-4 interactions
  ! as well as exclusion.  1-4 interactions will be marked with
  ! a "2"
  !
  ! normal interactions are marked with a "0"
  !
  ! IMPORTANT:  Do not zero the molecule_exclusions array in this subroutine,
  ! because it's possible that explicit exclusions have been read in from the
  ! topology file.
  !*****************************************
  subroutine generate_intramolecular_exclusions
    use global_variables
    integer :: i_molecule_type, i_bond, i_atom, j_atom, n_bonds, max_search
    integer, dimension(:), allocatable :: bond_trajectory
    integer, dimension(:,:), allocatable :: bond_lookup_table

    write(*,*) ""
    write(*,*) "automatically generating exclusions for intra-molecular non-bonded"
    write(*,*) "interactions that are within ", n_exclusions, " bonds away"
    write(*,*) ""

    ! if 1-4 are not excluded ( 3 bonds away ), we want to label these specially
    max_search = max( n_exclusions , 3 )

    ! this array keeps a list of all atoms along the bonding trajectory, so that we don't loop back over an atom
    allocate( bond_trajectory( max_search+1 ) )

    do i_molecule_type=1, n_molecule_type
       bond_trajectory=0

       ! make temporarary bond lookup table for this molecule type
       allocate( bond_lookup_table( molecule_type_data(i_molecule_type)%n_atom , molecule_type_data(i_molecule_type)%n_atom ) )
       bond_lookup_table=0
       do i_bond=1, size( molecule_type_data(i_molecule_type)%bond_list )
          i_atom = molecule_type_data(i_molecule_type)%bond_list(i_bond)%i_atom
          j_atom = molecule_type_data(i_molecule_type)%bond_list(i_bond)%j_atom
          bond_lookup_table(i_atom,j_atom)=1
          bond_lookup_table(j_atom,i_atom)=1
       end do

       do i_atom=1, molecule_type_data(i_molecule_type)%n_atom
          ! fill in self-exclusion
          molecule_type_data(i_molecule_type)%pair_exclusions(i_atom, i_atom) = 1

          ! this subroutine generates exclusions by searching over bonded neighbors recursively
          n_bonds=1
          bond_trajectory(1) = i_atom
          call search_bonds_recursive( i_molecule_type, i_atom, i_atom, n_bonds , bond_trajectory, bond_lookup_table, max_search )
       end do

       deallocate( bond_lookup_table )
    end do

  end subroutine generate_intramolecular_exclusions

  !******************************************
  ! this subroutine finds bonded neighbors of j_atom, and if those neighbors
  ! are n_bonds away from i_atom, then an exclusion is generated between that
  ! neighbor and i_atom
  !******************************************
  recursive subroutine search_bonds_recursive( i_molecule_type, i_atom, j_atom, n_bonds_in , bond_trajectory, bond_lookup_table, max_search )
    use global_variables
    integer, intent(in) :: i_molecule_type, i_atom, j_atom, n_bonds_in, max_search
    integer, dimension(:), intent(in) :: bond_trajectory
    integer, dimension(:,:), intent(in) :: bond_lookup_table

    integer, dimension(size(bond_trajectory)) :: update_bond_trajectory
    integer :: local_atom, n_bonds_out, flag

    update_bond_trajectory = bond_trajectory

    do local_atom=1, molecule_type_data(i_molecule_type)%n_atom
       ! if there's a bond between this atom and j_atom
       if ( bond_lookup_table(j_atom, local_atom) == 1 ) then
          ! make sure we are not doubling back on bonding trajectory
          flag = check_list( bond_trajectory , n_bonds_in , local_atom )
          if ( flag == 0 ) then
             if ( ( n_bonds_in == 3 ) .and. ( n_exclusions < 3 ) .and. ( molecule_type_data(i_molecule_type)%pair_exclusions(i_atom, local_atom) /= 1 ) ) then
                ! this is a 1-4 interaction, label this separately, unless 1-4 are excluded, or we
                ! have explicitly read this in as an exclusion
                molecule_type_data(i_molecule_type)%pair_exclusions(i_atom, local_atom) = 2
             else
                ! fill in exclusion for i_atom
                molecule_type_data(i_molecule_type)%pair_exclusions(i_atom, local_atom) = 1
             end if

             n_bonds_out = n_bonds_in + 1
             update_bond_trajectory(n_bonds_out)= local_atom

             if ( n_bonds_out <= max_search ) then
                call search_bonds_recursive( i_molecule_type, i_atom, local_atom, n_bonds_out , update_bond_trajectory, bond_lookup_table, max_search )
             end if
          end if
       end if
    end do

  end subroutine search_bonds_recursive

  !**********************************
  ! this function returns 0 if integer
  ! index is not present in the first
  ! length elements of the array,
  ! otherwise returns 1
  !*********************************
  integer function check_list( array , length, index )
    integer,dimension(:), intent(in) :: array
    integer, intent(in) :: length, index
    integer :: i
    check_list=0
    do i=1, length
       if ( array(i) == index ) then
          check_list=1
       endif
    end do

  end function check_list

  !******************************************
  !  This subroutine reads a topology file that
  !  specifies the bonds, angles, and dihedrals
  !  of every molecule
  !
  !  the topology file should also give force
  !  field parameters for all of these interactions
  !
  !  the format for the file should be the same as in GROMACS,
  !  with topology given for each molecule type under the
  !  [ moleculetype ]  heading,
  !  with the sections headed with [ bonds ], [ angles ], [ dihedrals ], 
  !  etc.. the force field parameter sections should be headed by
  !  [ bondtypes ] , [ angletypes ] , [ dihedraltypes ] , etc.
  !******************************************
  subroutine read_topology_file( file_h, ifile_topology )
    use global_variables
    integer, intent(in)    :: file_h
    character(*),intent(in):: ifile_topology

    integer :: flag, flag_eof, i_mole
    integer :: flag_bondtypes , flag_angletypes , flag_dihedraltypes 
    integer,dimension(MAX_N_MOLE_TYPE) :: flag_moleculetype
    integer,parameter :: max_param=20


    open(unit=file_h,file=ifile_topology,status="old")

    ! read the [ bondtypes ] , [ angletypes ] , [ dihedraltypes ] sections

    call read_file_find_heading( file_h, '[ bondtypes ]' , flag_bondtypes, flag_eof )
    if ( flag_eof == -1 ) goto 100  ! end of file
    call read_topology_bondtypes( file_h, flag_eof )
    if ( flag_eof == -1 ) goto 100  ! end of file

    call read_file_find_heading( file_h, '[ angletypes ]' , flag_angletypes, flag_eof )
    if ( flag_eof == -1 ) goto 100  ! end of file
    call read_topology_angletypes( file_h, flag_eof )
    if ( flag_eof == -1 ) goto 100  ! end of file

    call read_file_find_heading( file_h, '[ dihedraltypes ]' , flag_dihedraltypes, flag_eof )
    if ( flag_eof == -1 ) goto 100  ! end of file       
    call read_topology_dihedraltypes( file_h, flag_eof )
    if ( flag_eof == -1 ) goto 100  ! end of file


    flag_moleculetype=0
    ! now read the topology file until end, look for each [ moleculetype ] section,
    do 
       call read_file_find_heading( file_h, '[ moleculetype ]' , flag, flag_eof )
       if ( flag_eof == -1 ) Exit  ! end of file
       call read_topology_moleculetype( file_h, flag_eof, flag_moleculetype )
       if ( flag_eof == -1 ) Exit  ! end of file
    enddo

    close(file_h)


100 continue
    ! make sure we found bondtypes section
    if ( flag_bondtypes == 0 ) then
       write(*,*) ""
       write(*,*) "couldn't find '[ bondtypes ]' section in topology file!"
       write(*,*) ""
       stop
    end if
    ! make sure we found angletypes section
    if ( flag_angletypes == 0 ) then
       write(*,*) ""
       write(*,*) "couldn't find '[ angletypes ]' section in topology file!"
       write(*,*) ""
       stop
    end if
    ! make sure we found dihedraltypes section
    if ( flag_dihedraltypes == 0 ) then
       write(*,*) ""
       write(*,*) "couldn't find '[ dihedraltypes ]' section in topology file!"
       write(*,*) ""
       stop
    end if
    ! make sure we found moleculetype section for all molecule types
    do i_mole = 1 , n_molecule_type
       if ( flag_moleculetype(i_mole) == 0 ) then
          write(*,*) ""
          write(*,*) "couldn't find '[ moleculetype ]' section for moleculetype"
          write(*,*) i_mole, "  in topology file!"
          write(*,*) ""
          stop
       end if
    enddo


  end subroutine read_topology_file

  !***************************************
  ! This subroutine reads the bondtypes section
  ! of a topology file
  ! follows gromacs format, but the function
  ! type is ignored (assumed to be harmonic)
  ! also units are consistent with this code rather
  ! than the gromacs topology file, so units are
  ! angstrom, kJ/mol, etc.
  ! end of section should be denoted with blank line
  !***************************************
  subroutine read_topology_bondtypes( file_h , flag_eof )
    use global_variables
    integer, intent(in) :: file_h
    integer, intent(out) :: flag_eof

    integer :: nargs, flag
    integer,parameter :: max_param=20
    character(300) :: input_string
    character(20),dimension(max_param)::args
    character(MAX_ANAME) :: atomtype1, atomtype2
    real*8 :: b0, kb, D, beta
    integer :: index1, index2, bond_type

    ! initialize to zero as a check that we have read in specific bond types
    atype_bond_parameter = 0d0

    flag_eof = 0 
    do
       call read_topology_line( file_h , input_string , flag )
       ! if end of file
       if ( flag == -1 ) then
          flag_eof=-1
          exit
       end if
       ! if end of section
       if ( flag == 1 ) exit

       call parse(input_string," ",args,nargs)
       ! move spaces to end of string for atom types
       call trim_head( args(1) )
       call trim_head( args(2) )
       atomtype1 = args(1)(1:MAX_ANAME)
       atomtype2 = args(2)(1:MAX_ANAME)
       call trim_end(atomtype1)
       call trim_end(atomtype2)
       call atype_name_reverse_lookup( atomtype1, index1 )
       call atype_name_reverse_lookup( atomtype2, index2 )

       read(args(3),*) bond_type
       atype_bond_type( index1, index2 ) = bond_type
       atype_bond_type( index2, index1 ) = bond_type

       Select Case(bond_type)
       Case(1)
          ! harmonic bond type
          read(args(4),*) b0
          read(args(5),*) kb
          ! store bond parameters for this pair of atomtypes
          atype_bond_parameter( index1, index2, 1 ) = b0
          atype_bond_parameter( index1, index2, 2 ) = kb
          ! transpose
          atype_bond_parameter( index2, index1, 1 ) = b0
          atype_bond_parameter( index2, index1, 2 ) = kb
       Case(2)
          ! this is the GROMOS-96 fourth power potential
          read(args(4),*) b0
          read(args(5),*) kb
          atype_bond_parameter( index1, index2, 1 ) = b0
          atype_bond_parameter( index1, index2, 2 ) = kb
          ! transpose
          atype_bond_parameter( index2, index1, 1 ) = b0
          atype_bond_parameter( index2, index1, 2 ) = kb          
       Case(3)
          ! Morse potential
          read(args(4),*) D
          read(args(5),*) beta
          read(args(6),*) b0
          ! store bond parameters for this pair of atomtypes
          atype_bond_parameter( index1, index2, 1 ) = D
          atype_bond_parameter( index1, index2, 2 ) = beta
          atype_bond_parameter( index1, index2, 3 ) = b0
          ! transpose
          atype_bond_parameter( index2, index1, 1 ) = D
          atype_bond_parameter( index2, index1, 2 ) = beta
          atype_bond_parameter( index2, index1, 3 ) = b0
       End Select

    end do

  end subroutine read_topology_bondtypes

  !***************************************
  ! This subroutine reads the angletypes section
  ! of a topology file
  ! follows gromacs format
  ! also units are consistent with this code rather
  ! than the gromacs topology file, so units are
  ! angstrom, kJ/mol, etc.
  ! end of section should be denoted with blank line
  !***************************************
  subroutine read_topology_angletypes( file_h , flag_eof )
    use global_variables
    integer, intent(in) :: file_h
    integer, intent(out) :: flag_eof

    integer :: nargs, flag
    integer,parameter :: max_param=20
    character(300) :: input_string
    character(20),dimension(max_param)::args
    character(MAX_ANAME) :: atomtype1, atomtype2, atomtype3
    real*8 :: th0, cth
    integer :: index1, index2, index3, angle_type

    ! initialize to zero as a check that we have read in specific bond types
    atype_angle_parameter = 0d0

    flag_eof = 0 
    do
       call read_topology_line( file_h , input_string , flag )
       ! if end of file
       if ( flag == -1 ) then
          flag_eof=-1
          exit
       end if
       ! if end of section
       if ( flag == 1 ) exit

       call parse(input_string," ",args,nargs)
       ! move spaces to end of string for atom types
       call trim_head( args(1) )
       call trim_head( args(2) )
       call trim_head( args(3) )
       atomtype1 = args(1)(1:MAX_ANAME)
       atomtype2 = args(2)(1:MAX_ANAME)
       atomtype3 = args(3)(1:MAX_ANAME)
       call trim_end( atomtype1 )
       call trim_end( atomtype2 )
       call trim_end( atomtype3 )
       call atype_name_reverse_lookup( atomtype1, index1 )
       call atype_name_reverse_lookup( atomtype2, index2 )
       call atype_name_reverse_lookup( atomtype3, index3 )

       ! angle type
       read(args(4),*) angle_type
       atype_angle_type( index1, index2, index3 ) = angle_type
       atype_angle_type( index3, index2, index1 ) = angle_type

       ! parameter order is the same for both angle types, either the 
       ! harmonic or cosine angle potential
       read(args(5),*) th0
       read(args(6),*) cth

       ! convert to radians
       th0 = th0 * constants%pi / 180d0

       ! store angle parameters for these atomtypes
       atype_angle_parameter(index1,index2,index3,1) = th0
       atype_angle_parameter(index1,index2,index3,2) = cth
       ! transpose
       atype_angle_parameter(index3,index2,index1,1) = th0
       atype_angle_parameter(index3,index2,index1,2) = cth

    end do

  end subroutine read_topology_angletypes

  !***************************************
  ! This subroutine reads the dihedraltypes section
  ! of a topology file
  ! follows gromacs format
  ! but units are consistent with this code rather
  ! than the gromacs topology file, so units are
  ! angstrom, kJ/mol, etc.
  ! end of section should be denoted with blank line
  !***************************************
  subroutine read_topology_dihedraltypes( file_h, flag_eof )
    use global_variables
    integer, intent(in) :: file_h
    integer, intent(out) :: flag_eof

    integer :: nargs, flag
    integer,parameter :: max_param=20
    character(300) :: input_string
    character(20),dimension(max_param)::args
    character(MAX_ANAME) :: atomtype1, atomtype2, atomtype3, atomtype4
    real*8 :: xi0, kxi, n_mult, c0, c1, c2, c3, c4, c5
    integer :: index1, index2, index3,index4,  dihedral_type

    ! initialize to zero as a check that we have read in specific bond types
    atype_dihedral_parameter = 0d0

    flag_eof = 0 
    do
       call read_topology_line( file_h , input_string , flag )
       ! if end of file
       if ( flag == -1 ) then
          flag_eof=-1
          exit
       end if
       ! if end of section
       if ( flag == 1 ) exit

       call parse(input_string," ",args,nargs)
       ! move spaces to end of string for atom types
       call trim_head( args(1) )
       call trim_head( args(2) )
       call trim_head( args(3) )
       call trim_head( args(4) )
       atomtype1 = args(1)(1:MAX_ANAME)
       atomtype2 = args(2)(1:MAX_ANAME)
       atomtype3 = args(3)(1:MAX_ANAME)
       atomtype4 = args(4)(1:MAX_ANAME)
       call trim_end( atomtype1 )
       call trim_end( atomtype2 )
       call trim_end( atomtype3 )
       call trim_end( atomtype4 )
       call atype_name_reverse_lookup( atomtype1, index1 )
       call atype_name_reverse_lookup( atomtype2, index2 )
       call atype_name_reverse_lookup( atomtype3, index3 )
       call atype_name_reverse_lookup( atomtype4, index4 )

       ! dihedral type
       read(args(5),*) dihedral_type
       atype_dihedral_type( index1, index2, index3, index4 ) = dihedral_type
       atype_dihedral_type( index4, index3, index2, index1 ) = dihedral_type
       


       Select Case(dihedral_type)
       !For RB potential
       Case(3) 
       read(args(6),*) c0
       read(args(7),*) c1
       read(args(8),*) c2
       read(args(9),*) c3
       read(args(10),*) c4
       read(args(11),*) c5
       
       ! store dihedral parameters for these atomtypes
       atype_dihedral_parameter(index1,index2,index3,index4,1) = c0
       atype_dihedral_parameter(index1,index2,index3,index4,2) = c1
       atype_dihedral_parameter(index1,index2,index3,index4,3) = c2
       atype_dihedral_parameter(index1,index2,index3,index4,4) = c3
       atype_dihedral_parameter(index1,index2,index3,index4,5) = c4
       atype_dihedral_parameter(index1,index2,index3,index4,6) = c5
       ! transpose
       atype_dihedral_parameter(index4,index3,index2,index1,1) = c0
       atype_dihedral_parameter(index4,index3,index2,index1,2) = c1
       atype_dihedral_parameter(index4,index3,index2,index1,3) = c2
       atype_dihedral_parameter(index4,index3,index2,index1,4) = c3
       atype_dihedral_parameter(index4,index3,index2,index1,5) = c4
       atype_dihedral_parameter(index4,index3,index2,index1,6) = c5

       Case default
       read(args(6),*) xi0
       read(args(7),*) kxi

       ! convert to radians
       xi0 = xi0 * constants%pi / 180d0

       ! store dihedral parameters for these atomtypes
       atype_dihedral_parameter(index1,index2,index3,index4,1) = xi0
       atype_dihedral_parameter(index1,index2,index3,index4,2) = kxi
       ! transpose
       atype_dihedral_parameter(index4,index3,index2,index1,1) = xi0
       atype_dihedral_parameter(index4,index3,index2,index1,2) = kxi
       End Select

       ! for proper dihedral, we also need multiplicity
       Select Case(dihedral_type)
       Case(1)
          read(args(8),*) n_mult
          atype_dihedral_parameter(index1,index2,index3,index4,3) = n_mult
          atype_dihedral_parameter(index4,index3,index2,index1,3) = n_mult
       End Select
       
    end do


  end subroutine read_topology_dihedraltypes

  !*********************************
  ! this subroutine reads the topology information
  ! for a specific molecule type
  !
  ! we previously determine molecule types from the unique
  ! molecules specified in the .conf file, however this doesn't
  ! always work in an evb simulation, since both conjugate acid and base
  ! may not always be initially present (i.e., we might start with all-protonated
  ! sulfonic acid in .conf file).  Therefore, new molecule types may be present in 
  ! the .top file which weren't present in the .conf file
  !*********************************
  subroutine read_topology_moleculetype( file_h, flag_eof, flag_moleculetype )
    use global_variables
    integer, intent(in) :: file_h
    integer, intent(out) :: flag_eof 
    integer, dimension(:), intent(inout) :: flag_moleculetype

    integer :: nargs, flag, i_mole_type, ind, i_atom, j_atom, k_atom, l_atom, i_type, j_type, k_type , l_type
    integer :: i_bond, i_angle, i_dihedral, n_bond, n_angle, n_dihedral
    integer :: flag_atoms, flag_bonds, flag_angles, flag_dihedrals, flag_exclusions, flag_newmole
    real*8  :: mass_atom
    real*8, parameter :: small = 1D-6
    integer,parameter :: max_param=20
    character(300) :: input_string
    character(20),dimension(max_param)::args
    character(MAX_MNAME) :: i_molecule_name
    character(MAX_ANAME) :: i_atom_name

    flag_atoms = 0
    flag_bonds = 0
    flag_angles = 0
    flag_dihedrals = 0
    flag_exclusions = 0
    flag_eof = 0 

    ! *********** first get molecule name
    call read_topology_line( file_h , input_string , flag )
    ! if end of file
    if ( flag == -1 ) then
       flag_eof=-1
       return
    end if
    call parse(input_string," ",args,nargs)   
    ! first non-comment line should be molecule name
    call trim_head( args(1) )
    i_molecule_name=args(1)(1:MAX_MNAME)
    call trim_end( i_molecule_name )

    ! get molecule_type index
    call mtype_name_reverse_lookup( i_molecule_name, i_mole_type )

    ! if this molecule type wasn't present in the .conf file (i.e for evb simulation), then
    ! define a new molecule type
    if ( i_mole_type < 0 ) then
       n_molecule_type = n_molecule_type + 1
       i_mole_type = n_molecule_type
       molecule_type_data(i_mole_type)%mname = i_molecule_name

       flag_newmole = 1
       write(*,*) "...reading new molecule type ", i_molecule_name
       write(*,*) "from .top file.  This molecule type wasn't present in .gro file,"
       write(*,*) "which makes sense only for an evb simulation..."
       write(*,*) ""
    else
       flag_newmole = 0
    end if

    ! found this moleculetype, set flag
    flag_moleculetype(i_mole_type)=1

    ! ************** now get bond, angles, dihedrals information
    loop1: do

       ! exit loop when we have read all atoms, bonds, angles, dihedrals, exclusions sections
       if ( ( flag_atoms == 1 ) .and. ( flag_bonds == 1 ) .and. ( flag_angles == 1 ) .and. ( flag_dihedrals == 1 ) .and. ( flag_exclusions == 1 ) ) Exit loop1

       call read_topology_line( file_h , input_string , flag )
       ! if end of file
       if ( flag == -1 ) then
          flag_eof=1
          exit
       end if


       !***************** atoms section ********************
       ind=INDEX(input_string,'[ atoms ]')
       if ( ind .ne. 0 ) then
          flag_atoms=1

          ! if this is a new molecule type, first get the number of atoms so that we can allocate data arrays
          if ( flag_newmole == 1 ) then
             molecule_type_data(i_mole_type)%n_atom=0
             do
                call read_topology_line( file_h , input_string , flag )
                ! if end of file, exit outer loop
                if ( flag == -1 ) then
                   flag_eof=1
                   exit loop1
                end if
                ! if blank line, exit inner loop
                if ( flag == 1 ) Exit

                molecule_type_data(i_mole_type)%n_atom = molecule_type_data(i_mole_type)%n_atom + 1
             end do

             ! rewind file
             do i_atom=1,molecule_type_data(i_mole_type)%n_atom+1
                backspace(file_h)
             enddo

             ! now allocate data structures for this molecule type
             allocate( molecule_type_data(i_mole_type)%atom_type_index( molecule_type_data(i_mole_type)%n_atom ) )
             allocate( molecule_type_data(i_mole_type)%pair_exclusions( molecule_type_data(i_mole_type)%n_atom , molecule_type_data(i_mole_type)%n_atom ) )
             ! zero molecule_exclusions here
             molecule_type_data(i_mole_type)%pair_exclusions = 0
             ! for MS-EVB
             allocate(molecule_type_data(i_mole_type)%evb_reactive_basic_atoms( molecule_type_data(i_mole_type)%n_atom ))
             molecule_type_data(i_mole_type)%evb_reactive_basic_atoms=0


          end if


          loop2: do 
             call read_topology_line( file_h , input_string , flag )
             ! if end of file, exit outer loop
             if ( flag == -1 ) then
                flag_eof=1
                exit loop1
             end if
             ! if blank line, exit inner loop
             if ( flag == 1 ) Exit loop2

             call parse(input_string," ",args,nargs) 
             if ( nargs /= 3 ) stop "must have 3 arguments (#,name,mass) under [ atoms ] section in topology file!"
             ! 3 arguments, atom number, name, and mass
             read(args(1),*) i_atom
             read(args(2),*) i_atom_name
             read(args(3),*) mass_atom
             call trim_end( i_atom_name )

             ! if old molecule, make sure this atom type matches .conf file
             if ( flag_newmole == 0 ) then
                i_type = molecule_type_data(i_mole_type)%atom_type_index(i_atom)
                if ( i_atom_name .ne. atype_name(i_type) ) then
                   write(*,*) "for molecule ", i_mole_type, " atom order "
                   write(*,*) "in .top file does not match atom order in "
                   write(*,*) ".conf file"
                   stop
                end if
             else
                ! new molecule, fill in molecule types
                call atype_name_reverse_lookup( i_atom_name, i_type )
                molecule_type_data(i_mole_type)%atom_type_index(i_atom) = i_type
             endif

             if ( atype_mass(i_type) < 0d0 ) then
                ! fill in the mass of this atom type if it hasn't been filled in
                atype_mass(i_type) = mass_atom
             else
                ! make sure mass is same as stored mass for this atom type
                if ( abs(mass_atom - atype_mass(i_type) ) > small ) then
                   write(*,*) "mass for atom type ", atype_name(i_type)
                   write(*,*) "is inconsistent in topology file "
                   stop
                end if
             end if

          end do loop2

       end if


       !***************** bonds section ********************
       ind=INDEX(input_string,'[ bonds ]')
       if ( ind .ne. 0 ) then
          flag_bonds=1
          ! first figure out how many bonds for this molecule type
          n_bond=0
          loop3: do 
             call read_topology_line( file_h , input_string , flag )
             ! if end of file, exit outer loop
             if ( flag == -1 ) then
                flag_eof=1
                exit loop1
             end if
             ! if blank line, exit inner loop
             if ( flag == 1 ) Exit loop3

             n_bond = n_bond + 1
          end do loop3

          ! rewind file
          do i_bond=1,n_bond+1
             backspace(file_h)
          enddo
          ! allocate bond_list
          allocate( molecule_type_data(i_mole_type)%bond_list(n_bond) )

          ! now store bond list
          do i_bond=1, n_bond
             call read_topology_line( file_h , input_string , flag )
             call parse(input_string," ",args,nargs)   
             ! ignore function type here
             read(args(1),*) i_atom
             read(args(2),*) j_atom
             ! add to bondlist
             molecule_type_data(i_mole_type)%bond_list(i_bond)%i_atom = i_atom
             molecule_type_data(i_mole_type)%bond_list(i_bond)%j_atom = j_atom

             ! make sure we have bonded parameters for this pair
             ! check that force constant is non-zero
             i_type = molecule_type_data(i_mole_type)%atom_type_index(i_atom)
             j_type = molecule_type_data(i_mole_type)%atom_type_index(j_atom)

             if ( atype_bond_parameter(i_type,j_type,2) < small ) then
                write(*,*) ""
                write(*,*) "Must have reasonable force constant for bond potential"
                write(*,*) "between atomtypes ", atype_name(i_type), atype_name(j_type)
                write(*,*) ""
                stop
             end if

          end do
       end if

       !***************** angles section ********************
       ind=INDEX(input_string,'[ angles ]')
       if ( ind .ne. 0 ) then
          flag_angles=1

          ! first figure out how many angles for this molecule type
          n_angle=0
          loop4: do
             call read_topology_line( file_h , input_string , flag )
             ! if end of file, exit outer loop
             if ( flag == -1 ) then
                flag_eof=1
                exit loop1
             end if
             ! if blank line, exit inner loop
             if ( flag == 1 ) Exit loop4

             n_angle = n_angle + 1
          end do loop4

          ! rewind file
          do i_angle=1,n_angle+1
             backspace(file_h)
          enddo
          ! allocate angle_list
          allocate( molecule_type_data(i_mole_type)%angle_list(n_angle) )

          ! now store angle list
          do i_angle=1, n_angle
             call read_topology_line( file_h , input_string , flag )
             call parse(input_string," ",args,nargs)   
             ! ignore function type here
             read(args(1),*) i_atom
             read(args(2),*) j_atom
             read(args(3),*) k_atom
             ! add to anglelist
             molecule_type_data(i_mole_type)%angle_list(i_angle)%i_atom=i_atom
             molecule_type_data(i_mole_type)%angle_list(i_angle)%j_atom=j_atom
             molecule_type_data(i_mole_type)%angle_list(i_angle)%k_atom=k_atom

             ! make sure we have angle parameters for this triplet
             ! check that force constant is non-zero
             i_type = molecule_type_data(i_mole_type)%atom_type_index(i_atom)
             j_type = molecule_type_data(i_mole_type)%atom_type_index(j_atom)
             k_type = molecule_type_data(i_mole_type)%atom_type_index(k_atom)

             if ( atype_angle_parameter(i_type,j_type,k_type,2) < small ) then
                write(*,*) ""
                write(*,*) "Must have reasonable force constant for angle potential"
                write(*,*) "between atomtypes ", atype_name(i_type), atype_name(j_type), atype_name(k_type)
                write(*,*) ""
                stop
             end if

          end do
       end if

       !***************** dihedrals section ********************
       ind=INDEX(input_string,'[ dihedrals ]')
       if ( ind .ne. 0 ) then
          flag_dihedrals=1

          ! first figure out how many dihedrals for this molecule type
          n_dihedral=0
          loop5: do
             call read_topology_line( file_h , input_string , flag )
             ! if end of file, exit outer loop
             if ( flag == -1 ) then
                flag_eof=1
                exit loop5
             end if
             ! if blank line, exit inner loop
             if ( flag == 1 ) Exit loop5

             n_dihedral = n_dihedral + 1
          end do loop5

          ! rewind file
          do i_dihedral=1,n_dihedral+1
             backspace(file_h)
          enddo
          ! allocate dihedral_list
          allocate( molecule_type_data(i_mole_type)%dihedral_list(n_dihedral) )

          ! now store dihedral list
          do i_dihedral=1, n_dihedral
             call read_topology_line( file_h , input_string , flag )
             call parse(input_string," ",args,nargs)   
             ! ignore function type here
             read(args(1),*) i_atom
             read(args(2),*) j_atom
             read(args(3),*) k_atom
             read(args(4),*) l_atom
             ! add to dihedral list
             molecule_type_data(i_mole_type)%dihedral_list(i_dihedral)%i_atom=i_atom
             molecule_type_data(i_mole_type)%dihedral_list(i_dihedral)%j_atom=j_atom
             molecule_type_data(i_mole_type)%dihedral_list(i_dihedral)%k_atom=k_atom
             molecule_type_data(i_mole_type)%dihedral_list(i_dihedral)%l_atom=l_atom

             ! make sure we have dihedral parameters for this quartet
             ! check that force constant is non-zero
             i_type = molecule_type_data(i_mole_type)%atom_type_index(i_atom)
             j_type = molecule_type_data(i_mole_type)%atom_type_index(j_atom)
             k_type = molecule_type_data(i_mole_type)%atom_type_index(k_atom)
             l_type = molecule_type_data(i_mole_type)%atom_type_index(l_atom)

             Select Case(atype_dihedral_type( i_type , j_type, k_type, l_type ))
             Case(3)
                continue
             case default
                if ( atype_dihedral_parameter(i_type,j_type,k_type,l_type,2) < small ) then
                   write(*,*) ""
                   write(*,*) "Must have reasonable force constant for dihedral potential"
                   write(*,*) "between atomtypes ", atype_name(i_type), atype_name(j_type), atype_name(k_type), atype_name(l_type) 
                   write(*,*) ""
                   stop
                end if
             end Select

          end do
       end if

       !****************** exclusions section **********************
       ind=INDEX(input_string,'[ exclusions ]')
       if ( ind .ne. 0 ) then
          flag_exclusions=1
          loop6: do 
             call read_topology_line( file_h , input_string , flag )
             ! if end of file, exit outer loop
             if ( flag == -1 ) then
                flag_eof=1
                exit loop1
             end if
             ! if blank line, exit inner loop
             if ( flag == 1 ) Exit loop6
             call parse(input_string," ",args,nargs)   
             read(args(1),*) i_atom
             read(args(2),*) j_atom
             molecule_type_data(i_mole_type)%pair_exclusions(i_atom, j_atom)=1
             molecule_type_data(i_mole_type)%pair_exclusions(j_atom, i_atom)=1
          enddo loop6
       end if


    end do loop1




    ! make sure we found atoms, bonds, angles, dihedrals, and exclusions
    if ( flag_atoms == 0 ) then
       write(*,*) ""
       write(*,*) " couldn't find '[ atoms ]' section in topology file! "
       write(*,*) ""
       stop
    endif
    if ( flag_bonds == 0 ) then
       write(*,*) ""
       write(*,*) " couldn't find '[ bonds ]' section in topology file! "
       write(*,*) ""
       stop
    endif
    if ( flag_angles == 0 ) then
       write(*,*) ""
       write(*,*) " couldn't find '[ angles ]' section in topology file! "
       write(*,*) ""
       stop
    endif
    if ( flag_dihedrals == 0 ) then
       write(*,*) ""
       write(*,*) " couldn't find '[ dihedrals ]' section in topology file! "
       write(*,*) ""
       stop
    endif
    if ( flag_exclusions == 0 ) then
       write(*,*) ""
       write(*,*) " couldn't find '[ exclusions ]' section in topology file! "
       write(*,*) ""
       stop
    endif


  end subroutine read_topology_moleculetype



end module bonded_interactions
