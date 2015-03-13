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
  subroutine intra_molecular_energy_force( E_bond, E_angle, E_dihedral , force_atoms, n_mole , n_atom, xyz )
   use global_variables
    real*8, intent(out) :: E_bond, E_angle, E_dihedral
    real*8, dimension(:,:,:), intent(inout) :: force_atoms
    integer, intent(in) :: n_mole
    integer, dimension(:), intent(in) :: n_atom
    real*8, dimension(:,:,:),intent(in) :: xyz

    E_bond=0d0
    E_angle=0d0
    E_dihedral=0d0

    ! this force array is cumulative, so will be updated by the bond, angle, etc. subroutines
    ! at this point, it should contain contributions from non-bond terms

    ! bond energy, force
    call intra_molecular_bond_energy_force( E_bond, force_atoms, n_mole, n_atom, xyz )

    ! angle energy, force
    call intra_molecular_angle_energy_force( E_angle, force_atoms, n_mole, n_atom, xyz )

    ! dihedral energy force
       !****************************timing**************************************!
       if(debug .eq. 1) then
          call date_and_time(date,time)
          write(*,*) "", "dihedral started at", time
       endif
       !***********************************************************************!
    call intra_molecular_dihedral_energy_force( E_dihedral, force_atoms, n_mole, n_atom, xyz ) 
       !****************************timing**************************************!
       if(debug .eq. 1) then
          call date_and_time(date,time)
          write(*,*) "dihedral finished at", time
       endif
       !***********************************************************************!


  end subroutine intra_molecular_energy_force



  !*************************************
  ! this subroutine returns the intra-molecular
  ! bond energy of the entire system
  ! note there is no pbc image check, because
  ! as is convention in this program, molecules
  ! are not split up over periodic boundaries
  !*************************************
  subroutine intra_molecular_bond_energy_force( E_bond, force_atoms, n_mole, n_atom, xyz ) 
    use global_variables
    real*8, intent(out) :: E_bond
    real*8,dimension(:,:,:),intent(inout) :: force_atoms
    integer, intent(in) :: n_mole
    integer, dimension(:), intent(in) :: n_atom
    real*8, dimension(:,:,:),intent(in) :: xyz

    integer :: i_mole, i_mole_type,i_atom_type, j_atom_type, i_atom, j_atom
    real*8, dimension(3) :: r_ij, force_ij
    real*8 :: r_mag, E_bond_ij


    E_bond=0d0

    do i_mole=1,n_mole
       ! find the molecule type
       i_mole_type = molecule_index(i_mole)
       ! loop over bond list for this molecule
       do i_atom=1,n_atom(i_mole)
          do j_atom=i_atom+1,n_atom(i_mole)
             if ( molecule_bond_list(i_mole_type,i_atom,j_atom) == 1 ) then
                ! bond between these two atoms
                i_atom_type = atom_index(i_mole,i_atom)
                j_atom_type = atom_index(i_mole,j_atom)

                r_ij(:) = xyz(i_mole,i_atom,:) - xyz(i_mole,j_atom,:)
                r_mag = dsqrt( dot_product( r_ij , r_ij ) )

                call pairwise_bond_energy_force( E_bond_ij, force_ij, i_atom_type, j_atom_type, r_ij , r_mag )
                E_bond = E_bond + E_bond_ij 

                force_atoms(i_mole,i_atom,:) = force_atoms(i_mole,i_atom,:) + force_ij(:)
                force_atoms(i_mole,j_atom,:) = force_atoms(i_mole,j_atom,:) - force_ij(:)

             end if
          enddo
       enddo
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
  subroutine intra_molecular_angle_energy_force( E_angle, force_atoms, n_mole, n_atom, xyz ) 
    use global_variables
    real*8, intent(out) :: E_angle
    real*8,dimension(:,:,:),intent(inout) :: force_atoms
    integer, intent(in) :: n_mole
    integer, dimension(:), intent(in) :: n_atom
    real*8, dimension(:,:,:),intent(in) :: xyz

    integer :: i_mole, i_mole_type,i_atom_type, j_atom_type, k_atom_type, i_atom, j_atom, k_atom
    real*8, dimension(3) :: r_ij, r_kj, f_ij , f_kj
    real*8 :: E_angle_ijk

    E_angle=0d0

    do i_mole=1,n_mole
       ! find the molecule type
       i_mole_type = molecule_index(i_mole)
       ! loop over angle list for this molecule
       do i_atom=1,n_atom(i_mole)   ! outer atom
          do k_atom=i_atom+1,n_atom(i_mole)  ! outer atom index > i_atom
             do j_atom=1,n_atom(i_mole)  ! middle atom, loop over all indices
                if ( molecule_angle_list(i_mole_type,i_atom,j_atom,k_atom) == 1 ) then
                   i_atom_type = atom_index(i_mole,i_atom)
                   j_atom_type = atom_index(i_mole,j_atom)
                   k_atom_type = atom_index(i_mole,k_atom)

                   r_ij(:) = xyz(i_mole,i_atom,:) - xyz(i_mole,j_atom,:)
                   r_kj(:) = xyz(i_mole,k_atom,:) - xyz(i_mole,j_atom,:)

                   call trimer_angle_energy_force( E_angle_ijk, f_ij, f_kj, i_atom_type, j_atom_type, k_atom_type , r_ij, r_kj)

                   E_angle = E_angle + E_angle_ijk
                   force_atoms(i_mole,i_atom,:) = force_atoms(i_mole,i_atom,:) + f_ij(:)
                   force_atoms(i_mole,k_atom,:) = force_atoms(i_mole,k_atom,:) + f_kj(:)
                   force_atoms(i_mole,j_atom,:) = force_atoms(i_mole,j_atom,:) - f_ij(:) - f_kj(:)

                end if
             end do
          end do
       end do
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
       theta = acos( cosine )
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
  subroutine intra_molecular_dihedral_energy_force( E_dihedral, force_atoms, n_mole, n_atom, xyz )
    use global_variables
    real*8, intent(out) :: E_dihedral
    real*8,dimension(:,:,:),intent(inout) :: force_atoms
    integer, intent(in) :: n_mole
    integer, dimension(:), intent(in) :: n_atom
    real*8, dimension(:,:,:),intent(in) :: xyz

    integer :: i_mole, i_mole_type,i_atom_type, j_atom_type, k_atom_type, l_atom_type, i_atom, j_atom, k_atom, l_atom
    real*8,dimension(3) :: r_ji, r_kj, r_lk, f_ji, f_kj, f_lk
    real*8 :: E_dihedral_ijkl, E_dihedral_molecule

	! test
	integer :: i

    E_dihedral=0d0

    ! don't use connectivity information here to loop over dihedrals, because impropers might be defined using
    ! non connected atoms


    call OMP_SET_NUM_THREADS(n_threads)
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(xyz,n_mole,molecule_index,molecule_dihedral_flag, n_atom,molecule_dihedral_list,atom_index,atype_dihedral_type,atype_dihedral_parameter, force_atoms ) REDUCTION(+:E_dihedral)
    !$OMP DO SCHEDULE(DYNAMIC,n_threads)
    do i_mole=1,n_mole
       E_dihedral_molecule=0d0
       ! find the molecule type
       i_mole_type = molecule_index(i_mole)
       ! see if this molecule has any dihedrals
       if ( molecule_dihedral_flag( i_mole_type ) == 1 ) then
          ! loop over dihedral list for this molecule
          do i_atom=1,n_atom(i_mole)  ! outer atom                 
             do l_atom = i_atom+1, n_atom(i_mole)  ! outer atom, index greater than i_atom
                ! note we're not double counting angle terms here, even though j_atom and k_atom are
                ! looping over the same atoms, because the order matters for the dihedral potential
                do j_atom=1,n_atom(i_mole)  ! inner atom
                   do k_atom=1,n_atom(i_mole) ! inner atom

                      if ( molecule_dihedral_list(i_mole_type,i_atom,j_atom,k_atom,l_atom) == 1 ) then

                         i_atom_type = atom_index(i_mole,i_atom)
                         j_atom_type = atom_index(i_mole,j_atom)
                         k_atom_type = atom_index(i_mole,k_atom)
                         l_atom_type = atom_index(i_mole,l_atom)

                         r_ji(:) = xyz(i_mole,j_atom,:) - xyz(i_mole,i_atom,:)
                         r_kj(:) = xyz(i_mole,k_atom,:) - xyz(i_mole,j_atom,:)
                         r_lk(:) = xyz(i_mole,l_atom,:) - xyz(i_mole,k_atom,:)


                         call quartet_dihedral_energy_force( E_dihedral_ijkl, f_ji, f_kj, f_lk, i_atom_type, j_atom_type, k_atom_type , l_atom_type, r_ji, r_kj, r_lk )

                         E_dihedral_molecule = E_dihedral_molecule + E_dihedral_ijkl

                         force_atoms(i_mole,i_atom,:) = force_atoms(i_mole,i_atom,:) - f_ji(:)
                         force_atoms(i_mole,j_atom,:) = force_atoms(i_mole,j_atom,:) + f_ji(:) - f_kj(:)
                         force_atoms(i_mole,k_atom,:) = force_atoms(i_mole,k_atom,:) + f_kj(:) - f_lk(:)
                         force_atoms(i_mole,l_atom,:) = force_atoms(i_mole,l_atom,:) + f_lk(:)


!!$                   if ( (abs(f_ji(1)) > 10d4 ) .or. (abs(f_kj(1)) > 10d4 ) .or. (abs(f_lk(1)) > 10d4 ) ) then
!!$                	write(*,*) i_mole, i_atom, j_atom, k_atom, l_atom
!!$			do i=1, n_atom(i_mole)
!!$			write(*,*) i, xyz(i_mole,i,:)
!!$			enddo	
!!$		   write(*,*) "dihedral forces too big"
!!$                      stop
!!$                   end if

                      endif

                   enddo
                enddo
             enddo
          enddo
          E_dihedral = E_dihedral + E_dihedral_molecule
       end if
    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

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
    xi = acos( cosine )

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
          ! make sure shift is zero so that sine is zero
          if ( abs( xi0 ) < small ) then
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
       if ( xi > ( pi / 2d0 ) ) then
          xi = abs( xi - pi )
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

!!$       if ( abs(fac) > 10d4 ) then
!!$          write(*,*) "improper"
!!$          write(*,*) fac
!!$          write(*,*) xi, xi0
!!$          write(*,*) cosine
!!$          stop
!!$       end if

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

    End Select

  end subroutine quartet_dihedral_energy_force



  !*****************************************
  ! This subroutine generates an exclusions list
  ! "molecule_exclusions for every molecule type
  ! for intra-molecular non-bonded interactions
  ! based on the setting of n_excl global variable parameter
  ! atom pairs connected by n_excl bonds or less will be
  ! excluded
  ! exclusions will be marked with a "1"
  !
  ! in addition to exclusions, we may have a special treatment
  ! of 1-4 interactions, and therefore we need to identify 1-4 interactions
  ! as well as exclusion.  1-4 interactions will be marked with
  ! a "2"
  !
  ! normal interactions are marked with a "0"
  !*****************************************
  subroutine generate_intramolecular_exclusions
    use global_variables
    integer :: i_molecule_type, i_atom, j_atom, n_bonds, max_search
    integer, dimension(:), allocatable :: bond_trajectory

    write(*,*) ""
    write(*,*) "automatically generating exclusions for intra-molecular non-bonded"
    write(*,*) "interactions that are within ", n_excl, " bonds away"
    write(*,*) ""

    ! if 1-4 are not excluded ( 3 bonds away ), we want to label these specially
    max_search = max( n_excl , 3 )

    ! this array keeps a list of all atoms along the bonding trajectory, so that we don't loop back over an atom
    allocate( bond_trajectory( max_search+1 ) )
    bond_trajectory=0

    molecule_exclusions=0

    do i_molecule_type=1, n_molecule_type
       do i_atom=1, MAX_N_ATOM
          ! fill in self-exclusion
          molecule_exclusions(i_molecule_type, i_atom, i_atom) = 1

          if ( molecule_type( i_molecule_type, i_atom ) == ( MAX_N_ATOM + 1 ) ) then
             ! end of atoms in molecule
             exit
          end if
          ! this subroutine generates exclusions by searching over bonded neighbors
          ! recursively
            n_bonds=1
            bond_trajectory(1) = i_atom
            call search_bonds_recursive( i_molecule_type, i_atom, i_atom, n_bonds , bond_trajectory, max_search )
       end do
    end do


  end subroutine generate_intramolecular_exclusions



  !******************************************
  ! this subroutine finds bonded neighbors of j_atom, and if those neighbors
  ! are n_bonds away from i_atom, then an exclusion is generated between that
  ! neighbor and i_atom
  !******************************************
  recursive subroutine search_bonds_recursive( i_molecule_type, i_atom, j_atom, n_bonds_in , bond_trajectory, max_search )
    use global_variables
    integer, intent(in) :: i_molecule_type, i_atom, j_atom, n_bonds_in, max_search
    integer, dimension(:), intent(in) :: bond_trajectory
    integer, dimension(size(bond_trajectory)) :: update_bond_trajectory
    integer :: local_atom, n_bonds_out, flag

    update_bond_trajectory = bond_trajectory

    do local_atom=1, MAX_N_ATOM
       ! if there's a bond between this atom and j_atom
       if ( molecule_bond_list(i_molecule_type, j_atom, local_atom) == 1 ) then
          ! make sure we are not doubling back on bonding trajectory
          flag = check_list( bond_trajectory , n_bonds_in , local_atom )
          if ( flag == 0 ) then
             if ( ( n_bonds_in == 3 ) .and. ( n_excl < 3 ) ) then
                ! this is a 1-4 interaction, label this separately, unless 1-4 are excluded
                molecule_exclusions(i_molecule_type, i_atom, local_atom) = 2
             else
                ! fill in exclusion for i_atom
                molecule_exclusions(i_molecule_type, i_atom, local_atom) = 1
             end if

             n_bonds_out = n_bonds_in + 1
             update_bond_trajectory(n_bonds_out)= local_atom

             if ( n_bonds_out <= max_search ) then
                call search_bonds_recursive( i_molecule_type, i_atom, local_atom, n_bonds_out , update_bond_trajectory, max_search )
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
  subroutine read_topology_file( ifile_topology )
    use global_variables
    character(*),intent(in)::ifile_topology

    integer :: file_h = 16
    integer :: nargs, ind, flag, flag_eof, i_mole
    integer :: flag_bondtypes , flag_angletypes , flag_dihedraltypes 
    integer,dimension(MAX_N_MOLE_TYPE) :: flag_moleculetype
    integer,parameter :: max_param=20
    character(300) :: input_string
    character(20),dimension(max_param)::args


    flag_bondtypes=0;flag_angletypes=0;flag_dihedraltypes=0;flag_moleculetype=0;

    open(unit=file_h,file=ifile_topology,status="old")

    ! read the topology file until end, look for each [ moleculetype ] section,
    ! [ bondtypes ] , [ angletypes ] , [ dihedraltypes ] sections

    do
       call read_topology_line( file_h , input_string , flag )
       ! if end of file
       if ( flag == -1 ) exit


       ind=INDEX(input_string,'[ bondtypes ]')
       if ( ind .ne. 0 ) then
          ! bondtypes section
          flag_bondtypes=1
          call read_topology_bondtypes( file_h, flag_eof )
          if ( flag_eof == 1 ) Exit  ! end of file
       end if

       ind=INDEX(input_string,'[ angletypes ]')
       if ( ind .ne. 0 ) then
          ! angletypes section
          flag_angletypes=1
          call read_topology_angletypes( file_h, flag_eof )
          if ( flag_eof == 1 ) Exit  ! end of file
       end if

       !*********** need to code in dihedral types***********
       ind=INDEX(input_string,'[ dihedraltypes ]')
       if ( ind .ne. 0 ) then
          ! angletypes section
          flag_dihedraltypes=1
          call read_topology_dihedraltypes( file_h, flag_eof )
          if ( flag_eof == 1 ) Exit  ! end of file
       end if

       ind=INDEX(input_string,'[ moleculetype ]')
       if ( ind .ne. 0 ) then
          !***************** get moleculetype information
          call read_topology_moleculetype( file_h, flag_eof, flag_moleculetype )
          if ( flag_eof == 1 ) Exit  ! end of file
       end if

    enddo

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

    flag_eof = 0 
    do
       call read_topology_line( file_h , input_string , flag )
       ! if end of file
       if ( flag == -1 ) then
          flag_eof=1
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

    flag_eof = 0 
    do
       call read_topology_line( file_h , input_string , flag )
       ! if end of file
       if ( flag == -1 ) then
          flag_eof=1
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
       th0 = th0 * pi / 180d0

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
    real*8 :: xi0, kxi, n_mult
    integer :: index1, index2, index3,index4,  dihedral_type

    flag_eof = 0 
    do
       call read_topology_line( file_h , input_string , flag )
       ! if end of file
       if ( flag == -1 ) then
          flag_eof=1
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

       read(args(6),*) xi0
       read(args(7),*) kxi

       ! convert to radians
       xi0 = xi0 * pi / 180d0

       ! store dihedral parameters for these atomtypes
       atype_dihedral_parameter(index1,index2,index3,index4,1) = xi0
       atype_dihedral_parameter(index1,index2,index3,index4,2) = kxi
       ! transpose
       atype_dihedral_parameter(index4,index3,index2,index1,1) = xi0
       atype_dihedral_parameter(index4,index3,index2,index1,2) = kxi

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
  !*********************************
  subroutine read_topology_moleculetype( file_h, flag_eof, flag_moleculetype )
    use global_variables
    integer, intent(in) :: file_h
    integer, intent(out) :: flag_eof 
    integer, dimension(:), intent(inout) :: flag_moleculetype

    integer :: nargs, flag, i_mtype, i_mole_type, ind, i_atom, j_atom, k_atom, l_atom
    integer :: flag_bonds, flag_angles, flag_dihedrals
    integer,parameter :: max_param=20
    character(300) :: input_string
    character(20),dimension(max_param)::args
    character(MAX_MNAME) :: i_molecule_name

    flag_bonds = 0
    flag_angles = 0
    flag_dihedrals = 0
    flag_eof = 0 

    ! *********** first get molecule name
    call read_topology_line( file_h , input_string , flag )
    ! if end of file
    if ( flag == -1 ) then
       flag_eof=1
       return
    end if
    call parse(input_string," ",args,nargs)   
    ! first non-comment line should be molecule name
    call trim_head( args(1) )
    i_molecule_name=args(1)(1:MAX_MNAME)
    call trim_end( i_molecule_name )

    ! get molecule_type index
    call mtype_name_reverse_lookup( i_molecule_name, i_mole_type )
    ! initialize bondlist and anglelist and dihedrallist for this molecule type
    molecule_bond_list( i_mole_type,:,:) = 0
    molecule_angle_list( i_mole_type,:,:,:) = 0
    molecule_dihedral_list( i_mole_type,:,:,:,:) = 0
    molecule_dihedral_flag( i_mole_type ) = 0
    ! found this moleculetype, set flag
    flag_moleculetype(i_mole_type)=1

    ! ************** now get bond, angles, dihedrals information
    loop1: do

       ! exit loop when we have read all bonds, angles, dihedrals sections
       if ( ( flag_bonds == 1 ) .and. ( flag_angles == 1 ) .and. ( flag_dihedrals == 1 ) ) Exit loop1

       call read_topology_line( file_h , input_string , flag )
       ! if end of file
       if ( flag == -1 ) then
          flag_eof=1
          exit
       end if

       !***************** bonds section ********************
       ind=INDEX(input_string,'[ bonds ]')
       if ( ind .ne. 0 ) then
          flag_bonds=1
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
             ! ignore function type here
             read(args(1),*) i_atom
             read(args(2),*) j_atom
             ! add to bondlist
             molecule_bond_list(i_mole_type,i_atom,j_atom)=1
             molecule_bond_list(i_mole_type,j_atom,i_atom)=1

          end do loop2
       end if

       !***************** angles section ********************
       ind=INDEX(input_string,'[ angles ]')
       if ( ind .ne. 0 ) then
          flag_angles=1
          loop3: do 
             call read_topology_line( file_h , input_string , flag )
             ! if end of file, exit outer loop
             if ( flag == -1 ) then
                flag_eof=1
                exit loop1
             end if
             ! if blank line, exit inner loop
             if ( flag == 1 ) Exit loop3

             call parse(input_string," ",args,nargs)   
             ! ignore function type here
             read(args(1),*) i_atom
             read(args(2),*) j_atom
             read(args(3),*) k_atom
             ! add to anglelist
             molecule_angle_list(i_mole_type,i_atom,j_atom,k_atom)=1
             molecule_angle_list(i_mole_type,k_atom,j_atom,i_atom)=1

          end do loop3
       end if

       !***************** dihedrals section ********************
       ind=INDEX(input_string,'[ dihedrals ]')
       if ( ind .ne. 0 ) then
          flag_dihedrals=1
          loop4: do 
             call read_topology_line( file_h , input_string , flag )
             ! if end of file, exit outer loop
             if ( flag == -1 ) then
                flag_eof=1
                exit loop1
             end if
             ! if blank line, exit inner loop
             if ( flag == 1 ) Exit loop4

             call parse(input_string," ",args,nargs)   
             ! ignore function type here
             read(args(1),*) i_atom
             read(args(2),*) j_atom
             read(args(3),*) k_atom
             read(args(4),*) l_atom
             ! add to dihedral list
             molecule_dihedral_list(i_mole_type,i_atom,j_atom,k_atom,l_atom)=1
             molecule_dihedral_list(i_mole_type,l_atom,k_atom,j_atom,i_atom)=1
             ! if any dihedrals, set to 1
             molecule_dihedral_flag( i_mole_type ) = 1

          end do loop4
       end if


    end do loop1

    ! make sure we found bonds and angles and dihedrals
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

  end subroutine read_topology_moleculetype



  !******************************
  ! this subroutine returns the next non-comment
  ! line of a topology file.
  !
  ! lines that begin with ";" are taken as comments
  ! (gromacs convention)
  !
  ! settings of output flag:
  ! 0 : normal line
  ! 1 : blank line
  ! -1 : end of file
  !******************************
  subroutine read_topology_line( file_h , line , flag )
    integer, intent(in) :: file_h
    character(*), intent(out) :: line
    integer, intent(out) :: flag

    ! use gromacs comment convention
    character(5) :: comment=";"
    integer :: nargs, inputstatus, ind
    integer,parameter :: max_param=20
    character(20),dimension(max_param)::args

    flag=0
    do
       Read(file_h,'(A)',Iostat=inputstatus) line
       If(inputstatus < 0) then
          ! end of file
          flag=-1
          Exit
       end if

       ! check end of section
       call parse(line," ",args,nargs)
       if ( nargs == 0 ) then
          flag=1
          Exit
       end if
       ! check comment
       ind=INDEX(args(1),comment)
       ! if not a comment line, then return this string
       if ( ind == 0 ) exit
    end do

  end subroutine read_topology_line



end module bonded_interactions
