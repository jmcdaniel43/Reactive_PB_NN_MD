module pairwise_interaction
  use global_variables
  use routines
  implicit none

  !***********************************************************************
  ! This module contains subroutines to calculate non-electrostatic pairwise
  ! interaction energies such as Lennard-Jones or Buckingham interactions
  !
  !**********************************************************************

contains


  !***********************************************************************
  ! This subroutine computes lennard jones (buckingham) forces and energy
  ! It is implemented for generic solute, as atom index array stores the type
  ! of atom for each molecule, which is then used to recover LJ parameters in atype_lj_parameter array
  ! also, storing framework atoms in tot_n_mole - n_mole indices allows solute framework interactions to
  ! be calculated
  !  As implemented, some forces are added to framework atom indices of force array for free due to Newton's 3rd law,
  !  this is not important and does nothing.
  !
  ! we have removed some of the subroutine calls inside the loops and replaced these calls
  ! with explicit code to enhance performance
  !
  ! Either verlet list or naive looping can be used
  !***********************************************************************
  subroutine lennard_jones( lj_energy, lj_force, tot_n_mole, n_mole, n_atom, xyz, box, r_com)
    real*8,intent(out) :: lj_energy
    real*8,dimension(:,:,:),intent(out)::lj_force
    integer,intent(in)::n_mole,tot_n_mole
    integer,dimension(:),intent(in)::n_atom
    real*8,dimension(:,:,:),intent(in)::xyz
    real*8,dimension(:,:),intent(in)::box
    real*8,dimension(:,:),intent(in)::r_com
    integer::i_mole,i_atom, a_index
    real*8 :: E_intra_lj, lj_cutoff2
    real*8,dimension(:,:), allocatable :: atom_list_xyz, atom_list_force
    integer,dimension(:), allocatable  :: atom_list_atom_index


    lj_force=0d0

    ! inter-molecular contributions
    Select Case(use_verlet_list)
    Case("no")
       call  lennard_jones_no_verlet( lj_energy, lj_force, tot_n_mole, n_mole, n_atom, xyz, box, r_com)
    Case("yes")
       ! here we need to map molecular data structures to atomic data structures that are
       ! consistent with the verlet list
       ! as these arrays are in inner loops, put innermost index first
       allocate( atom_list_xyz(3,verlet_lj_atoms), atom_list_force(3,verlet_lj_atoms) , atom_list_atom_index(verlet_lj_atoms) )
       call map_molecule_atom_data_structures_3d_m_to_a( tot_n_mole, n_atom , verlet_molecule_map_lj, atom_list_xyz , xyz )
       call map_molecule_atom_data_structures_2d_m_to_a_int( tot_n_mole, n_atom , verlet_molecule_map_lj, atom_list_atom_index , atom_index )
       ! ****** lj energy and force subroutine
       call  lennard_jones_use_verlet( lj_energy, atom_list_force, atom_list_xyz , atom_list_atom_index, box )
       ! map force back to molecular data structure
       call map_molecule_atom_data_structures_3d_a_to_m( tot_n_mole, n_atom, verlet_molecule_map_lj , atom_list_force , lj_force )
       deallocate( atom_list_xyz , atom_list_force, atom_list_atom_index )
    End Select

    lj_cutoff2 = lj_cutoff ** 2
    ! now intra-molecular contributions
    do i_mole = 1, n_mole
       call intra_lennard_jones_energy_force( E_intra_lj, lj_force, i_mole, n_atom, molecule_index, atom_index,  xyz , lj_cutoff2 )
       lj_energy = lj_energy + E_intra_lj
    enddo


  end subroutine lennard_jones



  !***********************************************************************
  ! In this subroutine, the data structures have been changed from molecule based
  ! data structures, to atom_list based data structures.  This is because we only want
  ! to consider inter-molecular interactions between atoms that have non-zero
  ! lennard-jones parameters.  The indexing of atoms in the atom_list data structures
  ! should be consistent with the indexing of atoms in the verlet list
  !***********************************************************************
  subroutine lennard_jones_use_verlet( lj_energy,atom_list_force, atom_list_xyz , atom_list_atom_index, box )
    use omp_lib
    real*8,intent(out) :: lj_energy
    real*8,dimension(:,:), intent(out) :: atom_list_force
    real*8,dimension(:,:), intent(in)  :: atom_list_xyz
    integer,dimension(:), intent(in)   :: atom_list_atom_index
    real*8,dimension(:,:), intent(in)  :: box

    integer:: i_atom, j_atom, verlet_start, verlet_finish,i_index,atom_id1,atom_id2,i
    integer :: split_do
    real*8 :: lj_cutoff2, norm_dr2, norm_dr6, norm_dr12, term6, term12
    real*8,dimension(3)::f_ij,rij,shift,shift_direct, dr_com,dr_direct

    atom_list_force=0d0
    lj_energy=0d0

    lj_cutoff2 = lj_cutoff ** 2

    ! decide how to split the parallel section
    if (n_threads .eq. 1 ) then
       split_do = 1
    else
       split_do = verlet_lj_atoms/n_threads+1
    endif


    !**************************************** use Verlet list ****************************************************************
    call OMP_SET_NUM_THREADS(n_threads)
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(atom_list_xyz, atom_list_atom_index, verlet_lj_atoms, box,atype_lj_parameter,lj_cutoff2,split_do,xyz_to_box_transform,verlet_point_lj,verlet_neighbor_list_lj) REDUCTION(+:atom_list_force, lj_energy)
    !$OMP DO SCHEDULE(DYNAMIC,split_do)
    ! note that in Verlet list, neighbors ( on j_mole) are only stored if j_mole > i_mole , to avoid double counting
    do i_atom=1, verlet_lj_atoms

       atom_id1 = atom_list_atom_index(i_atom)
       ! get j_atom from Verlet list
       verlet_start = verlet_point_lj(i_atom)
       i_index = i_atom + 1
       verlet_finish = verlet_point_lj(i_index) - 1

       ! make sure there's at least one neighbor
       if ( verlet_finish .ge. verlet_start ) then
          do i_index= verlet_start, verlet_finish    
             j_atom = verlet_neighbor_list_lj(i_index)
             atom_id2 = atom_list_atom_index(j_atom)

             ! coordinates should be stored in this order, optimizing fortran
             ! column major memory retrieval
             rij = atom_list_xyz(:,i_atom) - atom_list_xyz(:,j_atom)

             ! shift for general box
!!$             dr_direct(:) = matmul( xyz_to_box_transform, rij )
!!$             do i=1,3
!!$                shift_direct(i) = dble(floor( dr_direct(i) + 0.5d0 ))
!!$             enddo
!!$             shift = matmul( shift_direct , box )
             ! shift for orthorhombic box
             do i=1,3
                shift(i) = box(i,i) * floor( rij(i) / box(i,i) + 0.5d0 )
             end do

             rij = rij - shift
             norm_dr2 = dot_product( rij, rij )              

             if ( norm_dr2 < lj_cutoff2 ) then
                norm_dr6 = norm_dr2 ** 3
                norm_dr12 = norm_dr6 ** 2
                ! at this stage, all lj parameters should be expressed as C12 and C6, even though they were read in as epsilon and sigma
                term12 = atype_lj_parameter(atom_id1,atom_id2,1) / norm_dr12
                term6 = atype_lj_parameter(atom_id1,atom_id2,2) / norm_dr6
                lj_energy = lj_energy + term12 - term6 
                f_ij = rij / norm_dr2 * ( 12d0 * term12  - 6d0 * term6 )

                ! forces should be stored in this order, optimizing fortran
                ! column major memory retrieval
                atom_list_force(:,i_atom) = atom_list_force(:,i_atom) + f_ij(:)
                atom_list_force(:,j_atom) = atom_list_force(:,j_atom) - f_ij(:)
             end if
          end do
       endif
    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL


  end  subroutine lennard_jones_use_verlet





  !***********************************************************************
  ! self-explanatory
  !***********************************************************************
  subroutine lennard_jones_no_verlet( lj_energy, lj_force, tot_n_mole, n_mole, n_atom, xyz, box, r_com)
    use omp_lib
    real*8,intent(out) :: lj_energy
    real*8,dimension(:,:,:),intent(out)::lj_force
    integer,intent(in)::n_mole,tot_n_mole
    integer,dimension(:),intent(in)::n_atom
    real*8,dimension(:,:,:),intent(in)::xyz
    real*8,dimension(:,:),intent(in)::box
    real*8,dimension(:,:),intent(in)::r_com
    integer::i_mole,j_mole,i_atom,j_atom,atom_id1,atom_id2,i
    integer :: split_do
    real*8 :: lj_cutoff2, norm_dr2, norm_dr6, norm_dr12, term6, term12
    real*8,dimension(3)::f_ij,rij,shift,shift_direct, dr_com,dr_direct

    lj_force=0d0
    lj_energy=0d0

    lj_cutoff2 = lj_cutoff ** 2

    ! decide how to split the parallel section
    if (n_threads .eq. 1 ) then
       split_do = 1
    else
       split_do = n_mole/n_threads+1
    endif


    !**************************************** no Verlet list, naive looping ***********************************************
    call OMP_SET_NUM_THREADS(n_threads)

    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(tot_n_mole,n_mole,r_com,box,n_atom,xyz,atom_index,atype_lj_parameter,lj_cutoff2,split_do,xyz_to_box_transform) REDUCTION(+:lj_force, lj_energy)
    !$OMP DO SCHEDULE(DYNAMIC,split_do)
    do i_mole=1,n_mole
       do j_mole=i_mole+1,tot_n_mole
          do i_atom=1,n_atom(i_mole)
             do j_atom=1,n_atom(j_mole)
                atom_id1 = atom_index(i_mole,i_atom)
                atom_id2 = atom_index(j_mole,j_atom)
                rij = xyz(i_mole,i_atom,:) - xyz(j_mole,j_atom,:)

                dr_direct(:) = matmul( xyz_to_box_transform, rij )
                do i=1,3
                   shift_direct(i) = dble(floor( dr_direct(i) + 0.5d0 ))
                enddo
                shift = matmul( shift_direct , box )
                rij = rij - shift
                norm_dr2 = dot_product( rij, rij )  
                if ( norm_dr2 < lj_cutoff2 ) then          
                   norm_dr6 = norm_dr2 ** 3
                   norm_dr12 = norm_dr6 ** 2
                   ! at this stage, all lj parameters should be expressed as C12 and C6, even though they were read in as epsilon and sigma
                   term12 = atype_lj_parameter(atom_id1,atom_id2,1) / norm_dr12
                   term6 = atype_lj_parameter(atom_id1,atom_id2,2) / norm_dr6
                   lj_energy = lj_energy + term12 - term6
                   f_ij = rij / norm_dr2 * ( 12d0 * term12  - 6d0 * term6 )
                   lj_force(i_mole,i_atom,:) = lj_force(i_mole,i_atom,:) + f_ij(:)
                   lj_force(j_mole,j_atom,:) = lj_force(j_mole,j_atom,:) - f_ij(:)

                end if
             end do
          end do
       enddo
    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL


  end  subroutine lennard_jones_no_verlet




  !*******************************************
  ! this subroutine computes the intra-molecular lennard jones force and energy contribution
  !*******************************************
  subroutine intra_lennard_jones_energy_force( E_intra_lj , lj_force, i_mole, n_atom, molecule_index_local, atom_index_local, xyz, lj_cutoff2 )
    use global_variables
    real*8 , intent(out) :: E_intra_lj
    real*8 , dimension(:,:,:),intent(inout) :: lj_force
    integer, intent(in) :: i_mole
    integer, dimension(:), intent(in) :: n_atom
    integer, dimension(:), intent(in) :: molecule_index_local
    integer, dimension(:,:), intent(in) :: atom_index_local
    real*8, dimension(:,:,:), intent(in) :: xyz
    real*8, intent(in) :: lj_cutoff2

    integer :: i_atom, j_atom, i_mole_type, atom_id1, atom_id2
    real*8, dimension(3) :: rij, f_ij
    real*8 :: norm_dr2, norm_dr6, norm_dr12, term12, term6, C6, C12

    E_intra_lj = 0d0

    i_mole_type = molecule_index_local(i_mole)

    do i_atom=1, n_atom(i_mole)
       do j_atom = i_atom + 1, n_atom(i_mole)
          ! check for exclusions 
          if ( molecule_exclusions( i_mole_type, i_atom, j_atom ) /= 1 ) then
             ! not excluded
             ! as always, molecules are not broken up over pbc, so no need to consider shift

             rij = xyz(i_mole,i_atom,:) - xyz(i_mole,j_atom,:)
             norm_dr2 = dot_product( rij, rij )              
             if ( norm_dr2 < lj_cutoff2 ) then

                norm_dr6 = norm_dr2 ** 3
                norm_dr12 = norm_dr6 ** 2

                atom_id1 = atom_index_local(i_mole,i_atom)
                atom_id2 = atom_index_local(i_mole,j_atom)

                ! at this stage, all lj parameters should be expressed as C12 and C6, even though they were read in as epsilon and sigma
                ! if this is a 1-4 interaction, take parameters from 1-4 interaction parameter array
                if ( molecule_exclusions( i_mole_type, i_atom , j_atom ) == 2 ) then
                   C12 = atype_lj_parameter_14(atom_id1,atom_id2,1)
                   C6  = atype_lj_parameter_14(atom_id1,atom_id2,2)
                else
                   C12 = atype_lj_parameter(atom_id1,atom_id2,1)
                   C6  = atype_lj_parameter(atom_id1,atom_id2,2)
                end if

                term12 = C12 / norm_dr12
                term6 =  C6 / norm_dr6

                E_intra_lj = E_intra_lj + term12 - term6
                f_ij = rij / norm_dr2 * ( 12d0 * term12  - 6d0 * term6 )
                lj_force(i_mole,i_atom,:) = lj_force(i_mole,i_atom,:) + f_ij(:)
                lj_force(i_mole,j_atom,:) = lj_force(i_mole,j_atom,:) - f_ij(:)


             end if
          end if
       enddo
    enddo


  end subroutine intra_lennard_jones_energy_force



end module pairwise_interaction
