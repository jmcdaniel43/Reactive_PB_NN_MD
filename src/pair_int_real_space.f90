module pair_int_real_space
  use routines
  implicit none

  !*******************************************************************
  ! REAL-SPACE INTERACTION SUBROUTINES
  !
  ! This module contains subroutines to compute energy and forces
  ! for pairwise real-space interactions using a Verlet neighbor list
  ! this covers both PME real space and VDWs real-space interactions
  !*******************************************************************

    !** this is data structure which contains all the necessary
    !** data to compute pairwise interactions for neighbors and should be
    !** allocated locally
    type pairwise_neighbor_data_type
      integer, dimension(:) , allocatable :: atom_index
      real*8 , dimension(:,:) , allocatable :: xyz
      real*8 , dimension(:,:) , allocatable :: dr
      real*8 , dimension(:,:) , allocatable :: f_ij
      real*8 , dimension(:)   , allocatable :: dr2
      ! these are force field parameters needed for interactions
      real*8 , dimension(:)   , allocatable :: qi_qj   ! charge * charge
      real*8 , dimension(:,:) , allocatable :: atype_vdw_parameter
    end type pairwise_neighbor_data_type

contains

    !****** these subroutines control allocatation/deallocation of pairwise_neighbor data structures
    subroutine allocate_pairwise_neighbor_data( pairwise_neighbor_data , n_neighbors , size_vdw_parameter )
       type(pairwise_neighbor_data_type), intent(inout) :: pairwise_neighbor_data
       integer , intent(in)   :: n_neighbors, size_vdw_parameter

       allocate( pairwise_neighbor_data%atom_index(n_neighbors) )
       allocate( pairwise_neighbor_data%xyz(3,n_neighbors) )
       allocate( pairwise_neighbor_data%dr(3,n_neighbors) )
       allocate( pairwise_neighbor_data%f_ij(3,n_neighbors) )
       allocate( pairwise_neighbor_data%dr2(n_neighbors) )
       allocate( pairwise_neighbor_data%qi_qj(n_neighbors) )
       allocate( pairwise_neighbor_data%atype_vdw_parameter(size_vdw_parameter , n_neighbors) )

       ! zero forces
       pairwise_neighbor_data%f_ij=0d0

    end subroutine allocate_pairwise_neighbor_data


    subroutine deallocate_pairwise_neighbor_data( pairwise_neighbor_data )
       type(pairwise_neighbor_data_type), intent(inout) ::  pairwise_neighbor_data
       deallocate( pairwise_neighbor_data%atom_index , pairwise_neighbor_data%xyz , pairwise_neighbor_data%dr , pairwise_neighbor_data%f_ij , pairwise_neighbor_data%dr2 , pairwise_neighbor_data%qi_qj ,  pairwise_neighbor_data%atype_vdw_parameter )
    end subroutine deallocate_pairwise_neighbor_data




  !***********************************************************************
  ! This calculates real-space, pairwise contributions to force and energy,
  ! including non-bonded VDWs interactions and real-space portion of 
  ! PME electrostatic energy and forces
  !***********************************************************************
  subroutine  real_space_energy_force( system_data, molecule_data, atom_data, verlet_list_data, PME_data )
    use global_variables
    use MKL_DFTI
    use routines
    use omp_lib
    implicit none

    Type(system_data_type),intent(inout)             :: system_data
    Type(molecule_data_type),dimension(:),intent(in) :: molecule_data
    Type(atom_data_type),intent(inout)               :: atom_data
    Type(verlet_list_data_type),intent(in)           :: verlet_list_data
    Type(PME_data_type), intent(in)                  :: PME_data

    !******** this is a local data structure with pointers that will be set
    ! to subarrays of atom_data arrays for the specific atoms in the molecule
    type(single_molecule_data_type) :: single_molecule_data

    integer :: i_mole
    

    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "real space energy and force calculation started at", time
    endif
    !***********************************************************************!


    !******************************** inter-molecular real-space interactions *******************************

    !*******   pme real space energy and force subroutine
    call pairwise_real_space_verlet( system_data , atom_data , verlet_list_data , PME_data )


    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "intra-molecular non-bonded interactions started at", time
    endif
    !***********************************************************************!


    !******************************** intra-molecular real-space interactions *******************************
    do i_mole=1,system_data%n_mole
      if( molecule_data(i_mole)%n_atom > 1) then

        ! set pointers for this data structure to target molecule
        call return_molecule_block( single_molecule_data , molecule_data(i_mole)%n_atom, molecule_data(i_mole)%atom_index , atom_xyz=atom_data%xyz , atom_force=atom_data%force, atom_charge=atom_data%charge, atom_type_index=atom_data%atom_type_index )

        ! here we explicitly pass energy and force data structures to be modified.  We allow this flexibility for calls within the MS-EVB subroutines in
        ! which case we want to use local data structures
        call intra_molecular_pairwise_energy_force( single_molecule_data%force, system_data%E_elec, system_data%E_vdw,  single_molecule_data ,  molecule_data(i_mole)%molecule_type_index , molecule_data(i_mole)%n_atom, PME_data )

       call dissociate_single_molecule_data(single_molecule_data)
       endif
    enddo



    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "real space energy and force calculation finished at", time
    endif
    !***********************************************************************!


  end subroutine  real_space_energy_force




  !***********************************************
  ! this subroutine calculates the real space part of 
  ! the pme energy using an atom-atom based verlet list
  !
  ! In this subroutine, the data structures have been changed from molecule based
  ! data structures, to atom_list based data structures.  This is because we only want
  ! to consider inter-molecular interactions between atoms that have non-zero
  ! charges.  The indexing of atoms in the atom_list data structures
  ! should be consistent with the indexing of atoms in the verlet list
  !***********************************************************************
  subroutine pairwise_real_space_verlet( system_data , atom_data , verlet_list_data , PME_data  )
    use global_variables
    use routines
    use omp_lib
    implicit none
    type(system_data_type), intent(inout) :: system_data
    type(atom_data_type), intent(inout)   :: atom_data
    type(verlet_list_data_type) , intent(in) :: verlet_list_data
    type(PME_data_type) , intent(in)         :: PME_data

    !****** note the storage dimensions of these arrays, consistent with atom_data%force
    real*8, dimension(3,size(atom_data%force(1,:))) :: local_force
    real*8, dimension(3,size(atom_data%force(1,:)),n_threads) :: temp_force

    Type(pairwise_neighbor_data_type) :: pairwise_neighbor_data_verlet , pairwise_neighbor_data_cutoff

    integer , dimension(:) , allocatable    :: cutoff_mask
    real*8, dimension(3,3)   :: box, xyz_to_box_transform
    real*8                   :: real_space_cutoff2, alpha_sqrt, erf_factor, E_elec, E_vdw, E_elec_local , E_vdw_local
    integer                  :: split_do, total_atoms, n_neighbors, n_cutoff, size_vdw_parameter
    integer :: i, i_atom, j_atom, thread_id, i_thread, verlet_start, verlet_finish, i_index
    real*8, dimension(3)     :: xyz_i_atom(3)


    E_elec=0d0
    E_vdw =0d0
    local_force = 0.D0
    temp_force= 0.D0
    ! real_space_cutoff is a global variable
    real_space_cutoff2 = real_space_cutoff ** 2

    ! store small data structures as local variable for convenience
    total_atoms          = system_data%total_atoms
    box                  = system_data%box
    xyz_to_box_transform = system_data%xyz_to_box_transform 
    alpha_sqrt           = PME_data%alpha_sqrt
    erf_factor           =  2.D0*alpha_sqrt/constants%pi_sqrt
    ! this is (maximum) number of parameters that we need for VDWs interaction
    size_vdw_parameter=size(atype_vdw_parameter(1,1,:))

    ! decide how to split the parallel section
    if (n_threads .eq. 1 ) then
       split_do = 1
    else
       split_do = total_atoms/n_threads+1
    endif


    !**************************************** use Verlet list ****************************************************************
    call OMP_SET_NUM_THREADS(n_threads)
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(atom_data, total_atoms, box, temp_force, real_space_cutoff2, split_do , xyz_to_box_transform, atype_vdw_parameter, size_vdw_parameter, erf_factor, alpha_sqrt,verlet_list_data, PME_data, constants) REDUCTION(+:E_elec,E_vdw)
    !$OMP CRITICAL
    local_force = 0.D0
    !$OMP END CRITICAL
    !$OMP DO SCHEDULE(DYNAMIC,split_do)
    ! note that in Verlet list, neighbors are only stored if j_atom > i_atom , for efficiency and to avoid double counting
    do i_atom=1 , total_atoms

       xyz_i_atom(:) = atom_data%xyz(:,i_atom)

       verlet_start = verlet_list_data%verlet_point(i_atom)
       i_index = i_atom + 1
       verlet_finish = verlet_list_data%verlet_point(i_index) - 1
       n_neighbors = verlet_finish - verlet_start + 1

       ! make sure there's at least one neighbor
       if ( n_neighbors > 0 ) then

          ! allocate data structure to store data for these neighbors
          call allocate_pairwise_neighbor_data( pairwise_neighbor_data_verlet , n_neighbors, size_vdw_parameter )
          
          ! group neighbor coordinates and parameters sequentially in array for vectorization
          call gather_neighbor_data( i_atom, verlet_start , verlet_finish , verlet_list_data%neighbor_list, pairwise_neighbor_data_verlet, atom_data%xyz, atom_data%charge, atom_data%atom_type_index, atype_vdw_parameter )

          ! this loop should now perfectly vectorize...
          ! $omp simd
          do j_atom = 1 , n_neighbors
             pairwise_neighbor_data_verlet%dr(:,j_atom) = xyz_i_atom(:) - pairwise_neighbor_data_verlet%xyz(:,j_atom)
          enddo

          ! PBC minimum image, this loop should now perfectly vectorize...
          do j_atom = 1 , n_neighbors

             ! shift for general box
             !             dr_direct(:) = matmul( xyz_to_box_transform, rij )
             !             do i=1,3
             !                shift_direct(i) = dble(floor( dr_direct(i) + 0.5d0 ))
             !             enddo
             !             shift = matmul( shift_direct , box )

             ! shift for orthorhombic box
             do i=1,3
                pairwise_neighbor_data_verlet%dr(i,j_atom) = pairwise_neighbor_data_verlet%dr(i,j_atom) -  box(i,i) * floor( pairwise_neighbor_data_verlet%dr(i,j_atom) / box(i,i) + 0.5d0 )
             end do
             pairwise_neighbor_data_verlet%dr2(j_atom) = dot_product( pairwise_neighbor_data_verlet%dr(:,j_atom), pairwise_neighbor_data_verlet%dr(:,j_atom) )
          enddo
          ! $omp end simd

          ! now check cutoff, and re-organize data structures
          allocate( cutoff_mask(n_neighbors) )
          cutoff_mask=0
          ! $omp simd
          do j_atom = 1 , n_neighbors
             if ( pairwise_neighbor_data_verlet%dr2(j_atom) < real_space_cutoff2 ) then
                cutoff_mask(j_atom) = 1
             endif
          enddo
          ! $omp end simd
          n_cutoff = sum(cutoff_mask)

          ! now allocate datastructure for atoms within cutoff distance
          ! allocate data structure to store data for these neighbors
          call allocate_pairwise_neighbor_data( pairwise_neighbor_data_cutoff , n_cutoff, size_vdw_parameter )
          
          ! transfer data from pairwise_neighbor_data_verlet to
          ! pairwise_neighbor_data_cutoff using cutoff_mask
          call apply_cutoff_mask( n_neighbors, pairwise_neighbor_data_cutoff, pairwise_neighbor_data_verlet, cutoff_mask )

          ! deallocate old arrays that we don't need anymore
          call deallocate_pairwise_neighbor_data( pairwise_neighbor_data_verlet )
          deallocate( cutoff_mask )


          ! now calculate pairwise forces and energies for this atom and its neighbors         
          ! all of these subroutines should vectorize...
          call pairwise_real_space_ewald( E_elec_local , pairwise_neighbor_data_cutoff%f_ij ,  pairwise_neighbor_data_cutoff%dr,  pairwise_neighbor_data_cutoff%dr2,  pairwise_neighbor_data_cutoff%qi_qj, erf_factor , alpha_sqrt, PME_data%erfc_table , PME_data%erfc_grid , PME_data%erfc_max, constants%conv_e2A_kJmol )  
          call pairwise_real_space_LJ( E_vdw_local , pairwise_neighbor_data_cutoff%f_ij ,  pairwise_neighbor_data_cutoff%dr, pairwise_neighbor_data_cutoff%dr2 ,  pairwise_neighbor_data_cutoff%atype_vdw_parameter )
 

          E_elec = E_elec + E_elec_local
          E_vdw  = E_vdw  + E_vdw_local


          ! running addition of forces...
          do i_index=1, n_cutoff
             j_atom = pairwise_neighbor_data_cutoff%atom_index(i_index)
             local_force(:,i_atom) =  local_force(:,i_atom) + pairwise_neighbor_data_cutoff%f_ij(:,i_index)
             local_force(:,j_atom) =  local_force(:,j_atom) - pairwise_neighbor_data_cutoff%f_ij(:,i_index) 
          enddo

          ! deallocate old arrays that we don't need anymore
          call deallocate_pairwise_neighbor_data( pairwise_neighbor_data_cutoff)
         
       end if
    end do

    !$OMP END DO NOWAIT
    thread_id = OMP_GET_THREAD_NUM()
    temp_force(:,:,thread_id+1)= local_force(:,:)
    !$OMP END PARALLEL


    do i_thread=1,n_threads
       atom_data%force = atom_data%force + temp_force(:,:,i_thread)
    enddo

    system_data%E_elec = system_data%E_elec + E_elec
    system_data%E_vdw = system_data%E_vdw + E_vdw


    contains



    !***** this subroutine gathers neighbor data from global atom_data array using the Verlet neighbor list
    subroutine gather_neighbor_data( i_atom, verlet_start , verlet_finish , neighbor_list, pairwise_neighbor_data, xyz, charge, atom_type_index, atype_vdw_parameter )
       integer, intent(in) :: i_atom , verlet_start, verlet_finish
       integer, dimension(:), intent(in) :: neighbor_list
       type(pairwise_neighbor_data_type),intent(inout) :: pairwise_neighbor_data
       real*8, dimension(:,:) , intent(in) :: xyz
       real*8, dimension(:)   , intent(in) :: charge
       integer, dimension(:) , intent(in)  :: atom_type_index
       real*8, dimension(:,:,:) , intent(in) :: atype_vdw_parameter

       integer :: j_verlet, j_atom, i_type, j_type, j_index
       real*8  :: q_i , q_j

       ! type of central atom i_atom
       i_type = atom_type_index(i_atom)
       ! charge of central atom i_atom
       q_i = charge(i_atom)

       j_index=1
       ! loop over Verlet neighbors
       do j_verlet=verlet_start, verlet_finish
           j_atom = neighbor_list(j_verlet)
           j_type = atom_type_index(j_atom)
           q_j = charge(j_atom)

           pairwise_neighbor_data%atom_index(j_index) = j_atom
           pairwise_neighbor_data%xyz(:,j_index) = xyz(:,j_atom)
           pairwise_neighbor_data%qi_qj(j_index) = q_i * q_j
           pairwise_neighbor_data%atype_vdw_parameter(:,j_index) = atype_vdw_parameter(i_type,j_type,:)

           j_index = j_index + 1
       end do

    end subroutine gather_neighbor_data

   

  end subroutine pairwise_real_space_verlet








  !*******************************************
  ! this subroutine computes intra-molecular
  ! real-space interactions, taking into account
  ! intra-molecular exclusions
  !
  ! the data structure 'single_molecule_data' contains pointers to
  ! the molecule block of associated global atom_type data structures, and so when
  ! we manipulate these quantities we are manipulating the global data 
  !
  ! We input/output local energy, force data structures, as this will be used in
  ! MS-EVB routines and we don't necessarily want to update the
  ! global force, energy data.
  !*******************************************
  subroutine intra_molecular_pairwise_energy_force( force_local, E_elec, E_vdw, single_molecule_data , i_mole_type , n_atom, PME_data  )
    use global_variables
    real*8, dimension(:,:), intent(inout) :: force_local
    real*8, intent(inout)                 :: E_elec, E_vdw
    type(single_molecule_data_type), intent(inout) :: single_molecule_data
    type(PME_data_type) , intent(in)    :: PME_data
    integer, intent(in)    ::   i_mole_type, n_atom    

    ! we use two data structures, one for excluded interactions and one for non-excluded interactions
    Type(pairwise_neighbor_data_type) :: pairwise_neighbor_data_excluded , pairwise_neighbor_data_nonexcluded

    integer :: i_atom, j_atom, i_type, j_type, i_index
    integer :: index_excluded, index_nonexcluded, n_excluded, n_nonexcluded, size_vdw_parameter
    real*8 :: alpha_sqrt , erf_factor, E_elec_local, E_vdw_local 

    alpha_sqrt           = PME_data%alpha_sqrt
    erf_factor           =  2.D0*alpha_sqrt/constants%pi_sqrt
    ! this is (maximum) number of parameters that we need for VDWs interaction
    size_vdw_parameter=size(atype_vdw_parameter(1,1,:))


    ! by convention, molecules are not broken up over pbc, so we don't need to
    ! consider shift or minimum image for these intra-molecular interactions

    do i_atom=1, n_atom

       n_excluded=0
       n_nonexcluded=0
       E_elec_local=0d0
       E_vdw_local=0d0

       ! get number of exclusions, non-exclusions between this and rest of atoms in molecule
       do j_atom = i_atom + 1, n_atom
          if ( molecule_type_data(i_mole_type)%pair_exclusions( i_atom, j_atom ) /= 1 ) then
             n_nonexcluded = n_nonexcluded + 1  ! not excluded
          else
             n_excluded = n_excluded + 1    ! excluded
          end if
       enddo

       ! allocate data structures for excluded/non-excluded interactions
       if ( n_excluded > 0 ) then
           call allocate_pairwise_neighbor_data( pairwise_neighbor_data_excluded , n_excluded, size_vdw_parameter )       
       endif
       if ( n_nonexcluded > 0 ) then
           call allocate_pairwise_neighbor_data( pairwise_neighbor_data_nonexcluded , n_nonexcluded, size_vdw_parameter )
       endif

       i_type = single_molecule_data%atom_type_index(i_atom)
       ! now loop over atoms and fill in data structures
       index_excluded=1
       index_nonexcluded=1
       do j_atom = i_atom + 1, n_atom
          j_type = single_molecule_data%atom_type_index(j_atom)
          if ( molecule_type_data(i_mole_type)%pair_exclusions( i_atom, j_atom ) /= 1 ) then  
              ! nonexcluded            
              pairwise_neighbor_data_nonexcluded%atom_index(index_nonexcluded) = j_atom
              pairwise_neighbor_data_nonexcluded%dr(:,index_nonexcluded) = single_molecule_data%xyz(:,i_atom) - single_molecule_data%xyz(:,j_atom)
              pairwise_neighbor_data_nonexcluded%qi_qj(index_nonexcluded)= single_molecule_data%charge(i_atom) * single_molecule_data%charge(j_atom)
              pairwise_neighbor_data_nonexcluded%dr2(index_nonexcluded) = dot_product( pairwise_neighbor_data_nonexcluded%dr(:,index_nonexcluded), pairwise_neighbor_data_nonexcluded%dr(:,index_nonexcluded) )
              ! at this stage, all lj parameters should be expressed as C12
              ! and C6, even though they were read in as epsilon and sigma
              ! if this is a 1-4 interaction, take parameters from 1-4
              ! interaction parameter array
              if ( molecule_type_data(i_mole_type)%pair_exclusions( i_atom , j_atom ) == 2 ) then
                 pairwise_neighbor_data_nonexcluded%atype_vdw_parameter(:,index_nonexcluded) = atype_vdw_parameter_14(i_type,j_type,:)
              else
                 pairwise_neighbor_data_nonexcluded%atype_vdw_parameter(:,index_nonexcluded) = atype_vdw_parameter(i_type,j_type,:)
              end if
              index_nonexcluded = index_nonexcluded + 1
           else
           ! excluded, don't need data structures for VDWs interaction, just electrostatics
              pairwise_neighbor_data_excluded%atom_index(index_excluded) = j_atom 
              pairwise_neighbor_data_excluded%dr(:,index_excluded) = single_molecule_data%xyz(:,i_atom) - single_molecule_data%xyz(:,j_atom)
              pairwise_neighbor_data_excluded%qi_qj(index_excluded)= single_molecule_data%charge(i_atom) * single_molecule_data%charge(j_atom)
              pairwise_neighbor_data_excluded%dr2(index_excluded) = dot_product( pairwise_neighbor_data_excluded%dr(:,index_excluded), pairwise_neighbor_data_excluded%dr(:,index_excluded) )
              index_excluded = index_excluded + 1
           endif
       enddo

       ! first remove intra-molecular PME contributions for excluded atoms
       if ( n_excluded > 0 ) then
           call intra_pme_exclusion( E_elec_local , pairwise_neighbor_data_excluded%f_ij ,  pairwise_neighbor_data_excluded%dr, pairwise_neighbor_data_excluded%dr2, pairwise_neighbor_data_excluded%qi_qj, erf_factor , alpha_sqrt , constants%conv_e2A_kJmol )
           E_elec = E_elec + E_elec_local
       endif
           
       ! now add non-excluded intra-molecular real-space electrostatic and VDWs interactions
       if ( n_nonexcluded > 0 ) then
           call pairwise_real_space_ewald( E_elec_local , pairwise_neighbor_data_nonexcluded%f_ij ,  pairwise_neighbor_data_nonexcluded%dr, pairwise_neighbor_data_nonexcluded%dr2,  pairwise_neighbor_data_nonexcluded%qi_qj, erf_factor , alpha_sqrt, PME_data%erfc_table , PME_data%erfc_grid , PME_data%erfc_max, constants%conv_e2A_kJmol )
           call pairwise_real_space_LJ( E_vdw_local , pairwise_neighbor_data_nonexcluded%f_ij ,  pairwise_neighbor_data_nonexcluded%dr, pairwise_neighbor_data_nonexcluded%dr2 , pairwise_neighbor_data_nonexcluded%atype_vdw_parameter )
           E_elec = E_elec + E_elec_local
           E_vdw  = E_vdw  + E_vdw_local
       endif

       ! add forces from exclusions...
       do i_index=1, n_excluded
           j_atom = pairwise_neighbor_data_excluded%atom_index(i_index)
           force_local(:,i_atom) =  force_local(:,i_atom) + pairwise_neighbor_data_excluded%f_ij(:,i_index)
           force_local(:,j_atom) =  force_local(:,j_atom) - pairwise_neighbor_data_excluded%f_ij(:,i_index)
       enddo

       ! add forces from non-excluded interactions...
       do i_index=1, n_nonexcluded
           j_atom = pairwise_neighbor_data_nonexcluded%atom_index(i_index)
           force_local(:,i_atom) = force_local(:,i_atom) + pairwise_neighbor_data_nonexcluded%f_ij(:,i_index)
           force_local(:,j_atom) = force_local(:,j_atom) - pairwise_neighbor_data_nonexcluded%f_ij(:,i_index)
       enddo


       ! deallocate data structures
       if ( n_nonexcluded > 0 ) then
          call deallocate_pairwise_neighbor_data(pairwise_neighbor_data_nonexcluded)
       endif
       if ( n_excluded > 0 ) then
          call deallocate_pairwise_neighbor_data(pairwise_neighbor_data_excluded)
       endif

    enddo


  end subroutine intra_molecular_pairwise_energy_force





    !*********** this subroutine uses the cutoff_mask to reorganize data
    !structures and only keep the neighbor atoms within the cutoff distance
    subroutine apply_cutoff_mask( n_neighbors, pairwise_neighbor_data_cutoff, pairwise_neighbor_data_verlet, cutoff_mask )
      integer,intent(in) :: n_neighbors
      type(pairwise_neighbor_data_type) , intent(inout) ::  pairwise_neighbor_data_cutoff
      type(pairwise_neighbor_data_type) ,  intent(in) ::  pairwise_neighbor_data_verlet
      integer, dimension(:), intent(in) :: cutoff_mask

      integer :: i_atom, i_index

      i_index=1
      do i_atom=1, n_neighbors
         if ( cutoff_mask(i_atom) == 1 ) then
            pairwise_neighbor_data_cutoff%dr(:,i_index) = pairwise_neighbor_data_verlet%dr(:,i_atom)
            pairwise_neighbor_data_cutoff%dr2(i_index) =  pairwise_neighbor_data_verlet%dr2(i_atom)
            pairwise_neighbor_data_cutoff%qi_qj(i_index) = pairwise_neighbor_data_verlet%qi_qj(i_atom)
            pairwise_neighbor_data_cutoff%atype_vdw_parameter(:,i_index) = pairwise_neighbor_data_verlet%atype_vdw_parameter(:,i_atom)
            pairwise_neighbor_data_cutoff%atom_index(i_index) = pairwise_neighbor_data_verlet%atom_index(i_atom)
         i_index = i_index + 1
         endif
      enddo

    end subroutine apply_cutoff_mask




 !********************************************
 ! this subroutine computes force and energy contributions from pairwise
 ! Lennard-Jones (LJ) interactions
 !********************************************
  subroutine pairwise_real_space_LJ( E_vdw , f_ij ,  dr, dr2, lj_parameters )
     real*8, intent(out) :: E_vdw
     real*8, dimension(:,:), intent(inout) :: f_ij
     real*8, dimension(:), intent(in)  :: dr2
     real*8, dimension(:,:), intent(in) :: dr
     real*8, dimension(:,:), intent(in) :: lj_parameters

     real*8, dimension(size(dr2)) :: dr6 , dr12
     integer :: i_atom

     ! $omp simd
     dr6 = dr2**3
     dr12 = dr6**2

     ! at this stage, all lj parameters should be expressed as C12 and C6, even though they were read in as epsilon and sigma
     E_vdw = sum( lj_parameters(1,:) / dr12(:) - lj_parameters(2,:) / dr6(:) )

     ! this should vectorize...
     do i_atom=1,size(dr2)
        f_ij(:,i_atom) = f_ij(:,i_atom) + dr(:,i_atom) / dr2(i_atom) * ( 12d0 * lj_parameters(1,i_atom) / dr12(i_atom) - 6d0 * lj_parameters(2,i_atom) / dr6(i_atom) )
     enddo

     ! $omp end simd

  end subroutine pairwise_real_space_LJ





  !**************************************************
  ! this subroutine computes the real space Ewald (PME) force and energy
  ! contribution for pairwise interactions
  !**************************************************
  subroutine pairwise_real_space_ewald( E_elec , f_ij ,  dr, dr2, qi_qj, erf_factor , alpha_sqrt, erfc_table , erfc_grid , erfc_max, conv_e2A_kJmol )
     real*8, intent(out) :: E_elec
     real*8, dimension(:,:), intent(inout) :: f_ij
     real*8, dimension(:), intent(in)  :: dr2 , qi_qj
     real*8, dimension(:,:), intent(in) :: dr
     real*8, dimension(:),  intent(in) :: erfc_table
     integer, intent(in)                :: erfc_grid
     real*8, intent(in)                 :: erfc_max, alpha_sqrt, erf_factor, conv_e2A_kJmol

     real*8, dimension(size(dr2)) :: dr_mag, erfc_value
     integer :: i_atom
  
     ! $omp simd

     dr_mag=sqrt(dr2)
     ! make array with values of the complementary error function applied to alpha_sqrt * distance
     erfc_value(:) = erfc_table(ceiling((dr_mag(:)*alpha_sqrt)/erfc_max*dble(erfc_grid)))

     E_elec = sum( qi_qj / dr_mag * erfc_value * conv_e2A_kJmol )  ! conversion is e^2/Angstrom to kJ/mol

     ! this should vectorize...
     do i_atom=1, size(dr2)
        f_ij(:,i_atom)  = f_ij(:,i_atom) + qi_qj(i_atom) / dr_mag(i_atom) * dr(:,i_atom) * ( erfc_value(i_atom) / dr2(i_atom) + erf_factor * exp(-(alpha_sqrt * dr_mag(i_atom)) **2) / dr_mag(i_atom) ) * conv_e2A_kJmol
     enddo

     ! $omp end simd

  end subroutine pairwise_real_space_ewald





  !************************************************
  ! this subroutine removes intra-molecular PME interactions
  ! from the reciprocal space sum
  !
  ! use explicit erfc since r can be small ( this is still in code from
  ! Drude oscillator treatment, probably don't need anymore without Drude
  ! oscillators)
  !
  ! Note that if we wish to subtract out the intra molecular interaction between two charges on top 
  ! of each other, (suppose drude oscillators are directly on top of their corresponding atoms)
  ! we need to subtract this contribution to the reciprocal space energy.  In the reciprocal
  ! space contribution, there is a factor of 1/2 for overcounting, so this contribution is
  !  1/2 * ( E12int + E21int ) = E12int , where E12int is charge 1 interacting with a gaussian
  ! of charge 2 right underneath it.  This energy is then
  ! 2 * q1 * q2 * ( alpha / pi ) ^ (1/2)  
  ! note that this is a factor of 2 greater than the Ewald self correction, because there are two
  ! interactions that we are accounting for, rather than one.
  !************************************************

  subroutine intra_pme_exclusion( E_elec_local , f_ij ,  dr,  dr2, qi_qj, erf_factor , alpha_sqrt , conv_e2A_kJmol )
    implicit none
    real*8, intent(out) :: E_elec_local
    real*8, dimension(:,:), intent(inout) :: f_ij 
    real*8, dimension(:,:), intent(in)  :: dr
    real*8, dimension(:), intent(in)    :: dr2
    real*8, dimension(:), intent(in)    :: qi_qj
    real*8, intent(in)                  :: erf_factor, alpha_sqrt, conv_e2A_kJmol

    real*8,parameter::small=1D-8
    real*8, dimension(size(dr2)) :: dr_mag, dr_mag3
    integer :: i_atom

    dr_mag=sqrt(dr2)
    dr_mag3=dr2*dr_mag

    E_elec_local=0d0

    ! remove charge-gaussion interaction, the if-statements are probably not
    ! necessary, problem occurs only if dr_mag is numerically zero, which only
    ! would happen if there are Drude oscillators that are sitting right on top
    ! of atom.  We leave this check in for safety, as intra-molecular
    ! interactions are not rate-limiting computation

    do i_atom=1, size(dr2)
       if( dr_mag(i_atom) < small ) then
          E_elec_local = E_elec_local - erf_factor * qi_qj(i_atom) * conv_e2A_kJmol  ! conversion is e^2/Angstrom to kJ/mol
          ! no contribution to force
       else
          E_elec_local = E_elec_local + qi_qj(i_atom) * (erfc(dr_mag(i_atom)*alpha_sqrt)-1.D0) / dr_mag(i_atom) * conv_e2A_kJmol
          f_ij(:,i_atom) = f_ij(:,i_atom) + qi_qj(i_atom) * dr(:,i_atom) * ( (erfc(dr_mag(i_atom)*alpha_sqrt)-1.D0) / dr_mag3(i_atom) + erf_factor*exp(-(dr_mag(i_atom)*alpha_sqrt)**2)/dr2(i_atom)) * conv_e2A_kJmol       
       endif
    enddo


  end subroutine intra_pme_exclusion





end module pair_int_real_space
