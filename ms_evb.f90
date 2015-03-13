module ms_evb
  use routines
  use global_variables

  !************************************
  ! this is the module to implement multi-state empirical valance bond force 
  ! fields.
  ! 
  ! currently, this is customized to implement the MS-EVB3 model of Voth and coworkers
  ! JPC B, 2008, 112, 467-482 (and errata, JPC B, 2008, 112, 7146)
  !
  ! essentially, this module is in charge of constructing the EVB Hamiltonian
  ! and diagonalizing this Hamiltonian to find the ground state adiabat, 
  ! and then getting the forces from the Hellman-Feynman theorem
  !
  ! To efficiently construct the diabats, some PME tricks need to be used
  ! (i.e. calculating the electrostatic potential on a grid)
  !************************************



  ! this maps the water and hydronium atom types after proton hop
  integer, dimension(MAX_N_ATOM_TYPE) :: map_atom_index

  ! based on the atom_type_index, this constructs a matrix detailing the potential evb reactivity
  ! note this matrix is NOT symmetric.  For example for hydrogen transfer, we are always transfering
  ! a proton from H3O+ to H2O, and so 
  integer, dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE) :: evb_reactivity_matrix  

  real*8, parameter :: evb_first_solvation_cutoff = 5d0
  real*8, parameter :: evb_reactive_pair_distance = 2.5d0
  integer, parameter :: evb_max_neighbors=10  ! this controls the dimensions of some of the evb data structures
  ! evb_max_states is the maximum dimension for evb-hamiltonian, and therefore should be set to the maximum number of diabats.  
  integer, parameter :: evb_max_states=50 

  ! maximum size of water chain for proton hopping.  Note this may need to be larger than the number
  ! of solvation shells for proton transfer.  For instance in a water hexamer, for a particular initial
  ! proton configuration, may need 6 hops to generate all diabats
  integer, parameter :: evb_max_chain=6  
  real*8 , dimension(evb_max_states,evb_max_states)  :: evb_hamiltonian
  integer, dimension(evb_max_states,evb_max_states)  :: evb_forces_lookup_index
  real*8, dimension(:,:,:,:), allocatable :: evb_forces_store
  ! this stores the information about the proton hops necessary to get from the principle diabat to the specified diabat
  ! data is stored in groups of 2 integers, specifying index of hydrogen atom proton donor water molecule, and 
  ! index of proton acceptor water molecule.  The index of proton donor water molecule is always known, because
  ! the principle hydronium is the first donor, and the last acceptor is the next donor in the chain
  integer, dimension(evb_max_states, 2 * evb_max_chain ) :: evb_diabat_proton_log

  ! this is a counter that is constantly updated with every recursive call to evb_conduct_proton_transfer_recursive subroutine
  integer :: diabat_index

contains

  !************************
  ! initialization stuff involved with MS-EVB simulation
  !***********************
  subroutine initialize_evb_simulation( n_mole, n_atom )
    integer,intent(in) :: n_mole
    integer, dimension(:),intent(in) :: n_atom   

    ! always match strings with spaces at end
    call trim_end(hydronium_proton_label)
    call trim_end(hydronium_oxygen_label)
    call trim_end(water_proton_label)
    call trim_end(water_oxygen_label)

    ! update hydronium molecule list
    call update_hydronium_molecule_index( n_mole, n_atom )

    ! check consistency
    call evb_consistency_checks( n_mole, n_atom )

    ! form reactivity matrix
    call construct_evb_reactivity_matrix


  end subroutine initialize_evb_simulation



  !************************
  ! this subroutine updates the hydronium_molecule_index array
  ! with the indices of hydronium molecules
  !************************
  subroutine update_hydronium_molecule_index( n_mole, n_atom )
    use global_variables
    integer,intent(in) :: n_mole
    integer, dimension(:),intent(in) :: n_atom   

    integer :: i_mole, i_atom, index, i_hydronium

    i_hydronium = 1
    do i_mole = 1 , n_mole
       ! see if molecule is hydronium 
       index = index_first_atom_type( i_mole , n_atom, hydronium_oxygen_label )
       if ( index > 0 ) then
          hydronium_molecule_index( i_hydronium ) = i_mole
          i_hydronium = i_hydronium + 1
       end if
    enddo

    n_hydronium_molecules = i_hydronium - 1

    ! need at least 1 hydronium
    if ( n_hydronium_molecules < 1 ) then
       stop "need at least one hydronium molecule!"
    endif

    ! can't have more than 1 hydronium
    if ( n_hydronium_molecules > 1 ) then
       stop "can't have more than 1 hydronium, see code comments"
       ! if more than one hydronium, need to tell subroutines
       ! which hydronium we're interested in
       ! in particular, need to change evb_create_diabat_data_structures
       ! to deal with more than one hydronium
    endif

  end subroutine update_hydronium_molecule_index



  !*************************
  ! we make use of some specific input requirements for a ms-evb simulation
  ! to allow easier use and manipulation of data structures.  This subroutine 
  ! checks for such requirements
  !*************************
  subroutine evb_consistency_checks( n_mole, n_atom )
    integer,intent(in) :: n_mole
    integer, dimension(:),intent(in) :: n_atom

    integer              :: i_mole,i_atom,i_type
    character(MAX_ANAME) :: atomname


    write(*,*) ""
    write(*,*) "For MS-EVB simulation, we require specific labels for hydronium"
    write(*,*) "and water atoms.  The required atomtypes are "
    write(*,*) "Hydronium proton ", hydronium_proton_label
    write(*,*) "Hydronium oxygen ", hydronium_oxygen_label
    write(*,*) "water proton ", water_proton_label
    write(*,*) "water oxygen ", water_oxygen_label
    write(*,*) ""

    ! We will be adding and removing hydrogen atoms to and from
    ! hydronium and water molecules during proton transfer.
    ! This requires manipulation of data structures.
    ! However we need all molecules of the same type to have the
    ! same ordering of atoms, so that the same moleculetype data
    ! structures can be used for all like molecules.
    ! Therefore, we require that all water and hydronium molecules
    ! have atoms ordered O, H , H , (H) , so that the removal/addition
    ! of a proton does not effect the atom ordering.

    do i_mole=1,n_mole
       do i_atom=1,n_atom(i_mole)
          i_type = atom_index(i_mole, i_atom)
          atomname = atype_name(i_type)
          if ( ( atomname .eq. hydronium_proton_label ) .and. ( i_atom == 1 ) ) then
             write(*,*) ""
             write(*,*) "hydronium molecule topology must have oxygen atom first"
             write(*,*) ""
             stop
          else if ( ( atomname .eq. hydronium_oxygen_label ) .and. ( i_atom /= 1 ) ) then
             write(*,*) ""
             write(*,*) "hydronium molecule topology must have oxygen atom first"
             write(*,*) ""
             stop
          end if
          if ( ( atomname .eq. water_proton_label ) .and. ( i_atom == 1 ) ) then
             write(*,*) ""
             write(*,*) "water molecule topology must have oxygen atom first"
             write(*,*) ""
             stop
          else if ( ( atomname .eq. water_oxygen_label ) .and. ( i_atom /= 1 ) ) then
             write(*,*) ""
             write(*,*) "water molecule topology must have oxygen atom first"
             write(*,*) ""
             stop
          end if
       end do
    end do


  end subroutine evb_consistency_checks



  !*************************
  ! this subroutine calculates the total force and energy
  ! for an ms-evb hamiltonian
  !
  ! adiabatic_force and adiabatic_energy are the output force and energy
  ! solved for using the MS-EVB Hamiltonian
  !
  ! upon a change in the identity of the hydronium ion, the following data structures are
  ! changed with the call to evb_create_diabat_data_structures subroutine
  ! local :: xyz, r_com, chg, n_atom, mass
  ! global :: hydronium_molecule_index, atom_index, molecule_index
  !*************************
  subroutine ms_evb_calculate_total_force_energy( adiabatic_force, adiabatic_energy, tot_n_mole, n_mole, n_atom, xyz, r_com, chg, mass, box, dfti_desc,dfti_desc_inv,log_file )
    use MKL_DFTI
    real*8,dimension(:,:,:),intent(out) :: adiabatic_force
    real*8 , intent(out)  :: adiabatic_energy
    integer, intent(in) :: tot_n_mole, n_mole
    integer, intent(inout), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: box
    real*8, intent(inout), dimension(:,:) :: r_com, chg, mass
    real*8, intent(inout), dimension(:,:,:) :: xyz
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv
    integer,intent(in)::log_file

    integer :: i_mole, i_mole_principle, i_mole_hydronium, new_i_mole_hydronium, principle_diabat
    ! temporary diabat data structures
    real*8, dimension(:,:,:), allocatable :: xyz_temp1
    real*8, dimension(:,:), allocatable :: r_com_temp1
    real*8, dimension(:,:), allocatable :: chg_temp1, mass_temp1
    integer, dimension(:), allocatable :: n_atom_temp1
    !****** this is junk variables to pass  *********
    integer, dimension(MAX_N_MOLE) :: n_atom_drude

    ! test
    integer :: i_atom, i_state, j_state, i
    real*8 ,dimension(:,:),allocatable :: hamiltonian
    real*8 , dimension(:,:,:,:,:), allocatable :: force_store
    real*8, dimension(3) :: dr
    real*8,dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: xyz_test
    real*8 :: dE, E_new

    ! not test
    n_atom_drude = n_atom
    ! make sure we are not using verlet list, as we need to update verlet list when
    ! changing diabats, and this isn't implemented yet
    Select Case(use_verlet_list)
       Case("yes")
          stop "can't run ms_evb with verlet list yet"
    end Select

!!$    !*********** test
!!$
!!$    i_mole_hydronium = hydronium_molecule_index( 1 )
!!$    call construct_evb_hamiltonian( i_mole_hydronium, tot_n_mole, n_mole, n_atom, xyz, r_com, chg, mass, box, dfti_desc,dfti_desc_inv,log_file )
!!$
!!$    allocate( hamiltonian(diabat_index, diabat_index), force_store(diabat_index, diabat_index, MAX_N_MOLE, MAX_N_ATOM,3 ) )
!!$    ! store initial hamiltonian
!!$    do i_state = 1 , diabat_index
!!$       do j_state = i_state , diabat_index
!!$          hamiltonian( i_state, j_state ) = evb_hamiltonian(i_state,j_state)
!!$          hamiltonian(j_state, i_state ) = hamiltonian(i_state, j_state )
!!$       end do
!!$    end do
!!$    ! store initial force
!!$    force_store=0d0
!!$    do i_state=1, diabat_index
!!$       do j_state = i_state, diabat_index
!!$          index = evb_forces_lookup_index(i_state,j_state)
!!$          if ( index > 0 ) then
!!$             ! add forces from this matrix element
!!$             do i_mole = 1, tot_n_mole
!!$                do i_atom = 1 , n_atom(i_mole)
!!$                   force_store(i_state, j_state,i_mole,i_atom,:) =  evb_forces_store(index, i_mole, i_atom, :)
!!$                enddo
!!$             enddo
!!$          end if
!!$       end do
!!$    end do
!!$
!!$
!!$    ! now test forces
!!$    do i_mole=1,n_mole
!!$       do i_atom=1, n_atom(i_mole)
!!$          write(*,*) "testing forces for atom ", i_mole, i_atom
!!$          do i=1,3
!!$             call random_number(dr(i))
!!$          enddo
!!$          dr = dr * 0.001d0
!!$          xyz_test = xyz
!!$          xyz_test(i_mole,i_atom,:) = xyz(i_mole,i_atom,:) + dr(:)
!!$          r_com(i_mole,:) = pos_com( xyz_test, i_mole, n_atom, mass )
!!$
!!$          call construct_evb_hamiltonian( i_mole_hydronium, tot_n_mole, n_mole, n_atom, xyz_test, r_com, chg, mass, box, dfti_desc,dfti_desc_inv,log_file )      
!!$          do i_state =1, diabat_index
!!$             do j_state = i_state, diabat_index
!!$                dE = -dot_product( force_store(i_state,j_state,i_mole,i_atom,:) , dr(:) )
!!$                E_new = hamiltonian(i_state, j_state) + dE
!!$                write(*,'(I5,I5,F12.6,F12.6,F12.6)') i_state, j_state, hamiltonian(i_state, j_state), E_new , evb_hamiltonian(i_state,j_state)
!!$                write(*,*) ( evb_diabat_proton_log( i_state , i ) , i=1,3 )
!!$             enddo
!!$          enddo
!!$       enddo
!!$    enddo
!!$    stop
!!$
!!$
!!$    !************* end test



    if ( n_hydronium_molecules > 1 ) then
       stop "can't have more than 1 hydronium"
    end if

    ! loop over hydronium molecules
    do i_mole=1, n_hydronium_molecules
       i_mole_hydronium = hydronium_molecule_index( i_mole )
       call construct_evb_hamiltonian( i_mole_hydronium, tot_n_mole, n_mole, n_atom, xyz, r_com, chg, mass, box, dfti_desc,dfti_desc_inv,log_file )

       call diagonalize_evb_hamiltonian( principle_diabat, new_i_mole_hydronium,  i_mole_hydronium, adiabatic_energy, adiabatic_force, tot_n_mole, n_atom )

       ! if principle hydronium molecule has changed, then change
       ! data structures.  Note we need to pass temporary arrays here, and then copy them back
       ! this is because the arrays passed in with intent(in), or intent(out) attribute will be passed in
       ! as pointers, so we can't pass (e.g.) n_atom, n_atom, without the first intent in array changing.
       ! which will cause a bug
       if ( new_i_mole_hydronium /= i_mole_hydronium ) then
          allocate( xyz_temp1(MAX_N_MOLE,MAX_N_ATOM,3), r_com_temp1(MAX_N_MOLE,3), chg_temp1(MAX_N_MOLE,MAX_N_ATOM), mass_temp1(MAX_N_MOLE,MAX_N_ATOM),n_atom_temp1(MAX_N_MOLE) )
          call evb_create_diabat_data_structures( xyz, xyz_temp1, r_com, r_com_temp1, chg, chg_temp1 , n_atom, n_atom_temp1, n_atom_drude, mass, mass_temp1, box, principle_diabat, i_mole_hydronium )
          xyz = xyz_temp1  
          r_com = r_com_temp1
          chg = chg_temp1
          mass = mass_temp1
          n_atom = n_atom_temp1
          deallocate( xyz_temp1, r_com_temp1, chg_temp1, mass_temp1, n_atom_temp1 )
       endif

    enddo

    ! deallocate stored force array
    deallocate( evb_forces_store )

  end subroutine ms_evb_calculate_total_force_energy




  !**************************
  ! this subroutine calculates the adiabatic energy
  ! and force for the ms-evb model by diagonalizing the
  ! evb hamiltonian
  !**************************
  subroutine diagonalize_evb_hamiltonian( principle_diabat, new_i_mole_principle, i_mole_principle, adiabatic_potential, force_atoms, tot_n_mole, n_atom )
    integer, intent(out) :: principle_diabat, new_i_mole_principle
    integer, intent(in)  :: i_mole_principle
    real*8 , intent(out) :: adiabatic_potential
    real*8, dimension(:,:,:), intent(out) :: force_atoms
    integer, intent(in) :: tot_n_mole
    integer, dimension(:) , intent(in) :: n_atom

    integer :: i_state, j_state, n_rot, ground_state_index, index, i_hop, i_mole, i_atom
    real*8 :: factor, state_coefficient
    real*8, dimension(:), allocatable  :: eigenvalues, ground_state_eigenvector
    real*8, dimension(:,:), allocatable :: hamiltonian, eigenvectors

    ! when the evb_hamiltonian was constructed, it was constructed as an upper triangular matrix
    ! fill in lower triangle.  Also we most probably have fewer diabats than the setting of
    ! evb_max_states, in which case the evb_hamiltonian is partially empty
    ! the number of total diabats is stored in the diabat_index variable

    allocate( eigenvalues(diabat_index) , ground_state_eigenvector(diabat_index), hamiltonian(diabat_index, diabat_index) , eigenvectors(diabat_index,diabat_index) )

    hamiltonian=0d0
!!$    write(*,*) "evb hamiltonian"
    do i_state = 1 , diabat_index
       do j_state = i_state , diabat_index
          hamiltonian( i_state, j_state ) = evb_hamiltonian(i_state,j_state)
          hamiltonian(j_state, i_state ) = hamiltonian(i_state, j_state )
       end do
         !write(*,'(I5,6F10.3)') i_state, hamiltonian( i_state, : )
!!$         write(*,*) i_state, hamiltonian( i_state, : )
    end do

    ! diagonalize
    call jacobi(hamiltonian,eigenvalues,eigenvectors,n_rot)

    ! find ground state
    ground_state_index=1
    adiabatic_potential = eigenvalues(1)
    do i_state = 2, diabat_index
       if ( eigenvalues(i_state) < adiabatic_potential ) then
          adiabatic_potential = eigenvalues(i_state)
          ground_state_index = i_state
       end if
    end do
    ground_state_eigenvector(:)=eigenvectors(:,ground_state_index)


    ! now get force using Hellmann-Feynman theorem
    force_atoms=0d0
    do i_state=1, diabat_index
       do j_state = i_state, diabat_index
          ! non-zero coupling here if evb_forces_lookup_index(i_state, j_state) > 0
          index = evb_forces_lookup_index(i_state,j_state)
          if ( index > 0 ) then
             if ( i_state /= j_state ) then
                ! factor of 2 to account for transpose coupling element
                factor=2d0 * ground_state_eigenvector(i_state) * ground_state_eigenvector(j_state)
             else
                factor= ground_state_eigenvector(i_state) * ground_state_eigenvector(j_state)
             end if

             ! add forces from this matrix element
             do i_mole = 1, tot_n_mole
                do i_atom = 1 , n_atom(i_mole)
                   force_atoms(i_mole,i_atom,:) = force_atoms(i_mole,i_atom,:) + factor * evb_forces_store(index, i_mole, i_atom, :)
                enddo
             enddo

          end if
       end do
    end do


    ! test
!!$    write(*,*) "eigenvalues"
!!$    !write(*,'(6F10.3)') eigenvalues
!!$    write(*,*) eigenvalues
!!$    write(*,*) "ground state"
!!$    !write(*,'(6F10.3)') ground_state_eigenvector
!!$    write(*,*) ground_state_eigenvector
!!$    write(*,*) "forces"
!!$    do i_mole=1, tot_n_mole
!!$       do i_atom = 1 , n_atom(i_mole)
!!$          write(*,*) i_mole, i_atom, force_atoms(i_mole,i_atom,:)
!!$       enddo
!!$    enddo
!!$    stop
    ! end test

    ! now figure out which hydronium molecule corresponds to the new principle diabat
    principle_diabat=1
    state_coefficient = abs(ground_state_eigenvector(1))
    do i_state=1, diabat_index
       if ( state_coefficient < abs(ground_state_eigenvector(i_state)) ) then
          state_coefficient = abs(ground_state_eigenvector(i_state))
          principle_diabat = i_state
       end if
    enddo

    ! now get hydronium molecule for this diabat
    new_i_mole_principle = i_mole_principle
    index = 1
    do i_hop=1, evb_max_chain
       if ( evb_diabat_proton_log( principle_diabat , index ) < 0 ) exit 
       new_i_mole_principle = evb_diabat_proton_log( principle_diabat , index + 1 )
       index = index + 2
    enddo


    deallocate( eigenvalues , ground_state_eigenvector, hamiltonian , eigenvectors )

  end subroutine diagonalize_evb_hamiltonian



  !************************************
  ! This subroutine loops through all ms evb states, and constructs the evb hamiltonian matrix
  !
  ! i_mole_principle is the molecule that defines the principle diabat, namely the index of the hydronium
  ! ion for the diabatic state with the largest ci**2 coefficient in the adiabat
  !
  ! as we construct the matrix elements for the evb-hamiltonian, we also calculate the forces on all
  ! the atoms in the system for that particular hamiltonian matrix element (including the diabatic 
  ! coupling matrix elements!  However, many of these are zero) .  Because this is a large number of
  ! forces for a large number of matrix elements, we need to be smart about allocating storage. 
  ! We therefore use two different arrays, the first to store the lookup table index for all non-zero
  ! matrix elements, and the second to store the forces for this particular non-zero hamitonian element
  ! note that we store the forces in the constant molecular topology of the principle diabat, regardless
  ! of the particular matrix element
  ! 
  ! WARNING:  THIS SUBROUTINE MODIFIES GLOBAL VARIABLES! 
  ! Global array atom_index is modified in order to allow for energy calculation of diabat
  ! when subroutine evb_create_diabat_data_structures is called
  ! This change can be restored by calling evb_restore_global_data_structures
  !************************************
  subroutine construct_evb_hamiltonian( i_mole_principle, tot_n_mole, n_mole, n_atom, xyz, r_com, chg, mass, box, dfti_desc,dfti_desc_inv,log_file )
    use MKL_DFTI
    integer, intent(in) :: i_mole_principle, tot_n_mole, n_mole
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: box,r_com, chg, mass
    real*8, intent(in), dimension(:,:,:) :: xyz
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv
    integer,intent(in)::log_file

    integer ::  diabat_index_donor,  i_mole_donor


    ! initialize proton hop log
    evb_diabat_proton_log=-1
    ! zero hamiltonian
    evb_hamiltonian=0d0

    ! initialize diabat_index
    diabat_index=1
    diabat_index_donor=1
    i_mole_donor = i_mole_principle

    ! construct and store the evb hamiltonian matrix element for the principle diabat.  This subroutine will update the evb_hamiltonian and evb_forces_store arrays.
    call evb_hamiltonian_elements_donor_acceptor( diabat_index_donor, diabat_index_donor, i_mole_principle, tot_n_mole, n_mole, n_atom, xyz, r_com, chg, mass, box, dfti_desc,dfti_desc_inv,log_file )


    ! this is a recursive subroutine that loops through all diabatic states using recursive proton hops
    call evb_conduct_proton_transfer_recursive( i_mole_donor, diabat_index_donor,  i_mole_principle, tot_n_mole, n_mole, n_atom, xyz, r_com, chg, mass, box, dfti_desc,dfti_desc_inv,log_file )

  end subroutine construct_evb_hamiltonian




  !***********************************************
  ! this subroutine constructs all diabatic states using recursive proton hops
  ! and subsequently calls the evb_hamiltonian_elements_donor_acceptor subroutine
  ! to calculate and store the evb hamiltonian matrix elements for that
  ! particular donor-acceptor state
  !
  ! note that we don't run into the problem that an acceptor molecule back transfers to its donor
  ! when the acceptor molecule becomes the next donor, because in the find_evb_reactive_neighbors
  ! subroutine we always input the original principle "xyz" geometry array
  !***********************************************
  recursive subroutine evb_conduct_proton_transfer_recursive( i_mole_donor, diabat_index_donor,  i_mole_principle, tot_n_mole, n_mole, n_atom, xyz, r_com, chg, mass, box, dfti_desc,dfti_desc_inv,log_file )
    use MKL_DFTI
    integer, intent(in) :: i_mole_donor, diabat_index_donor, i_mole_principle, tot_n_mole, n_mole
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: box,r_com, chg, mass
    real*8, intent(in), dimension(:,:,:) :: xyz
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv
    integer,intent(in)::log_file

    integer, dimension(evb_max_neighbors) :: evb_neighbor_list
    integer :: i_molecule, i_atom, i_mole_acceptor, i_hop, diabat_index_acceptor, atom_id, flag_reactive_atom, flag_cycle, index, count


    ! first see if we've reached the maximum limit for the proton hopping chain, figure this out by the number of stored elements in proton hopping log
    index = 1
    count = 0
    do i_hop=1, evb_max_chain
       if ( evb_diabat_proton_log( diabat_index_donor , index ) < 0 ) exit
       count = count + 1
       index = index + 2
    enddo

    !*******************************************
    ! this search for possible proton hops merits some comments
    ! at this stage in the code, all the data structures are preserved
    ! in their initial values corresponding to the principle diabat (i.e. before any proton hops)
    ! data structures are not changed to reflect the specific proton hop
    ! until the energy is evaluated (in subroutine  evb_hamiltonian_elements_donor_acceptor )
    !
    ! therefore, ALL protons (water, hydronium) should be viewed as reactive atoms, because
    ! the protons that are involved in secondary hops are currently in the water topology
    ! in the data structures.  Also, in the subroutine "find_evb_reactive_neighbors",
    ! a hop from a proton to a hydronium molecule should be considered, as this is necessary 
    ! to preserve topology invariance in e.g. a pentamer ring.  Note this will not create
    ! identical diabats from back transfers, because the initial protons will always
    ! belong the initial hydronium, since initial data structures are preserved here
    !
    ! if we're in a ring, and we have a complete proton transfer cycle around a ring, (proton is transfered to initial hydronium) 
    !we need to stop the transfer cycle, so we don't keep going.  This is done by using "flag_cycle"
    !*******************************************

    ! if not yet at max hops, find acceptor for this donor
    if ( count < evb_max_chain ) then

       ! loop over all hydrogen atoms for this donor
       do i_atom=1, n_atom(i_mole_donor)
          atom_id = atom_index(i_mole_donor,i_atom)
          flag_reactive_atom = check_atom_reactivity(atom_id)

          ! if reactive atom (proton), consider all possible diabats formed from the transfer of this proton
          if ( flag_reactive_atom == 1 ) then

!!$             write(*,*) "reactive", i_mole_donor, i_atom

             ! first identify molecules within first solvation shell to which we can transfer a proton
             call find_evb_reactive_neighbors( i_mole_donor, i_atom, evb_neighbor_list, tot_n_mole, n_mole, n_atom, xyz,r_com, box )

             ! loop over neighbors
             loop1:          do i_molecule = 1 , evb_max_neighbors
                if ( evb_neighbor_list(i_molecule) < 0 ) exit loop1

                ! starting a new diabat, update global diabat_index counter
                diabat_index = diabat_index + 1
                diabat_index_acceptor = diabat_index

                ! make sure we don't have too many states
                call evb_check_number_diabats                 

                i_mole_acceptor = evb_neighbor_list(i_molecule)

                ! check to see if this is a hydronium molecule (from i.e. cyclic transfer)
                flag_cycle = check_hydronium_molecule( i_mole_acceptor )

!!$                write(*,*) "acceptor" , i_mole_acceptor

                ! create proton hop log for this acceptor, copy previous proton transfer log from donor
                index = 1
                loop2:             do i_hop=1, evb_max_chain
                   if ( evb_diabat_proton_log( diabat_index_donor , index ) < 0 ) exit loop2
                   evb_diabat_proton_log(diabat_index_acceptor,index) = evb_diabat_proton_log(diabat_index_donor,index)
                   evb_diabat_proton_log(diabat_index_acceptor,index+1) = evb_diabat_proton_log(diabat_index_donor,index+1)
                   index = index + 2
                enddo loop2
                ! now add new proton transfer         
                evb_diabat_proton_log(diabat_index_acceptor,index) = i_atom
                evb_diabat_proton_log(diabat_index_acceptor,index+1) = i_mole_acceptor

!!$                write(*,*) "proton log"
!!$                write(*,*) evb_diabat_proton_log(diabat_index_acceptor,:)

                ! now construct the evb hamiltonian matrix elements and store them for this proton donor/acceptor pair.  This subroutine will update the evb_hamiltonian and evb_forces_store arrays.
                call evb_hamiltonian_elements_donor_acceptor( diabat_index_donor, diabat_index_acceptor, i_mole_principle, tot_n_mole, n_mole, n_atom, xyz, r_com, chg, mass, box, dfti_desc,dfti_desc_inv,log_file )

                ! now call this recursive subroutine, with the proton acceptor molecule as the new donor, unless this is the end of a cyclic transfer
                if ( flag_cycle < 1 ) then
                call evb_conduct_proton_transfer_recursive( i_mole_acceptor, diabat_index_acceptor,  i_mole_principle, tot_n_mole, n_mole, n_atom, xyz, r_com, chg, mass, box, dfti_desc,dfti_desc_inv,log_file )
                end if

             end do loop1 ! end loop over first neighbors

          end if ! end if on reactive atom for hydrogen atom
       end do  ! end loop over hydrogen atoms

    end if ! end if on less than max hops

  end subroutine evb_conduct_proton_transfer_recursive






  !***************************************
  ! this subroutine constructs the evb hamiltonian elements
  ! for a particular donor and acceptor pair
  ! this involves creating the necessary data structures,
  ! storing the necessary output, and restoring the
  ! modified data structures
  !****************************************
  subroutine evb_hamiltonian_elements_donor_acceptor( diabat_index_donor, diabat_index_acceptor, i_mole_principle, tot_n_mole, n_mole, n_atom, xyz, r_com, chg, mass, box, dfti_desc,dfti_desc_inv,log_file )
    use MKL_DFTI
    use total_energy_forces
    integer, intent(in) :: diabat_index_donor, diabat_index_acceptor, i_mole_principle, tot_n_mole, n_mole
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: box,r_com, chg, mass
    real*8, intent(in), dimension(:,:,:) :: xyz
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv
    integer,intent(in)::log_file

    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: force_atoms
    integer :: diabat_index
    real*8 :: potential, E_elec, E_bh , E_bond, E_angle, E_dihedral 

    ! temporary diabat data structures
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: xyz_temp1
    real*8, dimension(MAX_N_MOLE,3) :: r_com_temp1
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM) :: chg_temp1, mass_temp1
    integer, dimension(MAX_N_MOLE) :: n_atom_temp1

    !****** these are junk variables to pass to calculate_total_force_energy subroutine *********
    !****** as we don't implement ms-evb for polarizable force field
    real*8 :: E_elec_nopol, E_3body
    real*8,dimension(5) :: energy_components
    integer :: iteration
    integer, dimension(MAX_N_MOLE) :: n_atom_drude
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: xyz_drude

    n_atom_drude = n_atom
    !********************************************************************************************


    !********** if calculating the energy of the primary diabat
    if ( ( diabat_index_donor == 1 ) .and. ( diabat_index_acceptor == 1 ) ) then

       diabat_index = diabat_index_donor

       ! save global data structures, 1 indicates initialize
       call evb_restore_global_data_structures(1)

!!$       write(*,*) "diabat", diabat_index

       !*********************** energy and force of primary diabat ******************************!
       call calculate_total_force_energy( force_atoms,potential, E_elec,E_elec_nopol,E_bh, E_3body, E_bond, E_angle, E_dihedral, energy_components, iteration, tot_n_mole, n_mole, n_atom, n_atom_drude, r_com, xyz, chg, box, dfti_desc,dfti_desc_inv,log_file,xyz_drude)
       ! store the forces, the "1" is to initialize
       call evb_store_forces( i_mole_principle, diabat_index, diabat_index, tot_n_mole, n_atom, force_atoms, 1 )
       ! store the energy
       evb_hamiltonian( diabat_index , diabat_index ) = potential


    else
       ! *********** calculating proton acceptor diabat and coupling to donor

       ! create new data structures for this diabat
       ! WARNING, this subroutine modifies global data structures:: See comments at top of subroutine
       call evb_create_diabat_data_structures( xyz, xyz_temp1, r_com, r_com_temp1, chg, chg_temp1 , n_atom, n_atom_temp1, n_atom_drude, mass, mass_temp1, box, diabat_index_acceptor, i_mole_principle )


!!$       write(*,*) "diabat", diabat_index_acceptor

       ! energy and force of the proton acceptor diabat
       call calculate_total_force_energy( force_atoms,potential, E_elec,E_elec_nopol,E_bh, E_3body, E_bond, E_angle, E_dihedral, energy_components, iteration, tot_n_mole, n_mole, n_atom_temp1, n_atom_drude, r_com_temp1, xyz_temp1, chg_temp1, box, dfti_desc,dfti_desc_inv,log_file,xyz_drude)
       ! store the forces
       call evb_store_forces( i_mole_principle, diabat_index_acceptor, diabat_index_acceptor, tot_n_mole, n_atom_temp1, force_atoms )
       ! store the energy
       evb_hamiltonian( diabat_index_acceptor , diabat_index_acceptor ) = potential

       ! calculate off-diagonal diabatic coupling, energy and force
       call evb_diabatic_coupling( force_atoms, potential, diabat_index_acceptor, i_mole_principle, tot_n_mole, n_atom_temp1, r_com_temp1, xyz_temp1, chg_temp1, box)

       ! store the forces, note this coupling is between the current diabat, and the diabat_index_1neighbor (principle) diabat
       call evb_store_forces( i_mole_principle, diabat_index_donor, diabat_index_acceptor, tot_n_mole, n_atom_temp1, force_atoms )
       ! store the energy
       evb_hamiltonian( diabat_index_donor , diabat_index_acceptor ) = potential

       ! restore global atom_index data structure
       call evb_restore_global_data_structures  


    end if


  end subroutine evb_hamiltonian_elements_donor_acceptor




  !************************************
  ! this subroutine searches for evb reactive neighbors
  ! of reactive atom i_mole, i_atom, and puts them in neighbor list
  !
  ! output variable flag_cycle is given a value of "1" if we transfer a proton
  ! to a hydronium oxygen (in a cyclic transfer).  This is to prevent a new cycle of the original
  ! cyclic proton transfer
  !************************************
  subroutine find_evb_reactive_neighbors( i_mole, i_atom, neighbor_list, tot_n_mole, n_mole, n_atom, xyz, r_com, box ) 
    integer, intent(in) :: i_mole, i_atom,tot_n_mole, n_mole
    integer, dimension(:), intent(out) :: neighbor_list
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: box,r_com
    real*8, intent(in), dimension(:,:,:) :: xyz

    integer :: j_mole, j_atom, index, atom_id1, atom_id2
    real*8, dimension(3) :: r_com_i, r_com_j, dr_com, rij, shift

    atom_id1 = atom_index(i_mole,i_atom)
    index=1
    ! initialize list to negative value to signal end of neighbors
    neighbor_list=-1

    do j_mole = 1, tot_n_mole
       if ( i_mole /= j_mole ) then
          ! center of mass distance
          r_com_i(:) = r_com(i_mole,:)
          r_com_j(:) = r_com(j_mole,:)

          shift = pbc_shift( r_com_i, r_com_j, box, xyz_to_box_transform)
          dr_com = pbc_dr( r_com_i, r_com_j, shift)

          if ( dot_product( dr_com, dr_com ) < evb_first_solvation_cutoff **2 ) then
             ! loop over all atoms of principle molecule and target molecule
             ! see if atom pair is a potential reactivity pair
             loop2:          do j_atom = 1 , n_atom(j_mole)

                atom_id2 = atom_index(j_mole,j_atom)
                ! if reactive pair, evb_reactivity_matrix(atom_id1, atomid2) = 1, otherwise 0
                ! IMPORTANT-note that evb_reactivity_matrix is not symmetric, and therefore the order
                ! of atoms is important.  For proton transfer, we want to use the atom order
                ! atom(H3O+), atom(H2O), since the H,O matrix element should be 1, but the O, H should
                ! be zero


                if ( evb_reactivity_matrix( atom_id1, atom_id2 ) == 1 ) then
                   ! see if this reactivity pair is within the cutoff distance
                   rij = pbc_dr( xyz(i_mole,i_atom,:), xyz(j_mole,j_atom,:), shift ) ! Note for COM cutoff, shift values unchanged
                   if ( dot_product(rij,rij) < evb_reactive_pair_distance**2 ) then
                      ! this is a reactive pair
                      neighbor_list(index)=j_mole
                      index=index+1
                      exit loop2
                   endif
                endif
             enddo loop2
          end if
       end if
    end do


  end subroutine find_evb_reactive_neighbors







  !***********************************************
  ! WARNING:  THIS SUBROUTINE MODIFIES GLOBAL VARIABLES! 
  ! Global array atom_index is modified in order to allow for energy calculation of diabat
  ! This change can be restored by calling evb_restore_global_data_structures
  !
  ! the new data structures are created based on the proton hopping information contained in evb_diabat_proton_log
  !***********************************************
  subroutine evb_create_diabat_data_structures( xyz, xyz_temp, r_com, r_com_temp, chg, chg_temp , n_atom, n_atom_temp, n_atom_drude, mass, mass_temp, box, diabat_index, i_mole_principle )
    real*8, intent(in), dimension(:,:,:) :: xyz
    real*8, intent(in), dimension(:,:) :: r_com, chg, mass, box
    integer, intent(in), dimension(:) :: n_atom
    integer, intent(in) :: i_mole_principle, diabat_index
    real*8, intent(out), dimension(:,:,:) :: xyz_temp
    real*8, intent(out), dimension(:,:) :: r_com_temp, chg_temp, mass_temp
    integer, intent(out), dimension(:) :: n_atom_temp, n_atom_drude

    integer,save :: initialize=0
    integer :: i_hop, i_mole_donor, i_atom_donor, i_mole_acceptor, i_atom_acceptor, index, i_atom, molecule_index_temp

    if ( initialize == 0 ) then
       ! initialize atom_index mapping
       call initialize_atom_index_mapping
       initialize=1
    end if

    xyz_temp = xyz
    r_com_temp = r_com
    chg_temp = chg
    n_atom_temp = n_atom
    mass_temp = mass

!!$    write(*,*) "diabat ",diabat_index
!!$    write(*,*) evb_diabat_proton_log( diabat_index, :)

    ! first donor is always principle hydronium
    ! set i_mole_acceptor equal to this, because
    ! i_mole_donor will be copied from previous i_mole_acceptor
    i_mole_acceptor = i_mole_principle
    index = 1
    ! loop through all proton hops.  A negative integer signals the end of the chain
    loop1: do i_hop=1, evb_max_chain
        if ( evb_diabat_proton_log( diabat_index , index ) < 0 ) exit loop1
       i_mole_donor = i_mole_acceptor
       i_atom_donor = evb_diabat_proton_log( diabat_index , index )
       i_mole_acceptor = evb_diabat_proton_log( diabat_index , index + 1 )
       
       ! need to change this if more than one hydronium
       if ( n_hydronium_molecules > 1 ) stop "see code below"
       hydronium_molecule_index(1) = i_mole_acceptor

       !******************** first change topology information
       i_atom_acceptor = n_atom_temp(i_mole_acceptor) + 1
       xyz_temp(i_mole_acceptor,i_atom_acceptor,:) = xyz_temp(i_mole_donor,i_atom_donor,:)
       mass_temp(i_mole_acceptor,i_atom_acceptor) = mass_temp(i_mole_donor,i_atom_donor)
       n_atom_temp(i_mole_acceptor) = n_atom_temp(i_mole_acceptor) + 1

       ! here it's possible that the transfered proton could be split from the rest of the
       ! molecule by a periodic image. We need to fix this, as energy routines(bonds, angles, dihedrals)
       ! assume the molecule to not be split over pbc conditions
       call check_intra_molecular_shifts(i_mole_acceptor, n_atom_temp, xyz_temp, box)
       r_com_temp(i_mole_acceptor,:) = pos_com( xyz_temp, i_mole_acceptor, n_atom_temp , mass_temp )


       !**********************  now change force field information
       atom_index(i_mole_acceptor,i_atom_acceptor) = atom_index(i_mole_donor,i_atom_donor)
       ! don't overwrite the atom index we just filled in, hence -1
       do i_atom = 1 , n_atom_temp(i_mole_acceptor) - 1
          atom_index(i_mole_acceptor, i_atom) = map_atom_index(atom_index(i_mole_acceptor,i_atom))
       enddo

       do i_atom = 1 , n_atom_temp(i_mole_donor)
          atom_index(i_mole_donor, i_atom) = map_atom_index(atom_index(i_mole_donor,i_atom))
       enddo

       do i_atom = 1, n_atom_temp(i_mole_acceptor)
          chg_temp(i_mole_acceptor,i_atom) = atype_chg( atom_index(i_mole_acceptor, i_atom) )
       enddo

       do i_atom = 1, n_atom_temp(i_mole_donor)
          chg_temp(i_mole_donor,i_atom) = atype_chg( atom_index(i_mole_donor, i_atom) )
       enddo

       ! now we have to delete the transferred proton data from i_mole_donor, and shift the contents
       ! of the data structures for that molecule.  If the transferred proton was the last atom
       ! in the molecule's data structure, then we don't have to do anything because that data
       ! will not get looped over, since the n_atom array was decremented
       do i_atom= i_atom_donor + 1 , n_atom_temp(i_mole_donor)
          xyz_temp(i_mole_donor, i_atom - 1, :) = xyz_temp(i_mole_donor, i_atom , :)
          mass_temp(i_mole_donor, i_atom - 1) = mass_temp(i_mole_donor, i_atom )
          chg_temp(i_mole_donor, i_atom - 1) = chg_temp(i_mole_donor, i_atom )
          atom_index(i_mole_donor, i_atom - 1) = atom_index(i_mole_donor, i_atom )
       enddo

       ! now decrement n_atom_temp array for donor
       n_atom_temp(i_mole_donor) = n_atom_temp(i_mole_donor) - 1

       ! center of mass of donor
       r_com_temp(i_mole_donor,:) = pos_com( xyz_temp, i_mole_donor, n_atom_temp , mass_temp )

       ! molecule index
       molecule_index_temp = molecule_index(i_mole_acceptor)
       molecule_index(i_mole_acceptor) = molecule_index(i_mole_donor)
       molecule_index(i_mole_donor) = molecule_index_temp

       index = index + 2
    enddo loop1

    ! n_atom_drude should just be a copy here (not using drude oscillators)
    n_atom_drude = n_atom_temp

!!$    write(*,*) "test data structure"
!!$    write(*,*) i_mole_donor, i_mole_acceptor
!!$    do i_mole=1,4
!!$       do i_atom=1,n_atom_temp(i_mole)
!!$          write(*,*) i_mole, i_atom
!!$          write(*,*) atom_index(i_mole,i_atom)
!!$          write(*,*) i_mole, i_atom, xyz_temp(i_mole,i_atom,:)
!!$          write(*,*) chg_temp(i_mole,i_atom), mass_temp(i_mole,i_atom)
!!$       enddo
!!$    enddo
!!$    write(*,*) "end test data structure"
!!$    stop
          

  end subroutine evb_create_diabat_data_structures




  !*************************************************
  ! this subroutine calculates the diabatic coupling 
  ! hamiltonian matrix elements and their corresponding forces
  !
  ! the formula for these matrix elements is given by equations
  ! 10-12 of JPC B, 2008, 112, 467-482
  !
  ! we note the sum in equation 11 implies that no ewald sum was used to
  ! calculate the solvent dependent electrostatics
  ! We therefore do not use an ewald sum, and also do not use a cutoff
  ! to calculate this interaction.  Therefore this term should 
  ! be order(N) to evaluate
  !************************************************
  subroutine evb_diabatic_coupling( force_atoms, potential, diabat_index, i_mole_principle, tot_n_mole, n_atom, r_com, xyz, chg, box )
    real*8, intent(out), dimension(:,:,:) :: force_atoms
    real*8, intent(in), dimension(:,:,:) :: xyz
    real*8, intent(out)                 :: potential
    real*8, intent(in), dimension(:,:) :: r_com, chg, box
    integer, intent(in), dimension(:) :: n_atom
    integer, intent(in) :: i_mole_principle, diabat_index, tot_n_mole
    integer ::  i_hop, i_mole_donor, i_atom_donor, i_mole_acceptor, index, i_atom
    real*8 ::  Vex , A
    real*8, dimension(:,:,:), allocatable :: dVex
    real*8, dimension(3,3) :: dA

    ! test
    real*8,dimension(MAX_N_MOLE, MAX_N_ATOM,3) :: dVtest, xyz_test
    real*8,dimension(3) :: dr, force
    real*8,dimension(3,3) :: dAjunk
    real*8 :: A_new, A_dnew, V_new, V_dnew
    integer :: i_mole
    ! test

    potential=0d0
    force_atoms=0d0

    allocate( dVex( size(force_atoms(:,1,1)) , size(force_atoms(1,:,1)) , size(force_atoms(1,1,:)) ) )

    ! figure out proton donor and proton acceptor from evb_diabat_proton_log.
    ! we want the final donor and acceptor, which is the last proton transfer

    ! first donor is always principle hydronium
    ! set i_mole_acceptor equal to this, because
    ! i_mole_donor will be copied from previous i_mole_acceptor
    i_mole_acceptor = i_mole_principle
    index = 1
    ! loop through all proton hops.  A negative integer signals the end of the chain
    do i_hop=1, evb_max_chain
       if ( evb_diabat_proton_log( diabat_index , index ) < 0 ) exit
       i_mole_donor = i_mole_acceptor
       i_atom_donor = evb_diabat_proton_log( diabat_index , index )
       i_mole_acceptor = evb_diabat_proton_log( diabat_index , index + 1 )
       index = index + 2
    enddo


    ! first calculate solvent dependent prefactor, and its derivate w.r.t. atomic coordinates
    call evb_diabatic_coupling_electrostatics( Vex, dVex, i_mole_donor, i_mole_acceptor, tot_n_mole, n_atom, r_com, xyz , box, chg )


    ! now calculate geometry dependent scale factor and its derivatives
    ! three derivatives are contained in the output dA array, in order:
    ! O donor, O acceptor, H central Zundel Hydrogen
    call evb_diabatic_coupling_geometric( A , dA, i_mole_donor, i_mole_acceptor, n_atom, xyz , box )


!!$    write(*,*) "test Vex derivative"
!!$    do i_mole=1,tot_n_mole
!!$       do i_atom = 1, n_atom(i_mole)
!!$          xyz_test=xyz
!!$          dr=(/ 0.01 , 0.01, 0.01 /)
!!$          xyz_test(i_mole, i_atom,:) = xyz(i_mole,i_atom,:) + dr(:)
!!$
!!$          force(:) = dVex(i_mole,i_atom,:)
!!$          V_dnew = Vex + dot_product( force, dr )
!!$          call evb_diabatic_coupling_electrostatics( V_new, dVtest, i_mole_donor, i_mole_acceptor, tot_n_mole, n_atom, r_com, xyz_test , box, chg )
!!$          write(*,*) i_mole, i_atom, Vex, V_new, V_dnew
!!$       enddo
!!$    enddo
!!$    stop
    ! test


    ! now form evb diabatic coupling matrix element
    potential = ( Vconstij_evb + Vex ) * A


    ! now forces for this matrix element, first forces from geometric derivative
    ! O donor
    call match_first_atom_type( i_atom, i_mole_donor, n_atom, water_oxygen_label )
    force_atoms(i_mole_donor, i_atom,:) = - ( Vconstij_evb + Vex ) * dA(1,:)
    ! O acceptor
    call match_first_atom_type( i_atom, i_mole_acceptor, n_atom, hydronium_oxygen_label )
    force_atoms(i_mole_acceptor, i_atom,:) = - ( Vconstij_evb + Vex ) * dA(2,:)

    ! H central Zundel 
    force_atoms(i_mole_acceptor, n_atom(i_mole_acceptor),:) = - ( Vconstij_evb + Vex ) * dA(3,:)

    ! now forces from electrostatic derivative
    force_atoms(:,:,:) = force_atoms(:,:,:) - dVex(:,:,:) * A

    deallocate( dVex )

  end subroutine evb_diabatic_coupling



  !********************************
  ! this returns the geometric term for the
  ! diabatic coupling element, as well as its
  ! derivatives with respect to the positions
  ! of the oxygen atoms of the donor and acceptor,
  ! and the central hydrogen atom
  ! 
  ! dA(1,:) contains derivative w.r.t  O_donor coordinates
  ! dA(2,:) contains derivative w.r.t. O_acceptor coordinates
  ! dA(3,:) contains derivative w.r.t. H zundel coordinates
  !********************************
  subroutine evb_diabatic_coupling_geometric( A , dA, i_mole_donor, i_mole_acceptor, n_atom, xyz , box )
    real*8, intent(out) :: A
    real*8, dimension(:,:), intent(out) :: dA
    integer, intent(in)    :: i_mole_donor, i_mole_acceptor
    real*8, intent(in), dimension(:,:,:) :: xyz
    real*8, intent(in), dimension(:,:) :: box
    integer, intent(in), dimension(:) :: n_atom

    real*8 , dimension(3) :: r_O1, r_O2, r_H, r_OO , q, shift, r_ij
    real*8                :: r_OO_mag, q2_mag, q_mag, fac1, fac2, fac3 , dfac1, dfac2, dfac3
    integer :: i_atom


    ! first get distances, need r_OO and q, q is defined by
    ! q = ( rO1 + rO2 ) / 2 - rH
    ! to consider pbc, shift all atoms relative to the donor oxygen

    ! get O1, note we're using water_oxygen_label for the donor, because we've already transferred the proton
    call match_first_atom_type( i_atom, i_mole_donor, n_atom, water_oxygen_label )
    r_O1(:) = xyz(i_mole_donor, i_atom, :)

    ! get O2, see above comment
    call match_first_atom_type( i_atom, i_mole_acceptor, n_atom, hydronium_oxygen_label )
    r_O2(:) = xyz(i_mole_acceptor, i_atom, :)
    shift = pbc_shift( r_O1 , r_O2 , box , xyz_to_box_transform )
    r_ij = pbc_dr( r_O1 , r_O2 , shift )
    r_O2(:) = r_O1(:) + r_ij(:)

    ! H atom is last atom in acceptor
    r_H(:) = xyz(i_mole_acceptor, n_atom(i_mole_acceptor) , : )
    r_ij = pbc_dr( r_O1 , r_H , shift )
    r_H(:) = r_O1(:) + r_ij(:)

    r_OO = r_O1 - r_O2
    q = ( r_O1 + r_O2 ) / 2d0 - r_H

    r_OO_mag = dsqrt( dot_product( r_OO , r_OO ) )
    q2_mag   = dot_product( q , q )
    q_mag = dsqrt ( q2_mag )

    ! geometric factor is given as a product of three terms
    fac1 = exp( -gamma_evb * q2_mag )
    fac2 = 1d0 + Pu_evb * exp( -k_evb * ( r_OO_mag - DuOO_evb ) ** 2 )
    fac3 = 0.5d0 * ( 1d0 - tanh( beta_evb * ( r_OO_mag - RuOO0_evb ) ) ) + Pup_evb * exp( -alpha_evb * ( r_OO_mag - rlOO0_evb ) )

    ! all derivatives are relative to the general distance coordinate ( either q, r_OO )
    dfac1 = -gamma_evb * 2d0 * q_mag *  exp( -gamma_evb * q2_mag )
    dfac2 = Pu_evb * -k_evb * 2d0 * ( r_OO_mag - DuOO_evb ) * exp( -k_evb * ( r_OO_mag - DuOO_evb ) ** 2 )
    dfac3 = -0.5d0 * beta_evb / cosh( beta_evb * ( r_OO_mag - RuOO0_evb ) )**2 - Pup_evb * alpha_evb * exp( -alpha_evb * ( r_OO_mag - rlOO0_evb ) ) 


    ! geometric factor
    A = fac1 * fac2 * fac3

    ! derivative w.r.t  O_donor coordinates
    dA(1,:) = dfac1 * fac2 * fac3 * 0.5d0 * q(:) / q_mag
    dA(1,:) = dA(1,:) + fac1 * dfac2 * fac3 * r_OO(:) / r_OO_mag
    dA(1,:) = dA(1,:) + fac1 * fac2 * dfac3 * r_OO(:) / r_OO_mag

    ! derivative w.r.t. O_acceptor coordinates
    dA(2,:) = dfac1 * fac2 * fac3 * 0.5d0 * q(:) / q_mag
    dA(2,:) = dA(2,:) + fac1 * dfac2 * fac3 * -r_OO(:) / r_OO_mag
    dA(2,:) = dA(2,:) + fac1 * fac2 * dfac3 * -r_OO(:) / r_OO_mag  

    ! derivative w.r.t. H zundel coordinates
    dA(3,:) = dfac1 * fac2 * fac3 * -q(:) / q_mag


  end subroutine evb_diabatic_coupling_geometric






 !**********************************************
 ! this is the electrostatic part of the diabatic coupling
 ! it is composed of Coulombic interactions between
 ! the solvent water molecules and the H5O2 proton
 ! transfer transition state, with special charges
 ! for the H5O2 complex
 !**********************************************
  subroutine evb_diabatic_coupling_electrostatics( Vex, dVex, i_mole_donor, i_mole_acceptor, tot_n_mole, n_atom, r_com, xyz , box, chg )
    real*8, intent(out) :: Vex
    real*8, intent(out), dimension(:,:,:) :: dVex
    integer, intent(in) :: i_mole_donor, i_mole_acceptor, tot_n_mole
    real*8, intent(in), dimension(:,:,:) :: xyz
    real*8, intent(in), dimension(:,:) :: r_com, chg, box
    integer, intent(in), dimension(:) :: n_atom

    integer :: j_mole, i_atom, j_atom,  i_type
    real*8 :: q_i , q_j, r_mag 
    real*8, dimension(3) :: shift, r_ij, dV_ij

    Vex=0d0
    dVex=0d0

    ! this is a sum over coulomb interactions of all water molecules with the 7 atoms of the 
    ! H5O2+ Zundel complex


    ! first calculate interaction with donor water
    do i_atom=1, n_atom(i_mole_donor)
       ! figure out charge
       i_type = atom_index(i_mole_donor, i_atom )

       ! note we're using water_oxygen_label for the donor, because we've already transferred the proton
       if( atype_name( i_type ) .eq. water_proton_label ) then
          q_i = qHex_evb
       else if( atype_name( i_type ) .eq. water_oxygen_label ) then
          q_i = qOex_evb
       else
          write(*,*) ""
          write(*,*) " can't recognize atom type in evb_diabatic_coupling subroutine "
          write(*,*) ""
          stop
       end if

       do j_mole = 1, tot_n_mole  ! loop over solvent water molecules
          if ( ( j_mole /= i_mole_donor ) .and. ( j_mole /= i_mole_acceptor ) ) then
             shift = pbc_shift( r_com(i_mole_donor,:) , r_com(j_mole,:) , box , xyz_to_box_transform )
             do j_atom =1 , n_atom(j_mole)

                q_j = chg( j_mole, j_atom )
                r_ij = -pbc_dr( xyz(i_mole_donor, i_atom,:) , xyz(j_mole, j_atom, :) , shift )
                r_mag = dsqrt( dot_product( r_ij , r_ij ) )

                Vex = Vex + q_i * q_j / r_mag

                dV_ij = - q_i * q_j / r_mag ** 3 * r_ij
                dVex(i_mole_donor,i_atom,:) = dVex(i_mole_donor,i_atom,:) + dV_ij
                dVex(j_mole, j_atom, :) = dVex(j_mole, j_atom, :) - dV_ij
             enddo
          endif
       enddo
    enddo



    ! now calculate interaction with acceptor water
    do i_atom=1, n_atom(i_mole_acceptor)
       ! figure out charge, last atom in acceptor is central Zundel hydrogen which has a unique charge
       i_type = atom_index(i_mole_acceptor, i_atom )
       ! note we're using hydronium proton label for acceptor, because we've already transferred the proton
       if( atype_name( i_type ) .eq. hydronium_proton_label ) then
          if ( i_atom == n_atom(i_mole_acceptor) ) then
             ! central Zundel hydrogen
             q_i = qHexs_evb
          else
             q_i = qHex_evb
          end if
       else if( atype_name( i_type ) .eq. hydronium_oxygen_label ) then
          q_i = qOex_evb
       else
          write(*,*) ""
          write(*,*) " can't recognize atom type in evb_diabatic_coupling subroutine "
          write(*,*) ""
          stop
       end if


       do j_mole = 1, tot_n_mole  ! loop over solvent water molecules
          if ( ( j_mole /= i_mole_donor ) .and. ( j_mole /= i_mole_acceptor ) ) then
             shift = pbc_shift( r_com(i_mole_acceptor,:) , r_com(j_mole,:) , box , xyz_to_box_transform )
             do j_atom =1 , n_atom(j_mole)

                q_j = chg( j_mole, j_atom )
                r_ij = -pbc_dr( xyz(i_mole_acceptor, i_atom,:) , xyz(j_mole, j_atom, :) , shift )
                r_mag = dsqrt( dot_product( r_ij , r_ij ) )

                Vex = Vex + q_i * q_j / r_mag

                dV_ij = - q_i * q_j / r_mag ** 3 * r_ij
                dVex(i_mole_acceptor,i_atom,:) = dVex(i_mole_acceptor,i_atom,:) + dV_ij
                dVex(j_mole, j_atom, :) = dVex(j_mole, j_atom, :) - dV_ij
             enddo
          endif
       enddo
    enddo

    ! convert from e^2/A to kJ/mol
    Vex = Vex * 0.52914D0 * 627.51D0 * 4.184D0 
    dVex = dVex * 0.52914D0 * 627.51D0 * 4.184D0 


  end subroutine evb_diabatic_coupling_electrostatics



  !******************************
  ! WARNING:  THIS SUBROUTINE MODIFIES GLOBAL VARIABLES! 
  ! this subroutine saves the global data structure
  ! atom_index on its first call, and restores
  ! this data structure to the saved data on 
  ! subsequent calls
  !******************************
  subroutine evb_restore_global_data_structures( initialize )
    integer, optional :: initialize
    integer, dimension(MAX_N_MOLE,MAX_N_ATOM), save :: atom_index_store 
    integer, dimension(MAX_N_MOLE), save :: molecule_index_store
    integer, dimension( n_proton_max ), save :: hydronium_molecule_index_store

    if ( present(initialize) ) then
       atom_index_store = atom_index
       molecule_index_store = molecule_index
       hydronium_molecule_index_store = hydronium_molecule_index
    end if

    atom_index = atom_index_store
    molecule_index = molecule_index_store
    hydronium_molecule_index = hydronium_molecule_index_store

  end subroutine evb_restore_global_data_structures




!!$  subroutine ms_evb_diabat_electrostatic_energy( tot_n_mole, n_mole, n_atom, xyz, chg, box, r_com,dfti_desc,dfti_desc_inv,Q_change)
!!$
!!$    use MKL_DFTI
!!$    use omp_lib
!!$    integer, intent(in) :: tot_n_mole, n_mole,Q_change
!!$    integer, intent(in), dimension(:) :: n_atom
!!$    real*8, intent(in), dimension(:,:) :: chg, box,r_com
!!$    real*8, intent(in), dimension(:,:,:) :: xyz
!!$    TYPE(DFTI_DESCRIPTOR),pointer,intent(in):: dfti_desc,dfti_desc_inv
!!$    integer :: i_mole, i_atom, j_mole, j_atom
!!$
!!$
!!$  end subroutine ms_evb_diabat_electrostatic_energy




  !*********************************
  ! this subroutine stores the forces for the evb-hamitonian matrix elements
  ! by constructing a lookup table for the non-zero elements
  !
  ! note that each time this subroutine is called for a different diabat, the 
  ! force_atoms array will correspond to a different xyz topology (i.e. the protons
  ! will correspond to different molecules)
  !
  ! therefore, we need to map the force_atoms array back to the principle
  ! xyz topology of the principle diabat, and store all the forces
  ! in this same, consistent topology
  !********************************* 
  subroutine evb_store_forces( i_mole_principle, diabat_index1, diabat_index2, tot_n_mole, n_atom, force_atoms, initialize )
    integer, intent(in) :: i_mole_principle, diabat_index1, diabat_index2, tot_n_mole
    integer, dimension(:), intent(in) :: n_atom
    real*8,dimension(:,:,:), intent(in) :: force_atoms
    integer, intent(in),optional :: initialize

    integer, save :: store_index=1, size
    real*8, dimension(:,:,:), allocatable :: force_atoms_map
    integer :: i_mole, i_atom, i_hop, i_mole_donor, i_atom_donor, i_mole_acceptor, i_atom_acceptor, index, diabat_index

    ! initialize, allocate evb_force_store array.  Here we guess as to how large an array we need to allocate, which is based on the number of diabats plus the number of non-zero diabat coupling matrix elements.  The couplings are only non-zero for diabats that have hydronium on adjacent water molecules, so this corresponds to the number of donor hydrogen bonds per water molecule.  The number of water molecules is strictly less than the number of diabats. so 2 * n_diabats is a safe guess
    if ( present(initialize) ) then
       size = 3 * evb_max_states
       ! we call this subroutine with initialize present to reset the indices, so this array may already be allocated
       if ( .not. allocated(evb_forces_store) ) then
          allocate( evb_forces_store(size,tot_n_mole,MAX_N_ATOM,3) )
          evb_forces_store=0d0
       end if
       ! initialize evb_forces_lookup_index to negative value, so that we know which couplings have non-zero force
       evb_forces_lookup_index=-1
       ! reset this
       store_index=1
    end if

    ! store mapping
    evb_forces_lookup_index( diabat_index1 , diabat_index2 ) = store_index



    allocate( force_atoms_map(tot_n_mole,MAX_N_ATOM,3) )
    ! copy forces.  Note force_atoms array is typically of size (MAX_N_MOLE, MAX_N_ATOM,3) , with 
    ! MAX_N_MOLE > tot_n_mole, so we loop to copy elements
    do i_mole=1, tot_n_mole
       force_atoms_map( i_mole, : , : ) = force_atoms(i_mole, : , : )
    enddo

    if ( diabat_index2 > diabat_index1 ) then
       diabat_index = diabat_index2
    else
       diabat_index = diabat_index1
    end if

    ! now map force_atoms back to principle topology to construct force_atoms_map array
    if ( diabat_index > 1 ) then
       ! first donor is always principle hydronium
       ! set i_mole_acceptor equal to this, because
       ! i_mole_donor will be copied from previous i_mole_acceptor
       i_mole_acceptor = i_mole_principle
       index = 1
       ! loop through all proton hops.  A negative integer signals the end of the chain
       do i_hop=1, evb_max_chain
          if ( evb_diabat_proton_log( diabat_index , index ) < 0 ) exit
          i_mole_donor = i_mole_acceptor
          i_atom_donor = evb_diabat_proton_log( diabat_index , index )
          i_mole_acceptor = evb_diabat_proton_log( diabat_index , index + 1 )

          ! fix the forces for this donor.  Note that the donated proton is always the last
          ! atom of the acceptor molecule
          do i_atom = i_atom_donor, n_atom(i_mole_donor)
             ! shift forces up
             force_atoms_map( i_mole_donor, i_atom + 1 , : ) = force_atoms( i_mole_donor, i_atom  , : )
          enddo
          ! now copy force from acceptor molecule back to donor
          force_atoms_map( i_mole_donor, i_atom_donor , : ) = force_atoms( i_mole_acceptor , n_atom(i_mole_acceptor) , : )

          index = index + 2
          i_mole_donor = i_mole_acceptor
       end do
    end if

    evb_forces_store( store_index, : , : , : ) = force_atoms_map(: , : , :)

    store_index = store_index + 1

    ! check to make sure we're not overfilling array
    if ( store_index > size ) then
       write(*,*) " "
       write(*,*) " too many non-zero ms-evb matrix elements for allocated data structures "
       write(*,*) " please increase the size of evb_max_states parameter                   "
       write(*,*) ""
       stop
    endif

    deallocate( force_atoms_map )


  end subroutine evb_store_forces


  !****************************
  ! this subroutine initializes the mapping
  ! between water and hydronium ion atom types
  ! for use after a proton hop
  !****************************
  subroutine initialize_atom_index_mapping
    integer :: i_type, flag, H_hydronium_index, H_water_index, O_hydronium_index, O_water_index

    flag=0
    do i_type = 1, n_atom_type
       if ( atype_name(i_type ) .eq. hydronium_proton_label ) then
          H_hydronium_index=i_type
          flag = flag + 1
       else if( atype_name(i_type ) .eq. hydronium_oxygen_label ) then
          O_hydronium_index=i_type
          flag = flag + 1
       else if( atype_name(i_type ) .eq. water_proton_label ) then
          H_water_index=i_type
          flag = flag + 1
       else if( atype_name(i_type ) .eq. water_oxygen_label ) then
          O_water_index=i_type
          flag = flag + 1
       end if
    enddo

    ! make sure we found all atom types
    if ( flag /= 4 ) then
       write(*,*) ""
       write(*,*) "error finding hydronium and water atom types in "
       write(*,*) "initialize_atom_index_mapping of ms_evb module "
       write(*,*) ""
       stop
    end if

    ! now construct mapping

    map_atom_index(H_hydronium_index) = H_water_index
    map_atom_index(H_water_index) = H_hydronium_index
    map_atom_index(O_hydronium_index) = O_water_index
    map_atom_index(O_water_index) = O_hydronium_index

  end subroutine initialize_atom_index_mapping


  !*****************************
  ! this function checks whether the atom
  ! type is reactive, i.e can be transferred
  ! to another molecule
  !*****************************
  function check_atom_reactivity( atom_id )
    integer :: check_atom_reactivity
    integer, intent(in) :: atom_id

    ! all protons are reactive
    if ( ( atype_name(atom_id) .eq. hydronium_proton_label  ) .or. ( atype_name(atom_id) .eq. water_proton_label  ) )  then
       check_atom_reactivity=1
    else
       check_atom_reactivity=0
    end if

  end function check_atom_reactivity


  !*************************
  ! this function checks if the input molecule is 
  ! a hydronium molecule, using the global array
  ! hydronium_molecule_index
  ! returns 1 for hydronium molecule
  !*************************
  integer function check_hydronium_molecule(i_mole )
    integer, intent(in) :: i_mole
    integer :: index, j_mole
    check_hydronium_molecule=-1

    do j_mole=1,size(hydronium_molecule_index)
       if ( hydronium_molecule_index(j_mole) == i_mole ) then
          check_hydronium_molecule=1
       endif
    enddo

  end function check_hydronium_molecule



  !******************************
  ! this subroutine constructs the evb_reactivity_matrix
  ! which determines which atom pairs are reactive, i.e.
  ! where proton transfers can occur
  !
  ! note, this matrix will not be symmetric
  !
  ! see comments in subroutine 
  ! evb_conduct_proton_transfer_recursive
  ! concerning which atom types should be considered reactive
  ! we need to consider hydronium and water protons transfered
  ! to water oxygen, as well as water protons transfered to
  ! hydronium oxygen
  !******************************
  subroutine construct_evb_reactivity_matrix
    use global_variables

    integer :: i_atom, j_atom

    evb_reactivity_matrix=0

    do i_atom = 1, n_atom_type
       ! protons in first index
    if ( ( atype_name(i_atom) .eq. hydronium_proton_label  ) .or. ( atype_name(i_atom) .eq. water_proton_label  ) )  then
          do j_atom = 1 , n_atom_type
             ! proton transfer to water
             if ( atype_name(j_atom) .eq. water_oxygen_label  ) then
                evb_reactivity_matrix(i_atom,j_atom)=1
             end if

             ! proton transfer back to hydronium
             if ( ( atype_name(i_atom) .eq. water_proton_label  ) .and. ( atype_name(j_atom) .eq. hydronium_oxygen_label  ) ) then
                evb_reactivity_matrix(i_atom,j_atom)=1
             end if

          enddo
       endif
    enddo

  end subroutine construct_evb_reactivity_matrix



  !******************************
  ! this subroutine makes sure we don't have too many
  ! diabats for the setting of evb_max_states
  !******************************
  subroutine evb_check_number_diabats

    ! note these are shared variables within the module
    if ( diabat_index > evb_max_states ) then
       write(*,*) ""
       write(*,*) "Found more diabat states than the current setting of"
       write(*,*) "evb_max_states= ", evb_max_states
       write(*,*) "either increase the value of this parameter, or "
       write(*,*) "decrease the number of solvation shells to which"
       write(*,*) "a proton can hop.  This is controlled by the parameter"
       write(*,*) "'evb_max_chain', which has a current value of ", evb_max_chain
       write(*,*) ""
    end if

  end subroutine evb_check_number_diabats



end module ms_evb
