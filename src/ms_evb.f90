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
  ! To efficiently construct the diabats, some PME tricks need to be used:
  ! Q_grid and dQ_dr grid should be stored for principle diabat, and efficiently
  ! updated for new diabats
  !************************************



  real*8, parameter :: evb_first_solvation_cutoff = 5d0
  real*8, parameter :: evb_reactive_pair_distance = 2.5d0
  integer, parameter :: evb_max_neighbors=10  ! this controls the dimensions of some of the evb data structures
  ! evb_max_states is the maximum dimension for evb-hamiltonian, and therefore should be set to the maximum number of diabats.  
  integer, parameter :: evb_max_states=80

  ! maximum size of water chain for proton hopping.  Note this may need to be larger than the number
  ! of solvation shells for proton transfer.  For instance in a water hexamer, for a particular initial
  ! proton configuration, may need 6 hops to generate all diabats
  integer, parameter :: evb_max_chain=3 
  real*8 , dimension(evb_max_states,evb_max_states)  :: evb_hamiltonian
  integer, dimension(evb_max_states,evb_max_states)  :: evb_forces_lookup_index
  real*8, dimension(:,:,:), allocatable :: evb_forces_store


  ! this data structure is analogous to global atom_data_type, but is meant to
  ! be locally utilized for each diabat, so we don't use pointers here
  ! we don't really need velocity, mass arrays here, but we keep them for
  ! consistency with atom_data_type structure
  type atom_data_diabat_type
   real*8, dimension(:,:), allocatable :: xyz
   real*8, dimension(:,:), allocatable :: velocity
   real*8, dimension(:,:), allocatable :: force
   real*8, dimension(:),  allocatable  :: mass
   real*8, dimension(:),  allocatable  :: charge
   integer, dimension(:), allocatable  :: atom_type_index  ! this is index of atom_type to look up force field parameters for this atom
   character(MAX_ANAME)(:),allocatable :: aname
  end type atom_data_diabat_type



  !************************************ evb topology information ***********************
  ! this stores the information about the proton hops necessary to get from the principle diabat to the specified diabat
  ! data is stored in groups of 2 integers, specifying index of hydrogen atom proton donor water molecule, and 
  ! for proton hop of each diabat, all relevant information involving transferring proton is stored.
  ! Thus, the last index of this array, contains, in order, 1) index of donor molecule, 2) index of transferring proton, 3) index of
  ! basic atom of donor molecule to which donating proton is attached, 4) index of acceptor molecule, 5) index of basic atom of acceptor molecule
  ! we store these three atom indices because, in general, evb hamiltonian elements explicitly depend on these three atoms.
  integer, dimension(evb_max_states, evb_max_chain, 5 ) :: evb_diabat_proton_log
  !*****************************************************************************************************

  ! this stores the diabat index of the state that has a non-zero coupling hamiltonian matrix element with the particular diabat
  integer, dimension(evb_max_states) :: evb_diabat_coupling_matrix

  !*************** this array saves the Q_grid for each diabat.  Note, to read these from memory efficiently, store as Q_grid_diabats(i,j,k,i_diabat)
  real*8,dimension(:,:,:,:), allocatable :: Q_grid_diabats
  integer, dimension(evb_max_states) :: Q_grid_filled

  ! this is a counter that is constantly updated with every recursive call to evb_conduct_proton_transfer_recursive subroutine
  integer :: diabat_index

contains

  !************************
  ! initialization stuff involved with MS-EVB simulation
  !***********************
  subroutine initialize_evb_simulation( n_mole, file_io_data )
    integer, intent(in) :: n_mole
    type(file_io_data_type), intent(in) :: file_io_data

    integer :: ios
    character(200) :: line

    ! open file to print evb information.  If this trajectory is restarting, this file will be old
    Select Case(restart_trajectory)
    Case("yes")
       open( file_io_data%ofile_hop_file_h, file=file_io_data%ofile_hop , status='old' )
       ! prepare this file to be appended to
       do   
          read(file_io_data%ofile_hop_file_h, '(A)', iostat=ios ) line
          if ( ios < 0 ) Exit
       enddo
    case default
       open( file_io_data%ofile_hop_file_h, file=file_io_data%ofile_hop , status='new' )
    End Select

    ! get evb parameters and topology information from .top file
    call read_evb_parameters( file_io_data%ifile_top_file_h, file_io_data%ifile_top )
    call read_evb_topology( file_io_data%ifile_top_file_h, file_io_data%ifile_top )


    ! update hydronium molecule list
    call update_hydronium_molecule_index( n_mole )

    ! check consistency
    call evb_consistency_checks( )


    ! inform about ms evb parameters
    write(*,*) ""
    write(*,*) "*********************************************"
    write(*,*) "    Settings of ms-evb simulation parameters :  "
    write(*,*) "evb_max_chain:  ", evb_max_chain
    write(*,*) "evb_reactive_pair_distance:  ", evb_reactive_pair_distance
    write(*,*) "evb_max_states:  ", evb_max_states
    write(*,*) "*********************************************"
    write(*,*) ""

  end subroutine initialize_evb_simulation



  !************************
  ! this subroutine updates the hydronium_molecule_index array
  ! with the indices of hydronium molecules
  ! this is now general, and applies to any molecule
  ! that is in its acidic state
  !************************
  subroutine update_hydronium_molecule_index( n_mole )
    integer,intent(in) :: n_mole

    integer :: i_mole, i_type, i_hydronium

    i_hydronium = 1
    do i_mole = 1 , n_mole
       ! see if this molecule type contains an acidic hydrogen 
       i_type = molecule_index(i_mole)
       if ( evb_proton_index(i_type) > 0 ) then
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
       ! in particular, need to change evb_change_diabat_data_structure_topology
       ! to deal with more than one hydronium
       !
       ! also, we consider both acids and bases as possible acceptor and donor molecules
       ! if more than one hydronium, it is potentially possible (but i think unlikely) to
       ! consider a protonated acid as an acceptor for a second hydrogen.  This will crash the
       ! code currently, and will need to be fixed.
    endif

  end subroutine update_hydronium_molecule_index



  !*************************
  ! we make use of some specific input requirements for a ms-evb simulation
  ! to allow easier use and manipulation of data structures.  This subroutine 
  ! checks for such requirements
  !*************************
  subroutine evb_consistency_checks( )
    integer              :: i_mole,i_atom,i_type, flag_h
    character(MAX_ANAME) :: atomname


    ! We will be adding and removing hydrogen atoms to and from
    ! conjugate acids and bases
    ! This requires manipulation of data structures.
    ! The simplest way to do this is to require that 
    ! acidic protons are stored last

    do i_mole=1,n_molecule_type
       if ( evb_acid_molecule( i_mole ) == 1 ) then
          flag_h = 0
          loop2:  do i_atom=1,MAX_N_ATOM_TYPE
             if ( molecule_type(i_mole,i_atom) == MAX_N_ATOM_TYPE + 1 ) exit loop2
             if ( evb_reactive_protons( i_mole , i_atom ) == 1 ) flag_h=1
             if ( ( flag_h == 1 ) .and. ( evb_reactive_protons( i_mole , i_atom ) == 0 ) ) then
                write(*,*) "acidic protons must be defined last in the molecule topology"
                write(*,*) "this is not true for molecule ", molecule_type_name(i_mole)
                stop
             end if
          end do loop2
       end if
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
  ! changed with the call to evb_change_diabat_data_structure_topology subroutine
  ! global :: hydronium_molecule_index, atom_index, molecule_index
  !
  !*************************
  subroutine ms_evb_calculate_total_force_energy( system_data, molecule_data, atom_data, verlet_list_data, PME_data, file_io_data, n_output  )
    type(system_data_type), intent(inout)      :: system_data
    type(molecule_data_type), dimension(:), intent(inout)    :: molecule_data
    type(atom_data_type), intent(inout)        :: atom_data
    type(verlet_list_data_type), intent(inout) :: verlet_list_data
    type(PME_data_type)   , intent(inout)      :: PME_data
    type(file_io_data_type), intent(in)        :: file_io_data
    integer, intent(in)                        :: n_output

    integer :: i_mole_principle, i_mole_hydronium, new_i_mole_hydronium, principle_diabat
    integer   :: flag_junk

    if ( n_hydronium_molecules > 1 ) then
       stop "can't have more than 1 hydronium"
    end if


    ! in principle, we have this array if there is more than one excess proton. 
    ! right now, MS-EVB is implemented for only 1 proton

    i_mole_hydronium = hydronium_molecule_index(1)

    call construct_evb_hamiltonian( i_mole_hydronium,  system_data, molecule_data, atom_data, verlet_list_data, PME_data   )

    ! global energy and force, system_data%potential_energy and atom_data%force will be updated.  These are currently the energy and force of the
    ! principle diabat, but will be updated to the adiabatic energy and force of the diagonalized MS-EVB Hamilonian
    call diagonalize_evb_hamiltonian( principle_diabat, new_i_mole_hydronium,  i_mole_hydronium, system_data%potential_energy, atom_data%force, file_io_data, n_output )

    ! if principle hydronium molecule has changed, then change
    ! data structures.  Note we need to pass temporary arrays here, and then copy them back
    ! this is because the arrays passed in with intent(in), or intent(out) attribute will be passed in
    ! as pointers, so we can't pass (e.g.) n_atom, n_atom, without the first intent in array changing.
    ! which will cause a bug
    !
    ! note the forces currently are in the old hydronium molecule topology, these need to be changed to
    ! the new topology
    if ( new_i_mole_hydronium /= i_mole_hydronium ) then

       call evb_change_diabat_data_structure_topology( principle_diabat, i_mole_hydronium, system_data, molecule_data, atom_data, hydronium_molecule_index )

       !******* update verlet_list after a proton transfer
       call construct_verlet_list( verlet_list_data, atom_data, system_data%total_atoms, system_data%box, system_data%xyz_to_box_transform  )
       ! the "1" input to update_verlet_displacements signals to initialize the displacement array
       call update_verlet_displacements( system_data%total_atoms, atom_data%xyz, verlet_list_data , system_data%box, system_data%xyz_to_box_transform , flag_junk, 1 )

    endif

!!$    enddo

    ! deallocate stored force array
    deallocate( evb_forces_store )

    Select Case(electrostatic_type)
    Case("pme")
       deallocate( Q_grid_diabats )
    End Select


  end subroutine ms_evb_calculate_total_force_energy




  !**************************
  ! this subroutine calculates the adiabatic energy
  ! and force for the ms-evb model by diagonalizing the
  ! evb hamiltonian
  !**************************
  subroutine diagonalize_evb_hamiltonian( principle_diabat, new_i_mole_principle, i_mole_principle, adiabatic_potential, force_atoms, file_io_data, n_output )
    integer, intent(out) :: principle_diabat, new_i_mole_principle
    integer, intent(in)  :: i_mole_principle
    real*8 , intent(out) :: adiabatic_potential
    real*8, dimension(:,:), intent(out) :: force_atoms
    type(file_io_data_type),intent(in) :: file_io_data
    integer, intent(in) :: n_output

    integer :: i_state, j_state, n_rot, ground_state_index, index, i_hop
    real*8 :: factor, state_coefficient
    real*8, dimension(:), allocatable  :: eigenvalues, ground_state_eigenvector
    real*8, dimension(:,:), allocatable :: hamiltonian, eigenvectors

    integer, save :: i_step = 0

    ! when the evb_hamiltonian was constructed, it was constructed as an upper triangular matrix
    ! fill in lower triangle.  Also we most probably have fewer diabats than the setting of
    ! evb_max_states, in which case the evb_hamiltonian is partially empty
    ! the number of total diabats is stored in the diabat_index variable

    allocate( eigenvalues(diabat_index) , ground_state_eigenvector(diabat_index), hamiltonian(diabat_index, diabat_index) , eigenvectors(diabat_index,diabat_index) )


!!$    write(*,*) "hamiltonian"
    hamiltonian=0d0
    do i_state = 1 , diabat_index
       do j_state = i_state , diabat_index
          hamiltonian( i_state, j_state ) = evb_hamiltonian(i_state,j_state)
          hamiltonian(j_state, i_state ) = hamiltonian(i_state, j_state )
       end do
!!$       write(*,*) hamiltonian( i_state , : )
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

!!$    write(*,*) "ground state eigenvector"
!!$    write(*,*) ground_state_eigenvector

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

          force_atoms(:,:) = force_atoms(:,:) + factor * evb_forces_store(:,:,index)
          
          end if
       end do
    end do


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

    do i_hop=1, evb_max_chain
       if ( evb_diabat_proton_log( principle_diabat , i_hop, 1 ) < 0 ) exit 
       ! acceptor molecule is in 4th index of evb_diabat_proton_log
       new_i_mole_principle = evb_diabat_proton_log( principle_diabat , i_hop, 4 )
    enddo

    !**** print evb trajectory info if requested ******
    Select Case( print_ms_evb_data )
    Case("yes")
       if ( new_i_mole_principle /= i_mole_principle ) then
          write(file_io_data%ofile_hop_file_h,*) "step ", trajectory_step
          write(file_io_data%ofile_hop_file_h,*) "proton hop from ", i_mole_principle, " to ", new_i_mole_principle
       end if
       ! print
       if ( mod( i_step, n_output ) == 0 ) then
          ! if restart, don't print for initial configuration
          if ( ( restart_trajectory .eq. "no" ) .or. ( i_step > 0 ) ) then
             call print_evb_trajectory_data( ground_state_eigenvector , i_mole_principle, file_io_data )
          end if
       endif
    End Select

    ! update step in case of printing
    i_step = i_step + 1

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
  ! for every diabatic state, we will allocate new data structures for 
  !
  ! system_data
  ! atom_data
  ! molecule_data
  ! 
  !************************************
  subroutine construct_evb_hamiltonian( i_mole_principle, system_data, molecule_data, atom_data, verlet_list_data, PME_data )
    use total_energy_forces
    integer, intent(in) :: i_mole_principle
    type(system_data_type) , intent(inout)      :: system_data
    type(molecule_data_type), dimension(:), intent(in)       :: molecule_data
    type(atom_data_type)   , intent(inout)      :: atom_data
    type(verlet_list_data_type), intent(inout)  :: verlet_list_data
    type(PME_data_type) , intent(inout)         :: PME_data
    
    integer ::  diabat_index_donor,  i_mole_donor, store_index
    real*8  ::  E_ms_evb_repulsion, E_reference

    ! initialize proton hop log
    evb_diabat_proton_log=-1
    ! zero hamiltonian
    evb_hamiltonian=0d0

    ! initialize diabat_index
    diabat_index=1
    diabat_index_donor=1
    i_mole_donor = i_mole_principle


    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "ms-evb energy and force for principle diabat started at", time
    endif
    !***********************************************************************!

    !*********************** energy and force of primary diabat ******************************!
    call calculate_total_force_energy( system_data, molecule_data, atom_data, verlet_list_data, PME_data )

    ! need to add contribution from special evb repulsive interactions
    call ms_evb_intermolecular_repulsion( atom_data%force, E_ms_evb_repulsion, system_data, atom_data, molecule_data, hydronium_molecule_index )

    system_data%potential_energy = system_data%potential_energy + E_ms_evb_repulsion

    ! now add reference chemical energy for this adiabatic state
    call get_adiabatic_reference_energy( E_reference, molecule_data(i_mole_donor)%molecule_type_index )

    system_data%potential_energy = system_data%potential_energy + E_reference



    !****************************************************************************************

    ! store the forces, the "1" is to initialize
    store_index=1
    call evb_store_forces( i_mole_principle, diabat_index, diabat_index, system_data%n_mole, system_data%total_atoms, molecule_data, atom_data, store_index, 1 )
    ! store the energy
    evb_hamiltonian( diabat_index , diabat_index ) = system_data%potential_energy

    ! for diabat PME calculations...
    ! the "1" is to initialize the update Q_grid subroutine, this needs to be initialized because it stores the Q_grids for each diabat to save computation time.  Q_grid_junk is junk argument
    allocate(Q_grid_junk(PME_data%pme_grid,PME_data%pme_grid,PME_data%pme_grid))
    call ms_evb_diabat_lookup_Q_grid( Q_grid_junk, diabat_index, 1 )
    deallocate(Q_grid_junk)



    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "ms-evb energy and force for principle diabat finished at", time
    endif
    !***********************************************************************!

    ! this is a recursive subroutine that loops through and constructs all diabatic states using recursive proton hops
    ! the main purpose of this subroutine is to construct the array evb_diabat_proton_log, which contains all the information necessary to
    ! characterize a particular diabat, by storing the indices of all the consecutive proton donors and acceptors for each hop
    evb_diabat_coupling_matrix=-1 ! initialize to negative integer
    call evb_conduct_proton_transfer_recursive( i_mole_donor, diabat_index_donor, system_data , molecule_data, atom_data )


    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "ms-evb diabat calculations started at", time
    endif
    !***********************************************************************!


    ! now update energies for each diabat.  This takes care of everything except the reciprocal space component, which is done below
    call evb_hamiltonian_elements_donor_acceptor( i_mole_principle, system_data , molecule_data, atom_data, PME_data )


    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "ms-evb reciprocal space started at", time
    endif
    !***********************************************************************!

    ! we have stored the dQ_dr array from the principle diabat, so 
    ! we can use a fast update for the reciprocal space electrostatic energy of each diabat
    ! We parallelize this part, as it is expensive.  At this point, we should have 
    ! Q_grids stored for every diabat, and we have the dQ_dr grid
    call calculate_reciprocal_space_pme( i_mole_principle, system_data , molecule_data, atom_data, PME_data )


    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "ms-evb reciprocal space finished at", time
    endif
    !***********************************************************************!

  end subroutine construct_evb_hamiltonian




  !***********************************************
  ! this subroutine constructs all diabatic states using recursive proton hops
  !
  ! note that we don't run into the problem that an acceptor molecule back transfers to its donor
  ! when the acceptor molecule becomes the next donor, because in the find_evb_reactive_neighbors
  ! subroutine we always input the original principle "xyz" geometry array
  !***********************************************
  recursive subroutine evb_conduct_proton_transfer_recursive( i_mole_donor, diabat_index_donor,  i_mole_principle, system_data, molecule_data, atom_data )
    integer, intent(in) :: i_mole_donor, diabat_index_donor, i_mole_principle
    type(system_data_type), intent(in) :: system_data
    type(molecule_data_type), dimension(:), intent(in) :: molecule_data
    type(atom_data_type) ,  intent(in) :: atom_data

    ! for second index, first entry is molecule index, second is atom index
    ! in general, we could have ifferent proton transfers to different atoms of
    ! same molecule
    integer, dimension(evb_max_neighbors,2) :: evb_neighbor_list
    integer :: i_molecule, i_atom, j_atom,i_mole_acceptor, i_atom_acceptor, i_hop, diabat_index_acceptor, i_mole_type, flag_reactive_atom, flag_cycle, index, count


    ! first see if we've reached the maximum limit for the proton hopping chain, figure this out by the number of stored elements in proton hopping log
    count = 0
    do i_hop=1, evb_max_chain
       if ( evb_diabat_proton_log( diabat_index_donor , i_hop, 1 ) < 0 ) exit
       count = count + 1
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

       i_mole_type = molecule_data(i_mole_donor)%molecule_type_index
       ! loop over all hydrogen atoms for this donor
       do i_atom=1, molecule_data(i_mole_donor)%n_atom

          flag_reactive_atom = evb_reactive_protons(i_mole_type,i_atom)
          ! if reactive atom (proton), consider all possible diabats formed from the transfer of this proton
          if ( flag_reactive_atom == 1 ) then

             ! first identify molecules within first solvation shell to which we can transfer a proton
             call find_evb_reactive_neighbors( i_mole_donor, i_atom, evb_neighbor_list, system_data, molecule_data, atom_data )


             ! loop over neighbors
             loop1:          do i_molecule = 1 , evb_max_neighbors
                if ( evb_neighbor_list(i_molecule,1) < 0 ) exit loop1

                ! starting a new diabat, update global diabat_index counter
                diabat_index = diabat_index + 1
                diabat_index_acceptor = diabat_index

                ! make sure we don't have too many states
                call evb_check_number_diabats  

                ! store donor diabat index for this acceptor diabat
                evb_diabat_coupling_matrix(diabat_index) = diabat_index_donor              

                i_mole_acceptor = evb_neighbor_list(i_molecule,1)
                i_atom_acceptor = evb_neighbor_list(i_molecule,2) 

!!$                write(*,*) "diabat ", diabat_index
!!$                write(*,*) i_mole_donor, i_atom , i_mole_acceptor, i_atom_acceptor

                ! check to see if this is a hydronium molecule (from i.e. cyclic transfer)
                flag_cycle = check_hydronium_molecule( i_mole_acceptor )

                ! create proton hop log for this acceptor, copy previous proton transfer log from donor
                ! the proton log stores 1) index of donor molecule, 2) index of transferring proton, 
                ! 3) index of basic atom of donor molecule to which donating proton is attached, 4) index of acceptor molecule, 
                ! 5) index of basic atom of acceptor molecule
                loop2:             do i_hop=1, evb_max_chain
                   if ( evb_diabat_proton_log( diabat_index_donor , i_hop, 1 ) < 0 ) exit loop2
                   ! copy previous hop information from donor diabat
                   evb_diabat_proton_log(diabat_index_acceptor,i_hop,:) = evb_diabat_proton_log(diabat_index_donor,i_hop,:)
                enddo loop2

                ! now add new proton transfer information
                ! find basic bonded atom of donating proton
                call find_bonded_atom_hydrogen( i_mole_type, molecule_data(i_mole_donor)%n_atom, i_atom, j_atom )
                index=count+1
                evb_diabat_proton_log(diabat_index_acceptor,index,1) = i_mole_donor
                evb_diabat_proton_log(diabat_index_acceptor,index,2) = i_atom
                evb_diabat_proton_log(diabat_index_acceptor,index,3) = j_atom
                evb_diabat_proton_log(diabat_index_acceptor,index,4) = i_mole_acceptor
                evb_diabat_proton_log(diabat_index_acceptor,index,5) = i_atom_acceptor                

                ! now call this recursive subroutine, with the proton acceptor molecule as the new donor, unless this is the end of a cyclic transfer
                if ( flag_cycle < 1 ) then
                   call evb_conduct_proton_transfer_recursive( i_mole_acceptor, diabat_index_acceptor,  i_mole_principle, system_data, molecule_data, atom_data )
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
  ! storing the necessary output
  !
  ! for each diabat, we make local copies of 
  !  system_data   ==> system_data_diabat
  !  molecule_data ==> molecule_data_diabat
  !  atom_data     ==> molecule_data_diabat
  !
  !  where the new data structures take on the topology
  !  of the particular diabat
  !****************************************
  subroutine evb_hamiltonian_elements_donor_acceptor( i_mole_principle, system_data, molecule_data, atom_data, PME_data )
    integer, intent(in) :: i_mole_principle
    type(system_data_type),intent(in) :: system_data
    type(molecule_data_type),dimension(:),intent(in) :: molecule_data
    type(atom_data_type), intent(in)                 :: atom_data
    type(PME_data_type), intent(in)                  :: PME_data

    integer :: initialize, lookup_index, i_mole, i_atom, i_diabat, diabat_index_donor


    ! temporary diabat data structures
    ! note atom_data_diabat_type does not contain pointers, which is what we want, unlike atom_data_type
    type(atom_data_diabat_type)                             :: atom_data_diabat
    type(molecule_data_type   ),dimension(:), allocatable   :: molecule_data_diabat
    type(system_data_type     )                             :: system_data_diabat 

    integer, dimension( n_proton_max ), save :: hydronium_molecule_index_temp
    integer :: split_do


    ! decide how to split the parallel section
    if (n_threads .eq. 1 ) then
       split_do = 1
    else
       ! integer arithmetic
       if ( diabat_index > n_threads ) then
          split_do = diabat_index/n_threads
       else
          split_do = 1
       endif
    endif

    ! loop over diabats, excluding principle diabat
    call OMP_SET_NUM_THREADS(n_threads)
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(n_threads,split_do, system_data, molecule_data, atom_data, PME_data,  hydronium_molecule_index, i_mole_principle, diabat_index,  evb_hamiltonian, evb_diabat_coupling_matrix, store_index ) 
    !$OMP DO SCHEDULE(dynamic, split_do)
    do i_diabat = 2, diabat_index

       call evb_create_diabat_data_structures( atom_data_diabat, molecule_data_diabat, system_data_diabat, atom_data, molecule_data , system_data , hydronium_molecule_index_temp, hydronium_molecule_index )

       !************************** calculate diagonal matrix element energy and forces for this diabat
       ! important note, after call to ms_evb_diabat_force_energy, the topology of the data structures will be changed from the donor to acceptor topology
       call ms_evb_diabat_force_energy( system_data_diabat, i_diabat, i_mole_principle, atom_data_diabat, molecule_data_diabat, PME_data, hydronium_molecule_index_temp)

       ! store the forces, this needs to be in critical section, because forces will be stored using store_index, and then store_index will be incremented
       !$OMP CRITICAL
         call evb_store_forces( i_mole_principle, i_diabat, i_diabat, system_data_diabat%n_mole, system_data_diabat%total_atoms, molecule_data_diabat, atom_data_diabat, store_index )
       !$OMP END CRITICAL
       ! store the energy
       evb_hamiltonian( i_diabat , i_diabat ) = system_data_diabat%potential_energy

       !************************************* calculate off-diagonal diabatic coupling, energy and force
       ! here, the data structures should be in acceptor topology, after the call to ms_evb_diabat_force_energy
       ! note the atom_data_diabat%force array and system_data_diabat%potential_energy will be reinitialized to zero in this subroutine, before computing off-diagonal terms
       call evb_diabatic_coupling( i_diabat, i_mole_principle, system_data_diabat, atom_data_diabat, molecule_data_diabat )

       ! get donor diabat for which this coupling was calculated
       diabat_index_donor = evb_diabat_coupling_matrix(i_diabat)
       ! store the forces, note this coupling is between the current diabat, and the diabat_index_1neighbor (principle) diabat, again needs to be in critical section
       !$OMP CRITICAL
       call evb_store_forces( i_mole_principle, diabat_index_donor, i_diabat, system_data_diabat%n_mole, system_data_diabat%total_atoms, molecule_data_diabat, atom_data_diabat, store_index  )
       !$OMP END CRITICAL
       ! store the energy
       evb_hamiltonian( diabat_index_donor , i_diabat ) = system_data_diabat%potential_energy

    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL


  end subroutine evb_hamiltonian_elements_donor_acceptor




  !************************************
  ! this subroutine searches for evb reactive neighbors
  ! of reactive atom i_mole, i_atom, and puts them in neighbor list
  !
  ! output variable flag_cycle is given a value of "1" if we transfer a proton
  ! to a hydronium oxygen (in a cyclic transfer).  This is to prevent a new cycle of the original
  ! cyclic proton transfer
  !************************************
  subroutine find_evb_reactive_neighbors( i_mole, i_atom, neighbor_list, system_data, molecule_data, atom_data ) 
    integer, intent(in) :: i_mole, i_atom
    integer, dimension(:,:), intent(out) :: neighbor_list
    type(system_data_type), intent(in) :: system_data
    type(molecule_data_type), dimension(:),intent(in) :: molecule_data
    type(atom_data_type), intent(in) :: atom_data 

    !******** this is a local data structure with pointers that will be set
    ! to subarrays of atom_data arrays for the specific atoms in the molecule
    type(single_molecule_data_type) :: single_molecule_data_i, single_molecule_data_j

    integer :: j_mole, j_atom, index, atom_id1, atom_id2, j_mole_type
    real*8, dimension(3) :: r_com_i, r_com_j, dr_com, rij, shift


    ! set pointers for input molecule to atom_data data structure section
    call return_molecule_block( single_molecule_data_i , molecule_data(i_mole)%n_atom, molecule_data(i_mole)%atom_index, atom_xyz=atom_data%xyz )

    index=1
    ! initialize list to negative value to signal end of neighbors
    neighbor_list=-1

    do j_mole = 1, system_data%n_mole
       j_mole_type = molecule_data(j_mole)%molecule_type_index
       if ( i_mole /= j_mole ) then
          ! center of mass distance
          r_com_i(:) = molecule_data(i_mole)%r_com(:)
          r_com_j(:) = molecule_data(j_mole)%r_com(:)

          shift = pbc_shift( r_com_i, r_com_j, system_data%box, system_data%xyz_to_box_transform)
          dr_com = pbc_dr( r_com_i, r_com_j, shift)

          if ( dot_product( dr_com, dr_com ) < evb_first_solvation_cutoff **2 ) then
             ! loop over all atoms of  target molecule to find acceptor atoms
             ! in general, we could have multiple acceptor atoms for the same acceptor molecule
             ! and if this is the case these are stored as different neighbors
             loop2:          do j_atom = 1 , molecule_data(j_mole)%n_atom

                if ( evb_reactive_basic_atoms(j_mole_type, j_atom) == 1 ) then
                   ! see if this reactivity pair is within the cutoff distance
                   rij = pbc_dr( single_molecule_data_i%xyz(:,i_atom), single_molecule_data_j%xyz(:,j_atom), shift ) ! Note for COM cutoff, shift values unchanged
                   if ( dot_product(rij,rij) < evb_reactive_pair_distance**2 ) then
                      ! this is a reactive pair
                      neighbor_list(index,1)=j_mole
                      neighbor_list(index,2)=j_atom
                      index=index+1
                   endif
                endif
             enddo loop2
          end if
       end if
    end do


  end subroutine find_evb_reactive_neighbors



  !*********************************************
  ! this subroutine allocates diabat data structures
  ! and copies over data from primary diabat structures
  !*********************************************
  subroutine evb_create_diabat_data_structures( atom_data_diabat, molecule_data_diabat, system_data_diabat, atom_data, molecule_data , system_data, hydronium_molecule_index_local, hydronium_molecule_index )
    type(atom_data_diabat_type), intent(out) :: atom_data_diabat
    type(molecule_data_type),dimension(:),intent(out) :: molecule_data_diabat
    type(system_data_type), intent(out)      :: system_data_diabat
    type(atom_data_type), intent(in) :: atom_data_diabat
    type(molecule_data_type),dimension(:),intent(in) :: molecule_data_diabat
    type(system_data_type), intent(in)      :: system_data_diabat
    integer,dimension(:), intent(out)       :: hydronium_molecule_index_local
    integer,dimension(:), intent(in)        :: hydronium_molecule_index

    integer :: total_atoms, n_mole, i_mole

    total_atoms = system_data%total_atoms
    n_mole      = system_data%n_mole

    ! first allocate arrays in atom_data_diabat structure
    allocate( atom_data_diabat%xyz(3,total_atoms) )
    allocate( atom_data_diabat%force(3,total_atoms) )
    allocate( atom_data_diabat%charge(total_atoms) )
    allocate( atom_data_diabat%atom_type_index(total_atoms) )
    allocate( atom_data_diabat%aname(total_atoms) )

    ! now copy atom_data
    atom_data_diabat%xyz = atom_data%xyz
    atom_data_diabat%force = atom_data%force
    atom_data_diabat%charge = atom_data%charge
    atom_data_diabat%atom_type_index = atom_data%atom_type_index
    atom_data_diabat%aname = atom_data%aname

    ! now allocate molecule_data_diabat and copy molecule_data 
    allocate(molecule_data_diabat(n_mole))
    do i_mole=1,n_mole
       molecule_data_diabat(i_mole) = molecule_data(i_mole)
    enddo

    system_data_diabat = system_data
    hydronium_molecule_index_local = hydronium_molecule_index


  end subroutine evb_create_diabat_data_structures




  !***********************************************
  ! this subroutine modifies data structures to be consistent with
  ! the molecular topology of a particular diabat (diabat_index)
  !
  ! the new data structures are created based on the proton hopping information contained in evb_diabat_proton_log
  !***********************************************
  subroutine evb_change_diabat_data_structure_topology( diabat_index, i_mole_principle , system_data, molecule_data, atom_data, hydronium_molecule_index )
    integer, intent(in) :: i_mole_principle, diabat_index
    type(system_data_type), intent(in) :: system_data
    type(molecule_data_type), dimension(:), intent(inout) :: molecule_data
    type(atom_data_type),intent(inout) :: atom_data
    integer, dimension(:), intent(inout) :: hydronium_molecule_index

    integer :: i_hop, i_mole_donor, i_atom_donor, i_mole_acceptor, i_atom_acceptor, i_heavy_acceptor, i_atom


    ! first donor is always principle hydronium
    ! set i_mole_acceptor equal to this, because
    ! i_mole_donor will be copied from previous i_mole_acceptor
    i_mole_acceptor = i_mole_principle

    ! loop through all proton hops.  A negative integer signals the end of the chain
    loop1: do i_hop=1, evb_max_chain
       if ( evb_diabat_proton_log( diabat_index , i_hop, 1 ) < 0 ) exit loop1
       i_mole_donor = i_mole_acceptor
       i_atom_donor = evb_diabat_proton_log( diabat_index , i_hop, 2 )
       i_mole_acceptor = evb_diabat_proton_log( diabat_index , i_hop, 4 )
       i_heavy_acceptor = evb_diabat_proton_log( diabat_index , i_hop, 5 )
       i_atom_acceptor = n_atom_temp(i_mole_acceptor) + 1

       ! now update data structures for this proton hop from donor to acceptor
       call evb_change_data_structures_proton_transfer( i_mole_donor, i_atom_donor, i_mole_acceptor, i_atom_acceptor, i_heavy_acceptor, atom_data, molecule_data, system_data, hydronium_molecule_index )

    enddo loop1


  end subroutine evb_change_diabat_data_structure_topology





  !*************************************
  !
  ! this subroutine modifies atom_data_diabat and molecule_data_diabat
  ! data structures to account for the new molecular topology from a proton hop from
  ! i_mole_donor (i_atom_donor) to i_mole_acceptor (i_atom_acceptor)
  !
  !**************************************
  subroutine evb_change_data_structures_proton_transfer( i_mole_donor, i_atom_donor, i_mole_acceptor, i_atom_acceptor, i_heavy_acceptor, atom_data_diabat, molecule_data_diabat, system_data_diabat, hydronium_molecule_index_local )
    integer, intent(in) :: i_mole_donor, i_atom_donor, i_mole_acceptor, i_atom_acceptor, i_heavy_acceptor
    type(system_data_type) , intent(in)                   :: system_data_diabat
    type(atom_data_diabat_type), intent(inout)                   :: atom_data_diabat
    type(molecule_data_type), dimension(:), intent(inout) :: molecule_data_diabat
    integer, dimension(:), intent(inout) :: hydronium_molecule_index_local

    !******** this is a local data structure with pointers that will be set
    ! to subarrays of atom_data arrays for the specific atoms in the molecule
    type(single_molecule_data_type) :: single_molecule_data_donor, single_molecule_data_acceptor

    integer :: i_atom_donor_global, i_atom_acceptor_global
    integer :: i_atom, i, count, i_index, i_type


    ! need to change this if more than one hydronium
    if ( n_hydronium_molecules > 1 ) stop "see code below"
    hydronium_molecule_index_local(1) = i_mole_acceptor

      
    ! these are indices of atoms in global atom_data arrays
    i_atom_donor_global =   molecule_data_diabat(i_mole_donor)%atom_index(i_atom_donor)
    i_atom_acceptor_global = molecule_data_diabat(i_mole_acceptor)%atom_index(i_atom_acceptor)


    ! this subroutine shifts data structures for atom transfer.
    ! Here, we're transfering an H atom back from the donor to acceptor
    ! and modifying data structures for this proton transfer
    call shift_array_data_donor_acceptor_transfer( i_mole_donor, i_atom_donor_global, i_mole_acceptor, i_atom_acceptor_global, molecule_data_diabat, atom_data_diabat%xyz, atom_data_diabat%velocity, atom_data_diabat%force, atom_data_diabat%mass , atom_data_diabat%charge , atom_data_diabat%atom_type_index , atom_data_diabat%aname )

    ! set pointers for this data structure to target molecule
    ! molecule structure is for new, post proton transfer topology
    call return_molecule_block( single_molecule_data_donor , molecule_data_diabat(i_mole_donor)%n_atom, molecule_data_diabat(i_mole_donor)%atom_index, atom_xyz=atom_data_diabat%xyz, atom_mass=atom_data_diabat%mass, atom_charge=atom_data_diabat%charge, atom_type_index=atom_data_diabat%atom_type_index, atom_name=atom_data_diabat%aname )

    call return_molecule_block( single_molecule_data_acceptor , molecule_data_diabat(i_mole_acceptor)%n_atom, molecule_data_diabat(i_mole_acceptor)%atom_index, atom_xyz=atom_data_diabat%xyz, atom_mass=atom_data_diabat%mass, atom_charge=atom_data_diabat%charge, atom_type_index=atom_data_diabat%atom_type_index, atom_name=atom_data_diabat%aname )

    ! here it's possible that the transfered proton could be split from the rest of the
    ! molecule by a periodic image. We need to fix this, as energy routines(bonds, angles, dihedrals)
    ! assume the molecule to not be split over pbc conditions
    call make_molecule_whole( molecule_data_diabat(i_mole_acceptor)%n_atom, single_molecule_data_acceptor%xyz, system_data_diabat%box, system_data_diabat%xyz_to_box_transform)

    ! update center of mass for donor and acceptor after proton transfer
    molecule_data_diabat(i_mole_donor)%r_com = pos_com( single_molecule_data_donor%xyz, molecule_data_diabat(i_mole_donor)%n_atom , single_molecule_data_donor%mass )
    molecule_data_diabat(i_mole_acceptor)%r_com = pos_com( single_molecule_data_acceptor%xyz, molecule_data_diabat(i_mole_acceptor)%n_atom , single_molecule_data_acceptor%mass )


    ! modify atom_type_index of proton if we are transfering to a different base (not its conjugate base), we can't just copy the index
    call change_proton_index_proton_transfer( single_molecule_data_acceptor%atom_type_index, i_atom_acceptor , molecule_data_diabat(i_mole_acceptor)%molecule_type_index )
    ! single_molecule_data_acceptor%atom_type_index(i_atom_acceptor) should now be filled in with correct proton atom type

    ! now fill in all other acceptor atom types with conjugate base/acid atomtypes
    do i_atom = 1 , molecule_data_diabat(i_mole_acceptor)%n_atom
      i_type = single_molecule_data_acceptor%atom_type_index(i_atom) 
      if ( i_atom /= i_atom_acceptor ) then
          single_molecule_data_acceptor%atom_type_index(i_atom) = evb_conjugate_atom_index(i_type)
       end if
       ! fill in charges of new atom types of newly formed acid.  atype_chg is a global variable array
       single_molecule_data_acceptor%charge(i_atom) = atype_chg(i_type)
    enddo

    ! now map the heavy atom of the acceptor to it's specific atom type
    i_index = evb_conjugate_pairs( molecule_data_diabat(i_mole_acceptor)%molecule_type_index )
    single_molecule_data_acceptor%atom_type_index(i_heavy_acceptor) = evb_heavy_acid_index(i_index)

    ! now map the donor atom index
    do i_atom = 1 , molecule_data_diabat(i_mole_donor)%n_atom
       i_type = single_molecule_data_donor%atom_type_index(i_atom)
       single_molecule_data_donor%atom_type_index(i_atom) = evb_conjugate_atom_index(i_type)
       ! fill in charges of new atom types of newly formed base.  atype_chg is a global variable array
       single_molecule_data_donor%charge(i_atom) = atype_chg(i_type)
    enddo

    ! molecule index
    molecule_data_diabat(i_mole_donor)%molecule_type_index = evb_conjugate_pairs( molecule_data_diabat(i_mole_donor)%molecule_type_index )
    molecule_data_diabat(i_mole_acceptor)%molecule_type_index = evb_conjugate_pairs( molecule_data_diabat(i_mole_acceptor)%molecule_type_index )  

    ! may need to reorder data structures for acceptor molecule based on index of acceptor
    ! heavy atom to be consistent with molecule_type array
    call reorder_molecule_data_structures( molecule_data_diabat(i_mole_acceptor), single_molecule_data_acceptor )


  end subroutine evb_change_data_structures_proton_transfer



  !************************************************
  ! this subroutine reorders data structures
  ! if they are inconsistent with molecule_type array
  ! for example, if we have donated a proton to an
  ! R-SO3 group, the new O-H oxygen must be correctly ordered
  ! in the data structure array
  !************************************************
  subroutine reorder_molecule_data_structures( molecule_data_acceptor , single_molecule_data_acceptor )
    type(molecule_data_type) , intent(inout) :: molecule_data_acceptor
    type(single_molecule_data_type) , intent(inout) :: single_molecule_data_acceptor

    integer :: i_mole_type,i_atom,j_atom, index, count, i
    real*8 :: xyz_atom(3) , force_atom(3) , charge_atom, mass_atom, velocity_atom(3)
    integer :: atom_index_atom
    character(MAX_ANAME) :: aname_atom

    i_mole_type = molecule_data_acceptor%molecule_type_index

    ! loop over all atoms in molecule_type array
    loop1 : do i_atom=1, MAX_N_ATOM_TYPE
       if ( molecule_type(i_mole_type,i_atom) == MAX_N_ATOM_TYPE + 1 ) exit loop1

       ! if atomtypes don't correspond, then we need to shift data
       if ( molecule_type(i_mole_type,i_atom) /= single_molecule_data_acceptor%atom_type_index(i_atom) ) then

          ! loop through remaining atoms to find first occurance of this atomtype
          index=-1
          loop2: do j_atom = i_atom+1, molecule_data_acceptor%n_atom
             if ( molecule_type(i_mole_type,i_atom) == single_molecule_data_acceptor%atom_type_index(j_atom) ) then
                index = j_atom
                exit loop2
             endif
          enddo loop2

          if ( index == -1 ) stop "error in subroutine reorder_molecule_data_structures"

          ! now store data, then shift all data structures up
          xyz_atom(:) = single_molecule_data_acceptor%xyz(:,index)
          velocity_atom(:) = single_molecule_data_acceptor%velocity(:,index)
          force_atom(:) = single_molecule_data_acceptor%force(:,index)
          charge_atom = single_molecule_data_acceptor%charge(index)
          mass_atom = single_molecule_data_acceptor%mass(index)
          atom_index_atom = single_molecule_data_acceptor%atom_type_index(index)
          aname_atom = single_molecule_data_acceptor%aname(index)

          count=1
          do i = i_atom+1 , index
             ! start at end of array so we overwrite after we copy
             j_atom = index - count

             single_molecule_data_acceptor%xyz(:,j_atom+1) = single_molecule_data_acceptor%xyz(:,j_atom)
             single_molecule_data_acceptor%velocity(:,j_atom+1) = single_molecule_data_acceptor%velocity(:,j_atom)
             single_molecule_data_acceptor%force(:,j_atom+1) = single_molecule_data_acceptor%force(:,j_atom)
             single_molecule_data_acceptor%charge(j_atom+1) = single_molecule_data_acceptor%charge(j_atom)
             single_molecule_data_acceptor%mass(j_atom+1) = single_molecule_data_acceptor%mass(j_atom)
             single_molecule_data_acceptor%atom_type_index(j_atom+1) = single_molecule_data_acceptor%atom_type_index(j_atom)
             single_molecule_data_acceptor%aname(j_atom+1) = single_molecule_data_acceptor%aname(j_atom)            
             count=count+1
          enddo


          single_molecule_data_acceptor%xyz(:,i_atom) = xyz_atom(:)
          single_molecule_data_acceptor%velocity(:,i_atom) = velocity_atom(:)
          single_molecule_data_acceptor%force(:,i_atom) = force_atom(:)
          single_molecule_data_acceptor%charge(i_atom) =  charge_atom
          single_molecule_data_acceptor%mass(i_atom) = mass_atom
          single_molecule_data_acceptor%atom_type_index(i_atom) =atom_index_atom
          single_molecule_data_acceptor%aname(i_atom) = aname_atom

       end if
    end do loop1


  end subroutine reorder_molecule_data_structures




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
  subroutine evb_diabatic_coupling( diabat_index, i_mole_principle, system_data_diabat, atom_data_diabat, molecule_data_diabat )
    integer, intent(in) :: i_mole_principle, diabat_index
    type(system_data_type), intent(inout)  :: system_data_diabat
    type(atom_data_diabat_type), intent(inout)    :: atom_data_diabat
    type(molecule_data_type), dimension(:), intent(inout) :: molecule_data_diabat   

    !******** this is a local data structure with pointers that will be set
    ! to subarrays of atom_data arrays for the specific atoms in the molecule
    type(single_molecule_data_type) :: single_molecule_data_donor, single_molecule_data_acceptor

    integer ::  i_hop, i_mole_donor, i_proton_donor, i_atom_donor, i_mole_acceptor, i_atom_acceptor, i_atom, total_atoms
    real*8 ::  Vex , Vconstij, A
    real*8, dimension(:,:), allocatable :: dVex
    real*8, dimension(3,3) :: dA

    ! zero potential and force as 
    system_data_diabat%potential_energy=0d0
    atom_data_diabat%force=0d0

    total_atoms = system_data_diabat%total_atoms
    allocate( dVex( 3 , total_atoms ) )

    ! figure out proton donor and proton acceptor from evb_diabat_proton_log.
    ! we want the final donor and acceptor, which is the last proton transfer

    ! first donor is always principle hydronium
    ! set i_mole_acceptor equal to this, because
    ! i_mole_donor will be copied from previous i_mole_acceptor
    i_mole_acceptor = i_mole_principle

    !***** comment on named variables
    ! here, i_atom_donor refers to the heavy atom on the donor molecule that is attached to the proton
    ! i_proton_donor, is the proton, and i_atom_acceptor is the basic atom on the acceptor that
    ! is accepting the proton


    ! loop through all proton hops.  A negative integer signals the end of the chain
    do i_hop=1, evb_max_chain
       if ( evb_diabat_proton_log( diabat_index , i_hop, 1 ) < 0 ) exit
       i_mole_donor = i_mole_acceptor
       i_proton_donor = evb_diabat_proton_log( diabat_index , i_hop , 2 )
       i_mole_acceptor = evb_diabat_proton_log( diabat_index , i_hop, 4 )

    enddo



    ! first calculate solvent dependent prefactor, and its derivative w.r.t. atomic coordinates
    call evb_diabatic_coupling_electrostatics( Vex, dVex, i_mole_donor, i_mole_acceptor, system_data_diabat, atom_data_diabat, molecule_data_diabat )


    ! set pointers for acceptor and donor molecules to atom_data data structures
    call return_molecule_block( single_molecule_data_donor , molecule_data_diabat(i_mole_donor)%n_atom, molecule_data_diabat(i_mole_donor)%atom_index, atom_xyz=atom_data_diabat%xyz, atom_force=atom_data_diabat%force, atom_type_index=atom_data_diabat%atom_type_index )
    call return_molecule_block( single_molecule_data_acceptor , molecule_data_diabat(i_mole_acceptor)%n_atom, molecule_data_diabat(i_mole_acceptor)%atom_index, atom_xyz=atom_data_diabat%xyz, atom_force=atom_data_diabat%force, atom_type_index=atom_data_diabat%atom_type_index )

    ! now calculate geometry dependent scale factor and its derivatives
    ! three derivatives are contained in the output dA array, in order:
    ! i_atom_donor, i_atom_acceptor, i_proton_donor
    call evb_diabatic_coupling_geometric( A , dA, Vconstij , i_mole_donor, i_atom_donor, i_mole_acceptor, i_atom_acceptor, system_data_diabat, molecule_data_diabat, single_molecule_data_donor, single_molecule_data_acceptor )

    ! now form evb diabatic coupling matrix element
    system_data_diabat%potential_energy = ( Vconstij + Vex ) * A

    ! use pointers to molecular blocks of atomic_data_diabat%force data structure to fill in forces
    ! first forces from geometric derivative
    ! O donor
    single_molecule_data_donor%force(:,i_atom_donor) = - ( Vconstij + Vex ) * dA(:,1)
    ! O acceptor
    single_molecule_data_acceptor%force(:,i_atom_acceptor) = - ( Vconstij + Vex ) * dA(:,2)
    ! H central Zundel 
    single_molecule_data_acceptor%force(:,molecule_data_diabat(i_mole_acceptor)%n_atom) = - ( Vconstij + Vex ) * dA(:,3)

    ! now forces from electrostatic derivative, here we can directly update atomic atom_data_diabat%force array
    atom_data_diabat%force(:,:) = atom_data_diabat%force(:,:) - dVex(:,:) * A

    deallocate( dVex )

  end subroutine evb_diabatic_coupling



  !********************************
  ! this returns the geometric term for the
  ! diabatic coupling element, as well as its
  ! derivatives with respect to the positions
  ! of the oxygen atoms of the donor and acceptor,
  ! and the central hydrogen atom
  ! 
  ! dA(:,1) contains derivative w.r.t  O_donor coordinates
  ! dA(:,2) contains derivative w.r.t. O_acceptor coordinates
  ! dA(:,3) contains derivative w.r.t. H zundel coordinates
  !********************************
  subroutine evb_diabatic_coupling_geometric( A , dA, Vconstij, i_mole_donor, i_atom_donor, i_mole_acceptor,  i_atom_acceptor, system_data_diabat, molecule_data_diabat, single_molecule_data_donor, single_molecule_data_acceptor  )
    real*8, intent(out) :: A, Vconstij
    real*8, dimension(:,:), intent(out) :: dA
    integer, intent(in)    :: i_mole_donor, i_mole_acceptor
    integer, intent(out)  :: i_atom_donor, i_atom_acceptor
    type(system_data_type), intent(inout)  :: system_data_diabat
    type(molecule_data_type), dimension(:), intent(in) :: molecule_data_diabat
    type(single_molecule_data_type), intent(in) :: single_molecule_data_donor, single_molecule_data_acceptor

    real*8 , dimension(3) :: r_O1, r_O2, r_H, r_OO , q, shift, r_ij
    integer :: itype1, itype2, itype3, index, function_type


    ! get heavy atoms involved in proton transfer, note acceptor should currently be in acid topology
    call get_heavy_atom_transfer_base( i_atom_donor, molecule_data_diabat(i_mole_donor)%molecule_type_index )
    call get_heavy_atom_transfer_acid( i_atom_acceptor, molecule_data_diabat(i_mole_acceptor)%molecule_type_index )

    ! first get distances, need r_OO and q, q is defined by
    ! q = ( rO1 + rO2 ) / 2 - rH
    ! to consider pbc, shift all atoms relative to the donor oxygen

    r_O1(:) = single_molecule_data_donor%xyz(:, i_atom_donor)
    r_O2(:) = single_molecule_data_acceptor%xyz(:, i_atom_acceptor)

    shift = pbc_shift( r_O1 , r_O2 , system_data_diabat%box , system_data_diabat%xyz_to_box_transform )
    r_ij = pbc_dr( r_O1 , r_O2 , shift )
    r_O2(:) = r_O1(:) + r_ij(:)

    ! H atom is last atom in acceptor, because we've already transferred the proton
    r_H(:) = single_molecule_data_acceptor%xyz(:, molecule_data_diabat(i_mole_acceptor)%n_atom )
    r_ij = pbc_dr( r_O1 , r_H , shift )
    r_H(:) = r_O1(:) + r_ij(:)

    r_OO = r_O1 - r_O2
    q = ( r_O1 + r_O2 ) / 2d0 - r_H

    ! get parameters for this diabatic coupling term
    ! the diabatic coupling should be symmetric in atom types,
    ! so parameters shouldn't depend on whether we are considering
    ! the donor or acceptor topology
    itype1 = single_molecule_data_donor%atom_type_index(i_atom_donor)
    itype2 = single_molecule_data_acceptor%atom_type_index(i_atom_acceptor)
    itype3 = single_molecule_data_acceptor%atom_type_index(molecule_data_diabat(i_mole_acceptor)%n_atom)

    call get_index_atom_set( index, evb_diabat_coupling_interaction, (/itype1, itype2 , itype3/) )
    if ( index == -1 ) stop "couldn't find index in subroutine 'get_index_atom_set'"

    ! get coupling function type
    function_type = evb_diabat_coupling_type(index)

    call evb_diabatic_coupling_function( A , Vconstij , dA , function_type , evb_diabat_coupling_parameters(index,:) , q , r_OO )


    ! Note dA derivative array is returned as follows:
    ! dA(:,1) is derivative w.r.t. O_donor coordinates
    ! dA(:,2) is derivative w.r.t. O_acceptor coordinates
    ! dA(:,3) is derivative w.r.t. H zundel coordinates

  end subroutine evb_diabatic_coupling_geometric




  !**********************************************
  ! this subroutine evaluates different functional forms of the diabatic_coupling geometric
  ! component
  !**********************************************
  subroutine evb_diabatic_coupling_function( A, Vconstij, dA, function_type, function_parameters, q, r_OO )
    real*8, intent(out) :: A, Vconstij
    real*8, dimension(:,:), intent(out) :: dA
    integer, intent(in) :: function_type
    real*8, dimension(:), intent(in) :: function_parameters
    real*8, dimension(3), intent(in) :: q , r_OO

    real*8  :: r_OO_mag, q2_mag, q_mag, fac1, fac2, fac3 , dfac1, dfac2, dfac3

    ! diabatic coupling parameters
    real*8 :: gamma, P, k, D, beta, R0, Pprime, alpha, rl0    

    r_OO_mag = dsqrt( dot_product( r_OO , r_OO ) )
    q2_mag   = dot_product( q , q )
    q_mag = dsqrt ( q2_mag )

    Select Case(function_type)
    Case(1)
       !****************************** MS-EVB3 function type *************************

       Vconstij = function_parameters(1) ; gamma = function_parameters(2)
       P        = function_parameters(3) ; k     = function_parameters(4)
       D        = function_parameters(5) ; beta  = function_parameters(6)
       R0       = function_parameters(7) ; Pprime= function_parameters(8)
       alpha    = function_parameters(9) ; rl0   = function_parameters(10)

       ! geometric factor is given as a product of three terms
       fac1 = exp( -gamma * q2_mag )
       fac2 = 1d0 + P * exp( -k * ( r_OO_mag - D ) ** 2 )
       fac3 = 0.5d0 * ( 1d0 - tanh( beta * ( r_OO_mag - R0 ) ) ) + Pprime * exp( -alpha * ( r_OO_mag - rl0 ) )

       ! all derivatives are relative to the general distance coordinate ( either q, r_OO )
       dfac1 = -gamma * 2d0 * q_mag *  exp( -gamma * q2_mag )
       dfac2 = P * -k * 2d0 * ( r_OO_mag - D ) * exp( -k * ( r_OO_mag - D ) ** 2 )
       dfac3 = -0.5d0 * beta / cosh( beta * ( r_OO_mag - R0 ) )**2 - Pprime * alpha * exp( -alpha * ( r_OO_mag - rl0 ) ) 

       ! geometric factor
       A = fac1 * fac2 * fac3

       ! derivative w.r.t  O_donor coordinates
       dA(:,1) = dfac1 * fac2 * fac3 * 0.5d0 * q(:) / q_mag
       dA(:,1) = dA(:,1) + fac1 * dfac2 * fac3 * r_OO(:) / r_OO_mag
       dA(:,1) = dA(:,1) + fac1 * fac2 * dfac3 * r_OO(:) / r_OO_mag

       ! derivative w.r.t. O_acceptor coordinates
       dA(:,2) = dfac1 * fac2 * fac3 * 0.5d0 * q(:) / q_mag
       dA(:,2) = dA(:,2) + fac1 * dfac2 * fac3 * -r_OO(:) / r_OO_mag
       dA(:,2) = dA(:,2) + fac1 * fac2 * dfac3 * -r_OO(:) / r_OO_mag  

       ! derivative w.r.t. H zundel coordinates
       dA(:,3) = dfac1 * fac2 * fac3 * -q(:) / q_mag


    Case(2)
       !************************** product of 2 gaussians *******************
       Vconstij = function_parameters(1)
       gamma = function_parameters(2)
       k = function_parameters(3)
       D = function_parameters(4)

       ! gaussian 1
       fac1 = exp( -gamma * q2_mag )
       ! gaussian 2
       fac2 = exp( -k * ( r_OO_mag - D ) ** 2 )

       ! all derivatives are relative to the general distance coordinate ( either q, r_OO )
       dfac1 = -gamma * 2d0 * q_mag *  exp( -gamma * q2_mag )
       dfac2 = -k * 2d0 * ( r_OO_mag - D ) * exp( -k * ( r_OO_mag - D ) ** 2 )

       ! geometric factor
       A = fac1 * fac2

       ! derivative w.r.t  O_donor coordinates
       dA(:,1) = dfac1 * fac2 * 0.5d0 * q(:) / q_mag
       dA(:,1) = dA(:,1) + fac1 * dfac2 * r_OO(:) / r_OO_mag

       ! derivative w.r.t. O_acceptor coordinates
       dA(:,2) = dfac1 * fac2 * 0.5d0 * q(:) / q_mag
       dA(:,2) = dA(:,2) + fac1 * dfac2 * -r_OO(:) / r_OO_mag

       ! derivative w.r.t. H zundel coordinates
       dA(:,3) = dfac1 * fac2 * -q(:) / q_mag

    End Select


  end subroutine evb_diabatic_coupling_function






  !**********************************************
  ! this is the electrostatic part of the diabatic coupling
  ! it is composed of Coulombic interactions between
  ! the solvent water molecules and the H5O2 proton
  ! transfer transition state, with special charges
  ! for the H5O2 complex.  This has been generalized to
  ! other reactive molecules.
  !**********************************************
  subroutine evb_diabatic_coupling_electrostatics( Vex, dVex, i_mole_donor, i_mole_acceptor, system_data_diabat, atom_data_diabat, molecule_data_diabat )
    real*8, intent(out) :: Vex
    real*8, intent(out), dimension(:,:) :: dVex
    integer, intent(in) :: i_mole_donor, i_mole_acceptor
    type(system_data_type), intent(in) :: system_data_diabat
    type(atom_data_diabat_type), intent(in)   :: atom_data_diabat
    type(molecule_data_type), dimension(:), intent(in) :: molecule_data_diabat

    !******** this is a local data structure with pointers that will be set
    ! to subarrays of atom_data arrays for the specific atoms in the molecule
    type(single_molecule_data_type) :: single_molecule_data_donor, single_molecule_data_acceptor

    integer, dimension(:), allocatable :: screen_list, molecule_list
    integer :: j_mole, i_atom, j_atom,  i_type, i_mole_type, j_mole_type
    real*8 :: q_i , q_j, r_mag , q_exchange_transfer
    real*8, dimension(3) :: shift, shiftd, shifta, r_ij, dV_ij, r_com_zundel, xyz_donor, xyz_acceptor, dr(3)

    Vex=0d0
    dVex=0d0


    ! set pointers for this data structure to target molecule
    ! note we set force pointers of single_molecule_data arrays to point to
    ! molecule block of dVex, as this is what we want to update

    call return_molecule_block( single_molecule_data_donor , molecule_data_diabat(i_mole_donor)%n_atom, molecule_data_diabat(i_mole_donor)%atom_index, atom_xyz=atom_data_diabat%xyz, atom_force=dVex, atom_mass=atom_data_diabat%mass, atom_charge=atom_data_diabat%charge, atom_type_index=atom_data_diabat%atom_type_index )
    call return_molecule_block( single_molecule_data_acceptor , molecule_data_diabat(i_mole_acceptor)%n_atom, molecule_data_diabat(i_mole_acceptor)%atom_index, atom_xyz=atom_data_diabat%xyz, atom_force=dVex, atom_mass=atom_data_diabat%mass, atom_charge=atom_data_diabat%charge, atom_type_index=atom_data_diabat%atom_type_index )

    ! get Zundel center of mass, note we use this later for pbc shifts.  Note that the donor and acceptor may be split over
    ! periodic boundary conditions, and so there may be a shift from the center of mass from the zundel to either the donor
    ! or acceptor.  these shifts are stored as shifta and shiftd

    ! therefore, to get the positions of the donor and acceptor at the zundel r_com location, use
    ! pbc_dr( r_com_zundel , single_molecule_data_donor%xyz(:,i_atom), shiftd ) and 
    ! pbc_dr( r_com_zundel , single_molecule_data_acceptor%xyz(:,i_atom), shifta )

    r_com_zundel = zundel_r_com( i_mole_donor, i_mole_acceptor, system_data_diabat, molecule_data_diabat , single_molecule_data_donor%mass, single_molecule_data_acceptor%mass, shiftd, shifta ) 

    ! create_screening list for donor/acceptor atoms which we don't include in interactions
    allocate(molecule_list(2))
    molecule_list(1) = i_mole_donor; molecule_list(2) = i_mole_acceptor
    call create_atom_screen_list( screen_list , molecule_list , molecule_data_diabat )
    deallocate(molecule_list)

    ! this is a sum over coulomb interactions of all water molecules with the 7 atoms of the 
    ! H5O2+ Zundel complex

    ! first calculate interaction with donor
    ! donor is currently in its basic topology, as proton has been transferred
    do i_atom=1, molecule_data_diabat(i_mole_donor)%n_atom
       ! figure out charge
       i_type = single_molecule_data_donor%atom_type_index(i_atom)
       q_i = evb_exchange_charge_atomic( i_type )

       ! note in this subroutine, PBC boundary conditions are taken with respect to the whole Zundel complex, rather
       ! than the donor and acceptor water molecules individually.  This is so that this energy is invarient to switching the
       ! donor and acceptor, which may not be the case otherwise, since the centers of mass of the molecules would change, and therefore
       ! there's no guarantee that the same minimum images would be taken in both cases

       ! get donor coordinate relative to zundel
       dr(:)  =   pbc_dr( r_com_zundel(:) , single_molecule_data_donor%xyz(:,i_atom), shiftd )
       xyz_donor = r_com_zundel + dr

       do j_atom = 1, system_data_diabat%total_atoms  ! loop over all solvent atoms
       ! screen out donor/acceptor interactions which we don't include
          if ( screen_list(j_atom) == 1 ) then

             shift = pbc_shift( r_com_zundel(:) , atom_data_diabat%xyz(:,j_atom) , system_data_diabat%box , system_data_diabat%xyz_to_box_transform )
             q_j = atom_data_diabat%charge( j_atom )
             r_ij = -pbc_dr( xyz_donor(:) , atom_data_diabat%xyz(:,j_atom) , shift )
             r_mag = dsqrt( dot_product( r_ij , r_ij ) )

             Vex = Vex + q_i * q_j / r_mag * constants%conv_e2A_kJmol     ! convert from e^2/A to kJ/mol
             dV_ij = - q_i * q_j / r_mag ** 3 * r_ij * constants%conv_e2A_kJmol     ! convert from e^2/A to kJ/mol
                
             ! here we need to use pointer to reference molecule block of atom array.  This was set earlier..
             single_molecule_data_donor%force(:,i_atom) = single_molecule_data_donor%force(:,i_atom) + dV_ij(:)
             ! here we can directly access the atomic data structure
             dVex(:, j_atom) = dVex(:, j_atom) - dV_ij
          endif
       enddo
    enddo

    ! the transfering proton exchange charge depends on the identity of the donor and acceptor molecules
    i_mole_type = molecule_data_diabat(i_mole_acceptor)%molecule_type_index
    j_mole_type = molecule_data_diabat(i_mole_donor)%molecule_type_index
    q_exchange_transfer = evb_exchange_charge_proton( i_mole_type, j_mole_type )


    ! now calculate interaction with acceptor
    do i_atom=1, molecule_data_diabat(i_mole_acceptor)%n_atom
       ! figure out charge, last atom in acceptor is transferring proton
       i_type = molecule_data_diabat(i_mole_acceptor)%atom_index( i_atom )
       ! note we're using acid topology for acceptor, because we've already transferred the proton
       if ( i_atom == molecule_data_diabat(i_mole_acceptor)%n_atom ) then
          ! transferring proton
          q_i = q_exchange_transfer
       else
          q_i = evb_exchange_charge_atomic( i_type )
       end if

       ! see above comment about use of zundel center of mass for PBC shift

       ! get acceptor coordinate relative to zundel
       dr(:)  =   pbc_dr( r_com_zundel(:) , single_molecule_data_acceptor%xyz(:,i_atom), shifta )
       xyz_acceptor = r_com_zundel + dr

       do j_atom = 1, system_data_diabat%total_atoms  ! loop over all solvent atoms
       ! screen out donor/acceptor interactions which we don't include
          if ( screen_list(j_atom) == 1 ) then

             shift = pbc_shift( r_com_zundel(:) , atom_data_diabat%xyz(:,j_atom), system_data_diabat%box , system_data_diabat%xyz_to_box_transform )
             q_j = atom_data_diabat%charge( j_atom )
             r_ij = -pbc_dr( xyz_acceptor(:) , atom_data_diabat%xyz(:,j_atom) , shift )
             r_mag = dsqrt( dot_product( r_ij , r_ij ) )

             Vex = Vex + q_i * q_j / r_mag * constants%conv_e2A_kJmol     !convert from e^2/A to kJ/mol
             dV_ij = - q_i * q_j / r_mag ** 3 * r_ij * constants%conv_e2A_kJmol ! convert from e^2/A to kJ/mol

             ! here we need to use pointer to reference molecule block of atom array.  This was set earlier..
             single_molecule_data_acceptor%force(:,i_atom) = single_molecule_data_acceptor%force(:,i_atom) + dV_ij(:)
             ! here we can directly access the atomic data structure
             dVex(:, j_atom) = dVex(:, j_atom) - dV_ij

          endif
       enddo
    enddo



  end subroutine evb_diabatic_coupling_electrostatics






  !*********************************************
  ! this subroutine calculates the total energy and force for a 
  ! particular diabat.  Note, a call to this subroutine essentially
  ! takes the place of a call to calculate_total_force_energy subroutine,
  ! but this calculation is much cheaper, because it only calculates the
  ! updated force and energy terms from the last proton transfer, thus
  ! requiring order(N) real-space operations
  !
  ! data structures:
  !    system_data_diabat
  !    atom_data_diabat
  !    molecule_data_diabat
  !
  !  are local data structures for each diabat topology, and contain all necessary
  !  data to compute energy and force of system
  !*********************************************
  subroutine ms_evb_diabat_force_energy( system_data_diabat, i_diabat, i_mole_principle, atom_data_diabat, molecule_data_diabat, PME_data, hydronium_molecule_index_diabat)
    use pme_routines
    integer, intent(in) :: i_diabat, i_mole_principle
    type(system_data_type), intent(inout) :: system_data_diabat
    type(atom_data_diabat_type), intent(inout) :: atom_data_diabat
    type(molecule_data_type), dimension(:),intent(inout) :: molecule_data_diabat
    type(PME_data_type), intent(in)            :: PME_data
    integer, intent(inout) , dimension(:) :: hydronium_molecule_index_diabat

    !******** this is a local data structure with pointers that will be set
    ! to subarrays of atom_data arrays for the specific atoms in the molecule
    type(single_molecule_data_type) :: single_molecule_data_donor, single_molecule_data_acceptor

    ! this is local Q_grid for diabat
    real*8, dimension(:,:,:),allocatable :: Q_grid_local

    real*8, dimension(:,:), allocatable :: xyz_scale, d_force_atoms
    integer :: total_atoms

    integer :: lookup_index, i_hop,i_mole, i_atom, i_mole_acceptor, i_atom_acceptor, i_mole_donor, i_atom_donor, i_heavy_donor, i_heavy_acceptor, i_proton_index
    real*8  :: dE_donor_diabat_intra, dE_donor_diabat_real_space, dE_acceptor_diabat_intra, dE_acceptor_diabat_real_space, E_acceptor_ms_evb_repulsion, E_donor_ms_evb_repulsion
    real*8  :: E_reference_donor, E_reference_acceptor
    real*8,dimension(3,3) :: kk

    allocate(Q_grid_local(PME_data%pme_grid,PME_data%pme_grid,PME_data%pme_grid))
    ! set Q_grid_local to Q_grid global variable which should have Q_grid stored for the principle diabat
    Q_grid_local = Q_grid

    total_atoms = system_data_diabat%total_atoms
    allocate( d_force_atoms(3,total_atoms) )

    !******************************** loop over proton hops relative to principle diabat.  for each proton hop, reorganize data structures in the new acceptor topology, and calculate
    !******************************** energy differences between the donor and acceptor topology

    ! first donor is always principle hydronium
    ! set i_mole_acceptor equal to this, because
    ! i_mole_donor will be copied from previous i_mole_acceptor
    i_mole_acceptor = i_mole_principle

    ! loop through all proton hops.  A negative integer signals the end of the chain
    do i_hop=1, evb_max_chain
       if ( evb_diabat_proton_log( i_diabat , i_hop, 1 ) < 0 ) exit
       i_mole_donor = i_mole_acceptor
       i_atom_donor = evb_diabat_proton_log( i_diabat , i_hop, 2 )
       i_mole_acceptor = evb_diabat_proton_log( i_diabat , i_hop, 4 )
       i_heavy_acceptor = evb_diabat_proton_log( i_diabat , i_hop, 5 )
       i_atom_acceptor = n_atom_temp(i_mole_acceptor) + 1


       ! **********************   forces and energy updates from donor topology
       d_force_atoms=0d0

       ! first get energy and forces from donor and acceptor intra-molecular bond,angle, etc interactions in donor topology
       call ms_evb_diabat_force_energy_update_intra( d_force_atoms, dE_donor_diabat_intra, i_mole_donor, i_mole_acceptor, system_data_diabat, atom_data_diabat, molecule_data_diabat)
       ! real-space non-bonded energy and forces in donor topology--note reciprocal space PME will be computed later
       call ms_evb_diabat_force_energy_update_real_space( d_force_atoms, dE_donor_diabat_real_space, i_mole_donor, i_mole_acceptor, system_data_diabat, atom_data_diabat, molecule_data_diabat, PME_data)
       ! special ms-evb repulsion terms
       call ms_evb_intermolecular_repulsion( d_force_atoms , E_donor_ms_evb_repulsion, system_data_diabat, atom_data_diabat, molecule_data_diabat, hydronium_molecule_index_diabat )
       ! reference chemical energy for this adiabatic state
       call get_adiabatic_reference_energy( E_reference_donor, i_mole_donor, molecule_data(i_mole_donor)%molecule_type_index )


       ! now, while we still have data structures in the donor diabat topology, subtract the donor and acceptor contribution
       ! the updated Q_grids will be stored for use in the separate reciprocal space part
       call return_molecule_block( single_molecule_data_donor , molecule_data_diabat(i_mole_donor)%n_atom, molecule_data_diabat(i_mole_donor)%atom_index, atom_xyz=atom_data_diabat%xyz, atom_charge=atom_data_diabat%charge )
       call return_molecule_block( single_molecule_data_acceptor , molecule_data_diabat(i_mole_acceptor)%n_atom, molecule_data_diabat(i_mole_acceptor)%atom_index, atom_xyz=atom_data_diabat%xyz, atom_charge=atom_data_diabat%charge )

       call construct_reciprocal_lattice_vector(kk,system_data_diabat%box)
       !********** donor
       allocate(xyz_scale(3,molecule_data_diabat(i_mole_donor)%n_atom)   
       call create_scaled_direct_coordinates( xyz_scale, single_molecule_data_donor%xyz, molecule_data_diabat(i_mole_donor)%n_atom, kk, PME_data%pme_grid )
       ! the -1 signals to subtract these contributions rather than add them
       call modify_Q_grid( Q_grid_local, single_molecule_data_donor%charge, xyz_scale, molecule_data_diabat(i_mole_donor)%n_atom, PME_data%pme_grid, PME_data%spline_order, PME_data%spline_grid, -1 )  ! donor
       deallocate(xyz_scale)

       !********** acceptor
       allocate(xyz_scale(3,molecule_data_diabat(i_mole_acceptor)%n_atom)
       call create_scaled_direct_coordinates( xyz_scale, single_molecule_data_acceptor%xyz, molecule_data_diabat(i_mole_acceptor)%n_atom, kk,PME_data%pme_grid )
       ! the -1 signals to subtract these contributions rather than add them
       call modify_Q_grid( Q_grid_local, single_molecule_data_acceptor%charge, xyz_scale, molecule_data_diabat(i_mole_acceptor)%n_atom, PME_data%pme_grid, PME_data%spline_order, PME_data%spline_grid, -1 ) ! acceptor
       deallocate(xyz_scale)

       !************************************* subtract the forces from the donor topology
       atom_data_diabat%force = atom_data_diabat%force - d_force_atoms


       ! ********************** Now switch to acceptor topology, and compute force and energy updates from acceptor topology
       d_force_atoms=0d0
       ! Change data structures to acceptor topology 
       call evb_change_data_structures_proton_transfer( i_mole_donor, i_atom_donor, i_mole_acceptor, i_atom_acceptor, i_heavy_acceptor, atom_data_diabat, molecule_data_diabat, system_data_diabat, hydronium_molecule_index_diabat )

       ! get energy and forces from donor and acceptor intra-molecular bond,angle, etc interactions in acceptor topology
       call ms_evb_diabat_force_energy_update_intra( d_force_atoms, dE_acceptor_diabat_intra, i_mole_donor, i_mole_acceptor, system_data_diabat, atom_data_diabat, molecule_data_diabat)
       ! real-space non-bonded energy and forces in acceptor topology--note reciprocal space PME will be computed later
       call ms_evb_diabat_force_energy_update_real_space( d_force_atoms, dE_acceptor_diabat_real_space, i_mole_donor, i_mole_acceptor, system_data_diabat, atom_data_diabat, molecule_data_diabat, PME_data)
       ! special ms-evb repulsion terms
       call ms_evb_intermolecular_repulsion( d_force_atoms , E_acceptor_ms_evb_repulsion, system_data_diabat, atom_data_diabat, molecule_data_diabat, hydronium_molecule_index_diabat )
       ! reference chemical energy for this adiabatic state
       call get_adiabatic_reference_energy( E_reference_acceptor, i_mole_acceptor, molecule_data(i_mole_acceptor)%molecule_type_index )


       !********** donor
       allocate(xyz_scale(3,molecule_data_diabat(i_mole_donor)%n_atom)
       call create_scaled_direct_coordinates( xyz_scale, single_molecule_data_donor%xyz, molecule_data_diabat(i_mole_donor)%n_atom, kk, PME_data%pme_grid )
       ! the +1 signals to add these contributions
       call modify_Q_grid( Q_grid_local, single_molecule_data_donor%charge, xyz_scale, molecule_data_diabat(i_mole_donor)%n_atom, PME_data%pme_grid, PME_data%spline_order, PME_data%spline_grid, 1 )  ! donor
       deallocate(xyz_scale)

       !********** acceptor
       allocate(xyz_scale(3,molecule_data_diabat(i_mole_acceptor)%n_atom)
       call create_scaled_direct_coordinates( xyz_scale, single_molecule_data_acceptor%xyz, molecule_data_diabat(i_mole_acceptor)%n_atom, kk,PME_data%pme_grid )
       ! the +1 signals to add these contributions
       call modify_Q_grid( Q_grid_local, single_molecule_data_acceptor%charge, xyz_scale, molecule_data_diabat(i_mole_acceptor)%n_atom, PME_data%pme_grid, PME_data%spline_order, PME_data%spline_grid, 1 ) ! acceptor
       deallocate(xyz_scale)

      !************************************* add the forces from the acceptor topology
       atom_data_diabat%force = atom_data_diabat%force + d_force_atoms

       ! output energy
       system_data_diabat%potential_energy = system_data_diabat%potential_energy + E_reference_acceptor + dE_acceptor_diabat_intra + dE_acceptor_diabat_real_space + E_acceptor_ms_evb_repulsion - E_reference_donor - dE_donor_diabat_intra - dE_donor_diabat_real_space - E_donor_ms_evb_repulsion

    end do   ! end loop over proton hops

    ! store this Q_grid for future use
    Q_grid_diabats(:,:,:,i_diabat) = Q_grid_local
    Q_grid_filled(i_diabat)=1
    deallocate(Q_grid_local)



  end subroutine ms_evb_diabat_force_energy




  !**********************************
  ! this subroutine calculates the real-space force and energy
  ! for molecules i_mole_donor and i_mole_acceptor interacting
  ! with all the other molecules in the system
  !**********************************
  subroutine ms_evb_diabat_force_energy_update_real_space( force_atoms, dE_real_space, i_mole_donor, i_mole_acceptor, system_data_diabat, atom_data_diabat, molecule_data_diabat, PME_data )
    use pair_int_real_space
    real*8, dimension(:,:),intent(inout) :: force_atoms
    real*8, intent(out)                  :: dE_real_space
    integer, intent(in) :: i_mole_donor, i_mole_acceptor
    type(system_data_type), intent(in)                 :: system_data_diabat
    type(atom_data_diabat_type), intent(in)            :: atom_data_diabat
    type(molecule_data_type), dimension(:), intent(in) :: molecule_data_diabat    
    type(PME_data_type), intent(in)                    :: PME_data

    !******** this is a local data structure with pointers that will be set
    ! to subarrays of atom_data arrays for the specific atoms in the molecule
    type(single_molecule_data_type) :: single_molecule_data_donor, single_molecule_data_acceptor

    ! this is data structure used in pairwise interactions
    type(pairwise_neighbor_data_type) :: pairwise_atom_data, pairwise_neighbor_data_cutoff

    integer , dimension(:) , allocatable    :: cutoff_mask
    real*8  :: real_space_cutoff2, alpha_sqrt, erf_factor, E_elec_local, E_vdw_local
    integer :: i_atom, j_atom, i_index, j_index, total_atoms, n_cutoff, size_lj_parameter
    integer :: i , j
    integer, dimension(:), allocatable :: screen_list, molecule_list 

    !************** NOTE **********************
    ! we are using both pointers as well as explicit calls to the force_atoms
    ! array to update forces.  For inter-molecular interactions, we explicitly
    ! update the force_atoms array, but for intra-molecular interactions,
    ! we need to reference a molecular data structure, so we create pointers
    ! in the single_molecule_data structures to the force_atoms array, and
    ! utilize these pointers instead
    !******************************************


    ! set pointers for donor and acceptor molecules to atom_data data structure
    ! note here we set single_molecule_data pointers to local force array
    ! 'force_atoms', as we don't want to update the force in the atom_data arrays

    call return_molecule_block( single_molecule_data_donor , molecule_data_diabat(i_mole_donor)%n_atom, molecule_data_diabat(i_mole_donor)%atom_index, atom_xyz=atom_data_diabat%xyz, atom_force=force_atoms, atom_charge=atom_data_diabat%charge, atom_type_index=atom_data_diabat%atom_type_index )
    call return_molecule_block( single_molecule_data_acceptor , molecule_data_diabat(i_mole_acceptor)%n_atom, molecule_data_diabat(i_mole_acceptor)%atom_index, atom_xyz=atom_data_diabat%xyz, atom_force=force_atoms, atom_charge=atom_data_diabat%charge, atom_type_index=atom_data_diabat%atom_type_index )

    ! define local variables for convenience
    total_atoms = system_data_diabat%total_atoms
    alpha_sqrt  = PME_data%alpha_sqrt
    erf_factor  = PME_data%erf_factor

    ! this is (maximum) number of parameters that we need for VDWs interaction
    size_lj_parameter=size(atype_lj_parameter(1,1,:))

    ! this is for screening out donor-donor, and acceptor-acceptor intra-molecular interactions
    allocate(screen_list(total_atoms))

    ! allocate pairwise atom data for interactions
    call allocate_pairwise_neighbor_data( pairwise_atom_data , total_atoms , size_lj_parameter )

    allocate( cutoff_mask(total_atoms) )

    dE_real_space = 0d0
    real_space_cutoff2 = real_space_cutoff ** 2

    ! ***************************************
    ! ************************ Donor molecule first
    ! ***************************************
    ! loop over atoms in donor
    do i_index = 1, molecule_data_diabat(i_mole_donor)%n_atom
       i_atom  =  molecule_data_diabat(i_mole_donor)%atom_index(i_index)

       ! to vectorize, here we loop over all atoms, including atoms within the
       ! molecule.  We will adjust for exclusions before computing forces             
       do j_atom=1,total_atoms
             ! compute pair displacement
             pairwise_atom_data%dr(:,j_atom) = single_molecule_data_donor%xyz(:,i_index) - atom_data_diabat%xyz(:,j_atom)
       enddo
       ! now minimum image
       do j_atom=1,total_atoms     
             ! shift for orthorhombic box
             do i=1,3
                pairwise_atom_data%dr(i,j_atom) = pairwise_atom_data%dr(i,j_atom) -  system_data_diabat%box(i,i) * floor( pairwise_atom_data%dr(i,j_atom) / system_data_diabat%box(i,i) + 0.5d0 )
             end do
             pairwise_atom_data%dr2(j_atom) = dot_product( pairwise_atom_data%dr(:,j_atom), pairwise_atom_data%dr(:,j_atom) )
       enddo

       ! screen out donor atoms
       allocate(molecule_list(1))
       molecule_list(1) = i_mole_donor
       call create_atom_screen_list( screen_list , molecule_list , molecule_data_diabat ) 
       deallocate(molecule_list)

       ! now check cutoff
       cutoff_mask=0
       n_cutoff=0
       do j_atom = 1 , total_atoms
          if ( screen_list(j_atom) == 1 ) then
             if ( pairwise_atom_data%dr2(j_atom) < real_space_cutoff2 ) then
                cutoff_mask(j_atom) = 1
                n_cutoff = n_cutoff + 1
             endif
          endif
       enddo
       ! now allocate datastructure for atoms within cutoff distance
       ! allocate data structure to store data for these neighbors
       call allocate_pairwise_neighbor_data( pairwise_neighbor_data_cutoff , n_cutoff, size_lj_parameter )
       ! transfer data from pairwise_neighbor_data_verlet to pairwise_neighbor_data_cutoff using cutoff_mask
       call apply_cutoff_mask( total_atoms, pairwise_neighbor_data_cutoff, pairwise_atom_data, cutoff_mask )
       ! deallocate old arrays that we don't need anymore
       call deallocate_pairwise_neighbor_data( pairwise_atom_data )

       ! now calculate pairwise forces and energies for this atom and its neighbors
       ! all of these subroutines should vectorize...
       call pairwise_real_space_ewald( E_elec_local , pairwise_neighbor_data_cutoff%f_ij ,  pairwise_neighbor_data_cutoff%dr, pairwise_neighbor_data_cutoff%dr2,  pairwise_neighbor_data_cutoff%qi_qj, erf_factor , alpha_sqrt, PME_data%erfc_table , PME_data%erfc_grid , PME_data%erfc_max, constants%conv_e2A_kJmol )
       call pairwise_real_space_LJ( E_vdw_local , pairwise_neighbor_data_cutoff%f_ij ,  pairwise_neighbor_data_cutoff%dr, pairwise_neighbor_data_cutoff%dr2 , pairwise_neighbor_data_cutoff%atype_lj_parameter )

       dE_real_space = dE_real_space + E_elec_local + E_vdw_local

       ! running addition of forces...
       do j_index=1, n_cutoff
             j_atom = pairwise_neighbor_data_cutoff%atom_index(j_index)
             force_atoms(:,i_atom) =  force_atoms(:,i_atom) + pairwise_neighbor_data_cutoff%f_ij(:,j_index)
             force_atoms(:,j_atom) =  force_atoms(:,j_atom) - pairwise_neighbor_data_cutoff%f_ij(:,j_index)
       enddo

       ! deallocate old arrays that we don't need anymore
       call deallocate_pairwise_neighbor_data( pairwise_neighbor_data_cutoff)

    enddo  ! end loop over donor atoms



    ! ***************************************
    ! ************************ now Acceptor molecule
    ! ***************************************
    ! loop over atoms in acceptor
    do i_index = 1, molecule_data_diabat(i_mole_acceptor)%n_atom
       i_atom  =  molecule_data_diabat(i_mole_acceptor)%atom_index(i_index)

       ! to vectorize, here we loop over all atoms, including atoms within the
       ! molecule.  We will adjust for exclusions before computing forces
       do j_atom=1,total_atoms
             ! compute pair displacement
             pairwise_atom_data%dr(:,j_atom) = single_molecule_data_acceptor%xyz(:,i_index) - atom_data_diabat%xyz(:,j_atom)
       enddo
       ! now minimum image
       do j_atom=1,total_atoms
             ! shift for orthorhombic box
             do i=1,3
                pairwise_atom_data%dr(i,j_atom) = pairwise_atom_data%dr(i,j_atom) -  system_data_diabat%box(i,i) * floor( pairwise_atom_data%dr(i,j_atom) / system_data_diabat%box(i,i) + 0.5d0 )
             end do
             pairwise_atom_data%dr2(j_atom) = dot_product( pairwise_atom_data%dr(:,j_atom), pairwise_atom_data%dr(:,j_atom) )
       enddo


      ! screen out acceptor atoms,
      ! also don't count acceptor-donor interaction twice, it has been counted above
       allocate(molecule_list(2))
       molecule_list(1) = i_mole_donor ; molecule_list(2) = i_mole_acceptor
       call create_atom_screen_list( screen_list , molecule_list , molecule_data_diabat )
       deallocate(molecule_list)


       ! now check cutoff
       cutoff_mask=0
       n_cutoff=0
       do j_atom = 1 , total_atoms
          if ( screen_list( j_atom ) == 1 ) then
             if ( pairwise_atom_data%dr2(j_atom) < real_space_cutoff2 ) then
                cutoff_mask(j_atom) = 1
                n_cutoff = n_cutoff + 1
             endif
          end if
       enddo
       ! now allocate datastructure for atoms within cutoff distance
       ! allocate data structure to store data for these neighbors
       call allocate_pairwise_neighbor_data( pairwise_neighbor_data_cutoff , n_cutoff, size_lj_parameter )
       ! transfer data from pairwise_neighbor_data_verlet to
       ! pairwise_neighbor_data_cutoff using cutoff_mask
       call apply_cutoff_mask( total_atoms, pairwise_neighbor_data_cutoff, pairwise_atom_data, cutoff_mask )
       ! deallocate old arrays that we don't need anymore
       call deallocate_pairwise_neighbor_data( pairwise_atom_data )

       ! now calculate pairwise forces and energies for this atom and its neighbors
       ! all of these subroutines should vectorize...
       call pairwise_real_space_ewald( E_elec_local , pairwise_neighbor_data_cutoff%f_ij ,  pairwise_neighbor_data_cutoff%dr, pairwise_neighbor_data_cutoff%dr2,  pairwise_neighbor_data_cutoff%qi_qj, erf_factor , alpha_sqrt, PME_data%erfc_table , PME_data%erfc_grid , PME_data%erfc_max, constants%conv_e2A_kJmol )
       call pairwise_real_space_LJ( E_vdw_local , pairwise_neighbor_data_cutoff%f_ij ,  pairwise_neighbor_data_cutoff%dr, pairwise_neighbor_data_cutoff%dr2 , pairwise_neighbor_data_cutoff%atype_lj_parameter )

       dE_real_space = dE_real_space + E_elec_local + E_vdw_local

       ! running addition of forces...
       do j_index=1, n_cutoff
             j_atom = pairwise_neighbor_data_cutoff%atom_index(j_index)
             force_atoms(:,i_atom) =  force_atoms(:,i_atom) + pairwise_neighbor_data_cutoff%f_ij(:,j_index)
             force_atoms(:,j_atom) =  force_atoms(:,j_atom) - pairwise_neighbor_data_cutoff%f_ij(:,j_index)
       enddo

       ! deallocate old arrays that we don't need anymore
       call deallocate_pairwise_neighbor_data( pairwise_neighbor_data_cutoff)

    enddo  ! end loop over acceptor atoms



    !**********************************  Now intra-molecular interactions *******************

    !************* donor
    E_elec_local = 0d0; E_vdw_local = 0d0
    call intra_molecular_pairwise_energy_force( single_molecule_data_donor%force, E_elec_local, E_vdw_local, system_data_diabat, single_molecule_data_donor ,  molecule_data_diabat(i_mole_donor)%molecule_type_index , molecule_data(i_mole_donor)%n_atom, PME_data )
    dE_real_space = dE_real_space + E_elec_local + E_vdw_local

    !************ acceptor
    E_elec_local = 0d0; E_vdw_local = 0d0
    call intra_molecular_pairwise_energy_force( single_molecule_data_acceptor%force, E_elec_local, E_vdw_local, system_data_diabat, single_molecule_data_acceptor , molecule_data_diabat(i_mole_acceptor)%molecule_type_index , molecule_data(i_mole_acceptor)%n_atom, PME_data )
    dE_real_space = dE_real_space + E_elec_local + E_vdw_local



  end subroutine ms_evb_diabat_force_energy_update_real_space





  !**********************************
  ! this subroutine calculates the intra-molecular bond and angle force and energy
  ! for molecules i_mole_donor and i_mole_acceptor
  !**********************************
  subroutine ms_evb_diabat_force_energy_update_intra( force_atoms, E_intra, i_mole_donor, i_mole_acceptor, system_data_diabat, atom_data_diabat, molecule_data_diabat)
    use bonded_interactions
    real*8, dimension(:,:),intent(inout) :: force_atoms
    real*8, intent(out)                  :: E_intra
    integer, intent(in) :: i_mole_donor, i_mole_acceptor
    type(system_data_type) , intent(inout) :: system_data_diabat
    type(molecule_data_type), dimension(:), intent(inout) :: molecule_data_diabat
    type(atom_data_diabat_type), intent(inout) :: atom_data_diabat

    !******** this is a local data structure with pointers that will be set
    ! to subarrays of atom_data arrays for the specific atoms in the molecule
    type(single_molecule_data_type) :: single_molecule_data_donor, single_molecule_data_acceptor

    integer :: i_donor_type, i_acceptor_type
    real*8  ::  E_local

    E_intra=0d0

    ! this is used to look up bonded force field parameters for this molecule
    i_donor_type = molecule_data_diabat(i_mole_donor)%molecule_type_index
    i_acceptor_type = molecule_data_diabat(i_mole_acceptor)%molecule_type_index

    ! set pointers for single_molecule_data to target molecule subarrays of atom_data
    ! set force pointers in single_molecule_data arrays to point to local force_atoms array
    call return_molecule_block( single_molecule_data_donor , molecule_data_diabat(i_mole_donor)%n_atom, molecule_data_diabat(i_mole_donor)%atom_index, atom_xyz=atom_data_diabat%xyz, atom_force=force_atoms, atom_type_index=atom_data_diabat%atom_type_index )
    call return_molecule_block( single_molecule_data_acceptor , molecule_data_diabat(i_mole_acceptor)%n_atom, molecule_data_diabat(i_mole_acceptor)%atom_index,  atom_xyz=atom_data_diabat%xyz, atom_force=force_atoms, atom_type_index=atom_data_diabat%atom_type_index )

    !****************** bond terms ***********************************
    ! donor
    call intra_molecular_bond_energy_force( E_local, single_molecule_data_donor, i_donor_type,  molecule_data_diabat(i_mole_donor)%n_atom )
    E_intra = E_intra + E_local
    ! acceptor
    call intra_molecular_bond_energy_force( E_local, single_molecule_data_acceptor, i_acceptor_type,  molecule_data_diabat(i_mole_acceptor)%n_atom  )
    E_intra = E_intra + E_local

    !******************** angle terms ****************************
    ! donor
    call intra_molecular_angle_energy_force( E_local, single_molecule_data_donor, i_donor_type,  molecule_data_diabat(i_mole_donor)%n_atom  )
    E_intra = E_intra + E_local
    ! acceptor
    call intra_molecular_angle_energy_force( E_local, single_molecule_data_acceptor, i_acceptor_type,  molecule_data_diabat(i_mole_acceptor)%n_atom  )
    E_intra = E_intra + E_local

    !******************* dihedral terms ***************************
    ! donor
    call intra_molecular_dihedral_energy_force( E_local, single_molecule_data_donor, i_donor_type,  molecule_data_diabat(i_mole_donor)%n_atom  ) 
    E_intra = E_intra + E_local
    ! acceptor
    call intra_molecular_dihedral_energy_force( E_local, single_molecule_data_acceptor, i_acceptor_type,  molecule_data_diabat(i_mole_acceptor)%n_atom  ) 
    E_intra = E_intra + E_local

  end subroutine ms_evb_diabat_force_energy_update_intra







  !************************************
  ! this subroutine calculates the reciprocal_space_pme
  ! energy and force for each diabat, using the dQ_dr stored grid from
  ! the principle diabat, along with the stored Q_grids from
  ! each diabat
  !***********************************
  subroutine calculate_reciprocal_space_pme( i_mole_principle, system_data , molecule_data, atom_data, PME_data  )
    use MKL_DFTI
    use omp_lib
    implicit none
    integer, intent(in) :: i_mole_principle
    type(system_data_type), intent(in) :: system_data
    type(molecule_data_type),dimension(:),intent(in) :: molecule_data
    type(atom_data_type), intent(in)   :: atom_data
    type(PME_data_type),  intent(in)   :: PME_data

    TYPE(DFTI_DESCRIPTOR),pointer :: dfti_desc_local,dfti_desc_inv_local
    integer :: index_store, index, i_mole, i_atom, i, j, i_diabat, K, status, length(3)
    real*8 :: kk(3,3), force_temp(3), E_recip_local, dE_recip
    real*8,dimension(:,:,:), allocatable :: Q_local, theta_conv_Q_local
    complex*16,dimension(:,:,:),allocatable::FQ
    real*8,dimension(:), allocatable::q_1r
    complex*16,dimension(:), allocatable::q_1d
    real*8,dimension(:,:), allocatable :: pme_force_recip_diabat
    integer :: split_do

    K=PME_data%pme_grid
    ! note three dimensions
    length=PME_data%pme_grid

    call construct_reciprocal_lattice_vector(kk, system_data%box)


    !*************************************************************
    !                             IMPORTANT:: 
    ! for the shared dfti descriptors that we use in the parallel code, below, we must set
    ! DFTI_NUMBER_OF_USER_THREADS to be equal to n_threads,
    ! as n_threads are using the same shared DFTI descriptors in the below parallel code. 
    ! If DFTI_NUMBER_OF_USER_THREADS is not set to the appropriate value, parallel code below will not work correctly
    !************************************************************

    status=DftiCreateDescriptor(dfti_desc_local, DFTI_DOUBLE, DFTI_COMPLEX, 3, length)
    status=DftiSetValue(dfti_desc_local, DFTI_NUMBER_OF_USER_THREADS, n_threads)
    status=DftiCommitDescriptor(dfti_desc_local)
    ! don't scale back transform because pick up factor of K^3 from convolution
    status=DftiCreateDescriptor(dfti_desc_inv_local, DFTI_DOUBLE, DFTI_COMPLEX, 3, length)
    status=DftiSetValue(dfti_desc_inv_local, DFTI_NUMBER_OF_USER_THREADS, n_threads)
    status = DftiCommitDescriptor(dfti_desc_inv_local)


    !***************************************************************
    ! Note this code below may look confusing, but essentially we are just
    ! doing reciprocal space pme, using stored data structures
    ! therefore, see pme_reciprocal_space_energy_force for comments on this algorithm
    !***************************************************************

    ! decide how to split the parallel section
    if (n_threads .eq. 1 ) then
       split_do = 1
    else
       ! integer arithmetic
       split_do = diabat_index/n_threads
    endif

    !**** parallelize over diabats
    call OMP_SET_NUM_THREADS(n_threads)
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(n_threads,split_do,theta_conv_Q, CB, dQ_dr,dQ_dr_index,Q_grid_diabats,system_data,atom_data,molecule_data,PME_data,K,kk,i_mole_principle, diabat_index,evb_forces_lookup_index,evb_hamiltonian,evb_forces_store,dfti_desc_local,dfti_desc_inv_local) 
    !$OMP CRITICAL
    allocate( FQ(K,K,K), Q_local(K,K,K), theta_conv_Q_local(K,K,K), q_1r(K**3), q_1d(K**3), pme_force_recip_diabat(3,system_data%total_atoms) )
    !$OMP END CRITICAL
    !$OMP DO SCHEDULE(dynamic, split_do)
    do i_diabat=2,diabat_index  ! skip the first diabat, as that is the principle diabat, and the force is correct

       ! Q grid for this diabat should already be stored
       Q_local(:,:,:) = Q_grid_diabats(:,:,:,i_diabat)
       q_1r=RESHAPE(Q_local, (/K**3/) )
       q_1d=cmplx(q_1r,0.,16)

       ! take FT
       status=DftiComputeForward(dfti_desc_local, q_1d)
       FQ=RESHAPE(q_1d, (/K,K,K/) )

       !multiply B*C*F(Q)
       FQ=FQ*CB

       ! take Finv
       q_1d=RESHAPE(FQ,(/K**3/) )
       status = DftiComputeBackward(dfti_desc_inv_local, q_1d)
       FQ=RESHAPE(q_1d, (/K,K,K/) )

       theta_conv_Q_local=dble(FQ)

       ! reciprocal space energy for this diabat
       E_recip_local = 0.5D0*sum((Q_local*theta_conv_Q_local))*constants%conv_e2A_kJmol
       ! subtract reciprocal space energy of principle diabat
       dE_recip =  E_recip_local - PME_data%E_recip 

       ! now force
       pme_force_recip_diabat=0d0

       do i_atom=1, system_data%total_atoms
          do j=1, size(dQ_dr(1,:,1))    
             pme_force_recip_diabat(:,i_atom) = pme_force_recip_diabat(:,i_atom) + dQ_dr(:,j,i_atom) * theta_conv_Q_local(dQ_dr_index(1,j,i_atom), dQ_dr_index(2,j,i_atom),dQ_dr_index(3,j,i_atom))
          enddo
       enddo

       ! dQ_dr is stored in scaled coordinates, convert to general coordinates
       do i_atom=1, system_data%total_atoms
          force_temp=0d0
          do i=1,3
             do j=1,3
                force_temp(i) = force_temp(i) - dble(K) * kk(j,i) * pme_force_recip_diabat(j,i_atom)
             enddo
          enddo
          pme_force_recip_diabat(:,i_atom) = force_temp(:)
       enddo

       ! now subtract off the reciprocal space electrostatic force of the principle diabat
       pme_force_recip_diabat(:,i_atom) = pme_force_recip_diabat(:,i_atom) - PME_data%force_recip


       ! now get contribution to force from changes in dQ_dr for this diabat
       ! note here data structures should be in principle diabat topology
       call update_reciprocal_space_force_dQ_dr( pme_force_recip_diabat, theta_conv_Q_local, i_diabat, i_mole_principle, system_data, atom_data, molecule_data, PME_data, kk )


       ! add the energy and force differences to the stored energy and force arrays
       index = evb_forces_lookup_index(i_diabat,i_diabat)
       evb_hamiltonian(i_diabat,i_diabat) = evb_hamiltonian(i_diabat,i_diabat) + dE_recip
       evb_forces_store(:,:,index) = evb_forces_store(:,:,index) + pme_force_recip_diabat(:, :)

    enddo

    !$OMP END DO NOWAIT
    deallocate( FQ, Q_local, q_1r, q_1d,theta_conv_Q_local, pme_force_recip_diabat )
    !$OMP END PARALLEL

    status=DftiFreeDescriptor(dfti_desc_local)
    status=DftiFreeDescriptor(dfti_desc_inv_local)

  end subroutine calculate_reciprocal_space_pme




  !*********************************************
  ! this subroutine updates the reciprocal space force
  ! for a diabat, by using the actual theta_conv_Q array
  ! for the diabat, and updating the changes
  ! ro the dQ_dr grid from the principle diabat
  !*********************************************
  subroutine update_reciprocal_space_force_dQ_dr( pme_force_recip_diabat, theta_conv_Q_local, i_diabat, i_mole_principle, system_data, atom_data, molecule_data, PME_data, kk )
    use pme_routines
    implicit none
    real*8,dimension(:,:),intent(inout) :: pme_force_recip_diabat
    real*8,dimension(:,:,:),intent(in)   :: theta_conv_Q_local
    integer, intent(in) :: i_diabat, i_mole_principle
    type(system_data_type), intent(in) :: system_data
    type(atom_data_type),   intent(in) :: atom_data
    type(molecule_data_type), dimension(:), intent(in) :: molecule_data
    real*8,dimension(:,:),intent(in) :: kk  

    ! these are temporary local arrays that we modify for local diabat topology
    ! note that atom_data_diabat_type does not contain pointers, which is what we want, unlike atom_data_type
    type(system_data_type)             :: system_data_diabat
    type(atom_data_diabat_type)               :: atom_data_diabat
    type(molecule_data_type), dimension(:), allocatable :: molecule_data_diabat

    !******** this is a local data structure with pointers that will be set
    ! to subarrays of atom_data arrays for the specific atoms in the molecule
    type(single_molecule_data_type) :: single_molecule_data_donor , single_molecule_data_acceptor

    real*8, dimension(:,:),allocatable :: xyz_scale, force
    integer, dimension(:), allocatable :: hydronium_molecule_index_temp
    integer :: i_mole, i_mole_donor, i_atom_donor, i_mole_acceptor, i_atom_acceptor, i_heavy_acceptor, i_hop, i_atom, i, count, flag_update=1


    allocate( molecule_data_diabat(system_data%n_mole) )
    allocate( hydronium_molecule_index_temp(size(hydronium_molecule_index)) )

    ! all of these data structures, and also pme_force_recip_diabat should initially be in principle diabat data structures
    system_data_diabat = system_data
    atom_data_diabat   = atom_data
    do i_mole=1, system_data%n_mole
       molecule_data_diabat(i_mole) =  molecule_data(i_mole)    
    enddo

    ! NOTE, here set atom_data_diabat%force = pme_force_recip_diabat, as we need
    ! to use the atom_data_diabat data structure when we change topology
    atom_data_diabat%force = pme_force_recip_diabat 

    hydronium_molecule_index_temp = hydronium_molecule_index


    ! ****** loop through proton hops starting from principle diabat
    ! first donor is always principle hydronium
    ! set i_mole_acceptor equal to this, because
    ! i_mole_donor will be copied from previous i_mole_acceptor
    i_mole_acceptor = i_mole_principle

    ! loop through all proton hops.  A negative integer signals the end of the chain
    loop1: do i_hop=1, evb_max_chain
       if ( evb_diabat_proton_log( i_diabat , i_hop, 1 ) < 0 ) exit loop1
       i_mole_donor = i_mole_acceptor
       i_atom_donor = evb_diabat_proton_log( i_diabat , i_hop, 2 )
       i_mole_acceptor = evb_diabat_proton_log( i_diabat , i_hop, 4 )
       i_heavy_acceptor = evb_diabat_proton_log( i_diabat , i_hop, 5 )
       i_atom_acceptor = molecule_data_diabat(i_mole_acceptor)%n_atom + 1

       !*************************************** subtract the contribution of the donor as an h3o+ molecule
       ! set pointers for this data structure to target molecule.
       call return_molecule_block( single_molecule_data_donor , molecule_data_diabat(i_mole_donor)%n_atom, molecule_data_diabat(i_mole_donor)%atom_index, atom_xyz=atom_data_diabat%xyz, atom_force=atom_data_diabat%force, atom_charge=atom_data_diabat%charge )
       call return_molecule_block( single_molecule_data_acceptor , molecule_data_diabat(i_mole_acceptor)%n_atom, molecule_data_diabat(i_mole_acceptor)%atom_index, atom_xyz=atom_data_diabat%xyz, atom_force=atom_data_diabat%force, atom_charge=atom_data_diabat%charge )

       allocate( force(3,molecule_data_diabat(i_mole_donor)%n_atom) , xyz_scale(3,molecule_data_diabat(i_mole_donor)%n_atom) )
       call create_scaled_direct_coordinates_molecule( xyz_scale, single_molecule_data_donor%xyz, molecule_data_diabat(i_mole_donor)%n_atom, kk, PME_data%pme_grid )
       ! call derivative_grid_Q routine with flag, this signifies don't store dQ_dr
       do i_atom=1, molecule_data_diabat(i_mole_donor)%n_atom
          call derivative_grid_Q(force,theta_conv_Q_local, single_molecule_data_donor%charge, xyz_scale,i_atom, PME_data%pme_grid, system_data_diabat%box, PME_data%spline_order, kk, flag_update)
          ! subtract this contribution, need to use molecular pointers here
          ! which point to appropriate block of pme_force_recip_diabat array
          single_molecule_data_donor%force(:,i_atom) = single_molecule_data_donor%force(:,i_atom) - force(:)
       enddo
       deallocate( force, xyz_scale )

       ! subtract the contribution of the acceptor as an h2o molecule
       allocate( force(3,molecule_data_diabat(i_mole_acceptor)%n_atom) , xyz_scale(3,molecule_data_diabat(i_mole_acceptor)%n_atom) )
       call create_scaled_direct_coordinates_molecule(xyz_scale,single_molecule_data_acceptor%xyz, molecule_data_diabat(i_mole_acceptor)%n_atom, kk, PME_data%pme_grid )
       ! call derivative_grid_Q routine with flag, this signifies don't store dQ_dr
       do i_atom=1, molecule_data_diabat(i_mole_acceptor)%n_atom
          call derivative_grid_Q(force,theta_conv_Q_local, single_molecule_data_acceptor%charge, xyz_scale,i_atom, PME_data%pme_grid, system_data_diabat%box, PME_data%spline_order, kk, flag_update)
          ! subtract this contribution, need to use molecular pointers here
          ! which point to appropriate block of pme_force_recip_diabat array
          single_molecule_data_acceptor%force(:,i_atom) = single_molecule_data_acceptor%force(:,i_atom) - force(:)
       enddo
       deallocate( force, xyz_scale )

       !******************** now change data structure topology for proton transfer
       call evb_change_data_structures_proton_transfer( i_mole_donor, i_atom_donor, i_mole_acceptor, i_atom_acceptor, i_heavy_acceptor, atom_data_diabat, molecule_data_diabat, system_data_diabat, hydronium_molecule_index_temp )

       ! now generate new pointers to molecule data structures
       call return_molecule_block( single_molecule_data_donor , molecule_data_diabat(i_mole_donor)%n_atom, molecule_data_diabat(i_mole_donor)%atom_index, atom_xyz=atom_data_diabat%xyz, atom_force=atom_data_diabat%force, atom_charge=atom_data_diabat%charge )
       call return_molecule_block( single_molecule_data_acceptor , molecule_data_diabat(i_mole_acceptor)%n_atom, molecule_data_diabat(i_mole_acceptor)%atom_index, atom_xyz=atom_data_diabat%xyz, atom_force=atom_data_diabat%force, atom_charge=atom_data_diabat%charge )


       !******************** now contribution to force from acceptor diabat data structure topology


       ! add the contribution of the donor as an h2o molecule

       allocate( force(3,molecule_data_diabat(i_mole_donor)%n_atom) , xyz_scale(3,molecule_data_diabat(i_mole_donor)%n_atom) )
       call create_scaled_direct_coordinates_molecule(xyz_scale,single_molecule_data_donor%xyz, molecule_data_diabat(i_mole_donor)%n_atom, kk, PME_data%pme_grid )
       ! call derivative_grid_Q routine with flag, this signifies don't store dQ_dr
       do i_atom=1, molecule_data_diabat(i_mole_donor)%n_atom
          call derivative_grid_Q(force,theta_conv_Q_local, single_molecule_data_donor%charge, xyz_scale,i_atom, PME_data%pme_grid, system_data_diabat%box, PME_data%spline_order, kk, flag_update)
          ! add this contribution, need to use molecular pointers here
          ! which point to appropriate block of pme_force_recip_diabat array
          single_molecule_data_donor%force(:,i_atom) = single_molecule_data_donor%force(:,i_atom) + force(:)
       enddo
       deallocate( force, xyz_scale )

       ! add the contribution of the acceptor as an h2o molecule
       allocate( force(3,molecule_data_diabat(i_mole_acceptor)%n_atom) , xyz_scale(3,molecule_data_diabat(i_mole_acceptor)%n_atom) )
       call create_scaled_direct_coordinates_molecule(xyz_scale,single_molecule_data_acceptor%xyz, molecule_data_diabat(i_mole_acceptor)%n_atom, kk, PME_data%pme_grid )
       ! call derivative_grid_Q routine with flag, this signifies don't store dQ_dr
       do i_atom=1, molecule_data_diabat(i_mole_acceptor)%n_atom
          call derivative_grid_Q(force,theta_conv_Q_local, single_molecule_data_acceptor%charge, xyz_scale,i_atom, PME_data%pme_grid, system_data_diabat%box, PME_data%spline_order, kk, flag_update)
          ! add this contribution, need to use molecular pointers here
          ! which point to appropriate block of pme_force_recip_diabat array
          single_molecule_data_acceptor%force(:,i_atom) = single_molecule_data_acceptor%force(:,i_atom) + force(:)
       enddo
       deallocate( force, xyz_scale )

    enddo loop1


    !***************************************** now map forces back to principle diabat topology using recursive subroutine
    ! output pme_force_recip_diabat, note atom_data_diabat%force is not a  pointer...
    pme_force_recip_diabat = atom_data_diabat%force

    i_hop=1  ! i_hop is updated in recursive calls, start at 1
    ! molecule_data_diabat topology will be changed, but that's ok as we're done with it here
    call map_diabat_force_to_principle_recursive( i_diabat, i_hop, i_mole_principle, molecule_data_diabat, pme_force_recip_diabat )


  end subroutine update_reciprocal_space_force_dQ_dr




  !**************************************
  !  This subroutine calculates the non-Coulombic, non-Lennard Jones,
  !  intermolecular repulsion between hydronium and water molecules
  !  see JPC B, 2008, 112, 467-482 (and errata, JPC B, 2008, 112, 7146)
  !
  ! input "i_proton" and "i_heavy" give the atom indices of
  ! the acidic proton and it's connected basic atom of the donor
  ! which are used to calculate the repulsive interactions
  !**************************************
  subroutine ms_evb_intermolecular_repulsion( force_atoms_local, E_ms_evb_repulsion, system_data_diabat, molecule_data_diabat, atom_data_diabat, hydronium_molecule_index_local )
    real*8, dimension(:,:), intent(inout) :: force_atoms_local
    real*8, intent(out)  :: E_ms_evb_repulsion
    integer, dimension(:), intent(in)     :: hydronium_molecule_index_local
    type(system_data_type), intent(in)    :: system_data_diabat
    type(molecule_data_type), dimension(:), intent(in)   :: molecule_data_diabat
    type(atom_data_diabat_type), intent(inout)   :: atom_data_diabat

    integer :: i_mole, i_mole_hydronium 

    E_ms_evb_repulsion = 0d0

    ! loop over hydronium molecules
    ! as implemented, this should only be 1 hydronium!
    do i_mole=1, n_hydronium_molecules
       i_mole_hydronium = hydronium_molecule_index_local( i_mole )

       ! compute three-atom-special evb repulsion
       call ms_evb_three_atom_repulsion( force_atoms_local, E_ms_evb_repulsion, i_mole_hydronium, system_data_diabat, molecule_data_diabat, atom_data_diabat )

       ! now compute Born-Mayer terms
       call ms_evb_born_mayer( force_atoms_local, E_ms_evb_repulsion, i_mole_hydronium, system_data_diabat, molecule_data_diabat, atom_data_diabat  )

    end do

  end subroutine ms_evb_intermolecular_repulsion




  !**************************************
  !  This subroutine calculates the non-Coulombic, non-Lennard Jones,
  !  intermolecular repulsion between hydronium and water molecules
  !  see JPC B, 2008, 112, 467-482 (and errata, JPC B, 2008, 112, 7146)
  !
  ! input "i_proton" and "i_heavy" give the atom indices of
  ! the acidic proton and it's connected basic atom of the donor
  ! which are used to calculate the repulsive interactions
  !**************************************
  subroutine ms_evb_three_atom_repulsion( force_atoms_local, E_ms_evb_repulsion, i_mole_hydronium, system_data_diabat, molecule_data_diabat, atom_data_diabat )
    real*8, dimension(:,:), intent(inout) :: force_atoms_local
    real*8, intent(inout)  :: E_ms_evb_repulsion
    integer, intent(in)  :: i_mole_hydronium
    type(system_data_type), intent(in)      :: system_data_diabat
    type(molecule_data_type), dimension(:), intent(in)    :: molecule_data_diabat
    type(atom_data_type) , intent(in)    :: atom_data_diabat

    !******** this is a local data structure with pointers that will be set
    ! to subarrays of atom_data arrays for the specific atoms in the molecule
    type(single_molecule_data_type) :: single_molecule_data_i , single_molecule_data_j

    integer :: j_mole, i_atom, j_atom, index, i_type, i_type_H, i_type_heavy, j_type, i_atom_oxygen, i_heavy, index1
    real*8, dimension(3)  :: shift, rij, rij_O, q, fij
    real*8  :: r_ij, sum, q2, r_OO, fac_OO, fac_OH, exp_q, switch_OO, dswitch_OO
    ! interaction parameters
    ! these are for the heavy-atom/heavy-atom interaction
    real*8 ::  B, bl, d0_heavy, blprime, rs_heavy, rc_heavy


    ! set pointers for this data structure to hydronium molecule
    ! set force pointer to local force_atoms_local array
    call return_molecule_block( single_molecule_data_i , molecule_data(i_mole_hydronium)%n_atom, molecule_data(i_mole_hydronium)%atom_index, atom_xyz=atom_data_diabat%xyz, atom_force=force_atoms_local, atom_type_index=atom_data_diabat%atom_type_index )

    ! index of acid hydrogen type
    i_type_H = single_molecule_data_i%atom_type_index(molecule_data(i_mole_hydronium)%n_atom)

    call get_heavy_atom_transfer_acid( i_heavy , molecule_data(i_mole_hydronium)%molecule_type_index )

    ! index of acid heavy atom type
    i_type_heavy = single_molecule_data_i%atom_type_index(i_heavy)

    ! loop over solvent molecules
    do j_mole=1, system_data%n_mole
       if ( j_mole /= i_mole_hydronium ) then

          ! set pointers for this data structure to j_mole
          call return_molecule_block( single_molecule_data_j , molecule_data(j_mole)%n_atom, molecule_data(j_mole)%atom_index, atom_xyz=atom_data_diabat%xyz, atom_force=force_atoms_local, atom_type_index=atom_data_diabat%atom_type_index )

          ! see if any of the atoms of this solvent molecule are involved in the repulsive interaction
          do j_atom = 1 , molecule_data(j_mole)%n_atom

             j_type = single_molecule_data_j%atom_type_index(j_atom)
             call get_index_atom_set( index1, evb_donor_acceptor_interaction, (/j_type, i_type_heavy, i_type_H/) )

             if ( index1 > 0 ) then
                ! get interaction parameters
                B        = evb_donor_acceptor_parameters(index1,1) ; bl       = evb_donor_acceptor_parameters(index1,2)
                d0_heavy = evb_donor_acceptor_parameters(index1,3) ; blprime  = evb_donor_acceptor_parameters(index1,4)
                rs_heavy = evb_donor_acceptor_parameters(index1,5) ; rc_heavy = evb_donor_acceptor_parameters(index1,6)

                shift = pbc_shift( single_molecule_data_i%xyz(:,i_heavy), single_molecule_data_j%xyz(:,j_atom), system_data%box, system_data%xyz_to_box_transform)
                rij_O = -pbc_dr( single_molecule_data_i%xyz(:,i_heavy), single_molecule_data_j%xyz(:,j_atom), shift )
                r_OO = dsqrt( dot_product( rij_O , rij_O ) )

                ! switching function
                call ms_evb_repulsive_switch( switch_OO, dswitch_OO, r_OO, rs_heavy , rc_heavy )
                fac_OO = B * exp(-bl *( r_OO - d0_heavy ))

                sum=0d0
                ! in this loop, calculate the 
                ! sum over q coordinates, defined by (rO + rO)/2 - rH, for the oxygen-oxygen repulsion
                do i_atom=1, molecule_data(i_mole_hydronium)%n_atom

                   i_type = single_molecule_data_i%atom_type_index(i_atom)
                   if ( i_type == i_type_H ) then

                      ! hydronium proton, the negative sign creates displacement from water to hydronium atom
                      rij =  -pbc_dr( single_molecule_data_i%xyz(:,i_atom), single_molecule_data_j%xyz(:,j_atom), shift )
                      r_ij = dsqrt( dot_product( rij, rij ) )

                      ! THIS IS CORRECT equation. See errata, JPC B, 2008, 112, 7146.  Original paper
                      ! had R_HO distance instead of q for the oxygen-oxygen repulsion
                      q = ( 2d0 * single_molecule_data_j%xyz(:,j_atom) + rij_O ) / 2d0 - single_molecule_data_j%xyz(:,j_atom) + rij )
                      q2 = dot_product(q,q)
                      exp_q = exp(-blprime * q2)
                      sum = sum + exp_q

                      ! note extra negative sign for dq/drH
                      single_molecule_data_i%force(:,i_atom) = single_molecule_data_i%force(:,i_atom) + switch_OO * fac_OO * exp_q * -blprime * 2d0 * q
                      single_molecule_data_i%force(:,i_heavy) = single_molecule_data_i%force(:,i_heavy) + switch_OO * fac_OO * exp_q * blprime * q
                      single_molecule_data_j%force(:,j_atom) = single_molecule_data_j%force(:,j_atom) + switch_OO * fac_OO * exp_q * blprime * q

                   end if
                end do
                ! now oxygen-oxygen interaction, which depends on hydronium hydrogen positions    
                E_ms_evb_repulsion = E_ms_evb_repulsion + switch_OO * fac_OO * sum
                ! now derivatives of switch_OO * fac_OO
                fij = rij_O / r_OO * fac_OO * sum * ( switch_OO * bl  - dswitch_OO )

                single_molecule_data_i%force(:,i_heavy) = single_molecule_data_i%force(:,i_heavy) + fij
                single_molecule_data_j%force(:,j_atom) = single_molecule_data_j%force(:,j_atom)   - fij          

             end if
          end do
       end if
    end do



  end subroutine ms_evb_three_atom_repulsion



  !********************************
  ! these are general Born-Mayer interactions between
  ! atoms on the proton-donor and acceptor
  !********************************
  subroutine ms_evb_born_mayer( force_atoms_local, E_ms_evb_repulsion, i_mole_hydronium, system_data_diabat, molecule_data_diabat, atom_data_diabat )
    real*8, dimension(:,:), intent(inout) :: force_atoms_local
    real*8, intent(inout)  :: E_ms_evb_repulsion
    integer, intent(in)  :: i_mole_hydronium
    type(system_data_type), intent(in)      :: system_data_diabat
    type(molecule_data_type), dimension(:), intent(in)    :: molecule_data_diabat
    type(atom_data_diabat_type) , intent(in)    :: atom_data_diabat

    !******** this is a local data structure with pointers that will be set
    ! to subarrays of atom_data arrays for the specific atoms in the molecule
    type(single_molecule_data_type) :: single_molecule_data_i , single_molecule_data_j

    integer ::  j_mole, i_atom, j_atom, i_type, j_type,  index1
    real*8, dimension(3)  :: shift, rij, fij
    real*8  :: r_ij, fac_OH, switch_HO,dswitch_HO
    ! interaction parameters
    ! these are for the hydrogen/heavy-atom interaction
    real*8 ::  C, cl, d0_hyd , rs_hyd , rc_hyd


    ! set pointers for this data structure to hydronium molecule
    ! set force pointers to local force_atoms_local data structure
    call return_molecule_block( single_molecule_data_i , molecule_data(i_mole_hydronium)%n_atom, molecule_data(i_mole_hydronium)%atom_index, atom_xyz=atom_data_diabat%xyz, atom_force=force_atoms_local, atom_type_index=atom_data_diabat%atom_type_index )

    ! loop over atoms in proton donor
    do i_atom = 1 , molecule_data(i_mole_hydronium)%n_atom
       i_type = single_molecule_data_i%atom_type_index(i_atom)

       ! loop over solvent molecules
       do j_mole=1, system_data%n_mole
          if ( j_mole /= i_mole_hydronium ) then

             ! set pointers for this data structure to target molecule
             call return_molecule_block( single_molecule_data_j , atom_data , molecule_data(j_mole)%n_atom, molecule_data(j_mole)%atom_index )

             ! see if any of the atoms of this solvent molecule are involved in the repulsive interaction
             do j_atom = 1 , molecule_data(j_mole)%n_atom

                j_type = single_molecule_data_j%atom_type_index(j_atom)
                call get_index_atom_set( index1, evb_proton_acceptor_interaction, (/j_type, i_type/) )
                if ( index1 > 0 ) then

                   C       = evb_proton_acceptor_parameters(index1,1) ; cl       = evb_proton_acceptor_parameters(index1,2) 
                   d0_hyd  = evb_proton_acceptor_parameters(index1,3) ; rs_hyd   = evb_proton_acceptor_parameters(index1,4) 
                   rc_hyd  = evb_proton_acceptor_parameters(index1,5) 

                   shift = pbc_shift( single_molecule_data_i%xyz(:,i_atom), single_molecule_data_j%xyz(:,j_atom), system_data%box, system_data%xyz_to_box_transform)

                   ! hydronium proton, the negative sign creates displacement from water to hydronium atom
                   rij = -pbc_dr( single_molecule_data_i%xyz(:,i_atom), single_molecule_data_j%xyz(:,j_atom), shift ) 
                   r_ij = dsqrt( dot_product( rij, rij ) )

                   fac_OH = C * exp( -cl * ( r_ij - d0_hyd ) )
                   ! switching function
                   call ms_evb_repulsive_switch( switch_HO, dswitch_HO, r_ij, rs_hyd , rc_hyd )

                   fij = rij / r_ij * fac_OH * ( switch_HO * cl  - dswitch_HO )
                   E_ms_evb_repulsion = E_ms_evb_repulsion + switch_HO * fac_OH

                   single_molecule_data_i%force(:,i_atom) = single_molecule_data_i%force(:,i_atom) + fij(:)
                   single_molecule_data_j%force(:,j_atom) = single_molecule_data_j%force(:,j_atom) - fij(:)

                end if
             enddo
          endif
       enddo
    enddo


  end subroutine ms_evb_born_mayer




  !************************************
  ! this is switching function (and derivative) for repulsive hydronium-water
  ! interactions
  !************************************
  subroutine ms_evb_repulsive_switch( switch, dswitch, r, rs , rc )
    real*8, intent(in) :: r, rs, rc
    real*8, intent(out) :: switch, dswitch
    real*8 :: term1, term2

    ! switch = 0 if r > rcutoff
    switch = 0d0
    dswitch = 0d0

    if ( r < rc ) then
       if (  r < rs ) then
          switch = 1d0
       else
          term1 = (r - rs)**2 / ( rc - rs )**3
          term2 = 3d0*rc - rs - 2d0*r
          switch = 1d0 - term1 * term2
          dswitch = -2d0 * (r - rs) * term2 / ( rc - rs )**3 + 2d0 * term1
       end if
    end if

  end subroutine ms_evb_repulsive_switch



  !**************************************
  ! this subroutine returns the Q_grid for
  ! a particular diabat, for which it is
  ! assumed that the Q_grid of the diabat
  ! has already been calculated and stored
  !
  ! note that to optimize memory access, we store Q_grid_diabats as 
  ! Q_grid_diabats(i,j,k,i_diabat)
  !*************************************
  subroutine ms_evb_diabat_lookup_Q_grid( Q_grid_local, diabat_index, initialize )
    real*8,dimension(:,:,:),intent(out) :: Q_grid_local
    integer,intent(in)   ::  diabat_index
    integer, intent(in), optional :: initialize

    ! if initialize
    if ( present(initialize) ) then
       if ( .not. allocated(Q_grid_diabats) ) then
          allocate( Q_grid_diabats(pme_grid,pme_grid,pme_grid,evb_max_states) )
       endif
       ! no need to zero, this takes time...
       !Q_grid_diabats=0d0
       Q_grid_filled=0
       ! at this point, global array Q_grid should contain Q_grid of principle diabat
       Q_grid_diabats(:,:,:,1) = Q_grid(:,:,:)
       Q_grid_filled(1)=1
       return
    end if

    ! otherwise, return appropriate Q_grid
    if ( Q_grid_filled(diabat_index) == 0 ) then
       stop "trying to return a Q_grid that hasn't been computed for specified diabat"
    endif

    Q_grid_local(:,:,:) = Q_grid_diabats(:,:,:,diabat_index)

  end subroutine ms_evb_diabat_lookup_Q_grid





  !*********************************
  ! this subroutine stores the forces for the evb-hamitonian matrix elements
  ! by constructing a lookup table for the non-zero elements
  !
  ! note that each time this subroutine is called for a different diabat, the 
  ! molecule_data and atom_data structures will correspond to a different topology (i.e. the protons
  ! will correspond to different molecules)
  !
  ! we use labels molecule_data_diabat and atom_data_diabat to signify this
  !
  ! therefore, we need to map the force_atoms array back to the principle
  ! xyz topology of the principle diabat, and store all the forces
  ! in this same, consistent topology
  !
  ! we map the forces back to the principle diabat using a recursive subroutine
  ! that shifts arrays
  !********************************* 
  subroutine evb_store_forces( i_mole_principle, diabat_index1, diabat_index2, n_mole, total_atoms, molecule_data_diabat, atom_data_diabat, store_index, initialize )
    integer, intent(in) :: i_mole_principle, diabat_index1, diabat_index2, n_mole, total_atoms
    type(molecule_data_type),dimension(:), intent(in)  :: molecule_data_diabat
    type(atom_data_diabat_type) , intent(in)                  :: atom_data_diabat
    integer, intent(inout)  :: store_index
    integer, intent(in),optional :: initialize

    ! this is local temporary array , which will be modified in call to 'map_diabat_array_to_principle_recursive'
    ! we define local data structure because we don't want to modify input molecule_data_diabat structure
    type(molecule_data_type),dimension(:), allocatable :: molecule_data_local

    integer, save :: size
    integer :: diabat_index, i_mole, i_hop

    ! initialize, allocate evb_force_store array.  Here we guess as to how large an array we need to allocate, which is based on the number of diabats plus the number of non-zero diabat coupling matrix elements.  The couplings are only non-zero for diabats that have hydronium on adjacent water molecules, so this corresponds to the number of donor hydrogen bonds per water molecule.  The number of water molecules is strictly less than the number of diabats. so 2 * n_diabats is a safe guess
    if ( present(initialize) ) then
       size = 3 * evb_max_states
       ! we call this subroutine with initialize present to reset the indices, so this array may already be allocated
       if ( .not. allocated(evb_forces_store) ) then
          ! this is a big array, so allocate it for efficient looping for column major memory storage
          allocate( evb_forces_store(3,total_atoms,size) )
       end if
       ! initialize evb_forces_lookup_index to negative value, so that we know which couplings have non-zero force
       evb_forces_lookup_index=-1
    end if

    ! store mapping
    evb_forces_lookup_index( diabat_index1 , diabat_index2 ) = store_index


    ! assign forces, however this is currently in specific diabat storage, we
    ! will rearrange to primary diabat topology
    evb_forces_store(:,:, store_index) = atom_data_diabat%force(:,:)

    if ( diabat_index2 > diabat_index1 ) then
       diabat_index = diabat_index2
    else
       diabat_index = diabat_index1
    end if

    if ( diabat_index > 1 ) then
       ! set up temporary data structures needed for backmapping
       allocate( molecule_data_local(n_mole) )
       do i_mole=1,n_mole
         allocate(molecule_data_local%atom_index(i_mole)%array(size(molecule_data_diabat(i_mole)%atom_index)))
         molecule_data_local(i_mole)%atom_index = molecule_data_diabat(i_mole)%atom_index
         molecule_data_local(i_mole)%n_atom = molecule_data_diabat(i_mole)%n_atom
       enddo
       i_hop=1  ! i_hop is updated in recursive calls, start at 1
       ! now map forces back to principle diabat topology using recursive subroutine
       call map_diabat_force_to_principle_recursive( diabat_index, i_hop, i_mole_principle, molecule_data_local, evb_forces_store(:,:, store_index) )

    end if

    store_index = store_index + 1

    ! check to make sure we're not overfilling array
    if ( store_index > size ) then
       write(*,*) " "
       write(*,*) " too many non-zero ms-evb matrix elements for allocated data structures "
       write(*,*) " please increase the size of evb_max_states parameter                   "
       write(*,*) ""
       stop
    endif


  end subroutine evb_store_forces



  !************************
  ! this subroutine maps the force data structure from an arbitrary diabatic topology
  ! back to the principle diabat topology
  ! this has to be recursive, because we have to shift arrays back using
  ! reverse order in which they were changed, and the evb_diabat_proton_log that
  ! tells us how the diabat topology was generated is in order of sequential
  ! changes from the principle diabat
  !
  ! the input force_diabat array that we're mapping back should be of structure
  ! force_diabat(:,i_atom) 
  ! so that we're mapping the 2nd dimension
  !
  ! we don't want to change the molecule_data_diabat data structure
  ! as we back-map, so we input temporary molecule_data_local
  ! data structure
  !*************************
  recursive subroutine map_diabat_force_to_principle_recursive( diabat_index, i_hop, i_mole_principle, molecule_data_local, force_diabat )
    integer, intent(in) :: diabat_index, i_mole_principle
    integer, intent(in) :: i_hop
    type(molecule_data_type), dimension(:), intent(inout) :: molecule_data_local
    real*8, dimension(:,:), intent(inout) :: force_diabat

    integer :: i_hop_local, i_mole, i_atom, i_mole_donor, i_atom_donor, i_mole_acceptor, i_atom_acceptor
    integer :: i_atom_donor_global, i_atom_acceptor_global, i_index, n_elements

    ! first donor is always principle hydronium
    ! set i_mole_acceptor equal to this, because
    ! i_mole_donor will be copied from previous i_mole_acceptor
    i_mole_acceptor = i_mole_principle

    ! this recursion essentially does the following do loop
    ! loop through all proton hops.  A negative integer signals the end of the chain
    ! do i_hop=1, evb_max_chain
    !    if ( evb_diabat_proton_log( diabat_index , i_hop, 1 ) < 0 ) exit

    if ( evb_diabat_proton_log( diabat_index , i_hop, 1 ) < 0 ) then
       ! end the recursive loop, and start back-mapping data structure in reverse order
       return
       
    else
       i_mole_donor = i_mole_acceptor
       i_atom_donor = evb_diabat_proton_log( diabat_index , i_hop, 2 )
       i_mole_acceptor = evb_diabat_proton_log( diabat_index , i_hop, 4 )

       ! call recursive subroutine 
       i_hop_local = i_hop + 1
       call subroutine map_diabat_array_to_principle_recursive( diabat_index, i_hop_local, i_mole_acceptor, molecule_data_local, force_diabat )

       ! this code evaluates only after we're done with recursion and looping back
       ! these are indices of atoms in global atom_data arrays
       i_atom_donor_global = molecule_data_local(i_mole_donor)%atom_index(i_atom_donor)
       i_atom_acceptor_global = molecule_data_local(i_mole_acceptor)%atom_index(i_atom_acceptor)

       ! this subroutine shifts data structures for atom transfer.
       ! Here, we're transfering the atom back from the acceptor to the donor to
       ! recreate the original data topology of the principle diabat
       call shift_array_data_donor_acceptor_transfer( i_mole_acceptor, i_atom_acceptor_global, i_mole_donor , i_atom_donor_global , molecule_data_local , array2dr1=force_diabat )

    end if

  end subroutine map_diabat_force_to_principle_recursive




  !**************************
  ! this subroutine shifts input data structures for a single proton transfer
  ! from donor to acceptor molecule
  !
  ! this subroutine may be used to either change data structures for a proton
  ! hop ( donor ====> acceptor), or restore data structures when backmapping
  ! from a specific diabat topology ( acceptor ====> donor ).  
  !
  ! we use a local data structure 'molecule_data_local' for flexibility on 
  ! whether or not to update a temporary or global data structure
  !
  ! to account for this versatility, we use generic labels 'transfer_from' , 
  ! 'transfer_to' to explicitly specific how we're reordering data structures
  !
  ! arrays that are updated are input as optional arguments, to allow
  ! versatility in only updating some of them
  !
  ! for example, array2dr1 should be set to the first 2d real array that needs updated
  !**************************
  subroutine shift_array_data_donor_acceptor_transfer( i_mole_transfer_from, i_atom_transfer_from_global, i_mole_transfer_to, i_atom_transfer_to_global, molecule_data_local, array2dr1, array2dr2, array2dr3, array1dr1 , array1dr2 , array1di , array1dc )
    integer, intent(in)  :: i_mole_transfer_from, i_mole_transfer_to, i_atom_transfer_from_global, i_atom_transfer_to_global
    type(molecule_data_type), dimension(:), intent(inout) :: molecule_data_local
    real*8, dimension(:,:), intent(inout), optional :: array2dr1, array2dr2, array2dr3
    real*8, dimension(:)  , intent(inout), optional :: array1dr1, array1dr2
    integer, dimension(:) , intent(inout), optional :: array1di
    character(*), dimension(:), intent(inout),optional :: array1dc

    real*8 :: array_element1(3),array_element2(3),array_element3(3), elementr1, elementr2
    integer :: elementi
    character(MAX_ANAME) :: elementc
    
    integer :: i_atom, n_elements, i_index, shift

    ! i_atom_transfer_from_global and i_atom_transfer_to_global are absolute indices of atoms in global array



           !update data structures that are present
           if ( present( array2dr1 ) ) then
               array_element1(:) = array2dr1(:,i_atom_transfer_from_global)
           end if
           if ( present( array2dr2 ) ) then
               array_element2(:) = array2dr2(:,i_atom_transfer_from_global)
           end if
           if ( present( array2dr3 ) ) then
               array_element3(:) = array2dr3(:,i_atom_transfer_from_global)
           end if
           if ( present( array1dr1 ) ) then
               elementr1 = array1dr1(i_atom_transfer_from_global)
           end if
           if ( present( array1dr2 ) ) then
               elementr2 = array1dr2(i_atom_transfer_from_global)
           end if
           if ( present( array1di ) ) then
               elementi     = array1di(i_atom_transfer_from_global)
           end if
           if ( present( array1dc ) ) then
               elementc     = array1dc(i_atom_transfer_from_global)
           end if

           !********** loop through data arrays and shift elements up or down  *********

           n_elements = abs( i_atom_transfer_to_global - i_atom_transfer_from_global )
           ! if i_atom_transfer_from_global < i_atom_transfer_to_global, the
           ! following loop would look like
           !do i_atom = i_atom_transfer_from_global, i_atom_transfer_to_global-1
           do i_index=1, n_elements

               ! if i_atom_transfer_from_global is before i_atom_transfer_to_global, shift data up in array
               if ( i_atom_transfer_from_global < i_atom_transfer_to_global ) then
                  i_atom = i_atom_transfer_from_global + i_index - 1
                  shift=1
               else
                  ! else acceptor is after donor, shift data down in array
                  ! copy from last to first
                  i_atom = i_atom_transfer_from_global - i_index + 1
                  shift=-1
               endif

               if ( present( array2dr1 ) ) then
                  array2dr1(:,i_atom) = array2dr1(:,i_atom+shift)
               end if
               if ( present( array2dr2 ) ) then
                  array2dr2(:,i_atom) = array2dr2(:,i_atom+shift)
               end if
               if ( present( array2dr3 ) ) then
                  array2dr3(:,i_atom) = array2dr3(:,i_atom+shift)
               end if
               if ( present( array1dr1 ) ) then
                  array1dr1(i_atom) = array1dr1(i_atom+shift)
               end if
               if ( present( array1dr2 ) ) then
                  array1dr2(i_atom) = array1dr2(i_atom+shift)
               end if
               if ( present( array1di ) ) then
                  array1di(i_atom)  = array1di(i_atom+shift)
               end if
               if ( present( array1dc ) ) then
                  array1dc(i_atom)  = array1dc(i_atom+shift)
               end if

           enddo

           ! now fill in H atom that was transferred
           if ( present( array2dr1 ) ) then
               array2dr1(:,i_atom_transfer_to_global) = array_element1(:)
           end if
           if ( present( array2dr2 ) ) then
               array2dr2(:,i_atom_transfer_to_global) = array_element2(:)
           end if
           if ( present( array2dr3 ) ) then
               array2dr3(:,i_atom_transfer_to_global) = array_element3(:)
           end if
           if ( present( array1dr1 ) ) then
               array1dr1(i_atom_transfer_to_global)   = elementr1
           end if
           if ( present( array1dr2 ) ) then
               array1dr2(i_atom_transfer_to_global)   = elementr2
           end if
           if ( present( array1di ) ) then
               array1di(i_atom_transfer_to_global)    = elementi
           end if
           if ( present( array1dc ) ) then
               array1dc(i_atom_transfer_to_global)    = elementc
           end if



         if ( i_atom_transfer_from_global < i_atom_transfer_to_global )  then
           ! here we've shifted elements up in array
           ! shift atom_indexes down one for all molecules beyond i_mole_transfer_from
           do i_mole = i_mole_transfer_from + 1, i_mole_transfer_to
              do i_atom = 1, molecule_data_local(i_mole)%n_atom
                 molecule_data_local(i_mole)%atom_index(i_atom) = molecule_data_local(i_mole)%atom_index(i_atom) - 1
              enddo
           enddo
         else
           ! here we've shifted elements down in array
           ! shift atom_indexes up one for all molecules after i_mole_transfer_to
           do i_mole = i_mole_transfer_to + 1 , i_mole_transfer_from
             do i_atom = 1, molecule_data_local(i_mole)%n_atom
                 molecule_data_local(i_mole)%atom_index(i_atom) = molecule_data_local(i_mole)%atom_index(i_atom) + 1
             enddo
           enddo
         end if

         ! now add back index of donated proton, that the array size is n_atom+1
         ! to account for proton transfers
         molecule_data_local(i_mole_transfer_to)%atom_index(molecule_data_local(i_mole_transfer_to)%n_atom + 1) = molecule_data_local(i_mole_transfer_to)%atom_index(molecule_data_local(i_mole_transfer_to)%n_atom)+ 1


       ! now update number of atoms in donor/acceptor
       molecule_data_local(i_mole_transfer_to)%n_atom = molecule_data_local(i_mole_transfer_to)%n_atom + 1
       molecule_data_local(i_mole_transfer_from)%n_atom = molecule_data_local(i_mole_transfer_from)%n_atom - 1



  end subroutine shift_array_data_donor_acceptor_transfer





  !**************************
  ! this subroutine retrieves the reference energy (chemical energy) of
  ! each adiabatic state.  Currently, this is stored for each acid, 
  ! assuming the contribution from the basic molecule is included.  If 
  ! the acid does not uniquely determine the system state, this will have
  ! to be generalized
  !**************************
  subroutine get_adiabatic_reference_energy( E_reference, i_type )
    real*8,intent(out) :: E_reference
    integer, intent(in) :: i_type

    E_reference = evb_reference_energy( i_type )

  end subroutine get_adiabatic_reference_energy



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


  !*********************************
  ! this subroutine returns the index of the basic heavy atom involved
  ! in a proton transfer
  ! this is a bit tricky.  Consider sulfonic acid.  For the base, there
  ! are three equivalent heavy atoms, so we can't distinguish the atom that
  ! was involved in the proton transfer
  ! however, this base is formed from a proton transfer, and the heavy atom
  ! types are not rearranged from a acid-to-base transfer (only rearranged
  ! for a base-to-acid transfer).  Therefore, we can get the index of the 
  ! heavy atom involved in the transfer from the index of the corresponding
  ! atom in the acid
  !*********************************
  subroutine get_heavy_atom_transfer_base( i_atom_base, i_type_base )
    integer, intent(out) :: i_atom_base
    integer, intent(in) :: i_type_base

    integer :: i_type_base, i_type_acid, i_type_heavy, i_atom

    ! conjugate acid molecule type
    i_type_acid = evb_conjugate_pairs(i_type_base)
    ! heavy atom of acid bonded to proton
    i_type_heavy = evb_heavy_acid_index(i_type_acid)

    ! now search the molecule type array of the acid to find the index of this atom type
    i_atom_base=-1
    do i_atom=1, MAX_N_ATOM
       if ( molecule_type(i_type_acid,i_atom) == MAX_N_ATOM_TYPE + 1 ) exit
       if ( molecule_type(i_type_acid,i_atom) == i_type_heavy ) then
          i_atom_base = i_atom
          exit
       endif
    enddo

    ! make sure we've found atom type
    if ( i_atom_base == -1 ) stop "couldn't find atom type in subroutine get_heavy_atom_transfer_base"

  end subroutine get_heavy_atom_transfer_base




  !*********************************
  ! this subroutine returns the index of the acidic heavy atom involved
  ! in a proton transfer
  !************************************
  subroutine get_heavy_atom_transfer_acid( i_atom_acid, i_type_acid )
    integer, intent(out) :: i_atom_acid
    integer, intent(in) :: i_type_acid

    integer :: i_type_heavy, i_atom

    ! heavy atom of acid
    i_type_heavy = evb_heavy_acid_index(i_type_acid)

    ! now search the molecule type array of the acid to find the index of this atom type
    i_atom_acid=-1
    do i_atom=1, MAX_N_ATOM
       if ( molecule_type(i_type_acid,i_atom) == MAX_N_ATOM_TYPE + 1 ) exit
       if ( molecule_type(i_type_acid,i_atom) == i_type_heavy ) then
          i_atom_acid = i_atom
          exit
       endif
    enddo

    ! make sure we've found atom type
    if ( i_atom_acid == -1 ) stop "couldn't find atom type in subroutine get_heavy_atom_transfer_acid"

  end subroutine get_heavy_atom_transfer_acid



  !****************************
  ! this function returns the center of mass of the
  ! zundel species, as well as donor shifts (shiftd)
  ! and acceptor shifts (shifta) relative to this center
  ! of mass
  !***************************
  function zundel_r_com( i_mole_donor, i_mole_acceptor, system_data_diabat, molecule_data_diabat , mass_donor, mass_acceptor, shiftd, shifta )
    real*8,dimension(3) :: zundel_r_com
    integer, intent(in) :: i_mole_donor, i_mole_acceptor
    type(system_data_type), intent(in) :: system_data_diabat
    type(molecule_data_type), dimension(:), intent(in) :: molecule_data_diabat
    real*8, dimension(:), intent(in)   :: mass_donor, mass_acceptor

    real*8, dimension(3), intent(out) :: shiftd, shifta

    real*8  :: total_mass_donor, total_mass_acceptor
    real*8,dimension(3) :: rda, shift, r_com_a
    integer :: i_atom

    total_mass_donor=0d0
    total_mass_acceptor=0d0

    do i_atom =1, molecule_data_diabat(i_mole_donor)%n_atom
       total_mass_donor = total_mass_donor + mass_donor(i_atom)
    enddo

    do i_atom =1, molecule_data_diabat(i_mole_acceptor)%n_atom
       total_mass_acceptor = total_mass_acceptor + mass_acceptor(i_atom)
    enddo

    ! donor and acceptor may be broken up over pbc
    shift(:) = pbc_shift( molecule_data_diabat(i_mole_donor)%r_com ,  molecule_data_diabat(i_mole_acceptor)%r_com, system_data_diabat%box , system_data_diabat%xyz_to_box_transform )
    rda(:)  =   pbc_dr( molecule_data_diabat(i_mole_donor)%r_com , molecule_data_diabat(i_mole_acceptor)%r_com, shift )

    ! shifted acceptor position
    r_com_a = r_com_donor + rda
    zundel_r_com(:) = ( total_mass_donor * molecule_data_diabat(i_mole_donor)%r_com + total_mass_acceptor * r_com_a(:) ) / ( total_mass_donor + total_mass_acceptor )

    ! shift donor is zero by above definition
    shiftd=0d0
    shifta=shift  

  end function zundel_r_com




  !**************************************
  ! this subroutine changes the atom index of
  ! the transferring proton to the acceptor molecule type
  ! when this subroutine is called, the acceptor molecule is still in its
  ! basic topology, so molecule_index_local will return the index of the base
  !
  ! input atom_type_index should be pointer to data structure for acceptor molecule
  !**************************************
  subroutine  change_proton_index_proton_transfer( atom_type_index, i_atom_acceptor , i_base_type )
    integer, intent(in) :: i_atom_acceptor, i_base_type
    integer, dimension(:),intent(inout) :: atom_type_index

    integer :: i_acid_type

    ! get conjugate acid type
    i_acid_type = evb_conjugate_pairs( i_base_type )

    ! fill in proton atom index of acceptor with appropriate proton index
    atom_type_index(i_atom_acceptor) = evb_proton_index( i_acid_type )

  end subroutine change_proton_index_proton_transfer




  !*********************************
  ! here we return a screen list to screen out atoms of molecules
  ! that we don't want to consider in pairwise interactions
  ! we input an array full of indices of molecules, and all atoms
  ! of these molecules have screening elements set to zero
  !*********************************
  subroutine create_atom_screen_list( screen_list , molecule_list , molecule_data_diabat ) 
    integer, dimension(:), intent(out) :: screen_list
    integer, dimension(:), intent(in)  :: molecule_list
    type(molecule_data_type),dimension(:), intent(in) :: molecule_data_diabat

    integer :: i , i_mole, i_index, i_atom

    ! initialize to 1
    screen_list = 1

    do i=1, size(molecule_list)
       i_mole = molecule_list(i)
        
       do i_index=1, molecule_data_diabat(i_mole)%n_atom
          i_atom = molecule_data_diabat(i_mole)%atom_index(i_index)
          screen_list(i_atom) = 0
       end do
    end do 

  end subroutine create_atom_screen_list



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
       stop
    end if

  end subroutine evb_check_number_diabats





  !******************************
  ! this subroutine prints the evb trajectory information
  !******************************
  subroutine print_evb_trajectory_data( ground_state_eigenvector , i_mole_principle,  file_io_data )
    real*8,dimension(:), intent(in) :: ground_state_eigenvector
    integer, intent(in) :: i_mole_principle
    type(file_io_data_type) :: file_io_data

    integer :: i_state, i_mole_donor, i_mole_acceptor, i_atom_donor, i_hop,  shell
    real*8  :: coefficient

    write( file_io_data%ofile_log_file_h, * ) "number of diabat states : ", diabat_index
    write( file_io_data%ofile_log_file_h, * ) "diabat state    hydronium molecule   evb coefficient  solvation shell"

    do i_state = 1 , diabat_index

       coefficient = ground_state_eigenvector(i_state) ** 2

       !*************** find hydronium molecule for each diabat

       i_mole_acceptor = i_mole_principle

       shell = 0
       ! loop through all proton hops.  A negative integer signals the end of the chain
       loop1: do i_hop=1, evb_max_chain
          if ( evb_diabat_proton_log( i_state , i_hop, 1 ) < 0 ) exit loop1
          i_mole_donor = i_mole_acceptor
          i_atom_donor = evb_diabat_proton_log( i_state , i_hop, 2 )
          i_mole_acceptor = evb_diabat_proton_log( i_state , i_hop, 4 )
          shell = shell + 1
       enddo loop1

       write( file_io_data%ofile_log_file_h, '(I5,I10,F14.6,I5)' ) i_state, i_mole_acceptor, coefficient, shell

    enddo


  end subroutine print_evb_trajectory_data



  !************************************
  ! this subroutine reads evb topology information
  ! from topology file.  This specifies the conjugate acids and bases,
  ! as well as the transfering protons and acceptor atoms of each molecule
  ! and the mapping between atom types of donor and acceptor molecules
  !************************************
  subroutine read_evb_topology( file_h, ifile_top )
    integer, intent(in)      :: file_h
    character(*), intent(in) :: ifile_top

    integer :: flag_eof, flag
    integer :: flag_evb_topology, flag_evb_pairs, flag_acid_reactive_protons, flag_base_reactive_protons, flag_acid_acceptor_atoms, flag_base_acceptor_atoms, flag_conjugate_atoms
    integer :: nargs, itype1, itype2, atype1, atype2, index1, count
    integer,parameter :: max_param=20
    character(300) :: input_string
    character(20),dimension(max_param)::args
    character(MAX_ANAME) :: atomtype1, atomtype2
    character(MAX_MNAME) :: moletype1, moletype2

    open(unit=file_h,file=ifile_top,status="old")

    call read_file_find_heading( file_h, '[ evb_topology ]' , flag_evb_topology, flag_eof )
    if ( flag_eof == -1 ) goto 100  ! end of file

    ! zero arrays
    evb_acid_molecule = 0
    evb_basic_molecule = 0
    evb_conjugate_pairs = 0
    evb_conjugate_atom_index = 0 
    evb_proton_index = 0

    ! read the topology file until end, look for each [ evb_pairs ] section,
    count=0
    do 
       call read_file_find_heading( file_h, '[ evb_pairs ]' , flag_evb_pairs, flag_eof )
       if ( flag_eof == -1 ) Exit  ! end of file
       count=count+1
       call write_ms_evb_io_info( '[ evb_pairs ]' )
       call read_topology_line( file_h , input_string , flag )
       call parse(input_string," ",args,nargs)
       if ( nargs /= 4 ) stop "must have 4 arguments under [ evb_pairs ] section of topology file"
       read(args(1),*) moletype1
       read(args(2),*) moletype2
       read(args(3),*) atomtype1
       read(args(4),*) atomtype2
       call trim_end(moletype1)
       call trim_end(moletype2)
       call trim_end(atomtype1)
       call trim_end(atomtype2)
       call mtype_name_reverse_lookup( moletype1, itype1 )
       call mtype_name_reverse_lookup( moletype2, itype2 )
       call atype_name_reverse_lookup( atomtype1, atype1 )
       call atype_name_reverse_lookup( atomtype2, atype2 )
       evb_acid_molecule(itype1) = 1
       evb_basic_molecule(itype2) = 1
       evb_conjugate_pairs(itype1)=itype2
       evb_conjugate_pairs(itype2)=itype1
       evb_proton_index(itype1) = atype1
       evb_heavy_acid_index(itype1) = atype2


       ! ****************for each evb_pairs section that we find, get topology information
       ! first , reactive protons section
       call read_file_find_heading( file_h, '[ acid_reactive_protons ]' , flag_acid_reactive_protons, flag_eof )    
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
          if ( nargs /= 2 ) stop "must have 2 arguments in 'acid_reactive_protons' section of [ evb_pairs ] section of topology file"
          read(args(1),*) index1
          read(args(2),*) evb_reactive_protons( itype1, index1 )
       enddo

       call read_file_find_heading( file_h, '[ base_reactive_protons ]' , flag_base_reactive_protons, flag_eof )    
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
          if ( nargs /= 2 ) stop "must have 2 arguments in 'base_reactive_protons' section of [ evb_pairs ] section of topology file"
          read(args(1),*) index1
          read(args(2),*) evb_reactive_protons( itype2, index1 )
       enddo

       !***************** now reactive acceptor atoms section
       call read_file_find_heading( file_h, '[ acid_acceptor_atoms ]' , flag_acid_acceptor_atoms, flag_eof )    
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
          if ( nargs /= 2 ) stop "must have 2 arguments in 'acid_acceptor_atoms' section of [ evb_pairs ] section of topology file"
          read(args(1),*) index1
          read(args(2),*) evb_reactive_basic_atoms( itype1, index1 )
       enddo

       call read_file_find_heading( file_h, '[ base_acceptor_atoms ]' , flag_base_acceptor_atoms, flag_eof )    
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
          if ( nargs /= 2 ) stop "must have 2 arguments in 'base_acceptor_atoms' section of [ evb_pairs ] section of topology file"
          read(args(1),*) index1
          read(args(2),*) evb_reactive_basic_atoms( itype2, index1 )
       enddo

       !***************************** now conjugate atoms section
       call read_file_find_heading( file_h, '[ conjugate_atoms ]' , flag_conjugate_atoms, flag_eof )    
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
          if ( nargs /= 2 ) stop "must have 2 arguments in 'conjugate_atoms' section of [ evb_pairs ] section of topology file"
          read(args(1),*) atomtype1
          read(args(2),*) atomtype2
          call trim_end(atomtype1)
          call trim_end(atomtype2)
          call atype_name_reverse_lookup( atomtype1, itype1 )
          call atype_name_reverse_lookup( atomtype2, itype2 )
          evb_conjugate_atom_index( itype1 ) = itype2
          evb_conjugate_atom_index( itype2 ) = itype1
       enddo

       ! make sure we found all sections for this evb pair
       call check_evb_read( flag_acid_reactive_protons, '[ acid_reactive_protons ]' )
       call check_evb_read( flag_base_reactive_protons, '[ base_reactive_protons ]' )
       call check_evb_read( flag_acid_acceptor_atoms, '[ acid_acceptor_atoms ]' )
       call check_evb_read( flag_base_acceptor_atoms, '[ base_acceptor_atoms ]' )
       call check_evb_read( flag_conjugate_atoms, '[ conjugate_atoms ]' )

    enddo

    close(file_h)

    ! make sure we found at least 1 evb_pairs section

    if ( count == 0 ) then
       stop "need at least one evb_pairs section in topology file!"
    endif


100 continue

  end subroutine read_evb_topology





  !************************************
  ! this subroutine reads evb_parameters grouped by various
  ! heading sections from topology file, and stores them in
  ! global data arrays
  !************************************
  subroutine read_evb_parameters( file_h, ifile_top )
    character(*), intent(in) :: ifile_top
    integer, intent(in)      :: file_h

    integer :: flag_eof
    integer :: flag_evb_parameters, flag_donor_acceptor, flag_proton_acceptor, flag_anisotropic, flag_geometry_factor, flag_exchange_charge_atomic, flag_exchange_charge_proton, flag_reference_energy
    integer :: nargs, flag, count, itype1, itype2, itype3
    integer,parameter :: max_param=20
    character(300) :: input_string
    character(20),dimension(max_param)::args
    character(MAX_ANAME) :: atomtype1, atomtype2, atomtype3
    character(MAX_MNAME) :: moletype1, moletype2

    open(unit=file_h,file=ifile_top,status="old")

    ! zero arrays
    evb_reference_energy=0d0
    evb_donor_acceptor_interaction = 0
    evb_proton_acceptor_interaction = 0
    evb_diabat_coupling_interaction = 0

    !*********** evb parameters section 
    call read_file_find_heading( file_h, '[ evb_parameters ]' , flag_evb_parameters, flag_eof )
    if ( flag_eof == -1 ) goto 100  ! end of file
    call read_file_find_heading( file_h, '[ reference_energy ]' , flag_reference_energy, flag_eof )
    if ( flag_eof == -1 ) goto 100  ! end of file

    !************** reference energy
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
       if ( nargs /= 2 ) stop "must have 2 arguments in 'reference_energy' section of [ evb_parameters ] section of topology file"
       read(args(1),*) moletype1
       call trim_end(moletype1)
       call mtype_name_reverse_lookup( moletype1, itype1 )
       read(args(2),*) evb_reference_energy(itype1)
    enddo
    if ( flag_eof == -1 ) goto 100  ! end of file


    call read_file_find_heading( file_h, '[ adiabat_non_bond ]' , flag_evb_parameters, flag_eof )
    if ( flag_eof == -1 ) goto 100  ! end of file
    call read_file_find_heading( file_h, '[ donor_acceptor ]' , flag_donor_acceptor, flag_eof )
    if ( flag_eof == -1 ) goto 100  ! end of file

    ! ************************** donor-acceptor non-bond parameters
    call write_ms_evb_io_info( '[ donor_acceptor ]' )
    flag_eof=0
    count=1

    do
       call read_topology_line( file_h , input_string , flag )
       ! if end of file
       if ( flag == -1 ) then
          flag_eof=-1
          exit
       end if
       ! if end of section
       if ( flag == 1 ) exit

       if ( count > max_interaction_type ) stop "please increase value of 'max_interaction_type'"
       call parse(input_string," ",args,nargs)
       if ( nargs /= 9 ) stop "must have 9 arguments in 'donor_acceptor' section of [ evb_parameters ] section of topology file"
       read(args(1),*) atomtype1
       read(args(2),*) atomtype2
       read(args(3),*) atomtype3
       call trim_end(atomtype1)
       call trim_end(atomtype2)
       call trim_end(atomtype3)
       call atype_name_reverse_lookup( atomtype1, itype1 )
       call atype_name_reverse_lookup( atomtype2, itype2 )
       call atype_name_reverse_lookup( atomtype3, itype3 )
       ! fill in index for this interaction
       evb_donor_acceptor_interaction(count,1) = itype1
       evb_donor_acceptor_interaction(count,2) = itype2
       evb_donor_acceptor_interaction(count,3) = itype3
       read(args(4),*) evb_donor_acceptor_parameters(count,1)
       read(args(5),*) evb_donor_acceptor_parameters(count,2) 
       read(args(6),*) evb_donor_acceptor_parameters(count,3) 
       read(args(7),*) evb_donor_acceptor_parameters(count,4) 
       read(args(8),*) evb_donor_acceptor_parameters(count,5) 
       read(args(9),*) evb_donor_acceptor_parameters(count,6) 
       count = count + 1
    enddo
    if ( flag_eof == -1 ) goto 100  ! end of file


    !**************************  proton-acceptor non-bond parameters
    call read_file_find_heading( file_h, '[ proton_acceptor ]' , flag_proton_acceptor, flag_eof )
    if ( flag_eof == -1 ) goto 100  ! end of file
    call write_ms_evb_io_info( '[ proton_acceptor ]' )
    flag_eof=0
    count=1
    do
       call read_topology_line( file_h , input_string , flag )
       ! if end of file
       if ( flag == -1 ) then
          flag_eof=-1
          exit
       end if
       ! if end of section
       if ( flag == 1 ) exit

       if ( count > max_interaction_type ) stop "please increase value of 'max_interaction_type'"
       call parse(input_string," ",args,nargs)
       if ( nargs /= 7 ) stop "must have 7 arguments in 'proton_acceptor' section of [ evb_parameters ] section of topology file"
       read(args(1),*) atomtype1
       read(args(2),*) atomtype2
       call trim_end(atomtype1)
       call trim_end(atomtype2)
       call atype_name_reverse_lookup( atomtype1, itype1 )
       call atype_name_reverse_lookup( atomtype2, itype2 )
       ! fill in index for this interaction
       evb_proton_acceptor_interaction(count,1) = itype1
       evb_proton_acceptor_interaction(count,2) = itype2       
       read(args(3),*) evb_proton_acceptor_parameters(count,1)
       read(args(4),*) evb_proton_acceptor_parameters(count,2) 
       read(args(5),*) evb_proton_acceptor_parameters(count,3) 
       read(args(6),*) evb_proton_acceptor_parameters(count,4) 
       read(args(7),*) evb_proton_acceptor_parameters(count,5) 
       count = count + 1
    enddo
    if ( flag_eof == -1 ) goto 100  ! end of file


    ! ************************** diabat-coupling parameters
    call read_file_find_heading( file_h, '[ diabat_coupling ]' , flag_evb_parameters, flag_eof )
    if ( flag_eof == -1 ) goto 100  ! end of file
    call read_file_find_heading( file_h, '[ geometry_factor ]' , flag_geometry_factor, flag_eof )
    if ( flag_eof == -1 ) goto 100  ! end of file
    call write_ms_evb_io_info( '[ geometry_factor ]' )
    flag_eof=0
    count=1
    do
       call read_topology_line( file_h , input_string , flag )
       ! if end of file
       if ( flag == -1 ) then
          flag_eof=-1
          exit
       end if
       ! if end of section
       if ( flag == 1 ) exit

       if ( count > max_interaction_type ) stop "please increase value of 'max_interaction_type'"
       call parse(input_string," ",args,nargs)
       ! should have 4 arguments here, 3 atomtypes and functiontype
       if ( nargs /= 4 ) stop "should have 4 arguments present in 'diabat_coupling' section of [ evb_parameters ] section of topology file, 3 atomtypes and functiontype"
       read(args(1),*) atomtype1
       read(args(2),*) atomtype2
       read(args(3),*) atomtype3
       read(args(4),*) evb_diabat_coupling_type(count) 

       call trim_end(atomtype1)
       call trim_end(atomtype2)
       call trim_end(atomtype3)
       call atype_name_reverse_lookup( atomtype1, itype1 )
       call atype_name_reverse_lookup( atomtype2, itype2 )
       call atype_name_reverse_lookup( atomtype3, itype3 )
       ! fill in index for this interaction
       evb_diabat_coupling_interaction(count,1)=itype1
       evb_diabat_coupling_interaction(count,2)=itype2
       evb_diabat_coupling_interaction(count,3)=itype3

       ! now read parameters, based on functiontype used for coupling
       call read_topology_line( file_h , input_string , flag )
       call parse(input_string," ",args,nargs)

       Select Case( evb_diabat_coupling_type(count) )
       Case(1)
          ! MS-EVB3 function type:  10 parameters
          if ( nargs /= 10 ) stop "diabat_coupling function type of '1' corresponds to MS-EVB3 form, need 10 parameters"  
          read(args(1),*) evb_diabat_coupling_parameters(count,1)
          read(args(2),*) evb_diabat_coupling_parameters(count,2) 
          read(args(3),*) evb_diabat_coupling_parameters(count,3) 
          read(args(4),*) evb_diabat_coupling_parameters(count,4) 
          read(args(5),*) evb_diabat_coupling_parameters(count,5) 
          read(args(6),*) evb_diabat_coupling_parameters(count,6) 
          read(args(7),*) evb_diabat_coupling_parameters(count,7) 
          read(args(8),*) evb_diabat_coupling_parameters(count,8) 
          read(args(9),*) evb_diabat_coupling_parameters(count,9) 
          read(args(10),*) evb_diabat_coupling_parameters(count,10) 
       Case(2)
          ! product of Gaussians, need 4 parameters
          if ( nargs /= 4 ) stop "diabat_coupling function type of '2' is product of 2 gaussians, need 4 parameters (amplitude, exponent, exponent, center) " 
          read(args(1),*) evb_diabat_coupling_parameters(count,1)
          read(args(2),*) evb_diabat_coupling_parameters(count,2) 
          read(args(3),*) evb_diabat_coupling_parameters(count,3)
          read(args(4),*) evb_diabat_coupling_parameters(count,4)
       case default
          stop "cannot recoginize diabat_coupling function type!"
       End Select



       count = count + 1
    enddo
    if ( flag_eof == -1 ) goto 100  ! end of file


    ! ************************** exchange charges atomic
    call read_file_find_heading( file_h, '[ exchange_charge_atomic ]' , flag_exchange_charge_atomic, flag_eof )
    if ( flag_eof == -1 ) goto 100  ! end of file

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
       if ( nargs /= 2 ) stop "must have 2 arguments in 'exchange_charge_atomic' section of [ evb_parameters ] section of topology file"
       read(args(1),*) atomtype1
       call trim_end(atomtype1)
       call atype_name_reverse_lookup( atomtype1, itype1 )
       read(args(2),*) evb_exchange_charge_atomic( itype1 )
    enddo

    ! ************************** exchange charges donating proton
    call read_file_find_heading( file_h, '[ exchange_charge_proton ]' , flag_exchange_charge_proton, flag_eof )
    if ( flag_eof == -1 ) goto 100  ! end of file

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
       if ( nargs /= 3 ) stop "must have 3 arguments in 'exchange_charge_proton' section of [ evb_parameters ] section of topology file"
       read(args(1),*) moletype1
       read(args(2),*) moletype2
       call trim_end(moletype1)
       call trim_end(moletype2)
       call mtype_name_reverse_lookup( moletype1, itype1 )
       call mtype_name_reverse_lookup( moletype2, itype2 )
       read(args(3),*) evb_exchange_charge_proton( itype1, itype2 )
       evb_exchange_charge_proton( itype2, itype1 ) = evb_exchange_charge_proton( itype1, itype2 )
    enddo


    close(file_h)


100 continue

    ! make sure we found all parameters
    call check_evb_read( flag_evb_parameters, '[ evb_parameters ]' )
    call check_evb_read( flag_reference_energy, '[ reference_energy ]' )
    call check_evb_read( flag_donor_acceptor, '[ donor_acceptor ]' )
    call check_evb_read( flag_proton_acceptor, '[ proton_acceptor ]' )
!    call check_evb_read( flag_anisotropic, '[ anisotropic ]' )
    call check_evb_read( flag_geometry_factor, '[ geometry_factor ]' )
    call check_evb_read( flag_exchange_charge_atomic, '[ exchange_charge_atomic ]' )
    call check_evb_read( flag_exchange_charge_proton, '[ exchange_charge_proton ]' )


  end subroutine read_evb_parameters



  !***************************************************
  ! this subroutine writes i/o info to standard output based on
  ! the input string
  !**************************************************
  subroutine write_ms_evb_io_info( string )
    character(*),intent(in) :: string

    Select Case(string)
    Case( '[ donor_acceptor ]' )
       write(*,*) ""
       write(*,*) "...reading in '[ donor_acceptor ]' section of .top file for ms-evb simulation"
       write(*,*) " this section should specify three atom types, and the"
       write(*,*) " B   b    d0    bprime    rs     rc parameters"
       write(*,*) " atom types should be ordered as 1) heavy-atom acceptor, 2) heavy-atom donor,"
       write(*,*) " 3) proton "
       write(*,*) ""
    Case( '[ proton_acceptor ]' )
       write(*,*) ""
       write(*,*) "...reading in '[ proton_acceptor ]' section of .top file for ms-evb simulation"
       write(*,*) " this section should specify two atom types, and the"
       write(*,*) " C    c    d0    rs     rc parameters"
       write(*,*) " atom types should be ordered as 1) heavy-atom acceptor, 2) proton"
       write(*,*) "" 
    Case( '[ geometry_factor ]' )
       write(*,*) ""
       write(*,*) "...reading in '[ geometry_factor ]' section of .top file for ms-evb simulation"
       write(*,*) " this section should specify three atom types, and the"
       write(*,*) " Vconst gamma   P   k  D  beta  R0  P'  alpha   r0   parameters"
       write(*,*) " atom types should be ordered as 1) heavy-atom acceptor, 2) heavy-atom donor,"
       write(*,*) " 3) proton "
       write(*,*) ""      
    Case( '[ evb_pairs ]' )
       write(*,*) ""
       write(*,*) "...reading in '[ evb_pairs ]' section of .top file for ms-evb simulation"
       write(*,*) " three arguments should be given here, 1) name of conjugate acid,"
       write(*,*) "2) name of conjugate base, and 3) name of acidic proton"
       write(*,*) ""  
    end Select


  end subroutine write_ms_evb_io_info




  !*****************************************************
  ! this subroutine stops the code and writes an error is
  ! input flag is 0, meaning corresponding section heading wasn't found
  !*****************************************************
  subroutine check_evb_read( flag, string )
    integer, intent(in) :: flag
    character(*),intent(in) :: string

    if ( flag == 0 ) then
       write(*,*) ""
       write(*,*) "couldn't find ", string, " section in topology file!"
       write(*,*) ""
       stop
    end if

  end subroutine check_evb_read


end module ms_evb
