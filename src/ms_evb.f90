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
  real*8, dimension(:,:,:,:), allocatable :: evb_forces_store


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
       ! in particular, need to change evb_create_diabat_data_structures
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
  ! changed with the call to evb_create_diabat_data_structures subroutine
  ! local :: xyz, r_com, chg, n_atom, mass
  ! global :: hydronium_molecule_index, atom_index, molecule_index
  !
  ! additionally, if velocity array is present (vel_atom), then the velocities are updated to the correct topology given a proton transfer
  !*************************
  subroutine ms_evb_calculate_total_force_energy( adiabatic_force, adiabatic_energy, tot_n_mole, n_mole, n_atom, xyz, r_com, chg, mass, box, dfti_desc,dfti_desc_inv,log_file, vel_atom )
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
    real*8,intent(inout),dimension(:,:,:), optional :: vel_atom

    integer :: i_mole, i_mole_principle, i_mole_hydronium, new_i_mole_hydronium, principle_diabat
    ! temporary diabat data structures
    real*8, dimension(:,:,:), allocatable :: xyz_temp1, force_temp1, vel_atom_temp1
    real*8, dimension(:,:), allocatable :: r_com_temp1
    real*8, dimension(:,:), allocatable :: chg_temp1, mass_temp1
    integer, dimension(:), allocatable :: n_atom_temp1
    !****** this is junk variables to pass  *********
    integer, dimension(MAX_N_MOLE) :: n_atom_drude
    integer   :: flag_junk

    n_atom_drude = n_atom


    if ( n_hydronium_molecules > 1 ) then
       stop "can't have more than 1 hydronium"
    end if


!!$    ! loop over hydronium molecules
!!$    do i_mole=1, n_hydronium_molecules
    i_mole = 1
    i_mole_hydronium = hydronium_molecule_index( i_mole )
    call construct_evb_hamiltonian( i_mole_hydronium, tot_n_mole, n_mole, n_atom, xyz, r_com, chg, mass, box, dfti_desc,dfti_desc_inv,log_file )

    call diagonalize_evb_hamiltonian( principle_diabat, new_i_mole_hydronium,  i_mole_hydronium, adiabatic_energy, adiabatic_force, tot_n_mole, n_atom, log_file )

    ! if principle hydronium molecule has changed, then change
    ! data structures.  Note we need to pass temporary arrays here, and then copy them back
    ! this is because the arrays passed in with intent(in), or intent(out) attribute will be passed in
    ! as pointers, so we can't pass (e.g.) n_atom, n_atom, without the first intent in array changing.
    ! which will cause a bug
    !
    ! note the forces currently are in the old hydronium molecule topology, these need to be changed to
    ! the new topology
    if ( new_i_mole_hydronium /= i_mole_hydronium ) then

       allocate( xyz_temp1(MAX_N_MOLE,MAX_N_ATOM,3), force_temp1(MAX_N_MOLE,MAX_N_ATOM,3), r_com_temp1(MAX_N_MOLE,3), chg_temp1(MAX_N_MOLE,MAX_N_ATOM), mass_temp1(MAX_N_MOLE,MAX_N_ATOM),n_atom_temp1(MAX_N_MOLE) )
       ! if velocity was present, update this
       if ( present(vel_atom) ) then
          allocate( vel_atom_temp1(MAX_N_MOLE,MAX_N_ATOM,3) )
          call evb_create_diabat_data_structures( xyz, xyz_temp1, r_com, r_com_temp1, chg, chg_temp1 , n_atom, n_atom_temp1, n_atom_drude, mass, mass_temp1, box, principle_diabat, i_mole_hydronium, adiabatic_force, force_temp1, vel_atom, vel_atom_temp1 )
          vel_atom = vel_atom_temp1
          deallocate( vel_atom_temp1 )
       else
          call evb_create_diabat_data_structures( xyz, xyz_temp1, r_com, r_com_temp1, chg, chg_temp1 , n_atom, n_atom_temp1, n_atom_drude, mass, mass_temp1, box, principle_diabat, i_mole_hydronium, adiabatic_force, force_temp1 )
       endif
       xyz = xyz_temp1  
       r_com = r_com_temp1
       chg = chg_temp1
       mass = mass_temp1
       n_atom = n_atom_temp1
       adiabatic_force = force_temp1
       deallocate( xyz_temp1, force_temp1, r_com_temp1, chg_temp1, mass_temp1, n_atom_temp1 )

       !******* if we're using verlet neighbor lists, we need to update these after a proton transfer
          ! generate lennard jones verlet atom index
          call generate_verlet_atom_index( verlet_lj_atoms, verlet_molecule_map_lj, tot_n_mole, n_atom )
          ! generate electrostatic verlet atom index
          call generate_verlet_atom_index( verlet_elec_atoms, verlet_molecule_map_elec, tot_n_mole, n_atom, chg )
          call construct_verlet_list(tot_n_mole, n_mole, n_atom, xyz, box )  
          ! the "1" input to update_verlet_displacements signals to initialize the displacement array
          call update_verlet_displacements( n_mole, n_atom, xyz, box, flag_junk, 1 )

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
  subroutine diagonalize_evb_hamiltonian( principle_diabat, new_i_mole_principle, i_mole_principle, adiabatic_potential, force_atoms, tot_n_mole, n_atom, log_file )
    integer, intent(out) :: principle_diabat, new_i_mole_principle
    integer, intent(in)  :: i_mole_principle
    real*8 , intent(out) :: adiabatic_potential
    real*8, dimension(:,:,:), intent(out) :: force_atoms
    integer, intent(in) :: tot_n_mole
    integer, dimension(:) , intent(in) :: n_atom
    integer,intent(in)::log_file

    integer :: i_state, j_state, n_rot, ground_state_index, index, i_hop, i_mole, i_atom
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

             ! add forces from this matrix element
             do i_mole = 1, tot_n_mole
                do i_atom = 1 , n_atom(i_mole)
                   force_atoms(i_mole,i_atom,:) = force_atoms(i_mole,i_atom,:) + factor * evb_forces_store(:,i_atom,i_mole,index)
                enddo
             enddo

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
          write(ofile_hop_file_h,*) "step ", trajectory_step
          write(ofile_hop_file_h,*) "proton hop from ", i_mole_principle, " to ", new_i_mole_principle
       end if
       ! print
       if ( mod( i_step, n_output ) == 0 ) then
          ! if restart, don't print for initial configuration
          if ( ( restart_trajectory .eq. "no" ) .or. ( i_step > 0 ) ) then
             call print_evb_trajectory_data( ground_state_eigenvector , i_mole_principle, log_file )
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
  ! here we assume pme is being used to calculate electrostatics
  !************************************
  subroutine construct_evb_hamiltonian( i_mole_principle, tot_n_mole, n_mole, n_atom, xyz, r_com, chg, mass, box, dfti_desc,dfti_desc_inv,log_file )
    use MKL_DFTI
    use total_energy_forces
    integer, intent(in) :: i_mole_principle, tot_n_mole, n_mole
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: box,r_com, chg, mass
    real*8, intent(in), dimension(:,:,:) :: xyz
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv
    integer,intent(in)::log_file

    integer ::  diabat_index_donor,  i_mole_donor, store_index
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: force_atoms
    real*8  :: potential, E_ms_evb_repulsion, E_reference

    !****** these are junk variables to pass to calculate_total_force_energy subroutine *********
    !****** as we don't implement ms-evb for polarizable force field
    real*8 :: E_elec, E_bh, E_bond, E_angle, E_dihedral, E_elec_nopol, E_3body
    integer :: iteration
    integer, dimension(MAX_N_MOLE) :: n_atom_drude
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: xyz_drude
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: force_atoms_junk
    real*8, dimension(:,:,:), allocatable :: Q_grid_junk
    real*8     :: E_junk

    n_atom_drude = n_atom
    !********************************************************************************************

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
    call calculate_total_force_energy( force_atoms,potential, E_elec,E_elec_nopol,E_bh, E_3body, E_bond, E_angle, E_dihedral, iteration, tot_n_mole, n_mole, n_atom, n_atom_drude, r_com, xyz, chg, box, dfti_desc,dfti_desc_inv,log_file,xyz_drude)

    ! need to add contribution from special evb repulsive interactions
    call ms_evb_intermolecular_repulsion( force_atoms , E_ms_evb_repulsion, n_mole  , n_atom , xyz, hydronium_molecule_index, atom_index, molecule_index, box )

!!$    write(*,*) "msevb potential"
!!$    write(*,*) potential, E_ms_evb_repulsion

    potential = potential + E_ms_evb_repulsion

    ! now add reference chemical energy for this adiabatic state
    call get_adiabatic_reference_energy( E_reference, i_mole_donor, molecule_index )
    potential = potential + E_reference

!!$    write(*,*) "diabat potential"
!!$    write(*,*) potential


    !****************************************************************************************

    ! store the forces, the "1" is to initialize
    store_index=1
    call evb_store_forces( i_mole_principle, diabat_index, diabat_index, tot_n_mole, n_atom, force_atoms, store_index, 1 )
    ! store the energy
    evb_hamiltonian( diabat_index , diabat_index ) = potential

    Select Case(electrostatic_type)
    Case("pme")
       ! the "1" is to initialize the update Q_grid subroutine, this needs to be initialized because it stores the Q_grids for each diabat to save computation time.  Q_grid_junk is junk argument
       allocate(Q_grid_junk(pme_grid,pme_grid,pme_grid))
       call ms_evb_diabat_lookup_Q_grid( Q_grid_junk, diabat_index, 1 )
       deallocate(Q_grid_junk)
    End Select



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
    call evb_conduct_proton_transfer_recursive( i_mole_donor, diabat_index_donor,  i_mole_principle, tot_n_mole, n_mole, n_atom, xyz, r_com, chg, mass, box, log_file )


    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "ms-evb diabat calculations started at", time
    endif
    !***********************************************************************!


    ! now update energies for each diabat.  This takes care of everyting except the recirocal space component, which is done below
    call evb_hamiltonian_elements_donor_acceptor( potential, force_atoms, i_mole_principle, tot_n_mole, n_mole, n_atom, xyz, r_com, chg, mass, box, store_index, log_file )


    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "ms-evb reciprocal space started at", time
    endif
    !***********************************************************************!

    Select Case(electrostatic_type)
    Case("pme")
       ! we have stored the dQ_dr array from the principle diabat, so 
       ! we can use a fast update for the reciprocal space electrostatic energy of each diabat
       ! We parallelize this part, as it is expensive.  At this point, we should have 
       ! Q_grids stored for every diabat, and we have the dQ_dr grid
       call calculate_reciprocal_space_pme( i_mole_principle, tot_n_mole, n_atom, xyz, r_com, chg, box )
    End Select


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
  recursive subroutine evb_conduct_proton_transfer_recursive( i_mole_donor, diabat_index_donor,  i_mole_principle, tot_n_mole, n_mole, n_atom, xyz, r_com, chg, mass, box, log_file )
    integer, intent(in) :: i_mole_donor, diabat_index_donor, i_mole_principle, tot_n_mole, n_mole
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: box,r_com, chg, mass
    real*8, intent(in), dimension(:,:,:) :: xyz
    integer,intent(in)::log_file

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

       i_mole_type = molecule_index( i_mole_donor )
       ! loop over all hydrogen atoms for this donor
       do i_atom=1, n_atom(i_mole_donor)

          flag_reactive_atom = evb_reactive_protons(i_mole_type,i_atom)
          ! if reactive atom (proton), consider all possible diabats formed from the transfer of this proton
          if ( flag_reactive_atom == 1 ) then

             ! first identify molecules within first solvation shell to which we can transfer a proton
             call find_evb_reactive_neighbors( i_mole_donor, i_atom, evb_neighbor_list, tot_n_mole, n_mole, n_atom, xyz,r_com, box )


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
                call find_bonded_atom_hydrogen( i_mole_type, n_atom(i_mole_donor), i_atom, j_atom )
                index=count+1
                evb_diabat_proton_log(diabat_index_acceptor,index,1) = i_mole_donor
                evb_diabat_proton_log(diabat_index_acceptor,index,2) = i_atom
                evb_diabat_proton_log(diabat_index_acceptor,index,3) = j_atom
                evb_diabat_proton_log(diabat_index_acceptor,index,4) = i_mole_acceptor
                evb_diabat_proton_log(diabat_index_acceptor,index,5) = i_atom_acceptor                

                ! now call this recursive subroutine, with the proton acceptor molecule as the new donor, unless this is the end of a cyclic transfer
                if ( flag_cycle < 1 ) then
                   call evb_conduct_proton_transfer_recursive( i_mole_acceptor, diabat_index_acceptor,  i_mole_principle, tot_n_mole, n_mole, n_atom, xyz, r_com, chg, mass, box, log_file )
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
  !
  ! potential_ref and force_atoms_ref are the potential and force for the initial diabat, stored in that topology
  !****************************************
  subroutine evb_hamiltonian_elements_donor_acceptor( potential_ref, force_atoms_ref, i_mole_principle, tot_n_mole, n_mole, n_atom, xyz, r_com, chg, mass, box,  store_index, log_file )
    real*8, intent(in) :: potential_ref
    real*8, dimension(:,:,:), intent(in) :: force_atoms_ref
    integer, intent(in) :: i_mole_principle, tot_n_mole, n_mole
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: box,r_com, chg, mass
    real*8, intent(in), dimension(:,:,:) :: xyz
    integer, intent(inout) :: store_index
    integer,intent(in)::log_file

    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: force_atoms
    integer :: initialize, lookup_index, i_mole, i_atom, i_diabat, diabat_index_donor
    real*8 :: potential, E_elec, E_bh , E_bond, E_angle, E_dihedral 

    ! temporary diabat data structures
    real*8   :: potential_temp
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: xyz_temp, force_atoms_temp
    real*8, dimension(MAX_N_MOLE,3) :: r_com_temp
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM) :: chg_temp, mass_temp
    integer, dimension(MAX_N_MOLE) :: n_atom_temp
    integer, dimension(MAX_N_MOLE,MAX_N_ATOM)  :: atom_index_temp 
    integer, dimension(MAX_N_MOLE)    :: molecule_index_temp
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
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(n_threads,split_do,potential_ref, force_atoms_ref, xyz, r_com, chg, mass, box, tot_n_mole, n_mole, n_atom, atom_index, molecule_index,  hydronium_molecule_index, i_mole_principle, diabat_index,  evb_hamiltonian, evb_diabat_coupling_matrix, store_index ) 
    !$OMP DO SCHEDULE(dynamic, split_do)
    do i_diabat = 2, diabat_index

       ! copy data structures
       potential_temp = potential_ref
       force_atoms_temp = force_atoms_ref
       xyz_temp = xyz
       r_com_temp = r_com
       chg_temp = chg
       mass_temp = mass
       n_atom_temp = n_atom
       atom_index_temp = atom_index
       molecule_index_temp = molecule_index
       hydronium_molecule_index_temp = hydronium_molecule_index



       !************************** calculate diagonal matrix element energy and forces for this diabat
       ! important note, after call to ms_evb_diabat_force_energy, the topology of the data structures will be changed from the donor to acceptor topology
       call ms_evb_diabat_force_energy( force_atoms_temp, potential_temp, i_diabat, i_mole_principle, tot_n_mole, n_mole, n_atom_temp, mass_temp, r_com_temp, xyz_temp, chg_temp, atom_index_temp, molecule_index_temp, hydronium_molecule_index_temp, box, log_file)


       ! store the forces, this needs to be in critical section, because forces will be stored using store_index, and then store_index will be incremented
       !$OMP CRITICAL
       call evb_store_forces( i_mole_principle, i_diabat, i_diabat, tot_n_mole, n_atom_temp, force_atoms_temp, store_index )
       !$OMP END CRITICAL
       ! store the energy
       evb_hamiltonian( i_diabat , i_diabat ) = potential_temp

       !************************************* calculate off-diagonal diabatic coupling, energy and force
       ! here, the data structures should be in acceptor topology, after the call to ms_evb_diabat_force_energy
       call evb_diabatic_coupling( force_atoms_temp, potential_temp, i_diabat, i_mole_principle, tot_n_mole, n_atom_temp, r_com_temp, xyz_temp, chg_temp, molecule_index_temp, atom_index_temp, box)

       ! get donor diabat for which this coupling was calculated
       diabat_index_donor = evb_diabat_coupling_matrix(i_diabat)
       ! store the forces, note this coupling is between the current diabat, and the diabat_index_1neighbor (principle) diabat, again needs to be in critical section
       !$OMP CRITICAL
       call evb_store_forces( i_mole_principle, diabat_index_donor, i_diabat, tot_n_mole, n_atom_temp, force_atoms_temp, store_index )
       !$OMP END CRITICAL
       ! store the energy
       evb_hamiltonian( diabat_index_donor , i_diabat ) = potential_temp

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
  subroutine find_evb_reactive_neighbors( i_mole, i_atom, neighbor_list, tot_n_mole, n_mole, n_atom, xyz, r_com, box ) 
    integer, intent(in) :: i_mole, i_atom,tot_n_mole, n_mole
    integer, dimension(:,:), intent(out) :: neighbor_list
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: box,r_com
    real*8, intent(in), dimension(:,:,:) :: xyz

    integer :: j_mole, j_atom, index, atom_id1, atom_id2, j_mole_type
    real*8, dimension(3) :: r_com_i, r_com_j, dr_com, rij, shift

    index=1
    ! initialize list to negative value to signal end of neighbors
    neighbor_list=-1

    do j_mole = 1, tot_n_mole
       j_mole_type = molecule_index(j_mole)
       if ( i_mole /= j_mole ) then
          ! center of mass distance
          r_com_i(:) = r_com(i_mole,:)
          r_com_j(:) = r_com(j_mole,:)

          shift = pbc_shift( r_com_i, r_com_j, box, xyz_to_box_transform)
          dr_com = pbc_dr( r_com_i, r_com_j, shift)

          if ( dot_product( dr_com, dr_com ) < evb_first_solvation_cutoff **2 ) then
             ! loop over all atoms of  target molecule to find acceptor atoms
             ! in general, we could have multiple acceptor atoms for the same acceptor molecule
             ! and if this is the case these are stored as different neighbors
             loop2:          do j_atom = 1 , n_atom(j_mole)

                if ( evb_reactive_basic_atoms(j_mole_type, j_atom) == 1 ) then
                   ! see if this reactivity pair is within the cutoff distance
                   rij = pbc_dr( xyz(i_mole,i_atom,:), xyz(j_mole,j_atom,:), shift ) ! Note for COM cutoff, shift values unchanged
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




  !***********************************************
  ! this subroutine creates both local data structures ( *_temp )
  ! and modifies global data structures to be consistent with
  ! the molecular topology of a particular diabat (diabat_index)
  !
  ! WARNING:  THIS SUBROUTINE MODIFIES GLOBAL VARIABLES! 
  ! Global arrays atom_index, molecule_index, and hydronium_molecule_index are modified
  !
  !
  ! the new data structures are created based on the proton hopping information contained in evb_diabat_proton_log
  !
  ! if optional velocity data structure is present, then a new velocity data structure will be formed, analogously to the new xyz data structure
  !***********************************************
  subroutine evb_create_diabat_data_structures( xyz, xyz_temp, r_com, r_com_temp, chg, chg_temp , n_atom, n_atom_temp, n_atom_drude, mass, mass_temp, box, diabat_index, i_mole_principle , force , force_temp, velocity, velocity_temp )
    real*8, intent(in), dimension(:,:,:) :: xyz
    real*8, intent(in), dimension(:,:) :: r_com, chg, mass, box
    integer, intent(in), dimension(:) :: n_atom
    integer, intent(in) :: i_mole_principle, diabat_index
    real*8, intent(out), dimension(:,:,:) :: xyz_temp
    real*8, intent(out), dimension(:,:) :: r_com_temp, chg_temp, mass_temp
    integer, intent(out), dimension(:) :: n_atom_temp, n_atom_drude
    real*8, intent(in), dimension(:,:,:) :: force
    real*8, intent(out), dimension(:,:,:) :: force_temp
    real*8, intent(in), dimension(:,:,:), optional :: velocity
    real*8, intent(out), dimension(:,:,:), optional :: velocity_temp

    integer :: i_hop, i_mole_donor, i_atom_donor, i_mole_acceptor, i_atom_acceptor, i_heavy_acceptor, i_atom

    xyz_temp = xyz
    r_com_temp = r_com
    chg_temp = chg
    n_atom_temp = n_atom
    mass_temp = mass
    force_temp = force

    if ( present(velocity) ) then
       velocity_temp = velocity
    end if


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
       ! ******************* Important Note:  Passing global variables atom_index, molecule_index, hydronium_molecule_index to subroutine to be modified !!!! ***********************
       if ( present( velocity ) ) then
          call evb_change_data_structures_proton_transfer( i_mole_donor, i_atom_donor, i_mole_acceptor, i_atom_acceptor, i_heavy_acceptor, xyz_temp , r_com_temp, chg_temp, n_atom_temp, mass_temp, box, force_temp, atom_index, molecule_index, hydronium_molecule_index, velocity_temp )
       else
          call evb_change_data_structures_proton_transfer( i_mole_donor, i_atom_donor, i_mole_acceptor, i_atom_acceptor, i_heavy_acceptor, xyz_temp , r_com_temp, chg_temp, n_atom_temp, mass_temp, box, force_temp, atom_index, molecule_index, hydronium_molecule_index )
       end if

    enddo loop1

    ! n_atom_drude should just be a copy here (not using drude oscillators)
    n_atom_drude = n_atom_temp


  end subroutine evb_create_diabat_data_structures





  !*************************************
  ! this subroutine modifies  data structures
  ! to account for the new molecular topology from a proton hop from
  ! i_mole_donor (i_atom_donor) to i_mole_acceptor (i_atom_acceptor)
  !
  !**************************************
  subroutine evb_change_data_structures_proton_transfer( i_mole_donor, i_atom_donor, i_mole_acceptor, i_atom_acceptor, i_heavy_acceptor, xyz , r_com, chg, n_atom, mass, box, force, atom_index_local, molecule_index_local, hydronium_molecule_index_local, velocity )
    integer, intent(in) :: i_mole_donor, i_atom_donor, i_mole_acceptor, i_atom_acceptor, i_heavy_acceptor
    real*8, dimension(:,:,:),intent(inout) :: xyz, force
    real*8, dimension(:,:), intent(inout)  :: r_com, chg, mass
    integer, dimension(:),intent(inout)  :: n_atom
    real*8, dimension(:,:) , intent(in)    :: box
    integer, dimension(:,:) , intent(inout)  :: atom_index_local
    integer, dimension(:), intent(inout) :: molecule_index_local, hydronium_molecule_index_local
    real*8, dimension(:,:,:),intent(inout), optional :: velocity

    integer :: i_atom, i, count, i_index


    ! need to change this if more than one hydronium
    if ( n_hydronium_molecules > 1 ) stop "see code below"
    hydronium_molecule_index_local(1) = i_mole_acceptor

    !******************** first change topology information for acceptor
    n_atom(i_mole_acceptor) = n_atom(i_mole_acceptor) + 1
    ! if i_atom_acceptor is not the last atom, we need to shift the data structures up
    count=1
    do i = i_atom_acceptor + 1 , n_atom(i_mole_acceptor)
       ! start at end of array so we overwrite after we copy
       i_atom = n_atom(i_mole_acceptor) - count
       xyz(i_mole_acceptor,i_atom+1,:) = xyz(i_mole_acceptor,i_atom,:)
       mass(i_mole_acceptor,i_atom+1) = mass(i_mole_acceptor,i_atom)
       force(i_mole_acceptor,i_atom+1,:) = force(i_mole_acceptor,i_atom,:)
       if ( present(velocity) ) then
          velocity(i_mole_acceptor,i_atom+1,:) = velocity(i_mole_acceptor,i_atom,:)
       end if
       count=count+1
    enddo
    ! now copy donated proton information
    xyz(i_mole_acceptor,i_atom_acceptor,:) = xyz(i_mole_donor,i_atom_donor,:)
    mass(i_mole_acceptor,i_atom_acceptor) = mass(i_mole_donor,i_atom_donor)
    force(i_mole_acceptor,i_atom_acceptor,:) = force(i_mole_donor,i_atom_donor,:)
    if ( present(velocity) ) then
       velocity(i_mole_acceptor,i_atom_acceptor,:) = velocity(i_mole_donor,i_atom_donor,:)
    endif

    ! here it's possible that the transfered proton could be split from the rest of the
    ! molecule by a periodic image. We need to fix this, as energy routines(bonds, angles, dihedrals)
    ! assume the molecule to not be split over pbc conditions
    call check_intra_molecular_shifts(i_mole_acceptor, n_atom, xyz, box)
    r_com(i_mole_acceptor,:) = pos_com( xyz, i_mole_acceptor, n_atom , mass )


    !**********************  now change force field information

    ! first shift atom_index array up for acceptor molecule to make room for donated proton
    count=1
    do i = i_atom_acceptor + 1 , n_atom(i_mole_acceptor)
       ! start at end of array so we overwrite after we copy
       i_atom = n_atom(i_mole_acceptor) - count
       atom_index_local(i_mole_acceptor,i_atom+1) = atom_index_local(i_mole_acceptor,i_atom)
       count=count+1
    enddo

    ! copy atom_index of donated proton to acceptor molecule.  Note if we are transferring proton from an acid
    ! to a different base (not its conjugate base), we can't just copy the index
    call change_proton_index_proton_transfer( i_mole_acceptor, i_atom_acceptor , molecule_index_local, atom_index_local )
    ! atom_index_local(i_mole_acceptor,i_atom_acceptor) should now be filled in with correct proton atom type


    ! don't overwrite the atom index we just filled in
    do i_atom = 1 , n_atom(i_mole_acceptor)
       if ( i_atom /= i_atom_acceptor ) then
          atom_index_local(i_mole_acceptor, i_atom) = evb_conjugate_atom_index(atom_index_local(i_mole_acceptor,i_atom))
       end if
    enddo

    ! now map the heavy atom of the acceptor to it's specific atom type
    i_index = evb_conjugate_pairs( molecule_index_local(i_mole_acceptor) )
    atom_index_local(i_mole_acceptor, i_heavy_acceptor) = evb_heavy_acid_index(i_index)

    ! now map the donor atom index, note the proton doesn't have a conjugate atom, so don't map this, as this data will be removed
    ! anyway below
    do i_atom = 1 , n_atom(i_mole_donor)
       if ( i_atom /= i_atom_donor ) then
          atom_index_local(i_mole_donor, i_atom) = evb_conjugate_atom_index(atom_index_local(i_mole_donor,i_atom))
       end if
    enddo


    do i_atom = 1, n_atom(i_mole_acceptor)
       chg(i_mole_acceptor,i_atom) = atype_chg( atom_index_local(i_mole_acceptor, i_atom) )
    enddo

    do i_atom = 1, n_atom(i_mole_donor)
       chg(i_mole_donor,i_atom) = atype_chg( atom_index_local(i_mole_donor, i_atom) )
    enddo

    ! now we have to delete the transferred proton data from i_mole_donor, and shift the contents
    ! of the data structures for that molecule.  If the transferred proton was the last atom
    ! in the molecule's data structure, then we don't have to do anything because that data
    ! will not get looped over, since the n_atom will be decremented
    do i_atom= i_atom_donor + 1 , n_atom(i_mole_donor)
       xyz(i_mole_donor, i_atom - 1, :) = xyz(i_mole_donor, i_atom , :)
       mass(i_mole_donor, i_atom - 1) = mass(i_mole_donor, i_atom )
       chg(i_mole_donor, i_atom - 1) = chg(i_mole_donor, i_atom )
       atom_index_local(i_mole_donor, i_atom - 1) = atom_index_local(i_mole_donor, i_atom )
       force(i_mole_donor, i_atom - 1, :) = force(i_mole_donor, i_atom , :)          
       if ( present(velocity) ) then
          velocity(i_mole_donor, i_atom - 1, :) = velocity(i_mole_donor, i_atom , :)          
       end if
    enddo

    ! now decrement n_atom_temp array for donor
    n_atom(i_mole_donor) = n_atom(i_mole_donor) - 1

    ! center of mass of donor
    r_com(i_mole_donor,:) = pos_com( xyz, i_mole_donor, n_atom , mass )

    ! molecule index
    molecule_index_local(i_mole_acceptor) = evb_conjugate_pairs( molecule_index_local(i_mole_acceptor) )
    molecule_index_local(i_mole_donor) = evb_conjugate_pairs( molecule_index_local(i_mole_donor) )

    ! may need to reorder data structures for acceptor molecule based on index of acceptor
    ! heavy atom to be consistent with molecule_type array
    if ( present(velocity) ) then
       call reorder_molecule_data_structures( i_mole_acceptor, n_atom, xyz, chg, mass, force, atom_index_local, molecule_index_local, velocity )
    else
       call reorder_molecule_data_structures( i_mole_acceptor, n_atom, xyz, chg, mass, force, atom_index_local, molecule_index_local )
    endif


  end subroutine evb_change_data_structures_proton_transfer



  !************************************************
  ! this subroutine reorders data structures
  ! if they are inconsistent with molecule_type array
  ! for example, if we have donated a proton to an
  ! R-SO3 group, the new O-H oxygen must be correctly ordered
  ! in the data structure array
  !************************************************
  subroutine reorder_molecule_data_structures( i_mole_acceptor, n_atom, xyz, chg, mass, force, atom_index_local, molecule_index_local, velocity )
    integer, intent(in) :: i_mole_acceptor
    integer, dimension(:), intent(in) :: n_atom
    real*8,dimension(:,:,:), intent(inout) :: xyz, force
    real*8,dimension(:,:), intent(inout) :: chg, mass
    integer,dimension(:,:),intent(inout) :: atom_index_local
    integer, dimension(:),intent(in) :: molecule_index_local
    real*8,dimension(:,:,:),intent(inout),optional :: velocity

    integer :: i_mole_type,i_atom,j_atom, index, count, i
    real*8 :: xyz_atom(3) , force_atom(3) , chg_atom, mass_atom, atom_index_atom, velocity_atom(3)

    i_mole_type = molecule_index_local( i_mole_acceptor )

    ! loop over all atoms in molecule_index_local array
    loop1 : do i_atom=1, MAX_N_ATOM_TYPE
       if ( molecule_type(i_mole_type,i_atom) == MAX_N_ATOM_TYPE + 1 ) exit loop1

       ! if atomtypes don't correspond, then we need to shift data
       if ( molecule_type(i_mole_type,i_atom) /= atom_index_local(i_mole_acceptor,i_atom) ) then

          ! loop through remaining atoms to find first occurance of this atomtype
          index=-1
          loop2: do j_atom = i_atom+1, n_atom(i_mole_acceptor)
             if ( molecule_type(i_mole_type,i_atom) == atom_index_local(i_mole_acceptor,j_atom) ) then
                index = j_atom
                exit loop2
             endif
          enddo loop2

          if ( index == -1 ) stop "error in subroutine reorder_molecule_data_structures"

          ! now store data, then shift all data structures up
          xyz_atom(:) = xyz(i_mole_acceptor,index,:)
          force_atom(:) = force(i_mole_acceptor,index,:)
          chg_atom = chg(i_mole_acceptor,index)
          mass_atom = mass(i_mole_acceptor,index)
          atom_index_atom = atom_index_local(i_mole_acceptor,index)
          if (present(velocity)) then
             velocity_atom(:) = velocity(i_mole_acceptor,index,:)
          end if

          count=1
          do i = i_atom+1 , index
             ! start at end of array so we overwrite after we copy
             j_atom = index - count
             xyz(i_mole_acceptor,j_atom+1,:) = xyz(i_mole_acceptor,j_atom,:)
             force(i_mole_acceptor,j_atom+1,:) = force(i_mole_acceptor,j_atom,:)
             chg(i_mole_acceptor,j_atom+1) = chg(i_mole_acceptor,j_atom)
             mass(i_mole_acceptor,j_atom+1) = mass(i_mole_acceptor,j_atom)
             atom_index_local(i_mole_acceptor,j_atom+1) = atom_index_local(i_mole_acceptor,j_atom)
             if ( present(velocity) ) then
                velocity(i_mole_acceptor,j_atom+1,:) = velocity(i_mole_acceptor,j_atom,:)
             end if
             count=count+1
          enddo

          ! now fill in i_atom data locations
          xyz(i_mole_acceptor,i_atom,:) = xyz_atom(:)
          force(i_mole_acceptor,i_atom,:) = force_atom(:)
          chg(i_mole_acceptor,i_atom) = chg_atom
          mass(i_mole_acceptor,i_atom) = mass_atom
          atom_index_local(i_mole_acceptor,i_atom) = atom_index_atom
          if ( present(velocity) ) then
             velocity(i_mole_acceptor,i_atom,:) = velocity_atom(:)
          end if


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
  subroutine evb_diabatic_coupling( force_atoms, potential, diabat_index, i_mole_principle, tot_n_mole, n_atom, r_com, xyz, chg, molecule_index_local, atom_index_local, box )
    real*8, intent(out), dimension(:,:,:) :: force_atoms
    real*8, intent(in), dimension(:,:,:) :: xyz
    real*8, intent(out)                 :: potential
    real*8, intent(in), dimension(:,:) :: r_com, chg, box
    integer, intent(in), dimension(:) :: n_atom
    integer, intent(in) :: i_mole_principle, diabat_index, tot_n_mole
    integer, intent(in), dimension(:)  :: molecule_index_local
    integer, intent(in), dimension(:,:) :: atom_index_local
    integer ::  i_hop, i_mole_donor, i_proton_donor, i_atom_donor, i_mole_acceptor, i_atom_acceptor, i_atom
    real*8 ::  Vex , Vconstij, A
    real*8, dimension(:,:,:), allocatable :: dVex
    real*8, dimension(3,3) :: dA


    potential=0d0
    force_atoms=0d0

    allocate( dVex( size(force_atoms(:,1,1)) , size(force_atoms(1,:,1)) , size(force_atoms(1,1,:)) ) )

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
    call evb_diabatic_coupling_electrostatics( Vex, dVex, i_mole_donor, i_mole_acceptor, tot_n_mole, n_atom, r_com, xyz , box, molecule_index_local, atom_index_local, chg )

    ! now calculate geometry dependent scale factor and its derivatives
    ! three derivatives are contained in the output dA array, in order:
    ! i_atom_donor, i_atom_acceptor, i_proton_donor
    call evb_diabatic_coupling_geometric( A , dA, Vconstij , i_mole_donor, i_atom_donor, i_mole_acceptor, i_atom_acceptor, n_atom, xyz , molecule_index_local, atom_index_local, box )

    ! now form evb diabatic coupling matrix element
    potential = ( Vconstij + Vex ) * A

    ! now forces for this matrix element, first forces from geometric derivative
    ! O donor
    force_atoms(i_mole_donor, i_atom_donor,:) = - ( Vconstij + Vex ) * dA(1,:)
    ! O acceptor
    force_atoms(i_mole_acceptor, i_atom_acceptor,:) = - ( Vconstij + Vex ) * dA(2,:)
    ! H central Zundel 
    force_atoms(i_mole_acceptor, n_atom(i_mole_acceptor),:) = - ( Vconstij + Vex ) * dA(3,:)

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
  subroutine evb_diabatic_coupling_geometric( A , dA, Vconstij, i_mole_donor, i_atom_donor, i_mole_acceptor,  i_atom_acceptor, n_atom, xyz , molecule_index_local, atom_index_local, box )
    real*8, intent(out) :: A, Vconstij
    real*8, dimension(:,:), intent(out) :: dA
    integer, intent(in)    :: i_mole_donor, i_mole_acceptor
    integer, intent(out)  :: i_atom_donor, i_atom_acceptor
    integer, dimension(:), intent(in) :: n_atom
    real*8, intent(in), dimension(:,:,:) :: xyz
    real*8, intent(in), dimension(:,:) :: box
    integer, dimension(:), intent(in) :: molecule_index_local
    integer, dimension(:,:), intent(in) :: atom_index_local

    real*8 , dimension(3) :: r_O1, r_O2, r_H, r_OO , q, shift, r_ij
    integer :: itype1, itype2, itype3, index, function_type


    ! get heavy atoms involved in proton transfer, note acceptor should currently be in
    ! acid topology
    call get_heavy_atom_transfer_base( i_atom_donor, i_mole_donor, molecule_index_local )
    call get_heavy_atom_transfer_acid( i_atom_acceptor, i_mole_acceptor, molecule_index_local )

    ! test
!!$    if ( i_mole_donor == 1 ) then
!!$      write(*,*) "geometric diabatic coupling for donor molecule, atom "
!!$      write(*,*) i_mole_donor, i_atom_donor
!!$      write(*,*) "to acceptor atom ", i_atom_acceptor
!!$    endif
    ! test


    ! first get distances, need r_OO and q, q is defined by
    ! q = ( rO1 + rO2 ) / 2 - rH
    ! to consider pbc, shift all atoms relative to the donor oxygen

    r_O1(:) = xyz(i_mole_donor, i_atom_donor, :)
    r_O2(:) = xyz(i_mole_acceptor, i_atom_acceptor, :)

    shift = pbc_shift( r_O1 , r_O2 , box , xyz_to_box_transform )
    r_ij = pbc_dr( r_O1 , r_O2 , shift )
    r_O2(:) = r_O1(:) + r_ij(:)

    ! H atom is last atom in acceptor, because we've already transferred the proton
    r_H(:) = xyz(i_mole_acceptor, n_atom(i_mole_acceptor) , : )
    r_ij = pbc_dr( r_O1 , r_H , shift )
    r_H(:) = r_O1(:) + r_ij(:)

    r_OO = r_O1 - r_O2
    q = ( r_O1 + r_O2 ) / 2d0 - r_H

    ! get parameters for this diabatic coupling term
    ! the diabatic coupling should be symmetric in atom types,
    ! so parameters shouldn't depend on whether we are considering
    ! the donor or acceptor topology
    itype1 = atom_index_local(i_mole_donor,i_atom_donor)
    itype2 = atom_index_local(i_mole_acceptor,i_atom_acceptor)
    itype3 = atom_index_local(i_mole_acceptor, n_atom(i_mole_acceptor))

    call get_index_atom_set( index, evb_diabat_coupling_interaction, (/itype1, itype2 , itype3/) )
    if ( index == -1 ) stop "couldn't find index in subroutine 'get_index_atom_set'"

    ! get coupling function type
    function_type = evb_diabat_coupling_type(index)

    call evb_diabatic_coupling_function( A , Vconstij , dA , function_type , evb_diabat_coupling_parameters(index,:) , q , r_OO )

    ! test
!!$    if ( i_mole_donor == 1 ) then
!!$    write(*,*) "q" , sqrt(dot_product(q,q))
!!$    write(*,*) "rOO", sqrt(dot_product(r_OO,r_OO))
!!$    write(*,*) "coupling", A
!!$    end if

    ! test

    ! Note dA derivative array is returned as follows:
    ! dA(1,:) is derivative w.r.t. O_donor coordinates
    ! dA(2,:) is derivative w.r.t. O_acceptor coordinates
    ! dA(3,:) is derivative w.r.t. H zundel coordinates

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
       dA(1,:) = dfac1 * fac2 * fac3 * 0.5d0 * q(:) / q_mag
       dA(1,:) = dA(1,:) + fac1 * dfac2 * fac3 * r_OO(:) / r_OO_mag
       dA(1,:) = dA(1,:) + fac1 * fac2 * dfac3 * r_OO(:) / r_OO_mag

       ! derivative w.r.t. O_acceptor coordinates
       dA(2,:) = dfac1 * fac2 * fac3 * 0.5d0 * q(:) / q_mag
       dA(2,:) = dA(2,:) + fac1 * dfac2 * fac3 * -r_OO(:) / r_OO_mag
       dA(2,:) = dA(2,:) + fac1 * fac2 * dfac3 * -r_OO(:) / r_OO_mag  

       ! derivative w.r.t. H zundel coordinates
       dA(3,:) = dfac1 * fac2 * fac3 * -q(:) / q_mag


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
       dA(1,:) = dfac1 * fac2 * 0.5d0 * q(:) / q_mag
       dA(1,:) = dA(1,:) + fac1 * dfac2 * r_OO(:) / r_OO_mag

       ! derivative w.r.t. O_acceptor coordinates
       dA(2,:) = dfac1 * fac2 * 0.5d0 * q(:) / q_mag
       dA(2,:) = dA(2,:) + fac1 * dfac2 * -r_OO(:) / r_OO_mag

       ! derivative w.r.t. H zundel coordinates
       dA(3,:) = dfac1 * fac2 * -q(:) / q_mag

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
  subroutine evb_diabatic_coupling_electrostatics( Vex, dVex, i_mole_donor, i_mole_acceptor, tot_n_mole, n_atom, r_com, xyz , box, molecule_index_local, atom_index_local, chg )
    real*8, intent(out) :: Vex
    real*8, intent(out), dimension(:,:,:) :: dVex
    integer, intent(in) :: i_mole_donor, i_mole_acceptor, tot_n_mole
    real*8, intent(in), dimension(:,:,:) :: xyz
    real*8, intent(in), dimension(:,:) :: r_com, chg, box
    integer, intent(in), dimension(:)  :: molecule_index_local
    integer, intent(in) , dimension(:,:) :: atom_index_local
    integer, intent(in), dimension(:) :: n_atom

    integer :: j_mole, i_atom, j_atom,  i_type, i_mole_type, j_mole_type
    real*8 :: q_i , q_j, r_mag , q_exchange_transfer
    real*8, dimension(3) :: shift, shiftd, shifta, r_ij, dV_ij, r_com_zundel, xyz_donor, xyz_acceptor, dr(3)

    Vex=0d0
    dVex=0d0

    ! get Zundel center of mass, note we use this later for pbc shifts.  Note that the donor and acceptor may be split over
    ! periodic boundary conditions, and so there may be a shift from the center of mass from the zundel to either the donor
    ! or acceptor.  these shifts are stored as shifta and shiftd

    ! therefore, to get the positions of the donor and acceptor at the zundel r_com location, use
    ! pbc_dr( r_com_zundel , xyz(i_mole_donor,i_atom,:), shiftd ) and 
    ! pbc_dr( r_com_zundel , xyz(i_mole_acceptor,i_atom,:), shifta )

    r_com_zundel = zundel_r_com( i_mole_donor, i_mole_acceptor, r_com(i_mole_donor,:) , r_com(i_mole_acceptor,:) , shiftd, shifta , box,  n_atom, atom_index_local )

    ! this is a sum over coulomb interactions of all water molecules with the 7 atoms of the 
    ! H5O2+ Zundel complex


    ! first calculate interaction with donor
    ! donor is currently in its basic topology, as proton has been transferred
    do i_atom=1, n_atom(i_mole_donor)
       ! figure out charge
       i_type = atom_index_local(i_mole_donor, i_atom )
       q_i = evb_exchange_charge_atomic( i_type )

       ! note in this subroutine, PBC boundary conditions are taken with respect to the whole Zundel complex, rather
       ! than the donor and acceptor water molecules individually.  This is so that this energy is invarient to switching the
       ! donor and acceptor, which may not be the case otherwise, since the centers of mass of the molecules would change, and therefore
       ! there's no guarantee that the same minimum images would be taken in both cases

       ! get donor coordinate relative to zundel
       dr(:)  =   pbc_dr( r_com_zundel(:) , xyz(i_mole_donor,i_atom,:), shiftd )
       xyz_donor = r_com_zundel + dr

       do j_mole = 1, tot_n_mole  ! loop over solvent water molecules
          if ( ( j_mole /= i_mole_donor ) .and. ( j_mole /= i_mole_acceptor ) ) then
             shift = pbc_shift( r_com_zundel(:) , r_com(j_mole,:) , box , xyz_to_box_transform )
             do j_atom =1 , n_atom(j_mole)

                q_j = chg( j_mole, j_atom )
                r_ij = -pbc_dr( xyz_donor(:) , xyz(j_mole, j_atom, :) , shift )
                r_mag = dsqrt( dot_product( r_ij , r_ij ) )

                Vex = Vex + q_i * q_j / r_mag

                dV_ij = - q_i * q_j / r_mag ** 3 * r_ij
                dVex(i_mole_donor,i_atom,:) = dVex(i_mole_donor,i_atom,:) + dV_ij
                dVex(j_mole, j_atom, :) = dVex(j_mole, j_atom, :) - dV_ij
             enddo
          endif
       enddo
    enddo

    ! the transfering proton exchange charge depends on the identity of the donor and acceptor molecules
    i_mole_type = molecule_index_local(i_mole_acceptor)
    j_mole_type = molecule_index_local(i_mole_donor)
    q_exchange_transfer = evb_exchange_charge_proton( i_mole_type, j_mole_type )


    ! now calculate interaction with acceptor
    do i_atom=1, n_atom(i_mole_acceptor)
       ! figure out charge, last atom in acceptor is transferring proton
       i_type = atom_index_local(i_mole_acceptor, i_atom )
       ! note we're using acid topology for acceptor, because we've already transferred the proton
       if ( i_atom == n_atom(i_mole_acceptor) ) then
          ! transferring proton
          q_i = q_exchange_transfer
       else
          q_i = evb_exchange_charge_atomic( i_type )
       end if

       ! see above comment about use of zundel center of mass for PBC shift

       ! get acceptor coordinate relative to zundel
       dr(:)  =   pbc_dr( r_com_zundel(:) , xyz(i_mole_acceptor,i_atom,:), shifta )
       xyz_acceptor = r_com_zundel + dr

       do j_mole = 1, tot_n_mole  ! loop over solvent water molecules
          if ( ( j_mole /= i_mole_donor ) .and. ( j_mole /= i_mole_acceptor ) ) then
             shift = pbc_shift( r_com_zundel , r_com(j_mole,:) , box , xyz_to_box_transform )
             do j_atom =1 , n_atom(j_mole)

                q_j = chg( j_mole, j_atom )
                r_ij = -pbc_dr( xyz_acceptor(:) , xyz(j_mole, j_atom, :) , shift )
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






  !*********************************************
  ! this subroutine calculates the total energy and force for a 
  ! particular diabat.  Note, a call to this subroutine essentially
  ! takes the place of a call to calculate_total_force_energy subroutinem,
  ! but this calculation is much cheaper, because it only calculates the
  ! updated force and energy terms from the last proton transfer, thus
  ! requiring order(N) real-space operations
  !
  ! we use the notation xyz_temp, n_atom_temp, etc. for the data structures to reflect that
  ! such data structures will be topologically reorganized for each particular diabat.
  !
  ! IMPORTANT NOTE:  This subroutine is designed to be called from a parallel section, and thus we do not want to rely on global variables
  ! atom_index, molecule_index, etc, as such variables are different for each diabat.  Therefore, we instead use local variable arrays
  ! atom_index_temp, molecule_index_temp, etc.
  !*********************************************
  subroutine ms_evb_diabat_force_energy( force_atoms_temp, potential_temp, i_diabat, i_mole_principle, tot_n_mole, n_mole, n_atom_temp, mass_temp, r_com_temp, xyz_temp, chg_temp, atom_index_temp, molecule_index_temp, hydronium_molecule_index_temp, box, log_file)
    use pme_routines
    real*8, dimension(:,:,:),intent(inout) :: force_atoms_temp
    real*8, intent(inout)                  :: potential_temp
    integer, intent(in) :: i_diabat, i_mole_principle, tot_n_mole, n_mole
    integer, intent(inout), dimension(:) :: n_atom_temp
    real*8, intent(inout), dimension(:,:) :: r_com_temp, chg_temp, mass_temp
    real*8, intent(in), dimension(:,:) :: box
    real*8, intent(inout), dimension(:,:,:) :: xyz_temp
    integer, intent(inout) , dimension(:,:) :: atom_index_temp
    integer, intent(inout) , dimension(:) :: molecule_index_temp, hydronium_molecule_index_temp
    integer,intent(in)::log_file

    real*8, dimension(:,:,:),allocatable :: Q_grid_local
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: xyz_scale, d_force_atoms
    real*8,dimension(:,:),allocatable :: force_donor, force_acceptor
    integer :: lookup_index, i_hop,i_mole, i_atom, i_mole_acceptor, i_atom_acceptor, i_mole_donor, i_atom_donor, i_heavy_donor, i_heavy_acceptor, i_proton_index
    real*8  :: dE_donor_diabat_intra, dE_donor_diabat_lj, dE_donor_diabat_elec, dE_acceptor_diabat_intra, dE_acceptor_diabat_lj, dE_acceptor_diabat_elec, E_acceptor_ms_evb_repulsion, E_donor_ms_evb_repulsion
    real*8  :: E_reference_donor, E_reference_acceptor
    real*8,dimension(3,3) :: kk

    Select Case(electrostatic_type)
    Case("pme")
       allocate(Q_grid_local(pme_grid,pme_grid,pme_grid))
       ! set Q_grid_local to Q_grid global variable which should have Q_grid stored for the principle diabat
       Q_grid_local = Q_grid
    End Select


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
       call ms_evb_diabat_force_energy_update_intra( d_force_atoms,dE_donor_diabat_intra, i_mole_donor, i_mole_acceptor, n_atom_temp, xyz_temp, atom_index_temp, molecule_index_temp, box)
       ! now get energy and forces from donor and acceptor lj interactions (note O(N) interactions) in donor topology
       call ms_evb_diabat_force_energy_update_lj( d_force_atoms, dE_donor_diabat_lj, i_mole_donor, i_mole_acceptor, tot_n_mole, n_atom_temp, xyz_temp, molecule_index_temp , atom_index_temp, hydronium_molecule_index_temp, box)
       ! special ms-evb repulsion terms
       call ms_evb_intermolecular_repulsion( d_force_atoms , E_donor_ms_evb_repulsion, tot_n_mole  , n_atom_temp , xyz_temp, hydronium_molecule_index_temp, atom_index_temp, molecule_index_temp, box )
       ! calculate electrostatic energy
       call ms_evb_diabat_force_energy_update_electrostatic( d_force_atoms, dE_donor_diabat_elec, i_mole_donor, i_mole_acceptor, tot_n_mole, n_mole, n_atom_temp, xyz_temp, chg_temp, box, r_com_temp, molecule_index_temp  )
       ! reference chemical energy for this adiabatic state
       call get_adiabatic_reference_energy( E_reference_donor, i_mole_donor, molecule_index_temp )


       Select Case(electrostatic_type)
       Case("pme")
          ! now, while we still have data structures in the donor diabat topology, subtract the donor and acceptor contribution
          ! the updated Q_grids will be stored for use in the separate reciprocal space part
          call construct_reciprocal_lattice_vector(kk,box)
          call create_scaled_direct_coordinates( xyz_scale, xyz_temp, tot_n_mole, n_atom_temp, kk, pme_grid )
          ! the -1 signals to subtract these contributions rather than add them
          call modify_Q_grid( Q_grid_local, i_mole_donor, chg_temp, xyz_scale, n_atom_temp, pme_grid, spline_order, -1 )  ! donor
          call modify_Q_grid( Q_grid_local, i_mole_acceptor, chg_temp, xyz_scale, n_atom_temp, pme_grid, spline_order, -1 )  ! acceptor
       End Select
       ! subtract the forces from the donor topology
       force_atoms_temp = force_atoms_temp - d_force_atoms


       ! **********************   forces and energy updates from acceptor topology
       d_force_atoms=0d0
       ! Change data structures to acceptor topology 
       call evb_change_data_structures_proton_transfer( i_mole_donor, i_atom_donor, i_mole_acceptor, i_atom_acceptor, i_heavy_acceptor, xyz_temp , r_com_temp, chg_temp, n_atom_temp, mass_temp, box, force_atoms_temp, atom_index_temp, molecule_index_temp, hydronium_molecule_index_temp )

       ! get energy and forces from donor and acceptor intra-molecular bond,angle, etc interactions in acceptor topology
       call ms_evb_diabat_force_energy_update_intra( d_force_atoms, dE_acceptor_diabat_intra, i_mole_donor, i_mole_acceptor, n_atom_temp, xyz_temp, atom_index_temp, molecule_index_temp, box)
       ! now get energy and forces from donor and acceptor lj interactions (note O(N) interactions) in acceptor topology
       call ms_evb_diabat_force_energy_update_lj( d_force_atoms, dE_acceptor_diabat_lj, i_mole_donor, i_mole_acceptor, tot_n_mole, n_atom_temp, xyz_temp, molecule_index_temp, atom_index_temp, hydronium_molecule_index_temp, box)
       ! special ms-evb repulsion terms
       call ms_evb_intermolecular_repulsion( d_force_atoms , E_acceptor_ms_evb_repulsion, tot_n_mole  , n_atom_temp , xyz_temp, hydronium_molecule_index_temp, atom_index_temp, molecule_index_temp, box )
       ! now electrostatics
       call ms_evb_diabat_force_energy_update_electrostatic( d_force_atoms, dE_acceptor_diabat_elec, i_mole_donor, i_mole_acceptor, tot_n_mole, n_mole, n_atom_temp, xyz_temp, chg_temp, box, r_com_temp, molecule_index_temp )
       ! reference chemical energy for this adiabatic state
       call get_adiabatic_reference_energy( E_reference_acceptor, i_mole_acceptor, molecule_index_temp )

       Select Case(electrostatic_type)
       Case("pme")
          ! update Q_grid with donor and acceptor molecules in acceptor topology
          call create_scaled_direct_coordinates( xyz_scale, xyz_temp, tot_n_mole, n_atom_temp, kk, pme_grid )
          ! the +1 signals to add these contributions
          call modify_Q_grid( Q_grid_local, i_mole_donor, chg_temp, xyz_scale, n_atom_temp, pme_grid, spline_order, 1 )  ! donor
          call modify_Q_grid( Q_grid_local, i_mole_acceptor, chg_temp, xyz_scale, n_atom_temp, pme_grid, spline_order, 1 )  ! acceptor
       End Select

       ! output force
       force_atoms_temp = force_atoms_temp + d_force_atoms
       ! output energy
       potential_temp = potential_temp + E_reference_acceptor + dE_acceptor_diabat_intra + dE_acceptor_diabat_lj + E_acceptor_ms_evb_repulsion + dE_acceptor_diabat_elec - E_reference_donor - dE_donor_diabat_intra - dE_donor_diabat_lj - E_donor_ms_evb_repulsion - dE_donor_diabat_elec


    end do   ! end loop over proton hops

    Select Case(electrostatic_type)
    Case("pme")
       ! store this Q_grid for future use
       Q_grid_diabats(:,:,:,i_diabat) = Q_grid_local
       Q_grid_filled(i_diabat)=1
       deallocate(Q_grid_local)
    End Select




  end subroutine ms_evb_diabat_force_energy




  !**********************************
  ! this subroutine calculates the lj force and energy
  ! for molecules i_mole_donor and i_mole_acceptor interacting
  ! with all the other molecules in the system
  !**********************************
  subroutine ms_evb_diabat_force_energy_update_lj( force_atoms, E_lj, i_mole_donor, i_mole_acceptor, n_mole, n_atom, xyz, molecule_index_temp, atom_index_temp, hydronium_molecule_index_temp, box)
    use pairwise_interaction
    real*8, dimension(:,:,:),intent(inout) :: force_atoms
    real*8, intent(out)                  :: E_lj
    integer, intent(in) :: i_mole_donor, i_mole_acceptor, n_mole
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: box
    real*8, intent(in), dimension(:,:,:) :: xyz
    integer, intent(in), dimension(:) :: molecule_index_temp
    integer, intent(in), dimension(:,:) :: atom_index_temp
    integer, intent(in), dimension(:)  :: hydronium_molecule_index_temp

    integer :: i_mole, j_mole, i_atom, j_atom, atom_id1, atom_id2, i
    real*8  :: rij(3), shift(3), dr_direct(3), shift_direct(3), norm_dr2,norm_dr6,norm_dr12, term12, term6, f_ij(3),lj_cutoff2, E_intra_lj


    E_lj=0d0
    lj_cutoff2 = lj_cutoff ** 2

    i_mole = i_mole_donor
    do j_mole=1,n_mole
       if ( j_mole /= i_mole ) then
          do i_atom=1,n_atom(i_mole)
             do j_atom=1,n_atom(j_mole)
                atom_id1 = atom_index_temp(i_mole,i_atom)
                atom_id2 = atom_index_temp(j_mole,j_atom)
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
                   E_lj = E_lj + term12 - term6
                   f_ij = rij / norm_dr2 * ( 12d0 * term12  - 6d0 * term6 )

                   force_atoms(i_mole,i_atom,:) = force_atoms(i_mole,i_atom,:) + f_ij(:)
                   force_atoms(j_mole,j_atom,:) = force_atoms(j_mole,j_atom,:) - f_ij(:)

                end if
             end do
          end do
       end if
    enddo

    ! intra-molecular lj energy for donor
    call intra_lennard_jones_energy_force( E_intra_lj, force_atoms, i_mole, n_atom, molecule_index_temp, atom_index_temp, xyz , lj_cutoff2 )
    E_lj = E_lj + E_intra_lj


    i_mole = i_mole_acceptor
    do j_mole=1,n_mole
       ! don't count acceptor donor interaction twice, it has been counted above
       if ( (j_mole /= i_mole ) .and. (j_mole /= i_mole_donor ) ) then
          do i_atom=1,n_atom(i_mole)
             do j_atom=1,n_atom(j_mole)
                atom_id1 = atom_index_temp(i_mole,i_atom)
                atom_id2 = atom_index_temp(j_mole,j_atom)
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
                   E_lj = E_lj + term12 - term6
                   f_ij = rij / norm_dr2 * ( 12d0 * term12  - 6d0 * term6 )
                   force_atoms(i_mole,i_atom,:) = force_atoms(i_mole,i_atom,:) + f_ij(:)
                   force_atoms(j_mole,j_atom,:) = force_atoms(j_mole,j_atom,:) - f_ij(:)

                end if
             end do
          end do
       end if
    enddo

    ! intra-molecular lj energy for acceptor
    call intra_lennard_jones_energy_force( E_intra_lj, force_atoms, i_mole, n_atom, molecule_index_temp, atom_index_temp, xyz , lj_cutoff2 )
    E_lj = E_lj + E_intra_lj

  end subroutine ms_evb_diabat_force_energy_update_lj





  !**********************************
  ! this subroutine calculates the intra-molecular bond and angle force and energy
  ! for molecules i_mole_donor and i_mole_acceptor
  !**********************************
  subroutine ms_evb_diabat_force_energy_update_intra( force_atoms, E_intra, i_mole_donor, i_mole_acceptor, n_atom, xyz, atom_index_temp, molecule_index_temp, box)
    use bonded_interactions
    real*8, dimension(:,:,:),intent(inout) :: force_atoms
    real*8, intent(out)                  :: E_intra
    integer, intent(in) :: i_mole_donor, i_mole_acceptor
    integer, intent(in), dimension(:) :: n_atom
    integer, intent(in), dimension(:,:) :: atom_index_temp
    integer, intent(in), dimension(:)  :: molecule_index_temp
    real*8, intent(in), dimension(:,:) :: box
    real*8, intent(in), dimension(:,:,:) :: xyz

    integer ::  i_mole
    real*8  ::  E_dihedral_mole

    E_intra=0d0


    !****************** bond terms ***********************************
    i_mole = i_mole_donor
    call intra_molecular_bond_energy_force( E_intra, i_mole, force_atoms, n_atom, molecule_index_temp, atom_index_temp, xyz )
    i_mole = i_mole_acceptor
    call intra_molecular_bond_energy_force( E_intra, i_mole, force_atoms, n_atom, molecule_index_temp, atom_index_temp, xyz )

    !******************** angle terms ****************************
    i_mole = i_mole_donor
    call intra_molecular_angle_energy_force( E_intra, i_mole, force_atoms, n_atom, molecule_index_temp, atom_index_temp, xyz )
    i_mole = i_mole_acceptor
    call intra_molecular_angle_energy_force( E_intra, i_mole, force_atoms, n_atom, molecule_index_temp, atom_index_temp, xyz )

    !******************* dihedral terms ***************************
    i_mole = i_mole_donor
    call intra_molecular_dihedral_energy_force( E_dihedral_mole, i_mole, force_atoms, n_atom, molecule_index_temp, atom_index_temp, xyz ) 
    E_intra = E_intra + E_dihedral_mole

    i_mole = i_mole_acceptor
    call intra_molecular_dihedral_energy_force( E_dihedral_mole, i_mole, force_atoms, n_atom, molecule_index_temp, atom_index_temp, xyz ) 
    E_intra = E_intra + E_dihedral_mole

  end subroutine ms_evb_diabat_force_energy_update_intra






  !**********************************
  ! this subroutine calculates the real-space electrostatic force and energy
  ! for molecules i_mole_donor and i_mole_acceptor interacting
  ! with all the other molecules in the system
  ! 
  ! the reciprocal space component is calculated elsewhere
  !**********************************
  subroutine ms_evb_diabat_force_energy_update_electrostatic( force_atoms, E_elec, i_mole_donor, i_mole_acceptor, tot_n_mole, n_mole, n_atom, xyz, chg, box, r_com, molecule_index_local )
    real*8,dimension(:,:,:),intent(inout) :: force_atoms
    real*8,intent(out)   :: E_elec
    integer, intent(in) :: i_mole_donor, i_mole_acceptor, tot_n_mole, n_mole
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: chg, box,r_com
    real*8, intent(in), dimension(:,:,:) :: xyz
    integer, intent(in), dimension(:) :: molecule_index_local
    real*8 , dimension(size(force_atoms(:,1,1)),size(force_atoms(1,:,1)),size(force_atoms(1,1,:))) :: force_elec
    real*8 :: E_elec_real

    force_elec = 0d0
    ! calculate real space contribution
    call ms_evb_diabat_force_energy_update_realspace_elec( force_elec, E_elec_real, i_mole_donor, i_mole_acceptor, tot_n_mole, n_atom, xyz, chg, box, r_com, molecule_index_local )

    ! note we don't need to add Ewald_self energy here, because the self energy is the same
    ! for all the diabats, and so it just cancels in the update calculation
    E_elec =  E_elec_real  * 0.52914D0 * 627.51D0 * 4.184D0 ! convert from e^2/A to kJ/mol
    force_elec = force_elec * 0.52914D0 * 627.51D0 * 4.184D0 ! convert from e^2/A to kJ/mol
    force_atoms = force_atoms + force_elec

  end subroutine ms_evb_diabat_force_energy_update_electrostatic




  !**********************************
  ! this subroutine calculates the real space electrostatic force and energy
  ! for molecules i_mole_donor and i_mole_acceptor interacting
  ! with all the other molecules in the system
  !
  ! note for a pme, calculation, intra-molecular reciprocal space interactions
  ! for the acceptor and donor molecules need to be subtracted here
  !*********************************
  subroutine ms_evb_diabat_force_energy_update_realspace_elec( force_atoms, E_elec_real, i_mole_donor, i_mole_acceptor, n_mole, n_atom, xyz, chg, box, r_com, molecule_index_local )
    use pme_routines
    use electrostatic
    real*8,dimension(:,:,:),intent(inout) :: force_atoms
    real*8,intent(out)   :: E_elec_real
    integer, intent(in) :: i_mole_donor, i_mole_acceptor, n_mole
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: chg, box,r_com
    real*8, intent(in), dimension(:,:,:) :: xyz
    integer, intent(in), dimension(:) :: molecule_index_local

    integer :: i,i_mole, j_mole, i_atom, j_atom
    real*8 :: cutoff2, rij(3), f_ij(3), dr_direct(3), shift_direct(3), shift(3), erfc_value, coulomb, factor, norm_dr, norm_dr2, x, E_intra

    E_elec_real=0d0

    Select Case(electrostatic_type)
    Case("pme")
       cutoff2 = ewald_cutoff ** 2
    Case("cutoff")
       cutoff2 = Electro_cutoff ** 2
    end Select

    factor = 2.D0*alpha_sqrt/pi_sqrt

    ! donor
    i_mole = i_mole_donor
    do j_mole=1,n_mole
       if ( j_mole /= i_mole ) then
          do i_atom=1, n_atom(i_mole)
             do j_atom=1,n_atom(j_mole)

                rij = xyz(i_mole,i_atom,:) - xyz(j_mole,j_atom,:)
                dr_direct(:) = matmul( xyz_to_box_transform, rij )
                do i=1,3
                   shift_direct(i) = dble(floor( dr_direct(i) + 0.5d0 ))
                enddo
                shift = matmul( shift_direct , box )
                rij = rij - shift
                norm_dr2 = dot_product( rij, rij )
                if ( norm_dr2 < cutoff2 ) then
                   norm_dr = dsqrt( norm_dr2 )
                   ! pme
                   Select Case(electrostatic_type)
                   Case("pme")
                      x = norm_dr*alpha_sqrt
                      if ( x < erfc_max ) then
                         erfc_value = erfc_table(ceiling(x/erfc_max*dble(erfc_grid)))
                         coulomb = chg(i_mole,i_atom) * chg(j_mole,j_atom)  / norm_dr
                         E_elec_real = E_elec_real + erfc_value * coulomb
                         f_ij = coulomb * rij *(  erfc_value / norm_dr2 + factor * exp(-x**2) / norm_dr)
                         force_atoms(i_mole,i_atom,:) = force_atoms(i_mole,i_atom,:) + f_ij(:)
                         force_atoms(j_mole,j_atom,:) = force_atoms(j_mole,j_atom,:) - f_ij(:)
                      endif
                   Case("cutoff")
                      coulomb = chg(i_mole,i_atom) * chg(j_mole,j_atom)  / norm_dr
                      E_elec_real = E_elec_real + coulomb
                      f_ij = coulomb * rij / norm_dr2
                      force_atoms(i_mole,i_atom,:) = force_atoms(i_mole,i_atom,:) + f_ij(:)
                      force_atoms(j_mole,j_atom,:) = force_atoms(j_mole,j_atom,:) - f_ij(:)
                   End Select
                end if
             enddo
          enddo
       end if
    enddo

    ! acceptor
    i_mole = i_mole_acceptor
    do j_mole=1,n_mole
       ! we have already calculated donor acceptor real space interaction above, so don't calculate it again here
       if ( (j_mole /= i_mole ) .and. (j_mole /= i_mole_donor ) ) then
          do i_atom=1, n_atom(i_mole)
             do j_atom=1,n_atom(j_mole)

                rij = xyz(i_mole,i_atom,:) - xyz(j_mole,j_atom,:)
                dr_direct(:) = matmul( xyz_to_box_transform, rij )
                do i=1,3
                   shift_direct(i) = dble(floor( dr_direct(i) + 0.5d0 ))
                enddo
                shift = matmul( shift_direct , box )
                rij = rij - shift
                norm_dr2 = dot_product( rij, rij )
                if ( norm_dr2 < cutoff2 ) then
                   norm_dr = dsqrt( norm_dr2 )
                   ! pme
                   Select Case(electrostatic_type)
                   Case("pme")
                      x = norm_dr*alpha_sqrt
                      if ( x < erfc_max ) then
                         erfc_value = erfc_table(ceiling(x/erfc_max*dble(erfc_grid)))
                         coulomb = chg(i_mole,i_atom) * chg(j_mole,j_atom)  / norm_dr
                         E_elec_real = E_elec_real + erfc_value * coulomb
                         f_ij = coulomb * rij *(  erfc_value / norm_dr2 + factor * exp(-x**2) / norm_dr)
                         force_atoms(i_mole,i_atom,:) = force_atoms(i_mole,i_atom,:) + f_ij(:)
                         force_atoms(j_mole,j_atom,:) = force_atoms(j_mole,j_atom,:) - f_ij(:)
                      endif
                   Case("cutoff")
                      coulomb = chg(i_mole,i_atom) * chg(j_mole,j_atom)  / norm_dr
                      E_elec_real = E_elec_real + coulomb
                      f_ij = coulomb * rij / norm_dr2
                      force_atoms(i_mole,i_atom,:) = force_atoms(i_mole,i_atom,:) + f_ij(:)
                      force_atoms(j_mole,j_atom,:) = force_atoms(j_mole,j_atom,:) - f_ij(:)
                   End Select
                end if
             enddo
          enddo
       end if
    enddo

    ! now intra-molecular component.  Here, this subtracts off the reciprocal space intra-molecular interactions
    ! note these are different for the diabats, so we need to include this term here

    ! acceptor
    i_mole = i_mole_acceptor
    do i_atom=1,n_atom(i_mole)-1
       do j_atom=i_atom+1,n_atom(i_mole)
          Select Case(electrostatic_type)
          Case("pme") 
             call intra_pme_energy(E_intra,xyz,chg,i_mole,i_atom,j_atom,n_atom, molecule_index_local)
             call intra_pme_force(f_ij,xyz,chg,i_mole,n_atom,i_atom,j_atom, molecule_index_local)
          End Select
          E_elec_real = E_elec_real + E_intra
          force_atoms(i_mole,i_atom,:) = force_atoms(i_mole,i_atom,:) + f_ij(:)
          force_atoms(i_mole,j_atom,:) = force_atoms(i_mole,j_atom,:) - f_ij(:)       
       enddo
    enddo
    ! donor
    i_mole = i_mole_donor
    do i_atom=1,n_atom(i_mole)-1
       do j_atom=i_atom+1,n_atom(i_mole)
          Select Case(electrostatic_type)
          Case("pme") 
             call intra_pme_energy(E_intra,xyz,chg,i_mole,i_atom,j_atom,n_atom, molecule_index_local)
             call intra_pme_force(f_ij,xyz,chg,i_mole,n_atom,i_atom,j_atom, molecule_index_local)
          End Select
          E_elec_real = E_elec_real + E_intra
          force_atoms(i_mole,i_atom,:) = force_atoms(i_mole,i_atom,:) + f_ij(:)
          force_atoms(i_mole,j_atom,:) = force_atoms(i_mole,j_atom,:) - f_ij(:)       
       enddo
    enddo



  end subroutine ms_evb_diabat_force_energy_update_realspace_elec



  !************************************
  ! this subroutine calculates the reciprocal_space_pme
  ! energy and force for each diabat, using the dQ_dr stored grid from
  ! the principle diabat, along with the stored Q_grids from
  ! each diabat
  !***********************************
  subroutine calculate_reciprocal_space_pme( i_mole_principle, tot_n_mole, n_atom, xyz, r_com, chg, box  )
    use MKL_DFTI
    use omp_lib
    implicit none
    integer, intent(in) :: i_mole_principle, tot_n_mole
    integer,dimension(:),intent(in) :: n_atom
    real*8,dimension(:,:,:),intent(in) :: xyz
    real*8,dimension(:,:),intent(in) :: r_com
    real*8,dimension(:,:),intent(in) :: chg
    real*8,dimension(:,:),intent(in) :: box

    TYPE(DFTI_DESCRIPTOR),pointer :: dfti_desc_local,dfti_desc_inv_local
    integer :: index_store, index, i_mole, i_atom, i, j, i_diabat, n, K, status, length(3)
    real*8 :: kk(3,3), force_temp(3), E_recip_local, dE_recip
    real*8,dimension(:,:,:), allocatable :: Q_local, theta_conv_Q_local
    complex*16,dimension(:,:,:),allocatable::FQ
    real*8,dimension(:), allocatable::q_1r
    complex*16,dimension(:), allocatable::q_1d
    real*8,dimension(:,:), allocatable :: atom_list_force
    real*8,dimension(:,:,:), allocatable :: pme_force_recip_diabat
    integer,dimension(MAX_N_MOLE,MAX_N_ATOM) :: atom_molecule_map
    integer :: total_atoms_list, split_do

    n=spline_order
    K=pme_grid
    ! note three dimensions
    length=pme_grid

    call construct_reciprocal_lattice_vector(kk, box)


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

    ! note that if we are using a verlet list, we can use the verlet list as the atom-molecule map
    ! here, we generate a new map, in case we are not using a verlet list
    call generate_verlet_atom_index( total_atoms_list, atom_molecule_map, tot_n_mole, n_atom, chg )

    !***************************************************************
    ! Note this code below may look confusing, but essentially we are just
    ! doing reciprocal space pme, using stored data structures
    ! therefore, see pme_force_recip or pme_recip for comments on this algorithm
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
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(n_threads,split_do,theta_conv_Q, CB, dQ_dr,dQ_dr_index,Q_grid_diabats,xyz,chg,K,box,n,kk,spline_order, tot_n_mole, n_atom, atom_list_force_recip, E_recip, total_atoms_list, atom_molecule_map, i_mole_principle, diabat_index,evb_forces_lookup_index,evb_hamiltonian,evb_forces_store,dfti_desc_local,dfti_desc_inv_local) 
    !$OMP CRITICAL
    allocate( FQ(K,K,K), Q_local(K,K,K), theta_conv_Q_local(K,K,K), q_1r(K**3), q_1d(K**3), atom_list_force(3,total_atoms_list),pme_force_recip_diabat( tot_n_mole, MAX_N_ATOM, 3 ) )
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
       E_recip_local = 0.5D0*sum((Q_local*theta_conv_Q_local))
       ! subtract reciprocal space energy for principle diabat, and convert energies (E_recip is a stored global variable)
       dE_recip =  E_recip_local - E_recip 

       ! now force
       atom_list_force=0d0

       do index_store=1,size(dQ_dr(1,1,:))
          do j=1, size(dQ_dr(1,:,1))    
             atom_list_force(:,index_store) = atom_list_force(:,index_store) + dQ_dr(:,j,index_store) * theta_conv_Q_local(dQ_dr_index(1,j,index_store), dQ_dr_index(2,j,index_store),dQ_dr_index(3,j,index_store))
          enddo
       enddo

       ! dQ_dr is stored in scaled coordinates, convert to general coordinates
       do index_store=1,size(dQ_dr(1,1,:)) 
          force_temp=0d0
          do i=1,3
             do j=1,3
                force_temp(i) = force_temp(i) - dble(K) * kk(j,i) * atom_list_force(j,index_store)
             enddo
          enddo
          atom_list_force(:,index_store) = force_temp(:)
       enddo

       ! now subtract off the reciprocal space electrostatic force of the principle diabat ( atom_list_force_recip is a stored global data array )

       atom_list_force=atom_list_force - atom_list_force_recip

       ! convert forces back to molecular data structures
       call map_molecule_atom_data_structures_3d_a_to_m( tot_n_mole, n_atom, atom_molecule_map , atom_list_force , pme_force_recip_diabat )


       ! now get contribution to force from changes in dQ_dr for this diabat
       ! note here force should be in principle diabat data structure here, as xyz, n_atom datastructures are in principle diabat structure
       call update_reciprocal_space_force_dQ_dr( pme_force_recip_diabat, theta_conv_Q_local, i_diabat, i_mole_principle, tot_n_mole, n_atom, xyz, chg, kk, box, K, n )

       ! now convert units for energies and force
       dE_recip = dE_recip * 0.52914D0 * 627.51D0 * 4.184D0 ! convert from e^2/A to kJ/mol
       pme_force_recip_diabat = pme_force_recip_diabat * 0.52914D0 * 627.51D0 * 4.184D0 ! convert from e^2/A to kJ/mol


       ! add the energy and force differences to the stored energy and force arrays
       index = evb_forces_lookup_index(i_diabat,i_diabat)
       evb_hamiltonian(i_diabat,i_diabat) = evb_hamiltonian(i_diabat,i_diabat) + dE_recip
       do i_mole=1,tot_n_mole
          do i_atom=1,n_atom(i_mole)
             evb_forces_store(:,i_atom,i_mole,index) = evb_forces_store(:,i_atom,i_mole,index) + pme_force_recip_diabat(i_mole, i_atom, :)
          enddo
       enddo

    enddo

    !$OMP END DO NOWAIT
    deallocate( FQ, Q_local, q_1r, q_1d,theta_conv_Q_local, atom_list_force, pme_force_recip_diabat )
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
  subroutine update_reciprocal_space_force_dQ_dr( pme_force_recip_diabat, theta_conv_Q_local, i_diabat, i_mole_principle, tot_n_mole, n_atom, xyz, chg, kk, box, K, n )
    use pme_routines
    implicit none
    real*8,dimension(:,:,:),intent(inout) :: pme_force_recip_diabat
    real*8,dimension(:,:,:),intent(in)   :: theta_conv_Q_local
    integer, intent(in) :: i_diabat, i_mole_principle, tot_n_mole
    integer,dimension(:),intent(in) :: n_atom
    real*8,dimension(:,:,:),intent(in) :: xyz
    real*8,dimension(:,:),intent(in) :: chg
    real*8,dimension(:,:),intent(in) :: kk,box   
    integer, intent(in)  :: K,n

    real*8, dimension(3)   :: force
    real*8, dimension(4,3) :: xyz_scale ! 4 dimensions store either h3o, or h2o
    real*8, dimension(3,4) :: xyz_scale_transpose  ! this is just to input in reverse storage order
    real*8,dimension(:,:,:) , allocatable :: xyz_temp
    real*8,dimension(:,:), allocatable :: chg_temp
    integer, dimension(:), allocatable :: n_atom_temp
    integer, dimension(:,:), allocatable :: atom_index_temp
    integer, dimension(:), allocatable :: molecule_index_temp
    integer, dimension(:), allocatable :: hydronium_molecule_index_temp

    integer :: i_mole_donor, i_atom_donor, i_mole_acceptor, i_atom_acceptor, i_heavy_acceptor, i_hop, i_atom, i, count, flag_update=1

    ! these are junk variables to pass to subroutine evb_change_data_structures_proton_transfer
    real*8,dimension(size(xyz(:,1,1)),size(xyz(1,:,1))) :: mass_junk
    real*8,dimension(size(xyz(:,1,1)),3) :: r_com_junk

    mass_junk=0d0
    r_com_junk=0d0

    allocate( xyz_temp(size(xyz(:,1,1)),size(xyz(1,:,1)),size(xyz(1,1,:))) , chg_temp(size(chg(:,1)),size(chg(1,:))) , n_atom_temp(size(n_atom)) , atom_index_temp(size(atom_index(:,1)),size(atom_index(1,:))) )
    allocate( molecule_index_temp(size(molecule_index)), hydronium_molecule_index_temp(size(hydronium_molecule_index)) )

    ! all of these data structures, and also pme_force_recip_diabat should initially be in principle diabat data structures
    xyz_temp = xyz
    n_atom_temp = n_atom
    chg_temp = chg
    atom_index_temp = atom_index
    molecule_index_temp = molecule_index
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
       i_atom_acceptor = n_atom_temp(i_mole_acceptor) + 1

       !*************************************** subtract the contribution of the donor as an h3o+ molecule
       call create_scaled_direct_coordinates_molecule(xyz_scale,xyz_temp(i_mole_donor,:,:), n_atom_temp(i_mole_donor), kk, K )
       do i=1,3
          xyz_scale_transpose(i,:) = xyz_scale(:,i)
       enddo
       ! call derivative_grid_Q routine with flag, this signifies don't store dQ_dr
       do i_atom=1, n_atom_temp(i_mole_donor)
          call derivative_grid_Q(force,theta_conv_Q_local,chg_temp(i_mole_donor,:),xyz_scale_transpose,i_atom, K, box, n, kk, flag_update)
          ! subtract this contribution
          pme_force_recip_diabat(i_mole_donor,i_atom,:) = pme_force_recip_diabat(i_mole_donor,i_atom,:) - force(:) 
       enddo

       ! subtract the contribution of the acceptor as an h2o molecule
       call create_scaled_direct_coordinates_molecule(xyz_scale,xyz_temp(i_mole_acceptor,:,:), n_atom_temp(i_mole_acceptor), kk, K )
       do i=1,3
          xyz_scale_transpose(i,:) = xyz_scale(:,i)
       enddo
       ! call derivative_grid_Q routine with flag, this signifies don't store dQ_dr
       do i_atom=1, n_atom_temp(i_mole_acceptor)
          call derivative_grid_Q(force,theta_conv_Q_local,chg_temp(i_mole_acceptor,:),xyz_scale_transpose,i_atom, K, box, n, kk, flag_update)
          ! subtract this contribution
          pme_force_recip_diabat(i_mole_acceptor,i_atom,:) = pme_force_recip_diabat(i_mole_acceptor,i_atom,:) - force(:)
       enddo



       !******************** first change topology information for acceptor
       call evb_change_data_structures_proton_transfer( i_mole_donor, i_atom_donor, i_mole_acceptor, i_atom_acceptor, i_heavy_acceptor, xyz_temp , r_com_junk, chg_temp, n_atom_temp, mass_junk, box, pme_force_recip_diabat, atom_index_temp, molecule_index_temp, hydronium_molecule_index_temp )



       !******************** now contribution to force from acceptor diabat data structure topology

       ! add the contribution of the donor as an h2o molecule
       call create_scaled_direct_coordinates_molecule(xyz_scale,xyz_temp(i_mole_donor,:,:), n_atom_temp(i_mole_donor), kk, K )
       do i=1,3
          xyz_scale_transpose(i,:) = xyz_scale(:,i)
       enddo
       ! call derivative_grid_Q routine with flag, this signifies don't store dQ_dr
       do i_atom=1, n_atom_temp(i_mole_donor)
          call derivative_grid_Q(force,theta_conv_Q_local,chg_temp(i_mole_donor,:),xyz_scale_transpose,i_atom, K, box, n, kk, flag_update)
          ! add this contribution
          pme_force_recip_diabat(i_mole_donor,i_atom,:) = pme_force_recip_diabat(i_mole_donor,i_atom,:) + force(:)
       enddo

       ! add the contribution of the acceptor as an h3o+ molecule
       call create_scaled_direct_coordinates_molecule(xyz_scale,xyz_temp(i_mole_acceptor,:,:), n_atom_temp(i_mole_acceptor), kk, K )
       do i=1,3
          xyz_scale_transpose(i,:) = xyz_scale(:,i)
       enddo
       ! call derivative_grid_Q routine with flag, this signifies don't store dQ_dr
       do i_atom=1, n_atom_temp(i_mole_acceptor)
          call derivative_grid_Q(force,theta_conv_Q_local,chg_temp(i_mole_acceptor,:),xyz_scale_transpose,i_atom, K, box, n, kk, flag_update)
          ! add this contribution
          pme_force_recip_diabat(i_mole_acceptor,i_atom,:) = pme_force_recip_diabat(i_mole_acceptor,i_atom,:) + force(:)
       enddo

    enddo loop1


    !***************************************** now map forces back to principle diabat topology
    ! note here n_atom array should be in acceptor diabat topology, that's why we use n_atom_temp

    i_mole_acceptor = i_mole_principle
    ! loop through all proton hops.  A negative integer signals the end of the chain
    do i_hop=1, evb_max_chain
       if ( evb_diabat_proton_log( i_diabat , i_hop, 1 ) < 0 ) exit
       i_mole_donor = i_mole_acceptor
       i_atom_donor = evb_diabat_proton_log( i_diabat , i_hop , 2 )
       i_mole_acceptor = evb_diabat_proton_log( i_diabat , i_hop, 4 )

       ! fix the forces for this donor.  Note that the donated proton is always the last
       ! atom of the acceptor molecule
       count=0 ! this value is right, different from one above
       do i = i_atom_donor, n_atom_temp(i_mole_donor)
          ! start at end of array so we overwrite after we copy
          i_atom = n_atom_temp(i_mole_donor) - count
          ! shift forces up
          pme_force_recip_diabat( i_mole_donor, i_atom + 1 , : ) = pme_force_recip_diabat( i_mole_donor, i_atom  , : )
          count=count+1
       enddo
       ! now copy force from acceptor molecule back to donor
       pme_force_recip_diabat( i_mole_donor, i_atom_donor , : ) = pme_force_recip_diabat( i_mole_acceptor , n_atom_temp(i_mole_acceptor) , : )

    end do

    deallocate( xyz_temp , chg_temp , n_atom_temp , atom_index_temp )


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
  subroutine ms_evb_intermolecular_repulsion( force, E_ms_evb_repulsion, n_mole , n_atom , xyz, hydronium_molecule_index_local, atom_index_local, molecule_index_local, box )
    real*8,dimension(:,:,:),intent(inout) :: force
    real*8, intent(out)  :: E_ms_evb_repulsion
    integer, intent(in)  :: n_mole
    integer, dimension(:),intent(in) :: n_atom
    real*8, dimension(:,:,:), intent(in) :: xyz
    integer, dimension(:), intent(in) :: hydronium_molecule_index_local
    integer, dimension(:,:), intent(in) :: atom_index_local
    integer, dimension(:), intent(in) :: molecule_index_local
    real*8, dimension(:,:), intent(in)  :: box

    integer :: i_mole, i_mole_hydronium 

    E_ms_evb_repulsion = 0d0
    ! loop over hydronium molecules
    do i_mole=1, n_hydronium_molecules
       i_mole_hydronium = hydronium_molecule_index_local( i_mole )
       ! compute three-atom-special evb repulsion
       call ms_evb_three_atom_repulsion( force, E_ms_evb_repulsion, n_mole , n_atom , xyz, i_mole_hydronium, atom_index_local, molecule_index_local, box )
       ! now compute Born-Mayer terms
       call ms_evb_born_mayer( force, E_ms_evb_repulsion, n_mole , n_atom , xyz, i_mole_hydronium, atom_index_local, molecule_index_local, box )
       ! compute anisotropic morse interactions
       !call ms_evb_anisotropic_morse( force, E_ms_evb_repulsion, n_mole , n_atom , xyz, i_mole_hydronium, atom_index_local, molecule_index_local, box )
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
  subroutine ms_evb_three_atom_repulsion( force, E_ms_evb_repulsion, n_mole , n_atom , xyz, i_mole_hydronium, atom_index_local, molecule_index_local, box )
    real*8,dimension(:,:,:),intent(inout) :: force
    real*8, intent(inout)  :: E_ms_evb_repulsion
    integer, intent(in)  :: n_mole, i_mole_hydronium
    integer, dimension(:),intent(in) :: n_atom
    real*8, dimension(:,:,:), intent(in) :: xyz
    integer, dimension(:,:), intent(in) :: atom_index_local
    integer, dimension(:), intent(in) :: molecule_index_local
    real*8, dimension(:,:), intent(in)  :: box

    integer :: j_mole, i_atom, j_atom, index, i_type, i_type_H, i_type_heavy, j_type, i_atom_oxygen, i_heavy, index1
    real*8, dimension(3)  :: shift, rij, rij_O, q, fij
    real*8  :: r_ij, sum, q2, r_OO, fac_OO, fac_OH, exp_q, switch_OO, dswitch_OO
    ! interaction parameters
    ! these are for the heavy-atom/heavy-atom interaction
    real*8 ::  B, bl, d0_heavy, blprime, rs_heavy, rc_heavy


    i_type_H = atom_index_local(i_mole_hydronium,n_atom(i_mole_hydronium))
    call get_heavy_atom_transfer_acid( i_heavy , i_mole_hydronium, molecule_index_local )
    i_type_heavy = atom_index_local( i_mole_hydronium, i_heavy )

    ! loop over solvent molecules
    do j_mole=1, n_mole
       if ( j_mole /= i_mole_hydronium ) then
          ! see if any of the atoms of this solvent molecule are involved in the repulsive interaction
          do j_atom = 1 , n_atom(j_mole)
             j_type = atom_index_local(j_mole,j_atom)
             call get_index_atom_set( index1, evb_donor_acceptor_interaction, (/j_type, i_type_heavy, i_type_H/) )
             if ( index1 > 0 ) then
                ! get interaction parameters
                B        = evb_donor_acceptor_parameters(index1,1) ; bl       = evb_donor_acceptor_parameters(index1,2)
                d0_heavy = evb_donor_acceptor_parameters(index1,3) ; blprime  = evb_donor_acceptor_parameters(index1,4)
                rs_heavy = evb_donor_acceptor_parameters(index1,5) ; rc_heavy = evb_donor_acceptor_parameters(index1,6)

                shift = pbc_shift( xyz(i_mole_hydronium,i_heavy,:), xyz(j_mole,j_atom,:), box, xyz_to_box_transform)
                rij_O = -pbc_dr( xyz(i_mole_hydronium,i_heavy,:), xyz(j_mole,j_atom,:), shift )
                r_OO = dsqrt( dot_product( rij_O , rij_O ) )

                ! switching function
                call ms_evb_repulsive_switch( switch_OO, dswitch_OO, r_OO, rs_heavy , rc_heavy )
                fac_OO = B * exp(-bl *( r_OO - d0_heavy ))

                sum=0d0
                ! in this loop, calculate the 
                ! sum over q coordinates, defined by (rO + rO)/2 - rH, for the oxygen-oxygen repulsion
                do i_atom=1, n_atom(i_mole_hydronium)
                   i_type = atom_index_local(i_mole_hydronium,i_atom)
                   if ( i_type == i_type_H ) then

                      ! hydronium proton, the negative sign creates displacement from water to hydronium atom
                      rij = -pbc_dr( xyz(i_mole_hydronium,i_atom,:), xyz(j_mole,j_atom,:), shift ) 
                      r_ij = dsqrt( dot_product( rij, rij ) )

                      ! THIS IS CORRECT equation. See errata, JPC B, 2008, 112, 7146.  Original paper
                      ! had R_HO distance instead of q for the oxygen-oxygen repulsion
                      q = ( 2d0 * xyz(j_mole,j_atom,:) + rij_O ) / 2d0 - ( xyz(j_mole,j_atom,:) + rij )
                      q2 = dot_product(q,q)
                      exp_q = exp(-blprime * q2)
                      sum = sum + exp_q

                      ! note extra negative sign for dq/drH
                      force(i_mole_hydronium, i_atom,:) = force(i_mole_hydronium, i_atom,:) + switch_OO * fac_OO * exp_q * -blprime * 2d0 * q 
                      force(i_mole_hydronium, i_heavy,:) = force(i_mole_hydronium, i_heavy,:) + switch_OO * fac_OO * exp_q * blprime * q
                      force(j_mole, j_atom,:) = force(j_mole,j_atom,:) + switch_OO * fac_OO * exp_q * blprime * q

                   end if
                end do
                ! now oxygen-oxygen interaction, which depends on hydronium hydrogen positions    
                E_ms_evb_repulsion = E_ms_evb_repulsion + switch_OO * fac_OO * sum
                ! now derivatives of switch_OO * fac_OO
                fij = rij_O / r_OO * fac_OO * sum * ( switch_OO * bl  - dswitch_OO )
                force(i_mole_hydronium, i_heavy,:) = force(i_mole_hydronium, i_heavy,:) + fij             
                force(j_mole, j_atom,:) = force(j_mole,j_atom,:) - fij
             end if
          end do
       end if
    end do



  end subroutine ms_evb_three_atom_repulsion



  !********************************
  ! these are general Born-Mayer interactions between
  ! atoms on the proton-donor and acceptor
  !********************************
  subroutine ms_evb_born_mayer( force, E_ms_evb_repulsion, n_mole , n_atom , xyz, i_mole_hydronium, atom_index_local, molecule_index_local, box )
    real*8,dimension(:,:,:),intent(inout) :: force
    real*8, intent(inout)  :: E_ms_evb_repulsion
    integer, intent(in)  :: n_mole, i_mole_hydronium
    integer, dimension(:),intent(in) :: n_atom
    real*8, dimension(:,:,:), intent(in) :: xyz
    integer, dimension(:,:), intent(in) :: atom_index_local
    integer, dimension(:), intent(in) :: molecule_index_local
    real*8, dimension(:,:), intent(in)  :: box

    integer ::  j_mole, i_atom, j_atom, i_type, j_type,  index1
    real*8, dimension(3)  :: shift, rij, fij
    real*8  :: r_ij, fac_OH, switch_HO,dswitch_HO
    ! interaction parameters
    ! these are for the hydrogen/heavy-atom interaction
    real*8 ::  C, cl, d0_hyd , rs_hyd , rc_hyd


    ! loop over atoms in proton donor
    do i_atom = 1 , n_atom(i_mole_hydronium)
       i_type = atom_index_local(i_mole_hydronium,i_atom)

       ! loop over solvent molecules
       do j_mole=1, n_mole
          if ( j_mole /= i_mole_hydronium ) then
             ! see if any of the atoms of this solvent molecule are involved in the repulsive interaction
             do j_atom = 1 , n_atom(j_mole)
                j_type = atom_index_local(j_mole,j_atom)
                call get_index_atom_set( index1, evb_proton_acceptor_interaction, (/j_type, i_type/) )
                if ( index1 > 0 ) then

                   C       = evb_proton_acceptor_parameters(index1,1) ; cl       = evb_proton_acceptor_parameters(index1,2) 
                   d0_hyd  = evb_proton_acceptor_parameters(index1,3) ; rs_hyd   = evb_proton_acceptor_parameters(index1,4) 
                   rc_hyd  = evb_proton_acceptor_parameters(index1,5) 

                   shift = pbc_shift( xyz(i_mole_hydronium,i_atom,:), xyz(j_mole,j_atom,:), box, xyz_to_box_transform)

                   ! hydronium proton, the negative sign creates displacement from water to hydronium atom
                   rij = -pbc_dr( xyz(i_mole_hydronium,i_atom,:), xyz(j_mole,j_atom,:), shift ) 
                   r_ij = dsqrt( dot_product( rij, rij ) )

                   fac_OH = C * exp( -cl * ( r_ij - d0_hyd ) )
                   ! switching function
                   call ms_evb_repulsive_switch( switch_HO, dswitch_HO, r_ij, rs_hyd , rc_hyd )

                   fij = rij / r_ij * fac_OH * ( switch_HO * cl  - dswitch_HO )
                   E_ms_evb_repulsion = E_ms_evb_repulsion + switch_HO * fac_OH
                   force(i_mole_hydronium, i_atom,:) = force(i_mole_hydronium, i_atom,:) + fij(:)
                   force(j_mole,j_atom,:) = force(j_mole,j_atom,:) - fij(:)

                end if
             enddo
          endif
       enddo
    enddo


  end subroutine ms_evb_born_mayer



  !********************************
  ! these are anisotropic morse interactions
  ! atoms on the proton-donor and acceptor
  !********************************
  subroutine ms_evb_anisotropic_morse( force, E_ms_evb_repulsion, n_mole , n_atom , xyz, i_mole_hydronium, atom_index_local, molecule_index_local, box )
    real*8,dimension(:,:,:),intent(inout) :: force
    real*8, intent(inout)  :: E_ms_evb_repulsion
    integer, intent(in)  :: n_mole, i_mole_hydronium
    integer, dimension(:),intent(in) :: n_atom
    real*8, dimension(:,:,:), intent(in) :: xyz
    integer, dimension(:,:), intent(in) :: atom_index_local
    integer, dimension(:), intent(in) :: molecule_index_local
    real*8, dimension(:,:), intent(in)  :: box

    integer ::  j_mole, i_atom, j_atom, k_atom, i_type, j_type, k_type, index1, j_mole_type
    real*8, dimension(3)  :: shift, rij, rjk, fij, fjk
    real*8  :: r_ij, r_jk, cosine, E_morse, morse_term1,  morse_term2
    ! interaction parameters
    real*8 ::  D0, D1, alpha, req


    ! loop over atoms in proton donor
    do i_atom = 1 , n_atom(i_mole_hydronium)
       i_type = atom_index_local(i_mole_hydronium,i_atom)

       ! loop over solvent molecules
       do j_mole=1, n_mole
          j_mole_type = molecule_index_local(j_mole)

          if ( j_mole /= i_mole_hydronium ) then
             ! see if any of the atoms of this solvent molecule are involved in the repulsive interaction
             do j_atom = 1 , n_atom(j_mole)
                j_type = atom_index_local(j_mole,j_atom)
                ! find axis defining atom. This atom must be bonded
                do k_atom = 1 , n_atom(j_mole)
                   k_type = atom_index_local(j_mole,k_atom)

                   if ( molecule_bond_list(j_mole_type,j_atom,k_atom) == 1 ) then

                      ! see if this is a defined interaction
                      call get_index_atom_set( index1, evb_anisotropic_morse_interaction, (/i_type, j_type, k_type/) )

                      if ( index1 > 0 ) then

                         D0 =    evb_anisotropic_morse_parameters(index1,1) ; D1  = evb_anisotropic_morse_parameters(index1,2) ;
                         alpha = evb_anisotropic_morse_parameters(index1,3) ; req = evb_anisotropic_morse_parameters(index1,4) ;

                         shift = pbc_shift( xyz(i_mole_hydronium,i_atom,:), xyz(j_mole,j_atom,:), box, xyz_to_box_transform)

                         ! hydronium proton, the negative sign creates displacement from water to hydronium atom
                         rij = -pbc_dr( xyz(i_mole_hydronium,i_atom,:), xyz(j_mole,j_atom,:), shift ) 
                         r_ij = dsqrt( dot_product( rij, rij ) )
                         ! this defines z-axis
                         rjk(:) = xyz(j_mole,j_atom,:) - xyz(j_mole,k_atom,:)
                         r_jk = dsqrt( dot_product( rjk, rjk ) )
                         cosine = dot_product( rij , rjk ) / r_ij / r_jk

                         morse_term1 = D0 + D1 * cosine
                         morse_term2 = 1d0 - ( 1d0 - exp(-alpha*(r_ij - req)) ) ** 2 

                         E_morse = morse_term1 * morse_term2
                         E_ms_evb_repulsion = E_ms_evb_repulsion + E_morse

                         ! from derivative of radial morse part
                         fij = 2d0 * morse_term1 * alpha * exp(-alpha*(r_ij - req)) * ( 1d0 - exp(-alpha*(r_ij - req)) )* rij(:) / r_ij
                         ! from derivative of anistropic coeff part

                         fij = fij - ( rjk / r_ij / r_jk - cosine * rij(:) / r_ij ** 2 ) * D1 * morse_term2
                         fjk =     - ( rij / r_ij / r_jk - cosine * rjk(:) / r_jk ** 2 ) * D1 * morse_term2

                         force(i_mole_hydronium, i_atom,:) = force(i_mole_hydronium, i_atom,:) + fij(:)
                         force(j_mole,j_atom,:) = force(j_mole,j_atom,:) - fij(:) + fjk(:)
                         force(j_mole,k_atom,:) = force(j_mole,k_atom,:) - fjk(:)

                      end if

                   end if

                enddo
             enddo
          endif
       enddo
    enddo


  end subroutine ms_evb_anisotropic_morse





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
  ! force_atoms array will correspond to a different xyz topology (i.e. the protons
  ! will correspond to different molecules)
  !
  ! therefore, we need to map the force_atoms array back to the principle
  ! xyz topology of the principle diabat, and store all the forces
  ! in this same, consistent topology
  !********************************* 
  subroutine evb_store_forces( i_mole_principle, diabat_index1, diabat_index2, tot_n_mole, n_atom, force_atoms, store_index, initialize )
    integer, intent(in) :: i_mole_principle, diabat_index1, diabat_index2, tot_n_mole
    integer, dimension(:), intent(in) :: n_atom
    real*8,dimension(:,:,:), intent(in) :: force_atoms
    integer, intent(inout)  :: store_index
    integer, intent(in),optional :: initialize

    integer, save :: size
    real*8, dimension(:,:,:), allocatable :: force_atoms_map
    integer :: i_mole, i_atom, i_hop, i_mole_donor, i_atom_donor, i_mole_acceptor, i_atom_acceptor, diabat_index

    ! initialize, allocate evb_force_store array.  Here we guess as to how large an array we need to allocate, which is based on the number of diabats plus the number of non-zero diabat coupling matrix elements.  The couplings are only non-zero for diabats that have hydronium on adjacent water molecules, so this corresponds to the number of donor hydrogen bonds per water molecule.  The number of water molecules is strictly less than the number of diabats. so 2 * n_diabats is a safe guess
    if ( present(initialize) ) then
       size = 3 * evb_max_states
       ! we call this subroutine with initialize present to reset the indices, so this array may already be allocated
       if ( .not. allocated(evb_forces_store) ) then
          ! this is a big array, so allocate it for efficient looping for column major memory storage
          allocate( evb_forces_store(3,MAX_N_ATOM,tot_n_mole,size) )
          ! no need in zeroing this array as it is expensive, since it's so big (most of it is excessive anyway)
          !evb_forces_store=0d0
       end if
       ! initialize evb_forces_lookup_index to negative value, so that we know which couplings have non-zero force
       evb_forces_lookup_index=-1
    end if

    ! store mapping
    evb_forces_lookup_index( diabat_index1 , diabat_index2 ) = store_index


    ! allocate this array in order for efficient memory retrieval
    allocate( force_atoms_map(3,MAX_N_ATOM,tot_n_mole) )
    ! copy forces.  Note force_atoms array is typically of size (MAX_N_MOLE, MAX_N_ATOM,3) , with 
    ! MAX_N_MOLE > tot_n_mole, so we loop to copy elements
    do i_mole=1, tot_n_mole
       do i_atom=1,n_atom(i_mole)
          force_atoms_map(:,i_atom, i_mole) = force_atoms(i_mole, i_atom , : )
       end do
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

       ! loop through all proton hops.  A negative integer signals the end of the chain
       do i_hop=1, evb_max_chain
          if ( evb_diabat_proton_log( diabat_index , i_hop, 1 ) < 0 ) exit
          i_mole_donor = i_mole_acceptor
          i_atom_donor = evb_diabat_proton_log( diabat_index , i_hop, 2 )
          i_mole_acceptor = evb_diabat_proton_log( diabat_index , i_hop, 4 )

          ! fix the forces for this donor.  Note that the donated proton is always the last
          ! atom of the acceptor molecule
          do i_atom = i_atom_donor, n_atom(i_mole_donor)
             ! shift forces up, not we're not overwriting here because different arrays,
             force_atoms_map( :, i_atom + 1 , i_mole_donor ) = force_atoms( i_mole_donor, i_atom  , : )
          enddo
          ! now copy force from acceptor molecule back to donor
          force_atoms_map( :, i_atom_donor , i_mole_donor ) = force_atoms( i_mole_acceptor , n_atom(i_mole_acceptor) , : )
       end do
    end if

    evb_forces_store( :, : , : , store_index) = force_atoms_map(:,:,:)

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



  !**************************
  ! this subroutine retrieves the reference energy (chemical energy) of
  ! each adiabatic state.  Currently, this is stored for each acid, 
  ! assuming the contribution from the basic molecule is included.  If 
  ! the acid does not uniquely determine the system state, this will have
  ! to be generalized
  !**************************
  subroutine get_adiabatic_reference_energy( E_reference, i_mole_donor, molecule_index_local )
    real*8,intent(out) :: E_reference
    integer, intent(in) :: i_mole_donor
    integer, dimension(:), intent(in) :: molecule_index_local

    integer :: i_type

    i_type = molecule_index_local( i_mole_donor )
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
  subroutine get_heavy_atom_transfer_base( i_atom_base, i_mole_base, molecule_index_local )
    integer, intent(out) :: i_atom_base
    integer, intent(in) :: i_mole_base
    integer, dimension(:), intent(in) :: molecule_index_local

    integer :: i_type_base, i_type_acid, i_type_heavy, i_atom

    i_type_base = molecule_index_local(i_mole_base)
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
  subroutine get_heavy_atom_transfer_acid( i_atom_acid, i_mole_acid, molecule_index_local )
    integer, intent(out) :: i_atom_acid
    integer, intent(in) :: i_mole_acid
    integer, dimension(:), intent(in) :: molecule_index_local

    integer :: i_type_heavy, i_type_acid, i_atom

    i_type_acid = molecule_index_local(i_mole_acid)
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
  function zundel_r_com( i_mole_donor, i_mole_acceptor, r_com_donor , r_com_acceptor , shiftd, shifta ,  box , n_atom, atom_index_local )
    real*8,dimension(3) :: zundel_r_com
    integer, intent(in) :: i_mole_donor, i_mole_acceptor
    real*8,dimension(3),intent(in) :: r_com_donor, r_com_acceptor
    real*8, dimension(3), intent(out) :: shiftd, shifta
    real*8,dimension(3,3),intent(in) :: box
    integer, dimension(:),intent(in) :: n_atom
    integer, dimension(:,:), intent(in) :: atom_index_local

    real*8 :: mass_donor, mass_acceptor
    real*8,dimension(3) :: rda, shift, r_com_a
    integer :: i_atom, i_type

    mass_donor=0d0
    mass_acceptor=0d0

    do i_atom =1, n_atom(i_mole_donor)
       i_type = atom_index_local(i_mole_donor,i_atom)
       mass_donor = mass_donor + atype_mass(i_type)
    enddo

    do i_atom =1, n_atom(i_mole_acceptor)
       i_type = atom_index_local(i_mole_acceptor,i_atom)
       mass_acceptor = mass_acceptor + atype_mass(i_type)
    enddo

    ! donor and acceptor may be broken up over pbc
    shift(:) = pbc_shift( r_com_donor , r_com_acceptor, box , xyz_to_box_transform )
    rda(:)  =   pbc_dr( r_com_donor , r_com_acceptor, shift )

    ! shifted acceptor position
    r_com_a = r_com_donor + rda
    zundel_r_com(:) = ( mass_donor * r_com_donor(:) + mass_acceptor * r_com_a(:) ) / ( mass_donor + mass_acceptor )

    ! shift donor is zero by above definition
    shiftd=0d0
    shifta=shift  

  end function zundel_r_com




  !**************************************
  ! this subroutine changes the atom index of
  ! the transferring proton to the acceptor molecule type
  ! when this subroutine is called, the acceptor molecule is still in its
  ! basic topology, so molecule_index_local will return the index of the base
  !**************************************
  subroutine  change_proton_index_proton_transfer( i_mole_acceptor, i_atom_acceptor , molecule_index_local, atom_index_local )
    integer, intent(in) :: i_mole_acceptor, i_atom_acceptor
    integer, dimension(:), intent(in)  :: molecule_index_local
    integer, dimension(:,:),intent(inout) :: atom_index_local

    integer :: i_base_type, i_acid_type

    ! this gives base type, see above comment
    i_base_type = molecule_index_local(i_mole_acceptor)
    ! get conjugate acid type
    i_acid_type = evb_conjugate_pairs( i_base_type )

    ! fill in proton atom index of acceptor with appropriate proton index
    atom_index_local(i_mole_acceptor,i_atom_acceptor) = evb_proton_index( i_acid_type )

  end subroutine change_proton_index_proton_transfer






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
  subroutine print_evb_trajectory_data( ground_state_eigenvector , i_mole_principle, log_file )
    real*8,dimension(:), intent(in) :: ground_state_eigenvector
    integer, intent(in) :: i_mole_principle, log_file

    integer :: i_state,i_mole_donor, i_mole_acceptor, i_atom_donor, i_hop,  shell
    real*8  :: coefficient

    write( log_file, * ) "number of diabat states : ", diabat_index
    write( log_file, * ) "diabat state    hydronium molecule   evb coefficient  solvation shell"

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

       write( log_file, '(I5,I10,F14.6,I5)' ) i_state, i_mole_acceptor, coefficient, shell

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


    !**************************  anisotropic morse interaction
!!$    call read_file_find_heading( file_h, '[ anisotropic ]' , flag_anisotropic, flag_eof )
!!$    if ( flag_eof == -1 ) goto 100  ! end of file
!!$    call write_ms_evb_io_info( '[ anisotropic ]' )
!!$    flag_eof=0
!!$    count=1
!!$    do
!!$       call read_topology_line( file_h , input_string , flag )
!!$       ! if end of file
!!$       if ( flag == -1 ) then
!!$          flag_eof=-1
!!$          exit
!!$       end if
!!$       ! if end of section
!!$       if ( flag == 1 ) exit
!!$
!!$
!!$       if ( count > max_interaction_type ) stop "please increase value of 'max_interaction_type'"
!!$       call parse(input_string," ",args,nargs)
!!$       if ( nargs /= 7 ) stop "must have 7 arguments in 'anisotropic' section of [ evb_parameters ] section of topology file"
!!$       read(args(1),*) atomtype1
!!$       read(args(2),*) atomtype2
!!$       read(args(3),*) atomtype3
!!$       call trim_end(atomtype1)
!!$       call trim_end(atomtype2)
!!$       call trim_end(atomtype3)
!!$       call atype_name_reverse_lookup( atomtype1, itype1 )
!!$       call atype_name_reverse_lookup( atomtype2, itype2 )
!!$       call atype_name_reverse_lookup( atomtype3, itype3 )
!!$       ! fill in index for this interaction
!!$       evb_anisotropic_morse_interaction(count,1) = itype1
!!$       evb_anisotropic_morse_interaction(count,2) = itype2  
!!$       evb_anisotropic_morse_interaction(count,3) = itype3 
!!$       read(args(4),*) evb_anisotropic_morse_parameters(count,1)
!!$       read(args(5),*) evb_anisotropic_morse_parameters(count,2)
!!$       read(args(6),*) evb_anisotropic_morse_parameters(count,3)
!!$       read(args(7),*) evb_anisotropic_morse_parameters(count,4)
!!$       count = count + 1
!!$    enddo
!!$    if ( flag_eof == -1 ) goto 100  ! end of file



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
