!***********************************************************************
!
!  This module contains various initialization subroutines
!
!***********************************************************************
module initialize_routines
  use routines
  use ms_evb
  implicit none

contains

  !*************************************************************************
  ! this routine  initializes the simulation
  !*************************************************************************
  subroutine initialize_simulation(system_data, molecule_data, atom_data, integrator_data, file_io_data, verlet_list_data, PME_data, xyz, velocity, force, mass, charge, atom_type_index, aname )
    use global_variables
    use pairwise_interaction
    use pme_routines
    use bonded_interactions
    Type(system_data_type),intent(inout)                :: system_data
    Type(molecule_data_type),dimension(:),intent(inout) :: molecule_data
    Type(atom_data_type),intent(inout)                  :: atom_data
    Type(integrator_data_type),intent(inout)            :: integrator_data
    Type(file_io_data_type),intent(inout)               :: file_io_data
    Type(verlet_list_data_type),intent(inout)           :: verlet_list_data
    Type(PME_data_type), intent(inout)                  :: PME_data

    real*8, dimension(:,:), intent(inout) :: xyz
    real*8, dimension(:,:), intent(inout) :: velocity
    real*8, dimension(:,:), intent(inout) :: force
    real*8, dimension(:),   intent(inout) :: mass
    real*8, dimension(:),   intent(inout) :: charge
    integer, dimension(:),  intent(inout) :: atom_type_index  
    character(*)(:), intent(inout)        :: aname

    !******** this is a local data structure with pointers that will be set
    ! to subarrays of atom_data arrays for the specific atoms in the molecule
    type(single_molecule_data_type) :: single_molecule_data

    integer, dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE) :: gen_cross_terms
    real*8,dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE,9) :: temp_lj_parameter
    integer:: i_mole, total_atoms, size


    !********** this is for generating random velocities for Maxwell-Boltzmann
    call initialize_random_seed

    ! **************************** read atom coordinates   ***********************************!

    ! first get number of molecules, allocate datastructures
    call read_gro_number_molecules( file_io_data%ifile_gro_file_h, file_io_data%ifile_gro, system_data%n_mole )

    ! allocate molecule_data array
    allocate( molecule_data(system_data%n_mole) )

    ! now read gro file
    ! this call will fill in data structures atom_data and molecule_data
    ! open gro file as we may need to scan to end if restarting...
    open( file_io_data%ifile_gro_file_h, file=file_io_data%ifile_gro, status='old' )

    ! if trajectory restart, scan to end of .gro trajectory file...
    Select Case( restart_trajectory )
    Case("yes")
       call scan_grofile_restart( file_io_data%ifile_gro_file_h, n_old_trajectory )
    Case default

    ! this will allocate aname, xyz arrays, and attach pointers from atom_data structure
    call read_gro( file_io_data%ifile_gro_file_h, system_data, molecule_data, atom_data, aname, xyz )

    ! ******************  Allocate remainder of atom_data target arrays, initialize, and attach pointers
    total_atoms = system_data%total_atoms
    allocate( velocity(3,total_atoms) )
    allocate( force(3,total_atoms) )
    allocate( mass(total_atoms) )
    allocate( charge(total_atoms) )
    allocate( atom_type_index(total_atoms) )

    velocity=0d0;force=0d0;mass=0d0;charge=0d0; atom_type_index=0
    atom_data%velocity=>velocity
    atom_data%force=>force
    atom_data%mass=>mass
    atom_data%charge=>charge
    atom_data%atom_type_index=>atom_type_index


    ! if restarting, read in velocities from restart file
    Select Case( restart_trajectory )
    Case("yes")
       call read_velocity_restart_checkpoint( file_io_data%ifile_velocity_file_h, file_io_data%ifile_velocity, atom_data%velocity, n_old_trajectory )
    End Select


    ! initialize the transformation matrix to box coordinates
    call initialize_non_orth_transform ( system_data%box, system_data%xyz_to_box_transform )

    !***************** make sure molecules are not broken up by a box translation, as this can happen in GROMACs, and could be using GROMACS output as input here
    do i_mole=1,system_data%n_mole
       ! set pointers for this data structure to target molecule
       ! we will be just using xyz coordinates here
       ! note we are changing global atom_data%xyz data structure with pointer !
       call return_molecule_block( single_molecule_data , atom_data , molecule_data(i_mole)%n_atom, molecule_data(i_mole)%atom_index )
       call fix_intra_molecular_shifts( molecule_data(i_mole)%n_atom, single_molecule_data%xyz , system_data%box,  system_data%xyz_to_box_transform )
    enddo


    !***************** make sure cutoff distances are fine for box size, and
    !that box type is supported
    call check_cutoffs_box( real_space_cutoff, verlet_list_data%verlet_cutoff, system_data%box )


    !************************************** get parameters****************************************************!
    ! get parameters for all atom types, assuming there are no more than MAX_N_ATOM_TYPE types of atoms total
    !*****************************************************************************************************!
    call read_param( file_io_data_type%ifile_ffpmt_file_h, file_io_data_type%ifile_ffpmt,gen_cross_terms)
    ! make sure we don't have too many atom types
    if ( n_atom_type > MAX_N_ATOM_TYPE ) then
       stop "number of atom types g.t. MAX_N_ATOM_TYPE.  Increase this value"
    endif


    ! ********************* generate cross terms and fill in parameter arrays****************************!
    ! fill in total parameter arrays that are stored in global_variables
    call gen_param( atom_data%aname, system_data%total_atoms, atom_data%charge, atom_data%atom_type_index, gen_cross_terms)

    ! now get special 1-4 interactions from parameter file if needed
    call read_generate_14_interaction_parameters( file_io_data%ifile_ffpmt )

    ! sort molecules into types
    call gen_molecule_type_data( system_data%n_mole, molecule_data, atom_data )

    ! initialize mass to negative values to make sure we read all of these in
    atype_mass=-1d0

    ! read in topology information, zero molecule_exclusions here, as it is possible that
    ! exclusions are being read in from topology file
    molecule_exclusions=0
    call read_topology_file( file_io_data%ifile_top_file_h, file_io_data%ifile_top )

    ! mass is read in from topology file, fill in mass for each molecule
    call fill_mass( system_data%total_atoms, atom_data%mass, atom_data%atom_type_index )


    ! center of mass
    call update_r_com( system_data%n_mole, molecule_data, atom_data )
    ! center of mass of molecule might be outside of box after calling subroutine check_intra_molecular_shifts, fix this
    do i_mole=1,n_mole
       ! set pointers for this data structure to target molecule
       ! we will be just using xyz coordinates here
       call return_molecule_block( single_molecule_data , atom_data , molecule_data(i_mole)%n_atom, molecule_data(i_mole)%atom_index )
       call shift_move(molecule_data(i_mole)%n_atom,single_molecule_data%xyz,molecule_data(i_mole)%r_com,system_data%box, system_data%xyz_to_box_transform)
    enddo


    ! generate exclusions, note that some exclusions may have already been explicity read in
    ! from topology file.
    call generate_intramolecular_exclusions


    !*********************************initialize verlet list
    call allocate_verlet_list( verlet_list_data, system_data%total_atoms, system_data%volume )
    call construct_verlet_list( verlet_list_data, atom_data, system_data%total_atoms, system_data%box, system_data%xyz_to_box_transform  )
    ! the "1" input to update_verlet_displacements signals to initialize the displacement array
    call update_verlet_displacements( system_data%total_atoms, atom_data%xyz, verlet_list_data , system_data%box, system_data%xyz_to_box_transform , flag_junk, 1 )


    ! initialize ms_evb simulation
    Select Case(ms_evb_simulation)
    Case("yes")
       call initialize_evb_simulation( system_data%n_mole, file_io_data )
       ! we are storing dQ_dr grid for ms-evb, allocate data arrays for this

       ! grid size:  each atom has spline_order^3 non-zero derivatives in the Q grid
       size = system_data%total_atoms * spline_order**3
       ! size could be big depending on the system, make sure we don't ask for too much memory
       if ( size > 10000000 ) then
	  write(*,*) "are you sure you want to store dQ_dr?  This requires allocating an array"
	  write(*,*) " bigger than 3 x ", size
          stop
       end if
       ! allocate column major
       size=spline_order**3
       allocate( dQ_dr(3,size,system_data%total_atoms) , dQ_dr_index(3,size,system_data%total_atoms) )
       allocate(atom_list_force_recip(3,system_data%total_atoms) )
    End Select


  end subroutine initialize_simulation



  !**************************************************!
  ! this subroutine initializes energy and force calculations
  ! and computes energy and force for initial configuration
  ! 
  ! 3/2/2015 i/o attributes changed by JGM
  ! a call to the multi-state empirical valence bond module
  ! ( ms_evb_calculate_total_force_energy ) may change data structures
  ! xyz, r_com, chg, n_atom, mass, if the initially assigned 
  ! hydronium molecule is unstable w.r.t. a proton transfer.
  ! therefore these data structures were given the intent(inout) attribute
  !
  !**************************************************!
  subroutine initialize_energy_force(log_file,box,n_mole,tot_n_mole,n_atom,n_atom_drude,xyz,r_com,mass,chg,dfti_desc,dfti_desc_inv,iteration,potential,E_elec_nopol,E_elec,E_bh,E_3body,E_bond,E_angle,E_dihedral,xyz_drude,force_atoms)
    use global_variables
    use electrostatic
    use pairwise_interaction
    use pme_routines
    use MKL_DFTI
    use total_energy_forces
    integer,intent(in)::log_file
    real*8,dimension(:,:),intent(inout)::chg
    integer,intent(in)::n_mole,tot_n_mole
    integer,dimension(:),intent(inout)::n_atom,n_atom_drude
    real*8,dimension(:,:,:),intent(inout)::xyz
    real*8,dimension(:,:),intent(in) :: box
    real*8,dimension(:,:),intent(inout) :: r_com,mass
    TYPE(DFTI_DESCRIPTOR), pointer,intent(out):: dfti_desc,dfti_desc_inv
    integer,intent(out) :: iteration
    real*8,intent(out) :: potential, E_elec_nopol, E_elec, E_bh, E_3body,E_bond,E_angle,E_dihedral
    real*8, dimension(:,:,:),intent(out) :: xyz_drude,force_atoms

    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: lj_force
    integer,dimension(MAX_N_MOLE,MAX_N_ATOM)::drude_atoms
    integer, dimension(MAX_N_MOLE) :: n_atom_junk
    integer:: i,j,i_mole,i_atom,length(3),status
    real*8 :: a(3), b(3), c(3), ka(3), kb(3), kc(3),kk(3,3),vol
    real*8 :: x, exp_x2
    integer :: ntype_solute, atom_id1


    !************************************* initialize ewald/pme ************************************************!

    ! note here that tot_chg was passed to this subroutine, so drude oscillator charges are accounted for
    ! if they are present.  Use n_atom_drude for Ewald_self, so that drude oscillators are included
    call update_Ewald_self( tot_n_mole, n_atom_drude, chg )


    Select Case(electrostatic_type)
    Case("pme")

       ! allocate arrays
       allocate(CB(pme_grid,pme_grid,pme_grid),Q_grid(pme_grid,pme_grid,pme_grid),theta_conv_Q(pme_grid,pme_grid,pme_grid))

       ! set up fourier transform descriptors
       length=pme_grid

       status=DftiCreateDescriptor(dfti_desc, DFTI_DOUBLE, DFTI_COMPLEX, 3, length)
       status=DftiCommitDescriptor(dfti_desc)
       ! don't scale back transform because pick up factor of K^3 from convolution
       status=DftiCreateDescriptor(dfti_desc_inv, DFTI_DOUBLE, DFTI_COMPLEX, 3, length)
       status = DftiCommitDescriptor(dfti_desc_inv)


!!!!!!!!!!!!!!! compute CB array 
       a(:) = box(1,:);b(:) = box(2,:);c(:) = box(3,:)
       vol = volume( a, b, c )
       call crossproduct( a, b, kc ); kc = kc /vol 
       call crossproduct( b, c, ka ); ka = ka /vol
       call crossproduct( c, a, kb ); kb = kb /vol
       kk(1,:)=ka(:);kk(2,:)=kb(:);kk(3,:)=kc(:)

       call CB_array(CB,alpha_sqrt,vol,pme_grid,kk,spline_order)

!!!!!!!!!!!!!!!!grid B_splines
       if (spline_order .eq. 6) then
          do i=1,spline_grid
             B6_spline(i)=B_spline(6./dble(spline_grid)*dble(i),6)
             B5_spline(i)=B_spline(5./dble(spline_grid)*dble(i),5)
          enddo
       else if (spline_order .eq. 4) then
          do i=1,spline_grid
             B4_spline(i)=B_spline(4./dble(spline_grid)*dble(i),4)
             B3_spline(i)=B_spline(3./dble(spline_grid)*dble(i),3)
          enddo
       else
          stop "requested spline order not implemented"
       endif

       Select case(grid_erfc)
       case("yes")
!!!!!!!!!!!!!!!!grid complementary error function
          do i=1,erfc_grid
             erfc_table(i)= erfc(erfc_max/dble(erfc_grid)*dble(i))
          enddo
       end select

    end select

    ! initialize factorial array
    call initialize_factorial

    !************ grid Tang_Toennies C6-C12 damping functions if requested
    Select Case(grid_Tang_Toennies)
    Case("yes")
       do i = 1, Tang_Toennies_grid         
          x = Tang_Toennies_max/dble(Tang_Toennies_grid)*dble(i)
          ! C6 is in 1st index, C8 in 2nd, C10 in 3rd, C12 in 4th
          Tang_Toennies_table(1,i) = Tang_Toennies_damp(x,6)
          Tang_Toennies_table(2,i) = Tang_Toennies_damp(x,8)
          Tang_Toennies_table(3,i) = Tang_Toennies_damp(x,10)
          Tang_Toennies_table(4,i) = Tang_Toennies_damp(x,12)

          ! now derivatives
          dTang_Toennies_table(1,i) = dTang_Toennies_damp(x,6)
          dTang_Toennies_table(2,i) = dTang_Toennies_damp(x,8)
          dTang_Toennies_table(3,i) = dTang_Toennies_damp(x,10)
          dTang_Toennies_table(4,i) = dTang_Toennies_damp(x,12)
       enddo

       ! if three-body dispersion, grid 3rd order damping function
       Select Case(three_body_dispersion)
       Case("yes")
          do i = 1, Tang_Toennies_grid         
             x = Tang_Toennies_max/dble(Tang_Toennies_grid)*dble(i)
             Tang_Toennies_table3(i) = Tang_Toennies_damp(x,3)
          enddo
       End Select

    End Select


    !***************** initialize pme grid algorithm if using pme for dispersion ********************!

    if ( pme_disp == "yes" ) then
       !     the grid for long range dispersion interaction
       !     Count number of atom types in solute molecule
       stop " atype_solute data structure has been removed, and code below is commented "
!!$       do ntype_solute = 1, n_atom_type
!!$          if ( atype_solute(ntype_solute) == 0 ) then
!!$             goto 999
!!$          endif
!!$       enddo
!!$999    continue
!!$       ntype_solute = ntype_solute - 1
!!$
!!$       allocate( lrdisp_pot(ntype_solute,pme_grid,pme_grid,pme_grid) )
!!$       !     construct long-range dispersion grid
!!$       lrdisp_pot = 0d0
!!$       do atom_id1 = 1, ntype_solute
!!$          call construct_recip_disp_grid(atom_id1,n_mole,tot_n_mole,n_atom,box,xyz,dfti_desc,dfti_desc_inv)
!!$       enddo
    endif


    !  set long range correction and shifts
    ! this is in case there is no framework. If there is a framework, disp_lrc is automatically set to zero for obvious reasons


    !**************************** Get initial forces and energy *************************************!
    Select Case(ms_evb_simulation)
    Case("yes")
       call ms_evb_calculate_total_force_energy( force_atoms, potential, tot_n_mole, n_mole, n_atom, xyz, r_com, chg, mass, box, dfti_desc,dfti_desc_inv,log_file )
    Case("no")
       call calculate_total_force_energy( force_atoms,potential, E_elec,E_elec_nopol,E_bh, E_3body, E_bond,E_angle,E_dihedral,iteration, tot_n_mole, n_mole, n_atom, n_atom_drude, r_com, xyz, chg, box, dfti_desc,dfti_desc_inv,log_file,xyz_drude)
    End Select
    !************************************************************************************************!


    !**************** since this is initial energy calculation, it's important that drude oscillators converged, stop otherwise
    Select Case(drude_simulation)
    Case(1)
       if ( iteration > iteration_limit ) then
          stop "Drude oscillators did not converge in the initial energy calculation.  Check that variables force_threshold and iteration_limit have consistent settings in glob_v.f90."
       endif
    End Select



  end subroutine initialize_energy_force




  !**************************************************!
  ! this subroutine reads atom_type parameters from input file for a non SAPT-based force field
  !*************************************************
  subroutine read_param( file_h, ifile_pmt, gen_cross_terms)
    use global_variables
    integer, intent(in)    :: file_h
    character(*),intent(in):: ifile_pmt
    integer,dimension(:,:),intent(out) :: gen_cross_terms

    integer::i_type, j_type, i_param, n_cross, inputstatus,i_search,j_search,ind,ind1,ind2
    character(MAX_ANAME):: atom1, atom2
    integer,parameter :: max_param=20
    real*8,parameter :: small=1D-6
    character(20),dimension(max_param) :: args
    real*8,dimension(6) :: store
    integer        :: nargs
    character(300) :: input_string
    character(25)::junk
    character(50)::line


    gen_cross_terms=0

    ! store has dimension 6, since atype_lj_parameter has dimension 6 ( changed to accomodate C12).  Zero these components that will not be used
    store=0d0

    open(unit=file_h,file=ifile_pmt,status="old")

    ! get solute atom types
    do
       Read(file_h,'(A)',Iostat=inputstatus) line
       If(inputstatus < 0) Exit
       ind=INDEX(line,'solute_species')
       IF(ind .NE. 0) Exit
    enddo
    Read(file_h,'(A)',Iostat=inputstatus) line
    read(file_h,*) n_atom_type

    do i_type=1, n_atom_type
       ! new input format, we assume lj force field (not buckingham), and read in an integer to decide if that atom type will be frozen
       read(file_h,'(A)') input_string
       call parse(input_string," ",args,nargs)
       if ( nargs /= 5 ) then
          write(*,*) ""
          write(*,*) "should have 5 input arguments under solute_species section "
          write(*,*) "in parameter file :: atype_name, charge, epsilon, sigma, atype_freeze "
          write(*,*) "please check input file "
          write(*,*) ""
          stop
       end if

       atype_name(i_type) = args(1)(1:MAX_ANAME)
       read(args(2),*) atype_chg(i_type)
       read(args(3),*) atype_lj_parameter(i_type,i_type,1)     
       read(args(4),*) atype_lj_parameter(i_type,i_type,2)
       read(args(5),*) atype_freeze(i_type)

       ! warn if freezing atom type
       if ( atype_freeze(i_type) == 1 ) then
          write(*,*) "NOTE: Atomtype ", atype_name(i_type)," will be frozen during the simulation"
       end if

       !      read(file_h,*) atype_name(i_type),atype_chg(i_type),atype_lj_parameter(i_type,i_type,1),atype_lj_parameter(i_type,i_type,2),atype_lj_parameter(i_type,i_type,3),atype_pol(i_type)

       ! move spaces for matching
       call trim_end( atype_name(i_type) )
    enddo


    ! set up lj_bkghm_index array, in this case all atom pairs use the same potential
    lj_bkghm_index = lj_bkghm


    ind=0
    ! find out whether to read cross terms or generate them
    do
       Read(file_h,'(A)',Iostat=inputstatus) line
       If(inputstatus < 0) Exit
       ind=INDEX(line,'cross_terms')
       IF(ind .NE. 0) Exit
    enddo


    if ( ind .ne. 0 ) then
       read(file_h,*) n_cross
       if(n_cross > 0) then
          do i_param=1, n_cross
             read(file_h,*) ind1,ind2, store(1),store(2),store(3)
             Select Case(lj_comb_rule)
             Case("opls")
                ! C12 goes first, this is read in as second parameter
                atype_lj_parameter(ind1,ind2,1)=store(2)
                atype_lj_parameter(ind1,ind2,2)=store(1)
                atype_lj_parameter(ind2,ind1,1)=store(2)
                atype_lj_parameter(ind2,ind1,2)=store(1)
             Case default
                atype_lj_parameter(ind1,ind2,:)=store(:)
                atype_lj_parameter(ind2,ind1,:)=store(:)
                ! make sure we have corret combination rule selected, sigma and epsilon should be << 1000
                if ( ( atype_lj_parameter(ind1,ind2,1) > 1000d0 ) .or. ( atype_lj_parameter(ind1,ind2,2) > 1000d0 ) ) then
                   write(*,*) "looks like combination rule should be opls.  Cross term parameters look "
                   write(*,*) "like C6 and C12 instead of epsilon and sigma based on their magnitudes  "
                   write(*,*) "please check consistency of cross terms and parameters "
                   stop
                end if

             end Select
             gen_cross_terms(ind1,ind2)=1
             gen_cross_terms(ind2,ind1)=1
          enddo
       endif
    endif

    ! finally make sure we don't have two of the same atom type defined

    do i_type=1, n_atom_type - 1
       do j_type = i_type + 1 , n_atom_type
          if ( atype_name(i_type) .eq. atype_name(j_type) ) then
             write(*,*) "atomic parameters defined more than once for atom type", atype_name(i_type)
             stop
          endif
       enddo
    enddo

    close(file_h)


  end subroutine read_param






  !******************************************************!
  ! this subroutine fills in parameter arrays in global variables
  ! to be used during the simulation for a non SAPT-based force field
  !
  ! for lennard jones force field, the final parameters will be
  ! C12 and C6.  No matter the force field or combination rules used, parameters
  ! will be read in as epsilon and sigma.  
  !    For standard (Lorentz-Berthelot) combination rules,
  ! cross parameters for epsilon and sigma will be generated first, and then these parameters
  ! will be turned into C12 and C6
  !    For Opls combination rules, C12, C6 atomic parameters will be generated first, and then
  ! cross parameters will be created from these parameters
  ! 
  !*****************************************************!
  subroutine gen_param(aname, total_atoms, charge, atom_type_index, gen_cross_terms)
    use global_variables
    character(*), dimension(:),intent(in) :: aname
    integer, intent(in)              :: total_atoms
    real*8, dimension(:),intent(out) :: charge
    integer, dimension(:), intent(in) :: atom_type
    integer,dimension(:,:),intent(in) :: gen_cross_terms

    integer:: i_param, j_param, i_atom,i_type,flag

    ! if opls force field, we need to create C12 and C6 atomic terms first
    Select Case( lj_bkghm )
    Case(2)
       Select Case( lj_comb_rule )
       Case("opls")
          do i_param=1,n_atom_type
             call gen_C12_C6_epsilon_sigma(atype_lj_parameter,i_param,i_param)
          enddo
       End Select
    End Select

    ! create cross term parameters

    do i_param=1, n_atom_type
       do j_param=1,n_atom_type
          if (gen_cross_terms(i_param,j_param) .eq. 0) then
             ! if these cross terms weren't given in input file
             call combination_rule_cross_terms(atype_lj_parameter,i_param,j_param,lj_bkghm,lj_comb_rule)
          else
             ! use given terms
             atype_lj_parameter(i_param,j_param,:) = atype_lj_parameter(i_param,j_param,:)
          endif
       enddo
    enddo

    ! if Lorentz-Berthelot combination rules were used, we now  need to create C12 and C6 atomic and cross terms
    Select Case( lj_bkghm )
    Case(2)
       Select Case( lj_comb_rule )
       Case("standard")
          do i_param=1,n_atom_type
             do j_param=1,n_atom_type
                call gen_C12_C6_epsilon_sigma(atype_lj_parameter,i_param,j_param)
             enddo
          enddo
       End Select
    End Select

    ! now create atom index array that links atoms to parameters
       do i_atom =1, total_atoms
          flag=0
          do  i_type=1, n_atom_type
             if (atype_name(i_type) .eq. aname(i_atom) ) then
                flag=1
                i_param = i_type
                exit 
             endif
          enddo
          ! make sure all atom types have parameters
          if (flag .eq. 0 ) then
             write(*,*) "atom type    ", aname(i_atom), "doesn't have force field parameters!"
             stop
          endif

          ! set index, chg, and polarizability
          atom_type_index(i_atom) = i_param
          charge(i_atom) = atype_chg(i_param)
       enddo

  end subroutine gen_param


  !***************************************************
  ! this subroutine generates cross terms for buckingham or lennard jones
  ! force fields using the desired combination rules
  !
  ! for Lennard Jones force field,
  ! Lorentz-Berthelot combination rules operate on sigma and epsilon,
  ! Opls combination rules operate on C6, C12
  !***************************************************
  subroutine combination_rule_cross_terms(atype_lj_parameter,i_param,j_param,lj_bkghm,lj_comb_rule,exp_only)
    real*8, dimension(:,:,:), intent(inout) :: atype_lj_parameter
    integer, intent(in) :: i_param,j_param
    integer, intent(in) :: lj_bkghm
    character(*), intent(in) :: lj_comb_rule
    integer, intent(in),optional :: exp_only

    real*8 :: e1, e2
    integer :: flag_exponent=0

    ! if buckingham force field, we may want separate subroutine to generate cross terms for coefficients if we are using an energy decomposed force field.  We may only use this routine to generate exponents
    Select Case(lj_bkghm)
    Case(1)
       if ( present(exp_only) ) then
          flag_exponent = 1
       endif
    End Select

    Select Case(lj_bkghm)
    Case(1)
       !************* buckingham force field
       Select Case(lj_comb_rule)
       Case("standard")
          ! all geometric combination rules
          atype_lj_parameter(i_param,j_param,1) = sqrt(atype_lj_parameter(i_param,i_param,1)*atype_lj_parameter(j_param,j_param,1))    ! A coefficient
          atype_lj_parameter(i_param,j_param,2) = sqrt(atype_lj_parameter(i_param,i_param,2)*atype_lj_parameter(j_param,j_param,2))    ! B parameter
          atype_lj_parameter(i_param,j_param,3) = sqrt(atype_lj_parameter(i_param,i_param,3)*atype_lj_parameter(j_param,j_param,3))    ! C parameter
       Case("ZIFFF")
          ! only do prefactors if exp_only input isn't present
          if ( flag_exponent .eq. 0 ) then
             atype_lj_parameter(i_param,j_param,1) = sqrt(atype_lj_parameter(i_param,i_param,1)*atype_lj_parameter(j_param,j_param,1))    ! A coefficient 
          endif
          ! exponent combination rules given by ZIF FF rules
          ! B exponent
          e1 = atype_lj_parameter(i_param,i_param,2)
          e2 = atype_lj_parameter(j_param,j_param,2)
          atype_lj_parameter(i_param,j_param,2) = (e1+e2) * (e1*e2)/(e1**2+e2**2) 
          atype_lj_parameter(i_param,j_param,3) = sqrt(atype_lj_parameter(i_param,i_param,3)*atype_lj_parameter(j_param,j_param,3))    ! C coefficient
       Case default
          stop "lj_comb_rule parameter isn't recognized for buckingham type force field.  Please use either 'standard' for geometric combination rules, or 'ZIFFF' for ZIF FF combination rules"
       End Select
    Case(2)
       !************* lennard jones force field
       Select Case(lj_comb_rule)
       Case("standard")
          ! Lorentz_Berthelot combination rules
          atype_lj_parameter(i_param,j_param,1) = sqrt(atype_lj_parameter(i_param,i_param,1)*atype_lj_parameter(j_param,j_param,1))    ! epsilon
          atype_lj_parameter(i_param,j_param,2) = (atype_lj_parameter(i_param,i_param,2) + atype_lj_parameter(j_param,j_param,2))/dble(2)   ! sigma
       Case("opls")
          ! opls combination rules
          atype_lj_parameter(i_param,j_param,1) = sqrt(atype_lj_parameter(i_param,i_param,1)*atype_lj_parameter(j_param,j_param,1))   ! C12        
          atype_lj_parameter(i_param,j_param,2) = sqrt(atype_lj_parameter(i_param,i_param,2)*atype_lj_parameter(j_param,j_param,2))   ! C6
       Case default
          stop "lj_comb_rule parameter isn't recognized for lennard jones force field.  Please use either 'standard' for Lorentz-Berthelot combination rules or 'opls' for opls combination rules."
       End Select
    End Select


  end subroutine combination_rule_cross_terms


  !**************************************************
  ! this subroutine creates C12 and C6 coefficients from epsilon and sigma
  ! parameters for a lennard jones force field
  !**************************************************
  subroutine gen_C12_C6_epsilon_sigma(atype_lj_parameter,i_param,j_param)
    real*8, dimension(:,:,:), intent(inout) :: atype_lj_parameter
    integer, intent(in)  :: i_param,j_param

    real*8 :: epsilon, sigma

    epsilon = atype_lj_parameter(i_param,j_param,1)
    sigma   = atype_lj_parameter(i_param,j_param,2)

    ! C12
    atype_lj_parameter(i_param,j_param,1) = 4d0*epsilon*sigma**12
    ! C6
    atype_lj_parameter(i_param,j_param,2) = 4d0*epsilon*sigma**6

  end subroutine gen_C12_C6_epsilon_sigma




  !**************************************************
  ! this subroutine reads in special 1-4 interaction parameters
  ! as used in the GROMOS force field
  !
  ! by default, the 1-4 interaction parameters are the
  ! same as the normal LJ cross terms
  !**************************************************
  subroutine read_generate_14_interaction_parameters( ifile_pmt )
    use global_variables
    character(*),intent(in)::ifile_pmt
    integer::file_h=15, n_pairs, i_pairs, index1, index2, inputstatus, ind
    character(100)::line
    character(MAX_ANAME) :: atomtype1, atomtype2
    real*8 :: C6, C12

    ! by default, use the standard lj parameters
    atype_lj_parameter_14 = atype_lj_parameter

    open(unit=file_h,file=ifile_pmt,status="old")

    do
       Read(file_h,'(A)',Iostat=inputstatus) line
       If(inputstatus < 0) Exit
       ind=INDEX(line,'pairtypes')
       IF(ind /= 0) Exit
    enddo

    ! if 1-4 pair section is present
    if ( ind /= 0 ) then
       read(file_h,*) n_pairs

       do i_pairs=1, n_pairs
          read(file_h,*) atomtype1, atomtype2, C6, C12
          call trim_end( atomtype1 )
          call trim_end( atomtype2 )
          call atype_name_reverse_lookup( atomtype1, index1 )
          call atype_name_reverse_lookup( atomtype2, index2 )

          write(*,*) "explicit C6, C12 parameters read in for 1-4 interaction between"
          write(*,*) "atoms ", atomtype1, " and ", atomtype2

          atype_lj_parameter_14( index1, index2 , 1 ) = C12
          atype_lj_parameter_14( index2, index1 , 1 ) = C12
          atype_lj_parameter_14( index1, index2 , 2 ) = C6
          atype_lj_parameter_14( index2, index1 , 2 ) = C6
       end do

    end if

    close(file_h)


  end subroutine read_generate_14_interaction_parameters



  !******************************************************!
  ! this subroutine generates molecule_type arrays for force field
  ! this subroutine sets global variables:: n_molecule_type, molecule_type
  !*****************************************************!
  subroutine gen_molecule_type_data( n_mole, molecule_data, atom_data )
    use global_variables
    integer, intent(in) :: n_mole
    type(molecule_data_type), dimension(:), intent(in) :: molecule_data
    type(atom_data_type) , dimension(:), intent(in)    :: atom_data

    !******** this is a local data structure with pointers that will be set
    ! to subarrays of atom_data arrays for the specific atoms in the molecule
    type(single_molecule_data_type) :: single_molecule_data
    integer:: i_mole,j_mole,i_atom, old_type,count, flag_name

    ! first molecule type is first molecule in input file
    i_mole = 1
    ! attach pointers to atom subarrays for 1st molecule
    call return_molecule_block( single_molecule_data , atom_data, molecule_data(i_mole)%n_atom, molecule_data(i_mole)%atom_index )
    n_molecule_type = 1
    do i_atom=1, molecule_data(i_mole)%n_atom
       molecule_type(n_molecule_type,i_atom) = single_molecule_data%atom_type_index(i_atom)
    enddo
    ! mark end
    if ( molecule_data(i_mole)%n_atom < MAX_N_ATOM_TYPE ) then
       molecule_type(n_molecule_type,molecule_data(i_mole)%n_atom+1) = MAX_N_ATOM_TYPE + 1
    endif

    molecule_type_name(n_molecule_type)=molecule_data(i_mole)%mname
    molecule_index(i_mole)= n_molecule_type

    ! loop over rest of molecules
    do i_mole =2, n_mole
    ! reattach pointers to atom subarrays for new molecule
    call return_molecule_block( single_molecule_data , atom_data, molecule_data(i_mole)%n_atom, molecule_data(i_mole)%atom_index )

       old_type =0
       flag_name=0
       ! loop over all stored molecule types to see if this solute molecule is a new type
       do j_mole =1, n_molecule_type
          ! see if same name
          if ( molecule_type_name(j_mole) .eq. molecule_data(i_mole)%mname ) then
             flag_name = 1
          end if

          ! now check to see if they have same number of atoms
          do i_atom =1, MAX_N_ATOM_TYPE
             if( molecule_type(j_mole,i_atom) .eq. MAX_N_ATOM_TYPE + 1) exit
             count = i_atom
          enddo
          if( molecule_data(i_mole)%n_atom .eq. count) then
             do i_atom=1, n_atom(i_mole)
                if ( single_molecule_data%atom_type_index(i_atom) .ne. molecule_type(j_mole,i_atom)) then
                   goto 100
                endif
             enddo
             ! make sure that we haven't just matched a smaller molecule to part of a larger molecule
             if ( molecule_data(i_mole)%n_atom < MAX_N_ATOM_TYPE ) then
                if ( molecule_type(j_mole,molecule_data(i_mole)%n_atom+1 ) .ne. ( MAX_N_ATOM_TYPE + 1 ) ) then
                   goto 100
                endif
             end if
             ! make sure these molecules have the same name to avoid confusion, note this is different on the name
             ! check at the beginning of the loop, because that checks to see if any previous molecule
             ! had the same name, not just this specific molecule
             if ( molecule_type_name(j_mole) .ne. molecule_data(i_mole)%mname ) then
                stop "two identical  molecules have different names!"
             endif

             ! here we have matched this molecule with an old molecule type.  Record the index
             molecule_index(i_mole) = j_mole
             old_type = 1
100          continue
          endif
       enddo
       ! if this is a new type of molecule, record it
       if (old_type .eq. 0) then
          ! make sure this molecule has a new name
          if ( flag_name == 1 ) then
             stop "can't use the same name for two different molecules!"
          endif

          n_molecule_type = n_molecule_type + 1

          ! make sure we don't have too many molecule types
          if ( n_molecule_type > MAX_N_MOLE_TYPE ) then
             write(*,*) ""
             write(*,*) "Too many different types of molecules in simulation"
             write(*,*) "Please increase setting of MAX_N_MOLE_TYPE"
             write(*,*) ""
             stop
          end if

          do i_atom=1,molecule_data(i_mole)%n_atom
             molecule_type(n_molecule_type,i_atom) = single_molecule_data%atom_type_index(i_atom)
          enddo
          ! mark end
          if ( molecule_data(i_mole)%n_atom < MAX_N_ATOM_TYPE ) then
             molecule_type(n_molecule_type,molecule_data(i_mole)%n_atom+1 ) = MAX_N_ATOM_TYPE + 1
          end if
          molecule_type_name(n_molecule_type) = molecule_data(i_mole)%mname
          molecule_index(i_mole)= n_molecule_type

       endif
    enddo


  end subroutine gen_molecule_type_data



  !***************************************
  ! this subroutine fills in mass data structure array
  ! using atom-type masses that are stored as global variables
  !***************************************
  subroutine fill_mass( total_atoms, mass, atom_type_index )
    use global_variables
    integer, intent(in) :: total_atoms
    real*8,dimension(:),intent(out) :: mass
    integer, intent(in) :: atom_type_index

    integer :: i_atom, i_type
    real*8  :: mass_atom

       do i_atom=1, total_atoms
          ! get atom type
          i_type = atom_type_index(i_atom)
          mass_atom = atype_mass(i_type)
          if ( mass_atom < 0d0 ) then
             write(*,*) "mass for atom type ", atype_name(i_type)
             write(*,*) "has not been read in from topology file"
             stop
          end if
          mass(i_atom) = mass_atom
       enddo

  end subroutine fill_mass




end module initialize_routines
