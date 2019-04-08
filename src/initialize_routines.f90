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
  subroutine initialize_simulation(system_data, molecule_data, atom_data, file_io_data, verlet_list_data, PME_data )
    use global_variables
    use pme_routines
    use bonded_interactions
    Type(system_data_type),intent(inout)                :: system_data
    Type(molecule_data_type),dimension(:),allocatable, intent(inout) :: molecule_data
    Type(atom_data_type),intent(inout)                  :: atom_data
    Type(file_io_data_type),intent(inout)               :: file_io_data
    Type(verlet_list_data_type),intent(inout)           :: verlet_list_data
    Type(PME_data_type), intent(inout)                  :: PME_data

    integer, dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE) :: gen_cross_terms
    integer:: total_atoms, size, flag_junk


    !********** this is for generating random velocities for Maxwell-Boltzmann
    call initialize_random_seed

    ! **************************** read atom coordinates   ***********************************!

    ! first get number of molecules, allocate datastructures
    call read_gro_number_molecules( file_io_data%ifile_gro_file_h, file_io_data%ifile_gro, system_data%n_mole )

    ! allocate molecule_data array
    allocate( molecule_data(system_data%n_mole) )

    ! this call will fill in data structures atom_data and molecule_data

    ! if trajectory restart, scan to end of .gro trajectory file...
    ! this will allocate atom_data%aname, atom_data%xyz arrays
    Select Case( restart_trajectory )
    Case("yes")
       open( file_io_data%ofile_traj_file_h, file=file_io_data%ofile_traj, status='old' )
       call scan_grofile_restart( file_io_data%ofile_traj_file_h, n_old_trajectory )
       call read_gro( file_io_data%ofile_traj_file_h, system_data, molecule_data, atom_data )
    Case default
    open( file_io_data%ifile_gro_file_h, file=file_io_data%ifile_gro, status='old' )
    call read_gro( file_io_data%ifile_gro_file_h, system_data, molecule_data, atom_data )
    close( file_io_data%ifile_gro_file_h )
    End Select
    ! ******************  Allocate remainder of atom_data arrays, initialize
    total_atoms = system_data%total_atoms
    allocate( atom_data%velocity(3,total_atoms) )
    allocate( atom_data%force(3,total_atoms) )
    allocate( atom_data%mass(total_atoms) )
    allocate( atom_data%charge(total_atoms) )
    allocate( atom_data%atom_type_index(total_atoms) )

    atom_data%velocity=0d0;atom_data%force=0d0;atom_data%mass=0d0;atom_data%charge=0d0; atom_data%atom_type_index=0


    ! if restarting, read in velocities from restart file
    Select Case( restart_trajectory )
    Case("yes")
       call read_velocity_restart_checkpoint( file_io_data%ifile_velocity_file_h, file_io_data%ifile_velocity, atom_data%velocity, n_old_trajectory )
    case default
    ! open velocity_checkpoint file for printing if we are checkpointing velocities
       Select Case(checkpoint_velocity)
       Case("yes")
            open( file_io_data%ifile_velocity_file_h, file=file_io_data%ifile_velocity, status='new' )
       End Select
    End Select

    ! system volume
    system_data%volume = volume( system_data%box )

    ! initialize the transformation matrix to box coordinates
    call initialize_non_orth_transform ( system_data%box, system_data%xyz_to_box_transform )

    !***************** make sure molecules are not broken up by a box translation, as this can happen in GROMACs, and could be using GROMACS output as input here
    call fix_intra_molecular_shifts( system_data%n_mole , molecule_data , atom_data , system_data%box,  system_data%xyz_to_box_transform )

    !***************** make sure cutoff distances are fine for box size, and
    !that box type is supported
    call check_cutoffs_box( real_space_cutoff, verlet_list_data%verlet_cutoff, system_data%box )


    !************************************** get parameters****************************************************!
    ! get parameters for all atom types, assuming there are no more than MAX_N_ATOM_TYPE types of atoms total
    !*****************************************************************************************************!
    call read_param( file_io_data%ifile_ffpmt_file_h, file_io_data%ifile_ffpmt,gen_cross_terms)
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

    ! read topology file
    call read_topology_file( file_io_data%ifile_top_file_h, file_io_data%ifile_top )

    ! mass is read in from topology file, fill in mass for each molecule
    call fill_mass( system_data%total_atoms, atom_data%mass, atom_data%atom_type_index )


    ! center of mass
    call update_r_com( system_data%n_mole, molecule_data, atom_data )

    ! center of mass of molecule might be outside of box after calling subroutine fix_intra_molecular_shifts, fix this
    call shift_molecules_into_box( system_data%n_mole , molecule_data , atom_data , system_data%box,  system_data%xyz_to_box_transform )


    ! generate exclusions, note that some exclusions may have already been explicity read in from topology file.
    call generate_intramolecular_exclusions


    !*********************************initialize verlet list
    call allocate_verlet_list( verlet_list_data, system_data%total_atoms, system_data%volume )
    call construct_verlet_list( verlet_list_data, atom_data, molecule_data, system_data%total_atoms, system_data%box, system_data%xyz_to_box_transform  )
    ! the "1" input to update_verlet_displacements signals to initialize the displacement array
    call update_verlet_displacements( system_data%total_atoms, atom_data%xyz, verlet_list_data , system_data%box, system_data%xyz_to_box_transform , flag_junk, 1 )


    ! initialize ms_evb simulation
    Select Case(ms_evb_simulation)
    Case("yes")
       call initialize_evb_simulation( system_data%n_mole, molecule_data, file_io_data )
       ! we are storing dQ_dr grid for ms-evb, allocate data arrays for this

       ! grid size:  each atom has spline_order^3 non-zero derivatives in the Q grid
       size = system_data%total_atoms * PME_data%spline_order**3
       ! size could be big depending on the system, make sure we don't ask for too much memory
       if ( size > 10000000 ) then
         write(*,*) "are you sure you want to store dQ_dr?  This requires allocating an array"
         write(*,*) " bigger than 3 x ", size
          stop
       end if
       ! allocate column major
       size=PME_data%spline_order**3
       allocate( PME_data%dQ_dr(3,size,system_data%total_atoms) , PME_data%dQ_dr_index(3,size,system_data%total_atoms) )

    End Select

    ! for storing PME reciprocal space forces
    allocate(PME_data%force_recip(3,system_data%total_atoms) )

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
  subroutine initialize_energy_force(system_data, molecule_data, atom_data, verlet_list_data, PME_data, file_io_data, integrator_data )
    use global_variables
    use pme_routines
    use MKL_DFTI
    use total_energy_forces
    Type(system_data_type),intent(inout)                :: system_data
    Type(molecule_data_type),dimension(:),intent(inout) :: molecule_data
    Type(atom_data_type),intent(inout)                  :: atom_data
    Type(verlet_list_data_type),intent(inout)           :: verlet_list_data
    Type(PME_data_type), intent(inout)                  :: PME_data
    Type(file_io_data_type), intent(in)                 :: file_io_data
    Type(integrator_data_type), intent(in)              :: integrator_data

    integer:: i,length(3),status, pme_grid
    real*8 :: a(3), b(3), c(3), ka(3), kb(3), kc(3),kk(3,3)
    real*8 :: x


    !************************************* initialize ewald/pme ************************************************!

    ! note here that tot_chg was passed to this subroutine, so drude oscillator charges are accounted for
    ! if they are present.  Use n_atom_drude for Ewald_self, so that drude oscillators are included
    call update_Ewald_self( system_data%total_atoms, atom_data, PME_data )

    pme_grid=PME_data%pme_grid

    ! allocate arrays
    allocate(PME_data%CB(pme_grid,pme_grid,pme_grid),PME_data%Q_grid(pme_grid,pme_grid,pme_grid),PME_data%theta_conv_Q(pme_grid,pme_grid,pme_grid))

    ! set up fourier transform descriptors
    length(:)= pme_grid
    status=DftiCreateDescriptor(PME_data%dfti_desc, DFTI_DOUBLE, DFTI_COMPLEX, 3, length)
    status=DftiCommitDescriptor(PME_data%dfti_desc)
    ! don't scale back transform because pick up factor of K^3 from convolution
    status=DftiCreateDescriptor(PME_data%dfti_desc_inv, DFTI_DOUBLE, DFTI_COMPLEX, 3, length)
    status=DftiCommitDescriptor(PME_data%dfti_desc_inv)

    ! initialize PME dependency on system volume
    call periodic_box_change(system_data, PME_data, verlet_list_data, atom_data, molecule_data)

    ! grid B_splines
       if (PME_data%spline_order .eq. 6) then
          allocate(PME_data%B6_spline(PME_data%spline_grid),PME_data%B5_spline(PME_data%spline_grid))
          do i=1,PME_data%spline_grid
             PME_data%B6_spline(i)=B_spline(6./dble(PME_data%spline_grid)*dble(i),6)
             PME_data%B5_spline(i)=B_spline(5./dble(PME_data%spline_grid)*dble(i),5)
          enddo
       else if (PME_data%spline_order .eq. 4) then
          allocate(PME_data%B4_spline(PME_data%spline_grid),PME_data%B3_spline(PME_data%spline_grid))
          do i=1,PME_data%spline_grid
             PME_data%B4_spline(i)=B_spline(4./dble(PME_data%spline_grid)*dble(i),4)
             PME_data%B3_spline(i)=B_spline(3./dble(PME_data%spline_grid)*dble(i),3)
          enddo
       else
          stop "requested spline order not implemented"
       endif

    !  allocate and grid complementary error function
    allocate(PME_data%erfc_table(PME_data%erfc_grid))
    ! make sure erfc_max is set big enough to cover distances up to the cutoff
    if ( PME_data%erfc_max < PME_data%alpha_sqrt * real_space_cutoff ) then
       write(*,*) "please increase setting of erfc_max.  Must have erfc_max > alpha_sqrt * cutoff "
       stop
    end if

    do i=1,PME_data%erfc_grid
       PME_data%erfc_table(i)= erfc(PME_data%erfc_max/dble(PME_data%erfc_grid)*dble(i))
    enddo

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
    End Select


    !**************************** Get initial forces and energy *************************************!
    Select Case(ms_evb_simulation)
    Case("yes")
       call ms_evb_calculate_total_force_energy( system_data, molecule_data, atom_data, verlet_list_data, PME_data, file_io_data,integrator_data%n_output )
    Case("no")
       call calculate_total_force_energy( system_data, molecule_data, atom_data, verlet_list_data, PME_data )
    End Select
    !************************************************************************************************!

  end subroutine initialize_energy_force




  !**************************************************!
  ! this subroutine reads atom_type parameters from input file for a non SAPT-based force field
  !*************************************************
  subroutine read_param( file_h, ifile_pmt, gen_cross_terms)
    use global_variables
    integer, intent(in)    :: file_h
    character(*),intent(in):: ifile_pmt
    integer,dimension(:,:),intent(out) :: gen_cross_terms

    integer::i_type, j_type, i_param, n_cross, inputstatus,ind1,ind2,ind3
    integer,parameter :: max_param=20
    real*8,parameter :: small=1D-6, exp_init=3d0
    character(20),dimension(max_param) :: args
    real*8,dimension(6) :: store
    integer        :: nargs
    character(300) :: input_string
    character(50)::line
    character(5)::c_name


    gen_cross_terms=0

    ! if we're not using SAPT-FF interactions, initialize to zero.  But
    ! initialize exponents to finite value to be safe for combination rules...
    atype_vdw_tmp=0d0
    atype_vdw_tmp(:,:,5) = exp_init

    ! store has dimension 6, since atype_vdw_parameter has dimension 6 ( changed to accomodate C12).  Zero these components that will not be used
    store=0d0

    open(unit=file_h,file=ifile_pmt,status="old")

    ! READ FILE, look for three headings
    ! HEADING 1: 'solute_species'  -- standard coulomb and LJ interactions
    ! HEADING 2: 'custom_sapt_parameters' -- SAPT-FF force field section
    ! HEADING 3: 'cross_terms'  -- Explicit cross terms for LJ interactions
    ! HEADING 4: 'sapt_exclusions' -- We want water/hydronium to have SAPT
    ! interactions with all other species, except between themselves
    ! note that there could also be another section called 'pairtypes' which is
    ! for custom 1-4 interactions, but this is read by a different subroutine
    do
       Read(file_h,'(A)',Iostat=inputstatus) line
       If(inputstatus < 0) Exit
       ind1=INDEX(line,'solute_species')
       ind2=INDEX(line,'custom_sapt_parameters')
       ind3=INDEX(line,'cross_terms')

       ! ******** solute_species section ************
       IF(ind1 .NE. 0) then

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
              read(args(3),*) atype_vdw_parameter(i_type,i_type,1)     
              read(args(4),*) atype_vdw_parameter(i_type,i_type,2)
              read(args(5),*) atype_freeze(i_type)

              ! warn if freezing atom type
              if ( atype_freeze(i_type) == 1 ) then
                 write(*,*) "NOTE: Atomtype ", atype_name(i_type)," will be frozen during the simulation"
              end if

              ! move spaces for matching
              call trim_end( atype_name(i_type) )
          enddo

       !******* SAPT-FF section
       ELSE IF(ind2 .NE. 0) then
          Read(file_h,'(A)',Iostat=inputstatus) line
          do i_type=1, n_atom_type
             read(file_h,'(A)') input_string
             call parse(input_string," ",args,nargs)
             if ( nargs /= 10 ) then
                write(*,*) ""
                write(*,*) "should have 10 input arguments under custom_sapt_parameters"
                write(*,*) "in parameter file :: atype_name, 4 A parameters,1 B "
                write(*,*) "and 4 C parameters "
                write(*,*) "please check input file"
                stop
             end if
             c_name = args(1)(1:MAX_ANAME)
             read(args(2),*) atype_vdw_tmp(i_type,i_type,1)
             read(args(3),*) atype_vdw_tmp(i_type,i_type,2)
             read(args(4),*) atype_vdw_tmp(i_type,i_type,3)
             read(args(5),*) atype_vdw_tmp(i_type,i_type,4)
             read(args(6),*) atype_vdw_tmp(i_type,i_type,5)
             read(args(7),*) atype_vdw_tmp(i_type,i_type,6)
             read(args(8),*) atype_vdw_tmp(i_type,i_type,7)
             read(args(9),*) atype_vdw_tmp(i_type,i_type,8)
             read(args(10),*) atype_vdw_tmp(i_type,i_type,9)
          enddo

       !******* Explicit Cross terms 
       ELSE IF ( ind3 .ne. 0 ) then
          read(file_h,*) n_cross
          if(n_cross > 0) then
             do i_param=1, n_cross
                read(file_h,*) i_type,j_type, store(1),store(2),store(3)
                Select Case(lj_comb_rule)
                Case("opls")
                   ! C12 goes first, this is read in as second parameter
                   atype_vdw_parameter(i_type,j_type,1)=store(2)
                   atype_vdw_parameter(i_type,j_type,2)=store(1)
                   atype_vdw_parameter(j_type,i_type,1)=store(2)
                   atype_vdw_parameter(j_type,i_type,2)=store(1)
                Case default
                   atype_vdw_parameter(i_type,j_type,:)=store(:)
                   atype_vdw_parameter(j_type,i_type,:)=store(:)
                   ! make sure we have corret combination rule selected, sigma and epsilon should be << 1000
                   if ( ( atype_vdw_parameter(i_type,j_type,1) > 1000d0 ) .or. ( atype_vdw_parameter(i_type,j_type,2) > 1000d0 ) ) then
                      write(*,*) "looks like combination rule should be opls.  Cross term parameters look "
                      write(*,*) "like C6 and C12 instead of epsilon and sigma based on their magnitudes  "
                      write(*,*) "please check consistency of cross terms and parameters "
                      stop
                   end if
                end Select
                gen_cross_terms(i_type,j_type)=1
                gen_cross_terms(j_type,i_type)=1
             enddo
          endif

      END IF
   
    enddo

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
    real*8, dimension(:),intent(inout) :: charge
    integer, dimension(:), intent(inout) :: atom_type_index
    integer,dimension(:,:),intent(in) :: gen_cross_terms

    integer:: i_param, j_param, i_atom,i_type,flag
    real*8, parameter :: small=1D-6
    
    ! if opls force field, we need to create C12 and C6 atomic terms first
    Select Case( lj_comb_rule )
       Case("opls")
          do i_param=1,n_atom_type
             call gen_C12_C6_epsilon_sigma(atype_vdw_parameter,i_param,i_param)
          enddo
    End Select
    ! create cross terms first
    do i_param=1, n_atom_type
       do j_param=1, n_atom_type
        if (i_param .NE. j_param) then
          if (gen_cross_terms(i_param,j_param) .eq. 0) then

            if (atype_vdw_parameter(i_param,i_param,1) > small .and. atype_vdw_parameter(j_param,j_param,1) > small) then 
              ! this is Lennard-Jones interaction
              atype_vdw_type(i_param,j_param)=0
              call combination_rule_cross_terms(atype_vdw_parameter,i_param,j_param,lj_comb_rule, atype_vdw_type)
            else if ( atype_vdw_tmp(i_param,i_param,5) > small .and. atype_vdw_tmp(j_param,j_param,5) > small ) then
              ! this is SAPT-FF interaction
              atype_vdw_type(i_param,j_param)=1
              call combination_rule_cross_terms(atype_vdw_parameter,i_param,j_param,lj_comb_rule, atype_vdw_type, atype_vdw_tmp)
            else
              ! no interaction, set atype_vdw_type to -1
              atype_vdw_type(i_param,j_param)=-1
              ! call cross-term combination rule to zero parameters, even though
              ! we shouldn't be using ...
              call combination_rule_cross_terms(atype_vdw_parameter,i_param,j_param,lj_comb_rule, atype_vdw_type)
            end if

          else
             ! use given terms, this is an LJ interaction
             atype_vdw_type(i_param,j_param)=0
          endif
        endif
       enddo
    enddo

    ! now diagonal terms
    do i_param=1, n_atom_type
      if (gen_cross_terms(i_param,i_param) .eq. 0) then
        if (atype_vdw_parameter(i_param,i_param,1) > small )  then
           ! LJ interaction    
           atype_vdw_type(i_param,i_param) = 0
           call combination_rule_cross_terms(atype_vdw_parameter,i_param,i_param,lj_comb_rule,atype_vdw_type)
        else if (atype_vdw_tmp(i_param,i_param,1) > small )  then
           ! SAPT-FF interaction
           atype_vdw_type(i_param,i_param) = 1
           call combination_rule_cross_terms(atype_vdw_parameter,i_param,i_param,lj_comb_rule,atype_vdw_type,atype_vdw_tmp)             
        else
           ! no interaction, set atype_vdw_type to -1
           atype_vdw_type(i_param,i_param)=-1
           ! call cross-term combination rule to zero parameters, even though
           ! we shouldn't be using ...
           call combination_rule_cross_terms(atype_vdw_parameter,i_param,i_param,lj_comb_rule, atype_vdw_type)
        end if
      else
          atype_vdw_type(i_param,i_param)=0
     endif
    enddo
    
    
    ! if Lorentz-Berthelot combination rules were used, we now  need to create C12 and C6 atomic and cross terms
    do i_param=1,n_atom_type
      do j_param=1,n_atom_type
        Select Case( atype_vdw_type(i_param,j_param) )
         Case(0)
          Select Case( lj_comb_rule )
           Case("standard")
                call gen_C12_C6_epsilon_sigma(atype_vdw_parameter,i_param,j_param)
          End Select
        End Select
      enddo
    enddo
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

       if (maxval(atype_vdw_tmp(:,:,5)) * real_space_cutoff > Tang_Toennies_max) then
          write(*,*) maxval(atype_vdw_tmp(:,:,5)), real_space_cutoff, Tang_Toennies_max
          write(*,*) "Table size for damping function is smaller than B_max*the real space cutoff"
          stop
       endif
  end subroutine gen_param


  !***************************************************
  ! this subroutine generates cross terms for buckingham or lennard jones
  ! force fields using the desired combination rules
  !
  ! for Lennard Jones force field,
  ! Lorentz-Berthelot combination rules operate on sigma and epsilon,
  ! Opls combination rules operate on C6, C12
  !***************************************************
  subroutine combination_rule_cross_terms(atype_vdw_parameter,i_param,j_param,lj_comb_rule, atype_vdw_type, atype_vdw_tmp)
    real*8, dimension(:,:,:), intent(inout) :: atype_vdw_parameter
    integer, dimension(:,:), intent(in) :: atype_vdw_type
    real*8, dimension(:,:,:), intent(in), optional :: atype_vdw_tmp
    integer, intent(in) :: i_param,j_param
    character(*), intent(in) :: lj_comb_rule
    real*8 :: a_ex, a_el, a_ind, a_dhf, b_tmp

    Select Case(atype_vdw_type(i_param, j_param))
    Case(1)
       !First handle all A terms and then combine them in the permanent
       !atype_vdw_parameter structure
       a_ex = sqrt(atype_vdw_tmp(i_param,i_param,1)*atype_vdw_tmp(j_param,j_param,1))
       a_el = sqrt(atype_vdw_tmp(i_param,i_param,2)*atype_vdw_tmp(j_param,j_param,2))
       a_ind = sqrt(atype_vdw_tmp(i_param,i_param,3)*atype_vdw_tmp(j_param,j_param,3))
       a_dhf = sqrt(atype_vdw_tmp(i_param,i_param,4)*atype_vdw_tmp(j_param,j_param,4))
       atype_vdw_parameter(i_param,j_param,1) = a_ex - a_el - a_ind - a_dhf
       !B term 
       b_tmp = (atype_vdw_tmp(i_param,i_param,5) + atype_vdw_tmp(j_param,j_param,5))*&
atype_vdw_tmp(i_param,i_param,5)*atype_vdw_tmp(j_param,j_param,5)/(atype_vdw_tmp(i_param,i_param,5)**2+&
atype_vdw_tmp(j_param,j_param,5)**2)
       atype_vdw_parameter(i_param,j_param,2)=b_tmp
       !C terms
       atype_vdw_parameter(i_param,j_param,3) = sqrt(atype_vdw_tmp(i_param,i_param,6)*atype_vdw_tmp(j_param,j_param,6))
       atype_vdw_parameter(i_param,j_param,4) = sqrt(atype_vdw_tmp(i_param,i_param,7)*atype_vdw_tmp(j_param,j_param,7))
       atype_vdw_parameter(i_param,j_param,5) = sqrt(atype_vdw_tmp(i_param,i_param,8)*atype_vdw_tmp(j_param,j_param,8))
       atype_vdw_parameter(i_param,j_param,6) = sqrt(atype_vdw_tmp(i_param,i_param,9)*atype_vdw_tmp(j_param,j_param,9))
    Case default
       !************* lennard jones force field
       Select Case(lj_comb_rule)
       Case("standard")
          ! Lorentz_Berthelot combination rules
          atype_vdw_parameter(i_param,j_param,1) = sqrt(atype_vdw_parameter(i_param,i_param,1)*atype_vdw_parameter(j_param,j_param,1))    ! epsilon
          atype_vdw_parameter(i_param,j_param,2) = (atype_vdw_parameter(i_param,i_param,2) + atype_vdw_parameter(j_param,j_param,2))/dble(2)   ! sigma
       Case("opls")
          ! opls combination rules
          atype_vdw_parameter(i_param,j_param,1) = sqrt(atype_vdw_parameter(i_param,i_param,1)*atype_vdw_parameter(j_param,j_param,1))   ! C12        
          atype_vdw_parameter(i_param,j_param,2) = sqrt(atype_vdw_parameter(i_param,i_param,2)*atype_vdw_parameter(j_param,j_param,2))   ! C6
       Case default
          stop "lj_comb_rule parameter isn't recognized for lennard jones force field.  Please use either 'standard' for Lorentz-Berthelot combination rules or 'opls' for opls combination rules."
       End Select
    End Select


  end subroutine combination_rule_cross_terms


  !**************************************************
  ! this subroutine creates C12 and C6 coefficients from epsilon and sigma
  ! parameters for a lennard jones force field
  !**************************************************
  subroutine gen_C12_C6_epsilon_sigma(atype_vdw_parameter,i_param,j_param)
    real*8, dimension(:,:,:), intent(inout) :: atype_vdw_parameter
    integer, intent(in)  :: i_param,j_param

    real*8 :: epsilon, sigma

    epsilon = atype_vdw_parameter(i_param,j_param,1)
    sigma   = atype_vdw_parameter(i_param,j_param,2)

    ! C12
    atype_vdw_parameter(i_param,j_param,1) = 4d0*epsilon*sigma**12
    ! C6
    atype_vdw_parameter(i_param,j_param,2) = 4d0*epsilon*sigma**6

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
    atype_vdw_parameter_14 = atype_vdw_parameter

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

          atype_vdw_parameter_14( index1, index2 , 1 ) = C12
          atype_vdw_parameter_14( index2, index1 , 1 ) = C12
          atype_vdw_parameter_14( index1, index2 , 2 ) = C6
          atype_vdw_parameter_14( index2, index1 , 2 ) = C6
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
    type(molecule_data_type), dimension(:), intent(inout) :: molecule_data
    type(atom_data_type),  intent(in)    :: atom_data

    !******** this is a local data structure with pointers that will be set
    ! to subarrays of atom_data arrays for the specific atoms in the molecule
    type(single_molecule_data_type) :: single_molecule_data
    integer:: i_mole,j_mole,i_atom, old_type, flag_name

    ! first molecule type is first molecule in input file
    i_mole = 1
    ! attach pointers to atom subarrays for 1st molecule
    call return_molecule_block( single_molecule_data , molecule_data(i_mole)%n_atom, molecule_data(i_mole)%atom_index, atom_type_index=atom_data%atom_type_index )
    n_molecule_type = 1

    ! allocate data structure for first molecule type
    molecule_type_data(n_molecule_type)%n_atom =  molecule_data(i_mole)%n_atom
    allocate( molecule_type_data(n_molecule_type)%atom_type_index(molecule_data(i_mole)%n_atom) )
    allocate( molecule_type_data(n_molecule_type)%pair_exclusions( molecule_data(i_mole)%n_atom , molecule_data(i_mole)%n_atom) )
    ! zero molecule_exclusions here, as it is possible that exclusions are being read in from topology file
     molecule_type_data(n_molecule_type)%pair_exclusions = 0
    ! for MS-EVB
    allocate( molecule_type_data(n_molecule_type)%evb_reactive_basic_atoms(molecule_data(i_mole)%n_atom) ) 
    molecule_type_data(n_molecule_type)%evb_reactive_basic_atoms=0

    do i_atom=1, molecule_data(i_mole)%n_atom
       molecule_type_data(n_molecule_type)%atom_type_index(i_atom) = single_molecule_data%atom_type_index(i_atom)
    enddo
    molecule_type_data(n_molecule_type)%mname = molecule_data(i_mole)%mname
    molecule_data(i_mole)%molecule_type_index = n_molecule_type

    ! loop over rest of molecules
    do i_mole =2, n_mole
    ! reattach pointers to atom subarrays for new molecule
    call return_molecule_block( single_molecule_data , molecule_data(i_mole)%n_atom, molecule_data(i_mole)%atom_index, atom_type_index=atom_data%atom_type_index )

       old_type =0
       flag_name=0
       ! loop over all stored molecule types to see if this solute molecule is a new type
       do j_mole =1, n_molecule_type
          ! see if same name
          if ( molecule_type_data(j_mole)%mname .eq. molecule_data(i_mole)%mname ) then
             flag_name = 1
          end if

          ! now check to see if they have same number of atoms
          if( molecule_data(i_mole)%n_atom .eq. molecule_type_data(j_mole)%n_atom ) then

             do i_atom=1, molecule_data(i_mole)%n_atom
                if ( single_molecule_data%atom_type_index(i_atom) .ne. molecule_type_data(j_mole)%atom_type_index(i_atom) ) then
                   goto 100
                endif
             enddo

             ! make sure these molecules have the same name to avoid confusion, note this is different on the name
             ! check at the beginning of the loop, because that checks to see if any previous molecule
             ! had the same name, not just this specific molecule
             if ( molecule_type_data(j_mole)%mname .ne. molecule_data(i_mole)%mname ) then
                stop "two identical  molecules have different names!"
             endif

             ! here we have matched this molecule with an old molecule type.  Record the index
             molecule_data(i_mole)%molecule_type_index = j_mole
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

          ! allocate data structure for this molecule type
          molecule_type_data(n_molecule_type)%n_atom =  molecule_data(i_mole)%n_atom
          allocate( molecule_type_data(n_molecule_type)%atom_type_index(molecule_data(i_mole)%n_atom) )
          allocate( molecule_type_data(n_molecule_type)%pair_exclusions( molecule_data(i_mole)%n_atom , molecule_data(i_mole)%n_atom) )
          ! zero molecule_exclusions here, as it is possible that exclusions are being read in from topology file
          molecule_type_data(n_molecule_type)%pair_exclusions = 0
          ! for MS-EVB
          allocate( molecule_type_data(n_molecule_type)%evb_reactive_basic_atoms(molecule_data(i_mole)%n_atom) )
          molecule_type_data(n_molecule_type)%evb_reactive_basic_atoms=0


          do i_atom=1, molecule_data(i_mole)%n_atom
             molecule_type_data(n_molecule_type)%atom_type_index(i_atom) = single_molecule_data%atom_type_index(i_atom)
          enddo
          molecule_type_data(n_molecule_type)%mname = molecule_data(i_mole)%mname
          molecule_data(i_mole)%molecule_type_index = n_molecule_type

       endif
       ! dissociate pointers
       call dissociate_single_molecule_data(single_molecule_data)
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
    integer, dimension(:),intent(in) :: atom_type_index

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
