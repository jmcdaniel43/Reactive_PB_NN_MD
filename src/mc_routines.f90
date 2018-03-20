!***********************************************************************
!
!  This module contains various subroutines used in the Monte Carlo simulation,
!  including data initialization subrouines, as well as subroutines that
!  (for a Lennard-Jones type force field) read in the input files
!
!***********************************************************************
module mc_routines
  use routines
  use ms_evb
  implicit none

contains

  !*************************************************************************
  ! this routine reads the input files for the solute, and if present, the framework
  ! and initializes the simulation
  !*************************************************************************
  subroutine initialize_simulation(box,n_mole,tot_n_mole,n_atom,n_atom_drude,xyz,mass,vel_atom, r_com,chg,chg_drude,tot_chg,ifile_conf, ifile_top, ifile_fconf,ifile_pmt, conf_file, traj_file, vel_file, ofile_traj)
    use global_variables
    use electrostatic
    use pairwise_interaction
    use pme_routines
    use bonded_interactions
    character(MAX_FN),intent(in) :: ifile_conf,ifile_top, ifile_fconf,ifile_pmt, ofile_traj
    integer, intent(in) :: conf_file, traj_file, vel_file
    integer,intent(out)::n_mole,tot_n_mole
    integer,dimension(:),intent(out)::n_atom,n_atom_drude
    real*8,dimension(:,:,:),intent(out)::xyz, vel_atom
    real*8,dimension(:,:),intent(out)::mass,chg,chg_drude,tot_chg
    real*8,dimension(:,:),intent(out)::box,r_com
    integer, dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE) :: gen_cross_terms
    real*8,dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE,9) :: temp_lj_parameter
    real*8,dimension(MAX_N_MOLE,MAX_N_ATOM) :: temp_chg
    integer,dimension(MAX_N_MOLE,MAX_N_ATOM) :: atom_molecule_map
    integer :: total_atoms_list
    ! these are temporary atom name and molecule name data structures,
    ! for general use in the program, global arrays atype_name, molecule_type_name will be used instead
    character(MAX_MNAME), dimension(MAX_N_MOLE) :: molecule_name
    character(MAX_ANAME), dimension(MAX_N_MOLE,MAX_N_ATOM) :: alist  

    integer:: i,i_mole,i_atom,n_type_s,index,is_overlap,size, flag_junk, flag_charge
    real*8 :: unit_cell(3,3),sum_chg
    real*8,dimension(3) :: inertia
    real*8,parameter :: framework_chg_thresh = 1D-5 ,solute_chg_thresh = 1D-6 , small_stop=1d0, small_warn=1d-1
    integer :: warn_inertia=0

    ! zero this, in case we aren't reading velocities from a velocity checkpoint file
    vel_atom=0d0

    ! initialize random seed, this will be overwritten if seed is read in in read_conf routine
    call initialize_random_seed

    !************************************** get parameters****************************************************!
    ! get parameters for all atom types, assuming there are no more than MAX_N_ATOM_TYPE types of atoms total
    ! remember no two atom names can conflict, not even those of framework and solute
    !*****************************************************************************************************!

    call read_param( ifile_pmt,gen_cross_terms)

    ! make sure we don't have too many atom types
    if ( n_atom_type > MAX_N_ATOM_TYPE ) then
       stop "number of atom types g.t. MAX_N_ATOM_TYPE.  Increase this value"
    endif


    ! **************************** coordinates and box size************************************************!
    ! zero n_atom array values that wont be used, as maxval(n_atom) will used to allocate an array later on
    n_atom=0
    n_atom_drude=0

    Select Case( restart_trajectory )
    Case("yes")
       ! if this is a trajectory restart, get previous coordinates and velocities
       call read_coordinates_restart_trajectory( traj_file, ofile_traj, n_mole, n_atom, xyz, molecule_name, alist, box, n_old_trajectory )
       call read_velocity_restart_checkpoint( vel_file , velocity_file, n_mole, n_atom, vel_atom, n_old_trajectory )
    Case default
       ! open file here, as files won't be open in read_conf
       open( conf_file, file=ifile_conf, status='old' )    
       call read_conf( conf_file, n_mole, n_atom, xyz, molecule_name,alist, box )
       close( conf_file )
    End Select

    ! if framework simulation,get framework coordinates
    if (framework_simulation .eq. 1) then
       call  read_fconf( ifile_fconf, n_mole, tot_n_mole, n_atom, xyz, temp_chg, flag_charge, alist, unit_cell )
       ! if flag_charge is set to 1, charges have been read in from this file, and overwrite
       ! those read in from the .pmt file.
       box = unit_cell
    else
       flag_charge=0
       tot_n_mole = n_mole
    endif

    ! initialize the transformation matrix to box coordinates
    call initialize_non_orth_transform ( box )

    !***************** make sure molecules are not broken up by a box translation, as this can happen
    ! in GROMACs, and could be using GROMACS output as input here
    do i_mole = 1, n_mole
       call check_intra_molecular_shifts(i_mole,n_atom,xyz,box)
    enddo

    !****************** check to see that molecules are not too close together
    is_overlap= check_move( xyz, tot_n_mole , n_mole, n_atom, box, r_com )


    !***************** make sure cutoff distances are fine for box size, and that box type is supported
    call check_cutoffs_box( box )

    ! ********************* generate cross terms and fill in parameter arrays****************************!
    ! fill in total parameter arrays that are stored in global_variables
    call gen_param(alist,tot_n_mole,n_atom, chg, chg_drude,gen_cross_terms)

    ! now get special 1-4 interactions from parameter file if needed
    call read_generate_14_interaction_parameters( ifile_pmt )

    ! sort molecules into types
    call gen_molecule_type_data( xyz, n_mole, n_atom, molecule_name, mass, r_com )

    ! initialize mass to negative values to make sure we read all of these in
    atype_mass=-1d0
    ! read in topology information, zero molecule_exclusions here, as it is possible that
    ! exclusions are being read in from topology file
    molecule_exclusions=0
    call read_topology_file( ifile_top )
    ! mass is read in from topology file, fill in mass for each molecule
    call fill_mass( tot_n_mole, n_atom, mass )

    ! center of mass
    call update_r_com( xyz, n_mole, n_atom, mass, r_com )
    ! if framework, update r_com for framework
    if (framework_simulation .eq. 1) then
       do i_mole = n_mole+1, tot_n_mole
          r_com(i_mole,:) = xyz(i_mole,1,:)
       end do
    endif
    ! center of mass of molecule might be outside of box after calling subroutine check_intra_molecular_shifts, fix this
    do i_mole=1,n_mole
       call shift_move(n_atom(i_mole),xyz(i_mole,:,:),r_com(i_mole,:),box)
    enddo

    ! generate exclusions, note that some exclusions may have already been explicity read in
    ! from topology file.
    call generate_intramolecular_exclusions


    if (electrostatic_type .ne. "none") then
       ! warn if any solute molecules are charged
!!$       do i_mole = 1, n_mole
!!$          sum_chg=0d0
!!$          do i_atom=1,n_atom(i_mole)
!!$             sum_chg= sum_chg + chg(i_mole,i_atom)
!!$          end do
!!$          if ( abs(sum_chg) > solute_chg_thresh ) then
!!$             write(*,*) "WARNING molecule ", i_mole," has total charge ", sum_chg
!!$          endif
!!$       enddo
       ! check to make sure there is net zero charge in the unit cell, if we are using electrostatics
       if (framework_simulation .eq. 1) then
          sum_chg=0d0
          do i_mole = n_mole+1, tot_n_mole
             sum_chg= sum_chg + chg(i_mole,1)
          end do
          if ( abs(sum_chg) > framework_chg_thresh ) then
             write(*,*) "total charge on framework is", sum_chg
             stop " unit cell is not neutral! "
          endif
       endif
    endif


    ! make sure to store framework charges, and framework atom numbers
    tot_chg = chg
    n_atom_drude = n_atom
    ! if simulation has drude oscillators, correct atom charges, and intialize mapping between
    ! drude oscillators and their corresponding atoms
    ! even if this is not a drude oscillator simulation, call this subroutine as it will
    ! fill drude_atom_map trivially for the atoms, which will be used in screening functions
    call initialize_drude_oscillators(tot_n_mole, n_mole, n_atom, n_atom_drude, tot_chg , chg_drude)


    ! initialize verlet list if we are going to use it
    Select Case(use_verlet_list)
    Case("yes")
       ! make sure MAX_N_MOLE is at least n_mole + 1, since we use verlet_point(n_mole+1) to
       ! find finish position of verlet list for n_mole
       if ( ( n_mole + 1 ) > MAX_N_MOLE ) then
          write(*,*) "for verlet list, MAX_N_MOLE must be at least as big as "
          write(*,*) "n_mole + 1"
          stop
       end if
       call allocate_verlet_list( tot_n_mole, n_atom, tot_chg, box )
       call construct_verlet_list(tot_n_mole, n_mole, n_atom, xyz, box )
       ! the "1" input to update_verlet_displacements signals to initialize the displacement array
       call update_verlet_displacements( n_mole, n_atom, xyz, box, flag_junk, 1 )
    End Select

    ! initialize ms_evb simulation
    Select Case(ms_evb_simulation)
    Case("yes")
       call initialize_evb_simulation( n_mole, n_atom, ifile_top )
       ! we are storing dQ_dr grid for ms-evb, allocate data arrays for this

       ! note that if we are using a verlet list, we could use verlet_elec_atoms here.
       ! we generate total_atoms_list to use here just in case we aren't using a verlet list
       call generate_verlet_atom_index( total_atoms_list, atom_molecule_map, tot_n_mole, n_atom, tot_chg )
       ! grid size:  each atom has spline_order^3 non-zero derivatives in the Q grid
       size = total_atoms_list * spline_order**3
       ! size could be big depending on the system, make sure we don't ask for too much memory
       if ( size > 10000000 ) then
	  write(*,*) "are you sure you want to store dQ_dr?  This requires allocating an array"
	  write(*,*) " bigger than 3 x ", size
          stop
       end if
       ! 5 indices for 2nd dimension of dQ_dr_index, for molecule index, atom index, and 3 Qgrid indices
       ! allocate column major
       !allocate( dQ_dr(3,size), dQ_dr_index(5,size) )
       size=spline_order**3
       allocate( dQ_dr(3,size,total_atoms_list) , dQ_dr_index(3,size,total_atoms_list) )
       allocate(atom_list_force_recip(3,total_atoms_list) )
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
          !         initialize g functions used in dispersion PME, by Kuang
          do i = 1, gfun_grid
             x = gfun_max/dble(gfun_grid)*dble(i)
             exp_x2 = dexp(-x**2)
             g6_table(i) = (x**4/2d0+x**2+1d0)*exp_x2
             g8_table(i) = (x**6+3d0*x**4+6d0*x**2+6d0)/6d0*exp_x2
             g10_table(i) = (x**8+4d0*x**6+12d0*x**4+24d0*x**2+24d0)/24d0*exp_x2
             g12_table(i) = (x**10+5d0*x**8+20d0*x**6+60d0*x**4+120d0*x**2+120d0)/120d0*exp_x2
             g14_table(i) = (x**12+6d0*x**10+30d0*x**8+120d0*x**6+360d0*x**4+720d0*x**2+720d0)/720d0*exp_x2
          enddo
       end select

!!!!!!!!!!!!!!!!!!pme will add framework reciprocal space electrostatic interactions, remove this here
       Ewald_framework = 0D0
       if(framework_simulation .eq. 1 ) then
          call remove_framework_elec_recip( tot_n_mole, n_mole, n_atom, xyz, chg, box, dfti_desc,dfti_desc_inv )
       endif

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
  ! this subroutine reads crystal framework coordinates
  ! for a GCMC simulation
  !
  ! code detects whether one line is read for box vectors
  ! for a cubic box, or three lines are read for a non-cubic box
  !
  ! this will also read charges if those are given after the coordinates.  This is detected
  ! by the number of arguments on the atom line
  !**************************************************!
  subroutine read_fconf( ifile_fconf, n_mole, tot_n_mole, n_atom, xyz, chg, flag_charge, alist, box )
    use routines
    use global_variables
    character(*),intent(in) :: ifile_fconf
    character(*), dimension(:,:),intent(inout) :: alist
    integer,intent(in)::n_mole
    integer,intent(out)::tot_n_mole
    integer,dimension(:),intent(inout)::n_atom
    real*8,dimension(:,:,:),intent(inout)::xyz
    real*8,dimension(:,:),intent(out) :: chg
    integer  , intent(out)         :: flag_charge
    real*8,dimension(3,3),intent(out)::box

    real*8 :: tmp1, tmp2, tmp3
    character(100) :: line
    ! give args enough dimensions to accomodate all the arguments ( > 5) 
    character(20),dimension(10)  :: args
    integer:: file_h=15, i_mole, n_frwk_atoms,i, ios, nargs


    ! need to modify this to account for new total_atoms variable
    write(*,*) "need to modify read_fconf subroutine to account for total_atoms variable"
    stop

    open(unit=file_h,file=ifile_fconf,status="old")
    read(file_h,*) n_frwk_atoms
    ! see how many arguments are given on the atom line
    read(file_h,'(A)') line
    call parse( line," ", args, nargs )
    backspace(file_h)

    Select Case( nargs )
    Case(4)
       flag_charge=0
    Case(5)
       flag_charge=1
    case default
       write(*,*) "unrecognized file format in file ", ifile_fconf
       stop
    End Select
    chg=0d0

    do i_mole=1,n_frwk_atoms
       Select Case( flag_charge )
       Case(0)
          ! no charges
          read(file_h,*) alist( n_mole+i_mole,1),xyz(n_mole+i_mole,1,1),xyz(n_mole+i_mole,1,2),xyz(n_mole+i_mole,1,3)
       Case(1)
          !charges
          read(file_h,*) alist( n_mole+i_mole,1), xyz(n_mole+i_mole,1,1),xyz(n_mole+i_mole,1,2),xyz(n_mole+i_mole,1,3), chg(n_mole+i_mole,1)
       end Select
       n_atom( n_mole+i_mole) = 1
       call trim_end( alist( n_mole+i_mole,1) )
    enddo

    read(file_h,*) tmp1, tmp2, tmp3
    read( file_h, *, iostat=ios ) box(2,1), box(2,2), box(2,3)
    if ( ios < 0 ) then ! end of the file
       box = 0.0d0
       box(1,1) = tmp1
       box(2,2) = tmp2
       box(3,3) = tmp3
    else
       box(1,1) = tmp1
       box(1,2) = tmp2
       box(1,3) = tmp3
       read( file_h, * ) box(3,1), box(3,2), box(3,3)
    endif


    close(file_h)

    tot_n_mole = n_mole + n_frwk_atoms

    if ( tot_n_mole > MAX_N_MOLE ) then
       stop "please increase value of parameter MAX_N_MOLE in glob_v.f90 file to accomodate the size of your system"
    endif

  end subroutine read_fconf


  !**************************************************!
  ! this subroutine reads atom_type parameters from input file for a non SAPT-based force field
  !*************************************************
  subroutine read_param( ifile_pmt, gen_cross_terms)
    use global_variables
    character(*),intent(in)::ifile_pmt
    integer,dimension(:,:),intent(out) :: gen_cross_terms

    integer::file_h=15, i_type, j_type, i_param, n_cross, inputstatus,i_search,j_search,ind,ind1,ind2
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

    Select Case(framework_simulation)
    Case(1)
       stop "this code has been commented"
!!$       ! get framework atom types
!!$       do
!!$          Read(file_h,'(A)',Iostat=inputstatus) line
!!$          If(inputstatus < 0) Exit
!!$          ind=INDEX(line,'framework_species')
!!$          IF(ind .NE. 0) Exit
!!$       enddo
!!$       Read(file_h,'(A)',Iostat=inputstatus) line
!!$       read(file_h,*) n_type_f
!!$       do j_mole=1, n_type_f
!!$          i_mole= n_type_s + j_mole
!!$          ! make sure we are not exceeding array length
!!$          if ( i_mole > MAX_N_ATOM ) then
!!$             stop "number of atom type g.t. MAX_N_ATOM.  Increase this value"
!!$          endif
!!$          read(file_h,*) atype_name(i_mole),atype_chg(i_mole),atype_lj_parameter(i_mole,i_mole,1),atype_lj_parameter(i_mole,i_mole,2),atype_lj_parameter(i_mole,i_mole,3),atype_pol(i_mole)
!!$       enddo
!!$
!!$       n_atom_type =n_type_f + n_type_s
    Case(0)
       ! no framework

    End Select


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
  subroutine gen_param(alist,tot_n_mole,n_atom, chg, chg_drude,gen_cross_terms)
    use global_variables
    character(*), dimension(:,:),intent(in) :: alist
    real*8, dimension(:,:),intent(out) :: chg, chg_drude
    integer,dimension(:,:),intent(in) :: gen_cross_terms
    integer,dimension(:),intent(in):: n_atom
    integer,intent(in)::tot_n_mole

    integer:: i_param, j_param, i_mole , i_atom, ind, i_type,flag
    real*8,parameter:: small = 1D-10

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
    do i_mole = 1, tot_n_mole
       do i_atom =1, n_atom(i_mole)
          flag=0
          do  i_type=1, n_atom_type
             if (atype_name(i_type) .eq. alist(i_mole,i_atom) ) then
                flag=1
                i_param = i_type
                exit 
             endif
          enddo
          ! make sure all atom types have parameters
          if (flag .eq. 0 ) then
             write(*,*) "atom type    ", alist(i_mole,i_atom), "doesn't have force field parameters!"
             stop
          endif

          ! set index, chg, and polarizability
          atom_index(i_mole,i_atom) = i_param
          chg(i_mole,i_atom) = atype_chg(i_param)
          chg_drude(i_mole,i_atom) = -sqrt(atype_pol(i_param)*springcon)

          ! make sure all polarizabilities are zero if drude_simulation = 0 , so that input is consistent
          Select Case(drude_simulation)
          Case(0)
             if ( abs( chg_drude(i_mole,i_atom) ) > small ) then
                write(*,*)  ""
                write(*,*)  " polarizability on atom ", i_atom , " on molecule ", i_mole
                write(*,*)  " is non-zero, yet drude_simulation is set to 0 "
                write(*,*)  " if you want to run a drude oscillator simulation, set "
                write(*,*)  " drude_simulation to 1 ! "
                write(*,*) ""
                stop
             endif
          End Select

       enddo
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


  !***************************************************
  !  This subroutine reads C9 coefficients for three body dispersion from
  !  parameter file
  !
  !  for N atom types, the C9 coefficients should be input as follows:
  !  There should be N + (N-1) + (N-2) ... lines (records) containing
  !  C9 parameters.  These parameters are read in as 
  !  loop (i=1, N)
  !    loop ( j=i,N)
  !      loop ( k=j,N)
  !
  !  where the parameters of the last loop are on one record,
  !  and the first two loops account for the number of lines
  !
  !  these parameters should be given in units of (kJ/mol)*A^9, consistent with 
  !  the rest of the parameters
  !***************************************************
  subroutine read_C9_three_body_dispersion_parameters(ifile_pmt)
    use global_variables
    character(*),intent(in) :: ifile_pmt

    character(50)::line
    integer :: i, j, inputstatus,ind, file_h=15, itype,jtype,ktype, ordered_index(3),store(3)

    open(unit=file_h,file=ifile_pmt,status="old")

    do
       Read(file_h,'(A)',Iostat=inputstatus) line
       If(inputstatus < 0) Exit
       ind=INDEX(line,'three body dispersion')
       IF(ind .NE. 0) Exit
    enddo

    ! if we didn't find the three body dispersion section in parameter file, then stop!
    if ( ind .eq. 0 ) then
       stop "couldn't find three body dispersion section in parameter input file"
    endif

    do itype = 1 , n_atom_type
       do jtype = itype, n_atom_type
          read(file_h,*) ( atype_3body_C9(itype,jtype,ktype),ktype=jtype, n_atom_type )
       enddo
    enddo

    ! now fill in the rest of the array

    do itype = 1 , n_atom_type
       do jtype = 1 , n_atom_type
          do ktype =1 , n_atom_type

             store(1)=itype; store(2)=jtype; store(3)=ktype
             ! need to arrange these indices in increasing order.  Probably better ways to do this but...
             do i=1,3
                ordered_index(i) = min( store(1), store(2), store(3) )

                ! now give this index a higher value so we can keep using min() function
                do j=1,3
                   if ( store(j) .eq. ordered_index(i) ) then
                      store(j) = n_atom_type + 1
                      exit
                   endif
                enddo
             enddo

             ! now fill in array with C9
             atype_3body_C9(itype,jtype,ktype) = atype_3body_C9(ordered_index(1), ordered_index(2), ordered_index(3) )

          enddo
       enddo
    enddo


  end subroutine read_C9_three_body_dispersion_parameters



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
  ! this subroutine generates molecule_type arrays to be used for particle insertion
  ! notice this only includes solute molecules, as this is all we need
  !
  ! this subroutine sets global variables:: n_molecule_type, molecule_type
  !*****************************************************!
  subroutine gen_molecule_type_data( xyz, n_mole, n_atom, molecule_name, mass, r_com )
    use global_variables
    integer, intent(in) :: n_mole
    integer, intent(in), dimension(:) :: n_atom
    character(*),dimension(:),intent(in) :: molecule_name
    real*8, intent(in), dimension(:,:) :: mass
    real*8, intent(in), dimension(:,:) :: r_com
    real*8, intent(in), dimension(:,:,:) :: xyz

    integer:: i_mole,j_mole,i_atom, old_type,count, flag_name


    ! first molecule type is first molecule in input file
    i_mole = 1
    n_molecule_type = 1
    do i_atom=1,n_atom(i_mole)
       molecule_type(n_molecule_type,i_atom) = atom_index(i_mole,i_atom)
    enddo
    ! mark end
    if ( n_atom(i_mole) < MAX_N_ATOM ) then
       molecule_type(n_molecule_type,n_atom(i_mole)+1) = MAX_N_ATOM_TYPE + 1
    end if

    molecule_type_name(n_molecule_type)=molecule_name(i_mole)
    molecule_index(i_mole)= n_molecule_type

    ! loop over all solute molecules
    do i_mole =2, n_mole
       old_type =0
       flag_name=0
       ! loop over all stored molecule types to see if this solute molecule is a new type
       do j_mole =1, n_molecule_type
          ! see if same name
          if ( molecule_type_name(j_mole) .eq. molecule_name(i_mole) ) then
             flag_name = 1
          end if

          ! now check to see if they have same number of atoms
          do i_atom =1, MAX_N_ATOM
             if( molecule_type(j_mole,i_atom) .eq. MAX_N_ATOM_TYPE + 1) exit
             count = i_atom
          enddo
          if( n_atom(i_mole) .eq. count) then
             do i_atom=1, n_atom(i_mole)
                if (atom_index(i_mole,i_atom) .ne. molecule_type(j_mole,i_atom)) then
                   goto 100
                endif
             enddo
             ! make sure that we haven't just matched a smaller molecule to part of a larger molecule
             if ( n_atom(i_mole) < MAX_N_ATOM ) then
                if ( molecule_type(j_mole,n_atom(i_mole)+1 ) .ne. ( MAX_N_ATOM_TYPE + 1 ) ) then
                   goto 100
                endif
             end if
             ! make sure these molecules have the same name to avoid confusion, note this is different on the name
             ! check at the beginning of the loop, because that checks to see if any previous molecule
             ! had the same name, not just this specific molecule
             if ( molecule_type_name(j_mole) .ne. molecule_name(i_mole) ) then
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

          do i_atom=1,n_atom(i_mole)
             molecule_type(n_molecule_type,i_atom) = atom_index(i_mole,i_atom)
          enddo
          ! mark end
          if ( n_atom(i_mole) < MAX_N_ATOM ) then
             molecule_type(n_molecule_type,n_atom(i_mole)+1) = MAX_N_ATOM_TYPE + 1
          end if
          molecule_type_name(n_molecule_type) = molecule_name(i_mole)
          molecule_index(i_mole)= n_molecule_type

       endif
    enddo


  end subroutine gen_molecule_type_data



  !***************************************
  ! this subroutine fills in mass data structure array
  ! using atom-type masses that are stored as global variables
  !***************************************
  subroutine fill_mass( tot_n_mole, n_atom, mass )
    use global_variables
    integer, intent(in) :: tot_n_mole
    integer,dimension(:),intent(in) :: n_atom
    real*8,dimension(:,:),intent(out) :: mass

    integer :: i_mole, i_atom, i_type
    real*8  :: mass_atom

    do i_mole=1,tot_n_mole
       do i_atom=1, n_atom(i_mole)
          ! get atom type
          i_type = atom_index(i_mole,i_atom)
          mass_atom = atype_mass(i_type)
          if ( mass_atom < 0d0 ) then
             write(*,*) "mass for atom type ", atype_name(i_type)
             write(*,*) "has not been read in from topology file"
             stop
          end if
          mass(i_mole,i_atom) = mass_atom
       enddo
    enddo

  end subroutine fill_mass




  !******************************************************!
  ! this subroutine calculates reciprocal space framework electrostatic energy so that it can be subtracted out
  ! in a GCMC simulation
  !*****************************************************!
  subroutine remove_framework_elec_recip( tot_n_mole, n_mole, n_atom, xyz, chg, box, dfti_desc,dfti_desc_inv )
    use global_variables
    use MKL_DFTI
    use pme_routines
    use electrostatic
    integer, intent(in) :: n_mole,tot_n_mole
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: chg, box
    real*8, intent(in), dimension(:,:,:) :: xyz
    TYPE(DFTI_DESCRIPTOR),pointer,intent(in):: dfti_desc,dfti_desc_inv

    real*8, dimension(:,:,:),allocatable :: xyz_framework
    integer,dimension(:),allocatable :: n_atom_framework
    real*8, dimension(:,:),allocatable :: chg_framework
    real*8     :: Ewald_self_store
    integer:: i_mole, n_framework

    n_framework = tot_n_mole - n_mole


    allocate(xyz_framework(n_framework,1,3),n_atom_framework(n_framework),chg_framework(n_framework,1))
    n_atom_framework = 1

    do i_mole =1, n_framework
       xyz_framework(i_mole,1,:) = xyz(n_mole+i_mole,1,:)
       chg_framework(i_mole,1) = chg(n_mole+i_mole, 1)
    enddo


    Ewald_framework = pme_recip( n_framework, n_atom_framework, xyz_framework, chg_framework, box, dfti_desc,dfti_desc_inv,1 )


    Ewald_self_store = Ewald_self
    call update_Ewald_self( n_framework, n_atom_framework, chg_framework )
    Ewald_framework = Ewald_framework + Ewald_self

    Ewald_self = Ewald_self_store

  end subroutine remove_framework_elec_recip



end module mc_routines
