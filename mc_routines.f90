!***********************************************************************
!
!  This module contains various subroutines used in the Monte Carlo simulation,
!  including data initialization subrouines, as well as subroutines that
!  (for a Lennard-Jones type force field) read in the input files
!
!  NOTE : For a simple Lennard-Jones force field, the input files and data setup
!  is relatively simple, and is carried out mainly by the subroutines
!  read_param and gen_param found in this module
!
!  for our SAPT-based force fields, the data setup is more complicated and is
!  carried out by subroutines found in the sapt_ff_routines module
!
!***********************************************************************
module mc_routines
  use routines
  use rigid_body
  use sapt_ff_routines
  use ms_evb
  implicit none

contains

  !*************************************************************************
  ! this routine reads the input files for the solute, and if present, the framework
  ! and initializes the simulation
  !*************************************************************************
  subroutine initialize_simulation(box,alist,n_mole,tot_n_mole,n_atom,n_atom_drude,xyz,quaternion,mass,r_com,chg,chg_drude,tot_chg,ifile_conf,ifile_fconf,ifile_pmt)
    use global_variables
    use electrostatic
    use pairwise_interaction
    use pme_routines
    character(MAX_FN),intent(in) :: ifile_conf,ifile_fconf,ifile_pmt
    character(MAX_ANAME), dimension(MAX_N_MOLE,MAX_N_ATOM),intent(out) :: alist
    integer,intent(out)::n_mole,tot_n_mole
    integer,dimension(:),intent(out)::n_atom,n_atom_drude
    real*8,dimension(:,:,:),intent(out)::xyz
    real*8,dimension(:,:), intent(out):: quaternion
    real*8,dimension(:,:),intent(out)::mass,chg,chg_drude,tot_chg
    real*8,dimension(:,:),intent(out)::box,r_com
    integer, dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE) :: gen_cross_terms
    real*8,dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE,9) :: temp_lj_parameter
    real*8,dimension(MAX_N_MOLE,MAX_N_ATOM) :: temp_chg
    integer:: i,i_mole,i_atom,n_type_s,index,is_overlap, flag_junk, flag_charge
    real*8 :: unit_cell(3,3),sum_chg
    real*8,dimension(3) :: inertia
    real*8,parameter :: framework_chg_thresh = 1D-5 ,solute_chg_thresh = 1D-6 , small_stop=1d0, small_warn=1d-1
    integer :: warn_inertia=0


    ! initialize random seed, this will be overwritten if seed is read in in read_conf routine
    call initialize_random_seed

    !************************************** get parameters****************************************************!
    ! get parameters for all atom types, assuming there are no more than MAX_N_ATOM_TYPE types of atoms total
    ! remember no two atom names can conflict, not even those of framework and solute
    !*****************************************************************************************************!

    ! if sapt-type force field, need to read in more complicated parameter file.  
    Select Case(sapt_type_ff)       ! this determines input file type
    Case("no")
       call read_param( ifile_pmt,gen_cross_terms)
    Case("yes")
       ! temp_lj_parameter array is used to temporarily store solute-solute parameters
       ! if these are different than solute-framework
       call read_param_decomp( ifile_pmt , temp_lj_parameter, n_type_s )
    end Select
    ! make sure we don't have too many atom types
    if ( n_atom_type > MAX_N_ATOM_TYPE ) then
       stop "number of atom types g.t. MAX_N_ATOM_TYPE.  Increase this value"
    endif


    ! **************************** coordinates and box size************************************************!
    ! zero n_atom array values that wont be used, as maxval(n_atom) will used to allocate an array later on
    n_atom=0
    n_atom_drude=0

    ! get solute coordinates
    call read_conf( ifile_conf, n_mole, n_atom, xyz, alist, box )
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
    if(is_overlap .eq. 1) then
       write(*,*) ""
       write(*,*) "In the initial configuration, there exists atoms on two or more"
       write(*,*) "molecules that are closer than ", too_close, " Angstroms."
       write(*,*) "We use a hard wall potential at this distance, as determined"
       write(*,*) "by the parameter 'too_close' , and therefore the initial configuration"
       write(*,*) "must not contain atom-atom contacts this close.  "
       write(*,*) ""
       write(*,*) "Please use a different initial configuration, or change the value of "
       write(*,*) "'too_close', only if you believe this is a reasonable separation distance"
       write(*,*) ""
       write(*,*) "(NOTE : the hard-wall potential is important in a Monte Carlo simulation"
       write(*,*) "to avoid moves into close-contact, divergent regions of a Buckingham potential)"
       write(*,*) ""
       stop
    endif


    !***************** make sure cutoff distances are fine for box size, and that box type is supported
    call check_cutoffs_box( box )

    ! ********************* generate cross terms and fill in parameter arrays****************************!
    Select Case(sapt_type_ff)
    Case("no")
       ! fill in total parameter arrays that are stored in global_variables
       call gen_param(alist,tot_n_mole,n_atom, chg, chg_drude,gen_cross_terms)
    Case("yes")
       call gen_param_decomp(alist,tot_n_mole,n_atom, chg, chg_drude , temp_lj_parameter , n_type_s)
    end Select


    !******************** fill in framework charges with charges from the .fconf file if those have been read in there
    Select Case(flag_charge)
    Case(1)
       write(*,*) ""
       write(*,*) "*******************************"
       write(*,*) "Found charges in .fconf file.  These charges will be"
       write(*,*) "used in simulation and corresponding charges in     "
       write(*,*) "*.pmt file will be ignored"
       write(*,*) "*******************************"
       write(*,*) ""

       do i_mole=n_mole+1, tot_n_mole
          chg(i_mole,:) = temp_chg(i_mole,:)
       enddo
    End Select

    ! now get special 1-4 interactions from parameter file if needed
    call read_generate_14_interaction_parameters( ifile_pmt )

    !************************ generate masses from element type *****************************************!

!!!!!!!!zero for drude oscillators
    mass=0.
    ! generate masses for solute and framework atoms, even though framework is rigid, we need the masses on these atoms if we are using a Feynmann-Hibbs quantum correction
    do i_mole=1,tot_n_mole
       call gen_mass( alist(i_mole,:), n_atom(i_mole), mass(i_mole,:) )
    enddo

    ! generate atom type masses
    call gen_mass( atype_name, n_atom_type, atype_mass )

    !  careful using massage if more than co2
    !  call massage( xyz, n_mole, n_atom, r_mole_com, mass ) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

    ! if code is being used for replica exchange, and there are zero molecules initially, molecule information has been read in 
    Select Case(replica_input)
    Case("yes")
       if ( n_mole .ne. 0 ) then
          call gen_molecule_type_data( xyz, n_mole, n_atom, mass, r_com )
       endif
    Case default
       call gen_molecule_type_data( xyz, n_mole, n_atom, mass, r_com )
    End Select

    Select Case( flexible_molecule_simulation )
    Case("no")
       ! generate quaternions for molecule orientations
       do i_mole = 1, n_mole
          call gen_inertia( xyz, i_mole, n_atom, mass ,r_com, inertia, quaternion(i_mole,:) )
          ! check inertia tensor
          index = molecule_index(i_mole)

          if ( ( abs( molecule_inertia(index,1) - inertia(1) ) > small_stop ) .or. ( abs( molecule_inertia(index,2) - inertia(2) ) > small_stop ) .or. ( abs( molecule_inertia(index,3) - inertia(3) ) > small_stop ) ) then
             write(*,*) "mismatch of inertia components between particular molecule and stored molecule type inertia tensor"
             write(*,*) "moments of inertia for molecule ", i_mole
             write(*,*) inertia
             write(*,*) "moments of inertia for corresponding molecule type ", index
             write(*,*) molecule_inertia(index,:)
!!$          stop
          else if( ( abs( molecule_inertia(index,1) - inertia(1) ) > small_warn ) .or. ( abs( molecule_inertia(index,2) - inertia(2) ) > small_warn ) .or. ( abs( molecule_inertia(index,3) - inertia(3) ) > small_warn ) ) then
             if ( warn_inertia .eq. 0 ) then
                write(*,*) "geometries of identical molecules may not be exactly equivalent as principle moments of inertia deviate by as much as ", small_warn, " (but deviate less than ", small_stop, ")"
                warn_inertia=1
             endif
          endif

       end do
    End Select


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
          call initialize_evb_simulation( n_mole, n_atom )
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
  subroutine initialize_energy_force(log_file,box,n_mole,tot_n_mole,n_atom,n_atom_drude,alist,xyz,r_com,mass,chg,dfti_desc,dfti_desc_inv,iteration,potential,E_elec_nopol,E_elec,E_bh,E_3body,E_bond,E_angle,E_dihedral,xyz_drude,force_atoms,energy_components)
    use global_variables
    use electrostatic
    use insertion_bias_routines
    use penetration_ff_routines
    use pairwise_interaction
    use pme_routines
    use MKL_DFTI
    use total_energy_forces
    integer,intent(in)::log_file
    real*8,dimension(:,:),intent(inout)::chg
    integer,intent(in)::n_mole,tot_n_mole
    integer,dimension(:),intent(inout)::n_atom,n_atom_drude
    character(*),dimension(:,:),intent(in) :: alist
    real*8,dimension(:,:,:),intent(inout)::xyz
    real*8,dimension(:,:),intent(in) :: box
    real*8,dimension(:,:),intent(inout) :: r_com,mass
    TYPE(DFTI_DESCRIPTOR), pointer,intent(out):: dfti_desc,dfti_desc_inv
    integer,intent(out) :: iteration
    real*8,intent(out) :: potential, E_elec_nopol, E_elec, E_bh, E_3body,E_bond,E_angle,E_dihedral
    real*8, dimension(:,:,:),intent(out) :: xyz_drude,force_atoms
    real*8,dimension(:),intent(out) :: energy_components

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

    ! initialize atype_damp_dispersion array
    call initialize_atype_damp_dispersion


    !************ grid ewald real space interactions
    Select Case(grid_ewald_realspace_interaction)
    Case("yes")
       call construct_ewald_realspace_grid
    End Select


    !***************** initialize pme grid algorithm if using pme for dispersion ********************!

    if ( pme_disp == "yes" ) then
       !     the grid for long range dispersion interaction
       !     Count number of atom types in solute molecule
       do ntype_solute = 1, n_atom_type
          if ( atype_solute(ntype_solute) == 0 ) then
             goto 999
          endif
       enddo
999    continue
       ntype_solute = ntype_solute - 1

       allocate( lrdisp_pot(ntype_solute,pme_grid,pme_grid,pme_grid) )
       !     construct long-range dispersion grid
       lrdisp_pot = 0d0

       do atom_id1 = 1, ntype_solute
          call construct_recip_disp_grid(atom_id1,n_mole,tot_n_mole,n_atom,box,xyz,dfti_desc,dfti_desc_inv)
       enddo


    endif


    ! if we are using a penetration force field, set that up
    ! NOTE this uses damping functions, so we must do this after damping function tables have been created
    Select Case(penetration_force_field)
    Case("yes")
       call initialize_penetration_ff_routines
    Case("no")
       ! we are not using a penetration force field, set all distance thresholds to zero
       ! so that this will never be called
       atype_penetration_ff(:,:,2)=0d0
    End Select


    !  set long range correction and shifts
    ! this is in case there is no framework. If there is a framework, disp_lrc is automatically set to zero for obvious reasons
    call update_disp_lrc_noshift(n_mole, n_atom,alist, box)
    call update_lennard_jones_shift(box)


    !**************************** Get initial forces and energy *************************************!
    Select Case(ms_evb_simulation)
    Case("yes")
       call ms_evb_calculate_total_force_energy( force_atoms, potential, tot_n_mole, n_mole, n_atom, xyz, r_com, chg, mass, box, dfti_desc,dfti_desc_inv,log_file )
    Case("no")
       call calculate_total_force_energy( force_atoms,potential, E_elec,E_elec_nopol,E_bh, E_3body, E_bond,E_angle,E_dihedral,energy_components, iteration, tot_n_mole, n_mole, n_atom, n_atom_drude, r_com, xyz, chg, box, dfti_desc,dfti_desc_inv,log_file,xyz_drude)
    End Select
    !************************************************************************************************!


    !**************** since this is initial energy calculation, it's important that drude oscillators converged, stop otherwise
    Select Case(drude_simulation)
    Case(1)
       if ( iteration > iteration_limit ) then
          stop "Drude oscillators did not converge in the initial energy calculation.  Check that variables force_threshold and iteration_limit have consistent settings in glob_v.f90."
       endif
    End Select


    ! if puddle filling, intialize puddle_delta_E
    Select Case(puddle_filling_on)
    Case("yes")
       call puddle_fill_potential(n_mole,tot_n_mole,n_atom,n_atom_drude,xyz,box,r_com,chg,puddle_delta_E)
    end Select

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

    integer::file_h=15, i_mole,j_mole,n_type_s, n_type_f, n_cross, inputstatus,i_search,j_search,ind,ind1,ind2
    character(MAX_ANAME):: atom1, atom2
    real*8,dimension(6)::store
    character(25)::junk
    character(50)::line


    ! make sure the parameters seem reasonable for the type of force field
    call check_parameter_force_field_consistency( ifile_pmt )

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
    read(file_h,*) n_type_s

    do i_mole=1, n_type_s
       read(file_h,*) atype_name(i_mole),atype_chg(i_mole),atype_lj_parameter(i_mole,i_mole,1),atype_lj_parameter(i_mole,i_mole,2),atype_lj_parameter(i_mole,i_mole,3),atype_pol(i_mole)
       ! move spaces for matching
       call trim_end( atype_name(i_mole) )
    enddo

    Select Case(framework_simulation)
    Case(1)
       ! get framework atom types
       do
          Read(file_h,'(A)',Iostat=inputstatus) line
          If(inputstatus < 0) Exit
          ind=INDEX(line,'framework_species')
          IF(ind .NE. 0) Exit
       enddo
       Read(file_h,'(A)',Iostat=inputstatus) line
       read(file_h,*) n_type_f
       do j_mole=1, n_type_f
          i_mole= n_type_s + j_mole
          ! make sure we are not exceeding array length
          if ( i_mole > MAX_N_ATOM ) then
             stop "number of atom type g.t. MAX_N_ATOM.  Increase this value"
          endif
          read(file_h,*) atype_name(i_mole),atype_chg(i_mole),atype_lj_parameter(i_mole,i_mole,1),atype_lj_parameter(i_mole,i_mole,2),atype_lj_parameter(i_mole,i_mole,3),atype_pol(i_mole)
       enddo

       n_atom_type =n_type_f + n_type_s

    Case(0)
       ! no framework
       n_atom_type = n_type_s

    End Select

    ! fill in atype_solute array
    do i_mole=1,n_atom_type
       if ( i_mole .le. n_type_s ) then
          atype_solute(i_mole) = 1
       else
          atype_solute(i_mole) = 0
       endif
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
          do i_mole=1, n_cross
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
             end Select
             gen_cross_terms(ind1,ind2)=1
             gen_cross_terms(ind2,ind1)=1
          enddo
       endif
    endif

    ! finally make sure we don't have two of the same atom type defined

    do i_mole=1, n_atom_type - 1
       do j_mole = i_mole + 1 , n_atom_type
          if ( atype_name(i_mole) .eq. atype_name(j_mole) ) then
             write(*,*) "atomic parameters defined more than once for atom type", atype_name(i_mole)
             stop
          endif
       enddo
    enddo

    close(file_h)

  end subroutine read_param



  !********************************************
  ! this subroutine checks the consistency between
  ! the requested force field type, and the parameter
  ! file being read in.  If any inconsistencies are found,
  ! the code will stop with an error message
  !
  ! For example, if a buckingham force field is intended,
  ! but the 3 paramters after the charge are 
  ! xxx.x  xxx.x   0.0
  ! then this is probably an inconsistency, since
  ! that would be the input for a lennard jones potential
  ! ( epsilon, sigma, junk )
  !*********************************************
  subroutine check_parameter_force_field_consistency( ifile_pmt )
    use global_variables
    character(*),intent(in)::ifile_pmt

    real*8,parameter :: small = 1D-3
    integer::file_h=15, i_atom,i,n_type_s, inputstatus, ind, nargs
    integer,parameter :: max_param=20
    integer,parameter :: lj_bkghm_param_number=6
    character(MAX_ANAME):: atomname
    character(300) :: input_string
    character(50)::line
    character(50),dimension(max_param)::args
    real*8,dimension(5):: testparam



    open(unit=file_h,file=ifile_pmt,status="old")

    ! get solute atom types
    do
       Read(file_h,'(A)',Iostat=inputstatus) line
       If(inputstatus < 0) Exit
       ind=INDEX(line,'solute_species')
       IF(ind .NE. 0) Exit
    enddo
    Read(file_h,'(A)',Iostat=inputstatus) line
    read(file_h,*) n_type_s


    ! make sure each atom-type line has correct number of parameters
    do i_atom=1, n_type_s
       read(file_h,'(A300)') input_string
       call parse(input_string," ",args,nargs)

       ! if wrong number of input parameters, stop
       if ( nargs .ne. lj_bkghm_param_number ) then
          write(*,*) ""
          write(*,'(A56,I3)') "Incorrect number of force field parameters for atomtype", i_atom
          write(*,'(A4,A40,A10)') "In ", ifile_pmt , "input file."
          write(*,'(A17,I3,A25)') "There should be ", lj_bkghm_param_number, " force field parameters"
          write(*,*) "for a standard lennard-jones or buckingham force field "
          write(*,*) "Please see the code manual for further details"
          write(*,*) ""
          stop
       endif

       ! make sure Lennard-Jones/Buckingham force field functional form has been correctly selected
       backspace(file_h)
       read(file_h,*) atomname, (testparam(i),i=1,5)

       Select Case(lj_bkghm)
       Case(1)
          ! Buckingham potential, 4th argument (3rd after charge) should be non-zero and positive
          if ( testparam(4) < small ) then
             write(*,*) ""
             write(*,'(A4,A40,A12)') "In ", ifile_pmt, " input file"
             write(*,'(A70,I3,A5)') "the 4th parameter (3rd after the charge) has a value <= 0.0 for the ", i_atom, " atom."
             write(*,*) "For a Buckingham potential, this should be a positive value for the "
             write(*,*) "C6 coefficient.", "If a lennard jones force field was intended, please"
             write(*,*) "set variable 'lj_bkghm' equal to 2 in the simulation_parameters.pmt input file"
             write(*,*) ""
             stop
          endif
       Case(2)
          if ( testparam(4) > small ) then
             ! LJ potential, 4th argument should be zero
             write(*,*) ""
             write(*,'(A4,A40,A12)') "In ", ifile_pmt, " input file"
             write(*,'(A70,I3,A5)') "the 4th parameter (3rd after the charge) has a value > 0.0 for the ", i_atom, " atom."
             write(*,*) "For a LJ potential, this parameter is meaningless and should be set to zero"
             write(*,*) "for transparency. If a Buckingham force field was intended, please set "
             write(*,*) "variable 'lj_bkghm' equal to 1 in the simulation_parameters.pmt input file"
             write(*,*) ""
             stop
          endif
       End Select

    enddo


    close(file_h)

  end subroutine check_parameter_force_field_consistency





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
  subroutine gen_molecule_type_data( xyz, n_mole, n_atom, mass, r_com )
    use global_variables
    integer, intent(in) :: n_mole
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: mass
    real*8, intent(in), dimension(:,:) :: r_com
    real*8, intent(in), dimension(:,:,:) :: xyz

    integer:: i_mole,j_mole,i_atom, old_type,count, flag_name
    real*8,dimension(3)::inertia


    molecule_mass=0d0

    ! first molecule type is first molecule in input file
    i_mole = 1
    n_molecule_type = 1
    do i_atom=1,n_atom(i_mole)
       molecule_type(n_molecule_type,i_atom) = atom_index(i_mole,i_atom)
    enddo
    ! mark end
    if ( n_atom(i_mole) < MAX_N_ATOM ) then
       molecule_type(n_molecule_type,n_atom(i_mole)+1) = MAX_N_ATOM + 1
    end if

    molecule_type_name(n_molecule_type)=molecule_name(i_mole)
    molecule_index(i_mole)= n_molecule_type

    ! molecule mass
    do i_atom = 1, n_atom( i_mole )
       molecule_mass(n_molecule_type)= molecule_mass(n_molecule_type) + mass(i_mole,i_atom)
    end do

    ! update r_mole_com array in global variables, stored in principle axis coordinates
    ! as well as molecule_inertia
    call update_relative_atom_pos( xyz, 1, 1, n_atom, mass,r_com )

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
             if( molecule_type(j_mole,i_atom) .eq. MAX_N_ATOM + 1) exit
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
                if ( molecule_type(j_mole,n_atom(i_mole)+1 ) .ne. ( MAX_N_ATOM + 1 ) ) then
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
             molecule_type(n_molecule_type,n_atom(i_mole)+1) = MAX_N_ATOM + 1
          end if
          molecule_type_name(n_molecule_type) = molecule_name(i_mole)
          molecule_index(i_mole)= n_molecule_type

          ! molecule mass
          do i_atom = 1, n_atom( i_mole )
             molecule_mass( n_molecule_type )= molecule_mass( n_molecule_type ) + mass(i_mole,i_atom)
          end do

          ! update r_mole_com array in global variables, stored in principle axis coordinates
          ! as well as molecule_inertia tensor
          call update_relative_atom_pos( xyz, i_mole, n_molecule_type, n_atom, mass,r_com )

       endif
    enddo


  end subroutine gen_molecule_type_data



  !******************************************************!
  ! this subroutine updates data arrays for a newly inserted particle
  ! in a GCMC simulation
  ! 
  ! this subroutine changes global variable:
  ! integer, dimension(MAX_N_MOLE,MAX_N_ATOM):: atom_index   
  !*****************************************************!
  subroutine gen_insertion_data(new_n_mole,new_alist,new_n_atom,new_n_atom_drude,new_r_com,new_xyz,new_quaternion,new_chg,new_mass,ins_type)
    use global_variables
    use egc_routines
    integer, intent(in)  :: new_n_mole, ins_type
    character(*), dimension(:,:),intent(inout) :: new_alist
    integer,dimension(:),intent(inout)::new_n_atom,new_n_atom_drude
    real*8,dimension(:,:),intent(in)::new_r_com
    real*8, dimension(:,:,:),intent(inout) :: new_xyz
    real*8,dimension(:,:),intent(inout)::new_quaternion
    real*8,dimension(:,:),intent(inout)::new_mass,new_chg

    integer:: i_atom, t_atom, n_atom_ins, index_drude
    real*8,dimension(3):: unit
    real*8,parameter::small=1D-6


    ! figure out how many atoms there are in the inserted molecule
    do i_atom =1, MAX_N_ATOM
       if( molecule_type(ins_type,i_atom) .eq. MAX_N_ATOM + 1) exit
       n_atom_ins = i_atom
    enddo
    new_n_atom(new_n_mole) = n_atom_ins
    ! see subroutine initialize_drude_oscillators for how to index drude oscillators
    index_drude = n_atom_ins + 1

    molecule_name(new_n_mole) = molecule_type_name(ins_type)

    ! transfer parameters from stored atom_type arrays to system arrays
    new_chg(new_n_mole,:) = 0.
    do i_atom =1, n_atom_ins
       t_atom = molecule_type(ins_type,i_atom)
       new_alist(new_n_mole,i_atom) = atype_name(t_atom)
       new_mass(new_n_mole,i_atom) = atype_mass(t_atom)
       atom_index(new_n_mole,i_atom) = t_atom

       new_chg(new_n_mole,i_atom) = atype_chg(t_atom)

       ! need to update drude_atom_map regardless of whether drude oscillators are being used
       ! map drude oscillators consistent with subroutine initialize_drude_oscillators

       ! trivial map from atom to atom
       drude_atom_map(new_n_mole,i_atom) = i_atom

       ! if drude simulation correct atomic charges and store drude charges
       if ( abs( atype_pol(t_atom)  ) > small )   then   
          drude_atom_map(new_n_mole,index_drude) = i_atom          

          new_chg(new_n_mole,i_atom) = atype_chg(t_atom)+ sqrt(atype_pol(t_atom)*springcon)
          new_chg(new_n_mole,index_drude) = - sqrt(atype_pol(t_atom)*springcon)

          index_drude = index_drude + 1
       endif

    enddo
    ! total number of atoms and drude oscillators on this molecule
    new_n_atom_drude(new_n_mole) = index_drude - 1


    ! if egc ensemble, must create new atom types and scale parameters for this new molecule
    Select Case(select_ensemble)
    Case("egc")
       ! change partial molecule index
       partial_molecule_index = new_n_mole
       ! this is called for insertion, so increment = 1
       call adjust_partial_molecule_parameters(1,new_n_atom,new_chg)
    End Select


    ! pick random orientation and create new molecular orientation
    Select Case(molecule_shape(ins_type))
    Case(2)
       ! this is just an atom
       new_xyz(new_n_mole,1,:) = new_r_com(new_n_mole,:)
    case default
       call gen_random_orientation(new_xyz,new_r_com,new_n_mole,ins_type,new_quaternion)
    End Select

  end subroutine gen_insertion_data



  !*********************************************************
  ! This subroutine generates an array which gives atom positions relative to cofm
  ! for non-linear molecules, these configurations are stored in their principle axis basis
  ! this changes array r_mole_com in global variables
  !
  ! this also fills in molecule_inertia array in global variables
  !*********************************************************
  subroutine update_relative_atom_pos( xyz, i_mole, i_index, n_atom, mass,r_com )
    use global_variables
    implicit none
    real*8, intent(in), dimension(:,:,:) :: xyz
    integer, intent(in) :: i_mole,i_index
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: mass
    real*8, intent(in), dimension(:,:) :: r_com

    integer :: i,j, i_atom, flag
    real*8,dimension(3)  :: delta_x,u_axis,inertia
    real*8::r_mag,r_magp
    real*8,parameter:: small = .1 , linear_criteria=.99
    real*8,dimension(4) :: quaternion
    real*8, dimension(3,3) :: a , a_trans

    ! if more than one atom, find distances to cofm
    if ( n_atom(i_mole) > 1 ) then

!!! check to see if the molecule is linear and create r_mole_com array
       do i_atom=1,n_atom(i_mole)
          delta_x(:)=xyz(i_mole,i_atom,:)-r_com(i_mole,:)
          if(dot_product(delta_x,delta_x) > small) exit
       enddo

       u_axis=delta_x/sqrt(dot_product(delta_x,delta_x))

       flag = 1
       do i_atom=1,n_atom(i_mole)
          delta_x(:)=xyz(i_mole,i_atom,:)-r_com(i_mole,:)

          r_mag = sqrt(dot_product(delta_x,delta_x))
          r_magp = dot_product(delta_x,u_axis)
          ! don't divide by zero
          if (r_mag > small ) then
             if ( ( abs(r_magp)/r_mag ) < linear_criteria ) then
                flag = 0
             endif
          endif
       enddo

       if ( flag == 1 ) then
          ! this is a linear molecule, generate moments of inertia, don't need quaternions,
          ! and don't need to store r_mole_com vector in any particular coordinate frame,
          ! so just use global axis
          molecule_shape(i_index) = 1
          call gen_inertia( xyz, i_mole, n_atom, mass ,r_com, inertia, quaternion)
          molecule_inertia(i_index,:) = inertia(:)

          do i_atom=1,n_atom(i_mole)
             delta_x(:)=xyz(i_mole,i_atom,:)-r_com(i_mole,:)
             r_mole_com(i_index,i_atom,:) = delta_x(:)
          enddo

       else
          ! this is a non-linear molecule, make sure we are not using rotations meant for linear molecules
          Select Case( hybrid_md_mc )
          Case("no")
             Select Case( rotate_method )
             Case(1)
                stop "cant use rotate_method 1 if there are non-linear molecules"
             End Select
          End Select
          ! also orientation bias is only implemented for linear molecules
          if (orientation_try .ne. 1 ) then
             stop " cannot use orientation bias for non-linear molecules"
          endif

          molecule_shape(i_index) = 0

          ! generate inertia tensor and using the quaternions from this, generate
          ! r_mole_com array in principle axis basis
          call gen_inertia( xyz, i_mole, n_atom, mass ,r_com, inertia, quaternion)
          molecule_inertia(i_index,:) = inertia(:)

          ! get transformation matrix to body coordinates from quaternions
          call quaternion_transformation_matrix(a,quaternion)



          ! store r_mole_com in principle axis coordinates
          do i_atom=1,n_atom(i_mole)
             delta_x(:)=xyz(i_mole,i_atom,:)-r_com(i_mole,:)
             delta_x = matmul(a,delta_x)
             r_mole_com(i_index,i_atom,:) = delta_x(:)
          enddo

       endif


       ! don't want to divide by zero here, if there is only one atom
    else
       molecule_shape(i_index) = 2
       r_mole_com(i_index,1,:) = 0d0
    endif

  end subroutine update_relative_atom_pos



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


  !************************************************************
  ! this subroutine is meant to be used to test certain parts of the code,
  ! so that by setting the global variable "debug" to different values, different
  ! output can be printed for different test jobs
  !************************************************************
  subroutine run_code_tests(log_file,box,n_mole,tot_n_mole,n_atom,n_atom_drude,alist,xyz,r_com,chg,dfti_desc,dfti_desc_inv)
    use global_variables
    use MKL_DFTI
    use electrostatic
    use explicit_three_body_interaction
    use insertion_bias_routines
    use penetration_ff_routines
    use pairwise_interaction
    use pme_routines
    use eq_drude
    use routines
    use total_energy_forces
    integer,intent(in)::log_file
    real*8,dimension(:,:),intent(in)::chg
    integer,intent(in)::n_mole,tot_n_mole
    integer,dimension(:),intent(in)::n_atom,n_atom_drude
    character(*),dimension(:,:),intent(in) :: alist
    real*8,dimension(:,:,:),intent(in)::xyz
    TYPE(DFTI_DESCRIPTOR),pointer,intent(in):: dfti_desc,dfti_desc_inv
    real*8,dimension(:,:),intent(in)::box,r_com

    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: force_atoms, xyz_drude_final
    integer,dimension(MAX_N_MOLE,MAX_N_ATOM)::drude_atoms
    real*8 :: E_elec, E_bh, E_elec_nopol, dr(3)
    integer :: i_mole, i_atom, i_drude, iteration


    !*******************debug==10******************************************!
    ! print electrostatic and lennard jones forces
    if( debug == 10 ) then
       ! electrostatic forces, drude oscillators
       if( drude_simulation .eq. 1 ) then
          call scf_drude(E_elec,E_elec_nopol,force_atoms,xyz,chg,r_com,box,tot_n_mole,n_mole,n_atom,n_atom_drude,iteration,xyz_drude_final,dfti_desc,dfti_desc_inv,log_file)
       else
          ! electrostatic forces, no drude oscillators
          do i_mole=1,n_mole
             do i_atom=1,n_atom(i_mole)
                drude_atoms(i_mole,i_atom)=i_atom
             enddo
          enddo
          call Electrostatic_energy_force( E_elec, force_atoms, drude_atoms,tot_n_mole, n_mole, n_atom, xyz, chg, box, r_com,dfti_desc,dfti_desc_inv)
       end if

       write(*,*) "************** Electrostatic Forces *********************"
       write(*,*) " i_mole , i_atom, Fx , Fy , Fz "
       do i_mole=1,n_mole
          do i_atom=1,n_atom(i_mole)
             write(*,*) i_mole, i_atom, force_atoms(i_mole,i_atom,:)
          enddo
       enddo

       ! lennard jones forces
       call lennard_jones( E_bh, force_atoms, tot_n_mole, n_mole, n_atom, xyz, box, r_com)

       write(*,*) "************** Lennard Jones Forces *********************"
       write(*,*) " i_mole , i_atom, Fx , Fy , Fz "
       do i_mole=1,n_mole
          do i_atom=1,n_atom(i_mole)
             write(*,*) i_mole, i_atom, force_atoms(i_mole,i_atom,:)
          enddo
       enddo

       !***********************************************************************!
    else if( debug == 11 ) then
       !*******************debug==11******************************************!
       ! print converged Drude oscillator positions
       if( drude_simulation .eq. 1 ) then
          call scf_drude(E_elec,E_elec_nopol,force_atoms,xyz,chg,r_com,box,tot_n_mole,n_mole,n_atom,n_atom_drude,iteration,xyz_drude_final,dfti_desc,dfti_desc_inv,log_file)
       else 
          stop "error debug was set to 11, this is only meaningful for a Drude oscillator simulation"
       endif

       write(*,*) "********** converged Drude oscillator positions (dr) **********"  
       write(*,*) " dr = xyzdrude - xyzatom "
       write(*,*) " i_mole, i_atom, dr(x) , dr(y) , dr(z) "

       do i_mole = 1, n_mole
          if ( n_atom_drude(i_mole) > n_atom(i_mole) ) then

             do i_drude = n_atom(i_mole)+1 , n_atom_drude(i_mole)
                ! use map to figure out which atom this oscillator is on
                i_atom = drude_atom_map(i_mole,i_drude)      

                dr(:) = xyz_drude_final(i_mole,i_drude,:) - xyz(i_mole,i_atom,:)

                write(*,*) i_mole, i_atom, dr
             enddo
          endif
       enddo
       !***********************************************************************!
    endif




    write(*,*) "A setting of debug > 2 signals the code to call the run_code_tests subroutine"
    write(*,*) "which will calculate and print various information based on the setting of debug"
    write(*,*) "this is meant to be used by test jobs.  The code stops after this call"
    stop


  end subroutine run_code_tests

end module mc_routines
