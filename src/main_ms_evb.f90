!**************************************************************************************
! This is the head source file for the MS-EVB MD program
! This code was written by Jesse McDaniel 2015-, while working for Professor Arun Yethiraj
!
! This code is based off of subroutines and modules from a Monte Carlo program written by Jesse McDaniel
! and Kuang Yu, 2010-2014, while working for Professor JR Schmidt.
!
! The manual for the Monte-Carlo code, examples of how to run it, and parameters for our
! SAPT-based force fields can be downloaded from the Schmidt group website
!
!              schmidt.chem.wisc.edu
!
!
!**************************************************************************************
program main_ms_evb
  use routines
  use global_variables
  use simulation_parameters
  use total_energy_forces
  use md_integration
  use initialize_routines
  implicit none

  !**************************************** these data structures are defined in the global_variables module
  Type(system_data_type)                               :: system_data
  Type(molecule_data_type), dimension(:), allocatable  :: molecule_data
  Type(atom_data_type)                                 :: atom_data  
  Type(integrator_data_type)                           :: integrator_data
  Type(file_io_data_type)                              :: file_io_data
  Type(verlet_list_data_type)                          :: verlet_list_data
  Type(PME_data_type)                                  :: PME_data
 
  !**************************************** target variables for atom_data_type pointers
   real*8, dimension(:,:), allocatable :: xyz
   real*8, dimension(:,:), allocatable :: velocity
   real*8, dimension(:,:), allocatable :: force
   real*8, dimension(:), allocatable   :: mass
   real*8, dimension(:), allocatable   :: charge
   integer, dimension(:), allocatable  :: atom_type_index
   character(MAX_ANAME)(:),allocatable :: aname


  !***** Local variables
  integer :: i_step, i, j, sample_vel=1, iteration
 

  call sort_input_files( file_io_data )
  call check_restart_trajectory( file_io_data )

  call read_simulation_parameters( file_io_data%ifile_simpmt_file_h, file_io_data%ifile_simpmt, system_data, integrator_data, verlet_list_data, PME_data )
  call initialize_simulation( system_data, molecule_data, atom_data, integrator_data, file_io_data, verlet_list_data, PME_data, xyz, velocity, force, mass, charge, atom_type_index, aname )


  !*** open output files if this is not a continuation run.  If this is a continuation, these files will
  ! already be open, and prepared to be appended to, in the check_restart_trajectory subroutine
  Select Case( restart_trajectory )
  Case("no")
     open( file_io_data%ofile_traj_file_h, file=file_io_data%ofile_traj, status='new' )
     open( file_io_data%ofile_log_file_h,  file=file_io_data%ofile_log, status='new' )
     Select Case( checkpoint_velocity )
     Case("yes")
        ! status may be new, or append to old
        open( file_io_data%ifile_velocity_file_h, file=file_io_data%ifile_velocity, status ='new')
     End Select
  end Select

  ! assume orthorhombic box
  if ( ( abs(box(1,2)) > 10D-6 ) .or. ( abs(box(1,3)) > 10D-6 ) .or. ( abs(box(2,1)) > 10D-6 ) .or.  ( abs(box(2,3)) > 10D-6 ) .or. ( abs(box(3,1)) > 10D-6 ) .or. ( abs(box(3,2)) > 10D-6 ) ) then
     write(*,*)  ""
     write(*,*)  "code has been modified to assume orthorhombic box"
     write(*,*)  "to change this, fix commented sections in energy calculation subroutines"
     write(*,*)  ""
     stop
  end if
  !***********

  !**********
  ! verlet list size
  write(*,*) " verlet list size parameter verlet_safe is set to ", safe_verlet
  write(*,*) " too big a value could slow simulation down "
  write(*,*) ""
  !************


  call initialize_energy_force( system_data, molecule_data, atom_data, verlet_list_data, PME_data )


  Select Case( restart_trajectory )
  Case("no")
     call sample_atomic_velocities(n_mole, n_atom, mass, vel_atom )
  End Select

  call calculate_kinetic_energy( new_KE, n_mole,n_atom,mass,vel_atom, constants%conv_kJmol_ang2ps2gmol )
  vol= volume( box(1,:),box(2,:),box(3,:) )

  ! only print here if not a restart
  Select Case( restart_trajectory )
  Case("no")
     call print_run_param( log_file , n_mole, vol )
     call print_result( traj_file, log_file, n_mole, n_atom, xyz, box, 0, E_elec,E_elec_nopol,E_bh, E_3body, E_bond, E_angle, E_dihedral,potential,new_KE,vol )
  End Select



  do i_step = 1, n_step

     trajectory_step = n_old_trajectory + i_step

     call  mc_sample(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,force_atoms,vel_atom,mass,r_com,potential, E_elec,E_elec_nopol,E_bh,E_3body,E_bond,E_angle,E_dihedral,tot_chg,dfti_desc,dfti_desc_inv,log_file)

     if ( mod( i_step, n_output ) == 0 ) then
        call calculate_kinetic_energy( new_KE, n_mole,n_atom,mass,vel_atom, constants%conv_kJmol_ang2ps2gmol )
        vol= volume( box(1,:),box(2,:),box(3,:) )
        call print_result( traj_file, log_file, n_mole, n_atom, xyz, box, trajectory_step, E_elec,E_elec_nopol,E_bh, E_3body, E_bond, E_angle, E_dihedral, potential,new_KE,vol )
     endif

     ! checkpoint velocities
     Select Case( checkpoint_velocity )
     Case("yes")
        if ( mod( i_step, n_step_velocity ) == 0 ) then    
           call print_velocities_checkpoint( trajectory_step, vel_file, tot_n_mole, n_atom, vel_atom )
        endif
     End Select

  end do


end program main_ms_evb
