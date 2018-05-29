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
   real*8, dimension(:,:), allocatable, target :: xyz
   real*8, dimension(:,:), allocatable, target :: velocity
   real*8, dimension(:,:), allocatable, target :: force
   real*8, dimension(:), allocatable, target   :: mass
   real*8, dimension(:), allocatable, target   :: charge
   integer, dimension(:), allocatable, target  :: atom_type_index
   character(MAX_ANAME),dimension(:),allocatable, target :: aname


  !***** Local variables
  integer :: i_step
 

  !********** initialize constants in global variables
  call initialize_constants( file_io_data , verlet_list_data , PME_data )

  call sort_input_files( file_io_data )
  call check_restart_trajectory( file_io_data )

  call read_simulation_parameters( file_io_data%ifile_simpmt_file_h, file_io_data%ifile_simpmt, system_data, integrator_data, verlet_list_data, PME_data )
  call initialize_simulation( system_data, molecule_data, atom_data, file_io_data, verlet_list_data, PME_data, xyz, velocity, force, mass, charge, atom_type_index, aname )


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
  if ( ( abs(system_data%box(1,2)) > 10D-6 ) .or. ( abs(system_data%box(1,3)) > 10D-6 ) .or. ( abs(system_data%box(2,1)) > 10D-6 ) .or.  ( abs(system_data%box(2,3)) > 10D-6 ) .or. ( abs(system_data%box(3,1)) > 10D-6 ) .or. ( abs(system_data%box(3,2)) > 10D-6 ) ) then
     write(*,*)  ""
     write(*,*)  "code has been modified to assume orthorhombic box"
     write(*,*)  "to change this, fix commented sections in energy calculation subroutines"
     write(*,*)  ""
     stop
  end if
  !***********

  !**********
  ! verlet list size
  write(*,*) " verlet list size parameter verlet_safe is set to ", verlet_list_data%safe_verlet
  write(*,*) " too big a value could slow simulation down "
  write(*,*) ""
  !************


  call initialize_energy_force( system_data, molecule_data, atom_data, verlet_list_data, PME_data, file_io_data, integrator_data )


  Select Case( restart_trajectory )
  Case("no")
     call sample_atomic_velocities(system_data%n_mole, system_data%total_atoms, system_data%temperature, molecule_data, atom_data )
  End Select

  call calculate_kinetic_energy( system_data%kinetic_energy, system_data%total_atoms, atom_data%mass, atom_data%velocity, constants%conv_kJmol_ang2ps2gmol )

  ! only print here if not a restart
  Select Case( restart_trajectory )
  Case("no")
     i_step = 0
     call print_simulation_info( file_io_data%ofile_log_file_h , system_data , integrator_data, verlet_list_data, PME_data )
     call print_step( file_io_data%ofile_traj_file_h, file_io_data%ofile_log_file_h , i_step, system_data, integrator_data, molecule_data, atom_data )
  End Select


  do i_step = 1, integrator_data%n_step

     ! this is global variable which may be used for printing in ms-evb
     trajectory_step = n_old_trajectory + i_step

     call  mc_sample( system_data , molecule_data , atom_data, integrator_data, verlet_list_data, PME_data, file_io_data )

     if ( mod( i_step, integrator_data%n_output ) == 0 ) then
        call calculate_kinetic_energy( system_data%kinetic_energy, system_data%total_atoms, atom_data%mass, atom_data%velocity, constants%conv_kJmol_ang2ps2gmol )
        call print_step( file_io_data%ofile_traj_file_h, file_io_data%ofile_log_file_h , trajectory_step, system_data, integrator_data, molecule_data, atom_data )
     endif

     ! checkpoint velocities
     Select Case( checkpoint_velocity )
     Case("yes")
        if ( mod( i_step, n_step_velocity ) == 0 ) then    
           call print_velocities_checkpoint( file_io_data%ifile_velocity_file_h, i_step , integrator_data%delta_t, system_data , molecule_data , atom_data  )
        endif
     End Select

  end do


end program main_ms_evb
