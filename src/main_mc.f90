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
program mc
  use routines
  use global_variables
  use simulation_parameters
  use pairwise_interaction
  use MKL_DFTI
  use total_energy_forces
  use ms_evb_fit
  use bonded_interactions
  use sampling
  use mc_routines
  implicit none
  !************************** I/O files *************************************************
  ! note the ofile_hop file name is stored in global_variables
  character(MAX_FN) :: ifile_conf, ifile_traj, ifile_fconf, ifile_ffpmt, ifile_top, ifile_simpmt, ofile_traj, ofile_log
  !**************************    SIMULATION DATA    *************************************************
  !
  ! n_mole and tot_n_mole only differ in a simulation in which there is a fixed crystal lattice (GCMC adsorption simulation for example). In the
  ! case they differ, n_mole is the total number of solute molecules in the system, and ( tot_n_mole - n_mole ) is the number of fixed crystal lattice
  ! atoms in the simulation.  Note we treat each of these lattice atoms as an individual molecule.
  ! in all other cases, n_mole = tot_n_mole = total number of molecules in the system
  integer :: n_mole, tot_n_mole   
  !
  ! n_atom stores the number of atoms on each molecule.  n_atom_drude stores the total number of atoms + Drude oscillators on each molecule
  integer, dimension(MAX_N_MOLE) :: n_atom, n_atom_drude
  !
  ! r_com stores the center of mass position of each molecule. vel_com and ang_vel store the COM velocity, and rigid-body angular velocity of each molecule respectively
  real*8, dimension(MAX_N_MOLE,3) :: r_com
  !
  ! for a flexible molecule simulation, vel_atom stores the atomic velocities
  real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: vel_atom
  !
  !
  ! xyz stores the coordinates of every atom in the simulation, xyz_drude stores the coordinates of drude oscillators, and force_atoms stores the forces on each atom in the simulation
  real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: xyz,xyz_drude,force_atoms
  !
  ! chg stores the STATIC charges on all atoms, chg_drude stores the drude oscillator charges, and tot_chg stores the STATIC + POLARIZATION (due to drude oscillators) charges on each 
  ! atom.  Note we will normally just use the tot_chg array
  ! mass stores the masses of all atoms in the simulation
  real*8, dimension(MAX_N_MOLE,MAX_N_ATOM) :: chg, chg_drude,tot_chg,mass
  !
  ! box stores the box vectors of the system
  real*8, dimension(3,3) :: box
  !  
  !
  ! E_elec_nopol stores the total coulombic energy of the system not including polarization ( also, if present fixed crystal-crystal reciprocal space contributions are subtracted out and 
  ! therefore are not included).  E_elec stores the total coulombic energy of the system including Drude oscillator contributions, E_bh stores the total Lennard-Jones or Buckingham energy
  ! of the system, and E_3body stores the explicit 3body energy (not induction) of the system
  real*8 :: E_elec_nopol, E_elec, E_bh,potential, E_3body
  real*8 :: E_bond, E_angle, E_dihedral
  !*********************************************************************************************

  integer :: i_step
  real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: xyz_drude_in
  TYPE(DFTI_DESCRIPTOR), POINTER :: dfti_desc,dfti_desc_inv
  real*8 :: vol,new_KE
  integer :: sample_vel=1, i, j, status,is_overlap=0,count=0,iteration,n_files,index
  !************** file_handles**********
  ! note the file handle for the evb_hop_file is stored in global variables
  integer :: conf_file=99, traj_file=98, log_file=97,vel_file=96
  !*************************************
  real*8,parameter:: conv_fac=100.       ! converts kJ/mol to A^2/ps^2*g/mol



  call sort_input_files( n_files, ifile_conf, ifile_traj, ifile_fconf, ifile_ffpmt, ifile_top, ifile_simpmt, ofile_traj, ofile_log, ofile_hop )
  call check_restart_trajectory( ofile_traj , ofile_log, traj_file, log_file, vel_file )

  call read_simulation_parameters( ifile_simpmt )
  call initialize_simulation(box,n_mole,tot_n_mole,n_atom,n_atom_drude,xyz,mass,vel_atom, r_com,chg,chg_drude,tot_chg,ifile_conf,ifile_top, ifile_fconf,ifile_ffpmt, conf_file, traj_file, vel_file, ofile_traj)


  !*** open output files if this is not a continuation run.  If this is a continuation, these files will
  ! already be open, and prepared to be appended to, in the check_restart_trajectory subroutine
  Select Case( restart_trajectory )
  Case("no")
     open( traj_file, file=ofile_traj, status='new' )
     open( log_file,  file=ofile_log, status='new' )
     Select Case( checkpoint_velocity )
     Case("yes")
        ! status may be new, or append to old
        open( vel_file, file=velocity_file, status ='new')
     End Select
  end Select


  !**** write about code comments
  write(*,*) "this code has been streamlined for simple lj force fields"
  write(*,*) "in particular, parts of pme, and pair_int energy and force"
  write(*,*) "routines have been commented to speed code up"
  write(*,*) ""
  write(*,*) " NOTE beta spline grid size has been decreased to 100000 to speed up "
  write(*,*) "reciprocal space force calculation.         "
  write(*,*) ""
  ! assume orthorhombic box
  if ( ( abs(box(1,2)) > 10D-6 ) .or. ( abs(box(1,3)) > 10D-6 ) .or. ( abs(box(2,1)) > 10D-6 ) .or.  ( abs(box(2,3)) > 10D-6 ) .or. ( abs(box(3,1)) > 10D-6 ) .or. ( abs(box(3,2)) > 10D-6 ) ) then
     write(*,*)  ""
     write(*,*)  "code has been modified to assume orthorhombic box"
     write(*,*)  " to change this, fix commented sections in subroutines"
     write(*,*)  " pme_real_space_use_verlet, lennard_jones_use_verlet, and"
     write(*,*)  " construct_verlet_list_grid "
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

!!$  write(*,*) "****************************************************"
!!$  write(*,*) " WARNING:  evb_reactive_pair_distance temporarily set to 0.5 !!! "
!!$  write(*,*) "****************************************************"


  call initialize_energy_force(log_file,box,n_mole,tot_n_mole,n_atom,n_atom_drude,xyz,r_com,mass,tot_chg,dfti_desc,dfti_desc_inv,iteration,potential,E_elec_nopol,E_elec,E_bh,E_3body,E_bond,E_angle,E_dihedral,xyz_drude,force_atoms)


  !*** test
!!$  write(*,*) "forces"
!!$  do i=1,n_mole
!!$     do j=1,n_atom(i)
!!$        write(*,*) i, j, force_atoms(i,j,:)
!!$     enddo
!!$  enddo
!!$  stop
  !****

  !************ if parameter fitting, call that module
  Select Case(ms_evb_parameter_fit)
  Case("yes")

     call fit_ms_evb_parameters( ifile_traj , log_file, n_mole, n_atom, tot_chg , dfti_desc , dfti_desc_inv , box ) 
     stop
  end Select
  !*****************************

  Select Case( restart_trajectory )
  Case("no")
     call sample_atomic_velocities(n_mole, n_atom, mass, vel_atom )
  End Select

  call calculate_kinetic_energy( new_KE, n_mole,n_atom,mass,vel_atom,  conv_fac )
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
        call calculate_kinetic_energy( new_KE, n_mole,n_atom,mass,vel_atom,  conv_fac )
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


end program mc
