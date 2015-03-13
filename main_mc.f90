!**************************************************************************************
! This is the head source file for the Monte Carlo program.
!
! This code was written by Jesse McDaniel and Kuang Yu, 2010-2013, while working 
! for Professor JR Schmidt.
!
! The manual for this code, examples of how to run it, and parameters for our
! SAPT-based force fields can be downloaded from the group website
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
  use eq_drude
  use total_energy_forces
  use bonded_interactions
  use sampling
  use mc_routines
  use insertion_bias_routines
  implicit none
  !************************** I/O files *************************************************
  character(MAX_FN) :: ifile_conf, ifile_fconf, ifile_ffpmt, ifile_top, ifile_simpmt, ofile_traj, ofile_log
  !
  !**************************    SIMULATION DATA    *************************************************
  character(MAX_ANAME), dimension(MAX_N_MOLE,MAX_N_ATOM) :: alist   ! this array stores the atom names for every molecule in the simulation
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
  real*8, dimension(MAX_N_MOLE,3) :: r_com,vel_com,ang_vel
  !
  ! for a flexible molecule simulation, vel_atom stores the atomic velocities
  real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: vel_atom
  !
  ! As we treat all molecules as rigid bodies, the array quaternion stores the quaternions of each molecule which gives the orientation of the rigid body
  real*8, dimension(MAX_N_MOLE,4) :: quaternion
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
  ! energy components stores the energy decomposition of the configuration, in order:  exchange, electrostatic, induction, dhf, dispersion
  real*8,dimension(5) :: energy_components     
  !
  ! E_elec_nopol stores the total coulombic energy of the system not including polarization ( also, if present fixed crystal-crystal reciprocal space contributions are subtracted out and 
  ! therefore are not included).  E_elec stores the total coulombic energy of the system including Drude oscillator contributions, E_bh stores the total Lennard-Jones or Buckingham energy
  ! of the system, and E_3body stores the explicit 3body energy (not induction) of the system
  real*8 :: E_elec_nopol, E_elec, E_bh,potential, E_3body
  real*8 :: E_bond, E_angle, E_dihedral
  !*********************************************************************************************

  character(3) :: energy_decomposition_store
  integer :: i_step, i_mole, i_atom,accept_trans,accept_rot
  real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: force_atoms_junk   ! this is for passing to scf drude at bottom of program during energy decomposition
  real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: xyz_drude_in
  TYPE(DFTI_DESCRIPTOR), POINTER :: dfti_desc,dfti_desc_inv
  real*8 :: vol,new_KE,inertia_local, E_elec_junk, E_bh_junk
  real*8 :: i_acc=0.0,i_try=0.0,i_acc_p=0.0,i_try_p=0.0,i_try_v=0.0,i_acc_v=0.0,acc_pct
  real*8 :: i_acc_old=0.0,i_try_old=0.0 , i_acc_v_old=0.0,i_try_v_old=0.0
  integer :: traj_file=98, log_file=97, sample_vel=1, i, j, status,is_overlap=0,count=0,iteration,drude_p,n_files,index
  real*8,parameter:: conv_fac=100.       ! converts kJ/mol to A^2/ps^2*g/mol


  call sort_input_files( n_files, ifile_conf, ifile_fconf, ifile_ffpmt, ifile_top, ifile_simpmt, ofile_traj, ofile_log )

  open( traj_file, file=ofile_traj, status='new' )
  open( log_file,  file=ofile_log, status='new' )

  call read_simulation_parameters( ifile_simpmt )

  call initialize_simulation(box,alist,n_mole,tot_n_mole,n_atom,n_atom_drude,xyz,quaternion,mass,r_com,chg,chg_drude,tot_chg,ifile_conf,ifile_fconf,ifile_ffpmt)



  !**** write about code comments
  write(*,*) "this code has been streamlined for simple lj force fields"
  write(*,*) "in particular, parts of pme, and pair_int energy and force"
  write(*,*) "routines have been commented to speed code up"
  write(*,*) "The following array size parameters have been decreased"
  write(*,*) "decreased dramatically to save memory:"
  write(*,*) "gfun_grid , Tang_Toennies_grid,ewald_realspace_interaction_grid_size "
  write(*,*) ""
  write(*,*) " NOTE beta spline grid size has been decreased to 10000 to speed up "
  write(*,*) "reciprocal space force calculation.  This leads to errors of ~0.1%  "
  write(*,*) "in reciprocal space force compared to grid size of 1000000          "
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


  ! if flexible simulation, read topology file, and generate intra-molecular exclusions
  Select Case(flexible_molecule_simulation)
     Case("yes")
        call read_topology_file( ifile_top )
        call generate_intramolecular_exclusions
  End Select


  call initialize_energy_force(log_file,box,n_mole,tot_n_mole,n_atom,n_atom_drude,alist,xyz,r_com,mass,tot_chg,dfti_desc,dfti_desc_inv,iteration,potential,E_elec_nopol,E_elec,E_bh,E_3body,E_bond,E_angle,E_dihedral,xyz_drude,force_atoms,energy_components)



  !****************** run code tests if requested *************************!
  if(debug > 2) then
  call run_code_tests(log_file,box,n_mole,tot_n_mole,n_atom,n_atom_drude,alist,xyz,r_com,tot_chg,dfti_desc,dfti_desc_inv)
  endif
  !************************************************************************


  ! store energy decomposition setting, as the global variable will be changed back and forth, since decomposition only needs to be done while printing output
  energy_decomposition_store = energy_decomposition

  ! decide whether we are running an mc simulation or just griding energy
  Select Case(energy_grid_run)
  Case("yes")
     Select Case(framework_simulation)
     Case(0)
        stop "must have framework for energy grid run"
     Case(1)
        call energy_grid(box,alist,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,mass,r_com,tot_chg,dfti_desc,dfti_desc_inv,log_file)
        stop
     End Select
  End Select


  call calculate_kinetic_energy( new_KE, n_mole,n_atom,mass,vel_atom, vel_com, ang_vel, conv_fac )
  vol= volume( box(1,:),box(2,:),box(3,:) )
  call print_run_param( log_file , n_mole, vol )
  call print_result( traj_file, log_file, n_mole, n_atom, alist, xyz, box, 0, energy_components,E_elec,E_elec_nopol,E_bh, E_3body, E_bond, E_angle, E_dihedral,potential,new_KE,vol )



  ! don't waste time with energy decomposition during run until printing
  energy_decomposition = "no"


  do i_step = 1, n_step

     call  mc_sample(box,iteration,alist,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,quaternion,force_atoms,vel_com,ang_vel,vel_atom,mass,r_com,energy_components,potential, E_elec,E_elec_nopol,E_bh,E_3body,E_bond,E_angle,E_dihedral,tot_chg,i_acc,i_try,i_acc_p,i_try_p,i_try_v,i_acc_v,dfti_desc,dfti_desc_inv,log_file,sample_vel)

     if ( mod( i_step, n_output ) == 0 ) then


        ! for expanded grand canonical ensemble, can only sample if terminal coupling value
        if ( (select_ensemble .ne. "egc") .or. ( partial_molecule_coupling .eq. 0 ) ) then

           ! do energy decomposition if requested
           energy_decomposition = energy_decomposition_store
           Select Case(energy_decomposition)
           Case("yes")

              Select Case(drude_simulation)
              Case(1)
                 ! here we need to call scf_drude to get E_elec_nopol, as this isn't computed if energy_decomposition is turned off
                 ! no need to overwrite force or energy, so pass junk array to scf_drude 
                 call scf_drude(E_elec_junk,E_elec_nopol,force_atoms_junk,xyz,tot_chg,r_com,box,tot_n_mole,n_mole,n_atom,n_atom_drude,iteration,xyz_drude,dfti_desc,dfti_desc_inv,log_file,xyz_drude)
              Case(0)
                 E_elec_nopol = E_elec
              End Select

              ! no need to overwrite energy here, calling this routine because energy_components isn't computed if energy_decomposition is turned off
              !call lennard_jones( E_bh_junk, force_atoms_junk, energy_components, tot_n_mole, n_mole, n_atom, xyz, box, r_com )

           end select

           call calculate_kinetic_energy( new_KE, n_mole,n_atom,mass,vel_atom, vel_com, ang_vel, conv_fac )
           vol= volume( box(1,:),box(2,:),box(3,:) )

           call print_result( traj_file, log_file, n_mole, n_atom, alist, xyz, box, i_step, energy_components,E_elec,E_elec_nopol,E_bh, E_3body, E_bond, E_angle, E_dihedral, potential,new_KE,vol )


           ! turn energy_decomposition back off
           energy_decomposition = "no"

        endif
     endif

!!$     if( mod( i_step, n_update_stepsize ) == 0 ) then
!!$
!!$        ! if single molecule moves, need to decouple translation/rotation steps for appropriate update
!!$        Select case (hybrid_md_mc)
!!$        case("yes")
!!$           if ( select_ensemble .ne. "md" ) then
!!$              call update_stepsize_trans_rot ( i_acc, i_try, i_acc_old, i_try_old )
!!$           endif
!!$        case("no")
!!$           call decouple_trans_rot_accept(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,mass,r_com,energy_components,potential,E_elec,E_elec_nopol,E_bh,E_3body,tot_chg,i_acc,i_try,dfti_desc,dfti_desc_inv,log_file,accept_trans,accept_rot)
!!$           call update_stepsize_trans_rot ( i_acc, i_try, i_acc_old, i_try_old, accept_trans, accept_rot )
!!$        End Select
!!$
!!$        ! update volume step size if npt ensemble
!!$        Select case (select_ensemble)
!!$        Case("npt")
!!$           call update_stepsize_volume ( i_acc_v, i_try_v, i_acc_v_old, i_try_v_old )
!!$        End Select
!!$
!!$     end if
  end do

  Select case (hybrid_md_mc)
  case("yes")
     write(log_file,*) 'Final time step is (ps):', delta_t1
  case("no")
     write(log_file,*) 'Final translation step:', single_mol_step
     Select case (rotate_method)
     case(1)
        write(log_file,*) 'Final rotation step:', single_mol_rotate
     case(2)
        write(log_file,*) 'Final euler angle steps (phi,psi,theta):', phi_step, psi_step,theta_step
     end select
  end select


  Select case (select_ensemble)
  Case("uvt")
     write(log_file,*) 'Acceptance ratio for insertions/deletions: ', i_acc_p / i_try_p
  Case("egc")
     write(log_file,*) 'Acceptance ratio for insertions/deletions: ', i_acc_p / i_try_p
  Case("npt")
     write(log_file,*) 'Final volume step size (in moves of ln(volume) ) : ', vol_step
     write(log_file,*) 'Acceptance ratio for volume moves: ', i_acc_v / i_try_v
  End Select

  write(log_file,*) 'Acceptance ratio for position moves is : ', i_acc/i_try

  ! for replica exchange runs, write current random seed to output file so that next simulation can continue with this seed
  Select Case(replica_input)
  Case("yes")
     write(traj_file,*) "random number seed"
     call random_seed(get = seed_replica)
     write(traj_file,*) seed_replica
  End Select

end program mc
