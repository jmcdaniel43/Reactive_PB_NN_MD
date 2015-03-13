!**************************************************************
! This module only has one subroutine, which is pretty self-explanatory:  See below
!**************************************************************
module total_energy_forces
  use global_variables
  use MKL_DFTI
  use electrostatic
  use pairwise_interaction
  use explicit_three_body_interaction
  use pme_routines
  use eq_drude
  use bonded_interactions
  use routines

  
contains

!***********************************************************
!  This subroutine collects all the possible energy and force subroutine calls, so that
!  the total energy (and if desired force) is given by a call to the subroutine
!  
!  As total forces are only currently used for hybrid MD/MC moves, we need the total energy anyway
!  with these forces, and therefore whenever forces are calculated, so is the total energy
!
!  All energies can be computed using cutoffs or Ewald (PME) sums, including long range R^-6, R^-8 etc 
!  for forces, electrostatics can be computed with either cutoff or PME, while the pairwise 
!  Lennard-Jones or Buckingham terms can only be computed with cutoff (As forces are used for MD/MC
!  moves, it is fine to use a different Hamiltonian for forces and energy)
!
!  Note that forces will not be calculated for three-body-dispersion terms, so if three-body-dispersion
!  is used in a hybrid MC/MD simulation, this will also lead to discrepancies between actual Hamiltonian
!  and sampling Hamiltonian
!  
!  Variable I/O:
!
!  force_atoms:  This array contains the total force on the atoms. Note this is only computed for atoms 
!                with index up to n_mole, as atoms with indices from n_mole to tot_n_mole are 
!                assumed fixed
!  
!  potential:    This is the total potential energy of the system, with contributions from electrostatic
!                interactions (including Drude oscillators, E_elec; LJ type terms, E_bh; any explicit 
!                three body interactions such as three body dispersion, E_3body;
!                and any LJ long range correction energy, E_disp_lrc.
!
!  E_elec:       This is the total electrostatic energy ( including any Drude oscillator energies ! )
!                for solute-solute and solute-framework ( framework meaning fixed atoms)
!                interactions.  Any reciprocal space framework framework-framework interactions
!                have been subtracted out
!  
!  E_elec_nopol: This is the same as E_elec, but does not include Drude oscillator energies
!
!  E_bh:         This is all remaining pairwise additive energy ( either lennard-jones or buckingham
!                type potentials mostly, but can include penetration force field terms, or possibly
!                Feynman-Hibbs corrections ).
!
!  E_3body       This is explicit three body contributions to the energy, such as three body dispersion.
!                Note that while the induction energy is many-body as the Drude oscillators are 
!                optimized self-consistently, this isn't included in E_3body because after the 
!                oscillators are optimized, the energy is explicitly pairwise additive.
!
!  energy_components:
!                This array contains the energy decomposition explicitly included in the force field
!                The best way to understand each element is by looking at the "print_result" subroutine
!                in the rout.f90 file
!
!  iteration:    This gives the number of iterations used to converge the Drude oscillators
!
!  tot_n_mole:   This is the total number of molecules in the system, including fixed, framework atoms.
!                Note that every fixed, framework atom is treated as its own molecule
!
!  n_mole:       This is the number of solute molecules in the simulation.  Obviously for a neat fluid,
!                n_mole will equal tot_n_mole
!
!  n_atom:       This array gives the number of atoms for each molecule, not including drude oscillators
!
!  n_atom_drude: This array gives the number of atoms plus drude oscillators for each molecule
!
!  r_com:        This array contains the center of mass for each molecule in the system
!
!  xyz:          This array gives the cartesian coordinates of every atom in the system
!
!  chg:          This array gives the total charges of every atom in the system. IMPORTANT::  These
!                total charges include the contribution from any Drude oscillator attached to the atom. 
!                therefore if an atom has a Drude oscillator attached to it with a charge -q(drude),
!                this array will hold a charge of q(tot) = q(fixed) + q(drude) for that atom.
!                Also, the Drude oscillator charges are stored in this array, after all the actual
!                atoms.  So for example if molecule 5 has 3 atoms, then the charges of the 3 atoms
!                are stored in chg(5,1),chg(5,2), and chg(5,3) and the charges of the corresponding 
!                Drude oscillators are stored in chg(5,4),chg(5,5), and chg(5,6) respectively
!
!  box:          This array contains the three cartesian components (2nd index) of each of the three
!                box vectors (1st index).
!
!  dfti_desc, dfti_desc_inv:
!                These pointers are using in the MKL Fast Fourier Transform routines
!
!  log_file:     This is the file handle of one of the output files
!
!  xyz_drude_final:    
!                This is the final position of the optimized Drude oscillators
!
!  xyz_drude_initial:
!                This is (optionally), an initial guess for the positions of the Drude oscillators. 
!
!  r_com_old, xyz_old, E_elec_old, E_bh_old, energy_components_old, i_move:
!                These variables correspond to their base variables, but are for pre-moved 
!                configurations.  They are input and used if only one molecule was moved (non-hybrid
!                moves), so that energy can be updated, rather than completely recalculated. 
!                i_move is the index of the molecule that was moved. IMPORTANT::  If these arguments
!                are present, force will not be calculated, as it is not needed for single molecule
!                moves
!              
!
!***********************************************************
  subroutine calculate_total_force_energy( force_atoms, potential, E_elec,E_elec_nopol, E_bh, E_3body, E_bond, E_angle, E_dihedral,energy_components, iteration, tot_n_mole, n_mole, n_atom, n_atom_drude, r_com, xyz, chg, box, dfti_desc,dfti_desc_inv,log_file,xyz_drude_final, xyz_drude_initial, r_com_old, xyz_old, E_elec_old,E_bh_old, energy_components_old, i_move)
    real*8, dimension(:,:,:),intent(out) :: force_atoms
    real*8,intent(out) :: potential, E_elec, E_elec_nopol, E_bh, E_3body, E_bond, E_angle, E_dihedral
    real*8,dimension(:),intent(out) :: energy_components
    integer,intent(out) :: iteration
    integer,intent(in)::  tot_n_mole, n_mole
    integer,dimension(:),intent(in)::  n_atom, n_atom_drude
    real*8,dimension(:,:),intent(in)::  r_com
    real*8,dimension(:,:,:),intent(in)::  xyz
    real*8,dimension(:,:),intent(in)::  chg
    real*8,dimension(:,:),intent(in)::  box
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv
    integer,intent(in)::log_file
    real*8, dimension(:,:,:),intent(out) :: xyz_drude_final
    real*8, dimension(:,:,:),intent(in),optional :: xyz_drude_initial
    real*8,dimension(:,:),intent(in),optional::  r_com_old
    real*8,dimension(:,:,:),intent(in),optional::  xyz_old
    real*8,intent(in),optional :: E_elec_old, E_bh_old
    real*8,dimension(:),intent(in),optional :: energy_components_old
    integer,intent(in),optional :: i_move

    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: lj_force
    integer,dimension(MAX_N_MOLE,MAX_N_ATOM)::drude_atoms
    integer :: i_mole, i_atom, flag_verlet_list, flag_junk
    real*8 :: dE_bh_old, dE_bh_new
    real*8,dimension(size(energy_components)) :: dE_components_old, dE_components_new
    character(3) :: calculate_force

    E_bond=0d0; E_angle=0d0; E_dihedral=0d0
    energy_components=0d0

    ! zero forces.  It is important to zero All force elements, because for ms-evb we are constantly changing topology
    force_atoms=0d0

    !  Decide whether to calculate forces.  Forces are needed for hybrid MD/MC type moves, but are not needed for single molecule moves.
    !  The presence (or absence) of certain optional arguments tells us what kind of moves are being done

    if ( present(r_com_old) ) then
       ! make sure all other corresponding optional arguments are present
       if ( present(xyz_old) .and. present(E_elec_old) .and. present(E_bh_old) .and. present(energy_components_old) .and. present(i_move) ) then
          calculate_force="no"
       else
          stop "use of optional arguments is inconsistent in subroutine calculate_total_force_energy"
       endif
    else
       calculate_force="yes"
    endif


    !************************* check Verlet list, make sure it's updated ***********************
    Select Case(use_verlet_list)
    Case("yes")
       call update_verlet_displacements( n_mole, n_atom, xyz, box, flag_verlet_list )     
       ! if flag_verlet_list=0, verlet list is fine, if it's 1, we need to update verlet list, if 2, we need to reallocate
       Select Case(flag_verlet_list)
       Case(1)
          call construct_verlet_list(tot_n_mole, n_mole, n_atom, xyz, box )  
          ! the "1" input to update_verlet_displacements signals to initialize the displacement array
          call update_verlet_displacements( n_mole, n_atom, xyz, box, flag_junk, 1 )
       Case(2)
          call allocate_verlet_list( tot_n_mole, n_atom, chg, box )
          call construct_verlet_list(tot_n_mole, n_mole, n_atom, xyz, box )
          call update_verlet_displacements( n_mole, n_atom, xyz, box, flag_junk, 1 )
       End Select
    End Select
    !***************************************************************************************


    Select Case(calculate_force)
    Case("yes")
       !**************************** Hybrid MD/MC, therefore calculate forces ********************!

       !*********** Electrostatics *************

!!!if there are drude oscillators, use scf_drude to optimize positions and compute energy
!!!!!!!!!get energy of new config
       if( drude_simulation .eq. 1 ) then
          if ( present(xyz_drude_initial) ) then
             call scf_drude(E_elec,E_elec_nopol,force_atoms,xyz,chg,r_com,box,tot_n_mole,n_mole,n_atom,n_atom_drude,iteration,xyz_drude_final,dfti_desc,dfti_desc_inv,log_file,xyz_drude_initial)
          else
             call scf_drude(E_elec,E_elec_nopol,force_atoms,xyz,chg,r_com,box,tot_n_mole,n_mole,n_atom,n_atom_drude,iteration,xyz_drude_final,dfti_desc,dfti_desc_inv,log_file)
          endif
       else
!!!! no drude oscillators
          iteration = 0
          xyz_drude_final=0d0

          do i_mole=1,n_mole
             do i_atom=1,n_atom(i_mole)
                drude_atoms(i_mole,i_atom)=i_atom
             enddo
          enddo

          call Electrostatic_energy_force( E_elec, force_atoms, drude_atoms,tot_n_mole, n_mole, n_atom, xyz, chg, box, r_com,dfti_desc,dfti_desc_inv)

          ! this variable is meaningless without oscillators, but just assign a value so compiler doesn't get angry
          E_elec_nopol = E_elec

       endif


       !*** Lennard Jones/ Buckingham terms ****

       !****************************timing**************************************!
       if(debug .eq. 1) then
          call date_and_time(date,time)
          write(*,*) "lj force and energy started at", time
       endif
       !***********************************************************************!


       call lennard_jones( E_bh, lj_force, tot_n_mole, n_mole, n_atom, xyz, box, r_com)

       !****************************timing**************************************!
       if(debug .eq. 1) then
          call date_and_time(date,time)
          write(*,*) "lj force and energy finished at", time
       endif
       !***********************************************************************!

       Select Case(pme_disp)
       Case("yes")
          ! recalculate lj energy with pme if requested
          call tot_disp_usegrid( E_bh,energy_components, tot_n_mole, n_mole, n_atom, xyz, box, r_com )
       End Select

       force_atoms = force_atoms + lj_force


    Case("no")
       !**************************** single MC moves, therefore don't calculate forces ********************!

       !*********** Electrostatics *************

       !  Drude oscillators shouldn't be used with single molecule moves, as this is very slow
       !  Therefore, the code should probably never enter this if statement
       !  ( A warning will be issued if Drude oscillators are used here)
       if( drude_simulation .eq. 1 ) then
          if ( present(xyz_drude_initial) ) then
             call scf_drude(E_elec,E_elec_nopol,force_atoms,xyz,chg,r_com,box,tot_n_mole,n_mole,n_atom,n_atom_drude,iteration,xyz_drude_final,dfti_desc,dfti_desc_inv,log_file,xyz_drude_initial)
          else
             call scf_drude(E_elec,E_elec_nopol,force_atoms,xyz,chg,r_com,box,tot_n_mole,n_mole,n_atom,n_atom_drude,iteration,xyz_drude_final,dfti_desc,dfti_desc_inv,log_file)
          endif
       else
          !  no drude oscillators, here xyz_old, etc are used as arguments, so that fast energy update can be used if electrostatic cutoff is used
          iteration = 0
          xyz_drude_final=0d0

          call Electrostatic_energy_force( E_elec, force_atoms, drude_atoms,tot_n_mole, n_mole, n_atom, xyz, chg, box, r_com,dfti_desc,dfti_desc_inv)
          E_elec_nopol = E_elec
       endif


       !*** Lennard Jones/ Buckingham terms ****

       ! code will stop if pme dispersion is requested with non-hybrid MD/MC moves, so don't consider that here

       call update_lennard_jones_ins ( dE_bh_old, dE_components_old, tot_n_mole, n_atom, xyz_old, i_move, box, r_com_old)
       call update_lennard_jones_ins ( dE_bh_new, dE_components_new, tot_n_mole, n_atom, xyz, i_move, box, r_com)

       energy_components = energy_components_old - dE_components_old + dE_components_new
       E_bh = E_bh_old - dE_bh_old + dE_bh_new


    End Select

    !******************************three-body-dispersion**************************************
    Select Case(three_body_dispersion)
    Case("yes")

       !****************************timing**************************************!
       if(debug .eq. 1) then
          call date_and_time(date,time)
          write(*,*) "link list started at", time
       endif
       !***********************************************************************!

       ! construct linklist every time
       call construct_linklist(box,tot_n_mole,n_mole,r_com)

       !****************************timing**************************************!
       if(debug .eq. 1) then
          call date_and_time(date,time)
          write(*,*) "link list finished at", time
       endif
       !***********************************************************************!


       ! calculate energy
       call calculate_three_body_interaction(E_3body, xyz, r_com, n_atom, box)


       !****************************timing**************************************!
       if(debug .eq. 1) then
          call date_and_time(date,time)
          write(*,*) "three body finished at", time
       endif
       !***********************************************************************!

    Case("no")
       E_3body = 0d0
    End Select
    !*****************************************************************************************



    !***************************** intra-molecular bond, angle energy and forces *****************
    Select Case(flexible_molecule_simulation)
    Case("yes")
       !****************************timing**************************************!
       if(debug .eq. 1) then
          call date_and_time(date,time)
          write(*,*) "intra_molecular force and energy started at", time
       endif
       !***********************************************************************!
       call intra_molecular_energy_force( E_bond, E_angle, E_dihedral, force_atoms, n_mole, n_atom, xyz )
       !****************************timing**************************************!
       if(debug .eq. 1) then
          call date_and_time(date,time)
          write(*,*) "intra_molecular force and energy finished at", time
       endif
       !***********************************************************************!
    End Select

    ! total potential energy of system:  Note that E_disp_lrc is not computed here, it is stored in
    ! global_variables and only recomputed whenever volume or particle number changes, as it
    ! doesn't depend on atomic positions


    potential = E_elec + E_bh + E_3body + E_disp_lrc + E_bond + E_angle + E_dihedral

  end subroutine calculate_total_force_energy





  !**************************************
  ! this subroutine calculates the total kinetic energy
  ! of the system, either given by the atomic velocites for a flexible
  ! simulation, or by the center-of-mass velocities and angular velocities
  ! for rigid bodies
  !**************************************
  subroutine calculate_kinetic_energy( KE, n_mole,n_atom, mass, vel_atom, vel_com, ang_vel, conv_fac )
    use global_variables
    real*8,intent(out) :: KE
    integer,intent(in) :: n_mole
    integer,dimension(:),intent(in) :: n_atom
    real*8,dimension(:,:),intent(in) :: mass
    real*8,dimension(:,:,:), intent(in) :: vel_atom
    real*8,dimension(:,:), intent(in) :: vel_com, ang_vel
    real*8,intent(in) :: conv_fac

    integer :: i_mole, i_atom,index

    KE=0d0

    select case(flexible_molecule_simulation)
       Case("yes")
          ! flexible molecules
          do i_mole=1,n_mole
             do i_atom=1,n_atom(i_mole)
                KE = KE + 0.5d0 * mass(i_mole,i_atom) * dot_product(vel_atom(i_mole,i_atom,:),vel_atom(i_mole,i_atom,:))/ conv_fac
             end do
          end do

       Case("no")
          ! rigid bodies
           do i_mole=1,n_mole
              index = molecule_index(i_mole)
              KE = KE + 0.5d0 * molecule_mass(index) * dot_product(vel_com(i_mole,:),vel_com(i_mole,:))/conv_fac
              call rotational_kinetic_energy(KE,i_mole,ang_vel, conv_fac)
           enddo

    end Select


  end subroutine calculate_kinetic_energy


  !******************************************************
  ! this subroutine generates the rotational kinetic energy
  ! for linear and non-linear molecules.  Here we have to separate
  ! the cases, as ang_velocity is in body coordinate system for non-linear
  ! molecules, but is not for linear molecules
  !
  ! units: inertia (g/mol*Ang^2), ang_vel(ps^-1)
  ! conv_fac converts kJ/mol to A^2/ps^2*g/mol
  !****************************************************
  subroutine rotational_kinetic_energy(KE,i_mole,ang_vel, conv_fac )
    use global_variables
    integer,intent(in) :: i_mole
    real*8,dimension(:,:),intent(in) :: ang_vel
    real*8, intent(in) :: conv_fac
    real*8,intent(inout) :: KE

    integer  :: i,index
    real*8   :: inertia_local
    real*8,parameter :: small=1d-3

    index = molecule_index(i_mole)

    ! here ang_velocity is not in principle moment basis for linear molecules, so careful with
    ! kinetic energy
    Select Case(molecule_shape(index))
    Case(0)
       ! non-linear molecule, ang_velocity is in principle moment basis
       do i=1,3
          KE = KE + .5*molecule_inertia(index,i)* ang_vel(i_mole,i)**2 /conv_fac
       enddo
    Case(1)
       ! linear molecule, angular velocity is perpendicular to axis, but in global coordinates
       do i=1,3
          if ( molecule_inertia(index,i) > small ) then
             inertia_local = molecule_inertia(index,i)
             exit
          endif
       enddo
       KE = KE + .5* inertia_local * dot_product(ang_vel(i_mole,:),ang_vel(i_mole,:))/conv_fac
    End Select


  end subroutine rotational_kinetic_energy



end module total_energy_forces
