!**************************************************************
! This module only has one subroutine, which is pretty self-explanatory:  See below
!**************************************************************
module total_energy_forces
  use global_variables
  use MKL_DFTI
  use pairwise_interaction
  use pme_routines
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
  !
  !  E_bh:         This is all remaining pairwise additive energy ( either lennard-jones or buckingham
  !                type potentials mostly, but can possible include Feynman-Hibbs corrections ).
  !
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
  !
  !
  !              
  !
  !***********************************************************
  subroutine calculate_total_force_energy( force_atoms, potential, E_elec,  E_bh, E_bond, E_angle, E_dihedral, iteration, tot_n_mole, n_mole, n_atom, r_com, xyz, chg, box, dfti_desc,dfti_desc_inv,log_file)
    real*8, dimension(:,:,:),intent(out) :: force_atoms
    real*8,intent(out) :: potential, E_elec, E_bh, E_bond, E_angle, E_dihedral
    integer,intent(out) :: iteration
    integer,intent(in)::  tot_n_mole, n_mole
    integer,dimension(:),intent(in)::  n_atom
    real*8,dimension(:,:),intent(in)::  r_com
    real*8,dimension(:,:,:),intent(in)::  xyz
    real*8,dimension(:,:),intent(in)::  chg
    real*8,dimension(:,:),intent(in)::  box
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv
    integer,intent(in)::log_file
    real*8, dimension(:,:,:),intent(out) :: xyz_drude_final

    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: lj_force
    integer,dimension(MAX_N_MOLE,MAX_N_ATOM)::drude_atoms
    integer :: i_mole, i_atom, flag_verlet_list, flag_junk
    real*8 :: dE_bh_old, dE_bh_new

    E_bond=0d0; E_angle=0d0; E_dihedral=0d0

    ! zero forces.  It is important to zero All force elements, because for ms-evb we are constantly changing topology
    force_atoms=0d0

 

    !************************* check Verlet list, make sure it's updated ***********************
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
    !***************************************************************************************


       !*********** Electrostatics *************

!!!! no drude oscillators
          iteration = 0
          xyz_drude_final=0d0

          do i_mole=1,n_mole
             do i_atom=1,n_atom(i_mole)
                drude_atoms(i_mole,i_atom)=i_atom
             enddo
          enddo

          call pme_energy_force(E_elec, force_atoms,drude_atoms,tot_n_mole,n_mole, n_atom, xyz, chg, box, r_com, dfti_desc,dfti_desc_inv)

          ! this variable is meaningless without oscillators, but just assign a value so compiler doesn't get angry
          E_elec_nopol = E_elec


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

       force_atoms = force_atoms + lj_force


    !***************************** intra-molecular bond, angle energy and forces *****************

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

    ! total potential energy of system:  Note that E_disp_lrc is not computed here, it is stored in
    ! global_variables and only recomputed whenever volume or particle number changes, as it
    ! doesn't depend on atomic positions


    potential = E_elec + E_bh + E_disp_lrc + E_bond + E_angle + E_dihedral


!!$    write(*,*) "elec, bh , bond, angle, dihedral"
!!$    write(*,*) E_elec , E_bh , E_bond , E_angle, E_dihedral

  end subroutine calculate_total_force_energy





  !**************************************
  ! this subroutine calculates the total kinetic energy
  ! of the system
  !**************************************
  subroutine calculate_kinetic_energy( KE, n_mole,n_atom, mass, vel_atom, conv_fac )
    use global_variables
    real*8,intent(out) :: KE
    integer,intent(in) :: n_mole
    integer,dimension(:),intent(in) :: n_atom
    real*8,dimension(:,:),intent(in) :: mass
    real*8,dimension(:,:,:), intent(in) :: vel_atom
    real*8,intent(in) :: conv_fac

    integer :: i_mole, i_atom,index

    KE=0d0
    ! flexible molecules
    do i_mole=1,n_mole
       do i_atom=1,n_atom(i_mole)
          KE = KE + 0.5d0 * mass(i_mole,i_atom) * dot_product(vel_atom(i_mole,i_atom,:),vel_atom(i_mole,i_atom,:))/ conv_fac
       end do
    end do

  end subroutine calculate_kinetic_energy




end module total_energy_forces
