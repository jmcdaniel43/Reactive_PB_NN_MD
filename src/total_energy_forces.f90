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
  !***********************************************************
  subroutine calculate_total_force_energy(system_data, molecule_data, atom_data, verlet_list_data, PME_data )
    Type(system_data_type),intent(inout)                :: system_data
    Type(molecule_data_type),dimension(:),intent(inout) :: molecule_data
    Type(atom_data_type),intent(inout)                  :: atom_data
    Type(verlet_list_data_type),intent(inout)           :: verlet_list_data
    Type(PME_data_type), intent(inout)                  :: PME_data

    integer :: flag_verlet_list, flag_junk


    ! zero forces and all energy components.  It is important to zero All force elements, because for ms-evb we are constantly changing topology
    atom_data%force=0d0
    system_data%potential_energy=0d0
    system_data%E_elec = 0d0
    system_data%E_vdw  = 0d0
    system_data%E_bond = 0d0
    system_data%E_angle= 0d0
    system_data%E_dihedral=0d0


    !************************* check Verlet list, make sure it's updated ***********************
       call update_verlet_displacements( system_data%total_atoms, atom_data%xyz, verlet_list_data, system_data%box, system_data%xyz_to_box_transform, flag_verlet_list )     
       ! if flag_verlet_list=0, verlet list is fine, if it's 1, we need to update verlet list
       Select Case(flag_verlet_list)
       Case(1)
          call construct_verlet_list( verlet_list_data, atom_data, system_data%total_atoms, system_data%box, system_data%xyz_to_box_transform  )
          ! the "1" input to update_verlet_displacements signals to initialize the displacement array
          call update_verlet_displacements( system_data%total_atoms, atom_data%xyz, verlet_list_data , system_data%box, system_data%xyz_to_box_transform , flag_junk, 1 )
       End Select
    !***************************************************************************************


    !*********** Electrostatics *************
    call pme_energy_force( system_data, molecule_data, atom_data, verlet_list_data, PME_data )


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
