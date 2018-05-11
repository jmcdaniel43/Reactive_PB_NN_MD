!********************************************************************
! This module contains all of the subroutines that control the generation and acceptance
! of Monte Carlo moves.  The subroutine mc_sample is thus the "workhorse" of the program,
! as it is called whenever any Monte Carlo move is attempted
!********************************************************************

module md_integration
  use routines
  implicit none

contains


  !*************************************************************************
  !  this subroutine controls what ensemble to run, and calls appropriate subroutines
  !  currently implemented ensembles are nvt, uvt, egc (expanded grand canonical)
  !*************************************************************************

  subroutine mc_sample(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,force_atoms,vel_atom,mass,r_com,potential,E_elec,E_elec_nopol,E_bh,E_3body,E_bond, E_angle,E_dihedral,chg,dfti_desc,dfti_desc_inv,log_file)
    use global_variables
    use MKL_DFTI
    use conj_grad
    integer,intent(inout)::n_mole,tot_n_mole
    integer,intent(in)::log_file
    integer,dimension(:),intent(inout)::n_atom,n_atom_drude
    integer,intent(out)::iteration
    real*8,dimension(:,:,:),intent(inout)::xyz,xyz_drude,force_atoms, vel_atom
    real*8,dimension(:,:),intent(inout)::mass,chg
    real*8,dimension(:,:),intent(inout)::box
    real*8,dimension(:,:),intent(inout)::r_com
    real*8,intent(inout)::potential,E_elec,E_elec_nopol,E_bh,E_3body, E_bond, E_angle, E_dihedral
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv

    real*8:: randnum, pct_v_move
    integer:: i_mole,i_atom,i,md_flag
    integer,save :: warn_egc=0 , warn_gcmc_single=0, initialize=0


    Select Case (select_ensemble)
    Case("cg")
       !******************* conjugate gradient energy minimization ************
       call conjugate_gradient_minimize(box,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,mass,potential,chg,dfti_desc,dfti_desc_inv,log_file)
    Case("stp")
       !******************* steepest descent energy minimization ***************
       call steepest_descent_minimize(box,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,mass,potential,chg,dfti_desc,dfti_desc_inv,log_file)
    Case("md")
       !*********************** this is nve md ******************************!
       call md_integrate_atomic(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,force_atoms,vel_atom,mass,r_com,potential,E_elec,E_elec_nopol,E_bh,E_3body,E_bond, E_angle,E_dihedral,chg,dfti_desc,dfti_desc_inv,log_file)
    case default
       stop " select_ensemble option not recognized in gcmc_step"
    end select

  end subroutine mc_sample



  !*************************************************************************
  ! samples velocities for all atoms from Maxwell-Boltzmann
  ! velocities are in Angstrom/ps . angular velocities are in ps^-1
  !*************************************************************************
  subroutine sample_atomic_velocities(n_mole, n_atom, mass, vel_atom )
    use global_variables
    integer,intent(in)::n_mole
    integer,dimension(:),intent(in) :: n_atom
    real*8,dimension(:,:),intent(in) :: mass
    real*8,dimension(:,:,:),intent(out):: vel_atom

    integer :: i_mole, i_atom, n_tot, atom_id
    real*8,dimension(2)  :: vel
    real*8,parameter :: small=1D-3

    real*8,parameter:: conv_fac=100.       ! converts kJ/mol to A^2/ps^2*g/mol
    real*8 :: sum_KE, norm

    n_tot=0
    sum_KE=0d0

    ! zero velocity, in case we're freezing an atom
    vel_atom=0d0

    !**************** first pull velocities from maxwell-boltzmann distribution
    do i_mole=1,n_mole
       do i_atom=1, n_atom(i_mole)
          ! make sure mass is non-zero
          if ( mass(i_mole, i_atom) < small ) then
             write(*,*) "trying to assign velocity for atom ", i_atom
             write(*,*) "of molecule ", i_mole, " but mass is zero!"
             stop
          end if

          atom_id = atom_index(i_mole,i_atom)
          ! get velocity if atomtype isn't frozen
          if ( atype_freeze(atom_id) /= 1 ) then
             ! gaussian random numbers come in 2
             call max_boltz(vel,mass(i_mole,i_atom),temp)
             vel_atom(i_mole,i_atom,1)=vel(1)
             vel_atom(i_mole,i_atom,2)=vel(2)
             call max_boltz(vel,mass(i_mole,i_atom),temp)
             vel_atom(i_mole,i_atom,3)=vel(1)
          end if

       enddo
    enddo


    ! get rid of excess center of mass momentum of system
    call subtract_center_of_mass_momentum(n_mole, n_atom, mass, vel_atom )

    ! now rescale velocities to desired temperature
    do i_mole=1,n_mole
       do i_atom=1, n_atom(i_mole)
          atom_id = atom_index(i_mole,i_atom)
          ! add KE if atom isn't frozen
          if ( atype_freeze(atom_id) /= 1 ) then
             n_tot = n_tot + 1
             sum_KE = sum_KE + 0.5d0 * mass(i_mole,i_atom) * dot_product(vel_atom(i_mole,i_atom,:),vel_atom(i_mole,i_atom,:)) / conv_fac
          end if
       enddo
    enddo

    norm = 1.5d0 * 0.008314d0 * temp * dble(n_tot) / sum_KE
    vel_atom = vel_atom * sqrt(norm)


  end subroutine sample_atomic_velocities




  !*******************************************************
  ! this subroutine calculates the center of mass momentum of the system,
  ! and subtracts the total net per atom contribution from each atom's momentum,
  ! so that the net COM momentum is zero
  !*******************************************************
  subroutine subtract_center_of_mass_momentum(n_mole, n_atom, mass, vel_atom )
    use global_variables
    integer,intent(in)::n_mole
    integer,dimension(:),intent(in) :: n_atom
    real*8,dimension(:,:),intent(in) :: mass
    real*8,dimension(:,:,:),intent(inout):: vel_atom

    integer :: i_mole, i_atom, n_tot, atom_id
    real*8,dimension(3) :: rho_system, rho_excess

    n_tot=0
    rho_system=0d0

    !**************** calculate total COM momentum
    do i_mole=1,n_mole
       do i_atom=1, n_atom(i_mole)

          atom_id = atom_index(i_mole,i_atom)
          ! add momentum if atomtype isn't frozen
          if ( atype_freeze(atom_id) /= 1 ) then
             n_tot = n_tot + 1
             rho_system(:) = rho_system(:) + mass(i_mole,i_atom) * vel_atom(i_mole,i_atom,:)
          end if

       enddo
    enddo

    !*************  now subtract excess momentum from each atom, so that center of mass momentum is zero
    rho_excess(:) = rho_system(:) / dble(n_tot)

    do i_mole=1,n_mole
       do i_atom=1, n_atom(i_mole)
          atom_id = atom_index(i_mole,i_atom)
          ! change velocity if atom isn't frozen
          if ( atype_freeze(atom_id) /= 1 ) then
             vel_atom(i_mole,i_atom,:) = vel_atom(i_mole,i_atom,:) - rho_excess(:) / mass(i_mole,i_atom)
          end if
       enddo
    enddo

  end subroutine subtract_center_of_mass_momentum









  !************************************************************************
  ! this is the MD engine for flexible molecular simulations.  
  ! currently, this uses the velocity verlet algorithm to integrate Newton's equations
  !  velocity units are A/ps
  !
  ! NOTE we have already checked for non-zero atomic masses in the sample_atomic_velocities
  ! subroutine, and so here we don't worry about it
  !
  ! 3/2/2015 i/o attributes changed by JGM
  ! a call to the multi-state empirical valence bond module
  ! ( ms_evb_calculate_total_force_energy ) may change data structures
  ! chg, n_atom, mass, if a proton transfer occurs
  ! therefore these data structures were given the intent(inout) attribute
  !
  ! additionally, if the topology of the system changes due to a proton transfer in the ms-evb module,
  ! the forces and velocity data structures need to be updated consistently.  This is taken
  ! care of in the ms-evb module, and that's why the velocity array is passed
  !************************************************************************
  subroutine md_integrate_atomic(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,force_atoms,vel_atom,mass,r_com,potential,E_elec,E_elec_nopol,E_bh,E_3body,E_bond,E_angle,E_dihedral,chg,dfti_desc,dfti_desc_inv,log_file)
    use global_variables
    use MKL_DFTI
    use total_energy_forces
    use ms_evb
    integer,intent(in)::n_mole,tot_n_mole,log_file
    integer,dimension(:),intent(inout)::n_atom,n_atom_drude
    integer,intent(out)::iteration
    real*8,dimension(:,:,:),intent(inout):: vel_atom
    real*8,dimension(:,:,:),intent(inout)::xyz,xyz_drude,force_atoms
    real*8,dimension(:,:),intent(inout)::mass,chg
    real*8,dimension(:,:),intent(in)::box
    real*8,dimension(:,:),intent(inout)::r_com
    real*8,intent(inout)::potential,E_elec,E_elec_nopol,E_bh,E_3body, E_bond, E_angle, E_dihedral
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv

    integer :: i_mole, i_atom, i_drude, atom_id
    real*8,parameter:: conv_fac=100.       ! converts kJ/mol to A^2/ps^2*g/mol
    real*8,dimension(3)                       :: dr
    real*8   :: delta_t


    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "md integration step started at", time
    endif
    !***********************************************************************!


    !******************* Velocity Verlet Integrator


    ! first calculate velocities at delta_t / 2
    do i_mole = 1, n_mole
       do i_atom = 1, n_atom(i_mole)
          atom_id = atom_index(i_mole,i_atom)
          ! update position, velocity if atom isn't frozen
          if ( atype_freeze(atom_id) /= 1 ) then
             ! first calculate velocities at delta_t / 2
             ! here force is in kJ/mol*Ang^-1
             vel_atom(i_mole,i_atom,:) = vel_atom(i_mole,i_atom,:) + delta_t / 2d0 / mass(i_mole,i_atom) * force_atoms(i_mole,i_atom,:) * conv_fac
             ! now calculate new atomic coordinates at delta_t
             dr(:) = delta_t * vel_atom(i_mole,i_atom,:)
             xyz(i_mole,i_atom,:) = xyz(i_mole,i_atom,:) + dr(:)
          end if
       end do

       ! after updating atomic coordinates, calculate new center of mass of molecules
       r_com(i_mole,:) = pos_com( xyz, i_mole, n_atom, mass )
       ! translate molecules back into the box
       call shift_move(n_atom(i_mole),xyz(i_mole,:,:),r_com(i_mole,:),box)

    end do



    !**********************get total forces and energies*****************************!
    Select Case(ms_evb_simulation)
    Case("yes")
       ! note velocity array is passed, because if the topology changes with a proton transfer, this needs to be updated
       call ms_evb_calculate_total_force_energy( force_atoms, potential, tot_n_mole, n_mole, n_atom, xyz, r_com, chg, mass, box, dfti_desc,dfti_desc_inv,log_file, vel_atom )
    Case("no")
       call calculate_total_force_energy( force_atoms, potential,E_elec, E_elec_nopol, E_bh, E_3body, E_bond, E_angle, E_dihedral,iteration, tot_n_mole, n_mole, n_atom, n_atom_drude, r_com, xyz, chg, box, dfti_desc,dfti_desc_inv,log_file, xyz_drude)
    End Select
    !********************************************************************************!


    ! now get final velocities
    do i_mole = 1, n_mole
       do i_atom = 1, n_atom(i_mole)
          atom_id = atom_index(i_mole,i_atom)
          ! change velocity if atom isn't frozen
          if ( atype_freeze(atom_id) /= 1 ) then
             ! here force is in kJ/mol*Ang^-1
             vel_atom(i_mole,i_atom,:) = vel_atom(i_mole,i_atom,:) + delta_t /2d0 / mass(i_mole,i_atom) * force_atoms(i_mole,i_atom,:) * conv_fac

             if ( ( abs(force_atoms(i_mole,i_atom,1)) > 10d4 ) .or. ( abs(force_atoms(i_mole,i_atom,2)) > 10d4 ) .or. ( abs(force_atoms(i_mole,i_atom,3)) > 10d4 ) ) then
                write(*,*) "forces too big"
                write(*,*) i_mole, i_atom, atype_name(atom_index(i_mole,i_atom))
                write(*,*) force_atoms(i_mole,i_atom,:)
                stop
             end if
          end if
       end do
    enddo

    ! finally remove center of mass momentum.  This should be numerical noise, so shouldn't effect energy conservation
    call subtract_center_of_mass_momentum(n_mole, n_atom, mass, vel_atom )


    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "md integration step finished at", time
    endif
    !***********************************************************************!

  end subroutine md_integrate_atomic





end module md_integration
