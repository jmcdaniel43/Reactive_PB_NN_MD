!********************************************************************
! This module contains all of the subroutines that control the generation and acceptance
! of Monte Carlo moves.  The subroutine mc_sample is thus the "workhorse" of the program,
! as it is called whenever any Monte Carlo move is attempted
!********************************************************************

module sampling
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
       !**************  MD with flexible intra-molecular force field
       ! atomic velocities are now sampled in main_mc.f90
       !if ( initialize == 0 ) then
       !   call sample_atomic_velocities(n_mole, n_atom, mass, vel_atom )
       !   initialize = 1 
       !endif
       call md_integrate_atomic(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,force_atoms,vel_atom,mass,r_com,potential,E_elec,E_elec_nopol,E_bh,E_3body,E_bond, E_angle,E_dihedral,chg,dfti_desc,dfti_desc_inv,log_file)


    Case("nvt")
       call sample_atomic_velocities(n_mole, n_atom, mass, vel_atom )
       call hybridmd_move(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,force_atoms,vel_atom,mass,r_com,potential,E_elec,E_elec_nopol,E_bh,E_3body,chg,dfti_desc,dfti_desc_inv,log_file)
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
    ! displacement of drude oscillators
    real*8,dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: delta_drude
    real*8,dimension(3)                       :: dr
    real*8   :: delta_t


    !***** timestep ******
    delta_t=delta_t1


    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "md integration step started at", time
    endif
    !***********************************************************************!

    !********* first store displacements of drude oscillators
    if( drude_simulation .eq. 1 ) then
       do i_mole = 1, n_mole
          if ( n_atom_drude(i_mole) > n_atom(i_mole) ) then
             do i_drude = n_atom(i_mole)+1 , n_atom_drude(i_mole)
                ! use map to figure out which atom this oscillator is on
                i_atom = drude_atom_map(i_mole,i_drude)
                delta_drude(i_mole,i_atom,:) = xyz_drude(i_mole,i_drude,:) - xyz(i_mole,i_atom,:)
             enddo
          endif
       enddo
    endif

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


    ! update positions of drude oscillators based on atomic displacements
    if( drude_simulation .eq. 1 ) then
       do i_mole = 1, n_mole
          if ( n_atom_drude(i_mole) > n_atom(i_mole) ) then
             do i_drude = n_atom(i_mole)+1 , n_atom_drude(i_mole)
                ! use map to figure out which atom this oscillator is on
                i_atom = drude_atom_map(i_mole,i_drude)
                xyz_drude(i_mole,i_drude,:)= xyz(i_mole,i_atom,:) + delta_drude(i_mole,i_atom,:)
             enddo
          endif
       enddo
    endif


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





  !******************************************************************
  ! this is a hybrid MD/MC module for flexible molecules, for use
  ! with the ms-evb hamiltonian.  See comments in the (commented out) rigid-body
  ! hybrid md/mc subroutine for more details
  !******************************************************************
  subroutine hybridmd_move(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,force_atoms,vel_atom,mass,r_com,potential,E_elec,E_elec_nopol,E_bh,E_3body,chg,dfti_desc,dfti_desc_inv,log_file)
    use global_variables
    use MKL_DFTI
    use total_energy_forces
    integer,intent(in)::n_mole,tot_n_mole,log_file
    integer,dimension(:),intent(inout)::n_atom,n_atom_drude
    integer,intent(out)::iteration
    real*8,dimension(:,:,:),intent(inout):: vel_atom
    real*8,dimension(:,:,:),intent(inout)::xyz,xyz_drude,force_atoms
    real*8,dimension(:,:),intent(inout)::mass,chg
    real*8,dimension(:,:),intent(in)::box
    real*8,dimension(:,:),intent(inout)::r_com
    real*8,intent(inout)::potential,E_elec,E_elec_nopol,E_bh,E_3body
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv

    real*8,parameter:: conv_fac=100.       ! converts kJ/mol to A^2/ps^2*g/mol
    real*8 :: E_bond, E_angle, E_dihedral, p, potential_old, old_KE, new_KE, old_Tot_E, new_Tot_E
    !*************** stored data structures
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: xyz_store, force_atoms_store
    real*8, dimension(MAX_N_MOLE,3) :: r_com_store
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM) :: chg_store, mass_store
    integer, dimension(MAX_N_MOLE) :: n_atom_store
    !**** ms-evb stored variables
    integer, dimension(:,:), allocatable :: atom_index_store 
    integer, dimension(:), allocatable  :: molecule_index_store
    integer, dimension(:), allocatable :: hydronium_molecule_index_store

    !i_try=i_try+1.

    ! for ms-evb, we need to store the global arrays that are changed with a proton hop, in case we reject the move
    Select Case(ms_evb_simulation)
    Case("yes")
       allocate( atom_index_store(MAX_N_MOLE,MAX_N_ATOM), molecule_index_store(MAX_N_MOLE), hydronium_molecule_index_store( n_proton_max ) )
       atom_index_store = atom_index
       molecule_index_store = molecule_index
       hydronium_molecule_index_store = hydronium_molecule_index
    End Select
    ! store data structures
    ! note, don't need to store velocities, as these are resampled every time step
    xyz_store = xyz
    force_atoms_store = force_atoms
    r_com_store = r_com
    chg_store = chg
    mass_store = mass
    n_atom_store = n_atom
    potential_old = potential

    call calculate_kinetic_energy( old_KE, n_mole, n_atom,mass, vel_atom,  conv_fac )

    !*********************************  MD integration
    call md_integrate_atomic(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,force_atoms,vel_atom,mass,r_com,potential,E_elec,E_elec_nopol,E_bh,E_3body,E_bond,E_angle,E_dihedral,chg,dfti_desc,dfti_desc_inv,log_file)

    call calculate_kinetic_energy( new_KE, n_mole, n_atom,mass, vel_atom,  conv_fac )

    old_Tot_E=old_KE + potential_old
    new_Tot_E=new_KE + potential

    call random_number(p)

    if ( ( new_Tot_E < old_Tot_E ) .or. ( p < exp( -(new_Tot_E - old_Tot_E)/8.314*1000/temp ) )  ) then
       ! move is accepted! 
       !i_acc = i_acc + 1.

    else
       ! move was rejected, restore data structures
       xyz = xyz_store
       force_atoms = force_atoms_store
       r_com = r_com_store
       chg = chg_store
       mass = mass_store
       n_atom = n_atom_store
       potential = potential_old

       Select Case(ms_evb_simulation)
       Case("yes")
          atom_index = atom_index_store
          molecule_index = molecule_index_store
          hydronium_molecule_index = hydronium_molecule_index_store
       End Select
    endif


    Select Case(ms_evb_simulation)
    Case("yes")
       deallocate( atom_index_store , molecule_index_store , hydronium_molecule_index_store )
    End Select

  end subroutine hybridmd_move



  !*************************************************************************
  ! this subroutine controls the volume move part of an npt ensemble simulation
  ! it calls an internal subroutine to make the volume move and change molecule positions
  ! and then it does energy evaluations, decides whether to accept move, and updates variables
  !*************************************************************************
  subroutine npt_volume_move(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,force_atoms,mass,r_com,potential,E_elec,E_elec_nopol,E_bh,E_3body,chg,i_acc_v,i_try_v,dfti_desc,dfti_desc_inv,log_file)
    use global_variables
    use pairwise_interaction
    use MKL_DFTI
    use pme_routines
    use total_energy_forces
    integer,intent(in)::tot_n_mole,n_mole,log_file
    integer,intent(out)::iteration
    integer,dimension(:),intent(in)::n_atom,n_atom_drude
    real*8,intent(inout)::i_acc_v,i_try_v
    real*8,dimension(:,:,:),intent(inout)::xyz,xyz_drude,force_atoms
    real*8,dimension(:,:),intent(in)::mass,chg
    real*8,dimension(:,:),intent(inout)::box,r_com
    real*8,intent(inout)::potential,E_elec,E_elec_nopol,E_bh,E_3body
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv
    real*8,dimension(MAX_N_MOLE,3)::new_r_com
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: new_xyz,new_xyz_drude,xyz_drude_in,force_atoms_try,lj_force
    real*8,dimension(3,3)::new_box,kk, xyz_to_box_transform_store
    real*8,dimension(3)::a,b,c,ka,kb,kc
    real*8,dimension(pme_grid,pme_grid,pme_grid)::CB_store
    real*8::new_nrg,arg,old_vol,new_vol
    real*8 :: E_elec_try,E_elec_nopol_try, E_bh_try,E_3body_try,potential_try,p,accept_crit, E_disp_lrc_store
    integer :: i_mole, i_atom, i_drude, i, j, is_overlap
    ! these will not be used for rigid molecule simulation
    real*8 :: E_bond, E_angle, E_dihedral

    i_try_v=i_try_v+1.

    a(:) = box(1,:);b(:) = box(2,:);c(:) = box(3,:)
    old_vol = volume( a, b, c )
    call gen_vol_move(box,new_box,n_mole,n_atom,xyz,new_xyz,r_com,new_r_com,mass,vol_step,old_vol,new_vol)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! check to see that molecules are not too close together
    is_overlap= check_move( new_xyz, tot_n_mole,n_mole, n_atom, new_box, new_r_com )
    if(is_overlap .eq. 1) then
!!$       write(log_file,*) "*******************************************************"
!!$       write(log_file,*) " atoms too close, move rejected before energy calculation"
!!$       write(log_file,*) "*******************************************************"
       goto 500
    endif

    !**************** box size changed.  Need to update global data such as CB matrix, box transformation matrix, E_disp_lrc, etc for energy calculations
    ! make sure to save old values in case move is rejected
    Select Case(electrostatic_type)
    Case("pme")
       a(:) = new_box(1,:);b(:) = new_box(2,:);c(:) = new_box(3,:)
       call crossproduct( a, b, kc ); kc = kc /new_vol 
       call crossproduct( b, c, ka ); ka = ka /new_vol
       call crossproduct( c, a, kb ); kb = kb /new_vol
       kk(1,:)=ka(:);kk(2,:)=kb(:);kk(3,:)=kc(:)

!!!!!!!!!!!!store old CB array in case volume move is not accepted, and calculate new CB array
       CB_store=CB
       call CB_array(CB,alpha_sqrt,new_vol,pme_grid,kk,spline_order)

    end select

    xyz_to_box_transform_store = xyz_to_box_transform
    ! new transformation matrix to box coordinates
    call initialize_non_orth_transform ( new_box )

    !***************** make sure cutoff distances are fine for new box size
    call check_cutoffs_box( new_box )

    E_disp_lrc_store = E_disp_lrc
    ! need to update long range correction since volume has changed
    stop " below subroutine update_disp_lrc_noshift has been removed"
    !call update_disp_lrc_noshift(n_mole, n_atom, alist,new_box)

    !**************************************************************************************************************************



!!!if there are drude oscillators, use scf_drude to optimize positions and compute energy
    if( drude_simulation .eq. 1 ) then
!!!!!!!!!!!create xyz_drude_in array for initial positions of oscillators from new atom positions, and old optimized relative positions
       do i_mole = 1, n_mole
          if ( n_atom_drude(i_mole) > n_atom(i_mole) ) then
             do i_drude = n_atom(i_mole)+1 , n_atom_drude(i_mole)
                ! use map to figure out which atom this oscillator is on
                i_atom = drude_atom_map(i_mole,i_drude)

                xyz_drude_in(i_mole,i_drude,:)= new_xyz(i_mole,i_atom,:) + xyz_drude(i_mole,i_drude,:) - xyz(i_mole,i_atom,:)
             enddo
          endif
       enddo
    endif



    !**************************get energy of new config*************
    call calculate_total_force_energy( force_atoms_try, potential_try, E_elec_try, E_elec_nopol_try, E_bh_try, E_3body_try, E_bond,E_angle,E_dihedral,iteration, tot_n_mole, n_mole, n_atom, n_atom_drude, new_r_com, new_xyz, chg, new_box, dfti_desc,dfti_desc_inv,log_file, new_xyz_drude)
    !***************************************************************



!!!!!!!!!! if drude oscillators didn't converge, reject move
    if(iteration .gt. iteration_limit) then
       goto 400
    endif



    arg=-(1000./(temp*8.314))*((potential_try-potential)+pressure*p_conv*(new_vol-old_vol)-dble(n_mole+1)*log(new_vol/old_vol)/(1000./(temp*8.314)))
    accept_crit=min(1.,exp(arg))
    call random_number(p)
    if ( p < accept_crit) then
       !  move accepted! Fantastic! Now update data ...
       force_atoms = force_atoms_try 
       i_acc_v = i_acc_v + 1.
       xyz = new_xyz
       xyz_drude = new_xyz_drude
       box=new_box
       r_com=new_r_com
       E_elec = E_elec_try
       E_elec_nopol = E_elec_nopol_try
       E_bh = E_bh_try
       E_3body = E_3body_try
       potential = potential_try
    else
       ! move rejected :( , return old variables
400    continue
       Select Case(electrostatic_type)
       Case("pme")
          CB=CB_store
       End Select
       xyz_to_box_transform = xyz_to_box_transform_store
       E_disp_lrc = E_disp_lrc_store

    endif

500 continue

  end subroutine npt_volume_move



  !***********************************************************************
  ! This subroutine generates a volume move by varying ln(V)
  ! Currently, this is specific to a cubic box
  !***********************************************************************
  subroutine gen_vol_move(old_box,new_box,n_mole,n_atom,old_xyz,new_xyz,old_cofm,new_cofm,mass,vol_step,oldvol,newvol)
    integer,intent(in)::n_mole
    integer,dimension(:),intent(in)::n_atom
    real*8,dimension(:,:,:),intent(in)::old_xyz
    real*8,dimension(:,:),intent(in)::mass,old_box,old_cofm
    real*8,intent(in)::vol_step,oldvol
    real*8,intent(out)::newvol
    real*8,dimension(:,:,:),intent(out)::new_xyz
    real*8,dimension(:,:),intent(out)::new_box,new_cofm
    real*8,parameter::small=1D-4
    real*8,dimension(3)::mol_axis
    integer::i,j,k,count,flag
    integer, save :: check_box=0
    real*8::randnum,lnvol,projec

    Select Case( check_box )
    Case(0)
       ! make sure box is cubic
       flag=0
       if ( ( abs( old_box(1,1) - old_box(2,2) ) > small ) .or.  ( abs( old_box(1,1) - old_box(3,3) ) > small ) ) then
          flag = 1
       endif
       if ( (abs( old_box(1,2) ) > small ) .or. (abs( old_box(1,3) ) > small ) .or. (abs( old_box(2,1) ) > small ) .or. (abs( old_box(2,3) ) > small ) .or. (abs( old_box(3,1) ) > small ) .or. (abs( old_box(3,2) ) > small ) ) then
          flag=1
       endif

       if ( flag == 1 ) then
          stop "npt volume moves only implemented for a cubic box"
       endif
       check_box=1
    End Select

    ! change volume by steps in ln(volume)
    call random_number(randnum)


    lnvol=log(oldvol)+(randnum-.5)*vol_step
    newvol=exp(lnvol)

    new_box=0.0
    do i=1,3
       new_box(i,i)=newvol**(1./3.)
    enddo

!!!!!!!!!!!!!!!!now generate new coordinates for all atoms in all molecules relative to the new center of mass

    do i=1,n_mole
       do j=1,3
          new_cofm(i,j)=old_cofm(i,j)*new_box(j,j)/old_box(j,j)
       enddo
    enddo

    do i=1,n_mole
       do j=1,n_atom(i)
          new_xyz(i,j,:) = new_cofm(i,:) + ( old_xyz(i,j,:) - old_cofm(i,:) )
       enddo
    enddo

  end subroutine gen_vol_move



end module sampling
