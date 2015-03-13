!********************************************************************
! This module contains all of the subroutines that control the generation and acceptance
! of Monte Carlo moves.  The subroutine mc_sample is thus the "workhorse" of the program,
! as it is called whenever any Monte Carlo move is attempted
!********************************************************************

module sampling
  use routines
  use rigid_body
  implicit none

contains



  !*************************************************************************
  !  this subroutine controls what ensemble to run, and calls appropriate subroutines
  !  currently implemented ensembles are nvt, uvt, egc (expanded grand canonical)
  !*************************************************************************

  subroutine mc_sample(box,iteration,alist,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,quaternion,force_atoms,vel_com,ang_vel,vel_atom,mass,r_com,energy_components,potential,E_elec,E_elec_nopol,E_bh,E_3body,E_bond, E_angle,E_dihedral,chg,i_acc,i_try,i_acc_p,i_try_p,i_try_v,i_acc_v,dfti_desc,dfti_desc_inv,log_file,sample_vel)
    use global_variables
    use pairwise_interaction
    use MKL_DFTI
    use eq_drude
    use pme_routines
    use conj_grad
    character(*), dimension(:,:),intent(inout) :: alist
    integer,intent(inout)::n_mole,tot_n_mole
    integer,intent(in)::log_file
    integer,dimension(:),intent(inout)::n_atom,n_atom_drude
    integer,intent(out)::iteration
    integer,intent(inout)::sample_vel
    real*8,intent(inout)::i_acc,i_try,i_acc_p,i_try_p,i_try_v,i_acc_v
    real*8,dimension(:,:),intent(inout):: vel_com,ang_vel
    real*8,dimension(:,:,:),intent(inout)::xyz,xyz_drude,force_atoms, vel_atom
    real*8,dimension(:,:),intent(inout):: quaternion
    real*8,dimension(:,:),intent(inout)::mass,chg
    real*8,dimension(:,:),intent(inout)::box
    real*8,dimension(:,:),intent(inout)::r_com
    real*8,dimension(:),intent(inout):: energy_components
    real*8,intent(inout)::potential,E_elec,E_elec_nopol,E_bh,E_3body, E_bond, E_angle, E_dihedral
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv

    real*8:: randnum, pct_v_move
    integer:: i_mole,i_atom,i,md_flag
    integer,save :: warn_egc=0 , warn_gcmc_single=0, initialize=0


    Select Case (select_ensemble)
    Case("cg")
       !******************* conjugate gradient energy minimization ************
       call conjugate_gradient_minimize(box,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,mass,potential,chg,dfti_desc,dfti_desc_inv,log_file)
    Case("md")
       !*********************** this is nve md ******************************!
       Select Case(flexible_molecule_simulation)
       Case("yes")
          !**************  MD with flexible intra-molecular force field
          if ( initialize == 0 ) then
             call sample_atomic_velocities(n_mole, n_atom, temp, mass, vel_atom )
             initialize = 1 
          endif
          call md_integrate_atomic(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,force_atoms,vel_atom,mass,r_com,energy_components,potential,E_elec,E_elec_nopol,E_bh,E_3body,E_bond, E_angle,E_dihedral,chg,dfti_desc,dfti_desc_inv,log_file)

       Case("no")
          !************** rigid body MD
          if ( initialize == 0 ) then
             call hybrid_md_get_velocities(1, n_mole, temp, vel_com , ang_vel)
             initialize = 1 
          endif

          ! md_flag setting is meaningless, it just has to be present in call to hybridmd_move
          md_flag=1
          call hybridmd_move(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,quaternion,force_atoms,vel_com,ang_vel,mass,r_com,energy_components,potential,E_elec,E_elec_nopol,E_bh,E_3body,chg,i_acc,i_try,dfti_desc,dfti_desc_inv,log_file,sample_vel,md_flag)

       End Select
    Case("nvt")
       !**************************************nvt*******************************************************!
       Select Case (hybrid_md_mc)
       Case("yes")
          !****************************hybrid md/mc moves *************************************!
          call hybrid_md_get_velocities(sample_vel, n_mole, temp, vel_com , ang_vel)

          call hybridmd_move(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,quaternion,force_atoms,vel_com,ang_vel,mass,r_com,energy_components,potential,E_elec,E_elec_nopol,E_bh,E_3body,chg,i_acc,i_try,dfti_desc,dfti_desc_inv,log_file,sample_vel)

       Case("no")
          !*************************single particle moves******************************************!
          call mc_make_single_mole_move(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,mass,r_com,energy_components,potential,E_elec,E_elec_nopol,E_bh,E_3body,chg,i_acc,i_try,dfti_desc,dfti_desc_inv,log_file)

       case default
          stop "hybrid_md_mc option not recognized in mc_sample"
       end select

    Case("npt")
       !**************************************npt*******************************************************!
       Select Case( hybrid_md_mc )
       Case("yes")
          ! for hybrid md, a cycle is one hybrid step
          pct_v_move = 1d0/ (dble(cycles_per_vol_step) + 1d0 )
       Case("no")
          ! for single molecule moves, n_mole moves is a cycle
          pct_v_move = 1d0/ (dble(cycles_per_vol_step * n_mole ) + 1d0 )
       End Select

       call random_number(randnum)


       if(randnum .ge. pct_v_move) then
!!!!!!!!!!! moving particles
          Select Case (hybrid_md_mc)
          Case("yes")
             !****************************hybrid md/mc moves *************************************!
             call hybrid_md_get_velocities(sample_vel, n_mole, temp, vel_com , ang_vel )

             call hybridmd_move(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,quaternion,force_atoms,vel_com,ang_vel,mass,r_com,energy_components,potential,E_elec,E_elec_nopol,E_bh,E_3body,chg,i_acc,i_try,dfti_desc,dfti_desc_inv,log_file,sample_vel)

          Case("no")
             !****************************single molecule translations ****************************! 
             ! only allow single molecule moves if there are no drude oscillators
             Select Case( drude_simulation )
             Case (1)
                stop "cannot use single molecule moves in npt with drude oscillators"
             End Select
             call mc_make_single_mole_move(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,mass,r_com,energy_components,potential,E_elec,E_elec_nopol,E_bh,E_3body,chg,i_acc,i_try,dfti_desc,dfti_desc_inv,log_file)
          End Select
       else
!!!! making volume move
          call npt_volume_move(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,alist,xyz,xyz_drude,force_atoms,mass,r_com,energy_components,potential,E_elec,E_elec_nopol,E_bh,E_3body,chg,i_acc_v,i_try_v,dfti_desc,dfti_desc_inv,log_file)

       endif

    Case("uvt")
       !**************************************uvt******************************************************!
       call random_number(randnum)

       if(randnum .ge. gcmc_ins_rm) then
!!!!!!!!!!! moving particles
          Select Case (hybrid_md_mc)
          Case("yes")
             !****************************hybrid md/mc moves *************************************!
             call hybrid_md_get_velocities(sample_vel, n_mole, temp, vel_com , ang_vel )
             call hybridmd_move(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,quaternion,force_atoms,vel_com,ang_vel,mass,r_com,energy_components,potential,E_elec,E_elec_nopol,E_bh,E_3body,chg,i_acc,i_try,dfti_desc,dfti_desc_inv,log_file,sample_vel)

          Case("no")
             Select Case(translation_method)
             Case(1)
                !****************************single molecule translations ********************************! 
                ! only allow single molecule moves if there are no drude oscillators
                Select Case( drude_simulation )
                Case (1)
                   stop "cannot use single molecule moves in gcmc with drude oscillators"
                End Select
                call mc_make_single_mole_move(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,mass,r_com,energy_components,potential,E_elec,E_elec_nopol,E_bh,E_3body,chg,i_acc,i_try,dfti_desc,dfti_desc_inv,log_file)
             Case(2)    
                !**************************** lattice hops/hybrid combination*********************!
                ! warn if no hybrid moves are used with lattice hops, as lattice hops could get stuck on own
                Select Case( warn_gcmc_single )
                Case (0)
                   if (pct_singlemv_hybrid > .999) then
                      write(*,*) "WARNING using strictly lattice hops might result in molecules being stuck.  A better idea is to add some percentage of hybrid-md moves"
                   endif
                   warn_gcmc_single=1
                End Select
                ! decide whether to do hybrid moves or hops
                call random_number(randnum)
                if(randnum .ge. pct_singlemv_hybrid ) then   
                   call hybrid_md_get_velocities(sample_vel, n_mole, temp, vel_com , ang_vel )
                   call hybridmd_move(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,quaternion,force_atoms,vel_com,ang_vel,mass,r_com,energy_components,potential,E_elec,E_elec_nopol,E_bh,E_3body,chg,i_acc,i_try,dfti_desc,dfti_desc_inv,log_file,sample_vel)
                else
                   call mc_make_single_mole_move(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,mass,r_com,energy_components,potential,E_elec,E_elec_nopol,E_bh,E_3body,chg,i_acc,i_try,dfti_desc,dfti_desc_inv,log_file)
                endif

             End Select
          case default
             stop "hybrid_md_mc option not recognized in mc_sample"
          end select

       else
!!!!!!!!!!!!!!! inserting/deleting particles
          call grand_canonical_ins_rm(box,iteration,alist,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,quaternion,force_atoms,mass,r_com,energy_components,potential,E_elec,E_elec_nopol,E_bh,chg,i_acc_p,i_try_p,dfti_desc,dfti_desc_inv,log_file)
          ! make sure we haven't inserted a molecule that exceeds the maximum number of molecules limit
          if ( tot_n_mole > MAX_N_MOLE ) then
             stop "molecule just inserted exceeds the total number of molecules allowed.  Please increase the value of MAX_N_MOLE in glob_v.f90"
          endif
       endif

    Case("egc")
       !****************expanded grand canonical******************************************************!
       ! warn against overflowing number of atom types
       if ( warn_egc == 0 ) then
          write(*,*) " "
          write(*,*) " This is an expanded grand canonical simulation, note that an additional number of atom types up to the number of atoms"
          write(*,*) " in the largest molecule will be created.  Make sure MAX_N_ATOM is large enough to handle these additional atom types"
          write(*,*) " "
          warn_egc=1
       endif

       call random_number(randnum)

       if(randnum .ge. gcmc_ins_rm) then
!!!!!!!!!!! moving particles
          Select Case (hybrid_md_mc)
          Case("yes")
             !****************************hybrid md/mc moves *************************************!
             call hybrid_md_get_velocities(sample_vel, n_mole, temp, vel_com , ang_vel )
             call hybridmd_move(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,quaternion,force_atoms,vel_com,ang_vel,mass,r_com,energy_components,potential,E_elec,E_elec_nopol,E_bh,E_3body,chg,i_acc,i_try,dfti_desc,dfti_desc_inv,log_file,sample_vel)

          Case("no")
             stop "must choose hybrid_md_mc for expanded grand canonical ensemble"
          case default
             stop "hybrid_md_mc option not recognized in mc_sample"
          end select

       else
!!!!!!!!!!!!!!! changing coupling of molecule, this could correspond to insertion/deletion
          call expanded_grand_canonical_change_coupling(box,iteration,alist,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,quaternion,force_atoms,mass,r_com,energy_components,potential,E_elec,E_elec_nopol,E_bh,chg,i_acc_p,i_try_p,dfti_desc,dfti_desc_inv,log_file)
          ! make sure we haven't inserted a molecule that exceeds the maximum number of molecules limit
          if ( tot_n_mole > MAX_N_MOLE ) then
             stop "molecule just inserted exceeds the total number of molecules allowed.  Please increase the value of MAX_N_MOLE in glob_v.f90"
          endif
       endif


    case default
       stop " select_ensemble option not recognized in gcmc_step"
    end select

  end subroutine mc_sample


  !*************************************************************************
  ! this subroutine samples velocities for hybrid md/mc if sample_vel is set to 1
  ! see hybridmd_move for comments on setting sample_vel
  !
  ! velocities are in Angstrom/ps . angular velocities are in ps^-1
  !*************************************************************************
  subroutine hybrid_md_get_velocities(sample_vel, n_mole, temp, vel_com , ang_vel )
    use global_variables
    integer,intent(in)::sample_vel,n_mole
    real*8,dimension(:,:),intent(inout):: vel_com,ang_vel
    real*8, intent(in) :: temp

    integer :: i_mole, v_flag, i , index
    real*8,dimension(2)  :: vel
    real*8,dimension(3)  :: inertia_local
    real*8  :: mass_local
    real*8,parameter :: small=1D-3

!!!!!!!!!!! gaussian random numbers come in 2, so that's the whole business of v_flag
    if(sample_vel .eq. 1) then
       do i_mole=1,n_mole
          index = molecule_index(i_mole)
          mass_local = molecule_mass(index)
          v_flag=0
          do i=1,3
             if (v_flag.eq.0) then
                call max_boltz(vel,mass_local,temp)
                vel_com(i_mole,i)=vel(1)
                v_flag=1
             elseif(v_flag.eq.1) then
                vel_com(i_mole,i)=vel(2)
                v_flag=0
             endif
          enddo
       enddo
!!!!!!!!!! now sample angular velocities
       do i_mole=1,n_mole
          index = molecule_index(i_mole)
          inertia_local(:) = molecule_inertia(index,:)
          do i=1,3
             ! if inertia element is zero, set corresponding ang_velocity to zero
             if ( inertia_local(i) > small ) then
                call max_boltz(vel,inertia_local(i),temp)
                ang_vel(i_mole,i)=vel(1)
             else
                ang_vel(i_mole,i)= 0d0
             endif
          enddo
       enddo
    endif

  end subroutine hybrid_md_get_velocities


  !*************************************************************************
  ! samples velocities for all atoms
  ! velocities are in Angstrom/ps . angular velocities are in ps^-1
  !*************************************************************************
  subroutine sample_atomic_velocities(n_mole, n_atom, temp, mass, vel_atom )
    integer,intent(in)::n_mole
    integer,dimension(:),intent(in) :: n_atom
    real*8,dimension(:,:),intent(in) :: mass
    real*8,dimension(:,:,:),intent(out):: vel_atom
    real*8, intent(in) :: temp

    integer :: i_mole, i_atom
    real*8,dimension(2)  :: vel
    real*8,parameter :: small=1D-3

       do i_mole=1,n_mole
          do i_atom=1, n_atom(i_mole)

             ! make sure mass is non-zero
             if ( mass(i_mole, i_atom) < small ) then
                write(*,*) "trying to assign velocity for atom ", i_atom
                write(*,*) "of molecule ", i_mole, " but mass is zero!"
                stop
             end if
          ! gaussian random numbers come in 2
                call max_boltz(vel,mass(i_mole,i_atom),temp)
                vel_atom(i_mole,i_atom,1)=vel(1)
                vel_atom(i_mole,i_atom,2)=vel(2)
                call max_boltz(vel,mass(i_mole,i_atom),temp)
                vel_atom(i_mole,i_atom,3)=vel(1)
          enddo
       enddo

  end subroutine sample_atomic_velocities


  !************************************************************************
  ! this is the MD engine for flexible molecular simulations.  For rigid-body MD
  ! use subroutine hybridmd_move instead
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
  !************************************************************************
  subroutine md_integrate_atomic(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,force_atoms,vel_atom,mass,r_com,energy_components,potential,E_elec,E_elec_nopol,E_bh,E_3body,E_bond,E_angle,E_dihedral,chg,dfti_desc,dfti_desc_inv,log_file)
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
    real*8,dimension(:),intent(inout)::energy_components
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv

    integer :: i_mole, i_atom, i_drude
    real*8,parameter:: conv_fac=100.       ! converts kJ/mol to A^2/ps^2*g/mol
    ! displacement of drude oscillators
    real*8,dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: delta_drude
    real*8,dimension(3)                       :: dr
    real*8   :: delta_t

    !test
    integer :: flag
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3), save :: xyz_store, vel_store, force_store
    !test


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
          ! first calculate velocities at delta_t / 2
          ! here force is in kJ/mol*Ang^-1
          vel_atom(i_mole,i_atom,:) = vel_atom(i_mole,i_atom,:) + delta_t / 2d0 / mass(i_mole,i_atom) * force_atoms(i_mole,i_atom,:) * conv_fac
          ! now calculate new atomic coordinates at delta_t
          dr(:) = delta_t * vel_atom(i_mole,i_atom,:)
          xyz(i_mole,i_atom,:) = xyz(i_mole,i_atom,:) + dr(:)
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

    ! test
    flag=0
    do i_mole=1,n_mole
       do i_atom=1,n_atom(i_mole)
          if ( ( abs( xyz(i_mole,i_atom,1)) > 200d0 ) .or. ( abs( xyz(i_mole,i_atom,2)) > 200d0 ) .or. ( abs( xyz(i_mole,i_atom,1)) > 200d0 ) ) then
             flag=1
          endif
       enddo
    enddo
    if ( flag == 1 ) then
       write(*,*) "new coordinates"
       do i_mole=1,n_mole
          do i_atom=1,n_atom(i_mole)
             write(*,'(I5,I5,3F14.6)') i_mole, i_atom, xyz(i_mole,i_atom,:)
          enddo
       enddo
       write(*,*) "old coordinates"
       do i_mole=1,n_mole
          do i_atom=1,n_atom(i_mole)
             write(*,'(I5,I5,3F14.6)') i_mole, i_atom, xyz_store(i_mole,i_atom,:)
          enddo
       enddo
       write(*,*) "old velocities"
       do i_mole=1,n_mole
          do i_atom=1,n_atom(i_mole)
             write(*,'(I5,I5,3F14.6)') i_mole, i_atom, vel_store(i_mole,i_atom,:)
          enddo
       enddo
       write(*,*) "old forces"
       do i_mole=1,n_mole
          do i_atom=1,n_atom(i_mole)
             write(*,'(I5,I5,3F14.6)') i_mole, i_atom, force_store(i_mole,i_atom,:)
          enddo
       enddo
       stop
    endif
    ! test

    !**********************get total forces and energies*****************************!
    Select Case(ms_evb_simulation)
    Case("yes")
       call ms_evb_calculate_total_force_energy( force_atoms, potential, tot_n_mole, n_mole, n_atom, xyz, r_com, chg, mass, box, dfti_desc,dfti_desc_inv,log_file )
    Case("no")
       call calculate_total_force_energy( force_atoms, potential,E_elec, E_elec_nopol, E_bh, E_3body, E_bond, E_angle, E_dihedral,energy_components, iteration, tot_n_mole, n_mole, n_atom, n_atom_drude, r_com, xyz, chg, box, dfti_desc,dfti_desc_inv,log_file, xyz_drude, xyz_drude_initial=xyz_drude)
    End Select
    !********************************************************************************!


    ! now get final velocities
    do i_mole = 1, n_mole
       do i_atom = 1, n_atom(i_mole)
          ! here force is in kJ/mol*Ang^-1
          vel_atom(i_mole,i_atom,:) = vel_atom(i_mole,i_atom,:) + delta_t /2d0 / mass(i_mole,i_atom) * force_atoms(i_mole,i_atom,:) * conv_fac

          if ( ( abs(force_atoms(i_mole,i_atom,1)) > 10d4 ) .or. ( abs(force_atoms(i_mole,i_atom,2)) > 10d4 ) .or. ( abs(force_atoms(i_mole,i_atom,3)) > 10d4 ) ) then
             write(*,*) "forces too big"
             write(*,*) force_atoms(i_mole,i_atom,:)
             stop
          end if

       end do
    enddo

    ! test
    xyz_store = xyz
    force_store = force_atoms
    vel_store = vel_atom
    ! end test


    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "md integration step finished at", time
    endif
    !***********************************************************************!

  end subroutine md_integrate_atomic



  !*************************************************************************
  !  for an overview of hybrid MC/MD, see
  !  Duane et al. Phys Lett. B, 1987, 195, 216-222
  !
  !  this subroutine obtains new atom positions using a reversible hybrid md move
  !  in principle, any reversible integration will work
  !  velocity units are A/ps
  !  acceptance criteria includes kinetic energy as it must to obey detailed balance,
  !  as the probability of making the move depends on the probability of sampling the 
  !  given velocities
  !
  !  IMPORTANT NOTE.  There can be (and are for certain situations) two different Hamiltonians used here
  !  There is the real Hamiltonian of the system, which is used for the acceptance criteria.
  !  Then there is the Hamiltonian for generating the trial (MD) Markovian moves, this can in principle be any hamiltonian
  !  We will use two different Hamiltonians specifically for the case of PME dispersion, and Feynmann Hibbs correction, and three-body-dispersion
  !  since there is probably no reason to code in forces for these energy routines.  Therefore, if using either of these
  !  energy evaluations, forces will be calculated with the normal F_lennard_jones routine
  !
  !  REVERSIBLE INTEGRATION
  !  TRANSLATION:  velocity verlet algorithm
  !
  !  ROTATION:  Linear molecules, since we can always pick the principle axis along the
  !  current angular velocity vector, Euler's equations are decoupled, and can just use
  !  velocity verlet to integrate angular velocity in the same way as center of mass velocity
  !             Non Linear molecules.  Use reversible integration of quaternions as described in
  !  Matubayasi N. and Nakahara, M. J Chem. Phys. 110, 1999, 3291
  !
  !  this subroutine is currently for simulations with drude oscillators
  !  currently, velocites are resampled after every move (not just rejected ones) as
  !  it was found by trial runs that this is necessary for detailed balance
  !
  ! input array "chg" should be total charges,including drudes
  !*************************************************************************
  subroutine hybridmd_move(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,quaternion,force_atoms,vel_com,ang_vel,mass,r_com,energy_components,potential,E_elec,E_elec_nopol,E_bh,E_3body,chg,i_acc,i_try,dfti_desc,dfti_desc_inv,log_file,sample_vel,md_flag)
    use global_variables
    use insertion_bias_routines
    use MKL_DFTI
    use total_energy_forces
    integer,intent(in)::n_mole,tot_n_mole,log_file
    integer,dimension(:),intent(in)::n_atom,n_atom_drude
    integer,intent(out)::iteration
    integer,intent(inout)::sample_vel
    real*8,intent(inout)::i_acc,i_try
    real*8,dimension(:,:),intent(inout):: vel_com,ang_vel
    real*8,dimension(:,:,:),intent(inout)::xyz,xyz_drude,force_atoms
    real*8,dimension(:,:),intent(inout)::quaternion
    real*8,dimension(:,:),intent(in)::mass,chg
    real*8,dimension(:,:),intent(in)::box
    real*8,dimension(:,:),intent(inout)::r_com
    real*8,intent(inout)::potential,E_elec,E_elec_nopol,E_bh,E_3body
    real*8,dimension(:),intent(inout)::energy_components
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv
    integer, intent(in),optional :: md_flag
    real*8,dimension(MAX_N_MOLE,3)::new_r_com,force_com,torque,vel_store,ang_store
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: new_xyz,new_xyz_drude,xyz_drude_in,force_atoms_try
    real*8, dimension(MAX_N_MOLE,4):: new_quaternion
    real*8, dimension(3):: delta_com,t_q
    real*8,dimension(size(energy_components)) :: energy_components_try
    real*8 :: E_elec_try, E_elec_nopol_try,E_bh_try,E_3body_try,potential_try,p,old_KE,new_KE,old_Tot_E,new_Tot_E,delta_t,puddle_delta_E_try
    integer :: i_mole, i_atom, i_drude, is_overlap, index,i,testi=3
    real*8,parameter:: conv_fac=100.       ! converts kJ/mol to A^2/ps^2*g/mol
    ! these will not be used for rigid molecule simulation
    real*8 :: E_bond, E_angle, E_dihedral


    i_try=i_try+1.
    delta_t=delta_t1
    ! make sure framework info gets copied
    new_r_com = r_com
    new_xyz= xyz
    new_quaternion= quaternion

    vel_store = vel_com
    ang_store = ang_vel


    force_com=0.;torque=0.;old_KE=0.
!!!!!!!!!!!!!!!!!!!! calculate center of mass forces, and total torque on molecules and perform center of mass integration
!!!!!!!!!!!!!!!!!!!!!!! generate new velocities for delta_t/2, and new positions for delta_t
    do i_mole=1,n_mole
       do i_atom=1,n_atom(i_mole)
          force_com(i_mole,:)=force_com(i_mole,:) + force_atoms(i_mole,i_atom,:)
          delta_com(:)= xyz(i_mole,i_atom,:)-r_com(i_mole,:)
          call crossproduct(delta_com,force_atoms(i_mole,i_atom,:),t_q)
          torque(i_mole,:) = torque(i_mole,:) + t_q(:)
       enddo

       index = molecule_index(i_mole)

       ! mass is grams/mole, velocity is Angstroms/ps
       old_KE=old_KE + .5*molecule_mass(index)*dot_product(vel_com(i_mole,:),vel_com(i_mole,:))/conv_fac

       call rotational_kinetic_energy(old_KE,i_mole,ang_vel, conv_fac )


       ! this is velocity at delta_t/2 as given by velocity verlet
       ! here force is in kJ/mol*Ang^-1
       vel_com(i_mole,:)=vel_com(i_mole,:) + delta_t/2./molecule_mass(index) * force_com(i_mole,:)*conv_fac
       ! this is center of mass position at delta_t as given by velocity verlet
       new_r_com(i_mole,:)= r_com(i_mole,:) + delta_t * vel_com(i_mole,:)

       ! these are pre-rotatated positions at delta_t
       do i_atom=1, n_atom(i_mole)
          new_xyz(i_mole,i_atom,:)= xyz(i_mole,i_atom,:) + delta_t * vel_com(i_mole,:)
       enddo

       call shift_move(n_atom(i_mole),new_xyz(i_mole,:,:),new_r_com(i_mole,:),box)

       ! now rotate molecules to get final postions (quaternions) at delta_t, as well as 
       ! angular velocites at delta_t / 2
       if ( n_atom(i_mole) > 1 ) then
          call reversible_rotate_q_w(i_mole,new_xyz(i_mole,:,:),n_atom(i_mole),new_r_com(i_mole,:),new_quaternion(i_mole,:), ang_vel(i_mole,:),torque(i_mole,:),delta_t,conv_fac,sample_vel)
       endif

    enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! check to see that molecules are not too close together
    is_overlap= check_move( new_xyz, tot_n_mole,n_mole, n_atom, box, new_r_com )

    if(is_overlap .eq. 1) then
!!$       write(log_file,*) "*******************************************************"
!!$       write(log_file,*) " atoms too close, move rejected before energy calculation"
!!$       write(log_file,*) "*******************************************************"

       ! while we can reject move in Monte Carlo, if this is MD, this shouldn't happen
       if ( present(md_flag) ) then
       write(*,*) "overlap of molecules"
       stop
       endif

       goto 200
    endif

    xyz_drude_in=0d0
!!!!!!!!!!!create xyz_drude_in array for initial positions of oscillators from new atom positions, and old optimized relative positions
    if( drude_simulation .eq. 1 ) then
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

    !**********************get total forces and energies*****************************!
    call calculate_total_force_energy( force_atoms_try, potential_try,E_elec_try, E_elec_nopol_try, E_bh_try, E_3body_try, E_bond, E_angle, E_dihedral, energy_components_try, iteration, tot_n_mole, n_mole, n_atom, n_atom_drude, new_r_com, new_xyz, chg, box, dfti_desc,dfti_desc_inv,log_file, new_xyz_drude, xyz_drude_initial=xyz_drude_in)
    !********************************************************************************!



!!!!!!!!!! if drude oscillators didn't converge, reject move
    if(iteration .gt. iteration_limit) then
       goto 200
    endif



!!!!!!!!!!!!! update velocities 
    force_com=0.; torque=0.; new_KE=0.
    do i_mole=1,n_mole
       do i_atom=1,n_atom(i_mole)
          force_com(i_mole,:)=force_com(i_mole,:) + force_atoms_try(i_mole,i_atom,:)
          delta_com(:)= new_xyz(i_mole,i_atom,:)-new_r_com(i_mole,:)
          call crossproduct(delta_com,force_atoms_try(i_mole,i_atom,:),t_q)
          torque(i_mole,:) = torque(i_mole,:) + t_q(:)
       enddo

       index = molecule_index(i_mole)
       ! here force is in kJ/mol*Ang^-1
       vel_com(i_mole,:) = vel_com(i_mole,:) + delta_t/2./molecule_mass(index) * force_com(i_mole,:)*conv_fac

       ! now finish integration of angular velocities, which depend on final torques
       if ( n_atom(i_mole) > 1 ) then
          call reversible_rotate_w(i_mole,new_xyz(i_mole,:,:),n_atom(i_mole),new_r_com(i_mole,:),new_quaternion(i_mole,:), ang_vel(i_mole,:),torque(i_mole,:),delta_t,conv_fac,sample_vel)
       endif

       new_KE=new_KE + .5*molecule_mass(index)*dot_product(vel_com(i_mole,:),vel_com(i_mole,:))/conv_fac
       call rotational_kinetic_energy(new_KE,i_mole,ang_vel,conv_fac )

    enddo


    old_Tot_E=old_KE + potential
    new_Tot_E=new_KE + potential_try

    ! if puddle filling
    Select Case(puddle_filling_on)
    Case("yes")
       call puddle_fill_potential(n_mole,tot_n_mole,n_atom,n_atom_drude,new_xyz,box,new_r_com,chg,puddle_delta_E_try)
       old_Tot_E = old_Tot_E + puddle_delta_E
       new_Tot_E = new_Tot_E + puddle_delta_E_try
    end Select


    call random_number(p)


!!$    write(*,*) "md_flag", md_flag

    ! md_flag will be present if this is md run, and we are not doing monte carlo
    if ( ( new_Tot_E < old_Tot_E ) .or. ( p < exp( -(new_Tot_E - old_Tot_E)/8.314*1000/temp ) ) .or. ( present(md_flag) ) ) then
       ! move is accepted! Hurray!, now update variables...
!!$       write(log_file,*) "move accepted"

       i_acc = i_acc + 1.
       xyz= new_xyz
       xyz_drude=new_xyz_drude
       force_atoms=force_atoms_try
       quaternion=new_quaternion
       r_com=new_r_com
       E_elec = E_elec_try
       E_elec_nopol = E_elec_nopol_try
       E_bh = E_bh_try
       E_3body = E_3body_try

       Select Case(puddle_filling_on)
       Case("yes")
          puddle_delta_E = puddle_delta_E_try
       End Select

       potential = potential_try
       energy_components = energy_components_try
       sample_vel=1

    else
200    continue
!!!!!!!!!!!! move was rejected, resample velocities
       sample_vel=1
       vel_com=vel_store
       ang_vel=ang_store
    endif


  end subroutine hybridmd_move




  !*****************************************
  ! this subroutine carries out single molecule position moves, in contrast to 
  ! hybridmd_move, which moves all molecules simultaneously
  ! 
  ! because only moving single molecules, if there are not drude oscillators, could take advantage
  ! of fast ewald, or fast cutoff ( looping over N particles real space , only updating
  ! p(k) in reciprocal space, for particular particle).  But with drude oscillators we lose this
  ! advantage and must use pme
  ! therefore, since we are mainly interested in using drude oscillators, these fast methods aren't
  ! implemented
  !
  ! since our single molecule moves are expensive (i.e. not using O(N) electrostatic energy routine),
  ! we do translations and rotations separately.  Therefore, the step sizes of these are coupled, and
  ! we don't know how to appropriately update the step size with just the i_acc, and i_try variables
  ! therefore, we allow an optional argument "test_trans_rot" to be passed, and if this is present
  ! the subroutine will test either a translation or rotation (will not make the move), and will output
  ! the probability of accepting either of these move types
  !
  ! warn against using this with drude oscillators as it is slow
  !******************************************
  subroutine mc_make_single_mole_move(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,mass,r_com,energy_components,potential,E_elec,E_elec_nopol,E_bh,E_3body,chg,i_acc,i_try,dfti_desc,dfti_desc_inv,log_file,test_trans_rot)
    use global_variables
    use MKL_DFTI
    use total_energy_forces
    integer,intent(in)::n_mole,tot_n_mole,log_file
    integer,dimension(:),intent(in)::n_atom,n_atom_drude
    integer,intent(out)::iteration
    integer,intent(inout),optional:: test_trans_rot
    real*8,intent(inout)::i_acc,i_try
    real*8,dimension(:,:,:),intent(inout)::xyz,xyz_drude
    real*8,dimension(:,:),intent(in)::mass,chg
    real*8,dimension(:,:),intent(in)::box
    real*8,dimension(:,:),intent(inout)::r_com
    real*8,intent(inout)::potential,E_elec,E_elec_nopol,E_bh, E_3body
    real*8,dimension(:),intent(inout)::energy_components
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv

    real*8,dimension(MAX_N_MOLE,3)::new_r_com,r_com_old
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: force_atoms_junk 
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: new_xyz,xyz_old,new_xyz_drude,xyz_drude_in
    real*8 :: E_elec_try, E_elec_nopol_try,E_bh_try,E_3body_try,E_elec_old,E_bh_old,potential_try,p,rand
    real*8,dimension(size(energy_components)) :: energy_components_try, energy_components_old
    real*8 :: weight
    integer,save :: warn_drude = 0
    integer :: i_mole, i_atom, i_drude, new_mole, is_overlap,count
    real*8,parameter:: conv_fac=100D0,small = 1D-8
    ! these will not be used for rigid molecule simulation
    real*8 :: E_bond, E_angle, E_dihedral

    Select case (drude_simulation)
    case(1)
       if ( warn_drude .eq. 0 ) then
          write(*,*) "*****WARNING***** drude simulation with single molecule moves ===SLOW!!!!"
          warn_drude = 1
       endif
    end select

    !********* if there are no solute molecules, skip this routine
    if ( n_mole == 0 ) then
       goto 200
    endif


    !************** if we are just testing translation/rotation acceptance, do that here
    if ( present(test_trans_rot) ) then
       call single_molecule_move(xyz,r_com,box,new_xyz,new_r_com,n_mole,n_atom,new_mole,test_trans_rot)

    else

       i_try=i_try+1.

       ! decide whether we are doing single molecule displacements or lattice hops
       Select Case( translation_method )
       Case(1)
          call single_molecule_move(xyz,r_com,box,new_xyz,new_r_com,n_mole,n_atom,new_mole)
       Case(2)
          ! should only use lattice_hop if this is a framework simulation
          Select Case( framework_simulation)
          Case(0)
             stop "can't use single_molecule_lattice_hop in a non framework simulation"
          End Select
          call single_molecule_lattice_hop(xyz,r_com,box,new_xyz,new_r_com,tot_n_mole, n_mole,n_atom,new_mole,weight)
          ! reject move if weight =0, meaning we are attempting hop from illegal lattice site
          if (weight < small ) then
             goto 200
          endif
       Case default
          stop "translation_method setting not recognized, must be 1 or 2"
       End Select

    endif

    is_overlap = check_distance( new_xyz, tot_n_mole, n_atom, new_mole, new_xyz(new_mole,:,:), box, mass, new_r_com )
    if ( is_overlap .eq. 1) then
       goto 200
    endif



!!!if there are drude oscillators, use scf_drude to optimize positions and compute energy
    if( drude_simulation .eq. 1 ) then
!!!!!!!!!!!create xyz_drude_in array for initial positions of oscillators from new atom positions, and old optimized relative positions
       xyz_drude_in = xyz_drude
       if ( n_atom_drude(new_mole) > n_atom(new_mole) ) then
          do i_drude = n_atom(new_mole)+1 , n_atom_drude(new_mole)
             ! use map to figure out which atom this oscillator is on
             i_atom = drude_atom_map(new_mole,i_drude)
             xyz_drude_in(new_mole,i_drude,:)= new_xyz(new_mole,i_atom,:) + xyz_drude(new_mole,i_drude,:) - xyz(new_mole,i_atom,:)
          enddo
       endif
    endif

    ! use variable_old, so optional arguments in calculate_total_force_energy are completely clear
    r_com_old=r_com
    xyz_old=xyz
    E_elec_old=E_elec
    E_bh_old=E_bh
    energy_components_old=energy_components

    call calculate_total_force_energy( force_atoms_junk,potential_try,E_elec_try,E_elec_nopol_try,E_bh_try,E_3body_try,E_bond,E_angle,E_dihedral,energy_components_try, iteration, tot_n_mole, n_mole, n_atom, n_atom_drude, new_r_com, new_xyz, chg, box, dfti_desc,dfti_desc_inv,log_file,new_xyz_drude, xyz_drude_initial=xyz_drude_in,r_com_old=r_com_old,xyz_old=xyz_old,E_elec_old=E_elec_old,E_bh_old=E_bh_old,energy_components_old=energy_components_old, i_move = new_mole )

    if(iteration .gt. iteration_limit) then
       goto 200
    endif


    call random_number(p)

    !************* if we are just testing translation/rotation acceptance, test move, but dont accept it
    if ( present(test_trans_rot) ) then
       if ( potential_try < potential .or. p < exp( -(potential_try-potential)/8.314*1000d0/temp ) ) then
          test_trans_rot = -test_trans_rot
       endif
       goto 200
    endif

    if ( potential_try < potential .or. p < exp( -(potential_try-potential)/8.314*1000d0/temp ) ) then
       i_acc = i_acc + 1.
       xyz= new_xyz
       xyz_drude=new_xyz_drude
       r_com=new_r_com
       E_elec = E_elec_try
       E_elec_nopol = E_elec_nopol_try
       E_bh = E_bh_try
       E_3body = E_3body_try
       potential = potential_try
       energy_components = energy_components_try
    else
200    continue
!!!!!!!!!!!! move was rejected
    endif

  end subroutine mc_make_single_mole_move



  !************************************************************************
  ! this subroutine generates new coordinates for a single molecule translation/rotation move
  ! this is done with random displacements and rotations
  !
  ! if optional argument test_trans_rot is present and equals "1", only a translation will be done
  ! if it is present and equals "2", only a rotation will be done
  !************************************************************************
  subroutine single_molecule_move(xyz,r_com,box,new_xyz,new_r_com,n_mole,n_atom,newmol,test_trans_rot)
    use global_variables
    integer,intent(in)::n_mole
    integer,intent(out) :: newmol
    integer,dimension(:),intent(in)::n_atom
    integer,intent(in),optional:: test_trans_rot
    real*8,dimension(:,:,:),intent(in)::xyz
    real*8,dimension(:,:,:),intent(out)::new_xyz
    real*8,dimension(:,:),intent(in)::box
    real*8,dimension(:,:),intent(in)::r_com
    real*8,dimension(:,:),intent(out)::new_r_com

    real*8, dimension(3):: delta_x,rand3,unit,mol_axis
    real*8, dimension(4):: junk_quaternion
    real*8 :: rand
    real*8  :: projec,rel_x
    integer :: i,j,i_atom, count,generate_rotation
    real*8,parameter:: small=1D-4 
    ! for euler rotation
    real*8 :: a(3,3),ainv(3,3),phi,psi,theta

!!!!!!!!!!!!!!pick a molecule to move at random
    call random_number(rand)
    newmol=ceiling(rand*n_mole)
    new_r_com = r_com
    new_xyz = xyz


!!!!!!!!!!!!!!!!!translate center of mass
    call random_number(rand3)

    if ( present(test_trans_rot) ) then
       Select Case(test_trans_rot)
       Case(1)
          new_r_com(newmol,:) = new_r_com(newmol,:) + 2d0 * single_mol_step * (rand3(:) - .5d0)
          generate_rotation=0
       Case(2)
          generate_rotation=1
       End Select
    else
       new_r_com(newmol,:) = new_r_com(newmol,:) + 2d0 * single_mol_step * (rand3(:) - .5d0)
       generate_rotation=1
    endif

    if ( generate_rotation == 1 ) then
       Select case (rotate_method)
       case(1)
          !*********** these rotations are for linear molecules only.  If non-linear molecules are present,
          !*********** and this rotation option is on, code will be stopped before this point
          call random_unit_vector(unit)
!!!!!!!!!!!!!!!!determine the molecular axis by the first atom not located at the center of mass
          do i_atom = 1,n_atom(newmol)
             delta_x(:)=xyz(newmol,i_atom,:)-r_com(newmol,:)
             if(dot_product(delta_x,delta_x) > small) then
                mol_axis = delta_x
                exit
             endif
          enddo

          mol_axis=mol_axis/sqrt(dot_product(mol_axis,mol_axis))

!!!!!!!!!!!!! rotate molecule by adding the random unit vector times scale factor gamma to mol_axis

          mol_axis=mol_axis + single_mol_rotate * unit
          mol_axis=mol_axis/sqrt(dot_product(mol_axis,mol_axis))

          do i_atom = 1,n_atom(newmol)
             delta_x(:)=xyz(newmol,i_atom,:)-r_com(newmol,:)
             rel_x = sqrt(dot_product(delta_x,delta_x))
             projec=dot_product( delta_x , mol_axis)
             if(projec > 0) then
                new_xyz(newmol,i_atom,:) = rel_x * mol_axis(:) + new_r_com(newmol,:)
             else
                new_xyz(newmol,i_atom,:) =  -rel_x * mol_axis(:) + new_r_com(newmol,:)
             endif
          enddo

       case(2)
          call random_number(rand3)

          phi = rand3(1) * phi_step
          psi = rand3(2) * psi_step
          theta = rand3(3) * theta_step

          ! here quaternions are not needed since we are not using hybrid_md moves
          call euler_transformation_matrix(a,phi,psi,theta,junk_quaternion)
          ! the inverse is the transpose
          do i=1,3
             do j=1,3
                ainv(i,j) = a(j,i)
             enddo
          enddo

          ! for detailed balance, prob 1/2 rotate with a, prob 1/2 rotate with ainv
          call random_number(rand)

          if (rand .lt. .5 ) then
             do i_atom = 1,n_atom(newmol)
                delta_x(:)=xyz(newmol,i_atom,:)-r_com(newmol,:)
                delta_x = matmul(a,delta_x)
                new_xyz(newmol,i_atom,:) = delta_x(:) + new_r_com(newmol,:)
             enddo
          else
             do i_atom = 1,n_atom(newmol)
                delta_x(:)=xyz(newmol,i_atom,:)-r_com(newmol,:)
                delta_x = matmul(ainv,delta_x)
                new_xyz(newmol,i_atom,:) = delta_x(:) + new_r_com(newmol,:)
             enddo
          endif

       case default
          stop " unrecognized rotate_method "
       end select

    endif

    call shift_move(n_atom(newmol),new_xyz(newmol,:,:),new_r_com(newmol,:),box)

  end subroutine single_molecule_move




  !****************************************************
  !  this subroutine is to individually test translation and rotation moves
  !  for single molecule moves so as to decouple the acceptance ratio for each
  !
  !****************************************************
  subroutine decouple_trans_rot_accept(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,mass,r_com,energy_components,potential,E_elec,E_elec_nopol,E_bh,E_3body,chg,i_acc,i_try,dfti_desc,dfti_desc_inv,log_file,accept_trans,accept_rot)
    use global_variables
    use MKL_DFTI
    integer,intent(in)::n_mole,tot_n_mole,log_file
    integer,dimension(:),intent(in)::n_atom,n_atom_drude
    integer,intent(out)::iteration, accept_trans, accept_rot
    real*8,intent(inout)::i_acc,i_try
    real*8,dimension(:,:,:),intent(inout)::xyz,xyz_drude
    real*8,dimension(:,:),intent(in)::mass,chg
    real*8,dimension(:,:),intent(in)::box
    real*8,dimension(:,:),intent(inout)::r_com
    real*8,intent(inout)::potential,E_elec,E_elec_nopol,E_bh, E_3body
    real*8,dimension(:),intent(inout)::energy_components
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv

    integer  :: test_trans_rot

    ! first do translation, acceptance is signaled by whether flag test_trans_rot is set negative
    test_trans_rot = 1
    call mc_make_single_mole_move(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,mass,r_com,energy_components,potential,E_elec,E_elec_nopol,E_bh,E_3body,chg,i_acc,i_try,dfti_desc,dfti_desc_inv,log_file,test_trans_rot)
    if ( test_trans_rot == 1 ) then
       accept_trans = 0 
    elseif ( test_trans_rot == -1 ) then
       accept_trans = 1
    endif

    ! now do rotation
    test_trans_rot = 2
    call mc_make_single_mole_move(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,mass,r_com,energy_components,potential,E_elec,E_elec_nopol,E_bh,E_3body,chg,i_acc,i_try,dfti_desc,dfti_desc_inv,log_file,test_trans_rot)
    if ( test_trans_rot == 2 ) then
       accept_rot = 0 
    elseif ( test_trans_rot == -2 ) then
       accept_rot = 1
    endif

  end subroutine decouple_trans_rot_accept



  !************************************************************************
  ! this subroutine carries out lattice hops between cavities in a framework
  ! this is in place of random single particle displacements
  !************************************************************************
  subroutine single_molecule_lattice_hop(xyz,r_com,box,new_xyz,new_r_com,tot_n_mole,n_mole,n_atom,newmol,weight)
    use insertion_bias_routines
    use global_variables
    integer,intent(in)::n_mole,tot_n_mole
    integer,intent(out) :: newmol
    integer,dimension(:),intent(in)::n_atom
    real*8, intent(out) :: weight
    real*8,dimension(:,:,:),intent(in)::xyz
    real*8,dimension(:,:,:),intent(out)::new_xyz
    real*8,dimension(:,:),intent(in)::box
    real*8,dimension(:,:),intent(in)::r_com
    real*8,dimension(:,:),intent(out)::new_r_com

    real*8, dimension(3):: delta_x,rand3,unit,mol_axis
    real*8 :: rand,old_weight_fac,new_weight_fac
    real*8  :: projec,rel_x
    integer :: i,j,i_atom, count
    real*8,parameter:: small=1D-4 

!!!!!!!!!!!!!!pick a molecule to move at random
    call random_number(rand)
    newmol=ceiling(rand*n_mole)
    new_r_com = r_com
    new_xyz = xyz

    ! make sure energy bias is not being used since we want purely cavity bias
    Select Case(energy_bias)
    Case(1)
       stop "cannot use energy bias with single_molecule_lattice_hop"
    End Select

    ! check that this old location is in legal cavity
    call gcmc_cavity_bias(newmol,n_mole,tot_n_mole,n_atom,new_xyz,new_r_com,box,old_weight_fac,0)
    ! pick new lattice site for move
    call gcmc_cavity_bias(newmol,n_mole,tot_n_mole,n_atom,new_xyz,new_r_com,box,new_weight_fac,1)

    ! pick random orientation and create new molecular orientation
    call random_unit_vector(unit)

    do i_atom = 1,n_atom(newmol)
       delta_x(:)=new_xyz(newmol,i_atom,:)-new_r_com(newmol,:)
       if(dot_product(delta_x,delta_x) > small) then
          mol_axis = delta_x
          exit
       endif
    enddo

    mol_axis=mol_axis/sqrt(dot_product(mol_axis,mol_axis))

    do i_atom = 1,n_atom(newmol)
       delta_x(:)=new_xyz(newmol,i_atom,:)-new_r_com(newmol,:)
       rel_x = sqrt(dot_product(delta_x,delta_x))
       projec=dot_product( delta_x , mol_axis)
       if(projec > 0) then
          new_xyz(newmol,i_atom,:) = rel_x * unit(:) + new_r_com(newmol,:)
       else
          new_xyz(newmol,i_atom,:) =  -rel_x * unit(:) + new_r_com(newmol,:)
       endif
    enddo

    ! in order to obey detailed balance, in this case we have that the attempt probabilities are the same, as we are choosing legal lattice spots with the same probablity.  The only thing that we have to be careful of is that we are attempting a displacement from a legal lattice site.  If this is not true, old_weight_fac will be zero.  Therefore we need to output this value and reject these moves

    weight = old_weight_fac

  end subroutine single_molecule_lattice_hop





  !*************************************************************************
  ! this subroutine controls the volume move part of an npt ensemble simulation
  ! it calls an internal subroutine to make the volume move and change molecule positions
  ! and then it does energy evaluations, decides whether to accept move, and updates variables
  !*************************************************************************
  subroutine npt_volume_move(box,iteration,tot_n_mole,n_mole,n_atom,n_atom_drude,alist,xyz,xyz_drude,force_atoms,mass,r_com,energy_components,potential,E_elec,E_elec_nopol,E_bh,E_3body,chg,i_acc_v,i_try_v,dfti_desc,dfti_desc_inv,log_file)
    use global_variables
    use pairwise_interaction
    use MKL_DFTI
    use pme_routines
    use total_energy_forces
    integer,intent(in)::tot_n_mole,n_mole,log_file
    integer,intent(out)::iteration
    integer,dimension(:),intent(in)::n_atom,n_atom_drude
    character(*),dimension(:,:),intent(in) :: alist
    real*8,intent(inout)::i_acc_v,i_try_v
    real*8,dimension(:,:,:),intent(inout)::xyz,xyz_drude,force_atoms
    real*8,dimension(:,:),intent(in)::mass,chg
    real*8,dimension(:,:),intent(inout)::box,r_com
    real*8,intent(inout)::potential,E_elec,E_elec_nopol,E_bh,E_3body
    real*8,dimension(:),intent(inout)::energy_components
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv
    real*8,dimension(MAX_N_MOLE,3)::new_r_com
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: new_xyz,new_xyz_drude,xyz_drude_in,force_atoms_try,lj_force
    real*8,dimension(size(energy_components)) :: energy_components_try
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
    call update_disp_lrc_noshift(n_mole, n_atom, alist,new_box)

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
    call calculate_total_force_energy( force_atoms_try, potential_try, E_elec_try, E_elec_nopol_try, E_bh_try, E_3body_try, E_bond,E_angle,E_dihedral,energy_components_try, iteration, tot_n_mole, n_mole, n_atom, n_atom_drude, new_r_com, new_xyz, chg, new_box, dfti_desc,dfti_desc_inv,log_file, new_xyz_drude, xyz_drude_initial=xyz_drude_in )
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
       energy_components = energy_components_try
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

500    continue

  end subroutine npt_volume_move



  !***********************************************************************
  ! This subroutine generates a volume move by varying ln(V)
  ! Currently, this is specific to a cubic box
  !***********************************************************************
  subroutine gen_vol_move(old_box,new_box,n_mole,n_atom,old_xyz,new_xyz,old_cofm,new_cofm,mass,vol_step,oldvol,newvol)
    use global_variables
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



  !*************************************************************************
  ! this subroutine controls the insertion/deletion portion of a grand canonical simulation
  ! this subroutine will also be used for expanded grand canonical ensemble, in which case
  ! the optional argument insert_in will determine whether there is an insertion or deletion
  !
  ! because this subroutine is either adding or deleting molecules, it has to restructure the data
  ! arrays.  This is because, for computational efficiency, solute molecules are indexed from
  ! i_mole=1,n_mole, while fixed framework atoms are indexed from i_mole=n_mole+1, tot_n_mole
  ! therefore, if a solute molecule is deleted, all data arrays must be shifted back an index after the
  ! deleted molecule.  Also if a solute molecule is inserted, all framework atoms must be shift up one
  ! index, so that this new solute molecule can be inserted at the n_mole+1 index.
  !
  ! Because of this data array restructuring, it is not convenient to just call one energy 
  ! calculation routine such as "calculate_total_force_energy".  Therefore, the energy calculation
  ! is necessarily split up over the subroutine
  !
  ! explicit three body interactions such as three body dispersion cannot be used (the code will stop)
  ! in this ensemble, therefore we don't need to consider this energy contribution
  !
  !*************************************************************************
  subroutine grand_canonical_ins_rm(box,iteration,alist,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,quaternion,force_atoms,mass,r_com,energy_components,potential,E_elec,E_elec_nopol,E_bh,chg,i_acc_p,i_try_p,dfti_desc,dfti_desc_inv,log_file,insert_in,i_mole_in)
    use insertion_bias_routines
    use global_variables
    use egc_routines
    use pairwise_interaction
    use MKL_DFTI
    use eq_drude
    use electrostatic
    use pme_routines
    use mc_routines
    integer,intent(in)::log_file
    character(*), dimension(:,:),intent(inout) :: alist
    integer,intent(inout)::n_mole,tot_n_mole
    integer,intent(out)::iteration
    integer,dimension(:),intent(inout)::n_atom, n_atom_drude
    real*8,intent(inout)::i_acc_p,i_try_p
    real*8,dimension(:,:,:),intent(inout)::xyz,xyz_drude,force_atoms
    real*8,dimension(:,:),intent(inout)::quaternion
    real*8,dimension(:,:),intent(inout)::mass,chg
    real*8,dimension(:,:),intent(in)::box
    real*8,dimension(:,:),intent(inout)::r_com
    real*8,intent(inout)::potential,E_elec,E_elec_nopol,E_bh
    real*8,dimension(:),intent(inout) :: energy_components
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv
    integer,intent(in),optional :: insert_in,i_mole_in

    integer::new_n_mole,new_tot_n_mole
    character(MAX_ANAME), dimension(MAX_N_MOLE,MAX_N_ATOM) :: new_alist
    integer,dimension(MAX_N_MOLE)::new_n_atom,new_n_atom_drude,new_molecule_index,temp_molecule_index
    real*8,dimension(MAX_N_MOLE,3)::new_r_com
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: new_xyz,new_xyz_drude,xyz_drude_in,force_atoms_try,lj_force
    real*8,dimension(MAX_N_MOLE,4)::new_quaternion
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: junk_xyz
    real*8,dimension(MAX_N_MOLE,MAX_N_ATOM)::new_mass,new_chg,tmp_chg
    integer,dimension(MAX_N_MOLE,MAX_N_ATOM)::drude_atoms,store_atom_index,store_drude_atom_map
    real*8, dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE,6) :: store_atype_lj_parameter
    real*8,dimension(pme_grid,pme_grid,pme_grid)::CBtemp
    real*8,parameter::small=1D-4
    real*8,dimension(size(energy_components)) :: denergy_components,junk_energy_components
    real*8::arg,fac,or_fac,accept_crit,randnum,randnum3(3),unit(3),u_axis(3),disp(3),p
    real*8:: rosen_fac,cavity_fac,lambda,theta,egc_fac,puddle_delta_E_try,mass_local
    real*8 :: E_elec_try,E_elec_nopol_try, Ewald_self_store, E_bh_try,dE_bh,pot_try,vol
    integer :: i_mole, i_atom, i_drude, i, j, is_overlap,count,r_mole, insert, n_framework, i_framework,ins_type,n_mole_formal,partial_molecule_index_try,partial_molecule_index_store


    vol = volume ( box(1,:),box(2,:),box(3,:) )
    !****************************debug**************************************!
    if(debug .eq. 2) then
       call date_and_time(date,time)
       write(*,*) "gc_ins_rm subroutine started at", time
    endif
    !***********************************************************************!
    n_framework = tot_n_mole - n_mole
    is_overlap=0
    new_n_mole=n_mole
    new_tot_n_mole = tot_n_mole
    new_alist=alist
    new_n_atom=n_atom
    new_n_atom_drude = n_atom_drude
    new_r_com=r_com
    new_xyz=xyz
    new_xyz_drude=xyz_drude
    new_quaternion=quaternion
    new_mass=mass
    new_chg=chg
    new_molecule_index = molecule_index
    Ewald_self_store = Ewald_self

    ! partial_molecule_index will be used in egc routines, so we need to change it and store old value
    partial_molecule_index_store = partial_molecule_index

    ! need to change global_variable atom_index array, since this is what will be used in LJ routines
    ! same for drude_atom_map, as this will be used in electrostatic routines
    ! therefore, store the original to substitute back at end if move is rejected
    store_atom_index = atom_index
    store_drude_atom_map = drude_atom_map
    ! atype_lj_parameter will be changed if egc ensemble
    store_atype_lj_parameter = atype_lj_parameter


    ! if this subroutine is called by egc ensemble, insertion or deletion is controlled by that
    Select Case(select_ensemble)
    Case("egc")
       insert = insert_in
    Case("uvt")

       i_try_p=i_try_p+1

       !!decide whether to insert or remove a molecule, insert=0 means remove, =1 means insert

       call random_number(randnum)
       if(randnum .lt. .5 ) then
          insert=0
          if( n_mole .eq. 0) then
             ! important to reject this deletion at zero molecules, rather than just automatically insert at zero molecules, otherwise violate detailed balance
             goto 300
          endif
       else
          insert=1
       endif

    End Select

    ! a molecule is being removed
    if ( insert .eq. 0) then
       Select Case(select_ensemble)
       Case("egc")
          r_mole = i_mole_in
       Case("uvt")
          ! if grand canonical ensemble, randomly choose molecule, and condense arrays
          call random_number(randnum)
          r_mole=ceiling(randnum*real(n_mole))
       End Select

       !! we are rejecting insertion if selected molecule is too close, so do same for deletion for detailed balance
       is_overlap=check_distance( xyz, tot_n_mole, n_atom, r_mole, xyz(r_mole,:,:), box, mass, r_com )
       if(is_overlap .eq. 1 ) then
          goto 300
       endif

       Select Case(pme_disp)
       Case("no")
          ! calculate lj deletion energy before changing atom_index
          call update_lennard_jones_ins( dE_bh,denergy_components, tot_n_mole, n_atom, xyz, r_mole, box, r_com )
          dE_bh= -dE_bh
          denergy_components = -denergy_components
       End Select


       ! calculate any insertion bias before changing atom_index
       ! ****************** cavity/energy bias for deletion*****************************!
       Select case ( cavity_bias)
       case("yes")
          call gcmc_cavity_bias(r_mole,n_mole,tot_n_mole,n_atom,xyz,r_com,box,cavity_fac,insert)
          ! cavity fac defined as weight_cubelet*volume/cubelet_vol
          cavity_fac = cavity_fac * vol / cubelet_vol
       case("no")
          cavity_fac = 1d0
       case Default
          stop " cavity_bias selection not recognized in grand_canonical_ins_rm"
       end select
       !***************** orientation bias for deletion**********************************!
       if( orientation_try > 1 ) then
!!! must generate orientations for removed molecule to obey detailed balance, notice the rosenbluth factors are defined slightly differently than in Frenkel and Smit
!!! junk_xyz is passed so that same subroutine can be used for insertion, in which case corresponding array contains selected configuration
          call gcmc_orientation_bias(r_mole,n_mole,tot_n_mole,n_atom,n_atom_drude,xyz,junk_xyz,box,r_com,chg,rosen_fac,insert)  
       else
          rosen_fac= 1D0
       endif
       !*************************************************************************************!

       new_n_mole = new_n_mole - 1
       new_tot_n_mole = new_tot_n_mole - 1

       do i_mole=r_mole, new_tot_n_mole
          new_alist(i_mole,:)= new_alist(i_mole+1,:)
          new_n_atom(i_mole) = new_n_atom(i_mole+1)
          new_n_atom_drude(i_mole) = new_n_atom_drude(i_mole+1)
          new_r_com(i_mole,:) = new_r_com(i_mole+1,:)
          new_xyz(i_mole,:,:) = new_xyz(i_mole+1,:,:)
          new_xyz_drude(i_mole,:,:) = new_xyz_drude(i_mole+1,:,:)
          new_quaternion(i_mole,:) = new_quaternion(i_mole+1,:)
          new_chg(i_mole,:) = new_chg(i_mole+1,:)
          new_mass(i_mole,:) = new_mass(i_mole+1,:)
          new_molecule_index(i_mole) = new_molecule_index(i_mole+1)
          atom_index(i_mole,:)= atom_index(i_mole+1,:)
          drude_atom_map(i_mole,:) = drude_atom_map(i_mole+1,:)
       enddo



    else
       ! molecule is being inserted, move all framework atoms up one position, and put new molecule in new_n_mole position
       new_n_mole = new_n_mole + 1
       new_tot_n_mole = new_tot_n_mole + 1    
       if (n_framework > 0 ) then
          do i_mole =1, n_framework
             i_framework = new_tot_n_mole - i_mole
             new_alist(i_framework+1,:)= new_alist(i_framework,:)
             new_n_atom(i_framework+1) = new_n_atom(i_framework)
             new_n_atom_drude(i_framework+1) = new_n_atom_drude(i_framework)
             new_r_com(i_framework+1,:) = new_r_com(i_framework,:)
             new_xyz(i_framework+1,:,:) = new_xyz(i_framework,:,:)
             new_chg(i_framework+1,:) = new_chg(i_framework,:)
             atom_index(i_framework+1,:) = atom_index(i_framework,:)
             drude_atom_map(i_framework+1,:) = drude_atom_map(i_framework,:)
          enddo
       endif

!!! now select random position to insert center of mass of new molecule
       call random_number(randnum3)
       new_r_com(new_n_mole,:)= randnum3(1)*box(1,:)+randnum3(2)*box(2,:)+randnum3(3)*box(3,:)

!!! now select type of molecule to insert
       call random_number(randnum)
       ins_type= ceiling(real(n_molecule_type)* randnum)
       !! only reversible if 1 molecule type
       if ( n_molecule_type .gt. 1 ) then
          stop " more than one type of solute molecule.  Detailed balance is probably broken as selection of molecule to insert is inconsistent with selection of molecule to delete"
       endif

       Select Case(select_ensemble)
       Case("egc")
          ! for egc, set partial_molecule_index
          partial_molecule_index = new_n_mole
       End Select

       !! generate new data for inserted molecule
       new_molecule_index(new_n_mole) = ins_type
       call gen_insertion_data(new_n_mole,new_alist,new_n_atom,new_n_atom_drude,new_r_com,new_xyz,new_quaternion,new_chg,new_mass,ins_type)


       ! ****************** cavity/energy bias for insertion*****************************!
       Select case ( cavity_bias)
       case("yes")
          call gcmc_cavity_bias(new_n_mole,new_n_mole,new_tot_n_mole,new_n_atom,new_xyz,new_r_com,box,cavity_fac,insert)
          ! cavity fac defined as weight_cubelet*volume/cubelet_vol
          cavity_fac = cavity_fac * vol / cubelet_vol
       case("no")
          cavity_fac = 1d0
       case Default
          stop " cavity_bias selection not recognized in grand_canonical_ins_rm"
       end select
       !********************* orientation bias for insertion ****************************!
       if( orientation_try > 1 ) then
!!! generate biased orientaion of new molecule, notice the rosenbluth factors are defined slightly differently than in Frenkel and Smit
          call gcmc_orientation_bias(new_n_mole,new_n_mole,new_tot_n_mole,new_n_atom,new_n_atom_drude,new_xyz,new_xyz,box,new_r_com,new_chg,rosen_fac,insert)  
       else
          rosen_fac= 1D0
       endif
       !*************************************************************************************!

       !! for the newly inserted molecule, put drude oscillators on atoms for lack of better choice
         if ( new_n_atom_drude(new_n_mole) > new_n_atom(new_n_mole) ) then
             do i_drude = new_n_atom(new_n_mole)+1 , new_n_atom_drude(new_n_mole)
                ! use map to figure out which atom this oscillator is on
                i_atom = drude_atom_map(new_n_mole,i_drude)
                new_xyz_drude(new_n_mole,i_drude,:) = new_xyz(new_n_mole,i_atom,:)
             enddo
         endif

       !!check distance between newly inserted molecule and others
       is_overlap=check_distance( new_xyz, new_tot_n_mole, new_n_atom, new_n_mole, new_xyz(new_n_mole,:,:), box, new_mass, new_r_com )
       if(is_overlap .eq. 1 ) then
!!$           write(log_file,*) "*******************************************************"
!!$           write(log_file,*) " atoms too close, swap move rejected before energy calculation"
!!$           write(log_file,*) "*******************************************************"
          goto 300
       endif


       Select Case(pme_disp)
       Case("no")
          call update_lennard_jones_ins( dE_bh,denergy_components, new_tot_n_mole, new_n_atom, new_xyz, new_n_mole, box, new_r_com )
       End Select

    endif

    ! if pme,remember to update Ewald_self interaction, as there is a new particle number

    if( drude_simulation .eq. 1 ) then
       Select Case(electrostatic_type)
       Case("pme")
          call update_Ewald_self( new_tot_n_mole, new_n_atom_drude, new_chg )
       end Select

       xyz_drude_in=new_xyz_drude

       call scf_drude(E_elec_try,E_elec_nopol_try,force_atoms_try,new_xyz,new_chg,new_r_com,box,new_tot_n_mole,new_n_mole,new_n_atom,new_n_atom_drude,iteration,new_xyz_drude,dfti_desc,dfti_desc_inv,log_file,xyz_drude_in)

       if(iteration .gt. iteration_limit) then
          goto 300
       endif

    else
       do i_mole=1,new_n_mole
          do i_atom=1,new_n_atom(i_mole)
             drude_atoms(i_mole,i_atom)=i_atom
          enddo
       enddo
       Select Case(electrostatic_type)
       Case("pme")
          call update_Ewald_self( new_tot_n_mole, new_n_atom, new_chg )
       end Select

       call Electrostatic_energy_force( E_elec_try, force_atoms_try,drude_atoms,new_tot_n_mole, new_n_mole, new_n_atom, new_xyz, new_chg, box, new_r_com,dfti_desc,dfti_desc_inv)

       ! this variable is meaningless without oscillators, but assign a value so compiler doesn't get angry
       E_elec_nopol_try = E_elec_nopol

    endif



    Select Case(pme_disp)
    Case("yes")
       call tot_disp_usegrid( E_bh_try,junk_energy_components,new_tot_n_mole, new_n_mole, new_n_atom, new_xyz, box, new_r_com )
    Case("no")
       !!     for debugging update_lennard_jones_ins
        !   call lennard_jones( E_bh_try,junk_energy_components,new_tot_n_mole, new_n_mole, new_n_atom, new_xyz, box, new_r_com )

       E_bh_try= E_bh + dE_bh

    End Select

    call update_disp_lrc_noshift(new_n_mole, new_n_atom, alist, box,1)

    pot_try = E_elec_try + E_bh_try + E_disp_lrc_try

    !**********************************puddle filling, accepted insertion or deletion randomly pick new molecule for puddle hamiltonian
    Select Case(puddle_filling_on)
    Case("yes")
       if ( new_n_mole .eq. 0 ) then
          partial_molecule_index_try=0
          puddle_delta_E_try=0.
       else
          partial_molecule_index_store=partial_molecule_index
          call random_number(randnum)
          partial_molecule_index_try=ceiling(randnum*real(new_n_mole))
          partial_molecule_index=partial_molecule_index_try
          call puddle_fill_potential(new_n_mole,new_tot_n_mole,new_n_atom,new_n_atom_drude,new_xyz,box,new_r_com,new_chg,puddle_delta_E_try)
          partial_molecule_index=partial_molecule_index_store
          potential=potential+ puddle_delta_E
          pot_try = pot_try + puddle_delta_E_try
       endif
    end Select
    !*******************************************************************************************************************************

    !*********************************************************************
    ! commented lines involving theta, refer to the rotational contribution of a diatomic molecule to 
    ! the chemical potential.  If the chemical potential defined in global_variables which is used in the acceptance
    ! criteria includes these rotational contributions, then it is important to be consistent for the acceptance rule
    ! usually, the chemical potential in the input file will not include these contributions, so these lines remain commented
    ! in the code
    !**********************************************************************


    if(insert .eq. 0) then
       mass_local = molecule_mass(molecule_index(r_mole))
!       theta= 48.506D0/inertia(r_mole,1,1)
       lambda= 126.178D0/sqrt(2D0*pi*mass_local*8.314D0*temp)
!       fac= lambda**3 * n_mole / vol*theta/temp * exp(-1000./(temp*8.314)*chem_potential)
       fac= lambda**3 * dble(n_mole) / vol* exp(-1000./(temp*8.314)*chem_potential)
       or_fac = rosen_fac * dble(orientation_try)
       cavity_fac = cavity_fac
    else
       mass_local = molecule_mass(new_molecule_index(new_n_mole))
!       theta= 48.506D0/new_inertia(new_n_mole,1,1)
       lambda= 126.178D0/sqrt(2D0*pi*mass_local*8.314D0*temp)
!       fac= vol /(lambda**3 * new_n_mole)* temp/theta * exp(1000./(temp*8.314)*chem_potential)
       fac= vol /(lambda**3 * dble(new_n_mole))* exp(1000./(temp*8.314)*chem_potential)
       or_fac = 1D0/rosen_fac/dble(orientation_try)
       cavity_fac = 1d0/cavity_fac
    endif

    arg=-(1000./(temp*8.314))*(pot_try-potential)

    Select Case(select_ensemble)
    Case("uvt")
       accept_crit=min(1.,cavity_fac*or_fac*fac*exp(arg))
    Case("egc")
       Select Case(insert)
          ! n_mole_formal is the number of fully coupled molecules in the system, which doesn't
          ! depend on initial or final states in this case since we are only changing the partial coupling
       Case(0)
          temp_molecule_index=molecule_index
          n_mole_formal= new_n_mole
       Case(1)
          temp_molecule_index=new_molecule_index
          n_mole_formal= n_mole
       End Select
       ! get egc weight factor if using that ensemble
       ! use n_mole for either insertion/deletion, since for insertion we 
       ! formally don't have a new molecule until it's fully coupled
       call get_egc_fac(egc_fac,insert,box,temp_molecule_index,n_mole_formal)
       accept_crit=min(1.,egc_fac*exp(arg))
    End Select


    call random_number(p)
    if ( p < accept_crit) then

       ! move accepted. add lennard jones forces and update data
       call lennard_jones( E_bh_try, lj_force, new_tot_n_mole, new_n_mole, new_n_atom, new_xyz, box, new_r_com)
       n_mole=new_n_mole
       tot_n_mole = new_tot_n_mole
       alist=new_alist
       n_atom=new_n_atom
       n_atom_drude=new_n_atom_drude
       mass=new_mass
       chg=new_chg
       molecule_index=new_molecule_index
       xyz = new_xyz
       xyz_drude = new_xyz_drude
       quaternion = new_quaternion
       force_atoms = force_atoms_try + lj_force
       r_com = new_r_com
       E_elec = E_elec_try
       E_elec_nopol = E_elec_nopol_try
       E_bh = E_bh_try
       E_disp_lrc=E_disp_lrc_try

       Select Case(puddle_filling_on)
       Case("yes")
          puddle_delta_E=puddle_delta_E_try
       End Select

       potential = pot_try
       i_acc_p = i_acc_p + 1.

       ! update depends on if this subroutine is called by expanded ensemble
       Select Case(select_ensemble)
       Case("uvt")
          energy_components = energy_components + denergy_components
       Case("egc")
          call update_partial_molecule(insert,n_atom,n_atom_drude,chg,alist)
       end select


    else
300    continue
       atom_index = store_atom_index
       drude_atom_map = store_drude_atom_map
       Ewald_self = Ewald_self_store
       ! for failed insertion in egc ensemble, must restore some global variables
       Select Case(select_ensemble)
       Case("egc")
          Select Case(insert)
          Case(1)
             partial_molecule_index=partial_molecule_index_store
             n_atom_type = n_atom_type_store
             atype_lj_parameter = store_atype_lj_parameter
          End Select
       End Select
    endif


    !****************************debug**************************************!
    if(debug .eq. 2) then
       call date_and_time(date,time)
       write(*,*) "gc_ins_rm subroutine finished at", time
    endif
    !***********************************************************************!

  end subroutine grand_canonical_ins_rm



  !**********************************************************
  !  This subroutine implements coupling changes in the expanded grand canonical ensemble.  It basically just controls 
  !  the modification of the global parameter arrays when the coupling of a molecule is changed.  When molecule becomes either fully coupled
  !  or uncoupled, this corresponds to an insertion or deletion, and the grand_canonical_ins_rm subroutine can take care of this once the relevant 
  !  parameter arrays have been updated
  !
  !  references are : Escobedo F. A., Jour. Chem. Phys. 127, 174104 (2007) and Escobedo F. A., Jour. Chem. Phys. 105, 4391, 1996 
  !*************************************************************
  subroutine expanded_grand_canonical_change_coupling(box,iteration,alist,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude, quaternion, force_atoms,mass,r_com,energy_components,potential,E_elec,E_elec_nopol,E_bh,chg,i_acc_p,i_try_p,dfti_desc,dfti_desc_inv,log_file)
    use egc_routines
    use global_variables
    use MKL_DFTI
    use pairwise_interaction
    use eq_drude
    use electrostatic
    use pme_routines

    integer,intent(in)::log_file
    character(*), dimension(:,:),intent(inout) :: alist
    integer,intent(inout)::n_mole,tot_n_mole
    integer,intent(out)::iteration
    integer,dimension(:),intent(inout)::n_atom,n_atom_drude
    real*8,intent(inout)::i_acc_p,i_try_p
    real*8,dimension(:,:,:),intent(inout)::xyz,xyz_drude,force_atoms
    real*8,dimension(:,:),intent(inout):: quaternion
    real*8,dimension(:,:),intent(inout)::mass,chg
    real*8,dimension(:,:),intent(in)::box
    real*8,dimension(:,:),intent(inout)::r_com
    real*8,intent(inout):: potential,E_elec,E_elec_nopol,E_bh
    real*8,dimension(:),intent(inout) :: energy_components
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv

    real*8   :: randnum,p,arg,accept_crit
    integer  :: i,i_atom,i_drude,i_mole,increment,flag
    integer,save :: initialize=0
    real*8   :: pot_try,dE_bh_old,dE_bh_new,Ewald_self_store,E_bh_try,E_elec_try,E_elec_nopol_try,egc_fac
    integer,dimension(MAX_N_MOLE,MAX_N_ATOM)::drude_atoms
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: new_xyz_drude,xyz_drude_in,force_atoms_try,lj_force
    integer, dimension(MAX_N_MOLE,MAX_N_ATOM):: store_atom_index 
    real*8, dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE,6) :: store_atype_lj_parameter
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM) :: new_chg
    real*8,dimension(size(energy_components)) :: junk_components



    ! initialize weight, use naive choice here, also check that we have at least 3 steps
    ! also store n_atom_type before it is messed with
    if ( initialize .eq. 0 ) then
       do i=1,partial_coupling_steps
          step_weight(i) = dble(i)* coupling_increment
       enddo
       n_atom_type_store = n_atom_type
       initialize=1
       if ( partial_coupling_steps .lt. 3 ) then
          stop " partial_coupling_steps must be at least 3"
       endif
       ! make sure orientation bias/cavity bias aren't being used
       if ( ( orientation_try > 1 ) .or. (cavity_bias .eq. "yes" ) ) then
          stop "cavity/orientation bias not implemented for egc ensemble"
       endif
    endif


    ! store force field parameter global arrays, since these will have to be changed
    ! to compute energies, and will need to be changed back if move is rejected
    ! also store Ewald_self, as this will need replaced for rejected move
    store_atype_lj_parameter = atype_lj_parameter
    store_atom_index = atom_index
    Ewald_self_store = Ewald_self
    new_chg = chg
    new_xyz_drude=xyz_drude

    i_try_p=i_try_p+1

    ! our first step is to decide whether we are increasing or decreasing the coupling parameter
    call random_number(randnum)
    if(randnum .lt. .5 ) then
       increment=0
       if( n_mole .eq. 0) then
          ! important to reject this step for detailed balance, just as in gcmc
          goto 400
       endif
    else
       increment=1
    endif


    ! now we need to change global array data, this depends on the increment type, and also the coupling state that we are in.
    ! partial_molecule_coupling can be in range 0 to (partial_coupling_steps - 1)
    ! cases 0 and 1 are special, in the sense that 0 -> 1 adds a new molecule to all global arrays, while case 1 -> 0 removes the partial molecule from all global arrays
    ! these cases will use grand_canonical_ins_rm subroutine

    ! flag will tell us whether grand_canonical_ins_rm is being called
    flag=0

    Select Case(partial_molecule_coupling)
    Case(0)
       Select Case (increment)
       Case(1)
          ! partial molecule coupling 0 --> 1 is essentially an insertion
          ! modify global parameters in subroutine, partial_molecule_index will be set in subroutine
          call grand_canonical_ins_rm(box,iteration,alist,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,quaternion,force_atoms,mass,r_com,energy_components,potential,E_elec,E_elec_nopol,E_bh,chg,i_acc_p,i_try_p,dfti_desc,dfti_desc_inv,log_file,increment)
          flag=1
       Case(0)
          ! here we are starting to uncouple a new molecule, choose which molecule this is, then change parameters
          call random_number(randnum)
          partial_molecule_index = ceiling(randnum*real(n_mole))

          call update_lennard_jones_ins( dE_bh_old,junk_components, tot_n_mole, n_atom, xyz, partial_molecule_index, box, r_com )
          call adjust_partial_molecule_parameters(increment,n_atom,new_chg)
       end Select

    Case(1)
       Select Case (increment)
       Case(0)
          ! partial_molecule_coupling steps 1 --> 0 is essentially a deletion
          call grand_canonical_ins_rm(box,iteration,alist,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,quaternion,force_atoms,mass,r_com,energy_components,potential,E_elec,E_elec_nopol,E_bh,chg,i_acc_p,i_try_p,dfti_desc,dfti_desc_inv,log_file,increment,partial_molecule_index) 
          flag=1
       Case(1)
          ! just changing coupling
          call update_lennard_jones_ins( dE_bh_old,junk_components, tot_n_mole, n_atom, xyz, partial_molecule_index, box, r_com )
          call adjust_partial_molecule_parameters(increment,n_atom,new_chg)
       end select
    Case default
       ! we are done with the two special cases, here we just change coupling
       call update_lennard_jones_ins( dE_bh_old,junk_components, tot_n_mole, n_atom, xyz, partial_molecule_index, box, r_com )
       call adjust_partial_molecule_parameters(increment,n_atom,new_chg)

    end Select

    ! if grand_canonical_ins_rm was not called, we need to compute energy differences for this configuration

    Select Case(flag)
    Case(0)
       ! get lj energy with scaled parameters
       call update_lennard_jones_ins( dE_bh_new,junk_components, tot_n_mole, n_atom, xyz, partial_molecule_index, box, r_com )

       ! if pme,remember to update Ewald_self interaction, as charges on partial molecule have been rescaled
       if( drude_simulation .eq. 1 ) then
          Select Case(electrostatic_type)
          Case("pme")
             call update_Ewald_self( tot_n_mole, n_atom_drude, new_chg )
          end Select

          xyz_drude_in=new_xyz_drude

          call scf_drude(E_elec_try,E_elec_nopol_try,force_atoms_try,xyz,new_chg,r_com,box,tot_n_mole,n_mole,n_atom,n_atom_drude,iteration,new_xyz_drude,dfti_desc,dfti_desc_inv,log_file,xyz_drude_in)

          if(iteration .gt. iteration_limit) then
             goto 400
          endif

       else
          do i_mole=1,n_mole
             do i_atom=1,n_atom(i_mole)
                drude_atoms(i_mole,i_atom)=i_atom
             enddo
          enddo
          Select Case(electrostatic_type)
          Case("pme")
             call update_Ewald_self( tot_n_mole, n_atom, new_chg )
          end Select

          call Electrostatic_energy_force( E_elec_try, force_atoms_try,drude_atoms,tot_n_mole, n_mole, n_atom, xyz, new_chg, box, r_com,dfti_desc,dfti_desc_inv)
          E_elec_nopol_try = E_elec_try

       endif

       E_bh_try = E_bh + dE_bh_new - dE_bh_old
       pot_try = E_elec_try + E_bh_try


       ! get egc weight factor
       call get_egc_fac(egc_fac,increment,box,molecule_index,n_mole)

       arg=-(1000./(temp*8.314))*(pot_try-potential)
       accept_crit=min(1.,egc_fac*exp(arg))


       call random_number(p)
       if ( p < accept_crit) then
          !! add lennard jones forces and update data
          call lennard_jones( E_bh_try, lj_force, tot_n_mole, n_mole, n_atom, xyz, box, r_com)
          chg=new_chg
          xyz_drude = new_xyz_drude
          force_atoms = force_atoms_try + lj_force
          E_elec = E_elec_try
          E_elec_nopol = E_elec_nopol_try
          E_bh = E_bh_try

          potential = pot_try
          i_acc_p = i_acc_p + 1.

          ! update partial molecule information
          call update_partial_molecule(increment,n_atom,n_atom_drude,chg,alist)

       else
400       continue
       ! restore partial molecule information, atom_index may have been changed if start of decoupling
          Select Case(partial_molecule_coupling)
          Case(0)
             ! this is necessarily a decoupling, as an insertion would be treated in grand_canonical_ins_rm
             atom_index = store_atom_index
             n_atom_type = n_atom_type_store
          End Select
          atype_lj_parameter = store_atype_lj_parameter
          Ewald_self = Ewald_self_store

          ! 
       endif

    End Select  ! this ends select on flag, so all of this code was for a move evaluated internally, rather than evaluated in grand_canonical_ins_rm

  end subroutine expanded_grand_canonical_change_coupling


end module sampling
