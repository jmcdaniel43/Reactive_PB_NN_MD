module eq_drude

contains


!******************************************************************!
! this subroutine finds optimal drude oscillator positions
! for a given configuration of molecules
!
!  the number of Drude oscillators on each molecule (say i_mole) is determined by the difference
!  in the value of n_atom_drude(i_mole) and n_atom(i_mole), with the former storing the total number of 
!  atoms and drude oscillators on each molecule, and the later just storing the number of atoms
!
!  The atom that each drude oscillator is connected to can be accessed through the array
!  drude_atom_map
!
!  global variables used:
!  real*8,parameter::springcon (.1)
!  real*8,parameter::thole     (2.)
!  real*8,parameter:: force_threshold
!
!  work in energy units of e^2/A, which is unit in which electrostatics are calculated
!  input array "chg" should be total charges including drudes
!
! this subroutine uses a conjugate gradient method to find equilibrium 
! drude positions
! based off Lindan, P.J.D, Gillan,M.J., J. Phys.: Condens. Matter 5 (1993) 1019-1030
!******************************************************************!
  subroutine scf_drude(energy,energy_nopol,force_atoms,xyz,chg,r_com,box,tot_n_mole,n_mole,n_atom,n_atom_drude,iteration,xyz_drude,dfti_desc,dfti_desc_inv,log_file,xyz_drude_in)
    use pme_routines
    use electrostatic
    use global_variables
    use MKL_DFTI
    use routines
    real*8,intent(out)::energy,energy_nopol
    integer,intent(in)::tot_n_mole,n_mole,log_file
    integer,intent(out)::iteration
    integer,dimension(:),intent(in)::n_atom, n_atom_drude
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3),intent(in) :: xyz
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3),intent(out) :: xyz_drude,force_atoms
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM),intent(in) :: chg
    real*8, dimension(3,3),intent(in)::box
    real*8, dimension(MAX_N_MOLE,3),intent(in) :: r_com
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3),intent(in),optional :: xyz_drude_in
    real*8,dimension(MAX_N_MOLE,MAX_N_ATOM,3)::force,force_old,force_new,spring_force
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) ::tot_xyz
    integer,dimension(MAX_N_MOLE,MAX_N_ATOM)::drude_atoms
    real*8,dimension(MAX_N_MOLE,MAX_N_ATOM)::temp_chg
    real*8,dimension(MAX_N_MOLE,MAX_N_ATOM,3)::search_d
    integer::i,i_mole,i_atom,i_drude,converged,tot_drudes,flag,index
    real*8::Beta,lambda1,sum_f_old,sum_f_new,sum_lambda_n,sum_lambda_d,Ewald_self_store
    real*8,dimension(3)::delta_xyz,disp_xyz,unit
    real*8,parameter::drude_initial=.005
    real*8::springE
    real*8,parameter::drude_max=1.


   !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "drude minimization started at", time
    endif
    !***********************************************************************!


    iteration=0
    converged=0

    ! since in general we only loop over Drude oscillators, we need to set spring_force array to zero so that we can add it to atoms without drude oscillators at the end
    spring_force=0d0


    ! Number of Drude oscillators on each molecule is given by n_atom_drude(i_mole) - n_atom(i_mole)
    if(present(xyz_drude_in)) then
       xyz_drude = xyz_drude_in
    else
!!!!!!!!!!!!!! if first call to this subroutine, put drude oscillators on respective atoms, otherwise use passed array
!!!!!!!!!!!!!!! put drude oscillators slightly off atoms so that erf(x)/x doesn't blow up, also distribute oscillators randomly so no net dipole of system

       do i_mole=1,n_mole
          if ( n_atom_drude(i_mole) > n_atom(i_mole) ) then

             do i_drude = n_atom(i_mole)+1 , n_atom_drude(i_mole)
                ! use map to figure out which atom this oscillator is on
                i_atom = drude_atom_map(i_mole,i_drude)
                call random_unit_vector(unit)

                xyz_drude(i_mole,i_drude,:)= xyz(i_mole,i_atom,:) + drude_initial * unit

             enddo
          endif
       enddo
    endif


    ! make sure to include framework positions
    tot_xyz=xyz
    ! forces shouldn't be calculated on any atoms except oscillators
    drude_atoms=0

!!!!!!!!!!!!!! combine positions of atoms and oscillators into one array to pass to force subroutine
!!!!!!!!!!!!!! remember that chg on atom equals permanent chg - charge on drude oscillator
!!!!!!!!!!!!!! while looping, construct drude_atoms array, which tells force routine for which atoms to compute forces

    do i_mole=1,n_mole
       if ( n_atom_drude(i_mole) > n_atom(i_mole) ) then
          do i_drude = n_atom(i_mole)+1 , n_atom_drude(i_mole)

             tot_xyz(i_mole,i_drude,:)=xyz_drude(i_mole,i_drude,:)

             ! drude_atoms array should be filled, starting at the first index,  with indices of atoms for which forces should be computed
             index = i_drude - n_atom(i_mole)
             ! index should never be greater than n_atom(i_mole)
             if ( index > n_atom(i_mole) ) then
                stop "there is an error in indexing the drude oscillators !" 
             endif
             drude_atoms(i_mole,index)=i_drude
          enddo
       endif
    enddo



    call Electrostatic_energy_force( energy, force,drude_atoms,tot_n_mole, n_mole, n_atom_drude, tot_xyz, chg, box, r_com,dfti_desc,dfti_desc_inv)
!!!!!!!!!convert back from KJ/mol/A to e^2/A^2
    force=force/0.52914D0/627.51D0/4.184D0    


    do while(converged.eq.0)
       iteration=iteration+1
       converged=1

       if(iteration .eq. 1) then
          sum_lambda_n=0.D0
          lambda1 = 1.D0/springcon
          springE=0.D0
          sum_f_new=0.D0
          tot_drudes=0
!!!!!!!!!!!!if first iteration, search direction is in direction of force
!!!!!!!!!!add spring force,and first search direction is only determined by initial force
          do i_mole=1,n_mole
             if ( n_atom_drude(i_mole) > n_atom(i_mole) ) then
                do i_drude = n_atom(i_mole)+1 , n_atom_drude(i_mole)

                   ! use map to figure out which atom this oscillator is on
                   i_atom = drude_atom_map(i_mole,i_drude)

                   tot_drudes=tot_drudes+1
                   disp_xyz(:)=xyz_drude(i_mole,i_drude,:)-xyz(i_mole,i_atom,:)
!!!!!!!!!!total force including spring
                   ! remember, electrostatic force for this drude oscillator is mapped from
                   ! force array by the drude_atoms array that we constructed earlier
                   index = i_drude - n_atom(i_mole)
                   force(i_mole,i_drude,:)=force(i_mole,index,:)-springcon*disp_xyz(:)
                   sum_f_new = sum_f_new + dot_product(force(i_mole,i_drude,:),force(i_mole,i_drude,:))

!!!!!!!!! first trial positions determined by steepest descent
                   search_d(i_mole,i_drude,:)=force(i_mole,i_drude,:)
                   delta_xyz(:)= lambda1*search_d(i_mole,i_drude,:)
                   xyz_drude(i_mole,i_drude,:)=xyz_drude(i_mole,i_drude,:)+delta_xyz(:)

                   ! NOTE INDEX here:  spring_force array is only used to output total forces on atoms, therefore, use the atom index here
                   spring_force(i_mole,i_atom,:) = -springcon*disp_xyz(:)

                   tot_xyz(i_mole,i_drude,:)=xyz_drude(i_mole,i_drude,:)
                   sum_lambda_n=sum_lambda_n+dot_product(force(i_mole,i_drude,:),force(i_mole,i_drude,:))
                   disp_xyz(:)=xyz_drude(i_mole,i_drude,:)-xyz(i_mole,i_atom,:)
                   springE = springE + .5D0*springcon*dot_product(disp_xyz,disp_xyz)
                enddo
             endif
          enddo
!!!!!!!!!! check total net force, if rms force < threshold, converged
          sum_f_new = sum_f_new/dble(tot_drudes)

          if ( sum_f_new > force_threshold**2) then
             converged=0
          endif
!!!!!!!!!!!!!!!!!!if converged, lets get outta here
          if(converged.eq.1) then
             goto 100
          endif
       else

!!!!!!!!!!if not the first iteration, search direction is determined using information from previous step
!!!!calculate Beta, which is contribution of previous direction
          sum_f_old=0.;sum_f_new=0.
          do i_mole=1,n_mole
             if ( n_atom_drude(i_mole) > n_atom(i_mole) ) then
                do i_drude = n_atom(i_mole)+1 , n_atom_drude(i_mole)
                   ! use map to figure out which atom this oscillator is on
                   i_atom = drude_atom_map(i_mole,i_drude)
                   disp_xyz(:)=xyz_drude(i_mole,i_drude,:)-xyz(i_mole,i_atom,:)

                   ! NOTE INDEX here:  spring_force array is only used to output total forces on atoms, therefore, use the atom index here
                   spring_force(i_mole,i_atom,:) = -springcon*disp_xyz(:)

                   ! remember, electrostatic force for this drude oscillator is mapped from
                   ! force array by the drude_atoms array that we constructed earlier
                   index = i_drude - n_atom(i_mole)
                   force(i_mole,i_drude,:)=force(i_mole,index,:)-springcon*disp_xyz(:)
                   sum_f_new = sum_f_new+dot_product(force(i_mole,i_drude,:),force(i_mole,i_drude,:))
                   sum_f_old = sum_f_old+dot_product(force_old(i_mole,i_drude,:),force_old(i_mole,i_drude,:))
                enddo
             endif
          enddo
          Beta=sum_f_new/sum_f_old
!!!!!!!!!! check total net force, if rms force < threshold, converged
          sum_f_new = sum_f_new / dble(tot_drudes)
          if ( sum_f_new > force_threshold**2) then
             converged=0
          endif
!!!!!!!!!!!!!!!!!!if converged, lets get outta here
          if(converged.eq.1) then
             goto 100
          endif

!!!!!!!!!calculate directions of moves, and lambda
          sum_lambda_n=0.D0;sum_lambda_d=0.D0
          do i_mole=1,n_mole
             if ( n_atom_drude(i_mole) > n_atom(i_mole) ) then
                do i_drude = n_atom(i_mole)+1 , n_atom_drude(i_mole)
!!!!!!!!!!!!!!!!now iterations,like conjugate gradient, use direction of previous move
                   search_d(i_mole,i_drude,:)=force(i_mole,i_drude,:)+Beta*search_d(i_mole,i_drude,:)
                   sum_lambda_n=sum_lambda_n+dot_product(force(i_mole,i_drude,:),search_d(i_mole,i_drude,:))
                   sum_lambda_d=sum_lambda_d+springcon*dot_product(search_d(i_mole,i_drude,:),search_d(i_mole,i_drude,:))
                enddo
             endif
          enddo

          lambda1=sum_lambda_n/sum_lambda_d
          springE=0.D0
!!!!!!!!!!!!!!calculate positions of oscillators for this value of lambda
          do i_mole=1,n_mole
             if ( n_atom_drude(i_mole) > n_atom(i_mole) ) then
                do i_drude = n_atom(i_mole)+1 , n_atom_drude(i_mole)
                   ! use map to figure out which atom this oscillator is on
                   i_atom = drude_atom_map(i_mole,i_drude)

                   delta_xyz(:) = lambda1 * search_d(i_mole,i_drude,:)
                   xyz_drude(i_mole,i_drude,:)=xyz_drude(i_mole,i_drude,:)+delta_xyz(:)
                   tot_xyz(i_mole,i_drude,:)=xyz_drude(i_mole,i_drude,:)
                   disp_xyz(:)=xyz_drude(i_mole,i_drude,:)-xyz(i_mole,i_atom,:)
                   springE=springE + .5D0*springcon*dot_product(disp_xyz,disp_xyz)
                enddo
             endif
          enddo

       endif


!!!!!!!!!!!!!!!!!!!!!now calculate forces at new drude positions
       call Electrostatic_energy_force( energy, force_new,drude_atoms,tot_n_mole, n_mole, n_atom_drude, tot_xyz, chg, box, r_com,dfti_desc,dfti_desc_inv)

       force_new=force_new/0.52914D0/627.51D0/4.184D0

       ! need to loop over all atoms here, as the force_old array needs the forces
       ! in the drude indices, for the conjugate gradient minimization
       ! and the force array needs the forces in the atom positions, as this is where
       ! the above loop is expecting them
       do i_mole=1,n_mole
          do i_drude = 1 , n_atom_drude(i_mole)
             force_old(i_mole,i_drude,:) = force(i_mole,i_drude,:)
             force(i_mole,i_drude,:) = force_new(i_mole,i_drude,:)
          enddo
       enddo

100    continue

!!!!!!!!!!!!!!!!!!!!! if drude oscillators are not converging

       if(iteration .gt. iteration_limit) then
          write(log_file,*) "*******************************************************"
          write(log_file,*) " drude oscillators did not converge during this move, so move was rejected"
          write(log_file,*) "*******************************************************"
          exit
       endif


    enddo


    !when we get here, iterations have either converged or passed the max limit


    flag=0
    do i_mole=1,n_mole
       if ( n_atom_drude(i_mole) > n_atom(i_mole) ) then
          do i_drude = n_atom(i_mole)+1 , n_atom_drude(i_mole)
             ! use map to figure out which atom this oscillator is on
             i_atom = drude_atom_map(i_mole,i_drude)

             disp_xyz(:)=xyz_drude(i_mole,i_drude,:)-xyz(i_mole,i_atom,:)
             if(dot_product(disp_xyz,disp_xyz) > drude_max) then
                flag=1
             endif
          enddo
       endif
    enddo

    if(flag .eq.1) then
       write(log_file,*) "*******************************************************"
       write(log_file,*) " there is a drude oscillator farther than", drude_max, "away from its atom"
       write(log_file,*) "*******************************************************"
    endif


    drude_atoms=0
    ! calculate forces on atoms to output for hybrid mc/md
    do i_mole=1,n_mole
       do i_atom=1,n_atom(i_mole)
          drude_atoms(i_mole,i_atom)=i_atom
       enddo
    enddo


!!!!!!!!!!!!!!!!!!!! if drude oscillators did not converge, calculate force without polarization
    if(iteration .gt. iteration_limit) then
       do i_mole=1,n_mole
          if ( n_atom_drude(i_mole) > n_atom(i_mole) ) then
             do i_drude = n_atom(i_mole)+1 , n_atom_drude(i_mole)
                ! use map to figure out which atom this oscillator is on
                i_atom = drude_atom_map(i_mole,i_drude)

                xyz_drude(i_mole,i_drude,:)=xyz(i_mole,i_atom,:)
                tot_xyz(i_mole,i_drude,:) = xyz(i_mole,i_atom,:)
             enddo
          endif
       enddo
    endif

    call Electrostatic_energy_force( energy, force_atoms,drude_atoms,tot_n_mole, n_mole, n_atom_drude, tot_xyz, chg, box, r_com,dfti_desc,dfti_desc_inv)


!!!!!!!!!!!!!!!!!!!!!if oscillators converged, add spring terms
    if(iteration .le. iteration_limit) then
       force_atoms=force_atoms - spring_force * 0.52914D0*627.51D0*4.184D0 
       energy=energy + springE * 0.52914D0*627.51D0*4.184D0 
    endif

    !! if energy decomposition is requested, call energy routine once more with atoms and static charges (no drudes) for energy breakdown
    !! need to change ewald self for pme
    Select case (energy_decomposition)
    Case("yes")
       temp_chg = chg
       do i_mole=1,n_mole
          if ( n_atom_drude(i_mole) > n_atom(i_mole) ) then
             do i_drude = n_atom(i_mole)+1 , n_atom_drude(i_mole)
                ! use map to figure out which atom this oscillator is on
                i_atom = drude_atom_map(i_mole,i_drude)

                temp_chg(i_mole,i_atom) = chg(i_mole,i_atom) + chg(i_mole, i_drude)
             enddo
          endif
       enddo

       Ewald_self_store = Ewald_self
       Select case (electrostatic_type)
       Case("pme")
          call update_Ewald_self( tot_n_mole, n_atom, temp_chg )   
       end Select

    call Electrostatic_energy_force( energy_nopol, force_atoms,drude_atoms,tot_n_mole, n_mole, n_atom_drude, tot_xyz, chg, box, r_com,dfti_desc,dfti_desc_inv)

       Ewald_self = Ewald_self_store

    end Select


   !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "drude minimization finished at", time
    endif
    !***********************************************************************!


  end subroutine scf_drude

end module eq_drude

