!***********************************************************
! this module contains routines related to cavity or energy bias
! as well as orientation bias
!**************************************************************

module insertion_bias_routines
use routines
contains


!*************************************************************
! this subroutine calculates an orientationally averaged boltzmann factor
! for a chosen number of cubelets in a unit cell for a potential
! lattice site model at low loadings
!*************************************************************
  subroutine energy_grid(box,alist,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,mass,r_com,chg,dfti_desc,dfti_desc_inv,log_file)
    use egc_routines
    use global_variables
    use MKL_DFTI
    use pairwise_interaction
    use eq_drude
    use electrostatic
    use pme_routines
    use total_energy_forces

    integer,intent(in)::log_file
    character(*), dimension(:,:),intent(in) :: alist
    integer,intent(in)::n_mole,tot_n_mole
    integer,dimension(:),intent(in)::n_atom,n_atom_drude
    real*8,dimension(:,:,:),intent(in)::xyz
    real*8,dimension(:,:),intent(in)::mass,chg
    real*8,dimension(:,:),intent(in)::box
    real*8,dimension(:,:),intent(in)::r_com
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv
    integer,dimension(3)  :: cav_grid
    real*8,dimension(3,3) :: cubelet_vec
    real*8,dimension(cav_grid_a,cav_grid_b,cav_grid_c,3) :: cavity_bias_grid_center
    real*8,dimension(cav_grid_a,cav_grid_b,cav_grid_c) :: boltzmann_config, min_energy_config

    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: xyz_config,force_atoms, xyz_drude
    real*8,dimension(MAX_N_MOLE,3)     ::   new_r_com
    real*8, dimension(orientation_try,MAX_N_ATOM,3):: config_try
    real*8, dimension(orientation_try):: Energy_orientation,dE_elec,dE_bh


    ! other variables
    integer           :: i,j,k,ix,iy,iz,i_mole,j_mole,i_atom,j_atom,j_drude,count,i_orient,is_overlap,iteration
    real*8            :: E_ins_min,vol
    real*8,dimension(3) :: randvec
    real*8,dimension(5) :: energy_components
    real*8,parameter :: small = 1D-4, max_rep=100.
    integer :: orientation_low_e,count_avg
    real*8 :: E_elec, E_elec_nopol, E_bh, E_3body, potential
   ! these will not be used for rigid molecule simulation
    real*8 :: E_bond, E_angle, E_dihedral


    vol = volume ( box(1,:),box(2,:),box(3,:) )

    ! first collapse cavity grid size parameters from global_variables into a local array for use
    cav_grid(1) = cav_grid_a; cav_grid(2) = cav_grid_b; cav_grid(3) = cav_grid_c; 
    ! find cubelet vectors
    do i = 1,3
       cubelet_vec(i,:) = box(i,:) / dble (cav_grid(i))
    enddo
    ! cubelet  volume
    cubelet_vol = volume ( cubelet_vec(1,:),cubelet_vec(2,:),cubelet_vec(3,:) )


    ! assuming one molecule in the input file for energy grid
    if ( n_mole .ne. 1 ) then
       stop "must have one and only one solute molecule for energy grid"
    endif
   
    ! pme will always be used for electrostatics
    write(*,*) "pme is always used for electrostatics in energy grid (unless 'none' option is used)"

    ! ************************ write output******************************************
    write( log_file, *) "*******run parameters***************************"
    write( log_file, *) " code version = ", code_version
    write( log_file, *) " this is an energy grid run"
    write( log_file, *) "lj_bkghm ", lj_bkghm
    write( log_file, *) "lj cutoff", lj_cutoff
    write( log_file, *) "lj_shift ", lj_shift
    write( log_file, *) "lj_lrc ", lj_lrc
    write( log_file, *) "electrostatics","   ", electrostatic_type
    write( log_file, *) "electrostatic cutoff", ewald_cutoff, Electro_cutoff
    write( log_file, *) "grid error function","   ", grid_erfc
    write( log_file, *) "electrostatic screening", screen_type
    write( log_file, *) "drude simulation", drude_simulation
    write( log_file, *) "distance for assumed repulsive wall (A)", too_close
    write( log_file, *) "temperature", temp
    write( log_file, *) "volume", vol
    write( log_file, *) "cavity grid (x,y,z) " , cav_grid_a , cav_grid_b , cav_grid_c
    write( log_file, *) "number of orientations per cubelet ",orientation_try
    write( log_file, *) "cubelet_volume" , cubelet_vol

    i_mole = 1

    ! now loop over all cubelets, calculating their centers, and calculating energies if  there is no overlap
    ! center vectors are in cartesian basis, not unit cell basis, but these are usually the same

    boltzmann_config=0d0
    xyz_config = xyz
    new_r_com = r_com
    count=0
    do i = 1,cav_grid(1)
       do j=1,cav_grid(2)
          do k=1,cav_grid(3)

             cavity_bias_grid_center(i,j,k,:) = (dble(i) - .5d0) * cubelet_vec(1,:) + (dble(j) - .5d0) * cubelet_vec(2,:)+(dble(k) - .5d0) * cubelet_vec(3,:)   

             E_ins_min = max_rep

             do i_orient=1, orientation_try

                ! random position in grid cell
                call random_number(randvec)
                new_r_com(i_mole,:) = cavity_bias_grid_center(i,j,k,:) + (randvec(1)-0.5d0)*cubelet_vec(1,:) &
                     & + (randvec(2)-0.5d0)*cubelet_vec(2,:) + (randvec(3)-0.5d0)*cubelet_vec(3,:)

                ! using the first molecule
                do i_atom=1,n_atom(i_mole)
                   xyz_config(i_mole,i_atom,:) = xyz(i_mole,i_atom,:) - r_com(i_mole,:) + new_r_com(i_mole,:)
                enddo
                call random_rotation( xyz_config, new_r_com,n_atom, i_mole )

                ! check distance from nearest atom
                is_overlap=check_distance( xyz_config, tot_n_mole, n_atom, i_mole, xyz_config(i_mole,:,:), box, mass, new_r_com )
                if(is_overlap .eq. 1 ) then
                   ! if too close, set as max repulsive value
                   Energy_orientation(i_orient) = max_rep
                else

                   call calculate_total_force_energy( force_atoms,potential, E_elec,E_elec_nopol,E_bh, E_3body, E_bond, E_angle, E_dihedral, energy_components, iteration, tot_n_mole, n_mole, n_atom, n_atom_drude, new_r_com, xyz_config, chg, box, dfti_desc,dfti_desc_inv,log_file,xyz_drude)       

                   ! calculate total energy
                   Energy_orientation(i_orient) = potential
                   if (Energy_orientation(i_orient) < E_ins_min ) then
                      E_ins_min = Energy_orientation(i_orient)
                   endif
                endif

                boltzmann_config(i,j,k) = boltzmann_config(i,j,k) + dexp(-Energy_orientation(i_orient)*1000/8.314d0/temp)


             enddo


             min_energy_config(i,j,k)= E_ins_min

             boltzmann_config(i,j,k) = boltzmann_config(i,j,k) / dble(i_orient)

          enddo
       enddo
    enddo


    ! now print to output file

    do i = 1,cav_grid(1)
       do j=1,cav_grid(2)
          do k=1,cav_grid(3)

             write( log_file, *) " cavity index " , i,j,k
             write( log_file, *) cavity_bias_grid_center(i,j,k,:)
             write( log_file, *) "min energy (over orientations) ", min_energy_config(i,j,k)
             write( log_file, *) "orientationally averaged boltzmann factor for this cubelet ", boltzmann_config(i,j,k)
          enddo
       enddo
    enddo



  end subroutine energy_grid




!******************************************************
! this subroutine carries out cavity/energy bias insertions and deletions
!
! global variables that are specific to cavity bias
! integer,parameter :: cav_grid_a, cav_grid_b, cav_grid_c
! integer,parameter :: energy_bias
! real*8,parameter :: min_frmwk_dist
! real*8           :: cubelet_vol 
!
! global variables used but not changed
! real*8,parameter :: temp 
!
! global variables that will be changed only temporarily
! integer, dimension(MAX_N_MOLE,MAX_N_ATOM):: atom_index 
! 
! if energy bias =0, then this is naive cavity bias.  If energy bias =1, then 
! probabilities of inserting into allowed cubelets is not constant, but rather depends on a
! boltzmann factor computed with lennard_jones interactions of first solute atom type with framework
!
! saved variable init_cavity =0 will intiallize grid, and then will be set to 1
! also save cubelet vectors, these will only be needed locally
! and save cav_grid array which will be used locally
! saved array cavity_bias_grid_weight saves weights of cubelets
! saved array cavity_bias_grid_center saves vector from origin to center of cubelet
!
! references
! Mezei, Molecular Physics, 1980, vol 40, no 4, 901-906
! Snurr, R.Q. et al, JPC, 1993,97,13742-13752
!*******************************************************
  subroutine gcmc_cavity_bias(t_mole,n_mole,tot_n_mole,n_atom,xyz,r_com,box,weight_fac,insert)
    use global_variables
    use pairwise_interaction
    integer,intent(in)::t_mole,n_mole,tot_n_mole,insert
    integer,dimension(:),intent(in)::n_atom
    real*8,dimension(:,:,:),intent(inout)::xyz
    real*8,dimension(:,:),intent(in)::box
    real*8,dimension(:,:),intent(inout)::r_com
    real*8,intent(out):: weight_fac

    ! saved quantities
    integer, save     :: init_cavity = 0
    integer, save     :: eligible_cubes
    integer,dimension(3),save  :: cav_grid
    real*8,dimension(3,3),save :: cubelet_vec
    real*8,dimension(cav_grid_a,cav_grid_b,cav_grid_c),save :: cavity_bias_grid_weight
    real*8,dimension(cav_grid_a,cav_grid_b,cav_grid_c,3),save :: cavity_bias_grid_center
    integer,dimension(cav_grid_a*cav_grid_b*cav_grid_c,3),save:: eligible_index
    real*8,dimension(cav_grid_a*cav_grid_b*cav_grid_c),save :: eligible_sum_weight
    ! temporary variables for calculating energy bias
    real*8,dimension(size(xyz(:,1,1)),size(xyz(1,:,1)),3) :: temp_xyz
    real*8,dimension(size(r_com(:,1)),3) :: temp_r_com
    integer,dimension(size(n_atom))      :: temp_n_atom
    integer                              :: temp_n_mole, temp_tot_n_mole
    ! important, store global array atom_index
    integer,dimension(size(atom_index(:,1)),size(atom_index(1,:)))::temp_atom_index
    ! other variables
    integer           :: i,j,k,i_mole,i_atom,f_mole,count,i_select
    real*8            :: sum_lj,E_ins,weight,rand
    integer,dimension(3)  :: sel_index
    real*8,dimension(3) :: shift, dr_com,rand_vec,new_r_com,delta_com
    real*8,dimension(5) :: denergy_components
    real*8,parameter :: small = 1D-10

    Select case (init_cavity)
    Case(0)

       !*****************might be bug in energy bias, not fully tested************
       if (energy_bias .eq. 1 ) then
          write(*,*) "************************ WARNING ******************************"
          write(*,*) "energy bias has not been fully tested, there might be a bug in it"
       end if

       !*********************initialize and carry out cavity bias*******************************************************************************************************!  
       ! first collapse cavity grid size parameters from global_variables into a local array for use
       cav_grid(1) = cav_grid_a; cav_grid(2) = cav_grid_b; cav_grid(3) = cav_grid_c; 

       ! find cubelet vectors
       do i = 1,3
          cubelet_vec(i,:) = box(i,:) / dble (cav_grid(i))
       enddo

       ! cubelet  volume
       cubelet_vol = volume ( cubelet_vec(1,:),cubelet_vec(2,:),cubelet_vec(3,:) )

       ! if energy bias, must prepare global arrays for update_lennard_jones_ins subroutine, as this is not smart enough
       ! to only calculate solute-framework interactions ( and not solute-solute ) unless we modify some global arrays
       temp_atom_index = atom_index
       temp_n_mole = 1
       temp_tot_n_mole = tot_n_mole - n_mole + 1
       ! atom_index determines parameters for lj to use, shift all framework parameters down, and use first atom type for biased insertion
       atom_index(1,1) = 1
       do i = 1, temp_tot_n_mole
          ! only considering one atom per molecule, the test atom is the first molecule, all other framework atoms are molecules with 1 atom
          temp_n_atom(i) = 1
          i_mole = i + 1
          f_mole = i + n_mole
          atom_index(i_mole,1) = atom_index(f_mole,1)
          temp_xyz(i_mole,1,:) = xyz(f_mole,1,:)
          temp_r_com(i_mole,:) = r_com(f_mole,:)
       enddo

       ! now loop over all cubelets, calculating their centers, and determining if they are eligible for insertions and deletions
       ! if this is energy bias, weight factor depends on lennard jones based boltzmann factor
       ! center vectors are in cartesian basis, not unit cell basis, but these are usually the same

       count=0
       sum_lj=0d0
       do i = 1,cav_grid(1)
          do j=1,cav_grid(2)
             do k=1,cav_grid(3)
                weight = 1d0
                cavity_bias_grid_center(i,j,k,:) = (dble(i) - .5d0) * cubelet_vec(1,:) + (dble(j) - .5d0) * cubelet_vec(2,:)+(dble(k) - .5d0) * cubelet_vec(3,:)   

                ! search over framework atoms (each framework atom treated as a molecule, so use r_com) to test cubelet eligibility
                ! remember solute stored up to n_mole, framework n_mole+1 to tot_n_mole
                do i_mole = n_mole+1, tot_n_mole
                   shift = pbc_shift( r_com(i_mole,:), cavity_bias_grid_center(i,j,k,:), box , xyz_to_box_transform )
                   dr_com = pbc_dr( r_com(i_mole,:), cavity_bias_grid_center(i,j,k,:), shift )
                   if (dot_product( dr_com,dr_com ) < min_frmwk_dist**2) then
                      weight = 0d0
                   endif
                enddo

                if (weight > 0d0) then
                   count = count +1
                   ! this cubelet is in a "cavity" and is thus eligible for insertion/deletion
                   select case(energy_bias)
                   case(0)
                      eligible_sum_weight(count) = dble(count)
                   case(1)
                      ! for energy bias, weight of cubelet depends on boltzman factor with lj energy of test particle at center of cubelet
                      temp_xyz(1,1,:)= cavity_bias_grid_center(i,j,k,:)
                      temp_r_com(1,:)= cavity_bias_grid_center(i,j,k,:)
                      call update_lennard_jones_ins( E_ins,denergy_components, temp_tot_n_mole, temp_n_atom, temp_xyz, 1 , box, temp_r_com )

                      weight = exp (- E_ins * 1000d0/(temp*8.314d0) )        ! energy from lj is in KJ/mol
                      sum_lj = sum_lj + weight
                      eligible_sum_weight(count) = sum_lj
                   end select

                   ! this is an eligible cubelet, store this for possible insertions
                   eligible_index(count,1)=i; eligible_index(count,2)=j; eligible_index(count,3)=k;

                endif
                cavity_bias_grid_weight(i,j,k) = weight
             enddo
          enddo
       enddo

       eligible_cubes = count
       
       ! fix atom_index array
       atom_index = temp_atom_index

       ! normalize weights such that sum(weight) = 1
       if (energy_bias .eq. 1 ) then
          cavity_bias_grid_weight = cavity_bias_grid_weight/ sum_lj
          eligible_sum_weight = eligible_sum_weight / sum_lj
       else
          cavity_bias_grid_weight = cavity_bias_grid_weight/dble(count)
          eligible_sum_weight = eligible_sum_weight / dble(count)
       endif

       ! done with initialization, now do bias

       Select Case (insert)
       Case (1) 
          ! this is an insertion, pick cubelet according to weight distribution using binary search

          call pick_cubelet(sel_index,eligible_cubes,eligible_sum_weight,eligible_index)

          ! now we have our favorite cubelet. New center of mass is at a random location in this cubelet
          call random_number(rand_vec)
          new_r_com(:) = cavity_bias_grid_center(sel_index(1),sel_index(2),sel_index(3),:) + ( rand_vec(1)-.5d0 ) * cubelet_vec(1,:)+ ( rand_vec(2)-.5d0 ) * cubelet_vec(2,:)+ ( rand_vec(3)-.5d0 ) * cubelet_vec(3,:)
          ! now scoot the molecule over
          do i_atom =1, n_atom(t_mole)
             delta_com(:) = xyz(t_mole,i_atom,:) - r_com(t_mole,:)
             xyz(t_mole,i_atom,:) = delta_com(:) + new_r_com(:)
          enddo
          ! now update center of mass
          r_com(t_mole,:) = new_r_com(:)
          ! and output weight factor of cubelet
          weight_fac = cavity_bias_grid_weight(sel_index(1),sel_index(2),sel_index(3))

       Case (0)
          ! this is deletion
          ! we need to figure out which cubelet our targeted molecule is in
          ! this is easiest if our unit cell basis vectors are orthogonal and lie along cartesian basis

          if ( (box(1,2).lt. small) .and. (box(1,3).lt. small) .and. (box(2,1).lt. small) .and. (box(2,3).lt. small).and. (box(3,1).lt. small) .and. (box(3,2).lt. small) ) then
             do i=1,3
                sel_index(i) =  ceiling (r_com(t_mole,i)/box(i,i) * cav_grid(i) )
             enddo

          else
             stop "cavity bias is not implemented for desired box type.  Modify the source code or change boxes"
          endif

          ! output weight factor of this cubelet
          weight_fac = cavity_bias_grid_weight(sel_index(1),sel_index(2),sel_index(3))

       end select

       ! this is the end of the case (init_cavity=0) set this saved variable to one so we don't do the first part of code again

       init_cavity = 1


    Case(1)
       !*********************just carry out cavity bias, initialization has been done*******************************************************************************************************!  

  Select Case (insert)
       Case (1) 
          ! this is an insertion, pick cubelet according to weight distribution using binary search

          call pick_cubelet(sel_index,eligible_cubes,eligible_sum_weight,eligible_index)

            ! now we have our favorite cubelet. New center of mass is at a random location in this cubelet
            call random_number(rand_vec)
            new_r_com(:) = cavity_bias_grid_center(sel_index(1),sel_index(2),sel_index(3),:) + ( rand_vec(1)-.5d0 ) * cubelet_vec(1,:)+ ( rand_vec(2)-.5d0 ) * cubelet_vec(2,:)+ ( rand_vec(3)-.5d0 ) * cubelet_vec(3,:)
            ! now scoot the molecule over
            do i_atom =1, n_atom(t_mole)
               delta_com(:) = xyz(t_mole,i_atom,:) - r_com(t_mole,:)
               xyz(t_mole,i_atom,:) = delta_com(:) + new_r_com(:)
            enddo
            ! now update center of mass
            r_com(t_mole,:) = new_r_com(:)
            ! and output weight factor of cubelet
            weight_fac = cavity_bias_grid_weight(sel_index(1),sel_index(2),sel_index(3))

         Case (0)
            ! this is deletion
            ! we need to figure out which cubelet our targeted molecule is in
            ! this is easiest if our unit cell basis vectors are orthogonal and lie along cartesian basis

            if ( (box(1,2).lt. small) .and. (box(1,3).lt. small) .and. (box(2,1).lt. small) .and. (box(2,3).lt. small).and. (box(3,1).lt. small) .and. (box(3,2).lt. small) ) then
               do i=1,3
                  sel_index(i) =  ceiling (r_com(t_mole,i)/box(i,i) * cav_grid(i) )
               enddo

            else
               stop "cavity bias is not implemented for desired box type.  Modify the source code or change boxes"
            endif

            ! output weight factor of this cubelet
            weight_fac = cavity_bias_grid_weight(sel_index(1),sel_index(2),sel_index(3))

         end select
   

    end select


  end subroutine gcmc_cavity_bias


!*****************************************************************
! this subroutine picks a cubelet with probability corresponding to its weight
! and returns the indices of this cubelet
! this is basically a binary search routine
!*****************************************************************
  subroutine pick_cubelet(sel_index,max,weight,index)
    integer,dimension(:),intent(out) :: sel_index
    integer, intent(in) :: max
    real*8,dimension(:),intent(in) :: weight
    integer,dimension(:,:),intent(in) :: index

    integer  :: low, mid, high, i_select
    real*8,parameter :: bound_low = .99, bound_high = 1.01
    real*8           :: rand

    ! check to make sure that distribution is normalized and cube number is correct
    if ( ( bound_low > weight(max) ) .or. ( bound_high < weight(max) ) ) then
       stop "weight distribution for cublets is not normalized in cavity bias"
    endif

    call random_number(rand)

    if ( rand < weight(1) ) then
       i_select = 1
       goto 200
    endif

    low=1; high= max
    do 
       mid = low + (high - low ) / 2
       if ( (weight(mid) < rand ) .and. ( rand < weight(mid+1) ) ) then
          i_select = mid + 1
          exit
       else if ( weight(mid+1) < rand ) then
          low = mid + 1
       else
          high = mid - 1 
       endif
    enddo

200 continue

sel_index(:) = index(i_select,:)

end subroutine pick_cubelet



!***************************************************************************
! this subroutine carries out orientation_bias for insertion of linear molecules
! as based on Frenkel and Smit
!***************************************************************************
subroutine gcmc_orientation_bias(t_mole,n_mole,tot_n_mole,n_atom,n_atom_drude,xyz,new_xyz,box,r_com,chg,rosen_fac,insert)
  use global_variables
  use pairwise_interaction
  use electrostatic
  use pme_routines
  integer,intent(in)::t_mole,n_mole,tot_n_mole,insert
  integer,dimension(:),intent(in)::n_atom,n_atom_drude
  real*8,dimension(:,:,:),intent(in)::xyz
  real*8,dimension(:,:,:),intent(inout)::new_xyz
  real*8,dimension(:,:),intent(in)::box
  real*8,dimension(:,:),intent(in)::chg
  real*8,dimension(:,:),intent(in)::r_com
  real*8,intent(out):: rosen_fac

  real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: xyz_config
  real*8, dimension(orientation_try,MAX_N_ATOM,3):: config_try
  real*8, dimension(MAX_N_MOLE,MAX_N_ATOM) :: atom_chg
  real*8,parameter::small=1D-4
  real*8,dimension(3)::disp,u_axis,unit
  real*8,dimension(orientation_try)::dE_elec,dE_bh
  real*8,dimension(5) :: denergy_components
  integer::i,i_orient,i_atom,i_drude,i_mole, i_select
  real*8 :: sum,sum1,p,term
  real*8,parameter:: boltz_min= 1D-10, boltz_max = 1D8     

  do i_atom=1,n_atom(t_mole)
     disp(:)=xyz(t_mole,i_atom,:)-r_com(t_mole,:)
     if(dot_product(disp,disp) > small) exit
  enddo
  u_axis=disp/sqrt(dot_product(disp,disp))


  !! set original orientation
  config_try(1,:,:) = xyz(t_mole,:,:)

  do i_orient=2, orientation_try
     call random_unit_vector(unit)
     do i_atom=1,n_atom(t_mole)
        disp(:)=xyz(t_mole,i_atom,:)-r_com(t_mole,:)
        if(dot_product(disp,u_axis)>0.) then
           config_try(i_orient,i_atom,:)= sqrt(dot_product(disp,disp))*unit(:)+r_com(t_mole,:)
        else
           config_try(i_orient,i_atom,:)= -sqrt(dot_product(disp,disp))*unit(:)+r_com(t_mole,:)
        endif
     enddo
  enddo

  !! now do an energy calculation without polarization
  !! if drude oscillators are being used in this simulation, must
  !! add charges to get static charges on atoms

  atom_chg = chg
  do i_mole=1,tot_n_mole
     if ( n_atom_drude(i_mole) > n_atom(i_mole) ) then
        do i_drude = n_atom(i_mole)+1 , n_atom_drude(i_mole)
           ! use map to figure out which atom this oscillator is on
           i_atom = drude_atom_map(i_mole,i_drude)
           atom_chg(i_mole,i_atom)=chg(i_mole,i_atom)+chg(i_mole,i_drude)
        enddo
     endif
  enddo


  xyz_config(:,:,:)=xyz(:,:,:)
  do i_orient=1, orientation_try
     xyz_config(t_mole,:,:) = config_try(i_orient,:,:)

     !! only need energy differences
     Select Case(electrostatic_type)
     Case("none")
        dE_elec(i_orient) = 0d0
     case default
        call update_elec_cutoff_ins( dE_elec(i_orient), tot_n_mole, n_atom, atom_chg, xyz_config, box, t_mole, r_com )
     End Select

     call update_lennard_jones_ins( dE_bh(i_orient),denergy_components, tot_n_mole, n_atom, xyz_config, t_mole , box, r_com )

  enddo

  ! energies can be extremely divergent, so for numerical reasons, put bounds on boltzmann factors
  ! maximum factor is for buckingham, since exponenential does not cancel negatively differgent r^-6

  sum = 0.
  do i_orient=1, orientation_try
     term = exp(-(1000./(temp*8.314))* (dE_elec(i_orient)+dE_bh(i_orient)))
     if( term .lt. boltz_min) then
        term = boltz_min
     elseif (term .gt. boltz_max) then
        term = boltz_max
     endif
     sum = sum + term
  enddo

  if(insert .eq. 0 ) then             
     ! this is a deletion, so don't choose configuraion
     term = exp(-(1000./(temp*8.314))* (dE_elec(1)+dE_bh(1)))
     if( term .lt. boltz_min) then
        term = boltz_min
     elseif (term .gt. boltz_max) then
        term = boltz_max
     endif
     rosen_fac = term/sum

  else
     !! this is insertion, choose orientation based on boltzmann factor
     call random_number(p)
     sum1=0.
     do i_orient=1, orientation_try
        term = exp(-(1000./(temp*8.314))* (dE_elec(i_orient)+dE_bh(i_orient)))
        if( term .lt. boltz_min) then
           term = boltz_min
        elseif (term .gt. boltz_max) then
           term = boltz_max
        endif
        sum1 = sum1 + term/ sum
        if(p < sum1) then
           i_select = i_orient
           exit
        endif
     enddo


     new_xyz(t_mole,:,:)=config_try(i_select,:,:)

     rosen_fac= term / sum

  endif

end  subroutine gcmc_orientation_bias



!***************************************************************************
! this subroutine is for puddle filling, with the purpose of preventing molecules
! from getting trapped in potential wells
!***************************************************************************
subroutine puddle_fill_potential(n_mole,tot_n_mole,n_atom,n_atom_drude,xyz,box,r_com,chg,delta_E)
  use global_variables
  use pairwise_interaction
  use electrostatic
  use pme_routines
  integer,intent(in)::n_mole,tot_n_mole
  integer,dimension(:),intent(in)::n_atom,n_atom_drude
  real*8,dimension(:,:,:),intent(in)::xyz
  real*8,dimension(:,:),intent(in)::box
  real*8,dimension(:,:),intent(in)::chg
  real*8,dimension(:,:),intent(in)::r_com
  real*8,intent(out):: delta_E

  real*8  :: dE_bh,dE_elec,randnum
  real*8,dimension(5) :: junk_components
  real*8, dimension(MAX_N_MOLE,MAX_N_ATOM) :: atom_chg
  integer :: i_atom,i_drude,i_mole
  integer,save :: initialize_puddle=0

  Select Case( initialize_puddle )
  Case(0)
     ! initialize, check for conflicts, then randomly pick a molecule to puddle fill pair potential
     if( select_ensemble .ne. "uvt" ) then
        stop "can only use puddle filling in grand canonical ensemble"
     endif
     if (hybrid_md_mc .ne. "yes" ) then
        stop "can only use puddle filling with hybrid md_mc"
     endif
     if ( puddle_min > 0. ) then
        stop "if you want puddle_min to be > 0, you gotta change the source code"
     endif
     call random_number(randnum)
     partial_molecule_index=ceiling(randnum*real(n_mole))
     initialize_puddle=1
  End Select

  ! evaluate pairwise energy for the molecule for which we are filling the pairwise puddle.  
  ! this molecule is stored in partial_molecule_index
  call update_lennard_jones_ins( dE_bh, junk_components, tot_n_mole, n_atom, xyz, partial_molecule_index, box, r_com )

  !! now do an energy calculation without polarization
  !! if drude oscillators are being used in this simulation, must
  !! add charges to get static charges on atoms

  atom_chg = chg
  do i_mole=1,tot_n_mole
     if ( n_atom_drude(i_mole) > n_atom(i_mole) ) then
        do i_drude = n_atom(i_mole)+1 , n_atom_drude(i_mole)
           ! use map to figure out which atom this oscillator is on
           i_atom = drude_atom_map(i_mole,i_drude)
           atom_chg(i_mole,i_atom)=chg(i_mole,i_atom)+chg(i_mole,i_drude)
        enddo
     endif
  enddo


  Select Case(electrostatic_type)
  Case("none")
     dE_elec = 0d0
  case default
     call update_elec_cutoff_ins( dE_elec, tot_n_mole, n_atom, atom_chg, xyz, box, partial_molecule_index, r_com )
  End Select

  delta_E = dE_bh + dE_elec


  if ( delta_E > puddle_min ) then
     delta_E = 0d0
  else
     delta_E = puddle_min - delta_E
  endif


end subroutine puddle_fill_potential


!************************************
! this subroutine does a random rotation on a polyatomic molecule
! no quaternion information is given with this rotation
!***********************************
subroutine random_rotation(xyz,r_com,n_atom, i_mole)
  use rigid_body
  real*8,dimension(:,:,:),intent(inout) :: xyz
  real*8,dimension(:,:),intent(in) :: r_com
  integer,dimension(:),intent(in) :: n_atom
  integer, intent(in ) :: i_mole

 integer:: i_atom,i,j
    real*8 :: a(3,3),a_trans(3,3), r_vec(3)
    real*8 :: R1,R2,x0,y1,y2
    real*8 :: q0,q1,q2,q3
    real*8,dimension(4) :: quaternion
  real*8, parameter :: pi=3.14159265

    ! Shoemake algorithm for choosing quaternions for random rotation matrix

    call random_number(x0)
    call random_number(y1)
    call random_number(y2)
    y1=y1*2d0*pi
    y2=y2*2d0*pi

    R1=sqrt(1d0-x0)
    R2=sqrt(x0)

    q0=dcos(y2)*R2
    q1=dsin(y1)*R1
    q2=dcos(y1)*R1
    q3=dsin(y2)*R2

    quaternion(1)=q0 ;
    quaternion(2)=q1 ;
    quaternion(3)=q2 ;
    quaternion(4)=q3 ;

    call quaternion_transformation_matrix(a,quaternion )


    do i=1,3
       do j=1,3
          a_trans(i,j) = a(j,i)
       enddo
    enddo

    do i_atom=1, n_atom(i_mole)
       r_vec(:) = xyz(i_mole,i_atom,:) - r_com(i_mole,:)
       r_vec(:) = matmul( a_trans ,r_vec )
       xyz(i_mole,i_atom,:) = r_vec(:) + r_com(i_mole,:)
    enddo


  end subroutine random_rotation



end module insertion_bias_routines
