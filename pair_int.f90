module pairwise_interaction
  use global_variables
  use routines
  implicit none

  !***********************************************************************
  ! This module contains subroutines to calculate non-electrostatic pairwise
  ! interaction energies such as Lennard-Jones or Buckingham interactions
  !
  !**********************************************************************

contains


  !***********************************************************************
  ! This subroutine computes lennard jones (buckingham) forces and energy
  ! It is implemented for generic solute, as atom index array stores the type
  ! of atom for each molecule, which is then used to recover LJ parameters in atype_lj_parameter array
  ! also, storing framework atoms in tot_n_mole - n_mole indices allows solute framework interactions to
  ! be calculated
  !  As implemented, some forces are added to framework atom indices of force array for free due to Newton's 3rd law,
  !  this is not important and does nothing.
  !
  ! we have removed some of the subroutine calls inside the loops and replaced these calls
  ! with explicit code to enhance performance
  !
  ! Either verlet list or naive looping can be used
  !***********************************************************************
  subroutine lennard_jones( lj_energy, lj_force, tot_n_mole, n_mole, n_atom, xyz, box, r_com)
    real*8,intent(out) :: lj_energy
    real*8,dimension(:,:,:),intent(out)::lj_force
    integer,intent(in)::n_mole,tot_n_mole
    integer,dimension(:),intent(in)::n_atom
    real*8,dimension(:,:,:),intent(in)::xyz
    real*8,dimension(:,:),intent(in)::box
    real*8,dimension(:,:),intent(in)::r_com
    integer::i_mole,i_atom, a_index
    real*8 :: E_intra_lj, lj_cutoff2, E_ms_evb_repulsion
    real*8,dimension(:,:), allocatable :: atom_list_xyz, atom_list_force
    integer,dimension(:), allocatable  :: atom_list_atom_index


    lj_force=0d0

    ! inter-molecular contributions
    Select Case(use_verlet_list)
    Case("no")
       call  lennard_jones_no_verlet( lj_energy, lj_force, tot_n_mole, n_mole, n_atom, xyz, box, r_com)
    Case("yes")
       ! here we need to map molecular data structures to atomic data structures that are
       ! consistent with the verlet list
       allocate( atom_list_xyz(verlet_lj_atoms,3), atom_list_force(verlet_lj_atoms,3) , atom_list_atom_index(verlet_lj_atoms) )
       do i_mole = 1, tot_n_mole
          do i_atom = 1, n_atom(i_mole)
             ! get verlet atom index for this atom
             a_index = verlet_molecule_map_lj(i_mole,i_atom)
             if (  a_index > 0  ) then
                ! store atom information from molecular data arrays
                atom_list_xyz(a_index,:) = xyz(i_mole,i_atom,:)
                atom_list_atom_index(a_index) = atom_index(i_mole,i_atom)
             end if
          enddo
       enddo

       call  lennard_jones_use_verlet( lj_energy, atom_list_force, atom_list_xyz , atom_list_atom_index, box )
       ! now add forces to molecular data structures
       do i_mole = 1, tot_n_mole
          do i_atom = 1, n_atom(i_mole)
             ! get verlet atom index for this atom
             a_index = verlet_molecule_map_lj(i_mole,i_atom)
             if (  a_index > 0  ) then
                lj_force(i_mole,i_atom,:) = atom_list_force(a_index,:)
             end if
          enddo
       enddo

       deallocate( atom_list_xyz , atom_list_force, atom_list_atom_index )
    End Select


    lj_cutoff2 = lj_cutoff ** 2
    ! now intra-molecular contributions

    do i_mole = 1, n_mole
       call intra_lennard_jones_energy_force( E_intra_lj, lj_force, i_mole, n_atom, xyz , lj_cutoff2 )
       lj_energy = lj_energy + E_intra_lj
    enddo


    ! ms-evb hydronium-water repulsive terms.  These are the additional intermolecular interaction
    ! terms that are non-Coulombic, non-Lennard Jones
    Select Case(ms_evb_simulation)
    Case("yes")
       call ms_evb_intermolecular_repulsion( lj_force , E_ms_evb_repulsion, n_mole  , n_atom , xyz, box )
       lj_energy = lj_energy + E_ms_evb_repulsion
    End Select


  end subroutine lennard_jones



  !***********************************************************************
  ! In this subroutine, the data structures have been changed from molecule based
  ! data structures, to atom_list based data structures.  This is because we only want
  ! to consider inter-molecular interactions between atoms that have non-zero
  ! lennard-jones parameters.  The indexing of atoms in the atom_list data structures
  ! should be consistent with the indexing of atoms in the verlet list
  !***********************************************************************
  subroutine lennard_jones_use_verlet( lj_energy,atom_list_force, atom_list_xyz , atom_list_atom_index, box )
    use omp_lib
    real*8,intent(out) :: lj_energy
    real*8,dimension(:,:), intent(out) :: atom_list_force
    real*8,dimension(:,:), intent(in)  :: atom_list_xyz
    integer,dimension(:), intent(in)   :: atom_list_atom_index
    real*8,dimension(:,:), intent(in)  :: box

    integer:: i_atom, j_atom, verlet_start, verlet_finish,i_index,atom_id1,atom_id2,i
    integer :: split_do
    real*8 :: lj_cutoff2, norm_dr2, norm_dr6, norm_dr12, term6, term12
    real*8,dimension(3)::f_ij,rij,shift,shift_direct, dr_com,dr_direct


    atom_list_force=0d0
    lj_energy=0d0

    lj_cutoff2 = lj_cutoff ** 2

    ! decide how to split the parallel section
    if (n_threads .eq. 1 ) then
       split_do = 1
    else
       split_do = verlet_lj_atoms/n_threads+1
    endif


    !**************************************** use Verlet list ****************************************************************
    call OMP_SET_NUM_THREADS(n_threads)
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(atom_list_xyz, atom_list_atom_index, verlet_lj_atoms, box,atype_lj_parameter,lj_cutoff2,split_do,xyz_to_box_transform,verlet_point_lj,verlet_neighbor_list_lj,lj_shift_store) REDUCTION(+:atom_list_force, lj_energy)
    !$OMP DO SCHEDULE(DYNAMIC,split_do)
    ! note that in Verlet list, neighbors ( on j_mole) are only stored if j_mole > i_mole , to avoid double counting
    do i_atom=1, verlet_lj_atoms

       atom_id1 = atom_list_atom_index(i_atom)
       ! get j_atom from Verlet list
       verlet_start = verlet_point_lj(i_atom)
       i_index = i_atom + 1
       verlet_finish = verlet_point_lj(i_index) - 1

       ! make sure there's at least one neighbor
       if ( verlet_finish .ge. verlet_start ) then
          do i_index= verlet_start, verlet_finish    
             j_atom = verlet_neighbor_list_lj(i_index)
             atom_id2 = atom_list_atom_index(j_atom)

             rij = atom_list_xyz(i_atom,:) - atom_list_xyz(j_atom,:)

             ! shift for general box
!!$             dr_direct(:) = matmul( xyz_to_box_transform, rij )
!!$             do i=1,3
!!$                shift_direct(i) = dble(floor( dr_direct(i) + 0.5d0 ))
!!$             enddo
!!$             shift = matmul( shift_direct , box )
             ! shift for orthorhombic box
             do i=1,3
                shift(i) = box(i,i) * floor( rij(i) / box(i,i) + 0.5d0 )
             end do

             rij = rij - shift
             norm_dr2 = dot_product( rij, rij )              

             if ( norm_dr2 < lj_cutoff2 ) then
                norm_dr6 = norm_dr2 ** 3
                norm_dr12 = norm_dr6 ** 2
                ! at this stage, all lj parameters should be expressed as C12 and C6, even though they were read in as epsilon and sigma
                term12 = atype_lj_parameter(atom_id1,atom_id2,1) / norm_dr12
                term6 = atype_lj_parameter(atom_id1,atom_id2,2) / norm_dr6
                lj_energy = lj_energy + term12 - term6 - lj_shift_store(atom_id1,atom_id2)
                f_ij = rij / norm_dr2 * ( 12d0 * term12  - 6d0 * term6 )
                atom_list_force(i_atom,:) = atom_list_force(i_atom,:) + f_ij(:)
                atom_list_force(j_atom,:) = atom_list_force(j_atom,:) - f_ij(:)
             end if
          end do
       endif
    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL



  end  subroutine lennard_jones_use_verlet





  !***********************************************************************
  ! self-explanatory
  !***********************************************************************
  subroutine lennard_jones_no_verlet( lj_energy, lj_force, tot_n_mole, n_mole, n_atom, xyz, box, r_com)
    use omp_lib
    real*8,intent(out) :: lj_energy
    real*8,dimension(:,:,:),intent(out)::lj_force
    integer,intent(in)::n_mole,tot_n_mole
    integer,dimension(:),intent(in)::n_atom
    real*8,dimension(:,:,:),intent(in)::xyz
    real*8,dimension(:,:),intent(in)::box
    real*8,dimension(:,:),intent(in)::r_com
    integer::i_mole,j_mole,i_atom,j_atom,atom_id1,atom_id2,i
    integer :: split_do
    real*8 :: lj_cutoff2, norm_dr2, norm_dr6, norm_dr12, term6, term12
    real*8,dimension(3)::f_ij,rij,shift,shift_direct, dr_com,dr_direct

    lj_force=0d0
    lj_energy=0d0

    lj_cutoff2 = lj_cutoff ** 2

    ! decide how to split the parallel section
    if (n_threads .eq. 1 ) then
       split_do = 1
    else
       split_do = n_mole/n_threads+1
    endif


    !**************************************** no Verlet list, naive looping ***********************************************
    call OMP_SET_NUM_THREADS(n_threads)

    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(tot_n_mole,n_mole,r_com,box,n_atom,xyz,atom_index,atype_lj_parameter,atype_penetration_ff,lj_cutoff2,split_do,atype_solute,C8_10_dispersion_terms,xyz_to_box_transform,penetration_force_field,Feynmann_Hibbs_forces,atype_mass,temp,lj_bkghm_index,lj_shift_store) REDUCTION(+:lj_force, lj_energy)
    !$OMP DO SCHEDULE(DYNAMIC,split_do)
    do i_mole=1,n_mole
       do j_mole=i_mole+1,tot_n_mole
          do i_atom=1,n_atom(i_mole)
             do j_atom=1,n_atom(j_mole)
                atom_id1 = atom_index(i_mole,i_atom)
                atom_id2 = atom_index(j_mole,j_atom)
                rij = xyz(i_mole,i_atom,:) - xyz(j_mole,j_atom,:)

                dr_direct(:) = matmul( xyz_to_box_transform, rij )
                do i=1,3
                   shift_direct(i) = dble(floor( dr_direct(i) + 0.5d0 ))
                enddo
                shift = matmul( shift_direct , box )
                rij = rij - shift
                norm_dr2 = dot_product( rij, rij )  
                if ( norm_dr2 < lj_cutoff2 ) then          
                   norm_dr6 = norm_dr2 ** 3
                   norm_dr12 = norm_dr6 ** 2
                   ! at this stage, all lj parameters should be expressed as C12 and C6, even though they were read in as epsilon and sigma
                   term12 = atype_lj_parameter(atom_id1,atom_id2,1) / norm_dr12
                   term6 = atype_lj_parameter(atom_id1,atom_id2,2) / norm_dr6
                   lj_energy = lj_energy + term12 - term6
                   f_ij = rij / norm_dr2 * ( 12d0 * term12  - 6d0 * term6 )
                   lj_force(i_mole,i_atom,:) = lj_force(i_mole,i_atom,:) + f_ij(:)
                   lj_force(j_mole,j_atom,:) = lj_force(j_mole,j_atom,:) - f_ij(:)

                   ! only add shift if atom pair is closer than cutoff
                   if ( norm_dr2 < lj_cutoff2 ) then
                      lj_energy = lj_energy - lj_shift_store(atom_id1,atom_id2)
                   end if
                end if
             end do
          end do
       enddo
    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL


  end  subroutine lennard_jones_no_verlet






  !*******************************************
  ! this subroutine computes the intra-molecular lennard jones force and energy contribution
  !*******************************************
  subroutine intra_lennard_jones_energy_force( E_intra_lj , lj_force, i_mole, n_atom, xyz, lj_cutoff2 )
    use global_variables
    real*8 , intent(out) :: E_intra_lj
    real*8 , dimension(:,:,:),intent(inout) :: lj_force
    integer, intent(in) :: i_mole
    integer, dimension(:), intent(in) :: n_atom
    real*8, dimension(:,:,:), intent(in) :: xyz
    real*8, intent(in) :: lj_cutoff2

    integer :: i_atom, j_atom, i_mole_type, atom_id1, atom_id2
    real*8, dimension(3) :: rij, f_ij
    real*8 :: norm_dr2, norm_dr6, norm_dr12, term12, term6, C6, C12

    E_intra_lj = 0d0

    i_mole_type = molecule_index(i_mole)

    do i_atom=1, n_atom(i_mole)
       do j_atom = i_atom + 1, n_atom(i_mole)
          ! check for exclusions 
          if ( molecule_exclusions( i_mole_type, i_atom, j_atom ) /= 1 ) then
             ! not excluded
             ! as always, molecules are not broken up over pbc, so no need to consider shift

             rij = xyz(i_mole,i_atom,:) - xyz(i_mole,j_atom,:)
             norm_dr2 = dot_product( rij, rij )              
             if ( norm_dr2 < lj_cutoff2 ) then

                norm_dr6 = norm_dr2 ** 3
                norm_dr12 = norm_dr6 ** 2

                atom_id1 = atom_index(i_mole,i_atom)
                atom_id2 = atom_index(i_mole,j_atom)

                ! at this stage, all lj parameters should be expressed as C12 and C6, even though they were read in as epsilon and sigma
                ! if this is a 1-4 interaction, take parameters from 1-4 interaction parameter array
                if ( molecule_exclusions( i_mole_type, i_atom , j_atom ) == 2 ) then
                   C12 = atype_lj_parameter_14(atom_id1,atom_id2,1)
                   C6  = atype_lj_parameter_14(atom_id1,atom_id2,2)
                else
                   C12 = atype_lj_parameter(atom_id1,atom_id2,1)
                   C6  = atype_lj_parameter(atom_id1,atom_id2,2)
                end if

                term12 = C12 / norm_dr12
                term6 =  C6 / norm_dr6
                ! don't consider shift here
                E_intra_lj = E_intra_lj + term12 - term6
                f_ij = rij / norm_dr2 * ( 12d0 * term12  - 6d0 * term6 )
                lj_force(i_mole,i_atom,:) = lj_force(i_mole,i_atom,:) + f_ij(:)
                lj_force(i_mole,j_atom,:) = lj_force(i_mole,j_atom,:) - f_ij(:)


             end if
          end if
       enddo
    enddo


  end subroutine intra_lennard_jones_energy_force




  !***********************************************************************
  ! This subroutine computes lennard jones (buckingham) energy of particle insertion,
  ! only computing the energy between the tagged molecule newmol, and the rest of the atoms
  ! in the system.  This can be used for either particle-insertions/deletions, or MC moves
  ! where only a single molecule is displaced
  !***********************************************************************
  subroutine update_lennard_jones_ins( energy,energy_decomp,tot_n_mole, n_atom, xyz, newmol, box, r_com )
    integer,intent(in)::tot_n_mole,newmol
    integer,dimension(:),intent(in)::n_atom
    real*8,dimension(:,:,:),intent(in)::xyz
    real*8,dimension(:,:),intent(in)::box
    real*8,intent(out)::energy
    real*8,dimension(:),intent(out)::energy_decomp
    real*8,dimension(:,:),intent(in)::r_com
    integer::j_mole
    real*8::sum_short,sum_long,Eij_short,Eij_long,Eij_FH,E_FH
    real*8,dimension(size(energy_decomp)) :: energy_decomp_ij


    sum_short=0D0;sum_long=0D0;energy_decomp=0D0;E_FH=0d0

    ! since decompositon types are hardwired in here for 5 total types, make sure this is correct
    if (size(energy_decomp) .ne. 5) then
       stop " update_lennard_jones_ins is implemented for 5 energy decomposition terms"
    endif

    do j_mole=1,tot_n_mole
       if(newmol.ne.j_mole) then
          call lennard_jones_interaction_energy( Eij_short,Eij_long, Eij_FH, energy_decomp_ij,newmol, j_mole, n_atom, xyz, box, r_com )
          sum_short = sum_short + Eij_short
          sum_long = sum_long + Eij_long
          energy_decomp = energy_decomp + energy_decomp_ij
          E_FH = E_FH + Eij_FH
       endif
    enddo

    energy = sum_long + sum_short
    Select Case(Feynmann_Hibbs_correction)
    Case("yes")
       energy=energy+E_FH
    End Select

  end subroutine update_lennard_jones_ins




  !***********************************************************************
  ! This subroutine computes lennard jones (buckingham) interactions
  ! It is implemented for generic solute, as atom index array stores the type
  ! of atom for each molecule, which is then used to recover LJ parameters in atype_lj_parameter array
  ! also, storing framework atoms in tot_n_mole - n_mole indices allows solute framework interactions to
  ! be calculated
  ! if energy shift is turned on, lj_shift_store will hold these shifts.
  ! if energy shift is off, lj_shift_store contains zeros and adds nothing
  !
  ! this has been updated to include calculation of energy_decomposition
  ! indices of energy_decomposition are contributions from exchange, electrostatic, induction, dhf, and dispersion respectively
  !***********************************************************************
  subroutine lennard_jones_interaction_energy( Eij_short,Eij_long, Eij_FH, energy_decomp,i_mole, j_mole, n_atom, xyz, box, r_com )
    integer,intent(in):: i_mole, j_mole
    integer,dimension(:),intent(in)::n_atom
    real*8,dimension(:,:,:),intent(in)::xyz
    real*8,dimension(:,:),intent(in)::box
    real*8, dimension(:),intent(out) :: energy_decomp
    real*8,intent(out):: Eij_short, Eij_long, Eij_FH
    real*8,dimension(:,:),intent(in)::r_com
    integer::i_atom,j_atom,atom_id1,atom_id2
    real*8::sum_short,sum_long,r_ij,dist,lj_cutoff_use,Exponential,E_short,E_long,r6inv,r12inv,r2,E_local,r_mag
    real*8,dimension(3)::rij,shift,dr_com,r_com_i,r_com_j
    real*8 :: E_ms_evb_repulsion

    !lj_cutoff_use = min(box(1,1)/2.,box(2,2)/2.,box(3,3)/2.,lj_cutoff)
    lj_cutoff_use = lj_cutoff

    sum_short=0D0;sum_long=0D0;energy_decomp=0D0;Eij_FH=0d0

    r_com_i(:) = r_com(i_mole,:)
    r_com_j(:) = r_com(j_mole,:)
    shift = pbc_shift( r_com_i, r_com_j, box, xyz_to_box_transform)
    dr_com = pbc_dr( r_com_i, r_com_j, shift)


    if((abs(dr_com(1))<lj_cutoff_use).and.(abs(dr_com(2))<lj_cutoff_use).and.(abs(dr_com(3))<lj_cutoff_use)) then
       if ( dot_product( dr_com, dr_com ) < lj_cutoff_use**2 ) then

          do i_atom=1,n_atom(i_mole)
             do j_atom=1,n_atom(j_mole)
                rij = pbc_dr( xyz(i_mole,i_atom,:), xyz(j_mole,j_atom,:), shift ) ! Note for COM cutoff, shift values unchanged

                r2 = dot_product( rij, rij )
                r_ij= dsqrt(r2)

                atom_id1 = atom_index(i_mole,i_atom)
                atom_id2 = atom_index(j_mole,j_atom)

                ! now figure out whether interaction is lj or bkghm
                Select Case(lj_bkghm_index(atom_id1,atom_id2))
                Case(1)
                   !**********************************buckingham potential******************************************!
                   ! get long range component
                   E_long = eval_C6_10_dispersion(atom_id1,atom_id2,r_ij)
                   ! get short range component
                   call eval_short_range_bkghm(E_short, energy_decomp, atom_id1,atom_id2, r_ij)

                   Select Case (energy_decomposition)
                      !*****************no energy decompositon, one buckingham term************!
                   Case("no")
                      ! don't shift potential to values at cutoff (for molecule close to cutoff)
                      if( r_ij > lj_cutoff_use) then
                         sum_long = sum_long + min (E_short + E_long  - lj_shift_store(atom_id1,atom_id2) , 0d0)
                      else
                         sum_short=sum_short + E_short
                         sum_long=sum_long + E_long - lj_shift_store(atom_id1,atom_id2)
                      endif
                      !*****************energy decompositon, five buckingham terms************!
                   Case("yes")
                      ! don't shift potential to values at cutoff (for molecule close to cutoff)
                      if( r_ij > lj_cutoff_use) then
                         ! notice we might be overcounting E_short here in the energy_decomposition if lj_cutoff_use is so short that
                         ! the buckingham terms haven't gone to zero yet, but this is minor and it is not overcounted in the total potential
                         ! as we are not updating sum_short, so it doesn't effect the simulation
                         energy_decomp(5) = energy_decomp(5) + min (E_short + E_long  - lj_shift_store(atom_id1,atom_id2) , 0d0)
                         sum_long = sum_long + min (E_short + E_long  - lj_shift_store(atom_id1,atom_id2) , 0d0)
                      else
                         sum_short=sum_short + E_short
                         energy_decomp(5) = energy_decomp(5) + E_long - lj_shift_store(atom_id1,atom_id2)
                         sum_long=sum_long + E_long - lj_shift_store(atom_id1,atom_id2)
                      endif
                   end Select
                   !************************************************************************!

                Case(2)
                   !**********************************lennard jones potential******************************************!
                   r6inv = r2**(-3)
                   r12inv = r6inv**2
                   ! at this stage, all lj parameters should be expressed as C12 and C6, even though they were read in as epsilon and sigma
                   E_short = atype_lj_parameter(atom_id1,atom_id2,1) * r12inv
                   E_long =  - atype_lj_parameter(atom_id1,atom_id2,2) * r6inv
                   ! don't shift potential to values at cutoff (for molecule close to cutoff)
                   if( r_ij > lj_cutoff_use) then
                      sum_long = sum_long + min (E_short + E_long  - lj_shift_store(atom_id1,atom_id2) , 0d0)
                   else
                      sum_short=sum_short + E_short
                      sum_long=sum_long + E_long - lj_shift_store(atom_id1,atom_id2)
                   endif

                End Select


                Select Case(Feynmann_Hibbs_correction)
                Case("yes")
                   ! if adding Feynmann_Hibbs correction
                   call Feynmann_Hibbs_lj_correction(E_local, r_ij , atom_id1, atom_id2,i_mole,j_mole,n_atom )
                   Eij_FH =Eij_FH + E_local
                End Select
             enddo
          enddo

       endif
    endif

    Eij_short = sum_short
    Eij_long = sum_long


  end subroutine lennard_jones_interaction_energy







  !***********************************************************************
  ! This subroutine computes lennard jones (buckingham) forces for the molecule pair i_mole, j_mole
  ! It is implemented for generic solute, as atom index array stores the type
  ! of atom for each molecule, which is then used to recover LJ parameters in atype_lj_parameter array
  !***********************************************************************
  subroutine F_lennard_jones_interaction( lj_force_ij, i_mole, j_mole, n_atom, xyz, box, r_com)
    real*8,dimension(:,:,:),intent(out)::lj_force_ij
    integer,intent(in)::i_mole,j_mole
    integer,dimension(:),intent(in)::n_atom
    real*8,dimension(:,:,:),intent(in)::xyz
    real*8,dimension(:,:),intent(in)::box
    real*8,dimension(:,:),intent(in)::r_com
    integer::i_atom,j_atom,atom_id1,atom_id2
    real*8::norm_dr,lj_cutoff_use,Exponent,r6inv, r12inv,r2
    real*8 :: small_thresh=1D-4
    real*8,dimension(3)::r_com_i,r_com_j
    real*8,dimension(3)::shift,dr_com,r_ij,f_ij,f_ij_lr

    lj_force_ij=0d0
    ! lj_cutoff_use = min(box(1,1)/2.,box(2,2)/2.,box(3,3)/2.,lj_cutoff)
    lj_cutoff_use = lj_cutoff


    r_com_i(:) = r_com( i_mole, : ) 
    r_com_j(:) = r_com( j_mole, : )

    shift = pbc_shift( r_com_i, r_com_j, box , xyz_to_box_transform )
    dr_com = pbc_dr( r_com_i, r_com_j, shift )
    if((abs(dr_com(1))<lj_cutoff_use).and.(abs(dr_com(2))<lj_cutoff_use).and.(abs(dr_com(3))<lj_cutoff_use)) then
       if ( dot_product( dr_com, dr_com ) < lj_cutoff_use**2 ) then
          do i_atom=1,n_atom(i_mole)
             do j_atom=1,n_atom(j_mole)
                r_ij = -pbc_dr( xyz(i_mole,i_atom,:), xyz(j_mole,j_atom,:), shift ) ! Note for COM cutoff, shift values unchanged
                r2=dot_product(r_ij,r_ij)
                norm_dr = sqrt(r2)
                atom_id1 = atom_index(i_mole,i_atom)
                atom_id2 = atom_index(j_mole,j_atom)

                ! now figure out whether interaction is lj or bkghm
                Select Case(lj_bkghm_index(atom_id1,atom_id2))
                Case(1)
                   !**********************************buckingham potential******************************************!
                   f_ij = eval_short_range_force(atom_id1,atom_id2,r_ij,norm_dr)
                   f_ij_lr = eval_C6_10_force(atom_id1,atom_id2,r_ij,norm_dr)
                   f_ij = f_ij + f_ij_lr

                Case(2)
                   !**********************************lennard jones potential******************************************!
                   ! don't waste time if no interaction
                   if ( (atype_lj_parameter(atom_id1,atom_id2,1) < small_thresh) .and. (atype_lj_parameter(atom_id1,atom_id2,2) < small_thresh)  ) then
                      f_ij =0d0
                   else

                      r6inv = r2**(-3)
                      r12inv = r6inv**2
                      ! at this stage, all lj parameters should be expressed as C12 and C6, even though they were read in as epsilon and sigma
                      f_ij = atype_lj_parameter(atom_id1,atom_id2,1) * 12d0 * r12inv * r_ij/r2
                      f_ij = f_ij - atype_lj_parameter(atom_id1,atom_id2,2) * 6d0 * r6inv * r_ij/r2
                   endif

                End Select

                ! store forces for i_mole in 1st index of lj_force_ij array, and forces for j_mole in 2nd index
                lj_force_ij(1,i_atom,:)=lj_force_ij(1,i_atom,:) + f_ij(:)
                lj_force_ij(2,j_atom,:)=lj_force_ij(2,j_atom,:) - f_ij(:)

                Select Case(Feynmann_Hibbs_forces)
                Case("yes")
                   ! if adding Feynmann_Hibbs forces
                   call Feynmann_Hibbs_lj_forces(f_ij, r_ij , atom_id1, atom_id2,i_mole,j_mole,n_atom )
                   lj_force_ij(1,i_atom,:)=lj_force_ij(1,i_atom,:) + f_ij(:)
                   lj_force_ij(2,j_atom,:)=lj_force_ij(2,j_atom,:) - f_ij(:)                               
                End Select

             enddo
          enddo


       endif
    end if


  end subroutine F_lennard_jones_interaction





  !*************************************************
  ! this subroutine evaluates the short range part of a buckingham pair wise potential
  ! the main reason for a separate routine here is if a separate penetration force field is 
  ! being used.  Note there is no select case for penetration_force_field.  Rather,
  ! if a penetration force field is not used, atype_penetration(:,:,2) will be set to zero
  ! for all atoms,which is sufficient
  !*************************************************
  subroutine eval_short_range_bkghm(E_short, energy_decomp,atom_id1,atom_id2,r_ij)
    integer,intent(in) :: atom_id1, atom_id2
    real*8,intent(in) :: r_ij
    real*8,intent(out) :: E_short
    real*8,dimension(:),intent(inout) :: energy_decomp

    real*8 :: Exponential,sum, energy_test=1d-8
    integer :: i_type

    ! determine whether pairwise distance is beyond cutoff
    if ( r_ij > atype_penetration_ff(atom_id1,atom_id2,2) ) then
       ! decomposed force field, check energy decomposition
       Select Case(energy_decomposition)
       Case("no")
          E_short =  atype_lj_parameter(atom_id1,atom_id2,1)*exp(-atype_lj_parameter(atom_id1,atom_id2,2)*r_ij)
       Case("yes")
          Exponential = exp(-atype_lj_parameter(atom_id1,atom_id2,2)*r_ij)
          sum=0d0
          do i_type =1, 5
             energy_decomp(i_type) = energy_decomp(i_type) + atype_bkghm_decomp(atom_id1,atom_id2,i_type) * Exponential
             sum = sum + atype_bkghm_decomp(atom_id1,atom_id2,i_type) * Exponential
          enddo
          E_short = atype_lj_parameter(atom_id1,atom_id2,1) * Exponential
          ! check consistency of energy decomp
          if ( abs(sum - E_short) > energy_test ) then
             write(*,*) "short range energy decomposition doesn't equal E_short for atom types ", atom_id1, atom_id2
             stop
          endif
       End Select
    else
       ! penetration force field, no contribution to energy decomposition
       E_short = atype_penetration_ff(atom_id1,atom_id2,1) * exp(-atype_penetration_ff(atom_id1,atom_id2,4) * atype_lj_parameter(atom_id1,atom_id2,2)*r_ij) - atype_penetration_ff(atom_id1,atom_id2,3)

    endif

  end subroutine eval_short_range_bkghm



  !*************************************************
  ! this function evaluates the short range buckingham force
  ! Note there is no select case for penetration_force_field.  Rather,
  ! if a penetration force field is not used, atype_penetration(:,:,2) will be set to zero
  ! for all atoms,which is sufficient
  !*************************************************
  function eval_short_range_force(atom_id1,atom_id2,r_ij,norm_dr)
    real*8,dimension(3) :: eval_short_range_force
    integer,intent(in) :: atom_id1,atom_id2
    real*8,dimension(3),intent(in) :: r_ij
    real*8,intent(in) :: norm_dr

    ! determine whether pairwise distance is beyond cutoff
    if ( norm_dr > atype_penetration_ff(atom_id1,atom_id2,2) ) then
       ! decomposed force field
       eval_short_range_force = atype_lj_parameter(atom_id1,atom_id2,1)*exp(-atype_lj_parameter(atom_id1,atom_id2,2)*norm_dr)*(atype_lj_parameter(atom_id1,atom_id2,2))*r_ij/norm_dr

    else
       ! penetration force field
       eval_short_range_force = atype_penetration_ff(atom_id1,atom_id2,1)*exp(-atype_penetration_ff(atom_id1,atom_id2,4) * atype_lj_parameter(atom_id1,atom_id2,2)*norm_dr)*(atype_penetration_ff(atom_id1,atom_id2,4) * atype_lj_parameter(atom_id1,atom_id2,2))*r_ij/norm_dr      

    endif

  end function eval_short_range_force


  !**************************************************
  ! this function returns long range dispersion energy for two input atoms
  ! if C8,C10 terms are used, all long range components are damped
  ! Note there is no select case for penetration_force_field.  Rather,
  ! if a penetration force field is not used, atype_penetration(:,:,2) will be set to zero
  ! for all atoms,which is sufficient
  !**************************************************
  function eval_C6_10_dispersion(atom_id1,atom_id2,r_ij)
    real*8 :: eval_C6_10_dispersion
    integer,intent(in) :: atom_id1,atom_id2
    real*8,intent(in) :: r_ij

    real*8 :: r_ij2,r_ij6,r_ij8,r_ij10,r_ij12,term6,term8,term10,term12

    ! determine whether pairwise distance is beyond cutoff
    if ( r_ij > atype_penetration_ff(atom_id1,atom_id2,2) ) then

       r_ij2 = r_ij**2
       r_ij6 = r_ij2**3
       ! determine whether we are using C8,C10 terms
       Select Case(C8_10_dispersion_terms)
       Case("no")
          eval_C6_10_dispersion= -atype_lj_parameter(atom_id1,atom_id2,3)/r_ij6
       Case("yes")
          r_ij8=r_ij6*r_ij2
          r_ij10=r_ij8*r_ij2
          term6= - C6_C10_damp(atom_id1,atom_id2,r_ij,6) * atype_lj_parameter(atom_id1,atom_id2,3)/r_ij6
          term8= - C6_C10_damp(atom_id1,atom_id2,r_ij,8) * atype_lj_parameter(atom_id1,atom_id2,4)/r_ij8
          term10= - C6_C10_damp(atom_id1,atom_id2,r_ij,10) * atype_lj_parameter(atom_id1,atom_id2,5)/r_ij10

          eval_C6_10_dispersion= term6  + term8 + term10

          ! add C12 dispersion if requested
          Select Case(C12_dispersion)
          Case("yes")
             r_ij12=r_ij10*r_ij2
             term12= - C6_C10_damp(atom_id1,atom_id2,r_ij,12) * atype_lj_parameter(atom_id1,atom_id2,6)/r_ij12

             eval_C6_10_dispersion=eval_C6_10_dispersion + term12
          End Select



       end select

    else
       eval_C6_10_dispersion=0d0

    endif

  end function eval_C6_10_dispersion


  !**********************************************
  ! this function returns long range contribution to the pairwise dispersion force.
  ! if C8,C10 terms are included, these derivative contribute as well
  ! Note there is no select case for penetration_force_field.  Rather,
  ! if a penetration force field is not used, atype_penetration(:,:,2) will be set to zero
  ! for all atoms,which is sufficient
  !**********************************************
  function eval_C6_10_force(atom_id1,atom_id2,r_vec,r_ij)
    real*8,dimension(3) :: eval_C6_10_force
    integer,intent(in) :: atom_id1,atom_id2
    real*8,dimension(3),intent(in) :: r_vec
    real*8,intent(in) :: r_ij

    real*8 :: r_ij2,r_ij6,r_ij8,r_ij10,r_ij12, fac6,fac8,fac10,fac12
    real*8,dimension(3) :: term6,term8,term10,term12

    ! determine whether pairwise distance is beyond cutoff
    if ( r_ij > atype_penetration_ff(atom_id1,atom_id2,2) ) then

       r_ij2 = r_ij**2
       r_ij6 = r_ij2**3
       ! determine whether we are using C8,C10 terms
       Select Case(C8_10_dispersion_terms)
       Case("no")
          eval_C6_10_force = -(atype_lj_parameter(atom_id1,atom_id2,3)/r_ij6)*(6d0 *r_vec/r_ij2)
       Case("yes")
          r_ij8=r_ij6*r_ij2
          r_ij10=r_ij8*r_ij2       
          fac6= atype_lj_parameter(atom_id1,atom_id2,3)/r_ij6
          fac8= atype_lj_parameter(atom_id1,atom_id2,4)/r_ij8
          fac10= atype_lj_parameter(atom_id1,atom_id2,5)/r_ij10

          term6 = -C6_C10_damp(atom_id1,atom_id2,r_ij,6) * fac6 * (6d0 *r_vec/r_ij2) + dC6_C10_damp(atom_id1,atom_id2,r_vec,r_ij,6) * fac6
          term8 = -C6_C10_damp(atom_id1,atom_id2,r_ij,8) * fac8 * (8d0 *r_vec/r_ij2) + dC6_C10_damp(atom_id1,atom_id2,r_vec,r_ij,8) * fac8
          term10 = -C6_C10_damp(atom_id1,atom_id2,r_ij,10) * fac10 * (10d0 *r_vec/r_ij2) + dC6_C10_damp(atom_id1,atom_id2,r_vec,r_ij,10) * fac10

          eval_C6_10_force = term6 + term8 + term10

          ! add C12 dispersion if requested
          Select Case(C12_dispersion)
          Case("yes")
             r_ij12=r_ij10*r_ij2 
             fac12= atype_lj_parameter(atom_id1,atom_id2,6)/r_ij12
             term12 = -C6_C10_damp(atom_id1,atom_id2,r_ij,12) * fac12 * (12d0 *r_vec/r_ij2) + dC6_C10_damp(atom_id1,atom_id2,r_vec,r_ij,12) * fac12

             eval_C6_10_force = eval_C6_10_force + term12
          End Select



       end select

    else
       eval_C6_10_force=0d0

    endif

  end function eval_C6_10_force


  !***********************************************************************
  !  this subroutine stores shift values for pairwise interactions
  !  in which the Lennard-Jones/Buckingham part of the potential is shifted
  !  such that it has a value of zero at the cutoff
  !***********************************************************************

  subroutine update_lennard_jones_shift(box)
    real*8,dimension(:,:),intent(in)::box

    integer:: i,j
    real*8:: r_ij

    ! r_ij = min(box(1,1)/2.,box(2,2)/2.,box(3,3)/2.,lj_cutoff)
    r_ij = lj_cutoff

    Select case(lj_shift)
    case(0)
       lj_shift_store = 0D0
    case(1)
       ! atom types will be created with expanded grand canonical ensemble
       if ( select_ensemble .eq. "egc" ) then
          stop " can't use shift with egc ensemble as atom types will be created"
       endif

       if ( lj_bkghm .eq. 1 )  then
          !! shift buckingham potential
          do i=1,n_atom_type
             do j=1,n_atom_type
                Select Case(C8_10_dispersion_terms)
                Case("no")
                   lj_shift_store(i,j) = atype_lj_parameter(i,j,1)*exp(-atype_lj_parameter(i,j,2)*r_ij) - atype_lj_parameter(i,j,3)/r_ij**6
                Case("yes")
                   ! assume that cutoff is long enough so that we dont have to worry about exponential part of potential or damping functions
                   lj_shift_store(i,j) =  - atype_lj_parameter(i,j,3)/r_ij**6 -  atype_lj_parameter(i,j,4)/r_ij**8 -  atype_lj_parameter(i,j,5)/r_ij**10
                   Select Case(C12_dispersion)
                   Case("yes")
                      lj_shift_store(i,j) = lj_shift_store(i,j) -  atype_lj_parameter(i,j,6)/r_ij**12
                   End Select
                End Select
             enddo
          enddo
       elseif ( lj_bkghm .eq. 2 )  then
          !! shift lennard jones potential
          do i=1,n_atom_type
             do j=1,n_atom_type
                ! at this point , lennard jones parameters are expressed as C12,C6, even though they were read in as epsilon and sigma
                lj_shift_store(i,j) =  -atype_lj_parameter(i,j,2)/ r_ij**6
             enddo
          enddo
       endif
    case default
       stop "lj_shift selection not recognized (should be 0 or 1)"
    end select


  end subroutine update_lennard_jones_shift



  !***********************************************************************
  !  this subroutine evaluates the long range correction to dispersion
  !  for the case where there is no shift
  !  if there is a framework present, this will automatically be set to zero
  !  Repulsive interactions are neglected for the lrc since these easily die off for 
  !  any reasonable cutoff
  !  also, screening functions for C6,C8,C10 are neglected for the lrc,
  !  since these are effectively 1 for exponents/damping parameters(~3.2 A-1)
  !  and any reasonable cutoff close to or greater than 10 Angstrom
  !
  !  C12 terms will not be included in long range correction
  !
  !  optional argument flag=1 tells subroutine to update global variable E_disp_lrc_try instead
  !  of E_disp_lrc
  !***********************************************************************
  subroutine update_disp_lrc_noshift(n_mole, n_atom, alist, box, flag) 
    integer,intent(in)::n_mole
    integer,dimension(:),intent(in)::n_atom
    character(*),dimension(:,:),intent(in) :: alist
    real*8,dimension(:,:),intent(in)::box
    integer,intent(in),optional :: flag
    real*8::vol,lj_cutoff_use, E_disp_lrc_internal
    integer::i,j,i_mole,i_atom,i_type,j_type,ind
    integer,save::initialize=0,flag_warn=0
    integer, dimension(MAX_N_ATOM),save :: n_atoms_type

    !******************** initialize n_atoms_type array ************
    Select Case(initialize)
    Case(0)
       n_atoms_type=0
       do i_mole=1, n_mole
          do i_atom=1,n_atom(i_mole)

             do  i_type=1, n_atom_type
                if (atype_name(i_type) .eq. alist(i_mole,i_atom) ) then
                   n_atoms_type(i_type) = n_atoms_type(i_type) + 1
                   exit 
                endif
             enddo
          enddo
       enddo
       initialize=1
    End Select

    E_disp_lrc_internal = 0d0

    Select Case ( lj_lrc)
    Case(1)
       ! only long range correction if not a gcmc framework simulation
       if(framework_simulation .eq. 0) then

          ! warn if C12 dispersion is on
          Select Case( flag_warn )
          Case(0)
             Select Case( C12_dispersion )
             Case("yes")
                write(*,*) "NOTE: The dispersion long range correction will contain no contribution"
                write(*,*) "from C12 terms"
                write(*,*) ""
             End Select
             flag_warn=1
          End Select

          vol = volume( box(1,:) , box(2,:) , box(3,:) )
          !lj_cutoff_use = min(box(1,1)/2.,box(2,2)/2.,box(3,3)/2.,lj_cutoff)
          lj_cutoff_use = lj_cutoff

          !******************** calculate contribution to long range correction for each pair of atom types

          do i_type=1, n_atom_type
             do j_type=i_type,n_atom_type


                if( lj_bkghm .eq. 1 ) then
                   !**************************buckingham potential*************************************!
                   ! determine whether we are using C8,C10 terms
                   Select Case(C8_10_dispersion_terms)
                   Case("no")
                      ! if i_type equals j_type, we need to correct for overcounting
                      if ( i_type == j_type ) then
                         E_disp_lrc_internal = E_disp_lrc_internal -(2./3.)*pi*n_atoms_type(i_type) * n_atoms_type(j_type) /vol* atype_lj_parameter(i_type,j_type,3)/(lj_cutoff_use**3)
                      else
                         E_disp_lrc_internal = E_disp_lrc_internal -(4./3.)*pi*n_atoms_type(i_type) * n_atoms_type(j_type) /vol* atype_lj_parameter(i_type,j_type,3)/(lj_cutoff_use**3)
                      endif
                   Case("yes")
                      ! if i_type equals j_type, we need to correct for overcounting
                      if ( i_type == j_type ) then
                         E_disp_lrc_internal = E_disp_lrc_internal - pi*n_atoms_type(i_type) * n_atoms_type(j_type) /vol* ( ( 2./3.) * atype_lj_parameter(i_type,j_type,3)/(lj_cutoff_use**3) + (2./5.) * atype_lj_parameter(i_type,j_type,4)/(lj_cutoff_use**5) + (2./7.) * atype_lj_parameter(i_type,j_type,5)/(lj_cutoff_use**7) )
                      else
                         E_disp_lrc_internal = E_disp_lrc_internal - 2d0 * pi*n_atoms_type(i_type) * n_atoms_type(j_type) /vol* ( ( 2./3.) * atype_lj_parameter(i_type,j_type,3)/(lj_cutoff_use**3) + (2./5.) * atype_lj_parameter(i_type,j_type,4)/(lj_cutoff_use**5) + (2./7.) * atype_lj_parameter(i_type,j_type,5)/(lj_cutoff_use**7) )
                      endif
                   End Select

                elseif( lj_bkghm .eq. 2 ) then
                   !***********************lennard jones potential**************************************!
                   ! at this point, lennard jones parameters are expressed as C12,C6, even though they were read in as epsilon and sigma
                   ! if i_type equals j_type, we need to correct for overcounting
                   if ( i_type == j_type ) then
                      E_disp_lrc_internal = E_disp_lrc_internal -(2./3.)*pi * n_atoms_type(i_type) * n_atoms_type(j_type) /vol * atype_lj_parameter(i_type,j_type,2) /(lj_cutoff_use**3)
                   else
                      E_disp_lrc_internal = E_disp_lrc_internal -(4./3.)*pi * n_atoms_type(i_type) * n_atoms_type(j_type) /vol * atype_lj_parameter(i_type,j_type,2) /(lj_cutoff_use**3)
                   endif

                endif

             enddo
          enddo

       endif

    end select

    if ( present(flag) ) then
       if ( flag == 1 ) then
          E_disp_lrc_try = E_disp_lrc_internal
       else
          E_disp_lrc = E_disp_lrc_internal  
       endif
    else
       E_disp_lrc = E_disp_lrc_internal
    end if


  end subroutine update_disp_lrc_noshift



  !************************************************
  ! this subroutine evaluates the Feynmann_Hibbs quantum correction
  ! for either a lj or bkghm potential.
  ! for dummy atoms with zero mass, located at the center of the molecule, the correct treatment 
  ! is to use the mass of the molecule for the reduced mass.  This is what is done in 
  ! generate_reduced_mass subroutine
  !
  !  Feynmann Hibbs effective potential for pairwise Uij is
  !  Uij(eff) = Uij + (hbar^2 ) / (24 *kT * mass_red ) * ( d2Uij/dr2 + (2/r) * dU/dr )
  !  
  !  see 
  !  JPC B, 2008, 112, 11708
  !  and
  !  Feynmann, Statistical Mechanics: A Set of Lectures (Equation 3.81 )
  !************************************************
  subroutine Feynmann_Hibbs_lj_correction(E_local, r_mag , atom_id1, atom_id2,i_mole,j_mole,n_atom )
    real*8, intent(out) :: E_local
    real*8,intent(in) :: r_mag
    integer,intent(in) :: atom_id1, atom_id2,i_mole,j_mole
    integer,dimension(:),intent(in) :: n_atom

    real*8  :: A,B,C,F,dF,d2F,C12, C6,fac, mass_reduced,r2,r2inv, r6inv, r8inv, r14inv,r_ij(14), FH_short,FH_long
    real*8,parameter :: hbar_mole_sq=0.403280544d0  ! this is (hbar*NA)^2 ( in units of Kg^2 m^2*A^2/s^2 )
    integer :: n_terms,i,C_index,R_index,R_index1,R_index2

    call generate_reduced_mass(mass_reduced, atom_id1, atom_id2,i_mole,j_mole,n_atom)

    ! fac = hbar^2/(24*kT * mass_red ) ; since distance is in angstroms, this should be in A^2
    ! use (hbar*mol)^2 in Kg^2*m^2*A^2/s^2, then k in J/K/mol, mass in Kg
    fac= hbar_mole_sq / ( 24d0 * 8.31446d0 * temp * mass_reduced / 1000d0 )


    Select Case( lj_bkghm_index(atom_id1,atom_id2) )
    Case(1)
       ! ********************************************bkghm*********************************
       ! first , short range, exponential part
       A=atype_lj_parameter(atom_id1,atom_id2,1)
       B=atype_lj_parameter(atom_id1,atom_id2,2)
       FH_short =  fac * A*(  B**2 * exp(-B*r_mag) + 2d0 / r_mag * -B * exp(-B*r_mag) )

       r_ij(2) = r_mag**2
       r_ij(6) = r_ij(2)**3
       r_ij(7) = r_ij(6) * r_mag
       r_ij(8) = r_ij(7) * r_mag
       ! determine whether we are using C8,C10 terms
       Select Case(C8_10_dispersion_terms)
       Case("no")
          ! just C6 terms, no damping
          ! d2Uij/dr2 contribution
          FH_long =  fac * ( 42d0 * atype_lj_parameter(atom_id1,atom_id2,3)/r_ij(8) )
          ! (2/r) * dU/dr contribution
          FH_long =  FH_long + fac * ( - 12d0 * atype_lj_parameter(atom_id1,atom_id2,3)/r_ij(8) )
       Case("yes")
          r_ij(9) = r_ij(8) * r_mag
          r_ij(10) = r_ij(9) * r_mag
          r_ij(11) = r_ij(10) * r_mag
          r_ij(12) = r_ij(11) * r_mag
          r_ij(13) = r_ij(12) * r_mag
          r_ij(14) = r_ij(13) * r_mag

          ! see if we're doing C12 dispersion
          Select Case(C12_dispersion)
          Case("yes")
             n_terms=4
          Case("no")
             n_terms=3
          End Select

          FH_long=0d0

          do i=1,n_terms
             C_index = i+2
             R_index = 6 + 2*(i-1)
             R_index1 = R_index + 1
             R_index2 = R_index1 + 1
             F=C6_C10_damp(atom_id1,atom_id2,r_mag,R_index)
             dF=dC6_C10_damp_scalar1(atom_id1,atom_id2,r_mag,R_index)
             d2F=dC6_C10_damp_scalar2(atom_id1,atom_id2,r_mag,R_index)
             C= atype_lj_parameter(atom_id1,atom_id2,C_index)

             ! d2Uij/dr2 contribution
             FH_long = FH_long + d2F*C/r_ij(R_index) -2d0*dF*dble(R_index)*C/r_ij(R_index1) + F*dble(R_index)*dble(R_index1)*C/r_ij(R_index2)
             ! (2/r) * dU/dr contribution
             FH_long = FH_long + 2d0 * ( dF * C / r_ij(R_index1) - F * dble(R_index) * C / r_ij(R_index2) )

          enddo

          FH_long= fac * FH_long

       end select

       E_local = FH_short - FH_long

    Case(2)
       ! ****************************************************  lj **********************************
       ! at this stage, all lj parameters should be expressed as C12 and C6, even though they were read in as epsilon and sigma
       r2inv = r_mag**(-2)
       r6inv = r2inv**(3)
       r8inv = r6inv * r2inv
       r14inv = r6inv * r8inv

       C12 = atype_lj_parameter(atom_id1,atom_id2,1)
       C6 = atype_lj_parameter(atom_id1,atom_id2,2)

       ! d2Uij/dr2 contribution
       E_local=fac * ( 156d0 * C12 * r14inv - 42d0 * C6 * r8inv )
       ! (2/r) * dU/dr contribution
       E_local = E_local + fac * 2d0 * ( -12d0 * C12 * r14inv + 6d0 * C6 * r8inv )

    End Select

  end subroutine Feynmann_Hibbs_lj_correction




  !************************************************
  ! this subroutine gives forces due to the Feynmann_Hibbs effective potential
  !
  ! As of now, it is only coded in for bkghm potential (lj_bkghm=1 ), and the code is set to stop
  ! if Feynmann_Hibbs_forces='yes' and lj_bkghm=2
  !*************************************************
  subroutine Feynmann_Hibbs_lj_forces(f_ij, rvec , atom_id1, atom_id2,i_mole,j_mole,n_atom )
    real*8,dimension(3), intent(out) :: f_ij
    real*8,dimension(3),intent(in) :: rvec
    integer,intent(in) :: atom_id1, atom_id2,i_mole,j_mole
    integer,dimension(:),intent(in) :: n_atom

    real*8  :: A,B,C,F,dF,d2F,d3F,fac, mass_reduced,rmag, r2,r3,r_ij(15), FH_short(3),FH_long(3),r_unit(3)
    real*8,parameter :: hbar_mole_sq=0.403280544d0  ! this is (hbar*NA)^2 ( in units of Kg^2 m^2*A^2/s^2 )
    integer :: n_terms,i,C_index,R_index,R_index1,R_index2,R_index3

    call generate_reduced_mass(mass_reduced, atom_id1, atom_id2,i_mole,j_mole,n_atom)

    ! fac = hbar^2/(24*kT * mass_red ) ; since distance is in angstroms, this should be in A^2
    ! use (hbar*mol)^2 in Kg^2*m^2*A^2/s^2, then k in J/K/mol, mass in Kg
    fac= hbar_mole_sq / ( 24d0 * 8.31446d0 * temp * mass_reduced / 1000d0 )


    r2 = dot_product( rvec, rvec)
    rmag = sqrt(r2)
    r_unit = rvec / rmag

    Select Case( lj_bkghm_index(atom_id1,atom_id2) )
    Case(1)
       ! ********************************************bkghm*********************************
       ! first , short range, exponential part
       A=atype_lj_parameter(atom_id1,atom_id2,1)
       B=atype_lj_parameter(atom_id1,atom_id2,2)
       FH_short =  r_unit * fac * A*(  -B**3 * exp(-B*rmag) + 2d0 / r2 * B * exp(-B*rmag) + 2d0 / rmag * B**2 * exp(-B*rmag) )


       r_ij(2) = r2
       r_ij(6) = r_ij(2)**3
       r_ij(7) = r_ij(6) * rmag
       r_ij(8) = r_ij(7) * rmag
       r_ij(9) = r_ij(8) * rmag
       ! determine whether we are using C8,C10 terms
       Select Case(C8_10_dispersion_terms)
       Case("no")
          ! just C6 terms, no damping
          ! gradient (d2Uij/dr2 ) contribution
          FH_long =  -r_unit * fac * ( 336d0 * atype_lj_parameter(atom_id1,atom_id2,3)/r_ij(9) )
          ! gradient ( (2/r) * dU/dr ) contribution
          FH_long =  FH_long + r_unit * fac * (  96d0 * atype_lj_parameter(atom_id1,atom_id2,3)/r_ij(9) )
       Case("yes")
          r_ij(10) = r_ij(9) * rmag
          r_ij(11) = r_ij(10) * rmag
          r_ij(12) = r_ij(11) * rmag
          r_ij(13) = r_ij(12) * rmag
          r_ij(14) = r_ij(13) * rmag
          r_ij(15) = r_ij(14) * rmag

          ! see if we're doing C12 dispersion
          Select Case(C12_dispersion)
          Case("yes")
             n_terms=4
          Case("no")
             n_terms=3
          End Select

          FH_long=0d0

          do i=1,n_terms
             C_index = i+2
             R_index = 6 + 2*(i-1)
             R_index1 = R_index + 1
             R_index2 = R_index1 + 1
             R_index3 = R_index2 + 1

             F=C6_C10_damp(atom_id1,atom_id2,rmag,R_index)
             dF=dC6_C10_damp_scalar1(atom_id1,atom_id2,rmag,R_index)
             d2F=dC6_C10_damp_scalar2(atom_id1,atom_id2,rmag,R_index)
             d3F=dC6_C10_damp_scalar3(atom_id1,atom_id2,rmag,R_index)
             C= atype_lj_parameter(atom_id1,atom_id2,C_index)

             ! gradient ( d2Uij/dr2 contribution )
             FH_long = FH_long + d3F*C/r_ij(R_index) - d2F*3d0*dble(R_index)*C/r_ij(R_index1) + 3d0*dF*dble(R_index)*dble(R_index1)*C/r_ij(R_index2) - F*dble(R_index)*dble(R_index1)*dble(R_index2)*C/r_ij(R_index3)
             ! gradient ( (2/r) * dU/dr contribution )
             FH_long = FH_long + 2d0*d2F*C/r_ij(R_index1)  - 2d0*(2d0*dble(R_index)+1d0)*dF*C/r_ij(R_index2) + F*2d0*dble(R_index)*dble(R_index2)*C/r_ij(R_index3) 

          enddo

          FH_long= r_unit * fac * FH_long

       end select

       f_ij = -FH_short + FH_long


    case default
       stop " lj_bkghm setting not recognized in Feynmann_Hibbs_lj_forces "
    End Select



  end subroutine Feynmann_Hibbs_lj_forces



  !**************************************
  !  This subroutine calculates the non-Coulombic, non-Lennard Jones,
  !  intermolecular repulsion between hydronium and water molecules
  !  see JPC B, 2008, 112, 467-482 (and errata, JPC B, 2008, 112, 7146)
  !**************************************
  subroutine ms_evb_intermolecular_repulsion( force, E_ms_evb_repulsion, n_mole , n_atom , xyz, box )
    real*8,dimension(:,:,:),intent(out) :: force
    real*8, intent(out)  :: E_ms_evb_repulsion
    integer, intent(in)  :: n_mole
    integer, dimension(:),intent(in) :: n_atom
    real*8, dimension(:,:,:), intent(in) :: xyz
    real*8, dimension(:,:), intent(in)  :: box

    integer :: i_mole, j_mole, i_mole_hydronium , j_mole_water, i_atom, j_atom, index, i_type, i_atom_oxygen
    real*8, dimension(3)  :: shift, rij, rij_O, q, fij, forceH2O_O
    real*8, dimension(4,3) :: forceH3O
    real*8  :: r_ij, sum, q2, r_OO, fac_OO, fac_OH, exp_q, switch_OO, dswitch_OO,switch_HO,dswitch_HO
    
        
!!$    write(*,*) "ms evb repulsion"

    E_ms_evb_repulsion = 0d0

    ! loop over hydronium molecules
    do i_mole=1, n_hydronium_molecules
       i_mole_hydronium = hydronium_molecule_index( i_mole )
       ! loop over water molecules
       do j_mole=1, n_mole
          ! these arrays are pairwise forces between two molecules
          forceH3O=0d0
          forceH2O_O=0d0
          ! see if this is a water molecule
          index = index_first_atom_type( j_mole , n_atom, water_oxygen_label )
          if ( index > 0 ) then
             j_mole_water = j_mole
             j_atom = index


             ! get hydronium oxygen
             i_atom_oxygen = index_first_atom_type( i_mole_hydronium, n_atom, hydronium_oxygen_label )

             shift = pbc_shift( xyz(i_mole_hydronium,i_atom_oxygen,:), xyz(j_mole_water,j_atom,:), box, xyz_to_box_transform)
             rij_O = -pbc_dr( xyz(i_mole_hydronium,i_atom_oxygen,:), xyz(j_mole_water,j_atom,:), shift )
             r_OO = dsqrt( dot_product( rij_O , rij_O ) )

             ! switching function
             call ms_evb_repulsive_switch( switch_OO, dswitch_OO, r_OO, rshiftOO_evb , rcutoffOO_evb )

             fac_OO = Bu_evb * exp(-bl_evb*( r_OO - dlOO0_evb ))

!!$             write(*,*) i_mole_hydronium, j_mole

             sum=0d0
             ! in this loop, calculate the hydronium hydrogen- water oxygen repulsion, as well as
             ! the sum over q coordinates, defined by (rO + rO)/2 - rH, for the oxygen-oxygen repulsion
             do i_atom=1, n_atom(i_mole_hydronium)
                if ( i_atom /= i_atom_oxygen ) then

                   ! hydronium proton, the negative sign creates displacement from water to hydronium atom
                   rij = -pbc_dr( xyz(i_mole_hydronium,i_atom,:), xyz(j_mole_water,j_atom,:), shift ) 
                   r_ij = dsqrt( dot_product( rij, rij ) )

                   fac_OH = Cu_evb * exp( -cl_evb * ( r_ij - dlOH0_evb ) )

                    ! switching function
                   call ms_evb_repulsive_switch( switch_HO, dswitch_HO, r_ij, rshiftOH_evb , rcutoffOH_evb )
                  
                   fij = rij / r_ij * fac_OH * ( switch_HO * cl_evb  - dswitch_HO )
                   E_ms_evb_repulsion = E_ms_evb_repulsion + switch_HO * fac_OH

!!$                   write(*,*) i_atom , r_ij, E_ms_evb_repulsion, switch_HO

                   forceH3O(i_atom,:) = forceH3O(i_atom,:) + fij(:)
                   forceH2O_O(:) = forceH2O_O(:) -fij(:)

                   ! THIS IS CORRECT equation. See errata, JPC B, 2008, 112, 7146.  Original paper
                   ! had R_HO distance instead of q for the oxygen-oxygen repulsion
                   q = ( 2d0 * xyz(j_mole_water,j_atom,:) + rij_O ) / 2d0 - ( xyz(j_mole_water,j_atom,:) + rij )
                   q2 = dot_product(q,q)
                   exp_q = exp(-blp_evb * q2)
                   sum = sum + exp_q

!!$                   write(*,*) "q2", q2
!!$                   write(*,*) q
                   

                ! note extra negative sign for dq/drH
                forceH3O(i_atom,:) = forceH3O(i_atom,:) + switch_OO * fac_OO * exp_q * -blp_evb * 2d0 * q     
                forceH3O(i_atom_oxygen,:) = forceH3O(i_atom_oxygen,:) + switch_OO * fac_OO * exp_q * blp_evb * q
                forceH2O_O(:) = forceH2O_O(:) + switch_OO * fac_OO * exp_q * blp_evb * q


                end if
             end do

             ! now oxygen-oxygen interaction, which depends on hydronium hydrogen positions    
             E_ms_evb_repulsion = E_ms_evb_repulsion + switch_OO * fac_OO * sum

!!$             write(*,*) "O-O", fac_OO, sum , switch_OO

             ! now derivatives of switch_OO * fac_OO
             fij = rij_O / r_OO * fac_OO * sum * ( switch_OO * bl_evb  - dswitch_OO )
             forceH3O(i_atom_oxygen,:) = forceH3O(i_atom_oxygen,:) + fij
             forceH2O_O(:) = forceH2O_O(:) - fij

             ! now add forces to global arrays
             do i_atom = 1, n_atom(i_mole_hydronium)
                force(i_mole_hydronium, i_atom,:) = force(i_mole_hydronium,i_atom,:) + forceH3O(i_atom,:)
             enddo
             force(j_mole_water,j_atom,:) = force(j_mole_water,j_atom,:) + forceH2O_O(:)
 

!!$             write(*,*) i_mole_hydronium, j_mole, E_ms_evb_repulsion

          end if
       end do
    end do

  end subroutine ms_evb_intermolecular_repulsion




  !************************************
  ! this is switching function (and derivative) for repulsive hydronium-water
  ! interactions
  !************************************
  subroutine ms_evb_repulsive_switch( switch, dswitch, r, rs , rc )
    real*8, intent(in) :: r, rs, rc
    real*8, intent(out) :: switch, dswitch
    real*8 :: term1, term2

    ! switch = 0 if r > rcutoff
    switch = 0d0
    dswitch = 0d0

    if ( r < rc ) then
       if (  r < rs ) then
          switch = 1d0
       else
          term1 = (r - rs)**2 / ( rc - rs )**3
          term2 = 3d0*rc - rs - 2d0*r
          switch = 1d0 - term1 * term2
          dswitch = -2d0 * (r - rs) * term2 / ( rc - rs )**3 + 2d0 * term1
       end if
    end if

  end subroutine ms_evb_repulsive_switch



  !**************************************
  !  This subroutine calculates the non-Coulombic, non-Lennard Jones,
  !  intermolecular repulsion forces between hydronium and water molecules
  !  see JPC B, 2008, 112, 467-482 (and errata, JPC B, 2008, 112, 7146)
  !**************************************
  subroutine ms_evb_intermolecular_repulsion_forces( force_i_mole, force_j_mole, i_mole , j_mole , n_atom , xyz , shift )
    real*8, dimension(:,:), intent(inout)  :: force_i_mole, force_j_mole
    integer, intent(in)  :: i_mole, j_mole
    integer, dimension(:),intent(in) :: n_atom
    real*8, dimension(:,:,:), intent(in) :: xyz
    real*8, dimension(:),intent(in) :: shift

    integer :: i_mole_hydronium , j_mole_water, i_atom, j_atom, index, i_type, i_atom_oxygen, flag
    real*8, dimension(3)  :: rij, rij_O, q, fij, forceH2O_O
    real*8, dimension(4,3) :: forceH3O
    real*8  :: r_ij, sum, q2, r_OO, fac_OO

    i_mole_hydronium=0
    ! see if either molecule is hydronium 
    index = index_first_atom_type( i_mole , n_atom, hydronium_oxygen_label )
    if ( index > 0 ) then
       i_mole_hydronium = i_mole
       j_mole_water = j_mole
       flag=1
    else
       index = index_first_atom_type( j_mole , n_atom, hydronium_oxygen_label )
       if ( index > 0 ) then
          i_mole_hydronium = j_mole
          j_mole_water = i_mole
          flag=-1
       end if
    end if


    ! if one of these molecules is hydronium
    if ( i_mole_hydronium > 0 ) then
       ! look for water oxygen atom
       j_atom = index_first_atom_type( j_mole_water , n_atom, water_oxygen_label )

       ! if this is water
       if ( j_atom > 0 ) then
          ! interaction of all hydronium atoms with water oxygen

          ! get hydronium oxygen
          i_atom_oxygen = index_first_atom_type( i_mole_hydronium, n_atom, hydronium_oxygen_label )
          rij_O = -pbc_dr( xyz(i_mole_hydronium,i_atom_oxygen,:), xyz(j_mole_water,j_atom,:), shift )
          r_OO = dsqrt( dot_product( rij_O , rij_O ) )
          fac_OO = Bu_evb * exp(-bl_evb*( r_OO - dlOO0_evb ))

          sum=0d0
          ! forces 
          do i_atom=1, n_atom(i_mole_hydronium)      
             if ( i_atom /= i_atom_oxygen ) then
                rij = -pbc_dr( xyz(i_mole_hydronium,i_atom,:), xyz(j_mole_water,j_atom,:), shift ) 
                r_ij = dsqrt( dot_product( rij, rij ) )
                fij = Cu_evb * cl_evb * exp( -cl_evb * ( r_ij - dlOH0_evb ) ) * rij / r_ij

                forceH3O(i_atom,:) = fij(:)
                forceH2O_O(:) = -fij(:)

                q = ( 2d0 * xyz(j_mole_water,j_atom,:) + rij_O ) / 2d0 - ( xyz(j_mole_water,j_atom,:) + rij )           
                q2 = dot_product(q,q)
                sum = sum + exp(-blp_evb * q2)
                ! note extra negative sign for dq/drH
                forceH3O(i_atom,:) = forceH3O(i_atom,:) + fac_OO * exp(-blp_evb * q2) * -blp_evb * 2d0 * q     
                forceH3O(i_atom_oxygen,:) = fac_OO * exp(-blp_evb * q2) * blp_evb * q
                forceH2O_O(:) = forceH2O_O(:) + fac_OO * exp(-blp_evb * q2) * blp_evb * q

             end if
          enddo

          fij = fac_OO * sum * -bl_evb * rij_O / r_OO
          forceH3O(i_atom_oxygen,:) = forceH3O(i_atom_oxygen,:) + fij
          forceH2O_O(:) = forceH2O_O(:) - fij


          ! now add the forces to the correct arrays
          if ( flag > 0 ) then
             do i_atom = 1, n_atom(i_mole)
                force_i_mole(i_atom,:) = force_i_mole(i_atom,:) + forceH3O(i_atom,:)
             enddo
             force_j_mole(j_atom,:) = force_j_mole(j_atom,:) + forceH2O_O(:)
          else
             ! the i,j's here are correct, even though they're mixed
             do i_atom = 1, n_atom(i_mole)
                force_j_mole(i_atom,:) = force_j_mole(i_atom,:) + forceH3O(i_atom,:)
             enddo
             force_i_mole(j_atom,:) = force_i_mole(j_atom,:) + forceH2O_O(:)
          end if

       end if

    end if



  end subroutine ms_evb_intermolecular_repulsion_forces


end module pairwise_interaction
