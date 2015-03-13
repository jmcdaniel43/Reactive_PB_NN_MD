module pme_routines
  use routines
  implicit none

  !*******************************************************************
  ! PARTICLE MESH EWALD SUBROUTINES
  !
  ! This module contains PME subroutines for energy and force
  ! pme routines use discrete fourier transforms in MKL library
  !
  !  reference for the pme algorithm is
  !  Essmann et al , J. Chem. Phys. 1995, 103, 8577-8593
  !*******************************************************************


contains





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This calculates electrostatic energy and forces
  ! using pme for reciprocal contribution
  !
  ! Global variables used but not changed
  ! integer, parameter :: spline_order,pme_grid
  ! real*8,dimension(pme_grid,pme_grid,pme_grid)::CB
  !
  ! Global variables changed
  ! real*8,dimension(pme_grid,pme_grid,pme_grid)::Q_grid,theta_conv_Q
  !
  !
  !  There is an input array target_atoms, which specifies the atoms that forces should be calculated for.
  !  This array has size n_mole, MAX_N_ATOM, and gives the atom labels of atoms for which forces are desired
  !  The first zero entry in the ith row (for the ith molecule) tells the code when to move to the next molecule
  !
  !  The output is an array pme_force(molecule,atom,direction), where ordering of atoms in each molecule corresponds
  !  to ordering in input array target_atoms, forces are derivatives w.r.t box coordinates, not scaled coordinates
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine  pme_energy_force(pme_energy, pme_force,target_atoms,tot_n_mole,n_mole, n_atom, xyz, chg, box, r_com, dfti_desc,dfti_desc_inv)
    use global_variables
    use MKL_DFTI
    use routines
    use omp_lib
    implicit none
    real*8, intent(out) :: pme_energy
    integer, intent(in) :: n_mole,tot_n_mole
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: chg, box, r_com
    integer,intent(in), dimension(:,:) :: target_atoms
    real*8,intent(out),dimension(:,:,:):: pme_force
    real*8, intent(in), dimension(:,:,:) :: xyz
    TYPE(DFTI_DESCRIPTOR),pointer,intent(in):: dfti_desc,dfti_desc_inv
    integer :: i, i_mole, i_atom, j_atom, a_index
    real*8 ::   Er, Ek, x, E_intra
    real*8, dimension(3) :: f_ij
    real*8,dimension(:,:), allocatable :: atom_list_xyz, atom_list_force
    real*8,dimension(:), allocatable  :: atom_list_chg


    pme_energy = 0d0
    pme_force = 0.D0

    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "real space pme calculation started at", time
    endif
    !***********************************************************************!


    !******************************** inter-molecular real-space interactions *******************************

    ! this works with or without framework atoms.  Without framework atoms, tot_n_mole = n_mole. With framework atoms, forces on solute from framework are included since tot_n_mole includes framework molecules (atoms)

    Select Case(use_verlet_list)
    Case("no")    
       call pme_real_space_no_verlet(Er, pme_force,tot_n_mole,n_mole, n_atom, xyz, chg, box, r_com)
    Case("yes")
       ! here we need to map molecular data structures to atomic data structures that are
       ! consistent with the verlet list
       allocate( atom_list_xyz(verlet_elec_atoms,3), atom_list_force(verlet_elec_atoms,3) , atom_list_chg(verlet_elec_atoms) )
       do i_mole = 1, tot_n_mole
          do i_atom = 1, n_atom(i_mole)
             ! get verlet atom index for this atom
             a_index = verlet_molecule_map_elec(i_mole,i_atom)
             if (  a_index > 0  ) then
                ! store atom information from molecular data arrays
                atom_list_xyz(a_index,:) = xyz(i_mole,i_atom,:)
                atom_list_chg(a_index) = chg(i_mole,i_atom)
             end if
          enddo
       enddo
       call pme_real_space_use_verlet(Er, atom_list_force,atom_list_xyz, atom_list_chg, box)
       ! now add forces to molecular data structures
       do i_mole = 1, tot_n_mole
          do i_atom = 1, n_atom(i_mole)
             ! get verlet atom index for this atom
             a_index = verlet_molecule_map_elec(i_mole,i_atom)
             if (  a_index > 0  ) then
                pme_force(i_mole,i_atom,:) = atom_list_force(a_index,:)
             end if
          enddo
       enddo

       deallocate( atom_list_xyz , atom_list_force, atom_list_chg )
    End Select



    !******************************** intra-molecular real-space interactions *******************************
    do i_mole=1,n_mole
       if(n_atom(i_mole) .gt. 1) then
          do i_atom=1,n_atom(i_mole)-1
             do j_atom=i_atom+1,n_atom(i_mole)
                call intra_pme_energy(E_intra,xyz,chg,i_mole,i_atom,j_atom,n_atom)
                Er = Er + E_intra
                call intra_pme_force(f_ij,xyz,chg,i_mole,n_atom,i_atom,j_atom)
                pme_force(i_mole,i_atom,:)=pme_force(i_mole,i_atom,:) + f_ij(:)
                pme_force(i_mole,j_atom,:)=pme_force(i_mole,j_atom,:) - f_ij(:)                         
             enddo
          enddo
       endif
    enddo



!!$    write(*,*) "pme force"
!!$    write(*,*) pme_force(1,3,:)
!!$    write(*,*) pme_force(12,1,:)
!!$    write(*,*) pme_force(61,2,:)
!!$    write(*,*) pme_force(121,3,:)
!!$    write(*,*) pme_force(256,1,:)
!!$    write(*,*) pme_force(319,2,:)
!!$    stop

    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "real space part of pme calculation finished at", time
    endif
    !***********************************************************************!

    !******************************** PME reciprocal space ***************************************!

!!!!!!!!!!!!!!!!!!!!!!!remember pme_force array will be reordered after call to pme_force_recip based on target_atoms
    call pme_force_recip(pme_force,target_atoms,tot_n_mole,n_mole, n_atom, xyz, chg, box, dfti_desc,dfti_desc_inv)
    pme_force=pme_force * 0.52914D0 * 627.51D0 * 4.184D0 ! convert from e^2/A^2 to kJ/mol/A

    Ek= pme_recip( tot_n_mole, n_atom, xyz, chg, box,dfti_desc,dfti_desc_inv, 0)
    pme_energy = Ek + Er - Ewald_framework + Ewald_self
    pme_energy = pme_energy * 0.52914D0 * 627.51D0 * 4.184D0 ! convert from e^2/A to kJ/mol

    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "pme calculation finished at", time
    endif
    !***********************************************************************!

  end subroutine  pme_energy_force




  !***********************************************
  ! this subroutine calculates the real space part of 
  ! the pme energy using an atom-atom based verlet list
  !
  ! In this subroutine, the data structures have been changed from molecule based
  ! data structures, to atom_list based data structures.  This is because we only want
  ! to consider inter-molecular interactions between atoms that have non-zero
  ! charges.  The indexing of atoms in the atom_list data structures
  ! should be consistent with the indexing of atoms in the verlet list
  !***********************************************************************
  subroutine pme_real_space_use_verlet(pme_real_space_energy, atom_list_force,atom_list_xyz, atom_list_chg, box)
    use global_variables
    use routines
    use omp_lib
    implicit none
    real*8, intent(out) :: pme_real_space_energy
    real*8, intent(in), dimension(:,:) ::  box
    real*8,intent(out),dimension(:,:):: atom_list_force
    real*8, intent(in), dimension(:,:) :: atom_list_xyz
    real*8, intent(in), dimension(:)   :: atom_list_chg

    integer :: i, i_atom, j_atom, thread_id, i_thread, verlet_start, verlet_finish, i_index
    real*8, dimension(3) :: rij, shift, shift_direct, dr_com, dr_direct, f_ij
    real*8, dimension(size(atom_list_force(:,1)),3) :: local_force
    real*8, dimension(n_threads,size(atom_list_force(:,1)),3) :: temp_force
    real*8 ::  norm_dr, norm_dr2, ewald_cutoff2, erfc_value, factor, coulomb, x, Er
    integer :: split_do

    Er=0d0
    atom_list_force = 0.D0
    local_force = 0.D0
    temp_force= 0.D0
    ewald_cutoff2 = ewald_cutoff ** 2

    ! decide how to split the parallel section
    if (n_threads .eq. 1 ) then
       split_do = 1
    else
       split_do = verlet_elec_atoms/n_threads+1
    endif

    !**************************************** use Verlet list ****************************************************************
    call OMP_SET_NUM_THREADS(n_threads)
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(atom_list_xyz, atom_list_chg, verlet_elec_atoms, box, temp_force,ewald_cutoff2, split_do , xyz_to_box_transform, alpha_sqrt,verlet_point_elec,verlet_neighbor_list_elec,erfc_table) REDUCTION(+:Er)
    !$OMP CRITICAL
    local_force = 0.D0
    !$OMP END CRITICAL
    !$OMP DO SCHEDULE(DYNAMIC,split_do)
    ! note that in Verlet list, neighbors (for molecule j_mole) are only stored if j_mole > i_mole , to avoid double counting
    do i_atom=1 , verlet_elec_atoms
       factor = 2.D0*alpha_sqrt/pi_sqrt

       ! get j_mole from Verlet list
       verlet_start = verlet_point_elec(i_atom)
       i_index = i_atom + 1
       verlet_finish = verlet_point_elec(i_index) - 1

       ! make sure there's at least one neighbor
       if ( verlet_finish .ge. verlet_start ) then
          do i_index= verlet_start, verlet_finish 
             j_atom=verlet_neighbor_list_elec(i_index)

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

             if ( norm_dr2 < ewald_cutoff2 ) then
                norm_dr = dsqrt( norm_dr2 )
                x = norm_dr*alpha_sqrt
                if ( x < erfc_max ) then
                   erfc_value = erfc_table(ceiling(x/erfc_max*dble(erfc_grid)))
                   coulomb = atom_list_chg(i_atom) * atom_list_chg(j_atom) / norm_dr
                   Er = Er + erfc_value * coulomb
                   f_ij = coulomb * rij *(  erfc_value / norm_dr2 + factor * exp(-x**2) / norm_dr)
                   local_force(i_atom,:) = local_force(i_atom,:) + f_ij(:)
                   local_force(j_atom,:) = local_force(j_atom,:) - f_ij(:)
                end if
             end if
          end do
       end if
    end do

    !$OMP END DO NOWAIT
    thread_id = OMP_GET_THREAD_NUM()
    temp_force(thread_id+1,:,:)= local_force(:,:)
    !$OMP END PARALLEL


    do i_thread=1,n_threads
       do i_atom=1, verlet_elec_atoms
          atom_list_force(i_atom,:) = atom_list_force(i_atom,:) + temp_force(i_thread,i_atom,:)
       enddo
    enddo

    pme_real_space_energy = Er

  end subroutine pme_real_space_use_verlet




  !***********************************************
  ! this subroutine is self-explanatory
  !***********************************************
  subroutine pme_real_space_no_verlet(pme_real_space_energy, pme_force,tot_n_mole,n_mole, n_atom, xyz, chg, box, r_com)
    use global_variables
    use routines
    use omp_lib
    implicit none
    real*8, intent(out) :: pme_real_space_energy
    integer, intent(in) :: n_mole,tot_n_mole
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: chg, box, r_com
    real*8,intent(out),dimension(:,:,:):: pme_force
    real*8, intent(in), dimension(:,:,:) :: xyz
    integer :: i, i_mole, i_atom, j_mole, j_atom, thread_id, i_thread
    real*8, dimension(3) :: rij, shift, shift_direct, dr_com, dr_direct, f_ij
    real*8, dimension(n_mole,maxval(n_atom),3) :: local_force
    real*8, dimension(n_threads,n_mole,maxval(n_atom),3) :: temp_force
    real*8 ::  norm_dr, norm_dr2, ewald_cutoff2, erfc_value, factor, coulomb, x, Er
    integer :: split_do

    Er = 0d0
    pme_force = 0.D0
    local_force = 0.D0
    temp_force= 0.D0
    ewald_cutoff2 = ewald_cutoff ** 2

    ! decide how to split the parallel section
    if (n_threads .eq. 1 ) then
       split_do = 1
    else
       split_do = n_mole/n_threads+1
    endif

    !**************************************** no Verlet list, naive looping ***********************************************
    call OMP_SET_NUM_THREADS(n_threads)
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(n_mole,tot_n_mole,r_com,box,n_atom,xyz,chg,temp_force,ewald_cutoff2,screen_type, thole, springcon, split_do , xyz_to_box_transform, alpha_sqrt,erfc_table) REDUCTION(+:Er)
    !$OMP CRITICAL
    local_force = 0.D0
    !$OMP END CRITICAL
    !$OMP DO SCHEDULE(DYNAMIC,split_do)
    do i_mole=1,n_mole
       factor = 2.D0*alpha_sqrt/pi_sqrt
       do j_mole = i_mole+1, tot_n_mole
          do i_atom = 1, n_atom( i_mole )
             do j_atom = 1, n_atom( j_mole )

                rij = xyz(i_mole,i_atom,:) - xyz(j_mole,j_atom,:)
                dr_direct(:) = matmul( xyz_to_box_transform, rij )
                do i=1,3
                   shift_direct(i) = dble(floor( dr_direct(i) + 0.5d0 ))
                enddo
                shift = matmul( shift_direct , box )
                rij = rij - shift
                norm_dr2 = dot_product( rij, rij )
                if ( norm_dr2 < ewald_cutoff2 ) then
                   norm_dr = dsqrt( norm_dr2 )
                   x = norm_dr*alpha_sqrt
                   if ( x < erfc_max ) then
                      erfc_value = erfc_table(ceiling(x/erfc_max*dble(erfc_grid)))
                      coulomb = chg(i_mole,i_atom) * chg(j_mole,j_atom)  / norm_dr
                      Er = Er + erfc_value * coulomb
                      f_ij = coulomb * rij *(  erfc_value / norm_dr2 + factor * exp(-x**2) / norm_dr)

                      ! these next lines are code which allow for real space screening function
                      ! we have commented them and added the above code instead to gain speed
                      ! f_ij = pairwise_ewald_real_space_force( norm_dr2, rij, i_mole, j_mole, i_atom, j_atom, chg )
                      ! remove reciprocal space term
                      !Er = Er + chg(i_mole,i_atom) * chg(j_mole,j_atom) * (erfc_g(norm_dr*alpha_sqrt)-1.D0) / norm_dr
                      ! add "screened" real space interaction
                      !Er = Er + chg(i_mole,i_atom) * chg(j_mole,j_atom) * screen(i_mole,j_mole,i_atom,j_atom,norm_dr) / norm_dr

                      local_force(i_mole,i_atom,:) = local_force(i_mole,i_atom,:) + f_ij(:)
                      local_force(j_mole,j_atom,:) = local_force(j_mole,j_atom,:) - f_ij(:)
                   end if
                end if
             enddo
          enddo
       end do
    enddo

    !$OMP END DO NOWAIT
    thread_id = OMP_GET_THREAD_NUM()
    temp_force(thread_id+1,:,:,:)= local_force(:,:,:)
    !$OMP END PARALLEL

    !************************************************ distribute real space force ******************************
    do i_thread=1,n_threads
       do i_mole=1, n_mole
          do i_atom=1,n_atom(i_mole)
             pme_force(i_mole,i_atom,:)= pme_force(i_mole,i_atom,:) + temp_force(i_thread,i_mole,i_atom,:)
          enddo
       enddo
    enddo

    pme_real_space_energy = Er

  end subroutine pme_real_space_no_verlet




  !**********************************************************
  ! This function calculates the pairwise real space electrostatic force in an Ewald sum
  ! First, the reciprocal space contribution (charge - Gaussian interaction) is subtracted off,
  ! and then the true real space interaction (screened) point charge Coulomb is added on
  !
  ! this function can either explicitly evaluate the pairwise interaction, or use a 
  ! look up table for each atomtype pair
  !
  ! to use the lookup table, set grid_ewald_realspace_interaction="yes"
  !**********************************************************
  function pairwise_ewald_real_space_force( norm_dr2, rij, i_mole, j_mole, i_atom, j_atom, chg )
    use global_variables
    use routines
    real*8,dimension(3) :: pairwise_ewald_real_space_force
    real*8, intent(in), dimension(:,:) :: chg
    integer, intent(in) :: i_mole, j_mole, i_atom, j_atom
    real*8,dimension(3),intent(in) :: rij
    real*8,intent(in) :: norm_dr2

    real*8 :: norm_dr, norm_dr3
    real*8 , dimension(3) :: f_ij
    integer :: p_atom1, p_atom2, atom_id1, atom_id2,index


    Select Case(grid_ewald_realspace_interaction)
    Case("yes")

       !********* find the mapping to the table array ************
       ! first, find parent atoms for drude oscillators ( or just returns atom if not a shell )
       p_atom1 = drude_atom_map(i_mole,i_atom)
       p_atom2 = drude_atom_map(j_mole,j_atom)
       ! now find atomtypes for these atoms
       atom_id1 = atom_index(i_mole,p_atom1)
       atom_id2 = atom_index(j_mole,p_atom2)
       ! now  map atomtype pair to the table array
       index = ewald_realspace_interaction_grid_map(atom_id1, atom_id2)

       ! remember, lookup tables are tabulated in terms of r^2
       ! force is equal to chg*chg*rvec*tablevalue
       pairwise_ewald_real_space_force = chg(i_mole,i_atom) * chg(j_mole,j_atom) * rij * ewald_realspace_interaction_grid(index,ceiling(norm_dr2/max_ewald_table *dble(ewald_realspace_interaction_grid_size)))

    Case default
       ! explicitly evaluate the interactions
       norm_dr = dsqrt( norm_dr2 )
       norm_dr3 = norm_dr2 * norm_dr

       ! screening functions are expensive, only use them up to a certain distance
       ! for an exponent=3.0 Ang^-1 ( which is very conservative ), the Tang-Toennies screening function at
       ! 7 Angstroms is 0.9999999833, so setting it equal to 1 at this distance, should introduce negligable error
       if ( norm_dr2 > screen_distance_max_sq ) then
          f_ij = chg(i_mole,i_atom) * chg(j_mole,j_atom) * rij *( erfc_g(norm_dr*alpha_sqrt) / norm_dr3 +2.D0*alpha_sqrt/pi_sqrt*exp(-(alpha_sqrt*norm_dr)**2)/norm_dr2)
       else
          !! subtract reciprocal space term
          f_ij = chg(i_mole,i_atom) * chg(j_mole,j_atom) * rij *(( erfc_g(norm_dr*alpha_sqrt)- 1D0) / norm_dr3 +2.D0*alpha_sqrt/pi_sqrt*exp(-(alpha_sqrt*norm_dr)**2)/norm_dr2)
          !! add "screened" real space interaction
          f_ij = f_ij + chg(i_mole,i_atom) * chg(j_mole,j_atom) * rij * screen(i_mole,j_mole,i_atom,j_atom,norm_dr) / norm_dr3
          f_ij = f_ij - chg(i_mole,i_atom) * chg(j_mole,j_atom) * d_screen(i_mole,j_mole,i_atom,j_atom,rij,norm_dr) / norm_dr
       endif

       pairwise_ewald_real_space_force = f_ij

    End Select



  end function pairwise_ewald_real_space_force





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! this is for the reciprocal contribution
  ! this function uses particle mesh ewald to compute electrostatic
  ! energy. It uses cardinal B splines, which can be differentiated for
  ! forces.  interfaced to MKL library for DFT
  ! based on Essmann paper, J. Chem. Phys. 103 (19) 1995
  ! definition of forward and backward dft in paper is reversed with definition
  ! in MKL library
  ! currently, this uses complex to complex FT, because real to complex
  ! wasn't implemented in 3D when this was written
  !
  !  Q_change is an input variable that determines whether the Q_array needs updated (1), or not (0)
  !  only input Q_change=0 if forces have already been calculated
  !  If Q_change is zero, use Q array and theta_conv_Q from glob_v
  !
  ! Global variables used but not changed
  ! integer, parameter :: spline_order,pme_grid
  ! real*8,dimension(pme_grid,pme_grid,pme_grid)::CB
  !
  ! Global variables changed
  ! real*8,dimension(pme_grid,pme_grid,pme_grid)::Q_grid,theta_conv_Q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function pme_recip( n_mole, n_atom, xyz, chg, box, dfti_desc,dfti_desc_inv,Q_change)
    use global_variables
    use MKL_DFTI
    implicit none
    integer, intent(in) :: n_mole,Q_change
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: chg, box
    real*8, intent(in), dimension(:,:,:) :: xyz
    TYPE(DFTI_DESCRIPTOR),pointer,intent(in):: dfti_desc,dfti_desc_inv
    real*8, dimension(:,:,:),allocatable::xyz_scale
    real*8 :: a(3), b(3), c(3), ka(3), kb(3), kc(3),kk(3,3)
    real*8 :: vol,scale
    integer::i,j,l,status,length(3),n,K  ! n is spline order, K is grid size
!!!!!!!!!!!!!!!!!!!!!!!even though arrays begin at index 1 for convenience, 1,1,1 corresponds to 0,0,0
    complex*16,dimension(pme_grid,pme_grid,pme_grid)::FQ
    real*8,dimension(pme_grid**3)::q_1r
    complex*16,dimension(pme_grid**3)::q_1d


    n=spline_order
    K=pme_grid

    length(:)=K
    a(:) = box(1,:)
    b(:) = box(2,:)
    c(:) = box(3,:)

    ! calculate the volume and the reciprocal vectors (notice no 2pi)
    vol = volume( a, b, c )
    call crossproduct( a, b, kc ); kc = kc /vol 
    call crossproduct( b, c, ka ); ka = ka /vol
    call crossproduct( c, a, kb ); kb = kb /vol
    kk(1,:)=ka(:);kk(2,:)=kb(:);kk(3,:)=kc(:)


    ! create scaled coordinates
    allocate(xyz_scale(n_mole,maxval(n_atom(:)),3))
    do i=1,n_mole
       do j=1,n_atom(i)
          do l=1,3
             xyz_scale(i,j,l)=dble(K)*dot_product(kk(l,:),xyz(i,j,:))
             ! if atoms are out of grid, shift them back in
             if (xyz_scale(i,j,l)<0.) then
                xyz_scale(i,j,l)=xyz_scale(i,j,l)+dble(K)
                ! else if(xyz_scale(i,j,l)>dble(K)) then ! changed by Jesse, 6/10/2014 see below line.  In rare cases, this could lead to out of bounds array calls when constructing Q
             else if(xyz_scale(i,j,l)>= dble(K)) then
                xyz_scale(i,j,l)=xyz_scale(i,j,l)-dble(K)
             endif
          enddo
       enddo
    enddo



!!!!!!!! if need to update Q_array
    if(Q_change.ne.0) then
       call grid_Q(Q_grid,chg,xyz_scale,n_atom,n_mole,K,n)


       q_1r=RESHAPE(Q_grid, (/K**3/) )
       q_1d=cmplx(q_1r,0.,16)

       status=DftiComputeForward(dfti_desc, q_1d)


!!!!! need Finv(B*C)convoluted w/ Q
!!!!! equals Finv(F(Finv(B*C) conv Q))
!!!!! equals K**3*Finv(B*C*F(Q))

       FQ=RESHAPE(q_1d, (/K,K,K/) )


!!!!!!!multiply B*C*F(Q)
       !FQ=FQ*cmplx(CB,0.,16)
       FQ=FQ*CB

!!!!!!! take Finv
       q_1d=RESHAPE(FQ,(/K**3/) )
!!$  scale = 1.0D0/dble(K)**3


       status = DftiComputeBackward(dfti_desc_inv, q_1d)
       FQ=RESHAPE(q_1d, (/K,K,K/) )
       theta_conv_Q=dble(FQ)

    endif

    pme_recip=.5D0*sum((Q_grid*theta_conv_Q))

    deallocate(xyz_scale)

  end function pme_recip




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This calculates electrostatic forces
  ! based on Essmann paper, J. Chem. Phys. 103 (19) 1995
  ! definition of forward and backward dft in paper is reversed with definition
  ! in MKL library
  !
  ! Global variables used but not changed
  ! integer, parameter :: spline_order,pme_grid
  ! real*8,dimension(pme_grid,pme_grid,pme_grid)::CB
  !
  ! Global variables changed
  ! real*8,dimension(pme_grid,pme_grid,pme_grid)::Q_grid,theta_conv_Q
  !
  !  There is an input array target_atoms, which specifies the atoms that forces should be calculated for.
  !  This array has size n_mole, MAX_N_ATOM, and gives the atom labels of atoms for which forces are desired
  !  The first zero entry in the ith row (for the ith molecule) tells the code when to move to the next molecule
  !
  !  The output is an array pme_force(molecule,atom,direction), where ordering of atoms in each molecule corresponds
  !  to ordering in input array target_atoms, forces are derivatives w.r.t box coordinates, not scaled coordinates
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine  pme_force_recip(pme_force,target_atoms,tot_n_mole,n_mole, n_atom, xyz, chg, box, dfti_desc,dfti_desc_inv)
    use global_variables
    use MKL_DFTI
    use omp_lib
    implicit none
    integer, intent(in) :: n_mole,tot_n_mole
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: chg, box
    integer,intent(in), dimension(:,:) :: target_atoms
    real*8,intent(inout),dimension(:,:,:):: pme_force
    real*8, intent(in), dimension(:,:,:) :: xyz
    TYPE(DFTI_DESCRIPTOR),pointer,intent(in):: dfti_desc,dfti_desc_inv
    real*8, dimension(:,:,:),allocatable::xyz_scale
    real*8 :: a(3), b(3), c(3), ka(3), kb(3), kc(3),kk(3,3)
    real*8 :: vol,scale,force(3)
    integer::i,j,l,i_atom,i_mole,status,length(3),n,K  ! n is spline order, K is grid size
!!!!!!!!!!!!!!!!!!!!!!!even though arrays begin at index 1 for convenience, 1,1,1 corresponds to 0,0,0
    complex*16,dimension(pme_grid,pme_grid,pme_grid)::FQ
    real*8,dimension(pme_grid**3)::q_1r
    complex*16,dimension(pme_grid**3)::q_1d
    integer:: split_do
    real*8 :: small=1D-6


    n=spline_order
    K=pme_grid

    length(:)=K
    a(:) = box(1,:)
    b(:) = box(2,:)
    c(:) = box(3,:)

    ! calculate the volume and the reciprocal vectors (notice no 2pi)
    vol = volume( a, b, c )
    call crossproduct( a, b, kc ); kc = kc /vol 
    call crossproduct( b, c, ka ); ka = ka /vol
    call crossproduct( c, a, kb ); kb = kb /vol
    kk(1,:)=ka(:);kk(2,:)=kb(:);kk(3,:)=kc(:)

    ! create scaled coordinates
    allocate(xyz_scale(tot_n_mole,maxval(n_atom(:)),3))

    do i=1,tot_n_mole
       do j=1,n_atom(i)
          do l=1,3
             xyz_scale(i,j,l)=dble(K)*dot_product(kk(l,:),xyz(i,j,:))
             ! if atoms are out of grid, shift them back in
             if (xyz_scale(i,j,l)<0.) then
                xyz_scale(i,j,l)=xyz_scale(i,j,l)+dble(K)
                !else if(xyz_scale(i,j,l)>dble(K)) then ! changed by Jesse, 6/10/2014 see below line.  In rare cases, this could lead to out of bounds array calls when constructing Q
             else if(xyz_scale(i,j,l) >= dble(K)) then
                xyz_scale(i,j,l)=xyz_scale(i,j,l)-dble(K)
             endif

             ! test
             if ( (xyz_scale(i,j,l) < 0d0 ) .or. (xyz_scale(i,j,l) > dble(K) ) ) then
                write(*,*) "scaled coordinates out of box"
                stop
             end if

          enddo
       enddo
    enddo

    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "Q_grid started at", time
    endif
    !************************************************************************!

    call grid_Q(Q_grid,chg,xyz_scale,n_atom,tot_n_mole,K,n)

    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "Q_grid finished at", time
    endif
    !************************************************************************!

    q_1r=RESHAPE(Q_grid, (/K**3/) )
    q_1d=cmplx(q_1r,0.,16)

    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "Fourier transforms started at", time
    endif
    !************************************************************************!

    status=DftiComputeForward(dfti_desc, q_1d)

!!!!! need Finv(B*C)convoluted w/ Q
!!!!! equals Finv(F(Finv(B*C) conv Q)
!!!!! equals K**3*Finv(B*C*F(Q))


    FQ=RESHAPE(q_1d, (/K,K,K/) )

    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "B*C*FQ started at", time
    endif
    !************************************************************************!

!!!!!!!multiply B*C*F(Q)
    ! FQ=FQ*cmplx(CB,0.,16)
    FQ=FQ*CB

    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "B*C*FQ finished at", time
    endif
    !************************************************************************!


!!!!!!! take Finv
    q_1d=RESHAPE(FQ,(/K**3/) )
!!$  scale = 1.0D0/dble(K)**3


    status = DftiComputeBackward(dfti_desc_inv, q_1d)

    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "Fourier transforms finished at", time
    endif
    !************************************************************************!

    FQ=RESHAPE(q_1d, (/K,K,K/) )
    theta_conv_Q=dble(FQ)


!!!!!!!!!! now compute forces for all desired atoms
!!!!!!!!!! remember, pme_force is input with forces on all atoms, so have to adjust these just to 
!!!!!!!!!! target atoms, so that you don't waste time taking derivatives

    ! forces are not calculated on framework, so loop over n_mole, not tot_n_mole

    ! decide how to split the parallel section
    if (n_threads .eq. 1 ) then
       split_do = 1
    else
       split_do = n_mole/n_threads+1
    endif

    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "Q grid derivatives started at", time
    endif
    !************************************************************************!


    call OMP_SET_NUM_THREADS(n_threads)
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(n_threads,n_mole,n_atom,target_atoms,theta_conv_Q,xyz_scale,chg,K,box,n,kk,pme_force,spline_order,split_do,small) 
    !$OMP DO SCHEDULE(dynamic, split_do)
    do i_mole=1,n_mole
       do i_atom=1,n_atom(i_mole)
          if(target_atoms(i_mole,i_atom).eq.0) then
             goto 100
!!$!!!!!!!!!!!!!move to next molecule
          endif

          if ( abs(chg(i_mole, target_atoms(i_mole,i_atom))) > small ) then
          call derivative_grid_Q(force,theta_conv_Q,chg,xyz_scale,target_atoms(i_mole,i_atom),i_mole,K,box,n,kk)
          pme_force(i_mole,i_atom,:)=pme_force(i_mole,target_atoms(i_mole,i_atom),:)+force(:)
          end if
       enddo
100    continue
    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    deallocate(xyz_scale)



  end subroutine pme_force_recip


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! this subroutine removes intra-molecular energy contributions
  ! and if drude oscillators are present, adds screened intra-
  ! molecular contribution
  !
  ! use explicit erfc since r can be very small
  !
  ! note if charge is zero on drude oscillator, pol =0, and thole_screen will blow up
  ! but should be fine because drude_p is 0 for 0 charge
  !
  ! Note that if we wish to subtract out the intra molecular interaction between two charges on top 
  ! of each other, (suppose drude oscillators are directly on top of their corresponding atoms)
  ! we need to subtract this contribution to the reciprocal space energy.  In the reciprocal
  ! space contribution, there is a factor of 1/2 for overcounting, so this contribution is
  !  1/2 * ( E12int + E21int ) = E12int , where E12int is charge 1 interacting with a gaussian
  ! of charge 2 right underneath it.  This energy is then
  ! 2 * q1 * q2 * ( alpha / pi ) ^ (1/2)  
  ! note that this is a factor of 2 greater than the Ewald self correction, because there are two
  ! interactions that we are accounting for, rather than one.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine intra_pme_energy(E_intra,xyz,chg,i_mole,i_atom,j_atom,n_atom)
    use global_variables
    implicit none
    real*8, intent(out):: E_intra
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: chg
    real*8, intent(in), dimension(:,:,:) :: xyz
    integer,intent(in)::i_mole,i_atom,j_atom

    real*8,dimension(3)::rij
    integer::i_drude,j_drude,n_pairs,sign_chg_i,sign_chg_j,flag_drudei,flag_drudej, flag_same_atom, i_mole_type
    real*8::norm_dr,pol1,pol2
    real*8,parameter::small=1D-8


!!!!!!!!!!!!!! subtract unscreened interaction with total charge
    rij(:) = xyz(i_mole,i_atom,:) - xyz(i_mole, j_atom,:)
    norm_dr = sqrt( dot_product( rij, rij ) )
    if(norm_dr < small) then
       E_intra= - 2d0 * chg(i_mole,i_atom) * chg(i_mole,j_atom) * alpha_sqrt/sqrt(pi)
    else
       E_intra = chg(i_mole,i_atom) * chg(i_mole,j_atom) * (erfc(norm_dr*alpha_sqrt)-1.D0) / norm_dr
    endif

    ! add real space interaction if no exclusion between these atoms
    i_mole_type = molecule_index(i_mole)
    ! check for exclusions 
    if ( molecule_exclusions( i_mole_type, i_atom, j_atom ) /= 1 ) then 
       E_intra = E_intra + chg(i_mole,i_atom) * chg(i_mole,j_atom) / norm_dr
    end if




!!!!!!!!!if oscillators present, get rid of intra-molecular fixed chg-chg interactions, but allow screened induced-dipole interactions
    call drude_atom_indices(i_mole,n_atom,i_atom,j_atom,i_drude,j_drude,sign_chg_i,sign_chg_j,flag_drudei,flag_drudej, flag_same_atom)

    ! we only consider intra molecular interactions between drude oscillators, or the corresponding drude oscillator charge on an atom
    ! also, we do not consider intra molecular interactions between a drude oscillator and it's corresponding atom


    ! flag_drudei tells us whether i_atom is either a drude oscillator, or has a drude_oscillator attached to it
    ! flag_drudej tells us whether j_atom is either a drude oscillator, or has a drude_oscillator attached to it
    ! flag_same_atom tells us whether j_atom is the drude oscillator on i_atom, or vice-versa

    if ( (flag_drudei == 1 ) .and. (flag_drudej == 1 ) .and. (flag_same_atom == 0 ) ) then
!!!!!!!!!!!!!!! add screened atom-atom, get atom chg from corresponding oscillator
       pol1=chg(i_mole,i_drude)**2/springcon; pol2=chg(i_mole,j_drude)**2/springcon
       E_intra = E_intra + dble(sign_chg_i*sign_chg_j)*chg(i_mole,i_drude) * chg(i_mole,j_drude) * thole_screen(pol1,pol2,norm_dr,thole)/ norm_dr
    endif


  end subroutine intra_pme_energy




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! this subroutine removes intra-molecular interactions for normal molecules
  ! and adds intra-molecular interactions for drude oscillators.  
  !
  ! use explicit erfc since r can be very small
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine intra_pme_force(f_ij,xyz,chg,i_mole,n_atom,i_atom,j_atom)
    use global_variables
    implicit none
    real*8, dimension(:), intent(out):: f_ij
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: chg
    real*8, intent(in), dimension(:,:,:) :: xyz
    integer,intent(in)::i_mole,i_atom,j_atom

    real*8,dimension(3)::rij
    integer::i_drude,j_drude,sign_chg_i,sign_chg_j,flag_drudei,flag_drudej, flag_same_atom, i_mole_type
    real*8::norm_dr,pol1,pol2
    real*8,parameter::small=1D-8


!!!!!!!!!!!!remove unscreened interaction
    rij(:) = xyz(i_mole,i_atom,:) - xyz(i_mole, j_atom,:)
    norm_dr = sqrt( dot_product( rij, rij ) )
    if(norm_dr > small) then
       f_ij= chg(i_mole,i_atom) * chg(i_mole,j_atom) * rij * ( (erfc(norm_dr*alpha_sqrt)-1.D0) / norm_dr**3 + (2.D0*alpha_sqrt/sqrt(pi))*exp(-(alpha_sqrt*norm_dr)**2)/norm_dr**2)
    else
       f_ij=0.D0
    endif


    ! add real space interaction if no exclusion between these atoms
    i_mole_type = molecule_index(i_mole)
    ! check for exclusions 
    if ( molecule_exclusions( i_mole_type, i_atom, j_atom ) /= 1 ) then 
       f_ij = f_ij + chg(i_mole,i_atom) * chg(i_mole,j_atom) * rij / norm_dr**3
    end if


!!!!!!!!!if oscillators present, get rid of intra-molecular fixed chg-chg interactions, but allow screened induced-dipole interactions
    call drude_atom_indices(i_mole,n_atom,i_atom,j_atom,i_drude,j_drude,sign_chg_i,sign_chg_j,flag_drudei,flag_drudej, flag_same_atom)

    ! we only consider intra molecular interactions between drude oscillators, or the corresponding drude oscillator charge on an atom
    ! also, we do not consider intra molecular interactions between a drude oscillator and it's corresponding atom

    ! flag_drudei tells us whether i_atom is either a drude oscillator, or has a drude_oscillator attached to it
    ! flag_drudej tells us whether j_atom is either a drude oscillator, or has a drude_oscillator attached to it
    ! flag_same_atom tells us whether j_atom is the drude oscillator on i_atom, or vice-versa

    if ( (flag_drudei == 1 ) .and. (flag_drudej == 1 ) .and. (flag_same_atom == 0 ) ) then
!!!!!!!!!!!!add screened atom-drude, get atom chg from corresponding oscillator
       pol1=chg(i_mole,i_drude)**2/springcon; pol2=chg(i_mole,j_drude)**2/springcon
       f_ij= f_ij + dble(sign_chg_i*sign_chg_j)*chg(i_mole,i_drude) * chg(i_mole,j_drude) * (rij * thole_screen(pol1,pol2,norm_dr,thole)/ norm_dr**3 - d_thole_screen(pol1,pol2,rij,thole) / norm_dr)
    endif



  end subroutine intra_pme_force



  !*****************************************************************
  ! This subroutine interpolates charges onto Q grid to be used in pme reciprocal space
  ! routines
  !*****************************************************************
  subroutine grid_Q(Q,chg,xyz,n_atom,n_mole,K,n)
    use global_variables
    use omp_lib
    integer, intent(in) :: n_mole
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: chg
    real*8, intent(in), dimension(:,:,:) :: xyz
    integer,intent(in)::K,n
    real*8,dimension(:,:,:),intent(out)::Q
    integer::i,j,k1,k2,k3,n1,n2,n3,nn1,nn2,nn3,nearpt(3)
    real*8::sum,chg_i
    real*8,dimension(3)::u,arg
    real*8 :: small=1D-6

    !test
!!$    integer :: p
!!$    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3), save :: xyz_store
    !test

    Q=0D0

    ! parameter spline_grid undeclared, but ok
    call OMP_SET_NUM_THREADS(n_threads)
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(n_threads,n_mole,n_atom,xyz,chg,n,B6_spline,B4_spline,K,spline_order,small) REDUCTION(+:Q)
    !$OMP DO SCHEDULE(dynamic, n_threads)
    !$
    do i=1,n_mole
       do j=1,n_atom(i)
          chg_i=chg(i,j)
          if ( abs(chg_i) > small ) then
             u=xyz(i,j,:)
             nearpt=floor(u)
             do k1=0,n
                n1=nearpt(1)-k1
                arg(1)=u(1)-dble(n1);
                do k2=0,n
                   n2=nearpt(2)-k2
                   arg(2)=u(2)-dble(n2)
                   do k3=0,n
                      n3=nearpt(3)-k3
                      arg(3)=u(3)-dble(n3);
!!$                 sum=chg_i*B_spline(arg(1),n)*B_spline(arg(2),n)*B_spline(arg(3),n)
                      if((arg(1)<dble(n)).and.(0.<arg(1)).and.(arg(2)<dble(n)).and.(0.<arg(2)).and.(arg(3)<dble(n)).and.(0.<arg(3))) then 
                         if(spline_order .eq.6) then
                            sum=chg_i*B6_spline(ceiling(arg(1)/6.D0*dble(spline_grid)))*B6_spline(ceiling(arg(2)/6.D0*dble(spline_grid)))*B6_spline(ceiling(arg(3)/6.D0*dble(spline_grid)))
                         else
                            sum=chg_i*B4_spline(ceiling(arg(1)/4.D0*dble(spline_grid)))*B4_spline(ceiling(arg(2)/4.D0*dble(spline_grid)))*B4_spline(ceiling(arg(3)/4.D0*dble(spline_grid)))   
                         endif
!!$                        sum=chg_i*interp_B_spline(arg(1),n)*interp_B_spline(arg(2),n)*interp_B_spline(arg(3),n)
                      else
                         sum=0.D0
                      endif
                      if(n1<0) then
                         nn1=n1+K
                      else
                         nn1=n1
                      endif
                      if(n2<0) then
                         nn2=n2+K
                      else
                         nn2=n2
                      endif
                      if(n3<0) then
                         nn3=n3+K
                      else
                         nn3=n3
                      endif

                      !test
!!$                      if ( ( nn1 < 0 ) .or. ( nn2 < 0 ) .or. ( nn3 < 0 ) ) then
!!$                         write(*,*) "bad coordinates"
!!$                         do p=1,n_atom(i)
!!$                            write(*,*) p, xyz_store(i,p,:)
!!$                         enddo
!!$                         stop
!!$                      end if
                      !test

                      Q(nn1+1,nn2+1,nn3+1)=Q(nn1+1,nn2+1,nn3+1)+sum
                   enddo
                enddo
             enddo
          end if
       enddo
    enddo

    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    !test
!!$    do i=1,n_mole
!!$       do j=1,n_atom(i)
!!$          xyz_store(i,j,:) = xyz(i,j,:) / dble(K) * 64.0383d0
!!$        enddo
!!$     enddo
    ! end test

  end subroutine grid_Q




!!!!!!!!!!!!!!!!!!!!!!!!
  ! this subroutine computes derivatives of grid Q with respect 
  ! to changes in position of i_atom in i_mole
  ! output is the force, which is computed by taking the product of
  ! the derivative Q array with FQ, which is the convoluted array
  !
  !
  ! notice that this is being called in a parallel section of code
!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine derivative_grid_Q(force,FQ,chg,xyz,i_atom,i_mole,K,box,n,kk)
    use global_variables
    integer,intent(in)::K,n
    real*8,dimension(3),intent(out)::force
    real*8,dimension(:,:,:),intent(in)::FQ
    integer, intent(in) :: i_mole,i_atom
    real*8, intent(in), dimension(:,:) :: chg
    real*8,intent(in),dimension(3,3)::box,kk
    real*8, intent(in), dimension(:,:,:) :: xyz
    integer::i,j,k1,k2,k3,n1,n2,n3,nn1,nn2,nn3,nearpt(3)
    integer::g1n(3),g1nmin(3),g1(3),g2(3)
    real*8::fac(3),chg_i
    real*8,dimension(3)::u,arg1,arg2,force_temp,box_length


    force=0.D0
    chg_i=chg(i_mole,i_atom)
    u=xyz(i_mole,i_atom,:)
    nearpt=floor(u)

       !!!!!!!!get bspline values from grid, therefore must check whether to see if function is being evaluated in non-zero domain, otherwise, can't call array
    !do k1=0,n
    ! if k1=n, arg1 >= n, so we can exclude that term in the sum, also because we exlude that term in the sum, arg1 <= n , so we dont need to check the limits of the bspline evaluation
    do k1=0,n-1
       n1=nearpt(1)-k1
       arg1(1)=u(1)-dble(n1);
       arg2(1) = arg1(1) - 1d0
       ! note arg1 > 0, so don't check for that.  for k1=n, arg1 > n
       !if((arg1(1)<real(n)).and.(0d0<arg1(1))) then
          !do k2=0,n
          do k2=0,n-1
             n2=nearpt(2)-k2
             arg1(2)=u(2)-dble(n2)
             arg2(2) = arg1(2) - 1d0
             !if((arg1(2)<real(n)).and.(0d0<arg1(2))) then
                !do k3=0,n
                 do k3=0,n-1
                   n3=nearpt(3)-k3
                   arg1(3)=u(3)-dble(n3);
                   arg2(3) = arg1(3) - 1d0
                   !if((arg1(3)<real(n)).and.(0d0<arg1(3))) then

!!!!!force in x direction
                      fac=0d0

                      ! don't use ceiling here, int() is faster and gives the same result in limit of converged grid
                      !g2=ceiling(arg2/real(n-1)*real(spline_grid))
                      !g1n=ceiling(arg1/real(n)*real(spline_grid))
                      !g1nmin = ceiling(arg1/real(n-1)*real(spline_grid))
                      g2=int(arg2/real(n-1)*real(spline_grid))
                      g1n=int(arg1/real(n)*real(spline_grid))
                      g1nmin = int(arg1/real(n-1)*real(spline_grid))

                      g1=g1n;g1(1)=g1nmin(1);

                      if(arg1(1)<real(n-1)) then 
                         if(spline_order.eq.6) then
                            fac(1)=chg_i*(B5_spline(g1(1))*B6_spline(g1(2))*B6_spline(g1(3)))
                         else
                            fac(1)=chg_i*(B3_spline(g1(1))*B4_spline(g1(2))*B4_spline(g1(3)))  
                         endif
                      endif
                      ! at this point, arg1(1) < n , so arg2(1) < n - 1 , so don't check for that
                      !if( (arg2(1)<real(n-1) ) .and.(0.<arg2(1)) ) then 
                      if(0.<arg2(1)) then
                         if(spline_order.eq.6) then
                            fac(1)=fac(1)+chg_i*(-B5_spline(g2(1))*B6_spline(g1(2))*B6_spline(g1(3)))
                         else
                            fac(1)=fac(1)+chg_i*(-B3_spline(g2(1))*B4_spline(g1(2))*B4_spline(g1(3)))  
                         endif
                      endif
!!!!!!force in y direction
                      g1=g1n;g1(2)=g1nmin(2)
                      if ( arg1(2)<real(n-1) ) then 
                         if(spline_order.eq.6) then
                            fac(2)=chg_i*(B5_spline(g1(2))*B6_spline(g1(1))*B6_spline(g1(3)))
                         else
                            fac(2)=chg_i*(B3_spline(g1(2))*B4_spline(g1(1))*B4_spline(g1(3)))  
                         endif
                      endif
                      if (0.<arg2(2)) then 
                         if(spline_order.eq.6) then
                            fac(2)=fac(2)+chg_i*(-B5_spline(g2(2))*B6_spline(g1(1))*B6_spline(g1(3)))
                         else
                            fac(2)=fac(2)+chg_i*(-B3_spline(g2(2))*B4_spline(g1(1))*B4_spline(g1(3)))  
                         endif
                      endif
!!!!!force in z direction
                      g1=g1n;g1(3)=g1nmin(3)
                      if(arg1(3)<real(n-1)) then 
                         if(spline_order.eq.6) then
                            fac(3)=chg_i*(B5_spline(g1(3))*B6_spline(g1(1))*B6_spline(g1(2)))
                         else
                            fac(3)=chg_i*(B3_spline(g1(3))*B4_spline(g1(1))*B4_spline(g1(2)))  
                         endif
                      endif
                      if (0.<arg2(3))  then 
                         if(spline_order.eq.6) then
                            fac(3)=fac(3)+chg_i*(-B5_spline(g2(3))*B6_spline(g1(1))*B6_spline(g1(2)))
                         else
                            fac(3)=fac(3)+chg_i*(-B3_spline(g2(3))*B4_spline(g1(1))*B4_spline(g1(2)))  
                         endif
                      endif


                      if(n1<0) then
                         nn1=n1+K
                      else
                         nn1=n1
                      endif
                      if(n2<0) then
                         nn2=n2+K
                      else
                         nn2=n2
                      endif
                      if(n3<0) then
                         nn3=n3+K
                      else
                         nn3=n3
                      endif
                      force=force + fac * FQ(nn1+1,nn2+1,nn3+1)
             !  dQdr(nn1+1,nn2+1,nn3+1,:)=fac(:)

!                   end if
                enddo
!              end if
          enddo
!          end if
    enddo


    !****************change to global coordinates*****************
    ! for a general box.  The array force contains derivatives with respect to the
    ! scaled grid coordinates.  We need to convert this to derivatives with respect to cartesian coords
    force_temp=0d0
    do i=1,3
       do j=1,3
          force_temp(i) = force_temp(i) - dble(K) * kk(j,i) * force(j)
       enddo
    enddo

    force = force_temp



  end subroutine derivative_grid_Q





  !***********************************************
  ! this function calculates B_splines which are used in pme as interpolating
  ! functions.  B_splines are calculated recursively, and therefore it's a good idea
  ! to grid them
  !************************************************
  real*8 function B_spline(u,n)
    real*8,intent(in)::u
    integer,intent(in)::n
    integer::i,j
    real,dimension(n-1,n-1)::mn
    real*8::ui

!!!!!!!! define m2 for n-1 values
    do i=1,n-1
       ui=u-dble(i-1)
       if((ui<0.).or.(ui>2.)) then
          mn(1,i)=0.D0
       else
          mn(1,i)=1.D0-abs(ui-1.D0)
       endif
    enddo

!!!!!!!!!!! define mj recursively for n-1-(j-1) values

    do j=2,n-1
       do i=1,n-j
          ui=u-dble(i-1)
          mn(j,i)=(ui/dble(j))*mn(j-1,i)+((dble(j+1)-ui)/dble(j))*mn(j-1,i+1)
       enddo
    enddo

    B_spline=mn(n-1,1)

  end function B_spline


  !***********************************************************
  ! This is a routine for reciprocal space pme calculation
  !***********************************************************
  subroutine CB_array(CB,alpha_sqrt,vol,K,kk,n)
    real*8,dimension(:,:,:),intent(out)::CB
    real*8,intent(in)::alpha_sqrt,vol
    real*8,dimension(3,3),intent(in)::kk
    integer,intent(in)::K,n
    real*8, parameter :: pi=3.14159265
    real*8,dimension(3)::mm
    real*8::mag
    integer::i,j,l,m1,m2,m3
    do i=0,K-1
       if(i>K/2) then
          m1=i-K
       else
          m1=i
       endif
       do j=0,K-1
          if(j>K/2) then
             m2=j-K
          else
             m2=j
          endif
          do l=0,K-1
             if(l>K/2) then
                m3=l-K
             else
                m3=l
             endif
             mm(:)=m1*kk(1,:)+m2*kk(2,:)+m3*kk(3,:)
             mag=dot_product(mm,mm)
             CB(i+1,j+1,l+1)=1./(vol*pi)*exp(-pi**2*mag/alpha_sqrt**2)/mag*bm_sq(i,n,K)*bm_sq(j,n,K)*bm_sq(l,n,K)
          enddo
       enddo
    enddo

    CB(1,1,1)=0.D0

  end subroutine CB_array


  !******************************************************
  ! this is needed in reciprocal space pme calculation
  !******************************************************
  function bm_sq(m,n,K)
    use global_variables
    real*8::bm_sq
    integer,intent(in)::m,n,K
    integer::i
    complex*16::bm,sum
    real*8::tmp

    sum=0.D0
    do i=0,n-2
       tmp=2.D0*pi*dble(m*i)/dble(K)
       sum=sum+B_spline(dble(i+1),n)*cmplx(cos(tmp),sin(tmp))
!!$     sum=sum+B6_spline(dble(i+1)/6.*dble(spline_grid))*cmplx(cos(tmp),sin(tmp))
    enddo


    bm=1.D0/sum
    bm_sq=dble(bm)**2+aimag(bm)**2

  end function bm_sq



  !********************************************************
  ! this subroutine calls particular type of electrostatic force calculation and returns force
  ! since force calculations are used for drude oscillators, and ewald is way too slow for this,
  ! ewald is not in here.
  !*********************************************************
  subroutine Electrostatic_energy_force( electrostatic_energy, force_atoms,target_atoms,tot_n_mole, n_mole, n_atom, xyz, chg, box, r_com,dfti_desc,dfti_desc_inv)
    use global_variables
    use MKL_DFTI
    use electrostatic
    implicit none
    real*8, intent(out)               :: electrostatic_energy
    real*8,dimension(:,:,:),intent(out)::force_atoms
    integer,dimension(:,:),intent(in) :: target_atoms
    integer, intent(in) :: tot_n_mole, n_mole
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: chg, box,r_com
    real*8, intent(in), dimension(:,:,:) :: xyz
    TYPE(DFTI_DESCRIPTOR),pointer,intent(in):: dfti_desc,dfti_desc_inv

    Select Case(electrostatic_type)
    Case("pme")
       call pme_energy_force(electrostatic_energy, force_atoms,target_atoms,tot_n_mole,n_mole, n_atom, xyz, chg, box, r_com, dfti_desc,dfti_desc_inv)
    Case("cutoff")
       call F_elec_cutoff(force_atoms ,target_atoms,tot_n_mole, n_mole, n_atom, chg, xyz, box, r_com) 
       electrostatic_energy = E_elec_cutoff( tot_n_mole ,n_mole, n_atom, chg, xyz, box, r_com )
    Case("none")
       force_atoms=0d0
       electrostatic_energy=0d0
    Case Default
       stop "requested electrostatic type is not recognized/implemented to compute forces"
    end Select

  end subroutine Electrostatic_energy_force





  !***************************** subroutines for PME dispersion grid calculations *****************************************
  ! written by kuang and modified by jesse 7/13/2012 to be general for non-orthogonal boxes and to not use link lists
  ! modified again 9/27/2012 to include C12 terms
  ! modified 12/6/2012 to allow for pme dispersion treatment of Feynmann-Hibbs correction
  !************************************************************************************************************************


  !****************************************************************************
  ! This subroutine construct the reciprocal PME potenatial grids
  ! Major output:
  ! real*8, dimension(:,:,:) :: recip_pot ! the reciprocal PME potential grid
  ! Note only works for cav_grid_a = cav_grid_b = cav_grid_c
  ! It can be generalized to the dispersion situation, where, U ~ 1/r^p
  ! In that case, chg should be Cp(ij), the dispersion coefficients
  !****************************************************************************
  subroutine recip_pme_grid( box,tot_n_mole,n_mole,n_atom,xyz,chg,dfti_desc,dfti_desc_inv,recip_pot,p)
    use global_variables
    use routines
    use MKL_DFTI
    implicit none 
    real*8, dimension(:,:), intent(in) :: box
    integer, intent(in) :: tot_n_mole, n_mole, p
    integer, dimension(:), intent(in) :: n_atom
    real*8, dimension(:,:,:), intent(in) :: xyz
    real*8, dimension(:,:), intent(in) :: chg
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv
    real*8, dimension(:,:,:), intent(out) :: recip_pot
    !   local variables
    real*8, dimension(pme_grid,pme_grid,pme_grid) :: Q
    complex*16,dimension(pme_grid,pme_grid,pme_grid) :: FQ, CBp
    integer :: n,K, i_mole, i_atom, idim, status, m1,m2,m3
    real*8 :: vol,a(3),b(3),c(3),ka(3),kb(3),kc(3),kk(3,3),mag,mm(3)
    real*8, dimension(:,:,:), allocatable :: xyz_scale
    real*8,dimension(pme_grid**3)::q_1r
    complex*16,dimension(pme_grid**3)::q_1d

    n=spline_order
    !   This pme calculation is controled by pme_grid
    K=pme_grid


    a(:) = box(1,:)
    b(:) = box(2,:)
    c(:) = box(3,:)

    ! calculate the volume and the reciprocal vectors (notice no 2pi)
    vol = volume( a, b, c )
    call crossproduct( a, b, kc ); kc = kc /vol
    call crossproduct( b, c, ka ); ka = ka /vol
    call crossproduct( c, a, kb ); kb = kb /vol
    kk(1,:)=ka(:);kk(2,:)=kb(:);kk(3,:)=kc(:)

    ! create scaled coordinates
    allocate(xyz_scale(tot_n_mole-n_mole,maxval(n_atom(n_mole+1:)),3))

    !   only loops over frame atoms
    do i_mole = n_mole+1, tot_n_mole
       do i_atom = 1, n_atom(i_mole)
          do idim = 1, 3
             xyz_scale(i_mole-n_mole,i_atom,idim) = dble(K) * dot_product(kk(idim,:),xyz(i_mole,i_atom,:))
             ! if atoms are out of grid, shift them back in
             if ( xyz_scale(i_mole-n_mole,i_atom,idim) < 0. ) then
                xyz_scale(i_mole-n_mole,i_atom,idim) = xyz_scale(i_mole-n_mole,i_atom,idim) + dble(K)
                !else if( xyz_scale(i_mole-n_mole,i_atom,idim) > dble(K) ) then ! changed by Jesse, 6/10/2014 see below line.  In rare cases, this could lead to out of bounds array calls when constructing Q
             else if( xyz_scale(i_mole-n_mole,i_atom,idim) >= dble(K) ) then
                xyz_scale(i_mole-n_mole,i_atom,idim) = xyz_scale(i_mole-n_mole,i_atom,idim) - dble(K)
             endif
          enddo
       enddo
    enddo

    !   Fourier transform of Q grid
    call grid_Q(Q,chg(n_mole+1:,:),xyz_scale,n_atom(n_mole+1:),tot_n_mole-n_mole,K,n)
    q_1r = RESHAPE(Q, (/K**3/) )
    q_1d = cmplx(q_1r,0.,16)
    status=DftiComputeForward(dfti_desc, q_1d)
    FQ = RESHAPE(q_1d, (/K,K,K/) )
    !   CB' matrixm, debug
    call CB_prime(CBp,alpha_sqrt,vol,K,kk,n,p)
    FQ = FQ * CBp
    !   inverse transform
    q_1d = RESHAPE(FQ, (/K**3/))
    status = DftiComputeBackward(dfti_desc_inv, q_1d)
    FQ=RESHAPE(q_1d, (/K,K,K/) )

    recip_pot = recip_pot + dble(FQ)

    return
  end subroutine recip_pme_grid


  !***************************************************************
  ! Modified from bm_sq, this is used to calculate potential grid
  !***************************************************************
  function bm(m,n,K)
    use global_variables
    integer,intent(in)::m,n,K
    integer::i
    complex*16::bm,sum
    real*8::tmp
    sum=0.D0
    do i=0,n-2
       tmp=2.D0*pi*dble(m*i)/dble(K)
       sum=sum+cmplx(B_spline(dble(i+1),n),0.0d0,16)*cmplx(cos(tmp),sin(tmp),16)
!!$     sum=sum+B6_spline(dble(i+1)/6.*dble(spline_grid))*cmplx(cos(tmp),sin(tmp))
    enddo
    tmp = 2.0d0*pi*dble((n-1)*m)/dble(K)
    bm = cmplx(cos(tmp),sin(tmp),16)/sum
  end function bm

  !**********************************************************************
  ! Here, the complex CB' matrix is:
  ! CB' = 1/pi/V*exp(-pi**2*m**2/alpha)/m**2 * b1(m1)*b2(m2)*b3(m3)
  ! compared to normal CB array, the B_spline functions are not squared
  !**********************************************************************
  subroutine CB_prime(CBp,alpha_sqrt,vol,K,kk,n,p)
    complex*16,dimension(:,:,:),intent(out)::CBp
    real*8,intent(in)::alpha_sqrt,vol
    real*8,dimension(3,3),intent(in)::kk
    integer,intent(in)::K,n,p
    real*8,dimension(3)::mm
    real*8, parameter :: pi=3.14159265
    real*8::mag, x, pi_sqrt, C12_temp
    integer::i,j,l,m1,m2,m3,split_do

    pi_sqrt= dsqrt(pi)

    CBp = 0.d0

    do i=0,K-1
       if(i>K/2) then
          m1=i-K
       else
          m1=i
       endif
       do j=0,K-1
          if(j>K/2) then
             m2=j-K
          else
             m2=j
          endif
          do l=0,K-1
             if(l>K/2) then
                m3=l-K
             else
                m3=l
             endif
             mm(:)=m1*kk(1,:)+m2*kk(2,:)+m3*kk(3,:)
             mag=dot_product(mm,mm)
             x = pi*dsqrt(mag)/alpha_sqrt

             !            Note, due to the different definition of FFT, we need the complex conjugate here.
             select case(p)
             case(1)
                CBp(i+1,j+1,l+1) = pi_sqrt*pi/vol/alpha_sqrt**2*exp(-x**2)/pi_sqrt/x**2*conjg( bm(i,n,K)*bm(j,n,K)*bm(l,n,K) )
                !               CBp(i+1,j+1,l+1) = 1.d0/(vol*pi)*exp(-pi**2*mag/alpha_sqrt**2)/mag*conjg( bm(i,n,K)*bm(j,n,K)*bm(l,n,K) )
             case(6)
                CBp(i+1,j+1,l+1) = pi_sqrt*pi/vol*alpha_sqrt**3* &
                     & ( (1d0-2d0*x**2)*exp(-x**2) + 2d0*pi_sqrt*x**3*erfc_g(x) ) / 3d0 &
                     & *conjg( bm(i,n,K)*bm(j,n,K)*bm(l,n,K) )
             case(8)
                CBp(i+1,j+1,l+1) = pi_sqrt*pi/vol*alpha_sqrt**5* &
                     & ( (3d0-2d0*x**2+4d0*x**4)*exp(-x**2) - 4d0*pi_sqrt*x**5*erfc_g(x) ) / 45d0 &
                     & *conjg( bm(i,n,K)*bm(j,n,K)*bm(l,n,K) )
             case(10)
                CBp(i+1,j+1,l+1) = pi_sqrt*pi/vol*alpha_sqrt**7* &
                     & ( (15d0-6d0*x**2+4d0*x**4-8d0*x**6)*exp(-x**2) + 8d0*pi_sqrt*x**7*erfc_g(x) ) / 1260d0 &
                     & *conjg( bm(i,n,K)*bm(j,n,K)*bm(l,n,K) )
             case(12)
                CBp(i+1,j+1,l+1) = pi_sqrt*pi/vol*alpha_sqrt**9* &
                     & ( (105d0-30d0*x**2+12d0*x**4-8d0*x**6+16d0*x**8)*exp(-x**2) - 16d0*pi_sqrt*x**9*erfc_g(x) ) / 56700d0 &
                     & *conjg( bm(i,n,K)*bm(j,n,K)*bm(l,n,K) )
             case(14)
                C12_temp=  ( (105d0-30d0*x**2+12d0*x**4-8d0*x**6+16d0*x**8)*exp(-x**2) - 16d0*pi_sqrt*x**9*erfc_g(x) ) / 56700d0 
                CBp(i+1,j+1,l+1) = pi_sqrt*pi/vol*alpha_sqrt**11 * ( 1/3960d0 * exp(-x**2) - 1/33d0 * x**2 * C12_temp ) *conjg( bm(i,n,K)*bm(j,n,K)*bm(l,n,K) )

             case default
                write(*,*) 'Error: Interactions other than r^-1, -6, -8, -10 -12 -14 are not supported!'
                stop
             end select
             !             CBp(i+1,j+1,l+1) = conjg( bm(i,n,K)*bm(j,n,K)*bm(l,n,K) )
          enddo
       enddo
    enddo


    if ( p == 1 ) then
       CBp(1,1,1) = cmplx(0.0d0,0.0d0,16)
    endif

  end subroutine CB_prime


  !**********************************************************************************************************************
  ! This subroutine construct the long-range dispersion interaction grid for a certain type of atom, which is specified
  ! by the variable atom_id1. It constructs the pseudo "charge" array based on dispersion coefficients, and feed it to 
  ! the reciprocal PME grid subroutine to generate the dispersion interaction grid.
  ! 
  ! if employing Feynmann-Hibbs correction, then derivatives of lower order terms get absorbed into higher order terms, to
  ! give "effective" C8,C10 , etc.
  !**********************************************************************************************************************
  subroutine construct_recip_disp_grid(atom_id1,n_mole,tot_n_mole,n_atom,box,xyz,dfti_desc,dfti_desc_inv)
    use global_variables
    use MKL_DFTI
    implicit none
    integer, intent(in) :: atom_id1, n_mole, tot_n_mole
    integer, dimension(:), intent(in) :: n_atom
    real*8, dimension(:,:), intent(in) :: box
    real*8, dimension(:,:,:), intent(in) :: xyz
    TYPE(DFTI_DESCRIPTOR),pointer,intent(in):: dfti_desc,dfti_desc_inv
    !   locale variables
    integer :: i_mole, atom_id2, i_atom , i_index, loc_index
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM) :: cpij
    real*8 :: mass_reduced, fac,mass_id1,mass_id2, C6_local
    real*8,parameter :: hbar_mole_sq=0.403280544d0  ! this is (hbar*NA)^2 ( in units of Kg^2 m^2*A^2/s^2 )

    ! if we are using Feynmann-Hibbs correction, we need mass of atomtype "atom_id1".  If this is a dummy atom, the mass
    ! is zero, and the correct thing to do is use the entire molecule mass if the dummy is at the center of the molecule
    ! usually we would call generate_reduced_mass to do this, but since we don't know what molecule this atom type is in, we cant call this routine
    ! figure it out here
    Select Case(Feynmann_Hibbs_correction)
    Case("yes")
       ! if atom mass is 0, figure out correct mass to use
       if ( atype_mass(atom_id1) < 0.01 ) then
          i_index=-1
          ! find a molecule containing this atom type
          do i_mole=1, n_mole
             do i_atom=1, n_atom(i_mole)
                if ( atom_index(i_mole,i_atom) == atom_id1 ) then
                   i_index = i_mole
                endif
             enddo
             if ( i_index > -1 ) exit
          enddo
          if ( i_index == -1 ) then
             stop "error in recognizing atom type in subroutine construct_recip_disp_grid "
          endif
          mass_id1=0d0
          do i_atom=1,n_atom(i_index)
             loc_index = atom_index(i_index,i_atom)
             mass_id1=mass_id1 + atype_mass(loc_index)
          enddo
       else
          mass_id1 = atype_mass(atom_id1)
       endif
    End Select

    !************************* now main part of subroutine ********************************


    lrdisp_pot(atom_id1,:,:,:) = 0.0d0
    do i_mole = 1, tot_n_mole
       if ( i_mole <= n_mole ) then ! solute molecules
          do i_atom = 1, n_atom(i_mole)
             cpij(i_mole,i_atom) = 1.0d0
          enddo
       else
          do i_atom = 1, n_atom(i_mole)
             atom_id2 = atom_index(i_mole,i_atom)
             if ( lj_bkghm_index(atom_id1,atom_id2) == 1 ) then
                cpij(i_mole,i_atom) = -1d0*atype_lj_parameter(atom_id1,atom_id2,3)
             else if ( lj_bkghm_index(atom_id1,atom_id2) == 2 ) then
                !       cpij(i_mole,i_atom) = -4d0*atype_lj_parameter(atom_id1,atom_id2,1)*atype_lj_parameter(atom_id1,atom_id2,2)**6
                ! working in C12, C6 now, not epsilon, sigma
                cpij(i_mole,i_atom) =-atype_lj_parameter(atom_id1,atom_id2,2)

             endif
          enddo
       endif
    enddo


    call recip_pme_grid( box,tot_n_mole,n_mole,n_atom,xyz,cpij,dfti_desc,dfti_desc_inv,lrdisp_pot(atom_id1,:,:,:),6)

    !    for either lj potential, or bkghm-6 potential
    !    add C8 terms for Feynmann_Hibbs correction even if not using C8_C10_dispersion
    ! note that this if statement should be fine for hybrid lj/bkghm simulation, i.e. lj_bkghm=3, because in that case we have to use C8_10_dispersion_terms='yes' 
    if (  ( C8_10_dispersion_terms == 'no' ) .or. ( lj_bkghm == 2 ) )  then
       Select Case(Feynmann_Hibbs_correction)
       Case("yes")
          do i_mole = 1, tot_n_mole
             if ( i_mole <= n_mole ) then ! solute molecules
                do i_atom = 1, n_atom(i_mole)
                   cpij(i_mole,i_atom) = 1.0d0
                enddo
             else
                do i_atom = 1, n_atom(i_mole)
                   atom_id2 = atom_index(i_mole,i_atom)
                   mass_id2 = atype_mass(atom_id2)
                   if ( mass_id2 < 0.01 ) then
                      stop "mass of zero on framework atoms in subroutine construct_recip_disp_grid "
                   endif
                   mass_reduced = mass_id1 * mass_id2 / (mass_id1 + mass_id2)
                   ! fac = hbar^2/(24*kT * mass_red ) ; since distance is in angstroms, this should be in A^2
                   ! use (hbar*mol)^2 in Kg^2*m^2*A^2/s^2, then k in J/K/mol, mass in Kg
                   fac= hbar_mole_sq / ( 24d0 * 8.31446d0 * temp * mass_reduced / 1000d0 )
                   ! C6 contribution to C8 effective potential
                   Select Case( lj_bkghm)
                   Case(1)
                      C6_local = atype_lj_parameter(atom_id1,atom_id2,3)
                   Case(2)
                      C6_local = atype_lj_parameter(atom_id1,atom_id2,2)
                   End Select

                   cpij(i_mole,i_atom) = -1d0 * fac * C6_local * ( 42d0 - 12d0 )   ! for clarity d2U/dr2 contribution, then 2/r*dU/dr contribution
                enddo
             endif
          enddo
          call recip_pme_grid( box,tot_n_mole,n_mole,n_atom,xyz,cpij,dfti_desc,dfti_desc_inv,lrdisp_pot(atom_id1,:,:,:),8)
       End Select
    endif


    if ( C8_10_dispersion_terms == 'yes' ) then
       !     C8 terms
       do i_mole = 1, tot_n_mole
          if ( i_mole <= n_mole ) then ! solute molecules
             do i_atom = 1, n_atom(i_mole)
                cpij(i_mole,i_atom) = 1.0d0
             enddo
          else
             do i_atom = 1, n_atom(i_mole)
                atom_id2 = atom_index(i_mole,i_atom)
                cpij(i_mole,i_atom) = -1d0*atype_lj_parameter(atom_id1,atom_id2,4)
                ! if Feynmann-Hibbs correction
                Select Case(Feynmann_Hibbs_correction)
                Case("yes")
                   mass_id2 = atype_mass(atom_id2)
                   if ( mass_id2 < 0.01 ) then
                      stop "mass of zero on framework atoms in subroutine construct_recip_disp_grid "
                   endif
                   mass_reduced = mass_id1 * mass_id2 / (mass_id1 + mass_id2)
                   ! fac = hbar^2/(24*kT * mass_red ) ; since distance is in angstroms, this should be in A^2
                   ! use (hbar*mol)^2 in Kg^2*m^2*A^2/s^2, then k in J/K/mol, mass in Kg
                   fac= hbar_mole_sq / ( 24d0 * 8.31446d0 * temp * mass_reduced / 1000d0 )
                   ! C6 contribution to C8 effective potential
                   cpij(i_mole,i_atom) = cpij(i_mole,i_atom) - 1d0 * fac * atype_lj_parameter(atom_id1,atom_id2,3) * ( 42d0 - 12d0 )   ! for clarity d2U/dr2 contribution, then 2/r*dU/dr contribution
                End Select
             enddo
          endif
       enddo
       call recip_pme_grid( box,tot_n_mole,n_mole,n_atom,xyz,cpij,dfti_desc,dfti_desc_inv,lrdisp_pot(atom_id1,:,:,:),8)
       !     C10 terms
       do i_mole = 1, tot_n_mole
          if ( i_mole <= n_mole ) then ! solute molecules
             do i_atom = 1, n_atom(i_mole)
                cpij(i_mole,i_atom) = 1d0
             enddo
          else
             do i_atom = 1, n_atom(i_mole)
                atom_id2 = atom_index(i_mole,i_atom)
                cpij(i_mole,i_atom) = -1d0*atype_lj_parameter(atom_id1,atom_id2,5)
                ! if Feynmann-Hibbs correction
                Select Case(Feynmann_Hibbs_correction)
                Case("yes")
                   mass_id2 = atype_mass(atom_id2)
                   mass_reduced = mass_id1 * mass_id2 / (mass_id1 + mass_id2)
                   fac= hbar_mole_sq / ( 24d0 * 8.31446d0 * temp * mass_reduced / 1000d0 )
                   ! C8 contribution to C10 effective potential
                   cpij(i_mole,i_atom) = cpij(i_mole,i_atom) - 1d0 * fac * atype_lj_parameter(atom_id1,atom_id2,4) * ( 72d0 - 16d0 )   ! for clarity d2U/dr2 contribution, then 2/r*dU/dr contribution
                End Select
             enddo
          endif
       enddo
       call recip_pme_grid( box,tot_n_mole,n_mole,n_atom,xyz,cpij,dfti_desc,dfti_desc_inv,lrdisp_pot(atom_id1,:,:,:),10)

       !       if we are using C12 terms, include them
       Select Case(C12_dispersion)
       Case("yes")
          !     C12 terms
          do i_mole = 1, tot_n_mole
             if ( i_mole <= n_mole ) then ! solute molecules
                do i_atom = 1, n_atom(i_mole)
                   cpij(i_mole,i_atom) = 1d0
                enddo
             else
                do i_atom = 1, n_atom(i_mole)
                   atom_id2 = atom_index(i_mole,i_atom)
                   cpij(i_mole,i_atom) = -1d0*atype_lj_parameter(atom_id1,atom_id2,6)
                   ! if Feynmann-Hibbs correction
                   Select Case(Feynmann_Hibbs_correction)
                   Case("yes")
                      mass_id2 = atype_mass(atom_id2)
                      mass_reduced = mass_id1 * mass_id2 / (mass_id1 + mass_id2)
                      fac= hbar_mole_sq / ( 24d0 * 8.31446d0 * temp * mass_reduced / 1000d0 )
                      ! C10 contribution to C12 effective potential
                      cpij(i_mole,i_atom) = cpij(i_mole,i_atom) - 1d0 * fac * atype_lj_parameter(atom_id1,atom_id2,5) * ( 110d0 - 20d0 )   ! for clarity d2U/dr2 contribution, then 2/r*dU/dr contribution
                   End Select
                enddo
             endif
          enddo
          call recip_pme_grid( box,tot_n_mole,n_mole,n_atom,xyz,cpij,dfti_desc,dfti_desc_inv,lrdisp_pot(atom_id1,:,:,:),12)


          ! C14 contribution from Feynmann_Hibbs correction, if we are using C12_dispersion='yes'
          Select Case(Feynmann_Hibbs_correction)
          Case("yes")
             do i_mole = 1, tot_n_mole
                if ( i_mole <= n_mole ) then ! solute molecules
                   do i_atom = 1, n_atom(i_mole)
                      cpij(i_mole,i_atom) = 1d0
                   enddo
                else
                   do i_atom = 1, n_atom(i_mole)
                      atom_id2 = atom_index(i_mole,i_atom)
                      mass_id2 = atype_mass(atom_id2)
                      mass_reduced = mass_id1 * mass_id2 / (mass_id1 + mass_id2)
                      fac= hbar_mole_sq / ( 24d0 * 8.31446d0 * temp * mass_reduced / 1000d0 )
                      ! C12 contribution to C14 effective potential
                      cpij(i_mole,i_atom) = - 1d0 * fac * atype_lj_parameter(atom_id1,atom_id2,6) * ( 156d0 - 24d0 )   ! for clarity d2U/dr2 contribution, then 2/r*dU/dr contribution
                   enddo
                endif
             enddo
             call recip_pme_grid( box,tot_n_mole,n_mole,n_atom,xyz,cpij,dfti_desc,dfti_desc_inv,lrdisp_pot(atom_id1,:,:,:),14)
          End Select

          ! ************* this is C12 effective potential from C10 dispersion parameters, calculated for setting C12_dispersion= 'no'
       Case("no")
          ! even if we are not explicitly using C12 terms, could be contribution from Feynmann-Hibbs correction
          Select Case(Feynmann_Hibbs_correction)
          Case("yes")
             do i_mole = 1, tot_n_mole
                if ( i_mole <= n_mole ) then ! solute molecules
                   do i_atom = 1, n_atom(i_mole)
                      cpij(i_mole,i_atom) = 1d0
                   enddo
                else
                   do i_atom = 1, n_atom(i_mole)
                      atom_id2 = atom_index(i_mole,i_atom)
                      mass_id2 = atype_mass(atom_id2)
                      mass_reduced = mass_id1 * mass_id2 / (mass_id1 + mass_id2)
                      fac= hbar_mole_sq / ( 24d0 * 8.31446d0 * temp * mass_reduced / 1000d0 )
                      ! C10 contribution to C12 effective potential
                      cpij(i_mole,i_atom) = - 1d0 * fac * atype_lj_parameter(atom_id1,atom_id2,5) * ( 110d0 - 20d0 )   ! for clarity d2U/dr2 contribution, then 2/r*dU/dr contribution
                   enddo
                endif
             enddo
             call recip_pme_grid( box,tot_n_mole,n_mole,n_atom,xyz,cpij,dfti_desc,dfti_desc_inv,lrdisp_pot(atom_id1,:,:,:),12)
          End Select
       End Select

    endif

  end subroutine construct_recip_disp_grid



  function disp_usegrid( box,i_mole,tot_n_mole,n_mole,n_atom,r_com,xyz_config,mass)
    use global_variables
    use routines
    use pairwise_interaction
    implicit none
    integer, intent(in) :: tot_n_mole, n_mole, i_mole
    integer, dimension(:), intent(in) :: n_atom
    real*8, dimension(:,:,:), intent(in) :: xyz_config
    real*8, dimension(:,:), intent(in) :: box, r_com, mass
    real*8 :: disp_usegrid
    !   local variables
    integer :: i_atom, j_mole, j_atom, k_atom,tmp_mole, tmp_atom, atom_id1, atom_id2,local_index
    real*8 :: ra,rb,rc, dra,drb,drc,dr, ewald_cutoff_use, pme_real, E_short, E_long, E_local
    real*8 :: vol,a(3),b(3),c(3),ka(3), kb(3), kc(3),kk(3,3),rr(3)
    real*8, dimension(5) :: energy_decomp
    integer :: ia,ib,ic, ia0,ib0,ic0, igrida,igridb,igridc, split_do,l,K
    real*8, dimension(-lgrg_order+1:lgrg_order,-lgrg_order+1:lgrg_order,-lgrg_order+1:lgrg_order) :: weight
    real*8, dimension(3) :: r_com_i, r_com_j, shift, dr_com,rij
    real*8 :: mass_reduced, fac,mass_id1,mass_id2, C6_local
    real*8,parameter :: hbar_mole_sq=0.403280544d0  ! this is (hbar*NA)^2 ( in units of Kg^2 m^2*A^2/s^2 )

    !   effective cutoff for real space calculations
    !  ewald_cutoff_use = min(box(1,1)/2.,box(2,2)/2.,box(3,3)/2.,lj_cutoff)
    ewald_cutoff_use = ewald_cutoff


    K=pme_grid
    a(:) = box(1,:)
    b(:) = box(2,:)
    c(:) = box(3,:)
    ! calculate the volume and the reciprocal vectors (notice no 2pi)
    vol = volume( a, b, c )
    call crossproduct( a, b, kc ); kc = kc /vol 
    call crossproduct( b, c, ka ); ka = ka /vol
    call crossproduct( c, a, kb ); kb = kb /vol
    kk(1,:)=ka(:);kk(2,:)=kb(:);kk(3,:)=kc(:)


    disp_usegrid = 0.0d0
    do i_atom = 1, n_atom(i_mole)
       atom_id1 = atom_index(i_mole, i_atom)
       do l=1,3
          rr(l) = dble(K)*dot_product(kk(l,:),xyz_config(i_mole,i_atom,:))
          ! if atoms are out of grid, shift them back in
          if (rr(l)<0.) then
             rr(l)=rr(l)+dble(K)
          else if(rr(l)>dble(K)) then
             rr(l)=rr(l)-dble(K)
          endif
       enddo
       ra = rr(1) ; rb = rr(2) ; rc = rr(3)

       weight = 1.0d0
       ia0 = floor(ra)
       ib0 = floor(rb)
       ic0 = floor(rc)
       do ia = -lgrg_order+1, lgrg_order
          weight( ia,:,: ) = weight( ia,:,: ) * wfun( ia - (ra-ia0) )
       enddo
       do ib = -lgrg_order+1, lgrg_order
          weight( :,ib,: ) = weight( :,ib,: ) * wfun( ib - (rb-ib0) )
       enddo
       do ic = -lgrg_order+1, lgrg_order
          weight( :,:,ic ) = weight( :,:,ic ) * wfun( ic - (rc-ic0) )
       enddo
       !     interpolate the reciprocal potential grid
       do ia = -lgrg_order+1, lgrg_order
          igrida = ia0 + ia
          igrida = igrida - floor(dble(igrida)/dble(pme_grid))*pme_grid
          igrida = igrida + 1
          do ib = -lgrg_order+1, lgrg_order
             igridb = ib0 + ib
             igridb = igridb - floor(dble(igridb)/dble(pme_grid))*pme_grid
             igridb = igridb + 1
             do ic = -lgrg_order+1, lgrg_order
                igridc = ic0 + ic
                igridc = igridc - floor(dble(igridc)/dble(pme_grid))*pme_grid
                igridc = igridc + 1
                disp_usegrid = disp_usegrid + weight(ia,ib,ic)*lrdisp_pot(atom_id1,igrida,igridb,igridc)
             enddo
          enddo
       enddo


       !     real space interactions, only between solute and framework
       r_com_i = r_com(i_mole, :)
       pme_real = 0.0d0


       do j_mole=n_mole+1, tot_n_mole
          r_com_j = r_com(j_mole, :)
          shift(:) = pbc_shift( r_com_i, r_com_j, box, xyz_to_box_transform )
          dr_com = pbc_dr( r_com_i, r_com_j, shift )
          if ( abs(dr_com(1)) < ewald_cutoff_use .and. abs(dr_com(2)) < ewald_cutoff_use &
               &   .and. abs(dr_com(3)) < ewald_cutoff_use ) then
             if ( dot_product( dr_com, dr_com ) < ewald_cutoff_use**2 ) then
                do j_atom = 1, n_atom( j_mole )
                   rij = pbc_dr( xyz_config(i_mole,i_atom,:), xyz_config(j_mole,j_atom,:), shift )
                   dr = dsqrt( dot_product( rij, rij ) )
                   !                 only inter-molecular interaction is possible here.
                   atom_id2 = atom_index(j_mole,j_atom)
                   if ( lj_bkghm_index(atom_id1,atom_id2) == 1 ) then ! Buckingham potential
                      ! short range interactions
                      call eval_short_range_bkghm(E_short, energy_decomp, atom_id1,atom_id2, dr)
                      pme_real = pme_real + E_short
                      ! C6 term
                      pme_real = pme_real - ( gfun_g(alpha_sqrt*dr,6) - 1d0 + C6_C10_damp(atom_id1,atom_id2,dr,6) ) &
                           & * atype_lj_parameter(atom_id1,atom_id2,3)/dr**6
                      ! C8 term
                      pme_real = pme_real - ( gfun_g(alpha_sqrt*dr,8) - 1d0 + C6_C10_damp(atom_id1,atom_id2,dr,8) ) &
                           & *atype_lj_parameter(atom_id1,atom_id2,4)/dr**8
                      ! C10 term
                      pme_real = pme_real - ( gfun_g(alpha_sqrt*dr,10) - 1d0 + C6_C10_damp(atom_id1,atom_id2,dr,10) ) &
                           & *atype_lj_parameter(atom_id1,atom_id2,5)/dr**10
                      ! if using C12 terms
                      Select Case(C12_dispersion)
                      Case("yes")
                         pme_real = pme_real - ( gfun_g(alpha_sqrt*dr,12) - 1d0 + C6_C10_damp(atom_id1,atom_id2,dr,12) ) &
                              & *atype_lj_parameter(atom_id1,atom_id2,6)/dr**12
                      End Select


                   else if ( lj_bkghm_index(atom_id1,atom_id2) == 2 ) then ! Lennard-Jones potential
                      ! working in C12, C6 now, not epsilon, sigma
                      ! E_short = 4D0*atype_lj_parameter(atom_id1,atom_id2,1) *(atype_lj_parameter(atom_id1,atom_id2,2)/dr)**12

                      E_short = atype_lj_parameter(atom_id1,atom_id2,1)/dr**12
                      pme_real = pme_real + E_short
                      !              pme_real = pme_real - gfun_g(alpha_sqrt*dr,6)* &
                      !                     & 4D0*atype_lj_parameter(atom_id1,atom_id2,1)*(atype_lj_parameter(atom_id1,atom_id2,2)/dr)**6
                      pme_real = pme_real - gfun_g(alpha_sqrt*dr,6)* atype_lj_parameter(atom_id1,atom_id2,2)/dr**6
                   else
                      stop 'lj_bkghm in disp_usegrid is unknown'
                   endif

                   ! **************************************Feynmann_Hibbs correction ************************
                   Select Case(Feynmann_Hibbs_correction)
                   Case("yes")
                      call generate_reduced_mass(mass_reduced,atom_id1,atom_id2,i_mole,j_mole,n_atom)
                      ! fac = hbar^2/(24*kT * mass_red ) ; since distance is in angstroms, this should be in A^2
                      ! use (hbar*mol)^2 in Kg^2*m^2*A^2/s^2, then k in J/K/mol, mass in Kg
                      fac= hbar_mole_sq / ( 24d0 * 8.31446d0 * temp * mass_reduced / 1000d0 )

                      ! first subtract off all reciprocal space terms
                      Select Case( lj_bkghm_index(atom_id1,atom_id2) )
                      Case(1)   
                         C6_local = atype_lj_parameter(atom_id1,atom_id2,3)
                      Case(2)
                         C6_local = atype_lj_parameter(atom_id1,atom_id2,2)
                      End Select

                      ! C8 effective from C6 term
                      pme_real = pme_real - ( gfun_g(alpha_sqrt*dr,8) - 1d0 ) *  fac * C6_local * ( 42d0 - 12d0 )/dr**8
                      Select Case( C8_10_dispersion_terms)
                      Case("yes")
                         ! C10 effective from C8 term
                         pme_real = pme_real - ( gfun_g(alpha_sqrt*dr,10) - 1d0 ) *  fac * atype_lj_parameter(atom_id1,atom_id2,4) * ( 72d0 - 16d0 )/dr**10            
                         ! C12 effective from C10 term
                         pme_real = pme_real - ( gfun_g(alpha_sqrt*dr,12) - 1d0 ) *  fac * atype_lj_parameter(atom_id1,atom_id2,5) * ( 110d0 - 20d0 )/dr**12
                         Select Case(C12_dispersion)
                         Case("yes")
                            ! C14 effective from C12 term
                            pme_real = pme_real - ( gfun_g(alpha_sqrt*dr,14) - 1d0 ) *  fac * atype_lj_parameter(atom_id1,atom_id2,6) * ( 156d0 - 24d0 )/dr**14
                         End Select
                      End Select

                      ! now add back screened real space contribution
                      call Feynmann_Hibbs_lj_correction(E_local, dr , atom_id1, atom_id2,i_mole,j_mole,n_atom )
                      pme_real = pme_real + E_local
                   End Select
                   ! **************************************************************************************

                enddo
             endif
          endif

       enddo

       disp_usegrid = disp_usegrid + pme_real
    enddo

    return
  end function disp_usegrid

  !******************************************************************************************************************************
  ! This subroutine calculates the total dispersion energy for systems contain n_mole solute molecules and framework. 
  ! It uses the pme dispersion grid for solute-framework interactions, but use simple cutoff for the solute-solute interactions.
  ! Important variables:
  ! real*8, intent(out) :: tot_energy  ! total dispersion energy, note frame-frame interactions are excluded
  !******************************************************************************************************************************
  subroutine tot_disp_usegrid( tot_energy,energy_decomp,tot_n_mole, n_mole, n_atom, xyz, box, cofm)
    use global_variables
    use routines
    use pairwise_interaction
    implicit none
    real*8, intent(out) :: tot_energy
    real*8, dimension(:), intent(out) :: energy_decomp
    integer, intent(in) :: tot_n_mole, n_mole
    integer, dimension(:), intent(in) :: n_atom
    real*8, dimension(:,:) :: box, cofm
    real*8, dimension(:,:,:) :: xyz
    !   local variables
    integer :: i_mole, i_atom, j_mole, j_atom, atom_id1, atom_id2, split_do
    real*8 :: energy_tmp, lj_cutss_use, dr, E_short, E_long, E_local
    real*8, dimension(3) :: shift, dr_com, drij
    real*8, dimension(1,1) :: mass
    real*8 :: tmp, E_sol_frame


    lj_cutss_use = lj_cutoff
    !   No energy decomposition is supported
    energy_decomp = 0d0
    tot_energy = 0d0
    E_sol_frame = 0d0
    !   Stupid dummy variable, nothing useful
    mass = 0.0d0

    if (n_threads .eq. 1 ) then
       split_do = 1
    else 
       split_do = (n_mole/n_threads)/2+1
    endif
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(n_mole,box,cofm,lj_cutss_use,n_atom,xyz,atom_index,atype_lj_parameter,lj_bkghm_index,xyz_to_box_transform) REDUCTION(+:tot_energy)
    !$OMP DO SCHEDULE(DYNAMIC,split_do)
    do i_mole = 1, n_mole
       do j_mole = i_mole+1, n_mole
          shift = pbc_shift( cofm(i_mole,:), cofm(j_mole,:), box, xyz_to_box_transform  )
          dr_com = pbc_dr( cofm(i_mole,:), cofm(j_mole,:), shift )
          if ( abs(dr_com(1))<lj_cutss_use .and. abs(dr_com(2))<lj_cutss_use .and. abs(dr_com(3))<lj_cutss_use ) then
             if ( dot_product( dr_com, dr_com ) < lj_cutss_use**2 ) then
                do i_atom = 1, n_atom(i_mole)
                   do j_atom = 1, n_atom(j_mole)
                      drij = pbc_dr( xyz(i_mole,i_atom,:), xyz(j_mole,j_atom,:), shift )
                      dr = dsqrt(dot_product(drij,drij))
                      atom_id1 = atom_index(i_mole,i_atom)
                      atom_id2 = atom_index(j_mole,j_atom)
                      ! Lennard-Jones
                      if (lj_bkghm_index(atom_id1,atom_id2) .eq. 2) then
                         ! working in C12,C6 now, not epsilon, sigma
                         !  E_short = 4D0 * atype_lj_parameter(atom_id1,atom_id2,1)*(atype_lj_parameter(atom_id1,atom_id2,2)/dr)**12
                         E_short = atype_lj_parameter(atom_id1,atom_id2,1)/dr**12
                         !  E_long = -4D0 * atype_lj_parameter(atom_id1,atom_id2,1)*(atype_lj_parameter(atom_id1,atom_id2,2)/dr)**6
                         E_long = -atype_lj_parameter(atom_id1,atom_id2,2)/dr**6
                         ! Buckingham
                      else
                         E_long = eval_C6_10_dispersion(atom_id1,atom_id2,dr)
                         call eval_short_range_bkghm(E_short, energy_decomp, atom_id1,atom_id2, dr)
                      endif
                      tot_energy = tot_energy + (E_short+E_long)

                      ! if Feynmann_Hibbs correction
                      Select Case(Feynmann_Hibbs_correction)
                      Case("yes")
                         call Feynmann_Hibbs_lj_correction(E_local, dr , atom_id1, atom_id2,i_mole,j_mole,n_atom )
                         tot_energy = tot_energy + E_local
                      End Select

                   enddo
                enddo
             endif
          endif
       enddo
    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL


    !   Looping over all solute molecules
    do i_mole = 1, n_mole
       !     solute-frame contribution, calculated with pme
       E_sol_frame = E_sol_frame + disp_usegrid( box,i_mole,tot_n_mole,n_mole,n_atom,cofm,xyz,mass )
    enddo



    tot_energy = tot_energy + E_sol_frame

    return
  end subroutine tot_disp_usegrid


end module pme_routines
