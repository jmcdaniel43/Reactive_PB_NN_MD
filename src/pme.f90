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





  !***********************************************************************
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
  !***********************************************************************
  subroutine  pme_energy_force(pme_energy, pme_force,target_atoms,tot_n_mole,n_mole, n_atom, xyz, chg, box, r_com, dfti_desc,dfti_desc_inv)
    use global_variables
    use MKL_DFTI
    use routines
    use omp_lib
    implicit none
    real*8, intent(out) :: pme_energy
    integer, intent(in) :: n_mole,tot_n_mole
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: box, r_com
    real*8, intent(in), dimension(:,:) :: chg 
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

    Case("yes")
       ! here we need to map molecular data structures to atomic data structures that are
       ! consistent with the verlet list
       ! as these arrays are in inner loops, put innermost index first
       allocate( atom_list_xyz(3,verlet_elec_atoms), atom_list_force(3,verlet_elec_atoms) , atom_list_chg(verlet_elec_atoms) )
       call map_molecule_atom_data_structures_3d_m_to_a( tot_n_mole, n_atom , verlet_molecule_map_elec, atom_list_xyz , xyz )
       call map_molecule_atom_data_structures_2d_m_to_a( tot_n_mole, n_atom , verlet_molecule_map_elec, atom_list_chg , chg ) 
       !*******   pme real space energy and force subroutine
       call pme_real_space_use_verlet(Er, atom_list_force,atom_list_xyz, atom_list_chg, box)
       ! convert forces back to molecular data structures
       call map_molecule_atom_data_structures_3d_a_to_m( tot_n_mole, n_atom, verlet_molecule_map_elec , atom_list_force , pme_force )
       deallocate( atom_list_xyz , atom_list_force, atom_list_chg )

    !******************************** intra-molecular real-space interactions *******************************
    do i_mole=1,n_mole
       if(n_atom(i_mole) .gt. 1) then
          do i_atom=1,n_atom(i_mole)-1
             do j_atom=i_atom+1,n_atom(i_mole)
                call intra_pme_energy(E_intra,xyz,chg,i_mole,i_atom,j_atom,n_atom, molecule_index)
                Er = Er + E_intra
                call intra_pme_force(f_ij,xyz,chg,i_mole,n_atom,i_atom,j_atom, molecule_index)
                pme_force(i_mole,i_atom,:)=pme_force(i_mole,i_atom,:) + f_ij(:)
                pme_force(i_mole,j_atom,:)=pme_force(i_mole,j_atom,:) - f_ij(:)                         
             enddo
          enddo
       endif
    enddo


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
    !****** note the storage dimensions of these arrays, consistent with atom_list_force
    real*8, dimension(3,size(atom_list_force(1,:))) :: local_force
    real*8, dimension(3,size(atom_list_force(1,:)),n_threads) :: temp_force

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

             ! coordinates should be stored in this order, optimizing fortran
             ! column major memory retrieval
             rij = atom_list_xyz(:,i_atom) - atom_list_xyz(:,j_atom)

             ! shift for general box
             !             dr_direct(:) = matmul( xyz_to_box_transform, rij )
             !             do i=1,3
             !                shift_direct(i) = dble(floor( dr_direct(i) + 0.5d0 ))
             !             enddo
             !             shift = matmul( shift_direct , box )
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

                   ! forces should be stored in this order, optimizing fortran
                   ! column major memory retrieval

                   local_force(:,i_atom) = local_force(:,i_atom) + f_ij(:)
                   local_force(:,j_atom) = local_force(:,j_atom) - f_ij(:)
                end if
             end if
          end do
       end if
    end do

    !$OMP END DO NOWAIT
    thread_id = OMP_GET_THREAD_NUM()
    temp_force(:,:,thread_id+1)= local_force(:,:)
    !$OMP END PARALLEL


    do i_thread=1,n_threads
       do i_atom=1, verlet_elec_atoms
          atom_list_force(:,i_atom) = atom_list_force(:,i_atom) + temp_force(:,i_atom,i_thread)
       enddo
    enddo

    pme_real_space_energy = Er

  end subroutine pme_real_space_use_verlet



  !****************************************************
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
  !*****************************************************

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
    real*8 :: kk(3,3)
    real*8 :: scale
    integer::i,j,l,status,n,K  ! n is spline order, K is grid size
    !*************even though arrays begin at index 1 for convenience, 1,1,1 corresponds to 0,0,0
    complex*16,dimension(:,:,:),allocatable::FQ
    real*8,dimension(:), allocatable::q_1r
    complex*16,dimension(:), allocatable::q_1d
    real*8, dimension(:,:), allocatable :: atom_list_xyz
    real*8, dimension(:), allocatable :: atom_list_chg
    integer,dimension(MAX_N_MOLE,MAX_N_ATOM) :: atom_molecule_map
    integer :: total_atoms_list


    !****************** if need to update Q_array
    if(Q_change.ne.0) then

       allocate(xyz_scale(n_mole,maxval(n_atom(:)),3),FQ(pme_grid,pme_grid,pme_grid), q_1r(pme_grid**3), q_1d(pme_grid**3))
       n=spline_order
       K=pme_grid
       call construct_reciprocal_lattice_vector(kk, box)
       ! create scaled coordinates
       call create_scaled_direct_coordinates(xyz_scale, xyz, n_mole, n_atom, kk, K)


       ! note that if we are using a verlet list, we can use the verlet list as the atom-molecule map
       ! here, we generate a new map, in case we are not using a verlet list
       call generate_verlet_atom_index( total_atoms_list, atom_molecule_map, n_mole, n_atom, chg )

       ! Q grid will be constructed using atomic data storage for efficiency, convert data structures to this format
       allocate( atom_list_xyz(3,total_atoms_list) , atom_list_chg(total_atoms_list) )
       call map_molecule_atom_data_structures_3d_m_to_a( n_mole, n_atom , atom_molecule_map, atom_list_xyz , xyz_scale )
       call map_molecule_atom_data_structures_2d_m_to_a( n_mole, n_atom , atom_molecule_map, atom_list_chg , chg ) 
       ! grid_Q
       call grid_Q(Q_grid,atom_list_chg,atom_list_xyz,K,n)
       deallocate( atom_list_xyz, atom_list_chg )

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
       !  scale = 1.0D0/dble(K)**3


       status = DftiComputeBackward(dfti_desc_inv, q_1d)
       FQ=RESHAPE(q_1d, (/K,K,K/) )
       theta_conv_Q=dble(FQ)

       deallocate( FQ, q_1r, q_1d, xyz_scale )
    endif

    pme_recip=.5D0*sum((Q_grid*theta_conv_Q))


    ! store the reciprocal space energy for ms-evb if needed
    Select Case(ms_evb_simulation)
    Case("yes")
       E_recip=pme_recip
    End Select

  end function pme_recip




  !************************************************
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
  !
  !*************************************************

  subroutine  pme_force_recip(pme_force,target_atoms,tot_n_mole,n_mole, n_atom, xyz, chg, box, dfti_desc,dfti_desc_inv)
    use global_variables
    use MKL_DFTI
    use omp_lib
    implicit none
    integer, intent(in) :: n_mole,tot_n_mole
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: box
    real*8, intent(in), dimension(:,:) :: chg    
    integer,intent(in), dimension(:,:) :: target_atoms
    real*8,intent(inout),dimension(:,:,:):: pme_force
    real*8, intent(in), dimension(:,:,:) :: xyz 
    TYPE(DFTI_DESCRIPTOR),pointer,intent(in):: dfti_desc,dfti_desc_inv
    real*8, dimension(:,:,:),allocatable::xyz_scale
    real*8 :: a(3), b(3), c(3), ka(3), kb(3), kc(3),kk(3,3)
    real*8 :: vol,scale,force(3)
    integer::i,j,l,i_atom,i_mole,status,n,K  ! n is spline order, K is grid size
!!!!!!!!!!!!!!!!!!!!!!!even though arrays begin at index 1 for convenience, 1,1,1 corresponds to 0,0,0
    complex*16,dimension(:,:,:),allocatable::FQ
    real*8,dimension(:), allocatable::q_1r
    complex*16,dimension(:), allocatable::q_1d
    real*8, dimension(:,:), allocatable :: atom_list_xyz, atom_list_force
    real*8, dimension(:), allocatable :: atom_list_chg
    integer,dimension(MAX_N_MOLE,MAX_N_ATOM) :: atom_molecule_map
    integer :: total_atoms_list
    integer:: split_do, index_store
    real*8 :: small=1D-6
    !test
    integer :: nq(3), save_nt
    real*8,dimension(3) :: force_temp
    real*8,dimension(3,MAX_N_ATOM,MAX_N_MOLE) :: test_force, xyz_scale1
    real*8,dimension(MAX_N_ATOM,MAX_N_MOLE) :: test_chg
    real*8,dimension(:,:),allocatable :: test_force1



    n=spline_order
    K=pme_grid

    call construct_reciprocal_lattice_vector(kk, box)
    ! create scaled coordinates
    allocate(xyz_scale(tot_n_mole,maxval(n_atom(:)),3))
    call create_scaled_direct_coordinates(xyz_scale, xyz, tot_n_mole, n_atom, kk, K)


    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "Q_grid started at", time
    endif
    !************************************************************************!


    ! note that if we are using a verlet list, we can use the verlet list as the atom-molecule map
    ! here, we generate a new map, in case we are not using a verlet list
    call generate_verlet_atom_index( total_atoms_list, atom_molecule_map, n_mole, n_atom, chg )

    ! Q grid will be constructed using atomic data storage for efficiency, convert data structures to this format
    allocate( atom_list_xyz(3,total_atoms_list) , atom_list_chg(total_atoms_list) )
    call map_molecule_atom_data_structures_3d_m_to_a( n_mole, n_atom , atom_molecule_map, atom_list_xyz , xyz_scale )
    call map_molecule_atom_data_structures_2d_m_to_a( n_mole, n_atom , atom_molecule_map, atom_list_chg , chg ) 
    ! grid_Q
    call grid_Q(Q_grid,atom_list_chg,atom_list_xyz,K,n)
    deallocate( atom_list_xyz, atom_list_chg )


    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "Q_grid finished at", time
    endif
    !************************************************************************!

    allocate( FQ(pme_grid,pme_grid,pme_grid), q_1r(pme_grid**3), q_1d(pme_grid**3) )

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
    !  scale = 1.0D0/dble(K)**3


    status = DftiComputeBackward(dfti_desc_inv, q_1d)

    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "Fourier transforms finished at", time
    endif
    !************************************************************************!

    FQ=RESHAPE(q_1d, (/K,K,K/) )
    theta_conv_Q=dble(FQ)

    deallocate( FQ, q_1r, q_1d )

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


    ! notice these atom-based data structures are stored in column major order for looping over atom indices
    allocate( atom_list_xyz(3,total_atoms_list), atom_list_force(3,total_atoms_list) , atom_list_chg(total_atoms_list) )
    call map_molecule_atom_data_structures_3d_m_to_a( tot_n_mole, n_atom , atom_molecule_map, atom_list_xyz , xyz_scale )
    call map_molecule_atom_data_structures_2d_m_to_a( tot_n_mole, n_atom , atom_molecule_map, atom_list_chg , chg )
    atom_list_force=0d0

    call OMP_SET_NUM_THREADS(n_threads)
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(n_threads,theta_conv_Q,atom_list_xyz,atom_list_chg,K,box,n,kk,atom_list_force,spline_order,split_do) 
    !$OMP DO SCHEDULE(dynamic, split_do)
    do i_atom=1,size(atom_list_force(1,:))
       call derivative_grid_Q(force,theta_conv_Q,atom_list_chg,atom_list_xyz,i_atom,K,box,n,kk)
       atom_list_force(:,i_atom)=atom_list_force(:,i_atom)+force(:)
    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "Q grid derivatives finished at", time
    endif
    !************************************************************************!

    ! store the reciprocal space force for ms-evb if needed
    Select Case(ms_evb_simulation)
    Case("yes")
       atom_list_force_recip=atom_list_force
    End Select

    ! map force back to molecule data structures
    call map_molecule_atom_data_structures_3d_a_to_m( tot_n_mole, n_atom, atom_molecule_map, atom_list_force, pme_force, "add" )

    deallocate(xyz_scale, atom_list_xyz, atom_list_force, atom_list_chg)

  end subroutine pme_force_recip





  !************************************************
  ! this subroutine removes intra-molecular energy contributions
  !
  ! use explicit erfc since r can be very small
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
  !************************************************

  subroutine intra_pme_energy(E_intra,xyz,chg,i_mole,i_atom,j_atom,n_atom, molecule_index_local)
    use global_variables
    implicit none
    real*8, intent(out):: E_intra
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: chg
    real*8, intent(in), dimension(:,:,:) :: xyz
    integer,intent(in)::i_mole,i_atom,j_atom
    integer, dimension(:), intent(in) :: molecule_index_local

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
    i_mole_type = molecule_index_local(i_mole)
    ! check for exclusions 
    if ( molecule_exclusions( i_mole_type, i_atom, j_atom ) /= 1 ) then 
       E_intra = E_intra + chg(i_mole,i_atom) * chg(i_mole,j_atom) / norm_dr
    end if


  end subroutine intra_pme_energy




  !******************************************
  ! this subroutine removes intra-molecular interactions
  !
  ! use explicit erfc since r can be very small
  !******************************************
  subroutine intra_pme_force(f_ij,xyz,chg,i_mole,n_atom,i_atom,j_atom, molecule_index_local)
    use global_variables
    implicit none
    real*8, dimension(:), intent(out):: f_ij
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: chg
    real*8, intent(in), dimension(:,:,:) :: xyz
    integer,intent(in)::i_mole,i_atom,j_atom
    integer, dimension(:), intent(in) :: molecule_index_local

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
    i_mole_type = molecule_index_local(i_mole)
    ! check for exclusions 
    if ( molecule_exclusions( i_mole_type, i_atom, j_atom ) /= 1 ) then 
       f_ij = f_ij + chg(i_mole,i_atom) * chg(i_mole,j_atom) * rij / norm_dr**3
    end if


  end subroutine intra_pme_force



  !*****************************************************************
  ! This subroutine interpolates charges onto Q grid to be used in pme reciprocal space
  ! routines
  !
  !  4/15/15 : Note, we have changed this subroutine so that it uses atom-based storage
  ! data arrays, rather than molecule-based storage data arrays.  Also, the atom-based data
  ! arrays have been stored so that outer-loop runs over the last array index, to take
  ! advantage of column-major memory storage in fortran
  ! Also, at this point it is assumed that only atoms with non-zero charge are included in
  ! data structures and thus we've removed the check on small charge
  ! we have also changed loop over kvectors to loop over last index of Q_grid first,
  ! to more efficiently pull from memory
  !*****************************************************************
  subroutine grid_Q(Q,chg,xyz,K,n)
    use global_variables
    use omp_lib
    real*8, intent(in), dimension(:) :: chg
    real*8, intent(in), dimension(:,:) :: xyz
    integer,intent(in)::K,n
    real*8,dimension(:,:,:),intent(out)::Q
    integer::i_atom,tot_atoms, k1,k2,k3,n1,n2,n3,nn1,nn2,nn3,nearpt(3),splindex(3)
    real*8::sum
    real*8,dimension(3)::u,arg
    integer :: split_do

    Q=0D0
    tot_atoms=size(xyz(1,:))

    ! decide how to split the parallel section
    if (n_threads .eq. 1 ) then
       split_do = 1
    else
       ! integer arithmetic
       split_do = tot_atoms/n_threads
    endif

    ! parameter spline_grid undeclared, but ok
    call OMP_SET_NUM_THREADS(n_threads)
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(split_do,xyz,chg,tot_atoms,n,B6_spline,B4_spline,K,spline_order) REDUCTION(+:Q)
    !$OMP DO SCHEDULE(dynamic, split_do)
    do i_atom=1,tot_atoms
       u=xyz(:,i_atom)
       nearpt=floor(u)

       ! loop over outer index of Q grid first, to more efficiently use memory
       ! only need to go to k=0,n-1, for k=n, arg > n, so don't consider this
       do k3=0,n-1
          n3=nearpt(3)-k3
          arg(3)=u(3)-dble(n3);
          ! shift index of array storage if < 0
          if(n3<0) then
             n3=n3+K
          endif
          do k2=0,n-1
             n2=nearpt(2)-k2
             arg(2)=u(2)-dble(n2)
             ! shift index of array storage if < 0
             if(n2<0) then
                n2=n2+K
             endif
             do k1=0,n-1
                n1=nearpt(1)-k1
                arg(1)=u(1)-dble(n1);
                ! shift index of array storage if < 0
                if(n1<0) then
                   n1=n1+K
                endif

                sum=0d0
                splindex = ceiling(arg/6.D0*dble(spline_grid))

                ! note 0<arg<n , so arg should always be within bounds of gridded spline
                if(spline_order .eq.6) then
                   sum=chg(i_atom)*B6_spline(splindex(1))*B6_spline(splindex(2))*B6_spline(splindex(3))
                else
                   sum=chg(i_atom)*B4_spline(splindex(1))*B4_spline(splindex(2))*B4_spline(splindex(3))   
                endif

                Q(n1+1,n2+1,n3+1)=Q(n1+1,n2+1,n3+1)+sum
             enddo
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

  end subroutine grid_Q





  !*********************************
  ! this subroutine either adds or subtracts the contribution of 
  ! i_mole to the Q_grid based on whether input parameter
  ! "operation" is positive or negative
  !
  ! see subroutine grid_Q for comments on algorithm, here we have deleted
  ! the redundant comments
  !
  ! JGM 4/15/15
  ! note the data structures used in this subroutine (xyz, chg) are in
  ! molecular storage format (i_mole, i_atom) in contrast to those in grid_Q subroutine
  !
  ! we have also changed loop over kvectors to loop over last index of Q_grid first,
  ! to more efficiently pull from memory
  !*********************************
  subroutine modify_Q_grid( Q , i_mole, chg, xyz, n_atom, K, n, operation)
    use global_variables
    integer, intent(in) :: i_mole
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: chg
    real*8, intent(in), dimension(:,:,:) :: xyz
    integer,intent(in)::K,n
    real*8,dimension(:,:,:),intent(inout)::Q
    integer, intent(in) :: operation
    integer::i,j,k1,k2,k3,n1,n2,n3,nn1,nn2,nn3,nearpt(3),splindex(3)
    real*8::sum,chg_i
    real*8,dimension(3)::u,arg
    real*8 :: small=1D-6
    i= i_mole
    do j=1,n_atom(i)
       chg_i=chg(i,j)
       if ( abs(chg_i) > small ) then
          u=xyz(i,j,:)
          nearpt=floor(u)
          ! loop over outer index of Q grid first, to more efficiently use memory
          do k3=0,n-1
             n3=nearpt(3)-k3
             arg(3)=u(3)-dble(n3);
             if(n3<0) then
                n3=n3+K
             endif
             do k2=0,n-1
                n2=nearpt(2)-k2
                arg(2)=u(2)-dble(n2)
                if(n2<0) then
                   n2=n2+K
                endif
                do k1=0,n-1
                   n1=nearpt(1)-k1
                   arg(1)=u(1)-dble(n1);
                   if(n1<0) then
                      n1=n1+K
                   endif
                   sum=0d0
                   splindex = ceiling(arg/6.D0*dble(spline_grid))
                   if(spline_order .eq.6) then
                      sum=chg_i*B6_spline(splindex(1))*B6_spline(splindex(2))*B6_spline(splindex(3))
                   else
                      sum=chg_i*B4_spline(splindex(1))*B4_spline(splindex(2))*B4_spline(splindex(3))      
                   endif
                   if ( operation < 0 ) then
                      Q(n1+1,n2+1,n3+1)=Q(n1+1,n2+1,n3+1)-sum
                   else
                      Q(n1+1,n2+1,n3+1)=Q(n1+1,n2+1,n3+1)+sum
                   end if
                enddo
             enddo
          enddo
       end if
    enddo
  end subroutine modify_Q_grid



  !**************************************************
  ! this subroutine computes derivatives of grid Q with respect 
  ! to changes in position of i_atom in i_mole
  ! output is the force, which is computed by taking the product of
  ! the derivative Q array with FQ, which is the convoluted array
  !
  !
  ! notice that this is being called in a parallel section of code
  !**************************************************
  subroutine derivative_grid_Q(force,FQ,chg,xyz,i_atom,K,box,n,kk,flag_update)
    use global_variables
    integer,intent(in)::K,n
    real*8,dimension(3),intent(out)::force
    real*8,dimension(:,:,:),intent(in)::FQ
    integer, intent(in) :: i_atom
    real*8, intent(in), dimension(:) :: chg
    real*8,intent(in),dimension(3,3)::box,kk
    real*8, intent(in), dimension(:,:) :: xyz
    integer, intent(in), optional :: flag_update
    integer::i,j,k1,k2,k3,n1,n2,n3,nearpt(3)
    integer::g1n(3),g1nmin(3),g1(3),g2(3)
    real*8::fac(3),chg_i
    real*8,dimension(3)::u,arg1,arg2,force_temp

    integer :: count
    count=1

    force=0.D0
    chg_i=chg(i_atom)
    ! coordinates should be stored this way for optimal memory retrieval
    u=xyz(:,i_atom)
    nearpt=floor(u)

    ! loop over K grid according to memory storage n3, n2, n1
    ! first part of fac
    do k3=0,n-1
       n3=nearpt(3)-k3
       arg1(3)=u(3)-dble(n3);
       arg2(3) = arg1(3) - 1d0
       ! shift index of array storage if < 0
       if(n3<0) then
          n3=n3+K
       endif
!!!!!!!!get bspline values from grid, therefore must check whether to see if function is being evaluated in non-zero domain, otherwise, can't call array
       ! note arg1 > 0, so don't check for that.  for k1=n, arg1 > n
       !if((arg1(1)<real(n)).and.(0d0<arg1(1))) then
       do k2=0,n-1
          n2=nearpt(2)-k2
          arg1(2)=u(2)-dble(n2)
          arg2(2) = arg1(2) - 1d0
          ! shift index of array storage if < 0
          if(n2<0) then
             n2=n2+K
          endif
          ! if k1=n, arg1 >= n, so we can exclude that term in the sum, also because we exlude that term in the sum, arg1 <= n , so we dont need to check the limits of the bspline evaluation
          do k1=0,n-1
             n1=nearpt(1)-k1
             arg1(1)=u(1)-dble(n1);
             arg2(1) = arg1(1) - 1d0
             ! shift index of array storage if < 0
             if(n1<0) then
                n1=n1+K
             end if
             fac=0d0
             ! use ceiling here, instead of int (which is faster), if we are not
             ! going to check that g1 > 0
             g2=ceiling(arg2/real(n-1)*real(spline_grid))
             g1n=ceiling(arg1/real(n)*real(spline_grid))
             g1nmin = ceiling(arg1/real(n-1)*real(spline_grid))


!!!!!force in x direction
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

             force=force + fac * FQ(n1+1,n2+1,n3+1)

	     ! if we're storing dQ_dr,
             if ( .not. present(flag_update) ) then
                Select Case(ms_evb_simulation)
                Case("yes")
                   !if ( index_store > size(dQ_dr(:,1)) ) then
                   !stop "size dQ_dr not allocated correctly"
                   !endif
                   dQ_dr(:,count,i_atom) = fac(:)
                   dQ_dr_index(1,count,i_atom) =n1+1
                   dQ_dr_index(2,count,i_atom) =n2+1
                   dQ_dr_index(3,count,i_atom) =n3+1
                   count = count + 1
                End Select
             endif
          enddo
       enddo
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



  !***************************************************************************
  ! This function update Ewald self interaction based on current configuration
  ! It updates:
  ! real*8 :: Ewald_self
  !
  ! Note that Ewald self corrects for ONE HALF of the interaction between a
  ! charge
  ! interacting with it's own Gaussian.  This is because in the Ewald reciprocal
  ! sum,
  ! There is a factor of one half to compensate for overcounting, and this
  ! factor
  ! therefore carries over for the self interaction
  !
  ! the total interaction of a charge with its Gaussian is 2 * (q1)^2 * ( alpha
  ! / pi ) ^ (1/2)
  ! so that's why the self correction is (q1)^2 * ( alpha / pi ) ^ (1/2)
  ! **************************************************************************
  subroutine update_Ewald_self( n_mole, n_atom, chg )
    use global_variables
    implicit none
    integer, intent(in) :: n_mole
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: chg
    integer :: i_mole, i_atom
    Ewald_self = 0.0
    do i_mole = 1, n_mole ! loop over all atoms
       do i_atom = 1, n_atom(i_mole)
          Ewald_self = Ewald_self - chg(i_mole,i_atom)**2
       end do
    end do
    Ewald_self = Ewald_self * alpha_sqrt/sqrt(pi)
  end subroutine update_Ewald_self


end module pme_routines
