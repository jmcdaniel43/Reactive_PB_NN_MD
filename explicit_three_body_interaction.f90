!**********************************************************************
! This module is all about evaluating explicit three-body interactions, such
! as Axilrod-Teller 3-body dispersion interactions, or 3-body exchange interactions
!
! Note that such interactions are explicit 3-body interactions, as the formulas for the
! energy contributions explicitly contain 3-bodies.
!
! Contrast this to 3-body induction, which although is a 3-body energy, once the induced
! dipoles are converged, the energy is written as a sum over 2-body terms, and the many-body
! effect is only implicitly realized in the final magnitudes of the induced dipoles
! hence 3-body induction is NOT treated in this module, but in the Drude-oscillator modules
!***********************************************************************


module explicit_three_body_interaction


contains

  !********************************************************************************
  ! This subroutine construct a link list for neighbour searching
  ! For algorithm details, see Allen&Tildesley P149
  ! 
  ! Global Variables changed: 
  ! integer, dimension(:,:,:) :: headlist  ! the head pointers for all cells
  ! integer, dimension(:,:) :: nslist      ! the link list for all the molecules
  !********************************************************************************
  subroutine construct_linklist(box,tot_n_mole,n_mole,r_com)
    use global_variables
    implicit none
    real*8, dimension(:,:), intent(in)::box 
    integer, intent(in) :: tot_n_mole,n_mole
    real*8, dimension(:,:), intent(in) :: r_com
    integer :: i_mole, ix, iy, iz, nx, ny, nz,i_end_mole
    integer, dimension(:,:,:), allocatable :: endlist
    real*8, dimension(3) :: rtmp

    integer , save :: initialize=0

    ! allocate headlist if this is first call
    Select Case(initialize)
    Case(0)
       allocate( headlist(na_nslist,nb_nslist,nc_nslist))
       initialize=1
    End Select

    allocate( endlist(na_nslist,nb_nslist,nc_nslist) )

    nx = na_nslist ; ny = nb_nslist ; nz = nc_nslist ;

    !   pointers all set to null
    headlist = 0
    nslist = 0
    endlist = 0

    ! here n_mole should be the same as tot_n_mole, as the code should have stopped if three body dispersion was requested during a framework simulation

    do i_mole = 1, n_mole
       rtmp = matmul(xyz_to_box_transform(:,:), r_com(i_mole,:))

       ix = floor(rtmp(1)*nx) + 1
       iy = floor(rtmp(2)*ny) + 1
       iz = floor(rtmp(3)*nz) + 1
       !       shift if it is out of the box
       ix = ix - floor(dble(ix-1)/nx)*nx
       iy = iy - floor(dble(iy-1)/ny)*ny
       iz = iz - floor(dble(iz-1)/nz)*nz
       if ( headlist(ix,iy,iz) == 0 ) then
          headlist(ix,iy,iz) = i_mole
          endlist(ix,iy,iz) = i_mole
       else
          i_end_mole = endlist(ix,iy,iz)
          nslist(i_end_mole) = i_mole
          endlist(ix,iy,iz) = i_mole
       endif
    enddo


    return
  end subroutine construct_linklist



  !************************************************
  !   This subroutine loops over all molecule trimers
  !   of a system within a certain cutoff, using a cell-linklist.
  !   this linklist should be updated after every move
  !  
  !
  !   the cutoff is enforced between all pairs in the trimer.  So for a potential trimer ABC, in order
  !   for the trimer to be considered for the interaction, we must have
  !   ( r(AB) < cutoff ) and ( r(AC) < cutoff ) and  ( r(BC) < cutoff )
  !**************************************************
  subroutine calculate_three_body_interaction(E_3B, xyz, r_com, n_atom, box)
    use omp_lib
    use global_variables
    use routines
    real*8, intent(out) :: E_3B
    real*8, dimension(:,:,:), intent(in) :: xyz
    real*8, dimension(:,:), intent(in) :: r_com
    integer, dimension(:), intent(in) :: n_atom
    real*8, dimension(:,:), intent(in) :: box

    real*8 :: ra,rb,rc, dra,drb,drc,dr, rka, rkb, rkc, r_ij, r_ik, r_jk, Eijk
    real*8, save ::three_body_cutoff_molecule
    integer,save :: flag_init=0
    real*8, dimension(3) :: r_com_i, r_com_j, r_com_k, shiftij, shiftik, shiftjk, dr_com
    integer :: i_mole, j_mole,k_mole, l_mole, ia1,ib1,ic1, ia2,ib2,ic2, ia3,ib3,ic3, igrida1,igridb1,igridc1,igrida2,igridb2,igridc2,igrida3,igridb3,igridc3, dia,dib,dic, count1, count2, count3


    Select Case(flag_init)
    Case(0)
       ! get three body cutoff for molecules
       call determine_three_body_cutoff_molecule( three_body_cutoff_molecule )
       flag_init=1
    End Select

    !*****************   real space cutoff for cell link list  *******************
    ! note this is the correct formula for a non-orthogonal box
    ! for orthogonal box, rcut * ra_recip_mag * Ngrid = rcut / ra * Ngrid
    ! where ra is the magnitude of a box vector, and ra_recip_mag is the magnitude of the
    ! corresponding reciprocal box vector
    ! however, for a non-orthogonal box, rcut * ra_recip_mag * Ngrid >= rcut / ra * Ngrid
    ! to cover all possible points within cutoff, we therefore want rcut * ra_recip_mag * Ngrid
    ! i.e. compare lines of code:
    !     rka = dsqrt(dot_product(xyz_to_box_transform(1,:),xyz_to_box_transform(1,:)))
    !     dia = floor(three_body_cutoff_molecule*rka*dble(na_nslist)) + 1
    !versus
    !      ra = dsqrt(dot_product(box(1,:),box(1,:)))
    !      dia = floor(three_body_cutoff_molecule / ra * dble(na_nslist)) + 1
    !
    !
    rka = dsqrt(dot_product(xyz_to_box_transform(1,:),xyz_to_box_transform(1,:)))
    rkb = dsqrt(dot_product(xyz_to_box_transform(2,:),xyz_to_box_transform(2,:)))
    rkc = dsqrt(dot_product(xyz_to_box_transform(3,:),xyz_to_box_transform(3,:)))
    dia = floor(three_body_cutoff_molecule*rka*dble(na_nslist)) + 1
    dib = floor(three_body_cutoff_molecule*rkb*dble(nb_nslist)) + 1
    dic = floor(three_body_cutoff_molecule*rkc*dble(nc_nslist)) + 1
    !
    !******************************************************************************


    E_3B=0d0


    ! loop over all link list grids

    call OMP_SET_NUM_THREADS(n_threads)
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(na_nslist,nb_nslist,nc_nslist,headlist,dia,dib,dic,r_com,box, xyz_to_box_transform,three_body_cutoff_molecule, xyz, n_atom, nslist ) REDUCTION(+:E_3B)
    !$OMP DO SCHEDULE(DYNAMIC,1)
    do ia1 = 1, na_nslist
       igrida1=ia1
       do ib1 = 1, nb_nslist
          igridb1=ib1
          do ic1 = 1, nc_nslist
             igridc1=ic1

             !*************** first molecule *******************

             ! find head molecule for this grid point

             i_mole = headlist(igrida1,igridb1,igridc1)
             count1=0
             ! loop over all molecules in this cell, link list ends with a zero
             do 
                if ( i_mole == 0 ) then
                   exit
                endif
                count1=count1+1
                ! watch for infinite loop
                if ( count1 > MAX_N_MOLE ) then
                   stop " error in cell link list in subroutine calculate_three_body_interaction!"
                endif


                ! loop over cells within cutoff distance from the first cell

                do ia2 = -dia, dia
                   igrida2 = igrida1 + ia2
                   igrida2 = igrida2 - floor(dble(igrida2-1)/dble(na_nslist))*na_nslist
                   do ib2 = -dib, dib
                      igridb2 = igridb1 + ib2
                      igridb2 = igridb2 - floor(dble(igridb2-1)/dble(nb_nslist))*nb_nslist
                      do ic2 = -dic, dic
                         igridc2 = igridc1 + ic2
                         igridc2 = igridc2 - floor(dble(igridc2-1)/dble(nc_nslist))*nc_nslist

                         !********************* second molecule
                         j_mole = headlist(igrida2, igridb2, igridc2)

                         count2=0
                         ! loop over all molecules in this cell, link list ends with a zero
                         do 
                            if ( j_mole == 0 ) then
                               exit
                            endif
                            count2=count2+1
                            ! watch for infinite loop
                            if ( count2 > MAX_N_MOLE ) then
                               stop " error in cell link list in subroutine calculate_three_body_interaction!"
                            endif

                            ! avoid double counting
                            if ( i_mole < j_mole ) then


                               ! loop over cells within cutoff distance from the second cell

                               do ia3 = -dia, dia
                                  igrida3 = igrida2 + ia3
                                  igrida3 = igrida3 - floor(dble(igrida3-1)/dble(na_nslist))*na_nslist
                                  do ib3 = -dib, dib
                                     igridb3 = igridb2 + ib3
                                     igridb3 = igridb3 - floor(dble(igridb3-1)/dble(nb_nslist))*nb_nslist
                                     do ic3 = -dic, dic
                                        igridc3 = igridc2 + ic3
                                        igridc3 = igridc3 - floor(dble(igridc3-1)/dble(nc_nslist))*nc_nslist



                                        !********************* third molecule
                                        k_mole = headlist(igrida3, igridb3, igridc3)

                                        count3=0
                                        ! loop over all molecules in this cell, link list ends with a zero
                                        do 
                                           if ( k_mole == 0 ) then
                                              exit
                                           endif
                                           count3=count3+1
                                           ! watch for infinite loop
                                           if ( count3 > MAX_N_MOLE ) then
                                              stop " error in cell link list in subroutine calculate_three_body_interaction!"
                                           endif

                                           ! avoid double counting
                                           if ( j_mole < k_mole )  then

                                              ! now check whether trimer satisfies cutoffs
                                              r_com_i(:) = r_com(i_mole,:)
                                              r_com_j(:) = r_com(j_mole,:)
                                              r_com_k(:) = r_com(k_mole,:)

                                              ! check i, j pair
                                              shiftij = pbc_shift( r_com_i, r_com_j, box, xyz_to_box_transform )
                                              dr_com = pbc_dr( r_com_i, r_com_j, shiftij )     
                                              if((abs(dr_com(1))<three_body_cutoff_molecule).and.(abs(dr_com(2))<three_body_cutoff_molecule).and.(abs(dr_com(3))<three_body_cutoff_molecule)) then
                                                 r_ij = dsqrt(dot_product(dr_com,dr_com) )
                                                 if ( r_ij < three_body_cutoff_molecule ) then

                                                    ! check i, k pair
                                                    shiftik = pbc_shift( r_com_i, r_com_k, box, xyz_to_box_transform)
                                                    dr_com = pbc_dr( r_com_i, r_com_k, shiftik)     
                                                    if((abs(dr_com(1))<three_body_cutoff_molecule).and.(abs(dr_com(2))<three_body_cutoff_molecule).and.(abs(dr_com(3))<three_body_cutoff_molecule)) then
                                                       r_ij = dsqrt(dot_product(dr_com,dr_com) )
                                                       if ( r_ij < three_body_cutoff_molecule ) then

                                                          ! check j, k pair
                                                          ! At this point we need to be careful checking the third molecule, we don't want to check different periodic images for the three pairs
                                                          ! therefore, rather than computing a shift using the j, k  molecules in the box, use the distance between the j, k images that were used
                                                          ! to compute closest ij, ik pairs, these images are given by shiftij, shiftik

                                                          shiftjk = shiftik - shiftij
                                                          dr_com = pbc_dr( r_com_j, r_com_k, shiftjk)     
                                                          if((abs(dr_com(1))<three_body_cutoff_molecule).and.(abs(dr_com(2))<three_body_cutoff_molecule).and.(abs(dr_com(3))<three_body_cutoff_molecule)) then
                                                             r_ij = dsqrt(dot_product(dr_com,dr_com) )
                                                             if ( r_ij < three_body_cutoff_molecule ) then

                                                                ! ************************** three molecule interaction
                                                                ! call another subroutine that parses these into atom-atom interactions
                                                                call calculate_atom_atom_trimer_interactions(Eijk, xyz, r_com, i_mole, j_mole , k_mole, n_atom, box )
                                                                E_3B = E_3B + Eijk     

                                                             endif
                                                          endif
                                                       endif
                                                    endif                     ! end checks on whether trimer satisfies cutoffs
                                                 endif
                                              endif

                                           endif               ! end if j_mole < k_mole


                                           l_mole = nslist(k_mole)
                                           k_mole = l_mole
                                        enddo ! end loop over link list for molecules in second grid cell

                                     enddo
                                  enddo   ! end loops over 3rd grid cell , for third molecule
                               enddo

                            endif  ! end if i_mole < j_mole


                            k_mole = nslist(j_mole)
                            j_mole = k_mole
                         enddo ! end loop over link list for molecules in second grid cell

                      enddo
                   enddo   ! end loops over 2nd grid cell , for second molecule
                enddo



                j_mole = nslist(i_mole)
                i_mole = j_mole
             enddo             ! end loop over link list for molecules in first grid cell


          enddo
       enddo    ! end loops over 1st grid cell , for first molecule
    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

  end subroutine calculate_three_body_interaction




  !****************************************
  ! this subroutine parses interactions between three molecules, into 
  ! interactions between all intermolecular atom trimers within these
  ! three molecules
  !
  ! the main reason for an extra subroutine here is that we allow
  ! for the possiblity of applying the cutoff either on a molecule-molecule cofm
  ! basis, or a atom-atom basis
  !
  ! NOTE:  This subroutine should only be called for trimers that satisfy (for a cubic box) xij < box_length/2
  !        for every pair ij, and every cartesian coordinate x.  If this isn't satisfied, then there isn't a 
  !        well defined trimer, because the geometry depends with respect to which  molecule we calculate shifts 
  !****************************************
  subroutine calculate_atom_atom_trimer_interactions(Eijk, xyz, r_com, i_mole, j_mole , k_mole, n_atom, box )
    use global_variables
    use routines
    real*8 , intent(out) :: Eijk
    real*8, dimension(:,:,:), intent(in) :: xyz
    real*8, dimension(:,:), intent(in) :: r_com
    integer , intent(in) :: i_mole, j_mole, k_mole
    integer, dimension(:), intent(in) :: n_atom
    real*8, dimension(:,:), intent(in) :: box

    integer :: i_atom, j_atom, k_atom, atom_id_i, atom_id_j, atom_id_k, consider_interaction
    real*8 :: rmag_ij, rmag_jk,rmag_ik,rmag_ij2, rmag_jk2, rmag_ik2, rmag_ij3, rmag_jk3, rmag_ik3, three_body_cutoff_atom, E3atom
    real*8 :: cos_alpha, cos_beta, cos_gamma
    real*8, dimension(3) :: shiftj, shiftk , rij, rjk,rik, r_com_i, r_com_j, r_com_k


   ! set three_body_cutoff atom to three_body_dispersion_cutoff
    ! if we are considering interactions based on molecule-molecule cofm distance only,
    !then we won't use this cutoff
    three_body_cutoff_atom = three_body_dispersion_cutoff


    Eijk=0d0

    ! get any pbc shifts for molecule j, k relative to i
    r_com_i(:) = r_com(i_mole,:)
    r_com_j(:) = r_com(j_mole,:)
    r_com_k(:) = r_com(k_mole,:)

    shiftj = pbc_shift( r_com_i, r_com_j, box, xyz_to_box_transform)
    shiftk = pbc_shift( r_com_i, r_com_k, box, xyz_to_box_transform)

    ! loop over all atoms of all three molecules
    ! note we label vectors with the notation, rij is the vector from atom i to atom j, etc.

    do i_atom=1, n_atom(i_mole)
       do j_atom=1, n_atom(j_mole)

          ! if considering atom-atom cutoffs, then we need atom-atom shifts
          if ( three_body_cutoff_type == 1 ) then
             shiftj = pbc_shift( xyz(i_mole,i_atom,:), xyz(j_mole,j_atom,:), box, xyz_to_box_transform)
          endif

          ! shift j_mole atoms relative to molecule i_mole
          rij(:) = pbc_dr( xyz(i_mole,i_atom,:) , xyz(j_mole,j_atom,:), shiftj(:) )

          rmag_ij2 = dot_product(rij,rij)
          rmag_ij= dsqrt(rmag_ij2)
          rmag_ij3 = rmag_ij * rmag_ij2

          do k_atom=1, n_atom(k_mole)


             ! if considering atom-atom cutoffs, then we need atom-atom shifts
             if ( three_body_cutoff_type == 1 ) then
                shiftk = pbc_shift( xyz(i_mole,i_atom,:), xyz(k_mole,k_atom,:), box, xyz_to_box_transform)   
             endif


             ! shift k_mole atoms relative to molecule i_mole
             rik(:) = pbc_dr( xyz(i_mole,i_atom,:) , xyz(k_mole,k_atom,:), shiftk(:) )

             rjk = rik - rij



             rmag_ik2 = dot_product(rik,rik)
             rmag_ik= dsqrt(rmag_ik2)
             rmag_ik3 = rmag_ik * rmag_ik2
             rmag_jk2 = dot_product(rjk,rjk)
             rmag_jk= dsqrt(rmag_jk2)
             rmag_jk3 = rmag_jk * rmag_jk2


             ! decide whether to compute this interaction
             Select Case( three_body_cutoff_type )
             Case(0)
                ! cutoff considered on a molecule-molecule basis, therefore compute every atom-atom interaction
                consider_interaction=1
             Case(1)
                !******* apply atom-atom cutoff
                if ( ( rmag_ik < three_body_cutoff_atom ) .and. ( rmag_jk < three_body_cutoff_atom ) .and. ( rmag_ij < three_body_cutoff_atom ) ) then
                   consider_interaction=1
                else
                   consider_interaction=0
                endif
             End Select


             if ( consider_interaction == 1 ) then

                ! we need these cosines for axilrod-teller

                ! cosines of angles of trimer triangle
                cos_alpha = dot_product(rij,rik)/ (rmag_ij * rmag_ik)
                cos_beta = dot_product(rik,rjk) / (rmag_jk * rmag_ik)
                ! note the sign here
                cos_gamma = dot_product(rij,rjk) / (rmag_ij * rmag_jk)
                cos_gamma = -cos_gamma


                ! get the f9 damping function
                atom_id_i = atom_index(i_mole,i_atom)
                atom_id_j = atom_index(j_mole,j_atom)
                atom_id_k = atom_index(k_mole,k_atom)

                ! always three body dispersion if this subroutine is called
                call calculate_axilrod_teller_3atom(E3atom, atom_id_i, atom_id_j, atom_id_k, rmag_ij, rmag_ik, rmag_jk, rmag_ij3, rmag_ik3, rmag_jk3, cos_alpha, cos_beta, cos_gamma )
                Eijk = Eijk + E3atom
                ! if three body exchange
                Select Case(three_body_exchange)
                Case("yes")
                   call calculate_three_body_exchange_3atom(E3atom,atom_id_i, atom_id_j, atom_id_k, rmag_ij, rmag_ik, rmag_jk  )
                   Eijk = Eijk + E3atom  
                End Select

             end if

          enddo
       enddo
    enddo


  end subroutine calculate_atom_atom_trimer_interactions





  !********************************************
  ! 
  ! This subroutine calculates the three body dispersion energy of 
  ! a three molecule trimer, assuming an Axilrod-Teller form for the three body
  ! interaction.  
  !
  ! the Axilrod-Teller form is given by
  !
  ! f(9)*C(9,abc) * (1 + 3cos(alpha)*cos(beta)*cos(gamma)) / ( R(ab)^3 * R(bc)^3 * R(ac) ^3 )
  !
  ! where C(9,abc) is the dispersion coefficient for atom triplet abc,
  ! alpha, beta, gamma are the angles of the triangle formed by the three atoms,
  ! and f(9) is a damping function, which is a product of three two-body Tang-Toennies functions
  !
  ! See ::  Wen and Beran, JCTC, 2011, 7, 3733-3742
  !         Axilrod and Teller, JCP, 1943, 11, 299-300
  !         Lotrich and Szalewicz, JCP, 1997, 106, 9688-9702
  !
  !********************************************
  subroutine calculate_axilrod_teller_3atom(E3atom, atom_id_i, atom_id_j, atom_id_k, rmag_ij, rmag_ik, rmag_jk, rmag_ij3, rmag_ik3, rmag_jk3, cos_alpha, cos_beta, cos_gamma )
    use global_variables
    use routines
    real*8 , intent(out) :: E3atom
    real*8, intent(in) :: rmag_ij, rmag_jk,rmag_ik,rmag_ij3, rmag_jk3, rmag_ik3,cos_alpha, cos_beta, cos_gamma
    integer, intent(in) :: atom_id_i, atom_id_j, atom_id_k

    real*8 ::  f9

    ! get the f9 damping function
    f9 = C3_damp(atom_id_i,atom_id_j,rmag_ij) * C3_damp(atom_id_i,atom_id_k,rmag_ik) * C3_damp(atom_id_j,atom_id_k,rmag_jk)

    E3atom = f9 * atype_3body_C9(atom_id_i, atom_id_j, atom_id_k ) * ( 1d0 + 3d0 * cos_alpha * cos_beta * cos_gamma ) / ( rmag_ij3 * rmag_ik3 * rmag_jk3 )

  end subroutine calculate_axilrod_teller_3atom



  !********************************************
  ! 
  ! This subroutine calculates the three body exchange energy of 
  ! a three molecule trimer, assuming an exponential form
  ! E = A*exp(-B(Rab + Rbc + Rca ) )
  !
  !
  !********************************************
  subroutine calculate_three_body_exchange_3atom(E3atom,atom_id_i, atom_id_j, atom_id_k, rmag_ij, rmag_ik, rmag_jk  )
    use global_variables
    real*8 , intent(out) :: E3atom
    real*8, intent(in) :: rmag_ij, rmag_jk,rmag_ik
    integer, intent(in) :: atom_id_i, atom_id_j, atom_id_k

     E3atom = atype_3body_exchange(atom_id_i, atom_id_j, atom_id_k,1) * exp(-atype_3body_exchange(atom_id_i, atom_id_j, atom_id_k,2) * ( rmag_ij + rmag_ik + rmag_jk ) )


  end subroutine calculate_three_body_exchange_3atom






  !**********************************************
  ! this subroutine outputs the three_body_cutoff that should be used
  ! on a molecule by molecule basis
  !
  ! if the variable three_body_cutoff_type == 0 , then all three body interactions
  ! are considered on a molecule by molecule basis (i.e. all atoms within a molecule are grouped)
  ! and then three_body_cutoff_molecule is just the cutoff that was input
  !
  ! however if three_body_cutoff_type == 1 , then interactions should be considered
  ! on an atom by atom basis, and therefore when checking molecule-molecule distances
  ! the cutoff should be bigger to allow for exterior atoms to potentially interact
  !
  ! currently, the base cutoff that is considered is variable "three_body_dispersion_cutoff"
  !***********************************************
  subroutine determine_three_body_cutoff_molecule( three_body_cutoff_molecule )
    use global_variables
    real*8, intent(out) :: three_body_cutoff_molecule

    real*8 :: max_dr, dr_com(3), dr_test
    integer :: i_type,i_atom, index_atom, index_molecule


    Select Case( three_body_cutoff_type )
    Case(0)
       ! cutoff is molecule by molecule basis, so use the input cutoff
       three_body_cutoff_molecule = three_body_dispersion_cutoff
    Case(1)
       ! here the final cutoff is on an atom by atom basis, using the input cutoff
       ! therefore molecule_cutoff = atom_cutoff + 2 * max(dr_atom_cofm)
       ! where max(dr_atom_cofm) is the maximum distance of any atom from the
       ! center of mass of its molecule

       max_dr = 0d0

       ! loop over all molecule types, and figure out what max_dr is

       do i_type=1, n_molecule_type

          ! loop over all atoms in this molecule type, end is marked with MAX_N_ATOM + 1
          do i_atom =1, MAX_N_ATOM
             if( molecule_type(i_type,i_atom) .eq. MAX_N_ATOM + 1) exit

             dr_com(:) = r_mole_com(i_type,i_atom,:)
             dr_test = dsqrt(dot_product(dr_com, dr_com ) )

             if ( max_dr < dr_test ) then
                max_dr = dr_test
                index_molecule = i_type
                index_atom = i_atom
             endif

          enddo

       enddo

       three_body_cutoff_molecule = three_body_dispersion_cutoff + 2d0 * max_dr

       write(*,*) "three_body_cutoff_molecule is set to ", three_body_cutoff_molecule
       write(*,*) "this is greater than the three_body_dispersion_cutoff for atoms ", three_body_dispersion_cutoff, " ,"
       write(*,*) "so as to encompass the largest atom-cofm distance for atomtype ", atype_name(molecule_type(index_molecule,index_atom)), "of molecule ", molecule_type_name(index_molecule)
       write(*,*) ""

    case default
       write(*,*) " couldn't recognize setting of variable 'three_body_cutoff_type' in subroutine"
       write(*,*) "'determine_three_body_cutoff_molecule'.  This variable should be set to either 0 or 1"
       stop 

    End Select

  end subroutine determine_three_body_cutoff_molecule


end module explicit_three_body_interaction
