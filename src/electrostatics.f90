module electrostatic
  use routines
  implicit none

  !*****************************************************************
  !  Despite the name, not all electrostatic subroutines are contained in this module.
  !  Because of the size of the subroutines, ===>
  !
  ! *******ALL PARTICLE MESH EWALD (PME) SUBROUTINES ARE LOCATED IN PME.F90***************
  !
  !  NOTE:  The ewald routines in this module are tested and work, but are never
  !  used in the main code as they are slow and pme is much faster
  !  they are kept here just for testing only
  !
  !  The energy routines here that are currently used are the electrostatic cutoff routines
  !  for the Particle Mesh Ewald (PME) routines, see pme.f90
  !
  !  All energy should be output in kJ/mol, length in A
  !*****************************************************************

contains

  !*********************************************************************
  ! This function calculates ewald electrostatic energy
  ! 
  ! NOTE:: This ewald implementation should probably never be used in simulation, since PME is
  ! much faster. This function was mainly written to test the PME subroutines.
  !
  !  Local variables:
  !    n_mole : number of molecule in the system
  !    n_atom : a list contains number of atoms in each molecule
  !       chg : charge lists
  !       xyz : cartesian coordinates list
  !       box : a 3*3 matrix, including a, b, c of the box
  !      mass : a 2-D list of atomic mass
  !
  ! Global variables this function updates and uses:
  !  integer :: totk
  !  character(3) :: ewald_volume_move    ! this tells code which data arrays to update
  !  complex*8, dimension(maxk) :: rou_k  ! Charge density in rec. space
  !  real*8, dimension(maxk,3) :: k_vec   ! All k vectors, no need to change if box size is not changed
  !  real*8, dimension(maxk) :: B_factor  ! B-factors associated with k_vec
  !*********************************************************************
  real*8 function ewald( n_mole, n_atom, xyz, chg, box, r_com )
    use global_variables
    implicit none
    integer, intent(in) :: n_mole
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: chg, box, r_com
    real*8, intent(in), dimension(:,:,:) :: xyz
    integer :: i, j, k, N, i_mole, i_atom, j_mole, j_atom
    real*8, dimension(3) :: ri, rj, rij, k1, k2, k3, ktot, shift, r_com_i, r_com_j, dr_com
    real*8 :: Ek, Er, qi, qj, norm_dr
    real*8 :: a(3), b(3), c(3), ka(3), kb(3), kc(3), k_now(3)
    real*8 :: vol, ksq, factor, tmp, rou_k_2

    a(:) = box(1,:)
    b(:) = box(2,:)
    c(:) = box(3,:)

    ! calculate the volume and the reciprocal vectors
    vol = volume( a, b, c )
    call crossproduct( a, b, kc ); kc = kc * 2*pi/vol 
    call crossproduct( b, c, ka ); ka = ka * 2*pi/vol
    call crossproduct( c, a, kb ); kb = kb * 2*pi/vol

    ! calculate the energy in real space
    Er = 0.0
    do i_mole = 1, n_mole
       do j_mole = i_mole, n_mole
          r_com_i(:) = r_com( i_mole, : ) 
          r_com_j(:) = r_com( j_mole, : )
          shift = pbc_shift( r_com_i, r_com_j, box , xyz_to_box_transform )
          dr_com = pbc_dr( r_com_i, r_com_j, shift )
          if ( abs(dr_com(1)) < ewald_cutoff .and. abs(dr_com(2)) < ewald_cutoff .and. abs(dr_com(3)) < ewald_cutoff ) then
             if ( dot_product( dr_com, dr_com ) < ewald_cutoff**2 ) then
                do i_atom = 1, n_atom( i_mole )
                   do j_atom = 1, n_atom( j_mole )
                      if ( i_mole /= j_mole) then  ! if it is not self interaction
                         rij = pbc_dr( xyz(i_mole,i_atom,:), xyz(j_mole,j_atom,:), shift ) ! Note for COM cutoff, shift values unchanged
                         norm_dr = sqrt( dot_product( rij, rij ) )
                         Er = Er + chg(i_mole,i_atom) * chg(j_mole,j_atom) * erfc(norm_dr*alpha_sqrt) / norm_dr
                      else if (i_mole == j_mole .and. i_atom < j_atom) then
                         rij = pbc_dr( xyz(i_mole,i_atom,:), xyz(j_mole,j_atom,:), shift )
                         norm_dr = sqrt( dot_product( rij, rij ) )
                         Er = Er + chg(i_mole,i_atom) * chg(j_mole,j_atom) * (erfc(norm_dr*alpha_sqrt)-1.0) / norm_dr
                      end if
                   end do
                end do
             end if
          end if
       end do
    end do

    ! set up all k vectors used in the program
    totk = 0
    do i = 0, kmax
       k1(:) = i * ka(:)
       if ( i == 0 ) then
          factor = 1.0
       else 
          factor = 2.0
       end if
       do j = -kmax, kmax
          k2(:) = j * kb(:)
          do k = -kmax, kmax
             k3(:) = k * kc(:)
             ktot(:) = k1(:)+ k2(:) + k3(:)
             ksq = dot_product( ktot, ktot )
             if ( ( i/=0 .or. j/=0 .or. k/=0 ) .and. ksq <= ksqmax ) then
                totk = totk + 1
                if ( totk > maxk ) stop 'Maxk is too small, set another value!'
                ! decide which global data arrays to update, depending on whether this is a volume move
                Select Case( ewald_volume_move )
                Case("no")
                   k_vec(totk,:) = ktot(:)
                   B_factor( totk ) = 2*pi/vol * exp( -ksq / 4 / alpha_sqrt / alpha_sqrt ) / ksq * factor
                Case("yes")
                   new_k_vec(totk,:) = ktot(:)
                   new_B_factor( totk ) = 2*pi/vol * exp( -ksq / 4 / alpha_sqrt / alpha_sqrt ) / ksq * factor
                case default
                   stop " setting of parameter ewald_volume_move is not recognized in ewald function!"
                End Select
             end if
          end do
       end do
    end do
    ! Energy in reciprocal space
    Ek = 0.0
    do k = 1, totk
       ! calculate rou(k), again which global data arrays to update depends on whether this is a volume move
       Select Case( ewald_volume_move )
       Case("no")
          rou_k(k) = 0.0
       Case("yes")
          rou_k_try(k) = 0.0
       End Select
       do i_mole = 1, n_mole  ! loop over all atoms
          do i_atom = 1, n_atom(i_mole)
             Select Case( ewald_volume_move )
             Case("no")
                tmp = dot_product(k_vec(k,:),xyz(i_mole,i_atom,:)) 
                rou_k(k) = rou_k(k) + chg(i_mole,i_atom) * cmplx( cos(tmp), sin(tmp) )
             Case("yes")
                tmp = dot_product(new_k_vec(k,:),xyz(i_mole,i_atom,:)) 
                rou_k_try(k) = rou_k_try(k) + chg(i_mole,i_atom) * cmplx( cos(tmp), sin(tmp) )
             End Select
          end do
       end do
       Select Case( ewald_volume_move )
       Case("no")
          rou_k_2 = real(rou_k(k))**2 + aimag(rou_k(k))**2
          Ek = Ek + B_factor(k) * rou_k_2
       Case("yes")
          rou_k_2 = real(rou_k_try(k))**2 + aimag(rou_k_try(k))**2
          Ek = Ek + new_B_factor(k) * rou_k_2
       End Select
    end do
    ! Calculate total energy
    !write(*,'(A,2F12.6)') 'Ewald:', Er * 0.52914 * 627.51 * 4.184, Ek * 0.52914 * 627.51 * 4.184
    ewald = Ek  + Er + Ewald_self
    ewald = ewald * 0.52914 * 627.51 * 4.184 ! convert from e^2/A to kJ/mol
  end function ewald


  !*********************************************************************
  ! This function return new Ewald energy based on previous caculation: Eold
  ! for a move in which molecule i_move is moved, but all other molecules remain in the 
  ! configurations used to calculate Eold
  ! 
  ! Global variables this function updates
  !  complex*8, dimension(maxk) :: rou_k_try
  !
  !  Global variables this function uses, but DOESNT update
  !  integer :: totk
  !  complex*8, dimension(maxk) :: rou_k  ! Charge density in rec. space
  !  real*8, dimension(maxk,3) :: k_vec   ! All k vectors, no need to change if box size is not changed
  !  real*8, dimension(maxk) :: B_factor  ! B-factors associated with k_vec
  !*********************************************************************
  real*8 function update_ewald( Eold, n_mole, n_atom, xyz, chg, box, i_move, xyz_new, mass, r_com )
    use global_variables
    implicit none
    real*8, intent(in) :: Eold
    integer, intent(in) :: n_mole, i_move
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:,:) :: xyz
    real*8, intent(in), dimension(:,:) :: chg, box, xyz_new, mass, r_com
    real*8 :: dEr,dEk, tmp, norm_dr, rou_k_2, rou_k_2_prev
    real*8, dimension(3) :: rij, shift, r_com_i, r_com_j, dr_com
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: xyz_tmp
    integer :: k, i_atom, j_mole, j_atom
    xyz_tmp = xyz
    xyz_tmp(i_move,:,:) = xyz_new(:,:)
    dEr = 0.0

    ! Real space energy
    do j_mole = 1, n_mole
       if ( j_mole /= i_move ) then
          r_com_j(:) = r_com( j_mole, : ) 

          r_com_i(:) = r_com( i_move, : )
          shift(:) = pbc_shift( r_com_i, r_com_j, box , xyz_to_box_transform )
          dr_com(:) = pbc_dr( r_com_i, r_com_j, shift )
          if ( abs(dr_com(1)) < ewald_cutoff .and. abs(dr_com(2)) < ewald_cutoff .and. abs(dr_com(3)) < ewald_cutoff ) then
             if ( dot_product( dr_com, dr_com ) < ewald_cutoff**2 ) then
                do i_atom = 1, n_atom( i_move )
                   do j_atom = 1, n_atom( j_mole )
                      rij = pbc_dr( xyz(i_move,i_atom,:), xyz(j_mole,j_atom,:), shift )
                      norm_dr = sqrt( dot_product( rij, rij ) )
                      dEr = dEr - chg(i_move,i_atom) * chg(j_mole,j_atom) * erfc(norm_dr*alpha_sqrt) / norm_dr
                   end do
                end do
             end if
          end if

          r_com_i(:) = pos_com( xyz_tmp, i_move, n_atom, mass )
          shift(:) = pbc_shift( r_com_i, r_com_j, box , xyz_to_box_transform )
          dr_com(:) = pbc_dr( r_com_i, r_com_j, shift )
          if ( abs(dr_com(1)) < ewald_cutoff .and. abs(dr_com(2)) < ewald_cutoff .and. abs(dr_com(3)) < ewald_cutoff ) then
             if ( dot_product( dr_com, dr_com ) < ewald_cutoff**2 ) then
                do i_atom = 1, n_atom( i_move )
                   do j_atom = 1, n_atom( j_mole )
                      rij = pbc_dr( xyz_tmp(i_move,i_atom,:), xyz_tmp(j_mole,j_atom,:), shift )
                      norm_dr = sqrt( dot_product( rij, rij ) )
                      dEr = dEr + chg(i_move,i_atom) * chg(j_mole,j_atom) * erfc(norm_dr*alpha_sqrt) / norm_dr
                   end do
                end do
             end if
          end if

       end if
    end do

    ! Reciprocal space
    dEk = 0.0
    rou_k_try(1:totk) = rou_k(1:totk)
    do k = 1, totk
       rou_k_2_prev = real(rou_k_try(k))**2 + aimag(rou_k_try(k))**2
       ! Update rou_k_try
       do i_atom = 1, n_atom(i_move)
          tmp = dot_product(k_vec(k,:),xyz(i_move,i_atom,:))
          rou_k_try(k) = rou_k_try(k) - chg(i_move,i_atom) * cmplx( cos(tmp), sin(tmp) ) ! Substract previous one
          tmp = dot_product(k_vec(k,:),xyz_tmp(i_move,i_atom,:))
          rou_k_try(k) = rou_k_try(k) + chg(i_move,i_atom) * cmplx( cos(tmp), sin(tmp) ) ! Add new one
       end do
       rou_k_2 = real(rou_k_try(k))**2 + aimag(rou_k_try(k))**2
       dEk = dEk + B_factor(k) * ( rou_k_2 - rou_k_2_prev )
    end do
    !write(*,'(A,2F12.6)') 'Update_ewald:', dEr * 0.52914 * 627.51 * 4.184, dEk * 0.52914 * 627.51 * 4.184
    update_ewald = Eold + (dEr+dEk) * 0.52914 * 627.51 * 4.184

  end function update_ewald

  !***************************************************************************
  ! This function update Ewald self interaction based on current configuration
  ! It updates:
  ! real*8 :: Ewald_self  
  !
  ! Note that Ewald self corrects for ONE HALF of the interaction between a charge
  ! interacting with it's own Gaussian.  This is because in the Ewald reciprocal sum,
  ! There is a factor of one half to compensate for overcounting, and this factor
  ! therefore carries over for the self interaction
  !
  ! the total interaction of a charge with its Gaussian is 2 * (q1)^2 * ( alpha / pi ) ^ (1/2)
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



  !*********************************************************************
  ! This function calculates electrostatic energy with cutoff
  ! It has been updated to include drude_oscillators, screened electrostatics, and
  ! framework atoms
  ! Global variables used:
  !   real*8 :: Electro_cutoff, springcon, thole
  !
  ! note that while other energy subroutines ( pme real space, lennard jones ), apply
  ! a cutoff on an atom-by-atom basis, as a molecule-molecule cutoff fails for large molecules,
  ! we can't do that here, because splitting up molecules essentially creates ions which 
  ! screws up the electrostatic energy.  In PME, it's okay, because we always have the
  ! reciprocal space contribution, but here its not.
  !*********************************************************************
  real*8 function E_elec_cutoff( tot_n_mole,n_mole, n_atom,chg, xyz, box, r_com ) 
    use global_variables
    use omp_lib
    implicit none
    integer, intent(in) :: n_mole,tot_n_mole
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: chg, box, r_com
    real*8, intent(in), dimension(:,:,:) :: xyz
    integer :: i_mole, i_atom, j_mole, j_atom
    real*8, dimension(3) :: dr, r_com_i, r_com_j, dr_com, shift
    real*8 :: norm_dr,Electro_cutoff_use,E_intra
    E_elec_cutoff = 0.0

    !Electro_cutoff_use = min(box(1,1)/2.,box(2,2)/2.,box(3,3)/2.,Electro_cutoff)
    Electro_cutoff_use = Electro_cutoff

    call OMP_SET_NUM_THREADS(n_threads)
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(tot_n_mole,n_mole,r_com,box,n_atom,xyz,chg,Electro_cutoff_use, screen_type, thole, springcon, xyz_to_box_transform, molecule_index) REDUCTION(+:E_elec_cutoff)
    !$OMP DO SCHEDULE(DYNAMIC,n_mole/8 + 1)
    do i_mole = 1, n_mole
       do j_mole = i_mole, tot_n_mole
          r_com_i(:) = r_com( i_mole, : )
          r_com_j(:) = r_com( j_mole, : )
          shift(:) = pbc_shift( r_com_i, r_com_j, box , xyz_to_box_transform )
          dr_com(:) = pbc_dr( r_com_i, r_com_j, shift )
          if ( dot_product( dr_com, dr_com ) < Electro_cutoff_use**2 ) then
             if ( i_mole /= j_mole ) then  ! if it is not self interaction
                do i_atom = 1, n_atom(i_mole)
                   do j_atom = 1, n_atom(j_mole)
                      dr = pbc_dr( xyz(i_mole,i_atom,:), xyz(j_mole,j_atom,:), shift )
                      norm_dr = sqrt( dot_product( dr, dr ) )
                      E_elec_cutoff = E_elec_cutoff + chg(i_mole,i_atom)*chg(j_mole,j_atom)* screen(i_mole,j_mole,i_atom,j_atom,norm_dr) / norm_dr
                   end do
                end do
             else if (i_mole == j_mole) then
                if(n_atom(i_mole) .gt. 1) then
                   do i_atom=1,n_atom(i_mole)-1
                      do j_atom=i_atom+1,n_atom(i_mole)  
                         call intra_elec_energy(E_intra,xyz,chg,i_mole,i_atom,j_atom,n_atom, molecule_index)
                         E_elec_cutoff = E_elec_cutoff + E_intra
                      enddo
                   enddo
                endif
             end if
          end if
       end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    E_elec_cutoff = E_elec_cutoff * 0.52914*627.51*4.184 ! Convert to kJ/mol
  end function E_elec_cutoff




  !**************************************************************
  ! This subroutine calculates electrostatic energy of insertion of a molecule using cutoff
  ! this subroutine should not be used with drude oscillators, for two reasons.
  ! The first is that the positions of drude oscillators should be solved for
  ! self consistently, so and Order(N) operation can't be carried out anyway
  ! The second reason is that there are no intra molecule interactions here
  !**************************************************************
  subroutine update_elec_cutoff_ins( Energy, tot_n_mole, n_atom, chg, xyz, box, i_move, r_com )
    use global_variables
    implicit none
    real*8, intent(out) :: Energy
    integer, intent(in) :: tot_n_mole, i_move
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: chg, box,r_com
    real*8, intent(in), dimension(:,:,:) :: xyz  
    real*8 :: norm_dr, Electro_cutoff_use
    real*8, dimension(3) :: r_com_i, r_com_j, dr_com, dr, shift
    integer :: i_atom, j_mole, j_atom


    ! this shouldn't be called until modified, since we've changed E_elec_cutoff to shift based on atom-atom basis
    write(*,*) "update_elec_cutoff_ins needs to be modified to be consistent with E_elec_cutoff which has been changed"
    write(*,*) "to use atom-atom PBC shifts"
    stop

    Energy=0.0
    !Electro_cutoff_use = min(box(1,1)/2.,box(2,2)/2.,box(3,3)/2.,Electro_cutoff)
    Electro_cutoff_use = Electro_cutoff

    do j_mole = 1, tot_n_mole
       if ( j_mole /= i_move ) then
          r_com_j(:) = r_com( j_mole, : )
          r_com_i(:) = r_com( i_move, : )
          shift(:) = pbc_shift( r_com_i, r_com_j, box , xyz_to_box_transform )
          dr_com(:) = pbc_dr( r_com_i, r_com_j, shift )
          if ( abs(dr_com(1)) < Electro_cutoff_use .and. abs(dr_com(2)) < Electro_cutoff_use .and. abs(dr_com(3)) < Electro_cutoff_use ) then
             if ( dot_product( dr_com, dr_com ) < Electro_cutoff_use**2 ) then
                do i_atom = 1, n_atom( i_move ) 
                   do j_atom = 1, n_atom( j_mole )
                      dr = pbc_dr( xyz(i_move,i_atom,:), xyz(j_mole,j_atom,:), shift )
                      norm_dr = sqrt( dot_product( dr, dr ) )
                      Energy = Energy + chg(i_move,i_atom)*chg(j_mole,j_atom)* screen(i_move,j_mole,i_atom,j_atom,norm_dr)/norm_dr
                   end do
                end do
             end if
          end if

       end if
    end do
    Energy = Energy * 0.52914*627.51*4.184 ! Convert to kJ/mol

  end subroutine update_elec_cutoff_ins





  !*********************************************************************
  ! This function calculates electrostatic forces on the i_atom for i_mole through an ewald sum
  !
  ! NOTE: the code for this function is very inefficient, but because it was only used to test the
  ! PME force code, and never used in simulations, no effort was taken to improve it
  !
  ! input re_calc_kvec if set to "1" makes routine recalculate the kvectors and Bfac, which is
  ! necessary if a volume move has been made and no previous call to either ewald, or ewald_force has been
  ! used to update these quantities
  !
  ! input re_calc_rouk if set to "1" makes routine recalculate rou_k, which is necessary if displacement moves
  ! have been made and no previous call to ewald_force, or ewald (or update_ewald) has been used to 
  ! update rou_k
  !*********************************************************************
  function ewald_force( n_mole, n_atom, xyz, chg, box, r_com, i_mole,i_atom, re_calc_kvec, re_calc_rouk)
    use global_variables
    implicit none
    real*8,dimension(3)::ewald_force
    integer, intent(in) :: n_mole,i_atom,i_mole,re_calc_kvec,re_calc_rouk
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: chg, box, r_com
    real*8, intent(in), dimension(:,:,:) :: xyz
    integer :: i, j, k, N, j_mole, j_atom
    real*8, dimension(3) :: a, b, c, ri, rj, rij, k1, k2, k3, ka, kb, kc,  ktot, shift, r_com_i, r_com_j, dr_com,Fk
    real*8 ::  qi, qj, norm_dr
    real*8 :: vol, ksq, factor, tmp


    ! calculate the force contribution in real space
    ewald_force = 0.0
    do j_mole = 1, n_mole
       r_com_i(:) = r_com( i_mole, : ) 
       r_com_j(:) = r_com( j_mole, : )
       shift = pbc_shift( r_com_i, r_com_j, box , xyz_to_box_transform )
       dr_com = pbc_dr( r_com_i, r_com_j, shift )
       if ( abs(dr_com(1)) < ewald_cutoff .and. abs(dr_com(2)) < ewald_cutoff .and. abs(dr_com(3)) < ewald_cutoff ) then
          if ( dot_product( dr_com, dr_com ) < ewald_cutoff**2 ) then
             do j_atom = 1, n_atom( j_mole )
                if ( i_mole /= j_mole) then  ! if it is not self interaction
                   rij = -pbc_dr( xyz(i_mole,i_atom,:), xyz(j_mole,j_atom,:), shift ) ! Note for COM cutoff, shift values unchanged
                   norm_dr = sqrt( dot_product( rij, rij ) )
                   ewald_force = ewald_force + chg(i_mole,i_atom) * chg(j_mole,j_atom) * rij * erfc(norm_dr*alpha_sqrt) / norm_dr**3
                   ewald_force = ewald_force + chg(i_mole,i_atom) * chg(j_mole,j_atom) *(2.*alpha_sqrt/sqrt(pi))*exp(-(alpha_sqrt*norm_dr)**2)* rij/norm_dr**2
                else if (i_mole == j_mole .and. i_atom /= j_atom) then
                   rij = -pbc_dr( xyz(i_mole,i_atom,:), xyz(j_mole,j_atom,:), shift ) ! Note for COM cutoff, shift values unchanged
                   norm_dr = sqrt( dot_product( rij, rij ) )
                   ewald_force = ewald_force + chg(i_mole,i_atom) * chg(j_mole,j_atom) * rij * (erfc(norm_dr*alpha_sqrt)-1.0) / norm_dr**3
                   ewald_force = ewald_force + chg(i_mole,i_atom) * chg(j_mole,j_atom) *(2.*alpha_sqrt/sqrt(pi))*exp(-(alpha_sqrt*norm_dr)**2)* rij/norm_dr**2
                end if
             end do
          end if
       end if
    end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  if re_calc_kvec is 1, need to recalculate kvec,Bfac
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(re_calc_kvec .eq.1) then

       a(:) = box(1,:)
       b(:) = box(2,:)
       c(:) = box(3,:)

       ! calculate the volume and the reciprocal vectors
       vol = volume( a, b, c )
       call crossproduct( a, b, kc ); kc = kc * 2*pi/vol 
       call crossproduct( b, c, ka ); ka = ka * 2*pi/vol
       call crossproduct( c, a, kb ); kb = kb * 2*pi/vol

       ! set up all k vectors using in the program
       totk = 0
       do i = 0, kmax
          k1(:) = i * ka(:)
          if ( i == 0 ) then
             factor = 1.0
          else 
             factor = 2.0
          end if
          do j = -kmax, kmax
             k2(:) = j * kb(:)
             do k = -kmax, kmax
                k3(:) = k * kc(:)
                ktot(:) = k1(:)+ k2(:) + k3(:)
                ksq = dot_product( ktot, ktot )
                if ( ( i/=0 .or. j/=0 .or. k/=0 ) .and. ksq <= ksqmax ) then
                   totk = totk + 1
                   if ( totk > maxk ) stop 'Maxk is too small, set another value!'
                   k_vec(totk,:) = ktot(:)
                   B_factor( totk ) = 2*pi/vol * exp( -ksq / 4 / alpha_sqrt / alpha_sqrt ) / ksq * factor
                end if
             end do
          end do
       end do

    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  if re_calc_rouk is 1, need to recalculate rouk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(re_calc_rouk.eq.1) then

       do k = 1, totk
          ! calculate rou(k)
          rou_k(k) = 0.0
          do j_mole = 1, n_mole  ! loop over all atoms
             do j_atom = 1, n_atom(j_mole)
                tmp = dot_product(k_vec(k,:),xyz(j_mole,j_atom,:)) 
                rou_k(k) = rou_k(k) + chg(j_mole,j_atom) * cmplx( cos(tmp), sin(tmp) )
             end do
          end do
       end do
    endif



    Fk = 0.0
    do k = 1, totk
       tmp = dot_product(k_vec(k,:),xyz(i_mole,i_atom,:)) 
       Fk(:) = Fk(:) + B_factor(k) * (-2.*aimag(rou_k(k))*chg(i_mole,i_atom)*k_vec(k,:)*cos(tmp)+2.*real(rou_k(k))*chg(i_mole,i_atom)*k_vec(k,:)*sin(tmp))
    end do

    ewald_force=(ewald_force+Fk)  * 0.52914 * 627.51 * 4.184 ! convert to kJ/mol

  end function ewald_force

  !*********************************************************************
  !  SEE COMMENTS TO ABOVE EWALD_FORCE FUNCTION
  !
  !  This function updates electrostatic forces on the i_atom for i_mole through an ewald sum
  !  for a change in position of i_move,for i_move.ne.i_mole.  If i_move eq i_mole , call ewald force instead
  !*********************************************************************
  function update_ewald_force( old_force,n_mole, n_atom, xyz, chg, box, mass, r_com, i_mole,i_atom,i_move,xyz_new)
    use global_variables
    implicit none
    real*8,dimension(3)::update_ewald_force
    real*8,dimension(3),intent(in)::old_force
    integer, intent(in) :: n_mole,i_atom,i_mole,i_move
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: chg, box, mass, r_com,xyz_new
    real*8, intent(in), dimension(:,:,:) :: xyz
    integer :: i, j, k, N, j_atom,j_mole
    real*8, dimension(3) :: ri, rj, rij, shift, r_com_i, r_com_j, dr_com,Fk,dE
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: xyz_tmp
    complex*8, dimension(maxk) ::rou_k_temp
    real*8 ::  qi, qj, norm_dr
    real*8 :: vol, ksq, factor, tmp

    xyz_tmp = xyz
    xyz_tmp(i_move,:,:) = xyz_new(:,:)

    j_mole=i_move
    ! calculate the force contribution in real space
    dE=0.
    r_com_i(:) = r_com( i_mole, : ) 
    r_com_j(:) = r_com( j_mole, : )
    shift = pbc_shift( r_com_i, r_com_j, box, xyz_to_box_transform )
    dr_com = pbc_dr( r_com_i, r_com_j, shift )
    if ( abs(dr_com(1)) < ewald_cutoff .and. abs(dr_com(2)) < ewald_cutoff .and. abs(dr_com(3)) < ewald_cutoff ) then
       if ( dot_product( dr_com, dr_com ) < ewald_cutoff**2 ) then
          do j_atom = 1, n_atom( j_mole )
             rij = -pbc_dr( xyz(i_mole,i_atom,:), xyz(j_mole,j_atom,:), shift ) ! Note for COM cutoff, shift values unchanged
             norm_dr = sqrt( dot_product( rij, rij ) )
             dE = dE - chg(i_mole,i_atom) * chg(j_mole,j_atom) * rij * erfc(norm_dr*alpha_sqrt) / norm_dr**3
             dE = dE - chg(i_mole,i_atom) * chg(j_mole,j_atom) *(2.*alpha_sqrt/sqrt(pi))*exp(-(alpha_sqrt*norm_dr)**2)* rij/norm_dr**2
          end do
       end if
    end if

    r_com_j(:) = pos_com( xyz_tmp, j_mole, n_atom, mass )
    shift = pbc_shift( r_com_i, r_com_j, box , xyz_to_box_transform )
    dr_com = pbc_dr( r_com_i, r_com_j, shift )
    if ( abs(dr_com(1)) < ewald_cutoff .and. abs(dr_com(2)) < ewald_cutoff .and. abs(dr_com(3)) < ewald_cutoff ) then
       if ( dot_product( dr_com, dr_com ) < ewald_cutoff**2 ) then
          do j_atom = 1, n_atom( j_mole )
             rij = -pbc_dr( xyz(i_mole,i_atom,:), xyz_tmp(j_mole,j_atom,:), shift ) ! Note for COM cutoff, shift values unchanged
             norm_dr = sqrt( dot_product( rij, rij ) )
             dE = dE + chg(i_mole,i_atom) * chg(j_mole,j_atom) * rij * erfc(norm_dr*alpha_sqrt) / norm_dr**3
             dE = dE + chg(i_mole,i_atom) * chg(j_mole,j_atom) *(2.*alpha_sqrt/sqrt(pi))*exp(-(alpha_sqrt*norm_dr)**2)* rij/norm_dr**2
          end do
       end if
    end if

    rou_k_temp=rou_k
    Fk = 0.0
    do k = 1, totk
       ! Update rou_k_temp
       do j_atom = 1, n_atom(j_mole)
          tmp = dot_product(k_vec(k,:),xyz(j_mole,j_atom,:))
          rou_k_temp(k) = rou_k_temp(k) - chg(j_mole,j_atom) * cmplx( cos(tmp), sin(tmp) ) ! Substract previous one
          tmp = dot_product(k_vec(k,:),xyz_tmp(j_mole,j_atom,:))
          rou_k_temp(k) = rou_k_temp(k) + chg(j_mole,j_atom) * cmplx( cos(tmp), sin(tmp) ) ! Add new one
       end do
       tmp = dot_product(k_vec(k,:),xyz(i_mole,i_atom,:)) 
       Fk(:) = Fk(:) + B_factor(k) * (-2.*(aimag(rou_k_temp(k))-aimag(rou_k(k)))*chg(i_mole,i_atom)*k_vec(k,:)*cos(tmp)+2.*(real(rou_k_temp(k))-real(rou_k(k)))*chg(i_mole,i_atom)*k_vec(k,:)*sin(tmp))
    end do

    update_ewald_force=old_force+ (dE+Fk) * 0.52914 * 627.51 * 4.184 ! convert to kJ/mol

  end function update_ewald_force




  !*********************************************************************
  ! This subroutine calculate electrostatic forces with cutoff
  ! It has been updated to include drude_oscillators, screened electrostatics, and
  ! framework atoms, and thus has been changed to a subroutine rather than a function
  ! Global variables used:
  !   real*8 :: Electro_cutoff
  !*********************************************************************
  subroutine F_elec_cutoff( elec_force,target_atoms,tot_n_mole, n_mole, n_atom, chg, xyz, box, r_com) 
    use global_variables
    use omp_lib
    implicit none
    real*8,dimension(:,:,:),intent(out):: elec_force
    integer, intent(in) :: n_mole, tot_n_mole
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: chg, box, r_com
    integer, dimension(:,:),intent(in):: target_atoms
    real*8, intent(in), dimension(:,:,:) :: xyz
    integer :: i_mole, i_atom, j_mole, j_atom, t_atom, thread_id, i_thread
    real*8, dimension(3) :: rij, r_com_i, r_com_j, dr_com, shift,f_ij
    real*8, dimension(n_mole,maxval(n_atom),3) :: local_force
    real*8, dimension(n_threads,n_mole,maxval(n_atom),3) :: temp_force
    real*8 :: norm_dr,Electro_cutoff_use

    elec_force = 0D0
    local_force = 0.D0
    temp_force= 0.D0
    !Electro_cutoff_use = min(box(1,1)/2.,box(2,2)/2.,box(3,3)/2.,Electro_cutoff)
    Electro_cutoff_use = Electro_cutoff

    call OMP_SET_NUM_THREADS(n_threads)
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(n_mole,tot_n_mole,target_atoms,r_com,box,n_atom,xyz,chg,temp_force,screen_type, thole, springcon, Electro_cutoff_use , xyz_to_box_transform, molecule_index)
    !$OMP CRITICAL
    local_force = 0.D0
    !$OMP END CRITICAL
    !$OMP DO SCHEDULE(DYNAMIC,n_mole/8+1)
    do i_mole = 1,n_mole
       r_com_i(:) = r_com( i_mole, : )
       do j_mole = i_mole, tot_n_mole
          r_com_j(:) = r_com( j_mole, : )
          shift(:) = pbc_shift( r_com_i, r_com_j, box , xyz_to_box_transform )
          dr_com(:) = pbc_dr( r_com_i, r_com_j, shift )
          if ( dot_product( dr_com, dr_com ) < Electro_cutoff_use**2 ) then
             if ( i_mole /= j_mole ) then  ! if it is not self interaction
                if( j_mole .le. n_mole) then  ! if we are considering solute-solute interaction, loop over all atoms
                   do i_atom = 1, n_atom( i_mole )
                      do j_atom = 1, n_atom( j_mole )
                         rij = -pbc_dr( xyz(i_mole,i_atom,:), xyz(j_mole,j_atom,:), shift )
                         norm_dr = sqrt( dot_product( rij, rij ) )
                         !! add "screened" real space interaction
                         f_ij = chg(i_mole,i_atom) * chg(j_mole,j_atom) * rij * screen(i_mole,j_mole,i_atom,j_atom,norm_dr) / norm_dr**3
                         f_ij = f_ij - chg(i_mole,i_atom) * chg(j_mole,j_atom) * d_screen(i_mole,j_mole,i_atom,j_atom,rij,norm_dr) / norm_dr
                         local_force(i_mole,i_atom,:) = local_force(i_mole,i_atom,:) + f_ij(:)
                         local_force(j_mole,j_atom,:) = local_force(j_mole,j_atom,:) - f_ij(:)
                      enddo
                   enddo
                else                    ! we are considering solute-framework interaction, only calculate forces on the atoms you care about for the solute
                   do i_atom = 1, n_atom( i_mole )
                      if(target_atoms(i_mole,i_atom).eq.0) then
                         goto 100
                      endif
                      t_atom = target_atoms(i_mole,i_atom)
                      do j_atom = 1, n_atom( j_mole )  ! n_atom should be one here, this is framework
                         rij = -pbc_dr( xyz(i_mole,t_atom,:), xyz(j_mole,j_atom,:), shift ) ! Note for COM cutoff, shift values unchanged
                         norm_dr = sqrt( dot_product( rij, rij ) )
                         !! add "screened" real space interaction
                         f_ij = chg(i_mole,t_atom) * chg(j_mole,j_atom) * rij * screen(i_mole,j_mole,t_atom,j_atom,norm_dr) / norm_dr**3
                         f_ij = f_ij - chg(i_mole,t_atom) * chg(j_mole,j_atom) * d_screen(i_mole,j_mole,t_atom,j_atom,rij,norm_dr) / norm_dr
                         local_force(i_mole,t_atom,:) = local_force(i_mole,t_atom,:) + f_ij(:)
                      enddo
                   enddo
100                continue
                endif
             elseif (i_mole == j_mole) then
                if(n_atom(i_mole) .gt. 1) then
                   do i_atom=1,n_atom(i_mole)-1
                      do j_atom=i_atom+1,n_atom(i_mole)
                         call intra_elec_force(f_ij,xyz,chg,i_mole,n_atom,i_atom,j_atom, molecule_index)
                         local_force(i_mole,i_atom,:)=local_force(i_mole,i_atom,:) + f_ij(:)
                         local_force(i_mole,j_atom,:)=local_force(i_mole,j_atom,:) - f_ij(:)                         
                      enddo
                   enddo
                endif
             end if
          endif
       end do
    enddo

    !$OMP END DO NOWAIT
    thread_id = OMP_GET_THREAD_NUM()
    temp_force(thread_id+1,:,:,:)= local_force(:,:,:)
    !$OMP END PARALLEL

    do i_thread=1,n_threads
       do i_mole=1, n_mole
          do i_atom=1,n_atom(i_mole)
             elec_force(i_mole,i_atom,:)= elec_force(i_mole,i_atom,:) + temp_force(i_thread,i_mole,i_atom,:)
          enddo
       enddo
    enddo


    elec_force = elec_force  * 0.52914*627.51*4.184 ! Convert to kJ/mol

    !! reorganize force array so that this is consistent with output of pme_force
    do i_mole=1,n_mole
       do i_atom=1,n_atom(i_mole)
          if(target_atoms(i_mole,i_atom).eq.0) then
             goto 200
          endif
          elec_force(i_mole,i_atom,:)=elec_force(i_mole,target_atoms(i_mole,i_atom),:)
       enddo
200    continue
    enddo

  end subroutine F_elec_cutoff






  !*************************************************************
  !  this subroutine adds intra molecular polarization (drude oscillator) electrostatic interactions
  !  it is meant to be used with Electrostatic cutoff routine, as no reciprocal space terms are subtracted
  !
  !  Ewald and pme electrostatics should use intra_pme_energy routine instead
  !*************************************************************
  subroutine intra_elec_energy(E_intra,xyz,chg,i_mole,i_atom,j_atom,n_atom, molecule_index_local)
    use global_variables
    implicit none
    real*8, intent(out):: E_intra
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: chg
    real*8, intent(in), dimension(:,:,:) :: xyz
    integer,intent(in)::i_mole,i_atom,j_atom
    integer, dimension(:), intent(in) :: molecule_index_local

    real*8,dimension(3)::rij
    integer::i_drude,j_drude,sign_chg_i,sign_chg_j,flag_drudei,flag_drudej, flag_same_atom, i_mole_type
    real*8::norm_dr,pol1,pol2
    real*8,parameter::small=1D-8



!!!!!!!!!if oscillators present, get rid of intra-molecular fixed chg-chg interactions, but allow screened induced-dipole interactions
    call drude_atom_indices(i_mole,n_atom,i_atom,j_atom,i_drude,j_drude,sign_chg_i,sign_chg_j,flag_drudei,flag_drudej, flag_same_atom)


    ! we only consider intra molecular interactions between drude oscillators, or the corresponding drude oscillator charge on an atom
    ! also, we do not consider intra molecular interactions between a drude oscillator and it's corresponding atom

    ! flag_drudei tells us whether i_atom is either a drude oscillator, or has a drude_oscillator attached to it
    ! flag_drudej tells us whether j_atom is either a drude oscillator, or has a drude_oscillator attached to it
    ! flag_same_atom tells us whether j_atom is the drude oscillator on i_atom



    if ( (flag_drudei == 1 ) .and. (flag_drudej == 1 ) .and. (flag_same_atom == 0 ) ) then
       rij(:) = xyz(i_mole,i_atom,:) - xyz(i_mole, j_atom,:)
       norm_dr = sqrt( dot_product( rij, rij ) )
!!!!!!!!!!!!!!! add screened atom-atom, get atom chg from corresponding oscillator
       pol1=chg(i_mole,i_drude)**2/springcon; pol2=chg(i_mole,j_drude)**2/springcon
       E_intra = dble(sign_chg_i*sign_chg_j)*chg(i_mole,i_drude) * chg(i_mole,j_drude) * thole_screen(pol1,pol2,norm_dr,thole)/ norm_dr
    else
       E_intra = 0D0
    endif

    ! add real space interaction if no exclusion between these atoms
    i_mole_type = molecule_index_local(i_mole)
    ! check for exclusions 
    if ( molecule_exclusions( i_mole_type, i_atom, j_atom ) /= 1 ) then 
       rij(:) = xyz(i_mole,i_atom,:) - xyz(i_mole, j_atom,:)
       norm_dr = sqrt( dot_product( rij, rij ) )
       E_intra = E_intra + chg(i_mole,i_atom) * chg(i_mole,j_atom) / norm_dr
    end if


  end subroutine intra_elec_energy



  !*************************************************************
  !  this subroutine adds intra molecular polarization (drude oscillator) electrostatic forces
  !  it is meant to be used with Electrostatic cutoff routine, as no reciprocal space terms are subtracted
  !
  !  Ewald and pme electrostatics should use intra_pme_energy routine instead
  !*************************************************************
  subroutine intra_elec_force(f_ij,xyz,chg,i_mole,n_atom,i_atom,j_atom, molecule_index_local)
    use global_variables
    implicit none
    real*8, dimension(:), intent(out):: f_ij
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: chg
    real*8, intent(in), dimension(:,:,:) :: xyz
    integer,intent(in)::i_mole,i_atom,j_atom
    integer, dimension(:), intent(in) :: molecule_index_local

    real*8,dimension(3)::rij
    integer::i_drude,j_drude,n_pairs,sign_chg_i,sign_chg_j,flag_drudei,flag_drudej, flag_same_atom, i_mole_type
    real*8::norm_dr,pol1,pol2
    real*8,parameter::small=1D-8



!!!!!!!!!if oscillators present, get rid of intra-molecular fixed chg-chg interactions, but allow screened induced-dipole interactions
    call drude_atom_indices(i_mole,n_atom,i_atom,j_atom,i_drude,j_drude,sign_chg_i,sign_chg_j,flag_drudei,flag_drudej, flag_same_atom)

    ! we only consider intra molecular interactions between drude oscillators, or the corresponding drude oscillator charge on an atom
    ! also, we do not consider intra molecular interactions between a drude oscillator and it's corresponding atom


    ! flag_drudei tells us whether i_atom is either a drude oscillator, or has a drude_oscillator attached to it
    ! flag_drudej tells us whether j_atom is either a drude oscillator, or has a drude_oscillator attached to it
    ! flag_same_atom tells us whether j_atom is the drude oscillator on i_atom, or vice-versa  



    if ( (flag_drudei == 1 ) .and. (flag_drudej == 1 ) .and. (flag_same_atom == 0 ) ) then
       rij(:) = xyz(i_mole,i_atom,:) - xyz(i_mole, j_atom,:)
       norm_dr = sqrt( dot_product( rij, rij ) )
!!!!!!!!!!!!add screened atom-drude, get atom chg from corresponding oscillator
       pol1=chg(i_mole,i_drude)**2/springcon; pol2=chg(i_mole,j_drude)**2/springcon
       f_ij= dble(sign_chg_i*sign_chg_j)*chg(i_mole,i_drude) * chg(i_mole,j_drude) * (rij * thole_screen(pol1,pol2,norm_dr,thole)/ norm_dr**3 - d_thole_screen(pol1,pol2,rij,thole) / norm_dr)
    else
       f_ij =0D0
    endif

    ! add real space interaction if no exclusion between these atoms
    i_mole_type = molecule_index_local(i_mole)
    ! check for exclusions 
    if ( molecule_exclusions( i_mole_type, i_atom, j_atom ) /= 1 ) then 
       rij(:) = xyz(i_mole,i_atom,:) - xyz(i_mole, j_atom,:)
       norm_dr = sqrt( dot_product( rij, rij ) )
       f_ij = f_ij + chg(i_mole,i_atom) * chg(i_mole,j_atom) * rij / norm_dr**3
    end if


  end subroutine intra_elec_force


end module electrostatic


