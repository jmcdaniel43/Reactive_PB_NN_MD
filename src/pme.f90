module pme_routines
  use routines
  implicit none

  !*******************************************************************
  ! PARTICLE MESH EWALD SUBROUTINES
  !
  ! This module contains PME subroutines for energy and force
  ! pme routines use discrete fourier transforms in MKL library
  !
  ! the PME real space interactions are computed in the module
  !            pair_int_real_space.f90
  !
  ! reference for the pme algorithm is
  ! Essmann et al , J. Chem. Phys. 1995, 103, 8577-8593
  !*******************************************************************

contains

  !************************************************
  ! This calculates recprocal-space, PME electrostatic energy and forces
  ! based on Essmann paper, J. Chem. Phys. 103 (19) 1995
  ! definition of forward and backward dft in paper is reversed with definition
  ! in MKL library
  !
  !*************************************************

  subroutine pme_reciprocal_space_energy_force( system_data, atom_data, PME_data )
    use global_variables
    use MKL_DFTI
    use omp_lib
    implicit none
    type(system_data_type), intent(inout)   :: system_data
    type(atom_data_type) , intent(inout)    :: atom_data
    type(PME_data_type)  , intent(inout)    :: PME_data 

    real*8, dimension(:,:),allocatable::xyz_scale
    real*8 :: pme_Erecip, force(3), kk(3,3)
    integer:: i_atom, total_atoms, status,  K
    ! arrays begin at index 1 for convenience, 1,1,1 corresponds to 0,0,0
    complex*16,dimension(:,:,:),allocatable::FQ
    real*8,dimension(:), allocatable::q_1r
    complex*16,dimension(:), allocatable::q_1d
    integer:: split_do

    ! define local variables for convenience
    K=PME_data%pme_grid
    total_atoms = system_data%total_atoms

    call construct_reciprocal_lattice_vector(kk, system_data%box)

    ! create scaled coordinates
    allocate(xyz_scale(3,total_atoms))
    call create_scaled_direct_coordinates(xyz_scale, atom_data%xyz, total_atoms, kk, K)

    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "Q_grid started at", time
    endif
    !************************************************************************!

    ! grid_Q, use scaled coordinates
    call grid_Q(PME_data%Q_grid,atom_data%charge,xyz_scale,PME_data)

    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "Q_grid finished at", time
    endif
    !************************************************************************!

    allocate( FQ(K,K,K), q_1r(K**3), q_1d(K**3) )

    q_1r=RESHAPE(PME_data%Q_grid, (/K**3/) )
    q_1d=cmplx(q_1r,0.,16)

    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "Fourier transforms started at", time
    endif
    !************************************************************************!

    status=DftiComputeForward(PME_data%dfti_desc, q_1d)

    !*** need Finv(B*C)convoluted w/ Q
    !*** equals Finv(F(Finv(B*C) conv Q)
    !*** equals K**3*Finv(B*C*F(Q))

    FQ=RESHAPE(q_1d, (/K,K,K/) )

    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "B*C*FQ started at", time
    endif
    !************************************************************************!

    !*multiply B*C*F(Q)
    FQ=FQ*PME_data%CB

    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "B*C*FQ finished at", time
    endif
    !************************************************************************!

    !** take Finv
    q_1d=RESHAPE(FQ,(/K**3/) )

    status = DftiComputeBackward(PME_data%dfti_desc_inv, q_1d)

    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "Fourier transforms finished at", time
    endif
    !************************************************************************!

    FQ=RESHAPE(q_1d, (/K,K,K/) )
    PME_data%theta_conv_Q=dble(FQ)

    deallocate( FQ, q_1r, q_1d )

    !****** PME reciprocal space energy
    pme_Erecip=.5D0*sum((PME_data%Q_grid*PME_data%theta_conv_Q))*constants%conv_e2A_kJmol   
    system_data%E_elec =  system_data%E_elec + pme_Erecip

    ! store the reciprocal space energy for ms-evb if needed
    Select Case(ms_evb_simulation)
    Case("yes")
       PME_data%E_recip=pme_Erecip
    End Select

    !*** now compute forces on atoms
    ! decide how to split the parallel section
    if (n_threads .eq. 1 ) then
       split_do = 1
    else
       split_do = total_atoms/n_threads+1
    endif

    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "Q grid derivatives started at", time
    endif
    !************************************************************************!

    ! zero reciprocal space forces
    PME_data%force_recip=0d0

    call OMP_SET_NUM_THREADS(n_threads)
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(n_threads,total_atoms,PME_data, xyz_scale, atom_data, K, system_data,kk,split_do,constants) 
    !$OMP DO SCHEDULE(dynamic, split_do)
    do i_atom=1,total_atoms
       ! the reason we pass theta_conv_Q separately is we want the flexibility
       ! to input a temporary array for this in the ms_evb code
       call derivative_grid_Q(force, PME_data%theta_conv_Q, atom_data%charge, xyz_scale, i_atom, PME_data,kk, constants%conv_e2A_kJmol )
       PME_data%force_recip(:,i_atom)=PME_data%force_recip(:,i_atom)+force(:)
    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "Q grid derivatives finished at", time
    endif
    !************************************************************************!

    ! we've used the data structure force_recip to store these reciprocal space PME forces as we will use these in MS-EVB
    atom_data%force = atom_data%force + PME_data%force_recip

    deallocate(xyz_scale)

  end subroutine pme_reciprocal_space_energy_force

  !*****************************************************************
  ! This subroutine interpolates charges onto Q grid to be used in pme reciprocal space routines
  !*****************************************************************
  subroutine grid_Q(Q,chg,xyz, PME_data)
    use global_variables
    use omp_lib
    real*8, intent(in), dimension(:) :: chg
    real*8, intent(in), dimension(:,:) :: xyz
    type(PME_data_type), intent(in)    :: PME_data
    real*8,dimension(:,:,:),intent(out)::Q
    integer::i_atom,tot_atoms, k1,k2,k3,n1,n2,n3,nearpt(3),splindex(3)
    real*8::sum
    real*8,dimension(3)::u,arg
    integer :: split_do

    ! define local variables for convenience
    integer :: K,spline_order, spline_grid

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
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(split_do,xyz,chg,tot_atoms,PME_data) REDUCTION(+:Q)
       ! local variables for convenience
       K=PME_data%pme_grid
       spline_order=PME_data%spline_order
       spline_grid=PME_data%spline_grid
    !$OMP DO SCHEDULE(dynamic, split_do)
    do i_atom=1,tot_atoms
       u=xyz(:,i_atom)
       nearpt=floor(u)
       
       ! loop over outer index of Q grid first, to more efficiently use memory
       ! only need to go to k=0,spline_order-1, for k=spline_order, arg > spline_order, so don't consider this
       do k3=0,spline_order-1
          n3=nearpt(3)-k3
          arg(3)=u(3)-dble(n3);
          ! shift index of array storage if < 0
          if(n3<0) then
             n3=n3+K
          endif
          do k2=0,spline_order-1
             n2=nearpt(2)-k2
             arg(2)=u(2)-dble(n2)
             ! shift index of array storage if < 0
             if(n2<0) then
                n2=n2+K
             endif
             do k1=0,spline_order-1
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
                   sum=chg(i_atom)*PME_data%B6_spline(splindex(1))*PME_data%B6_spline(splindex(2))*PME_data%B6_spline(splindex(3))
                else
                   sum=chg(i_atom)*PME_data%B4_spline(splindex(1))*PME_data%B4_spline(splindex(2))*PME_data%B4_spline(splindex(3))   
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
  !*********************************
  subroutine modify_Q_grid( Q , chg, xyz, n_atom, PME_data, operation)
    use global_variables
    integer, intent(in)              :: n_atom
    real*8, intent(in), dimension(:) :: chg
    real*8, intent(in), dimension(:,:) :: xyz
    type(PME_data_type) , intent(in)   :: PME_data
    real*8,dimension(:,:,:),intent(inout)::Q
    integer, intent(in) :: operation
    integer::j,k1,k2,k3,n1,n2,n3,nearpt(3),splindex(3)
    real*8::sum
    real*8,dimension(3)::u,arg
    real*8 :: small=1D-6
    ! define local variables for convenience
    integer :: K,spline_order, spline_grid

    ! local variables for convenience
    K=PME_data%pme_grid
    spline_order=PME_data%spline_order
    spline_grid=PME_data%spline_grid

    do j=1,n_atom
       if ( abs(chg(j)) > small ) then
          u=xyz(:,j)
          nearpt=floor(u)
          ! loop over outer index of Q grid first, to more efficiently use memory
          do k3=0,spline_order-1
             n3=nearpt(3)-k3
             arg(3)=u(3)-dble(n3);
             if(n3<0) then
                n3=n3+K
             endif
             do k2=0,spline_order-1
                n2=nearpt(2)-k2
                arg(2)=u(2)-dble(n2)
                if(n2<0) then
                   n2=n2+K
                endif
                do k1=0,spline_order-1
                   n1=nearpt(1)-k1
                   arg(1)=u(1)-dble(n1);
                   if(n1<0) then
                      n1=n1+K
                   endif
                   sum=0d0
                   splindex = ceiling(arg/6.D0*dble(spline_grid))
                   if(spline_order .eq.6) then
                      sum=chg(j)*PME_data%B6_spline(splindex(1))*PME_data%B6_spline(splindex(2))*PME_data%B6_spline(splindex(3))
                   else
                      sum=chg(j)*PME_data%B4_spline(splindex(1))*PME_data%B4_spline(splindex(2))*PME_data%B4_spline(splindex(3))      
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
  subroutine derivative_grid_Q(force,FQ,chg,xyz,i_atom, PME_data, kk, conv_e2A_kJmol, flag_update)
    use global_variables
    real*8,dimension(3),intent(out)::force
    real*8,dimension(:,:,:),intent(in)::FQ
    type(PME_data_type), intent(inout)   :: PME_data
    integer, intent(in) :: i_atom
    real*8, intent(in), dimension(:) :: chg
    real*8,intent(in),dimension(3,3)::kk
    real*8, intent(in), dimension(:,:) :: xyz
    real*8, intent(in)                 :: conv_e2A_kJmol
    integer, intent(in), optional :: flag_update

    integer::i,j,k1,k2,k3,n1,n2,n3,nearpt(3)
    integer::g1n(3),g1nmin(3),g1(3),g2(3)
    real*8::fac(3),chg_i
    real*8,dimension(3)::u,arg1,arg2,force_temp
    ! temporary local variables
    integer :: K, spline_order, spline_grid

    integer :: count
    count=1

    ! local variables
    K=PME_data%pme_grid
    spline_order=PME_data%spline_order
    spline_grid=PME_data%spline_grid

    force=0.D0
    chg_i=chg(i_atom)
    ! coordinates should be stored this way for optimal memory retrieval
    u=xyz(:,i_atom)
    nearpt=floor(u)

    ! loop over K grid according to memory storage n3, n2, n1
    ! first part of fac
    do k3=0,spline_order-1
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
       do k2=0,spline_order-1
          n2=nearpt(2)-k2
          arg1(2)=u(2)-dble(n2)
          arg2(2) = arg1(2) - 1d0
          ! shift index of array storage if < 0
          if(n2<0) then
             n2=n2+K
          endif
          ! if k1=n, arg1 >= n, so we can exclude that term in the sum, also because we exlude that term in the sum, arg1 <= n , so we dont need to check the limits of the bspline evaluation
          do k1=0,spline_order-1
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
             g2=ceiling(arg2/real(spline_order-1)*real(spline_grid))
             g1n=ceiling(arg1/real(spline_order)*real(spline_grid))
             g1nmin = ceiling(arg1/real(spline_order-1)*real(spline_grid))

!!!!!force in x direction
             g1=g1n;g1(1)=g1nmin(1);
             if(arg1(1)<real(spline_order-1)) then 
                if(spline_order.eq.6) then
                   fac(1)=chg_i*(PME_data%B5_spline(g1(1))*PME_data%B6_spline(g1(2))*PME_data%B6_spline(g1(3)))
                else
                   fac(1)=chg_i*(PME_data%B3_spline(g1(1))*PME_data%B4_spline(g1(2))*PME_data%B4_spline(g1(3)))  
                endif
             endif
             ! at this point, arg1(1) < n , so arg2(1) < n - 1 , so don't check for that
             !if( (arg2(1)<real(n-1) ) .and.(0.<arg2(1)) ) then 
             if(0.<arg2(1)) then
                if(spline_order.eq.6) then
                   fac(1)=fac(1)+chg_i*(-PME_data%B5_spline(g2(1))*PME_data%B6_spline(g1(2))*PME_data%B6_spline(g1(3)))
                else
                   fac(1)=fac(1)+chg_i*(-PME_data%B3_spline(g2(1))*PME_data%B4_spline(g1(2))*PME_data%B4_spline(g1(3)))  
                endif
             endif
!!!!!!force in y direction
             g1=g1n;g1(2)=g1nmin(2)
             if ( arg1(2)<real(spline_order-1) ) then 
                if(spline_order.eq.6) then
                   fac(2)=chg_i*(PME_data%B5_spline(g1(2))*PME_data%B6_spline(g1(1))*PME_data%B6_spline(g1(3)))
                else
                   fac(2)=chg_i*(PME_data%B3_spline(g1(2))*PME_data%B4_spline(g1(1))*PME_data%B4_spline(g1(3)))  
                endif
             endif
             if (0.<arg2(2)) then 
                if(spline_order.eq.6) then
                   fac(2)=fac(2)+chg_i*(-PME_data%B5_spline(g2(2))*PME_data%B6_spline(g1(1))*PME_data%B6_spline(g1(3)))
                else
                   fac(2)=fac(2)+chg_i*(-PME_data%B3_spline(g2(2))*PME_data%B4_spline(g1(1))*PME_data%B4_spline(g1(3)))  
                endif
             endif
!!!!!force in z direction
             g1=g1n;g1(3)=g1nmin(3)
             if(arg1(3)<real(spline_order-1)) then 
                if(spline_order.eq.6) then
                   fac(3)=chg_i*(PME_data%B5_spline(g1(3))*PME_data%B6_spline(g1(1))*PME_data%B6_spline(g1(2)))
                else
                   fac(3)=chg_i*(PME_data%B3_spline(g1(3))*PME_data%B4_spline(g1(1))*PME_data%B4_spline(g1(2)))  
                endif
             endif
             if (0.<arg2(3))  then 
                if(spline_order.eq.6) then
                   fac(3)=fac(3)+chg_i*(-PME_data%B5_spline(g2(3))*PME_data%B6_spline(g1(1))*PME_data%B6_spline(g1(2)))
                else
                   fac(3)=fac(3)+chg_i*(-PME_data%B3_spline(g2(3))*PME_data%B4_spline(g1(1))*PME_data%B4_spline(g1(2)))  
                endif
             endif

             force=force + fac * FQ(n1+1,n2+1,n3+1) * conv_e2A_kJmol  ! convert e^2/Ang to kJ/mol

	     ! if we're storing dQ_dr,
             if ( .not. present(flag_update) ) then
                Select Case(ms_evb_simulation)
                Case("yes")
                   PME_data%dQ_dr(:,count,i_atom) = fac(:)
                   PME_data%dQ_dr_index(1,count,i_atom) =n1+1
                   PME_data%dQ_dr_index(2,count,i_atom) =n2+1
                   PME_data%dQ_dr_index(3,count,i_atom) =n3+1
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
       tmp=2.D0*constants%pi*dble(m*i)/dble(K)
       sum=sum+B_spline(dble(i+1),n)*cmplx(cos(tmp),sin(tmp))
!!$     sum=sum+B6_spline(dble(i+1)/6.*dble(spline_grid))*cmplx(cos(tmp),sin(tmp))
    enddo
    bm=1.D0/sum
    bm_sq=dble(bm)**2+aimag(bm)**2

  end function bm_sq

  !***************************************************************************
  ! This function updates Ewald self interaction energy
  !
  ! Note that Ewald self corrects for ONE HALF of the interaction between a charge
  ! interacting with it's own Gaussian.  This is because in the Ewald reciprocal
  ! sum, there is a factor of one half to compensate for overcounting, and this
  ! factor therefore carries over for the self interaction
  !
  ! the total interaction of a charge with its Gaussian is 2 * (q1)^2 * ( alpha / pi ) ^ (1/2)
  ! so that's why the self correction is (q1)^2 * ( alpha / pi ) ^ (1/2)
  ! **************************************************************************
  subroutine update_Ewald_self( total_atoms , atom_data , PME_data )
    use global_variables
    implicit none
    integer, intent(in) :: total_atoms
    type(atom_data_type) :: atom_data
    type(PME_data_type) :: PME_data
    integer :: i_atom

    PME_data%Ewald_self = 0.0
    do i_atom = 1, total_atoms
       PME_data%Ewald_self = PME_data%Ewald_self - atom_data%charge(i_atom)**2
    end do
    PME_data%Ewald_self = PME_data%Ewald_self * PME_data%alpha_sqrt/constants%pi_sqrt

    ! convert from e^2/A to kJ/mol energy units
    PME_data%Ewald_self = PME_data%Ewald_self * constants%conv_e2A_kJmol

  end subroutine update_Ewald_self

  !***************************************************************************
  ! This function updates the CB_array in the PME_data structure due to volume
  ! changes
  !
  ! If a change to the periodic box is made at some point in the code, the PME
  ! data must be updated to maintain integrity of the electrostatics energy
  ! **************************************************************************
  subroutine periodic_box_change(system_data, PME_data, verlet_list_data, atom_data, molecule_data)
      use global_variables
      Type(system_data_type),intent(inout)                :: system_data
      Type(PME_data_type), intent(inout)                  :: PME_data
      type(verlet_list_data_type), intent(inout) :: verlet_list_data
      Type(atom_data_type),intent(inout)                  :: atom_data
      Type(molecule_data_type),dimension(:),intent(in)    :: molecule_data

      integer :: flag_junk

      real*8 :: a(3), b(3), c(3), ka(3), kb(3), kc(3),kk(3,3)

      system_data%volume = system_data%box(1,1)**3
      call initialize_non_orth_transform(system_data%box, system_data%xyz_to_box_transform)

      ! compute CB array
      a(:) = system_data%box(1,:); b(:) = system_data%box(2,:); c(:) = system_data%box(3,:)
      call crossproduct( a, b, kc ); kc = kc / system_data%volume
      call crossproduct( b, c, ka ); ka = ka / system_data%volume
      call crossproduct( c, a, kb ); kb = kb / system_data%volume
      kk(1,:)=ka(:);kk(2,:)=kb(:);kk(3,:)=kc(:)
      
      call CB_array(PME_data%CB,PME_data%alpha_sqrt,system_data%volume,PME_data%pme_grid,kk,PME_data%spline_order)

      call allocate_verlet_list(verlet_list_data, system_data%total_atoms, system_data%volume)
      call construct_verlet_list( verlet_list_data, atom_data, molecule_data, system_data%total_atoms, system_data%box, system_data%xyz_to_box_transform  )

      ! the "1" input to update_verlet_displacements signals to initialize the displacement array
      call update_verlet_displacements( system_data%total_atoms, atom_data%xyz, verlet_list_data , system_data%box, system_data%xyz_to_box_transform , flag_junk, 1 )

  end subroutine periodic_box_change

end module pme_routines
