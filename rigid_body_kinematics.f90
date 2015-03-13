module rigid_body
  use routines

  !***********************************
  ! these subroutines are for the manipulation of rigid polyatomic molecules
  ! The orientation of such molecules are described by quaternions, see either
  ! Allen and Tildesley "Computer Simulation of Liquids" Oxford Press, 1987
  ! or Evans, D.J. Molecular Physics, 1977, 34, 317
  ! the convention of euler angles is the same in both references, and is the same
  ! as Goldstein, Classical Mechanics, Addison-Wesley Press, 1980
  ! the correspondance between quaternion notation in these papers is
  ! q0-> Chi , q1-> Eta , q2-> -(Xi) , q3-> Zeta
  !
  ! we use the q0,q1,q2,q3 convention.
  !
  ! For a reversible, numerical integration scheme, we use the method of 
  ! Matubayasi N. and Nakahara M. J Chem. Phys. 110, 1999, 3291
  ! which has the same notation as Allen and Tildesley
  !
  ! we use the notation in Allen and Tildesley for the rotation matrix from space fixed (global) to body
  ! coordinates, namely the matrix "a" transforms a vector from space fixed to body coordinates
  ! whereas "a_trnspose" transforms from body coordinates to space coordinates
  !************************************



contains


  !************************************
  !  This subroutine carries out the steps for reversible integration of Euler's equations of 
  !  motion for a rigid body.  Velocity Verlet method is used for linear molecules, and the
  !  method of Matubayasi N. and Nakahara M. J Chem. Phys. 110, 1999, 3291 is used to 
  !  reversibly integrate the quaternions and angular velocities for non-linear molecules
  !
  ! note the suffix _q_w, indicates that both positions and velocities are being updated
  ! compared to subroutine reversible_rotate_w, which is the last stage in the reversible integration
  ! updating velocities and angular momentum using final forces and torques
  !
  ! make sure all data input arrays are declared intent(inout) since we are only passing portions of these
  !************************************
  subroutine reversible_rotate_q_w(i_mole,xyz,n_atom,r_com,quaternion,ang_vel,torque,delta_t,conv_fac,sample_vel)
    use global_variables
    implicit none
    integer,intent(in)::i_mole
    integer,intent(in)::n_atom,sample_vel
    real*8,intent(in)::delta_t,conv_fac
    real*8,dimension(:,:),intent(inout)::xyz
    real*8,dimension(:),intent(inout)::quaternion
    real*8,dimension(3),intent(inout)::r_com,torque
    real*8,dimension(3),intent(inout)::ang_vel

    integer::i_atom,index, linear_flag,j
    real*8,parameter::small=1D-4
    real*8,dimension(3)::delta_x,inertia_local
    real*8::Inertia
    real*8,dimension(3,3) :: a


    ! get inertia tensor, and determine if molecule is linear
    index = molecule_index(i_mole)
    linear_flag = molecule_shape(index)
    inertia_local(:) = molecule_inertia(index,:)

    !******* use different method for linear and non-linear molecules

    Select Case(linear_flag)
    Case(0)
       ! non-linear molecule
       ! get transformation matrix to body coordinates
       call quaternion_transformation_matrix(a,quaternion)

       call  rotate_nonlinear_q_w(xyz,n_atom,r_com,inertia_local,ang_vel,torque,delta_t,conv_fac,a,quaternion)

    Case(1)
       ! linear molecule
       do j=1,3
          if ( inertia_local(j) > small ) then
             Inertia = inertia_local(j)
             exit
          endif
       enddo

       call rotate_linear_q_w(xyz,n_atom,r_com,Inertia,ang_vel,torque,delta_t,conv_fac,sample_vel)

    End Select

    ! at this stage, positions, (quaternions) have been completely updated to delta_t timestep
    ! angular velocities have been updated to final auxiliary angular velocities, where the final
    ! angular velocities are an integration of the torques at delta_t with the auxiliary angular velocities


  end subroutine reversible_rotate_q_w


  !************************************
  ! This is the last stage in the reversible integration
  ! updating velocities and angular momentum using final forces and torques
  !
  ! make sure all data input arrays are declared intent(inout) since we are only passing portions of these arrays
  !************************************
  subroutine reversible_rotate_w(i_mole,xyz,n_atom,r_com,quaternion,ang_vel,torque_s,delta_t,conv_fac,sample_vel)
    use global_variables
    implicit none
    integer,intent(in)::i_mole
    integer,intent(in)::n_atom,sample_vel
    real*8,intent(in)::delta_t,conv_fac
    real*8,dimension(:,:),intent(inout)::xyz
    real*8,dimension(:),intent(inout)::quaternion
    real*8,dimension(3),intent(inout)::r_com,torque_s
    real*8,dimension(3),intent(inout)::ang_vel

    integer::i_atom,index, linear_flag,j,i
    real*8,parameter::small=1D-4
    real*8,dimension(3)::delta_x,inertia_local,torque_b
    real*8::Inertia
    real*8,dimension(3,3) :: a

    ! get inertia tensor, and determine if molecule is linear
    index = molecule_index(i_mole)
    linear_flag = molecule_shape(index)
    inertia_local(:) = molecule_inertia(index,:)

    !******* use different method for linear and non-linear molecules

    Select Case(linear_flag)
    Case(0)
       ! non-linear molecule
       ! get transformation matrix to body coordinates
       call quaternion_transformation_matrix(a,quaternion)

       ! torque needs to be in principle body frame.
       torque_b = matmul(a,torque_s)


       ! angular velocity is in principle body frame, so we have Lx=Ix * wx , etc.
       ! units: inertia (g/mol*Ang^2), ang_vel(ps^-1), delta_t (ps) , torque (kJ/mol)
       ! conv_fac converts kJ/mol to A^2/ps^2*g/mol
       do i=1,3
          ang_vel(i) = ang_vel(i) + delta_t / 2d0 * torque_b(i) / inertia_local(i) * conv_fac
       enddo


    Case(1)
       ! linear molecule
       do j=1,3
          if ( inertia_local(j) > small ) then
             Inertia = inertia_local(j)
             exit
          endif
       enddo
       !   update angular velocity with torque and time step, note that these vectors are in 
       !  global frame rather than molecular frame, but this is okay as Eulers equations can be
       !  transferred back to global frame as Inertia is a scalar
       ang_vel(:) = ang_vel(:) + delta_t/2./Inertia * torque_s(:)*conv_fac

    End Select


  end subroutine reversible_rotate_w



  !*************************************
  ! this subroutine updates angular velocity and quaternions for non-linear molecules
  ! see Matubayasi N. and Nakahara M. J Chem. Phys. 110, 1999, 3291
  !
  ! make sure all data input arrays are declared intent(inout) since we are only passing portions of these arrays
  !*************************************
  subroutine rotate_nonlinear_q_w(xyz,n_atom,r_com,inertia_local,ang_vel,torque_s,delta_t,conv_fac,a,quaternion)
    integer,intent(in)::n_atom
    real*8,intent(in)::delta_t,conv_fac
    real*8,dimension(:,:),intent(inout)::xyz
    real*8,dimension(:),intent(inout)::quaternion
    real*8,dimension(3),intent(inout)::r_com,torque_s,inertia_local
    real*8,dimension(3),intent(inout)::ang_vel
    real*8,dimension(3,3),intent(in) :: a

    integer  :: i,j,i_atom
    ! here aux_angmom_0 is initial, _t2 is at delta_t/2 , and _t is at delta_t
    real*8, dimension(3) :: torque_b, aux_angmom_0 , aux_angmom_t2 , aux_angmom_t,delta_x
    real*8  :: aux_angmom_Lz1,aux_angmom_Lz2, theta, w_mag
    real*8, dimension(4,4) :: unit_4, Aw_4, Mat_prop
    real*8, dimension(4) :: quaternion_old
    real*8, dimension(3,3) :: a_old, a_new, a_new_trans

    ! first we need angular velocities at time step delta_t/2

    !************** STEP 1
    ! in the first step, the auxiliary angular momentum is formed from the original angular momentum and is updated by the original torque.

    ! torque needs to be in principle body frame.
    torque_b = matmul(a,torque_s)

    ! angular velocity is in principle body frame, so we have Lx=Ix * wx , etc.
    ! units: inertia (g/mol*Ang^2), ang_vel(ps^-1), delta_t (ps) , torque (kJ/mol)
    ! conv_fac converts kJ/mol to A^2/ps^2*g/mol , so angular momentum is in g/mol*Ang^2/ps
    do i=1,3
       aux_angmom_0(i) = inertia_local(i) * ang_vel(i) + delta_t / 2d0 * torque_b(i) * conv_fac
    enddo


    !************** STEP 2
    ! the second step is a free rotation (delta_t/2) about Ly, generating Lx(delta_t/2) and Lz1
    theta = delta_t /2d0 * ( 1d0/inertia_local(3) - 1d0/inertia_local(2) ) * aux_angmom_0(2)
    aux_angmom_t2(1) = dcos(theta) * aux_angmom_0(1) + dsin(theta) * aux_angmom_0(3)

    aux_angmom_Lz1 = -dsin(theta) * aux_angmom_0(1) + dcos(theta) * aux_angmom_0(3)


    !************** STEP 3
    ! the third step is the first part of a free rotation about Lx (delta_t/2), generating Ly(delta_t/2) and Lz(delta_t/2)
    theta = delta_t /2d0 * ( 1d0/inertia_local(3) - 1d0/inertia_local(1) ) * aux_angmom_t2(1)
    aux_angmom_t2(2) = dcos(theta) * aux_angmom_0(2) - dsin(theta) * aux_angmom_Lz1
    aux_angmom_t2(3) = dsin(theta) * aux_angmom_0(2) + dcos(theta) * aux_angmom_Lz1

    ! ********* Now we have the angular momentum at time step delta_t/2, use these to propagate
    ! ********* quaternions to time step delta_t

    do i=1,3
       ang_vel(i) = aux_angmom_t2(i) / inertia_local(i)
    enddo

    unit_4=0d0
    do i=1,4
       unit_4(i,i) = 1d0
    enddo

    Aw_4(1,1) = 0d0 ; Aw_4(2,1) = ang_vel(1) ; Aw_4(3,1) = ang_vel(2) ; Aw_4(4,1) = ang_vel(3) ;
    Aw_4(1,2) = -ang_vel(1) ; Aw_4(2,2) = 0d0 ; Aw_4(3,2) = -ang_vel(3) ; Aw_4(4,2) = ang_vel(2) ;
    Aw_4(1,3) = -ang_vel(2) ; Aw_4(2,3) = ang_vel(3) ; Aw_4(3,3) = 0d0 ; Aw_4(4,3) = -ang_vel(1) ;
    Aw_4(1,4) = -ang_vel(3) ; Aw_4(2,4) = -ang_vel(2) ; Aw_4(3,4) = ang_vel(1) ; Aw_4(4,4) = 0d0 ;

    w_mag = sqrt( dot_product( ang_vel , ang_vel ) )
    theta = delta_t /2d0 * w_mag

    Mat_prop= dcos(theta) * unit_4 + dsin(theta)/w_mag * Aw_4

    quaternion_old = quaternion
    quaternion = matmul( Mat_prop , quaternion )

    ! now let's finish evolving the auxiliary angular momentum up to delta_t , 
    ! so that the final angular momentum can be created from this auxiliary angular momentum
    ! and the torques at the new positions

    !**************** STEP 4
    ! this next step is the continuation of the rotation about Lx (delta_t/2 , for a total of delta_t rotation ) generating Ly(delta_t) and Lz2
    theta = delta_t /2d0 * ( 1d0/inertia_local(3) - 1d0/inertia_local(1) ) * aux_angmom_t2(1)
    aux_angmom_t(2) = dcos(theta) * aux_angmom_t2(2) - dsin(theta) * aux_angmom_t2(3)
    aux_angmom_Lz2 =  dsin(theta) * aux_angmom_t2(2) + dcos(theta) * aux_angmom_t2(3)

    !**************** STEP 5
    ! this last step is a rotation along Ly by delta_t/2, creating final auxiliary angular momentum
    theta = delta_t /2d0 * ( 1d0/inertia_local(3) - 1d0/inertia_local(2) ) * aux_angmom_t(2)  
    aux_angmom_t(1) = dcos(theta) * aux_angmom_t2(1) + dsin(theta) * aux_angmom_Lz2
    aux_angmom_t(3) = -dsin(theta) * aux_angmom_t2(1) + dcos(theta) * aux_angmom_Lz2

    ! finally, change back to angular velocities
    do i=1,3
       ang_vel(i) = aux_angmom_t(i) / inertia_local(i)
    enddo

    ! update molecule coordinates for time step delta_t using updated quaternions at delta_t
    call quaternion_transformation_matrix(a_old,quaternion_old)
    call quaternion_transformation_matrix(a_new,quaternion)

    do i=1,3
       do j=1,3
          a_new_trans(i,j) = a_new(j,i)
       enddo
    enddo

    ! remember, matrix "a" goes from global to body coordinates, and "a_trans" is reverse

    do i_atom = 1, n_atom
       delta_x(:) = xyz(i_atom,:) - r_com(:)
       ! first transform back to body with consistent quaternion set
       delta_x = matmul(a_old, delta_x)
       ! now transform to global space coordinates with new quaternions
       delta_x = matmul(a_new_trans, delta_x)
       xyz(i_atom,:) = delta_x(:) + r_com(:)
    enddo


  end subroutine rotate_nonlinear_q_w


  !*************************************
  ! this subroutine updates angular velocity and rotates linear molecules
  !
  ! make sure all data input arrays are declared intent(inout) since we are only passing portions of these arrays
  !*************************************
  subroutine rotate_linear_q_w(xyz,n_atom,r_com,Inertia,ang_vel,torque,delta_t,conv_fac,sample_vel)
    implicit none
    integer,intent(in)::n_atom,sample_vel
    real*8,intent(in)::delta_t,conv_fac,Inertia
    real*8,dimension(:,:),intent(inout)::xyz
    real*8,dimension(3),intent(inout)::r_com,torque
    real*8,dimension(3),intent(inout)::ang_vel

    integer::i_atom
    integer,dimension(n_atom)::projection
    real*8,dimension(n_atom)::rel_x
    real*8,parameter::small=1D-4
    real*8,dimension(3)::delta_x,u_axis,u_wvec,unit,u_rotd
    real*8::theta


!!!!!!!!!!!!!!!!determine the molecular axis by the first atom not located at the center of mass
    do i_atom=1,n_atom
       delta_x(:)=xyz(i_atom,:)-r_com(:)
       if(dot_product(delta_x,delta_x) > small) exit
    enddo
    u_axis=delta_x/sqrt(dot_product(delta_x,delta_x))


!!!!!!!!!!!!!!!!!!!determine if atom projections on this axis are positive or negative
    do i_atom=1,n_atom
       delta_x(:)=xyz(i_atom,:)-r_com(:)
       rel_x(i_atom) = sqrt(dot_product(delta_x,delta_x))
       if(dot_product(delta_x,u_axis)>0.) then
          projection(i_atom)=1
       else
          projection(i_atom)=-1
       endif
    enddo

!!!!!!!!if velocities were sampled, rather than kept, then angular velocity vector needs to be adjusted
!!!!!!!!so that it is perpendicular to molecular axis (since torque must be perpendicular to axis)
    if(sample_vel.eq.1) then
       do
          call random_unit_vector(unit)
          call crossproduct(unit,u_axis,u_wvec)
          if(dot_product(u_wvec,u_wvec)>small) exit
       enddo
       u_wvec=u_wvec/sqrt(dot_product(u_wvec,u_wvec))
       ang_vel=u_wvec*sqrt(dot_product(ang_vel,ang_vel))
    endif

    !   update angular velocity with torque and time step, note that these vectors are in 
    !  global frame rather than molecular frame, but this is okay as Eulers equations can be
    !  transferred back to global frame as Inertia is a scalar

    ! units: inertia (g/mol*Ang^2), ang_vel(ps^-1), delta_t (ps) , torque (kJ/mol)
    ! conv_fac converts kJ/mol to A^2/ps^2*g/mol

    ang_vel = ang_vel + delta_t/2./Inertia * torque * conv_fac

    u_wvec = ang_vel/sqrt(dot_product(ang_vel,ang_vel))

!!!!!!!!!!!ang_vel is perpendicular to molecular axis, find rotation direction and rotate
    call crossproduct(u_wvec,u_axis,u_rotd)
    theta= delta_t*sqrt(dot_product(ang_vel,ang_vel))

    u_axis=cos(theta)*u_axis+sin(theta)*u_rotd

!!!!!!!!!!!now generate new lab coordinates for atoms
    do i_atom=1,n_atom
       xyz(i_atom,:)= r_com(:) + dble(projection(i_atom))*rel_x(i_atom)*u_axis(:)
    enddo


  end subroutine rotate_linear_q_w


  !********************************
  ! this subroutines generates principle moments of inertia for molecules
  ! in the process, the transformation matrix from global coordinates to 
  ! this principle coordinate system is found, and therefore the quaternions
  ! describing this configuration are generated
  !
  ! note we are only storing one set of principle moment of inertia vectors for
  ! each unique molecule type, so we need to make sure all same molecule type
  ! quaternions are consistent with this vector.  
  ! we do this by storing principle moments of inertia in order from smallest to 
  ! largest for each molecule, and generate quaternions consistent with this
  !********************************
  subroutine gen_inertia( xyz, i_mole, n_atom, mass ,r_com, inertia, quaternion)
    use global_variables
    implicit none
    real*8,dimension(:,:,:),intent(in)::xyz
    integer,intent(in)::i_mole
    integer,dimension(:),intent(in)::n_atom
    real*8,dimension(:,:),intent(in)::mass,r_com
    real*8,dimension(:),intent(inout)::inertia
    real*8,dimension(:),intent(inout):: quaternion

    integer::i_atom, i_mass,i_shape,index,i,j,nrot
    real*8::Ix, r_mag,r_mag2,r_vec(3),kron_delta
    real*8,parameter:: min=.1
    real*8,dimension(3,3) :: inertia_tensor,a_trans
    real*8,dimension(3) :: p_inertia
    real*8 :: small=1d-3
    real*8, dimension(4,4) :: Mat_trans
    real*8, dimension(4) :: adiag_4, quat_squared
    real*8 :: q0, q1, q2, q3
    real*8,parameter :: small_neg=-1d-8


    if ( n_atom(i_mole) > 1 ) then
       index = molecule_index( i_mole )
       i_shape = molecule_shape( index )

       Select Case(i_shape)
       Case(1)
          ! linear molecule; let internal coordinate z point along molecular axis, then Ix=Iy, Iz = 0
          Ix=0d0
          do  i_atom=1,n_atom(i_mole)
             r_vec(:) = xyz(i_mole,i_atom,:) - r_com(i_mole,:)
             r_mag = sqrt (dot_product(r_vec,r_vec))

             Ix= Ix + mass(i_mole,i_atom) * r_mag**2
          enddo

          inertia(:) = Ix
          inertia(3)=0d0

       Case(0)
          ! general, non-linear polyatomic molecule
          ! create full inertia tensor
          inertia_tensor=0d0
          do  i_atom=1,n_atom(i_mole)
             r_vec(:) = xyz(i_mole,i_atom,:) - r_com(i_mole,:) 
             r_mag2 = (dot_product(r_vec,r_vec))
             do i=1,3
                do j=1,3
                   if ( i == j ) then
                      kron_delta = 1d0
                   else
                      kron_delta = 0d0
                   endif

                   inertia_tensor(i,j) = inertia_tensor(i,j) + mass(i_mole,i_atom) * ( r_mag2 * kron_delta - r_vec(i) * r_vec(j) )

                enddo
             enddo
          enddo


          ! now diagonalize inertia tensor
          call jacobi(inertia_tensor,p_inertia, a_trans, nrot )


          ! now we need to order moments of inertia from smallest to largest, and also corresponding
          ! eigenvector array

          call order_moments_inertia(p_inertia,a_trans)
          inertia(:) = p_inertia(:)


          ! make sure all elements are non zero for non linear molecule
          do i = 1 ,3
             if ( inertia(i) < small ) then
                stop " princple inertia value is zero for a non-linear molecule, something's wrong"
             endif
          enddo

          ! make sure rotation matrix is physical and doesn't involve inversion

          if ( determinant_3(a_trans) < 0d0 ) then
             stop "determinant of rotation matrix is negative.  This is unphysical"
          endif

          ! here, the eigenvectors of the inertia tensor ( principle axis in global axis basis ) are the columns of a_trans.  Therefore, the matrix a_trans is a transformation from principle axis to global coordinates

          ! invert rotation matrix to get quaternions
          ! 
          ! (q0^2)        ( 1  1  1  1 ) (a11)
          ! (q1^2)  = 1/4 ( 1 -1 -1  1 ) (a22)
          ! (q2^2)        (-1  1 -1  1 ) (a33)
          ! (q3^2)        (-1 -1  1  1 ) ( 1 )

          Mat_trans(1,1) = 1d0 ;  Mat_trans(2,1) = 1d0 ; Mat_trans(3,1) = -1d0 ;       Mat_trans(4,1) = -1d0 ;
          Mat_trans(1,2) = 1d0 ;  Mat_trans(2,2) = -1d0 ; Mat_trans(3,2) = 1d0 ;       Mat_trans(4,2) = -1d0 ;
          Mat_trans(1,3) = 1d0 ;  Mat_trans(2,3) = -1d0 ; Mat_trans(3,3) = -1d0 ;       Mat_trans(4,3) = 1d0 ;
          Mat_trans(1,4) = 1d0 ;  Mat_trans(2,4) = 1d0 ; Mat_trans(3,4) =  1d0 ;       Mat_trans(4,4) =  1d0 ;
          Mat_trans = Mat_trans/4d0

          adiag_4(1) = a_trans(1,1); adiag_4(2) = a_trans(2,2); adiag_4(3) = a_trans(3,3); adiag_4(4) = 1d0
          quat_squared = matmul( Mat_trans , adiag_4)



          ! if any quaternion squared values are less than zero this is numerical issue, set to zero
          do i=1, size(quat_squared)
             if ( quat_squared(i) < 0d0 ) then
                if ( quat_squared(i) < small_neg ) then
                   ! this is too small!
                   stop "squared quaternion value is negative, more than numerical error"
                else
                quat_squared(i) = 0d0
                endif
             endif
             quaternion(i) = sqrt(quat_squared(i) )
          enddo

          ! now get signs
          ! since q and -q represent the same rotation, we can arbitrarily choose one of these elements to
          ! be positive.

          ! here, we have 4 options, to avoid rare cases where certain quaternion elements are zero
          q0=quaternion(1) ; q1=quaternion(2) ; q2=quaternion(3) ; q3=quaternion(4) ;

          if ( ( q0 .ge. q1 ) .and. ( q0 .ge. q2 ) .and. ( q0 .ge. q3 ) ) then
             ! if q0 is non-zero, set it positive
             quaternion(2) = sign( quaternion(2) , ( a_trans(3,2) - a_trans(2,3) ) )
             quaternion(3) = sign( quaternion(3) , ( a_trans(1,3) - a_trans(3,1) ) )
             quaternion(4) = sign( quaternion(4) , ( a_trans(2,1) - a_trans(1,2) ) )
          else if ( ( q1 .ge. q0 ) .and. ( q1 .ge. q2 ) .and. ( q1 .ge. q3 ) ) then
             quaternion(1) = sign( quaternion(1) , ( a_trans(3,2) - a_trans(2,3) ) )  
             quaternion(3) = sign( quaternion(3) , ( a_trans(2,1) + a_trans(1,2) ) )    
             quaternion(4) = sign( quaternion(4) , ( a_trans(1,3) + a_trans(3,1) ) )   
          else if ( ( q2 .ge. q0 ) .and. ( q2 .ge. q1 ) .and. ( q2 .ge. q3 ) ) then
             quaternion(1) = sign( quaternion(1) , ( a_trans(1,3) - a_trans(3,1) ) )
             quaternion(2) = sign( quaternion(2) , ( a_trans(2,1) + a_trans(1,2) ) )
             quaternion(4) = sign( quaternion(4) , ( a_trans(3,2) + a_trans(2,3) ) )
          else if ( ( q3 .ge. q0 ) .and. ( q3 .ge. q1 ) .and. ( q3 .ge. q2 ) ) then
             quaternion(1) = sign( quaternion(1) , ( a_trans(2,1) - a_trans(1,2) ) )
             quaternion(2) = sign( quaternion(2) , ( a_trans(1,3) + a_trans(3,1) ) )
             quaternion(3) = sign( quaternion(3) , ( a_trans(3,2) + a_trans(2,3) ) )
          else
             stop " something wrong in generating quaternion signs in gen_inertia routine"
          endif


       End Select   ! end select on linear molecule test

       ! inertia tensor is zero for single atom
    else
       inertia(:)=0d0
    endif              ! end if on number of atoms greater than 1


  end subroutine gen_inertia


  !**********************************
  ! this subroutine orders the principle moments of inertia
  ! and corresponding eigenvectors from smallest to largest
  !
  ! this corresponds to a 90 degree rotation about one or more of the
  ! principle axis, so one negative sign is introduced per rotation to the eigenvectors
  !**********************************
  subroutine order_moments_inertia(p_inertia,a_trans)
    real*8,dimension(:),intent(inout) :: p_inertia
    real*8,dimension(:,:),intent(inout) :: a_trans

    real*8,dimension(size(p_inertia)) :: p_inertia_temp
    real*8,dimension(size(a_trans(:,1)),size(a_trans(1,:))) :: a_trans_temp
    integer  :: flag, i , j, index1,index2,index3

    p_inertia_temp = p_inertia
    a_trans_temp = a_trans

    ! find smallest eigenvalue
    do i=1,3
       flag=1
       do j=1,3
          if ( i .ne. j) then
             if ( p_inertia(i) .ge. p_inertia(j) ) then
                flag=0
             endif
          endif
       enddo
       if ( flag == 1 ) then
          index1 = i
          exit
       endif
    enddo

    ! find next smallest eigenvalue
    do i=1,3
       if ( i .ne. index1 ) then
          flag=1
          do j=1,3
             if ( (i .ne. j) .and. ( j .ne. index1 ) ) then
                if ( p_inertia(i) .ge. p_inertia(j) ) then
                   flag=0
                endif
             endif
          enddo
          if ( flag == 1 ) then
             index2 = i
             exit
          endif
       endif
    enddo

    ! and last
    do i=1,3
       if ( ( i .ne. index1 ) .and. ( i.ne. index2 ) ) then
          index3 = i
          exit
       endif
    enddo

! now we have to rotate inertia tensor

    p_inertia(1) = p_inertia_temp(index1)
    p_inertia(2) = p_inertia_temp(index2)
    p_inertia(3) = p_inertia_temp(index3)

   ! first figure out which eigenvalue is going in first index,
   ! and rotate eigenvectors accordingly if necessary

    if ( index1 .ne. 1 ) then
       ! need to switch first eigenvalue with either 2nd or 3rd
       if ( index2 .eq. 1 ) then
          ! switch 1st and 2nd eigenvalue
          a_trans(:,index1) = a_trans_temp(:,index2)
          a_trans(:,index2) = -a_trans_temp(:,index1)
          index2=index1
          index1=1
          a_trans_temp=a_trans
       elseif( index3 .eq. 1 ) then
          ! switch 1st and 2nd eigenvalue
          a_trans(:,index1) = a_trans_temp(:,index3)
          a_trans(:,index3) = -a_trans_temp(:,index1)
          index3=index1
          index1=1
          a_trans_temp=a_trans
       endif
    endif

    ! now we have to consider 2nd and third columns

    if ( index2 .ne. 2 ) then
       ! need to switch 2nd and third
       a_trans(:,index2) = a_trans_temp(:,index3)
       a_trans(:,index3) = -a_trans_temp(:,index2)
       index2=index3
       index3=3
    endif


  end subroutine order_moments_inertia



  !*************************************
  ! this subroutine generates a random molecular orientation for linear and 
  ! nonlinear molecules
  !
  ! note that for a uniformly distributed random rotation matrix, one cannot just randomly 
  ! choose euler angles, this is incorrect
  ! 
  ! instead choose quaternions using Shoemake's algorithm
  !*************************************
  subroutine gen_random_orientation(new_xyz,new_r_com,new_n_mole,ins_type,new_quaternion)
    use global_variables
    integer, intent(in)  :: new_n_mole, ins_type
    real*8,dimension(:,:),intent(in)::new_r_com
    real*8, dimension(:,:,:),intent(inout) :: new_xyz  
    real*8, dimension(:,:), intent(inout) :: new_quaternion

    integer:: i_atom, n_atom_ins,i,j
    real*8 :: a(3,3),a_trans(3,3), r_vec(3)
    real*8 :: R1,R2,x0,y1,y2
    real*8 :: q0,q1,q2,q3

    ! figure out how many atoms there are in the inserted molecule
    do i_atom =1, MAX_N_ATOM
       if( molecule_type(ins_type,i_atom) .eq. MAX_N_ATOM + 1) exit
       n_atom_ins = i_atom
    enddo


    ! remember, r_mole_com is stored in principle axis coordinates for non-linear molecules
    ! so consistent rotation with these new quaternion values is to rotate r_mole_com vector
    ! back to space coordinates using a_transpose.
    ! for linear molecules, r_mole_com is in arbitrary frame, but this is okay since we
    ! aren't using quaternions for linear molecules and rotation is still random

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

    new_quaternion(new_n_mole,1)=q0 ;
    new_quaternion(new_n_mole,2)=q1 ;
    new_quaternion(new_n_mole,3)=q2 ;
    new_quaternion(new_n_mole,4)=q3 ;

    call quaternion_transformation_matrix(a,new_quaternion(new_n_mole,:) )


    do i=1,3
       do j=1,3
          a_trans(i,j) = a(j,i)
       enddo
    enddo

    do i_atom=1, n_atom_ins
       r_vec(:) = matmul( a_trans ,r_mole_com(ins_type,i_atom,:) )
       new_xyz(new_n_mole,i_atom,:) = r_vec(:) + new_r_com(new_n_mole,:)
    enddo


  end subroutine gen_random_orientation


  !************************************
  ! this subroutine creates a transformation matrix from the space fixed coordinate system
  ! to the body centered coordinate system given by the euler angles, phi, psi, theta
  ! input are phi, psi, and z, where z is cos(theta)
  ! it also generates the quaternions corresponding to these euler angles
  !***********************************
  subroutine euler_transformation_matrix(a,phi,psi,theta,quaternion)
    real*8,dimension(3,3),intent(out) :: a
    real*8,dimension(:), intent(inout) :: quaternion
    real*8,intent(in) :: phi,psi,theta
    real*8 :: cphi,sphi,cpsi,spsi,sthe,z
    real*8 :: cthe2,sthe2

    cphi = dcos(phi)
    sphi = dsin(phi)
    cpsi = dcos(psi)
    spsi = dsin(psi)
    z    = dcos(theta)
    sthe = dsqrt(1.d0-z**2)

    ! first generate the transformation matrix
    a(1,1) =  cpsi*cphi-z*sphi*spsi
    a(2,1) = -spsi*cphi-z*sphi*cpsi
    a(3,1) =  sthe*sphi
    a(1,2) =  cpsi*sphi+z*cphi*spsi
    a(2,2) = -spsi*sphi+z*cphi*cpsi
    a(3,2) = -sthe*cphi
    a(1,3) = sthe*spsi
    a(2,3) = sthe*cpsi
    a(3,3) = z

    ! now generate the corresponding quaternions for these euler angles
    cthe2 = dcos(theta/2d0)
    sthe2 = dsin(theta/2d0)

    quaternion(1) = cthe2* dcos(.5d0*(phi+psi) )
    quaternion(2) = sthe2* dcos(.5d0*(phi-psi) )
    quaternion(3) = sthe2* dsin(.5d0*(phi-psi) )
    quaternion(4) = cthe2* dsin(.5d0*(phi+psi) )


  end subroutine euler_transformation_matrix


  !************************************
  ! this subroutine creates a transformation matrix from the space fixed coordinate system
  ! to the body centered coordinate system given by the quaternions, q0,q1,q2,q3
  !***********************************
  subroutine quaternion_transformation_matrix(a,quaternion)
    real*8,dimension(3,3),intent(out) :: a
    real*8,dimension(:), intent(in) :: quaternion

    real*8 :: q0,q1,q2,q3

    q0 = quaternion(1) ;  q1 = quaternion(2) ; q2 = quaternion(3) ; q3 = quaternion(4) ;

    a(1,1) = q0**2 + q1**2 - q2**2 - q3**2
    a(2,1) = 2d0 * (q1*q2 - q0*q3)
    a(3,1) = 2d0 * (q1*q3 + q0*q2)
    a(1,2) = 2d0 * (q1*q2 + q0*q3)
    a(2,2) = q0**2 - q1**2 + q2**2 - q3**2
    a(3,2) = 2d0 * (q2*q3 - q0*q1)
    a(1,3) = 2d0 * (q1*q3 - q0*q2)
    a(2,3) = 2d0 * (q2*q3 + q0*q1)
    a(3,3) = q0**2 - q1**2 - q2**2 + q3**2

  end subroutine quaternion_transformation_matrix




end module rigid_body
