!*********************************************************************
!  This module is filled with miscellaneous subroutines that are used
!  extensively throughout the program, but seemingly have no better
!  place to reside.
!*********************************************************************

module routines
  implicit none

contains

  !********************************************************
  ! this subroutine sorts and stores the input files given as
  ! command line arguments, based on what
  ! kind of simulation is being run
  !********************************************************
  subroutine sort_input_files( n_files, ifile_conf, ifile_traj, ifile_fconf, ifile_ffpmt, ifile_top, ifile_simpmt , ofile_traj, ofile_log, ofile_hop_local )
    use global_variables
    integer, intent(out) :: n_files
    character(*),intent(out) :: ifile_conf, ifile_traj, ifile_fconf, ifile_ffpmt, ifile_top, ifile_simpmt , ofile_traj, ofile_log, ofile_hop_local

    integer :: N_base_files

    ! may not be using topology or fconf file in certain cases.  Initialize these with empty strings
    ifile_fconf=""
    ifile_top=""
    ofile_hop_local=""

    Select Case(ms_evb_simulation)
    Case("yes")
       ! 7 base files : ifile_conf, ifile_ffpmt, ifile_top, ifile_simpmt, ofile_traj, ofile_log, ofile_hop
       ! or 8, if fitting
       Select Case(ms_evb_parameter_fit)
       Case("yes")
          N_base_files = 8
       Case("no")
          N_base_files = 7
       End Select
    Case("no")
       ! 6 base files : ifile_conf, ifile_ffpmt, ifile_top, ifile_simpmt, ofile_traj, ofile_log
       N_base_files = 6
    End Select

    n_files = IARGC ()

    if (n_files /= N_base_files ) then
       stop " can't recognize number of i/o files ! "
    end if

!!$    if(n_files == N_base_files) then
!!$       framework_simulation=0
!!$    elseif(n_files == ( N_base_files + 1 ) ) then
!!$       framework_simulation=1
!!$    else
!!$       write(*,*) ""
!!$       write(*,*) "number of input files is not recognized"
!!$       write(*,*) "for rigid-body simulation, need 5 input files:"
!!$       write(*,*) "ifile_conf, ifile_ffpmt, ifile_simpmt, ofile_traj, ofile_log"
!!$       write(*,*) "for flexible molecule simulation, need additional topology file"
!!$       write(*,*) "for framework simulation, need additional .fconf file"
!!$       stop 
!!$    endif

    ! assume no framework_simulation for now
    framework_simulation=0

!!$    if (framework_simulation .eq.0) then
!!$       call getarg( 1, ifile_conf )
!!$       call getarg( 2, ifile_ffpmt )
!!$       call getarg( 3, ifile_top )
!!$       call getarg( 4, ifile_simpmt )
!!$       call getarg( 5, ofile_traj )
!!$       call getarg( 6, ofile_log )
!!$    else
!!$       call getarg( 1, ifile_conf)
!!$       call getarg( 2, ifile_fconf)
!!$       call getarg( 3, ifile_ffpmt )
!!$       call getarg( 4, ifile_top )
!!$       call getarg( 5, ifile_simpmt )
!!$       call getarg( 6, ofile_traj )
!!$       call getarg( 7, ofile_log )
!!$    endif

    Select Case(ms_evb_simulation)
    Case("yes")
    if ( ms_evb_parameter_fit .eq. "yes" ) then
       call getarg( 1, ifile_conf )
       call getarg( 2, ifile_traj )
       call getarg( 3, ifile_ffpmt )
       call getarg( 4, ifile_top )
       call getarg( 5, ifile_simpmt )
       call getarg( 6, ofile_traj )
       call getarg( 7, ofile_log )
       call getarg( 8, ofile_hop_local )
    else
       call getarg( 1, ifile_conf )
       call getarg( 2, ifile_ffpmt )
       call getarg( 3, ifile_top )
       call getarg( 4, ifile_simpmt )
       call getarg( 5, ofile_traj )
       call getarg( 6, ofile_log )
       call getarg( 7, ofile_hop_local )
    endif
    Case("no")
       call getarg( 1, ifile_conf )
       call getarg( 2, ifile_ffpmt )
       call getarg( 3, ifile_top )
       call getarg( 4, ifile_simpmt )
       call getarg( 5, ofile_traj )
       call getarg( 6, ofile_log )
    End Select

  end subroutine sort_input_files

  

  !**********************************************************************
  ! This subroutine checks if we are restarting a trajectory
  ! this is determined by the existence of both the velocity checkpoint file
  ! and the two coordinate and energy output files.  Also the last
  ! time step of the velocity checkpoint file has to match the last time step
  ! in the trajectory output file, otherwise we will print an error
  !**********************************************************************
  subroutine check_restart_trajectory( ofile_traj , ofile_log, traj_file, log_file, vel_file )
    use global_variables
    character(*),intent(in) :: ofile_traj , ofile_log
    ! file handles
    integer, intent(in) :: traj_file, log_file, vel_file
    logical  :: ofile_traj_exist, ofile_log_exist, velocity_file_exist

    integer :: flag, flag_eof, nargs, i_step_traj, i_step_vel
    character(400) :: line
    character(40),dimension(10)  :: args
    character(10) :: junk

    Inquire( file=ofile_traj, exist=ofile_traj_exist )
    Inquire( file=ofile_log , exist=ofile_log_exist  )    
    Inquire( file=velocity_file, exist=velocity_file_exist )


    if ( ofile_traj_exist .and. ofile_log_exist .and. velocity_file_exist ) then

       ! ***************** read ofile_traj ******************
       ! get last step from ofile_traj
       open( traj_file, file=ofile_traj, status='old' )
       ! read until last printed step
       do 
          call read_file_find_heading( traj_file, 'i_step' , flag, flag_eof )
          if ( flag_eof == -1 ) Exit  ! end of file
          backspace(traj_file)
          read(traj_file,'(A)') line
          call parse( line," ", args, nargs )
          if ( nargs /= 4 ) stop "error reading line arguments in trajectory file"
          read(args(3),'(I10)') i_step_traj
       enddo

       ! ****************** read velocity_file ***************
       ! get last step from ofile_traj
       open( vel_file, file=velocity_file, status='old' )
       ! read until last printed step
       do 
          call read_file_find_heading( vel_file, 'step' , flag, flag_eof )
          if ( flag_eof == -1 ) Exit  ! end of file
          backspace(vel_file)
          read(vel_file,*) junk, i_step_vel
       enddo


       !************** now make sure the last steps in the trajectory and velocity files are the same
       if ( ( i_step_vel == i_step_traj ) .and. ( i_step_vel > 0 ) ) then
          restart_trajectory="yes"
          n_old_trajectory= i_step_traj
       else
          write(*,*) ""
          write(*,*) " error restarting trajectory.  last step is not the same in the trajectory "
          write(*,*) " output and velocity checkpointing files "
          write(*,*) ""
          stop
       end if

       ! for the log file, go to end of file, to prepare to append output
       ! leave this file open
       open( log_file, file=ofile_log, status='old' )  
       do 
          ! read garbage
          call read_file_find_heading( log_file, 'asdfdasdfasdfasdf' , flag, flag_eof )
          if ( flag_eof == -1 ) Exit  ! end of file
       enddo

       ! close traj and vel files, these will be reopened, and set up to append to, when
       ! last coordinates and velocities are read
       close(vel_file)
       close(traj_file)

    else
       restart_trajectory="no"
       n_old_trajectory=0
    endif


  end subroutine check_restart_trajectory

  


  !**********************************************************************
  ! this subroutine reads coordinates from a previous trajectory, to restart
  ! a simulation
  !**********************************************************************
  subroutine read_coordinates_restart_trajectory( traj_file, ofile_traj, n_mole, n_atom, xyz, molecule_name, alist, box, n_old_trajectory )
   character(*),intent(in) :: ofile_traj
    ! file handle
    integer, intent(in) :: traj_file, n_old_trajectory    
    integer, intent(out) :: n_mole
    integer, dimension(:), intent(out) :: n_atom
    real*8, dimension(:,:,:),intent(out) :: xyz
    character(*),dimension(:),intent(out) :: molecule_name
    character(*),dimension(:,:),intent(out) :: alist
    real*8,dimension(:,:), intent(out) :: box

    integer :: flag, flag_eof, nargs, i_step_traj, atoms_local, i_mole, i_atom, i_mole_local, i_atom_local, i
    character(len(molecule_name(1))) :: mname
    character(len(alist(1,1))) :: aname
    character(400) :: line
    character(40),dimension(10)  :: args
    real*8, dimension(3) :: box_temp
    real*8, parameter :: small =1D-6
    
       open( traj_file, file=ofile_traj, status='old' )
       ! read until the final step
       do 
          call read_file_find_heading( traj_file, 'i_step' , flag, flag_eof )
          if ( flag_eof == -1 ) Exit  ! end of file
          backspace(traj_file)
          read(traj_file,'(A)') line
          call parse( line," ", args, nargs )
          if ( nargs /= 4 ) stop "error reading line arguments in trajectory file"
          read(args(3),'(I10)') i_step_traj

          ! if final step, exit this loop
          if ( i_step_traj == n_old_trajectory ) Exit

       enddo

       ! now we should be to the coordinates of the last step, store this trajectory
       call read_conf( traj_file, n_mole, n_atom, xyz, molecule_name,alist, box )

      ! leave file open.  We should be at the end of the trajectory file, and so we are ready to append
      ! the new trajectory output

  end subroutine read_coordinates_restart_trajectory




  !**********************************************************************
  ! this subroutine reads velocities from the velocity checkpoint file to restart
  ! a simulation
  !**********************************************************************
  subroutine read_velocity_restart_checkpoint( vel_file, velocity_file, n_mole, n_atom, vel_atom, n_old_trajectory )
   character(*),intent(in) :: velocity_file
    ! file handle
    integer, intent(in) :: vel_file    
    integer, intent(in) :: n_mole, n_old_trajectory
    integer, dimension(:), intent(in) :: n_atom
    real*8, dimension(:,:,:),intent(out) :: vel_atom

    integer :: flag, flag_eof, nargs, i_step_traj, atoms_local, i_mole, i_atom, i_mole_local, i_atom_local, i
    character(5) :: mname
    character(5) :: aname
    character(400) :: line
    character(10) :: junk
    character(40),dimension(10)  :: args

    
       open( vel_file, file=velocity_file, status='old' )
       ! read until the final step
       do 
          call read_file_find_heading( vel_file, 'step' , flag, flag_eof )
          if ( flag_eof == -1 ) Exit  ! end of file
          backspace(vel_file)
          read(vel_file,*) junk , i_step_traj
          ! if final step, exit this loop
          if ( i_step_traj == n_old_trajectory ) Exit
       enddo

       ! now we should be to the velocities of the last step
       do i_mole=1, n_mole
          do i_atom=1, n_atom(i_mole)
          read( vel_file, '(I5,2A5,I5,3F14.6)' ), i_mole_local, mname, aname, i_atom_local, vel_atom(i_mole,i_atom,:)
          enddo
       enddo

      ! leave file open.  We should be at the end of the trajectory file, and so we are ready to append
      ! the new trajectory output

  end subroutine read_velocity_restart_checkpoint




  !***********************************************************************
  ! This subroutine read input files ( *.conf file )
  ! It also assigns charges to each atoms, based on parameters in glob_v.f90 file
  !***********************************************************************
  subroutine read_conf( file_handle, n_mole, n_atom, xyz, molecule_name, alist, box )
    use global_variables
    implicit none
    integer, intent(in) :: file_handle
    character(MAX_MNAME), intent(out),dimension(MAX_N_MOLE) :: molecule_name
    character(MAX_ANAME), intent(out), dimension(MAX_N_MOLE,MAX_N_ATOM) :: alist
    integer, intent(out) :: n_mole
    integer, intent(out), dimension(MAX_N_MOLE) :: n_atom
    real*8, intent(out), dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: xyz
    real*8, intent(out), dimension(3,3) :: box
    integer :: ifile, i,j, n_line, i_mole, i_atom, i_mole_prev, junk, inputstatus, ind, flag, i_type, i_param, n_mole_atoms
    real*8, dimension(3) :: r_tmp
    character(MAX_MNAME) :: mname
    character(MAX_ANAME) :: aname
    character(50)::line
    real*8 :: tmp1, tmp2, tmp3
    integer :: ios

    read( file_handle, '(I)' ) n_line
    ! total atoms in simulation
    total_atoms = n_line
    n_mole = 0
    i_atom = 0
    if (n_line .gt. 0 ) then
       do i = 1, n_line
          read( file_handle, '(I5,2A5,I5,3F14.6)' ), i_mole, mname, aname, junk, r_tmp(1), r_tmp(2), r_tmp(3)
          call trim_end( aname )
          if ( i_mole /= n_mole ) then
             if ( n_mole > 0 ) then
                n_atom(n_mole) = i_atom
             end if
             n_mole = i_mole
             i_atom = 0
             call trim_end( mname )
             molecule_name(i_mole) = mname
          end if
          i_atom = i_atom + 1

          ! make sure we dont have too many atoms for this molecule
          if ( i_atom > MAX_N_ATOM ) then
             stop "please increase the value of MAX_N_ATOM in glob_v.f90 to accomodate the size of the solute molecules"
          endif

          xyz( i_mole, i_atom, : ) = r_tmp(:)
          alist( i_mole, i_atom ) = aname
       end do
       n_atom(n_mole) = i_atom
    endif

    read(file_handle,*) tmp1, tmp2, tmp3
    read(file_handle, *, iostat=ios ) box(2,1), box(2,2), box(2,3)
    if ( ios < 0 ) then ! end of the file
       box = 0.0d0
       box(1,1) = tmp1
       box(2,2) = tmp2
       box(3,3) = tmp3
    else
       box(1,1) = tmp1
       box(1,2) = tmp2
       box(1,3) = tmp3
       read( file_handle, * ) box(3,1), box(3,2), box(3,3)
    endif


    ! make sure we don't have too many molecules
    if ( n_mole > MAX_N_MOLE ) then
       stop "please increase the value of MAX_N_MOLE in glob_v.f90 to accomodate the number of molecules in the system"
    endif


  end subroutine read_conf



  !******************************
  ! this subroutine returns the next non-comment
  ! line of a topology file.
  !
  ! lines that begin with ";" are taken as comments
  ! (gromacs convention)
  !
  ! settings of output flag:
  ! 0 : normal line
  ! 1 : blank line
  ! -1 : end of file
  !******************************
  subroutine read_topology_line( file_h , line , flag )
    integer, intent(in) :: file_h
    character(*), intent(out) :: line
    integer, intent(out) :: flag

    ! use gromacs comment convention
    character(5) :: comment=";"
    integer :: nargs, inputstatus, ind
    integer,parameter :: max_param=30
    character(80),dimension(max_param)::args

    flag=0
    do
       Read(file_h,'(A)',Iostat=inputstatus) line
       If(inputstatus < 0) then
          ! end of file
          flag=-1
          Exit
       end if

       ! check end of section
       call parse(line," ",args,nargs)
       if ( nargs == 0 ) then
          flag=1
          Exit
       end if
       ! check comment
       ind=INDEX(args(1),comment)
       ! if not a comment line, then return this string
       if ( ind == 0 ) exit
    end do

  end subroutine read_topology_line



  !*************************************
  ! This subroutine reads a file until specified heading is found
  ! if heading is found, flag_found is set to "1", else it is "0"
  ! if end of file is reached, flag_eof is set to "-1"
  !*************************************
  subroutine read_file_find_heading( file_h, heading, flag_found, flag_eof )
    integer, intent(in) :: file_h
    character(*), intent(in) :: heading
    integer, intent(out) :: flag_found, flag_eof

    integer  :: ind, flag
    character(300) :: input_string

    flag_found=0
    flag_eof=0

    do
       call read_topology_line( file_h , input_string , flag )
       ! if end of file
       if ( flag == -1 ) then
          flag_eof = -1   
          exit
       end if

       ind=INDEX(input_string,heading)
       if ( ind .ne. 0 ) then
          flag_found=1
          exit
       end if
    enddo

  end subroutine read_file_find_heading




  !********************************************************
  ! This function find the center of mass of a molecule
  !********************************************************
  function pos_com( xyz, i_mole, n_atom, mass )
    implicit none
    real*8, intent(in), dimension(:,:,:) :: xyz
    real*8, intent(in), dimension(:,:) :: mass
    integer, intent(in) :: i_mole
    integer, intent(in), dimension(:) :: n_atom
    integer :: i_atom
    real*8, dimension(3) :: pos_com
    real*8 :: m_tot
    pos_com(:) = 0D0
    m_tot = 0D0

    do i_atom = 1, n_atom(i_mole)
       pos_com(:) = pos_com(:) + xyz(i_mole,i_atom,:) * mass(i_mole,i_atom)
       m_tot = m_tot + mass(i_mole,i_atom)
    end do
    pos_com(:) = pos_com(:) / m_tot

  end function pos_com


  !*********************************************************
  ! This function update all COM positions
  !*********************************************************
  subroutine update_r_com( xyz, n_mole, n_atom, mass, r_com )
    implicit none
    real*8, intent(in), dimension(:,:,:) :: xyz
    integer, intent(in) :: n_mole
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:) :: mass
    real*8, intent(out), dimension(:,:) :: r_com
    integer :: i_mole
    do i_mole = 1, n_mole
       r_com( i_mole, : ) = pos_com( xyz, i_mole, n_atom, mass )
    end do
  end subroutine update_r_com



  !**************************************************
  ! this subroutine creates the transformation matrix
  ! from cartesian to box coordinates
  !*************************************************
  subroutine initialize_non_orth_transform ( box )
    use global_variables
    real*8, dimension(3,3),intent(in) :: box

    real*8, dimension(3,3) :: temp_xyz
    real*8, dimension(3,3) :: junk
    integer :: i,j

    junk=0d0 
    junk(1,1)=1d0;junk(2,2)=1d0;junk(3,3)=1d0
    ! inverse of box vectors
    do i=1,3
       do j=1,3
          temp_xyz(i,j) = box(j,i)
       enddo
    enddo

    call gaussj(temp_xyz, junk)

    xyz_to_box_transform = temp_xyz

  end subroutine initialize_non_orth_transform


  !******************************************
  ! reciprocal lattice vector.  This is essentially the same
  ! as initialize_non_orth_transform subroutine, but we keep both
  ! in for compatibility with older code
  !******************************************
  subroutine construct_reciprocal_lattice_vector(kk,box)
    real*8,dimension(:,:),intent(out) :: kk
    real*8,dimension(:,:),intent(in) :: box

    real*8 :: a(3), b(3), c(3), ka(3), kb(3), kc(3), vol

    a(:) = box(1,:)
    b(:) = box(2,:)
    c(:) = box(3,:)

    ! calculate the volume and the reciprocal vectors (notice no 2pi)
    vol = volume( a, b, c )
    call crossproduct( a, b, kc ); kc = kc /vol 
    call crossproduct( b, c, ka ); ka = ka /vol
    call crossproduct( c, a, kb ); kb = kb /vol
    kk(1,:)=ka(:);kk(2,:)=kb(:);kk(3,:)=kc(:)

  end subroutine construct_reciprocal_lattice_vector



  !********************************************
  ! this subroutine creates direct coordinates, scaled
  ! by input integer "K" (pme_grid), using the
  ! reciprocal lattice vectors
  !********************************************
  subroutine create_scaled_direct_coordinates(xyz_scale, xyz, n_mole, n_atom, kk, K)
    real*8,dimension(:,:,:),intent(out) :: xyz_scale
    real*8,dimension(:,:,:),intent(in) :: xyz
    integer, intent(in) :: n_mole
    integer, dimension(:),intent(in) :: n_atom
    real*8,dimension(:,:),intent(in) :: kk
    integer, intent(in) :: K

    integer :: i_mole
    real*8,dimension(:,:),allocatable :: xyz_local

    allocate( xyz_local(size(xyz(1,:,1)),size(xyz(1,1,:))) )

    do i_mole=1,n_mole
       xyz_local(:,:) = xyz(i_mole,:,:)
       call create_scaled_direct_coordinates_molecule(xyz_scale(i_mole,:,:), xyz_local, n_atom(i_mole), kk, K)
    enddo

  end subroutine create_scaled_direct_coordinates


  !********************************************
  ! this subroutine creates direct coordinates for a particular molecule, scaled
  ! by input integer "K" (pme_grid), using the
  ! reciprocal lattice vectors
  !********************************************
  subroutine create_scaled_direct_coordinates_molecule(xyz_scale, xyz, n_atom, kk, K)
    real*8,dimension(:,:),intent(inout) :: xyz_scale
    real*8,dimension(:,:),intent(inout) :: xyz
    integer, intent(in) :: n_atom
    real*8,dimension(:,:),intent(in) :: kk
    integer, intent(in) :: K

    integer :: j,l
    real*8,parameter :: small=1D-6

    do j=1,n_atom
       do l=1,3
          xyz_scale(j,l)=dble(K)*dot_product(kk(l,:),xyz(j,:))
          ! if atoms are out of grid, shift them back in
          if (xyz_scale(j,l)<0.) then
             xyz_scale(j,l)=xyz_scale(j,l)+dble(K)
             ! else if(xyz_scale(j,l)>dble(K)) then ! changed by Jesse, 6/10/2014 see below line.  In rare cases, this could lead to out of bounds array calls when constructing Q
          else if(xyz_scale(j,l)>= dble(K)) then
             xyz_scale(j,l)=xyz_scale(j,l)-dble(K)
          endif
          ! make sure scaled coordinates are not numerically equal to zero, otherwise this will screw up Q grid routine
          if ( abs(xyz_scale(j,l)) < small ) then
             xyz_scale(j,l) = small
          end if
       enddo
    enddo

  end subroutine create_scaled_direct_coordinates_molecule



  !***************************************************************************
  ! This function calculate shift vector for PBC
  ! for a general treatment of non-orthogonal boxes
  !
  ! notice the floor argument for finding closest image is not valid in general
  ! for a non-orthogonal box for a cutoff of half the box length.
  ! however, there will always be a certain cutoff, smaller than half the box length
  ! for which it is valid.  Make sure this is being used
  !***************************************************************************
  function pbc_shift( r_com_i, r_com_j, box , xyz_to_box_transform )
    real*8, dimension(3) :: pbc_shift
    real*8, dimension(3),intent(in) :: r_com_i, r_com_j
    real*8, dimension(3,3),intent(in) :: box , xyz_to_box_transform
    real*8, dimension(3) :: dr_com, dr_box , shift
    integer :: i,j


    dr_com(:) = r_com_j(:) - r_com_i(:)

    dr_box = matmul ( xyz_to_box_transform , dr_com )

    do i = 1, 3
       shift(i) = floor ( dr_box(i) + 0.5d0 )
    enddo

    pbc_shift(:) = 0.0
    do i = 1, 3
       do j = 1 , 3
          pbc_shift(i) = pbc_shift(i) + shift(j)*box(j,i)
       enddo
    enddo

  end function pbc_shift



  !******************************************************************
  ! This function calculate dr between two atoms, 
  ! @ith shift vector calcuated by pbc_shift
  !******************************************************************
  function pbc_dr( ri, rj, shift )
    real*8, dimension(3) :: pbc_dr
    real*8, intent(in), dimension(3) :: shift, ri, rj
    pbc_dr(:) = rj(:) - ri(:) - shift(:)
  end function pbc_dr



  !**************************************************************************
  ! this function searches the error function grid, and picks out best value
  ! this assumes all arguments are positive
  !**************************************************************************
  function erfc_g(x)
    use global_variables
    real*8::erfc_g
    real*8,intent(in)::x

    Select case (grid_erfc)
    case("yes")
       if(x > erfc_max) then
          erfc_g=0D0
          ! careful about calling erfc_table(0), rather call erfc_table(1)
       elseif ( x > 0d0 ) then
          erfc_g=erfc_table(ceiling(x/erfc_max*dble(erfc_grid)))
       else
          !here x=0
          erfc_g = 1d0
       endif
    case("no")
       erfc_g = erfc(x)
    case default
       stop "grid_erfc option not recognized"
    end select

  end function erfc_g

  !**************************************************************************
  ! this function searches the B spline grids, and interpolates B_spline of a
  ! real number from the grids, the accuracy gained from using this subroutine is 
  ! probably not worth the extra time, but this subroutine is mainly for
  ! use when calculating forces numerically using pme energy routine
  !**************************************************************************
  function interp_B_spline(x,n)
    use global_variables
    real*8::interp_B_spline
    real*8,intent(in)::x
    integer,intent(in)::n
    real*8::grid_val
    integer::n_int


    n_int=floor(x*dble(spline_grid)/dble(n))
    grid_val=x*dble(spline_grid)/dble(n)
    if(n_int.eq.0) then
       if(n .eq. 6 ) then
          interp_B_spline= (grid_val-dble(n_int))* B6_spline(n_int+1)
       elseif(n .eq. 5) then
          interp_B_spline= (grid_val-dble(n_int))* B5_spline(n_int+1)
       elseif(n .eq. 4) then
          interp_B_spline= (grid_val-dble(n_int))* B4_spline(n_int+1)
       elseif(n .eq. 3) then
          interp_B_spline= (grid_val-dble(n_int))* B3_spline(n_int+1)
       endif
    else
       if(n .eq. 6 ) then
          interp_B_spline= (1. -(grid_val-dble(n_int)))*B6_spline(n_int) + (grid_val-dble(n_int))* B6_spline(n_int+1)
       elseif(n.eq.5) then
          interp_B_spline= (1. -(grid_val-dble(n_int)))*B5_spline(n_int) + (grid_val-dble(n_int))* B5_spline(n_int+1)
       elseif(n.eq.4) then
          interp_B_spline= (1. -(grid_val-dble(n_int)))*B4_spline(n_int) + (grid_val-dble(n_int))* B4_spline(n_int+1)
       elseif(n.eq.3) then
          interp_B_spline= (1. -(grid_val-dble(n_int)))*B3_spline(n_int) + (grid_val-dble(n_int))* B3_spline(n_int+1)
       endif
    endif


  end function interp_B_spline



  !***********************************
  ! this subroutine is used to find the index
  ! of the heavy atom to which a hydrogen atom
  ! is bonded.
  !***********************************
  subroutine find_bonded_atom_hydrogen( i_mole_type, n_atom , hydrogen_atom , heavy_atom )
    use global_variables
    integer, intent(in) :: i_mole_type, hydrogen_atom
    integer, intent(in) :: n_atom
    integer, intent(out) :: heavy_atom

    integer :: j_atom, count

    count=0
    do j_atom=1, n_atom
       if ( molecule_bond_list(i_mole_type,hydrogen_atom,j_atom) == 1 ) then
          heavy_atom = j_atom
          count = count + 1
       end if
    enddo

    ! make sure we found 1 and only 1 atom
    if ( count /= 1 ) then
       stop "error in subroutine find_bonded_atom_hydrogen"
    end if

  end subroutine find_bonded_atom_hydrogen



  !*****************************************
  ! this subroutine is used to find the index corresponding
  ! to a set of atom types.  For example, this is used
  ! in ms-evb routines, in which there are specialized interactions
  ! that contribute to the matrix elements.
  ! In this case, it is efficient to just enumerate all of these interactions,
  ! since there are few of them, and have a corresponding index describing each
  ! interaction, from which the parameters can be looked up
  !******************************************
  subroutine get_index_atom_set( index, lookup_array , itype_array )
    integer, intent(out) :: index
    integer, dimension(:,:), intent(in) :: lookup_array
    integer, dimension(:), intent(in)   :: itype_array

    integer :: i , j, flag

    if (size(lookup_array(1,:)) /= size(itype_array) ) then
       stop "sizes of input arrays are not consistent in subroutine 'get_index_atom_set'"
    endif

    index=-1
    do i=1, size(lookup_array(:,1))
       if ( lookup_array(i,1) == 0 ) exit
       flag=1
       do j=1, size(lookup_array(1,:))
          if ( lookup_array(i,j) /= itype_array(j) ) flag=0
       enddo
       if ( flag == 1 ) then
          index=i
          exit
       end if
    enddo

  end subroutine get_index_atom_set




  !****************************
  ! this subroutine finds the first matching index
  ! of an atomtype for a particular molecule
  !****************************
  subroutine match_first_atom_type( index, i_mole, n_atom, atom_type, atom_index_local )
    integer, intent(out) :: index
    integer, intent(in) :: i_mole
    integer, dimension(:), intent(in) :: n_atom
    character(*),intent(in) :: atom_type
    integer, dimension(:,:), intent(in) :: atom_index_local

    integer :: i_atom, i_type

    index=index_first_atom_type( i_mole, n_atom, atom_type, atom_index_local )

    if ( index < 0 ) then
       write(*,*) ""
       write(*,*) "Can't find atomtype ", atom_type
       write(*,*) "molecule with index ", i_mole
       write(*,*) ""
       stop
    end if

  end subroutine match_first_atom_type



  function index_first_atom_type( i_mole, n_atom, atom_type, atom_index_local )
    use global_variables
    integer             :: index_first_atom_type
    integer, intent(in) :: i_mole
    integer, dimension(:), intent(in) :: n_atom
    character(*),intent(in) :: atom_type
    integer, dimension(:,:), intent(in) :: atom_index_local

    integer :: i_atom, i_type

    index_first_atom_type=-1

    do i_atom=1, n_atom(i_mole)
       i_type = atom_index_local(i_mole, i_atom )

       if( atype_name( i_type ) .eq. atom_type ) then
          index_first_atom_type = i_atom
          exit
       end if
    enddo

  end function index_first_atom_type



  !*********************************************************
  ! this subroutine maps data structures from molecule storage
  ! ( i_mole, i_atom ), to absolute atom storage (i_atom ) 
  !
  ! note that we assume that we are going to loop over the atom_array
  ! with the atom index in the outer loop, and so we put this dimension
  ! last to optimize the memory retrival from column-major fortran convention
  !*********************************************************
  subroutine map_molecule_atom_data_structures_3d_m_to_a( n_mole, n_atom , molecule_map, atom_array , molecule_array )
    integer, intent(in) :: n_mole
    integer, dimension(:), intent(in) :: n_atom
    integer, dimension(:,:),intent(in) :: molecule_map
    real*8, dimension(:,:),intent(out) :: atom_array
    real*8, dimension(:,:,:),intent(in) :: molecule_array

    integer :: i_mole, i_atom, a_index

    ! atom array should be in reverse order as molecule array.  We assume we will be looping over atom array
    ! in inner loops, so we want first index to correspond to inner loop
    if ( size(atom_array(:,1)) /= size(molecule_array(1,1,:)) ) then
       stop " array sizes not consistent in map_molecule_atom_data_structures"
    endif

    do i_mole = 1, n_mole
       do i_atom = 1, n_atom(i_mole)
          ! atom index for this atom
          a_index = molecule_map(i_mole,i_atom)
          if (  a_index > 0  ) then
             atom_array(:,a_index) = molecule_array(i_mole,i_atom,:)
          end if
       enddo
    enddo

  end subroutine map_molecule_atom_data_structures_3d_m_to_a



  !*********************************************************
  ! this subroutine maps data structures from absolute atom storage (i_atom )
  ! to molecule storage ( i_mole, i_atom )
  !
  ! note that we assume that we are going to loop over the atom_array
  ! with the atom index in the outer loop, and so we put this dimension
  ! last to optimize the memory retrival from column-major fortran convention
  !*********************************************************
  subroutine map_molecule_atom_data_structures_3d_a_to_m( n_mole, n_atom , molecule_map, atom_array , molecule_array,  add )
    integer, intent(in) :: n_mole
    integer, dimension(:), intent(in) :: n_atom
    integer, dimension(:,:),intent(in) :: molecule_map
    real*8, dimension(:,:),intent(in) :: atom_array
    real*8, dimension(:,:,:),intent(inout) :: molecule_array
    character(3),intent(in), optional :: add

    integer :: i_mole, i_atom, a_index

    ! atom array should be in reverse order as molecule array.  We assume we will be looping over atom array
    ! in inner loops, so we want first index to correspond to inner loop
    if ( size(atom_array(:,1)) /= size(molecule_array(1,1,:)) ) then
       stop " array sizes not consistent in map_molecule_atom_data_structures"
    endif

    do i_mole = 1, n_mole
       do i_atom = 1, n_atom(i_mole)
          ! atom index for this atom
          a_index = molecule_map(i_mole,i_atom)

          if (  a_index > 0  ) then
             if ( present(add) ) then
                molecule_array(i_mole,i_atom,:) = molecule_array(i_mole,i_atom,:) + atom_array(:,a_index)
             else
                molecule_array(i_mole,i_atom,:) = atom_array(:,a_index)
             endif
          end if
       enddo
    enddo

  end subroutine map_molecule_atom_data_structures_3d_a_to_m



  !*********************************************************
  ! this subroutine maps data structures from molecule storage
  ! ( i_mole, i_atom ), to absolute atom storage (i_atom ) 
  !*********************************************************
  subroutine map_molecule_atom_data_structures_2d_m_to_a( n_mole, n_atom , molecule_map, atom_array , molecule_array )
    integer, intent(in) :: n_mole
    integer, dimension(:), intent(in) :: n_atom
    integer, dimension(:,:),intent(in) :: molecule_map
    real*8, dimension(:),intent(out) :: atom_array
    real*8, dimension(:,:),intent(in) :: molecule_array

    integer :: i_mole, i_atom, a_index

    ! atom array should be in reverse order as molecule array.  We assume we will be looping over atom array
    ! in inner loops, so we want first index to correspond to inner loop

    do i_mole = 1, n_mole
       do i_atom = 1, n_atom(i_mole)
          ! atom index for this atom
          a_index = molecule_map(i_mole,i_atom)
          if (  a_index > 0  ) then
             ! store atom information from molecular data arrays
             atom_array(a_index) = molecule_array(i_mole,i_atom)
          end if
       enddo
    enddo

  end subroutine map_molecule_atom_data_structures_2d_m_to_a


  !*********************************************************
  ! this is same as above subroutine , except works on integer arrays
  ! this subroutine maps data structures from molecule storage
  ! ( i_mole, i_atom ), to absolute atom storage (i_atom ) 
  !*********************************************************
  subroutine map_molecule_atom_data_structures_2d_m_to_a_int( n_mole, n_atom , molecule_map, atom_array , molecule_array )
    integer, intent(in) :: n_mole
    integer, dimension(:), intent(in) :: n_atom
    integer, dimension(:,:),intent(in) :: molecule_map
    integer, dimension(:),intent(out) :: atom_array
    integer, dimension(:,:),intent(in) :: molecule_array

    integer :: i_mole, i_atom, a_index

    ! atom array should be in reverse order as molecule array.  We assume we will be looping over atom array
    ! in inner loops, so we want first index to correspond to inner loop

    do i_mole = 1, n_mole
       do i_atom = 1, n_atom(i_mole)
          ! atom index for this atom
          a_index = molecule_map(i_mole,i_atom)
          if (  a_index > 0  ) then
             ! store atom information from molecular data arrays
             atom_array(a_index) = molecule_array(i_mole,i_atom)
          end if
       enddo
    enddo

  end subroutine map_molecule_atom_data_structures_2d_m_to_a_int




  !*****************************************************************
  ! this subroutine generates the reduced mass for an atom pair. 
  ! if considering a dummy atom with zero mass, at the midpoint of a molecule,
  ! use the total molecule mass for this atom, as this is what we want for the
  ! Feynmann-Hibbs correction
  !*****************************************************************
  subroutine generate_reduced_mass(mass_reduced,atom_id1,atom_id2,i_mole,j_mole,n_atom)
    use global_variables
    real*8,intent(out) :: mass_reduced
    integer,intent(in) :: atom_id1, atom_id2,i_mole,j_mole
    integer,dimension(:),intent(in) :: n_atom

    integer :: i_atom, j_atom, local_id1, local_id2
    real*8:: mass1,mass2

    ! if neither atom is a dummy atom, use those masses
    mass1=atype_mass(atom_id1)
    mass2=atype_mass(atom_id2)

    ! for dummy atoms, add up molecular mass
    if ( mass1 < 0.01 ) then
       mass1=0d0
       do i_atom=1,n_atom(i_mole)
          local_id1= atom_index(i_mole, i_atom)
          mass1=mass1 + atype_mass(local_id1)
       enddo
    endif

    if ( mass2 < 0.01 ) then
       mass2=0d0
       do j_atom=1,n_atom(j_mole)
          local_id2= atom_index(j_mole, j_atom)
          mass2=mass2 + atype_mass(local_id2)
       enddo
    endif


    ! make sure we now have non-zero masses
    if ( ( mass1 < 0.001 ) .or. ( mass2 < 0.001 ) ) then
       stop "mass of atom is equal to zero in generate_reduced_mass subroutine"
    endif

    mass_reduced = mass1*mass2 / ( mass1 + mass2 )

  end subroutine generate_reduced_mass



  !************************************************************************
  ! this subroutine seeds the fortran random number generator based on the system clock
  ! based on the gnu documentation
  !************************************************************************
  subroutine initialize_random_seed
    use global_variables
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed

    call random_seed(size = n)
    allocate(seed(n))
    call system_clock(count=clock)
    seed = clock + 37 * (/ (i-1,i=1,n) /)
    call random_seed(put = seed )

  end subroutine initialize_random_seed

  !*************************************************************************
  ! This subroutine generates a random unit vector
  ! from 'computer simulation of liquids', Allen,M.P.,Tildesley,D.J.
  !*************************************************************************
  subroutine random_unit_vector(unit)
    real*8,dimension(3),intent(out)::unit
    real*8::ransq,ran1,ran2,ranh
    real*8,dimension(3)::randnum3
    ransq=2.
    do while (ransq.ge.1)
       call random_number(randnum3)
       ran1=1.-2.*randnum3(1)
       ran2=1.-2.*randnum3(2)
       ransq=ran1*ran1+ran2*ran2
    enddo
    ranh=2.*sqrt(1.-ransq)
    unit(1)=ran1*ranh
    unit(2)=ran2*ranh
    unit(3)=(1.-2.*ransq)

  end subroutine random_unit_vector


  !*************************************************************************
  ! This subroutine moves all spaces at the beginning of a string to the end
  !*************************************************************************
  subroutine trim_head( aname )
    implicit none
    character(*), intent(inout) :: aname
    integer :: i, n, n_space, flag
    n = len( aname )
    n_space = 0
    flag = 0
    do i = 1, n
       if ( flag == 0 .and. aname(i:i) == ' ' ) then
          n_space = n_space + 1
       else if ( flag == 0 ) then
          flag = 1
          aname(i-n_space:i-n_space) = aname(i:i)
       else 
          aname(i-n_space:i-n_space) = aname(i:i)
       end if
    end do
    do i = n-n_space+1, n
       aname(i:i) = ' ' 
    end do
  end subroutine trim_head


  !*************************************************************************
  ! This subroutine moves all spaces at the beginning of a string to the end
  !*************************************************************************
  subroutine trim_end( aname )
    implicit none
    character(*), intent(inout) :: aname
    integer :: i, n, n_space, flag
    n = len( aname )
    n_space=0
    flag = 0
    do i = n, 1, -1
       if ( flag == 0 .and. aname(i:i) == ' ' ) then
          n_space = n_space + 1
       else if ( flag == 0 ) then
          flag = 1
          aname(i+n_space:i+n_space) = aname(i:i)
       else
          aname(i+n_space:i+n_space) = aname(i:i)
       end if
    end do
    do i = 1, n_space
       aname(i:i) = ' '
    end do
  end subroutine trim_end

  !*************************************************************************
  ! this subroutine prints information about the simulation
  !*************************************************************************
  subroutine print_run_param( log_file , n_mole, volume )
    use global_variables
    integer, intent(in) :: n_mole,log_file
    real*8, intent(in)  :: volume
    integer  :: i

    write( log_file, *) "*******run parameters***************************"
    write( log_file, *) "ensemble", "    ",select_ensemble
    write( log_file, *) "framework simulation", framework_simulation
    write( log_file, *) "lj_bkghm ", lj_bkghm
    write( log_file, *) "lj cutoff", lj_cutoff
    write( log_file, *) "verlet list ", use_verlet_list
    Select Case(use_verlet_list)
    Case("yes")
       write( log_file, *) "verlet_cutoff_elec ", verlet_cutoff_elec
       write( log_file, *) "verlet_cutoff_lj ", verlet_cutoff_lj
       write( log_file, *) "verlet_grid_based_construction ", verlet_grid_based_construction
    End Select
    write( log_file, *) "pme dispersion ", pme_disp
    write( log_file, *) "electrostatics","   ", electrostatic_type
    Select Case(electrostatic_type)
    Case("pme")
       write( log_file, *) "Beta Spline order ", spline_order
       write( log_file, *) "PME grid points ", pme_grid
       write( log_file, *) "alpha_sqrt (for pme) ", alpha_sqrt
    End Select
    write( log_file, *) "electrostatic cutoff", ewald_cutoff, Electro_cutoff
    write( log_file, *) "grid error function","   ", grid_erfc
    write( log_file, *) "electrostatic screening", screen_type
    write( log_file, *) "drude simulation", drude_simulation
    write( log_file, *) "distance for auto-rejection of move (A)", too_close
    write( log_file, *) "******special force field parameters************"
    write( log_file, *) "******parallel or serial ************"
    write( log_file, *) "number of threads  ", n_threads
    write( log_file, *) "******ensemble specific info********************"
    write( log_file, *) "temperature", temp
    write( log_file, *) "molecules", n_mole
    write( log_file, *) "volume", volume
    write( log_file, *) "time step (ps) ", delta_t1

    write( log_file, *) "*******force field*****************************"
    write( log_file, *) " name, mass, chg, lj_param(3), polarizability"
    do i=1, n_atom_type
       write( log_file, *) atype_name(i), atype_mass(i),atype_chg(i), atype_lj_parameter(i,i,1),atype_lj_parameter(i,i,2),atype_lj_parameter(i,i,3), atype_pol(i)
    enddo


  end  subroutine print_run_param


  !**************************************************************************
  ! This subroutine print results to trajectory file and log file
  !**************************************************************************
  subroutine print_result( traj_file, log_file, n_mole, n_atom, xyz, box, i_step, E_elec,E_elec_nopol,E_bh, E_3body, E_bond, E_angle, E_dihedral, potential,new_KE,vol )
    use global_variables
    integer, intent(in) :: n_mole, traj_file, log_file, i_step
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:,:) :: xyz
    real*8, intent(in), dimension(3,3) :: box
    real*8, intent(in) :: potential, E_elec, E_elec_nopol,E_bh,E_3body,E_bond, E_angle, E_dihedral, new_KE,vol
    integer :: i_mole, i_atom, i, j, m_index, a_index


    write( traj_file, * ) '-------------  i_step ', i_step, ' -----------------' ! hua li de fen ge xian

    write( traj_file, * ) total_atoms

    i = 1
    do i_mole = 1, n_mole
       m_index = molecule_index(i_mole)
       do i_atom = 1, n_atom( i_mole )
          a_index = atom_index(i_mole,i_atom)
          write( traj_file, '(I5,2A5,I5,3F14.6)' ) i_mole, molecule_type_name(m_index), atype_name(a_index), i, &
               & xyz(i_mole,i_atom,1), xyz(i_mole,i_atom,2), xyz(i_mole,i_atom,3)
          i = i + 1
       end do
    end do

    do j=1,3
       write( traj_file, '(3F15.8)' ) box(j,1), box(j,2), box(j,3)
    enddo

    Select Case(ms_evb_simulation)
    Case("yes")
       write( log_file, '(I9,3E16.6)' ) i_step, potential
       write( log_file, * )  '------------------------------'
    case default   
       write( log_file, * ) '   i_step  Electro_static      Buckingham      Bond      Angle    Dihedral  Potential'
       write( log_file, '(I9,6E16.6)' ) i_step, E_elec, E_bh, E_bond, E_angle, E_dihedral, potential
       write( log_file, * )  '------------------------------'
    end Select

    write(log_file,*) "long range correction", E_disp_lrc
    write(log_file,*) "KE",new_KE
    write(log_file,*) "molecules", n_mole, "volume", vol

  end subroutine print_result




  !*******************************************
  ! this subroutine print grofile format output
  !*******************************************
  subroutine print_gro_file( i_step , xyz , n_mole , n_atom , box , file_h )
    use global_variables
    integer, intent(in) :: i_step, n_mole, file_h
    real*8,dimension(:,:,:), intent(in) :: xyz
    integer, dimension(:), intent(in) :: n_atom
    real*8,dimension(:,:) , intent(in) :: box

    integer :: i_mole, i_atom, m_index, a_index, count
 
          write( file_h, *) "step ", i_step
          write( file_h, *) total_atoms
          count=1
          do i_mole =1 , n_mole
             m_index = molecule_index(i_mole)
             do i_atom=1,n_atom(i_mole)
                a_index = atom_index(i_mole,i_atom)
                ! print in nanometers for gro file
                write( file_h, '(I5,2A5,I5,3F8.3)' ) i_mole , molecule_type_name(m_index) , atype_name(a_index), count, xyz(i_mole,i_atom,1) / 10d0 , xyz(i_mole,i_atom,2) / 10d0 , xyz(i_mole,i_atom,3) / 10d0
                count = count + 1
             end do
         enddo

         ! this print obviously only works for orthorhombic boxes
         write( file_h, * ) box(1,1) / 10d0 , box(2,2) / 10d0  , box(3,3) / 10d0 

  end subroutine print_gro_file





  !*******************************************
  ! this subroutine prints atomic velocities to a velocity checkpoint file
  ! the primary purpose of this is for continuing a simulation
  !*******************************************
  subroutine print_velocities_checkpoint( i_step , file_h, n_mole , n_atom , vel_atom )
    use global_variables
    integer, intent(in) :: i_step, n_mole, file_h
    integer, dimension(:), intent(in) :: n_atom
    real*8,dimension(:,:,:), intent(in) :: vel_atom

    integer :: i_mole, i_atom, m_index, a_index
 
          write( file_h, *) "step ", i_step
          do i_mole =1 , n_mole
             m_index = molecule_index(i_mole)
             do i_atom=1,n_atom(i_mole)
                a_index = atom_index(i_mole,i_atom)
                write( file_h, '(I5,2A5,I5,3F14.6)' ) i_mole , molecule_type_name(m_index) , atype_name(a_index), i_atom, vel_atom(i_mole,i_atom,1) , vel_atom(i_mole,i_atom,2) , vel_atom(i_mole,i_atom,3) 
             end do
         enddo

  end subroutine print_velocities_checkpoint





  !***********************************************************************
  !  This subroutine updates volume move step size for npt ensemble
  !***********************************************************************
  subroutine update_stepsize_volume ( i_acc_v, i_try_v, i_acc_v_old, i_try_v_old )
    use global_variables
    real*8, intent(inout) :: i_acc_v, i_try_v, i_acc_v_old, i_try_v_old
    real*8   :: acc_pct

!!!!!!!!!!!!!!!!!!update step size, using acceptance percentage in last interval
    acc_pct=(i_acc_v-i_acc_v_old)/(i_try_v-i_try_v_old)

    if(acc_pct < target_acc_v) then
       vol_step = (1d0 - step_update_v ) * vol_step 
    else
       if( vol_step < max_vol_step ) then
          vol_step = (1d0 + step_update_v ) * vol_step   
       endif
    endif

    i_acc_v_old=i_acc_v; i_try_v_old=i_try_v

  end subroutine update_stepsize_volume




  !**************************************************************
  ! This subroutine checks to make sure that molecules are not broken up
  ! by an intra-molecular pbc shift.  It assumes that the atoms are listed
  ! according to relative location in the molecule, so that the shift
  ! is based on the preceding atom
  !**************************************************************
  subroutine check_intra_molecular_shifts(i_mole,n_atom,xyz,box)
    use global_variables
    integer, intent(in) :: i_mole
    integer, dimension(:), intent(in) :: n_atom
    real*8,dimension(:,:,:), intent(inout) :: xyz
    real*8,dimension(:,:), intent(in) :: box

    integer ::  i_atom, j_atom
    real*8,dimension(3) :: shift, drij
    real*8,parameter :: small=1D-6

    ! shift atoms in molecule with respect to the preceding atom
    if ( n_atom(i_mole) > 1 ) then
       do i_atom=2, n_atom(i_mole)
          j_atom = i_atom - 1
          shift(:) = pbc_shift( xyz(i_mole, j_atom,:), xyz(i_mole,i_atom,:) , box , xyz_to_box_transform )

          if ( ( abs( shift(1) ) > small ) .or. ( abs( shift(2)) > small ) .or. ( abs( shift(3) ) > small ) ) then
             drij = pbc_dr( xyz(i_mole,j_atom,:), xyz(i_mole,i_atom,:), shift(:) )
             xyz(i_mole,i_atom,:) = xyz(i_mole,j_atom,:) + drij(:)

          endif
       enddo
    endif

  end subroutine check_intra_molecular_shifts




  !***************************************************************
  ! This subroutine checks that cutoffs are set appropriately for minimum image
  ! convention, and are consistent with our method for finding images
  !***************************************************************
  subroutine check_cutoffs_box( box )
    use global_variables
    implicit none
    real*8, intent(in), dimension(:,:) :: box
    real*8,parameter :: small_test=1d-6
    real*8,parameter :: warn_verlet=1.0



    ! if box is non-cubic, it must satisfy a set of criteria for our method of finding images
    ! criteria is same as Gromacs, namely ay = az = bz = 0 ; ax > 0 , by > 0 , cz > 0 ,
    ! |bx|<.5*ax , |cx|<.5*ax , |cy|<.5*by

    ! see if this box is ok
    if ( ( abs(box(1,2)) > small_test ) .or.  ( abs(box(1,3)) > small_test ) .or. ( abs(box(2,3)) > small_test ) ) then
       stop " box requirements ay = az = bz = 0 are not satisfied "
    else if ( ( box(1,1) .le. 0 ) .or. ( box(2,2) .le. 0 ) .or. ( box(3,3) .le. 0 ) ) then
       stop " box requirements ax > 0 , by > 0 , cz > 0 are not satisfied "
    else if ( (abs(box(2,1)) .gt. .5 * box(1,1) ) .or. (abs(box(3,1)) .gt. .5 * box(1,1) ) .or. (abs(box(3,2)) .gt. .5 * box(2,2) ) ) then
       stop " box requirements |bx|<.5*ax , |cx|<.5*ax , |cy|<.5*by are not satisfied "
    endif

    ! now check cutoff requirements.  For our simple minimum image searching, these are stricter
    ! Rcutoff < .5 * min(ax,by,cz)
    if ( ewald_cutoff .ge. .5 * min( box(1,1), box(2,2), box(3,3) ) ) then
       stop "ewald_cutoff is too big for this box"
    endif
    if ( Electro_cutoff .ge. .5 * min( box(1,1), box(2,2), box(3,3) ) ) then
       stop "Electro_cutoff is too big for this box"
    endif
    if ( lj_cutoff .ge. .5 * min( box(1,1), box(2,2), box(3,3) ) ) then
       stop "lj_cutoff is too big for this box"
    endif

    ! if verlet list, check the verlet skin

    Select Case(use_verlet_list)
    Case("yes")
       if ( verlet_cutoff_lj .ge. .5 * min( box(1,1), box(2,2), box(3,3) ) ) then
          write(*,*)  "verlet skin width is too big for this box"
          write(*,*)  "please set verlet_cutoff to a value less than "
          write(*,*)  "half the box length "
          stop
       endif
       if ( verlet_cutoff_elec .ge. .5 * min( box(1,1), box(2,2), box(3,3) ) ) then
          write(*,*)  "verlet skin width is too big for this box"
          write(*,*)  "please set verlet_cutoff to a value less than "
          write(*,*)  "half the box length "
          stop
       endif

       ! now make sure that the verlet cutoff is longer than the cutoffs for which verlet is used

       if ( verlet_cutoff_lj .le. lj_cutoff ) then
          stop "verlet_cutoff must be set greater than lj_cutoff "
       elseif ( ( verlet_cutoff_lj - lj_cutoff ) .lt. warn_verlet ) then
          write(*,*) ""
          write(*,*) "WARNING:  verlet_cutoff is less than ", warn_verlet
          write(*,*) "Angstrom larger than lj_cutoff.  Are you sure you want"
          write(*,*) "such a thin verlet skin thickness?"
          write(*,*) ""
       endif

       if ( verlet_cutoff_elec .le. ewald_cutoff ) then
          stop "verlet_cutoff must be set greater than ewald_cutoff "
       elseif ( ( verlet_cutoff_elec - ewald_cutoff ) .lt. warn_verlet ) then
          write(*,*) ""
          write(*,*) "WARNING:  verlet_cutoff is less than ", warn_verlet
          write(*,*) "Angstrom larger than ewald_cutoff.  Are you sure you want"
          write(*,*) "such a thin verlet skin thickness?"
          write(*,*) ""
       endif

    End Select


  end subroutine check_cutoffs_box



  !*************************************************************************
  ! This function check the distance between the newly moved particle and other molecules
  !*************************************************************************
  integer function check_distance( xyz, n_mole, n_atom, i_move, xyz_new, box, mass, r_com )
    use global_variables
    implicit none
    integer, intent(in) :: n_mole, i_move
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:,:) :: xyz
    real*8, intent(in), dimension(:,:) :: box, xyz_new, mass, r_com
    integer :: i_atom, j_mole, j_atom
    real*8, dimension(3) :: shift, r_com_i, r_com_j, drij
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: xyz_tmp

    check_distance = 0
    xyz_tmp = xyz
    xyz_tmp(i_move,:,:) = xyz_new(:,:)
    r_com_i(:) = pos_com( xyz_tmp, i_move, n_atom, mass )
    do j_mole = 1, n_mole
       if ( j_mole /= i_move ) then
          r_com_j(:) = r_com( j_mole, : )
          shift(:) = pbc_shift( r_com_i, r_com_j, box , xyz_to_box_transform )
          do i_atom = 1, n_atom( i_move ) 
             do j_atom = 1, n_atom( j_mole )
                drij = pbc_dr( xyz_new(i_atom,:), xyz(j_mole,j_atom,:), shift(:) )

                if ( abs(drij(1)) < too_close .and. abs(drij(2)) < too_close .and. abs(drij(3)) < too_close ) then
                   if ( dot_product( drij, drij ) < too_close**2 ) then
                      check_distance = 1
                   end if
                end if
             end do
          end do
       end if
    end do

  end function check_distance


  !**********************************************************************
  ! this function checks distances between all atoms in the simulation to make sure they are
  ! not "too close" to each other, where "too_close" is a distance that represents a hard wall
  ! potential which is important to use with a Buckingham potential in an MC simulation
  !
  ! Note that this function is order (N^2), so is costly, and it would be better for speed  to check
  ! distances in the loops calculating energy or force.  However, we keep it separate for code
  ! structure and readability advantages
  !**********************************************************************
  integer function check_move( xyz, tot_n_mole, n_mole, n_atom, box, r_com )
    use omp_lib
    use global_variables
    implicit none
    integer, intent(in) :: tot_n_mole, n_mole
    integer, intent(in), dimension(:) :: n_atom
    real*8, intent(in), dimension(:,:,:) :: xyz
    real*8, intent(in), dimension(:,:) :: r_com,box
    integer :: i_mole, i_atom, j_mole, j_atom
    real*8, dimension(3) :: shift, r_com_i, r_com_j, drij
    integer :: split_do

    check_move = 0

    ! decide how to split the parallel section
    if (n_threads .eq. 1 ) then
       split_do = 1
    else
       split_do = n_mole/n_threads+1
    endif

    call OMP_SET_NUM_THREADS(n_threads)
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(tot_n_mole, n_mole,r_com,box,n_atom,xyz,split_do,xyz_to_box_transform, too_close, check_move)
    !$OMP DO SCHEDULE(DYNAMIC,split_do)
    do i_mole = 1, n_mole
       r_com_i(:) = r_com(i_mole,:)
       if ((i_mole .lt. n_mole).or.(n_mole < tot_n_mole)) then
          do j_mole = i_mole+1, tot_n_mole      
             r_com_j(:) = r_com( j_mole, : )
             shift(:) = pbc_shift( r_com_i, r_com_j, box ,xyz_to_box_transform )
             do i_atom = 1, n_atom( i_mole ) 
                do j_atom = 1, n_atom( j_mole )
                   drij = pbc_dr( xyz(i_mole,i_atom,:), xyz(j_mole,j_atom,:), shift(:) )
                   if ( abs(drij(1)) < too_close .and. abs(drij(2)) < too_close .and. abs(drij(3)) < too_close ) then
                      if ( dot_product( drij, drij ) < too_close**2 ) then
                         ! this probably doesn't need to be in critical section, but as it is very unlikely threads will be here at the same time
                         ! it most certainly won't slow anything down
                         !$OMP CRITICAL
                         check_move = 1
                         !$OMP END CRITICAL
                      end if
                   end if
                end do
             end do
          enddo
       end if
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    ! stop if atoms are too close
    if ( check_move == 1 ) then
       write(*,*) ""
       write(*,*) "In the initial configuration, there exists atoms on two or more"
       write(*,*) "molecules that are closer than ", too_close, " Angstroms. "
       write(*,*) ""
       write(*,*) "Please use a different initial configuration, or change the value of "
       write(*,*) "'too_close', only if you believe this is a reasonable separation distance"
       write(*,*) ""
       stop
    end if

  end function check_move



  !*****************************************************
  ! this subroutine translates molecule back into box if it
  ! has been moved out
  !*****************************************************
  subroutine shift_move(n_atom,xyz,r_com,box)
    use global_variables
    implicit none
    integer,intent(in)::n_atom
    real*8,dimension(:,:),intent(inout)::xyz
    real*8,dimension(:),intent(inout)::r_com
    real*8,dimension(:,:),intent(in)::box

    real*8,dimension(3,3)::unitbox
    real*8,dimension(3) :: dr_box,shift,temp_xyz
    integer::i,j,i_atom
    real*8::dist

    ! general box, transform coordinates to box vector
    dr_box = matmul ( xyz_to_box_transform , r_com )
    do i = 1, 3
       shift(i)=0d0
       if ( dr_box(i) < 0d0 ) then
          shift(i)=1d0
       elseif( dr_box(i) > 1d0 ) then
          shift(i)=-1d0
       endif
    enddo
    temp_xyz(:) = 0.0
    do i = 1, 3
       do j = 1 , 3
          temp_xyz(i) = temp_xyz(i) + shift(j)*box(j,i)
       enddo
    enddo
    r_com(:)=r_com(:)+temp_xyz(:)
    do i_atom = 1,n_atom
       xyz(i_atom,:)=xyz(i_atom,:) + temp_xyz(:)
    enddo


  end subroutine shift_move




  !*********************************************
  ! this subroutine allocates the size of the verlet list
  ! we have a subroutine to do this, because the verlet list can take up
  ! a significant amount of memory if not allocated intelligently, and 
  ! the required size can change if the density of molecules in the simulation changes
  !  Verlet list size should be about 4*pi*Rverlet^3 * Natoms^2 / 6 * volume
  !*********************************************
  subroutine allocate_verlet_list( tot_n_mole, n_atom, chg, box )
    use global_variables
    integer, intent(in) :: tot_n_mole
    integer, dimension(:),intent(in) :: n_atom
    real*8,dimension(:,:),intent(in) :: chg
    real*8,dimension(:,:),intent(in) :: box

    real*8 :: vol
    integer :: size_verlet_lj, size_verlet_elec
    real*8 :: size_lj , size_elec

    if ( allocated ( verlet_neighbor_list_lj ) ) then
       deallocate( verlet_neighbor_list_lj )
    endif
    if ( allocated ( verlet_neighbor_list_elec ) ) then
       deallocate( verlet_neighbor_list_elec )
    endif
    if ( allocated ( verlet_point_lj ) ) then
       deallocate( verlet_point_lj )
    endif
    if ( allocated ( verlet_point_elec ) ) then
       deallocate( verlet_point_elec )
    endif

    ! generate lennard jones verlet atom index for this atom
    call generate_verlet_atom_index( verlet_lj_atoms, verlet_molecule_map_lj, tot_n_mole, n_atom )
    ! generate electrostatic verlet atom index for this atom
    call generate_verlet_atom_index( verlet_elec_atoms, verlet_molecule_map_elec, tot_n_mole, n_atom, chg )

    ! size verlet point should be equal to the number of atoms
    ! with non-zero interaction (plus one)
    allocate( verlet_point_lj(verlet_lj_atoms+1), verlet_point_elec(verlet_elec_atoms+1) )

    vol= volume( box(1,:),box(2,:),box(3,:) )

    size_lj = 4d0 * pi * verlet_cutoff_lj**3 * dble(verlet_lj_atoms)**2 / 6d0 / vol
    size_elec = 4d0 * pi * verlet_cutoff_elec**3 * dble(verlet_elec_atoms)**2 / 6d0 / vol

    ! for huge box, value of temp could be zero if we just have gas phase dimer
    ! verlet list should be at least as big as number of molecules
    size_verlet_lj = max( verlet_lj_atoms , floor( size_lj ) )
    size_verlet_elec = max( verlet_elec_atoms , floor( size_elec ) )

    ! safe_verlet is factor that we multiply theoretically needed size of verlet list to be safe
    ! found in global_variables
    size_verlet_lj = floor( dble(size_verlet_lj) * safe_verlet )
    size_verlet_elec = floor( dble(size_verlet_elec) * safe_verlet )

    allocate( verlet_neighbor_list_lj(size_verlet_lj) , verlet_neighbor_list_elec(size_verlet_elec)  )

  end subroutine allocate_verlet_list



  !**********************************************
  ! this subroutine keeps track of atomic displacements since the
  ! last construction of the Verlet list
  ! when the displacements get too large, a flag is turned on to tell
  ! the program to reconstruct the verlet list
  !
  ! also this subroutine will check the density of the system, and make sure it's not changing much
  ! if it does change significantly, the verlet list will be reallocated to accommodate this
  ! since verlet list can only be used in nvt, npt ensemble (or else code will stop), only care
  ! about the volume change
  !
  ! OUTPUT: flag_verlet_list = 0, nothing needs to be done
  !         flag_verlet_list = 1, verlet list needs to be updated, but not reallocated
  !         flag_verlet_list = 2, verlet list needs to be reallocated and updated
  !*********************************************
  subroutine update_verlet_displacements( n_mole, n_atom, xyz, box, flag_verlet_list, flag_initialize )
    use global_variables
    integer,intent(in) :: n_mole
    integer,dimension(:),intent(in) :: n_atom
    real*8,dimension(:,:,:),intent(in) :: xyz
    real*8,dimension(:,:),intent(in) :: box
    integer, intent(out) :: flag_verlet_list
    integer, intent(in), optional :: flag_initialize

    real*8,dimension(MAX_N_MOLE,MAX_N_ATOM,3), save :: r_atom_last, r_verlet_displacement
    real*8, save  :: vol_last
    integer :: i_mole, i_atom
    real*8, dimension(3) :: r_i_new, r_i_old, shift, dr
    real*8  :: max_d1 , max_d2, max_temp, norm_dr , verlet_skin, vol, max_delta_vol

    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "verlet list check displacements started at", time
    endif
    !***********************************************************************!

    ! if initialization
    if ( present(flag_initialize) ) then
       r_atom_last = xyz
       vol_last= volume( box(1,:),box(2,:),box(3,:) )
       r_verlet_displacement=0d0
       flag_verlet_list = 0
    else
       ! first make sure volume hasn't changed significantly (npt ensemble), or else reallocate and reconstruct verlet list
       vol = volume( box(1,:),box(2,:),box(3,:) )
       ! maximum allowed volume change
       max_delta_vol = ( safe_verlet - 1d0 ) / 2d0
       if ( (abs ( vol - vol_last ) / vol_last ) > max_delta_vol ) then
          ! volume has changed too much, reallocate size of verlet_neighbor_list
          flag_verlet_list = 2
       else
          ! don't need to reallocate, but check if we need to update verlet list
          max_d1=0d0; max_d2=0d0
          ! add new displacements and check max displacement
          do i_mole =1,n_mole
             do i_atom=1,n_atom(i_mole)
                r_i_old(:) = r_atom_last( i_mole, i_atom,: ) 
                r_i_new(:) = xyz( i_mole, i_atom, : )

                shift = pbc_shift( r_i_old, r_i_new, box , xyz_to_box_transform )
                dr = pbc_dr( r_i_old, r_i_new, shift )
                r_verlet_displacement(i_mole,i_atom,:) = r_verlet_displacement(i_mole,i_atom,:) + dr(:)

                ! check whether this is either 1st or 2nd biggest displacement
                norm_dr = sqrt( dot_product( r_verlet_displacement(i_mole,i_atom,:),r_verlet_displacement(i_mole,i_atom,:)) )
                ! store biggest displacement in max_d1, 2nd biggest in max_d2
                if ( norm_dr > max_d2 ) then
                   max_d2 = norm_dr
                   ! if this is biggest displacement, switch positions
                   if ( max_d2 > max_d1 ) then
                      max_temp = max_d1
                      max_d1 = max_d2
                      max_d2 = max_temp
                   end if
                end if

             end do
          enddo

          ! now store these most recent center of mass positions
          r_atom_last = xyz

          ! now see whether we need to update the verlet list
          verlet_skin = verlet_thresh * min ( (verlet_cutoff_lj - lj_cutoff) , (verlet_cutoff_elec - ewald_cutoff) )

          if ( ( max_d1 + max_d2 ) > verlet_skin ) then
             flag_verlet_list = 1 
          else
             flag_verlet_list = 0
          endif

       end if

    end if

    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "verlet list check displacements finished at", time
    endif
    !***********************************************************************!

  end subroutine update_verlet_displacements




  !**********************************************
  ! this subroutine constructs the verlet neighbor list that is used for 
  ! looping over atom-atom (intermolecular) interactions
  ! 
  ! See "Computer Simulation of Liquids" by Allen and Tildesley, Oxford Science Publications 2009
  ! for details
  !
  ! basically we construct an array "verlet_neighbor_list" which holds a list of all the neighbor atoms
  ! the array "verlet_point" points to the index of verlet_neighbor_list for an appropriate atom,
  ! for instance the neighbors of molecule i_atom are found in indices of verlet_neighbor_list between
  ! verlet_point(i_atom) and verlet_point(i_atom+1)
  !
  ! note that this means we need to set the value of verlet_point(last_atom+1), so that we know the finish position for
  ! neighbors of last_atom
  !*********************************************
  subroutine construct_verlet_list(tot_n_mole, n_mole, n_atom, xyz, box )
    use global_variables
    integer,intent(in) :: n_mole, tot_n_mole
    integer,dimension(:),intent(in) :: n_atom
    real*8,dimension(:,:,:),intent(in) :: xyz
    real*8,dimension(:,:),intent(in) :: box

    integer :: i_mole, j_mole, i_atom, j_atom, count_atom, verlet_index, a_index, last_atom
    real*8,dimension(3) :: r_ij, shift, norm_dr
    real*8 :: verlet_cutoff_lj2, verlet_cutoff_elec2
    real*8,parameter :: small = 1d-6

    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "verlet list construction started at", time
    endif
    !***********************************************************************!

    ! need to change code if n_mole not equal tot_n_mole
    if ( n_mole /= tot_n_mole) then
       write(*,*) "verlet list currently not supported for n_mole /= tot_n_mole"
       write(*,*) " this is because in the energy and force subroutines, we loop"
       write(*,*) " over either verlet_lj_atoms or verlet_elec_atoms which are "
       write(*,*) " the total number of interacting atoms in the system, not   "
       write(*,*) " the number of solute atoms for which the verlet list is constructed"
       write(*,*) " this needs to be modfied"
       stop
    end if

    Select Case(verlet_grid_based_construction)
    Case("yes")
       ! lj verlet list
       call construct_verlet_list_grid(tot_n_mole, n_mole, n_atom, xyz, box, 1 )
       ! elec verlet list
       call construct_verlet_list_grid(tot_n_mole, n_mole, n_atom, xyz, box, 2 )
    Case("no")
       ! this constructs both elec and lj verlet list
       call construct_verlet_list_loop(tot_n_mole, n_mole, n_atom, xyz, box )
    End Select


    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "verlet list construction finished at", time
    endif
    !***********************************************************************!


  end subroutine construct_verlet_list




  !****************************************
  ! this subroutine constructs the verlet lists using
  ! a cell list
  ! if flag_lj = 1, constructs lj_verlet list
  ! if flag_lj = 1, constructs elec verlet list
  !****************************************
  subroutine construct_verlet_list_grid(tot_n_mole, n_mole, n_atom, xyz, box, flag_lj )
    use global_variables
    integer,intent(in) :: n_mole, tot_n_mole, flag_lj
    integer,dimension(:),intent(in) :: n_atom
    real*8,dimension(:,:,:),intent(in) :: xyz
    real*8,dimension(:,:),intent(in) :: box

    integer :: ix, iy, iz, nx, ny, nz,i_end_atom, i_atom_head
    integer :: dia, dib, dic, igrid1, igrid2, igrid3, ia, ib, ic
    real*8, dimension(3) :: rtmp

    real*8 :: rka, rkb, rkc
    integer :: i,i_mole, j_mole, i_atom, j_atom, k_atom,count_atom, verlet_index, a_index, last_atom, verlet_neighbor_list_size
    real*8,dimension(3) :: r_ij, shift, norm_dr, dr_direct, shift_direct
    real*8 :: verlet_cutoff2
    real*8,parameter :: small = 1d-6
    ! this is a neighbor list local to each grid cell
    ! nslist_head, stores the starting atom for the neighbor list
    ! of the cell that the index atom is in
    integer, dimension(:), allocatable  :: nslist_cell, index_molecule
    integer, dimension(:,:,:), allocatable ::  endlist_cell, headlist_cell
    ! we use this kind of like a hash, to look up the cell index of an atom
    ! in order to do this, we split the string
    character(10), dimension(:), allocatable :: cell_index_hash
    real*8,dimension(:,:), allocatable :: atom_list_xyz

    nx = na_nslist ; ny = nb_nslist ; nz = nc_nslist ;
    ! hash won't work if these are not 2 digit
    if ( ( nx < 10 ) .or. ( ny < 10 ) .or. ( nz < 10  ) .or.  ( nx > 99 ) .or. ( ny > 99 ) .or. ( nz > 99  )) then
       stop " na_nslist, nb_nslist, nc_nslist must be between 10 and 100 "
    end if

    Select Case( flag_lj )
    Case(1)
       verlet_cutoff2 = verlet_cutoff_lj**2
       verlet_neighbor_list_size = size(verlet_neighbor_list_lj)
    Case(2)
       verlet_cutoff2 = verlet_cutoff_elec**2
       verlet_neighbor_list_size = size(verlet_neighbor_list_elec)
    End Select

    verlet_index = 1
    count_atom = 1

    ! first construct the cell-based neighbor list
    Select Case( flag_lj )
    Case(1)
       allocate( atom_list_xyz(verlet_lj_atoms,3), nslist_cell(verlet_lj_atoms), cell_index_hash(verlet_lj_atoms) , index_molecule(verlet_lj_atoms), headlist_cell(na_nslist,nb_nslist,nc_nslist), endlist_cell(na_nslist,nb_nslist,nc_nslist)  )
    Case(2)
       allocate( atom_list_xyz(verlet_elec_atoms,3), nslist_cell(verlet_elec_atoms), cell_index_hash(verlet_elec_atoms) , index_molecule(verlet_elec_atoms), headlist_cell(na_nslist,nb_nslist,nc_nslist), endlist_cell(na_nslist,nb_nslist,nc_nslist)  )
    End Select

    !   pointers all set to null
    headlist_cell = 0
    endlist_cell = 0
    nslist_cell = 0


    do i_mole=1, tot_n_mole
       do i_atom=1,n_atom(i_mole)

          ! get verlet atom index for this atom
          Select Case( flag_lj )
          Case(1)
             a_index = verlet_molecule_map_lj(i_mole,i_atom)
          Case(2)
             a_index = verlet_molecule_map_elec(i_mole,i_atom)
          End Select

          if ( a_index > 0  ) then
             ! store atom information from molecular data arrays
             atom_list_xyz(a_index,:) = xyz(i_mole,i_atom,:)
             ! store molecule index, this will be used later to avoid putting intramolecular atoms on neighbor list
             index_molecule(a_index) = i_mole

             rtmp = matmul(xyz_to_box_transform(:,:), xyz(i_mole,i_atom,:) )

             ix = floor(rtmp(1)*nx) + 1
             iy = floor(rtmp(2)*ny) + 1
             iz = floor(rtmp(3)*nz) + 1
             !       shift if it is out of the box
             ix = ix - floor(dble(ix-1)/nx)*nx
             iy = iy - floor(dble(iy-1)/ny)*ny
             iz = iz - floor(dble(iz-1)/nz)*nz
             ! construct cell index hash for this atom, 1 means construct
             call cell_index_hash_io( cell_index_hash, a_index, ix,iy,iz, 1 )

             if ( headlist_cell(ix,iy,iz) == 0 ) then
                headlist_cell(ix,iy,iz) = a_index
                endlist_cell(ix,iy,iz) = a_index
             else
                i_end_atom = endlist_cell(ix,iy,iz)
                nslist_cell(i_end_atom) = a_index
                endlist_cell(ix,iy,iz) = a_index
             endif
          end if
       enddo
    enddo

    ! this is the number of adjacent cells to search
    ! see explicit_three_body_interactions for comments on this formula
    rka = dsqrt(dot_product(xyz_to_box_transform(1,:),xyz_to_box_transform(1,:)))
    rkb = dsqrt(dot_product(xyz_to_box_transform(2,:),xyz_to_box_transform(2,:)))
    rkc = dsqrt(dot_product(xyz_to_box_transform(3,:),xyz_to_box_transform(3,:)))
    ! get verlet atom index for this atom
    Select Case( flag_lj )
    Case(1)
       dia = floor(verlet_cutoff_lj*rka*dble(na_nslist)) + 1
       dib = floor(verlet_cutoff_lj*rkb*dble(nb_nslist)) + 1
       dic = floor(verlet_cutoff_lj*rkc*dble(nc_nslist)) + 1
    Case(2)
       dia = floor(verlet_cutoff_elec*rka*dble(na_nslist)) + 1
       dib = floor(verlet_cutoff_elec*rkb*dble(nb_nslist)) + 1
       dic = floor(verlet_cutoff_elec*rkc*dble(nc_nslist)) + 1
    End Select

    ! note dia should be less than 1/2 * na_nslist, and similarly for dib, dic
    ! otherwise code below can loop over the same cell twice if dia = 1/2 na_nslist,
    ! and this will cause a bug.  We didn't have this problem in the explicit_three_body_interactions
    ! routines, because the cutoff was always short.  Therefore, here we have to add a check on this.
    call check_neighbor_list_grid_cutoff_consistency( dia, dib , dic , na_nslist, nb_nslist , nc_nslist, rka, rkb, rkc )


    ! now form verlet neighbor list using this cell list
    do i_mole =1,n_mole
       do i_atom =1, n_atom(i_mole)
          ! get verlet atom index for this atom
          Select Case( flag_lj )
          Case(1)
             a_index = verlet_molecule_map_lj(i_mole,i_atom)
          Case(2)
             a_index = verlet_molecule_map_elec(i_mole,i_atom)
          End Select

          if ( a_index > 0  ) then
             last_atom = a_index
             ! set verlet point for this atom
             Select Case( flag_lj )
             Case(1)
                verlet_point_lj(a_index) = verlet_index
             Case(2)
                verlet_point_elec(a_index) = verlet_index
             end Select

             ! get cell index for this atom, 2 means read
             call cell_index_hash_io (  cell_index_hash, a_index, ix,iy,iz, 2 )

             ! loop over cells within cutoff distance from first cell
             do ia = -dia, dia
                igrid1 = ix + ia
                igrid1 = igrid1 - floor(dble(igrid1-1)/dble(na_nslist))*na_nslist
                do ib = -dib, dib
                   igrid2 = iy + ib
                   igrid2 = igrid2 - floor(dble(igrid2-1)/dble(nb_nslist))*nb_nslist
                   do ic = -dic, dic
                      igrid3 = iz + ic
                      igrid3 = igrid3 - floor(dble(igrid3-1)/dble(nc_nslist))*nc_nslist

                      !********************* second atom
                      j_atom = headlist_cell(igrid1, igrid2, igrid3)

                      ! loop over all atoms in this cell
                      do
                         if ( j_atom == 0 ) exit
                         ! only add neighbors with index > a_index to avoid double counting
                         if ( a_index < j_atom ) then
                            ! make sure these atoms aren't on the same molecule
                            j_mole = index_molecule(j_atom)
                            if ( i_mole /= j_mole) then
                               ! remember, we are using atom_list_xyz to look up coordinates here, because j_atom is an absolute atom index

                               r_ij = xyz(i_mole, i_atom, : ) - atom_list_xyz(j_atom, : )

                               ! shift for general box
!!$                               dr_direct(:) = matmul( xyz_to_box_transform, r_ij )
!!$                               do i=1,3
!!$                                  shift_direct(i) = dble(floor( dr_direct(i) + 0.5d0 ))
!!$                               enddo
!!$                               shift = matmul( shift_direct , box )
                               ! shift for orthorhombic box
                               do i=1,3
                                  shift(i) = box(i,i) * floor( r_ij(i) / box(i,i) + 0.5d0 )
                               end do

                               r_ij = r_ij - shift

                               if ( dot_product( r_ij, r_ij ) < verlet_cutoff2 ) then

                                  if ( verlet_index > verlet_neighbor_list_size ) then
                                     write(*,*) "please increase size of verlet neighbor list"
                                     stop
                                  end if

                                  ! add this atom to verlet list
                                  Select Case( flag_lj )
                                  Case(1)
                                     verlet_neighbor_list_lj( verlet_index ) = j_atom
                                  Case(2)
                                     verlet_neighbor_list_elec( verlet_index ) = j_atom
                                  end Select

                                  verlet_index = verlet_index + 1
                               end if
                            end if
                         end if
                         k_atom = nslist_cell(j_atom)
                         j_atom = k_atom
                      enddo

                   enddo
                enddo
             enddo
          end if
       enddo
    enddo

    if ( n_mole == tot_n_mole ) then
       ! last atom shouldn't have any neighbors
       Select Case( flag_lj )
       Case(1)
          verlet_point_lj(last_atom+1) = verlet_index
       Case(2)
          verlet_point_elec(last_atom+1) = verlet_index
       End Select
    end if

    deallocate( atom_list_xyz, nslist_cell, cell_index_hash, index_molecule, headlist_cell, endlist_cell )

  end subroutine construct_verlet_list_grid





  !***********************************
  ! this subroutine makes sure that the interaction cutoff
  ! is appropriate for the neighborlist grid search
  ! dia, dib, dic should be the number of grid cells in each dimension
  ! to search for neighbors, as determined by the cutoff
  ! basically, these values should be less than 1/2 the total grid cells
  ! so that our current algorithm doesn't lead to a bug
  !***********************************
  subroutine check_neighbor_list_grid_cutoff_consistency( dia, dib , dic , na_nslist, nb_nslist , nc_nslist, rka, rkb, rkc )
    integer, intent(in) :: dia, dib, dic, na_nslist, nb_nslist, nc_nslist
    real*8, intent(in) :: rka, rkb, rkc

    integer :: diam, dibm, dicm
    real*8 :: max_cutoff, cutx,cuty,cutz

    ! note integer arithmetic
    if ( (dia >= na_nslist/2) .or.  (dib >= nb_nslist/2) .or. (dic >= nc_nslist/2) ) then

       ! calculate maximum cutoff that can be used with current algorithm
       diam = na_nslist/2
       dibm = nb_nslist/2
       dicm = nc_nslist/2
       cutx = dble( diam - 1 ) / rka / dble(na_nslist)
       cuty = dble( dibm - 1 ) / rkb / dble(nb_nslist)
       cutz = dble( dicm - 1 ) / rkc / dble(nc_nslist)

       max_cutoff = min( cutx, cuty, cutz )

       write(*,*) ""
       write(*,*) "for current grid-based neighbor searching, number of grid cells"
       write(*,*) "to search in each dimension must be less than half the total number"
       write(*,*) "of grid cells.  To achieve this, decrease cutoff (either interaction"
       write(*,*) "or verlet) to less than ", max_cutoff
       stop

    end if

  end subroutine check_neighbor_list_grid_cutoff_consistency






  !*************************************
  ! this subroutine acts like a hash, by
  ! writing to a string, and then reading
  ! three integers from a string
  !
  ! input "1" is construch hash element
  ! input "2" is retrieve 3 indices
  !*************************************
  subroutine cell_index_hash_io( cell_index_hash, a_index, ix,iy,iz, flag )
    character(*),dimension(:), intent(inout) :: cell_index_hash
    integer, intent(inout) :: a_index, ix, iy, iz
    integer, intent(in) :: flag
    character(2) :: int1, int2, int3
    character(2),dimension(3) :: args
    integer :: nargs

    if ( flag == 1 ) then
       ! create hash
       write(int1,'(I2)') ix
       write(int2,'(I2)') iy
       write(int3,'(I2)') iz
       cell_index_hash(a_index) = int1 // "_" // int2 // "_" // int3
    else if (flag == 2 ) then
       ! read hash
       call parse( cell_index_hash(a_index) , "_", args , nargs )
       if ( nargs /= 3 ) then
          stop "error in cell_index_hash_io"
       end if
       read(args(1), '(I2)') ix
       read(args(2), '(I2)') iy
       read(args(3), '(I2)') iz
    end if

  end subroutine cell_index_hash_io


  !**************************************
  ! this constructs the verlet lists in a naive
  ! fashion, by looping over atoms
  !**************************************
  subroutine construct_verlet_list_loop(tot_n_mole, n_mole, n_atom, xyz, box )
    use global_variables
    integer,intent(in) :: n_mole, tot_n_mole
    integer,dimension(:),intent(in) :: n_atom
    real*8,dimension(:,:,:),intent(in) :: xyz
    real*8,dimension(:,:),intent(in) :: box

    integer :: i_mole, j_mole, i_atom, j_atom, count_atom, verlet_index, a_index, last_atom
    real*8,dimension(3) :: r_ij, shift, norm_dr
    real*8 :: verlet_cutoff_lj2, verlet_cutoff_elec2
    real*8,parameter :: small = 1d-6

    ! *************************** lj verlet list *************************
    verlet_cutoff_lj2 = verlet_cutoff_lj**2

    verlet_index = 1
    count_atom = 1

    do i_mole =1,n_mole
       do i_atom =1, n_atom(i_mole)
          ! get verlet atom index for this atom
          a_index = verlet_molecule_map_lj(i_mole,i_atom)
          if ( a_index > 0  ) then
             last_atom = a_index
             ! set verlet point for this atom
             verlet_point_lj(a_index) = verlet_index
             if((i_mole < n_mole).or.(n_mole < tot_n_mole)) then
                do j_mole=i_mole+1,tot_n_mole
                   do j_atom=1, n_atom(j_mole)

                      ! get verlet atom index for this atom
                      a_index = verlet_molecule_map_lj(j_mole,j_atom)
                      if ( a_index > 0  ) then

                         shift = pbc_shift( xyz(i_mole, i_atom, : ) , xyz(j_mole, j_atom, : ) , box , xyz_to_box_transform )
                         r_ij = pbc_dr( xyz(i_mole,i_atom,:), xyz(j_mole,j_atom,:), shift )

                         if ( dot_product( r_ij, r_ij ) < verlet_cutoff_lj2 ) then

                            if ( verlet_index > size(verlet_neighbor_list_lj) ) then
                               write(*,*) "please increase size of verlet neighbor list"
                               stop
                            end if

                            ! add this atom to verlet list
                            verlet_neighbor_list_lj( verlet_index ) = a_index
                            verlet_index = verlet_index + 1
                         endif
                      endif
                   enddo
                enddo
             endif
          end if
       enddo
    enddo


    if ( n_mole == tot_n_mole ) then
       ! last atom shouldn't have any neighbors
       verlet_point_lj(last_atom+1) = verlet_index
    end if



    !******************************** elec verlet list ************************************
    verlet_cutoff_elec2 = verlet_cutoff_elec**2

    verlet_index = 1
    count_atom = 1

    do i_mole =1,n_mole
       do i_atom =1, n_atom(i_mole)
          ! get verlet atom index for this atom
          a_index = verlet_molecule_map_elec(i_mole,i_atom)
          if ( a_index > 0  ) then
             last_atom = a_index
             ! set verlet point for this atom
             verlet_point_elec(a_index) = verlet_index
             if((i_mole < n_mole).or.(n_mole < tot_n_mole)) then
                do j_mole=i_mole+1,tot_n_mole
                   do j_atom=1, n_atom(j_mole)

                      ! get verlet atom index for this atom
                      a_index = verlet_molecule_map_elec(j_mole,j_atom)
                      if ( a_index > 0  ) then

                         shift = pbc_shift( xyz(i_mole, i_atom, : ) , xyz(j_mole, j_atom, : ) , box , xyz_to_box_transform )
                         r_ij = pbc_dr( xyz(i_mole,i_atom,:), xyz(j_mole,j_atom,:), shift )

                         if ( dot_product( r_ij, r_ij ) < verlet_cutoff_elec2 ) then

                            if ( verlet_index > size(verlet_neighbor_list_elec) ) then
                               write(*,*) "please increase size of verlet neighbor list"
                               stop
                            end if

                            ! add this atom to verlet list
                            verlet_neighbor_list_elec( verlet_index ) = a_index
                            verlet_index = verlet_index + 1
                         endif
                      endif
                   enddo
                enddo
             endif
          end if
       enddo
    enddo


    if ( n_mole == tot_n_mole ) then
       ! last atom shouldn't have any neighbors
       verlet_point_elec(last_atom+1) = verlet_index
    end if


  end subroutine construct_verlet_list_loop



  !****************************************
  ! this subroutine assigns a total atom index for a given atom on a molecule
  ! to be used in the verlet list. The atom indices are assigned by looping over
  ! all atoms of all molecules in order, and every atom is assigned the increment index
  ! except atoms that do not interact with the specific type of interaction
  ! if chg array is present, we consider electrostatic interactions
  ! otherwise lj interactions
  ! verlet_molecule_map_(lj,elec)
  !
  ! note this can be used to construct any atom_molecule_map, since the map array is
  ! explicitly passed as output
  !****************************************
  subroutine generate_verlet_atom_index( total_atoms_list, atom_molecule_map, n_mole, n_atom, chg )
    integer, intent(out) :: total_atoms_list
    integer, dimension(:,:), intent(out) :: atom_molecule_map
    integer, intent(in) :: n_mole
    integer,dimension(:),intent(in) :: n_atom
    real*8,dimension(:,:),optional,intent(in) :: chg

    real*8,parameter :: small = 1d-6
    integer :: i_mole, i_atom, a_index, atom_id1

    a_index=1
    if ( present(chg) ) then
       ! generating electrostatic atom map
       ! intialize
       atom_molecule_map=0
       do i_mole=1, n_mole
          do i_atom=1,n_atom(i_mole)

             ! see if this atom has charge
             if (  abs(chg(i_mole,i_atom)) > small )   then
                ! assign atom index
                atom_molecule_map(i_mole,i_atom) = a_index
                a_index = a_index + 1
             end if
          end do
       end do

       total_atoms_list = a_index - 1

    else
       ! generating lennard jones atom map
       ! intialize
       atom_molecule_map=0
       do i_mole=1, n_mole
          do i_atom=1,n_atom(i_mole)

             ! NOTE the below commented lines will give a bug if we have explicit cross interactions, and so better leave them commented
             ! just in case, and accept a bigger Verlet list size
             ! see if this atom has LJ interaction
             !  atom_id1 = atom_index(i_mole,i_atom)
             !  if ( ( atype_lj_parameter(atom_id1,atom_id1,1) > small ) .or. ( atype_lj_parameter(atom_id1,atom_id1,2) > small ) ) then
             ! assign atom index
             atom_molecule_map(i_mole,i_atom) = a_index
             a_index = a_index + 1
             !  end if
          end do
       end do

       total_atoms_list = a_index - 1

    end if

  end subroutine generate_verlet_atom_index




  !*******************************
  ! this subroutine acts like a hash for
  ! the global variable atype_name array.
  ! Rather than looking up the atom name
  ! using the integer index, it finds the
  ! atomtype index using the name
  !********************************
  subroutine atype_name_reverse_lookup( atomtype, index )
    use global_variables
    character(*), intent(in) :: atomtype
    integer, intent(out) :: index
    integer :: i_atom

    index=-1
    do i_atom=1, n_atom_type
       if ( atomtype .eq. atype_name(i_atom) ) then
          index = i_atom
          exit
       end if
    end do

    if ( index < 0 ) then
       write(*,*) ""
       write(*,*) "can't match atomtype ", atomtype 
       write(*,*) "to any defined atomtype "
       write(*,*) ""
       stop
    end if

  end subroutine atype_name_reverse_lookup


  !*******************************
  ! this subroutine acts like a hash for
  ! the global variable molecule_type_name array.
  ! Rather than looking up the molecule name
  ! using the integer index, it finds the
  ! molecule_type index using the name
  !********************************
  subroutine mtype_name_reverse_lookup( moleculetype, index )
    use global_variables
    character(*), intent(in) :: moleculetype
    integer, intent(out) :: index
    integer :: i_mole

    index=-1
    do i_mole=1, n_molecule_type

       if ( moleculetype .eq. molecule_type_name(i_mole) ) then
          index = i_mole
          exit
       end if
    end do

!!$    if ( index < 0 ) then
!!$       write(*,*) ""
!!$       write(*,*) "can't match moleculetype ", moleculetype 
!!$       write(*,*) "in topology file to any of the moleculetypes "
!!$       write(*,*) "in configuration file"
!!$       write(*,*) ""
!!$       stop
!!$    end if

  end subroutine mtype_name_reverse_lookup





  !****************************************
  ! initializes factorial array for damping functions
  !****************************************
  subroutine initialize_factorial
    use global_variables
    integer :: i
    factorial(1)=1d0
    do i=2,size(factorial)
       factorial(i) = dble(i) * factorial(i-1)
    enddo
  end subroutine initialize_factorial




  !********************************************************
  ! this function calls the Tang-Toennies damping function for R6,R8,R10,R12 dispersion interaction
  ! even though it is the same functional form as that for electrostatics, (just a different summation
  ! limit), it is simpler to put this in its own subroutine as the electrostatic screening function 
  ! routine is a bit messy thanks to the drude oscillators
  !********************************************************
  function C6_C10_damp(atom_id1,atom_id2,norm_dr,n)
    use global_variables
    real*8            :: C6_C10_damp
    integer,intent(in):: atom_id1,atom_id2,n
    real*8,intent(in) :: norm_dr

    real*8 :: lambda
    integer :: index

    ! buckingham exponent is the second parameter in atype_lj_parameter array, use this for screening
    lambda = atype_lj_parameter(atom_id1,atom_id2,2) * norm_dr

    ! decide whether we are using a table lookup for these interactions
    Select Case(grid_Tang_Toennies)
    Case("yes")
       if ( lambda > Tang_Toennies_max ) then
          C6_C10_damp = 1d0
       else
          ! C6 is in 1st index, C8 in 2nd, C10 in 3rd, C12 in 4th
          index = ( n - 6 ) / 2 + 1
          C6_C10_damp = Tang_Toennies_table(index,ceiling(lambda/Tang_Toennies_max*dble(Tang_Toennies_grid)))
       endif

    Case("no")
       C6_C10_damp = Tang_Toennies_damp(lambda,n)
    End Select

  end function C6_C10_damp


  !********************************************************
  ! this function is identical to the function C6_C10_damp , except
  ! it is used for 3rd order Tang-Toennies damping functions which are
  ! used to damp Axilrod-Teller interactions.  The reason that we don't combine these 
  ! two subroutines is that they are both called from inner loops, and we want to
  ! avoid as many case and if statements as possible.  Also, we always damp
  ! Axilrod teller interactions
  !********************************************************
  function C3_damp(atom_id1,atom_id2,norm_dr)
    use global_variables
    real*8            :: C3_damp
    integer,intent(in):: atom_id1,atom_id2
    real*8,intent(in) :: norm_dr
    real*8 :: lambda
    integer :: n

    ! third order damping
    n=3

    ! buckingham exponent is the second parameter in atype_lj_parameter array, use this for screening
    lambda = atype_lj_parameter(atom_id1,atom_id2,2) * norm_dr

    ! decide whether we are using a table lookup for these interactions
    Select Case(grid_Tang_Toennies)
    Case("yes")
       if ( lambda > Tang_Toennies_max ) then
          C3_damp = 1d0
       else
          C3_damp = Tang_Toennies_table3(ceiling(lambda/Tang_Toennies_max*dble(Tang_Toennies_grid)))
       endif
    Case("no")
       C3_damp = Tang_Toennies_damp(lambda,n)
    End Select

  end function C3_damp


  !******************************************
  ! This returns the Tang_Toennies damping function of order n for argument x
  !******************************************
  function Tang_Toennies_damp(x,n)
    use global_variables
    real*8  :: Tang_Toennies_damp
    real*8,intent(in) :: x
    integer,intent(in) :: n

    integer :: i
    real*8 :: sum, xn

    ! sum = 1 from 0 order term which should be included in damping
    sum=1d0
    xn=1d0
    do i=1,n
       xn=xn*x
       sum=sum + xn/factorial(i)
    enddo

    Tang_Toennies_damp = 1d0 - sum * exp(-x)

  end function Tang_Toennies_damp


  !******************************************
  ! This returns the derivative w.r.t x of Tang_Toennies damping function of order n for argument x
  !******************************************
  function dTang_Toennies_damp(x,n)
    use global_variables
    real*8  :: dTang_Toennies_damp
    real*8,intent(in) :: x
    integer,intent(in) :: n

    ! two terms here, mostly cancel to give 
    dTang_Toennies_damp = exp(-x  ) * x ** n  / factorial(n)

  end function dTang_Toennies_damp


  !********************************************************
  ! this function is the gradient of the Tang-Toennies damping function for R6,R8,R10,R12 dispersion interaction
  !********************************************************
  function dC6_C10_damp(atom_id1,atom_id2,r_vec,norm_dr,n)
    real*8,dimension(3) :: dC6_C10_damp
    integer,intent(in):: atom_id1,atom_id2,n
    real*8,dimension(3),intent(in) :: r_vec
    real*8,intent(in) :: norm_dr

    dC6_C10_damp = dC6_C10_damp_scalar1(atom_id1,atom_id2,norm_dr,n) * r_vec / norm_dr 

  end function dC6_C10_damp



  !********************************************************
  ! this function returns the first order derivative of the 
  ! Tang-Toennies damping function for R6,R8,R10,R12 dispersion interaction
  ! this can either use an explicit evaluation or table loop up
  !********************************************************
  function dC6_C10_damp_scalar1(atom_id1,atom_id2,norm_dr,n)
    use global_variables
    real*8 :: dC6_C10_damp_scalar1
    integer,intent(in):: atom_id1,atom_id2,n
    real*8,intent(in) :: norm_dr
    real*8 :: exponent, lambda
    integer :: index



    ! buckingham exponent is the second parameter in atype_lj_parameter array, use this for screening
    exponent = atype_lj_parameter(atom_id1,atom_id2,2)
    lambda = exponent * norm_dr


    ! decide whether we are using a table lookup for these interactions
    Select Case(grid_Tang_Toennies)
    Case("yes")
       if ( lambda > Tang_Toennies_max ) then
          dC6_C10_damp_scalar1 = 0d0
       else
          ! C6 is in 1st index, C8 in 2nd, C10 in 3rd, C12 in 4th
          index = ( n - 6 ) / 2 + 1
          ! note that dTang_Toennies_damp is derivative w.r.t. lambda, we want derivative w.r.t distance, so multiply by exponent
          dC6_C10_damp_scalar1 = exponent * dTang_Toennies_table(index,ceiling(lambda/Tang_Toennies_max*dble(Tang_Toennies_grid)))
       endif

    Case("no")
       ! note that dTang_Toennies_damp is derivative w.r.t. lambda, we want derivative w.r.t distance, so multiply by exponent
       dC6_C10_damp_scalar1 = exponent * dTang_Toennies_damp(lambda,n)
    End Select

  end function dC6_C10_damp_scalar1



  !********************************************************
  ! this function computes the second order derivative of the 
  ! Tang-Toennies damping function for R6,R8,R10,R12 dispersion interaction
  ! this is explicitly used in Feynmann_Hibbs correction
  !********************************************************
  function dC6_C10_damp_scalar2(atom_id1,atom_id2,norm_dr,n)
    use global_variables
    real*8 :: dC6_C10_damp_scalar2
    integer,intent(in):: atom_id1,atom_id2,n
    real*8,intent(in) :: norm_dr
    real*8 :: exponent, lambda
    integer :: m

    ! buckingham exponent is the second parameter in atype_lj_parameter array, use this for screening
    exponent = atype_lj_parameter(atom_id1,atom_id2,2)
    lambda = exponent * norm_dr

    m=n-1
    dC6_C10_damp_scalar2 =  -exponent * exp(-lambda  ) * lambda ** ( n + 1 ) / norm_dr / factorial(n) + exp(-lambda  ) * lambda ** ( n + 1 ) / norm_dr**2 / factorial(m)

  end function dC6_C10_damp_scalar2


  !********************************************************
  ! this function computes the third order derivative of the 
  ! Tang-Toennies damping function for R6,R8,R10,R12 dispersion interaction
  ! this is explicitly used in Feynmann_Hibbs correction
  !********************************************************
  function dC6_C10_damp_scalar3(atom_id1,atom_id2,norm_dr,n)
    use global_variables
    real*8 :: dC6_C10_damp_scalar3
    integer,intent(in):: atom_id1,atom_id2,n
    real*8,intent(in) :: norm_dr
    real*8 :: exponent, lambda
    integer :: m,l

    ! buckingham exponent is the second parameter in atype_lj_parameter array, use this for screening
    exponent = atype_lj_parameter(atom_id1,atom_id2,2)
    lambda = exponent * norm_dr

    m=n-1
    l=n-2
    dC6_C10_damp_scalar3 =  ( exponent**2 * exp(-lambda  ) * lambda ** ( n + 1 ) / norm_dr / factorial(n) &
         & - 2d0 *exponent * exp(-lambda  ) * lambda ** ( n + 1 ) / norm_dr**2 / factorial(m) &
         & + exp(-lambda  ) * lambda ** ( n + 1 ) / norm_dr**3 / factorial(l)  )

  end function dC6_C10_damp_scalar3




  !********************************************************
  ! this subroutine determines whether the input atoms i_atom, j_atom on
  ! i_mole are either drude oscillators, or have drude oscillators attached to them
  !
  ! this subroutine returns corresponding drude oscillator indexes, signs of charges,
  ! and flags.  
  !
  ! This subroutine is called by intra_pme_energy, intra_pme_force, intra_elec_energy
  ! and intra_elec_force subroutines, so see those subroutines for how i_drude, sign_chg_i, flag_drudei,
  ! etc. are used.
  !********************************************************
  subroutine drude_atom_indices(i_mole,n_atom,i_atom,j_atom,i_drude,j_drude,sign_chg_i,sign_chg_j,flag_drudei,flag_drudej, flag_same_atom)
    use global_variables
    integer,intent(in)::i_mole,i_atom,j_atom
    integer, intent(in), dimension(:) :: n_atom
    integer, intent(out):: i_drude,j_drude,sign_chg_i,sign_chg_j,flag_drudei,flag_drudej,flag_same_atom

    integer :: i_test, j_test

    flag_drudei=0
    flag_drudej=0
    flag_same_atom=0

    ! give i_drude, j_drude, sign_chg_i,sign_chg_j values so that code doesn't get mad at undefined
    ! values
    ! if these variables are not set in the if statements below, then the flags will be set such
    ! that they are not used
    i_drude=0;j_drude=0;sign_chg_i=0;sign_chg_j=0;


    ! determine whether i_atom,j_atom are atoms or shells, this is needed because only shell charge is used for intra-molecular 
    ! screened interactions, and this charge needs to be taken from the corresponding drude_oscillator

    ! make sure that j_atom is not the drude_oscillator attached to i_atom or vice-versa
    ! use flag_same_atom to check this

    ! use drude_atom_map, iff drude_atom_map(i_mole,i_atom) = i_atom, then this is not a drude oscillator

    if ( drude_atom_map(i_mole,i_atom ) .eq. i_atom ) then
       ! this is not a drude oscillator, see if this atom has a shell on it
       sign_chg_i=-1
       ! shell indices are always larger than atom indices, so loop from i_atom + 1
       if ( i_atom < n_atom(i_mole) ) then
          do i_test=i_atom+1, n_atom(i_mole)
             if ( drude_atom_map(i_mole,i_test ) .eq. i_atom ) then
                ! found drude oscillator
                i_drude = i_test
                flag_drudei = 1
             endif
          enddo
       endif
    else
       ! this is a drude oscillator
       i_drude = i_atom
       flag_drudei = 1
       sign_chg_i=1
       ! make sure this drude oscillator is not attached to j_atom
       if ( drude_atom_map(i_mole,i_atom ) .eq. j_atom ) then
          flag_same_atom=1
       endif
    endif


    if ( drude_atom_map(i_mole,j_atom ) .eq. j_atom ) then
       ! this is not a drude oscillator, see if this atom has a shell on it
       sign_chg_j=-1
       ! shell indices are always larger than atom indices, so loop from j_atom + 1
       if ( j_atom < n_atom(i_mole) ) then
          do j_test=j_atom+1, n_atom(i_mole)
             if ( drude_atom_map(i_mole,j_test ) .eq. j_atom ) then
                ! found drude oscillator
                j_drude = j_test
                flag_drudej = 1
             endif
          enddo
       endif
    else
       ! this is a drude oscillator
       j_drude = j_atom
       flag_drudej = 1
       sign_chg_j=1
       ! make sure this drude oscillator is not attached to i_atom
       if ( drude_atom_map(i_mole,j_atom ) .eq. i_atom ) then
          flag_same_atom=1
       endif
    endif


  end subroutine drude_atom_indices


  !****************************************************
  !
  !  IMPORTANT: This subroutine should always be called, even if this is not
  !             a Drude oscillator simulation, as the array drude_atom_map will be used
  !             in electrostatic screening functions for all atoms
  ! 
  !  Framework atoms are stored in locations from n_mole to tot_n_mole, and we do not
  !  use drude oscillators on the framework, but still initialize those parts of arrays
  !
  !  This subroutine initializes the drude oscillator data
  !  The tot_chg array is modified to reflect the drude oscillator charges
  !  drude_atom_map, stored as a global variable, is initialized
  !  which maps which atom each drude oscillator is bound to
  !  Also, n_atom_drude is also constructed, which stores how many total atoms and drudes
  !  are on each molecule (for looping over), compared to n_atom , which is just atoms
  !****************************************************
  subroutine initialize_drude_oscillators(tot_n_mole, n_mole, n_atom, n_atom_drude, tot_chg , chg_drude)
    use global_variables
    integer,intent(in) :: tot_n_mole, n_mole
    integer,dimension(:),intent(in) :: n_atom
    integer,dimension(:),intent(inout) :: n_atom_drude
    real*8,dimension(:,:),intent(inout) :: tot_chg
    real*8,dimension(:,:),intent(in) :: chg_drude

    real*8,parameter::small=1D-6

    integer :: i_mole, i_atom, index_drude

    ! first take care of framework atoms, (just loop over all atoms here), and 
    ! fill in trivial map for atoms
    do i_mole=1, tot_n_mole
       do i_atom=1,n_atom(i_mole)
          ! trivial map from atom to atom
          drude_atom_map(i_mole,i_atom) = i_atom
       enddo
    enddo

    ! now take care of any drude oscillators, note we use n_mole here, so that we're not
    ! looping over framework atoms, since these don't have oscillators

    do i_mole=1,n_mole
       ! this is the starting index for Drude oscillators for this molecule
       index_drude = n_atom(i_mole) + 1
       do i_atom=1,n_atom(i_mole)

          if ( abs( chg_drude(i_mole,i_atom) ) > small ) then
             ! there is a drude oscillator on this atom, index this drude oscillator
             ! in the next empty spot in drude_atom_map

             drude_atom_map(i_mole,index_drude) = i_atom

             ! fix tot_chg array for both atom, and drude
             tot_chg(i_mole,i_atom) = tot_chg(i_mole,i_atom) - chg_drude(i_mole,i_atom)
             tot_chg(i_mole,index_drude)=chg_drude(i_mole,i_atom)

             index_drude = index_drude + 1

             ! make sure arrays are big enough to support these drude oscillators
             if ( index_drude > MAX_N_ATOM ) then   
                write(*,*) " "
                write(*,*) " for Drude oscillator simulations, must have MAX_N_ATOM set to at least the"
                write(*,*) " maximum number of atoms and drude oscillators in the largest molecule."
                write(*,*) " Currently, this is not the case."
                write(*,*) " Please increase this value in the global variables file"
                stop 
             endif

          endif

       enddo
       ! total number of atoms and drude oscillators on this molecule
       n_atom_drude(i_mole) = index_drude - 1

    enddo


  end subroutine initialize_drude_oscillators




  !********************************************************
  ! this function is a Tang-Toennies screening function for coulombic interactions
  ! when it is turned on, and "1.0" when it is off.  Currently, it is meant to use exponents
  ! from buckingham potential as parameters.  For drude oscillators, it takes parameters from
  ! corresponding atoms
  ! 
  !  we use the map drude_atom_map, to link Drude oscillators
  ! to their corresponding atoms.  We also use this map trivially for atoms, and so it 
  ! should be constructed regardless of whether this is a drude oscillator simulation
  !*******************************************************
  function screen(i_mole,j_mole,i_atom,j_atom,norm_dr)
    use global_variables
    real*8           :: screen
    integer,intent(in):: i_mole, j_mole, i_atom, j_atom
    real*8,intent(in) :: norm_dr


    integer :: p_atom1, p_atom2, atom_id1,atom_id2
    real*8  :: lambda

    if(screen_type .eq. 0 ) then

       screen = 1D0

    else if (screen_type .eq. 1 ) then

       ! these atoms may or may not be drude oscillators, either way, we can use the drude_atom_map array
       ! which will index both drude oscillators as well as regular atoms

       ! parent atom for drude oscillator ( or just atom )
       p_atom1 = drude_atom_map(i_mole,i_atom)
       p_atom2 = drude_atom_map(j_mole,j_atom)

       atom_id1 = atom_index(i_mole,p_atom1)
       atom_id2 = atom_index(j_mole,p_atom2)


       ! buckingham exponent is the second parameter in atype_lj_parameter array, use this for screening
       lambda = atype_lj_parameter(atom_id1,atom_id2,2) * norm_dr
       screen= Tang_Toennies_damp(lambda,1)

    else
       stop "screen type not implemented"

    endif

  end function screen

  !********************************************************
  ! this function is the gradient of the above screening function
  ! see comments to function above
  !*******************************************************
  function d_screen(i_mole,j_mole,i_atom,j_atom,r_ij,r_mag)
    use global_variables
    real*8,dimension(3) :: d_screen
    integer,intent(in):: i_mole, j_mole, i_atom, j_atom
    real*8,dimension(3),intent(in) :: r_ij
    real*8, intent(in)     :: r_mag

    integer :: p_atom1, p_atom2,atom_id1,atom_id2
    real*8  :: lambda,exponent

    if(screen_type .eq. 0 ) then

       d_screen = 0D0

    else if (screen_type .eq. 1 ) then

       ! these atoms may or may not be drude oscillators, either way, we can use the drude_atom_map array
       ! which will index both drude oscillators as well as regular atoms

       ! parent atom for drude oscillator ( or just atom )
       p_atom1 = drude_atom_map(i_mole,i_atom)
       p_atom2 = drude_atom_map(j_mole,j_atom)

       atom_id1 = atom_index(i_mole,p_atom1)
       atom_id2 = atom_index(j_mole,p_atom2)

       ! buckingham exponent is the second parameter in atype_lj_parameter array, use this for screening
       exponent = atype_lj_parameter(atom_id1,atom_id2,2)
       lambda = exponent * r_mag

       d_screen = exponent * r_ij / r_mag * dTang_Toennies_damp(lambda,1)

    endif

  end function d_screen


  !*******************************************************
  ! this function is a screening function for intermolecular induced dipole 
  ! (drude oscillator) interactions
  ! for our N2 model, the center site has no drude oscillator, however the 
  ! code treats the whole molecule as either all atoms having drude oscillators
  ! or none having drude oscillators.  Therefore we have to consider the case with
  ! a zero charge drude oscillator, so don't divide by zero
  !*******************************************************
  function thole_screen(pol1,pol2,r_ij,thole)
    real*8::thole_screen
    real*8,intent(in)::r_ij,pol1,pol2,thole
    real*8::a
    real*8,parameter :: drude_chg_thresh=1d-8

    a=thole

    ! don't divide by zero
    if ( ( pol1 > drude_chg_thresh ) .and. ( pol2 > drude_chg_thresh ) ) then
       thole_screen=1.D0-(1.D0+(a*r_ij)/(2.D0*(pol1*pol2)**(1.D0/6.D0)))*exp(-a*r_ij/(pol1*pol2)**(1.D0/6.D0))
    else
       ! this interaction is with a zero charge, so will be zero, just set screen to 1
       thole_screen=1d0
    endif


  end function thole_screen


  function  d_thole_screen(pol1,pol2,xyz,thole)
    real*8,dimension(3)::d_thole_screen
    real*8,intent(in)::pol1,pol2,thole
    real*8,dimension(3),intent(in)::xyz
    real*8::r_ij,fac,Ex,a
    real*8,parameter :: drude_chg_thresh=1d-8

    a=thole

    ! don't divide by zero
    if ( ( pol1 > drude_chg_thresh ) .and. ( pol2 > drude_chg_thresh ) ) then

       r_ij=sqrt(dot_product(xyz,xyz))
       fac=a/(pol1*pol2)**(1.D0/6.D0)
       Ex=exp(-fac*r_ij)

       d_thole_screen=(xyz/r_ij)*(fac*(1.D0+fac*r_ij/2.D0)*Ex-(fac/2.D0)*Ex)

    else
       d_thole_screen=0d0
    endif


  end function d_thole_screen



  !********************************
  ! this subroutines generates velocity from
  ! a maxwell boltzmann distribution
  !
  ! notice that angular velocity can also be generated, just by
  ! inputing moment of inertia instead of mass
  !
  ! notice two velocities are generated and output per call
  ! assumes temp input in K, mass in g/mol
  ! velocity is output in A/ps
  ! ang velocity is output in ps^-1
  !********************************
  subroutine max_boltz(vel,mass,temp)
    implicit none
    real*8,dimension(2),intent(out)::vel
    real*8,intent(in)::mass,temp
    real*8::rsq
    real*8,parameter::small=1D-4

    ! this subroutine will be used for angular velociies, in which case moment of inertia is input
    ! instead of mass.  For an atom, moment of inertia is zero.  Watch out for this

    if ( mass > small ) then
       do
          call random_number(vel)
          vel=2.0*vel-1.0
          rsq=vel(1)**2+vel(2)**2
          if (rsq > 0.0 .and. rsq < 1.0) exit
       end do
       rsq=sqrt(-2.0*log(rsq)/rsq)
       vel=vel*rsq

       !now multiply by stdev of mb distribution
       vel=vel*sqrt(8.314*temp/(mass/1000.))
       !now convert to A/ps for velocity
       !if this is angular velocity, inertia(mass) is input in g/mol*Ang^2, so this converts to ps^-1
       vel=vel/100.

    else
       vel=0d0
    endif

  end subroutine max_boltz

  !***************************
  ! this computes the determinant of a 3x3 matrix
  !***************************
  real*8 function determinant_3(a)
    real*8,dimension(:,:),intent(in) :: a

    if ( ( size(a(:,1)) .ne. 3 ) .or. ( size(a(1,:)) .ne. 3  ) ) then
       stop "matrix input to function determinant_3 must be 3x3 ! "
    endif

    determinant_3 = a(1,1) * (a(2,2) * a(3,3) - a(3,2) * a(2,3) ) - a(1,2) * (a(2,1) * a(3,3) - a(2,3)*a(3,1) ) + a(1,3) * (a(2,1) * a(3,2) - a(3,1) * a(2,2) )

  end function determinant_3


  subroutine crossproduct( a,b,ans )
    implicit none
    real*8, intent(in) :: a(3),b(3)
    real*8, intent(out) :: ans(3)
    ans(1) = a(2)*b(3)-a(3)*b(2)
    ans(2) = -a(1)*b(3)+a(3)*b(1)
    ans(3) = a(1)*b(2)-a(2)*b(1)
  end subroutine crossproduct

  real function volume( a, b, c )
    implicit none
    real*8, intent(in) :: a(3), b(3), c(3)
    volume = a(1) * (b(2)*c(3)-b(3)*c(2)) - a(2) * (b(1)*c(3)-b(3)*c(1)) + a(3) * (b(1)*c(2)-b(2)*c(1))
    volume = abs(volume)
  end function volume


  !**************************************************
  ! subroutines for pme dispersion
  !**************************************************
  !**************************************************************************
  ! function searches g function grid, which is used in pme dispersion calc
  ! this assumes all arguments are positive
  !**************************************************************************
  function gfun_g(x,p)
    use global_variables
    real*8 :: gfun_g
    real*8, intent(in) :: x
    integer, intent(in) :: p

    Select case (grid_erfc)
    case("yes") 
       if ( x > gfun_max) then
          gfun_g = 0d0
       else
          select case(p)
          case(6)
             gfun_g = g6_table(ceiling(x/gfun_max*dble(gfun_grid)))
          case(8)
             gfun_g = g8_table(ceiling(x/gfun_max*dble(gfun_grid)))
          case(10) 
             gfun_g = g10_table(ceiling(x/gfun_max*dble(gfun_grid)))
          case(12) 
             gfun_g = g12_table(ceiling(x/gfun_max*dble(gfun_grid)))
          case(14) 
             gfun_g = g14_table(ceiling(x/gfun_max*dble(gfun_grid)))
          case default
             stop "gfun_g option not recognized"
          end select
       endif
    case("no")
       select case(p)
       case(6)
          gfun_g = (x**4/2d0+x**2+1d0)*dexp(-x**2)
       case(8)
          gfun_g = (x**6+3d0*x**4+6d0*x**2+6d0)/6d0*dexp(-x**2)
       case(10)
          gfun_g = (x**8+4d0*x**6+12d0*x**4+24d0*x**2+24d0)/24d0*dexp(-x**2)
       case(12)
          gfun_g = (x**10+5d0*x**8+20d0*x**6+60d0*x**4+120d0*x**2+120d0)/120d0*dexp(-x**2)
       case(14)
          gfun_g = (x**12+6d0*x**10+30d0*x**8+120d0*x**6+360d0*x**4+720d0*x**2+720d0)/720d0*dexp(-x**2)
       case default
          stop "gfun_g option not recognized"
       end select
    case default
       stop "gfun_g option not recognized"
    end select

  end function gfun_g

  function wfun(u)
    use global_variables
    implicit none
    real*8 :: wfun
    real*8, intent(in) :: u
    integer :: j, k, p

    p = lgrg_order
    if ( u <= -p .or. u >= p ) then
       wfun = 0.0d0
       return
    endif
    k = floor(u)
    j = -p
    wfun = 1.0d0
    do while ( j <= p-1 )
       if ( j /= k ) then
          wfun = wfun * dble(u+j-k)/dble(j-k)
       endif
       j = j + 1
    enddo

    return
  end function wfun



  !***************************************************
  ! this gauss-jordan rotation subroutine is copied from numerical recipes
  ! with a few changes made as to variable type definitions
  !***************************************************

  SUBROUTINE gaussj(a,b)
    IMPLICIT NONE
    REAL*8, DIMENSION(:,:), INTENT(INOUT) :: a,b
    INTEGER, DIMENSION(size(a,1)) :: ipiv,indxr,indxc
    LOGICAL, DIMENSION(size(a,1)) :: lpiv
    REAL*8 :: pivinv
    REAL*8, DIMENSION(size(a,1)) :: dumc
    INTEGER, TARGET :: irc(2)
    INTEGER :: i,l,n
    INTEGER, POINTER :: irow,icol

    if ( ( size(a,1) .eq. size(a,2) ) .and. ( size(a,1) .eq. size(b,1) ) ) then
       ! good
       n = size(a,1)
    else
       stop "check array sizes in gaussj subroutine"
    endif

    irow => irc(1)
    icol => irc(2)
    ipiv=0
    do i=1,n
       lpiv = (ipiv == 0)
       irc=maxloc(abs(a),outerand(lpiv,lpiv))
       ipiv(icol)=ipiv(icol)+1
       if (ipiv(icol) > 1) then
          stop "gaussj: singular matrix (1)"
       endif
       if (irow /= icol) then
          call swap2(a(irow,:),a(icol,:))
          call swap2(b(irow,:),b(icol,:))
       end if
       indxr(i)=irow
       indxc(i)=icol
       if (a(icol,icol) == 0.0) then
          stop "gaussj: singular matrix (2)"
       endif

       pivinv=1.0d0/a(icol,icol)
       a(icol,icol)=1.0
       a(icol,:)=a(icol,:)*pivinv
       b(icol,:)=b(icol,:)*pivinv
       dumc=a(:,icol)
       a(:,icol)=0.0
       a(icol,icol)=pivinv
       a(1:icol-1,:)=a(1:icol-1,:)-outerprod(dumc(1:icol-1),a(icol,:))
       b(1:icol-1,:)=b(1:icol-1,:)-outerprod(dumc(1:icol-1),b(icol,:))
       a(icol+1:,:)=a(icol+1:,:)-outerprod(dumc(icol+1:),a(icol,:))
       b(icol+1:,:)=b(icol+1:,:)-outerprod(dumc(icol+1:),b(icol,:))
    end do
    do l=n,1,-1
       call swap2(a(:,indxr(l)),a(:,indxc(l)))
    end do
  END SUBROUTINE gaussj



  !***************************************************
  ! this jacobi rotation subroutine is copied from numerical recipes
  ! with a few changes made as to variable type definitions
  !***************************************************
  SUBROUTINE jacobi(a,d,v,nrot)
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: nrot
    REAL*8, DIMENSION(:), INTENT(OUT) :: d
    REAL*8, DIMENSION(:,:), INTENT(INOUT) :: a
    REAL*8, DIMENSION(:,:), INTENT(OUT) :: v
    INTEGER :: i,ip,iq,n
    REAL*8 :: c,g,h,s,sm,t,tau,theta,tresh
    REAL*8, DIMENSION(size(d)) :: b,z


    ! check sizes
    if ( (size(a,1)==size(a,2) ) .and. (size(a,2)==size(d) ) .and. (size(d)==size(v,1) ) .and. (size(v,1)==size(v,2) ) ) then
       n=size(a,1)
    else
       stop "mismatch in array dimensions in jacobi"
    endif

    call unit_matrix(v(:,:))
    b(:)=get_diag(a(:,:))
    d(:)=b(:)
    z(:)=0.0
    nrot=0
    do i=1,50
       sm=sum(abs(a),mask=upper_triangle(n,n))
       if (sm == 0.0d0) RETURN
       tresh=merge(0.2d0 *sm/n**2,0.0d0, i < 4 )
       do ip=1,n-1
          do iq=ip+1,n
             g=100.0d0 *abs(a(ip,iq))
             if ((i > 4) .and. (abs(d(ip))+g == abs(d(ip))) &
                  .and. (abs(d(iq))+g == abs(d(iq)))) then
                a(ip,iq)=0.0
             else if (abs(a(ip,iq)) > tresh) then
                h=d(iq)-d(ip)
                if (abs(h)+g == abs(h)) then
                   t=a(ip,iq)/h
                else
                   theta=0.5d0 *h/a(ip,iq)
                   t=1.0d0/(abs(theta)+sqrt(1.0d0+theta**2))
                   if (theta < 0.0) t=-t
                end if
                c=1.0d0/sqrt(1+t**2)
                s=t*c
                tau=s/(1.0d0+c)
                h=t*a(ip,iq)
                z(ip)=z(ip)-h
                z(iq)=z(iq)+h
                d(ip)=d(ip)-h
                d(iq)=d(iq)+h
                a(ip,iq)=0.0
                call jrotate(a(1:ip-1,ip),a(1:ip-1,iq))
                call jrotate(a(ip,ip+1:iq-1),a(ip+1:iq-1,iq))
                call jrotate(a(ip,iq+1:n),a(iq,iq+1:n))
                call jrotate(v(:,ip),v(:,iq))
                nrot=nrot+1
             end if
          end do
       end do
       b(:)=b(:)+z(:)
       d(:)=b(:)
       z(:)=0.0
    end do

    stop 'too many iterations in jacobi'

  CONTAINS

    SUBROUTINE jrotate(a1,a2)
      REAL*8, DIMENSION(:), INTENT(INOUT) :: a1,a2
      REAL*8, DIMENSION(size(a1)) :: wk1
      wk1(:)=a1(:)
      a1(:)=a1(:)-s*(a2(:)+a1(:)*tau)
      a2(:)=a2(:)+s*(wk1(:)-a2(:)*tau)
    END SUBROUTINE jrotate

  END SUBROUTINE jacobi



  !************** some subroutines from numerical recipes, redefined here

  FUNCTION iminloc(arr)
    REAL*8, DIMENSION(:), INTENT(IN) :: arr
    INTEGER, DIMENSION(1) :: imin
    INTEGER :: iminloc
    imin=minloc(arr(:))
    iminloc=imin(1)
  END FUNCTION iminloc

  FUNCTION imaxloc(arr)
    REAL*8, DIMENSION(:), INTENT(IN) :: arr
    INTEGER :: imaxloc
    INTEGER, DIMENSION(1) :: imax
    imax=maxloc(arr(:))
    imaxloc=imax(1)
  END FUNCTION imaxloc

  SUBROUTINE swap(a,b)
    REAL*8, INTENT(INOUT) :: a,b
    REAL*8 :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap

  SUBROUTINE swap2(a,b)
    REAL*8,dimension(:), INTENT(INOUT) :: a,b
    REAL*8,dimension(size(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap2

  FUNCTION outerand(a,b)
    LOGICAL, DIMENSION(:), INTENT(IN) :: a,b
    LOGICAL, DIMENSION(size(a),size(b)) :: outerand
    outerand = spread(a,dim=2,ncopies=size(b)) .and. &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerand

  FUNCTION outerdiff_i(a,b)
    integer, DIMENSION(:), INTENT(IN) :: a,b
    integer, DIMENSION(size(a),size(b)) :: outerdiff_i
    outerdiff_i = spread(a,dim=2,ncopies=size(b)) - &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_i

  FUNCTION outerprod(a,b)
    REAL*8, DIMENSION(:), INTENT(IN) :: a,b
    REAL*8, DIMENSION(size(a),size(b)) :: outerprod
    outerprod = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerprod

  FUNCTION get_diag(mat)
    REAL*8, DIMENSION(:,:), INTENT(IN) :: mat
    REAL*8, DIMENSION(size(mat,1)) :: get_diag
    INTEGER :: j

    if ( size(mat,1) == size(mat,2) ) then
       j= size(mat,1)
    else
       stop "mismatch of array dimensions in get_diag"
    endif

    do j=1,size(mat,1)
       get_diag(j)=mat(j,j)
    end do
  END FUNCTION get_diag

  SUBROUTINE unit_matrix(mat)
    REAL*8, DIMENSION(:,:), INTENT(OUT) :: mat
    INTEGER :: i,n
    n=min(size(mat,1),size(mat,2))
    mat(:,:)=0.0d0
    do i=1,n
       mat(i,i)=1.0d0
    end do
  END SUBROUTINE unit_matrix

  FUNCTION upper_triangle(j,k,extra)
    INTEGER, INTENT(IN) :: j,k
    INTEGER, OPTIONAL, INTENT(IN) :: extra
    LOGICAL, DIMENSION(j,k) :: upper_triangle
    INTEGER :: n
    n=0
    if (present(extra)) n=extra
    upper_triangle=(outerdiff_i(arth_i(1,1,j),arth_i(1,1,k)) < n)
  END FUNCTION upper_triangle

  FUNCTION arth_i(first,increment,n)
    INTEGER,PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
    INTEGER, INTENT(IN) :: first,increment,n
    INTEGER, DIMENSION(n) :: arth_i
    INTEGER :: k,k2,temp
    if (n > 0) arth_i(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_i(k)=arth_i(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_i(k)=arth_i(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_i

  !*****************************************************







  !******************************************************
  !        Here are some string manipulation routines 
  !        written by Dr. George Benthien and taken from
  !        http://www.gbenthien.net/strings/index.html
  !******************************************************

  subroutine parse(str,delims,args,nargs)

    ! Parses the string 'str' into arguments args(1), ..., args(nargs) based on
    ! the delimiters contained in the string 'delims'. Preceding a delimiter in
    ! 'str' by a backslash (\) makes this particular instance not a delimiter.
    ! The integer output variable nargs contains the number of arguments found.

    character(len=*) :: str,delims
    character(len=len_trim(str)) :: strsav
    character(len=*),dimension(:) :: args
    integer,intent(out) ::nargs
    integer :: i,k,na,lenstr


    strsav=str
    call compact(str)
    na=size(args)
    do i=1,na
       args(i)=' '
    end do
    nargs=0
    lenstr=len_trim(str)
    if(lenstr==0) return
    k=0

    do
       if(len_trim(str) == 0) exit
       nargs=nargs+1
       call split(str,delims,args(nargs))
       call removebksl(args(nargs))
    end do
    str=strsav

  end subroutine parse



  subroutine compact(str)

    ! Converts multiple spaces and tabs to single spaces; deletes control characters;
    ! removes initial spaces.

    character(len=*):: str
    character(len=1):: ch
    character(len=len_trim(str)):: outstr
    integer :: i,k,ich,isp,lenstr

    str=adjustl(str)
    lenstr=len_trim(str)
    outstr=' '
    isp=0
    k=0

    do i=1,lenstr
       ch=str(i:i)
       ich=iachar(ch)

       select case(ich)

       case(9,32)     ! space or tab character
          if(isp==0) then
             k=k+1
             outstr(k:k)=' '
          end if
          isp=1

       case(33:)      ! not a space, quote, or control character
          k=k+1
          outstr(k:k)=ch
          isp=0

       end select

    end do

    str=adjustl(outstr)

  end subroutine compact


  subroutine split(str,delims,before,sep)

    ! Routine finds the first instance of a character from 'delims' in the
    ! the string 'str'. The characters before the found delimiter are
    ! output in 'before'. The characters after the found delimiter are
    ! output in 'str'. The optional output character 'sep' contains the 
    ! found delimiter. A delimiter in 'str' is treated like an ordinary 
    ! character if it is preceded by a backslash (\). If the backslash 
    ! character is desired in 'str', then precede it with another backslash.

    character(len=*) :: str,delims,before
    character,optional :: sep
    logical :: pres
    character :: ch,cha
    integer:: i,k,lenstr,ibsl,ipos,iposa

    pres=present(sep)
    str=adjustl(str)
    call compact(str)
    lenstr=len_trim(str)
    if(lenstr == 0) return        ! string str is empty
    k=0
    ibsl=0                        ! backslash initially inactive
    before=' '
    do i=1,lenstr
       ch=str(i:i)
       if(ibsl == 1) then          ! backslash active
          k=k+1
          before(k:k)=ch
          ibsl=0
          cycle
       end if
       if(ch == '\') then          ! backslash with backslash inactive
          k=k+1
          before(k:k)=ch
          ibsl=1
          cycle
       end if
       ipos=index(delims,ch)         
       if(ipos == 0) then          ! character is not a delimiter
          k=k+1
          before(k:k)=ch
          cycle
       end if
       if(ch /= ' ') then          ! character is a delimiter that is not a space
          str=str(i+1:)
          if(pres) sep=ch
          exit
       end if
       cha=str(i+1:i+1)            ! character is a space delimiter
       iposa=index(delims,cha)
       if(iposa > 0) then          ! next character is a delimiter
          str=str(i+2:)
          if(pres) sep=cha
          exit
       else
          str=str(i+1:)
          if(pres) sep=ch
          exit
       end if
    end do
    if(i >= lenstr) str=''
    str=adjustl(str)              ! remove initial spaces
    return

  end subroutine split

  !**********************************************************************

  subroutine removebksl(str)

    ! Removes backslash (\) characters. Double backslashes (\\) are replaced
    ! by a single backslash.

    character(len=*):: str
    character(len=1):: ch
    character(len=len_trim(str))::outstr
    integer :: i,k,ibsl,lenstr

    str=adjustl(str)
    lenstr=len_trim(str)
    outstr=' '
    k=0
    ibsl=0                        ! backslash initially inactive

    do i=1,lenstr
       ch=str(i:i)
       if(ibsl == 1) then          ! backslash active
          k=k+1
          outstr(k:k)=ch
          ibsl=0
          cycle
       end if
  if(ch == '\') then          ! backslash with backslash inactive
   ibsl=1
   cycle
  end if
  k=k+1
  outstr(k:k)=ch              ! non-backslash with backslash inactive
end do

str=adjustl(outstr)

end subroutine removebksl

!**********************************************************************

end module routines
