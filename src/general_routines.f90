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
  ! kind of simulation is being run
  !********************************************************
  subroutine sort_input_files( file_io_data )
    use global_variables
    type(file_io_data_type), intent(inout)  :: file_io_data

    call getarg( 1, file_io_data%ifile_gro )
    call getarg( 2, file_io_data%ifile_ffpmt )
    call getarg( 3, file_io_data%ifile_top )
    call getarg( 4, file_io_data%ifile_simpmt )
    call getarg( 5, file_io_data%ofile_traj )
    call getarg( 6, file_io_data%ofile_log )
    call getarg( 7, file_io_data%ofile_hop )

  end subroutine sort_input_files

  

  !**********************************************************************
  ! This subroutine checks if we are restarting a trajectory
  ! this is determined by the existence of both the velocity checkpoint file
  ! and the two coordinate and energy output files.  Also the last
  ! time step of the velocity checkpoint file has to match the last time step
  ! in the trajectory output file, otherwise we will print an error
  !**********************************************************************
  subroutine check_restart_trajectory( file_io_data )
    use global_variables
    type(file_io_data_type) :: file_io_data

    integer  :: traj_file, log_file, vel_file
    logical  :: ofile_traj_exist, ofile_log_exist, velocity_file_exist
    integer :: flag, flag_eof, nargs, i_step_traj, i_step_vel
    character(400) :: line
    character(40),dimension(10)  :: args
    character(10) :: junk

    Inquire( file=file_io_data%ofile_traj, exist=ofile_traj_exist )
    Inquire( file=file_io_data%ofile_log , exist=ofile_log_exist  )    
    Inquire( file=file_io_data%ifile_velocity, exist=velocity_file_exist )


    if ( ofile_traj_exist .and. ofile_log_exist .and. velocity_file_exist ) then

       ! local variables
       traj_file = file_io_data%ofile_traj_file_h
       log_file  = file_io_data%ofile_log_file_h
       vel_file  = file_io_data%ifile_velocity_file_h

       ! ***************** read ofile_traj ******************
       ! get last step from ofile_traj
       open( traj_file, file=file_io_data%ofile_traj, status='old' )
       ! read until last printed step
       do 
          call read_file_find_heading( traj_file, 'step' , flag, flag_eof )
          if ( flag_eof == -1 ) Exit  ! end of file
          backspace(traj_file)
          read(traj_file,'(A)') line
          call parse( line," ", args, nargs )
          if ( nargs /= 4 ) stop "error reading line arguments in trajectory file"
          read(args(2),'(I10)') i_step_traj
       enddo

       ! ****************** read velocity_file ***************
       ! get last step from ofile_traj
       open( vel_file, file=file_io_data%ifile_velocity, status='old' )
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
       open( log_file, file=file_io_data%ofile_log, status='old' )  
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
  ! this subroutine scans to the end of a .gro trajectory to restart a simulation
  !**********************************************************************
  subroutine scan_grofile_restart( traj_file,  n_old_trajectory )
    integer, intent(in) :: traj_file  ! file handle
    integer, intent(in) :: n_old_trajectory    

    integer :: flag, flag_eof, nargs, i_step_traj
    character(400) :: line
    character(40),dimension(10)  :: args
    
       ! read until the final step
       do 
          call read_file_find_heading( traj_file, 'step' , flag, flag_eof )
          if ( flag_eof == -1 ) Exit  ! end of file
          backspace(traj_file)
          read(traj_file,'(A)') line
          call parse( line," ", args, nargs )
          if ( nargs /= 4 ) stop "error reading line arguments in trajectory file"
          read(args(2),'(I10)') i_step_traj

          ! if final step, exit this loop
          if ( i_step_traj == n_old_trajectory ) Exit
       enddo

       ! backspace one more time in prep for coordinate read...
       backspace(traj_file)

  end subroutine scan_grofile_restart




  !**********************************************************************
  ! this subroutine reads velocities from the velocity checkpoint file to restart a simulation
  !**********************************************************************
  subroutine read_velocity_restart_checkpoint( vel_file, velocity_file, velocity, n_old_trajectory )
   character(*),intent(in) :: velocity_file
    ! file handle
    integer, intent(in) :: vel_file    
    real*8,dimension(:,:), intent(inout) :: velocity
    integer, intent(in) :: n_old_trajectory

    integer :: flag, flag_eof, i_step_traj, i_atom, i_mole_local, i_atom_local
    character(5) :: mname
    character(5) :: aname
    character(10) :: junk

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
       do i_atom=1, size( velocity(1,:) )
          read( vel_file, '(I5,2A5,I5,3F14.6)' ), i_mole_local, mname, aname, i_atom_local, velocity(:,i_atom)
       enddo

      ! leave file open.  We should be at the end of the trajectory file, and so we are ready to append
      ! the new trajectory output

  end subroutine read_velocity_restart_checkpoint



  !***********************************************************************
  ! this subroutine figures out the number of molecules in gro file to allocate
  ! molecule arrays
  !***********************************************************************
  subroutine read_gro_number_molecules( file_handle, grofile, n_mole )
   integer, intent(in) :: file_handle
   character(*), intent(in) :: grofile
   integer, intent(out)     :: n_mole

   integer :: total_atoms, i, i_mole, junk
   character(5) :: mname, aname
   character(100) :: line
   real*8  :: r_tmp(3)

   open( file_handle, file=grofile, status='old' )

   read( file_handle, * ) line
   read( file_handle, '(I)' ) total_atoms 

   n_mole = 0
   do i = 1, total_atoms
      read( file_handle, '(I5,2A5,I5,3F8.3)' ), i_mole, mname, aname, junk, r_tmp(1), r_tmp(2), r_tmp(3)
      call trim_end( aname )
      if ( i_mole /= n_mole ) then
          n_mole = n_mole + 1
      end if
   end do
  
  close( file_handle )

  end subroutine read_gro_number_molecules




  !***********************************************************************
  ! this subroutine reads a .gro file and stores the atomic names and coordinates
  ! this will allocate aname, xyz arrays, and attach pointers from atom_data structure
  !***********************************************************************
  subroutine read_gro( file_handle, system_data, molecule_data, atom_data )
   use global_variables
   integer, intent(in) :: file_handle
   type(system_data_type), intent(inout) :: system_data
   type(molecule_data_type), dimension(:), intent(inout) :: molecule_data
   type(atom_data_type) , intent(inout)  :: atom_data

   integer :: total_atoms
   integer :: i,j, i_mole, i_atom, i_mole_prev, i_start, junk, nargs, inputstatus
   character(len(molecule_data(1)%mname)) :: mname
   real*8, dimension(3) :: r_tmp
   real*8, dimension(3,3) :: box
   character(40),dimension(9)  :: args
   character(400)::line

   read( file_handle, * ) line
   read( file_handle, '(I)' ) total_atoms

   system_data%total_atoms = total_atoms  ! store total atoms in simulation

   ! Allocate arrays (Fortran column major) and setup pointers. 
   allocate( atom_data%xyz(3,total_atoms), atom_data%aname(total_atoms) )

   i_mole_prev = 0
   i_atom = 0
   do i = 1, total_atoms
        read( file_handle, '(I5,2A5,I5,3F8.3)' ), i_mole, mname, atom_data%aname(i), junk, r_tmp(1), r_tmp(2), r_tmp(3)
        call trim_end( atom_data%aname(i) )
        if ( i_mole /= i_mole_prev ) then
             if ( i_mole_prev > 0 ) then
                molecule_data(i_mole_prev)%n_atom = i_atom
                ! IMPORTANT:  Allocate atom_index to be 1 unit bigger than
                ! n_atom, in case we protonate molecule in ms-evb
                allocate( molecule_data(i_mole_prev)%atom_index( i_atom + 1 ) )
                ! now fill in atom_index mapping array
                i_start = i - i_atom
                do j = 1 , i_atom
                   molecule_data(i_mole_prev)%atom_index( j ) = i_start + j - 1
                enddo                
             end if
             i_mole_prev = i_mole
             i_atom = 0
             call trim_end( mname )
             molecule_data(i_mole_prev)%mname = mname
        end if
        i_atom = i_atom + 1

        ! Gro file has coordinates in nm. This code uses angstrom
        ! convert nm to angstrom...
        atom_data%xyz(:,i) = r_tmp(:) * 10d0
   end do
   
   ! now fill in last molecule
   molecule_data(i_mole)%n_atom = i_atom
   allocate( molecule_data(i_mole)%atom_index( i_atom + 1 ) ) ! as previous, allocate 1 unit bigger...
   i_start = total_atoms - i_atom
   do j = 1 , i_atom
       molecule_data(i_mole)%atom_index( j ) = i_start + j
   enddo
   call trim_end( mname )
   molecule_data(i_mole)%mname = mname


   ! now get box.  Initialize to zero in case orthorhombic..
   box=0d0

    ! now box
    Read(file_handle,'(A)',Iostat=inputstatus) line
    ! if its orthogonal, 3 arguments, if not 9 arguments
    call parse(line," ",args,nargs)

    Select Case(nargs)
    Case(3)
       read(args(1),*) box(1,1)
       read(args(2),*) box(2,2)
       read(args(3),*) box(3,3)
    Case(9)
       read(args(1),*) box(1,1)
       read(args(2),*) box(2,2)
       read(args(3),*) box(3,3)
       read(args(4),*) box(1,2)
       read(args(5),*) box(1,3)
       read(args(6),*) box(2,1)
       read(args(7),*) box(2,3)
       read(args(8),*) box(3,1)
       read(args(9),*) box(3,2)
    case default
       stop "error reading box in read_trajectory_snapshot subroutine"
    End Select


    ! convert to angstroms
    box(:,:) = box(:,:) * 10d0

    system_data%box=box

  close( file_handle )


  end subroutine read_gro



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
  function pos_com( xyz, n_atom, mass )
    implicit none
    real*8, intent(in), dimension(:,:) :: xyz
    real*8, intent(in), dimension(:) :: mass
    integer, intent(in) :: n_atom
    integer :: i_atom
    real*8, dimension(3) :: pos_com
    real*8 :: m_tot
    pos_com(:) = 0D0
    m_tot = 0D0

    do i_atom = 1, n_atom
       pos_com(:) = pos_com(:) + xyz(:,i_atom) * mass(i_atom)
       m_tot = m_tot + mass(i_atom)
    end do
    pos_com(:) = pos_com(:) / m_tot

  end function pos_com


  !*********************************************************
  ! This function update all COM positions
  !*********************************************************
  subroutine update_r_com( n_mole, molecule_data , atom_data )
    use global_variables
    implicit none
    integer, intent(in)                      :: n_mole
    type(molecule_data_type), dimension(:), intent(inout)  :: molecule_data
    type(atom_data_type), intent(inout)      :: atom_data

    !******** this is a local data structure with pointers that will be set
    ! to subarrays of atom_data arrays for the specific atoms in the molecule
    type(single_molecule_data_type) :: single_molecule_data

    integer :: i_mole

    do i_mole = 1, n_mole
       ! set pointers for this data structure to target molecule
       call return_molecule_block( single_molecule_data , molecule_data(i_mole)%n_atom, molecule_data(i_mole)%atom_index, atom_xyz=atom_data%xyz, atom_mass=atom_data%mass )
       molecule_data(i_mole)%r_com(:) = pos_com( single_molecule_data%xyz, molecule_data(i_mole)%n_atom, single_molecule_data%mass )
    end do

  call dissociate_single_molecule_data(single_molecule_data)
  end subroutine update_r_com



  !**************************************************
  ! this subroutine creates the transformation matrix
  ! from cartesian to box coordinates
  !*************************************************
  subroutine initialize_non_orth_transform ( box, xyz_to_box_transform )
    real*8, dimension(3,3),intent(in) :: box
    real*8, dimension(3,3),intent(out) :: xyz_to_box_transform  

    real*8, dimension(3,3) :: temp_xyz
    real*8, dimension(3,3) :: junk
    integer :: i,j

    junk=0d0 
    junk(1,1)=1d0;junk(2,2)=1d0;junk(3,3)=1d0
    ! transpose of box vectors
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
    vol = volume( box )
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
  subroutine create_scaled_direct_coordinates(xyz_scale, xyz, n_atom, kk, K)
    real*8,dimension(:,:),intent(inout):: xyz_scale
    real*8,dimension(:,:),intent(in) :: xyz
    integer, intent(in)              :: n_atom
    real*8,dimension(:,:),intent(in) :: kk
    integer, intent(in) :: K

    integer :: i_atom,l
    real*8,parameter :: small=1D-6

    xyz_scale=0d0
    do i_atom=1,n_atom
       do l=1,3
          xyz_scale(l,i_atom)=dble(K)*dot_product(kk(l,:),xyz(:,i_atom))
          ! if atoms are out of grid, shift them back in
          if (xyz_scale(l,i_atom)<0d0) then
             xyz_scale(l,i_atom)=xyz_scale(l,i_atom)+dble(K)
          else if(xyz_scale(l,i_atom)>= dble(K)) then
             xyz_scale(l,i_atom)=xyz_scale(l,i_atom)-dble(K)
          endif
          ! make sure scaled coordinates are not numerically equal to zero, otherwise this will screw up Q grid routine
          if ( abs(xyz_scale(l,i_atom)) < small ) then
             xyz_scale(l,i_atom) = small
          end if
       enddo
    enddo

  end subroutine create_scaled_direct_coordinates





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




  !*********************************
  ! this function returns single_molecule_data_type structure
  ! which sets pointers to an array block of atoms that constitutes a molecule
  ! This local datastructure of pointers should be destroyed after use
  !
  ! while these pointers are generally set to point to the atom_data structure,
  ! we allow flexibility so that pointers can be set to temporary local data
  ! structures.  
  !
  ! we only set pointers for data structures that are input.
  ! all input data structures should be atomic data structures
  ! that have size total_atoms
  !*********************************
  subroutine return_molecule_block( single_molecule_data, n_atom, atom_index, atom_xyz, atom_velocity, atom_force, atom_mass, atom_charge, atom_type_index, atom_name )
    use global_variables
    type(single_molecule_data_type), intent(inout) :: single_molecule_data
    integer, intent(in)                          :: n_atom
    integer, dimension(:), intent(in)            :: atom_index
    !**** these are input atomic data structures for which we set pointers:  Not
    ! all of these may be needed, so not all of them will be input
    real*8, dimension(:,:), intent(in),target, optional :: atom_xyz
    real*8, dimension(:,:), intent(in),target, optional :: atom_velocity
    real*8, dimension(:,:), intent(in),target, optional :: atom_force
    real*8, dimension(:), intent(in),target, optional :: atom_mass
    real*8, dimension(:), intent(in),target, optional :: atom_charge
    integer, dimension(:), intent(in),target, optional :: atom_type_index
    character(*), dimension(:), intent(in),target, optional :: atom_name

    integer :: low_index, high_index

    ! get index of first and last atoms of molecule in global atom array
    low_index = atom_index(1)
    high_index = atom_index(n_atom)

    call dissociate_single_molecule_data(single_molecule_data)

    ! set pointers to molecule block of input atomic data structures that are present
    ! atomic xyz data
    if ( present(atom_xyz) ) then
       single_molecule_data%xyz=>atom_xyz(:,low_index:high_index)
    endif
    ! atomic velocity data
    if ( present(atom_velocity) ) then
       single_molecule_data%velocity=>atom_velocity(:,low_index:high_index)
    endif
    ! atomic force data
    if ( present(atom_force) ) then
       single_molecule_data%force=>atom_force(:,low_index:high_index)
    endif
    ! atomic mass data
    if ( present(atom_mass) ) then
       single_molecule_data%mass=>atom_mass(low_index:high_index)
    endif
    ! atomic charge data
    if ( present(atom_charge) ) then
       single_molecule_data%charge=>atom_charge(low_index:high_index)
    endif
    ! atom type index data
    if ( present(atom_type_index) ) then
       single_molecule_data%atom_type_index=>atom_type_index(low_index:high_index)
    endif
    ! atomic name data
    if ( present(atom_name) ) then
       single_molecule_data%aname=>atom_name(low_index:high_index)
    endif


  end subroutine return_molecule_block



  subroutine dissociate_single_molecule_data(single_molecule_data)
    use global_variables
    type(single_molecule_data_type), intent(inout) :: single_molecule_data

   ! nullify all pointers
    nullify( single_molecule_data%xyz )
    nullify( single_molecule_data%velocity )
    nullify( single_molecule_data%force )
    nullify( single_molecule_data%mass )
    nullify( single_molecule_data%charge )
    nullify( single_molecule_data%atom_type_index )
    nullify( single_molecule_data%aname )


  end subroutine dissociate_single_molecule_data



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
  ! to the log output file
  !*************************************************************************
  subroutine print_simulation_info( log_file , system_data, integrator_data, verlet_list_data, PME_data )
    use global_variables
    integer, intent(in) :: log_file
    type(system_data_type) :: system_data
    type(integrator_data_type) :: integrator_data
    type(verlet_list_data_type) :: verlet_list_data
    type(PME_data_type)    :: PME_data

    write( log_file, *) "********** Simulation Parameters **********"
    write( log_file, *) "ensemble                        ", integrator_data%ensemble
    write( log_file, *) "real-space cutoff               ", real_space_cutoff
    write( log_file, *) "verlet cutoff                   ", verlet_list_data%verlet_cutoff
    write( log_file, *) "PME beta spline order           ", PME_data%spline_order
    write( log_file, *) "PME grid size                   ", PME_data%pme_grid
    write( log_file, *) "PME Gaussian width (alpha_sqrt) ", PME_data%alpha_sqrt
    write( log_file, *) "number of threads               ", n_threads
    write( log_file, *) "temperature                     ", system_data%temperature
    write( log_file, *) "total number of atoms           ", system_data%total_atoms
    write( log_file, *) "volume                          ", system_data%volume
    write( log_file, *) "time step (ps)                  ", integrator_data%delta_t
    ! print tabulated functions
    write( log_file, *)  ""
    write( log_file, *) "********** Tabulated Functions (numerical accuracy)  **********"
    write( log_file, *) "Bspline table size              ", PME_data%spline_grid
    write( log_file, *) "erfc    table size              ", PME_data%erfc_grid
    write( log_file, *) "erfc    max value               ", PME_data%erfc_max
    ! print ms-evb parameters 
    Select Case(ms_evb_simulation)
    Case("yes")
       write( log_file, *) ""
       write( log_file, *) "********** MS-EVB Settings **********"
       write( log_file, *) "evb_max_chain                   ", evb_max_chain
       write( log_file, *) "evb_reactive_pair_distance      ", evb_reactive_pair_distance
       write( log_file, *) ""
    End Select

  end subroutine print_simulation_info


  !**************************************************************************
  ! This subroutine print trajectory and energy info of step to trajectory file and log file
  !**************************************************************************
  subroutine print_step( traj_file, log_file, i_step, system_data, integrator_data, molecule_data, atom_data )
    use global_variables
    integer, intent(in) :: traj_file, log_file, i_step
    type(system_data_type) , intent(in)    :: system_data
    type(integrator_data_type) , intent(in) :: integrator_data
    type(molecule_data_type), dimension(:), intent(in)   :: molecule_data
    type(atom_data_type)   , intent(in)    :: atom_data

    real*8 :: time_step

   time_step = dble(i_step) * integrator_data%delta_t

   ! first print frame to trajectory .gro file
   call print_gro_file( traj_file, i_step , time_step , system_data , molecule_data , atom_data ) 

   write( log_file, * ) 'i_step , time(ps), potential energy (kJ/mol), kinetic energy (kJ/mol)'
   write( log_file, '(I9,F10.3,2E16.6)' ) i_step, time_step, system_data%potential_energy, system_data%kinetic_energy

    Select Case(ms_evb_simulation)
    Case("yes")
       write( log_file, * )  '------------------------------'
    case default   
       write( log_file, * ) 'Electrostatic ,   VDWs ,   Bond   ,   Angle  ,  Dihedral'
       write( log_file, '(5E16.6)' ) system_data%E_elec, system_data%E_vdw, system_data%E_bond, system_data%E_angle, system_data%E_dihedral
       write( log_file, * )  '------------------------------'
    end Select

  end subroutine print_step




  !*******************************************
  ! this subroutine print grofile format output
  !*******************************************
  subroutine print_gro_file( grofile_h, i_step , time_step , system_data , molecule_data , atom_data )
    use global_variables
    integer, intent(in) :: grofile_h, i_step
    real*8, intent(in)  :: time_step
    type(system_data_type), intent(in)   :: system_data
    type(molecule_data_type), dimension(:), intent(in) :: molecule_data
    type(atom_data_type) , intent(in)    :: atom_data
    
    !******** this is a local data structure with pointers that will be set
    ! to subarrays of atom_data arrays for the specific atoms in the molecule
    type(single_molecule_data_type) :: single_molecule_data

    real*8  :: xyz(3), box(3,3)
    integer :: i_mole, i_atom, i_count

    write( grofile_h, *) "step ", i_step, "time(ps)", time_step
    write( grofile_h, *) system_data%total_atoms
    
    i_count=1
    do i_mole =1 , system_data%n_mole

       ! set pointers for this data structure to target molecule
       call return_molecule_block( single_molecule_data , molecule_data(i_mole)%n_atom, molecule_data(i_mole)%atom_index, atom_xyz=atom_data%xyz, atom_name=atom_data%aname )

       do i_atom=1, molecule_data(i_mole)%n_atom
           ! print in nanometers for gro file
           xyz(:) = single_molecule_data%xyz(:,i_atom) / 10d0
           write( grofile_h, '(I5,2A5,I5,3F8.3)' ) i_mole , molecule_data(i_mole)%mname , single_molecule_data%aname(i_atom), i_count, xyz
           i_count = i_count + 1
       enddo

    enddo

    ! print box
    box = system_data%box / 10d0
    write( grofile_h, '(9F7.4)' ) box(1,1) , box(2,2) , box(3,3) , box(1,2) , box(1,3) , box(2,1) , box(2,3) , box(3,1) , box(3,2)

  call dissociate_single_molecule_data(single_molecule_data)

  end subroutine print_gro_file





  !*******************************************
  ! this subroutine prints atomic velocities to a velocity checkpoint file
  ! the primary purpose of this is for continuing a simulation
  !*******************************************
  subroutine print_velocities_checkpoint( velfile_h, i_step , delta_t, system_data , molecule_data , atom_data )
    use global_variables
    integer, intent(in)      :: velfile_h, i_step
    real*8,  intent(in)      :: delta_t
    type(system_data_type), intent(in)   :: system_data
    type(molecule_data_type), dimension(:), intent(in) :: molecule_data
    type(atom_data_type) , intent(in)    :: atom_data

    !******** this is a local data structure with pointers that will be set
    ! to subarrays of atom_data arrays for the specific atoms in the molecule
    type(single_molecule_data_type) :: single_molecule_data

    real*8  :: time_step
    integer :: i_mole, i_atom
 
    time_step = dble(i_step) * delta_t

    write( velfile_h, *) "step ", i_step

    do i_mole =1 , system_data%n_mole
       ! set pointers for this data structure to target molecule
       call return_molecule_block( single_molecule_data , molecule_data(i_mole)%n_atom, molecule_data(i_mole)%atom_index, atom_velocity=atom_data%velocity, atom_name=atom_data%aname )
       do i_atom=1, molecule_data(i_mole)%n_atom
           write( velfile_h, '(I5,2A5,I5,3F14.6)' ) i_mole , molecule_data(i_mole)%mname , single_molecule_data%aname(i_atom), i_atom, single_molecule_data%velocity(:,i_atom)
       enddo
    enddo

  call dissociate_single_molecule_data(single_molecule_data)

  end subroutine print_velocities_checkpoint



  !**************************************************************
  ! This subroutine checks to make sure that molecules are not broken up
  ! by an intra-molecular pbc shift.  It assumes that the atoms are listed
  ! according to relative location in the molecule, so that the shift
  ! is based on the preceding atom
  !**************************************************************
  subroutine fix_intra_molecular_shifts( n_mole, molecule_data , atom_data, box, xyz_to_box_transform  )
    use global_variables
    integer, intent(in) :: n_mole
    type(molecule_data_type),dimension(:), intent(in) :: molecule_data
    type(atom_data_type), intent(inout)  :: atom_data
    real*8, dimension(:,:), intent(in)   :: box, xyz_to_box_transform

    !******** this is a local data structure with pointers that will be set
    ! to subarrays of atom_data arrays for the specific atoms in the molecule
    type(single_molecule_data_type) :: single_molecule_data

    integer :: i_mole

    ! loop over molecules
    do i_mole = 1 , n_mole

       ! set pointers for this data structure to target molecule
       ! we will be just using xyz coordinates here
       ! note we are changing global atom_data%xyz data structure with pointer !
       call return_molecule_block( single_molecule_data , molecule_data(i_mole)%n_atom, molecule_data(i_mole)%atom_index, atom_xyz=atom_data%xyz )
       call make_molecule_whole( molecule_data(i_mole)%n_atom , single_molecule_data%xyz , box, xyz_to_box_transform )

    enddo

  call dissociate_single_molecule_data(single_molecule_data)

  end subroutine fix_intra_molecular_shifts


  

  !******************** 
  ! this subroutine removes any pbc shifts for atoms within a molecule
  ! input xyz array should be pointer to coordinates of particular molecule
  !********************
  subroutine make_molecule_whole( n_atom , xyz, box, xyz_to_box_transform )
    integer, intent(in) :: n_atom
    real*8, dimension(:,:), intent(inout) :: xyz
    real*8, dimension(:,:), intent(in) :: box, xyz_to_box_transform

    integer ::  i_atom, j_atom
    real*8,dimension(3) :: shift, drij
    real*8,parameter :: small=1D-6

       if ( n_atom > 1 ) then
       ! shift atoms in molecule with respect to the preceding atom
          do i_atom=2, n_atom
             j_atom = i_atom - 1
             shift(:) = pbc_shift( xyz(:, j_atom), xyz(:, i_atom) , box , xyz_to_box_transform )
             if ( ( abs( shift(1) ) > small ) .or. ( abs( shift(2)) > small ) .or. ( abs( shift(3) ) > small ) ) then
                drij = pbc_dr( xyz(:,j_atom), xyz(:,i_atom), shift(:) )
                xyz(:,i_atom) = xyz(:,j_atom) + drij(:)
             endif
          enddo
       endif

  end subroutine make_molecule_whole





  !***************************************************************
  ! This subroutine checks that cutoffs are set appropriately for minimum image
  ! convention, and are consistent with our method for finding images
  !***************************************************************
  subroutine check_cutoffs_box( real_space_cutoff, verlet_cutoff, box )
    implicit none
    real*8, intent(in)                 :: real_space_cutoff, verlet_cutoff
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
    if ( real_space_cutoff .ge. .5 * min( box(1,1), box(2,2), box(3,3) ) ) then
       stop "real_space_cutoff is too big for this box"
    endif

    ! check the verlet skin
       if ( verlet_cutoff .ge. .5 * min( box(1,1), box(2,2), box(3,3) ) ) then
          write(*,*)  "verlet skin width is too big for this box"
          write(*,*)  "please set verlet_cutoff to a value less than "
          write(*,*)  "half the box length "
          stop
       endif

       ! now make sure that the verlet cutoff is longer than the cutoffs for which verlet is used

       if ( verlet_cutoff .le. real_space_cutoff ) then
          stop "verlet_cutoff must be set greater than real_space_cutoff "
       elseif ( ( verlet_cutoff - real_space_cutoff ) .lt. warn_verlet ) then
          write(*,*) ""
          write(*,*) "WARNING:  verlet_cutoff is less than ", warn_verlet
          write(*,*) "Angstrom larger than real_space_cutoff.  Are you sure you want"
          write(*,*) "such a thin verlet skin thickness?"
          write(*,*) ""
       endif

  end subroutine check_cutoffs_box




  !*****************************************************
  ! this subroutine translates molecules back into box if it
  ! has been moved out
  !*****************************************************
 subroutine shift_molecules_into_box( n_mole, molecule_data , atom_data, box, xyz_to_box_transform  )
    use global_variables
    integer, intent(in) :: n_mole
    type(molecule_data_type), dimension(:), intent(inout) :: molecule_data
    type(atom_data_type), intent(inout)  :: atom_data
    real*8, dimension(:,:), intent(in)   :: box, xyz_to_box_transform

    !******** this is a local data structure with pointers that will be set
    ! to subarrays of atom_data arrays for the specific atoms in the molecule
    type(single_molecule_data_type) :: single_molecule_data

    real*8,dimension(3) :: dr_box,shift,temp_xyz
    integer:: i_mole, i_atom, i, j

    ! loop over molecules
    do i_mole = 1 , n_mole

       ! set pointers for this data structure to target molecule
       ! we will be just using xyz coordinates here
       ! note we are changing global atom_data%xyz data structure with pointer !
       call return_molecule_block( single_molecule_data , molecule_data(i_mole)%n_atom, molecule_data(i_mole)%atom_index, atom_xyz=atom_data%xyz )

       ! general box, transform coordinates to box vector
       dr_box = matmul ( xyz_to_box_transform , molecule_data(i_mole)%r_com )
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

       ! shift center of mass back into box
       molecule_data(i_mole)%r_com(:)=molecule_data(i_mole)%r_com(:)+temp_xyz(:)

       ! shift atom positions
       do i_atom = 1, molecule_data(i_mole)%n_atom
          single_molecule_data%xyz(:,i_atom)=single_molecule_data%xyz(:,i_atom) + temp_xyz(:)
       enddo

    enddo  ! end loop over molecules

  call dissociate_single_molecule_data(single_molecule_data)

  end subroutine shift_molecules_into_box




  !*********************************************
  ! this subroutine allocates the size of the verlet list
  ! we have a subroutine to do this, because the verlet list can take up
  ! a significant amount of memory if not allocated intelligently, and 
  ! the required size can change if the density of molecules in the simulation changes
  !  Verlet list size should be about 4*pi*Rverlet^3 * Natoms^2 / 6 * volume
  !*********************************************
  subroutine allocate_verlet_list( verlet_list_data, total_atoms, volume )
    use global_variables
    type(verlet_list_data_type), intent(inout) :: verlet_list_data
    integer, intent(in) :: total_atoms
    real*8,  intent(in) :: volume

    integer :: size_verlet
    real*8  :: pi

    ! here we are allocating array data structures in global variables
    if ( allocated ( verlet_list_data%neighbor_list ) ) then
       deallocate( verlet_list_data%neighbor_list )
    endif

    if ( allocated ( verlet_list_data%verlet_point ) ) then
       deallocate( verlet_list_data%verlet_point )
    endif

    ! size verlet point should be equal to the number of atoms plus 1, where the
    ! last entry signals end of list
    allocate( verlet_list_data%verlet_point(total_atoms+1) )

    pi = constants%pi
    ! note size_verlet is integer, this will evaluate to integer...
    size_verlet = floor(4d0 * pi * verlet_list_data%verlet_cutoff**3 * dble(total_atoms)**2 / 6d0 / volume)

    ! for huge box, value of temp could be zero if we just have gas phase dimer
    ! verlet list should be at least as big as number of molecules
    size_verlet = max( total_atoms , size_verlet )

    ! safe_verlet is factor that we multiply theoretically needed size of verlet list to be safe
    ! found in global_variables
    size_verlet = floor( dble(size_verlet) *  verlet_list_data%safe_verlet )

    allocate( verlet_list_data%neighbor_list(size_verlet) )


  end subroutine allocate_verlet_list



  !**********************************************
  ! this subroutine keeps track of atomic displacements since the
  ! last construction of the Verlet list
  ! when the displacements get too large, a flag is turned on to tell
  ! the program to reconstruct the verlet list
  !
  ! OUTPUT: flag_verlet_list = 0, nothing needs to be done
  !         flag_verlet_list = 1, verlet list needs to be updated, but not reallocated
  !         flag_verlet_list = 2, verlet list needs to be reallocated and updated
  !*********************************************
  subroutine update_verlet_displacements( total_atoms, xyz, verlet_list_data, box, xyz_to_box_transform, flag_verlet_list, flag_initialize )
    use global_variables
    integer,intent(in) :: total_atoms
    type(verlet_list_data_type), intent(inout) :: verlet_list_data
    real*8,dimension(:,:),intent(in) :: xyz
    real*8,dimension(:,:),intent(in) :: box, xyz_to_box_transform
    integer, intent(out) :: flag_verlet_list
    integer, intent(in), optional :: flag_initialize

    integer :: i_atom
    real*8, dimension(3) :: r_i_new, r_i_old, shift, dr
    real*8  :: max_d1 , max_d2, max_temp, norm_dr , verlet_skin

    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "verlet list check displacements started at", time
    endif
    !***********************************************************************!

    ! if initialization
    if ( present(flag_initialize) ) then
      if ( .not. allocated( verlet_list_data%verlet_xyz_store ) ) then
         ! allocate verlet arrays to store old positions, displacements
         allocate( verlet_list_data%verlet_xyz_store(3,total_atoms) )
         allocate( verlet_list_data%verlet_displacement_store(3,total_atoms) )
      endif      

      verlet_list_data%verlet_xyz_store = xyz
      verlet_list_data%verlet_displacement_store=0d0
      flag_verlet_list = 0
    else
          max_d1=0d0; max_d2=0d0
          ! add new displacements and check max displacement
          do i_atom=1,total_atoms
                r_i_old(:) = verlet_list_data%verlet_xyz_store(:,i_atom)
                r_i_new(:) = xyz( :, i_atom )

                shift = pbc_shift( r_i_old, r_i_new, box , xyz_to_box_transform )
                dr = pbc_dr( r_i_old, r_i_new, shift )

                verlet_list_data%verlet_displacement_store(:,i_atom) = verlet_list_data%verlet_displacement_store(:,i_atom) + dr(:)

                dr(:) = verlet_list_data%verlet_displacement_store(:,i_atom)
                ! check whether this is either 1st or 2nd biggest displacement
                norm_dr = sqrt( dot_product(dr, dr) )
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
          enddo

          verlet_list_data%verlet_xyz_store = xyz

          ! now see whether we need to update the verlet list
          verlet_skin = verlet_list_data%verlet_thresh * (verlet_list_data%verlet_cutoff - real_space_cutoff)

          if ( ( max_d1 + max_d2 ) > verlet_skin ) then
             flag_verlet_list = 1 
          else
             flag_verlet_list = 0
          endif

    end if

    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "verlet list check displacements finished at", time
    endif
    !***********************************************************************!

  end subroutine update_verlet_displacements




  !**********************************************
  ! this subroutine constructs the verlet neighbor list using a cell list
  ! 
  ! See "Computer Simulation of Liquids" by Allen and Tildesley, Oxford Science Publications 2009
  ! for details
  !
  ! We construct an array "verlet_neighbor_list" which holds a list of all the neighbor atoms
  ! the array "verlet_point" points to the index of verlet_neighbor_list for an appropriate atom,
  ! for instance the neighbors of molecule i_atom are found in indices of verlet_neighbor_list between
  ! verlet_point(i_atom) and verlet_point(i_atom+1)
  !
  ! note that this means we need to set the value of verlet_point(last_atom+1), so that we know the finish position for
  ! neighbors of last_atom
  !*********************************************
  subroutine construct_verlet_list(verlet_list_data, atom_data, molecule_data, total_atoms, box, xyz_to_box_transform) 
    use global_variables
    type(verlet_list_data_type), intent(inout) :: verlet_list_data
    type(atom_data_type), intent(in)           :: atom_data
    type(molecule_data_type),dimension(:), intent(in) :: molecule_data
    integer,intent(in)          :: total_atoms
    real*8,dimension(:,:),intent(in) :: box, xyz_to_box_transform

    integer :: ix, iy, iz, nx, ny, nz,i_end_atom
    integer :: dia, dib, dic, igrid1, igrid2, igrid3, ia, ib, ic
    real*8, dimension(3) :: rtmp

    real*8 :: rka, rkb, rkc
    integer :: i,i_mole, j_mole, n_mole, i_atom, j_atom, k_atom, verlet_index, a_index, last_atom, verlet_neighbor_list_size
    real*8,dimension(3) :: r_ij, shift
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


    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "verlet list construction started at", time
    endif
    !***********************************************************************!

    n_mole = size(molecule_data)

    nx = verlet_list_data%na_nslist ; ny = verlet_list_data%nb_nslist ; nz = verlet_list_data%nc_nslist ;
    ! hash won't work if these are not 2 digit
    if ( ( nx < 10 ) .or. ( ny < 10 ) .or. ( nz < 10  ) .or.  ( nx > 99 ) .or. ( ny > 99 ) .or. ( nz > 99  )) then
       stop " na_nslist, nb_nslist, nc_nslist must be between 10 and 100 "
    end if

    verlet_cutoff2 = verlet_list_data%verlet_cutoff**2
    verlet_neighbor_list_size = size(verlet_list_data%neighbor_list)
    verlet_index = 1

   ! first construct the cell-based neighbor list
   allocate( nslist_cell(total_atoms), cell_index_hash(total_atoms) , index_molecule(total_atoms), headlist_cell(nx,ny,nz), endlist_cell(nx,ny,nz)  )

    !   pointers all set to null
    headlist_cell = 0
    endlist_cell = 0
    nslist_cell = 0


    do i_mole=1, n_mole
       do i_atom=1, molecule_data(i_mole)%n_atom

             a_index = molecule_data(i_mole)%atom_index(i_atom)
             ! store molecule index, this will be used later to avoid putting intramolecular atoms on neighbor list
             index_molecule(a_index) = i_mole

             rtmp = matmul(xyz_to_box_transform(:,:), atom_data%xyz(:,a_index) )

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
       enddo
    enddo

    ! this is the number of adjacent cells to search
    rka = dsqrt(dot_product(xyz_to_box_transform(1,:),xyz_to_box_transform(1,:)))
    rkb = dsqrt(dot_product(xyz_to_box_transform(2,:),xyz_to_box_transform(2,:)))
    rkc = dsqrt(dot_product(xyz_to_box_transform(3,:),xyz_to_box_transform(3,:)))

    dia = floor(verlet_list_data%verlet_cutoff*rka*dble(verlet_list_data%na_nslist)) + 1
    dib = floor(verlet_list_data%verlet_cutoff*rkb*dble(verlet_list_data%nb_nslist)) + 1
    dic = floor(verlet_list_data%verlet_cutoff*rkc*dble(verlet_list_data%nc_nslist)) + 1

    ! note dia should be less than 1/2 * na_nslist, and similarly for dib, dic
    ! otherwise code below can loop over the same cell twice if dia = 1/2 na_nslist,
    ! and this will cause a bug. 
    call check_neighbor_list_grid_cutoff_consistency( dia, dib , dic , verlet_list_data%na_nslist, verlet_list_data%nb_nslist , verlet_list_data%nc_nslist, rka, rkb, rkc )


    ! now form verlet neighbor list using this cell list
    do i_mole =1,n_mole
       do i_atom =1, molecule_data(i_mole)%n_atom
          a_index = molecule_data(i_mole)%atom_index(i_atom)
          last_atom = a_index
          ! set verlet point for this atom
          verlet_list_data%verlet_point(a_index) = verlet_index

          ! get cell index for this atom, 2 means read
          call cell_index_hash_io (  cell_index_hash, a_index, ix,iy,iz, 2 )

             ! loop over cells within cutoff distance from first cell
             do ia = -dia, dia
                igrid1 = ix + ia
                igrid1 = igrid1 - floor(dble(igrid1-1)/dble(verlet_list_data%na_nslist))*verlet_list_data%na_nslist
                do ib = -dib, dib
                   igrid2 = iy + ib
                   igrid2 = igrid2 - floor(dble(igrid2-1)/dble(verlet_list_data%nb_nslist))*verlet_list_data%nb_nslist
                   do ic = -dic, dic
                      igrid3 = iz + ic
                      igrid3 = igrid3 - floor(dble(igrid3-1)/dble(verlet_list_data%nc_nslist))*verlet_list_data%nc_nslist

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
                               ! a_index and j_atom are both global atom indices
                               r_ij = atom_data%xyz(:,a_index) - atom_data%xyz(:,j_atom)

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
                                  verlet_list_data%neighbor_list( verlet_index ) = j_atom
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
       enddo
    enddo

    ! last atom shouldn't have any neighbors
    verlet_list_data%verlet_point(total_atoms+1) = verlet_index

    deallocate( nslist_cell, cell_index_hash, index_molecule, headlist_cell, endlist_cell )


    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "verlet list construction finished at", time
    endif
    !***********************************************************************!


  end subroutine construct_verlet_list





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
  !********************************************************
  function C6_C10_damp(atom_id1,atom_id2,norm_dr,n)
    use global_variables
    real*8            :: C6_C10_damp
    integer,intent(in):: atom_id1,atom_id2,n
    real*8,intent(in) :: norm_dr

    real*8 :: lambda
    integer :: index

    ! buckingham exponent is the second parameter in atype_vdw_parameter array, use this for screening
    lambda = atype_vdw_parameter(atom_id1,atom_id2,2) * norm_dr

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
  function grad_C6_C10_damp(atom_id1,atom_id2,r_vec,norm_dr,n)
    real*8,dimension(3) :: grad_C6_C10_damp
    integer,intent(in):: atom_id1,atom_id2,n
    real*8,dimension(3),intent(in) :: r_vec
    real*8,intent(in) :: norm_dr

    grad_C6_C10_damp = dC6_C10_damp(atom_id1,atom_id2,norm_dr,n) * r_vec / norm_dr 

  end function grad_C6_C10_damp



  !********************************************************
  ! this function returns the first order derivative of the 
  ! Tang-Toennies damping function for R6,R8,R10,R12 dispersion interaction
  ! this can either use an explicit evaluation or table loop up
  !********************************************************
  function dC6_C10_damp(atom_id1,atom_id2,norm_dr,n)
    use global_variables
    real*8 :: dC6_C10_damp
    integer,intent(in):: atom_id1,atom_id2,n
    real*8,intent(in) :: norm_dr
    real*8 :: exponent, lambda
    integer :: index



    ! buckingham exponent is the second parameter in atype_vdw_parameter array, use this for screening
    exponent = atype_vdw_parameter(atom_id1,atom_id2,2)
    lambda = exponent * norm_dr


    ! decide whether we are using a table lookup for these interactions
    Select Case(grid_Tang_Toennies)
    Case("yes")
       if ( lambda > Tang_Toennies_max ) then
          dC6_C10_damp = 0d0
       else
          ! C6 is in 1st index, C8 in 2nd, C10 in 3rd, C12 in 4th
          index = ( n - 6 ) / 2 + 1
          ! note that dTang_Toennies_damp is derivative w.r.t. lambda, we want derivative w.r.t distance, so multiply by exponent
          dC6_C10_damp = exponent * dTang_Toennies_table(index,ceiling(lambda/Tang_Toennies_max*dble(Tang_Toennies_grid)))
       endif

    Case("no")
       ! note that dTang_Toennies_damp is derivative w.r.t. lambda, we want derivative w.r.t distance, so multiply by exponent
       dC6_C10_damp = exponent * dTang_Toennies_damp(lambda,n)
    End Select

  end function dC6_C10_damp





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
  subroutine max_boltz(vel,mass,temperature,kB)
    implicit none
    real*8,dimension(2),intent(out)::vel
    real*8,intent(in)::mass,temperature,kB

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
       vel=vel*sqrt(kB*1000d0*temperature/(mass/1000d0))  ! kB is kJ/mol/K convert to J, convert mass to kg/mol
       !now convert from m/s to A/ps
       vel=vel/100d0

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

  real function volume( box )
    implicit none
    real*8, dimension(3,3), intent(in) :: box
    real*8    :: a(3), b(3), c(3)
   
    a(:) = box(1,:)
    b(:) = box(2,:)
    c(:) = box(3,:)   

    volume = a(1) * (b(2)*c(3)-b(3)*c(2)) - a(2) * (b(1)*c(3)-b(3)*c(1)) + a(3) * (b(1)*c(2)-b(2)*c(1))
    volume = abs(volume)
  end function volume



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
