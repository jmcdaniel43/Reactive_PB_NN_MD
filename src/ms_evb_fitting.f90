!*****************************************
! this module fits parameters for ms-evb potential
! to ab initio computed energies
!*****************************************

module ms_evb_fit
  use ms_evb
  use routines
  use global_variables
  implicit none
  real*8, dimension(:,:,:,:), allocatable :: xyz_configs
  real*8, dimension(:), allocatable :: energy_configs

  ! these arrays are for amoeba
  real*8, dimension(:,:), allocatable :: amoeba_param



contains

  subroutine fit_ms_evb_parameters( ifile_traj, log_file ,n_mole, n_atom, chg, dfti_desc, dfti_desc_inv , box )
    use total_energy_forces
    use mc_routines
    use MKL_DFTI
    character(*) , intent(in) :: ifile_traj
    integer, intent(in)  :: log_file , n_mole
    integer, dimension(:), intent(in) :: n_atom
    real*8, dimension(:,:), intent(in) :: chg
    real*8, dimension(:,:), intent(in)  :: box
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv

    integer :: i_config, n_config, tot_n_mole, i_param, n_param, n_param_total, n_param_junk, iter,i
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: force_atoms
    real*8  :: potential
    ! need mass, and r_com to input to ms_evb routines, also copy data to temporary data structures,
    ! as ms_evb will change data structures if there is a proton hop, and we don't want this
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM)  :: mass_temp, r_com_temp, chg_temp
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: xyz_temp
    integer, dimension(MAX_N_MOLE) :: n_atom_temp
    ! store atom_index and molecule_index arrays as ms_evb subroutine may change these with a proton hop, and we need to change them back
    integer, dimension(MAX_N_MOLE,MAX_N_ATOM)  :: atom_index_store
    integer, dimension(MAX_N_MOLE)  :: molecule_index_store

    ! test
    real*8 :: E_elec , E_elec_nopol , E_bh, E_3body, E_bond, E_angle, E_dihedral
    integer :: iteration
    integer, dimension(MAX_N_MOLE) :: n_atom_drude
    real*8 , dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: xyz_drude
    ! test

    ! amoeba routines
    real*8, dimension(:) , allocatable :: y
    real*8, parameter      :: ftol = 1D-6


    ! store atom and molecule index arrays
    atom_index_store = atom_index
    molecule_index_store = molecule_index


    tot_n_mole = n_mole
    ! read energies and coordinates from input trajectory
    call read_input_trajectory( ifile_traj, n_config, n_mole, n_atom, box )

    ! construct mass array
    call fill_mass( tot_n_mole, n_atom, mass_temp )

    do i_config=1, n_config
       call update_r_com( xyz_configs(i_config,:,:,:), n_mole, n_atom, mass_temp, r_com_temp )
       xyz_temp(:,:,:) = xyz_configs(i_config,:,:,:)
       chg_temp = chg
       n_atom_temp = n_atom
       call ms_evb_calculate_total_force_energy( force_atoms, potential, tot_n_mole, n_mole, n_atom_temp, xyz_temp, r_com_temp, chg_temp, mass_temp, box, dfti_desc,dfti_desc_inv,log_file )

       ! test
!!$       n_atom_drude = n_atom
!!$       xyz_drude = xyz_temp
!!$       call calculate_total_force_energy( force_atoms,potential, E_elec, E_elec_nopol, E_bh, E_3body, E_bond, E_angle, E_dihedral, iteration, tot_n_mole, n_mole, n_atom, n_atom_drude, r_com_temp, xyz_temp, chg_temp, box, dfti_desc,dfti_desc_inv,log_file,xyz_drude)
       ! test


       write(*,*) i_config, potential
       ! now we need to fix global data structures if they were altered by a proton hop
       atom_index = atom_index_store
       molecule_index = molecule_index_store
       call update_hydronium_molecule_index( n_mole, n_atom, atom_index )

    end do

stop


    ! figure out how many parameters we're fitting, here n_param is equal to the number of force field parameters
    call amoeba_parameter_mapping( n_param )
    ! now we add one more parameter, which fits a constant shift relative to the ab initio energies.
    !n_param_total = n_param+1
    n_param_total = n_param

    ! allocate amoeba arrays
    allocate( amoeba_param( n_param_total + 1 , n_param_total ) , y(n_param_total+1) )
    ! contruct intial guess for amoeba
    call construct_amoeba_fitting_parameters( amoeba_param )

    ! get energies of these initial parameters
    do i_param=1, n_param_total + 1
       write(*,*) "************ i_param ****************** "
       write(*,*) i_param

       y(i_param) = kaisq( amoeba_param(i_param,:), n_mole, n_atom , chg , mass_temp, box, dfti_desc , dfti_desc_inv, log_file )

!!$       write(*,*) amoeba_param(i_param,:)
!!$       do i=2,4
!!$          write(*,*) i
!!$          write(*,*) evb_proton_acceptor_parameters(i,:)
!!$       enddo

    enddo

!!$stop

    write(*,*) "starting amoeba...."
    call amoeba(amoeba_param,y,ftol,iter, n_mole, n_atom, chg , mass_temp, box, dfti_desc , dfti_desc_inv, log_file  )

    write(*,*) "finished amoeba...."
    write(*,*) "results from amoeba"
    do i_param=1, n_param_total + 1
       write(*,*) "parameter set ", i_param
       write(*,*) "kaisq ", y(i_param)
       write(*,*) "parameters "
       write(*,*) amoeba_param(i_param,:)
    enddo


    ! write energy
    do i_config=1, n_config
       call update_r_com( xyz_configs(i_config,:,:,:), n_mole, n_atom, mass_temp, r_com_temp )
       xyz_temp(:,:,:) = xyz_configs(i_config,:,:,:)
       chg_temp = chg
       n_atom_temp = n_atom
       call ms_evb_calculate_total_force_energy( force_atoms, potential, tot_n_mole, n_mole, n_atom_temp, xyz_temp, r_com_temp, chg_temp, mass_temp, box, dfti_desc,dfti_desc_inv,log_file )


       ! test
!!$       n_atom_drude = n_atom
!!$       xyz_drude = xyz_temp
!!$       call calculate_total_force_energy( force_atoms,potential, E_elec, E_elec_nopol, E_bh, E_3body, E_bond, E_angle, E_dihedral, iteration, tot_n_mole, n_mole, n_atom, n_atom_drude, r_com_temp, xyz_temp, chg_temp, box, dfti_desc,dfti_desc_inv,log_file,xyz_drude)
       ! test

       write(*,*) i_config, potential
       ! now we need to fix global data structures if they were altered by a proton hop
       atom_index = atom_index_store
       molecule_index = molecule_index_store
       call update_hydronium_molecule_index( n_mole, n_atom, atom_index )

    end do



  end subroutine fit_ms_evb_parameters




  !*****************************
  ! this function returns the kaisq
  ! for a particular parameter set
  ! based on the input configurations
  ! and energies, stored in the
  ! arrays xyz_configs and energy_configs
  !*****************************
  function kaisq( p, n_mole, n_atom , chg , mass, box, dfti_desc , dfti_desc_inv, log_file )
    use total_energy_forces
    use MKL_DFTI
    real*8  :: kaisq
    real*8, dimension(:), intent(inout) :: p
    integer, intent(in) :: n_mole
    integer, dimension(:), intent(in) :: n_atom
    real*8, dimension(:,:) , intent(in) :: chg, mass
    real*8, dimension(:,:) , intent(in) :: box
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv    
    integer, intent(in)   :: log_file

    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: force_atoms
    real*8  :: potential
    ! need mass, and r_com to input to ms_evb routines, also copy data to temporary data structures,
    ! as ms_evb will change data structures if there is a proton hop, and we don't want this
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM)  :: r_com_temp, chg_temp, mass_temp
    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: xyz_temp
    integer, dimension(MAX_N_MOLE) :: n_atom_temp
    ! store atom_index and molecule_index arrays as ms_evb subroutine may change these with a proton hop, and we need to change them back
    integer, dimension(MAX_N_MOLE,MAX_N_ATOM)  :: atom_index_store
    integer, dimension(MAX_N_MOLE)  :: molecule_index_store

    integer :: i_config, n_config, n_param, n_param_total
    real*8  :: e_shift

    ! test
    real*8 :: E_elec , E_elec_nopol , E_bh, E_3body, E_bond, E_angle, E_dihedral
    integer :: iteration
    integer, dimension(MAX_N_MOLE) :: n_atom_drude
    real*8 , dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: xyz_drude
    ! test

    ! map these parameters to global variables
    call amoeba_parameter_mapping( n_param, p , 0 ) 

    n_param_total=size(p)
    !e_shift = p(n_param_total)
    e_shift=0d0

    ! store atom and molecule index arrays
    atom_index_store = atom_index
    molecule_index_store = molecule_index

    n_config = size(xyz_configs(:,1,1,1))

    kaisq=0d0
    do i_config=1, n_config
       call update_r_com( xyz_configs(i_config,:,:,:), n_mole, n_atom, mass, r_com_temp )
       xyz_temp(:,:,:) = xyz_configs(i_config,:,:,:)
       chg_temp = chg
       mass_temp = mass
       n_atom_temp = n_atom
       call ms_evb_calculate_total_force_energy( force_atoms, potential, n_mole, n_mole, n_atom_temp, xyz_temp, r_com_temp, chg_temp, mass_temp, box, dfti_desc,dfti_desc_inv,log_file )


       ! test
!!$       n_atom_drude = n_atom
!!$       xyz_drude = xyz_temp
!!$       call calculate_total_force_energy( force_atoms,potential, E_elec, E_elec_nopol, E_bh, E_3body, E_bond, E_angle, E_dihedral, iteration, n_mole, n_mole, n_atom, n_atom_drude, r_com_temp, xyz_temp, chg_temp, box, dfti_desc,dfti_desc_inv,log_file,xyz_drude)
       ! test

!!$       write(*,*) i_config, potential

       kaisq = kaisq + weight( energy_configs(i_config) ) * ( potential - energy_configs(i_config) + e_shift ) ** 2
       ! now we need to fix global data structures if they were altered by a proton hop
       atom_index = atom_index_store
       molecule_index = molecule_index_store
       call update_hydronium_molecule_index( n_mole, n_atom, atom_index )
    end do


  contains

    !***************************************************
    ! this is a weight function (fermi-dirac distribution)
    ! that weights the contribution of each configuration to the
    ! kaisq fitting function based on its total energy
    !***************************************************
    function weight(energy)
      real*8::weight
      real*8,intent(in)::energy

      real*8,parameter::Eff_mu=2000d0
      real*8,parameter::Eff_kt=1d0


      weight = 1d0 / (exp((energy-Eff_mu)/Eff_kt)+1d0)

    end function weight


  end function kaisq




  !******************************
  ! this subroutine constructs an initial
  ! parameter array to input to amoeba.
  ! currently, this is done by using the 
  ! initial parameters, and changing parameters
  ! individually by X pct
  !********************************
  subroutine construct_amoeba_fitting_parameters( p2 )
    real*8, dimension(:,:) , intent(inout) :: p2

    integer :: n_param, n_param_total, i, j
    real*8,parameter  :: factor1=2d0 , factor2=5d0 , factor3=10d0
    call amoeba_parameter_mapping( n_param , p2(1,:) , 1 )

    ! last parameter is shift relative to ab initio energies.  Set it to something small, 1 KJ/mol
    n_param_total = size(p2(1,:))
    !p2(1,n_param_total) = 1d0

    ! first copy
    do i=1, n_param_total
       j=i+1
       p2(j,:) = p2(1,:)
    enddo

    ! now make outer limit, this is a bit tricky, do this intelligently
    do i=1, n_param_total
       j=i+1
       ! if negative, make up to 2 times its absolute value
       if ( p2(1,i) < 0d0 ) then          
          p2(j,i) = abs( p2(1,i)) * factor1
       else if( p2(1,i) < 1.2d0 ) then
          ! if small, say less than 1.2, the parameter is probably exponent, let it be 5 times its initial value
          p2(j,i) = p2(1,i) * factor2
       else
          ! 
          p2(j,i) = p2(1,i) * factor1
       endif
    enddo

  end subroutine construct_amoeba_fitting_parameters




  !**********************************
  ! this subroutine controls the mapping
  ! between amoeba fitting parameters
  ! and global variable data arrays
  !
  ! IMPORTANT :: This subroutine changes
  ! global variable data structures
  !**********************************
  subroutine amoeba_parameter_mapping( n_param , p0, flag_fill )
    integer, intent(out) :: n_param
    real*8,dimension(:), optional, intent(inout) :: p0
    integer, optional, intent(in) :: flag_fill


    ! these data structures are hard-coded flags specifying which parameters are being fit
    ! if the parameter is being fit, the array contains an index specifying the mapping from
    ! the amoeba parameter array.  Each parameter may be mapped to more than one global variable,
    ! this is done when specific index is contained multiple times within the arrays.

    integer, save, dimension(max_interaction_type,6)  :: evb_donor_acceptor_parameters_fitflag
    integer, save, dimension(max_interaction_type,5)  :: evb_proton_acceptor_parameters_fitflag
    integer, save, dimension(max_interaction_type,10) :: evb_diabat_coupling_parameters_fitflag
    integer, save, dimension(MAX_N_ATOM_TYPE)    :: evb_exchange_charge_atomic_fitflag
    integer, save, dimension(MAX_N_MOLE_TYPE,MAX_N_MOLE_TYPE)  :: evb_exchange_charge_proton_fitflag
    integer, save, dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE,2) :: evb_atype_lj_parameter_fitflag
    integer, save, dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE,2) :: evb_atype_angle_parameter_fitflag
    integer, save, dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE,3) :: evb_atype_dihedral_parameter_fitflag

    integer       :: i_atom, j_atom, k_atom, l_atom, i, index_temp
    integer, save :: index
    integer, save :: initialize=0

    Select Case( initialize )
    Case(0)
       ! intialize to 0
       evb_donor_acceptor_parameters_fitflag=0
       evb_proton_acceptor_parameters_fitflag=0
       evb_diabat_coupling_parameters_fitflag=0
       evb_exchange_charge_atomic_fitflag=0
       evb_exchange_charge_proton_fitflag=0
       evb_atype_lj_parameter_fitflag=0

       !***************** hard code in parameters to be fitted
       index=1
       !  OW     O_ah     H_a  donor_acceptor repulsion, don't fit rs, rc
!!$       do i=1, 4
!!$          evb_donor_acceptor_parameters_fitflag(2,i) = index
!!$          index = index + 1
!!$       enddo
!!$       !  O_b    O_h3o    H_h3o  donor_acceptor repulsion, don't fit rs, rc
!!$       do i=1, 4
!!$          evb_donor_acceptor_parameters_fitflag(3,i) = index
!!$          index = index + 1
!!$       enddo
       !  O_b     H_h3o  proton_acceptor repulsion, don't fit rs, rc
       do i=1,2
          evb_proton_acceptor_parameters_fitflag(2,i) = index
          index = index + 1
       enddo
       !  O_b     O_h3o  proton_acceptor repulsion, only fit prefactor, exponent
       do i=1,2
          evb_proton_acceptor_parameters_fitflag(3,i) = index
          index = index + 1
       enddo
!!$
!!$!******************************* fitting lj C6 parameter
!!$          evb_atype_lj_parameter_fitflag(11,8,2) = index
!!$          evb_atype_lj_parameter_fitflag(8,11,2) = index
!!$          index = index + 1
!!$
!!$!*******************************

!****************************** fitting angle and dihedral parameters
!!$          evb_atype_angle_parameter_fitflag(6,7,8,2) = index
!!$          evb_atype_angle_parameter_fitflag(8,7,6,2) = index
!!$          index = index + 1
!!$          evb_atype_angle_parameter_fitflag(8,7,8,2) = index
!!$          index = index + 1
!!$
!!$          evb_atype_dihedral_parameter_fitflag(7,6,8,8,2) = index
!!$          index = index + 1
!!$          evb_atype_dihedral_parameter_fitflag(7,8,8,8,2) = index
!!$          index = index + 1
!************************************************************

!!$       !  S_b     O_h3o  proton_acceptor repulsion, only fit prefactor, exponent
       do i=1,1
          evb_proton_acceptor_parameters_fitflag(4,i) = index
          index = index + 1
       enddo
!!$       do i=1,1
!!$          evb_proton_acceptor_parameters_fitflag(5,i) = index
!!$          index = index + 1
!!$       enddo
!!$       do i=1,1
!!$          evb_proton_acceptor_parameters_fitflag(6,i) = index
!!$          index = index + 1
!!$       enddo
!!$       !  OW      O_ah  proton_acceptor repulsion, only fit prefactor, exponent
!!$       do i=1,1
!!$          evb_proton_acceptor_parameters_fitflag(5,i) = index
!!$          evb_proton_acceptor_parameters_fitflag(6,i) = index
!!$          index = index + 1
!!$       enddo
!!$       do i=1,1    
!!$          evb_proton_acceptor_parameters_fitflag(7,i) = index
!!$          index = index + 1
!!$       enddo
!!$       do i=1,1    
!!$          evb_proton_acceptor_parameters_fitflag(8,i) = index
!!$          index = index + 1
!!$       enddo
!!$       do i=1,1     
!!$          evb_proton_acceptor_parameters_fitflag(7,i) = index
!!$          index = index + 1
!!$       enddo
       ! diabatic coupling, same parameters for  O_b    O_h3o    H_h3o and OW    O_ah     H_a
!!$       do i=1,4
!!$          evb_diabat_coupling_parameters_fitflag(2,i) = index
!!$          evb_diabat_coupling_parameters_fitflag(3,i) = index
!!$          index = index + 1
!!$       enddo
       index = index - 1
       initialize=1
    End Select

    ! output number of parameters
    n_param = index


    ! if parameter array "p" is present, we are either filling it with initial values, or filling global variables with new parameters
    if ( present(p0)) then
       if ( flag_fill > 0 ) then
          ! filling p with initial parameters
          !  OW     O_ah     H_a  donor_acceptor repulsion, don't fit rs, rc
          do i=1, 4
             if ( evb_donor_acceptor_parameters_fitflag(2,i) > 0 ) then
                index_temp = evb_donor_acceptor_parameters_fitflag(2,i)
                p0(index_temp) = evb_donor_acceptor_parameters(2,i)
             end if
          enddo
          !  O_b    O_h3o    H_h3o  donor_acceptor repulsion, don't fit rs, rc
          do i=1, 4
             if ( evb_donor_acceptor_parameters_fitflag(3,i) > 0 ) then
                index_temp = evb_donor_acceptor_parameters_fitflag(3,i)
                p0(index_temp) = evb_donor_acceptor_parameters(3,i)
             end if
          enddo
          !  O_b     H_h3o  proton_acceptor repulsion, don't fit rs, rc
          do i=1, 3
             if ( evb_proton_acceptor_parameters_fitflag(2,i) > 0 ) then
                index_temp = evb_proton_acceptor_parameters_fitflag(2,i)
                p0(index_temp) = evb_proton_acceptor_parameters(2,i)
             end if
          enddo
          !  OW      H_a   proton_acceptor repulsion, don't fit rs, rc
          do i=1, 3
             if ( evb_proton_acceptor_parameters_fitflag(3,i) > 0 ) then
                index_temp = evb_proton_acceptor_parameters_fitflag(3,i)
                p0(index_temp) = evb_proton_acceptor_parameters(3,i)
             end if
          enddo
          do i=1, 3
             if ( evb_proton_acceptor_parameters_fitflag(4,i) > 0 ) then
                index_temp = evb_proton_acceptor_parameters_fitflag(4,i)
                p0(index_temp) = evb_proton_acceptor_parameters(4,i)
             end if
          enddo
          do i=1, 3
             if ( evb_proton_acceptor_parameters_fitflag(5,i) > 0 ) then
                index_temp = evb_proton_acceptor_parameters_fitflag(5,i)
                p0(index_temp) = evb_proton_acceptor_parameters(5,i)
             end if
          enddo
          do i=1, 3
             if ( evb_proton_acceptor_parameters_fitflag(7,i) > 0 ) then
                index_temp = evb_proton_acceptor_parameters_fitflag(7,i)
                p0(index_temp) = evb_proton_acceptor_parameters(7,i)
             end if
          enddo
          do i=1, 3
             if ( evb_proton_acceptor_parameters_fitflag(8,i) > 0 ) then
                index_temp = evb_proton_acceptor_parameters_fitflag(8,i)
                p0(index_temp) = evb_proton_acceptor_parameters(8,i)
             end if
          enddo
          !  diabatic coupling
          do i=1, 10
             if ( evb_diabat_coupling_parameters_fitflag(2,i) > 0 ) then
                index_temp = evb_diabat_coupling_parameters_fitflag(2,i)
                p0(index_temp) = evb_diabat_coupling_parameters(2,i)
             end if
          enddo

          !*************************** lj parameters
          do i_atom = 1, MAX_N_ATOM_TYPE
             do j_atom = 1 , MAX_N_ATOM_TYPE
                if ( evb_atype_lj_parameter_fitflag(i_atom,j_atom,2) > 0 ) then
                  index_temp = evb_atype_lj_parameter_fitflag(i_atom,j_atom,2)
                   p0(index_temp) = atype_lj_parameter(i_atom,j_atom,2)
                endif
             enddo
         enddo
         ! ******************** angle
         do i_atom =1 , MAX_N_ATOM_TYPE
            do j_atom =1 , MAX_N_ATOM_TYPE
               do k_atom = 1 , MAX_N_ATOM_TYPE
                  if ( evb_atype_angle_parameter_fitflag(i_atom,j_atom,k_atom,2) > 0 ) then
                     index_temp = evb_atype_angle_parameter_fitflag(i_atom,j_atom,k_atom,2)
                     p0(index_temp) = atype_angle_parameter(i_atom,j_atom,k_atom,2)
                  endif
               enddo
            enddo
        enddo

        ! ******************** dihedral
         do i_atom =1 , MAX_N_ATOM_TYPE
            do j_atom =1 , MAX_N_ATOM_TYPE
               do k_atom = 1 , MAX_N_ATOM_TYPE
                  do l_atom = 1 , MAX_N_ATOM_TYPE
                     if ( evb_atype_dihedral_parameter_fitflag(i_atom,j_atom,k_atom,l_atom,2) > 0 ) then
                        index_temp = evb_atype_dihedral_parameter_fitflag(i_atom,j_atom,k_atom,l_atom,2)
                        p0(index_temp) = atype_dihedral_parameter(i_atom,j_atom,k_atom,l_atom,2)
                     end if
                  enddo
               enddo
            enddo
        enddo
         
       else
          ! filling global data structures with parameter values
          !  OW     O_ah     H_a  donor_acceptor repulsion, don't fit rs, rc
          do i=1, 4
             if ( evb_donor_acceptor_parameters_fitflag(2,i) > 0 ) then
                index_temp = evb_donor_acceptor_parameters_fitflag(2,i)
                ! only coefficient can be negative
                if ( i > 1 ) then
                   evb_donor_acceptor_parameters(2,i) = abs(p0(index_temp))
                else
                   evb_donor_acceptor_parameters(2,i) = p0(index_temp) 
                endif
             end if
          enddo
          !  O_b    O_h3o    H_h3o  donor_acceptor repulsion, don't fit rs, rc
          do i=1, 4
             if ( evb_donor_acceptor_parameters_fitflag(3,i) > 0 ) then
                index_temp = evb_donor_acceptor_parameters_fitflag(3,i)
                ! only coefficient can be negative
                if ( i > 1 ) then
                   evb_donor_acceptor_parameters(3,i) = abs(p0(index_temp))
                else
                   evb_donor_acceptor_parameters(3,i) = p0(index_temp)
                endif
             end if
          enddo
          !  O_b     H_h3o  proton_acceptor repulsion, don't fit rs, rc
          do i=1, 3
             if ( evb_proton_acceptor_parameters_fitflag(2,i) > 0 ) then
                index_temp = evb_proton_acceptor_parameters_fitflag(2,i)
                ! only coefficient can be negative
                if ( i > 1 ) then
                   evb_proton_acceptor_parameters(2,i) = abs(p0(index_temp)) 
                else
                   evb_proton_acceptor_parameters(2,i) = p0(index_temp) 
                end if
             end if
          enddo
          !  OW      H_a   proton_acceptor repulsion, don't fit rs, rc
          do i=1, 3
             if ( evb_proton_acceptor_parameters_fitflag(3,i) > 0 ) then
                index_temp = evb_proton_acceptor_parameters_fitflag(3,i)

                ! only coefficient can be negative
                if ( i > 1 ) then
                   evb_proton_acceptor_parameters(3,i) = abs(p0(index_temp))
                else
                   evb_proton_acceptor_parameters(3,i) = p0(index_temp)
                end if
             end if
          enddo
          do i=1, 3
             if ( evb_proton_acceptor_parameters_fitflag(4,i) > 0 ) then
                index_temp = evb_proton_acceptor_parameters_fitflag(4,i)

                ! only coefficient can be negative
                if ( i > 1 ) then
                   evb_proton_acceptor_parameters(4,i) = abs(p0(index_temp))
                else
                   evb_proton_acceptor_parameters(4,i) = p0(index_temp)
                end if
             end if
          enddo
          do i=1, 3
             if ( evb_proton_acceptor_parameters_fitflag(5,i) > 0 ) then
                index_temp = evb_proton_acceptor_parameters_fitflag(5,i)
                ! only coefficient can be negative
                if ( i > 1 ) then
                   evb_proton_acceptor_parameters(5,i) = abs(p0(index_temp))
                   evb_proton_acceptor_parameters(6,i) = abs(p0(index_temp))
                else
                   evb_proton_acceptor_parameters(5,i) = p0(index_temp)
                   evb_proton_acceptor_parameters(6,i) = p0(index_temp)
                end if
             end if
          enddo
          do i=1, 3
             if ( evb_proton_acceptor_parameters_fitflag(7,i) > 0 ) then
                index_temp = evb_proton_acceptor_parameters_fitflag(7,i)
                ! only coefficient can be negative
                if ( i > 1 ) then
                   evb_proton_acceptor_parameters(7,i) = abs(p0(index_temp))
                else
                   evb_proton_acceptor_parameters(7,i) = p0(index_temp)
                end if
             end if
          enddo
          do i=1, 3
             if ( evb_proton_acceptor_parameters_fitflag(8,i) > 0 ) then
                index_temp = evb_proton_acceptor_parameters_fitflag(8,i)
                ! only coefficient can be negative
                if ( i > 1 ) then
                   evb_proton_acceptor_parameters(8,i) = abs(p0(index_temp))
                else
                   evb_proton_acceptor_parameters(8,i) = p0(index_temp)
                end if
             end if
          enddo
          !  diabatic coupling
          do i=1, 10
             if ( evb_diabat_coupling_parameters_fitflag(2,i) > 0 ) then
                index_temp = evb_diabat_coupling_parameters_fitflag(2,i)
                ! only coefficient can be negative
                if ( i > 1 ) then
                   evb_diabat_coupling_parameters(2,i) = abs(p0(index_temp)) 
                   evb_diabat_coupling_parameters(3,i) = abs(p0(index_temp)) 
                else
                   evb_diabat_coupling_parameters(2,i) = p0(index_temp) 
                   evb_diabat_coupling_parameters(3,i) = p0(index_temp) 
                endif
             end if
          enddo

          !*************************** lj parameters
          do i_atom = 1, MAX_N_ATOM_TYPE
             do j_atom = 1 , MAX_N_ATOM_TYPE
                if ( evb_atype_lj_parameter_fitflag(i_atom,j_atom,2) > 0 ) then
                  index_temp = evb_atype_lj_parameter_fitflag(i_atom,j_atom,2)
                  atype_lj_parameter(i_atom,j_atom,2) = abs(p0(index_temp))
                endif
             enddo
         enddo

         ! ******************** angle
         do i_atom =1 , MAX_N_ATOM_TYPE
            do j_atom =1 , MAX_N_ATOM_TYPE
               do k_atom = 1 , MAX_N_ATOM_TYPE
                  if ( evb_atype_angle_parameter_fitflag(i_atom,j_atom,k_atom,2) > 0 ) then
                     index_temp = evb_atype_angle_parameter_fitflag(i_atom,j_atom,k_atom,2)
                     atype_angle_parameter(i_atom,j_atom,k_atom,2) = abs(p0(index_temp))
                  endif
               enddo
            enddo
        enddo

        ! ******************** dihedral
         do i_atom =1 , MAX_N_ATOM_TYPE
            do j_atom =1 , MAX_N_ATOM_TYPE
               do k_atom = 1 , MAX_N_ATOM_TYPE
                  do l_atom = 1 , MAX_N_ATOM_TYPE
                     if ( evb_atype_dihedral_parameter_fitflag(i_atom,j_atom,k_atom,l_atom,2) > 0 ) then
                        index_temp = evb_atype_dihedral_parameter_fitflag(i_atom,j_atom,k_atom,l_atom,2)
                        atype_dihedral_parameter(i_atom,j_atom,k_atom,l_atom,2) = abs(p0(index_temp))
                     end if
                  enddo
               enddo
            enddo
        enddo

       endif
    end if


  end subroutine amoeba_parameter_mapping




  !********************************
  ! this subroutine reads in a fitting input file,
  ! which contains energies and configurations for fitting
  !********************************
  subroutine read_input_trajectory( ifile_traj , n_config, n_mole, n_atom, box )
    character(*), intent(in) :: ifile_traj
    integer, intent(out) :: n_config
    integer, intent(in) :: n_mole
    integer, dimension(:), intent(in) :: n_atom
    real*8, dimension(:,:) , intent(in) :: box

    integer :: file_h=16, count, atoms, ios, i_mole, i_atom, i_config, i ,j
    real*8  :: box1 , box2 , box3
    real*8, parameter :: small=1D-6
    character(MAX_MNAME) :: mname
    character(MAX_ANAME) :: aname
    character(200) :: line
    character(20)  :: junk, string

    open (unit=file_h,file=ifile_traj,status="old")
    ! find number of configs


    count=0
    do 
       ! energy
       read( file_h, '(A)', iostat=ios ) line
       if ( ios < 0 ) exit ! end of the file
       count = count + 1
       ! number of atoms
       read( file_h, * ) atoms
       if ( atoms /= total_atoms ) stop "total number of atoms is inconsistent in subroutine read_input_trajectory"
       ! coordinates
       do i_mole=1,n_mole
          do i_atom=1, n_atom(i_mole)
             read( file_h, '(A)' ) line 
          enddo
       enddo
       ! box
       read( file_h, '(A)' ) line 
    enddo
    close(file_h)

    n_config = count

    allocate( xyz_configs(n_config, MAX_N_MOLE, MAX_N_ATOM, 3), energy_configs(n_config) )

    ! now get trajectory info
    open (unit=file_h,file=ifile_traj,status="old")

    do i_config=1, n_config
       ! energy
       read( file_h, *  ) junk, string
       read( string, * ) energy_configs(i_config)
       read( file_h, '(I5)' ) atoms

       ! coordinates
       do i_mole=1,n_mole
          do i_atom=1, n_atom(i_mole)
             read( file_h, '(I5,2A5,I5,3F14.6)' ), i, mname, aname, j, xyz_configs(i_config,i_mole,i_atom,:)
             call trim_end( mname )
             call trim_end( aname )
             ! check data
             if ( ( i /= i_mole ) .or. ( mname .ne. molecule_type_name(molecule_index(i_mole)) ) .or. ( aname .ne. atype_name( atom_index( i_mole, i_atom ) ) ) ) then
                write(*,*) "input data inconsistency in subroutine read_input_trajectory"
                stop
             end if
          enddo
       enddo
       ! box
       read( file_h, * ) box1 , box2, box3
       if (  ( abs( box1 - box(1,1) ) > small ) .or. ( abs( box2 - box(2,2) ) > small ) .or. ( abs( box3 - box(3,3) ) > small ) ) stop "inconsistent box vectors in subroutine read_input_trajectory"         

    enddo
    close(file_h)

  end subroutine read_input_trajectory




  !************************************************
  ! This is the amoeba fitting routine modified from
  ! numerical recipes
  !***********************************************

  SUBROUTINE amoeba(p,y,ftol,iter, n_mole, n_atom, chg, mass, box, dfti_desc , dfti_desc_inv, log_file)
    use MKL_DFTI
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: iter
    REAL*8, INTENT(IN) :: ftol
    REAL*8, DIMENSION(:), INTENT(INOUT) :: y
    REAL*8, DIMENSION(:,:), INTENT(INOUT) :: p
    integer, intent(in) :: n_mole
    integer, dimension(:), intent(in) :: n_atom
    real*8, dimension(:,:) , intent(in) :: chg, mass
    real*8, dimension(:,:) , intent(in) :: box
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv    
    integer, intent(in)   :: log_file


    INTEGER, PARAMETER :: ITMAX=50000
    REAL*8, PARAMETER :: TINY=1.0e-10
    INTEGER :: ihi,ndim
    REAL*8, DIMENSION(size(p,2)) :: psum
    call amoeba_private
  CONTAINS
    !BL
    SUBROUTINE amoeba_private
      IMPLICIT NONE
      INTEGER :: i,ilo,inhi
      REAL*8 :: rtol,ysave,ytry,ytmp
      ndim=size(p,2)
      if ( (size(p,2) .ne. size(p,1)-1 ) .or. (size(p,2) .ne. size(y)-1) ) then
         stop "error in input array sizes in amoeba"
      endif


      iter=0
      psum(:)=sum(p(:,:),dim=1)
      do
         write(*,*) "iteration ",iter," in amoeba..."
         ilo=iminloc(y(:))
         ihi=imaxloc(y(:))
         ytmp=y(ihi)
         y(ihi)=y(ilo)
         inhi=imaxloc(y(:))
         y(ihi)=ytmp
         rtol=2.0d0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
         if (rtol < ftol) then
            call swap(y(1),y(ilo))
            call swap2(p(1,:),p(ilo,:))
            RETURN
         end if
         if (iter >= ITMAX) then
            stop "ITMAX exceeded in amoeba"
         endif
         ytry=amotry(-1.0d0)
         iter=iter+1
         if (ytry <= y(ilo)) then
            ytry=amotry(2.0d0)
            iter=iter+1
         else if (ytry >= y(inhi)) then
            ysave=y(ihi)
            ytry=amotry(0.5d0)
            iter=iter+1
            if (ytry >= ysave) then
               p(:,:)=0.5d0*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
               do i=1,ndim+1
                  if (i /= ilo) y(i)= kaisq( p(i,:), n_mole, n_atom , chg , mass, box, dfti_desc , dfti_desc_inv, log_file )
               end do
               iter=iter+ndim
               psum(:)=sum(p(:,:),dim=1)
            end if
         end if
      end do
    END SUBROUTINE amoeba_private
    !BL
    FUNCTION amotry(fac)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: fac
      REAL*8 :: amotry
      REAL*8 :: fac1,fac2,ytry
      REAL*8, DIMENSION(size(p,2)) :: ptry
      fac1=(1.0d0-fac)/ndim
      fac2=fac1-fac
      ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
      ytry=kaisq( ptry, n_mole, n_atom , chg , mass, box, dfti_desc , dfti_desc_inv, log_file )
      if (ytry < y(ihi)) then
         y(ihi)=ytry
         psum(:)=psum(:)-p(ihi,:)+ptry(:)
         p(ihi,:)=ptry(:)
      end if
      amotry=ytry
    END FUNCTION amotry
  END SUBROUTINE amoeba



end module ms_evb_fit
