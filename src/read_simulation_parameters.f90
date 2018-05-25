!**********************************************************
!  This module controls reading in the simulation parameters from
!  the input file, as well as assigning default values to those
!  which are not read in. 
!
!  We have tried to include as many checks as possible on the consistency
!  of input parameters, such that the code will stop if we think there
!  are inconsistent settings.
!
!**********************************************************

module simulation_parameters
implicit none

contains
  !*******************************
  !  this subroutine reads general simulation run parameters out of 
  !  the appropriate input file
  !  nothing fancy here, we separate parameters which are read as strings into
  !  "Simulation Methodology" group, and parameters which are read as
  !  either integers or real numbers into "Simulation Parameters" group
  !*******************************

  subroutine read_simulation_parameters( file_h, ifile_simpmt, system_data, integrator_data, verlet_list_data, PME_data )
    use global_variables
    integer, intent(in)      :: file_h
    character(*), intent(in) :: ifile_simpmt
    type(system_data_type), intent(inout)   :: system_data
    type(integrator_data_type), intent(inout) :: integrator_data
    type(verlet_list_data_type), intent(inout) :: verlet_list_data
    type(PME_data_type), intent(inout)  :: PME_data

    integer :: inputstatus,ind,ind1
    character(10) :: param_string
    character(30) :: param_name,line
    real*8 :: param_number
    integer:: flag_ensemble=0, flag_n_step=0, flag_n_output=0, flag_temperature=0
    integer:: flag_delta_t=0, flag_real_space_cutoff=0, flag_n_threads=0, flag_lj_bkghm=0,flag_lj_comb_rule=0,flag_lj_comb_rule2=0
    integer:: flag_na_nslist=0, flag_nb_nslist=0, flag_nc_nslist=0, flag_verlet_cutoff=0, flag_debug=0, flag_grid_Tang_Toennies=0, flag_alpha_sqrt=0, flag_pme_grid=0, flag_spline_order=0, flag_checkpoint_velocity=0


    open(unit=file_h,file=ifile_simpmt,status="old")

    ! read Simulation_Methodology input parameters
    do
       Read(file_h,'(A)',Iostat=inputstatus) line
       If(inputstatus < 0) Exit
       ind=INDEX(line,'Simulation Methodology')
       IF(ind .NE. 0) Exit
    enddo

    do 
       Read(file_h,* ,Iostat=inputstatus) param_name, param_string 
       If(inputstatus < 0) Exit
       ind=INDEX(param_name,'Simulation')
       ind1=INDEX(param_string,'Param')
       IF( (ind .NE. 0) .and. (ind1 .NE. 0) ) Exit

       Select Case (param_name)
       Case("ensemble")
          integrator_data%ensemble = param_string
          flag_ensemble=1
       Case("lj_comb_rule")
          lj_comb_rule = param_string
          flag_lj_comb_rule=1
       Case("lj_comb_rule2")
          lj_comb_rule2 = param_string
          flag_lj_comb_rule2=1
       Case("grid_Tang_Toennies")
          grid_Tang_Toennies = param_string
          flag_grid_Tang_Toennies=1
       End Select

    enddo

    ! now read Simulation_Parameters, parameters

    do 
       Read(file_h,* ,Iostat=inputstatus) param_name, param_number 
       If(inputstatus < 0) Exit

       Select Case (param_name)
       Case("n_step")
          integrator_data%n_step = NINT( param_number )
          flag_n_step=1
       Case("n_output")
          integrator_data%n_output = NINT( param_number )
          flag_n_output=1
       Case("checkpoint_velocity")
          n_step_velocity = NINT( param_number )
          flag_checkpoint_velocity=1
       Case("temperature")
          system_data%temperature = param_number
          flag_temperature=1
       Case("delta_t")
          integrator_data%delta_t = param_number
          flag_delta_t = 1
       Case("lj_bkghm")
          lj_bkghm = param_number
          flag_lj_bkghm = 1
       Case("real_space_cutoff")
          real_space_cutoff = param_number
          flag_real_space_cutoff = 1  
       Case('na_nslist')
          verlet_list_data%na_nslist = NINT( param_number )
          flag_na_nslist = 1
       Case('nb_nslist')
          verlet_list_data%nb_nslist = NINT( param_number )
          flag_nb_nslist = 1
       Case('nc_nslist')
          verlet_list_data%nc_nslist = NINT( param_number )
          flag_nc_nslist = 1
       Case('verlet_cutoff')
          verlet_list_data%verlet_cutoff = param_number 
          flag_verlet_cutoff = 1
       Case("alpha_sqrt")
          PME_data%alpha_sqrt = param_number
          flag_alpha_sqrt = 1
       Case("pme_grid")
          PME_data%pme_grid = NINT( param_number )
          flag_pme_grid = 1
       Case("spline_order")
          PME_data%spline_order = NINT( param_number )
          flag_spline_order = 1
       Case("n_threads")
          n_threads = NINT( param_number )
          flag_n_threads = 1
       Case("debug")
          debug = NINT( param_number )
          flag_debug = 1
       End Select

    enddo

    Close(file_h)





    ! ********************************************* REQUIRED VARIABLES have no default values *******************************************
    ! Simulation_Methodology section
    if( flag_lj_bkghm .eq. 0 ) then
       stop "variable lj_bkghm has not been given a value, set to either '1' for buckingham or '2' for lennard jones (or 3 for a hybrid treatment)"
    endif
    if( flag_n_step .eq. 0 ) then
       stop "variable n_step has not been given a value under the Simulation Parameters section in simulation parameter input file"
    endif
    if( flag_n_output .eq. 0 ) then
       stop "variable n_output has not been given a value under the Simulation Parameters section in simulation parameter input file"
    endif
    if( flag_temperature .eq. 0 ) then
       stop "variable temperature has not been given a value under the Simulation Parameters section in simulation parameter input file"
    endif
    if( flag_real_space_cutoff .eq. 0 ) then
       stop "variable real_space_cutoff has not been given a value under the Simulation Parameters section in simulation parameter input file"
    endif
    if( flag_n_threads .eq. 0 ) then
       stop "variable n_threads has not been given a value under the Simulation Parameters section in simulation parameter input file"
    endif
    if( flag_delta_t .eq. 0 ) then
       stop "variable delta_t has not been given a value under the Simulation Parameters section in simulation parameter input file"
    endif
    if ( flag_verlet_cutoff == 0 ) then
       stop "variable verlet_cutoff has not been given a value under the Simulation Parameters section in simulation parameter input file"         
    endif


    ! ******************************************* DEFAULT VALUES for some variables ****************************************************
    ! if these variables are not present in the simulation parameters input file, we try to set to standard settings
    ! note, these default settings will still be subject to consistency checks against the required variables further down

    if( flag_ensemble .eq. 0 ) then
       integrator_data%ensemble = 'NVE'
    endif

    if ( flag_alpha_sqrt .eq. 0 ) then
       PME_data%alpha_sqrt= 0.3d0 ! in A^-1 , reasonable default value for Gaussian width parameter
       write(*,*) "Default Gaussian width parameter for reciprocal space Ewald sum"
       write(*,*) "set to 0.3 inverse Angstroms.  Please make sure this is reasonable "
       write(*,*) "for the system size"
       write(*,*) ""
    endif
    if ( flag_pme_grid .eq. 0 ) then
       PME_data%pme_grid= 60 ! pme_grid default value
       write(*,*) "Default PME reciprocal space grid size set to 60"
       write(*,*) "Please make sure this is reasonable for the system size"
       write(*,*) ""
    endif
    if ( flag_spline_order .eq. 0 ) then
       PME_data%spline_order= 6 ! beta spline order default value
       write(*,*) "PME will use 6 order beta splines by default"
       write(*,*) ""
    endif

    if ( flag_checkpoint_velocity .eq. 0 ) then
       checkpoint_velocity="no"
       ! make sure we're not continuing a trajectory
       Select Case(restart_trajectory)
          Case("yes")
             write(*,*) ""
             write(*,*) " if continuing a trajectory, must have the number of steps for velocity"
             write(*,*) " checkpointing set in the simulation parameters file "
             write(*,*) ""
             stop
       end Select
       
    else
       checkpoint_velocity="yes"
       write(*,*) ""
       write(*,*) " will write atomic velocities to checkpoint velocity_file "
       write(*,*) "every ", n_step_velocity , " steps"
       write(*,*) ""
   endif

    if ( flag_grid_Tang_Toennies .eq. 0 ) then             ! grid Tang_Toennies set to yes for default
       grid_Tang_Toennies="yes"
    end if

    if( flag_debug .eq. 0 ) then                         ! default value for debug is 0
       debug=0
    endif



    !******************************************************************************************
    !
    !   From here until the end of the subroutine, all code is for checking consistency of parameter
    !   settings to make sure we don't find any inconsistencies.
    !
    !******************************************************************************************

    Select Case(lj_bkghm)
    Case(2)
       ! first check combination rule
       if( flag_lj_comb_rule .eq. 0 ) then 
          stop "variable lj_comb_rule has not been given a value.  For lj force field, Set to either 'opls' for opls combination rule, or 'standard' for Lorentz-Berthelot combination rule"
       endif

    case default
       stop " please set lj_bkghm variable to  '2' for a LJ potential "
    End Select


    ! using grid based construction of verlet list, make sure we've read in grid size
    if ( ( flag_na_nslist == 0 ) .or. ( flag_nb_nslist == 0 ) .or. ( flag_nc_nslist == 0 ) ) then
       write(*,*) "must input grid size settings, na_nslist, nb_nslist, nc_nslist, if"
       write(*,*) "using grid-based verlet-list construction "
       stop
    endif

  end subroutine read_simulation_parameters

end module simulation_parameters
