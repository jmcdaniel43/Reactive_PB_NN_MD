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
  subroutine read_simulation_parameters( ifile_simpmt )
    use global_variables
    character(*), intent(in) :: ifile_simpmt
    integer :: file_h = 10 , inputstatus,ind,ind1
    character(10) :: param_string
    character(30) :: param_name,line
    real*8 :: param_number
    integer:: flag_select_ensemble=0, flag_electrostatic_type=0, flag_n_step=0, flag_n_output=0, flag_temperature=0, flag_too_close=0
    integer:: flag_delta_t=0, flag_max_delta_t=0, flag_lj_cutoff=0, flag_ewald_cutoff=0, flag_Electro_cutoff=0, flag_screen_type=0, flag_springconstant=0, flag_thole=0, flag_drude_simulation=0, flag_n_threads=0, flag_cycles_per_vol_step=0,flag_vol_step=0,flag_max_vol_step=0,flag_pressure=0, flag_lj_bkghm=0,flag_lj_comb_rule=0,flag_lj_comb_rule2=0,flag_pme_disp=0,flag_Feynmann_Hibbs_correction=0,flag_Feynmann_Hibbs_forces=0
    integer:: flag_three_body_dispersion=0, flag_na_nslist=0, flag_nb_nslist=0, flag_nc_nslist=0, flag_three_body_dispersion_cutoff=0, flag_three_body_cutoff_type=0,flag_three_body_exchange=0,flag_use_verlet_list=0, flag_verlet_cutoff_lj=0,flag_verlet_cutoff_elec=0, flag_debug=0, flag_grid_Tang_Toennies=0, flag_alpha_sqrt=0, flag_pme_grid=0, flag_spline_order=0, flag_checkpoint_velocity=0


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
       Case("select_ensemble")
          select_ensemble = param_string
          flag_select_ensemble=1
          Select Case(select_ensemble)
          Case("npt")
             !******** don't run npt with a framework
             if ( framework_simulation == 1 ) then
                stop " cannot run npt ensemble with framework atoms!"
             endif
          End Select
       Case("verlet_list")
          use_verlet_list = param_string
          flag_use_verlet_list=1
       Case("electrostatic_type")
          electrostatic_type = param_string
          flag_electrostatic_type=1
       Case("lj_comb_rule")
          lj_comb_rule = param_string
          flag_lj_comb_rule=1
       Case("lj_comb_rule2")
          lj_comb_rule2 = param_string
          flag_lj_comb_rule2=1
       Case("pme_disp")
          pme_disp = param_string
          flag_pme_disp=1
       Case('three_body_dispersion')
          three_body_dispersion = param_string
          flag_three_body_dispersion=1
       Case('three_body_exchange')
          three_body_exchange = param_string
          flag_three_body_exchange=1
       Case("Feynmann_Hibbs_correction")
          Feynmann_Hibbs_correction = param_string
          flag_Feynmann_Hibbs_correction=1
       Case("Feynmann_Hibbs_forces")
          Feynmann_Hibbs_forces = param_string
          flag_Feynmann_Hibbs_forces=1
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
          n_step = NINT( param_number )
          flag_n_step=1
       Case("n_output")
          n_output = NINT( param_number )
          flag_n_output=1
       Case("checkpoint_velocity")
          n_step_velocity = NINT( param_number )
          flag_checkpoint_velocity=1
       Case("temperature")
          temp = param_number
          flag_temperature=1
       Case("cycles_per_vol_step")
          cycles_per_vol_step = NINT(param_number)
          flag_cycles_per_vol_step=1
       Case("too_close")
          too_close = param_number
          flag_too_close = 1
       Case("delta_t")
          delta_t1 = param_number
          flag_delta_t = 1
       Case("vol_step")
          vol_step = param_number
          flag_vol_step = 1
       Case("max_delta_t")
          max_delta_t1 = param_number
          flag_max_delta_t = 1
       Case("max_vol_step")
          max_vol_step = param_number
          flag_max_vol_step = 1  
       Case("pressure")
          pressure = param_number
          flag_pressure = 1
       Case("lj_bkghm")
          lj_bkghm = param_number
          flag_lj_bkghm = 1
       Case("lj_cutoff")
          lj_cutoff = param_number
          flag_lj_cutoff = 1  
       Case('three_body_dispersion_cutoff')
          three_body_dispersion_cutoff = param_number
          flag_three_body_dispersion_cutoff = 1
       Case('three_body_cutoff_type')
          three_body_cutoff_type = NINT( param_number )
          flag_three_body_cutoff_type = 1
       Case('na_nslist')
          na_nslist = NINT( param_number )
          flag_na_nslist = 1
       Case('nb_nslist')
          nb_nslist = NINT( param_number )
          flag_nb_nslist = 1
       Case('nc_nslist')
          nc_nslist = NINT( param_number )
          flag_nc_nslist = 1
       Case('verlet_cutoff_lj')
          verlet_cutoff_lj = param_number 
          flag_verlet_cutoff_lj = 1
       Case('verlet_cutoff_elec')
          verlet_cutoff_elec = param_number 
          flag_verlet_cutoff_elec = 1
       Case("ewald_cutoff")
          ewald_cutoff = param_number
          flag_ewald_cutoff = 1
       Case("Electro_cutoff")
          Electro_cutoff = param_number
          flag_Electro_cutoff = 1
       Case("alpha_sqrt")
          alpha_sqrt = param_number
          flag_alpha_sqrt = 1
       Case("pme_grid")
          pme_grid = NINT( param_number )
          flag_pme_grid = 1
       Case("spline_order")
          spline_order = NINT( param_number )
          flag_spline_order = 1
       Case("screen_type")
          screen_type = NINT( param_number )
          if ( (screen_type .NE. 0 ) .and. (screen_type .NE. 1 ) ) then
             stop "screen_type setting not recognized, should be 0 or 1"
          endif
          flag_screen_type = 1
       Case("springconstant")
          springcon = param_number * 1.889725989**3  ! this conversion factor is a.u.(e^2/Bohr^3) to e^2/A^3
          flag_springconstant = 1
       Case("thole")
          thole = param_number
          flag_thole = 1
       Case("drude_simulation")
          drude_simulation = NINT( param_number )
          if ( (drude_simulation .NE. 0) .and. (drude_simulation .NE. 1) ) then
             stop "drude_simulation setting not recognized, should be 0 or 1"
          endif
          flag_drude_simulation = 1 
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
    if( flag_select_ensemble .eq. 0 ) then
       stop "variable select_ensemble has not been given a value under the Simulation Methodology section in simulation parameter input file"
    endif
    if( flag_lj_bkghm .eq. 0 ) then
       stop "variable lj_bkghm has not been given a value, set to either '1' for buckingham or '2' for lennard jones (or 3 for a hybrid treatment)"
    endif
    if( flag_electrostatic_type .eq. 0 ) then
       stop "variable electrostatic_type has not been given a value under the Simulation Methodology section in simulation parameter input file"
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
    if( flag_too_close .eq. 0 ) then
       stop "variable too_close has not been given a value under the Simulation Parameters section in simulation parameter input file"
    endif
    if( flag_lj_cutoff .eq. 0 ) then
       stop "variable lj_cutoff has not been given a value under the Simulation Parameters section in simulation parameter input file"
    endif
    if( flag_Electro_cutoff .eq. 0 ) then
       stop "variable Electro_cutoff has not been given a value under the Simulation Parameters section in simulation parameter input file"
    endif
    if( flag_drude_simulation .eq. 0 ) then
       stop "variable drude_simulation has not been given a value under the Simulation Parameters section in simulation parameter input file"
    endif
    if( flag_n_threads .eq. 0 ) then
       stop "variable n_threads has not been given a value under the Simulation Parameters section in simulation parameter input file"
    endif


    !******************************* check ensemble specific required variables *******************************
    Select Case(select_ensemble)
    Case("npt")  
       if( flag_cycles_per_vol_step .eq. 0 ) then
          stop "variable cycles_per_vol_step has not been given a value under the Simulation Parameters section in simulation parameter input file"
       endif
       if( flag_vol_step .eq. 0 ) then
          stop "variable vol_step has not been given a value under the Simulation Parameters section in simulation parameter input file"
       endif
       if( flag_max_vol_step .eq. 0 ) then
          stop "variable max_vol_step has not been given a value under the Simulation Parameters section in simulation parameter input file"
       endif
    End Select
    if( flag_delta_t .eq. 0 ) then
       stop "variable delta_t has not been given a value under the Simulation Parameters section in simulation parameter input file"
    endif







    ! ******************************************* DEFAULT VALUES for some variables ****************************************************
    ! if these variables are not present in the simulation parameters input file, we try to set to standard settings
    ! note, these default settings will still be subject to consistency checks against the required variables further down

    if ( flag_alpha_sqrt .eq. 0 ) then
       alpha_sqrt= 0.3d0 ! in A^-1 , reasonable default value for Gaussian width parameter
       if ( electrostatic_type .eq. "pme" ) then
          write(*,*) "Default Gaussian width parameter for reciprocal space Ewald sum"
          write(*,*) "set to 0.3 inverse Angstroms.  Please make sure this is reasonable "
          write(*,*) "for the system size"
          write(*,*) ""
       endif
    endif
    if ( flag_pme_grid .eq. 0 ) then
       pme_grid= 60 ! pme_grid default value
       if ( electrostatic_type .eq. "pme" ) then
          write(*,*) "Default PME reciprocal space grid size set to 60"
          write(*,*) "Please make sure this is reasonable for the system size"
          write(*,*) ""
       endif
    endif
    if ( flag_spline_order .eq. 0 ) then
       spline_order= 6 ! beta spline order default value
       if ( electrostatic_type .eq. "pme" ) then
          write(*,*) "PME will use 6 order beta splines by default"
          write(*,*) ""
       endif
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
       write(*,*) " will write atomic velocities to file ", velocity_file
       write(*,*) "every ", n_step_velocity , " steps"
       write(*,*) ""
   endif

    if( flag_pme_disp .eq. 0 ) then                         ! pme treatment of solute-framework dispersion interactions, assume no
       pme_disp="no"
    endif
    if ( flag_three_body_dispersion .eq. 0 ) then           ! three-body, Axilrod-Teller dispersion interactions, assume no
       three_body_dispersion="no"
    endif
    if ( flag_three_body_exchange .eq. 0 ) then             ! three-body exchange, assume no
       three_body_exchange="no"
    end if
    if ( flag_grid_Tang_Toennies .eq. 0 ) then             ! grid Tang_Toennies set to yes for default
       grid_Tang_Toennies="yes"
    end if
    if( flag_ewald_cutoff .eq. 0 ) then         ! Ewald cutoff, don't need it if not using Ewald (PME)
       Select Case(electrostatic_type)
       Case("cutoff")
          ! don't need ewald cutoff here, set to zero
          ewald_cutoff=0d0
       Case("none")
          ! don't need ewald cutoff here, set to zero
          ewald_cutoff=0d0
       case default
          stop "variable ewald_cutoff has not been given a value under the Simulation Parameters section in simulation parameter input file"
       end Select
    endif
    if( flag_springconstant .eq. 0 ) then      ! spring constant for Drude oscillators, default to arbitrary value if not using Drude oscillators
       Select Case(drude_simulation)
       Case(0)
          ! don't need springconstant here, set to 0.1 au. arbitrarily
          springcon = 0.1 * 1.889725989**3  
       case default
          stop "variable springconstant has not been given a value under the Simulation Parameters section in simulation parameter input file"
       end select
    endif
    if( flag_thole .eq. 0 ) then              ! thole parameter for Drude oscillators, default to arbitrary value if not using Drude oscillators
       Select Case(drude_simulation)
       Case(0)
          ! don't need thole here, set to 2.0 arbitrarily
          thole = 2d0
       case default
          stop "variable thole has not been given a value under the Simulation Parameters section in simulation parameter input file"
       end select
    endif

    if( flag_use_verlet_list .eq. 0 ) then                         ! verlet list, don't use unless this is requested
       use_verlet_list="no"
    endif

    if( flag_debug .eq. 0 ) then                         ! default value for debug is 0
       debug=0
    endif






    !******************************************************************************************
    !
    !   From here until the end of the subroutine, all code is for checking consistency of parameter
    !   settings to make sure we don't find any inconsistencies.
    !
    !   By inconsistencies, we mean one of two things, either:
    !
    !   1) a combination of parameter settings that are logically inconsistent, and thus present confusion as to
    !   what actually the simulation is doing
    !           example: combined setting of electrostatic_type='none', and drude_simulation=1 .  ===> doesn't make
    !           sense to allow for Drude oscillators without electrostatics. In this case, the code would (probably)
    !           run fine, and run a Drude oscillator simulation in zero field (so the oscillators would be right on 
    !           their atoms), but we prefer not to do this, and stop the code for logical consistency reasons
    !
    !   2) a combination of parameter settings that the code is not meant to handle, and thus would likely produce
    !   a bug in some subroutine
    !           example: using pme treatment of dispersion ( pme_disp ='yes' ), with a non pme treatment of electrostatics
    !           ( electrostatic_type != 'pme' ) will produce bugs in the code, because pme dispersion subroutines share some
    !           data structures with pme electrostatics, and these data structures will not be initiallized if code is not 
    !           using pme electrostatics.  Obviously we could fix this, but choose not to because it's not very reasonable
    !           to use this parameter combination
    !
    !******************************************************************************************




    !**************************************************************** FORCE FIELD TYPE ***********************************************************************
    ! check these parameter consistencies

    Select Case(lj_bkghm)
    Case(2)
       ! ***************************************** check consistency of parameters for lennard jones potential *************************************************

       ! first check combination rule
       if( flag_lj_comb_rule .eq. 0 ) then 
          stop "variable lj_comb_rule has not been given a value.  For lj force field, Set to either 'opls' for opls combination rule, or 'standard' for Lorentz-Berthelot combination rule"
       endif

    case default
       stop " please set lj_bkghm variable to  '2' for a LJ potential "
    End Select




    ! cannot run drude simulation with no electrostatics!
    Select Case(electrostatic_type)
    Case("none")
       if (drude_simulation .eq.1 ) then
          stop "can't have drude_simulation without electrostatics, gonna have to change the laws of physics!"
       endif
    end select

    ! if pme_disp = yes, make sure we are using it for the correct simulation type
    Select Case(pme_disp)
    Case("yes")
       if ( select_ensemble .ne. "uvt" ) then
          stop "can only use pme dispersion in uvt ensemble"
       endif
       ! since pme dispersion is being used with hybrid_md_mc simulation, and we calculate forces using a cutoff,
       ! make user aware that we are using a different Hamiltonian for sampling
       write(*,*) ""
       write(*,*) "pme dispersion is being used to calculate energy.  However, to calculate forces"
       write(*,*) "for hybrid MD/MC moves, a cutoff is used for dispersion.  Therefore, we are "
       write(*,*) "using a different Hamiltonian for sampling.  This is ok as detailed balance"
       write(*,*) "still holds"
       write(*,*) ""

       ! can only use pme_disp with pme_electrostatics, or else there will be problems with things
       ! not being initialized, and the dispersion energy will be very wrong
       ! didn't see any need to fix this, because no idea why you would want to run pme_dispersion
       ! without pme_electrostatics
       if ( electrostatic_type .ne. "pme" ) then
          stop "cannot run pme_dispersion without pme_electrostatics !! "
       endif
    End Select


    ! if we are using verlet lists, make sure we've read in a cutoff
    Select Case(use_verlet_list)
    Case("yes")
       if ( flag_verlet_cutoff_lj == 0 ) then
          write(*,*) "verlet neighbor list was requested, but no value was input for the "
          write(*,*) "lj verlet skin thickness.  Please set this using the setting          "
          write(*,*) "`verlet_cutoff_lj =  number'       "
          stop
       endif
       if ( flag_verlet_cutoff_elec == 0 ) then
          write(*,*) "verlet neighbor list was requested, but no value was input for the "
          write(*,*) "elec verlet skin thickness.  Please set this using the setting          "
          write(*,*) "`verlet_cutoff_elec =  number'       "
          stop
       endif
       ! if we are using grid based construction of verlet list, make sure we've read in grid size
       Select Case(verlet_grid_based_construction)
       Case("yes")
          if ( ( flag_na_nslist == 0 ) .or. ( flag_nb_nslist == 0 ) .or. ( flag_nc_nslist == 0 ) ) then
             write(*,*) "must input grid size settings, na_nslist, nb_nslist, nc_nslist, if"
             write(*,*) "using grid-based verlet-list construction "
             stop
          endif
       End Select
       ! also only use with md, nvt, npt, as inserting molecules in uvt, etc. complicates updating
       ! the verlet list.  Plust we probably don't need to use it for uvt anyway
       if ( ( select_ensemble .ne. "md" ) .and. ( select_ensemble .ne. "nvt" ) .and. (select_ensemble .ne. "npt" ) ) then
          stop "can only use Verlet list for nvt or npt ensembles"
       endif
    End Select



    ! if using three body dispersion, make sure we are using it for the correct simulation type
    Select Case(three_body_dispersion)
    Case("yes")
       ! use with only nvt, npt
       if ( ( select_ensemble .ne. "nvt" ) .and. (select_ensemble .ne. "npt" ) ) then
          stop "will only calculate three_body_dispersion in nvt or npt ensembles"
       endif
       ! only use for neat fluids
       if ( framework_simulation == 1 ) then
          stop "three_body_dispersion can only be calculated for neat fluids!"
       endif

       ! now make sure all three_body_dispersion parameters have been read in
       if ( ( flag_na_nslist == 0 ) .or. ( flag_nb_nslist == 0 ) .or. ( flag_nc_nslist == 0 ) .or. ( flag_three_body_dispersion_cutoff == 0 ) ) then
          stop " For three body dispersion, must set parameters :  na_nslist , nb_nslist , nc_nslist, three_body_dispersion cutoff"
       endif

       ! need na_nslist, etc. parameters to construct verlet list
       Select Case(use_verlet_list)
       Case("yes")
          if ( ( flag_na_nslist == 0 ) .or. ( flag_nb_nslist == 0 ) .or. ( flag_nc_nslist == 0 )  ) then 
             write(*,*) "please input parameters na_nslist, nb_nslist, nc_nslist , which"
             write(*,*) "are used to determine resolution for which to grid box for     "
             write(*,*) "use in neighbor list construction "
             stop
          endif
       End Select

       ! if it isn't specified whether three body cutoff should apply to atoms or molecules, choose molecules
       if ( flag_three_body_cutoff_type == 0 ) then
          three_body_cutoff_type = 0
          write(*,*) "input parameter three_body_cutoff_type wasn't specified, so by default the"
          write(*,*) "three body dispersion cutoff will be applied on a molecule by molecule basis,"
          write(*,*) "grouping all atoms in the molecule"
          write(*,*) ""
       endif

    End Select


    Select Case(screen_type)
    Case(1)
       ! Tang-Toennies form, make sure buckingham potential is being used
       if (lj_bkghm .ne. 1 ) then
          stop "Tang_Toennies screening requested for a non-buckingham type potential"
       endif
    End Select



    ! can only use three-body exchange if using three-body-dispersion
    Select Case(three_body_exchange)
    Case("yes")
       if ( three_body_dispersion .ne. "yes" ) then
          stop " can only use three_body_exchange in conjunction with three_body_dispersion"
       endif
       ! will use the same cutoff and cell list as three body dispersion
       write(*,*) "for three body exchange interactions, code will use the same"
       write(*,*) "cutoff and cell list as for three body dispersion"
       write(*,*) ""
    End Select

    ! for Feynmann_Hibbs_correction, set to 'no' for default
    if( flag_Feynmann_Hibbs_correction .eq. 0 ) then
       Feynmann_Hibbs_correction="no"
    endif
    ! for Feynmann_Hibbs_forces, set to 'no' for default
    if( flag_Feynmann_Hibbs_forces .eq. 0 ) then
       Feynmann_Hibbs_forces="no"
    endif


  end subroutine read_simulation_parameters

end module simulation_parameters
