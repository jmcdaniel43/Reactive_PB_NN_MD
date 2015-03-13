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
    integer:: flag_select_ensemble=0,flag_sapt_type_ff=0, flag_solute_cross_parameter_set=0, flag_C8_10_dispersion_terms=0, flag_hybrid_md_mc=0, flag_electrostatic_type=0, flag_n_step=0, flag_n_output=0, flag_n_update_stepsize=0, flag_temperature=0, flag_too_close=0
    integer:: flag_delta_t=0, flag_max_delta_t=0, flag_target_acceptance=0, flag_lj_cutoff=0, flag_ewald_cutoff=0, flag_Electro_cutoff=0, flag_screen_type=0, flag_springconstant=0, flag_thole=0, flag_drude_simulation=0, flag_n_threads=0, flag_gcmc_ins_rm=0, flag_chemical_potential=0,flag_cycles_per_vol_step=0,flag_target_acceptance_v=0,flag_vol_step=0,flag_max_vol_step=0,flag_pressure=0,flag_dispersion_lrc=0,flag_damp_solute_solute_dispersion=0, flag_lj_bkghm=0,flag_lj_comb_rule=0,flag_lj_comb_rule2=0,flag_pme_disp=0,flag_C12_dispersion=0,flag_Feynmann_Hibbs_correction=0,flag_Feynmann_Hibbs_forces=0,flag_penetration_force_field=0
    integer:: flag_three_body_dispersion=0, flag_na_nslist=0, flag_nb_nslist=0, flag_nc_nslist=0, flag_three_body_dispersion_cutoff=0, flag_three_body_cutoff_type=0,flag_three_body_exchange=0,flag_orientation_try=0, flag_use_verlet_list=0, flag_verlet_cutoff_lj=0,flag_verlet_cutoff_elec=0, flag_debug=0, flag_grid_Tang_Toennies=0, flag_alpha_sqrt=0, flag_pme_grid=0, flag_spline_order=0


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
          Case("egc")
             !******** don't run egc with puddle filling
             if ( puddle_filling_on .eq. "yes" ) then
                stop " cannot run egc ensemble with puddle filling"
             endif
          End Select
       Case("sapt_type_ff")
          sapt_type_ff = param_string
          flag_sapt_type_ff=1
       Case("solute_cross_parameter_set")
          solute_cross_parameter_set = param_string
          flag_solute_cross_parameter_set=1
       Case("C8_10_dispersion_terms")
          C8_10_dispersion_terms = param_string
          flag_C8_10_dispersion_terms=1
       Case("C12_dispersion")
          C12_dispersion = param_string
          flag_C12_dispersion=1
       Case("damp_solute_solute_dispersion")
          damp_solute_solute_dispersion = param_string
          flag_damp_solute_solute_dispersion=1
       Case("hybrid_md_mc")
          hybrid_md_mc = param_string
          flag_hybrid_md_mc=1
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
       Case("penetration_force_field")
          penetration_force_field = param_string
          flag_penetration_force_field=1
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
       Case("n_update_stepsize")
          n_update_stepsize = NINT( param_number )
          flag_n_update_stepsize=1
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
       Case("target_acceptance")
          target_acc = param_number
          flag_target_acceptance = 1
       Case("target_acceptance_v")
          target_acc_v = param_number
          flag_target_acceptance_v = 1
       Case("pressure")
          pressure = param_number
          flag_pressure = 1
       Case("lj_bkghm")
          lj_bkghm = param_number
          flag_lj_bkghm = 1
       Case("lj_cutoff")
          lj_cutoff = param_number
          flag_lj_cutoff = 1
       Case("dispersion_lrc")
          lj_lrc = NINT(param_number)
          flag_dispersion_lrc = 1
          ! make sure we're using long range correction appropriately
          Select Case(lj_lrc)
          Case(1)
             if( ( select_ensemble .ne. "nvt" ) .and. ( select_ensemble .ne. "npt" ) ) then
                stop "can only use dispersion long range correction for nvt or npt ensembles"
             endif
             if (lj_shift .eq. 1 ) then
                stop " Cannot use long range corection with shifted potential"
             endif
             if (framework_simulation .eq. 1 ) then
                write(*,*) "long range correction is automatically set to zero for a framework simulation"
             endif
          End Select
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
       Case("gcmc_ins_rm")
          gcmc_ins_rm = param_number
          flag_gcmc_ins_rm = 1
       Case("orientation_try")
          orientation_try = NINT( param_number )
          flag_orientation_try = 1
       Case("chemical_potential")
          chem_potential = param_number
          flag_chemical_potential = 1
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
    if( flag_sapt_type_ff .eq. 0 ) then
       ! in older code, we had the parameter 'energy_decomposition' control both whether we were using a sapt-type force field, as well
       ! as whether we were calculating an explicit energy decompostion.  Therefore, the parameter  energy_decomposition was read in from the simulation_parameters input file
       ! now, we use the parameter 'sapt_type_ff' to signal whether a sapt-type force field is being used, so this is now the variable to read in , and the value
       ! of energy_decomposition will be automatically assigned by the code (see below)
       write(*,*) ""
       write(*,*) "variable sapt_type_ff has not been given a value under the Simulation Methodology section in simulation parameter input file"
       write(*,*) "note that this input parameter used to be called 'energy_decomposition', so if the parameter 'energy_decomposition' is "
       write(*,*) "assigned a value in the simulation parameter input file, and this is meant to determine whether a sapt-type force field"
       write(*,*) "is being used, please change the name of this parameter to 'sapt_type_ff' in the input file"
       write(*,*) ""
       stop
    endif
    if( flag_hybrid_md_mc .eq. 0 ) then
       stop "variable hybrid_md_mc has not been given a value under the Simulation Methodology section in simulation parameter input file"
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
    if( flag_n_update_stepsize .eq. 0 ) then
       stop "variable n_update_stepsize has not been given a value under the Simulation Parameters section in simulation parameter input file"
    endif
    if( flag_temperature .eq. 0 ) then
       stop "variable temperature has not been given a value under the Simulation Parameters section in simulation parameter input file"
    endif
    if( flag_too_close .eq. 0 ) then
       stop "variable too_close has not been given a value under the Simulation Parameters section in simulation parameter input file"
    endif
    if( flag_target_acceptance .eq. 0 ) then
       stop "variable target_acceptance has not been given a value under the Simulation Parameters section in simulation parameter input file"
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
       !*********************** grand canonical (uvt) **********************
    Case("uvt")
       if( flag_gcmc_ins_rm .eq. 0 ) then
          stop "variable gcmc_ins_rm has not been given a value under the Simulation Parameters section in simulation parameter input file"
       endif
       if( flag_chemical_potential .eq. 0 ) then
          stop "variable chemical_potential has not been given a value under the Simulation Parameters section in simulation parameter input file"
       endif
       !*********************** isothermal,isobaric (npt) **********************
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
       if( flag_target_acceptance_v .eq. 0 ) then
          stop "variable target_acceptance_v has not been given a value under the Simulation Parameters section in simulation parameter input file"
       endif
       if( flag_pressure .eq. 0 ) then
          stop "variable pressure has not been given a value under the Simulation Parameters section in simulation parameter input file"
       endif
    End Select


    !********************************** check move specific required variables
    Select Case(hybrid_md_mc)
    Case("yes")
       if( flag_delta_t .eq. 0 ) then
          stop "variable delta_t has not been given a value under the Simulation Parameters section in simulation parameter input file"
       endif
       if( flag_max_delta_t .eq. 0 ) then
          stop "variable max_delta_t has not been given a value under the Simulation Parameters section in simulation parameter input file"
       endif
    End Select






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

    if( flag_C8_10_dispersion_terms .eq. 0 ) then           ! C6-C10 treatment of dispersion, importance of this depends on lj_bkghm treatment
       Select Case(lj_bkghm)
       Case(2)
          ! lj potential, set C8_10_dispersion_terms ="no"
          C8_10_dispersion_terms = "no"
       case default
          Select Case( sapt_type_ff )
          Case("yes")
             write(*,*) "variable C8_10_dispersion_terms has not been given a value under the Simulation Methodology section"
             stop
          case("no")
             ! bkghm potential, set C8_10_dispersion_terms = "no" if not using SAPT ff, otherwise stop and require this as input
             write(*,*) ""
             write(*,*) "variable C8_10_dispersion_terms has not been given a value under the Simulation Methodology section"
             write(*,*) "in the simulation parameter input file.  We assign a default value of C8_10_dispersion_terms='no', "
             write(*,*) "however this setting may create a conflict with a different parameter setting "
             C8_10_dispersion_terms = "no"
          end Select
       end Select
    endif
    if( flag_C12_dispersion .eq. 0 ) then                   ! C12 dispersion parameters, assume no unless C8_10 is set to yes, then require this
       Select Case(C8_10_dispersion_terms)
       Case("yes")
          write(*,*) "variable C12_dispersion has not been given a value under the Simulation Methodology section"
          stop         
       Case("no")
          C12_dispersion = "no"
       End Select
    end if
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
    if( flag_dispersion_lrc .eq. 0 ) then  ! long range correction to van-der-waals cutoff distance, assume no in certain cases
       Select Case(framework_simulation)
       Case(1)
          lj_lrc = 0
       case default
          stop "variable dispersion_lrc has not been given a value under the Simulation Parameters section in simulation parameter input file"
       end select
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
    if( flag_screen_type .eq. 0 ) then          ! screen type set to zero, unless we're using a sapt-type force field, and then we usually use this, so don't set to default value in that case
       Select Case( sapt_type_ff )
       Case("yes")
          stop "variable screen_type has not been given a value under the Simulation Parameters section in simulation parameter input file"
       case default
          screen_type = 0
       End Select
    end if
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
    if( flag_orientation_try .eq. 0 ) then     ! # of orientations for orientation bias particle insertions/deletions, default, orientation_try=1, means don't use orientation bias
       orientation_try = 1
    else
       if ( orientation_try > 1 ) then
          write(*,*) ""
          write(*,*) "parameter orientation_try > 1, therefore an orientation bias algorithm will"
          write(*,*) "be used for more efficient insertions/deletions in gcmc"
          write(*,*) ""
       endif
    endif

    if( flag_solute_cross_parameter_set .eq. 0 ) then           ! this variable is for a framework simulation, if not a framework simulation, set to "no"
       Select Case(framework_simulation)
       Case(0)
          solute_cross_parameter_set="no"
       case default
          stop "variable cross_parameter_set has not been given a value under the Simulation Methodology section in simulation parameter input file"
       End Select
    end if

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
    Case(1)
       ! **************************************** check consistency of parameters for buckingham potential ****************************************************

       ! first check combination rule
       if( flag_lj_comb_rule .eq. 0 ) then 
          stop "variable lj_comb_rule has not been given a value.  For bkghm force field, Set to either 'standard' for all geometric combination rules, or 'ZIFFF' for ZIF FF exponent combination rule"
       else
          Select Case(sapt_type_ff)
          Case("yes")
             ! here lj_comb_rule parameter must be ZIFFF
             if ( lj_comb_rule .ne. "ZIFFF" ) then
                stop "for sapt-type force field, lj_comb_rule parameter must be set to 'ZIFFF' "
             endif
          case default
             if ( (lj_comb_rule .ne. "ZIFFF" ) .and. (lj_comb_rule .ne. "standard" ) )  then 
                write(*,*) ""
                write(*,*) "For a Buckingham force field, variable lj_comb_rule must be set to either "
                write(*,*) "'standard' for all geometric combination rules, or 'ZIFFF' for ZIF FF exponent combination rule"
                write(*,*) ""
                stop
             endif
          End Select
       endif


    Case(2)
       ! ***************************************** check consistency of parameters for lennard jones potential *************************************************

       ! first check combination rule
       if( flag_lj_comb_rule .eq. 0 ) then 
          stop "variable lj_comb_rule has not been given a value.  For lj force field, Set to either 'opls' for opls combination rule, or 'standard' for Lorentz-Berthelot combination rule"
       else
          if ( (lj_comb_rule .ne. "opls" ) .and. (lj_comb_rule .ne. "standard" ) )  then 
             write(*,*) ""
             write(*,*) "For a lennard jones force field, variable lj_comb_rule must be set to either "
             write(*,*) "'standard' for Lorentz-Berthelot combination rules, or 'opls' for opls combination rule"
             write(*,*) ""
             stop
          endif
       endif

       ! sapt-type force field makes no sense for lennard jones potential
       Select Case(sapt_type_ff)
       Case("yes")
          stop "sapt_type_ff='yes' makes no sense with lennard jones potential!"
       End Select

       ! C8,C10, etc terms make no sense for lennard jones potential
       Select Case(C8_10_dispersion_terms)
       Case("yes")
          stop "C8,C10,etc dispersion terms make no sense with lennard jones potential!"
       End Select


    Case(3)
       ! ***************************************** check consistency of parameters for hybrid lennard jones/buckingham  potential *******************************
       ! Here we are fairly strict on settings, as this is a special case

       ! first check combination rule, need two of them here, one for lj potential (lj_comb_rule2) and a second (lj_comb_rule) for the buckingham terms
       if( flag_lj_comb_rule .eq. 0 ) then 
          write(*,*) ""
          write(*,*) "variable lj_comb_rule has not been given a value. For a hybrid lj/bkghm potential (lj_bkghm=3)"
          write(*,*) "'lj_comb_rule' applies to the bkghm portion of the force field. Must be set to 'ZIFFF' in this case" 
          write(*,*) ""
          stop    
       else
          ! here lj_comb_rule parameter must be ZIFFF
          if ( lj_comb_rule .ne. "ZIFFF" ) then
             stop "for lj_bkghm=3, lj_comb_rule parameter must be set to 'ZIFFF' "
          endif
       endif

       if( flag_lj_comb_rule2 .eq. 0 ) then 
          write(*,*) ""
          write(*,*) "variable lj_comb_rule2 has not been given a value. For a hybrid lj/bkghm potential (lj_bkghm=3)"
          write(*,*) "'lj_comb_rule2' applies to the lj portion of the force field. Set 'standard' for Lorentz-Berthelot combination rule"
          write(*,*) ""
          stop 
       else
          if ( lj_comb_rule2 .ne. "standard" )  then 
             write(*,*) ""
             write(*,*) "For a hybrid lennard jones/bkghm force field, variable lj_comb_rule2 must be set to "
             write(*,*) "'standard' for Lorentz-Berthelot combination rule"
             write(*,*) ""
             stop
          endif
       endif

       ! must be using lj_bkghm=3 for a framework simulation
       Select Case(framework_simulation)
       Case(0) 
          stop "can only use lj_bkghm=3 for a rigid framework simulation"
       End Select

       ! for the case of lj_bkghm=3, must use sapt_type_ff with C8_10_dispersion and C12_dispersion, and also damp_solute_solute_dispersion must be set to 'no' since solute-solute is lj potential
       Select Case(sapt_type_ff)
       Case("yes")
          continue
       case default
          stop "must use setting sapt_type_ff='yes' for lj_bkghm=3"
       End Select

       Select Case(C8_10_dispersion_terms)
       Case("yes")
          continue
       case default
          stop "must use setting C8_10_dispersion_terms='yes' for lj_bkghm=3"
       End Select

       Select Case(C12_dispersion)
       Case("yes")
          continue
       case default
          stop "must use setting C12_dispersion_terms='yes' for lj_bkghm=3"
       End Select

       Select Case(damp_solute_solute_dispersion)
       Case("no")
          continue
       Case default
          stop "must use setting damp_solute_solute_dispersion='no' for lj_bkghm=3"
       End Select

       ! also, explicit_cross_term_exponents must be set to no, and solute_cross_parameter set must be set to "yes" for reading in parameters for lj_bkghm=3

       Select Case(explicit_cross_term_exponents)
       Case("no")
          continue
       case default
          stop "must use setting explicit_cross_term_exponents='no'  for lj_bkghm=3"
       end Select

       Select Case(solute_cross_parameter_set)
       Case("yes")
          continue
       case default
          stop "must use setting solute_cross_parameter_set='yes'  for lj_bkghm=3"
       end Select

       ! cannot use long range correction or shift with lj_bkghm=3
       Select Case(lj_shift)
       Case(0)
          continue
       Case default
          stop "must use setting lj_shift=0 for lj_bkghm=3"
       end Select

       Select Case(lj_lrc)
       Case(0)
          continue
       case default
          stop "must use setting lj_lrc=0 for lj_bkghm=3"
       end Select

       ! don't use with egc ensemble, because in egc ensemble we scale the ff parameters, and right now code isn't set up to scale both bkghm and lj interactions simulataneously.  Could easily fix this, but we never use egc ensemble anyway...

       Select Case(select_ensemble)
       Case("egc")
          stop "cannot use setting lj_bkghm=3 with egc ensemble"
       End Select


    case default
       stop " please set lj_bkghm variable to either '1' for a Buckingham potential or '2' for a LJ potential (or 3 for a hybrid treatment)"
    End Select





    !**************************************************************** SAPT-TYPE FORCE FIELDS  ***********************************************************************
    !****************  SAPT type force field requires many specific run parameters
    !****************  this section is all about checking consistency between input parameters
    !****************************************************************************************************************************************************************
    Select Case(sapt_type_ff)     
    Case("yes")
       ! give explicit energy decomposition unless either we are using pme_dispersion , or a hybrid lj/bkghm force field
       if ( ( pme_disp .eq. "yes" ) .or. ( lj_bkghm .eq. 3 ) ) then
          energy_decomposition = "no"
       else
          energy_decomposition = "yes"
       end if

       ! warn if using a sapt_type_ff without drude_oscillators
       if (drude_simulation .ne. 1 ) then
          write(*,*) "Running a sapt_type force field simulation without drude oscillators.  Are you sure you want to do this?"
          write(*,*) ""
          !       stop "Do you really want a sapt_type force field without drude oscillators??? Gonna have to change the source code"
       endif
       ! don't run sapt_type force field  without electrostatics
       select case(sapt_type_ff)
       case("none")
          stop "can't run sapt_type force field without electrostatics"
       end select

       ! sapt_type_ff only makes sense for buckingham ( allow lj_bkghm=1 (buckingham) or lj_bkghm=3 (lj/buckingham combination)
       if ( (lj_bkghm .eq. 1 ) .or. (lj_bkghm .eq. 3 ) ) then
          continue
       else
          stop "sapt_type_ff='yes' only makes sense for buckingham potential"
       endif


       ! must use C8_10_dispersion terms if we're reading in extra solute parameter set
       Select Case(solute_cross_parameter_set)
       Case("yes")
          Select Case(C8_10_dispersion_terms)
          Case("no")
             stop "must use C8_10_dispersion terms if reading in an extra solute parameter set"
          End Select
       End Select


       Select Case(solute_cross_parameter_set)
       Case("yes")
          ! make sure this is a framework simulation, otherwise solute_cross_parameter_set="yes" doesn't make sense
          Select Case(framework_simulation)
          Case(0)
             stop "setting 'cross_parameter_set=yes' doesn't makes sense for a non-rigid framework simulation"
          End Select
       End Select

       if( flag_damp_solute_solute_dispersion .eq. 0 ) then
          ! if this parameter isn't given, and if we are not using C8_C10_dispersion terms, assume no damping
          ! for C6 terms
          Select Case(C8_10_dispersion_terms)
          Case("no")
             damp_solute_solute_dispersion="no"
          case default
             stop "variable damp_solute_solute_dispersion has not been given a value under the Simulation Methodology section in simulation parameter input file"
          end select
       else
          ! damp_solute_solute_dispersion has been read in, make sure there is no damping for C6 treatment
          Select Case(C8_10_dispersion_terms)
          Case("no")
             if ( damp_solute_solute_dispersion .eq. "yes" ) then
                damp_solute_solute_dispersion= "no"
                write(*,*) "dispersion damping will not be used for a C6 treatment of dispersion"
             endif
          End Select
       endif




    Case("no")
       ! don't have explicit energy decomposition without sapt force field
       energy_decomposition = "no"

       Select Case(C8_10_dispersion_terms)
       Case("yes")
          write(*,*) "setting of variable C8_10_dispersion_terms will be overidden and set to 'no' as this is not a sapt_type_ff simulation"
       End Select
       C8_10_dispersion_terms="no"




       if ( lj_bkghm .eq. 1 ) then
          ! warn if ZIFFF combination rule is specified without sapt-type force field
          if ( lj_comb_rule .eq. "ZIFFF" ) then
             write(*,*) "Are you sure you want to be using ZIF FF combination without a sapt-type  force field?"
          endif
       endif


       if( flag_solute_cross_parameter_set .eq. 0 ) then
          solute_cross_parameter_set="no"
       else
          Select Case(solute_cross_parameter_set)
          Case("yes")
             write(*,*) "setting of variable solute_cross_parameter_set will be overidden and set to 'no' as this is not a sapt_type_ff simulation"
          End Select
          solute_cross_parameter_set="no"
       endif


       if( flag_damp_solute_solute_dispersion .eq. 0 ) then
          damp_solute_solute_dispersion="no"
       else
          Select Case(damp_solute_solute_dispersion)
          Case("yes")
             write(*,*) "setting of variable damp_solute_solute_dispersion will be overidden and set to 'no' as this is not a sapt_type_ff simulation"
          End Select
          damp_solute_solute_dispersion="no"
       endif




    case default
       stop "sapt_type_ff must be set to either 'yes' or 'no'"



    End Select



    !*************************************************************************



    ! cannot run drude simulation with no electrostatics!
    Select Case(electrostatic_type)
    Case("none")
       if (drude_simulation .eq.1 ) then
          stop "can't have drude_simulation without electrostatics, gonna have to change the laws of physics!"
       endif
    end select


    Select Case(solute_cross_parameter_set)
    Case("yes")
       ! make sure this is a framework simulation, otherwise solute_cross_parameter_set="yes" doesn't make sense
       Select Case(framework_simulation)
       Case(0)
          stop "setting 'cross_parameter_set=yes' doesn't makes sense for a non-rigid framework simulation"
       End Select
    End Select


    ! if pme_disp = yes, make sure we are using it for the correct simulation type
    Select Case(pme_disp)
    Case("yes")
       if ( select_ensemble .ne. "uvt" ) then
          stop "can only use pme dispersion in uvt ensemble"
       endif
       if ( hybrid_md_mc .ne. "yes" ) then
          stop "can only use pme dispersion with hybrid_md_mc simulation"
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
       ! also only use with md, nvt, npt, as inserting molecules in uvt, etc. complicates updating
       ! the verlet list.  Plust we probably don't need to use it for uvt anyway
       if ( ( select_ensemble .ne. "md" ) .and. ( select_ensemble .ne. "nvt" ) .and. (select_ensemble .ne. "npt" ) ) then
          stop "can only use Verlet list for nvt or npt ensembles"
       endif
       ! also only use for hybrid md/mc at this point
       if ( hybrid_md_mc .ne. "yes" ) then
          stop "can only use Verlet list for hybrid_md_mc simulation"
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
       ! sapt_type_ff must be set to "yes", as C9 coefficients will be read in using that module
       if ( sapt_type_ff .ne. "yes" ) then
          stop "sapt_type_ff must be set to 'yes' when using three body dispersion"
       endif

       ! damping is always used for 3body dispersion
       if ( damp_solute_solute_dispersion .ne. "yes" ) then
          write(*,*) "although damp_solute_solute_dispersion was set to 'no' ( no two body dispersion damping )"
          write(*,*) "damping will be used for three body dispersion interactions"
       endif

       ! warn that forces are not computed for three body dispersion, so hybrid moves are generated without this component of the Hamiltonian
       if ( hybrid_md_mc .eq. "yes" ) then
          write(*,*) "NOTE:  Forces are not computed for three body dispersion interactions,"
          write(*,*) "so Hamiltonian for hybrid MD sampling will be different from real Hamiltonian"
          write(*,*) "(i.e. missing three-body-dispersion).  This does not violate detailed"
          write(*,*) "balance so is ok"
          write(*,*) ""
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


    Select Case(C12_dispersion)
    Case("yes")
       if ( C8_10_dispersion_terms .ne. "yes" ) then
          stop "must use C8_C10_dispersion_terms='yes' in order to use C12 terms"
       endif
       ! only use C12 dispersion without explicit cross term exponents, if reading in two solute sets.  NOTE for one solute set, explicit_cross_term_exponents must be set to yes
       Select Case( solute_cross_parameter_set )
       Case("yes")
          if ( explicit_cross_term_exponents .ne. "no" ) then
             stop "if using C12_dispersion, explicit_cross_term_exponents must be set to 'no' in glob_v.f90"
          endif
       end select
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




    ! for penetration_force_field, set to 'no' for default
    if( flag_penetration_force_field .eq. 0 ) then
       penetration_force_field="no"
    else
       Select Case(penetration_force_field)
       Case("yes")
          ! warn that energy components will not add up to total energy with penetration ff
          write(*,*) ""
          write(*,*) "Using a piecewise force field to describe charge penetration"
          write(*,*) "The short range part of this piece wise force field is not included"
          write(*,*) "in any of the decomposed energy components, and therefore the "
          write(*,*) "printed energy components will not add up to the total potential"
          write(*,*) "This discrepancy is equal to the short range piecewise penetration"
          write(*,*) "force field energy."
          write(*,*) ""

          Select Case(sapt_type_ff)
          Case("no")
             stop "must use sapt_type force field with penetration force field"
          End Select
          ! use of C8,C10, etc terms implies buckingham force field, so we don't have to check that as well
          Select Case(C8_10_dispersion_terms)
          Case("no")
             stop "must use C8_10_dispersion_terms with penetration force field"
          End Select
          ! cannot use this with expanded ensemble as global arrays are not implemented for changing atom type
          Select Case(select_ensemble)
          Case("egc")
             stop "cannot use penetration force field with egc ensemble"
          End Select
       End Select
    endif

    ! for Feynmann_Hibbs_correction, set to 'no' for default
    if( flag_Feynmann_Hibbs_correction .eq. 0 ) then
       Feynmann_Hibbs_correction="no"
    else
       Select Case(Feynmann_Hibbs_correction)
       Case("yes")
          Write(*,*) "Feynmann Hibbs correction has been requested.  This will be used for the LJ/Buckingham part of the potential."
          Select Case( penetration_force_field)
          Case("yes")
             stop "Feynmann Hibbs correction is not implemented for use with penetration force fields"
          End Select
          Select Case( lj_shift)
          Case(1)
             stop "Feynmann Hibbs correction is not implemented for use with shifted potentials"
          End Select
          Select Case( lj_lrc)
          Case(1)
             write(*,*) "There will be no Feynmann Hibbs correction in the dispersion long range correction"
          End Select
       End Select
    endif
    ! for Feynmann_Hibbs_forces, set to 'no' for default
    if( flag_Feynmann_Hibbs_forces .eq. 0 ) then
       Feynmann_Hibbs_forces="no"
    else
       Select Case(Feynmann_Hibbs_forces)
       Case("yes")
          Select Case(Feynmann_Hibbs_correction)
          Case("no")
             stop " Can't use Feynmann_Hibbs_forces without Feynmann_Hibbs_correction='yes' ! "
          End Select
          ! Feynmann_Hibbs forces not coded in for lj potential
          Select Case(lj_bkghm)
          Case(1)
             continue
          case default
             stop " Feynmann_Hibbs forces are not coded in for lj potential"
          End Select
          Write(*,*) "Feynmann Hibbs forces has been requested.  The forces from the Feynmann-Hibbs correction will be used for hybrid MD/MC moves."
       End Select
    endif
    ! if using Feynman-Hibbs correction with hybrid MD/MC moves, note if not computing Feynman-Hibbs forces for these moves
    if ( ( Feynmann_Hibbs_correction .eq. "yes" ) .and. ( Feynmann_Hibbs_forces .ne. "yes" ) .and. ( hybrid_md_mc .eq. "yes" ) ) then
       write(*,*) "NOTE:  While Feynman-Hibbs correction is being used for energy, this correction"
       write(*,*) "is not being used for forces in hybrid MD/MC moves, and thus we are using a"
       write(*,*) "different Hamiltonian for sampling.  This is okay, as detailed balance is not violated"
    endif



  end subroutine read_simulation_parameters

end module simulation_parameters
