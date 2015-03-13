!*****************************************************************************
! this module is all about reading in and setting up the data structures for a 
! SAPT-based force field.
!
! for a SAPT force field, there is an explicit energy decomposition, in which the buckingham terms
! are separated into exchange, electrostatic,induction, and possible dispersion contributions
! also , the electrostatic calculation will separate electrostatic and induction energies,
! basically by calculating energy with and without drude oscillators
!
! the subroutines in this  module are for reading an input file and filling in global arrays with 
! parameters for such a simulation
!
!*****************************************************************************

module sapt_ff_routines
  use routines
  implicit none

  integer :: flag_solute_param

contains
  !*******************************
  ! this subroutine reads parameter file for a sapt_based force field format
  ! which has an explicit energy decomposition
  !
  ! the structure of the input file depends on a lot of different factors, such as
  ! whether there is a different force field for solute-solute and solute-framework 
  ! interactions (in GCMC), whether C8,C10 terms are given, etc.
  !
  ! temp_lj_parameter is used to temporarily store solute-solute parameters if these are 
  ! different than solute-framework parameters
  !
  ! note that this subroutine has a lot of excessive options that will usually not be used
  ! They were added when we were testing different force field development options, and will
  ! probably not be used in the future
  !
  !   here are some of the old force field development ideas that we don't really use anymore
  !
  !   - we used to use an exponential charge penetration term in addition to C6,C8,C10 terms to model
  !    the dispersion energy.  That's why there are 5 indices for the last element of
  !    atype_bkghm_decomp array.  We don't use this term for dispersion anymore, but this option
  !    is still coded in, so now we just input zero's for these parameters
  !
  !   - exponents used to be read in as explicit cross terms, as there was no explicit combination rule
  !    this is probably unneccesary as now we will probably always use an explicit combination rule
  !
  ! 
  !   NOTE:  temp_lj_parameter array is used to temporarily store parameters, until all combination
  !   rules can be carried out and all final lj/bkghm parameters will be put into atype_bkghm_decomp
  !   and atype_lj_parameter arrays.  
  !   the indices used to store different parameters in temp_lj_parameter array depend on 
  !   the force field parameters being read in.
  !     -- if we are using solute_cross_parameter_set="yes", then solute-solute parameters will
  !       be stored in temp_lj_parameter(:,:,1-8) as Aexch, Aelec, Aind, Adhf, Adisp, C6, C8, C10
  !       respectively
  !     -- if we are using solute_cross_parameter_set="no", this means that Aexch, Aelec, Aind, Adhf
  !        are the same for solute-solute and solute-framework parameters.  However, dispersion is
  !        different, as we always use C6 only in "SYM" , but we always use C6,C8,C10 for 
  !        solute-framework interactions.  In this case, solute-solute dispersion parameters will be
  !        stored as:  Adisp --> temp_lj_parameter(:,:,1), C6 --> temp_lj_parameter(:,:,3), 
  !         C8 --> temp_lj_parameter(:,:,4),  C10 --> temp_lj_parameter(:,:,5), 
  !     -- if explicit_cross_term_exponents = "no", then we use temp_lj_parameter(:,:,9) array 
  !        elements to store solute-solute exponents
  !
  !********************************
  subroutine read_param_decomp( ifile_pmt , temp_lj_parameter , n_type_s)
    use global_variables
    character(*),intent(in)::ifile_pmt
    real*8,dimension(:,:,:),intent(out) :: temp_lj_parameter
    integer,intent(out) :: n_type_s

    integer::file_h=15, i_mole,j_mole, n_type_f, n_cross, inputstatus,i_search,j_search,n_penetration,i_type,j_type,i,j
    integer::ind=0,ind_dhf=0,ind1=0,ind2=0,ind_frame=0
    character(MAX_ANAME):: atom1, atom2
    real*8,dimension(3)::store
    real*8,dimension(:),allocatable::temp_cross
    character(25)::junk
    character(2):: junk_label
    character(50)::line

    open(unit=file_h,file=ifile_pmt,status="old")



    !****************************************************************************************************************************************************************
    !**********
    !**********                                                          solute parameters
    !**********
    !****************************************************************************************************************************************************************

    ! get solute atom types
    do
       Read(file_h,'(A)',Iostat=inputstatus) line
       If(inputstatus < 0) Exit
       ind=INDEX(line,'solute_species')
       IF(ind .NE. 0) Exit
    enddo
    Read(file_h,'(A)',Iostat=inputstatus) line
    read(file_h,*) n_type_s

    ! zero dispersion buckingham elements, in case we are not using these for dispersion, also setting the dhf parameter to zero lets us know whether we have read in explicit cross terms
    atype_bkghm_decomp(:,:,:)=0d0
    ! zero dispersion coefficients, if we are using C12 terms, need these to be zeroed as default
    atype_lj_parameter(:,:,:)=0d0


    !*********************************** decide how many sets of solute parameters we are reading in (one or two) *************************
    Select Case( solute_cross_parameter_set )
    Case("yes")

       call read_solute_parameters_two_set(file_h,n_type_s,C12_dispersion,explicit_cross_term_exponents,atype_name,atype_chg,atype_bkghm_decomp,atype_lj_parameter,temp_lj_parameter, atype_pol, lj_bkghm)

    Case ("no")
       !************************************** here we are only reading in one set of solute parameters ********************************************************

       ! input is determined by whether or not this is a simulation with C8 and C10 dispersion terms
       Select Case( C8_10_dispersion_terms)
       Case('yes')

          call read_solute_parameters_one_set(file_h,n_type_s,C12_dispersion,explicit_cross_term_exponents,atype_name,atype_chg,atype_bkghm_decomp,atype_lj_parameter,atype_pol)

          ! allow for possibly different dispersion parameters
          ind=0
          do
             Read(file_h,'(A)',Iostat=inputstatus) line
             If(inputstatus < 0) Exit
             ind_dhf=INDEX(line,'solute dhf cross terms')
             IF(ind_dhf .NE. 0) Exit
             ind=INDEX(line,'solute-solute-dispersion')
             IF(ind .NE. 0) Exit
             ind1=INDEX(line,'solute-solute')
             IF(ind1 .NE. 0) Exit
             ind_frame=INDEX(line,'framework_species')
             IF(ind_frame .NE. 0) Exit
          enddo

          if ( ind .ne. 0 ) then

             ! make sure we're not using C12 dispersion here
             Select Case(C12_dispersion)
             Case("yes")
                stop "cannot use C12_dispersion='yes' with settings of solute_cross_parameter_set='no' but reading in a set of dispersion terms under 'solute-solute-dispersion' heading"
             End Select

             ! get solute-solute terms and store in temp_lj_parameter array
             ! 5 indices in this array, fill with A_disp, empty, C6,C8,C10
             do i_mole=1, n_type_s
                read(file_h,*) junk_label, temp_lj_parameter(i_mole,i_mole,1),temp_lj_parameter(i_mole,i_mole,3),temp_lj_parameter(i_mole,i_mole,4),temp_lj_parameter(i_mole,i_mole,5)
             enddo
             flag_solute_param=1
          else
             flag_solute_param=0
          endif

       Case('no')
          ! read in order for each atom, chg, A_exch, A_elec, A_induc,A_dhf,C6, pol
          do i_mole=1, n_type_s
             read(file_h,*) atype_name(i_mole),atype_chg(i_mole),atype_bkghm_decomp(i_mole,i_mole,1),atype_bkghm_decomp(i_mole,i_mole,2),atype_bkghm_decomp(i_mole,i_mole,3),atype_bkghm_decomp(i_mole,i_mole,4),atype_lj_parameter(i_mole,i_mole,3),atype_pol(i_mole)
          enddo
       case default
          stop "C8_10_dispersion_terms setting not recognized.  Must be yes or no."
       end select

    case default
       stop "solute_cross_parameter_set setting not recognized.  Must be yes or no. "
    end Select


    ! read in explicit solute cross-term exponents and/or explicit dhf cross terms if they are there
    call read_solute_explicit_cross_exponents_and_dhf(file_h,n_type_s,explicit_cross_term_exponents,atype_bkghm_decomp,atype_lj_parameter,temp_lj_parameter, framework_simulation,lj_bkghm, solute_cross_parameter_set, ind, ind1, ind_dhf, ind_frame)






    !****************************************************************************************************************************************************************
    !**********
    !**********                                                          framework parameters
    !**********
    !****************************************************************************************************************************************************************


    Select Case(framework_simulation)
    Case(1)
       ! if ind_frame .ne. 0, then we are already at this section
       if ( ind_frame .eq. 0 ) then
          do
             Read(file_h,'(A)',Iostat=inputstatus) line
             If(inputstatus < 0) Exit
             ind=INDEX(line,'framework_species')
             IF(ind .NE. 0) Exit
          enddo
       endif

       ! get framework atom types
       Read(file_h,'(A)',Iostat=inputstatus) line
       read(file_h,*) n_type_f
       ! input is determined by whether or not this is a simulation with C8 and C10 dispersion terms
       Select Case( C8_10_dispersion_terms)
       Case('yes')
          ! read in order for each atom, chg, A_exch, A_elec, A_induc,A_dhf,A_disp,C6,C8,C10, pol
          do j_mole=1, n_type_f
             i_mole= n_type_s + j_mole
             ! make sure we are not exceeding array length
             if ( i_mole > MAX_N_ATOM_TYPE ) then
                stop "number of atom type g.t. MAX_N_ATOM_TYPE .  Increase this value"
             endif
             ! code will stop if using explicit_cross_term_exponents="no" with options other than
             ! C8_10_dispersion_terms="yes", solute_cross_parameter_set="yes"
             Select Case( explicit_cross_term_exponents )
             Case("yes")
                ! cannot use C12 coeff here
                Select Case(C12_dispersion)
                Case("yes")
                   stop "cannot use C12 dispersion with explicit_cross_term_exponents='yes' for a framework simulation"
                End Select
                read(file_h,*) atype_name(i_mole),atype_chg(i_mole),atype_bkghm_decomp(i_mole,i_mole,1),atype_bkghm_decomp(i_mole,i_mole,2),atype_bkghm_decomp(i_mole,i_mole,3),atype_bkghm_decomp(i_mole,i_mole,4),atype_bkghm_decomp(i_mole,i_mole,5),atype_lj_parameter(i_mole,i_mole,3),atype_lj_parameter(i_mole,i_mole,4),atype_lj_parameter(i_mole,i_mole,5),atype_pol(i_mole)
             Case("no")
                ! decide whether we are reading in C12 coeff
                Select Case(C12_dispersion)
                Case("no")
                   read(file_h,*) atype_name(i_mole),atype_chg(i_mole),atype_bkghm_decomp(i_mole,i_mole,1),atype_bkghm_decomp(i_mole,i_mole,2),atype_bkghm_decomp(i_mole,i_mole,3),atype_bkghm_decomp(i_mole,i_mole,4),atype_bkghm_decomp(i_mole,i_mole,5),atype_lj_parameter(i_mole,i_mole,2), atype_lj_parameter(i_mole,i_mole,3),atype_lj_parameter(i_mole,i_mole,4),atype_lj_parameter(i_mole,i_mole,5),atype_pol(i_mole)
                Case("yes")
                   read(file_h,*) atype_name(i_mole),atype_chg(i_mole),atype_bkghm_decomp(i_mole,i_mole,1),atype_bkghm_decomp(i_mole,i_mole,2),atype_bkghm_decomp(i_mole,i_mole,3),atype_bkghm_decomp(i_mole,i_mole,4),atype_bkghm_decomp(i_mole,i_mole,5),atype_lj_parameter(i_mole,i_mole,2), atype_lj_parameter(i_mole,i_mole,3),atype_lj_parameter(i_mole,i_mole,4),atype_lj_parameter(i_mole,i_mole,5),atype_lj_parameter(i_mole,i_mole,6),atype_pol(i_mole)
                End Select
             End Select
          enddo
       Case('no')
          ! read in order for each atom, chg, A_exch, A_elec, A_induc,A_dhf,C6, pol
          do j_mole=1, n_type_f
             i_mole= n_type_s + j_mole
             ! make sure we are not exceeding array length
             if ( i_mole > MAX_N_ATOM_TYPE ) then
                stop "number of atom type g.t. MAX_N_ATOM_TYPE.  Increase this value"
             endif
             read(file_h,*) atype_name(i_mole),atype_chg(i_mole),atype_bkghm_decomp(i_mole,i_mole,1),atype_bkghm_decomp(i_mole,i_mole,2),atype_bkghm_decomp(i_mole,i_mole,3),atype_bkghm_decomp(i_mole,i_mole,4),atype_lj_parameter(i_mole,i_mole,3),atype_pol(i_mole)
          enddo
       end select

       n_atom_type =n_type_f + n_type_s

       ! if reading explicit cross term exponents
       Select Case( explicit_cross_term_exponents )
       Case("yes")
          do
             Read(file_h,'(A)',Iostat=inputstatus) line
             If(inputstatus < 0) Exit
             ind=INDEX(line,'solute framework')
             IF(ind .NE. 0) Exit
          enddo

          ! get solute-framework exponents
          n_cross = n_type_s * n_type_f
          allocate( temp_cross( n_cross ) )
          read(file_h, *) (temp_cross(ind) , ind=1,n_cross)
          ind =1
          do i_mole=1,n_type_s
             do ind1=1, n_type_f
                j_mole= n_type_s + ind1
                atype_lj_parameter(i_mole,j_mole,2) = temp_cross(ind)
                ind = ind + 1
             enddo
          enddo
          deallocate( temp_cross )

       End Select

    Case(0)  ! no framework here
       n_atom_type =n_type_s

    End Select

    !************************************************************************

    ! fill in atype_solute array
    do i_mole=1,n_atom_type
       if ( i_mole .le. n_type_s ) then
          atype_solute(i_mole) = 1
       else
          atype_solute(i_mole) = 0
       endif
    enddo

    ! set up lj_bkghm_index array to tell code which atom-pairs use lj or bkghm potentials
    Select Case(lj_bkghm)
    Case(3)
       ! solute-solute interactions are lj, solute-framework are sapt based
       do i_type=1,n_atom_type
          do j_type=1,n_atom_type
             if ( ( i_type .le. n_type_s ) .and. ( j_type .le. n_type_s ) ) then
                lj_bkghm_index(i_type,j_type) = 2
             else
                lj_bkghm_index(i_type,j_type) = 1
             endif
          end do
       end do

    case default
       ! here we only have one type of potential for everything
       lj_bkghm_index = lj_bkghm
    End Select



    ! if using additional penetration force field, read it in
    Select Case(penetration_force_field)
    Case("yes")
       close(file_h)
       open(unit=file_h,file=ifile_pmt,status="old")
       do
          Read(file_h,'(A)',Iostat=inputstatus) line
          If(inputstatus < 0) Exit
          ind=INDEX(line,'penetration force field')
          IF(ind .NE. 0) Exit
       enddo
       if ( ind == 0 ) then
          stop " cant find penetration force field section in parameter file "
       endif
       read(file_h,*) penetration_force_field_threshold
       read(file_h,*) n_penetration
       ! set switch distances to zero for all pairs not used
       atype_penetration_ff(:,:,2) = 0d0
       do i =1,n_penetration
          read(file_h,*) i_type, j_type, (atype_penetration_ff(i_type,j_type,j),j=1,4)
          atype_penetration_ff(j_type,i_type,:)=atype_penetration_ff(i_type,j_type,:)
       enddo
    end Select


    close(file_h)

    ! if using three body dispersion, get C9 coefficients
    Select Case(three_body_dispersion)
    Case("yes")
       call read_C9_three_body_dispersion_parameters(ifile_pmt)
    End Select

    ! if using three body exchange, get those parameters
    Select Case(three_body_exchange)
    Case("yes")
       call read_three_body_exchange_parameters(ifile_pmt)
    end Select

    ! finally make sure we don't have two of the same atom type defined

    do i_mole=1, n_atom_type - 1
       do j_mole = i_mole + 1 , n_atom_type
          if ( atype_name(i_mole) .eq. atype_name(j_mole) ) then
             write(*,*) "atomic parameters defined more than once for atom type", atype_name(i_mole)
             stop
          endif
       enddo
    enddo


  end subroutine read_param_decomp


  !**************************************
  ! this subroutine generates parameters for a simulation using sapt-based force field
  ! temporary array temp_lj_parameter is used to input solute-solute specific parameters (C6,etc.),
  ! whereas atype_lj_parameter array has solute parameters for solute framework interactions.
  ! as all diagonal and cross terms will be generated and filled into atype_lj_parameter array,
  ! temp_lj_parameter array will not be used after this subroutine.
  !
  ! note for C6,C8,C10 cross terms, if only nonzero C6 terms on hydrogen, C8 cross terms will be generated between
  ! hydrogen and heavy atoms using C8cross=sqrt(C6H*C10nonH).  This is appropriate as C6 is purely dipole polarizability
  ! and C10 is purely quadrupolar polarizability and so C8 will be dipole-quadrupole polarizability 
  !**************************************
subroutine   gen_param_decomp(alist,tot_n_mole,n_atom, chg, chg_drude , temp_lj_parameter , n_type_s)
  use global_variables
  character(*), dimension(:,:),intent(in) :: alist
  real*8, dimension(:,:),intent(out) :: chg, chg_drude
  integer,dimension(:),intent(in):: n_atom
  integer,intent(in)::tot_n_mole,n_type_s
  real*8,dimension(:,:,:),intent(inout) :: temp_lj_parameter

  integer:: i,i_param, j_param, i_mole , i_atom, ind, i_type,flag, lj_bkghm_local
  real*8 :: i_sign,j_sign
  real*8,parameter:: small = 1D-10, small_dhf=1D-2

  ! there's a chance here that lj_bkghm = 3 for a hybrid lj/bkghm force field.  The subroutine combination_rule_cross_terms only recognizes settings of lj_bkghm=1 or 2.
  ! since this subroutine only calles combination_rule_cross_terms subroutine for solute-framework interactions, and solute-framework interactions are always buckingham type
  ! (lj_bkghm=1), create a local variable lj_bkghm_local=1 and use this in the subroutine call
  Select Case(lj_bkghm)
  Case(1)
     continue
  Case(3)
     continue
  Case default
     stop "cannot recognize the setting of variable lj_bkghm in subroutine gen_param_decomp"
  End Select
  lj_bkghm_local=1


  !**********************************************************************************************************
  !
  !                 create cross terms for decomposed buckingham parameters
  !
  !**********************************************************************************************************

  ! this is  bkghm potential, use all geometric combination rules
  ! generate bkghm coeff cross terms, 1st index is exchange and is positive, 2nd,3rd,and 5th indices are electrostatic ,induction, and dispersion, and are negative, 4th index is dhf and can be
  ! either positive or negative
  do i_param=1, n_atom_type
     do j_param=1,n_atom_type
        do ind =1,5
           Select case(ind)
           case(1)
              atype_bkghm_decomp(i_param,j_param,ind) = sqrt(abs( atype_bkghm_decomp(i_param,i_param,ind)* atype_bkghm_decomp(j_param,j_param,ind)))
           case(2:3)
              atype_bkghm_decomp(i_param,j_param,ind) = -sqrt(abs( atype_bkghm_decomp(i_param,i_param,ind)* atype_bkghm_decomp(j_param,j_param,ind)))
           case(4)
              ! if cross terms have not been read in
              if ( abs(atype_bkghm_decomp(i_param,j_param,ind)) < small_dhf ) then
                 ! don't use combination rule for same atom, otherwise might mess up sign
                 if (i_param .eq. j_param) then
                    atype_bkghm_decomp(i_param,j_param,ind) = atype_bkghm_decomp(i_param,j_param,ind)
                 else
                    call construct_dhf_cross_terms(atype_bkghm_decomp(i_param,j_param,ind),atype_bkghm_decomp(i_param,i_param,ind),atype_bkghm_decomp(j_param,j_param,ind),dhf_combination_rule)
                 endif
              endif
           case(5)           
              atype_bkghm_decomp(i_param,j_param,ind) = -sqrt(abs( atype_bkghm_decomp(i_param,i_param,ind)* atype_bkghm_decomp(j_param,j_param,ind)))
           end Select
        enddo

     enddo
  enddo


  ! first fill in solute-framework parameters, as these are already partially stored in atype_lj_parameter array
  do i_param=1, n_atom_type
     do j_param=n_type_s + 1,n_atom_type
        call solute_framework_cross_terms(i_param,j_param,atype_bkghm_decomp,atype_lj_parameter,atype_name,C8_10_dispersion_terms,n_type_s,C12_dispersion)
     enddo
  enddo


  ! create exponent cross terms if combination rule is requested
  Select Case(explicit_cross_term_exponents)
  Case("no")
     do i_param=1, n_atom_type
        do j_param=n_type_s + 1,n_atom_type     
           ! see above comment as to setting of lj_bkghm_local
           call combination_rule_cross_terms(atype_lj_parameter,i_param,j_param,lj_bkghm_local,lj_comb_rule,1)
        enddo
     enddo
     ! construct solute-solute exponents if we are NOT using a hybrid lj/bkghm force field (i.e. lj_bkghm != 3 )
     Select Case(lj_bkghm)
     Case(3)
        continue  ! do nothing
     case default
        ! if we are using a single set of solute exponents for all cross terms, then temp_lj_parameter(1,1,9) will be set to -101 to signal this
        if ( temp_lj_parameter(1,1,9) > -100d0) then
           ! now fill solute-solute_exponents back in.  Remember these are stored in temp_lj_parameter(:,:,9)
           do i_param=1,n_type_s
              do j_param=i_param, n_type_s
                 atype_lj_parameter(i_param,j_param,2) = temp_lj_parameter(i_param,j_param,9)
              enddo
           enddo
           ! if we need to create solute-solute cross exponents
        else
           ! note that here, through previous consistency checks, we are sure that the parameter solute_cross_parameter_set='no'
           do i_param=1,n_type_s
              do j_param=i_param, n_type_s
                 ! see above comment as to setting of lj_bkghm_local
                 call combination_rule_cross_terms(atype_lj_parameter,i_param,j_param,lj_bkghm_local,lj_comb_rule,1) 
              enddo
           enddo
        endif
     End Select
  End Select

  ! now fill in transpose framework-solute terms
  do i_param= n_type_s + 1, n_atom_type
     do j_param= 1, n_atom_type
        atype_lj_parameter(i_param,j_param,:) = atype_lj_parameter(j_param,i_param,:)
     enddo
  enddo

  ! it is important that solute-solute cross terms are formed AFTER solute-framework cross terms, because atype_lj_parameter array will be filled in with solute-solute cross terms
  ! in appropriate locations, whereas before these diagonal indices held parameters used to generate solute-framework cross terms
  ! now solute-solute
  do i_param=1, n_type_s
     do j_param=1, n_type_s
        Select Case(lj_bkghm)
        Case(3)
           ! this is a hybrid lj/bkghm potential, need lj_comb_rule2 combination rule for solute-solute cross terms
           call solute_solute_cross_terms(i_param,j_param,lj_bkghm,atype_bkghm_decomp,atype_lj_parameter,atype_name,C8_10_dispersion_terms,temp_lj_parameter,solute_cross_parameter_set,small_dhf,dhf_combination_rule,C12_dispersion,lj_comb_rule2)
        Case default
           call solute_solute_cross_terms(i_param,j_param,lj_bkghm,atype_bkghm_decomp,atype_lj_parameter,atype_name,C8_10_dispersion_terms,temp_lj_parameter,solute_cross_parameter_set,small_dhf,dhf_combination_rule,C12_dispersion)
        End Select
     enddo
  enddo

  ! if we were using a hybrid lj/bkghm force field, lj_bkghm=3, must use Lorentz-Bethelot combination rule
  ! now we need to create C12_C6 cross terms out of epsilon and sigma
  Select Case(lj_bkghm)
  Case(3)
     do i_param=1,n_type_s
        do j_param=1,n_type_s
           call gen_C12_C6_epsilon_sigma(atype_lj_parameter,i_param,j_param)
        enddo
     enddo

  case default

     ! exponential cross terms have already been read in,just fill in transpose terms
     do i_param=1, n_atom_type
        do j_param=i_param+1,n_atom_type
           atype_lj_parameter(j_param,i_param,2) = atype_lj_parameter(i_param,j_param,2)
        enddo
     enddo

  End Select


  ! now create atom index array that links atoms to parameters
  do i_mole = 1, tot_n_mole
     do i_atom =1, n_atom(i_mole)
        flag=0
        do  i_type=1, n_atom_type
           if (atype_name(i_type) .eq. alist(i_mole,i_atom) ) then
              flag=1
              i_param = i_type
              exit 
           endif
        enddo
        ! make sure all atom types have parameters
        if (flag .eq. 0 ) then
           write(*,*)  "atom type    ", alist(i_mole,i_atom), "doesn't have force field parameters!"
           stop
        endif

        ! set index, chg, and polarizability
        atom_index(i_mole,i_atom) = i_param
        chg(i_mole,i_atom) = atype_chg(i_param)
        chg_drude(i_mole,i_atom) = -sqrt(atype_pol(i_param)*springcon)

        ! make sure all polarizabilities are zero if drude_simulation = 0 , so that input is consistent
        Select Case(drude_simulation)
        Case(0)
           if ( abs( chg_drude(i_mole,i_atom) ) > small ) then
              write(*,*)  ""
              write(*,*)  " polarizability on atom ", i_atom , " on molecule ", i_mole
              write(*,*)  " is non-zero, yet drude_simulation is set to 0 "
              write(*,*)  " if you want to run a drude oscillator simulation, set "
              write(*,*)  " drude_simulation to 1 ! "
              write(*,*) ""
              stop
           endif
        End Select
     enddo
  enddo


end subroutine gen_param_decomp


  !***********************************************************
  !  this subroutine generates cross_terms for solute molecules, as these could be different from solute-framework terms
  !  parameters for this are stored in temp_lj_parameter array if C8,C10 terms are used for solute-framework
  !  if only C6 terms are used for solute-framework, then solute-solute terms are the same and are stored in atype_lj_parameter array
  !  only C6 terms used for solute
  !
  ! important, note that if we are reading in two complete sets of solute parameters, then 8 indices of temp_lj_parameter(:,:,i) array will be used
  ! the first 5 for bkghm coefficients, and the last three for C6,C8,C10
  ! however, if we are only reading in a different set of dispersion terms for the solute, then only 3 indices of temp_lj_parameter(:,:,i) array will be used
  ! with indices 3-5 used for C6-C10 respectively
  !
  ! NOTE that if lj_bkghm = 3, the solute-solute cross terms are a lj potential.
  !***********************************************************
subroutine solute_solute_cross_terms(i_param,j_param,lj_bkghm,atype_bkghm_decomp,atype_lj_parameter,atype_name,C8_10_dispersion_terms,temp_lj_parameter,solute_cross_parameter_set,small_dhf,dhf_combination_rule,C12_dispersion,lj_comb_rule2)
  integer,intent(in) :: i_param,j_param, lj_bkghm
  real*8,dimension(:,:,:),intent(inout) :: atype_bkghm_decomp
  real*8,dimension(:,:,:),intent(inout) :: temp_lj_parameter
  real*8,dimension(:,:,:),intent(inout) :: atype_lj_parameter
  character(*),dimension(:),intent(in) :: atype_name
  integer,intent(in) :: dhf_combination_rule
  character(*), intent(in)  :: C8_10_dispersion_terms
  character(*),intent(in) :: solute_cross_parameter_set,C12_dispersion
  real*8 ,intent(in) :: small_dhf
  character(*),intent(in),optional :: lj_comb_rule2
  real*8,parameter :: small=1D-5

  real*8 :: i_sign,j_sign,sum
  integer :: i,ind, lj_bkghm_local

  Select Case(lj_bkghm)
  Case(1)
     ! ***************************** For lj_bkghm = 1 , this is a sapt-type force field for solute-solute interactions

     ! if we are are using two solute force fields, one for solute-solute, another for solute-framework, we have to fill in all solute-solute parameters now
     Select Case (solute_cross_parameter_set)
     Case("yes")
        !***************************************************** separate solute-framework and solute-framework solute parameters, we have some work to do....
        do ind =1,5
           Select case(ind)
           case(1)
              atype_bkghm_decomp(i_param,j_param,ind) = sqrt(abs( temp_lj_parameter(i_param,i_param,ind) * temp_lj_parameter(j_param,j_param,ind)))
           case(2:3)
              atype_bkghm_decomp(i_param,j_param,ind) = -sqrt(abs( temp_lj_parameter(i_param,i_param,ind)* temp_lj_parameter(j_param,j_param,ind)))
           case(4)
              ! if cross terms have not been read in
              if ( abs(temp_lj_parameter(i_param,j_param,ind)) < small_dhf ) then
                 ! don't use combination rule for same atom, otherwise might mess up sign
                 if (i_param .eq. j_param) then
                    atype_bkghm_decomp(i_param,j_param,ind) = temp_lj_parameter(i_param,j_param,ind)
                 else
                    call construct_dhf_cross_terms(atype_bkghm_decomp(i_param,j_param,ind),temp_lj_parameter(i_param,i_param,ind),temp_lj_parameter(j_param,j_param,ind),dhf_combination_rule)
                 endif
              Else
                 ! here we have a non-zero diagonal or cross term that has already been read in.  Just transfer arrays.
                 atype_bkghm_decomp(i_param,j_param,ind) = temp_lj_parameter(i_param,j_param,ind)
              endif
           case(5)      
              if ( (abs(temp_lj_parameter(i_param,i_param,5)) > small ) .or. (abs(temp_lj_parameter(j_param,j_param,5)) > small ) ) then
                 write(*,*) "cannot use dispersion penetration for solute-solute interaction, check atom types ", i_param, j_param
                 stop
              endif
              atype_bkghm_decomp(i_param,j_param,ind) = 0d0
           end Select
        enddo
     End Select



     !************ fill in buckingham coefficient, this will be the sum of the first 4 atype_bkghm_decomp array elements, as there should be no dispersion penetration term
     sum=0d0
     do i=1,4
        sum=sum + atype_bkghm_decomp(i_param,j_param,i)
     enddo
     atype_lj_parameter(i_param,j_param,1) = sum

     ! make sure we are not using dispersion penetration term for solute, and zero this element in atype_bkghm_decomp as solute-framework parameters have already been generated
     ! and need this for energy decomposition
     ! note this parameter occurs at different indices of temp_lj_parameter array, depending on whether we are using one or two solute force fields
     Select Case (solute_cross_parameter_set)
     Case("no")
        Select Case( flag_solute_param )
        Case(1)
           if ( (abs(temp_lj_parameter(i_param,i_param,1)) > small ) .or. (abs(temp_lj_parameter(j_param,j_param,1)) > small ) ) then
              write(*,*) "cannot use dispersion penetration for solute-solute interaction, check atom types ", i_param, j_param
              stop
           endif
        Case(0)
           if ( (abs(atype_bkghm_decomp(i_param,i_param,5)) > small ) .or. (abs(atype_bkghm_decomp(j_param,j_param,5)) > small ) ) then
              write(*,*) "cannot use dispersion penetration for solute-solute interaction, check atom types ", i_param, j_param
              stop
           endif
        End Select
     End Select

     atype_bkghm_decomp(i_param,j_param,5)=0d0


     ! now do C6,C8,C10 dispersion coefficients.  If we are using two sets of solute parameters, then these are stored in temp_lj_parameter array
     Select Case(solute_cross_parameter_set)
     Case("no")
        ! here there is only one set of solute parameters
        Select Case(C8_10_dispersion_terms)
        Case("no")
           ! only C6 terms for framework, solute has only one set of parameters and this is already in atype_lj_parameter array
           atype_lj_parameter(i_param,j_param,3) = sqrt(atype_lj_parameter(i_param,i_param,3)*atype_lj_parameter(j_param,j_param,3))
           atype_lj_parameter(i_param,j_param,4) = 0d0
           atype_lj_parameter(i_param,j_param,5) = 0d0
        Case("yes")

           ! we need to determine where solute parameters are.  If 2nd set of dispersion parameters were read in, then use temp_lj_parameter array
           Select Case( flag_solute_param )
           Case(1)
              ! here there is no solute_cross_parameter_set, but flag_solute_param=1 means we have different solute terms for dispersion.
              ! these have been stored in temp_lj_parameter(:,:,i) array in indices 3-5 for C6-C10 respectively
              ! note no C12 terms should be used here, and code should have already stopped if C12_dispersion='yes' was set

              atype_lj_parameter(i_param,i_param,3) = temp_lj_parameter(i_param,i_param,3) ;
              atype_lj_parameter(i_param,i_param,4) = temp_lj_parameter(i_param,i_param,4) ;
              atype_lj_parameter(i_param,i_param,5) = temp_lj_parameter(i_param,i_param,5) ;
              atype_lj_parameter(j_param,j_param,3) = temp_lj_parameter(j_param,j_param,3) ;
              atype_lj_parameter(j_param,j_param,4) = temp_lj_parameter(j_param,j_param,4) ;
              atype_lj_parameter(j_param,j_param,5) = temp_lj_parameter(j_param,j_param,5) ;

              call gen_C6_C10_cross_terms(i_param,j_param,atype_name,atype_lj_parameter,C12_dispersion )

           Case(0)
              ! here no additional dispersion terms were read in.  use atype_lj_parameter array
              call gen_C6_C10_cross_terms(i_param,j_param,atype_name,atype_lj_parameter,C12_dispersion )

           End select

        end select

     Case("yes")
        ! here there are two sets of solute parameters
        ! note here that C6,C8,C10 are stored in 6,7,8 indices of temp_lj_parameter array since we have two solute parameter sets
        ! note that C12 terms cannot be used here for solute-solute interactions, this should have been warned if C12_dispersion='yes' was set
        ! make sure to zero these components as they will have nonzero values from framework cross terms

        atype_lj_parameter(i_param,i_param,3) = temp_lj_parameter(i_param,i_param,6) ;
        atype_lj_parameter(i_param,i_param,4) = temp_lj_parameter(i_param,i_param,7) ;
        atype_lj_parameter(i_param,i_param,5) = temp_lj_parameter(i_param,i_param,8) ;
        atype_lj_parameter(i_param,i_param,6) = 0d0
        atype_lj_parameter(j_param,j_param,3) = temp_lj_parameter(j_param,j_param,6) ;
        atype_lj_parameter(j_param,j_param,4) = temp_lj_parameter(j_param,j_param,7) ;
        atype_lj_parameter(j_param,j_param,5) = temp_lj_parameter(j_param,j_param,8) ;
        atype_lj_parameter(j_param,j_param,6) = 0d0

        call gen_C6_C10_cross_terms(i_param,j_param,atype_name,atype_lj_parameter,C12_dispersion )

     End Select



  Case(3)
     !****************************** with the setting lj_bkghm=3, this is a hybrid lj/bkghm simulation **********************
     ! there are automatically two different parameter sets, since there are two different force fields.  Therefore, copy parameters from temp_lj_parameter array
     atype_lj_parameter(i_param,i_param,1) = temp_lj_parameter(i_param,i_param,1)
     atype_lj_parameter(i_param,i_param,2) = temp_lj_parameter(i_param,i_param,2)
     atype_lj_parameter(i_param,i_param,3) = temp_lj_parameter(i_param,i_param,3)
     atype_lj_parameter(j_param,j_param,1) = temp_lj_parameter(j_param,j_param,1)
     atype_lj_parameter(j_param,j_param,2) = temp_lj_parameter(j_param,j_param,2)
     atype_lj_parameter(j_param,j_param,3) = temp_lj_parameter(j_param,j_param,3)

     ! lj_comb_rule2 must be present in this case
     if ( present(lj_comb_rule2) ) then
        continue
     else
        stop " lj_comb_rule2 argument must be present if lj_bkghm='yes' in call to solute_solute_cross_terms subroutine"
     endif


     ! subroutine combination_rule_cross_terms only recognizes setting of lj_bkghm = 1 or 2.  In this case for lj_bkghm = 3, the solute
     ! solute portion of the force field is a lj type force field, so set lj_bkghm_local =2 and pass that variable instead
     lj_bkghm_local = 2
     call combination_rule_cross_terms(atype_lj_parameter,i_param,j_param,lj_bkghm_local,lj_comb_rule2) 

  case default
     stop "cannot recognize the setting of variable lj_bkghm in subroutine solute_solute_cross_terms"
  End Select


end subroutine solute_solute_cross_terms


  !***********************************************************
  ! this subroutine generates cross terms for solute - framework atoms
  ! the only tricky part here is generating C6,C8,C10 terms if necessary
  ! see comments about gen_param_decomp subroutine, but in short,
  ! for Hydrogen-heavy atom cross terms, we generate C8 with sqrt (C6H*C10) if no C8,C10 are given (zero)
  !***********************************************************
  subroutine solute_framework_cross_terms(i_param,j_param,atype_bkghm_decomp,atype_lj_parameter,atype_name,C8_10_dispersion_terms,n_type_s,C12_dispersion)
    integer,intent(in) :: i_param,j_param,n_type_s
    real*8,dimension(:,:,:),intent(in) :: atype_bkghm_decomp
    real*8,dimension(:,:,:),intent(inout) :: atype_lj_parameter
    character(*),dimension(:),intent(in) :: atype_name
    character(*), intent(in)  :: C8_10_dispersion_terms, C12_dispersion

    character(len(atype_name)) :: label_i, label_j
    integer :: flag_Hi, flag_Hj

    ! first fill in buckingham coefficient, which will be the sum of the 5 atype_bkghm_decomp array elements
    ! if no dispersion penetration term, this term will already be zeroed

    atype_lj_parameter(i_param,j_param,1) = sum( atype_bkghm_decomp(i_param,j_param,:))

    ! check that the net buckingham coefficient is repulsive
    ! diagonal framework elements could possibly be negative, without creating negative cross terms, so don't check framework-framework terms
    if ( (i_param .le. n_type_s ) .or. (j_param .le. n_type_s ) ) then
       if ( atype_lj_parameter(i_param,j_param,1) < 0d0 ) then
          write(*,*) "negative net buckingham coefficient between atom ", atype_name(i_param) , " and atom ", atype_name(j_param) , "this is unphysical"
          stop
       endif
    endif

    Select Case(C8_10_dispersion_terms)
    Case("no")
       ! only C6 terms
       atype_lj_parameter(i_param,j_param,3) = sqrt(atype_lj_parameter(i_param,i_param,3)*atype_lj_parameter(j_param,j_param,3))
       atype_lj_parameter(i_param,j_param,4) = 0d0
       atype_lj_parameter(i_param,j_param,5) = 0d0

    Case("yes")
       ! C6,C8,C10  dispersion terms
       ! first make sure C6,C8,C10 are all positive, as this is not necessarily the case in camcasp calculations
       if ( (atype_lj_parameter(i_param,i_param,3) < 0d0 ) .or. (atype_lj_parameter(i_param,i_param,4) < 0d0 ) .or. (atype_lj_parameter(i_param,i_param,5) < 0d0 ) ) then
          write(*,*) "C6,C8,C10 terms are not all positive for atom ", atype_name(i_param)
          stop
       endif
       if ( (atype_lj_parameter(j_param,j_param,3) < 0d0 ) .or. (atype_lj_parameter(j_param,j_param,4) < 0d0 ) .or. (atype_lj_parameter(j_param,j_param,5) < 0d0 ) ) then
          write(*,*) "C6,C8,C10 terms are not all positive for atom ", atype_name(j_param)
          stop
       endif

       call gen_C6_C10_cross_terms(i_param,j_param,atype_name,atype_lj_parameter,C12_dispersion )


    case default                               ! this is for Case structure on C8_10_dispersion_terms
       stop " C8_10_dispersion_terms option not recognized.  Must be yes/no "
    end Select

  end  subroutine solute_framework_cross_terms




  !***************************************************************
  !  This subroutine is for creating C6, C8, C10 cross terms
  !  we use a subroutine for this because in general we use
  !  a different combination rule for hydrogen, generating C8 cross terms
  ! from hydrogen C6 and a heavy atom C10
  !
  !****************************************************************
  subroutine gen_C6_C10_cross_terms(i_param,j_param,atype_name,atype_lj_parameter, C12_dispersion )
    integer,intent(in) :: i_param,j_param
    real*8,dimension(:,:,:),intent(inout) :: atype_lj_parameter
    character(*),dimension(:),intent(in) :: atype_name
    character(*), intent(in) :: C12_dispersion

    character(len(atype_name)) :: label_i, label_j
    integer :: flag_Hi, flag_Hj
    real*8,parameter :: small=1D-5
    integer :: flag_warn=0

    ! if we are using C12 coefficients, only use simple combination rule regardless of hydrogen
    Select Case(C12_dispersion)
    Case("yes")
       atype_lj_parameter(i_param,j_param,3) = sqrt(atype_lj_parameter(i_param,i_param,3)*atype_lj_parameter(j_param,j_param,3))
       atype_lj_parameter(i_param,j_param,4) = sqrt(atype_lj_parameter(i_param,i_param,4)*atype_lj_parameter(j_param,j_param,4))
       atype_lj_parameter(i_param,j_param,5) = sqrt(atype_lj_parameter(i_param,i_param,5)*atype_lj_parameter(j_param,j_param,5))
       atype_lj_parameter(i_param,j_param,6) = sqrt(atype_lj_parameter(i_param,i_param,6)*atype_lj_parameter(j_param,j_param,6))
    Case("no")

       ! we need to figure out whether we have hydrogen atoms here, since these probably only have C6 coefficients
       label_i = atype_name(i_param)
       label_j = atype_name(j_param)

       if( (label_i(1:1) .eq. "H") .or. (label_i(1:1) .eq. "h") ) then
          flag_Hi = 1
       else
          flag_Hi = 0
       endif
       if( (label_j(1:1) .eq. "H") .or. (label_j(1:1) .eq. "h") ) then
          flag_Hj = 1
       else
          flag_Hj = 0
       endif

       Select Case(flag_Hi)
       Case(0)
          Select Case(flag_Hj)
          Case(0)
             ! two non hydrogen atoms
             atype_lj_parameter(i_param,j_param,3) = sqrt(atype_lj_parameter(i_param,i_param,3)*atype_lj_parameter(j_param,j_param,3))
             atype_lj_parameter(i_param,j_param,4) = sqrt(atype_lj_parameter(i_param,i_param,4)*atype_lj_parameter(j_param,j_param,4))
             atype_lj_parameter(i_param,j_param,5) = sqrt(atype_lj_parameter(i_param,i_param,5)*atype_lj_parameter(j_param,j_param,5))
          Case(1)
             ! one hydrogen atom, use regular combination rule if C8,C10 terms on hydrogen are not zero
             if ( (atype_lj_parameter(j_param,j_param,4) > small) .or. (atype_lj_parameter(j_param,j_param,5) > small) ) then
                atype_lj_parameter(i_param,j_param,3) = sqrt(atype_lj_parameter(i_param,i_param,3)*atype_lj_parameter(j_param,j_param,3))
                atype_lj_parameter(i_param,j_param,4) = sqrt(atype_lj_parameter(i_param,i_param,4)*atype_lj_parameter(j_param,j_param,4))
                atype_lj_parameter(i_param,j_param,5) = sqrt(atype_lj_parameter(i_param,i_param,5)*atype_lj_parameter(j_param,j_param,5))
             else
                ! use special combination rule for C8
                atype_lj_parameter(i_param,j_param,3) = sqrt(atype_lj_parameter(i_param,i_param,3)*atype_lj_parameter(j_param,j_param,3))    
                atype_lj_parameter(i_param,j_param,4) = sqrt(atype_lj_parameter(i_param,i_param,5)*atype_lj_parameter(j_param,j_param,3))
                atype_lj_parameter(i_param,j_param,5) = 0d0

                ! warn about this combination rule
                Select Case( flag_warn )
                Case(0)
                   write(*,*) "IMPORTANT NOTE:  Combination rule C8=sqrt(C6*C10) is being used for hydrogen-heavy atom combination"
                   flag_warn=1
                End Select

             endif
          End Select
       Case(1)
          Select Case(flag_Hj)
          Case(0)                      
             ! one hydrogen atom, use regular combination rule if C8,C10 terms on hydrogen are not zero
             if ( (atype_lj_parameter(i_param,i_param,4) > small) .or. (atype_lj_parameter(i_param,i_param,5) > small) ) then
                atype_lj_parameter(i_param,j_param,3) = sqrt(atype_lj_parameter(i_param,i_param,3)*atype_lj_parameter(j_param,j_param,3))
                atype_lj_parameter(i_param,j_param,4) = sqrt(atype_lj_parameter(i_param,i_param,4)*atype_lj_parameter(j_param,j_param,4))
                atype_lj_parameter(i_param,j_param,5) = sqrt(atype_lj_parameter(i_param,i_param,5)*atype_lj_parameter(j_param,j_param,5))
             else
                ! use special combination rule for C8
                atype_lj_parameter(i_param,j_param,3) = sqrt(atype_lj_parameter(i_param,i_param,3)*atype_lj_parameter(j_param,j_param,3))    
                atype_lj_parameter(i_param,j_param,4) = sqrt(atype_lj_parameter(i_param,i_param,3)*atype_lj_parameter(j_param,j_param,5))
                atype_lj_parameter(i_param,j_param,5) = 0d0

                ! warn about this combination rule
                Select Case( flag_warn )
                Case(0)
                   write(*,*) "IMPORTANT NOTE:  Combination rule C8=sqrt(C6*C10) is being used for hydrogen-heavy atom combination"
                   flag_warn=1
                End Select

             endif
          Case(1)
             ! both atoms are hydrogen, use regular combination rule
             atype_lj_parameter(i_param,j_param,3) = sqrt(atype_lj_parameter(i_param,i_param,3)*atype_lj_parameter(j_param,j_param,3))
             atype_lj_parameter(i_param,j_param,4) = sqrt(atype_lj_parameter(i_param,i_param,4)*atype_lj_parameter(j_param,j_param,4))
             atype_lj_parameter(i_param,j_param,5) = sqrt(atype_lj_parameter(i_param,i_param,5)*atype_lj_parameter(j_param,j_param,5))
          End Select
       End Select

    End Select

  end subroutine gen_C6_C10_cross_terms


  !***************************************************
  ! this subroutine generates cross terms for buckingham or lennard jones
  ! force fields using the desired combination rules
  !
  ! for Lennard Jones force field,
  ! Lorentz-Berthelot combination rules operate on sigma and epsilon,
  ! Opls combination rules operate on C6, C12
  !***************************************************
  subroutine combination_rule_cross_terms(atype_lj_parameter,i_param,j_param,lj_bkghm,lj_comb_rule,exp_only)
    real*8, dimension(:,:,:), intent(inout) :: atype_lj_parameter
    integer, intent(in) :: i_param,j_param
    integer, intent(in) :: lj_bkghm
    character(*), intent(in) :: lj_comb_rule
    integer, intent(in),optional :: exp_only

    real*8 :: e1, e2
    integer :: flag_exponent=0

    ! if buckingham force field, we may want separate subroutine to generate cross terms for coefficients if we are using an energy decomposed force field.  We may only use this routine to generate exponents
    Select Case(lj_bkghm)
    Case(1)
       if ( present(exp_only) ) then
          flag_exponent = 1
       endif
    End Select

    Select Case(lj_bkghm)
    Case(1)
       !************* buckingham force field
       Select Case(lj_comb_rule)
       Case("standard")
          ! all geometric combination rules
          atype_lj_parameter(i_param,j_param,1) = sqrt(atype_lj_parameter(i_param,i_param,1)*atype_lj_parameter(j_param,j_param,1))    ! A coefficient
          atype_lj_parameter(i_param,j_param,2) = sqrt(atype_lj_parameter(i_param,i_param,2)*atype_lj_parameter(j_param,j_param,2))    ! B parameter
          atype_lj_parameter(i_param,j_param,3) = sqrt(atype_lj_parameter(i_param,i_param,3)*atype_lj_parameter(j_param,j_param,3))    ! C parameter
       Case("ZIFFF")
          ! only do prefactors if exp_only input isn't present
          if ( flag_exponent .eq. 0 ) then
             atype_lj_parameter(i_param,j_param,1) = sqrt(atype_lj_parameter(i_param,i_param,1)*atype_lj_parameter(j_param,j_param,1))    ! A coefficient 
          endif
          ! exponent combination rules given by ZIF FF rules
          ! B exponent
          e1 = atype_lj_parameter(i_param,i_param,2)
          e2 = atype_lj_parameter(j_param,j_param,2)
          atype_lj_parameter(i_param,j_param,2) = (e1+e2) * (e1*e2)/(e1**2+e2**2) 
          atype_lj_parameter(i_param,j_param,3) = sqrt(atype_lj_parameter(i_param,i_param,3)*atype_lj_parameter(j_param,j_param,3))    ! C coefficient
       Case default
          stop "lj_comb_rule parameter isn't recognized for buckingham type force field.  Please use either 'standard' for geometric combination rules, or 'ZIFFF' for ZIF FF combination rules"
       End Select
    Case(2)
       !************* lennard jones force field
       Select Case(lj_comb_rule)
       Case("standard")
          ! Lorentz_Berthelot combination rules
          atype_lj_parameter(i_param,j_param,1) = sqrt(atype_lj_parameter(i_param,i_param,1)*atype_lj_parameter(j_param,j_param,1))    ! epsilon
          atype_lj_parameter(i_param,j_param,2) = (atype_lj_parameter(i_param,i_param,2) + atype_lj_parameter(j_param,j_param,2))/dble(2)   ! sigma
       Case("opls")
          ! opls combination rules
          atype_lj_parameter(i_param,j_param,1) = sqrt(atype_lj_parameter(i_param,i_param,1)*atype_lj_parameter(j_param,j_param,1))   ! C12        
          atype_lj_parameter(i_param,j_param,2) = sqrt(atype_lj_parameter(i_param,i_param,2)*atype_lj_parameter(j_param,j_param,2))   ! C6
       Case default
          stop "lj_comb_rule parameter isn't recognized for lennard jones force field.  Please use either 'standard' for Lorentz-Berthelot combination rules or 'opls' for opls combination rules."
       End Select
    End Select


  end subroutine combination_rule_cross_terms


  !**************************************************
  ! this subroutine creates C12 and C6 coefficients from epsilon and sigma
  ! parameters for a lennard jones force field
  !**************************************************
  subroutine gen_C12_C6_epsilon_sigma(atype_lj_parameter,i_param,j_param)
    real*8, dimension(:,:,:), intent(inout) :: atype_lj_parameter
    integer, intent(in)  :: i_param,j_param

    real*8 :: epsilon, sigma

    epsilon = atype_lj_parameter(i_param,j_param,1)
    sigma   = atype_lj_parameter(i_param,j_param,2)

    ! C12
    atype_lj_parameter(i_param,j_param,1) = 4d0*epsilon*sigma**12
    ! C6
    atype_lj_parameter(i_param,j_param,2) = 4d0*epsilon*sigma**6

  end subroutine gen_C12_C6_epsilon_sigma




  !*********************************************
  ! this subroutine constructs dhf cross terms based on desired combination rule
  !  combination rule "1" is as in JPC C, 2012, 116, 1892-1903
  !  combination rule "2" is Aij=positive iff Aii and Ajj positive
  !  combination rule "3" is Aij always negative
  !*********************************************
  subroutine construct_dhf_cross_terms(Aij,Aii,Ajj,dhf_comb_rule)
    real*8,intent(inout) :: Aij,Aii,Ajj
    integer, intent(in) :: dhf_comb_rule

    real*8,parameter:: small = 1D-10
    real*8 :: i_sign,j_sign

    ! if zero, dont waste time, and dont divide by zero

    if ( (abs(Aii)>small) .and. (abs(Ajj)>small) ) then
       Select Case (dhf_comb_rule)
       Case(1)
          i_sign = Aii/abs(Aii)
          j_sign = Ajj/abs(Ajj)
          Aij = i_sign * j_sign * dsqrt(abs(Aii)*abs(Ajj))
       Case(2)
          if ( ( Aii > 0d0 ) .and. ( Ajj > 0d0 ) ) then
             Aij = dsqrt(abs(Aii)*abs(Ajj))    
          else
             Aij = -dsqrt(abs(Aii)*abs(Ajj))   
          endif
       Case(3)
          Aij = -dsqrt(abs(Aii)*abs(Ajj))   
       End Select
    else
       Aij=0d0
    endif

  end subroutine construct_dhf_cross_terms




  !*********************************************
  ! this subroutine reads in force field parameters for a sapt-type force field, if there is only one set of 
  ! solute parameters being used, (i.e. we are not using a different set of force field parameters for solute-framework
  ! interactions, if this latter case is true, use subroutine read_solute_parameters_two_set )
  !*********************************************
  subroutine read_solute_parameters_one_set(file_h,n_type_s,C12_dispersion,explicit_cross_term_exponents,atype_name,atype_chg,atype_bkghm_decomp,atype_lj_parameter,atype_pol)
    integer,intent(in)::file_h,n_type_s
    character(*),intent(in) :: C12_dispersion,explicit_cross_term_exponents
    character(*),dimension(:),intent(inout) :: atype_name
    real*8,dimension(:),intent(inout) :: atype_chg,atype_pol
    real*8,dimension(:,:,:),intent(inout) :: atype_bkghm_decomp
    real*8,dimension(:,:,:),intent(inout) :: atype_lj_parameter

    integer :: i_mole,nargs
    integer,parameter :: max_param=20
    character(300) :: input_string
    character(50),dimension(max_param)::args

    ! input format is dictated by how many parameters are given for each atom type
    read(file_h,'(A300)') input_string

    call parse(input_string," ",args,nargs)

    backspace(file_h)

    ! this subroutine only gets called if C8_C10_dispersion='yes', and solute_cross_parameter_set='no'

    ! The arguments we are expecting in the parameter file are, in order
    ! atomname,  charge, Aexch, Aelec, Aind, Adhf, Adisp, B(maybe), C6, C8, C10, C12(maybe), alpha

    ! so if 13 arguments present, assume everything above is there, in that order
    ! if 11 arguments only, assume B, C12 are not present (ie, not using C12 for solute-solute, using explicit exponents for solute-solute
    ! if 12 arguments, determine whether B or C12 is present based on setting of C12_dispersion

    ! if only 11 arguments, exponents are not there, and there is only one solute parameter set, this means that explicit_cross_term_exponents must be set to 'yes'
    Select Case(nargs)
    Case(11)
       Select Case(explicit_cross_term_exponents)
       Case("no")
          stop "explicit_cross_term_exponents was set to 'no', but solute exponents were not found in solute parameter section!"
       End Select
    Case(12)
       ! if 12 arguments, check for the same problem
       Select Case(C12_dispersion)
       Case("yes")
          Select Case(explicit_cross_term_exponents)
          Case("no")
             stop "explicit_cross_term_exponents was set to 'no', but solute exponents were not found in solute parameter section!"
          End Select
       End Select
    End Select


    do i_mole=1, n_type_s
       Select Case(nargs)
       Case(13)
          ! read  atomname,  charge, Aexch, Aelec, Aind, Adhf, Adisp, B, C6, C8, C10, C12, alpha
          read(file_h,*) atype_name(i_mole),atype_chg(i_mole),atype_bkghm_decomp(i_mole,i_mole,1),atype_bkghm_decomp(i_mole,i_mole,2),atype_bkghm_decomp(i_mole,i_mole,3),atype_bkghm_decomp(i_mole,i_mole,4),atype_bkghm_decomp(i_mole,i_mole,5),atype_lj_parameter(i_mole,i_mole,2),atype_lj_parameter(i_mole,i_mole,3),atype_lj_parameter(i_mole,i_mole,4),atype_lj_parameter(i_mole,i_mole,5),atype_lj_parameter(i_mole,i_mole,6),atype_pol(i_mole)

       Case(11)
          ! read  atomname,  charge, Aexch, Aelec, Aind, Adhf, Adisp, C6, C8, C10, alpha
          read(file_h,*) atype_name(i_mole),atype_chg(i_mole),atype_bkghm_decomp(i_mole,i_mole,1),atype_bkghm_decomp(i_mole,i_mole,2),atype_bkghm_decomp(i_mole,i_mole,3),atype_bkghm_decomp(i_mole,i_mole,4),atype_bkghm_decomp(i_mole,i_mole,5),atype_lj_parameter(i_mole,i_mole,3),atype_lj_parameter(i_mole,i_mole,4),atype_lj_parameter(i_mole,i_mole,5),atype_pol(i_mole)        
          ! decide whether we are reading in C12 coeff

       Case(12)
          ! here we have to decide whether C12 dispersion is being read
          Select Case(C12_dispersion)
          Case("no")
             read(file_h,*) atype_name(i_mole),atype_chg(i_mole),atype_bkghm_decomp(i_mole,i_mole,1),atype_bkghm_decomp(i_mole,i_mole,2),atype_bkghm_decomp(i_mole,i_mole,3),atype_bkghm_decomp(i_mole,i_mole,4),atype_bkghm_decomp(i_mole,i_mole,5),atype_lj_parameter(i_mole,i_mole,2),atype_lj_parameter(i_mole,i_mole,3),atype_lj_parameter(i_mole,i_mole,4),atype_lj_parameter(i_mole,i_mole,5),atype_pol(i_mole)
          Case("yes")
             read(file_h,*) atype_name(i_mole),atype_chg(i_mole),atype_bkghm_decomp(i_mole,i_mole,1),atype_bkghm_decomp(i_mole,i_mole,2),atype_bkghm_decomp(i_mole,i_mole,3),atype_bkghm_decomp(i_mole,i_mole,4),atype_bkghm_decomp(i_mole,i_mole,5),atype_lj_parameter(i_mole,i_mole,3),atype_lj_parameter(i_mole,i_mole,4),atype_lj_parameter(i_mole,i_mole,5),atype_lj_parameter(i_mole,i_mole,6),atype_pol(i_mole)
          End Select

       case default
          stop "Not a valid number of arguments (parameters) for solute atoms under solute-species section!"
       End Select

    enddo


  end subroutine read_solute_parameters_one_set



  !*********************************************
  ! this subroutine reads in force field parameters for a sapt-type force field, if there are two sets of
  ! solute parameters being used, (i.e. we ARE using a different set of force field parameters for solute-framework interactions)
  !
  ! note if lj_bkghm = 3, we are using a hybrid lj/bkghm force field, with the lj ff for solute-solute interactions, and the bkghm force field for
  ! solute framework interactions
  !*********************************************
  subroutine read_solute_parameters_two_set(file_h,n_type_s,C12_dispersion,explicit_cross_term_exponents,atype_name,atype_chg,atype_bkghm_decomp,atype_lj_parameter,temp_lj_parameter, atype_pol,lj_bkghm)
    integer,intent(in)::file_h,n_type_s,lj_bkghm
    character(*),intent(in) :: C12_dispersion,explicit_cross_term_exponents
    character(*),dimension(:),intent(inout) :: atype_name
    real*8,dimension(:),intent(inout) :: atype_chg,atype_pol
    real*8,dimension(:,:,:),intent(inout) :: atype_bkghm_decomp
    real*8,dimension(:,:,:),intent(inout) :: temp_lj_parameter
    real*8,dimension(:,:,:),intent(inout) :: atype_lj_parameter

    integer:: i_mole, inputstatus, ind
    character(50)::line


    !********************************** here we are reading in two sets of solute parameters, one for solute-solute interaction, the other for solute framework interactions. charges and drude oscillators are same for both


    Select Case(lj_bkghm)
    Case(1)
       ! sapt-type bkghm force field for solute-solute interactions
       ! read in order for each atom, chg, A_exch, A_elec, A_induc,A_dhf,A_disp,C6,C8,C10, pol , put bkghm, C6,C8,C10 in temp_lj_parameter array since these are solute-solute terms
       do i_mole=1, n_type_s
          read(file_h,*) atype_name(i_mole),atype_chg(i_mole),temp_lj_parameter(i_mole,i_mole,1),temp_lj_parameter(i_mole,i_mole,2),temp_lj_parameter(i_mole,i_mole,3),temp_lj_parameter(i_mole,i_mole,4),temp_lj_parameter(i_mole,i_mole,5),temp_lj_parameter(i_mole,i_mole,6),temp_lj_parameter(i_mole,i_mole,7),temp_lj_parameter(i_mole,i_mole,8),atype_pol(i_mole)
       enddo

    Case(3)
       ! lennard-jones force field for solute-solute interactions
       ! read in order for each atom, chg, epsilon, sigma, junk, pol , put in temp_lj_parameter array since these are solute-solute terms
       do i_mole=1, n_type_s
          read(file_h,*) atype_name(i_mole),atype_chg(i_mole),temp_lj_parameter(i_mole,i_mole,1),temp_lj_parameter(i_mole,i_mole,2),temp_lj_parameter(i_mole,i_mole,3),atype_pol(i_mole)
       enddo
    Case default
       stop "unrecognized setting of variable lj_bkghm in read_solute_parameters_two_set subroutine"
    End Select

    ! warn here if C12_disperion="yes", that we do not read in C12 terms for solute-solute interactions
    Select Case(C12_dispersion)
    Case("yes")
       write(*,*) "Although C12 dispersion terms are being used to model solute-framework interactions, C12 terms will not be read in or used for solute-solute interactions"
    End Select

    do
       Read(file_h,'(A)',Iostat=inputstatus) line
       If(inputstatus < 0) Exit
       ind=INDEX(line,'framework cross terms')
       IF(ind .NE. 0) Exit
    enddo

    ! if we are not using explicit cross terms for exponents, read them in here
    Select Case(explicit_cross_term_exponents)
    Case("yes")
       do i_mole=1, n_type_s
          read(file_h,*) atype_name(i_mole),atype_bkghm_decomp(i_mole,i_mole,1),atype_bkghm_decomp(i_mole,i_mole,2),atype_bkghm_decomp(i_mole,i_mole,3),atype_bkghm_decomp(i_mole,i_mole,4),atype_bkghm_decomp(i_mole,i_mole,5),atype_lj_parameter(i_mole,i_mole,3),atype_lj_parameter(i_mole,i_mole,4),atype_lj_parameter(i_mole,i_mole,5)
       enddo
    Case("no")
       ! decide if C12 dispersion
       Select Case(C12_dispersion)
       Case("no")
          do i_mole=1, n_type_s
             read(file_h,*) atype_name(i_mole),atype_bkghm_decomp(i_mole,i_mole,1),atype_bkghm_decomp(i_mole,i_mole,2),atype_bkghm_decomp(i_mole,i_mole,3),atype_bkghm_decomp(i_mole,i_mole,4),atype_bkghm_decomp(i_mole,i_mole,5),atype_lj_parameter(i_mole,i_mole,2),atype_lj_parameter(i_mole,i_mole,3),atype_lj_parameter(i_mole,i_mole,4),atype_lj_parameter(i_mole,i_mole,5)
          enddo
       Case("yes")
          do i_mole=1, n_type_s
             read(file_h,*) atype_name(i_mole),atype_bkghm_decomp(i_mole,i_mole,1),atype_bkghm_decomp(i_mole,i_mole,2),atype_bkghm_decomp(i_mole,i_mole,3),atype_bkghm_decomp(i_mole,i_mole,4),atype_bkghm_decomp(i_mole,i_mole,5),atype_lj_parameter(i_mole,i_mole,2),atype_lj_parameter(i_mole,i_mole,3),atype_lj_parameter(i_mole,i_mole,4),atype_lj_parameter(i_mole,i_mole,5),atype_lj_parameter(i_mole,i_mole,6)
          enddo
       End Select
    End Select


  end subroutine read_solute_parameters_two_set





  !**************************************************
  ! This subroutine reads in explicit solute-solute cross terms
  ! for either exponents and/or dhf penetration terms
  !**************************************************
  subroutine read_solute_explicit_cross_exponents_and_dhf(file_h,n_type_s,explicit_cross_term_exponents,atype_bkghm_decomp,atype_lj_parameter,temp_lj_parameter, framework_simulation, lj_bkghm, solute_cross_parameter_set, ind, ind1, ind_dhf, ind_frame)
    integer,intent(in)::file_h, n_type_s, framework_simulation, lj_bkghm
    integer,intent(inout) :: ind, ind1, ind_dhf, ind_frame
    character(*),intent(in) :: explicit_cross_term_exponents, solute_cross_parameter_set
    real*8,dimension(:,:,:),intent(inout) :: atype_bkghm_decomp
    real*8,dimension(:,:,:),intent(inout) :: temp_lj_parameter
    real*8,dimension(:,:,:),intent(inout) :: atype_lj_parameter

    integer :: i_mole, j_mole, inputstatus, n_cross
    real*8,dimension(:),allocatable::temp_cross
    character(50)::line

    ! if ind1 = 1, means we have skipped to solute-solute-exponents
    ! if ind_frame =1, means we have skipped to framework_species section
    if ( ( ind1 .eq. 0 ) .and. ( ind_frame .eq. 0 ) ) then
       ! read dhf explicit cross terms if present
       do
          Read(file_h,'(A)',Iostat=inputstatus) line
          If(inputstatus < 0) Exit
          ind_dhf=INDEX(line,'solute dhf cross terms')
          IF(ind_dhf .NE. 0) Exit
          ind=INDEX(line,'solute-solute')
          IF(ind .NE. 0) Exit
          ind_frame=INDEX(line,'framework_species')
          IF(ind_frame .NE. 0) Exit
       enddo
    else
       ind_dhf = 0
    endif



    if ( ind_dhf .ne. 0 ) then
       !***************************** here we are reading in explict dhf cross terms, decide which arrays to put them in based on number of sets of solute parameters
       Select Case ( solute_cross_parameter_set )
       Case("no")
          ! one set of solute parameters, no need to store in temporary array
          do i_mole=1,n_type_s-1
             read(file_h,*) (atype_bkghm_decomp(i_mole,j_mole,4),j_mole=i_mole+1,n_type_s)
          enddo
          ! now transpose
          do i_mole=1,n_type_s-1
             do j_mole = i_mole+1,n_type_s
                atype_bkghm_decomp(j_mole,i_mole,4) = atype_bkghm_decomp(i_mole,j_mole,4)
             enddo
          enddo
       Case("yes")
          ! two sets of solute parameters, store in temporary array
          do i_mole=1,n_type_s-1
             read(file_h,*) (temp_lj_parameter(i_mole,j_mole,4),j_mole=i_mole+1,n_type_s)
          enddo
          ! now transpose
          do i_mole=1,n_type_s-1
             do j_mole = i_mole+1,n_type_s
                temp_lj_parameter(j_mole,i_mole,4) = temp_lj_parameter(i_mole,j_mole,4)
             enddo
          enddo

       End Select
       ! now go to next section
       do
          Read(file_h,'(A)',Iostat=inputstatus) line
          If(inputstatus < 0) Exit
          ind=INDEX(line,'solute-solute')
          IF(ind .NE. 0) Exit
          ind_frame=INDEX(line,'framework_species')
          IF(ind_frame .NE. 0) Exit
       enddo

    endif

    if ( ind .ne. 0 ) then
       ! get solute-solute exponents, these are on same record, so that's why we're reading them in this way
       ! if using explicit_cross_term_exponents="no", we will have read solute-framework exponents into 
       ! atype_lj_parameter(:,:,2) array, so until we generate cross-term exponents, store 
       ! solute-solute exponents in temp_lj_parameter(:,:,9) array slots
       n_cross = (n_type_s**2 - n_type_s)/2 + n_type_s
       allocate( temp_cross( n_cross ) )
       read(file_h, *) (temp_cross(ind) , ind=1,n_cross)

       Select Case(explicit_cross_term_exponents)
       Case("yes")
          ind =1
          do i_mole=1,n_type_s
             do j_mole=i_mole, n_type_s
                atype_lj_parameter(i_mole,j_mole,2) = temp_cross(ind)
                ind = ind + 1
             enddo
          enddo
       Case("no")
          ind =1
          do i_mole=1,n_type_s
             do j_mole=i_mole, n_type_s
                temp_lj_parameter(i_mole,j_mole,9) = temp_cross(ind)
                ind = ind + 1
             enddo
          enddo
       End Select

       deallocate( temp_cross )

    else

       ! these checks only apply if lj_bkghm != 3
       Select Case(lj_bkghm)
       Case(3)
          continue
       case default
          ! didn't find solute-solute explicit exponent section.  Make sure explicit exponents isn't set to yes
          Select Case(explicit_cross_term_exponents)
          Case("yes")
             stop " explicit_cross_term_exponents was set to 'yes', but a solute-solute exponents section was not found in the parameter file"
          End Select
          ! if solute_cross_parameter_set was set to 'yes', it is assumed that explicit cross term exponents would be read for solute-solute interactions
          Select Case(solute_cross_parameter_set)
          Case("yes")
             stop "Setting solute_cross_parameter_set='yes' assumes that explicit solute-solute exponents will be read in.  However, this section has not been found"
          End Select

          ! if we didn't find a solute-solute exponents section, make sure to form_solute_exponents,
          ! also make sure we skipped to framework_species section if this is a framework simulation
          ! set temp_lj_param(1,1,9) < -100 to signal that we need to form cross solute-solute exponents
          temp_lj_parameter(1,1,9) = -101d0

       End Select

       Select Case(framework_simulation)
       Case(1)
          if ( ind_frame .eq. 0 ) then
             stop "error reading parameter input file.  Solute-solute exponents section wasn't found, and neither was framework species section"
          endif
       End Select

    endif


  end subroutine read_solute_explicit_cross_exponents_and_dhf






  !***************************************************
  !  This subroutine reads C9 coefficients for three body dispersion from
  !  parameter file
  !
  !  for N atom types, the C9 coefficients should be input as follows:
  !  There should be N + (N-1) + (N-2) ... lines (records) containing
  !  C9 parameters.  These parameters are read in as 
  !  loop (i=1, N)
  !    loop ( j=i,N)
  !      loop ( k=j,N)
  !
  !  where the parameters of the last loop are on one record,
  !  and the first two loops account for the number of lines
  !
  !  these parameters should be given in units of (kJ/mol)*A^9, consistent with 
  !  the rest of the parameters
  !***************************************************
  subroutine read_C9_three_body_dispersion_parameters(ifile_pmt)
    use global_variables
    character(*),intent(in) :: ifile_pmt

    character(50)::line
    integer :: i, j, inputstatus,ind, file_h=15, itype,jtype,ktype, ordered_index(3),store(3)

    open(unit=file_h,file=ifile_pmt,status="old")

    do
       Read(file_h,'(A)',Iostat=inputstatus) line
       If(inputstatus < 0) Exit
       ind=INDEX(line,'three body dispersion')
       IF(ind .NE. 0) Exit
    enddo

    ! if we didn't find the three body dispersion section in parameter file, then stop!
    if ( ind .eq. 0 ) then
       stop "couldn't find three body dispersion section in parameter input file"
    endif

    do itype = 1 , n_atom_type
       do jtype = itype, n_atom_type
          read(file_h,*) ( atype_3body_C9(itype,jtype,ktype),ktype=jtype, n_atom_type )
       enddo
    enddo

    ! now fill in the rest of the array

    do itype = 1 , n_atom_type
       do jtype = 1 , n_atom_type
          do ktype =1 , n_atom_type

             store(1)=itype; store(2)=jtype; store(3)=ktype
             ! need to arrange these indices in increasing order.  Probably better ways to do this but...
             do i=1,3
                ordered_index(i) = min( store(1), store(2), store(3) )

                ! now give this index a higher value so we can keep using min() function
                do j=1,3
                   if ( store(j) .eq. ordered_index(i) ) then
                      store(j) = n_atom_type + 1
                      exit
                   endif
                enddo
             enddo

             ! now fill in array with C9
             atype_3body_C9(itype,jtype,ktype) = atype_3body_C9(ordered_index(1), ordered_index(2), ordered_index(3) )

          enddo
       enddo
    enddo


  end subroutine read_C9_three_body_dispersion_parameters




 !***************************************************
  !  This subroutine reads three body exchange parameters
  !
  !  first coefficients are given, then under these, exponents should be given
  !
  !  for N atom types, the parameters should be input as follows:
  !  There should be N + (N-1) + (N-2) ... lines (records) containing
  !  C9 parameters.  These parameters are read in as 
  !  loop (i=1, N)
  !    loop ( j=i,N)
  !      loop ( k=j,N)
  !
  !  where the parameters of the last loop are on one record,
  !  and the first two loops account for the number of lines
  !
  !  these parameters should be given in units of kJ/mol and Angstrom, consistent with 
  !  the rest of the parameters
  !***************************************************
  subroutine read_three_body_exchange_parameters(ifile_pmt)
    use global_variables
    character(*),intent(in) :: ifile_pmt

    character(50)::line
    integer :: i, j, inputstatus,ind, file_h=15, itype,jtype,ktype, ordered_index(3),store(3)

    open(unit=file_h,file=ifile_pmt,status="old")

    do
       Read(file_h,'(A)',Iostat=inputstatus) line
       If(inputstatus < 0) Exit
       ind=INDEX(line,'three body exchange')
       IF(ind .NE. 0) Exit
    enddo

    ! if we didn't find the three body exchange section in parameter file, then stop!
    if ( ind .eq. 0 ) then
       stop "couldn't find three body exchange section in parameter input file"
    endif

    do itype = 1 , n_atom_type
       do jtype = itype, n_atom_type
          read(file_h,*) ( atype_3body_exchange(itype,jtype,ktype,1),ktype=jtype, n_atom_type )
       enddo
    enddo
    do itype = 1 , n_atom_type
       do jtype = itype, n_atom_type
          read(file_h,*) ( atype_3body_exchange(itype,jtype,ktype,2),ktype=jtype, n_atom_type )
       enddo
    enddo

    ! now fill in the rest of the array

    do itype = 1 , n_atom_type
       do jtype = 1 , n_atom_type
          do ktype =1 , n_atom_type

             store(1)=itype; store(2)=jtype; store(3)=ktype
             ! need to arrange these indices in increasing order.  Probably better ways to do this but...
             do i=1,3
                ordered_index(i) = min( store(1), store(2), store(3) )

                ! now give this index a higher value so we can keep using min() function
                do j=1,3
                   if ( store(j) .eq. ordered_index(i) ) then
                      store(j) = n_atom_type + 1
                      exit
                   endif
                enddo
             enddo

             ! now fill in array with C9
             atype_3body_exchange(itype,jtype,ktype,:) = atype_3body_exchange(ordered_index(1), ordered_index(2), ordered_index(3),: )

          enddo
       enddo
    enddo


  end subroutine read_three_body_exchange_parameters



end module sapt_ff_routines
