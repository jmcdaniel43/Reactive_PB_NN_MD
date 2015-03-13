module egc_routines
use routines

contains

  !*****************************************************
  ! this subroutine controls the adjustment of the partial
  ! molecule parameters for an expanded_grand_canonical simulation
  ! 
  ! the incorporation of a partial molecule into the existing codes for
  ! energy calculations, etc, is done by creating additional atomtypes
  ! for this partial molecule.  For electrostatics, global static and drude
  ! charges can be updated.  For lennard jones, parameters are taken from global
  ! atom type arrays.  Therefore, new atom types will be added, and the global array 
  ! elements for the new atom types can be adjusted with the amount of coupling
  !
  ! notice we will not need to scale atype_bkghm_decomp array as this is just
  ! for energy decomposition and sampling is not implemented for intermediate scale factors
  !***************************************************** 
  subroutine adjust_partial_molecule_parameters(increment,n_atom,chg)
    use global_variables
    integer,intent(in) :: increment
    integer,dimension(:),intent(in)::n_atom
    real*8,dimension(:,:),intent(inout)::chg

    integer :: count,i_atom,j_atom,index,flag,i
    integer,dimension(MAX_N_ATOM):: temp_index
    real*8 :: scale_old, scale_new
    integer,save :: warn_energy_decomp=0

    ! we are not scaling atype_bkghm_parameter arrays, and therefore these arrays will only be correct at terminal coupling values.
    ! thus the energy decomposition will be incorrect throughout simulation, and in order to get correct decomposition full LJ routine
    ! must be called at sampling time when a terminal coupling value is realized
    Select Case(energy_decomposition)
    Case("yes")
       if (warn_energy_decomp .eq. 0 ) then
          write(*,*) "atype_bkghm_dcomp arrays are not scaled, and therefore full lj subroutine should be called after molecule is either fully inserted or deleted to get energy decomposition correct"
          warn_energy_decomp = 1
       endif
    End Select


    ! first we need to decide whether we are creating new atom types, or just adjusting the parameters on these new atom types

    Select Case(partial_molecule_coupling)
    Case(0)
       ! we need to create new atom types for insertion or deletion, so either value of increment is ok

       ! figure out how many new atom types to create (one for each unique atom in partial molecule)
       ! need at least one
       count=1
       index = atom_index(partial_molecule_index,1)
       temp_index(1)= index
       atype_lj_parameter(n_atom_type+1,:,:) = atype_lj_parameter(index,:,:)
       atype_lj_parameter(:,n_atom_type+1,:) = atype_lj_parameter(:,index,:)
       ! now update atom index for this first atom
       atom_index(partial_molecule_index,1) = n_atom_type+1

       if( n_atom(partial_molecule_index) > 1 ) then
          do i_atom = 2, n_atom(partial_molecule_index)
             index = atom_index(partial_molecule_index,i_atom)
             ! check if we need a new atom type
             flag=1
             do j_atom =1 , i_atom-1
                if( index == temp_index(j_atom)) then
                   ! we already made a new type for this
                   atom_index(partial_molecule_index,i_atom) = atom_index(partial_molecule_index,j_atom)
                   flag=0
                   exit
                endif
             enddo
             if ( flag .eq. 1 ) then
                ! need a new type
                count=count+1
                atype_lj_parameter(n_atom_type+count,:,:) = atype_lj_parameter(index,:,:)
                atype_lj_parameter(:,n_atom_type+count,:) = atype_lj_parameter(:,index,:)  
                ! now update atom index for this atom
                atom_index(partial_molecule_index,i_atom) = n_atom_type+count
                temp_index(i_atom) = index
             endif
          enddo
       endif

       ! now change n_atom_type to its new value
       n_atom_type = n_atom_type_store + count

    end Select


    ! rescale parameters, figure out ratio for this

    Select Case(partial_molecule_coupling)
    Case(0)
       Select Case(increment)
       Case(0)
          ! we are decoupling a completely inserted molecule
          scale_old = 1d0
          scale_new = coupling_increment * dble(partial_coupling_steps - 1 )
       Case(1)
          ! we are inserting a lowest coupled molecule
          scale_old = 1d0
          scale_new = coupling_increment * dble(partial_molecule_coupling + 1 )
       End Select
       ! check to make sure deletion isn't treated here
    Case(1)
       select case(increment)
       Case(0)
          stop "check conditions under which adjust_partial_molecule_parameters is being called"
       Case(1)
          ! normal increment
          scale_old = coupling_increment * dble(partial_molecule_coupling)
          scale_new = coupling_increment * dble(partial_molecule_coupling + 1)
       end select
    case default
       scale_old = coupling_increment * dble(partial_molecule_coupling)
       Select Case(increment)
       Case(1) 
          scale_new = coupling_increment * dble(partial_molecule_coupling + 1)
       Case(0)
          scale_new = coupling_increment * dble(partial_molecule_coupling - 1)
       End Select
    End Select

    ! now scale , first do atype arrays, we are only messing with array components from n_atom_type_store +1 through n_atom_type

    do i_atom = n_atom_type_store +1, n_atom_type
       ! be a bit careful here, for buckingham we don't want to scale exponents
       ! for lj, just scale epsilon, since this is simplest
       Select Case(lj_bkghm)
       Case(1)
          do i=1,size(atype_lj_parameter(1,1,:))
             if ( i .ne. 2 ) then
                ! this is buckingham,scale everything except the exponents which are stored in 2nd index
                atype_lj_parameter(i_atom,:,i) =   (scale_new/scale_old) * atype_lj_parameter(i_atom,:,i)
                ! it's ok if diagonal elements from n_atom_type + 1 on get scaled more than once, as these will never be used
                atype_lj_parameter(:,i_atom,i) =   (scale_new/scale_old) * atype_lj_parameter(:,i_atom,i)      
             end if
          enddo
       Case(2)
          ! this is lj, scale epsilon
          atype_lj_parameter(i_atom,:,1) =   (scale_new/scale_old) * atype_lj_parameter(i_atom,:,1)
          atype_lj_parameter(:,i_atom,1) =   (scale_new/scale_old) * atype_lj_parameter(:,i_atom,1)    
       end Select
    enddo

    ! now scale charges
    chg(partial_molecule_index,:) = (scale_new/scale_old) * chg(partial_molecule_index,:)


  end subroutine adjust_partial_molecule_parameters




!*************************************************************************************
! this subroutine updates global data after the acceptance of a coupling move in the 
! expanded grand canonical ensemble.
! 
! specifically, this subroutine needs to update chg array and relink the partial molecule
! with the static atype parameters if the molecule has just been fully coupled into the system
!**********************************************************************************8
  subroutine update_partial_molecule(increment,n_atom,n_atom_drude,chg,alist)
    use global_variables
    integer,intent(in) :: increment
    integer,dimension(:),intent(in)::n_atom,n_atom_drude
    real*8,dimension(:,:),intent(inout)::chg
    character(*), dimension(:,:),intent(in) :: alist

    integer :: i_atom,i_drude,i_type,i_param,ind
    real*8,parameter::small=1D-6

    ! consider special case where molecule has just been fully coupled into the system
    ! in this case, we can get rid of all extraneous atype parameters and relink the molecule
    ! with the static parameters
    Select Case(partial_molecule_coupling)
    Case(partial_coupling_steps-1)
       Select Case(increment)
       Case(1)
          ! zero the atype_lj_parameters that were changed
          do i_atom = n_atom_type_store +1, n_atom_type  
             atype_lj_parameter(i_atom,:,:)=0d0
             atype_lj_parameter(:,i_atom,:)=0d0
          enddo
          ! reset n_atom_type
          n_atom_type = n_atom_type_store
          ! set partial_molecule_coupling to zero, since we have completely inserted a molecule
          partial_molecule_coupling=0

          ! relink molecule with old parameters
          do i_atom=1,n_atom(partial_molecule_index)
             do  i_type=1, n_atom_type
                if ( atype_name(i_type) .eq. alist(partial_molecule_index,i_atom ) ) then
                   i_param = i_type
                   exit 
                endif
             enddo
             ! set index, chg, and polarizability
             atom_index(partial_molecule_index,i_atom) = i_param

             chg(partial_molecule_index,i_atom) = atype_chg(i_param)

          enddo

          ! if there are drude oscillators on this molecule, fix charges
          if ( n_atom_drude(partial_molecule_index) > n_atom(partial_molecule_index) ) then
             do i_drude = n_atom(partial_molecule_index)+1, n_atom_drude(partial_molecule_index)
                ! use map to figure out which atom this oscillator is on
                i_atom = drude_atom_map(partial_molecule_index,i_drude)
                i_param = atom_index(partial_molecule_index,i_atom)

                chg(partial_molecule_index,i_atom) = atype_chg(i_param)+sqrt(atype_pol(i_param)*springcon)
                chg(partial_molecule_index,i_drude) = -sqrt(atype_pol(i_param)*springcon)
             enddo
         endif

                
       Case(0)
          ! this is a move from partial_coupling = M -1 to M -2 where M is partial_coupling_steps
          ! only need to change partial_molecule_coupling here
          partial_molecule_coupling = partial_molecule_coupling - 1
       End Select
    Case(0)
       Select Case(increment)
          Case(0)
       ! here we have successfully decoupled a fully inserted molecule by one step, so adjust partial_molecule_coupling accordingly
         partial_molecule_coupling = partial_coupling_steps-1
          Case(1)
       ! here we have successfully inserted a lowest coupled molecule
         partial_molecule_coupling = partial_molecule_coupling + 1          
       End Select
     Case(1)
        Select Case(increment)
           Case(0)
        ! here we have completely removed a lowest coupled molecule.  Clean up global arrays
         ! zero the atype_lj_parameters that were changed
          do i_atom = n_atom_type_store +1, n_atom_type  
             atype_lj_parameter(i_atom,:,:)=0d0
             atype_lj_parameter(:,i_atom,:)=0d0
          enddo
          ! reset n_atom_type
          n_atom_type = n_atom_type_store
          ! set partial_molecule_coupling to zero, since we have completely removed a molecule
          partial_molecule_coupling=0   
          Case(1)
             ! here we are just increasing the coupling
          partial_molecule_coupling = partial_molecule_coupling + 1
          End Select
    Case default
       ! no more special cases
       Select Case(increment)
       Case(0)
          partial_molecule_coupling = partial_molecule_coupling - 1    
       Case(1)
          partial_molecule_coupling = partial_molecule_coupling + 1
       End Select

    End Select



  end subroutine update_partial_molecule
   

!***************************************************************
! this subroutine returns the exponential factor in the acceptance
! criteria for the expanded grand canonical ensemble
!
! this is a little confusing because the stage at which we formally add/remove
! molecules is different than when we practically remove these.
! we practically add/remove molecules when we call the grand_canonical_ins_rm
! subroutine, since this is when data arrays are changed.
!
! however, we only "count" the molecule as being part of the system when we
! make a transition from M-1 to 0 , and we "count" removing a molecule when we transition
! from 0 to M-1, even though these transitions don't physically remove any molecules from our system
!
! the definition of chemical potential is the same as in grand_canonical_ins_rm subroutine
! basically, we do not include rotational ideal gas contributions
!
! be careful at zero molecules ( here we modify the molecule number to always be one
! so that the logarithm doesn't blow up.  The way it is implemented should satisfy detailed balance
!**************************************************************
  subroutine get_egc_fac(egc_fac,increment,box,temp_molecule_index,n_mole)
    use global_variables
    real*8, intent(out) :: egc_fac
    integer,intent(in) :: increment
    real*8,dimension(:,:),intent(in)::box    
    integer,dimension(:),intent(in)::temp_molecule_index
    integer,intent(in)::n_mole

    real*8 :: vol,lambda,ideal_cont,phi_i, phi_f,beta,n_mole_use,mass_local

    ! check ensemble
    ! make sure that we never use zero molecules as an argument

    Select Case(select_ensemble)
    Case("uvt")
       egc_fac = 1d0
    Case("egc")
       vol = volume ( box(1,:),box(2,:),box(3,:) )
       mass_local = molecule_mass(temp_molecule_index(partial_molecule_index))

       lambda= 126.178D0/sqrt(2D0*pi*mass_local*8.314D0*temp)
       beta = 1000./(temp*8.314)
       ! get weight for step, consider special cases at end points
       Select Case(partial_molecule_coupling)
       Case(0)
          Select Case(increment)
          Case(0)
             ! decoupling a fully inserted molecule, here we are formally removing a molecule, notice that when we decouple, we lose a molecule
             n_mole_use = max(1,n_mole-1)
             phi_i = 1d0 * (beta*chem_potential + log(vol /(lambda**3 * dble(n_mole))) )
             phi_f = step_weight(partial_coupling_steps-1) * (beta*chem_potential + log(vol /(lambda**3 * dble(n_mole_use))) )

          Case(1)
             ! inserting a lowest coupled new molecule, here n_mole is the old number of molecules, passed from grand_canonical_ins_rm subroutine, since the new molecule doesn't count yet
             n_mole_use = max(1,n_mole)
             phi_i = 0d0
             phi_f = step_weight(partial_molecule_coupling+1) * (beta*chem_potential + log(vol /(lambda**3 * dble(n_mole_use))) )
          End Select
       Case(1)
          Select Case(increment)
          Case(0)
             ! removing a lowest coupled molecule, notice here n_mole doesn't count this lowest coupled molecule as passed by grand_canonical_ins_rm
             n_mole_use = max(1,n_mole)
             phi_i = step_weight(partial_molecule_coupling) * (beta*chem_potential + log(vol /(lambda**3 * dble(n_mole_use))) )
             phi_f = 0d0
          Case(1)
             ! increasing the coupling, note this coupled molecule doesn't count in total number
             n_mole_use = max(1,n_mole-1)
             phi_i = step_weight(partial_molecule_coupling)* (beta*chem_potential + log(vol /(lambda**3 * dble(n_mole_use))) )
             phi_f = step_weight(partial_molecule_coupling+1)* (beta*chem_potential + log(vol /(lambda**3 * dble(n_mole_use))) )
          End Select
       Case(partial_coupling_steps-1)
          Select Case(increment)
          Case(0)
             ! decreasing the coupling
             n_mole_use = max(1,n_mole-1)
             phi_i = step_weight(partial_molecule_coupling)* (beta*chem_potential + log(vol /(lambda**3 * dble(n_mole_use))) )
             phi_f = step_weight(partial_molecule_coupling-1)* (beta*chem_potential + log(vol /(lambda**3 * dble(n_mole_use))) )
          Case(1)
             ! fully coupling a molecule, this is when the molecule technically enters the system
             n_mole_use = max(1,n_mole-1)
             phi_i = step_weight(partial_coupling_steps-1)* (beta*chem_potential + log(vol /(lambda**3 * dble(n_mole_use))) )
             phi_f = 1d0 * (beta*chem_potential + log(vol /(lambda**3 * dble(n_mole))) )
          End Select
       Case default
          ! done with special cases
          n_mole_use = max(1,n_mole-1)
          Select Case(increment)
          Case(0)
             phi_i = step_weight(partial_molecule_coupling)* (beta*chem_potential + log(vol /(lambda**3 * dble(n_mole_use))) )
             phi_f = step_weight(partial_molecule_coupling-1)* (beta*chem_potential + log(vol /(lambda**3 * dble(n_mole_use))) )
          Case(1)
             phi_i = step_weight(partial_molecule_coupling)* (beta*chem_potential + log(vol /(lambda**3 * dble(n_mole_use))) )
             phi_f = step_weight(partial_molecule_coupling+1)* (beta*chem_potential + log(vol /(lambda**3 * dble(n_mole_use))) )
          End Select
       End Select

       egc_fac= exp(phi_f - phi_i)

    End Select

  end subroutine get_egc_fac

 end module egc_routines
