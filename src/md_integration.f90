
module md_integration
  use routines
  implicit none

contains

  !*************************************************************************
  !  this subroutine controls what ensemble to run, and calls appropriate subroutines
  !  currently, only NVE molecular dynamics is implemented
  !
  !  data structure molecule_data will be changed if MS-EVB simulation, and we
  !  have proton hop, hence intent(inout)
  !
  !  verlet_list will be changed in energy routines if needs update
  !*************************************************************************

  subroutine mc_sample( system_data , molecule_data , atom_data, integrator_data, verlet_list_data, PME_data, file_io_data )
    use global_variables
    type( system_data_type ) , intent(inout)    :: system_data
    type( molecule_data_type ), dimension(:), intent(inout) :: molecule_data
    type( atom_data_type )  , intent(inout)   :: atom_data
    type( integrator_data_type ) , intent(in) :: integrator_data
    type(verlet_list_data_type) , intent(inout)  :: verlet_list_data
    type(PME_data_type)     , intent(inout)      :: PME_data
    type(file_io_data_type) , intent(in)         :: file_io_data


    Select Case (integrator_data%ensemble)
    Case("NVE")
       !*********************** this is nve md ******************************!
       call md_integrate_atomic( system_data , molecule_data , atom_data, integrator_data, verlet_list_data, PME_data, file_io_data )
    case default
       stop "Currently, only NVE ensemble is implemented"
    end select

  end subroutine mc_sample



  !*************************************************************************
  ! samples velocities for all atoms from Maxwell-Boltzmann
  ! units of velocities are Angstrom/ps 
  !
  ! we allow the possibility of freezing atoms during a simulation, which is
  ! flagged using the atype_freeze data structure.  These frozen atoms will
  ! have zero velocity and won't contribute to the temperature
  !*************************************************************************
  subroutine sample_atomic_velocities(n_mole, n_atom, temperature, molecule_data, atom_data )
    use global_variables
    integer, intent(in)          :: n_mole,n_atom
    real*8,   intent(in)         :: temperature
    type(molecule_data_type), dimension(:),intent(in) :: molecule_data
    type(atom_data_type) , intent(inout)   :: atom_data

    ! allocate temporary arrays here for convenience
    real*8,dimension(:), allocatable :: mass
    integer,dimension(:), allocatable :: atom_type_index
     
    integer :: i_atom, i_type, n_tot
    real*8,dimension(2)  :: vel
    real*8,parameter :: small=1D-3
    real*8           :: conv_fac, kB
    real*8 :: sum_KE, norm

    conv_fac = constants%conv_kJmol_ang2ps2gmol   ! converts kJ/mol to A^2/ps^2*g/mol
    kB       = constants%boltzmann                ! 0.008314 kJ/mol/K

    allocate(mass(n_atom), atom_type_index(n_atom) )
    mass            = atom_data%mass
    atom_type_index = atom_data%atom_type_index
    
    n_tot=0
    sum_KE=0d0

    ! zero velocity, in case we're freezing an atom
    atom_data%velocity=0d0

    !**************** first pull velocities from maxwell-boltzmann distribution
       do i_atom=1, n_atom
          ! make sure mass is non-zero
          if ( mass(i_atom) < small ) then
             write(*,*) "trying to assign velocity for atom ", i_atom
             write(*,*) "but mass is zero!"
             stop
          end if

          i_type = atom_type_index(i_atom)
          ! get velocity if atomtype isn't frozen
          if ( atype_freeze(i_type) /= 1 ) then
             ! gaussian random numbers come in 2
             call max_boltz(vel,mass(i_atom),temperature, kB)
             atom_data%velocity(1,i_atom)=vel(1)
             atom_data%velocity(2,i_atom)=vel(2)
             call max_boltz(vel,mass(i_atom),temperature, kB)
             atom_data%velocity(3,i_atom)=vel(1)
          end if

       enddo

    ! get rid of excess center of mass momentum of system
    call subtract_center_of_mass_momentum(n_mole, molecule_data, atom_data )

    ! now rescale velocities to desired temperature
       do i_atom=1, n_atom
          i_type = atom_type_index(i_atom)
          ! add KE if atom isn't frozen
          if ( atype_freeze(i_type) /= 1 ) then
             n_tot = n_tot + 1
             sum_KE = sum_KE + 0.5d0 * mass(i_atom) * dot_product(atom_data%velocity(:,i_atom),atom_data%velocity(:,i_atom)) / conv_fac
          end if
       enddo

    norm = 1.5d0 * kB * temperature * dble(n_tot) / sum_KE
    atom_data%velocity = atom_data%velocity * sqrt(norm)

    deallocate(mass, atom_type_index )

  end subroutine sample_atomic_velocities




  !*******************************************************
  ! this subroutine calculates the center of mass momentum of the system,
  ! and subtracts the total net per atom contribution from each atom's momentum,
  ! so that the net COM momentum is zero
  !*******************************************************
  subroutine subtract_center_of_mass_momentum(n_mole, molecule_data, atom_data )
    use global_variables
    integer,intent(in)::n_mole
    type(molecule_data_type) , dimension(:), intent(in) :: molecule_data
    type(atom_data_type) , intent(inout)  :: atom_data

    !******** this is a local data structure with pointers that will be set
    ! to subarrays of atom_data arrays for the specific atoms in the molecule
    type(single_molecule_data_type) :: single_molecule_data


    integer :: i_mole, i_atom, i_type, n_tot
    real*8,dimension(3) :: rho_system, rho_excess

    n_tot=0
    rho_system=0d0

    !**************** calculate total COM momentum
    do i_mole=1,n_mole

       ! set pointers for this data structure to target molecule
       call return_molecule_block( single_molecule_data , molecule_data(i_mole)%n_atom, molecule_data(i_mole)%atom_index, atom_velocity=atom_data%velocity, atom_mass=atom_data%mass, atom_type_index=atom_data%atom_type_index )

       do i_atom=1, molecule_data(i_mole)%n_atom
          i_type = single_molecule_data%atom_type_index(i_atom)
          ! add momentum if atomtype isn't frozen
          if ( atype_freeze(i_type) /= 1 ) then
             n_tot = n_tot + 1
             rho_system(:) = rho_system(:) + single_molecule_data%mass(i_atom) * single_molecule_data%velocity(:,i_atom)
          end if
       enddo
    enddo

    !*************  now subtract excess momentum from each atom, so that center of mass momentum is zero
    rho_excess(:) = rho_system(:) / dble(n_tot)

    do i_mole=1,n_mole

       ! set pointers for this data structure to target molecule
       call return_molecule_block( single_molecule_data , molecule_data(i_mole)%n_atom, molecule_data(i_mole)%atom_index, atom_velocity=atom_data%velocity, atom_mass=atom_data%mass, atom_type_index=atom_data%atom_type_index )

       do i_atom=1, molecule_data(i_mole)%n_atom
          i_type = single_molecule_data%atom_type_index(i_atom)     
          ! change velocity if atom isn't frozen
          if ( atype_freeze(i_type) /= 1 ) then
             single_molecule_data%velocity(:,i_atom) = single_molecule_data%velocity(:,i_atom) - rho_excess(:) / single_molecule_data%mass(i_atom)
          end if
       enddo
    enddo

  end subroutine subtract_center_of_mass_momentum









  !************************************************************************
  ! this is the MD engine for atomistic molecular simulations.  
  ! currently, this uses the velocity verlet algorithm to integrate Newton's equations
  !  velocity units are A/ps
  !
  ! NOTE we have already checked for non-zero atomic masses in the sample_atomic_velocities
  ! subroutine, and so here we don't worry about it
  !
  !  data structure molecule_data will be changed if MS-EVB simulation, and we
  !  have proton hop, hence intent(inout)
  !
  !  verlet_list will be changed in energy routines if needs update
  !
  !************************************************************************
  subroutine md_integrate_atomic( system_data , molecule_data , atom_data, integrator_data, verlet_list_data, PME_data, file_io_data )
    use global_variables
    use total_energy_forces
    use ms_evb
    type( system_data_type ), intent(inout)     :: system_data
    type( molecule_data_type ), dimension(:), intent(inout)   :: molecule_data
    type( atom_data_type ) , intent(inout)      :: atom_data
    type( integrator_data_type ), intent(in)    :: integrator_data
    type(verlet_list_data_type), intent(inout)  :: verlet_list_data
    type(PME_data_type)     , intent(inout)     :: PME_data
    type(file_io_data_type) , intent(in)        :: file_io_data

    integer :: i_atom, i_type, total_atoms
    real*8  :: dt , conv_fac

    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "md integration step started at", time
    endif
    !***********************************************************************!

    ! define local variables for convenience
    dt = integrator_data%delta_t
    conv_fac = constants%conv_kJmol_ang2ps2gmol  ! converts kJ/mol to A^2/ps^2*g/mol
    total_atoms = system_data%total_atoms

    !******************* Velocity Verlet Integrator

    ! first calculate velocities at delta_t / 2
       do i_atom = 1, total_atoms
          i_type = atom_data%atom_type_index(i_atom)

          ! update position, velocity if atom isn't frozen
          if ( atype_freeze(i_type) /= 1 ) then
             ! first calculate velocities at delta_t / 2
             ! here force is in kJ/mol*Ang^-1
             atom_data%velocity(:,i_atom) = atom_data%velocity(:,i_atom) + dt / 2d0 / atom_data%mass(i_atom) * atom_data%force(:,i_atom) * conv_fac

             ! now calculate new atomic coordinates at delta_t
             atom_data%xyz(:,i_atom) = atom_data%xyz(:,i_atom) + atom_data%velocity(:,i_atom) * dt
          end if
       end do

       ! after updating atomic coordinates, calculate new center of mass of molecules
       call update_r_com( system_data%n_mole, molecule_data, atom_data )
       ! translate molecules back into the box if they have left
       call shift_molecules_into_box( system_data%n_mole , molecule_data , atom_data , system_data%box,  system_data%xyz_to_box_transform )


    !**********************get total forces and energies*****************************!
    Select Case(ms_evb_simulation)
    Case("yes")
       call ms_evb_calculate_total_force_energy( system_data, molecule_data, atom_data, verlet_list_data, PME_data, file_io_data, integrator_data%n_output )
    Case("no")
       call calculate_total_force_energy( system_data, molecule_data, atom_data, verlet_list_data, PME_data )
    End Select
    !********************************************************************************!


     ! now final velocities
    do i_atom = 1, total_atoms
        i_type = atom_data%atom_type_index(i_atom)

          ! update position, velocity if atom isn't frozen
          if ( atype_freeze(i_type) /= 1 ) then
             ! here force is in kJ/mol*Ang^-1
             atom_data%velocity(:,i_atom) = atom_data%velocity(:,i_atom) + dt / 2d0 / atom_data%mass(i_atom) * atom_data%force(:,i_atom) * conv_fac

             ! make sure forces aren't crazy
             if ( ( abs( atom_data%force(1,i_atom) ) > 10d4 ) .or. ( abs( atom_data%force(2,i_atom) ) > 10d4 ) .or. ( abs( atom_data%force(3,i_atom) ) > 10d4 ) ) then
                write(*,*) "force on atom ", i_atom , " is too big ", atom_data%force(:,i_atom)
                stop
             end if

          end if
    end do


    ! finally remove center of mass momentum.  This should be numerical noise, so shouldn't effect energy conservation
    call subtract_center_of_mass_momentum(system_data%n_mole, molecule_data, atom_data )

    !****************************timing**************************************!
    if(debug .eq. 1) then
       call date_and_time(date,time)
       write(*,*) "md integration step finished at", time
    endif
    !***********************************************************************!

  end subroutine md_integrate_atomic





end module md_integration
