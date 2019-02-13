!************************************************************
!
! This module contains global variables that are used in many subroutines
! throughout the  program
!
! With a few exceptions (such as the PME, Ewald charge grids), the global data
! in this module doesn't change once initialized.
!
! trajectory data that is changing during the course of the simulations (such as xyz coordinates)
! is normally not stored in this module
!
! The initialization of this global data occurs in multiple ways
!
! 1) some of this data is hard-coded in.  Any parameters that are hard coded in that the user might want to change are
!    placed at the top of the module, as labeled
!
! 2) some of the technical simulation parameter data is read in from the simulation_parameters.pmt input file
!
! 3) The force field parameters are stored here as global data, and are read in from the force field input file
!
! 4) Finally, some of the mapping arrays are constructed logically by the code.
!
!*************************************************************


module global_variables
use MKL_DFTI
implicit none

!***************************** NOTE THIS SECTION ***********************************************

  ! these are number of atom/molecule types in force field, control size of
  ! force field parameter arrays
  integer, parameter ::  MAX_N_MOLE_TYPE=10, MAX_N_ATOM_TYPE=25

  integer, parameter :: MAX_FN=100, MAX_ANAME=5, MAX_MNAME=5

! these variables determine whether to grid expensive functions in memory.  Code is much faster when these are set to 'yes'
  character(3)  :: grid_Tang_Toennies ! grid damping functions

  character(3), parameter  :: ms_evb_simulation="yes"   ! ms_evb
  character(3), parameter  :: print_ms_evb_data = "yes"  ! if yes, this will print extra evb trajectory info

!***********************************************************************************************




  !***********************************************************
  !                         MS-EVB3 parameters
  !
  ! distance thresholds for finding acceptor molecules
  real*8, parameter :: evb_first_solvation_cutoff = 5d0
  real*8, parameter :: evb_reactive_pair_distance = 2.5d0
  integer, parameter :: evb_max_neighbors=10  ! this controls the dimensions of some of the evb data structures
  !
  ! evb_max_states is the maximum dimension for evb-hamiltonian, and therefore
  ! should be set to the maximum number of diabats.
  integer, parameter :: evb_max_states=80
  !
  ! maximum size of water chain for proton hopping.  Note this may need to be larger than the number
  ! of solvation shells for proton transfer.  For instance in a water hexamer, for a particular initial
  ! proton configuration, may need 6 hops to generate all diabats
  integer, parameter :: evb_max_chain=3           ! 3 is a good choice for bulk water
  !
  integer , parameter :: n_proton_max =10
  ! this keeps track of the hydronium molecules
  integer        :: n_hydronium_molecules
  integer, dimension( n_proton_max ) :: hydronium_molecule_index
  !
  integer, parameter :: max_interaction_type=15
  !
  ! this array lists the atom types for each ms-evb donor-acceptor interaction that
  ! is considered, i.e. equation 7 in JPCB,2008,112,467-482 (see erratum for correct formula)
  ! order of atomtypes in this array is 1) heavy-atom acceptor, 2) heavy-atom donor, 3) proton
  integer, dimension(max_interaction_type,3) :: evb_donor_acceptor_interaction
  ! this array contains donor_acceptor ms-evb interaction parameters
  real*8, dimension(max_interaction_type,6) :: evb_donor_acceptor_parameters
  !
  !
  ! this array lists the atomtypes for each ms-evb proton_acceptor interaction
  ! that is considered, i.e. equation 8 in JPCB,2008,112,467-482
  ! order of atomtypes in this array is 1) heavy-atom acceptor, 2) proton
  integer, dimension(max_interaction_type,2) :: evb_proton_acceptor_interaction
  ! this array contains donor_acceptor ms-evb interaction parameters
  real*8, dimension(max_interaction_type,5) :: evb_proton_acceptor_parameters
  !
  ! this array lists the atom types for each diabat coupling element that
  ! is considered, i.e. equation 10 and 12 in JPCB,2008,112,467-482
  ! order of atomtypes in this array is 1) heavy-atom acceptor, 2) heavy-atom donor, 3) proton
  integer, dimension(max_interaction_type,3) :: evb_diabat_coupling_interaction
  ! this array contains donor_acceptor ms-evb interaction parameters
  real*8, dimension(max_interaction_type,10) :: evb_diabat_coupling_parameters
  ! this array specifies which diabatic coupling function type we're using
  integer,dimension(max_interaction_type)    :: evb_diabat_coupling_type   ! currently, "1" specifies the  ms-evb3 functional form, "2" specifies product of gaussians in q and ROO
  !
  ! this array contains exchange charges for the atomtypes
  real*8, dimension(MAX_N_ATOM_TYPE)  :: evb_exchange_charge_atomic
  ! this array contains exchange charges for the transferring proton between two molecules
  real*8, dimension(MAX_N_MOLE_TYPE,MAX_N_MOLE_TYPE)  :: evb_exchange_charge_proton
  !
  ! this array stores molecules with acidic protons, molecule type is denoted "1" if it is an acid
  integer, dimension(MAX_N_MOLE_TYPE) :: evb_acid_molecule
  ! this array stores proton acceptor molecules, molecule type is denoted "1" if it is a base
  integer, dimension(MAX_N_MOLE_TYPE) :: evb_basic_molecule  
  ! this array stores index of conjugate acid/base moleculetype
  integer, dimension(MAX_N_MOLE_TYPE)  :: evb_conjugate_pairs
  ! 
  ! this array stores conjugate atoms of conjugate acid/base.  For each atomtype, index of
  ! conjugate atomtype is stored
  integer, dimension(MAX_N_ATOM_TYPE) :: evb_conjugate_atom_index
  !
  ! this array stores the relative chemical energy of each molecule, which
  ! will be incorporated into the adiabatic hamiltonian matrix elements
  real*8, dimension(MAX_N_MOLE_TYPE) :: evb_reference_energy
  ! this array stores index of proton atom type for each acid type
  integer, dimension(MAX_N_MOLE_TYPE) :: evb_proton_index
  ! this array stores index of acidic heavy atom type for  each acid type
  integer, dimension(MAX_N_MOLE_TYPE) :: evb_heavy_acid_index
  !***********************************************************


  !************************ defined type to store general system information
  type system_data_type
   integer :: n_mole  ! this is number of molecules in system
   integer :: total_atoms ! total number of atoms in simulation
   real*8, dimension(3,3) :: box 
   real*8, dimension(3,3) :: xyz_to_box_transform  ! transformation matrix to box vectors
   real*8                 :: volume
   real*8                 :: temperature      ! temperature for generating Maxwell-Boltzmann velocities (no thermostat yet)
   real*8                 :: kinetic_energy   ! total system KE
   real*8                 :: potential_energy ! total system PE
   real*8                 :: E_elec           ! electrostatic energy
   real*8                 :: E_vdw            ! all non-bond energy that is not electrostatic
   real*8                 :: E_bond           ! bond energy
   real*8                 :: E_angle          ! angle energy
   real*8                 :: E_dihedral       ! dihedral energy
  endtype system_data_type


  !***************  defined data type for integrator
  type integrator_data_type
   character(3)           :: ensemble  ! currently NVE molecular dynamics, or minimization
   integer                :: n_step    ! number of integration steps
   integer                :: n_output  ! frequency to print trajectory
   real*8                 :: delta_t   ! time step for integration (ps)
  end type integrator_data_type


  !************************* defined type to store all atom information
  type atom_data_type
   real*8, dimension(:,:), allocatable :: xyz
   real*8, dimension(:,:), allocatable :: velocity
   real*8, dimension(:,:), allocatable :: force
   real*8, dimension(:), allocatable   :: mass
   real*8, dimension(:), allocatable   :: charge
   integer, dimension(:), allocatable  :: atom_type_index  ! this is index of atom_type to look up force field parameters for this atom
   character(MAX_ANAME),dimension(:), allocatable :: aname
  end type atom_data_type

 !************************* defined type to store all molecule information
 type molecule_data_type
   integer :: n_atom            ! number of atoms in this molecule
   integer :: molecule_type_index    ! index of molecule type in molecule_type array
   character(MAX_MNAME) :: mname
   ! IMPORTANT : Note that we allocate atom_index to be 1 unit bigger than
   ! number of atoms of molecule, in case proton is transferred to this molecule
   integer, dimension(:),allocatable :: atom_index  ! this stores index of atom in atom_data array
   real*8, dimension(3) :: r_com  ! this stores center of mass of the molecule
 end type molecule_data_type



 !********************* this is meant to be used as a temporary local data structure
 ! as these pointers are meant to reference molecular subsections of the arrays in the
 ! atom_data_type structure, we don't want them floating around, so they should
 ! be used and destroyed locally
 !************************
 type single_molecule_data_type
   real*8, dimension(:,:), pointer :: xyz
   real*8, dimension(:,:), pointer :: velocity
   real*8, dimension(:,:), pointer :: force
   real*8, dimension(:), pointer   :: mass
   real*8, dimension(:), pointer   :: charge
   integer, dimension(:), pointer  :: atom_type_index  ! this is index of atom_type to look up force field parameters for this atom
   character(MAX_ANAME),dimension(:), pointer :: aname
 end type single_molecule_data_type


  
 !************************* defined type to store all i/o file information
 type file_io_data_type
  character(MAX_FN) :: ifile_gro
  character(MAX_FN) :: ifile_ffpmt
  character(MAX_FN) :: ifile_top
  character(MAX_FN) :: ifile_simpmt
  character(MAX_FN) :: ifile_velocity ! use this if restarting a simulation
  character(MAX_FN) :: ofile_traj
  character(MAX_FN) :: ofile_log
  character(MAX_FN) :: ofile_hop
  character(MAX_FN) :: ofile_hamiltonian
  integer           :: ifile_gro_file_h  ! these are the file handles
  integer           :: ifile_ffpmt_file_h 
  integer           :: ifile_top_file_h 
  integer           :: ifile_simpmt_file_h 
  integer           :: ifile_velocity_file_h
  integer           :: ofile_traj_file_h 
  integer           :: ofile_log_file_h 
  integer           :: ofile_hop_file_h 
  integer           :: ofile_hamiltonian_file_h
 end type file_io_data_type


 !********************* defined type to store constants
 ! the values of these constants will be initialized in the
 ! initialize_constants subroutine within this module
 type constants_data_type
  real*8      :: pi
  real*8      :: pi_sqrt
  real*8      :: conv_kJmol_ang2ps2gmol       ! converts kJ/mol to A^2/ps^2*g/mol
  real*8      :: conv_e2A_kJmol               ! converts e^2/A to kJ/mol
  real*8      :: boltzmann                    ! kB, kJ/mol/K
  real*8      :: friction_coeff               ! Friction coefficient for Langevin integrator, 1/ps
 end type constants_data_type

 Type(constants_data_type) :: constants


 !********************* defined type for verlet list data
 type verlet_list_data_type
  ! These are atom-atom based Verlet list, based on intermolecular interactions.
  ! intra-molecular atom-atom interactions are not included in these lists
  integer,dimension(:), allocatable             :: neighbor_list
  integer,dimension(:), allocatable             :: verlet_point
  integer                                   :: verlet_atoms  ! number of atoms in verlet list
  real*8                                    :: verlet_cutoff
  real*8                                    :: safe_verlet 
  real*8                                    :: verlet_thresh  ! see subroutine update_verlet_displacements for use of this variable
  integer                                   :: na_nslist      ! grid dimensions for neighbor searching
  integer                                   :: nb_nslist
  integer                                   :: nc_nslist  
  real*8, dimension(:,:), allocatable           :: verlet_xyz_store
  real*8, dimension(:,:), allocatable           :: verlet_displacement_store
 end type verlet_list_data_type



 !********************** defined type for PME data structures
 type PME_data_type
  ! alpha_sqrt controls the width of the smearing Gaussians in the Ewald sum.  As always, changing this value
  ! determines the relative magnitudes of reciprocal and real space sums, and therefore appropriate real space cutoff and PME grid size
  real*8                               ::  alpha_sqrt   ! Gaussian width parameter in A^-1,
  integer                              ::  pme_grid     ! this is the PME grid size in each dimension
  integer                              ::  spline_order ! this is the order of Beta-splines used in the PME charge interpolation
  real*8                               ::  Ewald_self                ! units: e^2/A
  real*8,dimension(:,:,:),allocatable  ::  CB, Q_grid,theta_conv_Q
  real*8,dimension(:,:,:),allocatable  ::  dQ_dr
  integer,dimension(:,:,:),allocatable ::  dQ_dr_index
  real*8,dimension(:,:),allocatable    ::  force_recip
  TYPE(DFTI_DESCRIPTOR), POINTER       ::  dfti_desc,dfti_desc_inv  ! for MKL FFT
  real*8                               ::  E_recip
  integer                              ::  spline_grid
  integer                              ::  erfc_grid
  real*8                               ::  erfc_max        ! this is max value up to which erfc is grid
  real*8, dimension(:),allocatable     ::  B6_spline,B5_spline,B4_spline,B3_spline
  real*8, dimension(:),allocatable     ::  erfc_table
 end type PME_data_type



 !*********************** defined type to store bond pair indices
 type bond_list_type
  integer                              ::  i_atom, j_atom 
 end type bond_list_type

!*********************** defined type to store angle atom indices
 type angle_list_type
  integer                              ::  i_atom, j_atom, k_atom
 end type angle_list_type

!*********************** defined type to store dihedral atom indices
 type dihedral_list_type
  integer                              ::  i_atom, j_atom, k_atom, l_atom
 end type dihedral_list_type

  ! exclusions for intra-molecular interactions
  ! we could easily define molecule-type specific number of exclusions, if we
  ! need to do this, put n_exclusions in molecule_type_data data structure and define
  ! for each molecule before setting molecule type exclusions
  integer       :: n_exclusions   ! this is the number of bonds away to exclude intra-molecular interactions

 !*********************** defined type to store atom and connectivity info for each molecule type
 type molecule_type_data_type
   integer :: n_atom                  ! number of atoms in this molecule type
   character(MAX_MNAME) :: mname
   integer, dimension(:), allocatable  :: atom_type_index  ! this is index of atom_type to look up force field parameters for this atom
   type(bond_list_type), dimension(:), allocatable :: bond_list
   type(angle_list_type), dimension(:), allocatable :: angle_list
   type(dihedral_list_type), dimension(:), allocatable :: dihedral_list
   integer, dimension(:,:), allocatable  :: pair_exclusions  ! list of all intra-molecular exclusions, generated based on the setting of n_excl.  Excluded atom pairs are marked with a "1", non-exclusions with a "0", and special 1-4 interactions labeled with a "2"
   ! these arrays are for MS-EVB, storing reactive atom indices. "1" indicates reactive proton (or basic atom).  Note that acids and bases can have both reactive protons and reactive basic atoms
   integer, dimension(:),allocatable  :: evb_reactive_protons
   integer, dimension(:),allocatable  :: evb_reactive_basic_atoms
 end type molecule_type_data_type


 ! data structure of molecule types, we don't make this allocatable as it's
 ! impractical to figure out the number of unique molecule types without having a
 ! data structure to store molecule_type_information first
 ! also, for MS-EVB simulation, we could have more molecule types than are in
 ! the .gro input file, as we need types for conjugate acids/bases
 type(molecule_type_data_type), dimension(MAX_N_MOLE_TYPE) :: molecule_type_data



 !********************************************* global atomtype data structures for force field **************************************************************
  integer:: n_atom_type, n_molecule_type
  character(MAX_ANAME), dimension(MAX_N_ATOM_TYPE) :: atype_name
  real*8, dimension(MAX_N_ATOM_TYPE) :: atype_chg
  integer, dimension(MAX_N_ATOM_TYPE) :: atype_freeze  ! this flag is for freezing atom types
  real*8, dimension(MAX_N_ATOM_TYPE) :: atype_mass    
  real*8, dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE,6) :: atype_vdw_parameter        ! for SAPT-FF, store A,B,C6,C8,C10,C12, for LJ, store epsilon,sigma
  real*8, dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE,9) :: atype_vdw_tmp
  integer, dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE) :: atype_vdw_type
  real*8, dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE,6) :: atype_vdw_parameter_14        ! this is same as atype_vdw_parameter array, but stores special values for 1-4 interactions as used in GROMOS-45a3 force field
  integer, dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE) :: atype_sapt_exclusions !This structure identifies the H2O-H3O parameters that need to be zeroed out for SAPT interactions
  character(10)            :: lj_comb_rule            ! for lj, set to "opls" or "standard"== Lorentz-Berthelot, for bkghm, set to "standard" or "ZIFFF"
! intra-molecular bond, angle and dihedral data structures
  integer, dimension(MAX_N_ATOM_TYPE, MAX_N_ATOM_TYPE)  :: atype_bond_type  ! stores type of bond, 1=harmonic, 2=GROMOS-96, 3=Morse
  real*8, dimension(MAX_N_ATOM_TYPE, MAX_N_ATOM_TYPE,3) :: atype_bond_parameter ! stores intra-molecular bond parameters for pair of atomtypes. For harmonic bond, first parameter is b0 (angstroms), second is kb (Kj/mol)/angstrom^2 , for Morse potential, first parameter is D (kJ/mol) , second is beta(angstroms^-1), third is b0(angstroms)
  integer, dimension(MAX_N_ATOM_TYPE, MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE)  :: atype_angle_type ! stores type of angle, 1=harmonic, 2=GROMOS96
  real*8, dimension(MAX_N_ATOM_TYPE, MAX_N_ATOM_TYPE, MAX_N_ATOM_TYPE,2) :: atype_angle_parameter ! stores intra-molecular angle parameters for set of three atomtypes.  First parameter is th0 (degrees), second is cth (Kj/mol)/ degree^2
  integer, dimension(MAX_N_ATOM_TYPE, MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE)  :: atype_dihedral_type ! stores type of dihedral, 1=proper, 2=improper, 3=ryckaert-belleman
  real*8, dimension(MAX_N_ATOM_TYPE, MAX_N_ATOM_TYPE, MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE,6) :: atype_dihedral_parameter ! stores intra-molecular dihedral parameters, for improper:first parameter is xi0, second is kxi ; for proper, first is phi0, second is kphi, third is multiplicity ; for ryckaert-bellemans, first is C0, second is C1, third is C2, fourth is C3, fifth is C4 and sixth is C5 

!*******************************************************************************************************************************************


  ! C6-C12 dispersion damping function tables
  ! the value of the 6th order Tang Toennies damping function at x=30 is
  ! 0.99999988268 , 8th order is 0.999997953924096, 10th order is 0.9999776512,
  ! 12th order is 0.9998323
  ! This should be a fine max values, since at x=30, the C10/R^10, C12/R^12
  ! terms should be basically zero for the corresponding distance
  real*8,parameter:: Tang_Toennies_max=50d0  
  integer,parameter:: Tang_Toennies_grid=1000 !Tang_Toennies_grid=1000000
  real*8,dimension(4,Tang_Toennies_grid) :: Tang_Toennies_table, dTang_Toennies_table ! C6 is in 1st index, C8 in 2nd, C10 in 3rd, C12 in 4th
 

! if trajectory is being restarted, set yes
  character(3) :: restart_trajectory
  integer  :: trajectory_step  ! this is the absolute trajectory step, based on the number of steps in this run, plus the number of steps from the continuation file
  integer  :: n_old_trajectory ! this is the number of  steps in the old trajectory
! if checkpointing velocity
  character(3) :: checkpoint_velocity
  integer      :: n_step_velocity
! cutoff for real-space interactions (VDWs, PME)
  real*8  :: real_space_cutoff 
  real*8,dimension(12) :: factorial    ! this saves the first 12 factorials which are used in Tang-Toennies damping functions


  ! parameters for parallel openmp
  integer:: n_threads

  ! debug: 0= normal, 1= print force calculation time decomposition for pme, 2=general time decomposition, options > 10 are for testing certain parts of code.  See specific test jobs.
  integer :: debug
  character(8)::date
  character(10)::time




  contains

  !***************************************************
  ! this subroutine fills in constant parameters to defined
  ! data structures.  We would prefer to just assign the parameters
  ! in the data type definitions, but Fortran won't let you do this...
  !***************************************************
  subroutine initialize_constants( file_io_data , verlet_list_data , PME_data )
  type(file_io_data_type), intent(inout) :: file_io_data
  type(verlet_list_data_type), intent(inout) :: verlet_list_data
  type(PME_data_type), intent(inout)         :: PME_data

  ! first fill in constants which is a global variable in this module

  constants%pi=3.141592654d0
  constants%pi_sqrt=1.772453851d0
  constants%conv_kJmol_ang2ps2gmol = 100d0  ! converts kJ/mol to A^2/ps^2*g/mol
  constants%conv_e2A_kJmol = 1389.35465     ! converts e^2/A to kJ/mol
  constants%boltzmann = 0.008314462d0       ! kB, kJ/mol/K
  constants%friction_coeff = 0.1d0            ! Friction coefficient for the Langevin integrator, 1/ps

  ! fill in verlet list parameters
  verlet_list_data%safe_verlet= 1.2 ! need bigger for long lamallae 1.6
  verlet_list_data%verlet_thresh =1.2d0  ! see subroutine update_verlet_displacements for use of this variable

  ! fill in sizes of PME lookup tables
  PME_data%spline_grid=100000
  PME_data%erfc_grid=1000000      !1000000
  PME_data%erfc_max=10d0       ! this is max value up to which erfc is grid


  ! name of velocity checkpoint file
  file_io_data%ifile_velocity = "velocity_checkpoint"  ! use this if restarting a simulation  

  ! now fill in file handles of file_io_data
  file_io_data%ifile_gro_file_h=90
  file_io_data%ifile_ffpmt_file_h=91
  file_io_data%ifile_top_file_h=92
  file_io_data%ifile_simpmt_file_h=93
  file_io_data%ifile_velocity_file_h=94
  file_io_data%ofile_traj_file_h=95
  file_io_data%ofile_log_file_h=96
  file_io_data%ofile_hop_file_h=97
  file_io_data%ofile_hamiltonian_file_h=98


  end subroutine initialize_constants


end module
