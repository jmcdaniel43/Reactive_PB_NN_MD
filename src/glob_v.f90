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
  integer, parameter ::  MAX_N_MOLE_TYPE=10, MAX_N_ATOM_TYPE=20

  ! this is maximum number of atoms per molecule, which is used to set size of
  ! bond, angle, dihedral list arrays.  Probably we should get rid of this
  ! parameter and allocate necessary size at runtime...
  integer, parameter :: MAX_N_ATOM=30

  integer, parameter :: MAX_FN=100, MAX_ANAME=5, MAX_MNAME=5

! these variables determine whether to grid expensive functions in memory.  Code is much faster when these are set to 'yes'
  character(3)  :: grid_Tang_Toennies ! grid damping functions

  character(3), parameter  :: ms_evb_simulation="no"   ! ms_evb
  character(3), parameter  :: print_ms_evb_data = "yes"  ! if yes, this will print extra evb trajectory info

!***********************************************************************************************




  !***********************************************************
  !                         MS-EVB3 parameters
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
  ! these arrays store reactive atom indices. "1" indicates reactive proton (or basic
  ! atom).  Note that acids and bases can have both reactive protons and reactive basic atoms, so
  ! in general we need to include both acids and bases in both these data structures
  integer, dimension(MAX_N_MOLE_TYPE,MAX_N_ATOM)  :: evb_reactive_protons
  integer, dimension(MAX_N_MOLE_TYPE,MAX_N_ATOM)  :: evb_reactive_basic_atoms
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
   real*8, dimension(:,:), pointer :: xyz
   real*8, dimension(:,:), pointer :: velocity
   real*8, dimension(:,:), pointer :: force
   real*8, dimension(:), pointer   :: mass
   real*8, dimension(:), pointer   :: charge
   integer, dimension(:), pointer  :: atom_type_index  ! this is index of atom_type to look up force field parameters for this atom
   character(MAX_ANAME),dimension(:), pointer :: aname
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
  integer           :: ifile_gro_file_h  ! these are the file handles
  integer           :: ifile_ffpmt_file_h 
  integer           :: ifile_top_file_h 
  integer           :: ifile_simpmt_file_h 
  integer           :: ifile_velocity_file_h
  integer           :: ofile_traj_file_h 
  integer           :: ofile_log_file_h 
  integer           :: ofile_hop_file_h 
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
 end type constants_data_type

 Type(constants_data_type) :: constants


 !********************* defined type for verlet list data
 type verlet_list_data_type
  ! These are atom-atom based Verlet list, based on intermolecular interactions.
  ! intra-molecular atom-atom interactions are not included in these lists
  integer,dimension(:), pointer             :: neighbor_list
  integer,dimension(:), pointer             :: verlet_point
  integer                                   :: verlet_atoms  ! number of atoms in verlet list
  real*8                                    :: verlet_cutoff
  real*8                                    :: safe_verlet 
  real*8                                    :: verlet_thresh  ! see subroutine update_verlet_displacements for use of this variable
  integer                                   :: na_nslist      ! grid dimensions for neighbor searching
  integer                                   :: nb_nslist
  integer                                   :: nc_nslist  
  real*8, dimension(:,:), pointer           :: verlet_xyz_store
  real*8, dimension(:,:), pointer           :: verlet_displacement_store
 end type verlet_list_data_type


  !************************ Target variables for verlet_list_data_type pointers
  integer,dimension(:),allocatable, target :: neighbor_list
  integer,dimension(:),allocatable, target :: verlet_point
  real*8, dimension(:,:), allocatable, target :: verlet_xyz_store    ! store old positions to check for update
  real*8, dimension(:,:), allocatable, target :: verlet_displacement_store  !store old displacements to check for update


 !********************** defined type for PME data structures
 type PME_data_type
  ! alpha_sqrt controls the width of the smearing Gaussians in the Ewald sum.  As always, changing this value
  ! determines the relative magnitudes of reciprocal and real space sums, and therefore appropriate real space cutoff and PME grid size
  real*8                               ::  alpha_sqrt   ! Gaussian width parameter in A^-1,
  integer                              ::  pme_grid     ! this is the PME grid size in each dimension
  integer                              ::  spline_order ! this is the order of Beta-splines used in the PME charge interpolation
  real*8                               ::  Ewald_self                ! units: e^2/A
  real*8,dimension(:,:,:), pointer     ::  CB, Q_grid,theta_conv_Q
  real*8,dimension(:,:,:), pointer     ::  dQ_dr
  integer,dimension(:,:,:), pointer    ::  dQ_dr_index
  real*8,dimension(:,:) , pointer      ::  force_recip
  TYPE(DFTI_DESCRIPTOR), POINTER       ::  dfti_desc,dfti_desc_inv  ! for MKL FFT
  real*8                               ::  E_recip
  integer                              ::  spline_grid
  integer                              ::  erfc_grid
  real*8                               ::  erfc_max        ! this is max value up to which erfc is grid
  real*8, dimension(:), pointer        ::  B6_spline,B5_spline,B4_spline,B3_spline
  real*8, dimension(:), pointer        ::  erfc_table
 end type PME_data_type


 !************************ Target variables for PME_data_type pointers
  real*8,dimension(:,:,:),allocatable, target   ::  CB, Q_grid,theta_conv_Q
  real*8,dimension(:,:,:), allocatable, target  ::  dQ_dr
  integer,dimension(:,:,:), allocatable, target ::  dQ_dr_index
  real*8,dimension(:,:) , allocatable,   target ::  force_recip  ! stores reciprocal space PME forces for MS-EVB
  real*8,dimension(:), allocatable, target      ::  B6_spline,B5_spline,B4_spline,B3_spline
  real*8,dimension(:), allocatable, target      ::  erfc_table


 !********************************************* global data structures for force field **************************************************************
  integer:: n_atom_type, n_molecule_type
  character(MAX_ANAME), dimension(MAX_N_ATOM_TYPE) :: atype_name
  real*8, dimension(MAX_N_ATOM_TYPE) :: atype_chg
  integer, dimension(MAX_N_ATOM_TYPE) :: atype_freeze  ! this flag is for freezing atom types
  real*8, dimension(MAX_N_ATOM_TYPE) :: atype_mass    
  real*8, dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE,6) :: atype_vdw_parameter        ! for SAPT-FF, store A,B,C6,C8,C10,C12, for LJ, store epsilon,sigma
  real*8, dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE,6) :: atype_vdw_parameter_14        ! this is same as atype_vdw_parameter array, but stores special values for 1-4 interactions as used in GROMOS-45a3 force field
  integer, dimension(MAX_N_MOLE_TYPE,MAX_N_ATOM) :: molecule_type      ! this array contains indices for all types of solute molecules (not framework) end is marked with MAX_N_ATOM+1
  character(MAX_MNAME), dimension(MAX_N_MOLE_TYPE) :: molecule_type_name
  integer                 :: lj_bkghm                ! =1 for bkghm, =2 for lj
  character(10)            :: lj_comb_rule            ! for lj, set to "opls" or "standard"== Lorentz-Berthelot, for bkghm, set to "standard" or "ZIFFF"
 ! if lj_bkghm is set to 3, meaning we are using hybrid lj/bkghm force field, we need a 2nd combination rule.  The first "lj_comb_rule" will then be used for the bkghm force field for solute-framework interactions, and "lj_comb_rule2" will be used for the lj force field for solute-solute interactions
  character(10)            :: lj_comb_rule2            ! this will only be used if lj_bkghm=3
  integer, dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE) :: lj_bkghm_index ! this  maps which atom-atom types use a lj interaction, and which use a buckingham interaction, which is really only necessary for lj_bkghm=3, but is used always

! intra-molecular bond, angle and dihedral data structures
  integer, dimension(MAX_N_ATOM_TYPE, MAX_N_ATOM_TYPE)  :: atype_bond_type  ! stores type of bond, 1=harmonic, 2=GROMOS-96, 3=Morse
  real*8, dimension(MAX_N_ATOM_TYPE, MAX_N_ATOM_TYPE,3) :: atype_bond_parameter ! stores intra-molecular bond parameters for pair of atomtypes. For harmonic bond, first parameter is b0 (angstroms), second is kb (Kj/mol)/angstrom^2 , for Morse potential, first parameter is D (kJ/mol) , second is beta(angstroms^-1), third is b0(angstroms)
  integer, dimension(MAX_N_ATOM_TYPE, MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE)  :: atype_angle_type ! stores type of angle, 1=harmonic, 2=GROMOS96
  real*8, dimension(MAX_N_ATOM_TYPE, MAX_N_ATOM_TYPE, MAX_N_ATOM_TYPE,2) :: atype_angle_parameter ! stores intra-molecular angle parameters for set of three atomtypes.  First parameter is th0 (degrees), second is cth (Kj/mol)/ degree^2
  integer, dimension(MAX_N_ATOM_TYPE, MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE)  :: atype_dihedral_type ! stores type of dihedral, 1=proper, 2=improper
  real*8, dimension(MAX_N_ATOM_TYPE, MAX_N_ATOM_TYPE, MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE,3) :: atype_dihedral_parameter ! stores intra-molecular dihedral parameters, for improper:first parameter is xi0, second is kxi ; for proper, first is phi0, second is kphi, third is multiplicity
  integer,dimension(MAX_N_MOLE_TYPE, MAX_N_ATOM, MAX_N_ATOM)   :: molecule_bond_list  ! this stores the bond list for each molecule type.  If there is a bond between the 2nd and 3rd atoms in the 1st molecule type, then the entry (1,2,3)=1, otherwise (1,2,3)=0
  integer,dimension(MAX_N_MOLE_TYPE, MAX_N_ATOM, MAX_N_ATOM, MAX_N_ATOM)   :: molecule_angle_list  ! this stores the angle list for each molecule type.  data storage is same format as molecule_bond_list
  integer,dimension(MAX_N_MOLE_TYPE, MAX_N_ATOM, MAX_N_ATOM, MAX_N_ATOM, MAX_N_ATOM)   :: molecule_dihedral_list  ! this stores the dihedral list for each molecule type.  data storage is same format as molecule_bond_list
  integer,dimension(MAX_N_MOLE_TYPE) :: molecule_dihedral_flag  ! =1, if molecule has any dihedrals, =0 otherwise.  This is used to avoid searching over atoms in molecules (water) that don't have any dihedrals

  ! exclusions for intra-molecular interactions
  integer, parameter  :: n_excl = 2  ! this is the number of bonds away to exclude intra-molecular interactions
  integer, dimension(MAX_N_MOLE_TYPE, MAX_N_ATOM, MAX_N_ATOM) :: molecule_exclusions ! this has a list of all atom-atom exlusions for each molecule type.  This is generated based on the setting of n_excl.  Excluded atom pairs are marked with a "1", non-exclusions with a "0".  In addition, this array has a special labeling for 1-4 interactions, which are labeled with a "2"

!*******************************************************************************************************************************************


  ! C6-C12 dispersion damping function tables
  ! the value of the 6th order Tang Toennies damping function at x=30 is
  ! 0.99999988268 , 8th order is 0.999997953924096, 10th order is 0.9999776512,
  ! 12th order is 0.9998323
  ! This should be a fine max values, since at x=30, the C10/R^10, C12/R^12
  ! terms should be basically zero for the corresponding distance
  real*8,parameter:: Tang_Toennies_max=30d0  
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

  ! fill in verlet list parameters
  verlet_list_data%safe_verlet= 1.2 ! need bigger for long lamallae 1.6
  verlet_list_data%verlet_thresh =2d0  ! see subroutine update_verlet_displacements for use of this variable

  ! fill in sizes of PME lookup tables
  PME_data%spline_grid=100000
  PME_data%erfc_grid=1000000
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


  end subroutine initialize_constants


end module
