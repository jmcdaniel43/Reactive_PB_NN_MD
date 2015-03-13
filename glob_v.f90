!************************************************************
!
! This module contains global variables that are used in many subroutines
! throughout the Monte Carlo program
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
implicit none

!***************************** NOTE THIS SECTION ***********************************************
! here we organize all global variable parameters that are hard-coded in, that the user
! may want to adjust in certain cases

  ! these parameters control the size of data arrays.  Make sure MAX_N_MOLE > maximum number of molecules in the simulations
  ! and MAX_N_ATOM > maximum number of atoms per molecule
   integer, parameter :: MAX_FN=100, MAX_ANAME=5, MAX_MNAME=5, MAX_N_MOLE=4000, MAX_N_ATOM=30 , MAX_N_MOLE_TYPE=3, MAX_N_ATOM_TYPE=8

  ! these parameters control the Drude oscillator convergence:  Maximum number of iterations and convergence threshold
  integer,parameter::iteration_limit = 15    ! if scf iterations exceed this limit, reject move
  real*8,parameter:: force_threshold = 4.0477D-6  ! this is for rms force in e^2/A^2, equal to .1 KJ/mol/nm

  ! this shifts the vdw part of the force field to zero at the cutoff
  integer, parameter :: lj_shift  = 0       ! 0 = no shift, 1 = shift to zero at cutoff

! interpolation for PME dispersion grid.
  integer, parameter :: lgrg_order=2 ! if this equals 1, then the grid is interpolated linearly

! dhf combination rule:  As we used sort of a weird convention for the sign of the dhf cross terms in the beginning, but changed to an completely attractive dhf cross term convention when fitting virials, we need to implement a variety of combination rules
! here "1" is as in JPC C, 2012, 116, 1892-1903
! "2" is Aij = positive iff Aii and Ajj positive (difference between 1 and 2 is Aij negative in 2 if Aii and Ajj both negative
! "3" is always attractive
integer, parameter :: dhf_combination_rule=2

! these variables determine whether to grid expensive functions in memory.  Code is much faster when these are set to 'yes'
character(3)  :: grid_erfc = "yes"        ! grid error function for pme
character(3)  :: grid_Tang_Toennies ! grid damping functions

character(3)  :: ms_evb_simulation="no"   ! ms_evb

character(3)  :: flexible_molecule_simulation="yes"


!***********************************************************************************************




  !***********************************************************
  !                         MS-EVB3 parameters
  ! from JPC B, 2008, 112, 467-482 (and errata, JPC B, 2008, 112, 7146)
  ! "u" stands for upper case, "l" stands for lower case, "p" stands for prime
  ! "s" stands for star
  ! I've tried to make parameter names identical to that in the reference
  ! all parameters have a "_evb" suffice
  ! units are same as rest of the code kJ/mol, Angstrom, etc
  !
  real*8, parameter :: Bu_evb=47.1118977d0 , bl_evb=1.1d0 , blp_evb=2.12d0 , dlOO0_evb=2.4d0
  real*8, parameter :: Cu_evb=19.1274639d0 , cl_evb=2.1d0 , dlOH0_evb=1.0d0 , Vconstij_evb=-97.0151921d0
  real*8, parameter :: qOex_evb=-0.0895456d0 , qHex_evb=0.0252683d0 , qHexs_evb=0.078018d0 , gamma_evb=1.8302895d0
  real*8, parameter :: Pu_evb=0.232726d0 , k_evb=9.562153d0 , DuOO_evb=2.94d0 , beta_evb=6.0179066d0 
  real*8, parameter :: RuOO0_evb=3.1d0, Pup_evb=10.8831327d0 , alpha_evb=10.0380922d0, rlOO0_evb=1.8136426d0
  ! these are the rs, rc parameters in the switching function (eqn 9) on page 470
  ! of JPC B, 2008, 112, 467-482
  real*8, parameter :: rshiftOO_evb=2.85d0, rcutoffOO_evb=3.05d0, rshiftOH_evb=2.3d0, rcutoffOH_evb=2.5d0

  ! because of the special treatment of hydronium molecules, reactive protons, etc.,
  ! we require special labels to identify these atom types
  character(MAX_ANAME) :: hydronium_proton_label="H1", hydronium_oxygen_label="O1"
  character(MAX_ANAME) :: water_proton_label="H0", water_oxygen_label="O0"

  integer , parameter :: n_proton_max =10
  ! this keeps track of the hydronium molecules
  integer        :: n_hydronium_molecules
  integer, dimension( n_proton_max ) :: hydronium_molecule_index

  !***********************************************************




  real*8, parameter :: pi=3.14159265
  real*8, parameter :: pi_sqrt=1.7724538509d0
  real*8,dimension(12) :: factorial    ! this saves the first 12 factorials which are used in Tang-Toennies damping functions


! framework
  integer            :: framework_simulation            ! tells subroutines whether or not this simulation includes fixed framework atoms

! code version
  integer  :: code_version = 30

! if yes, this will not run a simulation, but calculate the interaction energy between a solute molecule and the framework over specified grid
 character(3) :: energy_grid_run = "no"    ! will just call energy_grid subroutine

! for puddle filling
character(3) :: puddle_filling_on = "no"
real*8  :: puddle_min = -15.
real*8  :: puddle_delta_E

! Feynmann Hibbs quantum correction for hydrogen
character(3) :: Feynmann_Hibbs_correction
character(3) :: Feynmann_Hibbs_forces

! this changes the input data structure, as now the code must be able to handle starting with zero molecules, also it reads in a random number seed.
  character(3) :: replica_input = "no"
  integer,dimension(:),allocatable :: seed_replica

! choose ensemble
  character(3) :: select_ensemble        ! nvt, npt, uvt, egc

! input file type is determined by whether or not a sapt_type force field is being used
 character(3) :: sapt_type_ff
 character(3) :: energy_decomposition     ! an explicit energy decomposition will be given if using a (purely) sapt-based force field ( and without pme dispersion )


! penetration force field for short distances , this will only be used with energy decomposition
character(3) :: penetration_force_field
character(3) :: switch_at_minimum = "no"
real*8 :: penetration_force_field_threshold

! if we want to use a different solute parameter set to create solute framework cross terms, we have to read this in separately
 character(3) :: solute_cross_parameter_set 

! if, energy decomposition, have option of using C8,C10 parameters for dispersion.  This also determines input file structure
 character(3) :: C8_10_dispersion_terms  ! if yes, damping for dispersion will be turned on, and buckingham coefficients for dispersion will be read in
 character(3) :: C12_dispersion

 ! three body force field flags
 character(3) :: three_body_dispersion
 character(3) :: three_body_exchange


! if using pme dispersion
 character(3) :: pme_disp
 real*8, dimension(:,:,:,:), allocatable :: lrdisp_pot

! explicitly read solute-framework exponents
 character(3) :: explicit_cross_term_exponents="no"

! flag tells code whether to damp solute-solute dispersion interactions
 character(3) :: damp_solute_solute_dispersion

! hybrid mc/md or regular mc?
character(3)  :: hybrid_md_mc         ! "no" for standard mc translation/rotation moves (very slow for drudes!), "yes" for hybrid_md_mc

!force field parameters
  integer:: n_atom_type_store  ! used for expanded GCMC where we have to introduce additional atom types
  integer:: n_atom_type, n_molecule_type
  character(MAX_ANAME), dimension(MAX_N_ATOM_TYPE) :: atype_name
  character(MAX_MNAME), dimension(MAX_N_MOLE) :: molecule_name
  real*8, dimension(MAX_N_ATOM_TYPE) :: atype_chg,atype_pol
  real*8, dimension(MAX_N_ATOM_TYPE) :: atype_mass    
  real*8, dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE,5) :: atype_bkghm_decomp                   ! this stores "A" type coefficients for Bkingham, index 1-5 are exch,elec,induc,dhf,dispersion, respectively
  real*8, dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE,6) :: atype_lj_parameter        ! for buckingham, store A,B,C6 (possibly C8,C10), for LJ, store epsilon,sigma ; index 6 is for possible C12 terms
  real*8, dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE,6) :: atype_lj_parameter_14        ! this is same as atype_lj_parameter array, but stores special values for 1-4 interactions as used in GROMOS-45a3 force field
  real*8, dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE) :: atype_3body_C9            ! this stores C9 coefficients, for 3-body dispersion
  real*8, dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE,2) :: atype_3body_exchange      ! this stores parameters, for 3-body exchange
  real*8, dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE,4) :: atype_penetration_ff      ! penetration force field, first index is bkghm coeff, 2 index is distance cutoff, 3 index is energy shift, 4 index is exponent scale factor
  integer, dimension(MAX_N_ATOM_TYPE) ::  atype_solute                      ! keeps track of whether this atom type is framework(0) or solute(1)
  integer, dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE) ::  atype_damp_dispersion   ! keeps track of whether we damp the dispersion interaction between these atom types
  integer, dimension(MAX_N_MOLE_TYPE,MAX_N_ATOM) :: molecule_type      ! this array contains indices for all types of solute molecules (not framework) end is marked with MAX_N_ATOM+1
  character(MAX_MNAME), dimension(MAX_N_MOLE_TYPE) :: molecule_type_name
!  real*8, dimension(MAX_N_MOLE,3,3) :: molecule_type_inertia 
  real*8, dimension(MAX_N_MOLE_TYPE,MAX_N_ATOM,3) :: r_mole_com     ! stores configuration for all molecule types. for non-linear molecules, configurations are stored in principle coordinate basis
  integer, dimension(MAX_N_MOLE_TYPE) :: molecule_shape      ! flag "0" for nonlinear, polyatomic, "1" for linear molecules, "2" for atom
  real*8,dimension(MAX_N_MOLE_TYPE,3) :: molecule_inertia
  real*8,dimension(MAX_N_MOLE_TYPE) :: molecule_mass
  integer, dimension(MAX_N_MOLE,MAX_N_ATOM):: atom_index            ! integer specifies where to look up atom parameters in lj_parameter array
  integer, dimension(MAX_N_MOLE):: molecule_index            ! integer indexes molecule type for all solute molecules in simulation
  integer, dimension(MAX_N_MOLE,MAX_N_ATOM):: drude_atom_map            ! this array maps which atom each Drude oscillator resides on
  integer                 :: lj_bkghm                ! =1 for bkghm, =2 for lj
  character(10)            :: lj_comb_rule            ! for lj, set to "opls" or "standard"== Lorentz-Berthelot, for bkghm, set to "standard" or "ZIFFF"
 ! if lj_bkghm is set to 3, meaning we are using hybrid lj/bkghm force field, we need a 2nd combination rule.  The first "lj_comb_rule" will then be used for the bkghm force field for solute-framework interactions, and "lj_comb_rule2" will be used for the lj force field for solute-solute interactions
  character(10)            :: lj_comb_rule2            ! this will only be used if lj_bkghm=3
  integer, dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE) :: lj_bkghm_index ! this  maps which atom-atom types use a lj interaction, and which use a buckingham interaction, which is really only necessary for lj_bkghm=3, but is used always

! intra-molecular bond, angle and dihedral data structures
! note the mapping for each molecule to the molecule type index
! used to access the particular molecular parameters is contained
! in the molecule_index array
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

  ! this are for gridding the ewald real space atom type interactions
  ! note that these are atom type specific because of the damping functions
  character(3) :: grid_ewald_realspace_interaction="no"
  integer,parameter :: ewald_realspace_interaction_grid_size = 100 ! ewald_realspace_interaction_grid_size = 1000000
  real*8, dimension(:,:), allocatable :: ewald_realspace_interaction_grid
  integer,dimension(MAX_N_ATOM_TYPE, MAX_N_ATOM_TYPE) :: ewald_realspace_interaction_grid_map
  real*8 :: ewald_realspace_table_pad = 10d0 ! this is padding beyond ewald_cutoff for table construction which is used because cutoff is applied on a molecule by molecule basis, and so atom-atom distances can be farther than this
  real*8 :: max_ewald_table



! mc parameters, specified here temporarilly
  integer   :: n_update_stepsize, n_output, n_step
  real*8 :: temp  
  integer :: cycles_per_vol_step  ! this is the average cycles ( 1 hybrid move, or nmole single particle_moves)/vol_move
  real*8 :: too_close         ! if any atoms on different molecules closer than this, reject move right away
  real*8 :: delta_t1           ! this is time step for integrator in hybridmc (ps)
  real*8 :: vol_step           ! this is (twice the maximum) percentage change in ln(volume) (.05 is reasonable)
  real*8 :: max_delta_t1 , max_vol_step , target_acc
  real*8 :: target_acc_v      ! desired percentage accepted moves
  real*8, parameter:: step_update=.1,step_update_v=.05  ! percentage by which to change step sizes
  real*8 :: pressure   ! pressure in bar
  real*8, parameter :: p_conv=(1D-03/16.6054)  ! bar to kj/mol/A^3
! Buckingham related parameters
  real*8  :: lj_cutoff 
  integer  :: lj_lrc          ! 0 = no lrc, 1 = lrc (lrc will never be computed if there is framework)
  real*8, dimension(MAX_N_ATOM_TYPE,MAX_N_ATOM_TYPE) :: lj_shift_store

! transformation matrix to box vectors
  real*8,dimension(3,3) :: xyz_to_box_transform

! These are atom-atom based Verlet list, based on intermolecular interactions.
! intra-molecular atom-atom interactions are not included in these lists
! because this is an atom-atom list, we need an array that maps each atom
! of the molecule to the atom in the list.  These arrays are
! verlet_molecule_map_lj and verlet_molecule_map_elec
! here we have separate verlet lists for lj and electrostatic subroutines,
! because some atoms may have zero charge, and others may have no lj sites,
! and we want to have different lists to exclude such atoms
  integer      :: total_atoms  ! total number of atoms in simulation
  character(3) :: use_verlet_list
  character(3),parameter :: verlet_grid_based_construction="yes"
  integer,dimension(:),allocatable :: verlet_neighbor_list_lj , verlet_neighbor_list_elec
  integer,dimension(MAX_N_MOLE,MAX_N_ATOM) :: verlet_molecule_map_lj, verlet_molecule_map_elec
  integer,dimension(:),allocatable :: verlet_point_lj, verlet_point_elec
  integer                       :: verlet_lj_atoms, verlet_elec_atoms ! these are the number of different atoms in each verlet list
  real*8       :: verlet_cutoff_lj , verlet_cutoff_elec
  real*8,parameter :: safe_verlet=1.2  ! this is factor that we multiply theoretically needed size of verlet list to be safe
  real*8,parameter :: verlet_thresh =2d0  ! see subroutine update_verlet_displacements for use of this variable

! single molecule moves
  real*8 :: pct_singlemv_hybrid = .9
  integer, parameter :: translation_method= 1 ! 1 for displacement, 2 for lattice hop
  integer, parameter :: rotate_method = 2  ! 1 for linear, 2 for general euler angles
  integer, parameter :: reset_trans_rot = 10 ! see subroutine update_stepsize_trans_rot
  real*8 :: single_mol_step =.2 , max_single_mol_step = 4d0
  real*8 :: single_mol_rotate = .2 ,max_single_mol_rotate = 5d0 ! for type 1 rotation
  real*8 :: phi_step = .1 , psi_step = .1 , theta_step = .05 ! for type 2 rotation
  real*8 :: max_phi_step = 2.6 , max_psi_step = 2.6 , max_theta_step = 1.3 ! for type 2 rotation

! parameters for Electrostatic (Ewald/Cutoff/pme)
  integer, parameter :: kmax=100, maxk=400000, maxn=1000
  real*8, parameter :: ksqmax=10. ! ksqmax is cutoff for k vectors
  real*8            :: ewald_cutoff
  real*8            :: Electro_cutoff
  integer           :: screen_type           ! screening for coulomb potential 0 = none, 1 = Tang-Toennies type
  real*8,parameter  :: screen_distance_max_sq = 50d0 ! this is the maximum distance (squared) for which the screening function will be used in the electrostatic force calculation, since it is expensive, and it is extremely close to 1 at this point
  ! PME parameters.  alpha_sqrt controls the width of the smearing Gaussians in the Ewald sum.  As always, changing this value
  ! determines the relative magnitudes of reciprocal and real space sums, and therefore appropriate real space cutoff and PME grid size
  real*8             ::  alpha_sqrt   ! Gaussian width parameter in A^-1,
  integer            ::  pme_grid     ! this is the PME grid size in each dimension
  integer            ::  spline_order ! this is the order of Beta-splines used in the PME charge interpolation

  character(10)     :: electrostatic_type     ! choose from pme, cutoff, or none. ewald is only implemented for testing
  ! global variables in Ewald calculations
  integer :: totk
  character(3) :: ewald_volume_move ! this flag signals which global variables the ewald function should update
  complex*16, dimension(maxk) :: rou_k, rou_k_try
  real*8, dimension(maxk,3) :: k_vec,new_k_vec
  real*8, dimension(maxk) :: B_factor,new_B_factor
  real*8 :: Ewald_self, Ewald_framework ! These two have different units: e^2/A
  real*8 :: E_disp_lrc,E_disp_lrc_try

  ! parameters for pme
  real*8,dimension(:,:,:),allocatable::CB,Q_grid,theta_conv_Q
  ! real*8,dimension(pme_grid,pme_grid,pme_grid)::Q_grid_try
  !real*8::pme_realE,pme_realE_try
  ! grid Bsplines,for order n, need n spline and n-1 spline grids
  ! B(i) corresponds to B(i*n/gridsize)
  integer,parameter::spline_grid=10000,erfc_grid=1000000, gfun_grid=100 !gfun_grid=500000
  real*8,parameter::erfc_max=10d0              ! this is max value up to which erfc is grid
  real*8,parameter::gfun_max=5.0d0
  real*8,dimension(spline_grid)::B6_spline,B5_spline,B4_spline,B3_spline
  real*8,dimension(erfc_grid):: erfc_table
  real*8,dimension(gfun_grid):: g6_table, g8_table, g10_table, g12_table, g14_table


  ! C6-C12 dispersion damping function tables
  ! the value of the 6th order Tang Toennies damping function at x=30 is 0.99999988268 , 8th order is 0.999997953924096, 10th order is 0.9999776512, 12th order is 0.9998323
  ! This should be a fine max values, since at x=30, the C10/R^10, C12/R^12 terms should be basically zero for the corresponding distance
  real*8,parameter:: Tang_Toennies_max=30d0  
  integer,parameter:: Tang_Toennies_grid=100 !Tang_Toennies_grid=1000000
  real*8,dimension(4,Tang_Toennies_grid) :: Tang_Toennies_table, dTang_Toennies_table ! C6 is in 1st index, C8 in 2nd, C10 in 3rd, C12 in 4th
  real*8,dimension(Tang_Toennies_grid) :: Tang_Toennies_table3  ! this is used for three-body dispersion damping functions

  ! parameters for drude oscillator
  real*8          :: springcon
  real*8          :: thole 
  integer        :: drude_simulation         ! this flag will let sampling routines know whether to call scf drude (=1 for drudes)


  ! parameters for three body dispersion
  real*8 :: three_body_dispersion_cutoff
  integer :: three_body_cutoff_type   ! this should be "0" if three body cutoff acts on molecules as a whole, or "1" if it is an atom by atom basis
  integer :: na_nslist, nb_nslist, nc_nslist   ! link list parameters
  integer, dimension(:,:,:), allocatable :: headlist
  integer, dimension(MAX_N_MOLE) :: nslist

  ! parameters for parallel openmp
  integer:: n_threads

  ! debug: 0= normal, 1= print force calculation time decomposition for pme, 2=general time decomposition, options > 10 are for testing certain parts of code.  See specific test jobs.
  integer :: debug
  character(8)::date
  character(10)::time

  ! Gibbs ensemble, need copies of certain variables for 2nd box (first box uses regular variables)
real*8 :: E_disp_lrc2,E_disp_lrc_try2
real*8,dimension(2) :: Ewald_self_gibbs
!real*8,dimension(pme_grid,pme_grid,pme_grid)::CB2,Q_grid2,theta_conv_Q2
real*8,parameter:: gibbs_vol_pct=0.2,gibbs_swap_pct=0.2        ! fraction of moves that change volume,swap particles
integer  :: orientation_try    ! number of orientation tries for particle insertion bias, =1 means unbiased

  ! grand canonical ensemble
real*8          :: gcmc_ins_rm        ! fraction of moves that attempt to insert or remove a particle, for expanded grand, this is attempts to change coupling
real*8          :: chem_potential
   
 ! expanded grand canonical ensemble
integer :: partial_molecule_index
integer :: partial_molecule_coupling = 0  ! our initially energy calculation will have all molecules at full coupling, so this should have an initial value of 0
integer,parameter :: partial_coupling_steps = 5
real*8  :: coupling_increment = 1d0/dble(partial_coupling_steps)
real*8,dimension(partial_coupling_steps) :: step_weight

! cavity-bias
character(3)      :: cavity_bias = "no"   ! yes or no
integer,parameter :: cav_grid_a = 26, cav_grid_b = 26 , cav_grid_c = 26
integer,parameter :: energy_bias = 0 ! 0= naive cavity bias, 1= lennard jones boltzmann factor
real*8,parameter :: min_frmwk_dist = 2d0                 ! min distance from center of cubelet to framework for cubelet to be considered for insertion
real*8           :: cubelet_vol                            ! calculated from box size and grid size

end module
