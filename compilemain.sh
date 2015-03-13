#!/bin/bash

# important:: use -openmp flag not only for explicitly parallel subroutines, but also for anything called by parallel code
# this flag enables automatic memory allocation of arrays, which is essential for subroutines and functions called by parallel code
# better to just use it globally

OMPINCLUDE=/usr/mpi/intel/openmpi-1.2.8/lib64/openmpi/
MKLLIB=/opt/intel/mkl/10.0.3.020/lib/em64t/

#OPT="-openmp -static -check bounds -check uninit -check format -warn declarations -traceback"  #-warn unused 

OPT="-openmp -static"

# note using the -openmp option will result in automatic memory allocation of all arrays, which for ifort could result in a stack overflow
# therefore if code segfaults on some machines, try ulimit -s unlimited

# -O0 option for parallel part
#/opt/intel/mkl/10.0.3.020/include/mkl_dfti.f90

ifort $OPT -c glob_v.f90 read_simulation_parameters.f90 general_routines.f90 rigid_body_kinematics.f90 -I$MKLLIB -I$OMPINCLUDE
ifort $OPT -c intra_bonded_interactions.f90 electrostatics.f90 pair_int.f90 -I$MKLLIB -I$OMPINCLUDE
ifort $OPT -c pme.f90 explicit_three_body_interaction.f90 -I$MKLLIB -I$OMPINCLUDE
ifort $OPT -c eq_drude.f90 total_energy_forces.f90 sapt_ff_routines.f90 expanded_grand_canonical_routines.f90 insertion_bias_routines.f90 penetration_ff_routines.f90 ms_evb.f90 mc_routines.f90 frprmn.f90 sampling.f90 main_mc.f90  -I$MKLLIB -I$OMPINCLUDE
ifort $OPT *.o -L$MKLLIB -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -L$OMPINCLUDE -libompitv -o main_mc

export OMP_STACKSIZE=20M
# -liomp5 -lpthread  use these if not using -openmp
