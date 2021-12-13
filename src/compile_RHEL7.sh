#!/bin/bash

#************* Use this compile script to compile on RHEL7 machines on PACE ***********************

# important:: use -openmp flag not only for explicitly parallel subroutines, but also for anything called by parallel code
# this flag enables automatic memory allocation of arrays, which is essential for subroutines and functions called by parallel code
# better to just use it globally

module load intel/19.0

COMPILER=mpiifort
MKLLIB="-L${MKLROOT}/lib/intel64"
MKLINC="-I./mkl_modules"
#FFLAGS="-O3 -qopenmp -mkl=sequential $MKLLIB -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -ldl $MKLINC"
FFLAGS="-O1 -qopenmp -check all -traceback -debug all -mkl=sequential $MKLLIB -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -ldl $MKLINC"

# use this set of flags for profiling function timings.  Need to not link with openmp, since profiler only works on serial code...
#FFLAGS="-O3 -profile-functions -mkl=sequential $MKLLIB -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -ldl $MKLINC"

echo MKLROOT: $MKLROOT
echo MKLLIB: $MKLLIB
echo FFLAGS: $FFLAGS

export KMP_AFFINITY=verbose,none

# note using the -openmp option will result in automatic memory allocation of all arrays, which for ifort could result in a stack overflow
# therefore if code segfaults on some machines, try ulimit -s unlimited

set -x # show the compilation steps

$COMPILER $FFLAGS -c glob_v.f90 read_simulation_parameters.f90 general_routines.f90
$COMPILER $FFLAGS -c intra_bonded_interactions.f90 pair_int_real_space.f90 pme.f90
$COMPILER $FFLAGS -c total_energy_forces.f90 ms_evb.f90 initialize_routines.f90 md_integration.f90 main_ms_evb.f90
$COMPILER $FFLAGS *.o $MKLLIB -o main_ms_evb

set +x # end
