#!/bin/bash

# important:: use -openmp flag not only for explicitly parallel subroutines, but also for anything called by parallel code
# this flag enables automatic memory allocation of arrays, which is essential for subroutines and functions called by parallel code
# better to just use it globally


# newer libraries, currently, these are not working on PACE...
#module load intel/15.0
#module load openmpi/1.8
#module load mkl/11.2
#OMPINCLUDE=/usr/local/pacerepov1/openmpi/1.8/intel-15.0/lib/openmpi/
#MKLLIB=/usr/local/pacerepov1/intel/mkl/11.2/lib/intel64/

# older intel libraries
module load intel/12.1.4
module load openmpi/1.6
module load mkl/10.0
OMPINCLUDE=/usr/local/pacerepov1/openmpi/1.6/intel-12.1.4/lib/openmpi/
MKLLIB=/usr/local/pacerepov1/intel/mkl/10.0.5.25/lib/em64t/


COMPILER=ifort
MKLMOD=./mkl_modules

# this first compiler option is for testing...
OPT="-O1 -g -openmp -static -check all -warn all -traceback -debug all"
#OPT="-vec-report3 -openmp -static"
#OPT="-O3 -qopenmp -static"
#OPT="-Ofast -openmp -static"  # note Ofast includes -xHost ISA specific vectorization in addition to -O3 level optimization

# note using the -openmp option will result in automatic memory allocation of all arrays, which for ifort could result in a stack overflow
# therefore if code segfaults on some machines, try ulimit -s unlimited

#/opt/intel/mkl/10.0.3.020/include/mkl_dfti.f90

$COMPILER $OPT -c glob_v.f90 read_simulation_parameters.f90 general_routines.f90 -I$MKLLIB -I$MKLMOD -I$OMPINCLUDE
$COMPILER $OPT -c intra_bonded_interactions.f90 pair_int_real_space.f90 pme.f90 -I$MKLLIB -I$MKLMOD -I$OMPINCLUDE
$COMPILER $OPT -c total_energy_forces.f90 ms_evb.f90 initialize_routines.f90 md_integration.f90 main_ms_evb.f90  -I$MKLLIB -I$MKLMOD -I$OMPINCLUDE
#$COMPILER $OPT *.o -L$MKLLIB -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -L$OMPINCLUDE -o main_ms_evb
$COMPILER $OPT *.o -L$MKLLIB -lmkl_cdft_core -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -L$OMPINCLUDE -lompi_dbg_msgq -o main_ms_evb
export OMP_STACKSIZE=20M
# -liomp5 -lpthread  use these if not using -openmp
