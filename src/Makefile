define modules =
	source /etc/profile.d/modules.sh
	module load intel/12.1.4
	module load openmpi/1.6
	module load mkl/10.0
endef

COMPILER=ifort
MKLMOD=./mkl_modules
OMPINCLUDE=/usr/local/pacerepov1/openmpi/1.6/intel-12.1.4/lib/openmpi/
MKLLIB=/usr/local/pacerepov1/intel/mkl/10.0.5.25/lib/em64t/

OPT=-openmp -static
#OPT=-Ofast -openmp -static
#OPT=-vec-report3 -openmp -static
#OPT=-O3 -qopenmp -static

# note using the -openmp option will result in automatic memory allocation of all arrays, which for ifort could result in a stack overflow
# therefore if code segfaults on some machines, try ulimit -s unlimited
#
# -liomp5 -lpthread  use these if not using -openmp

OBJ = glob_v.o read_simulation_parameters.o general_routines.o intra_bonded_interactions.o pair_int_real_space.o pme.o total_energy_forces.o ms_evb.o initialize_routines.o md_integration.o main_ms_evb.o

all: fast

debug: OPT += -O1 -g -check all -warn all -traceback -debug all
debug: main_ms_evb

# note fast includes -xHost ISA specific vectorization in addition to -O3 level optimization
fast: OPT += -Ofast
fast: main_ms_evb

%.o: %.f90 ; $(value modules)
	$(COMPILER) $(OPT) -c $< -I$(MKLLIB) -I$(MKLMOD) -I$(OMPINCLUDE)

main_ms_evb: $(OBJ) ; $(value modules)
	$(COMPILER) $(OPT) *.o -L$(MKLLIB) -lmkl_cdft_core -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -L$(OMPINCLUDE) -lompi_dbg_msgq -o main_ms_evb

.PHONY: clean modules

clean:
	rm *.o
	rm main_ms_evb

.ONESHELL:
