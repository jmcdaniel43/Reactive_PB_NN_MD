# Makefile created by mkmf $Id: mkmf,v 13.0 2006/04/10 21:20:01 fms Exp $ 

FC=ifort
LD=ifort
OMPINCLUDE=/usr/mpi/intel/openmpi-1.2.8/lib64/openmpi/
MKLLIB=/opt/intel/mkl/10.0.3.020/lib/em64t/
FFLAGS=-O3 -static -I${MKLLIB} -I${OMPINCLUDE} -openmp
LDFLAGS=-L${MKLLIB} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -L${OMPINCLUDE} -libompitv -static -O3

.DEFAULT:
	-touch $@
all: main_mc
glob_v.o: ./glob_v.f90
	$(FC) $(FFLAGS) -c	./glob_v.f90
sapt_ff_routines.o: ./sapt_ff_routines.f90 general_routines.o glob_v.o
	$(FC) $(FFLAGS) -c	./sapt_ff_routines.f90
general_routines.o: ./general_routines.f90 glob_v.o
	$(FC) $(FFLAGS) -c	./general_routines.f90
rigid_body_kinematics.o: ./rigid_body_kinematics.f90 general_routines.o glob_v.o
	$(FC) $(FFLAGS) -c	./rigid_body_kinematics.f90
read_simulation_parameters.o: ./read_simulation_parameters.f90 glob_v.o
	$(FC) $(FFLAGS) -c	./read_simulation_parameters.f90
insertion_bias_routines.o: ./insertion_bias_routines.f90 general_routines.o expanded_grand_canonical_routines.o glob_v.o total_energy_forces.o pair_int.o eq_drude.o electrostatics.o pme.o
	$(FC) $(FFLAGS) -c	./insertion_bias_routines.f90
main_mc.o: ./main_mc.f90 general_routines.o glob_v.o read_simulation_parameters.o pair_int.o total_energy_forces.o eq_drude.o sampling.o mc_routines.o insertion_bias_routines.o
	$(FC) $(FFLAGS) -c	./main_mc.f90
expanded_grand_canonical_routines.o: ./expanded_grand_canonical_routines.f90 general_routines.o glob_v.o
	$(FC) $(FFLAGS) -c	./expanded_grand_canonical_routines.f90
sampling.o: ./sampling.f90 general_routines.o rigid_body_kinematics.o glob_v.o total_energy_forces.o pair_int.o eq_drude.o pme.o insertion_bias_routines.o expanded_grand_canonical_routines.o electrostatics.o mc_routines.o frprmn.o ms_evb.o
	$(FC) $(FFLAGS) -c	./sampling.f90
total_energy_forces.o: ./total_energy_forces.f90 eq_drude.o pme.o explicit_three_body_interaction.o electrostatics.o pair_int.o glob_v.o intra_bonded_interactions.o
	$(FC) $(FFLAGS) -c	./total_energy_forces.f90
explicit_three_body_interaction.o: ./explicit_three_body_interaction.f90 glob_v.o general_routines.o
	$(FC) $(FFLAGS) -c	./explicit_three_body_interaction.f90
eq_drude.o: ./eq_drude.f90 pme.o electrostatics.o glob_v.o general_routines.o
	$(FC) $(FFLAGS) -c	./eq_drude.f90
mc_routines.o: ./mc_routines.f90 general_routines.o rigid_body_kinematics.o sapt_ff_routines.o glob_v.o total_energy_forces.o electrostatics.o pair_int.o pme.o penetration_ff_routines.o insertion_bias_routines.o eq_drude.o expanded_grand_canonical_routines.o ms_evb.o
	$(FC) $(FFLAGS) -c	./mc_routines.f90
pair_int.o: ./pair_int.f90 glob_v.o general_routines.o
	$(FC) $(FFLAGS) -c	./pair_int.f90
penetration_ff_routines.o: ./penetration_ff_routines.f90 general_routines.o glob_v.o
	$(FC) $(FFLAGS) -c	./penetration_ff_routines.f90
electrostatics.o: ./electrostatics.f90 general_routines.o glob_v.o
	$(FC) $(FFLAGS) -c	./electrostatics.f90
pme.o: ./pme.f90 general_routines.o glob_v.o electrostatics.o pair_int.o
	$(FC) $(FFLAGS) -c	./pme.f90
intra_bonded_interactions.o: ./intra_bonded_interactions.f90 glob_v.o general_routines.o
	$(FC) $(FFLAGS) -c	./intra_bonded_interactions.f90
frprmn.o: ./frprmn.f90 glob_v.o general_routines.o ms_evb.o
	$(FC) $(FFLAGS) -c	./frprmn.f90
ms_evb.o: ./ms_evb.f90 glob_v.o general_routines.o total_energy_forces.o
	$(FC) $(FFLAGS) -c	./ms_evb.f90
SRC = ./glob_v.f90 ./sapt_ff_routines.f90 ./general_routines.f90 ./rigid_body_kinematics.f90 ./read_simulation_parameters.f90 ./insertion_bias_routines.f90 ./main_mc.f90 ./expanded_grand_canonical_routines.f90 ./sampling.f90 ./eq_drude.f90 ./explicit_three_body_interaction.f90 ./total_energy_forces.f90 ./mc_routines.f90 ./pair_int.f90 ./penetration_ff_routines.f90 ./electrostatics.f90 ./pme.f90 ./intra_bonded_interactions.f90 ./ms_evb.f90 ./frprmn.o
OBJ = glob_v.o sapt_ff_routines.o general_routines.o rigid_body_kinematics.o read_simulation_parameters.o insertion_bias_routines.o main_mc.o expanded_grand_canonical_routines.o sampling.o eq_drude.o explicit_three_body_interaction.o total_energy_forces.o mc_routines.o pair_int.o penetration_ff_routines.o electrostatics.o pme.o intra_bonded_interactions.o ms_evb.o frprmn.o
clean: neat
	-rm -f .cppdefs $(OBJ) main_mc
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
main_mc: $(OBJ)
	$(LD) $(OBJ) -o main_mc $(LDFLAGS)
