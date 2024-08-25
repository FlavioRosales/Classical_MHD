#FC = h5pfc
FC = mpif90
#FC=gfortran
#FC=/opt/intel/bin/ifort

FFLAGS = -g -O3
#FFLAGS = -g

LNK = mpif90
#LNK = h5pfc
#LNK=gfortran
#LNK=/opt/intel/bin/ifort

LIB=usr/lib/x86_64-linux-gnu/hdf5/openmpi
INC=usr/include/hdf5/openmpi

OBJS = main.o initial_conditions.o initialize_variables.o input_data.o output_data.o solve_Poisson.o system_evolve.o overall_diagnosis.o

MODS = \
../tools/mpi.o ../tools/FFT_1D.o ../tools/fft_lib.o ../tools/integrals.o ../tools/ODE.o ../tools/strings_lib.o ../tools/interpolation.o \
../tools/finite_differences_1d.o ../tools/finite_differences.o ../tools/differential_operators.o \
../tools/periodic_finite_differences_1d.o ../tools/periodic_finite_differences.o  ../tools/periodic_differential_operators.o \
../IO/save_ASCII.o ../IO/read_data.o ../IO/hdf5_lib.o \
../system/mesh.o \
../system/hydrodinamics/hydro_base.o ../system/hydrodinamics/hydro_diagnostics.o ../system/hydrodinamics/hydro_initial_data.o ../system/hydrodinamics/hydro.o \
global_numbers.o

$(OBJS):	$(MODS)

main:	$(OBJS) $(MODS)
		$(LNK) $(FFLAGS) -L/$(LIB) -o main $(OBJS) $(MODS) -lhdf5 -lhdf5_fortran
	@ mkdir -p exe
	@ mv main exe
#	@ cp input.par xxx


.PHONY:	clean

clean:
	-rm -f *.o exe/main ../IO/*.o ../system/*.o ../system/MoL/*.o ../system/Schrodinger/*.o ../system/hydrodinamics/*.o ../tools/*.o *.mod *.o

%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@ -I/$(INC)
