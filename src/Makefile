FC=gfortran
OUTPREFIX=../bin/
FFLAGS=-lhdf5 -lhdf5_fortran
FDEBUG=-g -Wall -Wextra -fcheck=all -fbacktrace
FOPTIM=-O3 -fopenmp
SRC= config_m.f90 diff.f90 evolvePDE.f90 helpers.f90 ogpf.f90 precision.f90 sparse_matrices.f90 lusol/src/lusol.f90\
data.f90
OBJ=$(SRC:.f90=.o)
FFLAGS+=$(FOPTIM)
LIB=usr/lib/x86_64-linux-gnu/hdf5/serial
INC=usr/include/hdf5/serial
INC2=./lusol/src

TESTS = test_diff test_evolvePDE test_helpers test_lusol test_sparse_matrices test_data

%.o: %.f90
	$(FC) $(FFLAGS) -I/$(INC) -I$(INC2) -c $<

tests: $(TESTS)
$(TESTS): $(OBJ)
	$(FC) $(OBJ) $(FFLAGS) -L/$(LIB) -I/$(INC) -I$(INC2) -o $(OUTPREFIX)$@ tests/$@.f90

#build: $(OBJ)
#	$(FC) $(FFLAGS) -L/$(LIB) -I/$(INC) -I$(INC2) -o bin/test_diff $(OBJ) tests/test_diff.f90


.PHONY: clean
clean:
	rm -f *.o *.mod

.PHONY: realclean
realclean:
	rm -f *.o *.mod ../bin/* ogpf_temp_script.gp
	cd lusol && make clean && cd ..


helpers.o: precision.o
config_m.o: helpers.o precision.o ogpf.o
sparse_matrices.o: helpers.o precision.o lusol/src/lusol.o
diff.o: sparse_matrices.o precision.o lusol/src/lusol.o
evolvePDE.o: config_m.o helpers.o ogpf.o diff.o precision.o lusol/src/lusol.o sparse_matrices.o data.o
lusol/src/lusol.o:
	cd lusol && make && cd ..
data.o: precision.o
