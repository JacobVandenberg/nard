FC=gfortran
OUTPREFIX=../bin/
#MKLROOT=/home/jacob/intel/oneapi/mkl/latest
FFLAGS=-lhdf5 -lhdf5_fortran
FDEBUG=-g -Wall -Wextra -fcheck=all -fbacktrace
FOPTIM=-O3 -fopenmp
SRC= config_m.f90 diff.f90 evolvePDE.f90 helpers.f90 ogpf.f90 precision.f90 sparse_matrices.f90 lusol/src/lusol.f90\
data.f90
OBJ=$(SRC:.f90=.o)
FFLAGS+=$(FOPTIM)
LIB=-L/usr/lib/x86_64-linux-gnu/hdf5/serial

LIB+= -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
FFLAGS +=  -fdefault-integer-8  -m64  -I"${MKLROOT}/include"

INC=-I/usr/include/hdf5/serial -I./lusol/src

TESTS = test_mkl test_diff test_evolvePDE test_helpers test_lusol test_sparse_matrices test_data test_bigdata test_config_m
NARD_EXE = nard2D

%.o: %.f90
	$(FC) $(FFLAGS) $(LIB) $(INC) -c $<


tests: $(TESTS)
$(TESTS): $(OBJ) ../bin
	$(FC) $(OBJ) $(FFLAGS) $(LIB) $(INC) -o $(OUTPREFIX)$@ tests/$@.f90

../bin:
	mkdir ../bin

.PHONY: nard
nard: $(NARD_EXE)
$(NARD_EXE): $(OBJ) ../bin
	#$(FC) $(OBJ) user_config.o $(FFLAGS) $(LIB) $(INC) -o $(OUTPREFIX)$@ $@.f90
	#rm user_config.o user_config.mod user_config.f90

.PHONY: new
new: nard user_functions.o
	$(FC) $(OBJ) user_functions.o $(FFLAGS) $(LIB) $(INC) -o $(OUTPREFIX)nard2D_$(ARGS) nard2D.f90
	rm -f user_functions.o user_functions.mod user_functions.f90

.PHONY: clean
clean:
	rm -f *.o *.mod user_config.o user_config.mod user_config.f90

.PHONY: realclean
realclean:
	rm -f *.o *.mod ogpf_temp_script.gp ../bin/* user_config.o user_config.mod user_config.f90
	cd lusol && make clean && cd ..


helpers.o: precision.o
config_m.o: helpers.o precision.o ogpf.o data.o
sparse_matrices.o: helpers.o precision.o lusol/src/lusol.o
diff.o: sparse_matrices.o precision.o lusol/src/lusol.o
evolvePDE.o: config_m.o helpers.o ogpf.o diff.o precision.o lusol/src/lusol.o sparse_matrices.o data.o
lusol/src/lusol.o:
	cd lusol && make && cd ..
data.o: precision.o
user_config.o: $(OBJ)