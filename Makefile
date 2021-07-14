
code: test_diff test_helpers test_evolvePDE test_sparse_matrices test_lusol
	echo "Finished"
test_evolvePDE: evolvePDE.o config_m.o helpers.o ogpf.o diff.o precision.o sparse_matrices.o lusol/src/lusol.o tests/test_evolvePDE.f90
	gfortran -O3 -floop-parallelize-all -I./lusol/src -o bin/test_evolvePDE precision.o evolvePDE.o lusol/src/lusol.o config_m.o helpers.o ogpf.o diff.o sparse_matrices.o tests/test_evolvePDE.f90 -fopenmp
test_diff: diff.o ogpf.o sparse_matrices.o lusol/src/lusol.o precision.o helpers.o tests/test_diff.f90
	gfortran -O3 -floop-parallelize-all -I./lusol/src -o bin/test_diff diff.o precision.o ogpf.o lusol/src/lusol.o sparse_matrices.o helpers.o tests/test_diff.f90 -fopenmp
test_sparse_matrices: sparse_matrices.o lusol/src/lusol.o helpers.o precision.o ogpf.o diff.o tests/test_sparse_matrices.f90
	gfortran -O3 -floop-parallelize-all -I./lusol/src -o bin/test_sparse_matrices precision.o diff.o lusol/src/lusol.o ogpf.o sparse_matrices.o helpers.o tests/test_sparse_matrices.f90 -fopenmp
test_lusol: lusol/src/lusol.o lusol/src/lusol_precision.o precision.o helpers.o diff.o sparse_matrices.o ogpf.o
	gfortran -O3 -floop-parallelize-all -I./lusol/src -o bin/test_lusol helpers.o precision.o diff.o lusol/src/lusol.o sparse_matrices.o ogpf.o tests/test_lusol.f90 -fopenmp
sparse_matrices.o: helpers.o precision.o lusol/src/lusol.o sparse_matrices.f90
	gfortran -O3 -floop-parallelize-all -I./lusol/src -c helpers.o precision.o lusol/src/lusol.o sparse_matrices.f90 -fopenmp
test_helpers: helpers.o ogpf.o precision.o tests/test_helpers.f90
	gfortran -O3 -floop-parallelize-all -o bin/test_helpers helpers.o precision.o ogpf.o tests/test_helpers.f90 -fopenmp
evolvePDE.o: config_m.o helpers.o ogpf.o diff.o precision.o lusol/src/lusol.o sparse_matrices.o evolvePDE.f90
	gfortran -O3 -floop-parallelize-all -I./lusol/src -c config_m.o precision.o helpers.o lusol/src/lusol.o ogpf.o diff.o sparse_matrices.o evolvePDE.f90
config_m.o: helpers.o precision.o config_m.f90
	gfortran -O3 -floop-parallelize-all -c precision.o helpers.o config_m.f90
ogpf.o:
	gfortran -O3 -floop-parallelize-all -c ogpf.f90
diff.o: sparse_matrices.o precision.o lusol/src/lusol.o diff.f90
	gfortran -O3 -floop-parallelize-all -I./lusol/src -c precision.o lusol/src/lusol.o sparse_matrices.o diff.f90
helpers.o: precision.o helpers.f90
	gfortran -O3 -floop-parallelize-all -c precision.o helpers.f90 -fopenmp
precision.o: precision.f90
	gfortran -O3 -floop-parallelize-all -c precision.f90 -fopenmp
clean:
	rm -f *.o *.mod

