
code: test_diff test_helpers test_evolvePDE test_sparse_matrices
	echo "Finished"
test_evolvePDE: evolvePDE.o config_m.o helpers.o ogpf.o diff.o sparse_matrices.o tests/test_evolvePDE.f90
	gfortran -O3 -floop-parallelize-all -o bin/test_evolvePDE evolvePDE.o config_m.o helpers.o ogpf.o diff.o sparse_matrices.o tests/test_evolvePDE.f90 -fopenmp
test_diff: diff.o ogpf.o sparse_matrices.o helpers.o tests/test_diff.f90
	gfortran -O3 -floop-parallelize-all -o bin/test_diff diff.o ogpf.o sparse_matrices.o helpers.o tests/test_diff.f90 -fopenmp
test_sparse_matrices: sparse_matrices.o helpers.o tests/test_sparse_matrices.f90
	gfortran -O3 -floop-parallelize-all -o bin/test_sparse_matrices sparse_matrices.o helpers.o tests/test_sparse_matrices.f90 -fopenmp
sparse_matrices.o: helpers.o sparse_matrices.f90
	gfortran -O3 -floop-parallelize-all -c helpers.o sparse_matrices.f90 -fopenmp
test_helpers: helpers.o ogpf.o tests/test_helpers.f90
	gfortran -O3 -floop-parallelize-all -o bin/test_helpers helpers.o ogpf.o tests/test_helpers.f90 -fopenmp
evolvePDE.o: config_m.o helpers.o ogpf.o diff.o sparse_matrices.o evolvePDE.f90
	gfortran -O3 -floop-parallelize-all -c config_m.o helpers.o ogpf.o diff.o sparse_matrices.o evolvePDE.f90
config_m.o: helpers.o config_m.f90
	gfortran -O3 -floop-parallelize-all -c helpers.o config_m.f90
ogpf.o:
	gfortran -O3 -floop-parallelize-all -c ogpf.f90
diff.o: sparse_matrices.o diff.f90
	gfortran -O3 -floop-parallelize-all -c sparse_matrices.o diff.f90
helpers.o: helpers.f90
	gfortran -O3 -floop-parallelize-all -c helpers.f90 -fopenmp
clean:
	rm -f *.o *.mod

