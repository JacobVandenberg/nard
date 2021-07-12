
code: test_diff test_helpers
	echo "Finished"
test_diff: diff.o tests/test_diff.f90
	gfortran -o bin/test_diff diff.o tests/test_diff.f90
test_helpers: helpers.o diff.o ogpf.o tests/test_helpers.f90
	gfortran -o bin/test_helpers helpers.o diff.o ogpf.o tests/test_helpers.f90
ogpf.o:
	gfortran -c ogpf.f90
diff.o: diff.f90
	gfortran -c diff.f90
helpers.o: helpers.f90
	gfortran -c helpers.f90
clean:
	rm -f *.o *.mod

