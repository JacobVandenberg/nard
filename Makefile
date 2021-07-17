.DEFAULT: tests
tests:
	cd src && make tests && cd ..
.PHONY: clean
clean:
	cd src && make clean && cd ..
.PHONY: realclean
realclean:
	cd src && make realclean && cd ..