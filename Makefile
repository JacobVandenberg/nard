.DEFAULT: nard
.PHONY:nard
nard:
	cd src && make nard && cd ..
tests:
	cd src && make tests && cd ..
.PHONY: clean
clean:
	cd src && make clean && cd ..
.PHONY: realclean
realclean:
	cd src && make realclean && cd ..