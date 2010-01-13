all: 
	make _libgptk
	make _examples
	make _tests

_libgptk:
	@echo
	@echo "***************************************************************"
	@echo "Compiling Gaussian Process Toolkit library"
	@echo "***************************************************************"
	cd libgptk; make all

_tests:
	@echo
	@echo "***************************************************************"
	@echo "Compiling tests"
	@echo "***************************************************************"
	cd tests; make all


_examples:
	@echo
	@echo "***************************************************************"
	@echo "Compiling examples"
	@echo "***************************************************************"
	cd examples; make all

clean:
	cd libgptk;     make clean
	cd examples; make clean
	cd tests;    make clean