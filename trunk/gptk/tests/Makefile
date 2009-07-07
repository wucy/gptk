ROOT = ..

include $(ROOT)/Makefile.common

LIBS = $(LIBS_COMMON) -lgptk -L$(ROOT)/libgptk

EXE = testSequentialGP testCSVStream testGradientCovFunc 

all: testSequentialGP testCSVStream testGradientCovFunc testOptimisationSCG testGP \
     testSeparableGP testITPPext


# Compile sources into objects
.cpp.o: 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@	

# CSV input/output
testCSVStream: TestCSVStream.o Test.o
	$(CXX) -o testCSVStream TestCSVStream.o Test.o $(LIBS)
	@echo "** Successfully built testCSVStream."

# CSV input/output
testGradientCovFunc: TestGradientCovFunc.o Test.o
	$(CXX) -o testGradientCovFunc TestGradientCovFunc.o Test.o $(LIBS)
	@echo "** Successfully built testGradientCovFunc."

# Sequentiall Gaussian Process
testSequentialGP: TestSequentialGP.o Test.o
	$(CXX) -o testSequentialGP TestSequentialGP.o Test.o $(LIBS)
	@echo "** Successfully built testSequentialGP."

# SCG
testOptimisationSCG: TestOptimisationSCG.o Test.o
	$(CXX) -o testOptimisationSCG TestOptimisationSCG.o Test.o $(LIBS)
	@echo "** Successfully built testOptimisationSCG."

# Gaussian Process
testGP: TestGP.o Test.o
	$(CXX) -o testGP TestGP.o Test.o $(LIBS)
	@echo "** Successfully built testGP."

# Separable Gaussian Process
testSeparableGP: TestSeparableGP.o Test.o
	$(CXX) -o testSeparableGP TestSeparableGP.o Test.o $(LIBS)
	@echo "** Successfully built testSeparableGP."

# ITPPEXT
testITPPext: TestITPPext.o Test.o
	$(CXX) -o testITPPext TestITPPext.o Test.o $(LIBS)
	@echo "** Successfully built testITPPext."


clean:
	rm -f *.o $(EXE)

