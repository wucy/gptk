ROOT = ..

include $(ROOT)/Makefile.common

LIBS = $(LIBS_COMMON) -lgptk -L$(ROOT)/libgptk

EXE = testSequentialGP testCSVStream testGradientCovFunc 

all: testSequentialGP testCSVStream testGradientCovFunc


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

clean:
	rm -f *.o $(EXE)

