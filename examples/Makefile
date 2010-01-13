ROOT = ..

include $(ROOT)/Makefile.common

OBJS = Rosenbrock.o RosenbrockND.o
EXE = demOpt1 demOpt2 demGP1 demGP2 demSSGP1 demSSGP2
LIBS = $(LIBS_COMMON) -lgptk -L$(ROOT)/libgptk

all:	demOpt1 demOpt2 demGP1 demGP2 demSSGP1 demSSGP2


# Compile sources into objects
.cpp.o: 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@	

# Optimisation Demonstration 1
demOpt1: demOpt1.o $(OBJS) 
	$(CXX) -o demOpt1 demOpt1.o  $(OBJS) $(LIBS)

# Optimisation Demonstration 2
demOpt2:	demOpt2.o  $(OBJS)  
	$(CXX) -o demOpt2 demOpt2.o  $(OBJS) $(LIBS)

# Gaussian Process Demonstration 1
demGP1:	demGP1.o
	$(CXX) -o demGP1 demGP1.o $(LIBS)

# Gaussian Process Demonstration 2
demGP2:	demGP2.o
	$(CXX) -o demGP2 demGP2.o $(LIBS)

# Sequentiall Gaussian Process Demonstration 1
demSSGP1:	demSSGP1.o 
	$(CXX) -o demSSGP1 demSSGP1.o $(LIBS)

# Sequentiall Gaussian Process Demonstration 2
demSSGP2:	demSSGP2.o
	$(CXX) -o demSSGP2 demSSGP2.o $(LIBS)

# Sequentiall Gaussian Process Demonstration 2
demSSGPAlgo:	demSSGPAlgo.o
	$(CXX) -o demSSGPAlgo demSSGPAlgo.o $(LIBS)

# Sequentiall Gaussian Process - Helicopter data set
demSSGPHelicopter:	demSSGPHelicopter.o
	$(CXX) -o demSSGPHelicopter demSSGPHelicopter.o $(LIBS)


clean:
	rm -f *.o $(EXE)
