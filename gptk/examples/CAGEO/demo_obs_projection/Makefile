# Compiler options
CXX = g++ 
CXXFLAGS = `itpp-config --cflags` -O3 -Wall -fmessage-length=0 -g

# Root folder of PSGP installation
ROOT = ../../..

INCLUDES = -I$(ROOT)/lib/xercesc/include -I$(ROOT)/libgptk
LIBS = `itpp-config --libs` -L$(ROOT)/lib/xercesc/lib -L$(ROOT)/libgptk/


# Name of resulting program
PROGRAM_NAME = demo_obs_projection

all: $(PROGRAM_NAME)

$(PROGRAM_NAME): main.o
	$(CXX) $(CXXFLAGS) $(LIBS) main.o -lgptk -lxerces-c -o $(PROGRAM_NAME)

main.o: main.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c main.cpp -o main.o

clean: 
	rm main.o $(PROGRAM_NAME)
	