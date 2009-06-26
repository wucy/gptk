ROOT = ..
include $(ROOT)/Makefile.common

# Output library
LIBRARY = libgptk
LIBRARY_PATH = .
VERSION = 0.0.1

LIBS = $(LIBS_COMMON)

# All general purpose files
OBJS =		    covarianceFunctions/GaussianCF.o \
				covarianceFunctions/ExponentialCF.o \
				covarianceFunctions/MaternCF.o \
				covarianceFunctions/ConstantCF.o \
				covarianceFunctions/LinearCF.o \
				covarianceFunctions/CovarianceFunction.o \
				covarianceFunctions/WhiteNoiseCF.o \
				covarianceFunctions/SumCovarianceFunction.o \
				covarianceFunctions/Transform.o \
				covarianceFunctions/IdentityTransform.o \
				covarianceFunctions/LogTransform.o \
				covarianceFunctions/NegLogSigmoidTransform.o \
				covarianceFunctions/SigmoidTransform.o \
				likelihoodModels/LikelihoodType.o \
				likelihoodModels/AnalyticLikelihood.o \
				likelihoodModels/GaussianLikelihood.o \
				likelihoodModels/SamplingLikelihood.o \
				likelihoodModels/GaussianSampLikelihood.o \
				likelihoodModels/SamplingLikelihoodMathML.o \
				likelihoodModels/GaussianSampLikelihoodMathML.o \
				likelihoodModels/MathMLParser/ExpressionNode.o \
				likelihoodModels/MathMLParser/MathMLParser.o \
				optimisation/Optimisable.o \
				optimisation/ModelTrainer.o \
				optimisation/SCGModelTrainer.o \
				optimisation/GDModelTrainer.o \
				optimisation/QuasiNewtonModelTrainer.o \
				optimisation/CGModelTrainer.o \
				GraphPlotter/GraphPlotter.o \
				gaussianProcesses/GaussianProcess.o \
				gaussianProcesses/SequentialGP.o \
				gaussianProcesses/ForwardModel.o \
				io/csvstream.o

.cpp.o: 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@	

all: libgptk

libgptk: $(OBJS)	
# 	Compile shared library 
	gcc -shared -Wl,-soname,$(LIBRARY).so -o $(LIBRARY_PATH)/$(LIBRARY).so.$(VERSION) $(OBJS) $(LIBS)
#   Remove symbolic link if existing
	if [[ -f $(LIBRARY_PATH)/$(LIBRARY).so ]]; then rm $(LIBRARY_PATH)/$(LIBRARY).so; fi
#   Create symbolic link to library
	ln -s $(LIBRARY_PATH)/$(LIBRARY).so.$(VERSION) $(LIBRARY_PATH)/$(LIBRARY).so
	@echo "** Library successfully built."

clean:
	rm -f $(OBJS)
	rm  $(LIBRARY_PATH)/$(LIBRARY).so.$(VERSION) $(LIBRARY_PATH)/$(LIBRARY).so
	@echo "** Cleaning done."
