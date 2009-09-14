#ifndef ASTONGEOSTATS_H
#define ASTONGEOSTATS_H



#include "R.h"
#include "Rmath.h"

#include "itpp/itbase.h"
#include "itpp/itstat.h"

#include "gaussianProcesses/ForwardModel.h"
// #include "gaussianProcesses/SequentialGP.h"
#include "gaussianProcesses/PSGP.h"
#include "gaussianProcesses/GaussianProcess.h"
#include "optimisation/Optimisable.h"
#include "optimisation/SCGModelTrainer.h"
#include "optimisation/QuasiNewtonModelTrainer.h"
#include "optimisation/CGModelTrainer.h"
#include "optimisation/GDModelTrainer.h"
#include "optimisation/ModelTrainer.h"

#include "likelihoodModels/LikelihoodType.h"
#include "likelihoodModels/ExponentialSampLikelihood.h"
#include "likelihoodModels/GaussianLikelihood.h"

#include "covarianceFunctions/CovarianceFunction.h"
#include "covarianceFunctions/GaussianCF.h"
#include "covarianceFunctions/ExponentialCF.h"
#include "covarianceFunctions/WhiteNoiseCF.h"
#include "covarianceFunctions/SumCovarianceFunction.h"
#include "covarianceFunctions/LogTransform.h"
#include "covarianceFunctions/IdentityTransform.h"
#include "covarianceFunctions/NegLogSigmoidTransform.h"
#include "covarianceFunctions/ConstantCF.h"

#include "itppext/itppext.h"
#include "io/csvstream.h"

#include <cassert>
#include <string>
#include <algorithm>
#include <stdexcept>

// Maximum number of observations kept for param estimation
#define MAX_OBS 1000

// Maximum number of active points
#define MAX_ACTIVE_POINTS 400

// ID for invalid noise model (when using observation noise) 
#define INVALID_MODEL_NAME "INVALID_MODEL"

// Whether to use a GP instead of PSGP for parameter estimation
#define PARAMETER_ESTIMATION_USING_GP false

// Outer loops in parameter estimation for PSGP
#define PSGP_PARAM_ITERATIONS 5

// Inner loop (i.e. SCG iterations in each outer loop) for PSGP
#define PSGP_SCG_ITERATIONS 5



using namespace std;
using namespace itpp;

void learnParameters(int numObs, int numError, double *xData, double *yData, double *errorData, 
                     int numMetadata, int *errPtr, int *sensorPtr, char **metaDataTable, 
                     double *range, double *sill, double *nugget, double *bias, int model);

void makePredictions(int numObs, int numPred, int numError, double *xData, double *yData, 
                     double *errorData, double *xPred, int numMetadata, int *errPtr, int *sensorPtr, 
                     char **metaDataTable, double *meanPred, double *varPred, double *range, 
                     double *sill, double *nugget, double *bias, int model);

void parseMetadata(char** metadataTable, const ivec likelihoodModelIndexes,  
                   string modelNames[], vec modelParams[]);

void buildLikelihoodVector(string modelName[], vec modelParams[], ivec iLikelihoodModel, 
                           Vec<LikelihoodType*> &likelihoodModels);

int minusOne(int n);


#endif // ASTONGEOSTATS_H
