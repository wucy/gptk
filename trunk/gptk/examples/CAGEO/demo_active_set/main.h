#ifndef MAIN_H_
#define MAIN_H_

#include <iostream>

#include "itpp/itbase.h"
#include "itpp/itstat.h"

#include "io/csvstream.h"
#include "optimisation/SCGModelTrainer.h"

#include "gaussianProcesses/GaussianProcess.h"
#include "gaussianProcesses/PSGP.h"

#include "likelihoodModels/GaussianLikelihood.h"

#include "covarianceFunctions/SumCovarianceFunction.h"
#include "covarianceFunctions/GaussianCF.h"
#include "covarianceFunctions/WhiteNoiseCF.h"

#include "GraphPlotter/GraphPlotter.h"

using namespace std;
using namespace itpp;

mat computeLikelihoodProfile(PSGP &psgp, mat paramRanges);
mat computeLikelihoodProfileGP(GaussianProcess &psgp, mat paramRanges);


#endif /*MAIN_H_*/
