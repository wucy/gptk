#ifndef MAIN_H_
#define MAIN_H_

#include <iostream>

#include "itpp/itbase.h"
#include "itpp/itstat.h"

#include "itppext/itppext.h"
#include "io/csvstream.h"
#include "optimisation/SCGModelTrainer.h"

#include "gaussianProcesses/GaussianProcess.h"
#include "gaussianProcesses/PSGP.h"

#include "likelihoodModels/GaussianLikelihood.h"

#include "covarianceFunctions/SumCovarianceFunction.h"
#include "covarianceFunctions/GaussianCF.h"
#include "covarianceFunctions/ExponentialCF.h"
#include "covarianceFunctions/Exponential2CF.h"
#include "covarianceFunctions/Matern3CF.h"
#include "covarianceFunctions/Matern5CF.h"
#include "covarianceFunctions/NeuralNetCF.h"
#include "covarianceFunctions/ConstantCF.h"
#include "covarianceFunctions/WhiteNoiseCF.h"

using namespace std;
using namespace itpp;

/** 
 * Main routine - runs the experiment for a given number of active points
 */
void run(int n_active);

#endif /*MAIN_H_*/
