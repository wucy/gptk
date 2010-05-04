#ifndef MAIN_H_
#define MAIN_H_

#include <iostream>

#include <itpp/itbase.h>

#include "itppext/itppext.h"
#include "io/csvstream.h"
#include "optimisation/SCGModelTrainer.h"

#include "gaussian_processes/GaussianProcess.h"
#include "gaussian_processes/PSGP.h"

#include "likelihood_models/GaussianLikelihood.h"

#include "covariance_functions/SumCF.h"
#include "covariance_functions/GaussianCF.h"
#include "covariance_functions/ExponentialCF.h"
#include "covariance_functions/Matern3CF.h"
#include "covariance_functions/Matern5CF.h"
#include "covariance_functions/NeuralNetCF.h"
#include "covariance_functions/ConstantCF.h"
#include "covariance_functions/WhiteNoiseCF.h"

#include "design/MaxMinDesign.h"
using namespace std;
using namespace itpp;

/** 
 * Main routine - runs the experiment for a given number of active points
 */
void run(int n_active);

#endif /*MAIN_H_*/
