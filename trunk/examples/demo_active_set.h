#ifndef DEMO_ACTIVE_SET_H_
#define DEMO_ACTIVE_SET_H_

#include <iostream>

#include "itpp/itbase.h"
#include "itpp/itstat.h"

#include "io/csvstream.h"
#include "optimisation/SCGModelTrainer.h"

#include "gaussian_processes/GaussianProcess.h"
#include "gaussian_processes/PSGP.h"

#include "likelihood_models/GaussianLikelihood.h"

#include "covariance_functions/SumCF.h"
#include "covariance_functions/GaussianCF.h"
#include "covariance_functions/WhiteNoiseCF.h"

#include "plotting/GraphPlotter.h"

using namespace std;
using namespace itpp;

#endif
