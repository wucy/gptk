/*
 * demo_heterogeneous_noise.h
 *
 *  Created on: Apr 27, 2010
 *      Author: barillrl
 */

#ifndef DEMO_HETEROGENEOUS_NOISE_H_
#define DEMO_HETEROGENEOUS_NOISE_H_

#include <iostream>

#include "itpp/itbase.h"
#include "itpp/itstat.h"

#include "optimisation/SCGModelTrainer.h"

#include "io/csvstream.h"

#include "gaussian_processes/PSGP.h"
#include "gaussian_processes/GaussianProcess.h"

#include "likelihood_models/GaussianLikelihood.h"
#include "likelihood_models/GaussianSampLikelihood.h"
#include "likelihood_models/ExponentialSampLikelihood.h"

#include "covariance_functions/SumCF.h"
#include "covariance_functions/Matern3CF.h"
#include "covariance_functions/GaussianCF.h"
#include "covariance_functions/WhiteNoiseCF.h"

#include "plotting/GraphPlotter.h"

using namespace std;
using namespace itpp;


#endif /* DEMO_HETEROGENEOUS_NOISE_H_ */
