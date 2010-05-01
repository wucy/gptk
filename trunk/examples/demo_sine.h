/*
 * demo_sine.h
 *
 *  Created on: 1 May 2010
 *      Author: Remi Barillec <r.barillec@aston.ac.uk>
 */

#ifndef DEMO_SINE_H_
#define DEMO_SINE_H_

#include <iostream>
#include <itpp/itbase.h>

#include <itppext/itppext.h>
#include "plotting/GraphPlotter.h"
#include "covariance_functions/SumCF.h"
#include "covariance_functions/GaussianCF.h"
#include "covariance_functions/WhiteNoiseCF.h"
#include "gaussian_processes/PSGP.h"
#include "likelihood_models/GaussianLikelihood.h"
#include "optimisation/SCGModelTrainer.h"

using namespace std;
using namespace itpp;
using namespace itppext;

#endif /* DEMO_SINE_H_ */

