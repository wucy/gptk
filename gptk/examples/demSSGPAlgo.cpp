/***************************************************************************
 *   AstonGeostats, algorithms for low-rank geostatistical models          *
 *                                                                         *
 *   Copyright (C) Ben Ingram, 2008-2009                                   *
 *                                                                         *
 *   Remi Barillec, r.barillec@aston.ac.uk                                 *
 *   Ben Ingram, IngramBR@Aston.ac.uk                                      *
 *   Neural Computing Research Group,                                      *
 *   Aston University,                                                     *
 *   Aston Street, Aston Triangle,                                         *
 *   Birmingham. B4 7ET.                                                   *
 *   United Kingdom                                                        *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/



#include <iostream>
#include <cassert>
#include <ctime>

#include "itpp/itbase.h"
#include "itpp/itstat.h"

#include "optimisation/Optimisable.h"
#include "optimisation/ModelTrainer.h"
#include "optimisation/SCGModelTrainer.h"
#include "optimisation/CGModelTrainer.h"
#include "optimisation/GDModelTrainer.h"
#include "optimisation/QuasiNewtonModelTrainer.h"

#include "ForwardModel.h"
#include "optimisation/Optimisable.h"
#include "SequentialGP.h"

#include "likelihoodModels/GaussianLikelihood.h"

#include "covarianceFunctions/CovarianceFunction.h"
#include "covarianceFunctions/SumCovarianceFunction.h"
#include "covarianceFunctions/GaussianCF.h"
#include "covarianceFunctions/MaternCF.h"
#include "covarianceFunctions/ExponentialCF.h"
#include "covarianceFunctions/WhiteNoiseCF.h"
#include "covarianceFunctions/ConstantCF.h"
#include "covarianceFunctions/LinearCF.h"

#include "GraphPlotter/GraphPlotter.h"

#include "covarianceFunctions/NegLogSigmoidTransform.h"
#include "covarianceFunctions/SigmoidTransform.h"

#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>


using namespace std;
using namespace itpp;

/**
 * This demonstration program illustrates the functioning of the 
 * SSGP algorithm. Noisy observations are sampled from a Gaussian Process
 * with known parameters. We then display the posterior distribution as 
 * active points are added.  
 **/

// Program constants
#define RANDOM_SEED 13 

#define N_OBS 50		  // Number of observations (total)
#define N_ACTIVE 5        // Number of active points


int main()
{
	RNG_reset(RANDOM_SEED);

	/*------------------------------------------------------------------------*
	 * Generate observations by sampling from a GP
	 *------------------------------------------------------------------------*/
	double lengthScale = 1.1;
	double variance = 1.0;
	double nugget = 0.5;
	
	// Input points sampled uniformly on [0,10]
	vec X = zeros(N_OBS);
	for (int i=0; i<N_OBS; i++) X(i) = 10.0*i/N_OBS;

	// Covariance function
	GaussianCF covf1(lengthScale, variance);
	WhiteNoiseCF covf2(nugget); 
	CovarianceFunction *covf = new SumCovarianceFunction(covf1);
	((SumCovarianceFunction*) covf)->addCovarianceFunction(covf2);
	
	covf->displayCovarianceParameters();
	
	mat K = zeros(N_OBS,N_OBS); 
	covf->computeCovariance(K, X, X);

	vec Y =  transpose(chol(K))*randn(N_OBS);
	
	// Plot observed data
	GraphPlotter gplot = GraphPlotter();
	gplot.clearPlot();
	gplot.plotPoints(X, Y, "data", CROSS, BLUE);
	
	/*------------------------------------------------------------------------*
	 * 
	 *------------------------------------------------------------------------*/
	
	
}
