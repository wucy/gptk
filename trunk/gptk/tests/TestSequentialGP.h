#ifndef TESTSEQUENTIALGP_H_
#define TESTSEQUENTIALGP_H_

#include <iostream>

#include "Test.h"
#include "itpp/itbase.h"
#include "itpp/itstat.h"
#include "io/csvstream.h"

#include "covarianceFunctions/SumCovarianceFunction.h"
#include "covarianceFunctions/GaussianCF.h"
#include "covarianceFunctions/WhiteNoiseCF.h"

#include "likelihoodModels/GaussianLikelihood.h"

#include "gaussianProcesses/SequentialGP.h"

#include "GraphPlotter/GraphPlotter.h"



class TestSequentialGP : public Test
{
public:
	TestSequentialGP();
	virtual ~TestSequentialGP();
	
	static bool testNoisySine();
};

#endif /*TESTSEQUENTIALGP_H_*/
