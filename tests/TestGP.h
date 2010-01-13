#ifndef TESTGP_H_
#define TESTGP_H_


#include <iostream>

#include "Test.h"
#include "itpp/itbase.h"
#include "itpp/itstat.h"
#include "io/csvstream.h"

#include "covarianceFunctions/SumCovarianceFunction.h"
#include "covarianceFunctions/GaussianCF.h"
#include "covarianceFunctions/WhiteNoiseCF.h"

#include "optimisation/SCGModelTrainer.h"

#include "likelihoodModels/GaussianLikelihood.h"

#include "gaussianProcesses/GaussianProcess.h"

#include "GraphPlotter/GraphPlotter.h"

/**
 * Test class for standard Gaussian Process
 **/
class TestGP : public Test
{
public:
	TestGP();
	virtual ~TestGP();
	
	static bool testNoisySineLearnParams();
	
private:
	static void loadNoisySineData(vec &Xtrn, vec &Ytrn, vec &Xtst, vec &Ytst, vec &gpmean, vec &gpvar);
	static void plotResults(vec ssgpmean, vec ssgpvar, vec gpmean, vec gpvar, 
	                            vec Xtrn, vec Ytrn, vec Xtst, vec Ytst);
};

#endif /*TESTGP_H_*/
