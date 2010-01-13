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
#include "covarianceFunctions/ConstantCF.h"

#include "optimisation/SCGModelTrainer.h"
#include "optimisation/CGModelTrainer.h"
#include "optimisation/QuasiNewtonModelTrainer.h"

#include "design/MaxMinDesign.h"

#include "likelihoodModels/GaussianLikelihood.h"

#include "gaussianProcesses/SequentialGP.h"

#include "GraphPlotter/GraphPlotter.h"



class TestSequentialGP : public Test
{
public:
	TestSequentialGP();
	virtual ~TestSequentialGP();
	
	static bool testNoisySineFixedParams();
	static bool testNoisySineLearnParams();
	static bool testCheckGradient();
	
private:
	static void loadNoisySineData(vec &Xtrn, vec &Ytrn, vec &Xtst, vec &Ytst, vec &gpmean, vec &gpvar);
	static void plotResults(vec ssgpmean, vec ssgpvar, vec gpmean, vec gpvar, 
                            vec Xtrn, vec Ytrn, vec Xtst, vec Ytst, SequentialGP ssgp);
	static void plotOptLog(mat theta);
	
	static ivec  maxDistance(mat x, int n);
	static int  maxTotalDistance(mat x, vec xstar);
};

#endif /*TESTSEQUENTIALGP_H_*/
