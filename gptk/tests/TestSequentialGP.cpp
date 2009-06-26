#include "TestSequentialGP.h"

TestSequentialGP::TestSequentialGP()
{
	header = "Test for class SequentialGP";
	
	addTest(&testNoisySine, "SSGP prediction with full training set as active points \n(should give same prediction as standard GP)");
}

TestSequentialGP::~TestSequentialGP()
{
}

bool TestSequentialGP::testNoisySine() 
{
	mat M;
	vec Xtrn, Xtst, Ytrn, Ytst;
	mat Xtrnmat;
	vec gpmean, gpvar, ssgpmean, ssgpvar;
	

	//---------------------------------------------------------------
	// Load data
	//---------------------------------------------------------------
	
	csvstream csv;
	
	csv.read(M, "noisy_sine_tst.csv");            	// Test set
	Xtst = M.get_col(0);
	Ytst = M.get_col(1);
	
	csv.read(M, "noisy_sine_trn.csv");			  	// Training set
	Xtrn = M.get_col(0);
	Xtrnmat = Xtrn;
	Ytrn = M.get_col(1);
	
	csv.read(M, "noisy_sine_gp.csv");				// GP prediction
	gpmean = M.get_col(0);
	gpvar  = M.get_col(1);
	
	
	//---------------------------------------------------------------
	// Set up SSGP
	//---------------------------------------------------------------
	
	// Fixed hyperparameters
	double range = sqrt(3.2);                     
	double sill  = 1.4;
	double nugget = 0.01;

	// Covariance function: Gaussian + Nugget
	GaussianCF 	 g1(range, sill);            
	WhiteNoiseCF g2(nugget);
	SumCovarianceFunction gCmp(g1);
	gCmp.addCovarianceFunction(g2);
	
	// Initialise sequential GP
	SequentialGP ssgp(1, 1, Ytrn.size(), Xtrnmat, Ytrn, gCmp);        
	
	// Gaussian observation likelihood
	GaussianLikelihood gaussLik(nugget);
	
	
	//---------------------------------------------------------------
	// Predict using SSGP
	//---------------------------------------------------------------
	ssgp.computePosterior(gaussLik);
	ssgp.makePredictions(ssgpmean, ssgpvar, Xtst);
	
	//---------------------------------------------------------------
	// TODO: Compare with GP
	//---------------------------------------------------------------
		

	//---------------------------------------------------------------
	// Plot
	//---------------------------------------------------------------
		
	GraphPlotter gplot = GraphPlotter();
	gplot.clearPlot();
	
	// Plot ssgp
	gplot.plotPoints(Xtst, ssgpmean, "ssgp mean", LINE, BLUE);
	gplot.plotPoints(Xtst, ssgpmean + 2.0*sqrt(ssgpvar), "error bar", LINE, CYAN);
	gplot.plotPoints(Xtst, ssgpmean - 2.0*sqrt(ssgpvar), "error bar", LINE, CYAN);
	
	// Plot gp
	gplot.plotPoints(Xtst, gpmean, "gp mean", LINE, GREEN);
	gplot.plotPoints(Xtst, gpmean + 2.0*sqrt(gpvar), "error bar", LINE, GREEN);
	gplot.plotPoints(Xtst, gpmean - 2.0*sqrt(gpvar), "error bar", LINE, GREEN);

	// Plot true function and observations
	gplot.plotPoints(Xtrn, Ytrn, "training set", CROSS, RED);
	gplot.plotPoints(Xtst, Ytst, "true function", LINE, RED);
	
	
	return false;	
}


int main() {
	TestSequentialGP test;
	test.run();
}
