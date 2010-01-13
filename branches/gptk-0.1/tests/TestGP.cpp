#include "TestGP.h"

TestGP::TestGP()
{
	addTest(&testNoisySineLearnParams, "SSGP prediction with parameter estimation");
}

TestGP::~TestGP()
{
}

/**
* Load noisy sine data
**/
void TestGP::loadNoisySineData(vec &Xtrn, vec &Ytrn, vec &Xtst, vec &Ytst, 
		                                  vec &gpmean, vec &gpvar) 
{
	mat M;
	csvstream csv;

	csv.read(M, "noisy_sine_tst.csv");            	// Test set
	Xtst = M.get_col(0);
	Ytst = M.get_col(1);

	csv.read(M, "noisy_sine_trn.csv");			  	// Training set
	Xtrn = M.get_col(0);		
	Ytrn = M.get_col(1);

	csv.read(M, "noisy_sine_gp.csv");				// GP prediction
	gpmean = M.get_col(0);
	gpvar  = M.get_col(1);

}


/**
* Plot the mean and variance of SSGP and GP, along with
* training data and noiseless function.
**/
void TestGP::plotResults(vec gpmean, vec gpvar, vec truegpmean, vec truegpvar, 
		                           vec Xtrn, vec Ytrn, vec Xtst, vec Ytst)
{
	GraphPlotter gplot = GraphPlotter();
	gplot.clearPlot();

	// Plot ssgp
	gplot.plotPoints(Xtst, gpmean, "ssgp mean", LINE, BLUE);
	gplot.plotPoints(Xtst, gpmean + 2.0*sqrt(gpvar), "error bar", LINE, CYAN);
	gplot.plotPoints(Xtst, gpmean - 2.0*sqrt(gpvar), "error bar", LINE, CYAN);

	// Plot gp
	gplot.plotPoints(Xtst, truegpmean, "gp mean", LINE, GREEN);
	gplot.plotPoints(Xtst, truegpmean + 2.0*sqrt(truegpvar), "error bar", LINE, GREEN);
	gplot.plotPoints(Xtst, truegpmean - 2.0*sqrt(truegpvar), "error bar", LINE, GREEN);

	// Plot true function and observations
	gplot.plotPoints(Xtrn, Ytrn, "training set", CROSS, RED);
	gplot.plotPoints(Xtst, Ytst, "true function", LINE, RED);
	
}


/**
 * Test learning parameters with a GP
 **/
bool TestGP::testNoisySineLearnParams()
{
	vec Xtrn, Xtst, Ytrn, Ytst;
	vec truegpmean, truegpvar, gpmean, gpvar;
	double range  = 1.0;                     
	double sill   = 1.0;
	double nugget = 0.1;

	// Covariance function used for prediction
	GaussianCF predCovFunc(range, sill);

	// Load noisy sine data
	loadNoisySineData(Xtrn, Ytrn, Xtst, Ytst, truegpmean, truegpvar);

	// Covariance function: Gaussian + Nugget
	GaussianCF 	 g1(range, sill);            
	WhiteNoiseCF g2(nugget);
	SumCovarianceFunction gCmp(g1);
	gCmp.addCovarianceFunction(g2);

	// Initialise sequential GP
	mat Xtrnmat = Xtrn;
	GaussianProcess gp(1, 1, Xtrnmat, Ytrn, gCmp);

	// Learn parameters
	SCGModelTrainer gpTrainer(gp);

	// Make sure that checkGradient always returns the same answer
	cout << endl << endl << "*** Check that gradient is not changed by successive calls *** " << endl; 
	for (int i=0; i<0; i++)
		gpTrainer.checkGradient();

	// Check the gradient while moving away from the stored
	// parameters
	cout << endl << endl << "*** Check that gradient remains OK when moving away from stored parameters *** " << endl;

	vec x = gp.getParametersVector();
	vec eps = 0.1*ones(x.length());

	for (int i=0; i<10; i++) {
		cout << "X = " << endl;
		gpTrainer.checkGradient();
		gp.setParametersVector(x+i*eps);
	}

	// Learn parameters
	gpTrainer.Train(100);
	
	gCmp.displayCovarianceParameters();
	
	// Make prediction
	gp.makePredictions(gpmean, gpvar, Xtst, g1);
	cout << gpvar << endl;
	plotResults(gpmean, gpvar, truegpmean, truegpvar, Xtrn, Ytrn, Xtst, Ytst);
	
	return true;
}

/**
 * Main function
 * Simply creates an instance of this test suite and runs it
 */
int main() {
	TestGP test;
	test.run();
}
