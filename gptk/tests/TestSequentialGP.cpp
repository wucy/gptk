#include "TestSequentialGP.h"


TestSequentialGP::TestSequentialGP()
{
	header = "Test for class SequentialGP";
	
	// addTest(&testNoisySineFixedParams, "SSGP prediction with full training set as active points \n(should give same prediction as standard GP)");
	addTest(&testNoisySineLearnParams, "SSGP prediction with parameter estimation");
	// addTest(&testCheckGradient, "Gradient checks");
}


TestSequentialGP::~TestSequentialGP()
{
}


/** 
 * Load noisy sine data
 **/
void TestSequentialGP::loadNoisySineData(vec &Xtrn, vec &Ytrn, vec &Xtst, vec &Ytst, 
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
void TestSequentialGP::plotResults(vec ssgpmean, vec ssgpvar, vec gpmean, vec gpvar, 
		                           vec Xtrn, vec Ytrn, vec Xtst, vec Ytst, SequentialGP ssgp)
{
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
	
	// Plot active points
	vec activeX = (ssgp.getActiveSetLocations()).get_col(0);
	vec activeY = Ytrn(ssgp.getActiveSetIndices());
	gplot.plotPoints(activeX, activeY, "active points", CIRCLE, BLUE);
	cout << ssgp.getActiveSetIndices() << endl;
}


/**
 * Compares SSGP with standard GP on noisy sine data. Parameters
 * (length scale, process variance and noise variance) are set to
 * the true values.
 **/ 
bool TestSequentialGP::testNoisySineFixedParams() 
{
	vec Xtrn, Xtst, Ytrn, Ytst;
	vec gpmean, gpvar, ssgpmean, ssgpvar;
	int n_active = 7;
	double range  = 3.2;                     
	double sill   = 1.4;
	double nugget = 0.01;

	// Covariance function used for prediction
	GaussianCF predCovFunc(range, sill);
		
	// Load noisy sine data
	loadNoisySineData(Xtrn, Ytrn, Xtst, Ytst, gpmean, gpvar);
		
	// Covariance function: Gaussian + Nugget
	GaussianCF 	 g1(range, sill);            
	WhiteNoiseCF g2(nugget);
	SumCovarianceFunction gCmp(g1);
	gCmp.addCovarianceFunction(g2);

	// Initialise sequential GP
	mat Xtrnmat = Xtrn;
	SequentialGP ssgp(1, 1, n_active, Xtrnmat, Ytrn, gCmp);
		
	// Gaussian observation likelihood
	GaussianLikelihood gaussLik(nugget);
		
	//---------------------------------------------------------------
	// Predict using SSGP
	//---------------------------------------------------------------

	ssgp.computePosterior(gaussLik);

	ssgp.makePredictions(ssgpmean, ssgpvar, Xtst, predCovFunc);
	
	plotResults(ssgpmean, ssgpvar, gpmean, gpvar, Xtrn, Ytrn, Xtst, Ytst, ssgp);	
	
	return false;	
}


/**
 * Compares SSGP with standard GP on noisy sine data. Parameters
 * (length scale, process variance and noise variance) are estimated
 * from the data.
 **/ 
bool TestSequentialGP::testNoisySineLearnParams() 
{
	vec Xtrn, Xtst, Ytrn, Ytst;
	vec gpmean, gpvar, ssgpmean, ssgpvar;
	int n_active;
	double range  = 4.2;                     
	double sill   = 2.1;
	double nugget = 0.01;

	// Covariance function used for prediction
	GaussianCF predCovFunc(range, sill);
		
	// Load noisy sine data
	loadNoisySineData(Xtrn, Ytrn, Xtst, Ytst, gpmean, gpvar);
	n_active = 20; // Ytrn.size();
		
	// Covariance function: Gaussian + Nugget
	GaussianCF 	 g1(range, sill);            
	WhiteNoiseCF g2(nugget);
	SumCovarianceFunction gCmp(g1);
	gCmp.addCovarianceFunction(g2);

	// Initialise sequential GP
	mat Xtrnmat = Xtrn;
	SequentialGP ssgp(1, 1, n_active, Xtrnmat, Ytrn, gCmp);
		
	// Gaussian observation likelihood
	GaussianLikelihood gaussLik(nugget);
	ivec iActive = to_ivec(floor(linspace(0,Xtrnmat.rows()-1,n_active))); // maxDistance(Xtrnmat,n_active);
	
	ssgp.computePosteriorFixedActiveSet(gaussLik, iActive);
		
	// Learn parameters
	SCGModelTrainer gpTrainer(ssgp);
	ssgp.setLikelihoodType(UpperBound);

	gpTrainer.setCheckGradient(true);

	
	bvec optMask(3);
	optMask(0) = true;
	optMask(1) = true;
	optMask(2) = false;
	
	for (int i=0; i<5; i++) {
		gpTrainer.setOptimisationMask(optMask);
		gpTrainer.Train(100);
		gpTrainer.checkGradient();
		ssgp.resetPosterior();
		ssgp.computePosteriorFixedActiveSet(gaussLik,iActive);
	}
	
	gCmp.displayCovarianceParameters();
	
	// Predict using SSGP
	// ssgp.recomputePosterior();

	ssgp.makePredictions(ssgpmean, ssgpvar, Xtst, g1);
	
	plotResults(ssgpmean, ssgpvar, gpmean, gpvar, Xtrn, Ytrn, Xtst, Ytst, ssgp);	
	
	return false;	
}


/**
 * Performs several gradient checks
 * This test was used to identify an issue with the checkGradient() 
 * method in ModelTrainer.
 **/
bool TestSequentialGP::testCheckGradient()
{
	vec Xtrn, Xtst, Ytrn, Ytst;
	vec gpmean, gpvar, ssgpmean, ssgpvar;
	int n_active = 7;
	double range  = 3.2;                     
	double sill   = 1.0;
	double nugget = 0.1;

	// Covariance function used for prediction
	GaussianCF predCovFunc(range, sill);
		
	// Load noisy sine data
	loadNoisySineData(Xtrn, Ytrn, Xtst, Ytst, gpmean, gpvar);
		
	// Covariance function: Gaussian + Nugget
	GaussianCF 	 g1(range, sill);            
	WhiteNoiseCF g2(nugget);
	SumCovarianceFunction gCmp(g1);
	gCmp.addCovarianceFunction(g2);

	// Initialise sequential GP
	mat Xtrnmat = Xtrn;
	SequentialGP ssgp(1, 1, n_active, Xtrnmat, Ytrn, gCmp);
		
	// Gaussian observation likelihood
	GaussianLikelihood gaussLik(nugget);
	ssgp.computePosterior(gaussLik);
		
	// Learn parameters
	SCGModelTrainer gpTrainer(ssgp);
	ssgp.setLikelihoodType(UpperBound);
	
	// Make sure that checkGradient always returns the same answer
	cout << endl << endl << "*** Check that gradient is not changed by successive calls *** " << endl; 
	for (int i=0; i<10; i++)
		gpTrainer.checkGradient();

	// Check the gradient while moving away from the stored
	// parameters
	cout << endl << endl << "*** Check that gradient remains OK when moving away from stored parameters *** " << endl;
	for (int i=0; i<10; i++) {
		vec x = ssgp.getParametersVector();
		cout << "X = " << endl;
		gpTrainer.checkGradient();
		ssgp.setParametersVector(x+0.1*ones(x.length()));
	}

	
	return true;
}



/** THE FOLLOWING TWO METHODS ARE NOT VERY USEFUL - 
 *  I WAS TRYING TO FIND A DESIGN THAT MAXIMISES THE DISTANCE
 *  BETWEEN POINTS, BUT THIS GIVES A POOR DESIGN. NOT WHAT I
 *  MEAN TO DO. LEFT HERE FOR FURTHER IMPROVEMENT.
 */

/**
 * Set of points with maximum total distance (edges)
 */
ivec TestSequentialGP::maxDistance(mat x, int n)
{
	assert(n <= x.rows());
	
	mat xbest(n,x.cols());
	ivec best(n);
	vec  vrand = randu(x.rows());
	ivec irand = sort_index(vrand);
 
		
	// Initialise active set
	for (int i=0; i<n; i++) {
		best(i) = irand(i);
		xbest.set_row(i,x.get_row(irand(i)));
	}
	cout << xbest << endl;
	cout << best << endl;
	
	// Look for points that increase the distance
	
	for (int i=n; i<x.rows(); i++) {
		vec obs = x.get_row(i);
		int ibetter = maxTotalDistance(xbest, obs);
		
		if (ibetter >= 0) {   // Replace point ibetter with point i
			cout << ibetter << endl;
			best(ibetter) = i;
			xbest.set_row(ibetter, x.get_row(i));
		}
		cout << best << endl;
	}
	
	return best;
}

/** 
 * Return a vector of total edge distances for each set x where
 * the i-th point of x has been replaced with x*
 */
int TestSequentialGP::maxTotalDistance(mat x, vec xstar)
{
	int N = x.rows();
	vec D = zeros(N), Dtotal = zeros(N);
	vec Dstar = zeros(N);
	vec d, Dnew = zeros(N);
	
	// Compute distance matrix
	for (int i=0; i<N; i++) { 
		D(i) = 0.0;
		
		for (int j=0; j<N; j++) { 
			d = x.get_row(i) - x.get_row(j);
			D(i) += sum(abs(d));      // Total distance between x_i and the other x_j's
		}
		
		d = x.get_row(i) - xstar; 
		Dstar(i) = sum(abs(d));       // Distance between x_i and x*
	}
	
	double dstartot = sum(Dstar);
	
	// Compute total distance for each set in which x_i has been replaced
	// by x*
	for (int i=0; i<N; i++) {
		Dnew(i) = dstartot - Dstar(i) - D(i);
	}
		
	// Find index of point maximising the new total distance
	int imax = max_index(Dnew);
	
	// Make sure the new distance is greater than the previous one 
	if (Dnew(imax)>0) 
		return imax;
	else 
		return -1;
	
}




/**
 * Main function
 * Simply creates an instance of this test suite and runs it
 */
int main() {
	TestSequentialGP test;
	test.run();
}