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


void TestSequentialGP::plotOptLog(mat theta) 
{
    GraphPlotter gplot = GraphPlotter();
    gplot.clearPlot();

    vec iter = linspace(1,theta.rows(),theta.rows());
    
    // Plot ssgp
    gplot.plotPoints(iter, theta.get_col(0), "length scale", LINE, BLUE);
    gplot.plotPoints(iter, theta.get_col(1), "proc var", LINE, RED);
    gplot.plotPoints(iter, theta.get_col(2), "nugget", LINE, GREEN);
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
	int n_active;
	double range  = 3.2;                     
	double sill   = 1.4;
	double nugget = 0.01;

	// Covariance function used for prediction
	GaussianCF predCovFunc(range, sill);
		
	// Load noisy sine data
	loadNoisySineData(Xtrn, Ytrn, Xtst, Ytst, gpmean, gpvar);
	n_active = 10; // Ytrn.size();
	
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
		
	// Predict using SSGP
	ssgp.computePosterior(gaussLik);

	ssgp.makePredictions(ssgpmean, ssgpvar, Xtst, g1);
	
	plotResults(ssgpmean, ssgpvar, gpmean, gpvar, Xtrn, Ytrn, Xtst, Ytst, ssgp);	
	
	return false;	
}

/** -----------------------------------------------------------------------------------------------------
 * Compares SSGP with standard GP on noisy sine data. Parameters
 * (length scale, process variance and noise variance) are estimated
 * from the data.
 **/ 
bool TestSequentialGP::testNoisySineLearnParams() 
{
	vec Xtrn, Xtst, Ytrn, Ytst;
	vec gpmean, gpvar, ssgpmean, ssgpvar;
	int n_active;
	double range  = 3.2;                     
	double sill   = 1.4;
	double nugget = 0.01;

	// Covariance function used for prediction
	GaussianCF predCovFunc(range, sill);
		
	// Load noisy sine data
	loadNoisySineData(Xtrn, Ytrn, Xtst, Ytst, gpmean, gpvar);
	n_active =  10; // Ytrn.length();
		
	// Covariance function: Gaussian + Nugget
	GaussianCF 	 g1(range, sill);            
	WhiteNoiseCF g2(nugget);
	SumCovarianceFunction gCmp(g1);
	gCmp.addCovarianceFunction(g2);

	// Initialise sequential GP
	mat Xtrnmat = Xtrn;
	SequentialGP ssgp(1, 1, n_active, Xtrnmat, Ytrn, gCmp, 1);
		
	// Gaussian observation likelihood
	GaussianLikelihood gaussLik(nugget);

	// ivec iActive = to_ivec(floor(linspace(0,Xtrnmat.rows()-1,n_active)));
	MaxMinDesign design(100);
	ivec iActive = design.subsample(Xtrnmat, Ytrn, n_active);

	// ssgp.computePosteriorFixedActiveSet(gaussLik, iActive);
	ssgp.computePosterior(gaussLik);
		
	ssgp.makePredictions(ssgpmean, ssgpvar, Xtst, g1);
	    
	plotResults(ssgpmean, ssgpvar, gpmean, gpvar, Xtrn, Ytrn, Xtst, Ytst, ssgp);    
	
	return 0;
	
	// Learn parameters
	SCGModelTrainer gpTrainer(ssgp);
	ssgp.setLikelihoodType(UpperBound);

	/*
	// Compute likelihood profile for length scale and variance
	vec params = oldParams;
	mat likprof = zeros(100,100);
	cout << endl;
	for (int i=0; i<100; i++) {
	    params(0) = oldParams(0) - 1.0 + 0.02*(1+i);
	    for (int j=0; j<100; j++)
	    {
	        params(1) = oldParams(1) - 1.0 + 0.02*(1+j);
	        ssgp.setParametersVector(params);
	        likprof(i,j) = ssgp.objective();
	    }
	    cout << "\r" << i << flush;
	}        
	cout << endl;
	csvstream csv;
	csv.write(likprof, "likelihood_profile_params.csv");
	
	// Compute likelihood profile for length scale and variance
	    mat likprof2 = zeros(400,400);
	    cout << endl;
	    for (int i=0; i<400; i++) {
	        params(0) = -4.0 + 0.02*(1+i);
	        for (int j=0; j<400; j++)
	        {
	            params(1) = - 4.0 + 0.02*(1+j);
	            ssgp.setParametersVector(params);
	            likprof2(i,j) = ssgp.objective();
	        }
	        cout << "\r" << i << flush;
	    }        
	    cout << endl;
	    csv.write(likprof2, "likelihood_profile.csv");
	  */  
	
	// Plot before training
	// ssgp.makePredictions(ssgpmean, ssgpvar, Xtst, predCovFunc);
	// plotResults(ssgpmean, ssgpvar, gpmean, gpvar, Xtrn, Ytrn, Xtst, Ytst, ssgp);    
	gpTrainer.setCheckGradient(true);
	
	bvec optMask(3);
	
	int niter = 10;
	mat theta = zeros(niter+1,3);
	theta.set_row(0, ssgp.getParametersVector());

	//===============================================================
	// It is critical to limit the number optimisation steps
	// We don't want to reach full convergence, as the "optimal"
	// solution is conditionned on our current estimate of the 
	// C's and alpha's, and thus only optimal with respect to these.
	// We need to optimise a little, then reestimate the C's and 
	// alpha's, and repeat.
	//===============================================================
	for (int i=0; i<niter; i++) 
	{
	    // Optimise covariance parameters 
	    optMask(0) = true;
	    optMask(1) = true;
	    optMask(2) = false;
	    gpTrainer.setOptimisationMask(optMask);
	    
	    gpTrainer.Train(5);
	    gpTrainer.checkGradient();
	    
	    // Recompute basis vectors 
	    ssgp.resetPosterior();
	    // ssgp.computePosterior(gaussLik);
	    ssgp.computePosteriorFixedActiveSet(gaussLik,iActive);
	    // ssgp.recomputePosteriorFixedActiveSet(gaussLik);
	    // ssgp.recomputePosterior();
	    	    
	    
	        // Optimise noise parameters
	    /*
	    optMask(0) = false;
	    optMask(1) = false;
	    optMask(2) = true;
	    gpTrainer.setOptimisationMask(optMask);
	    gpTrainer.Train(5);
	    gpTrainer.checkGradient();
	    // g2.displayCovarianceParameters();

	    ssgp.resetPosterior();
	    ssgp.computePosterior(gaussLik);
	    // ssgp.computePosteriorFixedActiveSet(gaussLik,iActive);
	    */
	    
	    theta.set_row(i+1, ssgp.getParametersVector());
	}
	
	gCmp.displayCovarianceParameters();
	
	// Recompute posterior using estimated parameters
	// and smaller active set
	// ssgp.setActiveSetSize(10);
	// ssgp.resetPosterior();
	// ivec iActive = to_ivec(floor(linspace(0,Xtrnmat.rows()-1,10)));
	// ssgp.computePosteriorFixedActiveSet(gaussLik, iActive);
	// ssgp.computePosterior(gaussLik);
	
	
	
	// Recompute posterior using estimated parameters
	// and smaller active set
	// ssgp.setActiveSetSize(10);
	// ssgp.resetPosterior();
	// ivec iActive = to_ivec(floor(linspace(0,Xtrnmat.rows()-1,10)));
	// ssgp.computePosteriorFixedActiveSet(gaussLik, iActive);
	// ssgp.computePosterior(gaussLik);
	
	// Plot after training
    ssgp.makePredictions(ssgpmean, ssgpvar, Xtst, g1);
	
	plotResults(ssgpmean, ssgpvar, gpmean, gpvar, Xtrn, Ytrn, Xtst, Ytst, ssgp);
	plotOptLog(theta); 
	
	return false;	
}



/** -----------------------------------------------------------------------------------------------------
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
	for (int i=0; i<0; i++)
		gpTrainer.checkGradient();

	// Check the gradient while moving away from the stored
	// parameters
	cout << endl << endl << "*** Check that gradient remains OK when moving away from stored parameters *** " << endl;
	
	vec x = ssgp.getParametersVector();
	vec eps = 0.1*ones(x.length());
			
	for (int i=0; i<10; i++) {
		cout << "X = " << endl;
		gpTrainer.checkGradient();
		ssgp.setParametersVector(x+i*eps);
	}

	cout << log(3.2) << " " << log(1.4) << " " << log(0.01) << endl;
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
