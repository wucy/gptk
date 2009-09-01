#ifndef SPATIALEXAMPLE_H_
#define SPATIALEXAMPLE_H_

#include <iostream>
#include <string>

#include "itpp/itbase.h"
// #include "itpp/itstat.h"

#include "itppext/itppext.h"
#include "io/csvstream.h"

//#include <cassert>
// #include <algorithm>

#include "gaussianProcesses/SequentialGP.h"
#include "gaussianProcesses/GaussianProcess.h"
#include "optimisation/SCGModelTrainer.h"

#include "likelihoodModels/GaussianLikelihood.h"

#include "covarianceFunctions/GaussianCF.h"
#include "covarianceFunctions/ExponentialCF.h"
#include "covarianceFunctions/WhiteNoiseCF.h"
#include "covarianceFunctions/SumCovarianceFunction.h"

using namespace std;
using namespace itpp;

enum ParameterEstimationMethod { PARAM_ESTIM_GP, PARAM_ESTIM_PSGP, PARAM_ESTIM_NO_ESTIMATION, 
                                  PARAM_ESTIM_CUSTOM };

class SpatialExample
{
    
public:
    SpatialExample();
    SpatialExample(bool normaliseData);
	virtual ~SpatialExample();
	
	void run(string datafile, string predfile);
	
	// INPUT/OUTPUT
    bool loadData(string filename);
	bool saveResults(string filename_pred);
	bool saveResults(string filename_pred, string filename_params);

	// COVARIANCE FUNCTION/PARAMETER ESTIMATION
	bool initParameters();
	bool initCovarianceFunction(); 
	bool learnParameters();
	
	void setNumberActivePoints(int value);
	void setNumberOptimIterations(int value);
	void setNumberOuterLoops(int value);
	void setParameterEstimationMethod(ParameterEstimationMethod method);


	// PREDICTION
	void setPredictionLocations(mat Xpred);
	bool makePredictions();
	
	// CURRENT OPTIONS
	void setNumberSweeps(int value);
	void displayCurrentOptions(); 
	
protected:
    mat X, Xpred;            // Input (observed, predicted)
    vec y, ypred;            // Output (observed, predicted)
    vec v, vpred;            // Output variance (observed, predicted)
    
    mat M;                   // Data matrix
    vec Mmean, Mcovdiag;     // Data mean and diagonal covariance (used in normalisation)
    
    int n_active;            // Number of active points to use in PSGP
    int n_sweeps;            // Number of sweeps through the data in PSGP
    
    bool observedNoise;      // Do we have information about obs noise?
    bool normaliseData;      // Do we normalise the data?
    bool isSetPredictionLocations; // Have prediction locations been specified?
    
    double defaultNuggetRatio;  // Default nugget to variance ratio (to initialise nugget)
    double range, sill, nugget; // Covariance function's parameters

    // COVARIANCE FUNCTION
    CovarianceFunction *kernelCF;
    CovarianceFunction *nuggetCF;
    CovarianceFunction *covFunc;
    
    // LIKELIHOOD FUNCTION
    LikelihoodType *likFunction;
    
    string dataFilename;        // Name of data file (to read from)
    string predFilename;        // Name of prediction file (to write to)
    
    void setDefaults();
    
    void setNoiseModel();       // TODO
    void setDefaultNuggetRatio(double value);
        
    bool preProcess();          // This is a template method - for overriding
    bool postProcess();         // This is a template method - for overriding
    
    void uniformGrid(mat X, mat& grid);
    
    // PARAMETER ESTIMATION
    ParameterEstimationMethod paramEstimationMethod;
    int  n_optim_iterations;    // Number of iteration in parameter estimation
    int  n_outer_loops;         // Number of outer loops in optimisation of PSGP
    bool learnParametersGP();
    bool learnParametersPSGP();
    bool learnParametersCustom();
    
};

#endif /*SPATIALEXAMPLE_H_*/
