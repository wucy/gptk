/**
 * This demonstration program illustrates the functioning of the 
 * PSGP algorithm. The evolution of the posterior process is shown
 * as we iterate through the observations. 
 * 
 * (c) 2009, Remi Barillec <r.barillec@aston.ac.uk>  
 **/


#include <iostream>

#include "itpp/itbase.h"
#include "itpp/itstat.h"

#include "optimisation/SCGModelTrainer.h"

#include "gaussianProcesses/PSGP.h"

#include "likelihoodModels/GaussianLikelihood.h"

#include "covarianceFunctions/SumCovarianceFunction.h"
#include "covarianceFunctions/GaussianCF.h"
#include "covarianceFunctions/WhiteNoiseCF.h"

#include "GraphPlotter/GraphPlotter.h"

using namespace std;
using namespace itpp;


int main(void)
{
    vec Xtrn, Xtst, Ytrn, Ytst;
    vec psgpmean, psgpvar;

    GraphPlotter gplot = GraphPlotter();

    // Generate some data from a GP
    double range  = 2.0;               // The range or length scale of the GP
    double sill   = 2.0;               // The sill or variance of the GP
    double nugget = 0.02;              // The noise variance

    // Covariance function: Gaussian + Nugget
    GaussianCF   gaussianCovFunc(range, sill);            
    WhiteNoiseCF nuggetCovFunc(nugget);
    SumCovarianceFunction covFunc(gaussianCovFunc);
    covFunc.addCovarianceFunction(nuggetCovFunc);

    // Generate some data from a GP
    int n_points = 200;
    Xtst = linspace(-10,10,n_points);       // Sample the input randomly

    mat K = zeros(Xtst.length(), Xtst.length());
    gaussianCovFunc.computeSymmetric(K, Xtst);              // Covariance of the inputs
    Ytst = (chol(K+0.0001*eye(n_points))).transpose()*randn(n_points);   // Outputs

    // Training set - use a random subsample from the data
    int n_train = 64;
    ivec itrn = to_ivec(linspace(0,Xtst.length()-1,n_train));
    Xtrn = Xtst(itrn);
    Ytrn = Ytst(itrn) + sqrt(nugget)*randn(n_train);

    // Use a Gaussian likelihood model with fixed variance (set to 
    // the nugget variance) 
    GaussianLikelihood gaussLik(nugget);

    // Use a fixed set of active points/basis vectors
    int n_active = 8;    // Start with 8 active points
    ivec active_indices = to_ivec( linspace(0, n_train-1, n_active) );
    
    // Compute the PSGP posterior distribution under the Gaussian likelihood
    // after more and more observations have been seen. The active set remains
    // fixed.
    int subset_size = 1;
    
    for (int i=0; i<7; i++) 
    {
        cout << "Posterior process after seeing " << subset_size << " observation(s)" << endl;
        
        // Keep 1, 2, 4...64 observations from the training set
        // + the 8 active points (at the begining)
        ivec subset_indices = to_ivec( linspace(0,subset_size-1,subset_size) );
        
        mat X_subset = mat( Xtrn( concat( active_indices, subset_indices ) ) );
        vec Y_subset = Ytrn( concat( active_indices, subset_indices ) );

        // PSGP with fixed active set - note that the number of iterations
        // through the data with replacement is set to 0 to ensure 
        // the active set remains fixed. The active set comprises the first 8
        // observations (indices 0 to 7).
        PSGP psgp(X_subset, Y_subset, covFunc, n_active, 0, 1);
        psgp.setGammaTolerance(0.0);
        psgp.setActiveSet( ivec("0:7"), gaussLik );
        
        // Compute the posterior - the active set is now fixed
        psgp.computePosterior(gaussLik);
        
        // Note that we are not using the covariange function without the
        // nugget term for prediction.
        psgpmean = zeros(n_points);
        psgpvar  = zeros(n_points);
        psgp.makePredictions(psgpmean, psgpvar, mat(Xtst), gaussianCovFunc);

        // Plot PSGP mean and error bars
        gplot.clearPlot();
        gplot.plotPoints(Xtst, psgpmean, "psgp mean", LINE, BLUE);
        gplot.plotPoints(Xtst, psgpmean + 2.0*sqrt(psgpvar), "error bar", LINE, CYAN);
        gplot.plotPoints(Xtst, psgpmean - 2.0*sqrt(psgpvar), "error bar", LINE, CYAN);

        // Plot observations
        gplot.plotPoints(X_subset.get_col(0), Y_subset, "training set", CROSS, RED);  
        gplot.plotPoints(Xtst, Ytst, "true GP", LINE, RED);

        // Plot active points
        vec activeX = (psgp.getActiveSetLocations()).get_col(0);
        vec activeY = Y_subset(psgp.getActiveSetIndices());
        gplot.plotPoints(activeX, activeY, "active points", CIRCLE, BLUE);

        // Double the number of observations for the next iteration
        subset_size *= 2;
        psgp.setActiveSetSize(n_active);
        
        cout << "Press a key to continue" << endl;
        getchar();
    }


    cout << "End of demonstration" << endl;
    
    return 0;   
}
