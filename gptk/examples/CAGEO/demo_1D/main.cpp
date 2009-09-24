/**
 * This demonstration program illustrates the functioning of the 
 * PSGP algorithm. Noisy observations are sampled from a Gaussian Process
 * with known parameters. We then display the posterior distribution as 
 * the number of active points is increased.
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
    double nugget = 0.05;              // The noise variance

    // Covariance function: Gaussian + Nugget
    GaussianCF   gaussianCovFunc(range, sill);            
    WhiteNoiseCF nuggetCovFunc(nugget);
    SumCovarianceFunction covFunc(gaussianCovFunc);
    covFunc.addCovarianceFunction(nuggetCovFunc);

    // Generate some data from a GP
    int n_points = 200;
    Xtst = 15.0*(randu(n_points)-0.5);       // Sample the input randomly
    sort(Xtst);                              // Sort the input

    mat K = zeros(Xtst.length(), Xtst.length());
    covFunc.computeSymmetric(K, Xtst);              // Covariance of the inputs
    Ytst = (chol(K)).transpose()*randn(n_points);   // Outputs

    // Training set - use a random subsample from the data
    ivec itrn = to_ivec(linspace(0,Xtst.length()-1,40));
    Xtrn = Xtst(itrn);
    Ytrn = Ytst(itrn);
    mat Xtrnmat = Xtrn;

    // Initialise the PSGP
    int n_active = 8;    // Start with 8 active points

    PSGP psgp(Xtrnmat, Ytrn, covFunc, n_active);

    // Use a Gaussian likelihood model with fixed variance (set to 
    // the nugget variance) 
    GaussianLikelihood gaussLik(nugget);

    
    // Compute the PSGP posterior distribution under the Gaussian likelihood
    // for an increasing number of active points
    for (int i=0; i<5; i++) 
    {
        cout << "PSGP with " << n_active << " active points" << endl;

        // Recompute the posterior - the best active points are selected from
        // the 40 points in the training set
        psgp.resetPosterior();
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
        gplot.plotPoints(Xtrn, Ytrn, "training set", CROSS, RED);  

        // Plot active points
        vec activeX = (psgp.getActiveSetLocations()).get_col(0);
        vec activeY = Ytrn(psgp.getActiveSetIndices());
        gplot.plotPoints(activeX, activeY, "active points", CIRCLE, BLUE);

        // Increase the number of active point for next iteration
        n_active += 8;
        psgp.setActiveSetSize(n_active);

        cout << "Press a key to continue" << endl;
        getchar();
    }


    return 0;   
}
