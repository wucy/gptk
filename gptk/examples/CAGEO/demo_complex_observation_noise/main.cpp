/**
 * This demonstration program illustrates the functioning of the 
 * PSGP algorithm. Observations are sampled from a Gaussian Process
 * with known parameters and noise of different types is used to corrupt
 * the data. We compare results with and without taking into account the
 * multiple noise models.
 * 
 * (c) 2009, Remi Barillec <r.barillec@aston.ac.uk>  
 **/




#include <iostream>

#include "itpp/itbase.h"
#include "itpp/itstat.h"

#include "optimisation/SCGModelTrainer.h"

#include "gaussianProcesses/PSGP.h"

#include "likelihoodModels/GaussianLikelihood.h"
#include "likelihoodModels/ExponentialSampLikelihood.h"

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
    double nugget = 0.0001;            // The noise variance

    // Covariance function: Gaussian + Nugget
    GaussianCF   gaussianCovFunc(range, sill);            
    WhiteNoiseCF nuggetCovFunc(nugget);
    SumCovarianceFunction covFunc(gaussianCovFunc);
    covFunc.addCovarianceFunction(nuggetCovFunc);

    // Generate some data from a GP
    int n_points = 200;
    Xtst = 12.0*(randu(n_points)-0.5);       // Sample the input randomly
    sort(Xtst);                              // Sort the input

    mat K = zeros(Xtst.length(), Xtst.length());
    covFunc.computeSymmetric(K, Xtst);      // Covariance of the inputs
    Ytst = (chol(K)).transpose()*randn(n_points);   // Outputs

    // Training set - use a random subsample from the data
    int n_train = 60;     // 60 training points
    ivec itrn = to_ivec( linspace(0,Xtst.length()-1,n_train) );
    Xtrn = Xtst(itrn);
    Ytrn = Ytst(itrn);

    // Noise model 1 is gaussian(0,0.1)
    // Noise model 2 is exponential(1.0)
    // Noise model 3 is gaussian(0,0.4)
    double sigma1 = 0.01;
    double lambda = 1.0;
    double sigma2 = 0.1;
    vec noise1 = sqrt(sigma1)*randn(n_train/3);
    vec noise2 = lambda * exp( -lambda * randu(n_train/3) );
    vec noise3 = sqrt(sigma2)*randn(n_train/3);

    Ytrn += concat(noise1, noise2, noise3);
    mat Xtrnmat = Xtrn;

    // Initialise the PSGP
    int n_active = n_train;    

    gaussianCovFunc.setParameter(0, 0.5);
    nuggetCovFunc.setParameter(0, 2.0);
    PSGP psgp(Xtrnmat, Ytrn, covFunc, n_active);
    psgp.setGammaTolerance(1e-6);
    
    // Use a Gaussian likelihood model with fixed variance (set to 
    // the nugget variance) 
    GaussianLikelihood gaussLik(0.01);
    
    // Compute the PSGP posterior distribution under the Gaussian likelihood
    // Recompute the posterior - the best active points are selected from
    // the training set
    psgp.computePosterior(gaussLik);
    
    SCGModelTrainer scg(psgp);
    
    for(int i=0; i<5; i++)
    {
        scg.Train(5);
        psgp.resetPosterior();
        psgp.computePosterior(gaussLik);
    }
        
    // Note that we are not using the covariange function without the
    // nugget term for prediction.
    psgpmean = zeros(n_points);
    psgpvar  = zeros(n_points);
    psgp.makePredictions(psgpmean, psgpvar, mat(Xtst), gaussianCovFunc);

    
    // Recompute with multiple likelihood models
    Vec<LikelihoodType*> multiLik(n_train);
    ivec multiLikIndex = to_ivec(linspace(0,n_train-1,n_train));
    for (int i=0; i<n_train; i++)
    {
        if (i<n_train/3) 
            multiLik[i] = new GaussianLikelihood(sigma1);
        else if (i<2*n_train/3)
            multiLik[i] = new ExponentialSampLikelihood(lambda);
        else 
            multiLik[i] = new GaussianLikelihood(sigma2);
    }
    
    psgp.computePosterior(multiLikIndex, multiLik);
        
    for(int i=0; i<5; i++)
    {
        scg.Train(5);
        psgp.resetPosterior();
        psgp.computePosterior(multiLikIndex, multiLik);
    }
    
    
    // Plot PSGP mean and error bars
    gplot.plotPoints(Xtst, psgpmean, "psgp mean", LINE, BLUE);
    gplot.plotPoints(Xtst, psgpmean + 2.0*sqrt(psgpvar), "error bar", LINE, BLUE);
    gplot.plotPoints(Xtst, psgpmean - 2.0*sqrt(psgpvar), "error bar", LINE, BLUE);

    psgp.makePredictions(psgpmean, psgpvar, mat(Xtst), gaussianCovFunc);
    gplot.plotPoints(Xtst, psgpmean, "psgp mean", LINE, GREEN);
    gplot.plotPoints(Xtst, psgpmean + 2.0*sqrt(psgpvar), "error bar", LINE, GREEN);
    gplot.plotPoints(Xtst, psgpmean - 2.0*sqrt(psgpvar), "error bar", LINE, GREEN);
    
    // Plot observations
    gplot.plotPoints(Xtst, Ytst, "GP", LINE, RED);
    gplot.plotPoints(Xtrn, Ytrn, "Obs", CROSS, RED);
          
    // Plot active points
    vec activeX = (psgp.getActiveSetLocations()).get_col(0);
    vec activeY = Ytrn(psgp.getActiveSetIndices());
    gplot.plotPoints(activeX, activeY, "active points", CIRCLE, BLUE);

    return 0;   
}
