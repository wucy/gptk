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
#include "likelihoodModels/GaussianSampLikelihood.h"
#include "likelihoodModels/ExponentialSampLikelihood.h"

#include "covarianceFunctions/SumCovarianceFunction.h"
#include "covarianceFunctions/GaussianCF.h"
#include "covarianceFunctions/WhiteNoiseCF.h"

#include "GraphPlotter/GraphPlotter.h"

using namespace std;
using namespace itpp;

/**
 * Observation operator for 3rd noise type (the data, x, is not observed
 * directly. Instead, a quartic expression is observed.)  
 */
double obsOperator (double x) 
{ 
    return  pow(x/3.0, 4.0); 
}


int main(int argc, char* argv[])
{
    vec Xtrn, Xtst, Ytrn, Ytst;
    vec psgpmean, psgpvar;

    GraphPlotter gplot = GraphPlotter();

    // Generate some data from a GP
    double range  = 10.0;               // The range or length scale of the GP
    double sill   = 10.0;               // The sill or variance of the GP
    double nugget = 0.0001;            // The noise variance

    // Covariance function: Gaussian + Nugget
    GaussianCF   gaussianCovFunc(range, sill);            
    WhiteNoiseCF nuggetCovFunc(nugget);
    SumCovarianceFunction covFunc(gaussianCovFunc);
    covFunc.addCovarianceFunction(nuggetCovFunc);

    // Generate some data from a GP
    int n_points = 1000;
    Xtst = linspace(0,210,n_points);

    mat K = zeros(Xtst.length(), Xtst.length());
    gaussianCovFunc.computeSymmetric(K, Xtst);      // Covariance of the inputs
    Ytst = (chol(K+0.0001*eye(n_points))).transpose()*randn(n_points);   // Outputs

    // Training set - use a random subsample from the data
    int n_train = 200;
    ivec itrn = to_ivec( linspace(0,Xtst.length()-1,n_train) );
    Xtrn = Xtst(itrn);
    Ytrn = Ytst(itrn);

    // Noise model 1 is gaussian(0,sigma1)
    // Noise model 2 is exponential(lambda)
    // Noise model 3 is gaussian(0,sigma3) 
    double sigma1 = 2.0;
    double lambda = 0.7;
    double mu3 = 3.0;
    double sigma3 = 0.5;
    vec noise1 = sqrt(sigma1)*randn(70);
    vec noise2 = - log(randu(70))/lambda;
    vec noise3 = mu3 + sqrt(sigma3)*randn(60);  

    Ytrn.set_subvector(0, 69, Ytrn(0,69) + noise1);
    Ytrn.set_subvector(70, 139, Ytrn(70,139) + noise2); 
    
    // An observation operator is applied to the third segment: 1/sigma2 * y^4 + 2.0
    Ytrn.set_subvector(140, 199, apply_function( &obsOperator, Ytrn(140, 199) ) + noise3);   
    
    mat Xtrnmat = Xtrn;

    // Initialise the PSGP
    int n_active = 50;    

    // gaussianCovFunc.setParameter(0, 0.5);
    // nuggetCovFunc.setParameter(0, 2.0);
    PSGP psgp(Xtrnmat, Ytrn, covFunc, n_active, 1, 2);
    // psgp.setGammaTolerance(1e-6);
    
    // Use a Gaussian likelihood model with fixed variance (set to 
    // the empirical observation variance) 
    GaussianLikelihood gaussLik( itpp::variance( concat(noise1, noise2, noise3) ) );
    
    // Compute the PSGP posterior distribution under the Gaussian likelihood
    // Recompute the posterior - the best active points are selected from
    // the training set
    psgp.computePosterior(gaussLik);
    
    SCGModelTrainer scg(psgp);
    
    for(int i=0; i<0; i++)
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
    vec activeX = (psgp.getActiveSetLocations()).get_col(0);
    vec activeY = Ytrn(psgp.getActiveSetIndices());
        
    
    // Recompute with multiple likelihood models
    ivec multiLikIndex = to_ivec( concat( zeros(70), ones(70), 2*ones(60) ) );
    Vec<LikelihoodType*> multiLik(3);
    multiLik[0] = new GaussianLikelihood(sigma1);
    multiLik[1] = new ExponentialSampLikelihood(lambda);
    multiLik[2] = new GaussianSampLikelihood(mu3, sigma3, &obsOperator);
   
    psgp.resetPosterior();
    psgp.computePosterior(multiLikIndex, multiLik);
        
    for(int i=0; i<0; i++)
    {
        scg.Train(5);
        psgp.resetPosterior();
        psgp.computePosterior(multiLikIndex, multiLik);
    }
    
    
    // Plot observations
    gplot.plotPoints(Xtst, Ytst, "GP", LINE, RED);
    gplot.plotPoints(Xtrn, Ytrn, "Obs", CIRCLE, RED);
    
    // Plot PSGP mean and error bars (Gaussian likelihood)
    gplot.plotPoints(Xtst, psgpmean, "psgp mean", LINE, GREEN);
    gplot.plotPoints(Xtst, psgpmean + 2.0*sqrt(psgpvar), "error bar", LINE, GREEN);
    gplot.plotPoints(Xtst, psgpmean - 2.0*sqrt(psgpvar), "error bar", LINE, GREEN);
    gplot.plotPoints(activeX, activeY, "active points", CIRCLE, GREEN);
    
    // Predictions and active points for PSGP (Multiple likelihood)
    psgp.makePredictions(psgpmean, psgpvar, mat(Xtst), gaussianCovFunc);
    
    activeX = (psgp.getActiveSetLocations()).get_col(0);
    activeY = Ytrn(psgp.getActiveSetIndices());
        
    // Plot PSGP mean and error bars (Multiple likelihood)
    gplot.plotPoints(Xtst, psgpmean, "psgp mean", LINE, BLUE);
    gplot.plotPoints(Xtst, psgpmean + 2.0*sqrt(psgpvar), "error bar", LINE, BLUE);
    gplot.plotPoints(Xtst, psgpmean - 2.0*sqrt(psgpvar), "error bar", LINE, BLUE);
    gplot.plotPoints(activeX, activeY, "active points", CIRCLE, BLUE);

    
          
    return 0;
       
}
