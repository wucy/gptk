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

#include "io/csvstream.h"

#include "gaussianProcesses/PSGP.h"
#include "gaussianProcesses/GaussianProcess.h"

#include "likelihoodModels/GaussianLikelihood.h"
#include "likelihoodModels/GaussianSampLikelihood.h"
#include "likelihoodModels/ExponentialSampLikelihood.h"

#include "covarianceFunctions/SumCovarianceFunction.h"
#include "covarianceFunctions/Matern3CF.h"
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
    return 2.0*pow(x/3.0,3.0) - 2.0; // pow(x/3.0, 4.0) + 5.0; 
}


mat computeLikelihoodProfile(PSGP &psgp, mat paramRanges);
mat computeLikelihoodProfileGP(GaussianProcess &gp, mat paramRanges);

int main(int argc, char* argv[])
{
    vec Xtrn, Xtst, Ytrn, Ytst;

    GraphPlotter gplot = GraphPlotter();

    // Generate some data from a GP
    double range  = 20.0;               // The range or length scale of the GP
    double sill   = 10.0;               // The sill or variance of the GP
    double nugget = 0.02;               // The noise variance

    // Covariance function: Gaussian + Nugget
    Matern3CF   maternCovFuncGP(range, sill);            
    WhiteNoiseCF nuggetCovFuncGP(nugget);
    SumCovarianceFunction gpCovFunc(maternCovFuncGP);
    gpCovFunc.addCovarianceFunction(nuggetCovFuncGP);
    
    Matern3CF   maternCovFuncPSGP(range, sill);            
    WhiteNoiseCF nuggetCovFuncPSGP(nugget);
    SumCovarianceFunction psgpCovFunc(maternCovFuncPSGP);
    psgpCovFunc.addCovarianceFunction(nuggetCovFuncPSGP);
        

    // Generate some data from a GP
    Matern3CF dataCovFunc(range, sill);
    int n_points = 1000;
    Xtst = linspace(0,210,n_points);

    mat K = zeros(Xtst.length(), Xtst.length());
    dataCovFunc.computeSymmetric(K, Xtst);      // Covariance of the inputs
    Ytst = (chol(K)).transpose()*randn(n_points);   // Outputs

    // Training set - use a random subsample from the data
    int n_train = 200;
    ivec itrn = to_ivec( linspace(0,Xtst.length()-1,n_train) );
    Xtrn = Xtst(itrn);
    Ytrn = Ytst(itrn);

    // Noise model 1 is gaussian(0,sigma1)
    // Noise model 2 is exponential(lambda)
    // Noise model 3 is gaussian(0,sigma3) 
    double sigma1 = 1.0;
    double lambda = 0.7;
    double mu3 = 0.0;
    double sigma3 = 0.1;
    vec noise1 = sqrt(sigma1)*randn(70);
    vec noise2 = - log(randu(70))/lambda;
    vec noise3 = sqrt(sigma3)*randn(60);  

    Ytrn.set_subvector(0, 69, Ytrn(0,69) + noise1);
    Ytrn.set_subvector(70, 139, Ytrn(70,139) + noise2); 
    
    // An observation operator is applied to the third segment: 1/sigma2 * y^4 + 2.0
    Ytrn.set_subvector(140, 199, apply_function( &obsOperator, Ytrn(140, 199) ) + noise3);   
    
    mat Xtrnmat = Xtrn;

        
    //-------------------------------------------------------------------------
    // Regression using PSGP with multiple likelihood models
    
    // Initialise the PSGP
    int n_active = 50;
    
    //nuggetCovFuncPSGP.setParameter(0, itpp::variance( concat(noise1, noise2, noise3) ) );
    
    PSGP psgp(Xtrnmat, Ytrn, psgpCovFunc, n_active, 1, 1);
    psgp.setLikelihoodType(UpperBound);
    
    ivec multiLikIndex = to_ivec( concat( zeros(70), ones(70), 2*ones(60) ) );
    Vec<LikelihoodType*> multiLik(3);
    multiLik[0] = new GaussianLikelihood(sigma1);
    multiLik[1] = new ExponentialSampLikelihood(lambda);
    multiLik[2] = new GaussianSampLikelihood(mu3, sigma3, &obsOperator);
   
    // Compute the posterior - the best active points are selected from
    // the training set
    psgp.resetPosterior();
    psgp.computePosterior(multiLikIndex, multiLik);
        
    // Estimate the covariance function parameters
    SCGModelTrainer scg(psgp);   // SCG model trainer for parameter estimation
    bvec optMask(3);
    optMask(0) = true; optMask(1) = true; optMask(2) = true;
    scg.setOptimisationMask(optMask);
    for(int i=0; i<3; i++)
    {
        scg.Train(5);
        psgp.resetPosterior();
        psgp.computePosterior(multiLikIndex, multiLik);
    }
    psgpCovFunc.displayCovarianceParameters();
    cout << endl << endl;
    maternCovFuncPSGP.displayCovarianceParameters();
    
    cout << "------------------------------------------------------------" << endl;
    
    // Predictions and active points for PSGP (Multiple likelihood)
    vec psgpmean = zeros(n_points);
    vec psgpvar  = zeros(n_points);
    psgp.makePredictions(psgpmean, psgpvar, mat(Xtst), maternCovFuncPSGP);

    // Active set
    ivec iactive = psgp.getActiveSetIndices();
    mat Xactive = psgp.getActiveSetLocations();
    vec Yactive = psgp.getActiveSetObservations();
    
    //-------------------------------------------------------------------------
    // Regression using standard GP and Gaussian likelihood
    // Use a Gaussian likelihood model with fixed variance (set to 
    // the empirical observation variance) 
    nuggetCovFuncGP.setParameter(0, itpp::variance( concat(noise1, noise2, noise3) ) );
    
    GaussianProcess gp(1, 1, Xtrnmat, Ytrn, gpCovFunc);

    // Estimate the covariance function parameters
    SCGModelTrainer scg2(gp);
    scg2.setOptimisationMask(optMask);
    for(int i=0; i<3; i++)
    {
        scg2.Train(5);
    }

    // Note that we are not using the covariange function without the
    // nugget term for prediction.
    vec gpmean = zeros(n_points);
    vec gpvar  = zeros(n_points);
    gp.makePredictions(gpmean, gpvar, mat(Xtst), maternCovFuncGP);
    
    gpCovFunc.displayCovarianceParameters();

    /*
    //-------------------------------------------------------------------------
    // Plot observations
    gplot.plotPoints(Xtst, Ytst, "True function", LINE, RED);
    gplot.plotPoints(Xtrn, Ytrn, "Obs", CIRCLE, RED);

    // Plot GP mean and error bars (Gaussian likelihood)
    gplot.plotPoints(Xtst, gpmean, "gp mean", LINE, GREEN);
    gplot.plotPoints(Xtst, gpmean + 2.0*sqrt(gpvar), "gp error bar", LINE, GREEN);
    gplot.plotPoints(Xtst, gpmean - 2.0*sqrt(gpvar), "gp error bar", LINE, GREEN);


    // Plot PSGP mean and error bars (Multiple likelihood)
    gplot.plotPoints(Xtst, psgpmean, "psgp mean", LINE, BLUE);
    gplot.plotPoints(Xtst, psgpmean + 2.0*sqrt(psgpvar), "error bar", LINE, BLUE);
    gplot.plotPoints(Xtst, psgpmean - 2.0*sqrt(psgpvar), "error bar", LINE, BLUE);
    gplot.plotPoints(Xactive.get_col(0), Yactive, "active points", CIRCLE, BLUE);
*/
    
    //-------------------------------------------------------------------------
    // Save data for Matlab plotting

    csvstream csv;
    
    mat active_set(2, n_active);             // Active set
    active_set.set_row(0, Xactive.get_col(0));
    active_set.set_row(1, Yactive);
    csv.write(active_set, "demo_noise_active_set.csv");
    
    mat gpmv = zeros(2, n_points);           // Mean and variance of GP
    gpmv.set_row(0, gpmean);
    gpmv.set_row(1, gpvar);
    csv.write(gpmv, "demo_noise_gp.csv");
    
    mat psgpmv = zeros(2, n_points);         // Mean and variance of PSGP
    psgpmv.set_row(0, psgpmean);
    psgpmv.set_row(1, psgpvar);
    csv.write(psgpmv, "demo_noise_psgp.csv");
    
    mat test = zeros(2,n_points);            // True function
    test.set_row(0, Xtst);
    test.set_row(1, Ytst);
    csv.write(test, "demo_noise_test.csv");
    
    mat obs = zeros(2,n_train);              // Observations
    obs.set_row(0, Xtrn);
    obs.set_row(1, Ytrn);
    csv.write(obs, "demo_noise_obs.csv");
    
    cout << "Done" << endl;


    //-------------------------------------------------------------------------
    // LIKELIHOOD PROFILES
    
    // Parameter ranges (for likelihood profile)
    int n = 100;
    mat paramRanges(3,n);
    paramRanges.set_row(0, linspace(range*0.1, range*2.0, n));
    paramRanges.set_row(1, linspace(sill*0.01,  sill*2.0, n));
    paramRanges.set_row(2, linspace(nugget*0.1, nugget*10.0, n));

    // Likelihood profile for full GP
    mat profilesGP = computeLikelihoodProfileGP(gp, paramRanges);;
    csv.write(profilesGP, "proflik.complex_noise.gp.csv");
    csv.write(paramRanges, "proflik.complex_noise.params.csv");
    
    n_active = 25;
    psgp.setActiveSetSize(n_active);
    
    for (int i=0; i<4; i++) 
    {
        //---------------------------------------------------------------------
        cout << "PSGP with " << n_active << " active points" << endl;

        // Recompute the posterior - the best active points are selected from
        // the training set
        psgp.resetPosterior();
        psgp.computePosterior( multiLikIndex, multiLik );

        // Note that we are not using the covariange function without the
        // nugget term for prediction.
        psgpmean = zeros(n_points);
        psgpvar  = zeros(n_points);
        psgp.makePredictions(psgpmean, psgpvar, mat(Xtst), maternCovFuncPSGP);

        //---------------------------------------------------------------------
        // Compute likelihood profile
        mat profile = computeLikelihoodProfile(psgp, paramRanges);
        stringstream filename;
        filename << "proflik.complex_noise.psgp." << n_active << ".csv";
        csv.write(profile,filename.str());    // Save profile to CSV file
    

        //---------------------------------------------------------------------
        // Increase the number of active point for next iteration
        n_active *= 2;
        psgp.setActiveSetSize(n_active);

        cout << "Press a key to continue" << endl;
    }

    return 0;

}






/**
 * Compute likelihood profiles for PSGP for a given number
 * of active points
 */ 
mat computeLikelihoodProfile(PSGP &psgp, mat paramRanges)
{
    cout << "Computing likelihood profiles for PSGP" << endl;
    
    vec params = psgp.getParametersVector(); 
    int nparams = params.length();
    mat profiles(nparams,paramRanges.cols());
 
    for(int i=0; i<nparams; i++) 
    {
        // Compute likelihood at each point in parameter range
        for (int j=0; j<paramRanges.cols(); j++)
        {
            vec new_params = params;
            new_params(i) = log(paramRanges(i,j));   // Parameters are in log space by default
            psgp.setParametersVector(new_params);
            profiles(i,j) = psgp.objective();
        }
    }
    psgp.setParametersVector(params);
    
    return profiles;
}



/**
 * Compute likelihood profiles for a full GP
 */ 
mat computeLikelihoodProfileGP(GaussianProcess &gp, mat paramRanges)
{
    cout << "Computing likelihood profiles for GP" << endl;
    
    vec params = gp.getParametersVector(); 
    int nparams = params.length();
    mat profiles(nparams,paramRanges.cols());
 
    for(int i=0; i<nparams; i++) 
    {
        // Compute likelihood at each point in parameter range
        for (int j=0; j<paramRanges.cols(); j++)
        {
            vec new_params = params;
            new_params(i) = log(paramRanges(i,j));   // Parameters are in log space by default
            gp.setParametersVector(new_params);
            profiles(i,j) = gp.objective();
        }
    }
    gp.setParametersVector(params);
    
    return profiles;
}
