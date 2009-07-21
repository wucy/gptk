/********************************************************************
 * 
 * Spare Sequential Gaussian Process on the Helicopter scenario
 * 
 * Author: Remi Barillec <r.barillec@aston.ac.uk>
 * 
 * Copyright 2009,  R. Barillec, B. Ingram, D. Cornford
 *  
 *******************************************************************/

#include <iostream>
#include <cstring>

#include "itpp/itbase.h"
#include "itpp/itstat.h"
#include "io/csvstream.h"

#include "covarianceFunctions/SumCovarianceFunction.h"
#include "covarianceFunctions/GaussianCF.h"
#include "covarianceFunctions/WhiteNoiseCF.h"

#include "optimisation/SCGModelTrainer.h"
#include "optimisation/CGModelTrainer.h"
#include "optimisation/QuasiNewtonModelTrainer.h"

#include "design/GreedyMaxMinDesign.h"

#include "likelihoodModels/GaussianLikelihood.h"

#include "gaussianProcesses/SequentialGP.h"

/**
 * Print error message and abort program
 */
void die(string msg)
{
    cerr << "** Error: " << msg << endl;
    exit(1);
}


/**
 * Read in some data and split it into inputs
 * and outputs (assuming last column is output)
 **/  
int loadData(string filename, mat &X, vec &Y)
{
    csvstream csv;
    mat M;
    
    // Read data - return error if problem
    if (csv.read(M,filename)) return 1;
    
    // Extract inputs and outputs 
    int nin = M.cols() - 1;
    int nobs = M.rows();
    
    X = zeros(nobs, nin);
    
    for (int i=0; i<nin; i++) X.set_col(i,M.get_col(i));
    
    Y = M.get_col(nin);
    
    return 0;
}


/**
 * Main method
 **/
int main() 
{
    mat X, Xtst;
    vec Y, Ytst;
    vec ssgpmean, ssgpvar;
    int n_active;
    
    double range  = 0.5;                     
    double sill   = 1.0;
    double nugget = 0.01;
    
    // Read data in
    cout << "Loading data" << endl;

    if (loadData("helicopter_normalised.csv", X, Y)) {
        die("Failed to load data. Aborting.");
    }
    
    // Test set - uniform grid between input min and input max
    int nin = X.cols();
    int ngrid = 200;
    Xtst = zeros(ngrid*ngrid, nin);
    
    vec X1 = linspace(min(X.get_col(0)), max(X.get_col(0)), ngrid);
    vec X2 = linspace(min(X.get_col(1)), max(X.get_col(1)), ngrid);
    
    for (int i=0; i<ngrid; i++) 
    {
        for (int j=0; j<ngrid; j++) {
            Xtst(i*ngrid+j,0) = X1(i);
            Xtst(i*ngrid+j,1) = X2(j);
        }
    }
    
    n_active = 100; //Y.length();

    // Covariance function: Gaussian + Nugget
    GaussianCF   g1(range, sill);
    WhiteNoiseCF g2(nugget);
    SumCovarianceFunction gCmp(g1);
    gCmp.addCovarianceFunction(g2);

    // Initialise sequential GP
    cout << "Initialising SSGP" << endl;
    SequentialGP ssgp(nin, 1, n_active, X, Y, gCmp);

    // Gaussian observation likelihood
    GaussianLikelihood gaussLik(nugget);

    // Select active set
    cout << "Selecting active set" << endl;
    
    GreedyMaxMinDesign design;
    mat XY(X.rows(), X.cols()+1);
    XY.set_cols(0,X);
    XY.set_col(2,Y);
    
    
    ivec iActive = design.subsample(X, n_active);
    cout << iActive << endl;

    // Compute posterior
    cout << "Computing posterior" << endl;
    // ssgp.computePosterior(gaussLik);
    ssgp.computePosteriorFixedActiveSet(gaussLik, iActive);
    
    // Learn parameters
    SCGModelTrainer gpTrainer(ssgp);
    ssgp.setLikelihoodType(UpperBound);
    bvec optMask(3);
    
    gpTrainer.setCheckGradient(true);
        
    cout << "Optimising parameters" << endl;
    int niter = 5;
    for (int i=0; i<niter; i++) 
    {
        cout << endl << endl << "-- " << i+1 << "/" << niter;
        cout << " ------------------------------------------------------------------------------" << endl;
        
        // Optimise covariance parameters 
        optMask(0) = true;
        optMask(1) = true;
        optMask(2) = false;
        gpTrainer.setOptimisationMask(optMask);

        gpTrainer.Train(5);
        gpTrainer.checkGradient();
        // g2.displayCovarianceParameters();

        // Recompute basis vectors 
        ssgp.resetPosterior();
        // ssgp.computePosterior(gaussLik);
        ssgp.computePosteriorFixedActiveSet(gaussLik,iActive);

        // Optimise noise parameters
        optMask(0) = false;
        optMask(1) = false;
        optMask(2) = true;
        gpTrainer.setOptimisationMask(optMask);
        gpTrainer.Train(5);
        gpTrainer.checkGradient();
        // g2.displayCovarianceParameters();

        ssgp.resetPosterior();
        ssgp.computePosteriorFixedActiveSet(gaussLik,iActive);
        // ssgp.computePosterior(gaussLik);

    }
    
    gCmp.displayCovarianceParameters();
    
    cout << "Making prediction" << endl;
    ssgp.makePredictions(ssgpmean, ssgpvar, Xtst, g1);
    
    mat xActive = ssgp.getActiveSetLocations();
    // cout << xActive << endl;
    // xActive.set_size(ngrid, nin, true);
    
    cout << "Saving data" << endl;
    mat data = zeros(ngrid*ngrid, nin+4);
    
    data.set_cols(0,Xtst);
    data.set_col(nin, ssgpmean);                         // Mean prediction
    data.set_col(nin+1, ssgpvar);                        // Var prediction
    data.set_col(nin+2, xActive.get_col(0));                        // Active locations
    data.set_col(nin+3, xActive.get_col(1));                        // Active locations
    // data.set_col(nin+2+nin, Y(ssgp.getActiveSetIndices()));  // Active indices
    
    // Save prediction data
    csvstream csv;
    if (csv.write(data, "helicopter_pred.csv"))
        die("Could not write prediction data to file.");
    
    // Save ssgp parameters
    if (csv.write(ssgp.getParametersVector(), "helicopter_params.csv"))
        die("Could not write parameters to file.");
    
    
    cout << "End of helicopter example" << endl;
    
    return 0;
    
}

