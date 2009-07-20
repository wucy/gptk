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

#include "design/MaxMinDesign.h"

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
    /*
    double range  = 10.0;                     
    double sill   = 10000.0;
    double nugget = 5.0;
    */
    double range  = 5.0;                     
    double sill   = 2.0;
    double nugget = 0.01;
    
    // Read data in
    cout << "Loading data" << endl;
    /*
    if (loadData("helicopter_shifted.csv", X, Y)) {
        die("Failed to load data. Aborting.");
    }
    */
    
    if (loadData("gp2d_trn.csv", X, Y) || loadData("gp2d_tst.csv", Xtst, Ytst) ) {
        die("Failed to load data. Aborting.");
    }
    
    // Test set - uniform grid between input min and input max
    int nin = X.cols();
    cout << nin << endl;
    int ngrid = Xtst.rows();
    /*
    int ngrid = 100;
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
    */
    
    n_active = 60; //Y.length();

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

    // Predict using SSGP
    cout << "Computing posterior" << endl;
    // ssgp.computePosterior(gaussLik);
    MaxMinDesign design(100);
    ivec iActive = design.subsample(X, Y, n_active);
    ssgp.computePosteriorFixedActiveSet(gaussLik, iActive);
    
    // Learn parameters
    SCGModelTrainer gpTrainer(ssgp);
    ssgp.setLikelihoodType(UpperBound);
    bvec optMask(3);
    
    gpTrainer.setCheckGradient(true);
        
    cout << "Optimising parameters" << endl;
    for (int i=0; i<5; i++) 
    {
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
    xActive.set_size(ngrid, nin, true);
    
    cout << "Saving data" << endl;
    mat data = Xtst;
    data.set_size(ngrid, nin+3+nin, true);
    data.set_col(nin, ssgpmean);
    data.set_col(nin+1, ssgpvar);
    data.set_col(nin+2, xActive.get_col(0));
    data.set_col(nin+3, xActive.get_col(1));
    data.set_col(nin+4, Y(ssgp.getActiveSetIndices()));
    
    // Save prediction data
    csvstream csv;
    if (csv.write(data, "gp2d_pred.csv"))
        die("Could not write prediction data to file.");
    
    // Save ssgp parameters
    if (csv.write(ssgp.getParametersVector(), "gp2d_params.csv"))
            die("Could not write parameters to file.");
    
    
    cout << "End of example" << endl;
    return 0;
    
}

