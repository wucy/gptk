/**
 * This demonstration program illustrates the functioning of the 
 * PSGP algorithm. Noisy observations are sampled from a Gaussian Process
 * with known parameters. We then display the posterior distribution as 
 * the number of active points is increased.
 * 
 * (c) 2009, Remi Barillec <r.barillec@aston.ac.uk>  
 **/

#include "main.h"

int main(void)
{
    vec Xtrn, Xtst, Ytrn, Ytst;
    vec psgpmean, psgpvar, gpmean, gpvar;
    
    GraphPlotter gplot = GraphPlotter();
    csvstream csv;                           // Input/output to CSV file
    stringstream filename;                   // Filename
    
    // True GP parameters
    double range  = 5.0;               // The range or length scale of the GP
    double sill   = 3.0;               // The sill or variance of the GP
    double nugget = 0.01;              // The noise variance
    double noisevar = 0.01;
    
    // Parameter ranges (for likelihood profile)
    int n = 100;
    mat paramRanges(3,n);
    paramRanges.set_row(0, linspace(range*0.1, range*10.0, n));
    paramRanges.set_row(1, linspace(sill*0.1,  sill*10.0, n));
    paramRanges.set_row(2, linspace(nugget*0.1, nugget*10.0, n));
    
    // Covariance function: Gaussian + Nugget
    GaussianCF   gaussianCovFunc(range, sill);            
    WhiteNoiseCF nuggetCovFunc(nugget);
    SumCovarianceFunction covFunc(gaussianCovFunc);
    covFunc.addCovarianceFunction(nuggetCovFunc);

    //-------------------------------------------------------------------------
    // Generate some data from GP with known parameters
    int n_points = 400;
    Xtst = linspace(-40,40, n_points);
    mat K = zeros(Xtst.length(), Xtst.length());
    gaussianCovFunc.computeSymmetric(K, Xtst);              // Covariance of the inputs
    Ytst = (chol(K+0.0001*eye(n_points))).transpose()*randn(n_points);   // Outputs

    // Training set - use a random subsample from the data
    int n_train = 64;
    ivec itrn = to_ivec(linspace(0,Xtst.length()-1,n_train));
    Xtrn = Xtst(itrn);
    Ytrn = Ytst(itrn) + sqrt(noisevar)*randn(n_train);
    mat Xtrnmat = Xtrn;
    
    //-------------------------------------------------------------------------
    // Use a full GP to fit the data
    GaussianProcess gp(1, 1, Xtrnmat, Ytrn, covFunc);
    gpmean = zeros(n_points);
    gpvar = zeros(n_points);
    SCGModelTrainer scggp(gp);
    scggp.Train(15);
    gp.makePredictions(gpmean, gpvar, Xtst, gaussianCovFunc);
    
    // Likelihood profile for full GP
    mat profilesGP = computeLikelihoodProfileGP(gp, paramRanges);;
    csv.write(profilesGP, "proflik.1d.gp.csv");
    csv.write(paramRanges, "proflik.1d.params.csv");
    
    //-------------------------------------------------------------------------
    // Data structures for Matlab plotting
    int N = 4;
    mat active_set(2*N, n_train);
    
    mat gpmv = zeros(2, n_points);
    gpmv.set_row(0, gpmean);
    gpmv.set_row(1, gpvar);
    
    mat test = zeros(2,n_points);
    test.set_row(0, Xtst);
    test.set_row(1, Ytst);
    
    mat obs = zeros(2,n_train);
    obs.set_row(0, Xtrn);
    obs.set_row(1, Ytrn);
        
    mat psgpm = zeros(N, n_points);
    mat psgpv = zeros(N, n_points);
    
    //-------------------------------------------------------------------------
    // Initialise the PSGP
    int n_active = 8;    // Start with 8 active points

    covFunc.setParameter(0, range);
    covFunc.setParameter(1, sill);
    covFunc.setParameter(2, nugget);
    PSGP psgp(Xtrnmat, Ytrn, covFunc, n_active, 1, 5);
    psgp.setGammaTolerance(1e-10);
    

    // Use a Gaussian likelihood model with fixed variance (set to 
    // the nugget variance) 
    GaussianLikelihood gaussLik(nugget);
    psgp.computePosterior(gaussLik);
    
    // Compute the PSGP posterior distribution under the Gaussian likelihood
    // for an increasing number of active points
    for (int i=0; i<4; i++) 
    {
        //---------------------------------------------------------------------
        cout << "PSGP with " << n_active << " active points" << endl;

        // Recompute the posterior - the best active points are selected from
        // the 40 points in the training set
        for (int i=0; i<3; i++) {
            SCGModelTrainer scg(psgp);
            scg.Train(5);
            psgp.resetPosterior();
            psgp.computePosterior(gaussLik);
        }
        
        // Note that we are not using the covariange function without the
        // nugget term for prediction.
        psgpmean = zeros(n_points);
        psgpvar  = zeros(n_points);
        psgp.makePredictions(psgpmean, psgpvar, mat(Xtst), gaussianCovFunc);
/*
 * 
        // Plot PSGP mean and error bars
        gplot.clearPlot();
        gplot.plotPoints(Xtst, psgpmean, "psgp mean", LINE, BLUE);
        gplot.plotPoints(Xtst, psgpmean + 2.0*sqrt(psgpvar), "error bar", LINE, CYAN);
        gplot.plotPoints(Xtst, psgpmean - 2.0*sqrt(psgpvar), "error bar", LINE, CYAN);

        // Plot full GP for comparison
        gplot.plotPoints(Xtst, gpmean, "psgp mean", LINE, GREEN);
        gplot.plotPoints(Xtst, gpmean + 2.0*sqrt(gpvar), "error bar", LINE, GREEN);
        gplot.plotPoints(Xtst, gpmean - 2.0*sqrt(gpvar), "error bar", LINE, GREEN);
        
        // Plot observations and true function
        gplot.plotPoints(Xtrn, Ytrn, "training set", CROSS, RED);  
        gplot.plotPoints(Xtst, Ytst, "true function", LINE, RED);

        // Plot active points
        vec activeX = (psgp.getActiveSetLocations()).get_col(0);
        vec activeY = Ytrn(psgp.getActiveSetIndices());
        gplot.plotPoints(activeX, activeY, "active points", CIRCLE, BLUE);
*/
        
        //---------------------------------------------------------------------
        // Compute likelihood profile
        mat profile = computeLikelihoodProfile(psgp, paramRanges);
        filename << "proflik.1d.psgp." << n_active << ".csv";
        csv.write(profile,filename.str());    // Save profile to CSV file
        filename.str("");                          // Reset filename
        
        //---------------------------------------------------------------------
        // Store data for Matlab plotting
        psgpm.set_row(i, psgpmean);
        psgpv.set_row(i, psgpvar);
        ivec iactive = psgp.getActiveSetIndices();
        cout << "Saving " << iactive.length() << " active points" << endl;
        active_set.set_row( i, concat( Xtrn(iactive), zeros(n_train-iactive.length()) ) );
        active_set.set_row( N+i, concat( Ytrn(iactive), zeros(n_train-iactive.length())) );
        
        //---------------------------------------------------------------------
        // Increase the number of active point for next iteration
        n_active *= 2;
        psgp.setActiveSetSize(n_active);

        cout << "Press a key to continue" << endl;
        // getchar();
    }

    //-------------------------------------------------------------------------
    // Write data for Matlab plotting
    csv.write(active_set,"demo_set_active_set.csv");
    csv.write(gpmv,  "demo_set_gp.csv");
    csv.write(psgpm, "demo_set_psgp_mean.csv");
    csv.write(psgpv, "demo_set_psgp_var.csv");
    csv.write(obs,   "demo_set_obs.csv");
    csv.write(test,  "demo_set_test.csv");
    
    return 0;   
}



/**
 * Compute likelihood profiles for PSGP for a given number
 * of active points
 */ 
mat computeLikelihoodProfile(PSGP &psgp, mat paramRanges)
{
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
