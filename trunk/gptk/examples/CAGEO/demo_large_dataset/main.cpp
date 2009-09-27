/**
 * This demonstration program illustrates the PSGP method in real
 * situation. Data from a simulated radioactive event is interpolated.
 * 
 * (c) 2009, Remi Barillec <r.barillec@aston.ac.uk>  
 **/

#include "main.h"

int main(void)
{
    int n_active = 400;
    time_t tstart = clock();
    run_scenario6(n_active);
    cout << (clock() - tstart)/1e6 << " seconds elapsed" << endl;
    
    return 0;
}    

/**
 * Run the Scenario 6 interpolation exercice with n_active active points
 */
void run_scenario6(int n_active)
{
    vec psgpmean, psgpvar;
    
    csvstream csv;                           // Input/output from/to CSV file
    
    
    //-------------------------------------------------------------------------
    // LOAD AND EXTRACT DATA
    mat data;
    string filename("INTAMP_scenario6_1080.csv");
    cout << "Loading data" << endl;
    if (csv.read(data, filename.c_str()))
    {
        cout << "Problem encountered while reading file " << filename << endl;
        cout << "Aborting." << endl;
        return;
    }
    
    // Use 50% of the data for training and 50% for testing (validation)
    cout << "Separating trainint/test set" << endl;
    int n_obs = data.rows();
    int n_train = ceil(0.8*n_obs);
    int n_test  = n_obs - n_train;
    
    ivec irand = itppext::randperm(n_obs);
    mat datarand = data.get_rows(irand);                    // Shuffle the data points
    
    mat data_trn = datarand.get_rows(0, n_train-1);
    mat data_tst = datarand.get_rows(n_train, n_obs-1);
    
    mat Xtrn = data_trn.get_cols(0, 1);         // Training inputs
    vec Ytrn = data_trn.get_col(2);             // Training outputs
    mat Xtst = data_tst.get_cols(0, 1);         // Testing inputs
    vec Ytst = data_tst.get_col(2);             // Testing outputs

    
    //-------------------------------------------------------------------------
    // COVARIANCE FUNCTION
    cout << "Finding sensible starting parameters." << endl;
    
    // A quick first guess for the covariance parameters
    double x1min = min(data.get_col(0));
    double x1max = max(data.get_col(0));
    double x2min = min(data.get_col(1));
    double x2max = max(data.get_col(1));
        
    double r1 = abs( x1max - x1min );
    double r2 = abs( x2max - x2min );
    double range  = 0.25 * ((r1+r2) / 2.0);
    double sill   = abs(variance(Ytrn));    
    double nugget = 0.01*sill;
    
    // Uniform grid for prediction
    int grid_size = 100;
    double dr = min(r1, r2)/grid_size;
    int n_grid1 = ceil(r1/dr);
    int n_grid2 = ceil(r2/dr);
    mat Xgrid(n_grid1*n_grid2, 2);
    vec x1grid = linspace(x1min, x1max, n_grid1);
    vec x2grid = linspace(x2min, x2max, n_grid2);

    int k=0;
    for (int i=0; i<n_grid1; i++) {
        for (int j=0; j<n_grid2; j++) {
            Xgrid(k,0) = x1grid(i);
            Xgrid(k,1) = x2grid(j);
            k++;
        }
    }
    cout << "Making predictions on gridded data " << x1min << ":" << x1max << " x " << x2min << ":" << x2max << "(" << n_grid1 << "x" << n_grid2 << ")" << endl;

    
    // We use a mixture of several kernels for more flexibility
    cout << "Initialising covariance function" << endl;
    GaussianCF  kernel1(range, sill);    
    Matern3CF   kernel2(range, sill);    
    Matern5CF   kernel3(range, sill);    
    NeuralNetCF kernel4(range, sill);    
    // ConstantCF  bias(1e3);
    
    SumCovarianceFunction kernelCovFunc(kernel1);           // Put together the sum of various kernels
    kernelCovFunc.addCovarianceFunction(kernel2);
    kernelCovFunc.addCovarianceFunction(kernel3);
    kernelCovFunc.addCovarianceFunction(kernel4);
    // kernelCovFunc.addCovarianceFunction(bias);
    
    WhiteNoiseCF noise(nugget);                             // And add some white noise to make up the final
    SumCovarianceFunction covFunc(kernelCovFunc);           // covariance function
    covFunc.addCovarianceFunction(noise);
        
    covFunc.displayCovarianceParameters();
    
    //-------------------------------------------------------------------------
    // Initialise the PSGP
    cout << "Initialising PSGP with " << n_active << " active points" << endl;
    PSGP psgp(Xtrn, Ytrn, covFunc, n_active);
    
    // Use a Gaussian likelihood model with fixed variance (set to 
    // a small percentage of the nugget variance) 
    GaussianLikelihood gaussLik(nugget);
    
    cout << "Compute posterior under Gaussian likelihood" << endl;
    psgp.computePosterior(gaussLik);
    
    cout << "PSGP retained " << psgp.getSizeActiveSet() << " active sites" << endl;
    
    // Estimate parameters
    SCGModelTrainer scg(psgp);
    
    for (int i=0; i<3; i++) 
    {
        scg.Train(5);
        psgp.resetPosterior();
        psgp.computePosterior(gaussLik);
    }
    
    covFunc.displayCovarianceParameters();
    
    // Prediction on the test set
    cout << "Making predictions on test data (" << n_test << " locations)" << endl;
    vec tstmean = zeros(n_test);
    vec tstvar  = zeros(n_test);
    psgp.makePredictions(tstmean, tstvar, Xtst, kernelCovFunc);
        
    // Prediction on a uniform grid (square grid cells)
    vec gridmean = zeros(n_grid1*n_grid2);
    vec gridvar  = zeros(n_grid1*n_grid2);
    
    cout << "Making predictions on gridded data " << x1min << ":" << x1max << " x " << x2min << ":" << x2max << "(" << n_grid1 << "x" << n_grid2 << ")" << endl;
    psgp.makePredictions(gridmean, gridvar, Xgrid, kernelCovFunc);
    
    
    //---------------------------------------------------------------------
    // Store data for Matlab plotting
    //-------------------------------------------------------------------------
    cout << "Saving results" << endl;
    // PSGP results on grid
    mat psgpgrid = zeros(4, n_grid1*n_grid2);
    psgpgrid.set_row(0, Xgrid.get_col(0));
    psgpgrid.set_row(1, Xgrid.get_col(1));
    psgpgrid.set_row(2, gridmean);
    psgpgrid.set_row(3, gridvar);
    csv.write(psgpgrid, "demo_scenario6_psgp_grid.csv");
    
    // PSGP results on test set
    mat psgptest = zeros(5, n_test);
    psgpgrid.set_row(0, Xtst.get_col(0));
    psgpgrid.set_row(1, Xtst.get_col(1));
    psgptest.set_row(2, tstmean);
    psgptest.set_row(3, tstvar);
    psgptest.set_row(4, Ytst);
    csv.write(psgptest, "demo_scenario6_psgp_test.csv");
    
    cout << "Done" << endl;
}


