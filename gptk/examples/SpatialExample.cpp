#include "SpatialExample.h"

SpatialExample::SpatialExample()
{
    setDefaults();
}

SpatialExample::SpatialExample(bool _normaliseData)
{
    setDefaults();
    normaliseData = _normaliseData;
}

SpatialExample::~SpatialExample()
{
}

void SpatialExample::setDefaults()
{
    observedNoise = false;            // No estimate for observation noise
    normaliseData = true;             // Normalise the data by default
    isSetPredictionLocations = false; // Prediction locations haven't been set
    
    defaultNuggetRatio = 0.1;       // Nugget to Variance (sill) ratio
    
    paramEstimationMethod = PARAM_ESTIM_GP;   // Use GP to estimate parameters
    n_optim_iterations = 50;        // Number of iterations in parameter estimation
    
    n_active = 400;           // Number of active points (max)
    n_sweeps = 4;             // Number of sweeps through data (PSGP)
    n_outer_loops = 5;        // Number of outer loops (PSGP parameters estimation)
    
    predictionType = PREDICTION_CHUNKS;  // Split prediction domain by default
    predictionChunkSize = 1000;         // Default size for prediction chunks
}


/**
 * Set the maximum number of active points
 */
void SpatialExample::setNumberActivePoints(int value)
{
    n_active = value;
}

/**
 * Set the number of iterations in parameter optimisation
 */
void SpatialExample::setNumberOptimIterations(int value)
{
    n_optim_iterations = value;
}

/**
 * Set the number of outer loops in PSGP parameter estimation
 */
void SpatialExample::setNumberOuterLoops(int value)
{
    n_outer_loops = value;
}

/**
 * Set number of (data recycling) sweeps through data in PSGP
 */
void SpatialExample::setNumberSweeps(int value)
{
    assert(value >= 0);
    n_sweeps = value;
}


/**
 * Set the default nugget (proportion of process variance)
 */ 
void SpatialExample::setDefaultNuggetRatio(double value) 
{
    defaultNuggetRatio = value;
}

/**
 * Set prediction locations
 */
void SpatialExample::setPredictionLocations(mat _Xpred) 
{
    if (Xpred.rows() > 0 && Xpred.cols() == 2) {
        Xpred = _Xpred;
        isSetPredictionLocations = true;
    }
    else {
        cerr << "Error setting in setPredictionLocations(): invalid locations." << endl;
    }
}


/**
 * 
 */
void SpatialExample::setParameterEstimationMethod(ParameterEstimationMethod method)
{
    paramEstimationMethod = method;
}

/**
 * Set prediction type: 
 *   PREDICTION_FULL: single prediction at all locations
 *   PREDICTION_CHUNKS: sequential prediction at subsets of locations  
 */
void SpatialExample::setPredictionType(PredictionType type) 
{
    predictionType = type; 
}

/**
 * Load data from a CSV matrix file
 * Data file must have 3 or 4 columns, understood as:
 * location (1st-2nd cols), observation (3rd col.) and if provided, 
 * observation noise (4th col.)
 * 
 * This sets observedNoise to true if noise information is available
 * (i.e. the 4th column is provided) and to false otherwise.
 * 
 * Returns true if there is a problem loading the data, false
 * otherwise.
 */
bool SpatialExample::loadData(string filename)
{
    csvstream csv;
    mat M;
    
    // Read data - return error if problem
    if (csv.read(M,filename)) {
        cerr << "Error reading data file" << endl;
        return true;
    }
    
    // M.set_size(600,M.cols(),true);
    
    if (M.cols() < 3 || M.cols() > 4) {
        cerr << "Invalid number of columns. Data file should have 3 or 4 columns." << endl;
        return true;
    }
    
    // Normalise (if required)
    if (normaliseData) {
        itppext::normalise(M, Mmean, Mcovdiag);
    }
    
    // Extract locations and observations 
    X = zeros(M.rows(), 2);
    y = zeros(M.rows(), 1);
    
    X.set_col(0,M.get_col(0));
    X.set_col(1,M.get_col(1));
    y = M.get_col(2);
    
    // Observation noise 
    if (M.cols() == 4) {            // Observation noise provided
        v = M.get_col(3);
        observedNoise = true;
    } 
    else {                          // No observation noise 
        observedNoise = false;
    }
    
    return false;
}


/**
 * Pre-process the data
 * 
 * This does nothing but provide a means to do some preprocessing if
 * needed (simply overload).
 * 
 * Returns true if an error occured, false otherwise.
 */ 
bool SpatialExample::preProcess() 
{
    return false;
}

/**
 * Computes and sets the initial parameters for the PSGP
 */
bool SpatialExample::initParameters() 
{
    if (normaliseData)         // Sphered data - set parameters to sensible values 
    {    
        range = 0.5;
        sill  = 1.0;
    }
    else                       // Raw data - use simple heuristic
    {
        double r1 = abs(max(X.get_col(0)) - min(X.get_col(0)));
        double r2 = abs(max(X.get_col(1)) - min(X.get_col(1)));
        range  = 0.25 * ((r1+r2) / 2.0);
        sill   = abs(variance(y));
    }
    
    // Default nugget set to proportion of of process variance (default is 10%) 
    nugget = defaultNuggetRatio*sill;
    
    return false;
}


/**
 * Initialises the covariance function and sets its parameters
 * to that determined above
 */ 
bool SpatialExample::initCovarianceFunction()
{
    kernelCF = new ExponentialCF(range, sill);
    
    if (observedNoise) {        // Use fixed diagonal covariance matrix
        // nuggetCF = new DiagonalMatrixCF(v);
        nuggetCF = new WhiteNoiseCF(nugget);
    } 
    else {                      // Use white noise covariance function
        nuggetCF = new WhiteNoiseCF(nugget); 
    }
        
    covFunc = new SumCovarianceFunction(*kernelCF);
    ((SumCovarianceFunction*) covFunc)->addCovarianceFunction(*nuggetCF);
    
    return false;
}


/**
 * Learn PSGP parameters
 * 
 * 4 methods are provided:
 * - Use a GP to estimate the parameters (learnParametersGP)
 * - Use the PSGP to estimate the parameters (learnParametersPSGP)
 * - Use a custom parameter estimation method (learnParametersCustom)
 * - Do not optimise parameters
 */
bool SpatialExample::learnParameters() 
{
    switch(paramEstimationMethod) 
    {
    case PARAM_ESTIM_GP:
        // Estimate parameters using GP
        return learnParametersGP();
    
    case PARAM_ESTIM_PSGP:
        // Estimate parameters using PSGP
        return learnParametersPSGP();
    
    case PARAM_ESTIM_CUSTOM:
        // Estimate parameters using custom method
        // Requires overloading of learnParametersCustom()
        return learnParametersCustom();
        
    case PARAM_ESTIM_NO_ESTIMATION:
        // Do not estimate parameters
        cout << "No parameter estimation was required. Using initial parameter values." << endl;
        break;
    
    default:
        cerr << "Unknown parameter estimation mode." << endl;
    }
   
    return false;
}

/**
 * Use a standard GP to estimate parameters
 */
bool SpatialExample::learnParametersGP()
{
    GaussianProcess gp(2, 1, X, y, *covFunc);

    SCGModelTrainer gpTrainer(gp);

    gpTrainer.setAnalyticGradients(true);
    gpTrainer.setCheckGradient(true);
    gpTrainer.Train(n_optim_iterations);
    
    return false;
}

/**
 * Use the PSGP to estimate parameters
 * 
 * This relies on several parameters such as number of inner/outer loops,
 * etc. which can be set using the relevant modifier methods.
 * 
 * Note that this optimisation method can be numerically unstable. 
 */
bool SpatialExample::learnParametersPSGP() 
{
    // Gaussian Likelihood function
    likFunction = new GaussianLikelihood(nugget);
        
    // PSGP for prediction
    n_active = min(n_active, X.rows()); // Do not exceed number of obs.
        
    SequentialGP psgp(2, 1, n_active, X, y, *covFunc, n_sweeps);
        
    cout << "  Compute posterior" << endl;
    psgp.computePosterior(*likFunction);
    
    SCGModelTrainer gpTrainer(psgp);
    psgp.setLikelihoodType(UpperBound);

    gpTrainer.setCheckGradient(true);

    for (int i=0; i<n_outer_loops; i++) 
    {
        cout << endl << endl << "-- " << i+1 << "/" << n_outer_loops << endl;
        gpTrainer.Train(n_optim_iterations);
        psgp.resetPosterior();
        psgp.computePosterior(*likFunction);
    }
    return false;
}

/**
 * To use a custom parameter estimation method, override the following.
 * This should not be called and will issue a warning message. 
 */
bool SpatialExample::learnParametersCustom()
{
    cout << "**WARNING: You need to override SpatialExample::learnParametersCustom()." << endl;
    cout << "           SpatialExample::learnParametersCustom() is a template method and" << endl; 
    cout << "           does nothing as such." << endl;
    return false;
}


void SpatialExample::displayCurrentOptions()
{
    cout << endl << "Current options:" << endl;
    cout << "  Using normalised data: " << normaliseData << endl;
    cout << "  Observation noise provided: " << observedNoise << endl;
    cout << "  Using " << n_active << " active points / " << X.rows() << " observations" << endl;
    covFunc->displayCovarianceParameters(2);
    cout << endl;
}

bool SpatialExample::makePredictions()
{
    // If no prediction locations have been provided,
    // default to uniform grid on observed area.
    if (!isSetPredictionLocations) uniformGrid(X, Xpred);
        
    // Reinitialise predictive mean and variance
    ypred = zeros(Xpred.rows());
    vpred = zeros(Xpred.rows());
    
    // Gaussian Likelihood function
    // TODO: Make this generic so can be set to something else
    likFunction = new GaussianLikelihood(nugget);
    
    // PSGP for prediction
    n_active = min(n_active, X.rows()); // Do not exceed number of obs.
    
    SequentialGP psgp(2, 1, n_active, X, y, *covFunc, n_sweeps);
    
    cout << "  Compute posterior" << endl;
    psgp.computePosterior(*likFunction);
    
    bool predictionError = false;
    
    switch (predictionType) {
    case PREDICTION_FULL:
        predictionError = makePredictionsFull(psgp);
        break;
        
    case PREDICTION_CHUNKS:
        predictionError = makePredictionsChunks(psgp);
        break;
        
    default: 
        cerr << "Unknown prediction type" << endl;
        return true;
    }
    
    // Make sure everything went fine
    if (predictionError) {
        cerr << "An error occurred during prediction." << endl;
        cerr << "Results might not be correct." << endl;
        return true;
    }
    
    return false;
}


/**
 * Predict at all specified locations at once
 */
bool SpatialExample::makePredictionsFull(SequentialGP &psgp) 
{
    cout << "  Predict at " << Xpred.rows() << " locations" << endl;
    psgp.makePredictions(ypred, vpred, Xpred, *kernelCF);

    return false;
}


/**
 * Limited memory prediciton - the prediction domain is split
 * into chunks of predefined size and prediction is made at each
 * chunk sequentially.
 */
bool SpatialExample::makePredictionsChunks(SequentialGP &psgp) 
{
    int chunkStart = 0;                        
    int chunkSize  = predictionChunkSize;
    int chunkEnd   = chunkSize - 1;

    int numPred = Xpred.rows();
    
    // Adjust chunk size if default chunk bigger than preduction domain
    if(chunkEnd >= numPred)
    {
        chunkEnd = numPred - 1;
        chunkSize = chunkEnd - chunkStart + 1;
    }

    // Predict at each chunk
    while(chunkStart < numPred)
    {
        cout << "  Predicting from locations " << chunkStart+1 << " to " << chunkEnd+1 << endl;

        // Extract chunk
        mat XpredChunk = Xpred.get_rows(chunkStart, chunkEnd);

        // Predict at chunk
        vec predYChunk(chunkSize);
        vec predVarChunk(chunkSize);

        psgp.makePredictions(predYChunk, predVarChunk, XpredChunk, *kernelCF);

        ypred.replace_mid(chunkStart, predYChunk);
        vpred.replace_mid(chunkStart, predVarChunk);

        // Move to next chunk
        chunkStart = chunkEnd + 1;
        chunkEnd +=  chunkSize;

        // Adjust chunk size if taking us beyond the total number of predictive 
        // locations (typically, when last chunk is smaller than the others)
        if(chunkEnd >= numPred)
        {
            chunkEnd = numPred - 1;
            chunkSize = chunkEnd - chunkStart + 1;
        }
    }
    
    return false;
}

/**
 * Given a set of locations X, determine a rectangle uniform grid covering
 * all locations in X and having resolution equal to the minimum distance
 * between any 2 observations. The resulting grid is returned as a set of 
 * locations (2-column matrix).
 * If the resolution is 0 or ranges are empty, the empty matrix is returned.
 */
void SpatialExample::uniformGrid(mat X, mat &grid)
{
    // Determine max resolution (min distance between any 2 observations)
    double minDist2 = pow(X(0,0)-X(1,0),2) + pow(X(0,1)-X(1,1),2);
    double x1min, x1max, x2min, x2max;
    
    // Initialise min and max along each coordinate
    x1min = X(0,0);
    x1max = X(0,0);
    x2min = X(0,1);
    x2max = X(0,1);
        
    // Compute minimum distance between any 2 points
    for (int i=0; i<X.rows(); i++) {
        vec x = X.get_row(i);
        for (int j=i+1; j<X.rows(); j++) {
            vec x2 = X.get_row(j);
            double dist2 = pow(x(0)-x2(0),2) + pow(x(1)-x2(1),2);
            
            // Update minimum squared distance
            if (dist2 < minDist2) minDist2 = dist2;
        }
        
        // Update xmin, xmax
        if (x(0) < x1min) x1min = x(0);
        if (x(0) > x1max) x1max = x(0);
        if (x(1) < x2min) x2min = x(1);
        if (x(1) > x2max) x2max = x(1);
    }

    double res = sqrt(minDist2);
    
    // Make sure resolution is positive and greater than machine precision
    if (0.5*res == res || abs(res) != res)
    {
        cerr << "Error computing uniform grid: resolution below machine precision" << endl;
        grid = mat(0,0);
    }
    
    vec X1 = linspace(x1min, x1max, (int) ((x1max-x1min)/res) );
    vec X2 = linspace(x2min, x2max, (int) ((x2max-x2min)/res) );
            
    cout << "  Grid size: " << X1.length() << "x" << X2.length() << "x" << res << endl;
    
    // Make sure ranges are not empty
    if (X1.length() == 0 || X2.length() == 0) {
        cerr << "Error computing uniform grid: empty range" << endl;
        grid = mat(0,0);
    }
    
    // Compute grid locations
    grid = zeros(X1.length()*X2.length(),2);
    for (int i=0; i<X1.length(); i++) {
        for (int j=0; j<X2.length(); j++) {
            grid(i*X2.length()+j,0) = X1(i);
            grid(i*X2.length()+j,1) = X2(j);
        }
    }
   
}



/*
 * Save prediction results to file
 */
bool SpatialExample::saveResults(string filename)
{
    csvstream csv;
    mat M = zeros(Xpred.rows(), 3);
    
    M.set_cols(0,Xpred);             // Prediction locations
    M.set_col(2, ypred);             // Prediction mean 
    
    // If normalised data, project back to initial domain
    if (normaliseData) {
        itppext::denormalise(M, Mmean, Mcovdiag);
        vpred = Mcovdiag(2)*vpred;
    }
    
    M.set_size(M.rows(),M.cols()+1,true);
    M.set_col(M.cols()-1,vpred);
    
    if (csv.write(M, filename)) {
        cerr << "Could not write prediction data to file." << endl;
        return true;
    }
    
    return false;
}


/**
 * Wrapper to run all operations (default mode, can be overridden 
 * by calling the inner methods separately)
 */
void SpatialExample::run(string datafile, string predfile)
{
    // Load data file
    cout << "Loading data file " << datafile << endl;
    loadData(datafile);

    // Set initial parameters
    cout << "Set/Compute initial parameters" << endl;
    initParameters();
    
    // Set covariance function
    cout << "Set covariance function" << endl;
    initCovarianceFunction();
    displayCurrentOptions();
    
    // Learn parameters
    cout << "Learning covariance function parameters" << endl;
    learnParameters();
    displayCurrentOptions();
    
    // Make predictions
    cout << "Make predictions" << endl;
    makePredictions();
    
    // Save results
    cout << "Saving results" << endl;
    saveResults(predfile);
    
}

int main(int argc, char* argv[])
{
    SpatialExample ex;
    
    switch (argc) {
    case 3: 
        {
            ex.setParameterEstimationMethod(PARAM_ESTIM_PSGP);
            ex.setNumberSweeps(1);
            ex.setNumberOuterLoops(1);
            ex.setNumberOptimIterations(5);
            
            long starttime = clock();
            ex.run(argv[1], argv[2]);
            long endtime = clock();
            cout << (endtime - starttime)/1e6 << " seconds elapsed" << endl;
        }
        break;
    
    default:
        cerr << "Usage: " << endl;
        cerr << "  demoName datafile predfile" << endl;
        return 1;
    }
    
    return 0;
}







