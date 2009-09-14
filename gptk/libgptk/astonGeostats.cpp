#include "astonGeostats.h"
#include <iostream>

/**
 * Function to convert metadata indices (1 to N) to C indices (0 to N-1)
 */
int minusOne(int n) { return n-1; }

void learnParameters(int numObs, int numError, double *xData, double *yData, double *errorData, 
        int numMetadata, int *errPtr, int *sensorPtr, char **metaDataTable, double *range, double *sill, 
        double *nugget, double *bias, int model)
{

	vec rndNums = itpp::randu(numObs);
	ivec rndIdx = itpp::sort_index(rndNums);
	
	// Limit to MAX_OBS observation for parameter estimation
	// so that computations remain timely
	int minLim = min(numObs, MAX_OBS);

	mat X = itpp::zeros(minLim, 2);
	vec y = itpp::zeros(minLim);
	vec E = itpp::zeros(minLim);
	Vec<LikelihoodType *> noiseModels;

	// need to make this use a factory method from IT++
	// this would make things much more efficient
	for(int i=0; i < minLim; i++) {
		X(i, 0) = xData[rndIdx[i]];
		X(i, 1) = xData[rndIdx[i] + numObs];
		y(i) = yData[rndIdx[i]];
		E(i) = errorData[rndIdx[i]];
	}

	// some reasonable starting parameters
	double r1=abs(min(max(X, 1) - min(X, 1)));
	double r2=abs(min(max(X, 2) - min(X, 2)));

	*range = 0.25 * ((r1+r2) / 2.0);
	*sill = abs(variance(y));
	*nugget = 0.5 * *sill;
	*bias = abs(1.0 / mean(y));

	
	// hardcode a gaussian covariance function so far
	ExponentialCF covComp1(*range, *sill);
  	WhiteNoiseCF covComp2(*nugget);
	ConstantCF covComp3(*bias);

	SumCovarianceFunction covSum(covComp1);
	covSum.addCovarianceFunction(covComp2);
	covSum.addCovarianceFunction(covComp3);

	cout << "===Starting Parameters===========================" <<endl;
	covSum.displayCovarianceParameters();

	// there will be a little memory leak here - should sort this out at somepoint...
	covSum.setTransform(0, new LogTransform());
	covSum.setTransform(1, new LogTransform());
	covSum.setTransform(2, new LogTransform());
	covSum.setTransform(3, new LogTransform());

	if(PARAMETER_ESTIMATION_USING_GP) 
	{
	    // Use a GP to learn parameters
	    GaussianProcess gp(2, 1, X, y, covSum);
	
    	SCGModelTrainer gpTrainer(gp);
    
    	gpTrainer.setAnalyticGradients(true);
    	gpTrainer.setCheckGradient(true);
    	gpTrainer.Train(50);
	}
	else
	{
	    GaussianLikelihood* defaultLikelihood;
	    Vec<LikelihoodType*> likelihoodModels(numObs);        // Vector of likelihood models
	    ivec iLikelihoodModel(numObs);
	    
	    // PARSE METADATA IF AVAILABLE, OTHERWISE USE DEFAULT LIKELIHOOD MODEL
	    if(numMetadata == 0)
	    {
	        cout << "No noise model specified" << endl;
	        cout << "Defaulting to GAUSSIAN:" << (*nugget * 0.01) << endl;
	        defaultLikelihood = new GaussianLikelihood(*nugget * 0.01);
	    }
	    else
	    {
	        cout << "Noise models specified. Extracting from metadata table." << endl;

	        string    modelName[numObs];                          // Name of likelihood model for each obs
	        itpp::vec modelParams[numObs];                        // Params of likelihood model for each obs
	        
	        // Shift indexes (starting from 1 in data, from 0 in ITPP)
	        iLikelihoodModel = apply_function(minusOne, ivec(errPtr, numObs));

	        // Make sure we have same number of locations and error models
	        assert(iLikelihoodModel.length() == numMetadata);
	        assert(iLikelihoodModel.length() == X.rows());

	        // Extract the noise models and their parameters for all observations
	        parseMetadata(metaDataTable, iLikelihoodModel, modelName, modelParams);

	        // Build vector of likelihod models 
	        buildLikelihoodVector(modelName, modelParams, iLikelihoodModel, likelihoodModels);
	    }
	    
	    // LEARN PARAMETERS USING PSGP AND LIKELIHOOD MODEL(S) GENERATED ABOVE
	    PSGP psgp(X, y, covSum);
	            
	    if (numMetadata == 0) {
	        psgp.computePosterior(*defaultLikelihood);
	    }
	    else {
	        psgp.computePosterior(iLikelihoodModel, likelihoodModels);
	    }
	    
	    SCGModelTrainer gpTrainer(psgp);
	    gpTrainer.setAnalyticGradients(true);
	    gpTrainer.setCheckGradient(true);
	    
	    for (int i=0; i<PSGP_PARAM_ITERATIONS; i++)
	    {
	        gpTrainer.Train(PSGP_SCG_ITERATIONS);
	        // psgp.resetPosterior();
	        // psgp.computePosterior(*likFunction);
	        psgp.recomputePosterior();
	    }
	}
	
	
	
	cout << "===Finishing Parameters===========================" <<endl;
	covSum.displayCovarianceParameters();

	*range  = covSum.getParameter(0);
	*sill   = covSum.getParameter(1);
	*nugget = covSum.getParameter(2);
	*bias   = covSum.getParameter(3);
}





/**
 * Prediction using Projected Sequential Gaussian Processes at a set of specified locations
 * 
 * int numObs               Number of observations
 * int numPred              Number of predictions
 * int numError             Number of errors
 * double *xData            Array of observed locations
 * double *yData            Array of observed values
 * double *xPred            Array of prediction locations
 * double *errorData        Array of error values
 * int numMetadata          Size of metadata table (i.e. number of likelihood models)
 * int *errPtr              Array of likelihood model indexes (e.g. 1 per observation)
 * int *sensorPtr           Array of sensor indexes (e.g. 1 per observation)
 * char **metaDataTable,    Metadata array of strings, each string containing the description of a 
 *                          likelihood model in the form "<DISTRIBUTION_NAME>,<PARAM_1>,...,<PARAM_N>" 
 *
 * double *meanPred,        Array of prediction mean values                
 * double *varPred,         Array of prediction variances
 * 
 * double *range,           PSGP range parameter                
 * double *sill,            PSGP sill parameter
 * double *nugget,          PSGP nugget parameter
 * double *bias,            PSGP bias parameter 
 * int model                ? (Unused)
 */
void makePredictions(int numObs, int numPred, int numError, double *xData, double *yData,  double *xPred, 
                     double *errorData, int numMetadata, int *errPtr, int *sensorPtr, char **metaDataTable, 
                     double *meanPred, double *varPred, double *range, double *sill, double *nugget,
                     double *bias, int model)
{
    
	LikelihoodType *defaultLikelihood = new GaussianLikelihood(*nugget * 0.01);
	
	Vec<LikelihoodType*> likelihoodModels(numObs);        // Vector of likelihood models
	ivec iLikelihoodModel(numObs);                        // Index of likelihood model for each obs
	
	cout<< "Number of observations: " << numObs << endl;
	cout<< "Number of prediction:   " << numPred << endl;
	cout<< "Number of error values: " << numError << endl;
	cout<< "Size of metadata:       " << numMetadata << endl;
	
	// Convert arguments to IT++ vectors/matrices
	vec rndNums = itpp::randu(numObs);
	ivec rndIdx = itpp::sort_index(rndNums);

	// Why only use 500 obs? Speed issue?
	// int minLim = min(numObs, 500);

	// mat X = itpp::zeros(numObs,2);
	// vec y = itpp::zeros(numObs);
	// vec E = itpp::zeros(numObs);
	mat X = reshape(vec(xData,2*numObs),numObs,2);
	vec y = vec(yData,numObs);
	vec E = vec(errorData,numObs);
	
	// mat Xpred = itpp::zeros(numPred, 2);
	mat Xpred = reshape(vec(xPred,2*numPred),numPred,2);
	vec predY = itpp::zeros(numPred);
	vec predVar = itpp::zeros(numPred);
	
	/*
	// DEBUG: Dump data to file(s)
	csvstream csv;
	mat debug_data(X);
	debug_data.append_col(y);
	debug_data.append_col(E);
	debug_data.append_col(to_vec(ivec(errPtr,numObs)));
	debug_data.append_col(to_vec(ivec(sensorPtr,numObs)));
	csv.write(debug_data, "debug_data.csv");            // Input data
	csv.write(Xpred, "debug_pred.csv");        // Prediction locations
	vec debug_params(9);
	debug_params(0) = numObs;
	debug_params(1) = numPred;
	debug_params(2) = numError;
	debug_params(3) = numMetadata;
	debug_params(4) = *range;
	debug_params(5) = *sill; 
	debug_params(6) = *nugget;
	debug_params(7) = *bias;
	debug_params(8) = model;
	csv.write(debug_params,"debug_params.csv");

	ofstream f("debug_metadata.csv", ios::trunc); 
	for(int i=0; i<numMetadata; i++) {
	    f << metaDataTable[i] << endl;
	}
	f.close();
	*/
	
	// need to make this use a factory method from IT++
	// this would make things much more efficient
	/*
	for(int i=0; i < numObs; i++) {
	    X(i, 0) = xData[rndIdx[i]];
	    X(i, 1) = xData[rndIdx[i] + numObs];
	    // y(i) = yData[rndIdx[i]];
	    // E(i) = errorData[rndIdx[i]];
	}
	
	for(int i=0; i < numPred; i++) {
            Xpred(i, 0) = xPred[i];
            Xpred(i, 1) = xPred[i + numPred];
    }
	*/
	
	//-------------------------------------------------------------------------
	//
	// PARSE METADATA IF AVAILABLE, OTHERWISE USE DEFAULT LIKELIHOOD MODEL
	//
	//-------------------------------------------------------------------------
	if(numMetadata == 0)
	{
		cout << "No noise model specified" << endl;
		cout << "Defaulting to GAUSSIAN:" << (*nugget * 0.01) << endl;
		defaultLikelihood = new GaussianLikelihood(*nugget * 0.01);
	}
	else
	{
	    cout << "Noise models specified. Extracting from metadata table." << endl;
	 
	    string    modelName[numObs];                          // Name of likelihood model for each obs
	    itpp::vec modelParams[numObs];                        // Params of likelihood model for each obs
	    
	    // Shift indexes (starting from 1 in data, from 0 in ITPP)
	    iLikelihoodModel = apply_function(minusOne, ivec(errPtr, numObs));
	    
	    // Make sure we have same number of locations and error models
	    // assert(iLikelihoodModel.length() == numMetadata);
	    // assert(iLikelihoodModel.length() == X.rows());
	    
	    // Extract the noise models and their parameters for all observations
	    parseMetadata(metaDataTable, iLikelihoodModel, modelName, modelParams);
	    
	    // Build vector of likelihod models 
	    buildLikelihoodVector(modelName, modelParams, iLikelihoodModel, likelihoodModels);
	}
	
	
	//-------------------------------------------------------------------------
	//
	// INITIALISE PSGP AND COVARIANCE FUNCTION
	//
	//-------------------------------------------------------------------------
	// Hardcode a gaussian covariance function so far
	ExponentialCF covComp1(*range, *sill);
  	WhiteNoiseCF covComp2(*nugget);
	ConstantCF covComp3(*bias);

	// PSGP Covariance function
	SumCovarianceFunction covSum(covComp1);
	covSum.addCovarianceFunction(covComp2);
	covSum.addCovarianceFunction(covComp3);

	// PSGP Covariance function for prediction (same as above 
	// but without nugget term - noise free prediction)
	SumCovarianceFunction covSumPred(covComp1);
	covSumPred.addCovarianceFunction(covComp3);
	
	// Set number of active points to max(MAX_ACTIVE_POINTS, size(obs))
	int n_active = min(MAX_ACTIVE_POINTS, numObs);
	
	// SequentialGP ssgp(2,1,n_active,X,y,covSum,1);
	PSGP ssgp(X,y,covSum,n_active,1,1);

	cout << "Computing posterior..." << endl;
	if (numMetadata == 0) {
	    ssgp.computePosterior(*defaultLikelihood);
	}
	else {
	    ssgp.computePosterior(iLikelihoodModel, likelihoodModels);
	}
	cout << "Predicting..."<<endl;

	// should add an option here to select one or the other...
	if(0)
	{
	    ssgp.makePredictions(predY, predVar, Xpred, covSumPred);
	}
	else
	{
	    int startVal = 0;
	    int chunkSize = 1000;
	    int endVal = chunkSize - 1;

	    if(endVal > numPred)
	    {
	        endVal = numPred - 1;
	        chunkSize = endVal - startVal + 1;
	    }

	    while(startVal < numPred)
	    {
	        cout << "  Predicting chunk [" << startVal << ":" << endVal << "/" << numPred << "]" << endl;
	        
	        mat XpredChunk = Xpred.get_rows(startVal, endVal);

	        // predYChunk.set_size(chunkSize);
	        // predVarChunk.set_size(chunkSize);
	        vec predYChunk(chunkSize);
	        vec predVarChunk(chunkSize);

	        ssgp.makePredictions(predYChunk, predVarChunk, XpredChunk, covSumPred);

	        predY.replace_mid(startVal, predYChunk);
	        predVar.replace_mid(startVal, predVarChunk);
	        
	        startVal = endVal + 1;
	        endVal = endVal + chunkSize;
	        
	        if(endVal >= numPred)
	        {
	            endVal = numPred - 1;
	            chunkSize = endVal - startVal + 1;
	        }
	    }
	}

	// should use an IT++ factory for vectors
	for(int i=0; i < numPred; i++)
	{
		meanPred[i] = predY(i);
		varPred[i]  = predVar(i);
	}

	covSum.displayCovarianceParameters();

}


/**
 * Parse the metadata table.
 * 
 * PARAMETERS:
 * 
 *  metadataTable:  a array of C strings, each element being a textual
 *                  representation of a distribution (noise model) and
 *                  its parameters (e.g. "GAUSSIAN,0.0,1.0" for the normal
 *                  distribution)
 *  
 *  modelIndex:     a vector of indexes giving, for each observation, 
 *                  the index of the row in the metadata table that
 *                  contains the corresponding noise model data.
 * 
 * RETURNS:
 * 
 *  vecModelName:   a std::vector of distribution names (strings). These
 *                  correspond to the first part of the metadata (e.g.
 *                  "GAUSSIAN")
 * 
 *  vecModelParams: a std::vector of distribution parameters (strings).
 *                  These correspond to the parameters (as text) of the 
 *                  distribution (e.g. "0.0,1.0")
 * 
 *  modelIndex:     Observations without a valid noise model have their index
 *                  set to -1 in this vector, for further processing in the
 *                  calling function. 
 */  
void parseMetadata(char** metadataTable, const ivec modelIndex,  
                   string modelNames[], vec modelParams[])
{
    // For each observation
    for(int iObs = 0; iObs < modelIndex.length(); iObs++)
    {
        // modelSet(iObs) = false;
        
        // Get index of corresponding noise model
        int iNoiseModel = modelIndex(iObs);
        
        // Retrieve corresponding string in metadata table
        // Should be something like "GAUSSIAN,0.0,0.123", i.e. a list of string
        // tokens where the first is the distribution identifier and the following
        // are the values for the parameters (number of these depends on the distribution
        // type, e.g. Gaussian has 2 parameters: mean and variance)
        // cout << "metadataTable[" << iNoiseModel << "]: " << metadataTable[iNoiseModel] << endl;
        string strNoiseModel(metadataTable[iNoiseModel]);
        
        // Split string in 2 w.r.t. to first (non-trailing) comma delimiter
        std::vector<string> tokens;
        itppext::tokenise(strNoiseModel, tokens, ", ", 1);

        // Retrieve distribution name
        modelNames[iObs] = tokens[0];
                
        // Use IT++ vector initialisation from string - if this does not work, 
        // i.e. string cannot be converted to double, flag the noise model as
        // invalid by setting the modelIndex to -1
        
        vec distParams;
        
        try {
            modelParams[iObs] = vec(tokens[1]);
        }
        catch (std::runtime_error e) {
            cout << "** Error in metadata parsing for observation " << iObs << ":" << endl;
            cout << "   Invalid parameter string \"" << tokens[1] << "\"" << endl;
            cout << "   Parameter string must be a sequence of numeric values" << endl;
            cout << "   separated by commas, e.g. \"1.23,4,5.6,78.9\"" << endl;
         
            // Flag observation as not having a valid model
            modelNames[iObs] = INVALID_MODEL_NAME;
        }
    }
    
}

/**
 * Return a vector of likelihood models based on array of model names
 * and array of model parameters.
 */
void buildLikelihoodVector(string modelName[], vec modelParams[], ivec iLikelihoodModel, 
                           Vec<LikelihoodType*> &likelihoodModels)
{
    int numModels  = iLikelihoodModel.length();
    
    int countGaussianModels = 0;   // The number of Gaussian noise models correctly parsed
    double avgErrMean = 0.0, avgErrVar = 0.0;  // The average mean and variance of these models

    for(int iModel=0; iModel<numModels; iModel++) 
    {
        // Get parameter vector for current model
        vec params = modelParams[iModel];

        // If a valid model name/params exist for this obs, 
        // create likelihood model
        if (modelName[iModel] == "GAUSSIAN") 
        {
            // TODO: At the moment, the GaussianLikelihood takes only the variance parameter
            // Change to GaussianLikelihood(mean, variance)
            // Also, add support for LikelihoodType(vec params) to make it more generic
            likelihoodModels[iModel] = new GaussianLikelihood( params(1) );

            // TODO: uncomment line below once the mean is taken into account in GaussianLikelihood
            // avgMean += distParams(0);
            avgErrVar += modelParams[iModel](1);
            countGaussianModels++;
        }
        /* Use template below to add other noise models 
                ---
                else if (modelName[iObs] == "MY_MODEL")
                {
                    myLikelihoods.push_back( MyLikelihoodType( params(0), ..., params(n) ) );
                    likelihoodModels[iModel] = new MyLikelihoodType( params(0), ..., params(n) );
                }
                ---
         */
        else
        {
            cout << "Unrecognized observation noise model: " << modelName[iModel] << 
            "for observation " << iModel << endl;

            // Flag observation as having an invalid noise model 
            modelName[iModel] = INVALID_MODEL_NAME;
        }
    }

    // If we found some valid Gaussian noise models, infer parameters for
    // default likelihood from theirs (take average mean and average variance).
    // Should hopefuly be slightly better than our first guess.
    if (countGaussianModels > 1) {
        avgErrMean = avgErrMean/countGaussianModels;
        avgErrVar  = avgErrVar/countGaussianModels;

        cout << "Observations without a valid noise model are given a default" << endl;
        cout << "Gaussian noise model with (mean, variance) = (" << avgErrMean; 
        cout << "," << avgErrVar << ")" << endl;
    }


    // Attribute all observations with an invalid noise model the default 
    // (Gaussian) likelihood model
    for(int iModel=0; iModel<numModels; iModel++) 
    {
        if(modelName[iModel] == INVALID_MODEL_NAME) 
        {
            // likelihoodModels[iLikelihoodModel(iModel)] = defaultLikelihood;
            likelihoodModels[iLikelihoodModel(iModel)] = new GaussianLikelihood(avgErrVar);
        }
    }
}

