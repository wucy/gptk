
#include <stdlib.h>
#include <stdio.h>

#include "astonGeostats.h"
#include "R.h"
#include "Rmath.h"
#include "Rinternals.h"

#define NUM_VARIOGRAM_PARAMETERS 5

using namespace std;

extern "C" {
    SEXP predict(SEXP xData, SEXP yData, SEXP eData, SEXP xPred, SEXP psgpPar, 
                 SEXP errorIdx, SEXP sensorIdx, SEXP metaData)
    {

        SEXP meanResult;
        SEXP varResult;
        SEXP ans;
        
        // These are not used by PSGP anymore. Instead, we use the parameters 
        // stored in psgpParameters 
        // double range, sill, nugget, bias;
        int xDataLen, yDataLen, eDataLen, xPredLen;
        double *xDataPtr, *yDataPtr, *xPredPtr, *eDataPtr;
        
        // PSGP covariance function parameters
        double  *psgpParameters = REAL(psgpPar);
        
        int     *errorPtr, *sensorPtr;
        char    **metaDataTable;

        int metadataSize = length(metaData);
        // The metadata table, if provided, is terminated by an empty line
        // Need to remove it form the size
        if (metadataSize > 0) metadataSize -= 1;

        // there must be a more obvious way to get the dimensions of a matrix?
        xDataLen = length(xData) / 2;
        yDataLen = length(yData);
        xPredLen = length(xPred) / 2;
        eDataLen = length(eData);

        PROTECT(meanResult = allocVector(REALSXP, xPredLen));
        PROTECT(varResult = allocVector(REALSXP, xPredLen));
        PROTECT(ans = allocVector(VECSXP, 2));

        xPredPtr = REAL(xPred);
        xDataPtr = REAL(xData);
        yDataPtr = REAL(yData);
        eDataPtr = REAL(eData);
        errorPtr = INTEGER(errorIdx);
        sensorPtr= INTEGER(sensorIdx);

        // If metadata is specified, it should have one likelihood model per
        // observation
        assert(metadataSize == 0 || metadataSize == xDataLen);
        
        // we need to make a table of pointers to the sensor model strings
        metaDataTable = (char**)calloc(metadataSize, sizeof(char *));

        // casting from a const pointer to a non-const pointer is a bit
        // of nightmare, but I'm not really sure what the other options are?
        for (int i = 0; i < metadataSize; i++ )
        {
            metaDataTable[i] = const_cast<char*>(CHAR(STRING_ELT(VECTOR_ELT(metaData, i),0)));
        }

        /*
        makePredictions(xDataLen, xPredLen, eDataLen, xDataPtr, yDataPtr, xPredPtr, eDataPtr, 
                metadataSize, errorPtr, sensorPtr, metaDataTable,
                REAL(meanResult), REAL(varResult),
                &range, &sill, &nugget, &bias, model);
        */
        printf("Making predictions\n");
        makePredictions(xDataLen, xPredLen, eDataLen, xDataPtr, yDataPtr, xPredPtr, eDataPtr, 
                        metadataSize, errorPtr, sensorPtr, metaDataTable,
                        REAL(meanResult), REAL(varResult), psgpParameters);
        
        SET_VECTOR_ELT(ans, 0, meanResult);
        SET_VECTOR_ELT(ans, 1, varResult);

        UNPROTECT(3);
        return ans;
    }
}

extern "C" {
    SEXP estParam(SEXP xData, SEXP yData, SEXP eData, SEXP vario, SEXP errorIdx, SEXP sensorIdx, SEXP metaData)
    {
        // SEXP meanResult;
        // SEXP varResult;
        SEXP params;
        int xDataLen, eDataLen;
        double *xDataPtr, *yDataPtr, *eDataPtr, *varioPtr;
        int *errorPtr, *sensorPtr;
        // int model = 0;
        char **metaDataTable;
        

        int metadataSize = length(metaData);
        // The metadata table, if provided, is terminated by an empty line
        // Need to remove it form the size
        if (metadataSize > 0) metadataSize -= 1;
        
        // there must be a more obvious way to get the dimensions of a matrix?
        xDataLen = length(xData) / 2;
        eDataLen = length(eData);
        xDataPtr = REAL(xData);
        yDataPtr = REAL(yData);
        eDataPtr = REAL(eData);
        varioPtr = REAL(vario);
        errorPtr = INTEGER(errorIdx);
        sensorPtr= INTEGER(sensorIdx);

        // Check that all data arrays are coherent 
        assert(xDataLen == length(yData));
        assert(xDataLen == eDataLen);
        
        // If metadata is specified, it should have one likelihood model per
        // observation
        assert(metadataSize == 0 || metadataSize == xDataLen);
                
        // we need to make a table of pointers to the sensor model strings
        metaDataTable = (char**)calloc(metadataSize, sizeof(char *));

        // casting from a const pointer to a non-const pointer is a bit
        // of nightmare, but I'm not really sure what the other options are?
        for (int i = 0; i < metadataSize; i++ )
        {
            metaDataTable[i] = const_cast<char*>(CHAR(STRING_ELT(VECTOR_ELT(metaData, i),0)));
        }

        // Variogram and PSGP parameters are stored in the same array, which is
        // eventrally returned on completion of this function.
        // 
        // The first NUM_VARIOGRAM_PARAMETERS parameters are for the variogram
        // and the next NUM_PSGP_PARAMETERS (defined in astonGeostats.h) parameters
        // are for the PSGP covariance function.
        PROTECT(params = allocVector(REALSXP, NUM_PSGP_PARAMETERS));
        double* psgpParameters = REAL(params);
        UNPROTECT(1);
        
        // Copy current variogram parameters to parameter array
        memcpy(psgpParameters, varioPtr, NUM_VARIOGRAM_PARAMETERS * sizeof(double));
        
        /*
        printf("Variogram model:\n");
        printf("  type   = %d\n", (int) varioPtr[0]);
        printf("  range  = %f\n", varioPtr[1]);
        printf("  sill   = %f\n", varioPtr[2]);
        printf("  nugget = %f\n", varioPtr[3]);
        printf("  bias   = %f\n", varioPtr[4]);
        */
        
        
        // Estimate parameters.
        // This also updates the parameter values in psgpParameters 
        // and in variogramParameters (if the initial parameters were not
        // valid, i.e. negative...)
        learnParameters(xDataLen, eDataLen, xDataPtr, yDataPtr, eDataPtr,
                        metadataSize, errorPtr, sensorPtr, 
                        metaDataTable, psgpParameters);
        
        /*
        // Display resulting parameter vector
        for(int i=0; i<NUM_VARIOGRAM_PARAMETERS+NUM_PSGP_PARAMETERS; i++) {
            printf("Parameter %2d: %f\n", i, paramsPtr[i]);
        }
        */
        
        return params;
    }
}

