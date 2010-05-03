#ifndef NEURALNETCF_H_
#define NEURALNETCF_H_

#include "CovarianceFunction.h"

/**
 * Neural network covariance function
 * 
 * This covariance function is of the form
 *  
 *   s * asin[ r * (1+x'*y) / sqrt( (1+r+r*x'*x)*(1+r+r*y'*y) ) ]
 * 
 * where r = 1/(lengthscale^2) and s = variance is a scaling factor.
 * 
 * (C) 2009 Remi Barillec <r.barillec@aston.ac.uk>
 */
class NeuralNetCF : public CovarianceFunction
{
public:
	NeuralNetCF(double lengthscale, double variance, double offset);
	virtual ~NeuralNetCF();

	inline double computeElement(const vec& A, const vec& B) const;

	virtual void covarianceGradient(mat& grad, const int parameterNumber, const mat& X) const;

private:
    double &lengthScale;                // Length scale parameter
    double &offset;                     // Length scale of the offset
    double &variance;                   // Variance parameter
};

#endif /*NEURALNETCF_H_*/
