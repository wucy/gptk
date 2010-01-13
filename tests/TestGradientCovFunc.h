#ifndef TESTGRADIENTCOVFUNC_H_
#define TESTGRADIENTCOVFUNC_H_

#include "Test.h"
#include "covarianceFunctions/GaussianCF.h"
#include "covarianceFunctions/WhiteNoiseCF.h"
#include "covarianceFunctions/SumCovarianceFunction.h"

#include <itpp/itbase.h>
#include <itpp/stat/misc_stat.h>

using namespace std;
using namespace itpp;

class TestGradientCovFunc : public Test
{
public:
  TestGradientCovFunc();
  virtual ~TestGradientCovFunc();
  
  /**
   * Test gradient of Gaussian covariance function
   */
  static bool testGradientWhiteNoiseCF();
  
  /**
   * Test gradient of WhiteNoise covariance function
   */
  static bool testGradientGaussianCF();
    
  /**
   * Test gradient of sum of 2 covariance functions
   */ 
  static bool testGradientSumCF();

  /**
   * Computes the error between analytic gradient matrix and finite differences
   * estimate. For each parameter of the covariance function, the mean and variance 
   * of the error (between the elements of the two gradient matrices) are displayed.
   * Returns true if mean, var and max below fixed tolerance (1e-3)
   */
  static bool matGradCheck(CovarianceFunction *cf);

private:
  
  /**
   * Computes the gradient matrix using finite differences
   * Pass in a covariance function and its parameters
   */
  static mat gradFiniteDifferences(CovarianceFunction *cf, vec X, int paramIndex);
  
};

#endif /*TESTGRADIENTCOVFUNC_H_*/
