#include "TestGradientCovFunc.h"

TestGradientCovFunc::TestGradientCovFunc() 
{
  header = "Test set for gradient of covariance functions";
  addTest(&testGradientGaussianCF, "Gradient of Gaussian Covariance Function");
  addTest(&testGradientWhiteNoiseCF, "Gradient of Gaussian White Noise Covariance Function");
  addTest(&testGradientSumCF, "Gradient of Sum of covariance functions");
}

TestGradientCovFunc::~TestGradientCovFunc() {}
  
/**
 * Test gradient of Gaussian covariance function
 */
bool TestGradientCovFunc::testGradientGaussianCF()
{
  double lengthScale = 2.1;
  double variance    = 3.3;
  GaussianCF *cf = new GaussianCF(lengthScale, variance);
  return matGradCheck(cf);
}

/**
 * Test gradient of WhiteNoise covariance function
 */
bool TestGradientCovFunc::testGradientWhiteNoiseCF()
{
  double variance = 4.1;
  WhiteNoiseCF *cf = new WhiteNoiseCF(variance);
  return matGradCheck(cf);
}

/**
 * Test gradient of sum of 2 covariance functions
 */ 
bool TestGradientCovFunc::testGradientSumCF() {
	double lengthScale = 2.1;
	double variance    = 3.3;
	double nugget = 4.1;
	
	GaussianCF 		cf1(lengthScale, variance);
	WhiteNoiseCF 	cf2(variance);

	SumCovarianceFunction *cf = new SumCovarianceFunction(cf1);
	cf->addCovarianceFunction(cf2);
	
	return matGradCheck(cf);
}

/**
 * A multivariate version of Netlab's gradchek function.
 * 
 * This computes the derivative of the covariance function with 
 * respect to each parameter and compares with a finite difference
 * estimate. Since the derivative is a matrix, we look at the error
 * (per matrix element) between the analytic derivative and the finite
 * difference approximation. The mean, variance and max of this error
 * are displayed on the standard output.
 * 
 * The function returns true if the mean, error and variance are all
 * below some minimum value, for all parameters.    
 */
bool TestGradientCovFunc::matGradCheck(CovarianceFunction *cf)
{
  vec params = cf->getParameters();
  vec X = linspace(0,10,100);
  mat gradK(X.size(),X.size()), gradKfd(X.size(),X.size());
  double max_errmean = 0.0, max_errvar = 0.0, max_errmax = 0.0;
  double tolerance = 1e-3;
  
  // For each paramter, compute the partial derivative matrix
  // using gradient and finite difference. Compute the error 
  // between the two and print out the mean and variance of 
  // this error.
  printf("\nMATRIX GRADCHECK: Error between gradient and finite difference\n");
  printf("\nParam   Mean       Var        Max\n");
  for (int i=0; i<params.size(); i++) {
    cf->getParameterPartialDerivative(gradK, i, X);
    gradKfd = gradFiniteDifferences(cf,X,i);
        
    mat err = abs(gradK - gradKfd);
    double errmean = mean(err);
    double errvar  = sum(sum(pow(err-errmean,2)))/(err.cols()*err.rows()-1);
    double errmax  = max(max(err));
   
    // Update current max value for mean, var and max
    if (errmean > max_errmean) max_errmean = errmean;
    if (errvar  > max_errvar)  max_errvar  = errvar;
    if (errmax  > max_errmax)  max_errmax  = errmax;
    
    printf("%3d: %10.5f %10.5f %10.5f\n", i, errmean, errvar, errmax);
  }
  
  return ((max_errmean < tolerance) && (max_errvar < tolerance) && (max_errmax < tolerance));
  
} 

/**
 * Computes the finite difference gradient of the specified covariance function
 * with repect to parameter i. The covariance function is evaluated at point X, 
 * i.e. this returns dK(X,X)/d\theta_i.
 */
mat TestGradientCovFunc::gradFiniteDifferences(CovarianceFunction *cf, vec X, int i)
{
  double epsilon = 1e-6;
  mat Kplus(X.size(),X.size()), Kminus(X.size(),X.size());
  vec params = cf->getParameters();
  
  params(i) += epsilon;
  cf->setParameters(params);
  // cf->computeCovariance(Kplus,X,X);
  cf->computeSymmetric(Kplus,X);
  
  params(i) -= 2.0*epsilon;
  cf->setParameters(params);
  // cf->computeCovariance(Kminus,X,X);
  cf->computeSymmetric(Kminus,X);
  
  return (Kplus - Kminus)/(2.0*epsilon);
}


/**
 * Run the tests
 */
int main() {
  TestGradientCovFunc test;
  test.run();
}
