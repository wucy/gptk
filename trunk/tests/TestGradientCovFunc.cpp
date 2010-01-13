#include "TestGradientCovFunc.h"

TestGradientCovFunc::TestGradientCovFunc() 
{
  header = "Test set for gradient of covariance functions";
  addTest(&testGradientGaussianCF, "Gradient of Gaussian covariance function");
  addTest(&testGradientWhiteNoiseCF, "Gradient of Gaussian White Noise covariance function");
//  addTest(&testGradientConstantCF, "Gradient of Constant covariance functions");
//  addTest(&testGradientMatern3CF, "Gradient of Matern 3/2 covariance function");
//  addTest(&testGradientMatern5CF, "Gradient of Matern 5/2 covariance function");
//  addTest(&testGradientSumCF, "Gradient of Sum of covariance functions");
//  addTest(&testGradientNeuralNetCF, "Gradient of neural network covariance function");
  // addTest(&testGradientExponential2CF, "Gradient of anisotropic exponential covariance function");
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
 * Test gradient of Constant covariance function
 */
bool TestGradientCovFunc::testGradientConstantCF()
{
    double amplitude = 0.0123;
//  ConstantCF *cf = new ConstantCF(amplitude);
//  return matGradCheck(cf);
    return false;
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
	
//	GaussianCF 		cf1(lengthScale, variance);
//	WhiteNoiseCF 	cf2(variance);
//
//	SumCovarianceFunction *cf = new SumCovarianceFunction(cf1);
//	cf->addCovarianceFunction(cf2);
//
//	return matGradCheck(cf);
	return false;
}

/**
 * Test gradient of Matern 3/2 covariance function
 */
bool TestGradientCovFunc::testGradientMatern3CF()
{
    double lengthScale = 2.1;
    double variance = 3.7;
//    Matern3CF *cf = new Matern3CF(lengthScale, variance);
//
//    return matGradCheck(cf);
    return false;
}


/**
 * Test gradient of Matern 5/2 covariance function
 */
bool TestGradientCovFunc::testGradientMatern5CF()
{
    double lengthScale = 2.1;
    double variance = 3.7;
//    Matern5CF *cf = new Matern5CF(lengthScale, variance);
//
//    return matGradCheck(cf);
    return false;
}


/**
 * Test gradient of neural network covariance function
 */
bool TestGradientCovFunc::testGradientNeuralNetCF()
{
    double lengthScale = 2.1;
    double variance = 3.7;
    double offset = 1.2;
//    NeuralNetCF *cf = new NeuralNetCF(lengthScale, variance, offset);
//
//    return matGradCheck(cf);
    return false;
}
    

bool TestGradientCovFunc::testGradientExponential2CF()
{
    double ls_x = sqrt(2.1);
    double ls_y = sqrt(2.6);
    double ls_xy = -0.23;
    
    double variance = 3.7;
//    Exponential2CF *cf = new Exponential2CF(ls_x, ls_y, ls_xy, variance);
    
//    return matGradCheck(cf);
    return false;
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
  vec params = cf->getTransformedParameters();
  mat X = 10.0*randn(10,2);
  mat gradK(X.rows(),X.rows()), gradKfd(X.rows(),X.rows());
  double max_errmean = 0.0, max_errvar = 0.0, max_errmax = 0.0;
  double tolerance = 1e-3;
  
  // For each paramter, compute the partial derivative matrix
  // using gradient and finite difference. Compute the error 
  // between the two and print out the mean and variance of 
  // this error.
  printf("\n  MATRIX GRADCHECK: Error between gradient and finite difference\n");
  printf("\n  Param   Mean       Var        Max       Parameter name\n");
  for (int i=0; i<params.size(); i++) 
  {
    cf->covariance_gradient(gradK, i, X);
    gradKfd = gradFiniteDifferences(cf,X,i);
        
    // cout << "grad analytic = " << endl << gradK << endl;
    // cout << "grad fin.diff. = " << endl << gradKfd << endl;
    
    mat err = abs(gradK - gradKfd);
    double errmean = mean(err);
    double errvar  = sum(sum(pow(err-errmean,2)))/(err.cols()*err.rows()-1);
    double errmax  = max(max(err));
   
    // Update current max value for mean, var and max
    if (errmean > max_errmean) max_errmean = errmean;
    if (errvar  > max_errvar)  max_errvar  = errvar;
    if (errmax  > max_errmax)  max_errmax  = errmax;
    
    printf("  %3d: %10.5f %10.5f %10.5f   %s\n", i, errmean, errvar, errmax,
                                               (cf->getParameterName(i)).c_str());
  }
                    
  return ((max_errmean < tolerance) && (max_errvar < tolerance) && (max_errmax < tolerance));
  
} 


/**
 * Computes the finite difference gradient of the specified covariance function
 * with repect to parameter i. The covariance function is evaluated at point X, 
 * i.e. this returns dK(X,X)/d\theta_i.
 */
mat TestGradientCovFunc::gradFiniteDifferences(CovarianceFunction *cf, mat X, int i)
{
  double epsilon = 1e-6;
  mat Kplus(X.rows(),X.rows()), Kminus(X.rows(),X.rows());
  vec params = cf->getTransformedParameters();
  
  params(i) += epsilon;
  cf->setTransformedParameters(params);
  cf->covariance(Kplus,X);
  
  params(i) -= 2.0*epsilon;
  cf->setTransformedParameters(params);
  cf->covariance(Kminus,X);
  
  return (Kplus - Kminus)/(2.0*epsilon);
}


/**
 * Run the tests
 */
int main() {
  TestGradientCovFunc test;
  test.run();
}
