/*
 * TestOptimisationSCG.cpp
 *
 *  Created on: 27-Jun-2009
 *      Author: barillrl
 */

#include "TestOptimisationSCG.h"

QuadraticError::QuadraticError(double x0, double y0)
{
	x = x0;
	y = y0;
}

QuadraticError::~QuadraticError()
{

}

double QuadraticError::objective() const
{
	//double x = params(0);
	//double y = params(1);
	return (x-3.2)*(x-3.2)+(y-3.2)*(y-3.2);
}

vec QuadraticError::gradient() const
{
	//double x = params(0);
	//double y = params(1);

	vec g = zeros(2);
	g(0) = (2.0*x-6.4);
	g(1) = (2.0*y-6.4);

	return g;
}

vec QuadraticError::getParametersVector() const
{
	vec params(2);
	params(0) = x;
	params(1) = y;
	return params;
}


void QuadraticError::setParametersVector(const vec p)
{
	x = p(0);
	y = p(1);
}

TestOptimisationSCG::TestOptimisationSCG() {
	header = "Test for class SequentialGP";

	// addTest(&testQuadraticError, "Quadratic error function");
	addTest(&testRosenbrock, "Rosenbrock function");

}

TestOptimisationSCG::~TestOptimisationSCG() {
}

bool TestOptimisationSCG::testQuadraticError()
{
	QuadraticError err(10.1, 23.5);
	SCGModelTrainer scg(err);
	scg.Train(100);

	vec params = err.getParametersVector();
	return (params(0) - 3.2 < 1e-4) && (params(1) - 3.2 < 1e-4);
}


bool TestOptimisationSCG::testRosenbrock()
{
    int iterations = 100;
    cout.setf(ios::fixed);
    cout.precision(4);
    Rosenbrock rosen(-1.0, 1.0);
    vec startingValue = rosen.getParametersVector();
    SCGModelTrainer scgTrainer(rosen);
    rosen.setParametersVector(startingValue);
    cout << "Scaled Conjugate Gradient" << endl;
    scgTrainer.setAnalyticGradients(true);
    scgTrainer.Train(iterations);
    scgTrainer.Summary();
    vec finalParams = rosen.getParametersVector();
    cout << "Final Parameters: " << finalParams << endl << endl;
    return (abs(finalParams(0)-1.0) < 1e-4 && abs(finalParams(1)-1.0 < 1e-4));
}

int main() {
	TestOptimisationSCG test;
	test.run();
}
