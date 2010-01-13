/*
 * TestOptimisationSCG.h
 *
 *  Created on: 27-Jun-2009
 *      Author: barillrl
 */

#ifndef TESTOPTIMISATIONSCG_H_
#define TESTOPTIMISATIONSCG_H_

#include "Test.h"
#include "optimisation/Optimisable.h"
#include "optimisation/SCGModelTrainer.h"
#include <cassert>
#include "../examples/Rosenbrock.h"

class QuadraticError : public Optimisable
{
public:
	QuadraticError(double x0, double y0);
	~QuadraticError();

	double 	objective() const;
	vec 	gradient() const;
	vec 	getParametersVector() const;
	void 	setParametersVector(const vec p);

private:
	double x;
	double y;
};

class TestOptimisationSCG: public Test
{
public:
	TestOptimisationSCG();
	virtual ~TestOptimisationSCG();

	static bool testQuadraticError();
	static bool testRosenbrock();
};

#endif /* TESTOPTIMISATIONSCG_H_ */
