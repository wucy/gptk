#ifndef TESTGP_H_
#define TESTGP_H_

#include "Test.h"
#include "gaussianProcesses/GaussianProcess.h"

class TestGP : public Test
{
public:
	TestGP();
	virtual ~TestGP();
	
	static bool testGPNoisySineLearnParams();
};

#endif /*TESTGP_H_*/
