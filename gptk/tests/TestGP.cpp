#include "TestGP.h"

TestGP::TestGP()
{
	addTest(&testGPNoisySineLearnParams, "Noisy sine with parameter estimation");
}

TestGP::~TestGP()
{
}

bool TestGP::testGPNoisySineLearnParams() 
{
	
	