#ifndef TESTITPPEXT_H_
#define TESTITPPEXT_H_

#include "Test.h"
#include "itppext/itppext.h" 

using namespace itppext;

class TestITPPext : public Test
{
public:
	TestITPPext();
	virtual ~TestITPPext();
	
	static bool testTriangularPacking();
	static bool testConcatenation();
	static bool testNormalisation(); 
};

#endif /*TESTITPPEXT_H_*/
