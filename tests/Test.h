#ifndef TEST_H_
#define TEST_H_

#include <iostream>
#include <vector>
#include <string>
#include <itpp/itbase.h>
#include <cassert>

using namespace std;
using namespace itpp;

class Test
{
public:
	/**
	 * Constructor
	 * Pass in the description for the current test set
	 */
	Test();
	
	virtual ~Test();
			

	/**
	 * Adds a test to the list of tests
	 */ 
	void addTest(bool function(), string header);
	
	/**
	 * Test wrapper
	 * This simply runs the tests and outputs a failed/passed message
	 */
	void runTest( bool testFunction() );
	
	/**
	 * Runs all tests
	 */ 
	void run();

	
protected:
	string header;                         // Overall description of this test class
	vector<string> descriptions;           // Vector of test descriptions
	vector<bool(*)()> tests;       // Vector of (test) function pointers   
	
	
	/**
	 * Displays the overall header for this test set
	 */ 
	void printHeader();

	/**
	 * Displays the description for a single test
	 */
	void printTestDescription(string message);

};

#endif /*TEST_H_*/
