#ifndef CSVSTREAM_H_
#define CSVSTREAM_H_

// Uncomment following line to enable debugging messages
// #define DEBUG
#include "../debug.h"
#include <iostream>
#include <string>
#include <vector>
#include <itpp/itbase.h>


using namespace std;
using namespace itpp;

class csvstream
{
public:
	csvstream();
	virtual ~csvstream();
	
	/**
	 * Reads a matrix from a CSV file
	 * 
	 * @param matrix The matrix to which the data is read (does not need to be pre-allocated)
	 * @param filename The name of the CSV file to read 
	 **/
	int read(mat &matrix, const string filename);

	/**
	 * Writes a matrix to a CSV file
	 * 
	 * @param matrix The matrix containing the data
	 * @param filename The name of the CSV file to write to 
	 * @params decimals The number of retained decimals 
	 **/
	int write(const mat matrix, const string filename, int decimals = 5);

};

#endif /*CSVSTREAM_H_*/
