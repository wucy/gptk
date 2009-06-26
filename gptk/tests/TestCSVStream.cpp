#include "TestCSVStream.h"

TestCSVStream::TestCSVStream() : Test()
{
	header = "Test set for csvstream class";
	addTest(&testReadWrite, "Writing and reading a CSV file from a generic matrix.");
	addTest(&testReadWrite1col, "Writing and reading a CSV file from a 1 column matrix.");
	addTest(&testReadWrite1row, "Writing and reading a CSV file from a 1 row matrix.");
	addTest(&testReadWriteEmpty, "Writing and reading a CSV file from an empty matrix.");
}

TestCSVStream::~TestCSVStream()
{
}

bool TestCSVStream::readWrite(const mat M) {
	csvstream csv;
	mat P, Q;
	double precision = 5.0;   // Use 5 decimal places

	// CSV format has limited precision
	Q = floor(pow(10.0,precision)*M)/pow(10.0,precision);
	
	csv.write(Q,"csvstream_test.csv", (int) precision);
	csv.read(P, "csvstream_test.csv");
	
	return (P==Q);
}

bool TestCSVStream::testReadWrite() 
{
	mat M = randn(20,15);
	return readWrite(M);
}

bool TestCSVStream::testReadWrite1col() 
{
	mat M = randn(20,1);
	return readWrite(M);
}

bool TestCSVStream::testReadWrite1row() 
{
	mat M = randn(1,20);
	return readWrite(M);
}

bool TestCSVStream::testReadWriteEmpty() 
{
	mat M;
	M.set_size(0,0);
	return readWrite(M);
}


int main() {
	TestCSVStream test;
	test.run();
}
