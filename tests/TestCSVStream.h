#ifndef TESTCSVSTREAM_H_
#define TESTCSVSTREAM_H_

#include "Test.h"
#include <io/csvstream.h>

class TestCSVStream : public Test
{
public:
	TestCSVStream();
	virtual ~TestCSVStream();
	
	static bool testReadWrite();
	static bool testReadWrite1col();
	static bool testReadWrite1row();
	static bool testReadWriteEmpty();

private:
	static bool readWrite(const mat M);
	
};

#endif /*TESTCSVSTREAM_H_*/
