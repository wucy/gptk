#include "TestITPPext.h"
#include <stdexcept>

TestITPPext::TestITPPext()
{
    addTest(&testTriangularPacking, "Packing/unpacking of triangular matrices");
    addTest(&testConcatenation, "Concatenation functions");
    addTest(&testNormalisation, "Normalisation functions");
    addTest(&testTokenise, "Tokeniser");
}

TestITPPext::~TestITPPext()
{
}

bool TestITPPext::testTriangularPacking()
{
    int N = 100;
    mat M = 100.0*randu(N,N);
    
    // Lower triangularise matrix M
    for (int i=0; i<N; i++)
        for (int j=i; j<N; j++)
            M(i,j) = 0.0;
    
    // Lower triangular packing/unpacking
    vec l = ltr_vec(M); 
    mat L = ltr_mat(l);
   
    // Upper triangular packing/unpacking
    vec u = utr_vec(M.transpose());
    mat U = utr_mat(u);
    
    return (L==M && U==M.transpose());
    
}

bool TestITPPext::testConcatenation() 
{
    int ncols = 7;
    int nrows = 13;
    
    mat M = 100.0*randu(nrows, ncols);
    mat P = M;
    vec y = M.get_col(ncols-1);
    P.set_size(nrows, ncols-1, true);
     
    /*
    cout << P << endl;
    cout << y << endl;
    cout << concat_cols(P,y) << endl;
    cout << M << endl;
    */
    return (M == concat_cols(P,y));
}


bool TestITPPext::testNormalisation() 
{
    int ncols = 3;
    int nrows = 7;
    
    mat M = 100.0*randu(nrows, ncols);
    vec y = 100.0*randu(nrows);
    
    mat M2 = M;
    vec y2 = y;
    
    vec zmean, zcov;
    
    normalise(M2, zmean, zcov);
    denormalise(M2, zmean, zcov);
    
    // Test passed if absolute element error < 1e-6
    return (sum(sum(abs(M - M2)))/(ncols*nrows) < 1e-6);
}


bool TestITPPext::testTokenise() 
{
    string s[5] = {"1.1","2.2","3.3","4.456","7.89"};
    string s1 = "1.1 2.2 3.3 4.456 7.89";
    string s2 = ",,,1.1,,2.2,3.3,4.456,7.89";
    string s3 = "1.1";
    string s4 = "2.2";
    string s5 = "3.3,4.456,7.89";
    
    vector<string> tokens;
    
    // Tokenise with default delimiter (space)
    tokenise(s1, tokens);
    for(int i=0; i<tokens.size(); i++) {
        if (tokens[i] != s[i]) return false;
    }
    
    // Tokenise with custom delimiter (comma)
    tokens.clear();
    tokenise(s2, tokens,",");
    for(int i=0; i<tokens.size(); i++) {
        if (tokens[i] != s[i]) return false;
    }
    
    // Tokenise with stopping after first token
    tokens.clear();
    tokenise(s2, tokens, ",", 2);
    for(int i=0; i<tokens.size(); i++) cout << tokens[i] << endl;
    if (tokens[0] != s3 || tokens[1] != s4 || tokens[2] != s5) return false;
    
    // Try to crash ITPP
    try {
        vec a("asdlkjaslkdj,dsas,ad,asd,asd,asd");
    }
    catch (std::runtime_error e) {
        cerr << "Invalid parameter string. Should be a sequence of numeric values" << endl;
        cerr << "separated by commas, e.g. \"1.23,4,5.6,78.9\"" << endl;
    }
    return true;
}


int main() {
    TestITPPext test;
    test.run();
}
