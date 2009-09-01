#include "TestITPPext.h"

TestITPPext::TestITPPext()
{
    addTest(&testTriangularPacking, "Packing/unpacking of triangular matrices");
    addTest(&testConcatenation, "Concatenation functions");
    addTest(&testNormalisation, "Normalisation functions");
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


int main() {
    TestITPPext test;
    test.run();
}
