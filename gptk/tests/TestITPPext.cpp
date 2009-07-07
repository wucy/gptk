#include "TestITPPext.h"

TestITPPext::TestITPPext()
{
    addTest(&testTriangularPacking, "Packing/unpacking of triangular matrices");
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



int main() {
    TestITPPext test;
    test.run();
}
