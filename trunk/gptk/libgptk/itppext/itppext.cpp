#include "itppext.h"

namespace itppext {

/**
 * Returns the lower triangular elements (including diagonal) 
 * of a matrix, in row-major order.
 * 
 * @param M An NxN square matrix
 * @return The vector of N(N+1)/2 lower triangular elements of M
 * @see ltr_mat, utr_vec, utr_mat
 */   
vec ltr_vec(mat M) 
{
    int N = M.cols();
    assert(N==M.rows());
        
    int k = 0;
    vec v(N*(N+1)/2); 
    
    for (int i=0; i<N; i++) 
        for (int j=0; j<=i; j++) 
            v(k++) = M(i,j); 
    
    return v;
}

/**
 * Returns the upper triangular elements (including diagonal) 
 * of a matrix, in row-major order.
 * 
 * @param M An NxN square matrix
 * @return The vector of N(N+1)/2 upper triangular elements of M
 * @see utr_mat, ltr_vec, ltr_mat
 */   
vec utr_vec(mat M) 
{
    int N = M.cols();
    assert(N==M.rows());
    
    int k = 0;
    vec v(N*(N+1)/2); 
    
    for (int i=0; i<N; i++) 
        for (int j=i; j<N; j++) 
            v(k++) = M(i,j); 
    
    return v;
}

/**
 * Returns a lower triangular matrix where the elements are taken
 * from the argument vector in row-major order.
 * 
 * @param v A vector of N(N+1)/2 lower triangular elements
 * @return A NxN lower triangular matrix
 * @see ltr_vec, utr_vec, utr_mat
 */
mat ltr_mat(vec v)
{
    // Retrieve dimension of matrix
    int N =  floor(sqrt(2*v.length()));
    
    assert(N*(N+1)/2 = v.length());
    
    mat M = zeros(N,N);
    int k = 0;
    
    for (int i=0; i<N; i++) 
            for (int j=0; j<=i; j++) 
                M(i,j) = v(k++);  
    
    return M;
}

/**
 * Returns an upper triangular matrix where the elements are taken
 * from the argument vector in row-major order.
 * 
 * @param v A vector of N(N+1)/2 upper triangular elements
 * @return A NxN upper triangular matrix
 * @see ltr_vec, utr_vec, utr_mat
 */
mat utr_mat(vec v)
{
    // Retrieve dimension of matrix
    int N =  floor(sqrt(2*v.length()));
    
    assert(N*(N+1)/2 = v.length());
    
    mat M = zeros(N,N);
    int k = 0;
    
    for (int i=0; i<N; i++) 
            for (int j=i; j<N; j++) 
                M(i,j) = v(k++);  
    
    return M;
}

double cond(mat M, int p) {
    assert(M.rows() == M.cols());
    return itpp::norm(M,p)*itpp::norm(inv(M),p);
}

/**
 * Returns a random permutation of numbers between 0 and N-1
 */
ivec randperm(int n) {
	vec rndNums = randu(n);
	return sort_index(rndNums);
}

/**
 * Returns the vector of minimum elements from 2 vectors, i.e.
 * z(i) = min(u(i), v(i)).
 */
vec min(vec u, vec v) {
    assert(u.length() == v.length());
    
    vec z(u.length());
    
    for (int i=0; i<u.length(); i++) z(i) = std::min(u(i), v(i));
    
    return z;
}

} // END OF NAMESPACE ITPPEXT


