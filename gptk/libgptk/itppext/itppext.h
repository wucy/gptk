#ifndef ITPPEXT_H_
#define ITPPEXT_H_

/**
 * Extension to IT++ library
 */

#include <cassert>
#include <math.h>
#include <itpp/itbase.h>
#include <itpp/itstat.h>

using namespace itpp;

namespace itppext 
{

vec ltr_vec(mat M);     // Vector of lower triangular elements
vec utr_vec(mat M);     // Vector of upper triangular elements
mat ltr_mat(vec v);     // Lower triangular matrix
mat utr_mat(vec v);     // Upper triangular matrix

double cond(mat M, int p=2); // Condition number for matrix p-norm (1 or 2)

ivec randperm(int n);  // Random permutation of numbers between 0 and N-1

vec min(vec u, vec v); // Minimum elements from 2 vectors of equal length 

} // END OF NAMESPACE ITPPEXT


#endif /*ITPPEXT_H_*/
