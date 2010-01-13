//
// Generalized Rosenbrock Function for optimisation
//
// History:
//
// 19 Dec 2008 - BRI - First implementation.
//

#include "RosenbrockND.h"

using namespace std;

RosenbrockND::RosenbrockND(int dimensions) : dim(dimensions)
{
	x = randu(dimensions);
}

RosenbrockND::~RosenbrockND()
{
}

vec RosenbrockND::getParametersVector() const
{
	vec p(x);
	return p;
}

void RosenbrockND::setParametersVector(const vec p)
{
	assert(p.size() == dim);
	for(int i = 0 ; i < dim ; i++)
	{
		x(i) = p(i);
	}
}

double RosenbrockND::objective() const
{
	double sum = 0.0;
	for(int i = 0 ; i < (dim - 1) ; i++)
	{
		sum = sum + (100 * pow(x(i + 1) - pow(x(i),2),2) + pow((1.0 - x(i)),2));
	}
	return sum;
}

vec RosenbrockND::gradient() const
{
	cerr << "Analytic gradients not implemented" << endl;

	vec p;
	p.set_size(dim, false);

	return p;
}


