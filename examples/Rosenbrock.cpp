//
// Rosenbrock Function for optimisation
//
// History:
//
// 19 Dec 2008 - BRI - First implementation.
//

#include "Rosenbrock.h"

using namespace std;

Rosenbrock::Rosenbrock(double p1, double p2)
{
	x1 = p1;
	x2 = p2;
}

Rosenbrock::~Rosenbrock()
{
}

vec Rosenbrock::getParametersVector() const
{
	vec p;
	p.set_size(2, false);
	p[0] = x1;
	p[1] = x2;
	return p;
}

void Rosenbrock::setParametersVector(const vec p)
{
	x1 = p[0];
	x2 = p[1];
}

double Rosenbrock::objective() const
{
	return 100 * pow(x2 - pow(x1,2),2) + pow((1.0 - x1),2);
}

vec Rosenbrock::gradient() const
{
	vec p;
	p.set_size(2, false);
	p[0] = -400.0 * (x2 - pow(x1, 2.0)) * x1 - 2.0 * (1 - x1);
	p[1] = 200.0 * (x2 - pow(x1, 2.0));

	return p;
}


