/*
 * GaussianCF.h
 *
 *  Created on: 28 Oct 2009
 *      Author: barillrl
 */

#ifndef GAUSSIANCF_H_
#define GAUSSIANCF_H_

#include "StationaryCF.h"

class GaussianCF: public StationaryCF {
public:
    GaussianCF(double variance, double lengthScale);
    virtual ~GaussianCF();

protected:
    virtual double correlation(double sqDist) const;
    virtual double correlation_gradient(int parameterNumber, double sqDist) const;
};

#endif /* GAUSSIANCF_H_ */
