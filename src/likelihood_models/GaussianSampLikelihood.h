/***************************************************************************
 *   AstonGeostats, algorithms for low-rank geostatistical models          *
 *                                                                         *
 *   Copyright (C) Ben Ingram, 2008                                        *
 *                                                                         *
 *   Ben Ingram, IngramBR@Aston.ac.uk                                      *
 *   Neural Computing Research Group,                                      *
 *   Aston University,                                                     *
 *   Aston Street, Aston Triangle,                                         *
 *   Birmingham. B4 7ET.                                                   *
 *   United Kingdom                                                        *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef GAUSSIANSAMPLIKELIHOOD_H_
#define GAUSSIANSAMPLIKELIHOOD_H_

#include "LikelihoodType.h"
#include "SamplingLikelihood.h"

class GaussianSampLikelihood : public SamplingLikelihood
{
public:
	GaussianSampLikelihood(double Mean, double Variance);
	GaussianSampLikelihood(double Mean, double Variance, double (*transform)(double));
	virtual ~GaussianSampLikelihood();

	double updateCoefficients(double& K1, double& K2, double Observation, double ModelMean, double ModelVariance) const;
	void displayParameters() const;

protected:
	double likelihoodMean;
	double likelihoodVariance;
};
#endif /*GAUSSIANSAMPLIKELIHOOD_H_*/
