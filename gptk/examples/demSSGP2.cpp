/***************************************************************************
 *   AstonGeostats, algorithms for low-rank geostatistical models          *
 *                                                                         *
 *   Copyright (C) Ben Ingram, 2008-2009                                   *
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



#include <iostream>
#include <cassert>
#include <ctime>

#include "itpp/itbase.h"
#include "itpp/itstat.h"

#include "optimisation/Optimisable.h"
#include "optimisation/ModelTrainer.h"
#include "optimisation/SCGModelTrainer.h"
#include "optimisation/CGModelTrainer.h"
#include "optimisation/GDModelTrainer.h"
#include "optimisation/QuasiNewtonModelTrainer.h"

#include "gaussianProcesses/ForwardModel.h"
#include "optimisation/Optimisable.h"
#include "gaussianProcesses/SequentialGP.h"
#include "gaussianProcesses/GaussianProcess.h"

#include "likelihoodModels/GaussianLikelihood.h"
#include "likelihoodModels/GaussianSampLikelihoodMathML.h"

#include "covarianceFunctions/CovarianceFunction.h"
#include "covarianceFunctions/SumCovarianceFunction.h"
#include "covarianceFunctions/GaussianCF.h"
//#include "covarianceFunctions/MaternCF.h"
#include "covarianceFunctions/ExponentialCF.h"
#include "covarianceFunctions/WhiteNoiseCF.h"
#include "covarianceFunctions/ConstantCF.h"
#include "covarianceFunctions/LinearCF.h"


#include "covarianceFunctions/NegLogSigmoidTransform.h"
#include "covarianceFunctions/SigmoidTransform.h"

#include "GraphPlotter/GraphPlotter.h"

#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>


using namespace std;
using namespace itpp;

/**
 * Returns a matrix read from a comma separated values (CSV) file
 */
//int readCSV(mat &data, string filename) {
//
//	std::fstream predfile("noisysinePred.csv", ios::in);
//	
//	while(getline(predfile, csvLine) && (i < predobs))
//	{
//		string csvElement;
//		istringstream csvStream(csvLine);
//
//		getline(csvStream, csvElement, ',');
//		pX(i,0) = strtod(csvElement.c_str(),NULL);
//
//		i++;
//	}
//}

int main()
{
	const int obs = 100, predobs = 400;
	mat X = zeros(obs,1);
	vec y = zeros(obs);
	mat pX = zeros(predobs, 1);
	vec meanP = zeros(predobs);
	vec varP = zeros(predobs);

	int i = 0;

	ivec likIdx(obs);

	cout.setf(ios::fixed);
   	cout.precision(4);

   	// **
   	// Read in noisy sine data
   	// **
	double valRead;
	string csvLine;
	std::fstream inputfile("noisysine.csv", ios::in);
	while(getline(inputfile, csvLine) && (i < obs))
	{
		string csvElement;
		istringstream csvStream(csvLine);

		getline(csvStream, csvElement, ',');
		X(i,0) = strtod(csvElement.c_str(),NULL);

//		getline(csvStream, csvElement, ',');
//		X(i,1) = strtod(csvElement.c_str(),NULL);
		
		getline(csvStream, csvElement, ',');
		valRead = strtod(csvElement.c_str(),NULL);

		if(i < 50)
		{
			valRead = valRead * 0.25;
			likIdx(i) = 2;
		}
		else
		{
			valRead = (valRead * 0.3) + 5.0;

			cout << valRead << endl;

			likIdx(i) = 1;
		}


		y(i) = valRead;
//		getline(csvStream, csvElement, ',');

		i++;
	}

	cout << y << endl;

	i = 0;
	std::fstream predfile("noisysinePred.csv", ios::in);
	while(getline(predfile, csvLine) && (i < predobs))
	{
		string csvElement;
		istringstream csvStream(csvLine);

		getline(csvStream, csvElement, ',');
		pX(i,0) = strtod(csvElement.c_str(),NULL);

//		getline(csvStream, csvElement, ',');
//		pX(i,1) = strtod(csvElement.c_str(),NULL);

		i++;
	}


	
	cout <<endl;

	cout << "===============================================" << endl;
	
//	double range = 0.20 * (min(max(X, 1) - min(X, 1)));
//	double sill = variance(y);
//	double nugget = 0.4 * sill;

	// double range = 0.9574;
	double range = 0.1;
	double sill = 0.66;
	double nugget = 0.05;

	GaussianCF g1(range, sill);
	WhiteNoiseCF g2(nugget);
//	ConstantCF g3(10000000.0);
	SumCovarianceFunction gCmp(g1);
	gCmp.addCovarianceFunction(g2);
//	gCmp.addCovarianceFunction(g3);	

	gCmp.displayCovarianceParameters();

	cout << "Setting up SSGP" << endl;
	SequentialGP ssgp(1, 1, 40, X, y, gCmp);

	bvec optMask(3);
	optMask(0) = true;
	optMask(1) = true;
	optMask(2) = true;
//	optMask(3) = true;	
	SCGModelTrainer gpTrainer(ssgp);
	gpTrainer.setOptimisationMask(optMask);


	cout << "Choosing a likelihood" << endl;
//	GaussianLikelihood gaussLik(0.001);

	cout << "likelihood 1: "<< endl;
	GaussianSampLikelihoodMathML gaussLik1(0.0, 0.005, "<apply><plus/><apply><times/><ci>x</ci><cn>0.3</cn></apply><cn>5</cn></apply>");


	cout << "likelihood 2: "<< endl;
	GaussianSampLikelihoodMathML gaussLik2(0.0, 0.005, "<apply><times/><ci>x</ci><cn>0.25</cn></apply>");


	cout << "vector of likelihoods: "<< endl;
	Vec<LikelihoodType *> likMod(2);
	likMod(0) = &gaussLik1;
	likMod(1) = &gaussLik2;

	gpTrainer.setAnalyticGradients(false);

	ssgp.setLikelihoodType(UpperBound);
	

	
	// **
	// 5 iterations of data recycling
	// **
	for(int i=1; i < 5; i++)
	{
		cout << "Data recycling : #" << i << endl;

		cout << "Computing posterior" << endl;
//		ssgp.computePosterior(gaussLik);
		ssgp.computePosterior(likIdx, likMod);

		cout << "Optimising parameters" << endl;
		gpTrainer.Train(5);

		cout << "Recomputing posterior" << endl;
		ssgp.recomputePosterior();

		cout << "Resetting posterior" << endl;
		ssgp.resetPosterior();

		gCmp.displayCovarianceParameters();
	}
/*
//	ssgp.computePosterior(gaussLik);
	ssgp.computePosterior(likIdx, likMod);

cout << "-------------------" << endl;	
		gCmp.displayCovarianceParameters();
		ssgp.displayModelParameters();
cout << "===================" << endl;



	ssgp.makePredictions(meanP, varP, pX);

	vec xTest = pX.get_col(0);
	vec xTrain = X.get_col(0);

	vec upperEB = meanP + sqrt(varP) * 2;
	vec lowerEB = meanP - sqrt(varP) * 2;

	GraphPlotter gplot = GraphPlotter();
	gplot.clearPlot();
	gplot.plotPoints(xTest, meanP, "test", LINE, BLUE);
	gplot.plotPoints(xTrain, y, "train", CROSS, RED);
	gplot.plotPoints(xTest, upperEB, "error bar", LINE, CYAN);
	gplot.plotPoints(xTest, lowerEB, "error bar", LINE, CYAN);
*/

//cout << "yTest = " << meanP << endl;
//cout << "xTest = " << pX << endl;
//cout << "xTrain = " << X << endl;
//cout << "yTrain = " << y << endl;

}

