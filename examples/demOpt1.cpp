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

//
// Optimisation Example 1
// Example of finding the root of the Rosenbrock function using different
// model trainers
//
// History:
//
// 12 Dec 2008 - BRI - First implementation. 
//
 
#include <iostream>

#include "optimisation/Optimisable.h"
#include "optimisation/ModelTrainer.h"

#include "optimisation/SCGModelTrainer.h"
#include "optimisation/GDModelTrainer.h"
#include "optimisation/CGModelTrainer.h"
#include "optimisation/QuasiNewtonModelTrainer.h"

#include "Rosenbrock.h"

using namespace std;

int main()
{
	// how many iterations of the algorithm should we run
	int iterations = 30;
	
	// setup output precision
   cout.setf(ios::fixed);
   cout.precision(4);
   
   // create a Rosenbrock object with default parameters
	Rosenbrock rosen(-1.0, 1.0);
	vec startingValue = rosen.getParametersVector();

	// create 4 different model trainers
	SCGModelTrainer scgTrainer(rosen);
	GDModelTrainer gdTrainer(rosen);
	CGModelTrainer cgTrainer(rosen);
	QuasiNewtonModelTrainer qnTrainer(rosen);
 
	cout << "Quasi-Newton" << endl;
	qnTrainer.setAnalyticGradients(true);
	qnTrainer.Train(iterations);
	qnTrainer.Summary();
	cout << "Final Parameters: " << rosen.getParametersVector() << endl << endl;

	rosen.setParametersVector(startingValue);
	cout << "Conjugate Gradient" << endl;
	cgTrainer.setAnalyticGradients(true);
	cgTrainer.Train(iterations);
	cgTrainer.Summary();
	cout << "Final Parameters: " << rosen.getParametersVector() << endl << endl;	

	rosen.setParametersVector(startingValue);
	cout << "Scaled Conjugate Gradient" << endl;
	scgTrainer.setAnalyticGradients(true);
	scgTrainer.Train(iterations);
	scgTrainer.Summary();
	cout << "Final Parameters: " << rosen.getParametersVector() << endl << endl;	

	rosen.setParametersVector(startingValue);
 	cout << "Gradient Descent" << endl;
	gdTrainer.setAnalyticGradients(true);
	gdTrainer.Train(iterations);
	gdTrainer.Summary();
	cout << "Final Parameters: " << rosen.getParametersVector() << endl << endl;	



	
}
