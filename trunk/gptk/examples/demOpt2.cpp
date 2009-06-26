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
// Example of finding the root of a generalized Rosenbrock function using
// different model trainers.  Use numerical gradients.
//
// History:
//
// 19 Dec 2008 - BRI - First implementation. 
//

#include <iostream>

#include "optimisation/Optimisable.h"
#include "optimisation/ModelTrainer.h"
#include "optimisation/SCGModelTrainer.h"

#include "RosenbrockND.h"

using namespace std;

int main()
{
	// how many iterations of the algorithm should we run
	int iterations = 2000;
	
	// setup output precision
   cout.setf(ios::fixed);
   cout.precision(4);
   
   // create a Rosenbrock 500D object with default parameters
   // the starting value for this is automatically created in the constructor
	RosenbrockND rosen(1500);
	
	bvec optMask(5);
	// create a SCG model trainers
	SCGModelTrainer scgTrainer(rosen);

//	optMask(0) = true;
//	optMask(1) = true;
//	optMask(2) = true;
//	optMask(3) = false;	
//	optMask(4) = true;
//	optMask(5) = false;
//	optMask(6) = false;
//	optMask(7) = true;	
//	optMask(8) = true;
//	optMask(9) = true;

//	vec a = rosen.getParametersVector();
//	a(4) = 1.0;
//	a(3) = 1.0;
//	rosen.setParametersVector(a);


//	scgTrainer.setOptimisationMask(optMask);

	cout << "Parameters Before: " << rosen.getParametersVector() << endl;	
	cout << "Scaled Conjugate Gradient..." << endl;
	scgTrainer.setAnalyticGradients(false);
	scgTrainer.Train(iterations);
	cout << "Parameters After : " << rosen.getParametersVector() << endl;	
	scgTrainer.Summary();

	cout << "Finishing..." << endl;
	
}
