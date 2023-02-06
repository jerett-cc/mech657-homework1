//============================================================================
// Name        : hw1.cpp
// Author      : Jerett Cherry
// Version     :
// Copyright   : Your copyright notice
// Description : Generate reference solutions to two problems.
//============================================================================

#include <cctype>
#include <iostream>
#include <math.h>
#include <vector>
#include <cmath>
#include <functional>
#include <fstream>
#include <stdio.h>
#include <cassert>

#include "Quasi1Dflow.h"


int main() {

	{
	int meshSize = 300;
	double start = 0.;
	double end = 10.;
	double dx = (end-start)/meshSize;
	std::vector<double> mesh(meshSize+1);
	//fill the vector mesh with explicit points
	for(int i=0; i < meshSize+1; ++i){
		mesh[i] = start + dx*i;
	}
	//problem1 parameters
	double R = 287.;
	double gamma = 1.4;
	double totalTemp = 300.;
	double inletPressure = 100.;
	double s_star = 0.8;
	double shock = 11;//shock oustide of domain...

	Quasi1DFlow problem1(mesh, R, gamma, totalTemp, inletPressure, s_star, meshSize, 0.2, 11, 0);

	problem1.calculateMach();
	problem1.printMachTempDensityPressure("problem1");
	}

	{
	int meshSize = 200;
	double start = 0.;
	double end = 7.;
	double dx = (end-start)/meshSize;
	std::vector<double> mesh(meshSize+1);
	//fill the vector mesh with explicit points
	for(int i=0; i < meshSize+1; ++i){
		mesh[i] = start + dx*i;
	}
	//problem2 parameters
	double R = 287.;
	double gamma = 1.4;
	double totalTemp = 300.;
	double inletPressure = 100.;
	double s_star = 1.;
	double shock = 7.;

	Quasi1DFlow problem2_preshock(mesh, R, gamma, totalTemp, inletPressure, s_star, meshSize, 1.2, shock, 1);
	problem2_preshock.calculateMach();
	problem2_preshock.printMachTempDensityPressure("problem2");

	//TODO. recalculate initial variables for after the shock, and create a second quasiflow with these initial values and a mesh 3-10. 
	//then in matlab, you'll need to append the results from each problem.
	//potential problem: when you might need to compare to another algorithm, it may be tough to generate the same points.
	}
	return 0;
}
