//============================================================================
// Name        : hw1.cpp
// Author      : Jerett Cherry
// Version     :
// Copyright   : Your copyright notice
// Description : Generate reference solutions to two problems.
//============================================================================

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;


/*Outline of steps:
 * step1) set up all of the requisite variables needed for computation
 * step2) write a bisection method to solve for M
 * step3) write a S(x) function
 * step4) generate a vector of points to evaluate at
 * step5) for every point, calculate M
 * step5.1) recalculate the downstream S* and p0
 * step6) store the corresponding M in a vector
 * step7) print the vector x,y in a file
 *
 *
 */

struct problem1Data {

	double S, S_star, gamma, tolerance, temp, pressure, T_01, P_01, M;

};

void bisection(double S, double S_star, double gamma, double tolerance){

	// complete the bisection method for the equation
	// 0 = -S/S_star +1/M*(2/(gamma + 1))*pow((1+((gamma-1)/2)M^2), (gamma + 1)/(2(gamma-1)))

}

void temperature(){

}


int main() {
	problem1Data* data1;
	data1->S_star = 0.8;
	data1->P_01 = 100.;
	data1->gamma = 1.4;

	return 0;
}
