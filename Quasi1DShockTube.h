/*
 * Quasi1DShockTube.h
 *
 * class to solve a shock tube problem reference Fundamental Algorithms in Computational Fluid Dynamics by Pulliam and Zingg
 * exercise 3.3
 *
 *  Created on: Feb 7, 2023
 *      Author: Jerett Cherry
 */

#ifndef QUASI1DSHOCKTUBE_H_
#define QUASI1DSHOCKTUBE_H_


#include <cctype>
#include <iostream>
#include <math.h>
#include <vector>
#include <cmath>
#include <functional>
#include <fstream>
#include <stdio.h>
#include <cassert>

class Quasi1DShockTube {
	std::vector<double> X, Mach, Temperature, Pressure, Density;
	double shock_location, contact_surface_location, head_of_expansion_fan_location;
	double AL, AR, pressure_ratio_across_shock, pressure_L, pressure_R, R, total_temperature;
};






#endif /* QUASI1DSHOCKTUBE_H_ */
