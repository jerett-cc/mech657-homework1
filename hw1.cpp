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
#include "Quasi1DShockTube.h"


void recalculateAfterShock(const Quasi1DFlow &preShock){
	double ML = preShock.M[preShock.mesh_size-1];
		double gamma = preShock.gamma;
		double ML2 = std::pow(ML,2);
//	std::cout << ML2 << std::endl;
		double PL = preShock.Pressure[preShock.mesh_size-1];
		double m_inside = (2+(gamma-1)*ML2)/(2*gamma*ML2 - (gamma -1));
		double MR = std::sqrt(m_inside);
		double PR = PL*(2*gamma*ML2 - (gamma -1))/(gamma+1);
		double p_top_inside = ((gamma+1)*ML2/2)/(1 + (gamma-1)*ML2/2);
		double p_top_exponent = gamma/(gamma-1);
		double p_bottom_inside = 2*gamma*ML2/(gamma+1) - (gamma-1)/(gamma+1);
		double p_bottom_exponent = 1/(gamma-1);
		double new_inletPressure = preShock.P_01* std::pow(p_top_inside, p_top_exponent)/std::pow(p_bottom_inside, p_bottom_exponent);
		double L_total_density = preShock.P_01/preShock.R/preShock.T_01;
		double R_total_density = new_inletPressure/preShock.R/preShock.T_01;
		double a01 = std::sqrt(gamma*preShock.P_01/L_total_density);
		double a0R = std::sqrt(gamma*new_inletPressure/R_total_density);
		double pLaL = L_total_density*a01*std::pow(2/(gamma+1),(gamma+1)/(2*(gamma-1)));
		double pRaR = R_total_density*a0R*std::pow(2/(gamma+1),(gamma+1)/(2*(gamma-1)));
		double new_S_star = preShock.S_star*pLaL/pRaR;
		int meshSize1 = 200;
		double start1 = 7.;
		double end1 = 10.;
		double dx1 = (end1-start1)/preShock.mesh_size;
		std::vector<double> mesh1(preShock.mesh_size+1);
		//fill the vector mesh with explicit points
		for(int i=0; i < meshSize1+1; ++i){
			mesh1[i] = start1 + dx1*i;
		}

//		std::cout<< new_S_star<<std::endl;

		Quasi1DFlow problem2_aftershock(mesh1, preShock.R, gamma, preShock.T_01 , new_inletPressure, new_S_star, meshSize1, 0.2, 11, 0);

		problem2_aftershock.calculatePhysicalQuantities();
		problem2_aftershock.printMachTempDensityPressure("problem2_after_shock");

}

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
	double inletPressure = 100000.;
	double s_star = 0.8;
	double shock = 11;//shock oustide of domain...

	Quasi1DFlow problem1(mesh, R, gamma, totalTemp, inletPressure, s_star, meshSize, 0.2, shock, 0);

	problem1.calculateMach();
	problem1.calculatePhysicalQuantities();
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
	double inletPressure = 100000.;
	double s_star = 1.;
	double shock = 7.;

	Quasi1DFlow problem2_preshock(mesh, R, gamma, totalTemp, inletPressure, s_star, meshSize, 1.2, shock, 1);
	problem2_preshock.calculateMach();
	problem2_preshock.calculatePhysicalQuantities();
	problem2_preshock.printMachTempDensityPressure("problem2_before_shock");

	//TODO. recalculate initial variables for after the shock, and create a second quasiflow with these initial values and a mesh 3-10. 
	//then in matlab, you'll need to append the results from each problem.
	//potential problem: when you might need to compare to another algorithm, it may be tough to generate the same points.

	recalculateAfterShock(problem2_preshock);

	}

	//problem 3
	{
	int meshSize = 2000;
	double start = 0.;
	double end = 10.;
	double dx = (end-start)/meshSize;
	std::vector<double> mesh(meshSize+1);
	//fill the vector mesh with explicit points
	for(int i=0; i < meshSize+1; ++i){
		mesh[i] = start + dx*i;
	}

	double PL = 1e5;
	double densityL = 1;
	double PR = 1e4;
	double densityR = 0.125;
	double membraneLocation = 5.;
	double gamma = 1.4;
	double time = 6.1/1000.;

	std::cout << "Starting pressure Left: " <<PL << std::endl;
	Quasi1DShockTube shockTube;

		shockTube.runSolution(mesh, PL, densityL, PR, densityR, gamma, time, 1., membraneLocation);
		shockTube.calculatePressure();
		shockTube.calculateDensity();
		shockTube.calculateMach();
		shockTube.printMachTempDensityPressure("problem3");
	}
	
	

	return 0;
}
