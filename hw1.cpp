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


void break_line(){
	//pause
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

	double ML = problem2_preshock.M[meshSize-1];
	//std::cout<<ML<<std::endl;
	double ML2 = std::pow(ML,2);
std::cout << ML2 << std::endl;
	double PL = problem2_preshock.Pressure[end];
std::cout<<"PL "<<PL<<std::endl;
	double m_inside = (2+(gamma-1)*ML2)/(2*gamma*ML2 - (gamma -1));
std::cout<<"M inside "<<m_inside<<std::endl;
	double MR = std::sqrt(m_inside);
std::cout<<MR<<std::endl;
	double PR = PL*(2*gamma*ML2 - (gamma -1))/(gamma+1);
std::cout<<PR<<std::endl;
	double p_top_inside = ((gamma+1)*ML2/2)/(1 + (gamma-1)*ML2/2);
std::cout<<p_top_inside<<std::endl;
	double p_top_exponent = gamma/(gamma-1);
std::cout<<p_top_exponent<<std::endl;
	double p_bottom_inside = 2*gamma*ML2/(gamma+1) - (gamma-1)/(gamma+1);
std::cout<<p_bottom_inside<<std::endl;
	double p_bottom_exponent = 1/(gamma-1);
std::cout<<p_bottom_exponent<<std::endl;
	double new_inletPressure = inletPressure* std::pow(p_top_inside, p_top_exponent)/std::pow(p_bottom_inside, p_bottom_exponent);
std::cout<<new_inletPressure<<std::endl;
	double L_total_density = inletPressure/R/totalTemp;
std::cout<<L_total_density<<std::endl;
	double R_total_density = new_inletPressure/R/totalTemp;
std::cout<<R_total_density<<std::endl;
	double a01 = std::sqrt(gamma*inletPressure/L_total_density);
std::cout<<a01<<std::endl;
	double a0R = std::sqrt(gamma*new_inletPressure/R_total_density);
std::cout<<a0R<<std::endl;
	double pLaL = L_total_density*a01*std::pow(2/(gamma+1),(gamma+1)/(2*(gamma-1)));
std::cout<<pLaL<<std::endl;
	double pRaR = R_total_density*a0R*std::pow(2/(gamma+1),(gamma+1)/(2*(gamma-1)));
std::cout<<pRaR<<std::endl;
	double new_S_star = s_star*pLaL/pRaR;
std::cout<<new_S_star<<std::endl;
	int meshSize1 = 200;
	double start1 = 7.;
	double end1 = 10.;
	double dx1 = (end1-start1)/meshSize;
	std::vector<double> mesh1(meshSize+1);
	//fill the vector mesh with explicit points
	for(int i=0; i < meshSize1+1; ++i){
		mesh1[i] = start1 + dx1*i;
	}

	std::cout<< new_S_star<<std::endl;

	Quasi1DFlow problem2_aftershock(mesh1, R, gamma, totalTemp, new_inletPressure, new_S_star, meshSize1, 0.2, 11, 0);

	problem2_aftershock.calculatePhysicalQuantities();
	problem2_aftershock.printMachTempDensityPressure("problem2_after_shock");

	//std::vector<double> M = problem2_preshock.M.push_back(problem2_aftershock.M);
	//std::cout<< M.size<< std::endl;


	}

	
	

	return 0;
}
