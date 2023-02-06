#ifndef QUASI_1D_FLOW_H
#define QUASI_1D_FLOW_H

#include <cctype>
#include <iostream>
#include <math.h>
#include <vector>
#include <cmath>
#include <functional>
#include <fstream>
#include <stdio.h>
#include <cassert>

#define TOL 1e-14

class Quasi1DFlow {
public:

	double S_star, gamma, T_01, P_01, R, initial_mach_guess, shock_location;
	double start = 0.;
	double end = 10.;
	int mesh_size;
	bool shock_present;
	std::vector<double> M, S, X, Temperature, Pressure, Density;

	
	std::vector<double> channelWidthP1(std::vector<double> x);
	double nonlinearFunctionToSolveP1(double M, const int i);
	double nonlinearFunctionToSolveP1Deriv(double M, const int i);
	double newtonSolveP1(double guess, const int index);
	void calculatePhysicalQuantities();
	void calculateTemp();
	void calculateDensity();
	void calculateMach();
	void calculatePressure();
	double reCalculateValuesAfterShock(const int i);

	void printMachTempDensityPressure(std::string a_file_name);

	Quasi1DFlow(std::vector<double> x, double r, double gamm, double totalTemp, double inletPressure, double s_star, int meshSize, double guessMach, double s_loc, bool sp)
		{
			//initialize all problem vectors to x, initialize s to shape function
			M=x;
			S=Quasi1DFlow::channelWidthP1(x);
			X=x;
			Pressure = x;
			Temperature = x;
			Density = x;
			//initialize problem parameters
			R = r;
			T_01 = totalTemp;
			P_01 = inletPressure;
			initial_mach_guess = guessMach;
			gamma = gamm;
			mesh_size = meshSize;
			S_star = s_star;
			shock_location = s_loc;
			shock_present = sp;
		}
};

std::vector<double> Quasi1DFlow::channelWidthP1(std::vector<double> x){

	int size = x.size();
	std::vector<double> S(size);
	for (int i = 0; i < size; ++i)
	{
		if(x[i] <= 5){
			S[i] = 1 + 1.5*pow(1-x[i]/5,2);
		}
		else if(x[i] >5){
			S[i] = 1 + 0.5*pow(1-x[i]/5,2);
		}
	}
	return S;
}

double Quasi1DFlow::nonlinearFunctionToSolveP1(double M, const int i){
	assert(S.size() > 0);
	double inside = 2/(gamma +1) + (gamma - 1)/(gamma+1)*pow(M,2);
	double exponent = (gamma + 1)/(2*(gamma-1));
	assert(fabs(exponent - 3)<TOL);
	return pow(inside, exponent) - (S[i]*M)/(S_star);
}

double Quasi1DFlow::nonlinearFunctionToSolveP1Deriv(double M, const int i){
	double inside = 2/(gamma +1) + (gamma - 1)/(gamma+1)*pow(M,2);
	double exponent = (gamma + 1)/(2*(gamma-1));
	return M*pow(inside, exponent-1) - (S[i])/(S_star);
}

double Quasi1DFlow::newtonSolveP1(double guess, const int index){
	double val = 1;
	int max_it = 10;
	int it = 1;
	while(!(std::fabs(val - 0)<TOL) && it < max_it){
		guess = guess - Quasi1DFlow::nonlinearFunctionToSolveP1(guess, index)/Quasi1DFlow::nonlinearFunctionToSolveP1Deriv(guess, index);
		val = Quasi1DFlow::nonlinearFunctionToSolveP1(guess, index);
		++it;
		}
	return guess;
}

void Quasi1DFlow::calculateMach(){
	std::vector<double> solution(X.size());
	double val = initial_mach_guess;
	
	if(!shock_present)//subsonic and no shock, use initial mach guess for all newton solves
	{
		for (int i = 0; i < mesh_size+1; ++i)
		{
			//TODO make this chack more robust, maybe a flag? what if shock is outside of the region <0?
			// if (X[i]<shock_location && initial_mach_guess<1)
			// {
			// 	val = 2; // before shock location we expect supersonic speed, so fix initial guess.
			// }
			M[i] = Quasi1DFlow::newtonSolveP1(val, i);
		}
	}
	else if (shock_present)//transonic and shock, use an increasing initial guess, jumping back to subsonic after shock with RH conditions
	{
		for (int i = 0; i < mesh_size+1; ++i)
		{
			//as we approach mach 1, we expect to go faster, while a typical solve will decrease in speed.
			if (X[i]>5)
			{
				val = 1.1; //use a supersonic guess to pass the 1 mach barrier.
			}

			M[i] = Quasi1DFlow::newtonSolveP1(val, i);
			val = M[i];
		}
	}
}

void Quasi1DFlow::calculateTemp(){
	assert(Temperature.size() == M.size());
	for (int i = 0; i<M.size();++i)
	{
		Temperature[i] = T_01/(1 + (gamma - 1)/(2)*std::pow(M[i],2));
	}
}

void Quasi1DFlow::calculateDensity(){
	assert(Density.size() == M.size());
	for (int i = 0; i<M.size();++i)
	{
		Density[i] = Pressure[i]/(R*Temperature[i]);
	}
}

void Quasi1DFlow::calculatePressure(){
	assert(Pressure.size() == M.size());
	for (int i = 0; i<M.size();++i)
	{
		double inside = (1 + (gamma - 1)/(2)*std::pow(M[i],2));
		double exponent = -gamma/(gamma-1);
		Pressure[i] = P_01*std::pow(inside,exponent);
	}
}

void Quasi1DFlow::calculatePhysicalQuantities(){
	Quasi1DFlow::calculateMach();
	Quasi1DFlow::calculateTemp();
	Quasi1DFlow::calculatePressure();
	Quasi1DFlow::calculateDensity();
}

double Quasi1DFlow::reCalculateValuesAfterShock(const int i){
//need to calculate new mach on right, pressure right, and pressure initial right

//then, we need to calculate S^* using the equation
//			Snew = Sold*A/B

// where A = 

	double G = std::pow(2/(gamma+1), (gamma + 1)/(2*(gamma-1)));
	double M_2 = std::pow(M[i-1],2);
	double M_r = std::sqrt((2 + (gamma-1)*M_2)/(2*gamma*M_2 - (gamma - 1)));
	M[i] = M_r;

	//double P_l = P_01*(1 + (gamma - 1)/2*M_2);
	//double P_r = (2*gamma*M_2 - gamma + 1)/(gamma+1)*P_l;



	return 0.0;
}

void Quasi1DFlow::printMachTempDensityPressure(std::string a_file_name){
	Quasi1DFlow::calculatePhysicalQuantities();

	std::ofstream a_file;
		a_file.open(a_file_name+ "_mach.csv", std::ios::out | std::ios::trunc);

			for (int i = 0; i< mesh_size+1; ++i)
			{
				a_file <<  " " << X[i]  << ", " << std::to_string(M[i]) << std::endl;
			}
			a_file.close();

	std::ofstream a_file1;
		a_file1.open(a_file_name+ "_temp.csv", std::ios::out | std::ios::trunc);

			for (int i = 0; i< mesh_size+1; ++i)
			{
				a_file1 <<  " " << X[i]  << ", " << std::to_string(Temperature[i]) << std::endl;
			}
			a_file1.close();

		std::ofstream a_file2;
			a_file2.open(a_file_name+ "_density.csv", std::ios::out | std::ios::trunc);

				for (int i = 0; i< mesh_size+1; ++i)
				{
					a_file2 <<  " " << X[i]  << ", " << std::to_string(Density[i]) << std::endl;
				}
				a_file2.close();

		std::ofstream a_file3;
			a_file3.open(a_file_name+ "_pressure.csv", std::ios::out | std::ios::trunc);

				for (int i = 0; i< mesh_size+1; ++i)
				{
					a_file3 <<  " " << X[i]  << ", " << std::to_string(Pressure[i]) << std::endl;
				}
				a_file3.close();
}



#endif


















