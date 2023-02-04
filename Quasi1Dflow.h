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

	double S_star, gamma, T_01, P_01, R, initial_mach_guess;
	std::vector<double> M, S, X;
	int mesh_size;
	double start = 0.;
	double end = 10.;
	double dx = (end-start)/mesh_size;
	void initializedata(std::vector<double> x, double R, double gamma, double totalTemp, double inletPressure, double s_star, int meshSize, double guessMach);
	std::vector<double> channelWidthP1(std::vector<double> x);
	double nonlinearFunctionToSolveP1(double M, const int i);
	double nonlinearFunctionToSolveP1Deriv(double M, const int i);
	double newtonSolveP1(double guess, const int index);
	void calculateTemp();
	void calculateDensity();
	void calculateMach();

	void printMachTempDensity(std::string a_file_name);

	Quasi1DFlow(std::vector<double> x, double r, double gamm, double totalTemp, double inletPressure, double s_star, int meshSize, double guessMach)
		{
				M=x;
				S=Quasi1DFlow::channelWidthP1(x);
				X=x;
				R = r;
				T_01 = totalTemp;
				P_01 = inletPressure;
				initial_mach_guess = guessMach;
				gamma = gamm;


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
	while(!(fabs(val - 0)<TOL) && it < max_it){
		guess = guess - Quasi1DFlow::nonlinearFunctionToSolveP1(guess, index)/Quasi1DFlow::nonlinearFunctionToSolveP1Deriv(guess, index);
		val = Quasi1DFlow::nonlinearFunctionToSolveP1(guess, index);
		++it;
		}
	return guess;
}

void Quasi1DFlow::calculateMach(){
	std::vector<double> solution(X.size());\
	double val = initial_mach_guess;
	for (int i = 0; i < mesh_size+1; ++i)
		{
			solution[i] = Quasi1DFlow::newtonSolveP1(val, i);
			val = solution[i];
		}
}

void Quasi1DFlow::calculateTemp(){

}

void Quasi1DFlow::calculateDensity(){

}


void Quasi1DFlow::printMachTempDensity(std::string a_file_name){
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
			a_file1 <<  " " << X[i]  << ", " << std::to_string(M[i]) << std::endl;
		}
		a_file1.close();

		std::ofstream a_file2;
			a_file2.open(a_file_name+ "_density.csv", std::ios::out | std::ios::trunc);

			for (int i = 0; i< mesh_size+1; ++i)
			{
				a_file2 <<  " " << X[i]  << ", " << std::to_string(M[i]) << std::endl;
			}
			a_file2.close();

		std::ofstream a_file3;
					a_file3.open(a_file_name+ "_pressure.csv", std::ios::out | std::ios::trunc);

						for (int i = 0; i< mesh_size+1; ++i)
						{
							a_file3 <<  " " << X[i]  << ", " << std::to_string(M[i]) << std::endl;
						}
						a_file3.close();
}



#endif


















