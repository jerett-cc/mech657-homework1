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


struct problemData {

	double S_star, gamma, T_01, P_01, R;
	std::vector<double> M, S, X;
	int mesh_size;

	problemData(std::vector<double> s, std::vector<double> x, double s_star, double gamm, double t_01, double p_01, double r, int meshSize){
		T_01 = t_01;
		P_01 = p_01;
		X = x;
		S_star = s_star;
		gamma = gamm;
		M = x;
		R = r;
		mesh_size = meshSize;
        S = s;
	}
};

std::vector<double> S1(std::vector<double> x){

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

double f(const problemData &data, double M, const int i){

	double inside = 2/(data.gamma +1) + (data.gamma - 1)/(data.gamma+1)*pow(M,2);
	double exponent = (data.gamma + 1)/(2*(data.gamma-1));

	assert(fabs(exponent - 3)<TOL);

	return pow(inside, exponent) - (data.S[i]*M)/(data.S_star);
}

double fprime(const problemData &data, double M, const int i){
	
	double inside = 2/(data.gamma +1) + (data.gamma - 1)/(data.gamma+1)*pow(M,2);
	double exponent = (data.gamma + 1)/(2*(data.gamma-1));
		return M*pow(inside, exponent-1) - (data.S[i])/(data.S_star);
}


double newtonSolve1(problemData &data, double guess, const int index){

//
	double val = 1;
	int max_it = 10;
	int it = 1;
	while(!(fabs(val - 0)<TOL) && it < max_it){
		guess = guess - f(data, guess, index)/fprime(data, guess, index);
		val = f(data, guess, index);
		++it;
	}
	return guess;
}

void calculateTemperature(){



}


void printResults(problemData data, std::string a_file_name){
	std::ofstream a_file;
	a_file.open(a_file_name, std::ios::out | std::ios::trunc);

	for (int i = 0; i< data.mesh_size+1; ++i)
	{
		a_file <<  " " << data.X[i]  << ", " << std::to_string(data.M[i]) << std::endl;
	}
	a_file.close();
}

int main() {

	int meshSize = 300;
	double start = 0.;
	double end = 10.;

	double dx = (end-start)/meshSize;

	//problem parameters
	double R = 287.;
	double gamma = 1.4;
	double totalTemp = 300.;
	double inletPressure = 100.;
	double s_star = 0.8;



	std::vector<double> mesh(meshSize+1);
	std::vector<double> solution(mesh.size());

	//fill the vector mesh with explicit points
	for(int i=0; i < meshSize+1; ++i){
		mesh[i] = start + dx*i;
	}

	Quasi1DFlow problem1(mesh, R, gamma, totalTemp, inletPressure, s_star, meshSize, 0.2);

	problem1.calculateMach();
	problem1.printMachTempDensity("new solution.txt");

//	std::vector<double> S = S1(mesh);
//	assert(fabs(S[meshSize]-1.5)<TOL);
//
//	problemData data1(S, mesh, s_star, gamma, totalTemp, inletPressure, R, meshSize);
//
//	double initial_guess1 = 0.2;
//
//	//loop through all points and generate the solution with newton method
//	for (int i = 0; i < meshSize+1; ++i)
//	{
//		solution[i] = newtonSolve1(data1, initial_guess1, i);
//		initial_guess1 = solution[i];
//	}
//
//	data1.M = solution;
//	assert((data1.M[3] - solution[3])< TOL);
//
//	printResults(data1, "problem1_mach.txt");

	return 0;
}
