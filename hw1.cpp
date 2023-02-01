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
#include <functional>
#include <fstream>
#include <stdio.h>
#include <cassert>

#define TOL = 1e-16;

using namespace std;


/*Outline of steps:
 * step1) set up all of the requisite variables needed for computation++
 * step2) write a bisection method to solve for M++
 * step3) write a S(x) function++
 * step4) generate a vector of points to evaluate at++
 * step5) for every point, calculate M
 * step5.1) recalculate the downstream S* and p0
 * step6) store the corresponding M in a vector
 * step7) print the vector x,y in a file
 *
 *
 */

struct problemData {

	double S_star, gamma, temp, pressure, T_01, P_01, R;
	vector<double> M, S, X;
	int mesh_size;

	problemData(vector<double> x, double s_star, double gamm, double t_01, double p_01, double r, int meshSize){
		T_01 = t_01;
		P_01 = p_01;
		X = x;
		S_star = s_star;
		gamma = gamm;
		M = x;
		R = r;
		mesh_size = meshSize;
	}


};

vector<double> S1(vector<double> x){

	int size = x.size();
	vector<double> S(size);

	for(int i = 0; i < size; ++i){

		if(x[i] <= 5){
			S[i] = 1 + 1.5*pow(1-x[i]/5,2);
		}
		else if(x[i] >5){
			S[i] = 1 + 0.5*pow(1-x[i]/5,2);
		}

	}

	return S;

}

//make proper return type
double f(problemData data, double M, int index){
	double inside = 2/(data.gamma + 1) * (1 + (data.gamma - 1)/2*pow(M,2));
	double exponent = (data.gamma + 1)/(2*(data.gamma - 1));
	double val = 1/M*pow(inside,exponent) - data.S[index]/data.S_star;

	return val;

}

double fprime(problemData data, double M, int index){

	double inside = 2/(data.gamma + 1) * (1 + (data.gamma - 1)/2*pow(M, 2));
	double exponent = (data.gamma + 1)/(2*(data.gamma - 1));

	double val = -1/pow(M,2)*pow(inside, exponent) + 1/M*exponent*pow(inside, exponent - 1);

	return val;

}



double newtonSolve1(problemData data, double guess, int index){

	while(f(data, guess, index)> 1e-16){

		guess = guess + f(data, guess, index)/fprime(data, guess, index);

	}

	return guess;

}

void calculateTemperature(){



}

void printResults(problemData data, string a_file_name){


	std::ofstream a_file;
	a_file.open(a_file_name, std::ios::out | std::ios::trunc);



	for (int i = 0; i< data.mesh_size; ++i)
	{
		a_file <<  " " << data.X[i]  << "   " << std::to_string(data.M[i]) << std::endl;
	}
	a_file.close();
}



int main() {

	//generate the x vector we will use for this problem

	int meshSize = 3000;
	double start = 0.;
	double end = 10.;

	double dx = (end-start)/meshSize;

	//problem parameters
	double R = 287.;
	double gamma = 1.4;
	double totalTemp = 300.;
	double inletPressure = 100.;
	double s_star = 0.8;


	vector<double> S(meshSize);//used for newton method
	vector<double> mesh(meshSize);

	//fill the vector mesh with explicit points
	for(int i=0; i < meshSize; ++i){
		S[i] = 0;
		mesh[i] = start + dx*i;
	}



	cout << mesh[0] << " " << mesh[1500] << " " << mesh[meshSize-1] << "\n";

	problemData data1(mesh, s_star, gamma, totalTemp, inletPressure, R, meshSize);
	data1.S = S1(data1.X);

	double initial_guess = 0.;

	for(int i = 0; i < meshSize; ++i){

		data1.M[i] = newtonSolve1(data1, initial_guess, i);

		initial_guess = data1.M[i];

	}

	//printResults(data1, "problem1_results.txt");


	return 0;
}
