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

#define TOL 1e-10

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


	double S_star, gamma, T_01, P_01, R;
	vector<double> M, S, X;
	int mesh_size;

	problemData(vector<double> s, vector<double> x, double s_star, double gamm, double t_01, double p_01, double r, int meshSize){
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
double f1(problemData &data, double M, int index){
	
	//cout << "M in f:" << M << endl;
	assert(fabs(data.gamma - 1.4)<TOL);
	double a = 2/(data.gamma + 1);
	double b = 1 + (data.gamma - 1)/2*pow(M,2);
	double c = (data.gamma +1)/(2*(data.gamma -1));
	//cout<< a<< " "<<b<<" "<< c-3<<endl;
	assert(fabs(c - 3.)<TOL);
	double inside = a * b;
	assert(fabs(data.S_star - 0.8)<TOL);
	double val = 1/M*pow(inside,c) - data.S[index]/data.S_star;
	//cout << "f = " << val << "\n";
	return val;

}

double f1prime(problemData &data, double M, int index){
	double a = 2/(data.gamma + 1);
	double b = 1 + (data.gamma - 1)/2*pow(M,2);
	double c = (data.gamma +1)/(2*(data.gamma -1));

	double inside = a*b;
	
	double val = -1/pow(M,2)*pow(inside, c) + 2*pow(inside, c - 1);

	return val;
}


double newtonSolve1(problemData &data, double guess, const int index){

//
	double f = 1;
	int max_it = 100;
	int it = 1;
	while(!(fabs(f - 0)<TOL) && it < max_it){
		//cout << "Guess sent to f " << guess << " at position " << index <<  endl;
		//data.M[index] = guess;
		guess = guess - f1(data, guess, index)/f1prime(data, guess, index);
		f = f1(data, guess, index);
		cout << "position" << index << "failed to converge" << endl;
		++it;
	}

	//for(int i = 0; i < 1000; ++i){
			//cout << "Guess sent to f " << guess << " at position " << index <<  endl;
			//data.M[index] = guess;
	//		guess = guess - f1(data, guess, index)/f1prime(data, guess, index);
	//		f = f1(data, guess, index);
			//cout << "value of f = " << f << endl;
	//	}

	return guess;

}

void calculateTemperature(){



}


double bisectionMtd(problemData &data, double a, double b, const int index){

	assert(a < b);
	//assert(f1(data, a, index)*f1(data, b, index)<0);

	double midpoint = (a + b)/2;

	while(fabs(f1(data, midpoint, index))<TOL){

		if(f1(data, midpoint, index) * f1(data, a, index)<0){
			b = midpoint;
		}
		else if(f1(data, midpoint, index) * f1(data, b, index)<0){
			a = midpoint;
		}
		midpoint = (a + b)/2;
	}
	return midpoint;
}

void printResults(problemData data, string a_file_name){


	std::ofstream a_file;
	a_file.open(a_file_name, std::ios::out | std::ios::trunc);


	for (int i = 0; i< data.mesh_size+1; ++i)
	{
		a_file <<  " " << data.X[i]  << ", " << std::to_string(data.M[i]) << std::endl;
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
    	

	vector<double> mesh(meshSize+1);
	vector<double> solution(mesh.size());

	//fill the vector mesh with explicit points
	for(int i=0; i < meshSize+1; ++i){
		mesh[i] = start + dx*i;
	}

	vector<double> S = S1(mesh);
	assert(fabs(S[meshSize]-1.5)<TOL);

	//cout << mesh[0] << " " << mesh[1500] << " " << mesh[meshSize-1] << "\n";

	problemData data1(S, mesh, s_star, gamma, totalTemp, inletPressure, R, meshSize);
	problemData* dataPtr = &data1;

	double initial_guess1 = 0.2;
	double initial_guess2 = 0.1;

	//loop through all points and generate the solution with 
	for(int i = 0; i < meshSize+1; ++i){
		cout << "Solving at: "<< i << endl;

		//uncomment below to solve with a newton method
		solution[i] = newtonSolve1(data1, initial_guess1, i);
		initial_guess1 = solution[i];

		//uncomment below to produce the bisection method
		//solution[i] = bisectionMtd(data1, initial_guess2, initial_guess1, i);
		//initial_guess1 = initial_guess1 + (-0.35/5 * fabs(  data1.X[i]-5 ) +0.55 );
		//initial_guess2 = initial_guess2 + (-0.35/5 * fabs(  data1.X[i]-5 ) +0.55 );
		//
	}

	data1.M = solution;
	assert((data1.M[3] - solution[3])< TOL);

	printResults(data1, "problem1_results.txt");
	
	std::ofstream a_file;
          a_file.open("problem1 S vec.txt", std::ios::out | std::ios::trunc);
  
  
          for (int i = 0; i< data1.mesh_size+1; ++i)
          {
                  a_file <<  " " << data1.X[i]  << "   " << std::to_string(data1.S[i]) << std::endl;
          }
          a_file.close();

	return 0;
}
