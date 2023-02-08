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

#define TOL 1e-14
//TODO: need to initialize the vectors all to X.
class Quasi1DShockTube {
	public:
		std::vector<double> X, Mach, Temperature, Pressure, Density;
		double initial_shock_location, head_of_expansion_fan_location, tail_of_expansion_fan_location;
		double AL, AR, pressure_ratio_across_shock, pressure_L, pressure_R, R, total_temperature, gamma;
		double density_L , density_R, time;
		double density_left_of_shock, density_left_of_contact, pressure_left_of_shock, pressure_left_of_contact;
		double shock_velocity, contact_velocity, shock_position, contact_position;
		double fluid_velocity_left_of_shock, fluid_velocity_left_of_contact;

		double functionToSolve(double P);
		double derivativeFunctionToSolve(double P);

		void setInitialStates(double pR, double pL, double densityR, double densityL, double gamm, double t, double x0);
		void calculateSoundSpeedsLR();
		void newtonSolvePressureRatio(double initial);

		void calculateDensityLeftOfShock();
		void calculatePressureLeftOfContact();
		void calculatePressureLeftOfShock();

		void calculateFluidSpeedLeftOfShock();
		void calculateFluidSpeedLeftOfContact();

		void calculatePropagationSpeedOfContact();

		void calculateDensityLeftOfContact();

		void calculatePropagationSpeedOfShock();

		void calculatePositionOfExpansionTail();
		void calculatePositionOfExpansionHead();
		void calculatePositionOfContact();
		void calculatePositionOfShock();

		void calculateMach();
		//void calculateTemperature();
		void calculateDensity();
		void calculatePressure();

		void runSolution(std::vector<double> &x, double PL, double densityL, double PR, double densityR, double gamma, double t, double initialP_guess, double x0);
		void printMachTempDensityPressure(std::string a_file_name);

};



double Quasi1DShockTube::functionToSolve(double P){

	double alpha = (gamma+1)/(gamma-1);
	double beta = std::sqrt(2/(gamma*(gamma-1)));

	double term_1 = 2/(gamma-1)*AL/AR*(1-std::pow(pressure_R/pressure_L*P, (gamma-1)/(2*gamma)));
	double term_2 = beta*(P-1)/(std::sqrt(1+alpha*P));

	return term_1-term_2;
}

double Quasi1DShockTube::derivativeFunctionToSolve(double P){

	double alpha = (gamma+1)/(gamma-1);
	double beta = std::sqrt(2/(gamma*(gamma-1)));

	double term_1 = -(AL*pressure_R)/(AR*pressure_L*gamma)*std::pow(pressure_R/pressure_L*P, (gamma-1)/(2*gamma)-1);
	double term_2 = beta*(std::sqrt(1+alpha*P) -(P-1)/2/(std::sqrt(1 + alpha*P)))/(1+alpha*P);

	return term_1-term_2;
}

void Quasi1DShockTube::calculateSoundSpeedsLR(){

	AR = std::sqrt(gamma*pressure_R/(density_R));
	AL = std::sqrt(gamma*pressure_L/(density_L));

}

void Quasi1DShockTube::setInitialStates(double pR, double pL, double densityR, double densityL, double gamm, double t, double x0){
	pressure_R = pR;
	pressure_L = pL;
	density_R = densityR;
	density_L = densityL;
	gamma = gamm;
	time = t;
	initial_shock_location = x0;
}

void Quasi1DShockTube::newtonSolvePressureRatio(double initial){
	double result = initial;
	//std::cout << "Function To Solve value " <<  std::fabs(Quasi1DShockTube::functionToSolve(result)) << std::endl;
	while( std::fabs(Quasi1DShockTube::functionToSolve(result)) > TOL ){
		result = result - Quasi1DShockTube::functionToSolve(result)/Quasi1DShockTube::derivativeFunctionToSolve(result);
	}
	//std::cout << "Function To Solve value after iterating " <<  std::fabs(Quasi1DShockTube::functionToSolve(result)) << std::endl;
	pressure_ratio_across_shock = result;
	//std::cout << "P is: "<< result << std::endl;

}

void Quasi1DShockTube::calculateDensityLeftOfShock(){
	double alpha = (gamma+1)/(gamma-1);
	density_left_of_shock = density_R*(1+alpha*pressure_ratio_across_shock)/(alpha+pressure_ratio_across_shock);
}


void Quasi1DShockTube::calculatePressureLeftOfContact(){
	pressure_left_of_contact = pressure_left_of_shock;
}

void Quasi1DShockTube::calculatePressureLeftOfShock(){
	pressure_left_of_shock = pressure_R * pressure_ratio_across_shock;
	//std::cout << "Pressure Left of Shock " << pressure_left_of_shock << std::endl;
}


void Quasi1DShockTube::calculateDensityLeftOfContact(){
	density_left_of_contact = density_L*std::pow(pressure_left_of_contact/pressure_L, 1/gamma);
}

void Quasi1DShockTube::calculatePropagationSpeedOfContact(){
	contact_velocity = 2/(gamma-1)*AL*(1-std::pow(pressure_left_of_contact/pressure_L, (gamma-1)/(2*gamma)));
}

void Quasi1DShockTube::calculateFluidSpeedLeftOfShock(){
	fluid_velocity_left_of_shock = contact_velocity;
}

void Quasi1DShockTube::calculateFluidSpeedLeftOfContact(){
	fluid_velocity_left_of_contact = contact_velocity;
}

void Quasi1DShockTube::calculatePropagationSpeedOfShock(){
	shock_velocity = ((pressure_ratio_across_shock - 1)*std::pow(AR,2))/(gamma*fluid_velocity_left_of_shock);
	//std::cout << pressure_ratio_across_shock <<' ' <<AR <<' ' << gamma <<' ' << fluid_velocity_left_of_shock << std::endl;
}


void Quasi1DShockTube::calculatePositionOfExpansionTail(){
	tail_of_expansion_fan_location = initial_shock_location + (contact_velocity*(gamma+1)/2 - AL)*time;
}

void Quasi1DShockTube::calculatePositionOfExpansionHead(){
	tail_of_expansion_fan_location = initial_shock_location - AL*time;
}

void Quasi1DShockTube::calculatePositionOfContact(){
	contact_position = initial_shock_location + (contact_velocity)*time;
}

void Quasi1DShockTube::calculatePositionOfShock(){
	shock_position = initial_shock_location + (shock_velocity)*time;
}

void Quasi1DShockTube::calculatePressure(){

}

void Quasi1DShockTube::calculateMach(){


	for (unsigned int i = 0; i <X.size(); ++i)
	{
		if (X[i]<head_of_expansion_fan_location){
			Mach[i] = 0;
		}
		else if (X[i]>= head_of_expansion_fan_location & X[i]<tail_of_expansion_fan_location)
		{
			Mach[i] = -1.;
		}
		else if (X[i]>=tail_of_expansion_fan_location & X[i]<contact_position)
		{
			Mach[i] = 1.;
		}
		else if (X[i]>=contact_position & X[i]<shock_position)
		{
			Mach[i] = -1.;
		}
		else if (X[i]>shock_position)
		{
			Mach[i] = 0.;
		}
	}
	std::cout << "Finished Calculating Mach " << std::endl;
}

void Quasi1DShockTube::calculateDensity(){

}








































void Quasi1DShockTube::runSolution(std::vector<double>&x, double PL, double densityL, double PR, double densityR, double gamma, double time, double initialP_guess, double x0){

	X = x;
	Quasi1DShockTube::setInitialStates(PR, PL, densityR, densityL, gamma, time, x0);
	Quasi1DShockTube::calculateSoundSpeedsLR();
	Quasi1DShockTube::newtonSolvePressureRatio(initialP_guess);
	Quasi1DShockTube::calculatePressureLeftOfShock();
	Quasi1DShockTube::calculatePressureLeftOfContact();
	Quasi1DShockTube::calculatePropagationSpeedOfContact();
	Quasi1DShockTube::calculateFluidSpeedLeftOfShock();
	Quasi1DShockTube::calculateFluidSpeedLeftOfContact();
	Quasi1DShockTube::calculateDensityLeftOfContact();
	Quasi1DShockTube::calculatePropagationSpeedOfShock();
	Quasi1DShockTube::calculatePositionOfExpansionTail();
	Quasi1DShockTube::calculatePositionOfExpansionHead();
	Quasi1DShockTube::calculatePositionOfContact();
	Quasi1DShockTube::calculatePositionOfShock();



}


void Quasi1DShockTube::printMachTempDensityPressure(std::string a_file_name){
	//Quasi1DFlow::calculatePhysicalQuantities();

	std::ofstream a_file;
		a_file.open(a_file_name+ "_mach.csv", std::ios::out | std::ios::trunc);

			for (unsigned int i = 0; i< X.size()+1; ++i)
			{
				a_file <<  " " << X[i]  << ", " << std::to_string(Mach[i]) << std::endl;
			}
			a_file.close();

	std::ofstream a_file1;
		a_file1.open(a_file_name+ "_temp.csv", std::ios::out | std::ios::trunc);

			for (unsigned int i = 0; i< X.size()+1; ++i)
			{
				//a_file1 <<  " " << X[i]  << ", " << std::to_string(Temperature[i]) << std::endl;
			}
			a_file1.close();

		std::ofstream a_file2;
			a_file2.open(a_file_name+ "_density.csv", std::ios::out | std::ios::trunc);

				for (unsigned int i = 0; i< X.size()+1; ++i)
				{
					//a_file2 <<  " " << X[i]  << ", " << std::to_string(Density[i]) << std::endl;
				}
				a_file2.close();

		std::ofstream a_file3;
			a_file3.open(a_file_name+ "_pressure.csv", std::ios::out | std::ios::trunc);

				for (unsigned int i = 0; i< X.size()+1; ++i)
				{
					//a_file3 <<  " " << X[i]  << ", " << std::to_string(Pressure[i]) << std::endl;
				}
				a_file3.close();

}


#endif /* QUASI1DSHOCKTUBE_H_ */
