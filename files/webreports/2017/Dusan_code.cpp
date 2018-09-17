//
//  main.cpp
//  NewEuler
//
//  Software tool implementing Euler's method for simulating meteoric properties
//
//  Created by Marcell Dorian Kovacs on 2017. 08. 04.
//
//  Project: Modeling Meteors
//
//  Summer School of Science
//  Pozega, Croatia
//
//  Project: Modeling Meteors
//
//  Young scientists:
//      Chloé Udressy
//      Marcell Dorian Kovacs
//      Nensi Komljenovic
//
//  Project leader:
//      Dusan Pavlovic
//
//  Copyright © 2017. Interplanetary Dušanoids. All rights reserved.
//
//
//

#include <iostream>     // For input and output on the console
#include <fstream>      // For file I/O
#include <cmath>        // For mathematical functions

using namespace std;

ofstream MeteorOutput("Meteor2.txt");// Output file for simulation


const double ZeroDensity = 0.122;    // Density of surface air [kg/m^3]
const double MolarMass = 0.029;      // Molar mass of ideal gas athmosphere [kg]
const double GasConstant = 8.314;    // Universal Gas Constant [kg*m^2/s^2/K/mol]
const double Gravity = 9.81;         // Gravitational acceleration on Earth [m/s^2]
const double Temperature = 200.0;    // Temperature of the gas [Kelvin]
const double Lambda = 1.0;           // Heat transfer coefficient (Fraction of energy)
const double Zet = 6e6;              // Latent heat coefficient (Energy lost during heating)
const double Gamma = 1.0;            // Drag coefficient (Fraction of momentum)
const double Shape = 1.21;           // Shape coefficient (Geometrical property of a sphere)

const double h = 0.0001;             // dt [s], timestep, time "instant"

double AthmosphericDensity(double height) {
    
    // Value of athmospheric density at a certain height
    // Simplifications: we model the athmosphere as an isothermic, ideal gas
    // This formula is also called the simple athmospheric model or bariometric formula
    
    return ZeroDensity * exp( - (MolarMass * Gravity * height) / (GasConstant * Temperature));
                             
}
double tau(double velocity) {
    
    // Tau parameter of a meteor
    // Describes the amount of kinetic energy is radiated away in form of (visible) light
    // Connects the kinetic energy to light intensity coming from the body
    
    const double v = velocity * 0.001; // Velocity must be expressed in [km/s]
    
    // Coefficients
    
    const double  c0 =  -0.001365811756924721;
    const double  c1 =  4.03316302906481e-4;
    const double  c2 = -3.837923319944585e-5;
    const double  c3 = 1.605242334632151e-6;
    const double  c4 = -2.995326211721562e-8;
    const double  c5 = 2.405216587475041e-10;
    const double  c6 = -4.425555198678479e-13;
    const double  d1 = -0.1939106756278702;
    const double  d2 = 0.02036594095461284;
    const double  d3 = -0.00130047417574090;
    const double  d4 = 4.580467150003708e-5;
    const double  d5 = -7.505453593442128e-7;
    const double  d6 = 4.80068343434027e-9;
    
    const double emi = 7.668;
    
    const double zetav2 = ( c0 + c1 * pow(v, 1.0)+c2 * pow(v, 2.0)+ c3 * pow(v, 3.0) + c4 * pow(v, 4.0) + c5 * pow(v, 5.0) + c6 * pow(v,6.0) ) / ( 1 + d1 * pow(v,1.0) + d2 * pow(v, 2.0) + d3 * pow(v, 3.0) + d4 * pow(v,4.0 ) + d5 * pow(v, 5.0) + d6 * pow(v, 6.0) );
    
    // parameter "emi" is the mean energy of excitation divided by a "mu" parameter determined from spectral analysis
    // parameter zetav2 is the coefficient of excitation divided by the velocity squared
    // our program and project approximated the value as a rational polynomial, which is only applicable in the range of our velocitites
    
    
    return 2 * emi * zetav2;
}


// Functions of equations of motion (originally a differential equation system known as Single Body Theory)

double velo(double velocity, double mass, double height, double density, double z){
    
    // A single body's velocity at a certain instant with the relevant parameters
    // Coming from transfer of momentum
    
    return - Gamma * Shape * pow( density , -2.0 / 3.0 ) * AthmosphericDensity(  height ) * pow( mass , -1.0 / 3.0 ) * pow( velocity , 2.0);
    
}
double mass(double velocity, double mass, double height, double density, double z){
    
    // A single body's mass at a certain instant with the relevant parameters
    // Coming from transfer of energy
    
    return -( (Lambda * Shape) / (2.0 * Zet) ) * pow(density, -2.0 / 3.0) * AthmosphericDensity(height) * pow(mass, 2.0 / 3.0) * pow(velocity , 3.0);
}
double height(double velocity, double mass, double height, double density, double z){
    
    // A single body's height at a certain instant with the relevant parameters
    // Coming from geometry and velocity (the time integral of the "y" axis projection of velocity)
    return - velocity * cos( z );
}
double intensity(double velocity, double mass, double height, double density, double z) {
    
    // Intensity of light coming from a meteoroid at a certain instant
    // Masschange is the mass() function for the arguments, but I could not include a function call of mass() here
    
    const double masschange = - ( (Lambda * Shape) / (2.0 * Zet) ) * pow(density, -2.0 / 3.0) * AthmosphericDensity(height) * pow(mass, 2.0 / 3.0) * pow(velocity , 3.0);
    
    return - ( 1.0 / 2.0 ) * tau(velocity) * pow( velocity , 2.0 ) * masschange ;
}


void Euler(double InitVelo, double InitMass, double InitHeight, double InitDensity, double Initz) {
    
    // Implementation of Euler's method of approximation to solve the differential equation system
    
    double Velo = InitVelo;
    double Mass = InitMass;
    double Height = InitHeight;
    double Density = InitDensity;
    double z = Initz;
    double Intensity = intensity(Velo, Mass, Height, Density, z);
    
    double MaxIntensity = 0;
    double MaxHeight = 0;
    
    int i = 0;
    
    while (true) {
        
        //  In our simplified model, the lifetime of the meteor ends, when it mass decreases under a certain point or we reach a certain point in calculating intensity
        
        if ( Mass < 0.001*InitMass || Intensity < MaxIntensity / 2)
           break;
        
        Velo += h * velo(Velo, Mass, Height, Density, z);
        
        Mass += h * mass(Velo, Mass, Height, Density, z);
        
        Height += h * height(Velo, Mass, Height, Density, z);
        
        Intensity = intensity(Velo, Mass, Height, Density, z);
        
        // Selecting the maximum of intensity and the height where such maximum occurs (MaxHeight)
        
        if (Intensity > MaxIntensity) {
            
            MaxIntensity = Intensity;
            MaxHeight = Height;
            
        }
        
        ++i;
        
        // Commands for outputting the actual light curves
        
        // MeteorOutput <<  i * h  << "       " << Velo << "       " << Mass << "       " << Height << "       "<< Intensity << endl;
        
        // cout <<  i * h  << " " << Velo << " " << Mass << " " << Height << " "<< Intensity << " " << MaxIntensity <<  endl;
    }
    
    // Commands for outputting only the maximums of the light curves
    
    MeteorOutput/* << InitVelo << " " <<  InitMass << " " << InitDensity << " " << Initz << " " */<< MaxIntensity << " " << MaxHeight << endl;
    cout/* << InitVelo << " " <<  InitMass << " " << InitDensity << " " << Initz << " " */<< MaxIntensity << " " << MaxHeight << endl;

    return;
    
    
}

int main(int argc, const char * argv[]) {
    
    // A little tribute to Hello World! but in physics style

    cout << "Hello Universe!" << endl;
    
    // Simulation of light curve maximums and height maximums with the following initial parameters for a meteor
    
    double InitCosZ, InitMass, InitVelo, InitDensity;
    
    const double InitHeight = 2e5; // Initial height is a paramater for the Euler, but is constant over the simulation
    
    int counter = 0; // Tracking the number of meteors calculated
    
    // Looping over several meteors with various initial conditions, saving the set of initials and the output of Euler's method
    
     
    for(InitCosZ = 1.0; InitCosZ >= 0.173648; InitCosZ /= pow(10.0 , 0.08))
        
        for(InitMass = pow(10.0, -8.0); InitMass <= pow(10.0, -1.5); InitMass *= pow(10.0, 0.08))
            
            for(InitVelo = 11200.0; InitVelo <= 72800.0; InitVelo *= pow(10.0, 0.032))
        
                for(InitDensity = 300.0; InitDensity <= 3100.0; InitDensity *= pow(10.0, 0.1)) {
                    
                    ++counter;
                    
                        // Printing the initial conditions to the file
                    
                    cout << counter << " " << InitVelo << " " << InitMass << " " << InitDensity << " " << InitCosZ << " ";
                        
                    MeteorOutput << counter << ". " << InitVelo << " " << InitMass << " " << InitDensity << " " << InitCosZ << " ";
                    
                        // Printing the results of Euler's method
                    
                    Euler(InitVelo, InitMass, InitHeight, InitDensity, acos(InitCosZ));
                    
                }
    
    cout << " A total of " << counter << " meteoroids were simulated." << endl;
    
    // Return gracefully
    
    return 0;
}


