//
//  SWING.h
//  Final_Problem1
//
//  Created by Huachen Qin on 12/17/15.
//  Copyright © 2015 Huachen Qin. All rights reserved.
//

#ifndef SWING_h
#define SWING_h

#include "matrix.h"
#include <vector>
#include <assert.h>
#include <cmath>
#include <iostream>

using namespace std;

class SwingOption
{
public:
    SwingOption(double strike_, double expiry_, int num_of_exercises_, double r_, double sigma_, vector<double>& exercises_); // constructor
    ~SwingOption();   // destructor
    
    double GetExpiry() const;  // get time to maturity
    int GetNumExpiry() const;  // get number of possible exercise dates
    
    double pricing(double T, double Spot, int k, double exercise_time_1, double exercise_time_2);  // pricing option given the exercise date
    double general_pricing(double T, double Spot, int k);  // pricing option considering all possible combinations of exercise dates
    
    vector<double> coeff_ABC(double m, double n, double delta_t,double delta_s);  // the A_n_m, B_n_m, C_n_m coefficients in Crank Nicolson Method
    void print_grid(double time, double Spot, int k_, double exercise_time1, double exercise_time2);      // print grid data
    
    friend double Call_PayOff(double spot, double strike); // payoff calculator of Call option
    friend double Put_PayOff(double spot, double strike);  // payoff calculator of Put option
    
private:
    double strike;        // stike price K
    double expiry;        // time to maturity T
    int num_of_exercises;   // number of exercise opportunities k
    double r;             // interest rate r
    double sigma;         // volatility sigma
    vector<double> exercises;  // vector containing possible exercise time
};

#endif /* SWING_h */
