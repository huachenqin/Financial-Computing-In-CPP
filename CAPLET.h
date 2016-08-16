//
//  CAPLET.h
//  Final_Problem2
//
//  Created by Huachen Qin on 12/17/15.
//  Copyright © 2015 Huachen Qin. All rights reserved.
//

#ifndef CAPLET_h
#define CAPLET_h

#include <vector>
#include <cmath>

using namespace std;

class Caplet
{
public:
    Caplet(vector<double> zeroRate, double Speed_, double Volatility_, long Steps_, double Time_, double Strike_, double Notional_, double delta);   // constructor
    ~Caplet();                    // destructor
    double pricing();             // return the price of the caplet
    void print_payoff();          // print the cumulative payoff tree
    
private:
    double Speed;                // alpha
    double Volatility;           // volatility sigma
    long Steps;                  // total steps  N = 5 / 0.25 = 20
    double Time;                 // T = 5 in this problem
    double Strike;               // K = 0.01
    double Notional;             // Notional Principal N = 1000000
    double delta;                // delta = 0.25
    
    
    vector<vector<double> > TheTree;      // Hull White Tree, each node contains instantaneous rate between (t,t+delta)
    vector<vector<double>> BondTree;      // BondTree, each node contains price of the Bond price between (t,t+delta) that pays off 1 at time t+delta
    vector<vector<double>> LIBORTree;     // LIBOR tree, each node contains LIBOR rate calculated from bond price between (t,t+delta)
    vector<vector<double> > r;            // HUll White model Gamma
    vector<vector<double> > Q;            // Hull White model bond coefficient Q
    vector<vector<double>> cumulative_payoff;   // cumulative payoff tree

};
#endif /* CAPLET_h */
