//
//  main.cpp
//  Final_Problem2
//
//  Created by Huachen Qin on 12/17/15.
//  Copyright © 2015 Huachen Qin. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <vector>
#include "CAPLET.h"
#include "Spline.h"
using namespace std;

int main()
{
    //construct yield curve
    vector<double> term = { 0.0795, 0.2521, 0.5014, 0.9425, 1.9836, 3.0247, 4.9863, 6.9863, 9.8630, 29.877 };//Time
    vector<double> price = { 100, 100, 99.9600, 99.7700, 100.0003, 100.2659, 100.0688, 100.0440, 99.7240, 100.1336 };//Price
    vector<double> coupon = { 0, 0, 0, 0, 0.625, 0.875, 1.375, 1.750, 2.00, 2.875 };//Coupon
    vector<double> time;
    vector<double> rate;
    
    for (int i = 0; i < 4; ++i)
    {
        time.push_back(term[i]);
        rate.push_back(log(100 / price[i]) / term[i]);
    }
    vector<SplineSet> sp = spline(time, rate);
    
    for (int i = 4; i < static_cast<int>(term.size()); ++i)
    {
        vector<double> tt;
        double temp = term[i];
        int count = 0;
        while (temp>0.5) {
            temp -= 0.5;
            ++count;
        }
        for (int j = 0; j < count; ++j) {
            tt.push_back(temp);
            temp += 0.5;
        }
        tt.push_back(temp);
        double s = 0;
        for (int j = 0; j < static_cast<int>(tt.size()) - 1; ++j) {
            s += exp(-findRate(tt[j], sp)*tt[j])*coupon[i] / 2;
        }
        time.push_back(term[i]);
        rate.push_back(-(log((price[i] - s) / (coupon[i] / 2 + 100))) / term[i]);
        sp = spline(time, rate);//update curve
    }
    
    
    // set up caplet class parameters
    double sigma = 0.2, a = 1, Time =5.0, Strike = 0.01, Notional = pow(10,6), delta = 0.25;
    long Steps = 20;
    double dt = Time / Steps;
    vector<double> zeroRate(Steps+1);
    for (long i = 1; i <= Steps + 1; ++i){
        zeroRate.push_back(findRate(i*dt, sp));
    }
    
    // create caplet class
    Caplet caplet(zeroRate,a,sigma,Steps,Time,Strike,Notional,delta);
    
    cout<<"The price of the caplet is: "<<caplet.pricing()<<endl;
    
    cout<<"----------------------------------"<<endl;
    
    cout<<"Now print out detail data."<<endl;
    cout<<"In this question, jmax=1, hence the tree is at most have 3 branches. Now print out the cumulative tree as follows:"<<endl;
    caplet.print_payoff();
    return 0;
}