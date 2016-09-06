//
//  main.cpp
//  Final_Problem1
//
//  Created by Huachen Qin on 12/17/15.
//  Copyright © 2015 Huachen Qin. All rights reserved.
//

#include "SWING.h"

using namespace std;

int main(int argc, const char * argv[])
{
    vector<double> exercises = {1.0/3, 2.0/3, 1.0};
    cout<<"Three possible exercise time are: "<<exercises[0]<<", "<<exercises[1]<<", "<<exercises[2]<<endl;
    
    SwingOption option(50, 1, 3, 0.07, 0.2, exercises);    // create swing option class
    
    cout<<"-------------------------"<<endl;
    cout<<"swing option price v(0,50,2): "<<endl;
    cout<<"When exercise dates are: "<<exercises[0]<<" and "<<exercises[1]<<": "<<option.pricing(0, 50, 2, exercises[0], exercises[1])<<endl;
    cout<<"When exercise dates are: "<<exercises[0]<<" and "<<exercises[2]<<": "<<option.pricing(0, 50, 2, exercises[0], exercises[2])<<endl;
    cout<<"When exercise dates are: "<<exercises[1]<<" and "<<exercises[2]<<": "<<option.pricing(0, 50, 2, exercises[1], exercises[2])<<endl;
    cout<<"Combining all situations, the swing option price is: "<<option.general_pricing(0, 50, 2)<<endl;
    cout<<"-------------------------"<<endl;
    
    cout<<"Now print detail data:"<<endl;
    cout<<"When exercise dates are: "<<exercises[0]<<" and "<<exercises[1]<<": "<<endl;
    option.print_grid(0, 50, 2, exercises[0], exercises[1]);
    cout<<"-------------------------"<<endl;
    cout<<"When exercise dates are: "<<exercises[0]<<" and "<<exercises[2]<<": "<<endl;
    option.print_grid(0, 50, 2, exercises[0], exercises[2]);
}
