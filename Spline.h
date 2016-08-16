//
//  Spline.h
//  Final_Problem2
//
//  Created by Huachen Qin on 12/17/15.
//  Copyright © 2015 Huachen Qin. All rights reserved.
//

#ifndef Spline_h
#define Spline_h

#include <vector>

struct SplineSet{
    double a;
    double b;
    double c;
    double d;
    double x;
};

std::vector<SplineSet> spline(std::vector<double>&, std::vector<double>&);
double findRate(double, std::vector<SplineSet>);


#endif /* Spline_h */
