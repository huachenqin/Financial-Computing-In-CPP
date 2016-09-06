//
//  Spline.cpp
//  Final_Problem2
//
//  Created by Huachen Qin on 12/17/15.
//  Copyright © 2015 Huachen Qin. All rights reserved.
//

#include <stdio.h>

#include "Spline.h"
using namespace std;

vector<SplineSet> spline(vector<double> &x, vector<double> &y)
{
    int n = int(x.size()) - 1;
    vector<double> a;
    a.insert(a.begin(), y.begin(), y.end());
    vector<double> b(n);
    vector<double> d(n);
    vector<double> h(n);
    vector<double> alpha(n);
    alpha[0] = 0;
    
    for (int i = 0; i < n; ++i)
        h[i] = x[i + 1] - x[i];
    
    for (int i = 1; i < n; ++i)
        alpha[i] = (3 * (a[i + 1] - a[i]) / h[i] - 3 * (a[i] - a[i - 1]) / h[i - 1]);
    
    vector<double> c(n + 1);
    vector<double> l(n + 1);
    vector<double> mu(n + 1);
    vector<double> z(n + 1);
    l[0] = 1;
    mu[0] = 0;
    z[0] = 0;
    
    for (int i = 1; i < n; ++i)
    {
        l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }
    
    l[n] = 1;
    z[n] = 0;
    c[n] = 0;
    
    for (int j = n - 1; j >= 0; --j)
    {
        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
        d[j] = (c[j + 1] - c[j]) / 3 / h[j];
    }
    
    vector<SplineSet> output_set(n);
    for (int i = 0; i < n; ++i)
    {
        output_set[i].a = a[i];
        output_set[i].b = b[i];
        output_set[i].c = c[i];
        output_set[i].d = d[i];
        output_set[i].x = x[i];
    }
    return output_set;
}

double findRate(double t, vector<SplineSet> sp)
{
    double rate = 0;
    if (t <= sp[0].x)//extrapolation
    {
        return sp[0].a + sp[0].b*(t - sp[0].x) + sp[0].c*(t - sp[0].x)*(t - sp[0].x) + sp[0].d*(t - sp[0].x)*(t - sp[0].x)*(t - sp[0].x);
    } else if (t >= sp[sp.size() - 1].x)//extrapolation
    {
        int i = int(sp.size()) - 1;
        return sp[i].a + sp[i].b*(t - sp[i].x) + sp[i].c*(t - sp[i].x)*(t - sp[i].x) + sp[i].d*(t - sp[i].x)*(t - sp[i].x)*(t - sp[i].x);
    } else {
        for (unsigned int i = 0; i < sp.size(); i++)
        {
            if (sp[i].x >= t)//interpolation
            {
                --i;
                rate = sp[i].a + sp[i].b*(t - sp[i].x) + sp[i].c*(t - sp[i].x)*(t - sp[i].x) + sp[i].d*(t - sp[i].x)*(t - sp[i].x)*(t - sp[i].x);
                break;
            }
        }
    }
    return rate;
}