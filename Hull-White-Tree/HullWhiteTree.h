// HW11 problem1
// Huachen Qin

#ifndef HullWhiteTree_h
#define HullWhiteTree_h

#include<iostream>
#include<cmath>
#include<vector>

using namespace std;

class HWTree
{
private:
    double alpha;
    double sigma;        // volatility
    int N;               // total time step, in this case, N = 360
    double T;            // total time length, in this case, T = 30
    double delta_x;      // delta x
    double delta_t;      // delta_t = 1 / 12
    int jmax;
    vector<vector<double>> tree;      // the interest rate tree before calibration
    vector<vector<double>> r_tudi;    // the process of r_tudi
    vector<double> gamma;             // gamma
    vector<double> mt;                // m_delta
    vector<double> bond_p;            // true interest rate over 30 years
    vector<vector<double>> rt;                // the interest rate after calibration
    
public:
    HWTree();                         // constructor
    ~HWTree();                        // deconstructor
    void get_r_tudi();
    void get_gamma();
    void get_mt();
    void get_rt();
    vector<double> normal_prob(int step, int j);   // normal tree branch probability
    vector<double> top_prob(int j);                // top tree branch probability
    vector<double> bot_prob(int j);                // bottom tree branch probability
    void print_gamma();
    void print_tree();
};


#endif /* HullWhiteTree_h */
