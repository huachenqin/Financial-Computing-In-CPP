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
    double sigma;
    int N;
    double T;
    double delta_x;
    double delta_t;
    int jmax;
    vector<vector<double>> tree;
    vector<vector<double>> r_tudi;
    vector<double> gamma;
    vector<double> mt;
    vector<double> bond_p;
    vector<double> rt;
    
public:
    HWTree();
    ~HWTree();
    HWTree(double T, int N);
    void get_r_tudi();
    void get_gamma();
    void get_mt();
    void get_rt();
    vector<double> normal_prob(int step, int j);
    vector<double> top_prob(int j);
    vector<double> bot_prob(int j);
    void print_tree();
};


#endif /* HullWhiteTree_h */
