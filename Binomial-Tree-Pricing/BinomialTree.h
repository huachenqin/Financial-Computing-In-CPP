#ifndef BinomialTree_h
#define BinomialTree_h

#include<iostream>
#include<cmath>
#include<vector>
#include<assert.h>

using namespace std;

class CRRBTree
{
private:
    double option_price;          // call option price
    double S;                     // spot stock price
    double K;                     // strike price
    double r;                     // interest rate
    double sigma;                 // volatility
    double T;                     // time to maturity
    int N;                        // number of steps in binomial tree
    double delta_t;               // the time length of each step
    double discount_factor;       // discount factor for each step
    double u, pu;                 // upfactor/downfactor and corresponding probability
    double d, pd;
    vector<vector<double>> stockprice;    // stockprice vector
    vector<vector<double>> optionvalue;   // optionvalue vector

public:
    CRRBTree();
    CRRBTree(double S, double K, double T, double r, double sigma, int N);
    ~CRRBTree();
    double optionpricing();      // pricing without pruning
    double pruning_1();          // pruning method one
    double pruning_2();          // pruning method two
    double pruning_3();          // pruning method three
    double ratio_analysis(int N_);
    double error_analysis(int N_);      // error analysis
    double Richardson(int N_);          // Richardson Extropolation
};

class TianBTree
{
private:
    double option_price;      // option value
    double S;                 // spot stock price
    double K;                 // strike price
    double r;                // interest rate
    double sigma;            // volatility
    double T;                 // time to maturity
    int N;                    // number of total steps
    double delta_t;           // time length of each step
    double discount_factor;       // discount factor for each step
    double u, pu;               // upfactor/downfactor and corresponding probability
    double d, pd;
    vector<vector<double>> stockprice;
    vector<vector<double>> optionvalue;

public:
    TianBTree();
    TianBTree(double S, double K, double T, double r, double sigma, int N);
    ~TianBTree();
    double optionpricing();
    double pruning_1();
    double pruning_2();
    double pruning_3();
    double ratio_analysis(int N_);
    double error_analysis(int N_);
    double Richardson(int N_);
};

class TriTree
{
private:
    double option_price;   // option price
    double S;              // spot stock price
    double K;              // strike price
    double r;              // interest rate
    double sigma;          // volatility
    double T;              // time to maturity
    int N;                 // number of steps in trinomial tree
    double delta_t;        // time length of each step
    double discount_factor;    // discount factor for each step
    double u, pu;          // upfactor, downfactor, middle factor and corresponding factor
    double m, pm;
    double d, pd;
    vector<vector<double>> stockprice;
    vector<vector<double>> optionvalue;
    
public:
    TriTree();
    TriTree(double S, double K, double T, double r, double sigma, int N);
    ~TriTree();
    double optionpricing();
    double pruning_1();
    double pruning_2();
    double pruning_3();
    double ratio_analysis(int N_);
    double error_analysis(int N_);
    double Richardson(int N_);

};

double BlackScholesCall(double S, double K, double sigma, double tau, double r); // define the function _ Call
double auxillary_function(double x);                         // define the auxiliary function to calculate Φ(x)


#endif
